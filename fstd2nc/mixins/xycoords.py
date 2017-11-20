###############################################################################
# Copyright 2017 - Climate Research Division
#                  Environment and Climate Change Canada
#
# This file is part of the "fstd2nc" package.
#
# "fstd2nc" is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# "fstd2nc" is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with "fstd2nc".  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

from fstd2nc import _, info, warn, error
from fstd2nc.mixins import _Buffer_Base


# Modify gdll to handle supergrids.
# Simply loops over each subgrid, then re-stacks them together.
def gdll (gdid):
  from rpnpy.librmn.interp import ezgprm, ezget_subgridids, gdll
  import numpy as np
  grtyp = ezgprm(gdid)['grtyp']
  if grtyp != 'U':
    return gdll(gdid)
  lat = []
  lon = []
  for subgrid in ezget_subgridids(gdid):
    subll = gdll(subgrid)
    lat.append(subll['lat'])
    lon.append(subll['lon'])
  # Stack them together again, to match the shape of the variable.
  lat = np.concatenate(lat, axis=1)
  lon = np.concatenate(lon, axis=1)
  return {'lat':lat, 'lon':lon}

# Modify gdgaxes to handle supergrids.
# Simply loops over each subgrid, then re-stacks them together.
def gdgaxes (gdid):
  from rpnpy.librmn.interp import ezgprm, ezget_subgridids, gdgaxes
  import numpy as np
  grtyp = ezgprm(gdid)['grtyp']
  if grtyp != 'U':
    return gdgaxes(gdid)
  ax = []
  ay = []
  for subgrid in ezget_subgridids(gdid):
    subll = gdgaxes(subgrid)
    ax.append(subll['ax'])
    ay.append(subll['ay'])
  # Stack them together again, to match the shape of the variable.
  ax = np.concatenate(ax, axis=1)
  ay = np.concatenate(ay, axis=1)
  return {'ax':ax, 'ay':ay}



#################################################
# Mixin for handling lat/lon coordinates.

class _XYCoords (_Buffer_Base):
  # Special records that contain coordinate info.
  # We don't want to output these directly as variables, need to decode first.
  _xycoord_nomvars = ('^^','>>','^>')
  # These are extra copies of the lat/lon information.
  _xycoord_nomvars = _xycoord_nomvars + ('LA','LO')
  # Grids that can be read directly from '^^','>>' records, instead of going
  # through ezqkdef (in fact, may crash ezqkdef if you try decoding them).
  _direct_grids = ('X','Y','T','+')

  @classmethod
  def _cmdline_args (cls, parser):
    super(_XYCoords,cls)._cmdline_args(parser)
    parser.add_argument('--subgrid-axis', action='store_true', help=_('For data on supergrids, split the subgrids along a "subgrid" axis.  The default is to leave the subgrids stacked together as they are in the RPN file.'))

  def __init__ (self, *args, **kwargs):
    self._subgrid_axis = kwargs.pop('subgrid_axis',False)
    # Tell the decoder not to process horizontal records as variables.
    self._meta_records = self._meta_records + self._xycoord_nomvars
    super(_XYCoords,self).__init__(*args,**kwargs)
    # Variables must have an internally consistent horizontal grid.
    self._var_id = self._var_id + ('grtyp',)
    self._human_var_id = self._human_var_id + ('%(grtyp)s',)
    # Also, must have consistent igX records for a variable.
    if 'ig1' not in self._var_id:
      self._var_id = self._var_id + ('ig1','ig2','ig3','ig4')
      self._human_var_id = self._human_var_id + ('grid_%(ig1)s_%(ig2)s_%(ig3)s_%(ig4)s',)

  # Helper method - look up a coordinate record for the given variable.
  # Need this for manual lookup of 'X' grids, since ezqkdef doesn't support
  # them?
  def _find_coord (self, var, coordname):
    from fstd2nc.mixins import dtype_fst2numpy
    from rpnpy.librmn.fstd98 import fstlir
    header = fstlir (self._meta_funit, nomvar=coordname, ip1=var.atts['ig1'],
                           ip2=var.atts['ig2'], ip3=var.atts['ig3'], rank=3)
    if header is not None:
      # Override output dtype
      dtype = dtype_fst2numpy(header['datyp'],header['nbits'])
      header['d'] = header['d'].view(dtype)
      return header
    raise KeyError("Unable to find matching '%s' for '%s'"%(coordname,var.name))


  # Add horizontal coordinate info to the data stream.
  def _iter (self):
    from fstd2nc.mixins import _iter_type, _var_type, _modify_axes
    from collections import OrderedDict
    from rpnpy.librmn.interp import ezqkdef, EzscintError, ezget_nsubgrids
    import numpy as np

    # Scan through the data, and look for any use of horizontal coordinates.
    grids = OrderedDict()
    # Only output 1 copy of 1D coords (e.g. could have repetitions with
    # horizontal staggering.
    coords = set()
    for var in super(_XYCoords,self)._iter():
      # Don't touch derived variables.
      if not isinstance(var,_iter_type):
        yield var
        continue
      # Get grid parameters.
      ni = int(var.atts['ni'])
      nj = int(var.atts['nj'])
      grtyp = var.atts['grtyp']
      ig1 = int(var.atts['ig1'])
      ig2 = int(var.atts['ig2'])
      ig3 = int(var.atts['ig3'])
      ig4 = int(var.atts['ig4'])
      # Uniquely identify the grid for this variable.
      #
      # Use a looser identifier for timeseries data (ni/nj have different
      # meanings here (not grid-related), and could have multiple grtyp
      # values ('+','Y') that should share the same lat/lon info.
      if var.atts['typvar'].strip() == 'T':
        key = ('T',ig1,ig2)
      else:
        key = (grtyp,ni,nj,ig1,ig2,ig3,ig4)
      if grtyp in ('Y','+'): key = key[1:]
      # Check if we already defined this grid.
      if key not in grids:

        try:
          # Get basic information about this grid.
          # First, handle non-ezqkdef grids.
          if grtyp in self._direct_grids:
            lat = self._find_coord(var,'^^')['d'].squeeze(axis=2)
            lon = self._find_coord(var,'>>')['d'].squeeze(axis=2)
            ll = {'lat':lat, 'lon':lon}
          # Everything else should be handled by ezqkdef.
          else:
            gdid = ezqkdef (ni, nj, grtyp, ig1, ig2, ig3, ig4, self._meta_funit)
            ll = gdll(gdid)
        except (TypeError,EzscintError,KeyError):
          warn(_("Unable to get grid info for '%s'")%var.name)
          yield var
          continue


        # Find X/Y coordinates (if applicable).
        try:
          # Can't do this for direct grids (don't have a gdid defined).
          if grtyp in self._direct_grids: raise TypeError
          xycoords = gdgaxes(gdid)
          ax = xycoords['ax'].transpose()
          ay = xycoords['ay'].transpose()
          # Convert from degenerate 2D arrays to 1D arrays.
          ax = ax[0,:]
          ay = ay[:,0]
          xaxis = _var_type('x',{'axis':'X'},{'x':tuple(ax)},ax)
          yaxis = _var_type('y',{'axis':'Y'},{'y':tuple(ay)},ay)
        except (TypeError,EzscintError):
          # Can't get X/Y coords for this grid?
          xaxis = yaxis = None

        # Construct lat/lon fields.
        latarray = ll['lat'].transpose() # Switch from Fortran to C order.
        latatts = OrderedDict()
        latatts['long_name'] = 'latitude'
        latatts['standard_name'] = 'latitude'
        latatts['units'] = 'degrees_north'
        lonarray = ll['lon'].transpose() # Switch from Fortran to C order.
        lonatts = OrderedDict()
        lonatts['long_name'] = 'longitude'
        lonatts['standard_name'] = 'longitude'
        lonatts['units'] = 'degrees_east'

        # Case 1: lat/lon can be resolved into 1D Cartesian coordinates.
        # Calculate the mean lat/lon arrays in double precision.
        meanlat = np.mean(np.array(latarray,dtype=float),axis=1,keepdims=True)
        meanlon = np.mean(np.array(lonarray,dtype=float),axis=0,keepdims=True)
        if latarray.shape[1] > 1 and lonarray.shape[1] > 1 and np.allclose(latarray,meanlat) and np.allclose(lonarray,meanlon):
          # Reduce back to single precision for writing out.
          meanlat = np.array(meanlat,dtype=latarray.dtype).squeeze()
          meanlon = np.array(meanlon,dtype=lonarray.dtype).squeeze()
          # Ensure monotonicity of longitude field.
          # (gdll may sometimes wrap last longitude to zero).
          # Taken from old fstd_core.c code.
          if meanlon[-2] > meanlon[-3] and meanlon[-1] < meanlon[-2]:
            meanlon[-1] += 360.
          latarray = meanlat
          lonarray = meanlon
          lat = _var_type('lat',latatts,{'lat':tuple(latarray)},latarray)
          lon = _var_type('lon',lonatts,{'lon':tuple(lonarray)},lonarray)
          gridaxes = [('lat',tuple(latarray)),('lon',tuple(lonarray))]

        # Case 2: lat/lon are series of points.
        elif latarray.shape[0] == 1 and lonarray.shape[0] == 1 and ('i' in var.axes or 'station_id' in var.axes):
          latarray = latarray[0,:]
          lonarray = lonarray[0,:]
          # Special case for station data
          if 'station_id' in var.axes:
            gridaxes = [('station_id',var.axes['station_id'])]
          else:
            gridaxes = [('i',var.axes['i'])]
          lat = _var_type('lat',latatts,OrderedDict(gridaxes),latarray)
          lon = _var_type('lon',lonatts,OrderedDict(gridaxes),lonarray)

        # Case 3: General 2D lat/lon fields on X/Y coordinate system.
        elif xaxis is not None and yaxis is not None:
          gridaxes = [('y',tuple(yaxis.array)),('x',tuple(xaxis.array))]
          # Special case: have supergrid data, and the user wants to split it?
          if grtyp == 'U' and self._subgrid_axis:
            ngrids = ezget_nsubgrids(gdid)
            ny = len(yaxis.array)//ngrids
            yaxis.array = yaxis.array[:ny]
            yaxis.axes['y'] = tuple(yaxis.array)
            gridaxes = [('subgrid',tuple(range(ngrids))), ('y',tuple(yaxis.array)), ('x',tuple(xaxis.array))]
            latarray = latarray.reshape(ngrids,ny,-1)
            lonarray = lonarray.reshape(ngrids,ny,-1)
          if tuple(yaxis.axes.items()) not in coords:
            yield yaxis
            coords.add(tuple(yaxis.axes.items()))
          if tuple(xaxis.axes.items()) not in coords:
            yield xaxis
            coords.add(tuple(xaxis.axes.items()))
          lat = _var_type('lat',latatts,OrderedDict(gridaxes),latarray)
          lon = _var_type('lon',lonatts,OrderedDict(gridaxes),lonarray)

        # Case 4: General 2D lat/lon fields with no coordinate system.
        elif 'i' in var.axes and 'j' in var.axes:
          gridaxes = [('j',var.axes['j']),('i',var.axes['i'])]
          lat = _var_type('lat',latatts,OrderedDict(gridaxes),latarray)
          lon = _var_type('lon',lonatts,OrderedDict(gridaxes),lonarray)

        else:
          warn(_("Unhandled lat/lon coords for '%s'")%var.name)
          yield var
          continue

        # Sanity check on lat/lon - make sure we have something of the right size.
        if lat.array.shape == tuple(map(len,lat.axes.values())) and lon.array.shape == tuple(map(len,lon.axes.values())):
          yield lat
          yield lon
          grids[key] = gridaxes
        else:
          warn(_("Wrong shape of lat/lon for '%s'")%var.name)
          yield var
          continue

      gridaxes = grids[key]

      # Update the var's horizontal coordinates.
      if len(gridaxes) == 1:
        var.axes = _modify_axes(var.axes, i=gridaxes[0])
      elif len(gridaxes) == 2:
        var.axes = _modify_axes(var.axes, j=gridaxes[0], i=gridaxes[1])
      elif len(gridaxes) == 3:
        var.axes = _modify_axes(var.axes, k=gridaxes[0], j=gridaxes[1], i=gridaxes[2])
      else:
        warn(_("Unusual grid axes for '%s' - ignoring.")%var.name)

      # For 2D lat/lon, need to reference them as coordinates in order for
      # netCDF viewers to display the field properly.
      if 'lat' not in var.axes or 'lon' not in var.axes:
        var.atts['coordinates'] = 'lon lat'

      yield var

