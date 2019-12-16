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

from fstd2nc.stdout import _, info, warn, error
from fstd2nc.mixins import BufferBase


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

# Base class for grid mapping classes
class GridMap(object):
  # Mean radius of the Earth used in OGC CRS (coord. ref. system) descriptions
  # Note that Earth's mean radius = 6371000. m in librmn!
  _earth_radius = 6371229. 
  def __init__(self, grd):
#   grd is a dictionary returned by readGrid containing grid parameters 
    from collections import OrderedDict
    self._grd = grd
    self._name = ''
    self._atts = OrderedDict()
    self._xaxisatts = OrderedDict()
    self._yaxisatts = OrderedDict()
# Factory method that creates various types of grid mapping objects
  @classmethod
  def gen_gmap(cls, grd):
    import numpy as np
    grref = grd['grref'].upper() 
    if grref == 'E' :
      # Special case: 'E' grid is not actually rotated.
      if np.allclose(grd['lat0'],grd['rlat0']) and np.allclose(grd['lon0'],grd['rlon0']):
        return LatLon(grd)
      # Usual case: 'E' grid is rotated.
      return RotLatLon(grd)
    elif grref in ('A','B','G','L') :
      return LatLon(grd)
    elif grref in ['N','S'] :
      return PolarStereo(grd)
    else:
      raise ValueError('Grid mapping variable cannot be created for grids based on grid type %s!' \
                       % (grref))

class LatLon(GridMap):
  def __init__(self, *args, **kwargs):
    super(LatLon,self).__init__(*args,**kwargs)
    # Grid mapping variable name
    self._name = 'crs_latlon'
  def gen_gmapvar(self):
    from fstd2nc.mixins import _var_type
    import numpy as np
    self._atts['grid_mapping_name'] = 'latitude_longitude'
    self._atts['earth_radius'] = self._earth_radius
    # Grid mapping variable
    self.gmap = _var_type(self._name,self._atts,[],np.array(b""))
    return self.gmap     
  # Generate true latitudes and longitudes
  def gen_ll(self):    
    from collections import OrderedDict
    from rpnpy.librmn.all import gdll 
    from fstd2nc.mixins import _axis_type
    ll = gdll(self._grd['id'])
    self._lonarray = ll['lon'][:,0]
    self._lonatts = OrderedDict()
    self._lonatts['long_name'] = 'longitude'
    self._lonatts['standard_name'] = 'longitude'
    self._lonatts['units'] = 'degrees_east'
    self._latarray = ll['lat'][0,:]
    self._latatts = OrderedDict()
    self._latatts['long_name'] = 'latitude'
    self._latatts['standard_name'] = 'latitude'
    self._latatts['units'] = 'degrees_north'
    self.lon = _axis_type('lon',self._lonatts,self._lonarray)
    self.lat = _axis_type('lat',self._latatts,self._latarray)
    self.gridaxes = [self.lat,self.lon]
    # Ensure monotonicity of longitude field.
    # (gdll may sometimes wrap last longitude to zero).
    # Taken from old fstd_core.c code.
    if len(self._lonarray) >= 3 and self._lonarray[-2] > self._lonarray[-3] and self._lonarray[-1] < self._lonarray[-2]:
      self._lonarray[-1] += 360.
    return (self.gridaxes, self.lon, self.lat)
  # Generic gen_xyll function.
  def gen_xyll(self):
    gridaxes, lon, lat = self.gen_ll()
    return (None, None, gridaxes, lon, lat)

class RotLatLon(GridMap):
  def __init__(self, *args, **kwargs):
    from rpnpy.librmn.all import egrid_rll2ll, egrid_ll2rll
    super(RotLatLon,self).__init__(*args,**kwargs)
    # Grid mapping variable name
    self._name = 'rotated_pole'
    # Calculate geographical latitude and longitude of rotated grid's North Pole
    (__rlat_RNP, __rlon_RNP) = (90., 0.)
    (self._grid_north_pole_latitude,self._grid_north_pole_longitude) = \
       egrid_rll2ll(self._grd['xlat1'], self._grd['xlon1'], 
       self._grd['xlat2'], self._grd['xlon2'], __rlat_RNP, __rlon_RNP)
    # Calculate rotated-grid latitude and longitude of the geographical North Pole
    (__lat_NP, __lon_NP) = (90., 0.)
    (self._north_pole_grid_latitude,self._north_pole_grid_longitude) = \
       egrid_ll2rll(self._grd['xlat1'], self._grd['xlon1'], 
       self._grd['xlat2'], self._grd['xlon2'], __lat_NP, __lon_NP)
    self._ax = self._grd['ax'][:,0] 
    # Offset applied to bring rotated longitude in range [-180,180]
    # Done to avoid crossing the dateline and representation problems 
    # in some software (e.g. IDV, Panoply, Iris)
    if self._ax.max() > 180. :
      self._ax = self._ax - self._north_pole_grid_longitude 
    self._ay = self._grd['ay'][0,:]
  def gen_gmapvar(self):
    from fstd2nc.mixins import _var_type
    import numpy as np
    self._atts['grid_mapping_name'] = 'rotated_latitude_longitude'
    self._atts['earth_radius'] = self._earth_radius
    self._atts['grid_north_pole_latitude'] = self._grid_north_pole_latitude
    self._atts['grid_north_pole_longitude'] = self._grid_north_pole_longitude
#    self._atts['north_pole_grid_longitude'] = self._north_pole_grid_longitude 
#   Set the optional grid mapping parameter 'north_pole_grid_longitude' to 0 to avoid 
#   some problems, such as the conversion from netcdf to grib performed by some tools
#    self._atts['north_pole_grid_longitude'] = 0.
#    self._atts['longitude_of_prime_meridian'] = 0.
    # Grid mapping variable
    self.gmap = _var_type(self._name,self._atts,[],np.array(b""))
    return self.gmap
  # Generate latitudes and longitudes in rotated pole grid 
  # and true latitudes and longitudes
  def gen_xyll(self):
    from collections import OrderedDict
    from rpnpy.librmn.all import gdll
    from fstd2nc.mixins import _var_type, _axis_type
    self._xaxisatts['long_name'] = 'longitude in rotated pole grid'
    self._xaxisatts['standard_name'] = 'grid_longitude'
    self._xaxisatts['units'] = 'degrees'
    self._xaxisatts['axis'] = 'X'
    self._yaxisatts['long_name'] = 'latitude in rotated pole grid'
    self._yaxisatts['standard_name'] = 'grid_latitude'
    self._yaxisatts['units'] = 'degrees'
    self._yaxisatts['axis'] = 'Y'
    self.xaxis = _axis_type('rlon',self._xaxisatts,self._ax)
    self.yaxis = _axis_type('rlat',self._yaxisatts,self._ay)
    self.gridaxes = [self.yaxis,self.xaxis]
    ll = gdll(self._grd['id'])
    self._lonarray = ll['lon'].transpose() # Switch from Fortran to C order.
    self._lonatts = OrderedDict()
    self._lonatts['long_name'] = 'longitude'
    self._lonatts['standard_name'] = 'longitude'
    self._lonatts['units'] = 'degrees_east'
    self._latarray = ll['lat'].transpose() # Switch from Fortran to C order.
    self._latatts = OrderedDict()
    self._latatts['long_name'] = 'latitude'
    self._latatts['standard_name'] = 'latitude'
    self._latatts['units'] = 'degrees_north'
    self.lon = _var_type('lon',self._lonatts,self.gridaxes,self._lonarray)
    self.lat = _var_type('lat',self._latatts,self.gridaxes,self._latarray)
    return (self.xaxis, self.yaxis, self.gridaxes, self.lon, self.lat)


class PolarStereo(GridMap):
  # The standard parallel is fixed at 60 N as this is
  # the only standard parallel used in the RPN routines
  def __init__(self, *args, **kwargs):
    super(PolarStereo,self).__init__(*args,**kwargs)
    # Grid mapping variable name
    self._std_parallel = 60.
    self._name = 'polar_stereo'
    # Grid resolution (spacing) in projection plane (at standard parallel)
    self._res = self._grd['d60']    # metres
    self.xaxis = None
    self.yaxis = None
  @staticmethod
  def map_scale_factor(std_parallel_deg):
    # Calculate map scale factor from latitude of standard parallel
    # with formula found at:
    # https://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
    from math import fabs, sin, radians
    abs_sin = fabs(sin(radians(std_parallel_deg)))
    # value returned for standard parallel at 60 deg: 0.933012701892
    # value retrieved using mscale function (rmnlib): 0.9330124
    return (1. + abs_sin)/2. 
  def gen_gmapvar(self):
    import numpy as np
    from fstd2nc.mixins import _var_type
    self._atts['grid_mapping_name'] = 'polar_stereographic'
    self._atts['earth_radius'] = self._earth_radius
    if self._grd['north']:
      self._atts['latitude_of_projection_origin'] = 90.
    else:
      self._atts['latitude_of_projection_origin'] = -90.
    # Set central meridian so that easting and northing directions match those returned by RPN routines
    self._atts['straight_vertical_longitude_from_pole'] = -(self._grd['dgrw'] + 90.)
    # Set either 'standard_parallel' or 'scale_factor_at_projection_origin', not both!
#    self._atts['standard_parallel'] = 60.
    self._atts['scale_factor_at_projection_origin']= self.map_scale_factor(self._std_parallel)
    self._atts['resolution_at_standard_parallel'] = self._res
    self._gen_xyll()
    self._atts['false_easting'] = self._false_easting
    self._atts['false_northing'] = self._false_northing
    # Grid mapping variable
    self.gmap = _var_type(self._name,self._atts,[],np.array(b""))
    return self.gmap
  # Generate projection coordinates
  def _gen_xyll(self):  
    from collections import OrderedDict 
    import numpy as np
    from rpnpy.librmn.all import gdll, gdxyfll
    from fstd2nc.mixins import _var_type, _axis_type
    ll = gdll(self._grd['id'])
    self._lonarray = ll['lon'].transpose() # Switch from Fortran to C order.
    self._lonatts = OrderedDict()
    self._lonatts['long_name'] = 'longitude'
    self._lonatts['standard_name'] = 'longitude'
    self._lonatts['units'] = 'degrees_east'
    self._latarray = ll['lat'].transpose() # Switch from Fortran to C order.
    self._latatts = OrderedDict()
    self._latatts['long_name'] = 'latitude'
    self._latatts['standard_name'] = 'latitude'
    self._latatts['units'] = 'degrees_north'
    xy = gdxyfll(self._grd['id'],ll['lat'],ll['lon'])
    # Scale grid coordinates back to actual coordinates in projection plane   
    self._ax = ( np.rint(xy['x'][:,0]) - 1) * self._res   # metres
    self._ay = ( np.rint(xy['y'][0,:]) - 1) * self._res
    # Determine the false easting and northing from 
    # the coordinates of the pole and of grid point (1,1)
    _pole_pi = self._grd['pi']
    _pole_pj = self._grd['pj']
    self._false_easting =  self._ax[int(round(_pole_pi)) - 1] - self._ax[0]
    self._false_northing = self._ay[int(round(_pole_pj)) - 1] - self._ay[0]
    self._xaxisatts['long_name'] = 'x-coordinate of polar-stereographic projection'
    self._xaxisatts['standard_name'] = 'projection_x_coordinate'
    self._xaxisatts['units'] = 'm'
    self._xaxisatts['axis'] = 'X'
    self._yaxisatts['long_name'] = 'y-coordinate of polar-stereographic projection'
    self._yaxisatts['standard_name'] = 'projection_y_coordinate'
    self._yaxisatts['units'] = 'm'
    self._yaxisatts['axis'] = 'Y'
    self.xaxis = _axis_type('xc',self._xaxisatts,self._ax)
    self.yaxis = _axis_type('yc',self._yaxisatts,self._ay)
    self.gridaxes = [self.yaxis,self.xaxis]
    self.lon = _var_type('lon',self._lonatts,self.gridaxes,self._lonarray)
    self.lat = _var_type('lat',self._latatts,self.gridaxes,self._latarray)
    return (self._false_easting, self._false_northing, self.xaxis, self.yaxis, \
            self.gridaxes, self.lon, self.lat)
  def gen_xyll(self):
    if self.xaxis == None:
      self._gen_xyll()
    return (self.xaxis, self.yaxis, self.gridaxes, self.lon, self.lat)
      

#################################################
# Mixin for handling lat/lon coordinates.

class XYCoords (BufferBase):
  # Special records that contain coordinate info.
  # We don't want to output these directly as variables, need to decode first.
  _xycoord_nomvars = ('^^','>>','^>')
  # Grids that can be read directly from '^^','>>' records, instead of going
  # through ezqkdef (in fact, may crash ezqkdef if you try decoding them).
  _direct_grids = ('X','Y','T','+')

  @classmethod
  def _cmdline_args (cls, parser):
    super(XYCoords,cls)._cmdline_args(parser)
    parser.add_argument('--subgrid-axis', action='store_true', help=_('For data on supergrids, split the subgrids along a "subgrid" axis.  The default is to leave the subgrids stacked together as they are in the RPN file.'))
    parser.add_argument('--keep-LA-LO', action='store_true', help=_('Include LA and LO records, even if they appear to be redundant.'))

  def __init__ (self, *args, **kwargs):
    """
    subgrid_axis : bool, optional
        For data on supergrids, split the subgrids along a "subgrid" axis.
        The default is to leave the subgrids stacked together as they are in
        the RPN file.
    keep_LA_LO : bool, optional
        Include LA and LO records, even if they appear to be redundant.
    """
    self._subgrid_axis = kwargs.pop('subgrid_axis',False)
    self._keep_LA_LO = kwargs.pop('keep_LA_LO',False)
    # Tell the decoder not to process horizontal records as variables.
    self._meta_records = self._meta_records + self._xycoord_nomvars
    super(XYCoords,self).__init__(*args,**kwargs)
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
    # Special case for series data - match any of the lat/lon grids.
    if var.atts['grtyp'] in ('+','Y'):
      header = fstlir (self._meta_funit, nomvar=coordname, rank=3)
    else:
      header = fstlir (self._meta_funit, nomvar=coordname, ip1=var.atts['ig1'],
                             ip2=var.atts['ig2'], ip3=var.atts['ig3'], rank=3)
    if header is not None:
      # Override output dtype
      dtype = dtype_fst2numpy(header['datyp'],header['nbits'])
      header['d'] = header['d'].view(dtype)
      return header
    raise KeyError("Unable to find matching '%s' for '%s'"%(coordname,var.name))


  # Add horizontal coordinate info to the data.
  def _makevars (self):
    from fstd2nc.mixins import _iter_type, _var_type, _axis_type, _dim_type
    from collections import OrderedDict
    from rpnpy.librmn.interp import ezqkdef, EzscintError, ezget_nsubgrids
    from rpnpy.librmn.all import readGrid, RMNError
    import numpy as np

    # Scan through the data, and look for any use of horizontal coordinates.
    grids = OrderedDict()
    gridmaps = OrderedDict()
    lats = OrderedDict()
    lons = OrderedDict()
    # Only output 1 copy of 1D coords (e.g. could have repetitions with
    # horizontal staggering.
    coords = set()

    super(XYCoords,self)._makevars()

    # Make sure any LA/LO records get processed first, so we can apply them as
    # coordinates to other variables.
    varlist = self._varlist
    varlist = [v for v in varlist if v.name in ('LA','LO')] + [v for v in varlist if v.name not in ('LA','LO')]

    for var in varlist:
      # Don't touch variables with no horizontal grid.
      if all(a not in var.dims for a in ('i','j','station_id')):
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
      if var.atts.get('typvar','').strip() == 'T':
        key = ('T',ig1,ig2)
      else:
        key = (grtyp,ni,nj,ig1,ig2,ig3,ig4)
      if grtyp in ('Y','+'): key = key[1:]
      # Check if we already defined this grid.
      if key not in grids:

        lat = lon = xaxis = yaxis = None

        # Check if GridMap recognizes this grid.
        if grtyp not in self._direct_grids:
          try:
            grd = readGrid(self._meta_funit, var.atts.copy())
            gmap = GridMap.gen_gmap(grd)
            gmapvar = gmap.gen_gmapvar()
            gridmaps[key] = gmapvar
            (xaxis,yaxis,gridaxes,lon,lat) = gmap.gen_xyll()
          except (TypeError,EzscintError,KeyError,RMNError,ValueError):
            pass # Wasn't supported.

        # Otherwise, need to decode the information here.
        if lat is None or lon is None:

          latatts = OrderedDict()
          latatts['long_name'] = 'latitude'
          latatts['standard_name'] = 'latitude'
          latatts['units'] = 'degrees_north'
          lonatts = OrderedDict()
          lonatts['long_name'] = 'longitude'
          lonatts['standard_name'] = 'longitude'
          lonatts['units'] = 'degrees_east'

          latarray = lonarray = None
          try:
            # First, handle non-ezqkdef grids.
            if grtyp in self._direct_grids:
              latarray = self._find_coord(var,'^^')['d'].squeeze(axis=2)
              lonarray = self._find_coord(var,'>>')['d'].squeeze(axis=2)
            # Handle ezqkdef grids.
            else:
              gdid = ezqkdef (ni, nj, grtyp, ig1, ig2, ig3, ig4, self._meta_funit)
              ll = gdll(gdid)
              latarray = ll['lat']
              lonarray = ll['lon']
              xycoords = gdgaxes(gdid)
              ax = xycoords['ax'].transpose()
              ay = xycoords['ay'].transpose()
              # Convert from degenerate 2D arrays to 1D arrays.
              ax = ax[0,:]
              ay = ay[:,0]
              xaxis = _axis_type('x',{'axis':'X'},ax)
              yaxis = _axis_type('y',{'axis':'Y'},ay)

          except (TypeError,EzscintError,KeyError,RMNError,ValueError):
            pass

          # Check for LA/LO variables, and use these as the coordinates if
          # nothing else available.
          if latarray is None and var.name == 'LA':
            var.name = 'lat'
            var.atts.update(latatts)
            #grids[key] = list(var.axes)
            lats[key] = var
            continue
          if lonarray is None and var.name == 'LO':
            var.name = 'lon'
            var.atts.update(lonatts)
            grids[key] = list(var.axes)
            lons[key] = var
            continue

          if latarray is None or lonarray is None:
            warn(_("Unable to get lat/lon coordinates for '%s'")%var.name)
            continue

          # Construct lat/lon variables from latarray and lonarray.
          latarray = latarray.transpose() # Switch from Fortran to C order.
          lonarray = lonarray.transpose() # Switch from Fortran to C order.

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
            lat = _axis_type('lat',latatts,latarray)
            lon = _axis_type('lon',lonatts,lonarray)
            gridaxes = [lat,lon]

          # Case 2: lat/lon are series of points.
          elif latarray.shape[0] == 1 and lonarray.shape[0] == 1 and ('i' in var.dims or 'station_id' in var.dims):
            latarray = latarray[0,:]
            lonarray = lonarray[0,:]
            # Special case for station data
            station_id = var.getaxis('station_id')
            if station_id is not None:
              gridaxes = [station_id]
              # Subset the lat/lon to the stations that are actually found.
              # Assuming the station id (ip3) starts at 1.
              if isinstance(station_id,_axis_type):
                indices = np.array(station_id.array,dtype=int) - 1
                latarray = latarray[indices]
                lonarray = lonarray[indices]
            else:
              gridaxes = [var.getaxis('i')]
            lat = _var_type('lat',latatts,gridaxes,latarray)
            lon = _var_type('lon',lonatts,gridaxes,lonarray)

          # Case 3: General 2D lat/lon fields on X/Y coordinate system.
          elif xaxis is not None and yaxis is not None:
            gridaxes = [yaxis,xaxis]
            # Special case: have supergrid data, and the user wants to split it?
            if grtyp == 'U' and self._subgrid_axis:
              ngrids = ezget_nsubgrids(gdid)
              ny = len(yaxis.array)//ngrids
              yaxis.array = yaxis.array[:ny]
              subgrid = _dim_type('subgrid',ngrids)
              gridaxes = [subgrid,yaxis,xaxis]
              latarray = latarray.reshape(ngrids,ny,-1)
              lonarray = lonarray.reshape(ngrids,ny,-1)
            lat = _var_type('lat',latatts,gridaxes,latarray)
            lon = _var_type('lon',lonatts,gridaxes,lonarray)

          # Case 4: General 2D lat/lon fields with no coordinate system.
          elif 'i' in var.dims and 'j' in var.dims:
            gridaxes = [var.getaxis('j'),var.getaxis('i')]
            lat = _var_type('lat',latatts,gridaxes,latarray)
            lon = _var_type('lon',lonatts,gridaxes,lonarray)

        # --- End of lat/lon/xaxis/yaxis decoding.

        if lat is None or lon is None:
          warn(_("Unable to get lat/lon coordinates for '%s'")%var.name)
          continue

        # Sanity check on lat/lon - make sure we have something of the right size.
        if lat.array.shape == lat.shape and lon.array.shape == lon.shape:
          grids[key] = gridaxes
          lats[key] = lat
          lons[key] = lon
        else:
          warn(_("Wrong shape of lat/lon for '%s'")%var.name)
          continue
      
      # --- End of grid decoding.

      gridaxes = grids[key]
      lat = lats[key]
      lon = lons[key]

      # Update the var's horizontal coordinates.
      newaxes = []
      if len(gridaxes) == 1:
        newaxes = [('i',gridaxes[0])]
      elif len(gridaxes) == 2:
        newaxes = [('j',gridaxes[0]),('i',gridaxes[1])]
      elif len(gridaxes) == 3:
        newaxes = [('k',gridaxes[0]),('j',gridaxes[1]),('i',gridaxes[2])]
      else:
        warn(_("Unusual grid axes for '%s' - ignoring.")%var.name)
      dims = var.dims
      for oldname,newaxis in newaxes:
        if oldname in dims:
          var.axes[dims.index(oldname)] = newaxis

      # For 2D lat/lon, need to reference them as coordinates in order for
      # netCDF viewers to display the field properly.
      if 'lat' not in var.dims or 'lon' not in var.dims:
        coordinates = var.atts.get('coordinates',[])
        coordinates.extend([lon,lat])
        var.atts['coordinates'] = coordinates

      if key in gridmaps:
        var.atts['grid_mapping'] = gridmaps[key]

      # Throw out superfluous LA/LO variables, if lat/lon was already decoded.
      if not self._keep_LA_LO:
        if var.name == 'LA' and ('lat' in var.dims or lat in coordinates):
          var.name = None
        if var.name == 'LO' and ('lon' in var.dims or lon in coordinates):
          var.name = None

    self._varlist = [v for v in varlist if v.name is not None]

    # Special case - there are no data records, ONLY a pair of lat/lon
    # coordinates.  In this case, include the lat/lon as variables.
    if len(self._varlist) == 0 and len(lats) == 1 and len(lons) == 1:
      self._varlist = [list(lats.values())[0], list(lons.values())[0]]

