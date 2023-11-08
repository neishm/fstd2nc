###############################################################################
# Copyright 2017-2023 - Climate Research Division
#                       Environment and Climate Change Canada
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

# Define coordinateless dimensions for grid cell boundaries.
from fstd2nc.mixins import _dim_type
bnds2 = _dim_type('bnds',2)
bnds4 = _dim_type('nv',4)

# Helper function - compute bounds for the given array.
# Returned dimension is (n,2) where n is the length of the array.
# [:,0] indices are the lefthand bounds, [:,1] indices are the righthand.
def get_bounds (array, Min=None, Max=None, snap_minmax=False):
  import numpy as np
  bounds = np.empty((len(array),2),dtype=array.dtype)
  bounds[1:,0] = bounds[:-1,1] = (array[1:]+array[:-1])/2.0
  bounds[0,0] = array[0] - (bounds[0,1]-array[0])
  bounds[-1,1] = array[-1] + (array[-1]-bounds[-1,0])
  if Min is not None: bounds[0,0] = max(bounds[0,0],Min)
  if Max is not None: bounds[-1,1] = min(bounds[-1,1],Max)
  if snap_minmax:
    if bounds[0,0] - Min < array[0] - bounds[0,0]:
      bounds[0,0] = Min
    if Max - bounds[-1,1] < bounds[-1,1] - array[-1]:
      bounds[-1,1] = Max
  return bounds

# Helper function - compute lat/lon cell boundaries from x/y cell boundaries.
def get_ll_vertices (gdid, x, y):
  from rpnpy.librmn.all import gdllfxy
  import numpy as np
  assert x.shape[1] == y.shape[1] == 2
  xpts = np.empty((y.shape[0],x.shape[0],4),dtype=x.dtype)
  ypts = np.empty((y.shape[0],x.shape[0],4),dtype=x.dtype)
  xpts[:,:,0] = x[None,:,0]
  xpts[:,:,3] = x[None,:,0]
  xpts[:,:,1] = x[None,:,1]
  xpts[:,:,2] = x[None,:,1]
  ypts[:,:,0] = y[:,None,0]
  ypts[:,:,1] = y[:,None,0]
  ypts[:,:,2] = y[:,None,1]
  ypts[:,:,3] = y[:,None,1]
  ll = gdllfxy(gdid, xpts.flatten(), ypts.flatten())
  return ll['lon'].reshape(xpts.shape), ll['lat'].reshape(ypts.shape)
  #return xpts, ypts

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
  # This is the value used in librmn, and also corresponds to the Normal
  # Sphere approximation referenced in
  # https://proj.org/usage/ellipsoids.html#built-in-ellipsoid-definitions
  _earth_radius = 6370997.
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
  def gen_gmap(cls, grd, no_adjust_rlon=False):
    import numpy as np
    grref = grd['grref'].upper() 
    if grref == 'E' :
      # Special case: 'E' grid is not actually rotated.
      if np.allclose(grd['lat0'],grd['rlat0']) and np.allclose(grd['lon0'],grd['rlon0']):
        return LatLon(grd)
      # Usual case: 'E' grid is rotated.
      return RotLatLon(grd, no_adjust_rlon=no_adjust_rlon)
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
  def gen_ll(self,bounds=False):
    from collections import OrderedDict
    from rpnpy.librmn.all import gdll 
    from fstd2nc.mixins import _axis_type, _var_type
    import numpy as np
    ll = gdll(self._grd['id'])
    self._lonarray = ll['lon'][:,0]
    self._lonatts = OrderedDict()
    self._lonatts['long_name'] = 'longitude'
    self._lonatts['standard_name'] = 'longitude'
    self._lonatts['units'] = 'degrees_east'
    self._lonatts['axis'] = 'X'
    self._latarray = ll['lat'][0,:]
    self._latatts = OrderedDict()
    self._latatts['long_name'] = 'latitude'
    self._latatts['standard_name'] = 'latitude'
    self._latatts['units'] = 'degrees_north'
    self._latatts['axis'] = 'Y'
    self.lon = _axis_type('lon',self._lonatts,self._lonarray)
    self.lat = _axis_type('lat',self._latatts,self._latarray)
    self.gridaxes = [self.lat,self.lon]
    # Ensure monotonicity of longitude field.
    # (gdll may sometimes wrap last longitude to zero).
    # Taken from old fstd_core.c code.
    if len(self._lonarray) >= 3 and self._lonarray[-2] > self._lonarray[-3] and self._lonarray[-1] < self._lonarray[-2]:
      self._lonarray[-1] += 360.
    # Use the ax/ay values instead of derived lat/lon to maintain perfect
    # precision of the coordinates.
    # NOTE: Could probably remove the gdll function call (and monotonicity fix)
    # and use these values directly.
    if 'ax' in self._grd and 'ay' in self._grd:
      if np.allclose(self._grd['ax'][:,0],self._lonarray,atol=1e-5):
        self._lonarray[:] = self._grd['ax'][:,0]
      if np.allclose(self._grd['ay'][0,:],self._latarray,atol=1e-5):
        self._latarray[:] = self._grd['ay'][0,:]
    # Add lat/lon boundaries.
    if bounds:
      lon_bnds = get_bounds(self._lonarray)
      lon_bnds = _var_type("lon_bnds",{},[self.lon,bnds2],lon_bnds)
      self.lon.atts["bounds"] = lon_bnds
      lat_bnds = get_bounds(self._latarray,Min=-90,Max=90,snap_minmax=True)
      lat_bnds = _var_type("lat_bnds",{},[self.lat,bnds2],lat_bnds)
      self.lat.atts["bounds"] = lat_bnds
    return (self.gridaxes, self.lon, self.lat)
  # Generic gen_xyll function.
  def gen_xyll(self, bounds=False):
    gridaxes, lon, lat = self.gen_ll(bounds=bounds)
    return (None, None, gridaxes, lon, lat)

class RotLatLon(GridMap):
  def __init__(self, *args, **kwargs):
    from rpnpy.librmn.all import egrid_rll2ll, egrid_ll2rll
    adjust_rlon = not kwargs.pop('no_adjust_rlon')
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
    if adjust_rlon:
      orig_ax = self._ax
      ax_adjust = self._north_pole_grid_longitude
      self._ax = self._ax - self._north_pole_grid_longitude
      self._north_pole_grid_longitude = 0
      # Make sure rlon axis is still in range.
      if self._ax.max() >= 360.:
        self._ax -= 360.
        ax_adjust += 360.
      if self._ax.min() <= -180.:
        self._ax += 360.
        ax_adjust -= 360.
      # Compensate for rounding error when undoing the adjustment.
      # This should make the values bit-for-bit identical if writing *back*
      # to FST format.
      self._ax = orig_ax - ax_adjust
    self._ay = self._grd['ay'][0,:]
  def gen_gmapvar(self):
    from fstd2nc.mixins import _var_type
    import numpy as np
    self._atts['grid_mapping_name'] = 'rotated_latitude_longitude'
    self._atts['earth_radius'] = self._earth_radius
    self._atts['grid_north_pole_latitude'] = self._grid_north_pole_latitude
    self._atts['grid_north_pole_longitude'] = self._grid_north_pole_longitude
    self._atts['north_pole_grid_longitude'] = self._north_pole_grid_longitude 
#   Set the optional grid mapping parameter 'north_pole_grid_longitude' to 0 to avoid 
#   some problems, such as the conversion from netcdf to grib performed by some tools
#    self._atts['north_pole_grid_longitude'] = 0.
#    self._atts['longitude_of_prime_meridian'] = 0.
    # Grid mapping variable
    self.gmap = _var_type(self._name,self._atts,[],np.array(b""))
    return self.gmap
  # Generate latitudes and longitudes in rotated pole grid 
  # and true latitudes and longitudes
  def gen_xyll(self, bounds=False):
    from collections import OrderedDict
    from rpnpy.librmn.all import gdll
    from fstd2nc.mixins import _var_type, _axis_type
    import numpy as np
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
    # Get boundary info.
    if bounds:
      x_bnds = get_bounds(self._ax)
      x_bnds = _var_type('rlon_bnds',{},[self.xaxis,bnds2],x_bnds)
      self.xaxis.atts['bounds'] = x_bnds
      y_bnds = get_bounds(self._ay,Min=-90,Max=90,snap_minmax=True)
      y_bnds = _var_type('rlat_bnds',{},[self.yaxis,bnds2],y_bnds)
      self.yaxis.atts['bounds'] = y_bnds
      # Compute bounds as array indices (needed for coordinate transformation).
      ax_bnds = get_bounds(np.arange(1,len(self._ax)+1,dtype='float32'))
      ay_bnds = get_bounds(np.arange(1,len(self._ay)+1,dtype='float32'))
      # Compute lat/lon of cell boundaries.
      lon_bnds, lat_bnds = get_ll_vertices (self._grd['id'], ax_bnds, ay_bnds)
      lon_bnds = _var_type('lon_bnds',{},[self.yaxis,self.xaxis,bnds4],lon_bnds)
      self.lon.atts['bounds'] = lon_bnds
      lat_bnds = _var_type('lat_bnds',{},[self.yaxis,self.xaxis,bnds4],lat_bnds)
      self.lat.atts['bounds'] = lat_bnds
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
    if self._grd['north']:
      pole = gdxyfll (self._grd['id'], 90, 0)
    else:
      pole = gdxyfll (self._grd['id'], -90, 0)
    px = np.rint(pole['x'][0] - 1) * self._res
    py = np.rint(pole['y'][0] - 1) * self._res
    self._false_easting =  px - self._ax[0]
    self._false_northing = py - self._ay[0]
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
  def gen_xyll(self, bounds=False):
    from fstd2nc.mixins import _var_type
    import numpy as np
    if self.xaxis == None:
      self._gen_xyll()
    # Get boundary info.
    if bounds:
      x_bnds = get_bounds(self._ax)
      x_bnds = _var_type('xc_bnds',{},[self.xaxis,bnds2],x_bnds)
      self.xaxis.atts['bounds'] = x_bnds
      y_bnds = get_bounds(self._ay)
      y_bnds = _var_type('yc_bnds',{},[self.yaxis,bnds2],y_bnds)
      self.yaxis.atts['bounds'] = y_bnds
      # Compute bounds as array indices (needed for coordinate transformation).
      ax_bnds = get_bounds(np.arange(1,len(self._ax)+1,dtype='float32'))
      ay_bnds = get_bounds(np.arange(1,len(self._ay)+1,dtype='float32'))
      # Compute lat/lon of cell boundaries.
      lon_bnds, lat_bnds = get_ll_vertices (self._grd['id'], ax_bnds, ay_bnds)
      lon_bnds = _var_type('lon_bnds',{},[self.yaxis,self.xaxis,bnds4],lon_bnds)
      self.lon.atts['bounds'] = lon_bnds
      lat_bnds = _var_type('lat_bnds',{},[self.yaxis,self.xaxis,bnds4],lat_bnds)
      self.lat.atts['bounds'] = lat_bnds
    return (self.xaxis, self.yaxis, self.gridaxes, self.lon, self.lat)
      

#################################################
# Mixin for handling lat/lon coordinates.

class XYCoords (BufferBase):
  # Tell the decoder not to process horizontal records as variables.
  @classmethod
  def _meta_records (cls):
    return super(XYCoords,cls)._meta_records() + (b'^^',b'>>',b'^>')
  # Grids that can be read directly from '^^','>>' records, instead of going
  # through ezqkdef (in fact, may crash ezqkdef if you try decoding them).
  _direct_grids = ('X','Y','T','+','O','M')

  @classmethod
  def _cmdline_args (cls, parser):
    from argparse import SUPPRESS
    super(XYCoords,cls)._cmdline_args(parser)
    parser.add_argument('--subgrid-axis', action='store_true', help=_('For data on supergrids, split the subgrids along a "subgrid" axis.  The default is to leave the subgrids stacked together as they are in the RPN file.'))
    parser.add_argument('--keep-LA-LO', action='store_true', help=_('Include LA and LO records, even if they appear to be redundant.'))
    parser.add_argument('--no-adjust-rlon', action='store_true', help=_('For rotated grids, do NOT adjust rlon coordinate to keep the range in -180..180.  Allow the rlon value to be whatever librmn says it should be.'))
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--bounds', action='store_true', default=False, help=_('Include grid cell boundaries in the output.'))
    group.add_argument('--no-bounds', action='store_false', dest='bounds', help=SUPPRESS)

  def __init__ (self, *args, **kwargs):
    """
    subgrid_axis : bool, optional
        For data on supergrids, split the subgrids along a "subgrid" axis.
        The default is to leave the subgrids stacked together as they are in
        the RPN file.
    keep_LA_LO : bool, optional
        Include LA and LO records, even if they appear to be redundant.
    no_adjust_rlon : bool, optional
        For rotated grids, do NOT adjust rlon coordinate to keep the range
        in -180..180.  Allow the rlon value to be whatever librmn says it
        should be.
    bounds : bool, optional
        Include grid cell boundaries in the output.
    """
    self._subgrid_axis = kwargs.pop('subgrid_axis',False)
    self._keep_LA_LO = kwargs.pop('keep_LA_LO',False)
    self._no_adjust_rlon = kwargs.pop('no_adjust_rlon',False)
    self._bounds = kwargs.pop('bounds',False)
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
    from fstd2nc.mixins.fstd import dtype_fst2numpy
    # Special case for series data - match any of the lat/lon grids.
    if var.atts['grtyp'] in ('+','Y'):
      header = self._fstlir (nomvar=coordname)
      # Make sure this actually matches a grid of the correct shape.
      if header['ni'] != var.atts['ni'] or header['nj'] != var.atts['nj']:
        header = None
    else:
      header = self._fstlir (nomvar=coordname, ip1=var.atts['ig1'],
                             ip2=var.atts['ig2'], ip3=var.atts['ig3'])
    if header is not None:
      # Override output dtype
      dtype = dtype_fst2numpy(int(header['datyp']),int(header['nbits']))
      header['d'] = header['d'][:,:,None].view(dtype)
      return header
    raise KeyError("Unable to find matching '%s' for '%s'"%(coordname,var.name))


  # Add horizontal coordinate info to the data.
  def _makevars (self):
    from fstd2nc.mixins import _iter_type, _chunk_type, _var_type, _axis_type, _dim_type
    from collections import OrderedDict
    from rpnpy.librmn.interp import ezqkdef, EzscintError, ezget_nsubgrids
    from rpnpy.librmn.all import cxgaig, ezgdef_fmem, ezgdef_supergrid, ezqkdef, decodeGrid, RMNError
    import numpy as np

    # Save a copy of librmn grid ids, which might be useful for other mixins.
    self._gids = np.empty(self._nrecs,dtype=int)
    self._gids[:] = -1

    # Scan through the data, and look for any use of horizontal coordinates.
    grids = OrderedDict()
    gridmaps = OrderedDict()
    gids = OrderedDict()
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
      # For '#' grid variables, ignore ni,nj,ig3,ig4
      # (they are different for each tile).
      elif grtyp == '#':
        key = (grtyp,ig1,ig2)
      else:
        key = (grtyp,ni,nj,ig1,ig2,ig3,ig4)
      if grtyp in ('Y','+'): key = key[1:]
      # Check if we already defined this grid.
      if key not in grids:

        lat = lon = xaxis = yaxis = None

        # Check if GridMap recognizes this grid.
        if grtyp not in self._direct_grids:
          atts = var.atts.copy()
          # For '#' grid, extract full coordinates of parent grid.
          if grtyp == '#':
            match = (self._headers['ip1'] == atts['ig1']) & (self._headers['ip2'] == atts['ig2'])
            match_nj = np.where(match & (self._headers['nomvar'] == b'^^  '))[0]
            match_ni = np.where(match & (self._headers['nomvar'] == b'>>  '))[0]
            if len(match_nj) >= 1 and len(match_ni) >= 1:
              atts['nj'] = int(self._headers['nj'][match_nj[0]])
              atts['ni'] = int(self._headers['ni'][match_ni[0]])
              atts['ig3'] = 1
              atts['ig4'] = 1
          try:
            # Get grid object from librmn.
            # Avoid readGrid, because it requires its own access to the
            # files(s), which may not be readily available on disk.
            if grtyp in ('Z','#'):
              ref = self._find_coord(var,b'>>  ')
              grd = atts.copy()
              grd['grref'] = ref['grtyp']
              grd['ig1'] = int(ref['ig1'])
              grd['ig2'] = int(ref['ig2'])
              grd['ig3'] = int(ref['ig3'])
              grd['ig4'] = int(ref['ig4'])
              grd['ax'] = ref['d'].squeeze()
              grd['ay'] = self._find_coord(var,b'^^  ')['d'].squeeze()
              grd = ezgdef_fmem(grd)
            else:
              grd = ezqkdef (ni, nj, grtyp, ig1, ig2, ig3, ig4)
            gids[key] = grd
            grd = decodeGrid(grd)
            gmap = GridMap.gen_gmap(grd,no_adjust_rlon=self._no_adjust_rlon)
            gmapvar = gmap.gen_gmapvar()
            gridmaps[key] = gmapvar
            (xaxis,yaxis,gridaxes,lon,lat) = gmap.gen_xyll(bounds=self._bounds)
          except (TypeError,EzscintError,KeyError,RMNError,ValueError,AttributeError):
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
              latarray = self._find_coord(var,b'^^  ')['d'].squeeze(axis=2)
              lonarray = self._find_coord(var,b'>>  ')['d'].squeeze(axis=2)
            # Handle ezqkdef grids.
            else:
              ###
              # U grids aren't yet supported by in-memory ezqkdef
              # (requires ext
              if grtyp == 'U':
                ref = self._find_coord(var,b'^>  ')
                data = ref['d'].flatten()
                nsubgrids = int(data[2])
                subgrids = []
                subdata = data[5:]
                for i in range(nsubgrids):
                  sub_ni = int(subdata[0])
                  sub_nj = int(subdata[1])
                  # Loosely based on steps in Lire_enrUvercode1 of librmn.
                  sub_ig1, sub_ig2, sub_ig3, sub_ig4 = cxgaig('E',*subdata[6:10])
                  sub_ax = subdata[10:10+sub_ni]
                  sub_ay = subdata[10+sub_ni:10+sub_ni+sub_nj]
                  subgrids.append(ezgdef_fmem(sub_ni, sub_nj, 'Z', 'E', sub_ig1, sub_ig2, sub_ig3, sub_ig4, sub_ax, sub_ay))
                  subdata = subdata[10+sub_ni+sub_nj:]
                gdid = ezgdef_supergrid(ni, nj, 'U', 'F', 1, subgrids)
              else:
                #TODO: check if this case still gets triggered?
                # GridMap grids (N,S,A,B,L,G,Z,E) are already handled.
                # Direct grids (X,Y,T,+,O,M) are already handled.
                # Supergrids (U) are already handled.
                # What's left?
                gdid = ezqkdef (ni, nj, grtyp, ig1, ig2, ig3, ig4, 0)
              gids[key] = gdid
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
          if latarray.shape[0] > 1 and lonarray.shape[1] > 1 and np.allclose(latarray,meanlat) and np.allclose(lonarray,meanlon):
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
            # Add boundary information.
            if self._bounds:
              lon_bnds = get_bounds(lonarray)
              lon_bnds = _var_type('lon_bnds',{},[lon,bnds2],lon_bnds)
              lon.atts['bounds'] = lon_bnds
              lat_bnds = get_bounds(latarray)
              lat_bnds = _var_type('lat_bnds',{},[lat,bnds2],lat_bnds)
              lat.atts['bounds'] = lat_bnds
            gridaxes = [lat,lon]

          # Case 2: lat/lon are series of points.
          elif 1 in latarray.shape and 1 in lonarray.shape and ('i' in var.dims or 'station_id' in var.dims):
            latarray = latarray.squeeze()
            lonarray = lonarray.squeeze()
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

        # If lat/lon are 1D axes, then add an 'axis' attribute.
        if isinstance(lat,_axis_type):
          lat.atts = dict(lat.atts, axis='Y')
        if isinstance(lon,_axis_type):
          lon.atts = dict(lon.atts, axis='X')

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
      # Store librmn grid ids for possible use in other parts of code.
      if key in gids and isinstance(var,_iter_type):
        ind = var.record_id.flatten()
        ind = ind[ind>=0]
        if len(ind) > 0:
          self._gids[ind] = gids[key]

      # Update the var's horizontal coordinates.
      newaxes = []
      if len(gridaxes) == 1:
        newaxes = [('i',gridaxes[0])]
      elif len(gridaxes) == 2:
        newaxes = [('j',gridaxes[0]),('i',gridaxes[1])]
      elif len(gridaxes) == 3:
        newaxes = [('j',((gridaxes[0],gridaxes[1]))),('i',gridaxes[2])]
      else:
        warn(_("Unusual grid axes for '%s' - ignoring.")%var.name)
      dims = var.dims
      for oldname,newaxis in newaxes:
        if oldname in dims:
          ind = dims.index(oldname)
          # Special case for splitting an inner axes into multiple components
          # (such as adding a subgrid index).
          if isinstance(newaxis,tuple):
            var.axes = var.axes[:ind] + list(newaxis) + var.axes[ind+1:]
            dims = var.dims # Number of dimensions is updated.
          else:
            var.axes[ind] = newaxis

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

    # Only create distinct grid_mappings when they're actually distinct.
    # (may have duplicates because of different ni,nj from staggered grid.)
    gridmaps = dict()
    for var in varlist:
      if 'grid_mapping' in var.atts:
        gmapvar = var.atts['grid_mapping']
        key = tuple(sorted(gmapvar.atts.items()))
        if key not in gridmaps:
          gridmaps[key] = gmapvar
        gmapvar = gridmaps[key]
        var.atts['grid_mapping'] = gmapvar

    self._varlist = [v for v in varlist if v.name is not None]

    # Reassemble '#' grids
    tiles = OrderedDict()
    for var in self._varlist:
      if var.atts.get('grtyp',None) == '#':
        key = tuple([var.atts.get(n,None) for n in self._var_id if n not in ('ni','nj','ig3','ig4')])
        tiles.setdefault(key,[]).append(var)
    self._varlist = [v for v in self._varlist if v.atts.get('grtyp',None) != '#']
    for key, vars in tiles.items():
      #TODO
      # Get i/j coordinates.
      ig3_list = sorted(set([v.atts['ig3'] for v in vars]))
      ig4_list = sorted(set([v.atts['ig4'] for v in vars]))
      # Arrange the tiles in their natural 2D order.
      var_tiles = np.empty((len(ig4_list),len(ig3_list)), dtype='object')
      for v in vars:
        ind = (ig4_list.index(v.atts['ig4']), ig3_list.index(v.atts['ig3']))
        var_tiles[ind] = v
      chunks = [(1,)*n for n in vars[0].record_id.shape] + [tuple(v.atts['nj'] for v in var_tiles[:,0])] + [tuple(v.atts['ni'] for v in var_tiles[0,:])]
      record_id = np.empty(vars[0].shape[:-2]+var_tiles.shape, dtype=int)
      record_id[()] = -1
      for ind in np.ndindex(var_tiles.shape):
        record_id[...,ind] = var_tiles[ind].record_id
      var = _chunk_type (vars[0].name, vars[0].atts, vars[0].axes, vars[0].dtype, chunks, record_id)
      self._varlist.append(var)

    # Special case - there are no data records, ONLY a pair of lat/lon
    # coordinates.  In this case, include the lat/lon as variables.
    if len(self._varlist) == 0 and len(lats) == 1 and len(lons) == 1:
      self._varlist = [list(lats.values())[0], list(lons.values())[0]]
      # Add grid mapping info.
      self._varlist.extend(gridmaps.values())

  # Decode horizontal grid metadata from variables back into table.
  def _unmakevars (self):
    from fstd2nc.mixins import _iter_type, _chunk_type, _var_type, _axis_type, _dim_type
    import numpy as np
    import rpnpy.librmn.all as rmn
    from math import sin, cos, asin, atan2, pi, sqrt

    # Helper - add a coordinate record
    def add_coord (name,nj,ni,values,**atts):
      import dask.array as da
      dims = []
      if nj > 1: dims.append(_dim_type('j',nj))
      if ni > 1: dims.append(_dim_type('i',ni))
      # Need to carefully encode the arrays so they are scalar object containing the data.
      array = np.empty(1,dtype=object)
      array[0] = da.from_array(values, chunks=-1)
      array = array.squeeze()
      self._varlist.append(_iter_type(name,atts,dims,'float32',array))

    # Helper method - determine rotated lat/lon grid parameters from
    # 2D lat/lon fields.
    def get_rotated_grid_params (ax, ay, lat, lon):
      ni = len(ax)
      nj = len(ay)
      # Convert to radians.
      ax = ax.astype('float64') * pi/180
      ay = ay.astype('float64') * pi/180
      ax = np.repeat(ax.reshape(1,-1),nj,axis=0)
      ay = np.repeat(ay.reshape(-1,1),ni,axis=1)
      lat = lat.astype('float64') * pi/180
      lon = lon.astype('float64') * pi/180
      def cartesian (lat, lon):
        return [np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)]
      # Find transformation from grid coordinates to real coordinates.
      # Use least-squares fit of lat/lon fields.
      rotated_coords = np.array(cartesian(ay.flatten(),ax.flatten())).T
      latlon_coords = np.array(cartesian(lat.flatten(),lon.flatten())).T
      # Get transformation matrix.
      transform, res, rank, s = np.linalg.lstsq(rotated_coords, latlon_coords, rcond=None)
      transform = transform.T
      # Find coordinates of grid north pole.
      x, y, z = np.dot(transform,[0,0,1])
      gpole_lat = asin(z) * 180/pi
      gpole_lon = (atan2(y,x) * 180/pi) % 360
      # Find coordinates of north pole in grid space.
      x, y, z = np.dot(transform.T,[0,0,1])
      npole = (atan2(y,x) * 180/pi) % 360
      # Estimate xlat1, xlon1.
      x, y, z = np.dot(transform,[-1,0,0])
      xlat1 = asin(z) * 180/pi
      xlon1 = (atan2(y,x) * 180/pi) % 360
      # Round to encoding precision.
      xlat1_round = np.round(xlat1*40)/40
      xlon1_round = np.round(xlon1*40)/40
      x_round, y_round, z_round = cartesian(xlat1_round*pi/180,xlon1_round*pi/180)
      # Re-orthogonalize gpole_lat, gpole_lon if we corrected for an error
      # from rounding.
      # 2D lat/lon fields are in single precision, so it will introduce some
      # errors in the encoding.
      x, y, z = np.dot(transform,[0,0,1])
      corr = x*x_round + y*y_round + z*z_round
      if abs(corr) > 0:
        x -= corr*x_round
        y -= corr*y_round
        z -= corr*z_round
        mag = np.sqrt(x**2 + y**2 + z**2)
        x /= mag
        y /= mag
        z /= mag
        gpole_lat = asin(z) * 180/pi
        gpole_lon = (atan2(y,x) * 180/pi) % 360
      return dict(
        grid_north_pole_latitude=gpole_lat,
        grid_north_pole_longitude=gpole_lon,
        north_pole_grid_longitude=npole
      )

    # Helper method - find grid longitudes along rotated equator that can be
    # losslessly encoded in FST.
    def get_nice_xlats_xlons (gpole_lat, gpole_lon):
      xlons = np.linspace(0,359.975,14400)
      xlats =  np.arctan2(-cos(gpole_lon)*cos(gpole_lat)*np.cos(xlons*pi/180) - \
                        sin(gpole_lon)*cos(gpole_lat)*np.sin(xlons*pi/180),  \
                        sin(gpole_lat)
               ) * 180/pi
      eps = (xlats*40)%1
      eps[eps>0.5] -= 1
      eps = abs(eps/40)
      match = np.where(np.isclose(eps,0))[0]
      # If only two matches (on opposite sides?), then adjust criteria to find 
      # another point.
      # First, look for a 'close' integer coordinate.
      # Useful, for instance, if detecting a 'yin' grid where the grid
      # rotation was inferred from 2D stacked lat/lon fields (less precision in
      # inferred params).
      if len(match) == 2:
        eps = xlats%1
        eps[eps>0.5] -= 1
        eps = abs(eps)
        match = np.where(np.isclose(eps,0,atol=1e-6))[0]
      # Otherwise, use another point 90 degrees away.
      # Less encoding precision, but this is what seems to be used for 'yang'
      # grids.  So this code block is effectively a yang grid detector.
      if len(match) == 2:
        def find_orthogonal(ind):
          xlat = xlats[ind]
          xlon = xlons[ind]
          x = cos(xlon*pi/180)*cos(xlat*pi/180)
          y = sin(xlon*pi/180)*cos(xlat*pi/180)
          z = sin(xlat*pi/180)
          X = np.cos(xlons*pi/180)*np.cos(xlats*pi/180)
          Y = np.sin(xlons*pi/180)*np.cos(xlats*pi/180)
          Z = np.sin(xlats*pi/180)
          eps = abs(x*X+y*Y+z*Z)
          ind2 = np.argmin(eps)
          if ind2 < ind: ind2 += 7200
          return ind2
        match = np.array([match[0], find_orthogonal(match[0]), match[1], find_orthogonal(match[1])])
      return xlats[match], xlons[match]

    # Helper method - get grid centre from given attributes.
    def get_rotated_centre (gpole_lat, gpole_lon, npole):
      cosp = cos(npole*pi/180)
      sinp = sin(npole*pi/180)
      X = cosp * sin(gpole_lat) * cos(gpole_lon) + sinp * sin(gpole_lon)
      Y = cosp * sin(gpole_lat) * sin(gpole_lon) - sinp * cos(gpole_lon)
      Z = -cosp * cos(gpole_lat)
      xlat1 = asin(Z) * 180 / pi
      xlon1 = atan2(Y,X) * 180 / pi
      if xlon1 < 0: xlon1 += 360
      return xlat1, xlon1

    # Helper method - grid other grid reference point (east of centre).
    def get_xlat_xlon (gpole_lat, gpole_lon, npole):
      xlat1, xlon1 = get_rotated_centre (gpole_lat, gpole_lon, npole)
      xlats, xlons = get_nice_xlats_xlons (gpole_lat, gpole_lon)
      # Find index of the centre.
      where_centre = np.where(np.isclose(xlons, xlon1))[0]
      if len(where_centre) == 0:
        raise Exception("TODO: case where no nice value available")
      # Find another point just to the right.
      ind = (where_centre[0]+1) % len(xlons)
      xlat2 = xlats[ind]
      xlon2 = xlons[ind]
      return xlat1, xlon1, xlat2, xlon2

    # Helper method - get best xlat1/xlon1/xlat2/xlon2 rotation parameters
    # for the given CF rotation attributes.
    def get_rotation_params (atts, ax):
      gpole_lat = atts['grid_north_pole_latitude'] * pi / 180
      gpole_lon = atts['grid_north_pole_longitude'] * pi / 180
      npole = atts.get('north_pole_grid_longitude',0.)
      # Case 1: have non-zero npole, so no adjustment was done to grid
      # longitudes.
      if npole != 0.0:
        xlat1, xlon1, xlat2, xlon2 = get_xlat_xlon (gpole_lat, gpole_lon, npole)
        return xlat1, xlon1, xlat2, xlon2, 0.0
      # Case 2: npole is zero, so need to determine the adjustment.
      else:
        # Get all good possible values of xlat1, xlon1.
        # (encodable without precision loss).
        xlats, xlons = get_nice_xlats_xlons (gpole_lat, gpole_lon)
        # Find possible values for npole.
        # Expect that ax values should be centered over the grid.
        # (mid-point should be 180 degrees).
        npole_guess = 180-(ax[0]+ax[-1])/2
        # Find which xlon value would give us the closest match to this npole.
        # Check for grid rotation, adjust if necessary.
        # Determine location of geographic pole in model coordinates.
        X = -sin(gpole_lat)*cos(gpole_lon)*np.cos(xlats*pi/180)*np.cos(xlons*pi/180) \
          - sin(gpole_lat)*sin(gpole_lon)*np.cos(xlats*pi/180)*np.sin(xlons*pi/180) \
          + cos(gpole_lat)*np.sin(xlats*pi/180)
        Y = -sin(gpole_lon)*np.cos(xlats*pi/180)*np.cos(xlons*pi/180) \
          + cos(gpole_lon)*np.cos(xlats*pi/180)*np.sin(xlons*pi/180)
        check = np.arctan2(Y,X)*180/pi
        closeness = (npole_guess - check)%360
        closeness[closeness>=180] -= 360
        closeness = abs(closeness)
        ind = np.argmin(closeness)
        npole = check[ind]
        xlat1, xlon1, xlat2, xlon2 = get_xlat_xlon (gpole_lat, gpole_lon, npole)
        ax_adjust = npole
        if ax_adjust < 0: ax_adjust += 360
        return xlat1, xlon1, xlat2, xlon2, ax_adjust

    # Helper method - apply adjustment to ax while recovering as much
    # accuracy as possible.
    def adjust_ax (ax, ax_adjust):
      ax = np.array(ax,dtype='float64')
      if np.max(ax+ax_adjust) >= 360.0:
        ax_adjust -= 360.0
      ax += ax_adjust
      # Check if we can apply a correction.
      # For LAM coming from yin-yang grid, check for integer start/end.
      # For mass grid, should have integer values.
      if np.allclose(ax[0],round(ax[0])) and np.allclose(ax[-1],round(ax[-1])):
        x0 = round(ax[0])
        x1 = round(ax[-1])
        x = np.linspace(x0,x1,len(ax))
        if np.allclose(ax,x):
          return x.astype('float32')
      # For staggered grid, should be offset by dx/2?
      dx = ax[1]-ax[0]
      if np.allclose(ax[0]-dx/2,round(ax[0]-dx/2)) and np.allclose(ax[-1]-dx/2,round(ax[-1]-dx/2)):
        x0 = round(ax[0]-dx/2)
        x1 = round(ax[-1]-dx/2)
        dx = (x1-x0)/(len(ax)-1)
        x = np.linspace(x0+dx/2,x1+dx/2,len(ax))
        if np.allclose(ax,x):
          return x.astype('float32')
      # Check if some inner points are integers?
      for i in range(1,len(ax)//3):
        if np.allclose(ax[i],round(ax[i])) and np.allclose(ax[-i-1],round(ax[-i-1])):
          xa = round(ax[i])
          xb = round(ax[-i-1])
          x = np.zeros(len(ax),float)
          x[i:-i] = np.linspace(xa,xb,len(ax)-2*i)
          dx = (xb-xa) / (len(ax)-2*i-1)
          x[:i] = np.linspace(xa-i*dx,xa-dx,i)
          x[-i:] = np.linspace(xb+dx,xb+dx*i,i)
          if np.allclose(ax,x):
            return x.astype('float32')
      return ax.astype('float32')

    gauss_table = {}
    lgrid_table = {}
    zlgrid_table = {}
    zegrid_table = {}
    yygrid_table = {}
    stacked_yaxes = {}
    projection_table = {}
    # Pull out projection variables.
    for var in self._varlist:
      if 'grid_mapping' in var.atts:
        projection_table[var.atts['grid_mapping']] = None
    for var in self._varlist:
      if var.name in projection_table.keys():
        projection_table[var.name] = var
    self._varlist = [v for v in self._varlist if v.name not in projection_table.keys()]
    # Loop over all variables, encode grid info.
    for var_ind, var in enumerate(self._varlist):
      # Skip variables already processed into records.
      if isinstance(var, _iter_type): continue
      axis_codes = [a.atts.get('axis','') if isinstance(a,_axis_type) else None for a in var.axes]
      # Skip records without any geophysical connection.
      if 'X' not in axis_codes or 'Y' not in axis_codes: continue
      xind = axis_codes.index('X')
      yind = axis_codes.index('Y')
      # Detect lat/lon axes.
      xaxis = var.axes[xind]
      yaxis = var.axes[yind]
      ni = len(xaxis)
      nj = len(yaxis)

      # Transpose to expected order.
      order = [i for i in range(len(var.axes)) if i not in (xind,yind)]
      order = order + [yind, xind]
      if order != list(range(len(var.axes))):
        var.axes = [var.axes[i] for i in order]
        var.array = var.array.transpose(*order)
      # Identify inner axes.
      var.axes[-2] = _dim_type('j',nj)
      var.axes[-1] = _dim_type('i',ni)
      # Special case for subgrid axis
      # Flatten it out so the yin/yang grids are stacked together.
      if 'subgrid' in var.dims:
        # First transpose so the subgrid axis is just before y axis.
        i = var.dims.index('subgrid')
        order = list(range(i)) + list(range(i+1,len(var.dims)-2)) + [i,-2,-1]
        if order != list(range(len(var.axes))):
          var.axes = [var.axes[i] for i in order]
          var.array = var.array.transpose(*order)
        # Next, make a stacked y axis.
        if id(yaxis) not in stacked_yaxes:
          stacked_yaxes[id(yaxis)] = _axis_type('j',yaxis.atts,np.concatenate([yaxis.array,yaxis.array]))
        yaxis = stacked_yaxes[id(yaxis)]
        var.array = var.array.reshape(var.shape[:-3] + (var.shape[-2]*2,var.shape[-1]))
        var.axes = var.axes[:-3] + [yaxis] + var.axes[-1:]
        nj *= 2
      ###
      # Process the variable.
      outer_shape = [len(a) for a in var.axes[:-2]]
      record_id = np.zeros(outer_shape, dtype=object)
      for ind in np.ndindex(tuple(outer_shape)):
        record_id[ind] = var.array[ind]
      var = _iter_type (var.name, var.atts, var.axes, var.array.dtype, record_id)
      self._varlist[var_ind] = var

      # Find best grtyp to use.
      # Use the specified one, if it works for the data.
      grtyp = var.atts.get('grtyp',None)
      grref = var.atts.get('grref',None)

      # Case 1: data is on lat/lon grid (and no coordinate records needed)
      have_lon = xaxis.name == 'lon' or xaxis.atts.get('standard_name',None) == 'longitude'
      have_lat = yaxis.name == 'lat' or yaxis.atts.get('standard_name',None) == 'latitude'
      if have_lon and have_lat:
        # 'A' grid
        #TODO: hemispheric
        agrid_lat = (np.arange(nj)+0.5)/nj*180-90
        agrid_lon = np.arange(ni)/ni*360
        if np.allclose(yaxis.array,agrid_lat,atol=1e-5) and np.allclose(xaxis.array,agrid_lon,atol=1e-5) and grtyp not in ('L','Z'):
          var.atts.update(grtyp='A',ig1=0,ig2=0,ig3=0,ig4=0)
          continue
        # 'B' grid
        #TODO: hemispheric
        bgrid_lat = np.linspace(-90,90,nj)
        bgrid_lon = np.linspace(0,360,ni)
        if np.allclose(yaxis.array,bgrid_lat,atol=1e-5) and np.allclose(xaxis.array,bgrid_lon,atol=1e-5) and grtyp not in ('L','Z'):
          var.atts.update(grtyp='B',ig1=0,ig2=0,ig3=0,ig4=0)
          continue
        # 'G' grid
        #TODO: hemispheric
        if nj not in gauss_table:
          gid = rmn.defGrid_G(ni,nj)
          gauss_table[nj] = rmn.gdll(gid)['lat']
        if np.allclose(yaxis.array,gauss_table[nj],atol=1e-5) and np.allclose(xaxis.array,agrid_lon,atol=1e-5) and grtyp != 'Z':
          var.atts.update(grtyp='G',ig1=0,ig2=0,ig3=0,ig4=0)
          continue
        # 'L' grid
        lat0 = float(yaxis.array[0])
        lon0 = float(xaxis.array[0])
        dlat = float(np.mean(yaxis.array[1:] - yaxis.array[:-1]))
        dlon = float(np.mean(xaxis.array[1:] - xaxis.array[:-1]))
        lgrid_key = (lat0,lon0,dlat,dlon)
        if lgrid_key not in lgrid_table:
          lgrid_code = rmn.cxgaig('L',*lgrid_key)
          # Check if 'L' grid has sufficient resolution to capture this
          # set of lat/lon.
          if np.allclose(rmn.cigaxg('L',*lgrid_code),lgrid_key,atol=1e-5):
            # Also check if this contains a repeated longitude.
            # In which case, this is probably GEM output and should be on
            # 'Z' grid?
            if not np.allclose(xaxis.array[0]+360,xaxis.array[-1]):
              lgrid_table[lgrid_key] = lgrid_code
        lgrid_lat = np.linspace(yaxis.array[0],yaxis.array[-1],nj)
        lgrid_lon = np.linspace(xaxis.array[0],xaxis.array[-1],ni)
        if lgrid_key in lgrid_table and np.allclose(lgrid_lat,yaxis.array,atol=1e-5) and np.allclose(lgrid_lon,xaxis.array,atol=1e-5) and grtyp != 'Z':
          ig1, ig2, ig3, ig4 = lgrid_table[(lat0,lon0,dlat,dlon)]
          var.atts.update(grtyp='L',ig1=ig1,ig2=ig2,ig3=ig3,ig4=ig4)
          continue
        # 'Z' grid (with degenerate rotation)
        if lgrid_key not in zlgrid_table:
          grid = rmn.defGrid_ZL(ni,nj,lat0=lat0,lon0=lon0,dlat=dlat,dlon=dlon)
          # Add extra positional records.
          add_coord('>>',1,ni,xaxis.array,typvar='X',etiket='POSX',datyp=5,nbits=32,grtyp='L',ip1=grid['tag1'],ip2=grid['tag2'],ip3=grid['tag3'],ig1=grid['ig1ref'],ig2=grid['ig2ref'],ig3=grid['ig3ref'])
          add_coord('^^',nj,1,yaxis.array,typvar='X',etiket='POSY',datyp=5,nbits=32,grtyp='L',ip1=grid['tag1'],ip2=grid['tag2'],ip3=grid['tag3'],ig1=grid['ig1ref'],ig2=grid['ig2ref'],ig3=grid['ig3ref'])
          zlgrid_table[lgrid_key] = grid
        grid = zlgrid_table[lgrid_key] = grid
        var.atts.update(grtyp='Z',ig1=grid['tag1'],ig2=grid['tag2'],ig3=grid['tag3'])
        continue

      # Case 2: other standard projections.
      projection = var.atts.get('grid_mapping',None)
      projection = projection_table.get(projection,None)
      projection_name = getattr(projection,'atts',{}).get('grid_mapping_name',None)
      # Polar stereographic
      if projection_name == 'polar_stereographic':
        origin = projection.atts.get('latitude_of_projection_origin',None)
        if 'standard_parallel' in projection.atts:
          stdlat = projection.atts['standard_parallel']
        else:
          k = projection.atts['scale_factor_at_projection_origin']
          stdlat = abs(asin(2*k-1)) / pi * 180
        if not np.allclose(stdlat,60.,atol=1e-5):
          warn(_('Standard parallel must be 60 deg to encode polar stereographic projections.  Found %s instead.')%stdlat)
          continue
        if 'straight_vertical_longitude_from_pole' not in projection.atts:
          warn(_('Sterographic projection missing attribute "straight_vertical_longitude_from_pole"'))
          continue
        dgrw = -(projection.atts['straight_vertical_longitude_from_pole']+90)
        dgrw %= 360.0
        d60 = np.mean(xaxis.array[1:]-xaxis.array[:-1])
        east = projection.atts.get('false_easting',0)
        north = projection.atts.get('false_northing',0)
        px = xaxis.array[0] + east
        py = yaxis.array[0] + north
        pi = px / d60 + 1
        pj = py / d60 + 1
        # 'N' grid
        if origin == 90:
          ig1, ig2, ig3, ig4 = rmn.cxgaig('N',pi,pj,d60,dgrw)
        elif origin == -90:
          ig1, ig2, ig3, ig4 = rmn.cxgaig('S',pi,pj,d60,dgrw)
        else:
          warn(_("latitude_of_projection_origin must be 90 or -90.  Found %s")%origin)
          continue
        var.atts.update(grtyp='N',ig1=ig1,ig2=ig2,ig3=ig3,ig4=ig4)
        continue
      # Rotated lat/lon
      if projection_name == 'rotated_latitude_longitude':
        gpole_lat = projection.atts['grid_north_pole_latitude'] * pi / 180
        gpole_lon = projection.atts['grid_north_pole_longitude'] * pi / 180
        npole = projection.atts.get('north_pole_grid_longitude',0.)
        zegrid_key = (gpole_lat,gpole_lon,npole,id(xaxis),id(yaxis))
        if zegrid_key not in zegrid_table:
          xlat1, xlon1, xlat2, xlon2, ax_adjust = get_rotation_params(projection.atts, xaxis.array)
          ax = adjust_ax (xaxis.array, ax_adjust)
          gid = rmn.defGrid_ZEraxes(ax=ax, ay=yaxis.array, xlat1=xlat1, xlon1=xlon1, xlat2=xlat2, xlon2=xlon2)
          zegrid_table[zegrid_key] = gid

        # Set grid descriptors to link to coordinate records.
        if zegrid_key in zegrid_table:
          gid = zegrid_table[zegrid_key]
          var.atts.update(grtyp='Z',ig1=gid['tag1'],ig2=gid['tag2'],ig3=gid['tag3'])
          continue
      # Case 3: yin-yang grid.
      # Find the 2D lat/lon fields associated with this grid.
      lat = None; lon = None
      for coord in var.atts.get('coordinates',[]):
        if coord.name == 'lat': lat = coord
        if coord.name == 'lon': lon = coord
      if lat is not None and lon is not None and nj%2==0 and np.all(yaxis.array[:nj//2]==yaxis.array[nj//2:]):
        yygrid_key = (id(xaxis),id(yaxis),id(lat),id(lon))
        if yygrid_key not in yygrid_table:
          # Fix shape of lat/lon (if had subgrid axes)
          lat = lat.array.reshape(nj,ni)
          lon = lon.array.reshape(nj,ni)
          # Get params for yin grid.
          yin_params = get_rotated_grid_params(xaxis.array,yaxis.array[:nj//2],lat[:nj//2,:],lon[:nj//2,:])
          yang_params = get_rotated_grid_params(xaxis.array,yaxis.array[nj//2:],lat[nj//2:,:],lon[nj//2:,:])
          # Encode yin grid
          xlat1, xlon1, xlat2, xlon2, ax_adjust = get_rotation_params(yin_params, xaxis.array)
          ax = adjust_ax (xaxis.array, ax_adjust)
          yin_gid = rmn.defGrid_ZEraxes(ax=ax, ay=yaxis.array[:nj//2], xlat1=xlat1, xlon1=xlon1, xlat2=xlat2, xlon2=xlon2)
          yinsize=15+ni+nj//2
          # Encode yang grid
          xlat1, xlon1, xlat2, xlon2, ax_adjust = get_rotation_params(yang_params, xaxis.array)
          ax = adjust_ax (xaxis.array, ax_adjust)
          yang_gid = rmn.defGrid_ZEraxes(ax=ax, ay=yaxis.array[nj//2:], xlat1=xlat1, xlon1=xlon1, xlat2=xlat2, xlon2=xlon2)
          # Consruct supergrid (to generate unique tags)
          gid = rmn.ezgdef_supergrid(ni, nj//2, 'U', 'F', 1, (yin_gid['id'],yang_gid['id']))
          gid = rmn.decodeGrid(gid)
          yygrid_table[yygrid_key] = gid
        if yygrid_key in yygrid_table:
          gid = yygrid_table[yygrid_key]
          var.atts.update(grtyp='U',ig1=gid['tag1'],ig2=gid['tag2'],ig3=gid['tag3'])
          continue

    # Add coordinate records to the table.
    for grid in zegrid_table.values():
      add_coord('>>',1,grid['ni'],grid['ax'],typvar='X',etiket='POSX',datyp=5,nbits=32,grtyp='E',ip1=grid['tag1'],ip2=grid['tag2'],ip3=grid['tag3'],ig1=grid['ig1ref'],ig2=grid['ig2ref'],ig3=grid['ig3ref'],ig4=grid['ig4ref'])
      add_coord('^^',grid['nj'],1,grid['ay'],typvar='X',etiket='POSY',datyp=5,nbits=32,grtyp='E',ip1=grid['tag1'],ip2=grid['tag2'],ip3=grid['tag3'],ig1=grid['ig1ref'],ig2=grid['ig2ref'],ig3=grid['ig3ref'],ig4=grid['ig4ref'])
    for grid in yygrid_table.values():
      add_coord('^>',1,25+2*grid['ni']+2*grid['nj'],grid['axy'],typvar='X',etiket='POSXY',datyp=5,nbits=32,grtyp='F',ip1=grid['tag1'],ip2=grid['tag2'],ip3=grid['tag3'],ig1=grid['ig1ref'],ig2=grid['ig2ref'],ig3=grid['ig3ref'],ig4=grid['ig4ref'])


    # Set default grtyp if no better one found.
    for var in self._varlist:
      #TODO: determine appropriate grtyp
      if 'grtyp' not in var.atts:
        var.atts['grtyp'] = 'X'
    super(XYCoords,self)._unmakevars()

    # Fill in columns with default values.
    self._headers['grtyp'] = self._headers['grtyp'].astype('|S1')
    for key in 'ig1','ig2','ig3','ig4':
      if key not in self._headers.keys():
        self._headers[key] = np.ma.masked_all(self._nrecs,dtype='int32')
      if hasattr(self._headers[key],'mask'):
        self._headers[key] = self._headers[key].filled(0)
