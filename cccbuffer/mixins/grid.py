# Mixin for decoding lat/lon coordinates for CCC files.

from fstd2nc.mixins import BufferBase
from fstd2nc.stdout import _, info, warn, error
class Grid(BufferBase):

  # This routine handles metadata and coordinates for the data
  # (after it has been structured into multi-dimensional arrays).
  # Here we add latitude/longitude variables.
  def _makevars (self):
    import numpy as np
    from fstd2nc.mixins import _axis_type, _var_type
    from collections import OrderedDict
    super(Grid,self)._makevars()
    # Annotations for axes
    lat_atts = OrderedDict([('units','degrees_north'),('axis','Y'),('long_name','Latitude'),('standard_name','latitude')])
    lon_atts = OrderedDict([('units','degrees_east'),('axis','X'),('long_name','Longitude'),('standard_name','longitude')])
    # Keep track of lat/lons
    handled_lons = {}
    handled_lats = {}
    for var in self._varlist:

      # Drop degenerate axes?
      for axisname in ('ilat','ilg'):
        axis = var.getaxis(axisname)
        if axis is None: continue
        if len(axis) != 1: continue
        ind = var.dims.index(axisname)
        var.axes = var.axes[:ind] + var.axes[ind+1:]

      # Add latitude
      if var.atts['kind'] == 'GRID' and 'ilat' in var.dims:
        axisname = 'ilat'
      elif var.atts['kind'] == 'ZONL' and 'ilg' in var.dims:
        axisname = 'ilg'
      else:
        continue
      nlat = int(var.atts[axisname])
      khem = int(var.atts['khem'])
      key = (nlat,khem)
      if key not in handled_lats:
        # Use Python-RPN to get Gaussian latitudes.
        # Could replace with another library if more appropriate.
        import rpnpy.librmn.all as rmn
        gid = rmn.defGrid_G (10, nlat, glb=(khem==0), north=(khem==1))['id']
        handled_lats[key] = _axis_type('lat',lat_atts.copy(),rmn.gdll(gid)['lat'][0,:])
      lat = handled_lats[key]
      var.axes[var.dims.index(axisname)] = lat

      # Add longtiude
      if var.atts['kind'] != 'GRID': continue
      if 'ilg' in var.dims:
        nlon = int(var.atts['ilg'])
        key = int(nlon)
        if key not in handled_lons:
          handled_lons[key] = _axis_type('lon',lon_atts.copy(),np.linspace(0,360,nlon,dtype='float32'))
        lon = handled_lons[key]
        var.axes[var.dims.index('ilg')] = lon

