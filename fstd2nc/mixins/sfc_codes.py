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


#################################################
# Mixin for handling various surface field codes.
#
# Converts "arbitrary code" levels from RPN files into surface types, when
# they are known.
#
# Note: due to the limited information available in the RPN files, the codes
# must be determined from a hard-coded list of variable names.

# Codes from "aggregated" surface variables in rpnphy.
# Equivalent CF labels from:
# http://cfconventions.org/Data/area-type-table/current/build/area-type-table.html
sfc_agg_codes = {
  1: "ice_free_land",  # (soil)
  2: "land_ice",       # (glaciers, ice caps)
  3: "ice_free_sea",   # (lake water, sea water)
  4: "lake_ice_or_sea_ice", # (lake ice, sea ice))
  5: "all_area_types", # (aggregated)
  6: "urban",          # (urban)
}
reverse_sfc_agg_codes = dict((v,k) for k,v in sfc_agg_codes.items())
# List of "aggregated" variables in rpnphy
# To get the list, from the rpnphy git repository run:
#  git grep 'ON=.*nagg' | sed 's/.*ON=//;s/ .*//'
sfc_agg_nomvars = (
  'AL','BT','FC','B5','FT','FV','HF','H1','IO','J9','Z0',
  'ZT','WT','SD','7A','TRAF',
  'DQST','DQSZ','TDST','TDSZ','TJST','TJSZ','UDST','UDSZ',
  'VDST','VDSZ','RTSU','RTHD','DXSU','DXHD','GXSU','GXHD',
  'GTSU','GTHD','WBT','QSSU','QSSK','QLSK','QSRD','QLRD',
  'SWE',
)

# List of variables that use a soil depth coordinate.
# For output of SPS.
default_soil_depths = "0.05,0.1,0.2,0.4,1.0,2.0,3.0"
soil_depth_nomvars = ('WSOL', 'ISOL')

class Sfc_Codes (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    super(Sfc_Codes,cls)._cmdline_args(parser)
    parser.add_argument('--sfc-agg-vars', metavar=_('NAME,NAME,...'), help=_('Define additional surface aggregate fields.'))
    parser.add_argument('--soil-depths', default=default_soil_depths, help=_('Define custom depths for soil fields (%s).  Defaults are %%(default)s.')%(','.join(soil_depth_nomvars)))

  def __init__ (self, *args, **kwargs):
    """
    sfc_agg_vars : str or list, optional
        Define additional surface aggregate fields.
    soil_depths : str or list, optional
        Define custom depths for soil fields.
    """
    import numpy as np
    sfc_agg_vars = kwargs.pop('sfc_agg_vars',None)
    if sfc_agg_vars is None:
      sfc_agg_vars = []
    if isinstance(sfc_agg_vars,str):
      sfc_agg_vars = sfc_agg_vars.replace(',', ' ')
      sfc_agg_vars = sfc_agg_vars.split()
    self._sfc_agg_nomvars = sfc_agg_nomvars + tuple(sfc_agg_vars)
    soil_depths = kwargs.pop('soil_depths',default_soil_depths)
    try:
      if isinstance(soil_depths,str):
        soil_depths = map(float,soil_depths.split(','))
      self._soil_bounds = np.array((0.,) + tuple(soil_depths), dtype='float32')
    except ValueError:
      error(_("Unable to parse soil-depths parameter."))
    super(Sfc_Codes,self).__init__(*args,**kwargs)

  def _makevars (self):
    from fstd2nc.mixins import _var_type, _axis_type, _dim_type
    from collections import OrderedDict
    import numpy as np

    handled_agg_codes = dict()
    handled_soil_depth_codes = dict()

    super(Sfc_Codes,self)._makevars()
    for var in self._varlist:

      levels = var.getaxis('level')

      # Look for variables with surface codes.
      if var.atts.get('kind',None) != 3 or levels is None:
        continue

      codes = tuple(levels.array)
      coordinates = []

      # Handle surface type codes.
      if var.name in self._sfc_agg_nomvars:
        # Generate the list of surface types.
        if codes not in handled_agg_codes:
          codenames = tuple(sfc_agg_codes.get(code,"unknown") for code in codes)
          array = np.array(codenames,dtype=np.char.string_).view('|S1').reshape(len(codes),-1)
          atts = OrderedDict([('standard_name','area_type')])
          sfctype = _dim_type('sfctype',array.shape[0])
          sfctype_strlen = _dim_type('sfctype_strlen',array.shape[1])
          surface_type = _var_type("surface_type",atts,[sfctype,sfctype_strlen],array)
          handled_agg_codes[codes] = surface_type
        surface_type = handled_agg_codes[codes]
        # "Levels" are actually surface type ids.
        var.axes[var.dims.index('level')] = surface_type.getaxis('sfctype')
        # Add the area type to the list of auxiliary coordinates.
        coordinates.append(surface_type)

      # Handle soil depth codes.
      if var.name in soil_depth_nomvars and set(codes) <= set(range(1,len(self._soil_bounds))):
        # Generate the soil depths.
        if codes not in handled_soil_depth_codes:
          indices = np.asarray(codes,dtype='int32')
          depths = (self._soil_bounds[indices-1] + self._soil_bounds[indices]) / 2.
          atts = OrderedDict([('standard_name','depth'), ('long_name','soil depth'), ('axis', 'Z'), ('positive', 'down'), ('units', 'm')])
          soil_depth = _axis_type('soil_depth',atts,depths)
          bnds = _dim_type('bnds',2)
          atts = OrderedDict([('units','m'),('long_name','soil depth bounds')])
          soil_depth_bnds = _var_type('soil_depth_bnds',atts,[soil_depth,bnds],np.stack([self._soil_bounds[indices-1],self._soil_bounds[indices]]).T)
          soil_depth.atts['bounds'] = soil_depth_bnds
          handled_soil_depth_codes[codes] = soil_depth_bnds
        soil_depth_bnds = handled_soil_depth_codes[codes]
        var.axes[var.dims.index('level')] = soil_depth_bnds.getaxis('soil_depth')
      elif var.name in soil_depth_nomvars:
        warn(_("More than the expected number of soil depths were found.  No depth values will be encoded."))
        soil_depth = _dim_type('soil_depth',len(codes))
        var.axes[var.dims.index('level')] = soil_depth
      if len(coordinates) > 0:
        var.deps.extend(coordinates)


  def _unmakevars (self):
    from fstd2nc.mixins import _axis_type
    from fstd2nc.mixins.vcoords import encode_ip1
    import numpy as np

    # Re-encode surface types back to ip1 values.
    levels = {}
    for var in self._varlist:
      if 'sfctype' in var.dims:
        sfctype = var.getaxis('sfctype')
        if id(sfctype) not in levels:
          # Find the surface types.
          surface_type = None
          for v in self._varlist:
            if v.name != 'surface_type': continue
            if sfctype in v.axes:
              surface_type = v
              break
          if surface_type is not None:
            codes = [reverse_sfc_agg_codes[t.decode()] for t in surface_type.array]
            codes = np.array(codes,dtype='float32')
            kind = np.empty(len(codes),dtype='int32')
            kind[:] = 3
            ip1 = encode_ip1(kind,codes)
            #TODO: let vcoords mixin do the encoding, if it can be run after
            # this mixin.
            levels[id(sfctype)] = _axis_type('ip1',{},ip1)
        if id(sfctype) in levels:
          var.axes[var.dims.index('sfctype')] = levels[id(sfctype)]

    super(Sfc_Codes,self)._unmakevars()
