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
  4: "sea_ice",        # (lake ice, sea ice))
  5: "all_area_types", # (aggregated)
  6: "urban",          # (urban)
}
# List of "aggregated" variables in rpnphy
# To get the list, from the rpnphy git repository run:
#  git grep 'ON=.*nagg' | sed 's/.*ON=//;s/ .*//'
sfc_agg_nomvars = (
  'AL','BT','FC','B5','FT','FV','HF','H1','IO','J9','Z0',
  'ZT','WT','SD','7A',
)

class Sfc_Codes (BufferBase):

  def _makevars (self):
    from fstd2nc.mixins import _var_type, _axis_type, _dim_type
    from collections import OrderedDict
    import numpy as np

    handled_agg_codes = dict()

    super(Sfc_Codes,self)._makevars()
    for var in self._varlist:

      levels = var.getaxis('level')

      # Look for variables with surface codes.
      if var.atts.get('kind',None) != 3 or levels is None:
        continue

      codes = tuple(levels.array)
      coordinates = var.atts.get('coordinates',[])

      if var.name in sfc_agg_nomvars:
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

      if len(coordinates) > 0:
        var.atts['coordinates'] = coordinates



