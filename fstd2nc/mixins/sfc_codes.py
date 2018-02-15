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

# Codes from "aggregated" variables in rpnphy
# Equivalent CF labels from:
# http://cfconventions.org/Data/area-type-table/current/build/area-type-table.html
agg_codes = {
  1: "all_area_types", # (Aggregated)
  2: "sea_ice",        # (Sea Ice)
  3: "ice_free_sea",   # (Water)
  4: "land_ice",       # (Glacier)
  5: "ice_free_land",  # (Soil)
}
# List of "aggregated" variables in rpnphy
# To get the list, from the rpnphy git repository run:
#  git grep 'ON=.*nagg' | sed 's/.*ON=//;s/ .*//'
agg_nomvars = (
  'AL','BT','FC','B5','FT','FV','HF','H1','IO','J9','Z0',
  'ZT','WT','SD','7A',
)

class Sfc_Codes (BufferBase):

  def _iter (self):
    from fstd2nc.mixins import _iter_type, _var_type, _modify_axes
    from collections import OrderedDict
    import numpy as np
    from datetime import timedelta
    from rpnpy.librmn.fstd98 import fstlir

    handled_agg_codes = []

    for var in super(Sfc_Codes,self)._iter():

      # Look for variables with surface codes:
      if var.name in agg_nomvars and var.atts.get('kind',None) == 3 and 'level' in var.axes:
        codes = tuple(var.axes['level'])
        # Change the axis name so the vcoord mixin doesn't look at it.
        var.axes = _modify_axes(var.axes, level='area_type_id')
        # Add the area type to the list of auxiliary coordinates.
        coordinates = var.atts.get('coordinates','').split()
        coordinates.append('area_type')
        var.atts['coordinates'] = ' '.join(coordinates)
        # Generate the list of surface types.
        if codes not in handled_agg_codes:
          codenames = tuple(agg_codes.get(code,"unknown") for code in codes)
          array = np.array(codenames).view('|S1').reshape(len(codes),-1)
          atts = OrderedDict([('standard_name','area_type')])
          axes = OrderedDict([('area_type_id',codes),('area_type_strlen',tuple(range(array.shape[1])))])
          yield _var_type("area_type",atts,axes,array)
          handled_agg_codes.append(codes)

      yield var


