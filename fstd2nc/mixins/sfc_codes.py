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

  def _iter (self):
    from fstd2nc.mixins import _var_type, _modify_axes
    from collections import OrderedDict
    import numpy as np

    handled_agg_codes = []
    handled_level_codes = []

    for var in super(Sfc_Codes,self)._iter():

      # Look for variables with surface codes.
      if var.atts.get('kind',None) != 3 or 'level' not in var.axes:
        yield var
        continue

      codes = tuple(var.axes['level'])
      coordinates = var.atts.get('coordinates','').split()

      if var.name in sfc_agg_nomvars:
        # Change the axis name so the vcoord mixin doesn't look at it.
        var.axes = _modify_axes(var.axes, level='sfctype')
        # Add the area type to the list of auxiliary coordinates.
        coordinates.append('surface_type')
        # Generate the list of surface types.
        if codes not in handled_agg_codes:
          codenames = tuple(sfc_agg_codes.get(code,"unknown") for code in codes)
          array = np.array(codenames).view('|S1').reshape(len(codes),-1)
          atts = OrderedDict([('standard_name','area_type')])
          axes = OrderedDict([('sfctype',codes),('sfctype_strlen',tuple(range(array.shape[1])))])
          yield _var_type("surface_type",atts,axes,array)
          handled_agg_codes.append(codes)

      if len(coordinates) > 0:
        var.atts['coordinates'] = ' '.join(coordinates)

      yield var

