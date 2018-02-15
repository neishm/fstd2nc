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

# Codes from "aggregated" variables in rpnphy.
# Equivalent CF labels from:
# http://cfconventions.org/Data/area-type-table/current/build/area-type-table.html
agg_codes = {
  1: "ice_free_land",  # (Soil)
  2: "land_ice",       # (Glacier)
  3: "ice_free_sea",   # (Water)
  4: "sea_ice",        # (Sea Ice)
  5: "all_area_types", # (Aggregated)
  6: "urban",          # (Urban)
}
# List of "aggregated" variables in rpnphy
# To get the list, from the rpnphy git repository run:
#  git grep 'ON=.*nagg' | sed 's/.*ON=//;s/ .*//'
agg_nomvars = (
  'AL','BT','FC','B5','FT','FV','HF','H1','IO','J9','Z0',
  'ZT','WT','SD','7A',
)

# Codes for deep / superficial surface layers in rpnphy.
# Note: Can't find any standard CF encoding for this, feel free to update if
# you know of one!
# Order of codes based on rpnphy src/surface/glaciers.F90
surface_codes = {
  1: "superficial",
  2: "deep",
}
# List of deep/superficial variables in rpnphy
# To get the list, from the rpnphy git repository run:
#  git grep 'ON=.*row\*2 ' | sed 's/.*ON=//;s/ .*//'
surface_nomvars = (
  '2W','I9','I0','9A','I1','5A',
)

class Sfc_Codes (BufferBase):

  def _iter (self):
    from fstd2nc.mixins import _var_type, _modify_axes
    from collections import OrderedDict
    import numpy as np

    handled_agg_codes = []
    handled_surface_codes = []

    for var in super(Sfc_Codes,self)._iter():

      # Look for variables with surface codes.
      if var.atts.get('kind',None) != 3 or 'level' not in var.axes:
        yield var
        continue

      codes = tuple(var.axes['level'])
      coordinates = var.atts.get('coordinates','').split()

      if var.name in agg_nomvars:
        # Change the axis name so the vcoord mixin doesn't look at it.
        var.axes = _modify_axes(var.axes, level='area_type_id')
        # Add the area type to the list of auxiliary coordinates.
        coordinates.append('area_type')
        # Generate the list of surface types.
        if codes not in handled_agg_codes:
          codenames = tuple(agg_codes.get(code,"unknown") for code in codes)
          array = np.array(codenames).view('|S1').reshape(len(codes),-1)
          atts = OrderedDict([('standard_name','area_type')])
          axes = OrderedDict([('area_type_id',codes),('area_type_strlen',tuple(range(array.shape[1])))])
          yield _var_type("area_type",atts,axes,array)
          handled_agg_codes.append(codes)

      if var.name in surface_nomvars:
        # Change the axis name so the vcoord mixin doesn't look at it.
        var.axes = _modify_axes(var.axes, level='surface_id')
        # Add the layer type to the list of auxiliary coordinates.
        coordinates.append('surface')
        # Generate the list of surface types.
        if codes not in handled_surface_codes:
          codenames = tuple(surface_codes.get(code,"unknown") for code in codes)
          array = np.array(codenames).view('|S1').reshape(len(codes),-1)
          atts = OrderedDict([('long_name','surface_layer')])
          axes = OrderedDict([('surface_id',codes),('surface_strlen',tuple(range(array.shape[1])))])
          yield _var_type("surface",atts,axes,array)
          handled_surface_codes.append(codes)

      # Sea ice layers
      if var.name == 'I7':  # Sea ice temperature
        var.axes = _modify_axes(var.axes, level='sea_ice_layer')
        # Don't know any coordinates for these ice layers.

      if len(coordinates) > 0:
        var.atts['coordinates'] = ' '.join(coordinates)

      yield var


