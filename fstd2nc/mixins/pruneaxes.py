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
# Mixin for pruning unnecessary axes from the data.

class PruneAxes (BufferBase):
  def _iter (self):
    from fstd2nc.mixins import _var_type
    from collections import OrderedDict
    import numpy as np
    from datetime import timedelta
    from rpnpy.librmn.fstd98 import fstlir

    used_axes = dict()
    all_axes = dict()

    for var in super(PruneAxes,self)._iter():
      # Pull axes out of the stream.
      if isinstance(var,_var_type) and list(var.axes.keys()) == [var.name]:
        all_axes[(var.name,tuple(var.axes[var.name]))] = var
        #yield var
        continue
      # Provide axes only when they are needed by variables.
      for name, values in var.axes.items():
        key = (name,tuple(values))
        if key in all_axes and key not in used_axes:
          axis = all_axes[key]
          yield axis
          used_axes[key] = axis
      yield var
