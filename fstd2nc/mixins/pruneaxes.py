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
# Mixin for pruning duplicate axes.

class PruneAxes (BufferBase):
  def _makevars (self):
    from fstd2nc.mixins import _dim_type

    handled_axes = dict()

    super(PruneAxes,self)._makevars()

    # Check for identical axes.
    for axis, varlist in self._iter_axes(varlist=True):
      if isinstance(axis,_dim_type):
        key = (axis.name,len(axis))
      else:
        key = (axis.name,tuple(axis.array))
      # Use only one version of the axis.
      if key in handled_axes:
        for var in varlist:
          var.axes[var.dims.index(axis.name)] = handled_axes[key]
      else:
        handled_axes[key] = axis

