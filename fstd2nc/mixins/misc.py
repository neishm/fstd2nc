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
# Miscellaneous mixins that don't fit anywhere else.


#################################################
# Remove extraneous dimensions from the output.

class NoNK (BufferBase):
  def __init__ (self, *args, **kwargs):
    import numpy as np
    super(NoNK,self).__init__(*args, **kwargs)
    if 'nk' in self._headers:
      if np.all(self._headers['nk'] == 1):
        del self._headers['nk']
  def _makevars (self):
    super(NoNK,self)._makevars()
    for var in self._varlist:
      if 'k' in var.dims and len(var.getaxis('k')) == 1:
        var.axes.pop(var.dims.index('k'))
      if 'j' in var.dims and len(var.getaxis('j')) == 1:
        var.axes.pop(var.dims.index('j'))

