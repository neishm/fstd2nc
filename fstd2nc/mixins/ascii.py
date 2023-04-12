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
# Mixin for handling ASCII metadata records.
#
# For example, the META record from GenPhysX.
#

class ASCII (BufferBase):

  def _makevars (self):
    from fstd2nc.mixins import _var_type, _iter_type, _dim_type
    import numpy as np

    super(ASCII,self)._makevars()
    for i,var in enumerate(self._varlist):

      # Decode GenPhysX META record and MLDPn INFO record.
      if isinstance(var,_iter_type) and var.name in ("META","INFO") and var.dtype == np.uint32 and var.record_id.size == 1 and var.atts['nj'] == 1:
        array = self._read_record(var.record_id.flatten()[0]).flatten()
        dim = _dim_type(var.name+'_strlen',array.size)
        array = np.asarray(array,'uint8').view('|S1')
        self._varlist[i] = _var_type(var.name,{},[dim],array)

