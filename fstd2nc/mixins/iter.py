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
from functools import reduce

#################################################
# A user-friendly iterator for using the Buffer in other Python scripts.

# Lightweight wrapper for the data.
# Allows the user to load the data through np.asarray or by slicing it.
class _Array (object):
  # Set some common attributes for the object.
  def __init__ (self, buffer, var):
    from functools import reduce
    self._buffer = buffer
    self._record_id = var.record_id
    # Expected shape and type of the array.
    self.shape = tuple(var.shape)
    self.ndim = len(self.shape)
    self.size = reduce(int.__mul__, self.shape, 1)
    self.dtype = var.dtype
  def __getitem__ (self, key):
    import numpy as np
    # Coerce key into a tuple of slice objects.
    if not isinstance(key,tuple):
      if hasattr(key,'__len__'): key = tuple(key)
      else: key = (key,)
    if len(key) == 1 and hasattr(key[0],'__len__'):
      key = tuple(key[0])
    if Ellipsis in key:
      i = key.index(Ellipsis)
      key = key[:i] + (slice(None),)*(self.ndim-len(key)+1) + key[i+1:]
    key = key + (slice(None),)*(self.ndim-len(key))
    if len(key) > self.ndim:
      raise ValueError(("Too many dimensions for slicing."))
    final_shape = tuple(len(np.arange(n)[k]) for n,k in zip(self.shape,key) if not isinstance(k,int))
    data = np.ma.empty(final_shape, dtype=self.dtype)
    outer_ndim = self._record_id.ndim
    record_id = self._record_id.__getitem__(key[:outer_ndim])
    # Final shape of each record (in case we did any reshaping along ni,nj,nk
    # dimensions).
    record_shape = self.shape[self._record_id.ndim:]
    # Iterate of each record.
    for ind in np.ndindex(record_id.shape):
      r = int(record_id[ind])
      if r >= 0:
        data[ind] = self._buffer._fstluk(r)['d'].transpose().reshape(record_shape)[(Ellipsis,)+key[outer_ndim:]]
      else:
        data.mask = np.ma.getmaskarray(data)
        data.mask[ind] = True
    return data
  def __array__ (self):
    return self.__getitem__(())

class Iter (BufferBase):
  def __iter__ (self):
    """
    Processes the records into multidimensional variables.
    Iterates over (name, atts, axes, array) tuples.
    Note that array may not be a true numpy array (values are not yet loaded
    in memory).  To load the array, pass it to numpy.asarray().

    Deprecated - use .to_xarray() to get a multidimensional structure.
    """
    from warnings import warn
    warn ("Iterating over a Buffer is deprecated.  Use .to_xarray() to access the multidimensional data.", stacklevel=2)
    from fstd2nc.mixins import _iter_type, _var_type
    self._makevars()
    for var in self._iter_objects():
      if isinstance(var, _iter_type):
        array = _Array(self, var)
        var = _var_type (var.name, var.atts, var.axes, array)
      yield var

