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
import collections

#################################################
# Provide various external array interfaces for the FSTD data.

# Helper function to embed information about the preferred chunk order.
# Use as a wrapper when constructing dask Array objects.
def _preferred_chunk_order (group, index, array):
  return array

# Helper interface for ordering dask tasks based on FSTD record order.
# Might provide a small speed boost when the OS employs a read-ahead buffer.
# Based on dask.core.get_dependencies
class _RecordOrder (object):
  def __init__(self, dsk):
    self.dask = dsk
  def __call__ (self, arg):
    dsk = self.dask
    work = [arg]
    while work:
        new_work = []
        for w in work:
            typ = type(w)
            if typ is tuple and w and isinstance(w[0], collections.Callable):  # istask(w)
                if w[0] is _preferred_chunk_order:
                  return w[1], w[2]
                else:
                  new_work += w[1:]
            elif typ is list:
                new_work += w
            elif typ is dict:
                new_work += list(w.values())
            else:
                try:
                    if w in dsk:
                        new_work.append(dsk[w])
                except TypeError:  # not hashable
                    pass
        work = new_work

# Add a callback to dask to ensure FSTD records are read in a good order.
try:
  from dask.callbacks import Callback
  class _FSTD_Callback (Callback):
    def _start_state (self, dsk, state):
      ready = sorted(state['ready'][::-1],key=_RecordOrder(dsk))
      state['ready'] = ready[::-1]
  _FSTD_Callback().register()
  del Callback

except ImportError:
  pass


class Extern (BufferBase):

  # The interface for getting chunks into dask.
  def _read_chunk (self, rec_id, shape, dtype, cache={}):
    with self._lock:
      # Cache the last record read from this interface, in case it is
      # immediately requested again.
      # Can happen, for instance, if the same data is sliced in multiple
      # different ways.
      if rec_id not in cache:
        cache.clear()
        cache[rec_id] = self._fstluk(rec_id,dtype=dtype)['d'].transpose().reshape(shape)
      return cache[rec_id]

  def _iter_dask (self):
    """
    Iterate over all the variables, and convert to dask arrays.
    """
    from fstd2nc.mixins import _iter_type, _var_type
    from dask import array as da
    from dask.base import tokenize
    import numpy as np
    from itertools import product
    unique_token = tokenize(self._files,id(self))
    # Make a local copy of two columns from the header table.
    # This speeds up access to their elements in the inner loop.
    all_file_ids = np.array(self._headers['file_id'],copy=True)
    all_keys = np.array(self._headers['key'],copy=True)
    self._makevars()
    for var in self._iter_objects():
      if isinstance(var,_iter_type):
        name = var.name+"-"+unique_token
        ndim = len(var.axes)
        shape = var.shape
        ndim_outer = var.record_id.ndim
        ndim_inner = ndim - ndim_outer
        inner_zeros = (0,)*ndim_inner
        dsk = dict()
        chunk_shape = (1,)*ndim_outer+shape[ndim_outer:]
        for ind in product(*map(range,var.record_id.shape)):
          # Pad index with all dimensions (including inner ones).
          key = (name,) + ind + inner_zeros
          rec_id = int(var.record_id[ind])
          # Add this record as a chunk in the dask Array.
          # Also, specify the preferred order of reading the chunks within the
          # file.
          if rec_id >= 0:
            file_id = all_file_ids[rec_id]
            filename = self._files[file_id]
            rec_key = all_keys[rec_id]
            dsk[key] = (_preferred_chunk_order,filename,rec_key,(self._read_chunk, rec_id, chunk_shape, var.dtype))
          else:
            # Fill missing chunks with fill value or NaN.
            if hasattr(self,'_fill_value'):
              var.atts['_FillValue'] = self._fill_value
              dsk[key] = (np.full, chunk_shape, self._fill_value, var.dtype)
            else:
              dsk[key] = (np.full, chunk_shape, float('nan'), var.dtype)
        chunks = []
        for i in range(ndim_outer):
          chunks.append((1,)*shape[i])
        for i in range(ndim_outer,ndim):
          chunks.append((shape[i],))
        array = da.Array(dsk, name, chunks, var.dtype)
        var = _var_type(var.name,var.atts,var.axes,array)
      yield var

  def to_xarray (self):
    """
    Create an xarray interface for the RPN data.
    Requires the xarray and dask packages.
    """
    from collections import OrderedDict
    import xarray as xr
    out = OrderedDict()
    for var in self._iter_dask():
      if not hasattr(var,'array'): continue
      out[var.name] = xr.DataArray(data=var.array, dims=var.dims, name=var.name, attrs=var.atts)
      # Preserve chunking information for writing to netCDF4.
      if hasattr(var.array,'chunks'):
        chunk_shape = [c[0] for c in var.array.chunks]
        out[var.name].encoding['chunksizes'] = chunk_shape
        out[var.name].encoding['original_shape'] = out[var.name].shape

    # Construct the Dataset from all the variables.
    out = xr.Dataset(out)
    # Decode CF metadata
    out = xr.conventions.decode_cf(out)

    # Make the time dimension unlimited when writing to netCDF.
    out.encoding['unlimited_dims'] = ('time',)

    return out

  def to_iris (self):
    """
    Create an iris interface for the RPN data.
    Requires iris >= 2.0, xarray >= 0.10.3, and dask.
    Returns a CubeList object.
    """
    from iris.cube import CubeList
    out = []
    for var in self.to_xarray().data_vars.values():
      # Omit some problematic variables.
      if var.dtype == '|S1': continue
      # Need to clean up some unrecognized metadata.
      for coord in var.coords.values():
        # Remove units of 'level' (confuses cf_units).
        if coord.attrs.get('units',None) in ('level','sigma_level'):
          coord.attrs.pop('units')
        # Remove non-standard standard names.
        if coord.attrs.get('standard_name',None) == 'atmosphere_hybrid_sigma_ln_pressure_coordinate':
          coord.attrs.pop('standard_name')
      out.append(var.to_iris())
    return CubeList(out)

  def to_pygeode (self):
    """
    Create a pygeode interface for the RPN data.
    Requires pygeode >= 1.2.0, and xarray/dask.
    """
    from pygeode.ext_xarray import from_xarray
    data = self.to_xarray()
    return from_xarray(data)

# Workaround for recent xarray (>0.10.0) which changed the methods in the
# conventions module.
# Fixes an AttributeError when using to_pygeode().
try:
  from xarray.coding import times
  from xarray import conventions
  if not hasattr(conventions,'maybe_encode_datetime'):
    conventions.maybe_encode_datetime = times.CFDatetimeCoder().encode
  if not hasattr(conventions,'maybe_encode_timedelta'):
    conventions.maybe_encode_timedelta = times.CFTimedeltaCoder().encode
except (ImportError,AttributeError):
  pass
