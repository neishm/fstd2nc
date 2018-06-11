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
# Provide an xarray+dask interface for the FSTD data.

# Helper interface for ordering tasks based on FSTD record order.
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
            if typ is tuple and w and callable(w[0]):  # istask(w)
                if getattr(w[0],'__name__',None) == '_read_chunk':
                  return w[1]
                else:
                  new_work += w[1:]
            elif typ is list:
                new_work += w
            elif typ is dict:
                new_work += w.values()
            else:
                try:
                    if w in dsk:
                        work.append(dsk[w])
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


class XArray (BufferBase):

  # The interface for getting chunks into dask.
  def _read_chunk (self, rec_id, shape, dtype):
    return self._fstluk(rec_id,dtype=dtype)['d'].transpose().reshape(shape)

  def to_xarray (self):
    """
    Create an xarray interface for the RPN data.
    Requires the xarray and dask packages.
    """
    from fstd2nc.mixins import _iter_type, _var_type
    import xarray as xr
    from dask import array as da
    from dask.base import tokenize
    import numpy as np
    unique_token = tokenize(self._files,self._headers)
    out = dict()
    for var in self._iter():
      if isinstance(var,_var_type):
        array = var.array
      elif isinstance(var,_iter_type):
        name = var.name+"-"+unique_token
        ndim = len(var.axes)
        shape = tuple(map(len,var.axes.values()))
        ndim_outer = var.record_id.ndim
        ndim_inner = ndim - ndim_outer
        dsk = dict()
        for ind in np.ndindex(var.record_id.shape):
          # Pad index with all dimensions (including inner ones).
          key = (name,) + ind + (0,)*ndim_inner
          chunk_shape = (1,)*ndim_outer+shape[ndim_outer:]
          rec_id = var.record_id[ind]
          if rec_id >= 0:
            dsk[key] = (self._read_chunk, rec_id, chunk_shape, var.dtype)
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
      else:
        warn(_("Unhandled type %s - ignoring variable.")%type(var))
        continue

      out[var.name] = xr.DataArray(data=array, dims=tuple(var.axes.keys()), name=var.name, attrs=var.atts)

    # Construct the Dataset from all the variables.
    out = xr.Dataset(out)
    # Decode CF metadata
    out = xr.conventions.decode_cf(out)

    # Make the time dimension unlimited when writing to netCDF.
    out.encoding['unlimited_dims'] = ('time',)

    return out

