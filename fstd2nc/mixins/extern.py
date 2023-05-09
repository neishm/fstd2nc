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
try:
  from collections import Callable
except ImportError:  # Python 3.10
  from collections.abc import Callable

#################################################
# Provide various external array interfaces for the FSTD data.

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
            if typ is tuple and w and isinstance(w[0], Callable):  # istask(w)
                if w[0] is _read_block:
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
      # Try sorting by FSTD filename, offset (if applicable)
      try:
        ready = sorted(state['ready'][::-1],key=_RecordOrder(dsk))
        state['ready'] = ready[::-1]
      except TypeError: pass  # Not applicable for this graph.
  del Callback
except ImportError:
  pass

# Only turn on this callback if explicitly requested by the user.
# It may cause problems when there are large graphs being communicated
# to a dask cluster environment.
def force_strict_ordering ():
  _FSTD_Callback().register()

# Method for reading a block from a file.
def _read_block (filename, offset, length):
  import numpy as np
  # Scalar version first.
  if not hasattr(offset,'__len__'):
    with open(filename,'rb') as f:
      f.seek (offset,0)
      return np.fromfile(f,'B',length)
  # Vectorized version.
  out = []
  with open(filename,'rb') as f:
    for o, l in zip(offset, length):
      f.seek (o,0)
      out.append(np.fromfile(f,'B',l))
  return out


class ExternOutput (BufferBase):

  # Helper method to call decode routine using the given inputs.
  @classmethod
  def _dasked_decode (cls, *args):
    keys = args[:-1:2]
    values = args[1::2]
    # Scalar case (decoding single record)
    if not any(isinstance(v,list) for v in values):
      kwargs = dict(zip(keys,values))
      return cls._decode(**kwargs)
    # Vectorized case (decoding multiple records into a single fused array)
    nrec = [len(v) for v in values if isinstance(v,list)][0]
    values = [v if isinstance(v,list) else [v]*nrec for v in values]
    out = []
    for i in range(nrec):
      kwargs = dict((k,v[i]) for k,v in zip(keys,values))
      out.append(cls._decode(**kwargs))
    import numpy as np
    return np.concatenate(out)

  def _iter_dask (self, include_coords=True, fused=True):
    """
    Iterate over all the variables, and convert to dask arrays.
    """
    from fstd2nc.mixins import _iter_type, _chunk_type, _var_type
    from dask import array as da
    from dask.base import tokenize
    import numpy as np
    from itertools import product
    unique_token = tokenize(self._files,id(self))
    files = np.array(self._files, dtype=object)
    self._makevars()
    for var in self._iter_objects():
      if not include_coords:
        if var not in self._varlist:
          continue
      if not isinstance(var,(_iter_type,_chunk_type)):
        yield var
        continue
      name = var.name+"-"+unique_token
      ndim = len(var.axes)
      shape = var.shape
      # Convert _iter_type to more generic _chunk_type.
      if isinstance(var,_iter_type):
        ndim_outer = var.record_id.ndim
        ndim_inner = ndim - ndim_outer
        chunks = [(1,)*n for n in shape[:ndim_outer]] + [(n,) for n in shape[ndim_outer:]]
        record_id = var.record_id.reshape(shape[:ndim_outer] + (1,)*ndim_inner)
        var = _chunk_type (var.name, var.atts, var.axes, var.dtype, chunks, record_id)

      # Fuse the records to make larger (and more dask-friendly) chunks.
      # Disable fusing for dask-based arrays.
      # Could make it work, but would need more explicit looping over
      # records to see where dask arrays are available (slower?)
      if 'd' in self._headers and any(d is not None for d in self._headers['d'][var.record_id].flatten()):
        pass # Shut off fusion.
      elif fused:
        file_ids = self._headers['file_id'][var.record_id]
        # Find dimension to fuse.
        dim = var.record_id.ndim-1
        while dim >= 0 and len(var.chunks[dim]) == 1:
          dim -= 1
        # Check if there is something available to fuse.
        if dim >= 0:
          # Look at part of record to guess at number of records to chunk.
          # Can only chunk within a file.
          sample = file_ids[(0,)*dim].squeeze()
          dx = sum(sample==sample[0])
          # Check if this guessed chunk size works for the data shape.
          if var.record_id.shape[dim] % dx == 0:
            if dx > 1:
              dy = var.record_id.shape[dim] // dx
              # Check if we always stay within file bounds.
              shape = var.record_id.shape[:dim] + (dy, dx)
              check = file_ids.reshape(shape)
              if np.all(check[...,:] == check[...,0:1]):
                chunks = list(var.chunks)
                chunks[dim] = (dx,)*dy
                shape = var.record_id.shape[:dim] + (dy,) + var.record_id.shape[dim+1:] + (dx,)
                # Put chunks record ids in extra dimension at end.
                record_id = var.record_id.reshape(shape)
                var = _chunk_type (var.name, var.atts, var.axes, var.dtype, chunks, record_id)
              else:
                warn(_("Unable to fuse some variables."))
          else:
            warn(_("Unable to fuse some variables."))

      # Transform _chunk_type data from record indices to graphs.
      # For fused data, use 2D record_id array (chunk, ids).
      # For unfused data, use 1D record_id array.
      if var.record_id.ndim > len(var.chunks):
        record_id = var.record_id.reshape(-1,var.record_id.shape[-1])
      else:
        record_id = var.record_id.flatten()
      file_ids = self._headers['file_id'][record_id]
      file_ids[record_id<0] = -1
      # Assume same file within chunk.
      if file_ids.ndim > 1: file_ids = file_ids[:,0]
      nchunks = record_id.shape[0]
      args = []
      # Add data sources (primary data plus possible secondary like mask).
      for key, (addr_col, len_col, dname) in self._decoder_data:
        if addr_col not in self._headers: continue # Skip inactive arguments.
        fname = files[file_ids]
        addr = self._headers[addr_col][record_id]
        if addr.ndim > 1: addr = map(list,addr)
        length = self._headers[len_col][record_id]
        if length.ndim > 1: length = map(list,length)
        rb = zip([_read_block]*nchunks, fname, addr, length)
        # Skip the _read_blocks construct where we don't have file data.
        # (e.g. for missing records or where data coming from dask).
        rb = [_rb if fid >= 0 else None for _rb,fid in zip(rb,file_ids)]
        # Inject dask array data where it's available.
        if dname in self._headers:
          rb = [d if d is not None else _rb for d,_rb in zip(self._headers[dname][record_id].flatten(),rb)]
        args.extend([[key]*nchunks, rb])
      # Add extra arguments from columns.
      for key in self._decoder_extra_args:
        if key not in self._headers: continue  # Skip inactive arguments.
        values = self._headers[key][record_id]
        if values.ndim > 1:
          values = map(list,values)
          # For case where all values are identical, reduce to scalar.
          values = [v[0] if len(set(v)) == 1 else v for v in values]
        args.extend([[key]*nchunks, values])
      # Add scalar arguments.
      for key, value in self._decoder_scalar_args().items():
        args.extend([[key]*nchunks, [value]*nchunks])
      graphs = zip([self._dasked_decode]*nchunks, *args)
      array = np.empty(nchunks, dtype=object)
      array[:] = list(graphs)
      shape = tuple(len(c) for c in var.chunks)
      array = array.reshape(shape)

      # Create dask array from this info.
      chunk_indices = [np.cumsum((0,)+c)[:-1] for c in var.chunks]
      # Loop over all indices, generate dask graph.
      dsk = dict()
      for ind, chunk_coord, chunk_shape in zip(np.ndindex(array.shape), product(*chunk_indices), product(*var.chunks)):
        # Unique key for this graph member.
        key = (name,) + chunk_coord
        # Add this record as a chunk in the dask Array.
        graph = array[ind]
        if graph is not None:
          dsk[key] = (np.reshape, graph, chunk_shape)
        else:
          # Fill missing chunks with fill value or NaN.
          if hasattr(self,'_fill_value'):
            var.atts['_FillValue'] = self._fill_value
            dsk[key] = (np.full, chunk_shape, self._fill_value, var.dtype)
          else:
            dsk[key] = (np.full, chunk_shape, float('nan'), var.dtype)
      array = da.Array(dsk, name, var.chunks, var.dtype)
      var = _var_type(var.name,var.atts,var.axes,array)
      yield var

  def to_xarray (self, fused=True):
    """
    Create an xarray interface for the RPN data.
    Requires the xarray and dask packages.
    """
    from collections import OrderedDict
    import xarray as xr
    out = OrderedDict()
    for var in self._iter_dask(fused=fused):
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

    # Remove coordinates from encoding, so that xarray can determine the
    # appropriate values upon writing to netCDF4.
    for var in out.variables:
      out[var].encoding.pop('coordinates',None)

    return out

  def to_xarray_list (self, fused=True):
    """
    Similar to the to_xarray method, but returns a list of xarray Datasets,
    one for each variable, instead of a single Dataset object.
    Could be useful for case where the dimension names are non-unique, to avoid
    name clobbering (in conjunction with unique_names initialization option).
    E.g.,
    Buffer("filename",unique_names=False).to_xarray_list()
    """
    from collections import OrderedDict
    import xarray as xr
    from fstd2nc.mixins import _axis_type
    out_list = []
    for var in self._iter_dask(include_coords=False, fused=fused):
      if not hasattr(var,'array'): continue
      out = OrderedDict()
      out[var.name] = xr.DataArray(data=var.array, dims=var.dims, name=var.name, attrs=var.atts)
      for extra in self._iter_objects(var):
        if not hasattr(extra,'array'): continue
        out[extra.name] = xr.DataArray(data=extra.array, dims=extra.dims, name=extra.name, attrs=extra.atts)
      # Preserve chunking information for writing to netCDF4.
      if hasattr(var.array,'chunks'):
        chunk_shape = [c[0] for c in var.array.chunks]
        out[var.name].encoding['chunksizes'] = chunk_shape
        out[var.name].encoding['original_shape'] = out[var.name].shape
      # Construct the Dataset from the variable.
      out = xr.Dataset(out)
      # Decode CF metadata
      out = xr.conventions.decode_cf(out)
      # Make the time dimension unlimited when writing to netCDF.
      out.encoding['unlimited_dims'] = ('time',)
      out_list.append(out)

    return out_list

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
    _fix_to_pygeode()
    from pygeode.ext_xarray import from_xarray
    data = self.to_xarray()
    return from_xarray(data)

  def to_fstpy (self):
    """
    Create a table compatible with the fstpy module.
    Requires pandas and dask.
    """
    import pandas as pd
    import numpy as np
    from fstd2nc.extra import decode
    from dask import delayed
    import dask.array as da
    # Put all the header info into a dictionary.
    fields = ['nomvar', 'typvar', 'etiket', 'ni', 'nj', 'nk', 'dateo', 'ip1', 'ip2', 'ip3', 'deet', 'npas', 'datyp', 'nbits', 'grtyp', 'ig1', 'ig2', 'ig3', 'ig4', 'datev']
    table = dict()
    # Create a mask to exclude deleted / overwritten / unselected records.
    # Include all meta (coordinate) records in the output.
    mask = self._headers['selected'] | self._headers['ismeta']
    for field in fields:
      col = self._headers[field][mask]
      # Convert byte arrays to strings, which is what fstpy expects.
      if col.dtype.str.startswith('|S'):
        col = np.asarray(col,dtype=col.dtype.str.replace('|S','<U'))
      table[field] = col
    # Convert to pandas table.
    table = pd.DataFrame.from_dict(table)
    # Generate dask objects.
    # NOTE: using raw (records for this, no transformations / masking applied).
    if 'file_id' in self._headers and 'address' in self._headers and 'length' in self._headers:
      filenames = np.array(self._files+[None])
      filenames = filenames[self._headers['file_id']]
      arrays = map(delayed(_read_block), filenames, self._headers['address'], self._headers['length'])
      arrays = map(delayed(decode), arrays)
      shape = zip(self._headers['nj'], self._headers['ni'])
      arrays = map(delayed(np.reshape), arrays, shape)
      arrays = map(delayed(np.transpose), arrays)
      shape = zip(self._headers['ni'], self._headers['nj'])
      arrays = map(da.from_delayed, arrays, shape, self._headers['dtype'])
      arrays = list(arrays)
      d = np.empty(self._nrecs, dtype=object)
      isvalid = self._headers['file_id'] >= 0
      for i in range(self._nrecs):
        if isvalid[i]:
          d[i] = arrays[i]
    else:
      d = np.empty(self._nrecs, dtype=object)
    # Check if we already have dask arrays to use.
    if 'd' in self._headers:
      isvalid = self._headers['d'] != None
      for i in range(self._nrecs):
        if isvalid[i]:
          d[i] = self._headers['d'][i]
    table['d'] = d[mask]
    return table

# Workaround for recent xarray (>0.10.0) which changed the methods in the
# conventions module.
# Fixes an AttributeError when using to_pygeode().
def _fix_to_pygeode (fixed=[False]):
  if fixed[0] is True: return
  try:
    from xarray.coding import times
    from xarray import conventions
    if not hasattr(conventions,'maybe_encode_datetime'):
      conventions.maybe_encode_datetime = times.CFDatetimeCoder().encode
    if not hasattr(conventions,'maybe_encode_timedelta'):
      conventions.maybe_encode_timedelta = times.CFTimedeltaCoder().encode
  except (ImportError,AttributeError):
    pass
  fixed[0] = True

class ExternInput (BufferBase):
  @classmethod
  def from_fstpy (cls, table, **kwargs):
    import numpy as np
    if hasattr(table,'to_pandas'):
      table = table.to_pandas()
    # Construct the record header info from the table.
    fields = ['nomvar', 'typvar', 'etiket', 'ni', 'nj', 'nk', 'dateo', 'ip1', 'ip2', 'ip3', 'deet', 'npas', 'datyp', 'nbits', 'grtyp', 'ig1', 'ig2', 'ig3', 'ig4', 'datev']
    headers = {}
    for col in fields:
      headers[col] = table[col].values.copy()
    # Pad out string variables with spaces.
    headers['nomvar'] = np.asarray(headers['nomvar'], dtype='|S4')
    headers['typvar'] = np.asarray(headers['nomvar'], dtype='|S2')
    headers['etiket'] = np.asarray(headers['etiket'], dtype='|S12')
    headers['grtyp'] = np.asarray(headers['grtyp'], dtype='|S1')
    headers['nomvar'] = np.char.ljust(headers['nomvar'], 4, ' ')
    headers['typvar'] = np.char.ljust(headers['typvar'], 2, ' ')
    headers['etiket'] = np.char.ljust(headers['etiket'], 12, ' ')
    # Add other fields that may be needed.
    if 'dltf' not in headers:
      headers['dltf'] = np.zeros(len(headers['nomvar']), dtype='int32')
      # We don't have file address info, so mark this as 'None' in case
      # any subroutine is looking for this info.
      headers['address'] = np.empty(len(headers['nomvar']), dtype=object)
      headers['length'] = np.empty(len(headers['nomvar']), dtype=object)
    # Fake file id (just so netcdf mixin _quick_load function doesn't crash).
    headers['file_id'] = np.empty(len(headers['nomvar']), dtype='int32')
    headers['file_id'][:] = -1
    # Add in data.
    headers['d'] = np.empty(len(headers['nomvar']), dtype=object)
    headers['d'][:] = table['d']

    # Encapsulate this info in a structure.
    fake_buffer = cls.__new__(cls)
    fake_buffer._files = [None]
    fake_buffer._headers = headers

    # Initialize a Buffer object with this info.
    b = cls(fake_buffer, **kwargs)

    return b

