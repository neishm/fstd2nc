import fstd2nc
from xarray.backends import BackendEntrypoint, BackendArray

# Create a lazy array interface from the given index data.
class FSTDBackendArray(BackendArray):
  __slots__ = ("shape", "dtype", "_buffer", "_files", "_group", "_outer_dims","_chunk_shape")
  def __init__ (self, buffer, files, group):
    self._buffer = buffer
    self._files = files
    self._group = group
    template = group.variables['template']
    self.shape = template.shape
    self.dtype = template.dtype
    self._outer_dims = len(group.groups['data'].variables['address'].dimensions)
    self._chunk_shape = self.shape[self._outer_dims:]
  def __getitem__ (self, key):
    import xarray as xr
    return xr.core.indexing.explicit_indexing_adapter (
      key,
      self.shape,
      xr.core.indexing.IndexingSupport.OUTER,
      self._raw_indexing_method
    )
  def _raw_indexing_method (self, key):
    import numpy as np
    outer_key = key[:self._outer_dims]
    inner_key = key[self._outer_dims:]
    # Get the telemetry of the data requested.
    file_ids = self._group.variables['file_id'][outer_key]
    # Get arguments for data construction.
    args = {}
    # Scalar arguments
    for argname in self._group.variables:
      if argname in ('file_id','template'): continue
      args[argname] = self._group.variables[argname][outer_key]
    # Address / length arguments.
    for argname in self._group.groups:
      address = self._group.groups[argname].variables['address'][outer_key]
      length = self._group.groups[argname].variables['length'][outer_key]
      args[argname] = (address,length)
    # Loop over each file, build up the output.
    # Will use the existing dask interface for fetching the data, since it's
    # similar enough in terms of input arguments.
    out = None
    unique_file_ids = np.unique(file_ids)
    for file_id in unique_file_ids:
      selection = (file_ids == file_id)
      nrecs = np.sum(selection)
      # Set up the arguments going into the data loader.
      filename = self._files[file_id].flatten()
      filename = filename.view('|S%d'%len(filename))[0].decode()
      current_args = [filename]
      current_args.append(args['data'][0].flatten())
      current_args.append(args['data'][1].flatten())
      for argname, argval in args.items():
        if argname == 'data': continue
        current_args.append(argname)
        if isinstance(argval,tuple):
          current_args.append(argval[0])
          current_args.append(argval[1])
        else:
          current_args.append(argval)
      current_args.append((nrecs,)+self._chunk_shape)
      # Fetch the data.
      # Vectorized over multiple records per input file.
      data = self._buffer._dasked_read(*current_args)
      # Apply the inner slicing.
      data = data[(slice(None),)+inner_key]
      # Build up the output array.
      if out is None:
        out = np.empty(file_ids.shape+data.shape[1:], dtype=data.dtype)
      out[selection] = data
      return out

# Quick and dirty hack to get list of parameters that could be passed to
# to Buffer interface.
# Use the command-line arguments, converted to Python arguments.
# This is currently the only definitive list of arguments.
def _get_open_dataset_parameters ():
  import argparse
  parser = argparse.ArgumentParser()
  fstd2nc.Buffer._cmdline_args(parser)
  args = parser.parse_args([])
  return ('filename_or_obj', 'drop_variables', 'fused') + tuple(vars(args).keys())

# Register this as a backend for xarray, so FST files can be directly
# opened with xarray.open_dataset
class FSTDBackendEntrypoint(BackendEntrypoint):
  description = "Open FST files in xarray."
  open_dataset_parameters = _get_open_dataset_parameters()
  def open_dataset (self, filename_or_obj, drop_variables=None, fused=False, **kwargs):
    import fstd2nc
    if drop_variables is not None: kwargs['exclude'] = drop_variables
    return fstd2nc.Buffer(filename_or_obj, **kwargs).to_xarray(fused=fused)

  def guess_can_open (self, filename_or_obj):
    from fstd2nc.mixins import _expand_files
    try:
      infiles = _expand_files(filename_or_obj)
      infiles = list(zip(*infiles))[1]
      if len(infiles) == 0: return False
      # Check first matching file.
      with open(infiles[0],'rb') as f:
        magic = f.read(16)
      return len(magic) >= 16 and magic[12:] == b'STDR'
    except Exception:
      # If something unexpected happened, then assume it's not a valid FST file.
      return False

# Add a "to_fstd" shortcut to Dataset objects.
import xarray as xr
@xr.register_dataset_accessor("to_fstd")
class to_fstd_accessor:
  def __init__ (self, xarray_obj):
    self._obj = xarray_obj
  def __call__ (self, filename, append=False, **kwargs):
    import fstd2nc
    b = fstd2nc.Buffer.from_xarray(self._obj, **kwargs)
    b.to_fstd(filename, append=append)

