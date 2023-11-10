import fstd2nc
from xarray.backends import BackendEntrypoint, BackendArray

# Wrapper for managing index file.
# Allows the file reference to be pickled.
class IndexFile(object):
  __slots__ = ("_filename", "root")
  def __init__ (self, filename):
    self.__setstate__ (filename)
  def __getstate__ (self):
    return self._filename
  def __setstate__ (self, filename):
    import netCDF4 as nc
    self._filename = filename
    self.root = nc.Dataset(filename,'r')
    self.root.set_auto_mask(False)

# Create a lazy array interface from the given index data.
class FSTDBackendArray(BackendArray):
  __slots__ = ("shape", "dtype", "_buffer", "_name", "_source", "_root", "_group", "_outer_dims","_chunk_shape","_pickles")
  def __init__ (self, buffer, name, source, pickles):
    from pickle import loads
    self._buffer = buffer
    self._name = name
    self._source = source
    self._pickles = pickles
    self._finalize()
  def _finalize (self):
    self._root = self._source.root
    self._group = self._root.groups[self._name]
    template = self._group.variables['template']
    self.shape = template.shape
    self.dtype = template.dtype
    self._outer_dims = len(self._group.groups['data'].variables['address'].dimensions)
    self._chunk_shape = self.shape[self._outer_dims:]
    # Construct pickle objects if this is the first time decoding.
    if len(self._pickles) == 0:
      for varname in self._root.variables:
        if varname.endswith('_pickle'):
          struct = self._root.variables[varname][:]
          struct = [s.view('|S%d'%len(s)) for s in struct]
          struct = [loads(s) for s in struct]
          self._pickles[varname.rstrip('_pickle')] = struct
  def __getstate__ (self):
    state = dict()
    for k in ('_buffer','_name','_source','_pickles'):
      state[k] = getattr(self,k)
    return state
  def __setstate__ (self, state):
    for k, v in state.items():
      setattr(self,k,v)
    self._finalize()
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
    files = self._root.variables['files']
    file_ids = self._group.variables['file_id'][outer_key]
    # Get arguments for data construction.
    args = {}
    # Scalar arguments
    for argname in self._group.variables:
      if argname in ('file_id','template'): continue
      # Check for compound data.
      if argname.endswith('_lookup'):
        ind = self._group.variables[argname][outer_key]
        # For now, assume this type of data is invariant over the variable.
        ind0 = ind.flatten()[0]
        assert np.all(ind==ind0)
        argname = argname.rstrip('_lookup')
        args[argname] = self._pickles[argname][ind0]
        continue
      # Regular scalar data.
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
      filename = files[file_id].flatten()
      filename = filename.view('|S%d'%len(filename))[0].decode()
      current_args = [filename]
      current_args.append(args['data'][0][selection].flatten())
      current_args.append(args['data'][1][selection].flatten())
      for argname, argval in args.items():
        if argname == 'data': continue
        current_args.append(argname)
        # Compound data (single copy over all data)
        if argname in self._pickles:
          current_args.append(argval)
        # Address / length arguments (for data from file)?
        elif isinstance(argval,tuple):
          current_args.append(argval[0][selection].flatten())
          current_args.append(argval[1][selection].flatten())
        # Scalar arguments?
        else:
          # Decode bytes into boolean (from netcdf compatibility).
          if argval.dtype.name == 'uint8':
            argval = argval.astype('bool')
          shape = args['data'][0].shape
          argval = np.array(argval).reshape(shape)[selection].flatten()
          current_args.append(list(argval))
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

