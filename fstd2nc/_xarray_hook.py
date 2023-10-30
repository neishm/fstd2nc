import fstd2nc
from xarray.backends import BackendEntrypoint, BackendArray

# Create a lazy array interface from the given index data.
class FSTDBackendArray(BackendArray):
  __slots__ = ("shape", "dtype", "_buffer", "_files", "_args")
  def __init__ (self, buffer, files, group):
    self._buffer = buffer
    self._files = files
    template = group.variables['template']
    self.shape = template.shape
    self.dtype = template.dtype
  def __getitem__ (self, key):
    raise NotImplementedError("Can't read data yet!")
# Create a dataset from the given index file.
def open_index (buffer, indexfile):
  # Get dimensions / coordinates.
  import xarray as xr
  import netCDF4 as nc
  # Open with xarray to get axes / coordinates.
  ds = xr.open_dataset(indexfile)
  # Open with netCDF4 get get access to the group structures containing the
  # metadata.
  root = nc.Dataset(indexfile, 'r')
  files = root.variables['files']
  # Generate the variables.
  vardict = {}
  for varname, group in root.groups.items():
    # Get a lazy accessor for the data.
    array = FSTDBackendArray (buffer, files, group)
    array = xr.core.indexing.LazilyIndexedArray(array)
    # Define the preferred chunking for the data, based on how it's stored
    # on disk.
    preferred_chunks = {}
    all_dims = group.variables['template'].dimensions
    outer_dims = group.groups['data'].variables['address'].dimensions
    for dim in all_dims:
      if dim in outer_dims:
        preferred_chunks[dim] = 1
      else:
        preferred_chunks[dim] = len(root.dimensions[dim])
    encoding = {'preferred_chunks':preferred_chunks}
    # Annotate the variable with dimensions and other metadata.
    template = group.variables['template']
    vardict[varname] = xr.Variable (template.dimensions, array, attrs=template.__dict__, encoding=encoding)
  # Put these variables into the Dataset structure (where the coordinates are
  # already defined.
  ds = ds.merge(vardict)
  # Remove bookkeeping variables.
  del ds['files']
  del ds.attrs['version']
  return (ds)

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

