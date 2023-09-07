import fstd2nc
from xarray.backends import BackendEntrypoint

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
  def __call__ (self, filename, **kwargs):
    import fstd2nc
    b = fstd2nc.Buffer.from_xarray(self._obj, **kwargs)
    b.to_fstd(filename)

