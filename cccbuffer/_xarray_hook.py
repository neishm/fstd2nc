import cccbuffer
from xarray.backends import BackendEntrypoint

# Quick and dirty hack to get list of parameters that could be passed to
# to Buffer interface.
# Use the command-line arguments, converted to Python arguments.
# This is currently the only definitive list of arguments.
def _get_open_dataset_parameters ():
  import argparse
  parser = argparse.ArgumentParser()
  cccbuffer.Buffer._cmdline_args(parser)
  args = parser.parse_args([])
  return ('filename_or_obj', 'drop_variables', 'fused') + tuple(vars(args).keys())

# Register this as a backend for xarray, so CCC files can be directly
# opened with xarray.open_dataset
class CCCBackendEntrypoint(BackendEntrypoint):
  description = "Open CCC files in xarray."
  open_dataset_parameters = _get_open_dataset_parameters()
  def open_dataset (self, filename_or_obj, drop_variables=None, fused=False, **kwargs):
    if drop_variables is not None: kwargs['exclude'] = drop_variables
    return cccbuffer.Buffer(filename_or_obj, **kwargs).to_xarray(fused=fused)
