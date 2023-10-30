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
  return ('filename_or_obj', 'drop_variables', 'fused', 'batch', 'cachefile') + tuple(vars(args).keys())

# Helper method - given a list of files, produce the graphs.
def graph_maker (**kwargs):
  from threading import Lock
  _lock = Lock()
  def get_graphs (infiles, first=[True]):
    import fstd2nc
    # Only print warning messages for first batch of files, assume the
    # warnings will be the same for the rest of the files as well.
    with _lock:
      if first == [True]:
        first[0] = False
        return list(fstd2nc.Buffer(infiles, **kwargs)._iter_graph())
    fstd2nc.stdout.streams = ('error',)
    b = fstd2nc.Buffer(infiles, **kwargs)
    return list(b._iter_graph())
  return get_graphs

# Helper method - write the given graph information into a cache file.
def _write_graphs (f, ind, batch, graphs, concat_axis='time'):
  import numpy as np
  init = (ind == 0)

  # Extract the filesnames needed for this batch.
  filename_len = len(f.dimensions['filename_len'])
  allfiles = f.variables['files'][ind*batch:(ind+1)*batch,:].view('|S%d'%filename_len).flatten()

  # Write dimensions / coordinates.
  # Non-concatenated coordinates only need to be written once.
  for var, args in graphs:
    sl = [slice(None)]*len(var.axes)
    dtype = var.dtype if hasattr(var,'dtype') else var.array.dtype
    concatenate = False
    ###
    # Define the coordinates (if not already defined).
    for iaxis,axis in enumerate(var.axes):
      # Check for concatenation axis.
      if axis.name.startswith(concat_axis):
        concatenate = True
        dt = len(axis)
        sl[iaxis] = slice(ind*dt,ind*dt+dt)
        full_length = None
      else:
        full_length = len(axis)
      if axis.name not in f.dimensions:
        f.createDimension(axis.name, full_length)
    ###
    sl = tuple(sl)

    # Write the variable (direct case)
    if args is None:
      if var.name not in f.variables:
        v = f.createVariable(var.name,dtype,var.dims,zlib=True)
        v.setncatts(var.atts)
      if init or concatenate:
        f.variables[var.name][sl] = var.array
      continue

    # Write the address / length info for indirect case.
    if var.name not in f.groups:
      f.createGroup(var.name)

    g = f.groups[var.name]

    if 'template' not in g.variables:
      v = g.createVariable('template',dtype,var.dims)
      v.setncatts(var.atts)

    # Define outer axes (for storing address / length values).
    ndim_outer = var.record_id.ndim
    while var.record_id.shape[ndim_outer-1] == 1 and var.shape[ndim_outer-1] != 1:
      ndim_outer -= 1

    dims = var.dims[:ndim_outer]
    outer_shape = var.shape[:ndim_outer]
    outer_sl = sl[:ndim_outer]
    # Set file_id
    files = np.array(args[0],dtype='|S%d'%filename_len)
    if 'file_id' not in g.variables:
      g.createVariable('file_id','i4',dims,zlib=True,chunksizes=outer_shape)
    file_id = g.variables['file_id']
    file_id[outer_sl] = np.searchsorted(allfiles,files).reshape(outer_shape) + ind*batch
    # Set address / length arguments.
    i = 1
    while i < len(args):
      assert isinstance(args[i][0],str)
      label = args[i][0]
      # Scalar argument
      if i == len(args)-2 or isinstance(args[i+2],str):
        if label not in g.variables:
          g.createVariable(label,type(args[i+1][0]),dims,zlib=True,chunksizes=outer_shape)
        v = g.variables[label]
        v[outer_sl] = np.array(args[i+1]).reshape(outer_shape)
        i += 2
        continue
      # Address / length arguments
      if label not in g.groups:
        g.createGroup(label)
      al = g.groups[label]
      if 'address' not in al.variables:
        al.createVariable('address','i8',dims,zlib=True,chunksizes=outer_shape)
      address = al.variables['address']
      address[outer_sl] = np.array(list(args[i+1]),'int64').reshape(outer_shape)
      if 'length' not in al.variables:
        al.createVariable('length','i4',dims,zlib=True,chunksizes=outer_shape)
      length = al.variables['length']
      length[outer_sl] = np.array(list(args[i+2]),'int32').reshape(outer_shape)
      i += 3


# Register this as a backend for xarray, so FST files can be directly
# opened with xarray.open_dataset
class FSTDBackendEntrypoint(BackendEntrypoint):
  description = "Open FST files in xarray."
  open_dataset_parameters = _get_open_dataset_parameters()
  def open_dataset (self, filename_or_obj, drop_variables=None, fused=False, batch=None, cachefile=None, **kwargs):
    from fstd2nc.stdout import _, info, warn, error
    from fstd2nc.mixins import _expand_files, _FakeBar, _ProgressBar
    import fstd2nc
    import numpy as np
    import xarray as xr
    import netCDF4 as nc  # For on-disk storage of meta information.
    from multiprocessing.pool import ThreadPool
    if drop_variables is not None: kwargs['exclude'] = drop_variables
    # Simple case: no batching done.
    if batch is None:
      return fstd2nc.Buffer(filename_or_obj, **kwargs).to_xarray(fused=fused)
    # Complicated case: processing the files in batches.
    if fused is True:
      warn(_("'fused' option ignored for batch processing."))
      fused = False
    if cachefile is None:
      error(_("Need to specify 'cachefile' when using 'batch' argument."))
    # First, get full list of files.
    allfiles = _expand_files(filename_or_obj)
    allfiles = list(zip(*allfiles))[1]
    allfiles = sorted(allfiles)
    if len(allfiles) % batch != 0:
      error(_("'batch' size (%d) does not divide evenly into number of files (%s)")%(batch,len(allfiles)))
    nbatch = len(allfiles) // batch
    # Construct a progress bar encompassing the whole file set.
    if kwargs.get('progress',False) is True:
      kwargs['progress'] = _ProgressBar(_("Checking input files"), suffix='%(percent)d%% (%(index)d/%(max)d)', max=len(allfiles))

    # Define the cache file.
    f = nc.Dataset(cachefile,'w')
    f.setncattr('version',1)
    # Write file information.
    filename_len = max(map(len,allfiles))
    f.createDimension('filename_len',filename_len)
    nfiles_dim = f.createDimension('nfiles',len(allfiles))
    file_var = f.createVariable('files','|S1',('nfiles','filename_len'),zlib=True,chunksizes=(batch,filename_len))
    allfiles = np.array(allfiles,dtype='|S%d'%filename_len)
    file_var[:,:] = allfiles.reshape(-1,1).view('|S1')

    # Put the files into batches for processing.
    file_batches = allfiles.reshape(nbatch,batch)
    # Remember original I/O streams (will get mangled by write_graphs).
    orig_streams = fstd2nc.stdout.streams
    # Start a thread pool for processing the graphs in parallel.
    #NOTE: disabled - not much speed improvement with this approach.
    # To re-enable, uncomment the following line, and then change 'map'
    # to 'pool.imap'.
    #pool = ThreadPool(2)
    all_graphs = map(graph_maker(**kwargs), file_batches)
    # Iterate through the graphs from this pool, write to the cache file.
    for ind, graphs in enumerate(all_graphs):
      _write_graphs(f, ind, batch, graphs)
    # Restore original I/O streams
    fstd2nc.stdout.streams = orig_streams
    f.close()
    return xr.Dataset({})

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

