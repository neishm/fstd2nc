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

# Helper method - write the given graph information into a cache file.
def _write_graphs (f, ind, nbatch, graphs, concat_axis='time'):
  import numpy as np
  init = (ind == 0)
  # Define batch / template groups.
  if 'batch' not in f.groups:
    f.createGroup('batch').createDimension('batch',nbatch)
  if 'templates' not in f.groups:
    f.createGroup('templates')
  b = f.groups['batch']
  t = f.groups['templates']

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
        full_length = len(axis)*nbatch
      else:
        full_length = len(axis)
      if axis.name not in b.dimensions:
        b.createDimension(axis.name, len(axis))
      if axis.name not in t.dimensions:
        t.createDimension(axis.name, full_length)
    ###
    sl = tuple(sl)

    if var.name not in t.variables:
      v = t.createVariable(var.name,dtype,var.dims)
      v.setncatts(var.atts)

    # Write the variable (direct case)
    if args is None:
      if init or concatenate:
        t.variables[var.name][sl] = var.array
      continue

    # Write the address / length info for indirect case.
    if var.name not in b.groups:
      b.createGroup(var.name)

    g = b.groups[var.name]

    # Define outer axes (for storing address / length values).
    ndim_outer = var.record_id.ndim
    while var.record_id.shape[ndim_outer-1] == 1 and var.shape[ndim_outer-1] != 1:
      ndim_outer -= 1

    dims = ('batch',) + var.dims[:ndim_outer]
    shape = var.shape[:ndim_outer]
    # Set file_id
    filename_len = len(f.dimensions['filename_len'])
    allfiles = f.variables['files'][:,:].view('|S%d'%filename_len)[0]
    files = np.array(args[0],dtype='|S%d'%filename_len)
    if 'file_id' not in g.variables:
      g.createVariable('file_id','i4',dims,zlib=True,chunksizes=(1,)+shape)
    file_id = g.variables['file_id']
    s = np.searchsorted(allfiles,files)
    file_id[0,...] = np.searchsorted(allfiles,files).reshape(shape)
    #TODO
    continue

    # Set address/length info
    for i in range(1,len(args)-2,3):
      label = args[i][0]
      address = vargroup.createVariable(label+'.address','i8',dims,zlib=True)
      address[0,...] = np.array(list(args[i+1]),'int64').reshape(shape[1:])
      length = vargroup.createVariable(label+'.length','i4',dims,zlib=True)
      length[0,...] = np.array(list(args[i+2]),'int32').reshape(shape[1:])
    # Define a blank version of the final variable with the final structure.
    template = vargroup.createVariable('template',var.dtype,var.dims)


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
    file_var = f.createVariable('files','|S1',('nfiles','filename_len'),zlib=True)
    allfiles = np.array(allfiles,dtype='|S%d'%filename_len)
    file_var[:,:] = allfiles.reshape(-1,1).view('|S1')

    # Process first batch, get basic structure of the dataset.
    b = fstd2nc.Buffer(allfiles[0:batch], **kwargs)
    graphs = list(b._iter_graph())
    #TODO
    _write_graphs (f, 0, nbatch, graphs)
    f.close()
    return xr.Dataset({})
    # Only print warning messages for first batch of files, assume the
    # warnings will be the same for the rest of the files as well.
    # Continue processing the rest of the batches, spawning threads for getting
    # the variable structures.
    # Turn warnings messages back on.
    #TODO
    ind = 0
    pieces = []
    original_streams = fstd2nc.stdout.streams
    while ind < len(allfiles):
      b = fstd2nc.Buffer(allfiles[ind:ind+batch], **kwargs)
      #TODO: spawn thread here.
      graphs = list(b._iter_graph())
      pieces.append(graphs)
      if ind == 0: fstd2nc.stdout.streams = ('error',)
      ind += batch
    fstd2nc.stdout.streams = original_streams
    # Finalize the progress bar.
    if 'progress' in kwargs and hasattr(kwargs['progress'],'finish'):
      kwargs['progress'].finish()

    #TODO
    import xarray as xr
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

