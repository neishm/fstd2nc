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

# Method for reading a block from a file.
def _read_block (filename, offset, length):
  import numpy as np
  # Scalar version first.
  if not hasattr(offset,'__len__'):
    # Skip addresses that are -1 (indicates no data available).
    # E.g. for masked data, if no corresponding mask available
    if offset < 0: return None
    with open(filename,'rb') as f:
      f.seek (offset,0)
      return np.fromfile(f,'B',length)
  # Vectorized version.
  out = []
  with open(filename,'rb') as f:
    for o, l in zip(offset, length):
      # Skip addresses that are -1 (indicates no data available).
      # E.g. for masked data, if no corresponding mask available
      if o < 0:
        out.append(None)
        continue
      f.seek (o,0)
      out.append(np.fromfile(f,'B',l))
  return out

# Helper method - write the given graph information into an index file.
def _write_graphs (f, ind, batch, graphs, concat_axis='time'):
  import numpy as np
  from pickle import dumps
  init = (ind == 0)

  # Global collection of all compound data (pickled).
  pickles = {}

  # Extract the filesnames needed for this batch.
  filename_len = len(f.dimensions['filename_len'])
  allfiles = f.variables['files'][ind*batch:(ind+1)*batch,:].view('|S%d'%filename_len).flatten()

  # Write dimensions / coordinates.
  # Non-concatenated coordinates only need to be written once.
  for var, args in graphs:
    # Skip dimensions (encoded as part of variable encoding below).
    if not hasattr(var,'axes'): continue
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
    # Evaluate any map objects.
    args = [list(a) if isinstance(a,map) else a for a in args]
    while i < len(args):
      assert isinstance(args[i][0],str)
      label = args[i][0]
      # Scalar argument
      if i == len(args)-2 or isinstance(args[i+2][0],str):
        scalar = args[i+1]
        # Convert boolean flags to bytes, for netcdf compatibility.
        if scalar[0] in (True,False):
          scalar = np.array(scalar,dtype='B')
        # Handle compound data (nested tuples).
        # NOTE: assuming the order of these tuples is always consistent across
        # every batch of files.
        if isinstance(scalar[0],tuple):
          # Get the unique data structures.
          unique = pickles.setdefault(label,{})
          for s in scalar:
            if id(s) not in unique:
              unique[id(s)] = s
          ids, structs = zip(*unique.items())
          if label+'_lookup' not in g.variables:
            g.createVariable(label+'_lookup','i4',dims,zlib=True,chunksizes=outer_shape)
          g.variables[label+'_lookup'][outer_sl] = np.array([ids.index(id(s)) for s in scalar],'int32').reshape(outer_shape)
          i += 2
          continue
          # End of compound data case
        # Initialize the scalar variable?
        if label not in g.variables:
          g.createVariable(label,type(scalar[0]),dims,zlib=True,chunksizes=outer_shape)
        # Fill in the scalar values.
        v = g.variables[label]
        v[outer_sl] = np.array(scalar).reshape(outer_shape)
        i += 2
        continue
      # Handle arguments from address / length pairs.
      # Skip for case where there are no valid records.
      if np.all(np.array(args[i+1])==-1):
        i += 3
        continue
      if label not in g.groups:
        g.createGroup(label)
      al = g.groups[label]
      if 'address' not in al.variables:
        al.createVariable('address','i8',dims,zlib=True,chunksizes=outer_shape)
      address = al.variables['address']
      address[outer_sl] = np.array(args[i+1],'int64').reshape(outer_shape)
      if 'length' not in al.variables:
        al.createVariable('length','i4',dims,zlib=True,chunksizes=outer_shape)
      length = al.variables['length']
      length[outer_sl] = np.array(args[i+2],'int32').reshape(outer_shape)
      i += 3

  if len(pickles) == 0: return
  # Final encoding of compound arguments.
  if 'pickle' not in f.vltypes:
    f.createVLType('B','pickle')
  for label, unique in pickles.items():
    ids, structs = zip(*unique.items())
    structs = [np.array([dumps(s)]).view('B') for s in structs]
    if label+'_index' not in f.dimensions:
      f.createDimension(label+'_index',len(structs))
    if label+'_pickle' not in f.variables:
      f.createVariable(label+'_pickle',f.vltypes['pickle'],(label+'_index',),zlib=True)
      for ind,s in enumerate(structs):
        f.variables[label+'_pickle'][ind] = s


class ExternOutput (BufferBase):

  # Helper method to read, decode, and post-process the data in one go.
  @classmethod
  def _dasked_read (cls, *args):
    import numpy as np
    filename = args[0]
    shape = args[-1]
    # Extract arguments, read and decode any data arguments.
    kwargs = {}
    key = 'data'
    i = 1
    while i < len(args)-1:
      if isinstance(args[i],str):
        key = args[i]
        i = i + 1
        continue
      # Single argument - scalar or external array?
      if isinstance(args[i+1],str) or i == len(args)-2:
        kwargs[key] = args[i]
        i = i + 1
        continue
      # Dual arguments (address + length)?
      else:
        offset = args[i]
        length = args[i+1]
        ###
        # Scalar version first.
        if not hasattr(offset,'__len__'):
          # Skip addresses that are -1 (indicates no data available).
          # E.g. for masked data, if no corresponding mask available
          if offset < 0:
            kwargs[key] = None
            i = i + 2
            continue
          try:
            with open(filename,'rb') as f:
              f.seek (offset,0)
              data = np.fromfile(f,'B',length)
              data = cls._decode(data)
              kwargs[key] = data
          except Exception:
            warn(_("Cannot read from %s - file may be missing or damaged.")%filename)
            kwargs[key] = None
          i = i + 2
          continue
        # Vectorized version.
        kwargs[key] = []
        try:
          with open(filename,'rb') as f:
            for o, l in zip(offset, length):
              # Skip addresses that are -1 (indicates no data available).
              # E.g. for masked data, if no corresponding mask available
              if o < 0:
                kwargs[key].append(None)
                continue
              f.seek (o,0)
              data = np.fromfile(f,'B',l)
              data = cls._decode(data)
              kwargs[key].append(data)
        except Exception:
          warn(_("Cannot read from %s - file may be missing or damaged.")%filename)
          kwargs[key] = [None]*len(offset)
        i = i + 2
        continue
    # Post-processing
    # Scalar case:
    if not any(isinstance(v,list) for v in kwargs.values()):
      # Skip post-processing for missing data.
      if kwargs.get('data',None) is None:
        return np.full(shape,float('nan'))
      data = cls._postproc(**kwargs)
      return data.reshape(shape)
    # Vector case
    # First, broadcast any scalar arguments into vector arguments.
    nrec = [len(v) for v in kwargs.values() if isinstance(v,list)][0]
    out = []
    for i in range(nrec):
      kw = dict((k,v[i] if isinstance(v,list) else v) for k,v in kwargs.items())
      # Skip post-processing for missing data.
      if kw.get('data',None) is None:
        out.append(None)
      out.append(cls._postproc(**kw))
    # Fill in any missing data.
    subshapes = [o.shape for o in out if o is not None]
    if len(subshapes) == 0: return np.full(shape, float('nan'))
    subshape = subshapes[0]
    out = [np.full(subshape,float('nan'),'float32') if o is None else o for o in out]
    return np.array(out).reshape(shape)

  def _iter_graph (self, include_coords=True, fused=False):
    """
    Iterate over all the variables, and generate graphs for them.
    """
    from fstd2nc.mixins import _iter_type, _chunk_type
    import numpy as np
    files = np.array(self._files, dtype=object)
    self._makevars()
    for var in self._iter_objects():
      if not include_coords:
        if var not in self._varlist:
          continue
      if not isinstance(var,(_iter_type,_chunk_type)):
        yield var, None
        continue
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
          dx = np.sum(sample==sample[0])
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
      args = [files[file_ids]]
      # Add data sources (primary data plus possible secondary like mask).
      for key, (addr_col, len_col, dname) in self._decoder_data:
        if addr_col not in self._headers: continue # Skip inactive arguments.
        # Source data is coming from dask / numpy?
        if dname in self._headers:
          d = self._headers[dname][record_id]
          args.extend([[key]*nchunks, d])
        # Source data is from disk?
        else:
          fname = files[file_ids]
          addr = self._headers[addr_col][record_id]
          addr[record_id<0] = -1
          if addr.ndim > 1: addr = map(np.array,addr)
          else: addr = map(int,addr)
          length = self._headers[len_col][record_id].astype('int32')
          if length.ndim > 1: length = map(np.array,length)
          else: length = map(int,length)
          args.extend([[key]*nchunks, addr, length])
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
      yield var, args

  # Helper method - given a list of files, produce the graphs.
  @classmethod
  def _graph_maker (cls, **kwargs):
    def get_graphs (infiles, first=[True]):
      import fstd2nc
      # Only print warning messages for first batch of files, assume the
      # warnings will be the same for the rest of the files as well.
      if first == [True]:
        first[0] = False
      else:
        fstd2nc.stdout.streams = ('error',)
      # Catch any problems with dates out of bounds.
      # pandas can't handle dates in 2263 or later.
      try:
        b = cls(infiles, **kwargs)
        graphs = list(b._iter_graph())
      except Exception as e:
        if 'Out of bounds' in str(e):
          warn (_("Dates out of bounds for pandas routines.  Using alternative (slower) routines."))
          fstd2nc.mixins._use_pandas = lambda: False
          _kwargs = dict(kwargs, progress=False)
          b = cls(infiles, **_kwargs)
          graphs = list(b._iter_graph())
        else:
          raise
      return graphs
    return get_graphs

  def _iter_dask (self, include_coords=True, fused=False, graph_iterator=None):
    """
    Iterate over all the variables, and convert to dask arrays.
    """
    from fstd2nc.mixins import _var_type
    from itertools import product
    import numpy as np
    from dask.base import tokenize
    from dask import array as da
    unique_token = tokenize(self._files,id(self))
    if graph_iterator is None:
      graph_iterator = self._iter_graph (include_coords=include_coords, fused=fused)
    for var, args in graph_iterator:
      if args is None:
        yield var
        continue
      # Don't need to label the primary data argument (called 'data').
      if args[1][0] == "data":
        args = args[0:1] + args[2:]
      name = var.name+"-"+unique_token
      nchunks = len(args[0])
      # Add shape as final argument.
      args.append(product(*var.chunks))
      graphs = zip([self._dasked_read]*nchunks, *args)
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
          dsk[key] = graph
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

  def to_xarray (self, fused=False):
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
        out[var.name].encoding['preferred_chunks'] = dict(zip(var.dims,chunk_shape))
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

  def to_xarray_list (self, fused=False):
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
        out[var.name].encoding['preferred_chunks'] = dict(zip(var.dims,chunk_shape))
        out[var.name].encoding['original_shape'] = out[var.name].shape
      # Construct the Dataset from the variable.
      out = xr.Dataset(out)
      # Decode CF metadata
      out = xr.conventions.decode_cf(out)
      # Make the time dimension unlimited when writing to netCDF.
      out.encoding['unlimited_dims'] = ('time',)
      out_list.append(out)

    return out_list

  # Generate an index file for fast loading of a dataset.
  @classmethod
  def make_index (cls, filenames, indexfile, batch=None, **kwargs):
    from fstd2nc.stdout import _, info, warn, error
    from fstd2nc.mixins import _expand_files, _FakeBar, _ProgressBar
    import fstd2nc
    import numpy as np
    import xarray as xr
    import netCDF4 as nc  # For on-disk storage of meta information.
    import os
    # Need a consistent reference date across all batches.
    if 'reference_date' not in kwargs and batch is not None:
      kwargs['reference_date'] = '1900-01-01'
    # First, get full list of files.
    if isinstance(filenames,list):
      allfiles = filenames
    else:
      allfiles = _expand_files(filenames)
      allfiles = list(zip(*allfiles))[1]
      allfiles = sorted(allfiles)
    if batch is None:
      batch = len(allfiles)
    if len(allfiles) % batch != 0:
      error(_("'batch' size (%d) does not divide evenly into number of files (%s)")%(batch,len(allfiles)))
    nbatch = len(allfiles) // batch
    # Construct a progress bar encompassing the whole file set.
    if kwargs.get('progress',False) is True:
      kwargs['progress'] = _ProgressBar(_("Indexing the files"), suffix='%(percent)d%% (%(index)d/%(max)d)', max=len(allfiles))

    # Define the index file.
    if os.path.exists(indexfile):
      os.remove(indexfile)
    f = nc.Dataset(indexfile,'w')
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
    # Get the funtion for generating the graphs.
    # The reason for currying the function was originally to support calling
    # the graph maker from threads (using just the file_batches as the input
    # argument).  However, in practice this didn't give much of a speed
    # improvement.
    graph_maker = cls._graph_maker(**kwargs)
    # Iterate through the graphs, write to the index file.
    for ind, file_batch in enumerate(file_batches):
      graphs = graph_maker (file_batch)
      _write_graphs(f, ind, batch, graphs)
      f.sync()
    # Restore original I/O streams
    fstd2nc.stdout.streams = orig_streams
    # Finalize the progress bar.
    if 'progress' in kwargs:
      if hasattr(kwargs['progress'],'finish'):
        kwargs['progress'].finish()
    f.close()

  # Create a dataset from the given index file.
  @classmethod
  def open_index (cls, indexfile):
    import xarray as xr
    import netCDF4 as nc
    from fstd2nc._xarray_hook import IndexFile, FSTDBackendArray
    # Open with xarray to get axes / coordinates.
    ds = xr.open_dataset(indexfile)
    # Open with netCDF4 get get access to the group structures containing the
    # metadata.
    index = IndexFile(indexfile)
    # Generate the variables.
    vardict = {}
    pickles = {}
    for varname, group in index.root.groups.items():
      # Get a lazy accessor for the data.
      array = FSTDBackendArray (cls, varname, index, pickles)
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
          preferred_chunks[dim] = len(index.root.dimensions[dim])
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
    # Decode CF metadata
    ds = xr.conventions.decode_cf(ds)
    return (ds)

  @classmethod
  def from_xarray (cls, ds, **params):
    """
    Create a Buffer object from the given xarray object.
    """
    from fstd2nc.mixins import _var_type, _iter_type, _dim_type, _axis_type
    import numpy as np
    varlist = []
    # Handle FSTD parameters passed in the method call.
    ds = ds.copy()
    ds.attrs.update(params)
    # Handle FSTD parameters passed in as global attributes of the Dataset.
    for varname, var in ds.variables.items():
      var.attrs.update(ds.attrs)
    # Handle FSTD parameters set as variable attributes.
    for varname, var in ds.variables.items():
      var.encoding.update(fstd_attrs=var.attrs)
    # Collect dimensions and coordinates into a separate structure.
    dims = {}
    coords = {}
    for dimname,dimsize in ds.dims.items():
      if dimname in ds.variables:
        dims[dimname] = _axis_type(dimname,dict(ds[dimname].attrs),np.array(ds[dimname]))
      else:
        dims[dimname] = _dim_type(dimname,dimsize)
    for coordname, coord in ds.coords.items():
      if coordname in dims: continue  # Already counted as a dimension.
      coords[coordname] = _var_type(coordname, dict(coord.attrs), [dims[d] for d in coord.dims], np.array(coord))
    # Construct the varlist.
    for varname, var in ds.variables.items():
      if varname in dims: continue
      if varname in coords: continue
      encoded = _var_type(varname, dict(var.attrs), [dims[d] for d in var.dims], ds[varname].data)
      # Add coordinates.
      encoded.atts.setdefault('coordinates',[coords[coordname] for coordname in ds[varname].coords if coordname in coords])
      varlist.append(encoded)
    # Decode the varlist into a table.
    return cls._deserialize(varlist)

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
    fields = ['nomvar', 'typvar', 'etiket', 'ni', 'nj', 'nk', 'ip1', 'ip2', 'ip3', 'deet', 'npas', 'datyp', 'nbits', 'grtyp', 'ig1', 'ig2', 'ig3', 'ig4', 'datev']
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
          d[i] = self._headers['d'][i].T
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
    """
    Create a Buffer object from an fstpy table of records.

    Parameters
    ----------
    table : pandas DataFrame
        The table of records from fstpy.
    """
    import numpy as np
    if hasattr(table,'to_pandas'):
      table = table.to_pandas()
    # Construct the record header info from the table.
    fields = ['nomvar', 'typvar', 'etiket', 'ni', 'nj', 'nk', 'ip1', 'ip2', 'ip3', 'deet', 'npas', 'datyp', 'nbits', 'grtyp', 'ig1', 'ig2', 'ig3', 'ig4', 'datev']
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
    # Fake file id (just so _read_record function doesn't crash).
    headers['file_id'] = np.empty(len(headers['nomvar']), dtype='int32')
    headers['file_id'][:] = -1
    # Add in data.
    headers['d'] = np.empty(len(headers['nomvar']), dtype=object)
    headers['d'][:] = [d.T for d in table['d']]

    # Encapsulate this info in a structure.
    fake_buffer = cls.__new__(cls)
    fake_buffer._files = [None]
    fake_buffer._headers = headers

    # Initialize a Buffer object with this info.
    b = cls(fake_buffer, **kwargs)

    return b

