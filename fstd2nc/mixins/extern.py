###############################################################################
# Copyright 2017-2021 - Climate Research Division
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
import collections

#################################################
# Provide various external array interfaces for the FSTD data.

# Method for reading a block from a file.
def _read_block (filename, offset, length):
  import numpy as np
  print ("READ OFFSET %X"%offset)
  with open(filename,'rb') as f:
    f.seek (offset,0)
    return np.fromfile(f,'B',length)

# Open a file for raw binary access using dask.
def _open_binary (filename):
  from fstd2nc.extra import blocksize
  from os.path import getsize
  from dask import array as da
  # Get a good chunk size for splitting the binary file into dask chunks.
  # Use the filesystem block size unless it's too small.
  chunksize = max(blocksize(filename), 2**20)
  filesize = getsize(filename)
  dsk = dict()
  chunks = []
  for n, offset in enumerate(range(0,filesize,chunksize)):
    end = min(offset+chunksize,filesize)
    length = end - offset
    chunks.append(length)
    key = (filename,n)
    dsk[key] = (_read_block, filename, offset, length)
  array = da.Array(dsk, filename, [chunks], 'B')
  return array

# Open an FSTD file, return a table with delayed access to the data.
def _open_fstd (filename, table=None):
  from fstd2nc.extra import all_params, decode
  from fstd2nc.mixins import dtype_fst2numpy
  import pandas as pd
  from dask import delayed
  from dask import array as da
  import numpy as np
  decode = delayed(decode)
  if table is None:
    with open(filename) as f:
      table = all_params(f)
  table = pd.DataFrame(table)
  b = _open_binary(filename)
  dlist = np.zeros(table.shape[0], dtype='object')
  for i, swa, lng, ni, nj, datyp, nbits in table[['swa','lng','ni','nj','datyp','nbits']].itertuples():
    offset = swa*8-8
    length = lng*8
    raw = b[offset:offset+length]
    data = decode(raw)
    #TODO: speedup
    dtype = dtype_fst2numpy (datyp, nbits)
    shape = (ni, nj)
    d = da.from_delayed (data, shape, dtype)
    dlist[i] = d
  table['d'] = dlist
  return table

class ExternOutput (BufferBase):

  _read_chunk_cache = None

  # The interface for getting chunks into dask.
  def _read_chunk (self, rec_id, shape, dtype):
    with self._lock:
      # Cache the last record read from this interface, in case it is
      # immediately requested again.
      # Can happen, for instance, if the same data is sliced in multiple
      # different ways.
      if self._read_chunk_cache is None:
        self._read_chunk_cache = {}
      if rec_id not in self._read_chunk_cache:
        self._read_chunk_cache.clear()
        self._read_chunk_cache[rec_id] = self._fstluk(rec_id,dtype=dtype)['d'].transpose().reshape(shape)
      return self._read_chunk_cache[rec_id]

  def _iter_dask (self):
    """
    Iterate over all the variables, and convert to dask arrays.
    """
    from fstd2nc.mixins import _iter_type, _chunk_type, _var_type
    from dask import array as da
    from dask.base import tokenize
    import numpy as np
    from itertools import product
    _add_callback()
    unique_token = tokenize(self._files,id(self))
    # Make a local copy of two columns from the header table.
    # This speeds up access to their elements in the inner loop.
    all_file_ids = np.array(self._headers['file_id'],copy=True)
    all_keys = np.array(self._headers['key'],copy=True)
    ###
    # For science network installation, make sure the stack size of child
    # threads is large enough to accomodate workspace for librmn.
    from os.path import basename, splitext
    from rpnpy.librmn import librmn
    libname = basename(getattr(librmn,'_name',''))
    # Skip this step for patched versions of librmn (ending in -rpnpy) that
    # use the heap instead of the stack for workspace.
    if not splitext(libname)[0].endswith('-rpnpy'):
      import threading
      # Use the size of the largest record, plus 1MB just in case.
      # Use multiples of 4K for the size, as suggested in Python theading
      # documentation.
      # https://docs.python.org/3/library/threading.html#threading.stack_size
      # Allow for possibility of double precision values (64-bit)
      arraysize = np.max(8 * self._headers['ni'] * self._headers['nj'] * self._headers['nk'])
      stacksize = arraysize//4096*4096 + (1<<20)
      if stacksize > threading.stack_size():
        threading.stack_size(stacksize)
    ###
    self._makevars()
    for var in self._iter_objects():
      if not isinstance(var,(_iter_type,_chunk_type)):
        yield var
        continue
      name = var.name+"-"+unique_token
      ndim = len(var.axes)
      shape = var.shape
      # Convert _iter_type to more generic _chunk_type.
      if isinstance(var,_iter_type):
        chunks = {}
        ndim_outer = var.record_id.ndim
        ndim_inner = ndim - ndim_outer
        chunk_shape = shape[ndim_outer:]
        for ind in product(*map(range,var.record_id.shape)):
          rec_id = var.record_id[ind]
          ind = ind + tuple((0,dx) for dx in shape[ndim_outer:])
          chunks[ind] = rec_id
        var = _chunk_type (var.name, var.atts, var.axes, var.dtype, chunks, chunk_shape)
      # Convert _chunk_type to dask Array objects.
      if isinstance(var,_chunk_type):
        ndim_inner = len(var.chunksize)
        ndim_outer = ndim - ndim_inner
        # Get chunk dimensions.
        # First, size of single (untruncated) chunk, full indices.
        untruncated_chunksize = (1,)*(ndim-len(var.chunksize)) + var.chunksize
        # Next, breakdown of chunks along all variable dimensions.
        chunks = []
        chunk_indices = []
        for i in range(ndim):
          dx = untruncated_chunksize[i]
          ch = tuple(dx for j in range(dx,shape[i]+1,dx))
          if shape[i] % dx > 0:
            ch = ch + (shape[i] % dx, )
          chunks.append(ch)
          chunk_indices.append(range(len(ch)))
        # Loop over all indices, generate dask graph.
        dsk = dict()
        for ind, chunk_shape in zip(product(*chunk_indices), product(*chunks)):
          # Unique key for this graph member.
          key = (name,) + ind
          # Get record id.
          slices = [(i*dx,i*dx+res) for i,dx,res in zip(ind,untruncated_chunksize,chunk_shape)]
          slices[:ndim_outer] = ind[:ndim_outer]
          rec_id = var.chunks.get(tuple(slices),-1)
          # Add this record as a chunk in the dask Array.
          # Also, specify the preferred order of reading the chunks within the
          # file.
          if rec_id >= 0:
            file_id = all_file_ids[rec_id]
            filename = self._files[file_id]
            rec_key = all_keys[rec_id]
            # Special case: already have a dask object from external source.
            # (E.g., from fstpy)
            if hasattr(self, '_extern_table'):
              d = self._extern_table['d'].iloc[rec_id]
              if hasattr(d,'compute'):
                # Start a dask graph using the external dask array as the source.
                dsk[key] = (d.compute,)
                # Ensure the dask array includes the degenerate outer dimensions.
                # Otherwise get a runtime error if slicing is done on this.
                dsk[key] = (np.ravel, dsk[key], 'K')
                dsk[key] = (np.reshape, dsk[key], chunk_shape)
              else:  # Special case: have a numpy array in memory.
                dsk[key] = d
                dsk[key] = np.ravel(dsk[key], 'K')
                dsk[key] = np.reshape(dsk[key], chunk_shape)
            # Otherwise, construct one with our own dask wrapper.
            else:
              dsk[key] = (_preferred_chunk_order,filename,rec_key,(self._read_chunk, rec_id, chunk_shape, var.dtype))
          else:
            # Fill missing chunks with fill value or NaN.
            if hasattr(self,'_fill_value'):
              var.atts['_FillValue'] = self._fill_value
              dsk[key] = (np.full, chunk_shape, self._fill_value, var.dtype)
            else:
              dsk[key] = (np.full, chunk_shape, float('nan'), var.dtype)
        array = da.Array(dsk, name, chunks, var.dtype)
        var = _var_type(var.name,var.atts,var.axes,array)
      yield var

  def to_xarray (self):
    """
    Create an xarray interface for the RPN data.
    Requires the xarray and dask packages.
    """
    from collections import OrderedDict
    import xarray as xr
    out = OrderedDict()
    for var in self._iter_dask():
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

    return out

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
    from fstpy.dataframe import add_grid_column
    from fstpy.std_io import add_dask_column
    # Special case: our data is already from an fstpy table, not from an FSTD
    # file in our control.
    # E.g., if some smartass does Buffer.from_fstpy(df).to_fstpy()
    if hasattr(self, '_extern_table'):
      return self._extern_table
    # Put all the header info into a dictionary.
    fields = ['nomvar', 'typvar', 'etiket', 'ni', 'nj', 'nk', 'dateo', 'ip1', 'ip2', 'ip3', 'deet', 'npas', 'datyp', 'nbits', 'grtyp', 'ig1', 'ig2', 'ig3', 'ig4', 'datev']
    table = dict()
    # Create a mask to exclude deleted / overwritten records.
    mask = self._headers['dltf'] == 0
    for field in fields:
      col = self._headers[field][mask]
      # Convert byte arrays to strings, which is what fstpy expects.
      if col.dtype.str.startswith('|S'):
        col = np.asarray(col,dtype=col.dtype.str.replace('|S','<U'))
      table[field] = col
    # Convert to pandas table.
    table = pd.DataFrame.from_dict(table)
    # Add grid info.
    add_grid_column (table)
    # Temporarily insert some extra columns needed for the data.
    table['shape'] = list(zip(table['ni'],table['nj']))
    filenames = dict((i,f) for i,f in enumerate(self._files))
    table['path'] = pd.Series(self._headers['file_id'][mask]).map(filenames)
    table['key'] = (self._headers['key'][mask] << 10)
    # Generate dask objects
    #TODO: use our own, in case we modified the data?
    # (doesn't normally happen, but you never know...)
    # For instance could happen if interp is used.
    add_dask_column(table)
    # Clean up temporary columns and return.
    table.drop(columns=['shape','path','key'], inplace=True)
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
    import tempfile
    from os import path
    import rpnpy.librmn.all as rmn
    import numpy as np
    if hasattr(table,'to_pandas'):
      table = table.to_pandas()
    # Construct the record header info from the table.
    fields = ['nomvar', 'typvar', 'etiket', 'ni', 'nj', 'nk', 'dateo', 'ip1', 'ip2', 'ip3', 'deet', 'npas', 'datyp', 'nbits', 'grtyp', 'ig1', 'ig2', 'ig3', 'ig4', 'datev']
    headers = np.zeros(len(table), dtype=cls.__new__(cls,**kwargs)._headers_dtype)
    for col in fields:
      headers[col][:] = table[col]
    # Pad out string variables with spaces.
    headers['nomvar'] = np.char.ljust(headers['nomvar'], 4, ' ')
    headers['typvar'] = np.char.ljust(headers['typvar'], 2, ' ')
    headers['etiket'] = np.char.ljust(headers['etiket'], 12, ' ')
    # Generate temporary file with target grid info.
    try: # Python 3
      grid_tmpdir = tempfile.TemporaryDirectory()
      gridfile = path.join(grid_tmpdir.name,"grid.fst")
    except AttributeError: # Python 2 (no auto cleanup)
      grid_tmpdir = tempfile.mkdtemp()
      gridfile = path.join(grid_tmpdir,"grid.fst")
    # Write all grid records to a temporary file, so they are accessible to
    # the vgrid / librmn helper functions.
    iun = rmn.fstopenall(gridfile, rmn.FST_RW)
    for nomvar in ('!!','>>','^^','^>','!!5F'):
      for ind in np.where(table['nomvar'] == nomvar)[0]:
        rec = table.iloc[ind].to_dict()
        rec['d'] = np.asarray(rec['d'])
        rmn.fstecr(iun, rec)
    rmn.fstcloseall(iun)

    # Initialize the Buffer object with this info.
    b = cls(gridfile, header_cache={'__ROOT__'+gridfile:headers}, **kwargs)
    b._grid_tmpdir = grid_tmpdir  # Save tmpdir until cleanup.

    # Save the dataframe for reference.
    # Will need the dask objects for getting the data.
    b._extern_table = table

    return b

  # Handle external data sources.
  # Overrides the usual reading of data from a file.
  def _fstluk (self, rec_id, dtype=None, rank=None, dataArray=None):
    import numpy as np
    # Check if there is custom data enabled for this Buffer.
    if hasattr(self, '_extern_table'):
      # Make sure we are looking for something in our list of records.
      # (could be a key pointing into something else?)
      if not isinstance(rec_id,dict):
        # Extract the record info from the table.
        rec = self._extern_table.iloc[rec_id].to_dict()
        # Load the data (if delayed).
        rec['d'] = np.asarray(rec['d'])
        return rec
    # Otherwise, continue as usual.
    return super(ExternInput,self)._fstluk (rec_id, dtype, rank, dataArray)

