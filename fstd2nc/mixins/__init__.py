###############################################################################
# Copyright 2017 - Climate Research Division
#                  Environment and Climate Change Canada
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


#####################################################################
# Mixins for different features / behaviour for the conversions.

# This module contains the base class for the mixins, from which all others
# are derived.


from fstd2nc.stdout import _, info, warn, error

# Override default dtype for "binary" (datyp 0) data.
# This is probably float32 data, not integer.
# See (internal) mailing list discussion at http://internallists.cmc.ec.gc.ca/pipermail/python-rpn/2015-February/000010.html where this problem is discussed.
# Also see examples on (internal) wiki dealing with this type of data:
# https://wiki.cmc.ec.gc.ca/wiki/Python-RPN/2.0/examples#Example_4:_Plot_RPN_STD_field_on_an_X-grid
# https://wiki.cmc.ec.gc.ca/wiki/Talk:Python-RPN/2.0/examples#Plot_GIOPS_Forecast_Data_with_Basemap
def dtype_fst2numpy (datyp, nbits=None):
  from rpnpy.librmn.fstd98 import dtype_fst2numpy
  if datyp == 0:
    warn (_("Raw binary records detected.  The values may not be properly decoded if you're opening on a different platform."))
    datyp = 5
  return dtype_fst2numpy(datyp,nbits)

# Decorator for efficiently converting a scalar function to a vectorized
# function.
def vectorize (f):
  from functools import wraps
  try:
    from pandas import Series, unique
    @wraps(f)
    def vectorized_f (x):
      # If we're given a scalar value, then simply return it.
      if not hasattr(x,'__len__'):
        return f(x)
      # Get unique values
      inputs = unique(x)
      outputs = map(f,inputs)
      table = dict(zip(inputs,outputs))
      result = Series(x).map(table)
      return result.values
  except ImportError:
    def cached_f(x, cache={}):
      if x not in cache:
        cache[x] = f(x)
      return cache[x]
    @wraps(f)
    def vectorized_f (x):
      # If we're given a scalar value, then simply return it.
      if not hasattr(x,'__len__'):
        return cached_f(x)
      return map(cached_f,x)
  return vectorized_f


# The type of data returned by the Buffer iterator.
class _var_type (object):
  __slots__ = ('name','atts','axes','array')
  def __init__ (self, name, atts, axes, array):
    self.name = name
    self.atts = atts
    self.axes = axes
    self.array = array
  def __iter__ (self):
    return (getattr(self,n) for n in self.__slots__)

# An internal type used in _iter methods.
class _iter_type (object):
  __slots__ = ('name','atts','axes','dtype','record_id')
  def __init__ (self, name, atts, axes, dtype, record_id):
    self.name = name
    self.atts = atts
    self.axes = axes
    self.dtype = dtype
    self.record_id = record_id
  def __iter__ (self):
    return (getattr(self,n) for n in self.__slots__)

# Helper method - modify the axes of an array.
def _modify_axes (axes, **kwargs):
  from collections import OrderedDict
  new_axes = OrderedDict()
  for name,values in axes.items():
    # Check if modifications requested.
    if name in kwargs:
      # Can either modify just the axis name, or the name and values too.
      if isinstance(kwargs[name],str):
        name = kwargs[name]
      else:
        name,values = kwargs[name]
    new_axes[name] = values
  return new_axes

# Fake progress bar - does nothing.
class _FakeBar (object):
  def __init__ (self, *args, **kwargs): pass
  def iter(self, it):
    for i in it: yield i
  def next(self): pass
  def finish(self): pass

# Try importing progress module.
try:
  from progress.bar import IncrementalBar
  class _ProgressBar(IncrementalBar):
    # Define a custom ETA, which prints the elapsed time at the end.
    @property
    def myeta(self):
      from datetime import timedelta
      # Print elapsed time at the end.
      if self.index == self.max:
        return self.elapsed_td
      # Don't update the eta too quickly, or the estimate will jump around.
      if not hasattr(self,'_last_eta_update'):
        self._last_eta_update = self._ts
      if self._ts - self._last_eta_update < 1.0 and hasattr(self,'_last_eta'):
        return self._last_eta
      dt = self._ts - self.start_ts
      if dt == 0: return '??:??:??'
      speed = self.index / dt
      time_remaining = (self.max-self.index) / speed
      self._last_eta = timedelta(seconds=int(time_remaining))
      self._last_eta_update = self._ts
      return self._last_eta
    # Rate limit the message update so it only happens once a second.
    def update(self):
      if not hasattr(self,'_last_update'):
        self._last_update = self._ts
      if self.index < self.max:
        if self._ts - self._last_update < 0.1:
          return
      super(_ProgressBar,self).update()
      self._last_update = self._ts
  del IncrementalBar
except ImportError:
  _ProgressBar = _FakeBar


# Define a class for encoding / decoding FSTD data.
# Each step is placed in its own "mixin" class, to make it easier to patch in 
# new behaviour if more exotic FSTD files are encountered in the future.
class BufferBase (object):

  # Keep a reference to fstd98 so it's available during cleanup.
  try:
    from rpnpy.librmn import fstd98 as _fstd98
  except ImportError: pass

  # Names of records that should be kept separate (never grouped into
  # multidimensional arrays).
  _meta_records = ()

  # Attributes which could potentially be used as axes.
  _outer_axes = ()

  # Attributes which uniquely identify a variable.
  # Note: nomvar should always be the first attribute
  _var_id = ('nomvar','ni','nj','nk')

  # Similar to above, but a human-readable version of the id.
  # Could be used for appending suffixes to variable names to make them unique.
  # Uses string formatting operations on the variable metadata.
  _human_var_id = ('%(nomvar)s', '%(ni)sx%(nj)s', '%(nk)sL')

  # Field names and types for storing record headers in a structured array.
  _headers_dtype = [
     ('file_id','int32'),('key','int32'),
     ('dateo','int32'),('datev','int32'),('deet','int32'),('npas','int32'),
     ('ni','int32'),('nj','int32'),('nk','int32'),('nbits','int32'),
     ('datyp','int32'),('ip1','int32'),('ip2','int32'),('ip3','int32'),
     ('typvar','|S2'),('nomvar','|S4'),('etiket','|S12'),('grtyp','|S1'),
     ('ig1','int32'),('ig2','int32'),('ig3','int32'),('ig4','int32'),
     ('swa','int32'),('lng','int32'),('dltf','int32'),('ubc','int32'),
     ('xtra1','int32'),('xtra2','int32'),('xtra3','int32'),
  ]

  # Record parameters which should not be used as nc variable attributes.
  # (They're either internal to the file, or part of the data, not metadata).
  _ignore_atts = ('file_id','swa','lng','dltf','ubc','xtra1','xtra2','xtra3','key','shape','d')

  # Define any command-line arguments for reading FSTD files.
  @classmethod
  def _cmdline_args (cls, parser):
    from fstd2nc import __version__
    from argparse import SUPPRESS
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('--progress', action='store_true', help=_('Display a progress bar during the conversion.  Requires the "progress" module.'))
    parser.add_argument('--minimal-metadata', action='store_true', help=_("Don't include RPN record attributes and other internal information in the output metadata."))
    parser.add_argument('--ignore-typvar', action='store_true', help=_('Tells the converter to ignore the typvar when deciding if two records are part of the same field.  Default is to split the variable on different typvars.'))
    parser.add_argument('--ignore-etiket', action='store_true', help=_('Tells the converter to ignore the etiket when deciding if two records are part of the same field.  Default is to split the variable on different etikets.'))
    parser.add_argument('--no-quick-scan', action='store_true', help=SUPPRESS)
    #help=_('Don't read record headers from the raw librmn structures, instead call fstprm.  This can be much slower, but safer.')

  # Do some checks on the command-line arguments after parsing them.
  @classmethod
  def _check_args (cls, parser, args):
    return  # Nothing to check here.

  # Clean up a buffer (close any attached files, etc.)
  def __del__ (self):
    try:
      self._close()
      from rpnpy.librmn.fstd98 import fstcloseall
      fstcloseall(self._meta_funit)
    except (ImportError, AttributeError):
      pass  # May fail if Python is doing a final cleanup of everything.
            # Or, if buffer wasn't fully initialized yet.

  # Extract metadata from a particular header.
  def _get_header_atts (self, header):
    for n,v in header.items():
      if n in self._ignore_atts: continue
      if isinstance(v,str):
        v = v.strip()
      yield (n,v)

  # Open the specified file (by index)
  def _open (self, file_id):
    from rpnpy.librmn.base import fnom
    from rpnpy.librmn.fstd98 import fstouv
    from rpnpy.librmn.const import FST_RO
    from fstd2nc.extra import librmn
    opened_file_id = getattr(self,'_opened_file_id',-1)
    # Check if this file already opened.
    if opened_file_id == file_id:
      return self._opened_funit
    # Close any open files before continuing.
    self._close()
    filename = self._files[file_id]
    # Open the file.
    self._opened_file_id = file_id
    self._opened_funit = fnom(filename,FST_RO)
    fstouv(self._opened_funit,FST_RO)
    self._opened_librmn_index = librmn.file_index(self._opened_funit)
    return self._opened_funit

  # Close any currently opened file.
  def _close (self):
    from rpnpy.librmn.base import fclos
    from rpnpy.librmn.fstd98 import fstfrm
    opened_funit = getattr(self,'_opened_funit',-1)
    if opened_funit >= 0:
      fstfrm(opened_funit)
      fclos(opened_funit)
    self._opened_file_id = -1
    self._opened_funit = -1
    self._opened_librmn_index = -1
        

  ###############################################
  # Basic flow for reading data

  def __init__ (self, filename, progress=False, minimal_metadata=False, ignore_typvar=False, ignore_etiket=False, no_quick_scan=False):
    """
    Read raw records from FSTD files, into the buffer.
    Multiple files can be read simultaneously.
    """
    from rpnpy.librmn.fstd98 import fstnbr, fstinl, fstprm, fstopenall
    from rpnpy.librmn.const import FST_RO
    from fstd2nc.extra import maybeFST as isFST
    from collections import Counter
    import numpy as np
    from glob import glob, has_magic
    import os
    import warnings
    from threading import Lock

    # Define a lock for controlling threaded access to the RPN file(s).
    # This is necessary because we can only open a limited number of files
    # at a time, so need to control which ones are opened and available in
    # a thread-safe way.
    self._lock = Lock()

    # Set up a progress bar for scanning the input files.
    Bar = _ProgressBar if progress is True else _FakeBar

    self._minimal_metadata = minimal_metadata
    if not ignore_typvar:
      # Insert typvar value just after nomvar.
      self._var_id = self._var_id[0:1] + ('typvar',) + self._var_id[1:]
      self._human_var_id = self._human_var_id[0:1] + ('%(typvar)s',) + self._human_var_id[1:]
    if not ignore_etiket:
      # Insert etiket value just after nomvar.
      self._var_id = self._var_id[0:1] + ('etiket',) + self._var_id[1:]
      self._human_var_id = self._human_var_id[0:1] + ('%(etiket)s',) + self._human_var_id[1:]

    if isinstance(filename,str):
      infiles = [filename]
    else:
      infiles = list(filename)

    # Apply wildcard and directory expansion to filenames.
    expanded_infiles = []
    for infile in infiles:
      for f in sorted(glob(infile)) or [infile]:
        if os.path.isdir(f):
          for dirpath, dirnames, filenames in os.walk(f,followlinks=True):
            for filename in filenames:
              expanded_infiles.append((infile,os.path.join(dirpath,filename)))
        else:
          expanded_infiles.append((infile,f))

    # Inspect all input files, and extract the headers from valid RPN files.
    matches = Counter()
    headers = []
    self._files = []

    # Show a progress bar when there are multiple input files.
    if len(expanded_infiles) > 1:
      expanded_infiles = Bar(_("Inspecting input files"), suffix='%(percent)d%% (%(index)d/%(max)d)').iter(expanded_infiles)

    for infile, f in expanded_infiles:
      if not os.path.exists(f) or not isFST(f):
        matches[infile] += 0
        continue
      matches[infile] += 1

      # Read the headers from the file(s) and store the info in the table.
      filenum = len(self._files)
      self._files.append(f)
      funit = self._open(filenum)
      nrecs = fstnbr(funit)
      h = np.ma.empty(nrecs, dtype=self._headers_dtype)

      if no_quick_scan:
        keys = fstinl(funit)
        params = map(fstprm, keys)
        for i,prm in enumerate(params):
          for n,v in prm.items():
            if n in h.dtype.names:
              h[n][i] = v
      else:
        from fstd2nc.extra import all_params
        params = all_params(funit,out=h)
        keys = params['key']

      # Encode the keys without the file index info.
      h['key'] = keys
      h['key'] >>= 10
      # The file info will be an index into a separate file list.
      h['file_id'] = filenum

      headers.append(h)

    # Check if the input entries actually matched anything.
    for infile, count in matches.items():
      if count == 0:
        if os.path.isfile(infile):
          warn(_("'%s' is not an RPN standard file.")%infile)
        elif os.path.isdir(infile):
          warn(_("Directory '%s' does not contain any RPN standard files.")%infile)
        elif has_magic(infile):
          warn(_("No RPN standard files match '%s'.")%infile)
        elif not os.path.exists(infile):
          warn(_("'%s' does not exist.")%infile)
        else:
          warn(_("Problem with input file '%s'")%infile)

    nfiles = len(headers)
    if nfiles == 0:
      error(_("no input files found!"))
    info(_("Found %d RPN input file(s)"%nfiles))

    self._headers = np.concatenate(headers)


    # Find all unique meta (coordinate) records, and link a subset of files
    # that provide all unique metadata records.
    # This will make it easier to look up the meta records later.
    meta_mask = np.zeros(len(self._headers),dtype='bool')
    for meta_name in self._meta_records:
      meta_name = (meta_name+'   ')[:4]
      meta_mask |= (self._headers['nomvar'] == meta_name)
    meta_recids = np.where(meta_mask)[0]
    # Use the same unique parameters as regular variables.
    # Plus, ig1,ig2,ig3,ig4.
    # Suppress FutureWarning from numpy about doing this.  Probably benign...
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      meta_keys = self._headers.data[meta_mask][list(self._var_id)+['ig1','ig2','ig3','ig4']]
    meta_keys, ind = np.unique(meta_keys, return_index=True)
    meta_recids = meta_recids[ind]
    # Find the files that give these unique coord records.
    file_ids = sorted(set(self._headers['file_id'][meta_recids]))
    filenames = [self._files[f] for f in file_ids]
    if len(filenames) > 500:
      error(_("Holy crap, how many coordinates do you have???"))
    # If no coordinates found, just open the first file as a dummy file.
    # Less error-prone than checking if _meta_funit is defined every time
    # an FSTD function is called.
    if len(filenames) == 0:
      filenames = self._files[0:1]
    # Open these files and link them together
    self._meta_funit = fstopenall(filenames, FST_RO)


  # Internal iterator.

  # Normal version (without use of pandas).
  def _iter_slow (self):
    from collections import OrderedDict
    import numpy as np
    import warnings

    nrecs = len(self._headers)

    # Degenerate case: no data in buffer
    if nrecs == 0: return

    records = self._headers

    # Ignore deleted / invalidated records.
    deleted = (records['dltf'] == 1)
    if np.any(deleted):
      records = records[~deleted]
    header_indices = np.where(~deleted)[0]

    # Determine the variable identifiers.
    # First, extract the uniquely identifying information from the metadata.
    # Suppress FutureWarning from numpy about doing this.  Probably benign...
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      all_var_ids = records.data[list(self._var_id)]

    # Do a pre-processing step to remove ids that are identical to the one
    # immediately before it.
    # This is purely an optimization thing - the np.unique call later on is
    # O(n log n) so want to prune this array as much as possible beforehand.
    # TODO: Update this once numpy has an unsorted "unique"-like function that
    # can run more efficiently when there are relatively few unique elements.
    # (could get ~O(n) with a hash table).
    var_ids = all_var_ids
    flag = np.concatenate(([True], var_ids[1:] != var_ids[:-1]))
    var_ids = var_ids[flag]

    # Now, find the unique var_ids from this pruned list.
    var_ids = np.unique(var_ids)

    # Loop over each variable and construct the data & metadata.
    for var_id in var_ids:
      selection = (all_var_ids == var_id)
      var_records = records[selection]
      var_record_indices = np.where(selection)[0]
      nomvar = var_id['nomvar'].strip()
      # Ignore coordinate records.
      if nomvar in self._meta_records: continue

      # Get the metadata for each record.
      atts = OrderedDict()
      for n in records.dtype.names:
        if n in self._outer_axes or n in self._ignore_atts: continue
        v = var_records[n]
        # Remove missing values before continuing.
        v = np.ma.compressed(v)
        if len(v) == 0: continue
        # Only use attributes that are consistent across all variable records.
        if len(set(v)) > 1: continue
        v = v[0]
        # Use regular integers for numeric types.
        if np.can_cast(v.dtype,int):
          v = int(v)
        # Trim string attributes (remove whitespace padding).
        if isinstance(v,str): v = v.strip()
        atts[n] = v

      # Get the axis coordinates.
      axes = OrderedDict()
      for n in self._outer_axes:
        values = var_records[n]
        # Remove missing values before continuing.
        values = np.ma.compressed(values)
        # Ignore axes that have no actual coordinate values.
        if len(values) == 0: continue
        # Get all unique values (sorted).
        values = tuple(sorted(set(values)))
        axes[n] = values

      # Construct a multidimensional array to hold the record keys.
      record_id = np.empty(map(len,axes.values()), dtype='int32')

      # Assume missing data (nan) unless filled in later.
      record_id[()] = -1
      
      # Arrange the record keys in the appropriate locations.
      indices = []
      for n in axes.keys():
        u, ind = np.unique(var_records[n], return_inverse=True)
        indices.append(ind)
      record_id[indices] = header_indices[var_record_indices]

      # Check if we have full coverage along all axes.
      have_data = [k >= 0 for k in record_id.flatten()]
      if not np.all(have_data):
        warn (_("Missing some records for %s.")%nomvar)

      # Add dummy axes for the ni,nj,nk record dimensions.
      axes['k'] = tuple(range(var_id['nk']))
      axes['j'] = tuple(range(var_id['nj']))
      axes['i'] = tuple(range(var_id['ni']))

      # Determine the optimal data type to use.
      # First, find unique combos of datyp, nbits
      # (want to minimize calls to dtype_fst2numpy).
      datyp, nbits = zip(*np.unique(var_records.data[['datyp','nbits']]))
      datyp = map(int,datyp)
      nbits = map(int,nbits)
      dtype_list = map(dtype_fst2numpy, datyp, nbits)
      dtype = np.result_type(*dtype_list)

      var = _iter_type( name = nomvar, atts = atts,
                        axes = axes,
                        dtype = dtype,
                        record_id = record_id )
      yield var

      #TODO: Find a minimum set of partial coverages for the data.
      # (e.g., if we have surface-level output for some times, and 3D output
      # for other times).

  # Faster version of iterator (using pandas).
  def _iter_pandas (self):
    from collections import OrderedDict
    import numpy as np
    import pandas as pd
    import warnings

    nrecs = len(self._headers)

    # Degenerate case: no data in buffer
    if nrecs == 0: return

    # Convert records to a pandas DataFrame, which is faster to operate on.
    records = pd.DataFrame.from_records(self._headers)

    # Ignore deleted / invalidated records.
    records = records[records['dltf']==0]

    # Iterate over each variable.
    # Variables are defined by the entries in _var_id.
    for var_id, var_records in records.groupby(self._var_id):
      var_id = OrderedDict(zip(self._var_id, var_id))
      nomvar = var_id['nomvar'].strip()
      # Ignore coordinate records.
      if nomvar in self._meta_records: continue

      # Get the metadata, axes, and corresponding indices of each record.
      atts = OrderedDict()
      axes = OrderedDict()
      indices = []
      for n in records.columns:
        if n in self._ignore_atts: continue
        # Ignore columns which are masked out.
        # https://stackoverflow.com/questions/29530232/python-pandas-check-if-any-value-is-nan-in-dataframe
        if var_records[n].isnull().values.any(): continue
        # Get the unique values, in order.
        cat = pd.Categorical(var_records[n])
        # Is this column a coordinate?
        if n in self._outer_axes:
          axes[n] = tuple(cat.categories)
          indices.append(cat.codes)
        # Otherwise, does it have a consistent value?
        # If so, can add it to the metadata.
        elif len(cat.categories) == 1:
          v = cat[0]
          # Trim string attributes (remove whitespace padding).
          if isinstance(v,str): v = v.strip()
          # Use regular integers for numeric types.
          elif np.can_cast(v.dtype,int):
            v = int(v)
          atts[n] = v

      # Construct a multidimensional array to hold the record keys.
      record_id = np.empty(map(len,axes.values()), dtype='int32')

      # Assume missing data (nan) unless filled in later.
      record_id[()] = -1

      # Arrange the record keys in the appropriate locations.
      record_id[indices] = var_records.index

      # Check if we have full coverage along all axes.
      have_data = [k >= 0 for k in record_id.flatten()]
      if not np.all(have_data):
        warn (_("Missing some records for %s.")%nomvar)

      # Add dummy axes for the ni,nj,nk record dimensions.
      axes['k'] = tuple(range(var_id['nk']))
      axes['j'] = tuple(range(var_id['nj']))
      axes['i'] = tuple(range(var_id['ni']))

      # Determine the optimal data type to use.
      # First, find unique combos of datyp, nbits
      # (want to minimize calls to dtype_fst2numpy).
      x = var_records[['datyp','nbits']].drop_duplicates()
      datyp = map(int,x['datyp'])
      nbits = map(int,x['nbits'])
      dtype_list = map(dtype_fst2numpy, datyp, nbits)
      dtype = np.result_type(*dtype_list)

      var = _iter_type( name = nomvar, atts = atts,
                        axes = axes,
                        dtype = dtype,
                        record_id = record_id )
      yield var

      #TODO: Find a minimum set of partial coverages for the data.
      # (e.g., if we have surface-level output for some times, and 3D output
      # for other times).

  # Choose which method to iterate over the data
  # (depending on if pandas is installed).
  def _iter (self):
    try:
      import pandas as pd
      iterfunc = self._iter_pandas
    except ImportError:
      iterfunc = self._iter_slow
    for var in iterfunc():
      yield var

  # How to read the data.
  # Override fstluk so it takes a record index instead of a key/handle.
  def _fstluk (self, rec_id, dtype=None, rank=None, dataArray=None):
    from rpnpy.librmn.fstd98 import fstluk
    # Use local overrides for dtype.
    if dtype is None and isinstance(rec_id,dict):
      if 'datyp' in rec_id and 'nbits' in rec_id:
        dtype = dtype_fst2numpy(rec_id['datyp'],rec_id['nbits'])

    # Lock the opened file so another thread doesn't close it.
    with self._lock:
      if isinstance(rec_id,dict):
        # If we were given a dict, assume it's for a valid open file.
        key = rec_id['key']
      else:
        # Otherwise, need to open the file and get the proper key.
        file_id = self._headers['file_id'][rec_id]
        self._open(file_id)
        key = int(self._headers['key'][rec_id]<<10) + self._opened_librmn_index

      return fstluk (key, dtype, rank, dataArray)

  #
  ###############################################



