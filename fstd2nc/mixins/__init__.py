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


#####################################################################
# Mixins for different features / behaviour for the conversions.

# This module contains the base class for the mixins, from which all others
# are derived.


from fstd2nc.stdout import _, info, warn, error

# Indicate if pandas should be used or not.
_pandas_needed = False
def _use_pandas ():
  if not _pandas_needed: return False
  try:
    import pandas
    return True
  except ImportError:
    return False

# Decorator for efficiently converting a scalar function to a vectorized
# function.
def vectorize (f):
  from functools import wraps
  import numpy as np
  def cached_f(x, cache={}):
    if x not in cache:
      cache[x] = f(x)
    return cache[x]
  @wraps(f)
  def vectorized_f (x):
    mask = None
    if hasattr(x,'mask'):
      mask = np.ma.getmaskarray(x)
      n = len(x)
      x = x[~mask]
    if _use_pandas():
      from pandas import Series, unique
      # If we're given a scalar value, then simply return it.
      if not hasattr(x,'__len__'):
        return f(x)
      # Get unique values
      x = np.array(x,copy=True)
      inputs = unique(x)
      outputs = map(f,inputs)
      table = dict(zip(inputs,outputs))
      result = Series(x).map(table)
      result = result.values
    else:
      # If we're given a scalar value, then simply return it.
      if not hasattr(x,'__len__'):
        return cached_f(x)
      result = list(map(cached_f,x))
    if mask is not None:
      result = np.asarray(result)
      masked_result = np.ma.masked_all(n,dtype=result.dtype)
      masked_result[~mask] = result
      result = masked_result
    return result
  return vectorized_f


class _base_type (object):
  @property
  def shape(self):
    return tuple(map(len,self.axes))
  @property
  def dims(self):
    return tuple(a.name for a in self.axes)
  def getaxis(self,name):
    for a in self.axes:
      if a.name == name:
        return a
    return None
  def __iter__ (self):
    return (getattr(self,n) for n in self.__slots__)

# Data that's loaded in memory.
class _var_type (_base_type):
  __slots__ = ('name','atts','axes','array','deps')
  def __init__ (self, name, atts, axes, array):
    self.name = name
    self.atts = atts
    self.axes = list(axes)
    self.array = array
    self.deps = []

# Axis data.
class _axis_type (_base_type):
  __slots__ = ('name','atts','axes','array','deps')
  def __init__ (self, name, atts, array):
    self.name = name
    self.atts = atts
    self.axes = [self]
    self.array = array
    self.deps = []
  def __len__ (self):
    return len(self.array)

# Dimension (without coordinate values).
class _dim_type (object):
  __slots__ = ('name','length','deps')
  def __init__ (self, name, length):
    self.name = name
    self.length = length
    self.deps = []
  def __len__ (self):
    return self.length

# Data that resides in FST records on disk.
class _iter_type (_base_type):
  __slots__ = ('name','atts','axes','dtype','record_id','deps')
  def __init__ (self, name, atts, axes, dtype, record_id):
    self.name = name
    self.atts = atts
    self.axes = axes
    self.dtype = dtype
    self.record_id = record_id
    self.deps = []

# Similar to above, but more general.
# Can allow for multiple records along the inner axes.
class _chunk_type (_base_type):
  __slots__ = ('name','atts','axes','dtype','chunks','record_id','deps')
  def __init__ (self, name, atts, axes, dtype, chunks, record_id):
    self.name = name
    self.atts = atts
    self.axes = axes
    self.dtype = dtype
    self.chunks = chunks
    self.record_id = record_id
    self.deps = []

# Fake progress bar - does nothing.
class _FakeBar (object):
  def __init__ (self, *args, **kwargs): pass
  def iter(self, it):
    for i in it: yield i
  def __next__(self): pass
  def next(self): pass
  def finish(self): pass

# Try importing progress module.
try:
  import sys
  if sys.stdout.encoding is not None and sys.stdout.encoding.lower().startswith('utf'):
    from progress.bar import IncrementalBar
  # Handle environments that can't print unicode characters.
  else:
    from progress.bar import Bar as IncrementalBar
  del sys
  class _ProgressBar(IncrementalBar):
    # Also display (even if not a tty).
    # Could be a jupyter notebook, which can still benefit from a progress bar.
    check_tty = False
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

# Special version of dictionary that keeps track of unaccessed keys.
class AccessCountDict (object):
  def __init__ (self, **items):
    from collections import Counter
    self._dict = dict(items)
    self._counter = Counter()
  def __getitem__ (self, key):
    self._counter[key] += 1
    return self._dict[key]
  def __setitem__ (self, key, value):
    self._counter[key] += 1
    self._dict[key] = value
  def __getattr__ (self, name): return getattr(self._dict, name)
  def __str__ (self): return self._dict.__str__()
  def __repr__ (self): return self._dict.__repr__()

def _expand_files (filename):
  """
  Expands the given filename/directory/glob/Path/etc. into an explicit list
  of matching files.

  Parameters
  ----------
  filename : str or Path or iterable
      The thing to expand.

  Returns
  -------
  List of matches, each match a tuple of (pattern,[files])
  """
  from pathlib import Path
  from glob import glob
  import os
  if isinstance(filename,(str,Path)):
    infiles = [filename]
  else:
    infiles = list(filename)

  # Apply wildcard and directory expansion to filenames.
  expanded_infiles = []
  for infile in infiles:
    if isinstance(infile,Path):
      infile = str(infile)
    for f in sorted(glob(infile)) or [infile]:
      if os.path.isdir(f):
        for dirpath, dirnames, filenames in os.walk(f,followlinks=True):
          for filename in filenames:
            expanded_infiles.append((infile,os.path.join(dirpath,filename)))
      else:
        expanded_infiles.append((infile,f))
  return expanded_infiles


# Define a class for encoding / decoding the data.
# Each step is placed in its own "mixin" class, to make it easier to patch in 
# new behaviour if more exotic files are encountered in the future.
class BufferBase (object):

  # Names of records that should be kept separate (never grouped into
  # multidimensional arrays).
  @classmethod
  def _meta_records(cls):
    return ()
  # Other records that should also be kept separate, but only if they are
  # 1-dimensional.  If they are 2D, then they should be processed as a normal
  # variable.
  @classmethod
  def _maybe_meta_records(cls):
    return ()

  # Attributes which could potentially be used as outer axes.
  # The values from the attribute will become the axis values.
  _outer_axes = ()

  # Attributes which could be used as inner axes.
  # The values represent the *length* of the axis, which will have no
  # coordinate values by default.
  _inner_axes = ()

  # Attributes which could be used as auxiliary coordinates for the outer
  # axes.  The dictionary keys are the outer axis names, and the values are
  # a list of columns which can act as coordinates.
  from collections import OrderedDict
  _outer_coords = OrderedDict()
  del OrderedDict

  # Attributes which uniquely identify a variable.
  _var_id = ()

  # Similar to above, but a human-readable version of the id.
  # Could be used for appending suffixes to variable names to make them unique.
  # Uses string formatting operations on the variable metadata.
  _human_var_id = ()

  # Record parameters which should not be used as nc variable attributes.
  # (They're either internal to the file, or part of the data, not metadata).
  _ignore_atts = ('file_id','name','address','length','dtype','selected')

  # Header columns which contain extra data that should be passed to _decode.
  # Tuple of (offset, length, d) where:
  # - offset is the address in the file where the data is located.
  # - length is how much data to read.
  # - d is an alternate source (dask graph?) for case where data is not
  #   stored in a file.
  _decoder_data = (('data',('address','length','d')),)

  # Extra arguments to pull from columns of the table.
  _decoder_extra_args = ()

  # Extra (scalar) arguments needed for the decoder at runtime.
  def _decoder_scalar_args (self):
    return {}

  # Define any command-line arguments for reading FSTD files.
  @classmethod
  def _cmdline_args (cls, parser):
    from fstd2nc import __version__
    from argparse import SUPPRESS
    from sys import stdout
    parser.add_argument('--version', action='version', version=__version__, help=_("show program's version number and exit"))
    group = parser.add_mutually_exclusive_group()
    _('Display a progress bar during the conversion, if the "progress" module is installed.')
    group.add_argument('--progress', action='store_true', default=stdout.isatty(), help=SUPPRESS)
    group.add_argument('--no-progress', action='store_false', dest='progress', help=_('Disable the progress bar.'))
    parser.add_argument('--serial', action='store_true', help=_('Disables multithreading/multiprocessing.  Useful for resource-limited machines.'))
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--minimal-metadata', action='store_true', default=True, help=_("Don't include internal record attributes and other internal information in the output metadata.")+" "+_("This is the default behaviour."))
    group.add_argument('--internal-metadata','--rpnstd-metadata', action='store_false', dest='minimal_metadata', help=_("Include all internal record attributes in the output metadata."))
    group.add_argument('--metadata-list','--rpnstd-metadata-list', metavar='nomvar,...', help=_("Specify a minimal set of internal record attributes to include in the output file."))

  # Do some checks on the command-line arguments after parsing them.
  @classmethod
  def _check_args (cls, parser, args):
    return  # Nothing to check here.

  # Extract metadata from a particular header.
  def _get_header_atts (self, header):
    for n,v in header.items():
      if n in self._ignore_atts: continue
      # Python 3: convert bytes to str
      if isinstance(v,bytes):
        v = str(v.decode())
      if isinstance(v,str):
        v = v.strip()
      yield (n,v)

  # Control the pickling / unpickling of BufferBase objects.
  def __getstate__ (self):
    state = self.__dict__.copy()
    return state
  def __setstate__ (self, state):
    self.__dict__.update(state)


  ###############################################
  # Basic flow for reading data

  def __init__ (self, filename, progress=False, serial=False, **kwargs):
    """
    Read raw records from FSTD files, into the buffer.
    Multiple files can be read simultaneously.

    Parameters
    ----------
    filename : str or list or existing Buffer object
        The input files to convert.
    progress : bool, optional
        Display a progress bar during the conversion, if the "progress"
        module is installed.
    serial : bool, optional
        Disables multithreading/multiprocessing.  Useful for resource-limited
        machines.
    internal_metadata : bool, optional
        Include all internal record attributes in the output metadata.
    metadata_list : str or list, optional
        Specify a minimal set of internal record attributes to include in the
        output file.
    """
    from collections import Counter
    import numpy as np
    from glob import has_magic
    import os
    from multiprocessing import Pool
    try:
      from itertools import imap  # Python 2
    except ImportError:
      imap = map

    # Set up a progress bar for scanning the input files.
    Bar = _ProgressBar if progress is True else _FakeBar

    self._serial = serial

    # Detect if an existing Buffer object was provided.
    if hasattr(filename,'_headers') and hasattr(filename,'_files'):
      existing_buffer = filename
      filename = []
    else:
      existing_buffer = None

    expanded_infiles = _expand_files (filename)

    # How to handle internal metadata.
    internal_metadata = kwargs.pop('internal_metadata',None)
    # Check for legacy rpnstd_metadata argument
    if internal_metadata is None:
      internal_metadata = kwargs.pop('rpnstd_metadata',None)
    metadata_list = kwargs.pop('metadata_list',None)
    # Check for legacy rpnstd_metadata_list argument
    if metadata_list is None:
      metadata_list = kwargs.pop('rpnstd_metadata_list',None)
    minimal_metadata = kwargs.pop('minimal_metadata',None)

    # Set default for minimal_metadata
    if internal_metadata is not None:
      minimal_metadata = not internal_metadata
    if minimal_metadata is None:
      minimal_metadata = True
    # Set default for metadata_list
    if minimal_metadata is True and metadata_list is None:
      metadata_list = ''
    if isinstance(metadata_list,str):
      metadata_list = metadata_list.replace(',',' ')
      metadata_list = metadata_list.split()
    if hasattr(metadata_list,'__len__'):
      metadata_list = tuple(metadata_list)
    self._metadata_list = metadata_list

    # Should not have any unprocessed keyword arguments after this.
    if len(kwargs) > 0:
      error(_("Unexpected arguments: %s"%(kwargs.keys())))

    # Extract headers from the files.
    if hasattr(progress,'next'):
      bar = progress
    else:
      if len(expanded_infiles) == 1: Bar = _FakeBar
      bar = Bar(_("Inspecting input files"), suffix='%(percent)d%% (%(index)d/%(max)d)', max=len(expanded_infiles))

    if len(expanded_infiles) > 1 and not self._serial:
      with Pool() as p:
        headers = p.imap (self._raw_headers, [f for (infile,f) in expanded_infiles])
        headers = ([h,bar.next()][0] for h in headers)
        headers = list(headers) # Start scanning.
    else:
      headers = imap (self._raw_headers, [f for (infile,f) in expanded_infiles])
      headers = ([h,bar.next()][0] for h in headers)
      headers = list(headers) # Start scanning.

    # If we made this progress bar, then finish it here.
    # Otherwise, if it was passed in it could be part of a bigger workflow
    # so don't finish it yet.
    if not hasattr(progress,'next'):
      bar.finish()

    # Check which files had headers used, report on the results.
    matches = Counter()
    self._files = []
    file_ids = []
    for i, (infile, f) in enumerate(expanded_infiles):
      if headers[i] is not None:
        matches[infile] += 1
        filenum = len(self._files)
        self._files.append(f)
        file_ids.extend([filenum]*len(headers[i]))
      else:
        matches[infile] += 0

    # Decode all the headers
    headers = [h for h in headers if h is not None and len(h) > 0]
    # Remember indices in the headers (needed for reconstructing keys)
    indices = [range(len(h)) for h in headers]
    if len(headers) > 0:
      headers = np.concatenate(headers)
      indices = np.concatenate(indices)
      headers = self._decode_headers(headers)
      # Encode the keys without the file index info.
      headers['file_id'] = np.array(file_ids, dtype='int32')
      headers['selected'] = np.ones(len(file_ids),'bool')

    # Check if the input entries actually matched anything.
    for infile, count in matches.items():
      if count == 0:
        if os.path.isfile(infile):
          warn(_("'%s' is not %s.")%(infile,self._format_singular))
        elif os.path.isdir(infile):
          warn(_("Directory '%s' does not contain any %s.")%(infile,self._format_plural))
        elif has_magic(infile):
          warn(_("No %s match '%s'.")%(self._format_plural,infile))
        elif not os.path.exists(infile):
          warn(_("'%s' does not exist.")%infile)
        else:
          warn(_("Problem with input file '%s'")%infile)

    # Use existing Buffer object data, if that was provided instead of files.
    if existing_buffer is not None:
      self._files = list(existing_buffer._files)
      headers = {k:v.copy() for k,v in existing_buffer._headers.items()}

    nfiles = len(self._files)
    if nfiles == 0 and existing_buffer is None:
      error(_("no input files found!"))
    elif nfiles > 10:
      global _pandas_needed
      _pandas_needed = True
    info(_("Found %d %s"%(nfiles,self._format_plural)))

    self._headers = headers
    self._nrecs = max(len(self._headers[key]) for key in self._headers.keys())


  # Generate structured variables from the data records.

  # Normal version (without use of pandas).
  def _makevars_slow (self):
    from collections import OrderedDict
    import numpy as np
    import warnings
    from fstd2nc.extra import structured_array

    nrecs = self._nrecs

    # Degenerate case: no data in buffer
    if nrecs == 0: return

    records = structured_array(self._headers)

    # Ignore deleted / invalidated records, and coordinate records.
    valid = (records['selected'] == True)
    records = records[valid]
    header_indices = np.where(valid)[0]

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
    if len(var_ids) > 0:
      flag = np.concatenate(([True], var_ids[1:] != var_ids[:-1]))
      var_ids = var_ids[flag]

    # Now, find the unique var_ids from this pruned list.
    var_ids = np.unique(var_ids)

    # Keep track of axes that were generated
    known_axes = dict()

    # Keep track of any auxiliary coordinates that were generated.
    known_coords = dict()

    # Loop over each variable and construct the data & metadata.
    self._varlist = []
    for var_id in var_ids:
      selection = (all_var_ids == var_id)
      var_records = records[selection]
      var_record_indices = np.where(selection)[0]
      nomvar = var_id['name'].strip()
      nomvar = str(nomvar.decode()) # Python3: convert bytes to str.

      # Get the metadata for each record.
      atts = OrderedDict()
      for n in records.dtype.names:
        if n in self._outer_axes or n in self._outer_coords or n in self._ignore_atts:
          continue
        v = var_records[n]
        # Remove missing values before continuing.
        v = np.ma.compressed(v)
        if len(v) == 0: continue
        # Only use attributes that are consistent across all variable records.
        if len(set(v)) > 1: continue
        v = v[0]
        # Python3: convert bytes to str.
        if isinstance(v,bytes): v = str(v.decode())
        # Trim string attributes (remove whitespace padding).
        if isinstance(v,str): v = v.strip()
        # Use regular integers for numeric types.
        elif np.can_cast(v.dtype,int):
          v = int(v)
        atts[n] = v

      # Get the axes.
      axes = OrderedDict()
      for n in self._outer_axes:
        values = var_records[n]
        # Ignore axes that have no actual coordinate values.
        if len(np.ma.compressed(values)) == 0: continue
        # Remove missing values before continuing.
        values = np.ma.unique(values).data
        # Get all unique values (sorted).
        values = tuple(sorted(set(values)))
        if (n,values) not in known_axes:
          known_axes[(n,values)] = _axis_type(name = n, atts = OrderedDict(),
                                              array = np.array(values))
        axes[n] = known_axes[(n,values)]

      # Construct a multidimensional array to hold the record keys.
      record_id = np.empty(list(map(len,axes.values())), dtype='int32')

      # Assume missing data (nan) unless filled in later.
      record_id[()] = -1

      # Arrange the record keys in the appropriate locations.
      indices = []
      for n in axes.keys():
        u, ind = np.unique(var_records[n], return_inverse=True)
        indices.append(ind)
      record_id[tuple(indices)] = header_indices[var_record_indices]

      # Get the auxiliary coordinates.
      coords = []
      for n, coordaxes in self._outer_coords.items():
        # Get the axes for this coordinate.
        # Use the same order of columns as was used for the outer axes,
        # so we get the right coordinate order after sorting.
        coordaxes = OrderedDict((k,v) for k,v in axes.items() if k in coordaxes)
        # Sanity check - do we actually have any of the coordinate axes?
        if len(coordaxes) == 0: continue
        # Unique key for this coordinate
        key = (n,tuple(coordaxes.items()))
        # Arrange the coordinate values in the appropriate location.
        shape = list(map(len,list(coordaxes.values())))
        # Extract all values of the coordinate (including duplicates over
        # other axes).  Will determine which values to put in what order later.
        all_coord_values = np.ma.compressed(var_records[n])
        if len(all_coord_values) == 0: continue
        values = np.zeros(shape,dtype=all_coord_values.dtype)
        indices = []
        for k in coordaxes.keys():
          u, ind = np.unique(var_records[k], return_inverse=True)
          indices.append(ind)
        values[tuple(indices)] = all_coord_values
        if key not in known_coords:
          coord = _var_type (name = n, atts = OrderedDict(),
                           axes = list(coordaxes.values()), array = values )
          known_coords[key] = coord
        coords.append(known_coords[key])

      # Check if we have full coverage along all axes.
      have_data = (record_id.flatten() >= 0)
      if not np.all(have_data):
        warn (_("Missing some records for %s.")%nomvar)

      # Add dummy axes for the inner dimensions.
      for n in self._inner_axes:
        values = var_records[n]
        # Remove missing values before continuing.
        values = np.ma.compressed(values)
        # Ignore dimensions that are masked out or are inconsistent.
        # Note: inner dimensions should never be inconsistent.  They should
        # be encoded into the var_id fields to ensure that.
        if len(values) == 0: continue
        length = values[0]
        axes[n] = _dim_type(name=n, length = length)

      # Determine the optimal data type to use.
      dtype_list = map(np.dtype,var_records['dtype'])
      dtype = np.result_type(*dtype_list)

      var = _iter_type( name = nomvar, atts = atts,
                        axes = list(axes.values()),
                        dtype = dtype,
                        record_id = record_id )
      # Add auxiliary coordinate variables through silent dependency.
      # Don't want these in the 'coordinates' attribute because they seem to
      # cause issues with certain netCDF decoders.
      if len(coords) > 0:
        var.deps = coords
      self._varlist.append(var)

      #TODO: Find a minimum set of partial coverages for the data.
      # (e.g., if we have surface-level output for some times, and 3D output
      # for other times).

  # Faster version of iterator (using pandas).
  def _makevars_pandas (self):
    from collections import OrderedDict
    import numpy as np
    import pandas as pd
    import warnings

    nrecs = self._nrecs

    # Degenerate case: no data in buffer
    if nrecs == 0: return

    # Convert records to a pandas DataFrame, which is faster to operate on.
    records = pd.DataFrame.from_dict(self._headers)
    # Keep track of original dtypes (may need to re-cast later).
    original_dtypes = dict([(key,value.dtype) for key,value in self._headers.items()])

    # Get the subset of records destined for output.
    records = records[records['selected']==True]

    # Keep track of any axes that were generated.
    known_axes = dict()

    # Keep track of any auxiliary coordinates that were generated.
    known_coords = dict()

    # Iterate over each variable.
    # Variables are defined by the entries in _var_id.
    self._varlist = []
    for var_id, var_records in records.groupby(list(self._var_id)):
      var_id = OrderedDict(zip(self._var_id, var_id))
      nomvar = var_id['name'].strip()
      nomvar = str(nomvar.decode()) # Python3: convert bytes to str.

      # Get the attributes, axes, and corresponding indices of each record.
      atts = OrderedDict()
      axes = OrderedDict()
      indices = OrderedDict()
      coordnames = []
      coord_axes = OrderedDict()
      coord_indices = OrderedDict()
      for n in records.columns:
        if n in self._ignore_atts: continue
        # Ignore columns which are masked out.
        # https://stackoverflow.com/questions/29530232/python-pandas-check-if-any-value-is-nan-in-dataframe
        try:
          if np.isnan(var_records[n].values[0]): continue
        except TypeError:
          if var_records[n].values[0] is None: continue
        # Get the unique values, in order.
        # Coerce back to original dtype, since masked columns get upcasted to
        # float64 in pandas.DataFrame.from_records.
        try:
          column = var_records[n].astype(original_dtypes[n])
        except TypeError:
          # Some types may not be re-castable.
          # For instance, pandas < 0.23 can't convert between datetime64 with
          # different increments ([ns] and [s]).
          column = var_records[n]
        # Is this column an outer axis?
        if n in self._outer_axes:
          values, indices[n] = np.unique(column, return_inverse=True)
          values = tuple(values)
          if (n,values) not in known_axes:
            known_axes[(n,values)] = _axis_type(name = n, atts = OrderedDict(),
                                   array = np.array(values,dtype=column.dtype))
          axes[n] = known_axes[(n,values)]
          # Is this also an axis for an auxiliary coordinate?
          for coordname,coordaxes in self._outer_coords.items():
            if n in coordaxes:
              coordnames.append(coordname)
              coord_axes.setdefault(coordname,OrderedDict())[n] = axes[n]
              coord_indices.setdefault(coordname,OrderedDict())[n] = indices[n]
        # Otherwise, does it have a consistent value?
        # If so, can add it to the metadata.
        # Ignore outer coords, since the value is already encoded elsewhere.
        elif np.all(column.values == column.values[0]) and n not in self._outer_coords:
          try:
            v = column.values[0]
            # Python3: convert bytes to str.
            if isinstance(v,bytes): v = str(v.decode())
            # Trim string attributes (remove whitespace padding).
            if isinstance(v,str): v = v.strip()
            # Use regular integers for numeric types.
            elif np.can_cast(v,int):
              v = int(v)
            atts[n] = v
          except (TypeError,ValueError):
            pass

      # Recover the proper order for the outer axes.
      # Not necessarily the order of the header table columns.
      axes = OrderedDict((n,axes[n]) for n in self._outer_axes if n in axes)
      indices = tuple([indices[n] for n in self._outer_axes if n in indices])
      for coordname in coord_axes.keys():
        coord_axes[coordname] = OrderedDict((n,coord_axes[coordname][n]) for n in self._outer_axes if n in coord_axes[coordname])
        coord_indices[coordname] = [coord_indices[coordname][n] for n in self._outer_axes if n in coord_indices[coordname]]

      # Construct a multidimensional array to hold the record keys.
      record_id = np.empty(list(map(len,axes.values())), dtype='int32')

      # Assume missing data (nan) unless filled in later.
      record_id[()] = -1

      # Arrange the record keys in the appropriate locations.
      record_id[indices] = var_records.index

      # Get the auxiliary coordinates.
      coords = []
      for n in coordnames:
        # Ignore auxiliary coordinates which are masked out.
        try:
          if np.isnan(var_records[n].values[0]): continue
        except TypeError:
          if var_records[n].values[0] is None: continue
        # Unique key for this coordinate
        key = (n,tuple(coord_axes[n].items()))
        # Arrange the coordinate values in the appropriate location.
        shape = list(map(len,list(coord_axes[n].values())))
        values = np.zeros(shape,dtype=var_records[n].dtype)
        indices = tuple(coord_indices[n])
        values[indices] = var_records[n]
        if key not in known_coords:
          coord = _var_type (name = n, atts = OrderedDict(),
                           axes = list(coord_axes[n].values()),
                           array = values )
          known_coords[key] = coord
        coords.append(known_coords[key])



      # Check if we have full coverage along all axes.
      have_data = (record_id.flatten() >= 0)
      if not np.all(have_data):
        warn (_("Missing some records for %s.")%nomvar)

      # Add dummy axes for the inner dimensions.
      for n in self._inner_axes:
        # [Note: this code copied from above for dealing with masked data
        #  in pandas DataFrames]
        #
        # Ignore dimensions that are masked out or are inconsistent.
        # Note: inner dimensions should never be inconsistent.  They should
        # be encoded into the var_id fields to ensure that.
        try:
          if np.isnan(var_records[n].values[0]): continue
        except TypeError:
          if var_records[n].values[0] is None: continue
        # Get the unique values, in order.
        column = var_records[n].astype(original_dtypes[n])
        ###
        values = column.values
        # Ignore dimensions that are masked out or are inconsistent.
        # Note: inner dimensions should never be inconsistent.  They should
        # be encoded into the var_id fields to ensure that.
        if len(values) == 0: continue
        length = values[0]
        axes[n] = _dim_type(name=n, length = length)

      # Determine the optimal data type to use.
      dtype_list = map(np.dtype,var_records['dtype'])
      dtype = np.result_type(*dtype_list)

      var = _iter_type( name = nomvar, atts = atts,
                        axes = list(axes.values()),
                        dtype = dtype,
                        record_id = record_id )
      # Add auxiliary coordinate variables through silent dependency.
      # Don't want these in the 'coordinates' attribute because they seem to
      # cause issues with certain netCDF decoders.
      if len(coords) > 0:
        var.deps = coords
      self._varlist.append(var)

      #TODO: Find a minimum set of partial coverages for the data.
      # (e.g., if we have surface-level output for some times, and 3D output
      # for other times).

  # Choose which method to iterate over the data
  # (depending on if pandas is installed).
  def _makevars (self):
    if _use_pandas():
      self._makevars_pandas()
    else:
      self._makevars_slow()

  # Decode variable structures into a format-independent structure.
  def _serialize (self):
    self._makevars()
    return self._varlist

  # Re-encode format-independent info into a format-specific table.
  @classmethod
  def _deserialize (cls, data):
    # Encapsulate this info in a structure.
    # Adapted from 'fake_buffer' object created in from_fstpy method.
    # TODO: simpler way of making empty buffer?
    b = cls.__new__(cls)
    b._files = [None]
    b._varlist = data
    b._unmakevars()
    # Check for unused columns.
    for colname in b._headers.keys():
      if b._headers._counter[colname] == 0:
        warn(_("Unknown axis '%s'.  Encoding not complete.")%colname)
    # Convert headers to regular dictionary.
    b._headers = dict(b._headers)
    # Return the result.
    return b


  # Inverse operation - determine record structure from a varlist.
  def _unmakevars (self):
    from fstd2nc.mixins import _var_type, _iter_type, _dim_type, _axis_type
    import numpy as np
    # By the time we get to this base method, the mixins should have filtered
    # out all variables that shouldn't be in the record table.
    varlist = []
    for var in self._varlist:
      if isinstance(var,_iter_type):
        varlist.append(var)
      else:
        warn (_("Unable to encode %s.")%var.name)
    if len(varlist) == 0:
      error (_("Nothing to encode!"))
    self._varlist = varlist
    self._nrecs = sum(np.product(var.record_id.shape) if var.record_id.ndim > 0 else 1 for var in varlist)
    # Keep track of columns that we don't need to worry about.
    # (don't care if they get used or not).
    dontcare = set()
    # Start constructing header information.
    # The initial columns will be the outer axes of the variables.
    # plus the name and data columns, and inner dimensions.
    headers = dict()
    for var in self._varlist:
      # Add outer axes.
      for axis in var.axes[:var.record_id.ndim]:
        if axis.name in headers: continue
        headers[axis.name] = np.ma.masked_all(self._nrecs, dtype=axis.array.dtype)
      # Add coordinates as well.
      for coord in var.atts.get('coordinates',[]):
        if coord.name in headers: continue
        # It's ok if these values are never actually referenced.
        dontcare.add(coord.name)
        headers[coord.name] = np.ma.masked_all(self._nrecs, dtype=coord.array.dtype)
      # Add inner axes.
      for dimname in self._inner_axes:
        headers[dimname] = np.ma.ones(self._nrecs, dtype='int32')
      for axis in var.axes[var.record_id.ndim:]:
        if axis.name not in self._inner_axes:
          error (_("Unhandled inner axis %s.")%axis.name)
      # Add extra metadata.
      for attname, attval in var.atts.items():
        if attname in headers: continue
        # It's ok if these values are never actually referenced.
        dontcare.add(attname)
        # For simple structures, create array with appropriate dtype.
        # For complicated structures, use object dtype.
        sample = np.array(attval)
        if sample.ndim == 0:
          headers[attname] = np.ma.masked_all(self._nrecs, dtype=sample.dtype)
        else:
          headers[attname] = np.ma.masked_all(self._nrecs, dtype=object)

    namesize = max(len(var.name) for var in varlist)
    headers['name'] = np.ma.masked_all(self._nrecs, dtype='|S%d'%namesize)
    headers['d'] = np.empty(self._nrecs, dtype=object)
    headers['file_id'] = np.empty(self._nrecs,dtype='int32')
    headers['file_id'][:] = -1
    offset = 0
    for var in varlist:
      axes = var.axes[:var.record_id.ndim]
      # Add outer axes.
      for i,axis in enumerate(axes):
        shape = [1]*len(axes)
        shape[i] = len(axis)
        array = np.empty(var.record_id.shape, dtype=axis.array.dtype)
        array[()] = axis.array.reshape(shape)
        headers[axis.name][offset:offset+var.record_id.size] = array.flatten()
      # Add coordinate info as well (only if defined on outer axes).
      for coord in var.atts.get('coordinates',[]):
        if not all(axis in axes for axis in coord.axes): continue
        shape = [len(axis) if axis in coord.axes else 1 for axis in axes]
        array = np.empty(var.record_id.shape, dtype=coord.array.dtype)
        array[()] = coord.array.reshape(shape)
        headers[coord.name][offset:offset+var.record_id.size] = array.flatten()
      # Add inner axes.
      for axis in var.axes[var.record_id.ndim:]:
        headers[axis.name][offset:offset+var.record_id.size] = len(axis)
      # Add extra metadata.
      for attname, attval in var.atts.items():
        if attname in ('coordinates',): continue
        headers[attname][offset:offset+var.record_id.size] = attval
      headers['name'][offset:offset+var.record_id.size] = var.name
      headers['d'][offset:offset+var.record_id.size] = var.record_id.flatten()
      offset = offset + var.record_id.size

    # Use special dictionary to keep track of which columns were used.
    # It would be a problem if the columns were never processed!
    self._headers = AccessCountDict(**headers)
    # Pre-access certain columns, so they don't get flagged as unhandled.
    # (columns that aren't important).
    dontcare |= set(('name','d','file_id'))
    for col in dontcare:
      self._headers[col]
    # From here, the mixins will start processing the header info

  # Iterate over all unique axes found in the variables.
  # Requires _makevars() to have already been called.
  def _iter_axes (self, name=None, varlist=False):
    from collections import OrderedDict
    if not varlist:
      handled = set()
      for var in self._iter_objects():
        if not hasattr(var,'axes'): continue
        for axis in var.axes:
          if name is not None and axis.name != name: continue
          if id(axis) in handled: continue
          yield axis
          handled.add(id(axis))
    else:
      id_lookup = dict()
      output = OrderedDict()
      for var in self._iter_objects():
        if not hasattr(var,'axes'): continue
        for axis in var.axes:
          if name is not None and axis.name != name: continue
          if id(axis) not in id_lookup:
            id_lookup[id(axis)] = axis
          output.setdefault(id(axis),[]).append(var)
      for axis_id,varlist in output.items():
        yield id_lookup[axis_id],varlist


  # Iterate over all unique coordinates found in the variables.
  # Requires _makevars() to have already been called.
  def _iter_coords (self):
    handled = set()
    for var in self._varlist:
      for coord in var.atts.get('coordinates',[]):
        if id(coord) in handled: continue
        yield coord
        handled.add(id(coord))

  # Iterate over all data objects.
  # Requires _makevars() to have already been called.
  def _iter_objects (self, obj=None, handled=None):
    if obj is None:
      obj = self._varlist
    if handled is None:
      handled = set()

    if id(obj) in handled:
      return

    if isinstance(obj,(_iter_type,_chunk_type,_var_type,_axis_type,_dim_type)):
      yield obj
      handled.add(id(obj))

    if isinstance(obj,list):
      for var in obj:
        for o in self._iter_objects(var,handled):
          yield o
      return

    if isinstance(obj,dict):
      for key,value in obj.items():
        for o in self._iter_objects(key,handled):
          yield o
        for o in self._iter_objects(value,handled):
          yield o
      return

    if hasattr(obj,'axes'):
      for o in self._iter_objects(obj.axes,handled):
        yield o
    if hasattr(obj,'atts'):
      for o in self._iter_objects(obj.atts,handled):
        yield o
    if hasattr(obj,'deps'):
      for o in self._iter_objects(obj.deps,handled):
        yield o

  # How to decode the data from a raw binary array.
  @staticmethod
  def _decode (data):
    raise NotImplementedError("No decoder found.")

  # Stub for any extra processing needed for the data after it's decoded.
  @classmethod
  def _postproc (cls, data):
    return data

  # Shortcuts to header decoding functions.
  # Put into the class so they can potentially be overridden for other formats.
  @staticmethod
  def _decode_headers (headers):
    raise NotImplementedError("No decoder found.")
  @staticmethod
  def _raw_headers (filename):
    raise NotImplementedError("No decoder found.")

  # Shortcut for reading a record, given a record id.
  def _read_record (self, rec):
    import numpy as np
    kwargs = {}
    # Add file-based data.
    file_id = self._headers['file_id'][rec]
    if file_id >= 0:
      f = open(self._files[file_id], 'rb')
    else:
      f = None
    for key, (addr_key, len_key, d_key) in self._decoder_data:
      if addr_key not in self._headers: continue
      # Special case: have dask array to read.
      if d_key in self._headers:
        d = self._headers[d_key][rec]
        if d is not None:
          kwargs[key] = np.asarray(d)
          continue
      address = self._headers[addr_key][rec]
      length = self._headers[len_key][rec]
      if address == -1 or length == -1: continue
      f.seek(address,0)
      data = self._decode(np.fromfile(f,'B',length))
      kwargs[key] = data
    for key in self._decoder_extra_args:
      if key in self._headers:
        value = self._headers[key][rec]
        kwargs[key] = value
    kwargs.update(self._decoder_scalar_args())
    if f is not None:
      f.close()
    return self._postproc(**kwargs)


  #
  ###############################################



