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
      return result.values
    else:
      # If we're given a scalar value, then simply return it.
      if not hasattr(x,'__len__'):
        return cached_f(x)
      return list(map(cached_f,x))
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
  __slots__ = ('name','atts','axes','dtype','chunks','chunksize','deps')
  def __init__ (self, name, atts, axes, dtype, chunks, chunksize):
    self.name = name
    self.atts = atts
    self.axes = axes
    self.dtype = dtype
    self.chunks = chunks
    self.chunksize = chunksize
    self.deps = []
  def items(self):
    for key, value in self.chunks.items():
      key = tuple(k if isinstance(k,int) else slice(*k) for k in key)
      yield key, value
  def keys(self):
    for key, value in self.items():
      yield key
  def values(self):
    return self.chunks.values()

# Fake progress bar - does nothing.
class _FakeBar (object):
  def __init__ (self, *args, **kwargs): pass
  def iter(self, it):
    for i in it: yield i
  def __next__(self): pass
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


# Define a lock for controlling threaded access to the RPN file(s).
# This is necessary because we can only open a limited number of files
# at a time, so need to control which ones are opened and available in
# a thread-safe way.
from threading import RLock
_lock = RLock()
del RLock

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
  # Other records that should also be kept separate, but only if they are
  # 1-dimensional.  If they are 2D, then they should be processed as a normal
  # variable.
  _maybe_meta_records = ()

  # Attributes which could potentially be used as axes.
  _outer_axes = ()

  # Attributes which could be used as auxiliary coordinates for the outer
  # axes.  The dictionary keys are the outer axis names, and the values are
  # a list of columns which can act as coordinates.
  from collections import OrderedDict
  _outer_coords = OrderedDict()
  del OrderedDict

  # Attributes which uniquely identify a variable.
  # Note: nomvar should always be the first attribute
  _var_id = ('nomvar','ni','nj','nk')

  # Similar to above, but a human-readable version of the id.
  # Could be used for appending suffixes to variable names to make them unique.
  # Uses string formatting operations on the variable metadata.
  _human_var_id = ('%(nomvar)s', '%(ni)sx%(nj)s', '%(nk)sL')

  # Record parameters which should not be used as nc variable attributes.
  # (They're either internal to the file, or part of the data, not metadata).
  _ignore_atts = ('file_id','swa','lng','dltf','ubc','xtra1','xtra2','xtra3','key','shape','d')

  # Define any command-line arguments for reading FSTD files.
  @classmethod
  def _cmdline_args (cls, parser):
    from fstd2nc import __version__
    from argparse import SUPPRESS
    parser.add_argument('--version', action='version', version=__version__)
    group = parser.add_mutually_exclusive_group()
    _('Display a progress bar during the conversion, if the "progress" module is installed.')
    group.add_argument('--progress', action='store_true', default=True, help=SUPPRESS)
    group.add_argument('--no-progress', action='store_false', dest='progress', help=_('Disable the progress bar.'))
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--minimal-metadata', action='store_true', default=True, help=_("Don't include RPN record attributes and other internal information in the output metadata.")+" "+_("This is the default behaviour."))
    group.add_argument('--rpnstd-metadata', action='store_false', dest='minimal_metadata', help=_("Include all RPN record attributes in the output metadata."))
    group.add_argument('--rpnstd-metadata-list', metavar='nomvar,...', help=_("Specify a minimal set of RPN record attributes to include in the output file."))
    parser.add_argument('--ignore-typvar', action='store_true', help=_('Tells the converter to ignore the typvar when deciding if two records are part of the same field.  Default is to split the variable on different typvars.'))
    parser.add_argument('--ignore-etiket', action='store_true', help=_('Tells the converter to ignore the etiket when deciding if two records are part of the same field.  Default is to split the variable on different etikets.'))
    parser.add_argument('--serial', action='store_true', help=_('Disables multithreading/multiprocessing.  Useful for resource-limited machines.'))

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
      # Python 3: convert bytes to str
      if isinstance(v,bytes):
        v = str(v.decode())
      if isinstance(v,str):
        v = v.strip()
      yield (n,v)

  # Open the specified file (by index)
  def _open (self, file_id):
    from rpnpy.librmn.base import fnom
    from rpnpy.librmn.fstd98 import fstouv
    from rpnpy.librmn.const import FST_RO
    from rpnpy.librmn import librmn
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

  def __init__ (self, filename, _headers=None, progress=False, minimal_metadata=None, rpnstd_metadata=None, rpnstd_metadata_list=None, ignore_typvar=False, ignore_etiket=False, serial=False):
    """
    Read raw records from FSTD files, into the buffer.
    Multiple files can be read simultaneously.

    Parameters
    ----------
    filename : str or list
        The RPN standard file(s) to convert.
    progress : bool, optional
        Display a progress bar during the conversion, if the "progress"
        module is installed.
    rpnstd_metadata : bool, optional
        Include all RPN record attributes in the output metadata.
    rpnstd_metadata_list : str or list, optional
        Specify a minimal set of RPN record attributes to include in the
        output file.
    ignore_typvar : bool, optional
        Tells the converter to ignore the typvar when deciding if two
        records are part of the same field.  Default is to split the
        variable on different typvars.
    ignore_etiket : bool, optional
        Tells the converter to ignore the etiket when deciding if two
        records are part of the same field.  Default is to split the
        variable on different etikets.
    serial : bool, optional
        Disables multithreading/multiprocessing.  Useful for resource-limited
        machines.
    """
    from rpnpy.librmn.fstd98 import fstnbr, fstinl, fstprm, fstopenall
    from rpnpy.librmn.const import FST_RO
    from fstd2nc.extra import raw_headers, decode_headers
    from collections import Counter
    import numpy as np
    from glob import glob, has_magic
    import os
    import warnings
    from multiprocessing import Pool
    try:
      from itertools import imap  # Python 2
    except ImportError:
      imap = map

    # Set up lock for threading.
    # The same lock is shared for all Buffer objects, to synchronize access to
    # librmn.
    self._lock = _lock

    # Set up a progress bar for scanning the input files.
    Bar = _ProgressBar if progress is True else _FakeBar

    # Set default for minimal_metadata
    if rpnstd_metadata is not None:
      minimal_metadata = not rpnstd_metadata
    if minimal_metadata is None:
      minimal_metadata = True
    # Set default for rpnstd_metadata_list
    if minimal_metadata is True and rpnstd_metadata_list is None:
      rpnstd_metadata_list = ''
    if isinstance(rpnstd_metadata_list,str):
      rpnstd_metadata_list = rpnstd_metadata_list.replace(',',' ')
      rpnstd_metadata_list = rpnstd_metadata_list.split()
    if hasattr(rpnstd_metadata_list,'__len__'):
      rpnstd_metadata_list = tuple(rpnstd_metadata_list)
    self._rpnstd_metadata_list = rpnstd_metadata_list

    if not ignore_typvar:
      # Insert typvar value just after nomvar.
      self._var_id = self._var_id[0:1] + ('typvar',) + self._var_id[1:]
      self._human_var_id = self._human_var_id[0:1] + ('%(typvar)s',) + self._human_var_id[1:]
    if not ignore_etiket:
      # Insert etiket value just after nomvar.
      self._var_id = self._var_id[0:1] + ('etiket',) + self._var_id[1:]
      self._human_var_id = self._human_var_id[0:1] + ('%(etiket)s',) + self._human_var_id[1:]

    self._serial = serial

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

    # Extract headers from the files.
    if len(expanded_infiles) == 1: Bar = _FakeBar
    bar = Bar(_("Inspecting input files"), suffix='%(percent)d%% (%(index)d/%(max)d)', max=len(expanded_infiles))

    if len(expanded_infiles) > 1 and not self._serial:
      with Pool() as p:
        headers = p.imap (raw_headers, [f for (infile,f) in expanded_infiles])
        headers = bar.iter(headers)
        headers = list(headers) # Start scanning.
    else:
      headers = imap (raw_headers, [f for (infile,f) in expanded_infiles])
      headers = bar.iter(headers)
      headers = list(headers) # Start scanning.

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
        file_ids.extend([filenum]*(len(headers[i])//72))
      else:
        matches[infile] += 0

    # Decode all the headers
    headers = [h for h in headers if h is not None]
    # Remember indices in the headers (needed for reconstructing keys)
    indices = [range(len(h)//72) for h in headers]
    if len(headers) > 0:
      headers = np.concatenate(headers)
      indices = np.concatenate(indices)
      headers = decode_headers(headers)
      # Encode the keys without the file index info.
      headers['key'] = (indices % 256) | ((indices//256)<<9)
      headers['file_id'] = np.array(file_ids, dtype='int32')

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

    nfiles = len(self._files)
    if nfiles == 0:
      error(_("no input files found!"))
    elif nfiles > 10:
      global _pandas_needed
      _pandas_needed = True
    info(_("Found %d RPN input file(s)"%nfiles))

    # Add extra headers? (hack for generating Buffer objects from external
    # sources (not FSTD files)
    if _headers is not None:
      headers = dict(_headers)
      # Keep link to first file for metadata purposes below?
      headers['file_id'] = np.zeros(len(headers['nomvar']), dtype='int32')

    self._headers = headers
    self._nrecs = len(headers['nomvar'])

    # Find all unique meta (coordinate) records, and link a subset of files
    # that provide all unique metadata records.
    # This will make it easier to look up the meta records later.
    meta_mask = np.zeros(self._nrecs,dtype='bool')
    for meta_name in self._meta_records:
      meta_name = (meta_name+b'   ')[:4]
      meta_mask |= (self._headers['nomvar'] == meta_name)
    for meta_name in self._maybe_meta_records:
      meta_name = (meta_name+b'   ')[:4]
      meta_mask |= (self._headers['nomvar'] == meta_name) & ((self._headers['ni']==1)|(self._headers['nj']==1))

    # Store this metadata identification in case it's useful for some mixins.
    self._headers['ismeta'] = np.empty(self._nrecs, dtype='bool')
    self._headers['ismeta'][:] = meta_mask

    meta_recids = np.where(meta_mask)[0]
    # Use the same unique parameters as regular variables.
    # Plus, ig1,ig2,ig3,ig4.
    # Suppress FutureWarning from numpy about doing this.  Probably benign...
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      from fstd2nc.extra import structured_array
      headers = structured_array(self._headers)
      meta_keys = headers.data[meta_mask][list(self._var_id)+['ip1','ip2','ip3','ig1','ig2','ig3','ig4']]
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


  # Generate structured variables from FST records.

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
    valid = (records['dltf'] == 0) & (records['ismeta'] == 0)
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
      nomvar = var_id['nomvar'].strip()
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
        # Use regular integers for numeric types.
        if np.can_cast(v.dtype,int):
          v = int(v)
        # Python3: convert bytes to str.
        if isinstance(v,bytes): v = str(v.decode())
        # Trim string attributes (remove whitespace padding).
        if isinstance(v,str): v = v.strip()
        atts[n] = v

      # Get the axes.
      axes = OrderedDict()
      for n in self._outer_axes:
        values = var_records[n]
        # Remove missing values before continuing.
        values = np.ma.compressed(values)
        # Ignore axes that have no actual coordinate values.
        if len(values) == 0: continue
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
      if len(coords) > 0:
        atts['coordinates'] = coords

      # Check if we have full coverage along all axes.
      have_data = [k >= 0 for k in record_id.flatten()]
      if not np.all(have_data):
        warn (_("Missing some records for %s.")%nomvar)

      # Add dummy axes for the ni,nj,nk record dimensions.
      axes['k'] = _dim_type(name='k', length = int(var_id['nk']))
      axes['j'] = _dim_type(name='j', length = int(var_id['nj']))
      axes['i'] = _dim_type(name='i', length = int(var_id['ni']))

      # Determine the optimal data type to use.
      # First, find unique combos of datyp, nbits
      # (want to minimize calls to dtype_fst2numpy).
      datyp, nbits = zip(*np.unique(var_records.data[['datyp','nbits']]))
      datyp = map(int,datyp)
      nbits = map(int,nbits)
      dtype_list = map(dtype_fst2numpy, datyp, nbits)
      dtype = np.result_type(*dtype_list)

      var = _iter_type( name = nomvar, atts = atts,
                        axes = list(axes.values()),
                        dtype = dtype,
                        record_id = record_id )
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

    # Ignore deleted / invalidated records.
    records = records[records['dltf']==0]

    # Ignore coordinate records.
    records = records[records['ismeta']==0]

    # Keep track of any axes that were generated.
    known_axes = dict()

    # Keep track of any auxiliary coordinates that were generated.
    known_coords = dict()

    # Iterate over each variable.
    # Variables are defined by the entries in _var_id.
    self._varlist = []
    for var_id, var_records in records.groupby(list(self._var_id)):
      var_id = OrderedDict(zip(self._var_id, var_id))
      nomvar = var_id['nomvar'].strip()
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
        if var_records[n].isnull().values.any(): continue
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
        cat = pd.Categorical(column)
        # Is this column an outer axis?
        if n in self._outer_axes:
          values = tuple(cat.categories)
          if (n,values) not in known_axes:
            known_axes[(n,values)] = _axis_type(name = n, atts = OrderedDict(),
                                   array = np.array(values,dtype=column.dtype))
          axes[n] = known_axes[(n,values)]
          indices[n] = cat.codes
          # Is this also an axis for an auxiliary coordinate?
          for coordname,coordaxes in self._outer_coords.items():
            if n in coordaxes:
              coordnames.append(coordname)
              coord_axes.setdefault(coordname,OrderedDict())[n] = axes[n]
              coord_indices.setdefault(coordname,OrderedDict())[n] = cat.codes
        # Otherwise, does it have a consistent value?
        # If so, can add it to the metadata.
        # Ignore outer coords, since the value is already encoded elsewhere.
        elif len(cat.categories) == 1 and n not in self._outer_coords:
          try:
            v = cat[0]
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
        if var_records[n].isnull().values.any(): continue
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
      if len(coords) > 0:
        atts['coordinates'] = coords



      # Check if we have full coverage along all axes.
      have_data = [k >= 0 for k in record_id.flatten()]
      if not np.all(have_data):
        warn (_("Missing some records for %s.")%nomvar)

      # Add dummy axes for the ni,nj,nk record dimensions.
      axes['k'] = _dim_type('k',int(var_id['nk']))
      axes['j'] = _dim_type('j',int(var_id['nj']))
      axes['i'] = _dim_type('i',int(var_id['ni']))

      # Determine the optimal data type to use.
      # First, find unique combos of datyp, nbits
      # (want to minimize calls to dtype_fst2numpy).
      x = var_records[['datyp','nbits']].drop_duplicates()
      datyp = map(int,x['datyp'])
      nbits = map(int,x['nbits'])
      dtype_list = map(dtype_fst2numpy, datyp, nbits)
      dtype = np.result_type(*dtype_list)

      var = _iter_type( name = nomvar, atts = atts,
                        axes = list(axes.values()),
                        dtype = dtype,
                        record_id = record_id )
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

  # How to decode the data from a raw binary array.
  def _decode (self, data):
    from fstd2nc.extra import decode
    nbits = int(data[0x0b])
    datyp = int(data[0x13])
    dtype = dtype_fst2numpy(datyp, nbits)
    out = decode(data).view(dtype)
    return out

  #
  ###############################################



