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


#################################################
# Mixin for handling FSTD file data.
# Contains low-level details of extracting information for that format.

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
# Vectorized version, using datyp,bits packed into 64-bit integer.
from fstd2nc.mixins import vectorize
@vectorize
def packed_dtype_fst2numpy (datyp_nbits):
  datyp = int(datyp_nbits//(1<<32))
  nbits = int(datyp_nbits%(1<<32))
  return dtype_fst2numpy (datyp, nbits)
def fast_dtype_fst2numpy (datyp, nbits):
  import numpy as np
  args = np.array(datyp,'uint64')
  args <<= 32
  args += np.asarray(nbits,'uint64')
  return packed_dtype_fst2numpy(args)

# Define a lock for controlling threaded access to the RPN file(s).
# This is necessary because we can only open a limited number of files
# at a time, so need to control which ones are opened and available in
# a thread-safe way.
from threading import RLock
_lock = RLock()
del RLock


# Define a class for encoding / decoding FSTD data.
class FSTD (BufferBase):
  _format = _("RPN standard file")
  _format_singular = _("An RPN standard file")
  _format_plural = _("RPN standard file(s)")

  # Keep a reference to fstd98 so it's available during cleanup.
  try:
    from rpnpy.librmn import fstd98 as _fstd98
  except ImportError: pass

  # Define any command-line arguments for reading FSTD files.
  @classmethod
  def _cmdline_args (cls, parser):
    super(FSTD,cls)._cmdline_args (parser)
    parser.add_argument('--ignore-typvar', action='store_true', help=_('Tells the converter to ignore the typvar when deciding if two records are part of the same field.  Default is to split the variable on different typvars.'))
    parser.add_argument('--ignore-etiket', action='store_true', help=_('Tells the converter to ignore the etiket when deciding if two records are part of the same field.  Default is to split the variable on different etikets.'))

  # Clean up a buffer (close any attached files, etc.)
  def __del__ (self):
    try:
      from rpnpy.librmn.fstd98 import fstcloseall
      fstcloseall(self._meta_funit)
    except (ImportError, AttributeError):
      pass  # May fail if Python is doing a final cleanup of everything.
            # Or, if buffer wasn't fully initialized yet.

  # Open metadata funit (done once during init of object).
  def _open_meta_funit (self):
    from rpnpy.librmn.fstd98 import fstopenall
    from rpnpy.librmn.const import FST_RO
    import numpy as np
    import warnings
    meta_mask = self._headers['ismeta'][:]
    meta_recids = np.where(meta_mask)[0]
    _var_id = ('nomvar','ni','nj','nk')
    # Use the same unique parameters as regular variables.
    # Plus, ig1,ig2,ig3,ig4.
    # Suppress FutureWarning from numpy about doing this.  Probably benign...
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      from fstd2nc.extra import structured_array
      headers = dict()
      for key in list(_var_id)+['ip1','ip2','ip3','ig1','ig2','ig3','ig4']:
        headers[key] = self._headers[key]
      headers = structured_array(headers)
      meta_keys = headers.data[meta_mask][list(_var_id)+['ip1','ip2','ip3','ig1','ig2','ig3','ig4']]
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
    self._meta_filenames = filenames  # Store this for further hacking.
    self._meta_funit = fstopenall(filenames, FST_RO)
        
  # Control the pickling / unpickling of BufferBase objects.
  def __getstate__ (self):
    state = super(FSTD,self).__getstate__()
    state.pop('_lock',None)  # librmn lock will be defined upon unpickling.
    state.pop('_meta_funit',None) # Same with librmn file handles.
    return state
  def __setstate__ (self, state):
    super(FSTD,self).__setstate__(state)
    self._lock = _lock
    self._open_meta_funit()


  ###############################################
  # Basic flow for reading data

  def __init__ (self, *args, **kwargs):
    """
    ignore_typvar : bool, optional
        Tells the converter to ignore the typvar when deciding if two
        records are part of the same field.  Default is to split the
        variable on different typvars.
    ignore_etiket : bool, optional
        Tells the converter to ignore the etiket when deciding if two
        records are part of the same field.  Default is to split the
        variable on different etikets.
    """
    import numpy as np

    # Set up lock for threading.
    # The same lock is shared for all Buffer objects, to synchronize access to
    # librmn.
    self._lock = _lock

    self._inner_axes = ('k','j','i')

    # Note: name should always be the first attribute
    self._var_id = ('name','ni','nj','nk') + self._var_id
    self._human_var_id = ('%(name)s', '%(ni)sx%(nj)s', '%(nk)sL') + self._human_var_id
    self._ignore_atts = ('swa','lng','dltf','ubc','xtra1','xtra2','xtra3','key','shape','d','i','j','k','ismeta') + self._ignore_atts

    ignore_typvar = kwargs.pop('ignore_typvar',False)
    ignore_etiket = kwargs.pop('ignore_etiket',False)

    if not ignore_typvar:
      # Insert typvar value just after nomvar.
      self._var_id = self._var_id[0:1] + ('typvar',) + self._var_id[1:]
      self._human_var_id = self._human_var_id[0:1] + ('%(typvar)s',) + self._human_var_id[1:]
    if not ignore_etiket:
      # Insert etiket value just after nomvar.
      self._var_id = self._var_id[0:1] + ('etiket',) + self._var_id[1:]
      self._human_var_id = self._human_var_id[0:1] + ('%(etiket)s',) + self._human_var_id[1:]

    super(FSTD,self).__init__(*args,**kwargs)

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

    # Make metadata records available for FSTD-related mixins.
    self._open_meta_funit()

    # Aliases for iner dimensions
    self._headers['k'] = self._headers['nk']
    self._headers['j'] = self._headers['nj']
    self._headers['i'] = self._headers['ni']

    # Add some standard fields needed for the Buffer.
    self._headers['name'] = self._headers['nomvar']
    # These two fields may not exist for externally-sourced data
    # (such as from fstpy)
    if 'swa' in self._headers:
      self._headers['address'] = np.array(self._headers['swa'],int)*8-8
    if 'lng' in self._headers:
      self._headers['length'] = np.array(self._headers['lng'],int)*4
    self._headers['dtype'] = np.array(fast_dtype_fst2numpy(self._headers['datyp'],self._headers['nbits']))
    self._headers['selected'] = (self._headers['dltf']==0) & (self._headers['ismeta'] == False)

  # How to decode the data from a raw binary array.
  def _decode (self, data, unused):
    from fstd2nc.extra import decode
    nbits = int(data[0x0b])
    datyp = int(data[0x13])
    dtype = dtype_fst2numpy(datyp, nbits)
    out = decode(data).view(dtype)
    return out

  # Shortcuts to header decoding functions.
  # Put into the class so they can potentially be overridden for other formats.
  @staticmethod
  def _decode_headers (headers):
    from fstd2nc.extra import decode_headers
    return decode_headers(headers)
  @staticmethod
  def _raw_headers (filename):
    from fstd2nc.extra import raw_headers
    return raw_headers(filename)

  #
  ###############################################



