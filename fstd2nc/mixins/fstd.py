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
  import numpy as np
  datyp = int(datyp_nbits//(1<<32))
  nbits = int(datyp_nbits%(1<<32))
  dtype = dtype_fst2numpy (datyp, nbits)
  dtype = np.dtype(dtype)
  return (dtype.kind + str(dtype.itemsize)).encode()
def fast_dtype_fst2numpy (datyp, nbits):
  import numpy as np
  args = np.array(datyp,'uint64')
  args <<= 32
  args += np.asarray(nbits,'uint64')
  return packed_dtype_fst2numpy(args)

# Define a class for encoding / decoding FSTD data.
class FSTD (BufferBase):
  _format = _("RPN standard file")
  _format_singular = _("an RPN standard file")
  _format_plural = _("RPN standard file(s)")

  _inner_axes = ('j','i')

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

  # Helper method - find all records with the given criteria.
  # Mimics fstinl, but returns table indices instead of record handles.
  def _fstinl (self, **criteria):
    import numpy as np
    mask = True
    for k, v in criteria.items():
      mask &= (self._headers[k]==v)
    return np.where(mask)[0]

  # Helper method - get metadata of the given record.
  def _fstprm (self, rec):
    prm = dict((k,v[rec]) for k,v in self._headers.items())
    for k in ('typvar','nomvar','grtyp','etiket'):
      prm[k] = prm[k].decode()
    return prm

  # Helper method - read the specified record.
  # Mimics the behaviour of fstluk.
  def _fstluk (self, rec):
    prm = self._fstprm(rec)
    ni = prm['ni']
    nj = prm['nj']
    prm['d'] = self._read_record(rec).T.reshape(ni,nj)
    return prm

  # Helper method - read a record with the specified criteria.
  # Mimics the behaviour of fstlir.
  # Note: not the fastest method.  Should be used sparingly.
  def _fstlir (self, **criteria):
    recs = self._fstinl(**criteria)
    if len(recs) == 0: return None
    rec = recs[0]
    return self._fstluk(rec)



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

    # Note: name should always be the first attribute
    self._var_id = ('name','ni','nj') + self._var_id
    self._human_var_id = ('%(name)s', '%(ni)sx%(nj)s') + self._human_var_id
    self._ignore_atts = ('swa','lng','dltf','ubc','xtra1','xtra2','xtra3','i','j','k','nk','ismeta','d') + self._ignore_atts

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
    for meta_name in self._meta_records():
      meta_name = (meta_name+b'   ')[:4]
      meta_mask |= (self._headers['nomvar'] == meta_name)
    for meta_name in self._maybe_meta_records():
      meta_name = (meta_name+b'   ')[:4]
      meta_mask |= (self._headers['nomvar'] == meta_name) & ((self._headers['ni']==1)|(self._headers['nj']==1))

    # Store this metadata identification in case it's useful for some mixins.
    self._headers['ismeta'] = np.empty(self._nrecs, dtype='bool')
    self._headers['ismeta'][:] = meta_mask

    # Aliases for inner dimensions
    self._headers['j'] = self._headers['nj']
    self._headers['i'] = self._headers['ni']

    # Add some standard fields needed for the Buffer.
    self._headers['name'] = self._headers['nomvar']
    # These two fields may not exist for externally-sourced data
    # (such as from fstpy)
    if 'swa' in self._headers:
      self._headers['address'] = np.array(self._headers['swa'],int)*8-8
      del self._headers['swa']
    if 'lng' in self._headers:
      self._headers['length'] = np.array(self._headers['lng'],'int32')*4
      del self._headers['lng']
    self._headers['dtype'] = np.empty(self._nrecs, dtype='|S2')
    self._headers['dtype'][:] = np.array(fast_dtype_fst2numpy(self._headers['datyp'],self._headers['nbits']))
    self._headers['selected'] = (self._headers['dltf']==0) & (self._headers['ismeta'] == False)

  # How to decode the data from a raw binary array.
  @staticmethod
  def _decode (data):
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

  # Reconstructing FSTD records from external data.
  def _unmakevars (self):
    import numpy as np
    # Generate a table of records (with incomplete information).
    super(FSTD,self)._unmakevars()
    # Aliases for inner dimensions
    #self._headers['nk'] = self._headers['k']
    self._headers['nj'] = self._headers['j']
    self._headers['ni'] = self._headers['i']
    # Add other FSTD-related columns that are expected to be there.
    self._headers['nomvar'] = np.empty(self._nrecs,dtype='|S4')
    self._headers['nomvar'][:] = self._headers['name']
    self._headers['ismeta'] = np.zeros(self._nrecs,'bool')
    self._headers['ismeta'] |= np.isin(self._headers['nomvar'],self._meta_records())
    self._headers['ismeta'] |= np.isin(self._headers['nomvar'],self._maybe_meta_records()) & ((self._headers['ni']==1)|(self._headers['nj']==1))
    self._headers['selected'] = ~self._headers['ismeta']
    self._headers['address'] = np.empty(self._nrecs,dtype=int)
    self._headers['address'][:] = -1
    # Some columns may not have a specific value to put in the table, so put
    # in some placeholder value.
    def add_column (colname, dtype, default):
      if colname in self._headers.keys():
        self._headers[colname] = self._headers[colname].astype(dtype)
      else:
        self._headers[colname] = np.empty(self._nrecs,dtype=dtype)
        self._headers[colname][:] = default
      if hasattr(self._headers[colname],'mask'):
        self._headers[colname] = self._headers[colname].filled(default)

    add_column ('typvar', '|S2', default='P')
    add_column ('etiket', '|S12', default='')
    add_column ('datyp', 'int32', default=133)
    add_column ('nbits', 'int32', default=32)
    add_column ('ip3', 'int32', default=0)
    add_column ('deet', 'int32', default=60)
    add_column ('npas', 'int32', default=0)
    add_column ('nk', 'int32', default=1)

  # Define an entry point for writing records to a file.
  # Allows the logic to be modified by mixins (such as masks).
  @classmethod
  def _fstecr (cls, outfile, rec, **extra):
    import rpnpy.librmn.all as rmn
    # Make sure we pass in single-precision when expected.
    if rec.get('nbits',32) <= 32 and rec['d'].dtype.name.endswith('64'):
      rec['d'] = rec['d'].astype(rec['d'].dtype.name[:-2]+'32')
    # Disable compression on small records, to silence warning from
    # armn_compress32.
    if rec['ni'] <= 16 or rec['nj'] <= 16:
      rec['datyp'] %= 128
    rmn.fstecr(outfile, rec, **extra)

  def to_fstd (self, filename, append=False, rewrite=False):
    """
    Write data to an FSTD file.
    """
    from os.path import exists
    from os import remove
    import rpnpy.librmn.all as rmn
    import numpy as np
    if exists(filename) and not append and not rewrite: remove(filename)
    outfile = rmn.fstopenall(filename, rmn.FST_RW)
    for i in np.where(self._headers['selected'] | self._headers['ismeta'])[0]:
      rec = self._fstluk(i)
      # Ensure data is Fortran-contiguous for librmn.
      rec['d'] = np.ascontiguousarray(rec['d'].T).T
      self._fstecr(outfile, rec, rewrite=rewrite)
    rmn.fstcloseall(outfile)

  #
  ###############################################



