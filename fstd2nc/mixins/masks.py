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

# Disable automatic decoding of datyp+64 missing value codes.
# Easier to find the encoded missing values (maximum value in the field).
# TODO: Remove this once the missing value interface is wrapped in rpnpy.
try:
  import ctypes
  from rpnpy.librmn import librmn
  librmn.ForceMissingValueUsage(ctypes.c_int(0))
  del ctypes, librmn
except ImportError:  # rpnpy not available, so nothing to do here.
  pass


#################################################
# Mixin for handling masks.

class Masks (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    super(Masks,cls)._cmdline_args(parser)
    parser.add_argument('--fill-value', type=float, default=1e30, help=_("The fill value to use for masked (missing) data.  Gets stored as '_FillValue' attribute in the metadata.  Default is '%(default)s'."))

  def __init__ (self, *args, **kwargs):
    """
    fill_value : scalar, optional
        The fill value to use for masked (missing) data.  Gets stored as
        '_FillValue' attribute in the metadata.  Default is 1e30.
    """
    import numpy as np
    from fstd2nc.extra import structured_array
    from collections import OrderedDict

    self._fill_value = kwargs.pop('fill_value',1e30)
    self._ignore_atts = ('mask_offset','mask_rec') + self._ignore_atts
    super(Masks,self).__init__(*args,**kwargs)

    nomvar = self._headers['nomvar']
    typvar = self._headers['typvar']
    etiket = self._headers['etiket']
    datev = self._headers['datev']
    ip1 = self._headers['ip1']
    ip2 = self._headers['ip2']
    ip3 = self._headers['ip3']
    dltf = self._headers['dltf']
    uses_mask = np.array(typvar,dtype='|S2').view('|S1').reshape(-1,2)[:,1] == b'@'
    if not np.any(uses_mask): return

    # Remove all mask records from the table, they should not become variables
    # themselves.
    is_mask = (self._headers['typvar'] == b'@@')
    self._headers['selected'][is_mask] = False

    nrecs = len(self._headers['name'])

    # Easy case - masks appear immediately following the data.
    # Find all such pairings in the headers.
    has_mask = (nomvar[:-1] == nomvar[1:]) & (etiket[:-1] == etiket[1:]) & (datev[:-1] == datev[1:]) & (ip1[:-1] == ip1[1:]) & (ip2[:-1] == ip2[1:]) & (ip3[:-1] == ip3[1:]) & uses_mask[:-1] & uses_mask[1:] & (dltf[:-1] == dltf[1:])
    has_mask = np.where(has_mask)[0]
    # If data is from an FSTD file on disk, then try overlapping the read of
    # the data and mask.
    self._headers['length'][has_mask] = self._headers['address'][has_mask+1] + self._headers['length'][has_mask+1] - self._headers['address'][has_mask]
    # Add a column in the headers to keep track of where to find the mask data
    # from the combined fields.
    #TODO: check how this interacts with externally generated records.
    # (from_fstpy)
    try:
      self._headers['mask_offset'] = np.empty(nrecs,dtype='int32')
      self._headers['mask_offset'][:] = -1
      if len(has_mask) > 0:
        self._headers['mask_offset'][has_mask] = self._headers['address'][has_mask+1] - self._headers['address'][has_mask]
    except KeyError:
      pass  # Not our data from disk (maybe from fstpy?) so can't do that.

    # Harder case - masks are in the file, but in random order.
    unordered_mask = np.zeros(nrecs,dtype='bool')
    unordered_mask[(typvar==b'@@')&(dltf==0)] = True
    if len(has_mask) > 0:
      unordered_mask[has_mask+1] = False
    if not np.any(unordered_mask): return

    # Figure out how to pair up the data and mask.
    # Requires O(n log n) time, which is better than O(n^2) for naive lookup
    # on each record.
    keys = OrderedDict([('dltf',dltf),('nomvar',nomvar),('etiket',etiket),('datev',datev),('ip1',ip1),('ip2',ip2),('ip3',ip3),('typvar',typvar)])
    ind = np.argsort(structured_array(keys))
    # Find mask / data pairs.
    # (In this case, with the above sort, masks will appear before the data)
    # Use similar comparisons as with the easy case.
    nomvar = nomvar[ind]
    typvar = typvar[ind]
    etiket = etiket[ind]
    datev = datev[ind]
    ip1 = ip1[ind]
    ip2 = ip2[ind]
    ip3 = ip3[ind]
    dltf = dltf[ind]
    uses_mask = np.array(typvar,dtype='|S2').view('|S1').reshape(-1,2)[:,1] == b'@'
    has_mask = (nomvar[:-1] == nomvar[1:]) & (etiket[:-1] == etiket[1:]) & (datev[:-1] == datev[1:]) & (ip1[:-1] == ip1[1:]) & (ip2[:-1] == ip2[1:]) & (ip3[:-1] == ip3[1:]) & uses_mask[:-1] & uses_mask[1:] & (dltf[:-1] == dltf[1:])
    has_mask = np.where(has_mask)[0]
    has_mask += 1  # Data appears after mask in this case.
    rec_id = np.arange(nrecs)[ind]
    self._headers['mask_rec'] = np.empty(nrecs,'int32')
    self._headers['mask_rec'][:] = -1
    self._headers['mask_rec'][rec_id[has_mask]] = rec_id[has_mask-1]
    del nomvar, typvar, etiket, datev, ip1, ip2, ip3, dltf, uses_mask, has_mask

  # Apply the fill value to the data.
  def _makevars (self):
    from fstd2nc.mixins import _iter_type
    super(Masks,self)._makevars()
    for var in self._varlist:
      if not isinstance(var,_iter_type):
        continue
      # Look for typvars such as 'P@'.
      if var.atts.get('typvar','').endswith('@') or var.atts.get('datyp',0) & 64 == 64:
        try:
          var.atts['_FillValue'] = var.dtype.type(self._fill_value)
        except OverflowError:
          warn(_("Can't set fill value '%g' for %s.")%(self._fill_value,var.name))

  # Apply the mask data from raw binary array.
  def _decode (self, data, rec_id):
    import numpy as np
    prm = self._decode_headers(data[:72])
    prm = dict((k,v[0]) for k,v in prm.items())
    # Get first field.
    field1 = super(Masks,self)._decode(data, rec_id)
    # If this data is encoded by the datyp+64 flag, then mask out the
    # largest value, assuming the mask was generated by (max-min)*1.01
    # or something similar.
    # TODO: Use the proper interface once it's wrapped in rpnpy.
    if (prm['datyp'] & 64) == 64:
      mx = field1.max()
      field1 = np.where(field1==mx, self._fill_value, field1)
      return field1
    # If typvar doesn't end in '@', then nothing to do here.
    if prm['typvar'][1:] != b'@':
      return field1
    if 'mask_offset' in self._headers:
      mask_offset = self._headers['mask_offset'][rec_id]
      if mask_offset >= 0:
        field2 = super(Masks,self)._decode(data[mask_offset:], rec_id+1)
        return ((field1*(field2>0)) + self._fill_value * (field2==0)).astype(field1.dtype)
    if 'mask_rec' in self._headers:
      mask_rec = self._headers['mask_rec'][rec_id]
      if mask_rec >= 0:
        field2 = self._read_record(mask_rec)
        return ((field1*(field2>0)) + self._fill_value * (field2==0)).astype(field1.dtype)
    return field1
