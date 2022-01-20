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
    self._fill_value = kwargs.pop('fill_value',1e30)
    super(Masks,self).__init__(*args,**kwargs)
    # Remove all mask records from the table, they should not become variables
    # themselves.
    is_mask = (self._headers['typvar'] == b'@@')
    self._headers['dltf'][is_mask] = 1
    # Modify 'lng' to overlap the data and mask fields (where they are stored
    # consecutively).
    # This will allow the _decode method to apply the mask directly.
    nomvar = self._headers['nomvar']
    typvar = self._headers['typvar']
    etiket = self._headers['etiket']
    datev = self._headers['datev']
    ip1 = self._headers['ip1']
    ip2 = self._headers['ip2']
    ip3 = self._headers['ip3']
    swa = self._headers['swa']
    lng = self._headers['lng'].copy()
    uses_mask = np.array(typvar,dtype='|S2').view('|S1').reshape(-1,2)[:,1] == b'@'
    if np.sum(uses_mask) > 0:
      overlap = (nomvar[:-1] == nomvar[1:]) & (etiket[:-1] == etiket[1:]) & (datev[:-1] == datev[1:]) & (ip1[:-1] == ip1[1:]) & (ip2[:-1] == ip2[1:]) & (ip3[:-1] == ip3[1:]) & uses_mask[:-1] & uses_mask[1:]
      overlap = np.where(overlap)[0]
      lng[overlap] = swa[overlap+1] + lng[overlap+1] - swa[overlap]
      self._headers['lng'] = lng

  # Apply the fill value to the data.
  def _makevars (self):
    from fstd2nc.mixins import _iter_type
    super(Masks,self)._makevars()
    for var in self._varlist:
      if not isinstance(var,_iter_type):
        continue
      # Look for typvars such as 'P@'.
      if var.atts.get('typvar','').endswith('@'):
        var.atts['_FillValue'] = var.dtype.type(self._fill_value)

  # Apply the mask data
  def _fstluk (self, rec_id, dtype=None, rank=None, dataArray=None):
    import numpy as np
    from rpnpy.librmn.fstd98 import fstinf
    with self._lock:
      prm = super(Masks,self)._fstluk(rec_id, dtype, rank, dataArray)
      # If this data is not masked, or if this data *is* as mask, then just
      # return it.
      if not prm['typvar'].endswith('@'): return prm
      if prm['typvar'] == '@@' : return prm
      mask_key = fstinf(self._opened_funit, nomvar=prm['nomvar'], typvar = '@@',
                        datev=prm['datev'], etiket=prm['etiket'],
                        ip1 = prm['ip1'], ip2 = prm['ip2'], ip3 = prm['ip3'])
    if mask_key is not None:
      mask = super(Masks,self)._fstluk(mask_key, rank=rank)['d']
      prm['d'] *= mask
      prm['d'] += self._fill_value * (1-mask)
    return prm

  # Apply the mask data from raw binary array.
  def _decode (self, data):
    d = data.view('>i4')
    # Get first field.
    code1 = (d[12]&0x000FFF00)>>8
    field1 = super(Masks,self)._decode(data)
    # Find out where the next field should be.
    swa = d[1]
    offset = d[0]&(0x00FFFFFF)
    while offset < len(d)//2 and swa+offset != d[offset*2+1]:
      offset += 1
    else:
      # No extra fields found?
      if offset >= len(d)//2:
        return field1
    data = data.view('>i4')[offset*2:].view(data.dtype)
    d = data.view('>i4')
    # Get second field.
    code2 = (d[12]&0x000FFF00)>>8
    field2 = super(Masks,self)._decode(data)
    # Apply the mask.
    if code1 == 2080:
      return (field1*field2) + self._fill_value * (1-field1)
    elif code2 == 2080:
      return (field1*field2) + self._fill_value * (1-field2)
    else:
      # Don't know how to apply this typvar for masking purposes.
      return field1
