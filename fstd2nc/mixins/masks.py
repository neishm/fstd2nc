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

from fstd2nc.stdout import _, info, warn, error
from fstd2nc.mixins import Buffer_Base


#################################################
# Mixin for handling masks.

class Masks (Buffer_Base):
  @classmethod
  def _cmdline_args (cls, parser):
    super(Masks,cls)._cmdline_args(parser)
    parser.add_argument('--fill-value', type=float, default=1e30, help=_("The fill value to use for masked (missing) data.  Gets stored as '_FillValue' attribute in the metadata.  Default is '%(default)s'."))
  def __init__ (self, *args, **kwargs):
    self._fill_value = kwargs.pop('fill_value',1e30)
    super(Masks,self).__init__(*args,**kwargs)
    # Remove all mask records from the table, they should not become variables
    # themselves.
    is_mask = (self._headers['typvar'] == '@@')
    self._headers['dltf'][is_mask] = 1

  # Apply the fill value to the data.
  def _iter (self):
    from fstd2nc.mixins import _iter_type
    for var in super(Masks,self)._iter():
      if not isinstance(var,_iter_type):
        yield var
        continue
      # Look for typvars such as 'P@'.
      if var.atts.get('typvar','').endswith('@'):
        var.atts['_FillValue'] = var.dtype.type(self._fill_value)
      yield var

  # Apply the mask data
  def _fstluk (self, rec_id, dtype=None, rank=None, dataArray=None):
    import numpy as np
    from rpnpy.librmn.fstd98 import fstinf
    prm = super(Masks,self)._fstluk(rec_id, dtype, rank, dataArray)
    # If this data is not masked, or if this data *is* as mask, then just
    # return it.
    if not prm['typvar'].endswith('@'): return prm
    if prm['typvar'] == '@@' : return prm
    mask_key = fstinf(self._opened_funit, nomvar=prm['nomvar'], typvar = '@@',
                      datev=prm['datev'], etiket=prm['etiket'],
                      ip1 = prm['ip1'], ip2 = prm['ip2'], ip3 = prm['ip3'])
    if mask_key is not None:
      mask = self._fstluk(mask_key, rank=rank)['d']
      prm['d'] *= mask
      prm['d'] += self._fill_value * (1-mask)
    return prm

