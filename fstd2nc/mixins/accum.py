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
# Mixin for handling accumulated values.
# E.g. precipitation.

accum_nomvars = ('PR','PRMM','PB','PC','PY','PZ')

class Accum (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    super(Accum,cls)._cmdline_args(parser)
    parser.add_argument('--accum-vars', metavar=_('NAME,NAME,...'), help=_('Specify which fields to treat as accumulated quantities (using IP3 as accumulation period).'))

  def __init__ (self, *args, **kwargs):
    """
    accum_vars : str or list, optional
        Specify which fields to treat as accumulated quantities
        (using IP3 as accumulation period).
    """
    import numpy as np
    accum_vars = kwargs.pop('accum_vars',None)
    if accum_vars is None:
      accum_vars = []
    if isinstance(accum_vars,str):
      accum_vars = accum_vars.replace(',', ' ')
      accum_vars = accum_vars.split()
    self._accum_nomvars = accum_nomvars + tuple(accum_vars)

    self._outer_axes = ('accum',) + self._outer_axes

    super(Accum,self).__init__(*args,**kwargs)

    # Compute an accumulation period.
    accum = np.ma.asarray (self._headers['ip3'], dtype='float32')
    accum.mask = (np.isin(self._headers['grtyp'],(b'T',b'+'))) | (self._headers['ismeta']) | (self._headers['ip3'] >= 1000) | (self._headers['ip1'] > 0)
    self._headers['accum'] = accum

  def _makevars (self):
    from fstd2nc.mixins import _axis_type
    import numpy as np
    super(Accum,self)._makevars()
    # Remove degenerate accumumlation periods.
    for var in self._varlist:
      if 'accum' not in var.dims: continue
      ind = var.dims.index('accum')
      # Check if 'accum' axis should be stripped out.
      if var.shape[ind] == 1:
        ip3 = int(var.axes[ind].array[0])
        if var.name not in self._accum_nomvars and ip3 != 0:
          warn(_('Ignoring IP3 value %d for %s, of unknown purpose.  If this is an accumulated variable, you could specify it in the --accum-vars option.')%(ip3,var.name))
        # For single accumulation periods, don't add an 'accum' axis.
        # Instead, put the accumulation period in an attribute.
        # This is for backwards-compatibility with previous versions of this
        # package where this axis did not exist.
        if var.name in self._accum_nomvars:
          var.atts['accumulation_period'] = '%d hours'%ip3
        var.axes.pop(ind)
        if ind < var.record_id.ndim:
          var.record_id = var.record_id.squeeze(axis=ind)
      # Check if we have an 'accum' axis that we don't know what to do with.
      elif var.shape[ind] > 1 and var.name not in self._accum_nomvars:
        warn(_('Multiple IP3 values of unknown purpose found for %s.  If this is an accumulated variable, you could specify it in the --accum-vars option.')%var.name)
        axis = var.axes[ind]
        var.axes[ind] = _axis_type('ip3_values', {}, axis.array)
      # Annotate the 'accum' axis.
      else:
        var.axes[ind].atts.update(long_name='accumulation_period',units='hours')

  # Re-encode accumulation time back into ip3.
  def _unmakevars (self):
    from fstd2nc.mixins import _axis_type
    import numpy as np
    for var in self._varlist:
      # Catch generic ip3 axis, put it into 'accum' column in table for now.
      if 'ip3_values' in var.dims:
        ind = var.dims.index('ip3_values')
        axis = var.axes[ind]
        var.axes[ind] = _axis_type('accum', {}, axis.array)
      # Check for single accumulation period stored as an attribute.
      if 'accumulation_period' in var.atts:
        axis = _axis_type('accum', {}, np.array([int(var.atts['accumulation_period'].split()[0])]))
        var.axes = [axis] + var.axes
        var.record_id = var.record_id.reshape((1,)+var.record_id.shape)
    super (Accum,self)._unmakevars()
    if 'accum' not in self._headers.keys(): return
    # Put accumulation times into ip3.
    if 'ip3' not in self._headers.keys():
      self._headers['ip3'] = np.ma.maskedall(self._nrecs,dtype='uint32')
    # Convert accum to units of hours if it's a timedelta64.
    if self._headers['accum'].dtype == 'timedelta64[ns]':
      self._headers['accum'] = self._headers['accum'] // np.timedelta64(3600,'s')
    valid = ~np.ma.getmaskarray(self._headers['accum'])
    self._headers['ip3'][valid] = self._headers['accum'][valid]

