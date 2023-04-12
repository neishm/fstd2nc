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

# Helper method - decode name information (from bytes) into string.
from fstd2nc.mixins import vectorize
@vectorize
def to_string (name):
  return name.decode().strip()


#################################################
# Mixin for selecting particular fields.

class SelectVars (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    super(SelectVars,cls)._cmdline_args(parser)
    parser.add_argument('--vars', metavar='VAR1,VAR2,...', help=_('Comma-separated list of variables to convert.  By default, all variables are converted.'))

  def __init__ (self, *args, **kwargs):
    """
    vars : str or list, optional
        Comma-separated list of variables to convert.  By default, all
        variables are converted.
    """
    import numpy as np
    vars = kwargs.pop('vars',None)
    super(SelectVars,self).__init__(*args,**kwargs)

    if vars is None:
      return

    if isinstance(vars,str):
      vars = vars.replace(',', ' ')
      vars = vars.split()
    info (_('Will look for variables: ') + ' '.join(vars))
    select = np.zeros(self._nrecs,dtype='bool')
    missing = []
    names = to_string(self._headers['name'])
    names = np.array(names,object)
    for v in vars:
      f = (names == v)
      if not np.any(f):
        missing.append(v)
      select |= f
    if len(missing) > 0:
      warn(_('Unable to find variable(s): ') + ' '.join(missing))
    if not np.any(select):
      error(_('Nothing to convert.'))
    # Marked unselected variables.
    self._headers['selected'] = self._headers['selected'] & select

