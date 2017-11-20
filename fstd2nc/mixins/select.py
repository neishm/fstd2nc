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

from fstd2nc import _, info, warn, error
from fstd2nc.mixins import _Buffer_Base


#################################################
# Mixin for selecting particular fields.

class _SelectVars (_Buffer_Base):
  @classmethod
  def _cmdline_args (cls, parser):
    super(_SelectVars,cls)._cmdline_args(parser)
    parser.add_argument('--vars', metavar='VAR1,VAR2,...', help=_('Comma-separated list of variables to convert.  By default, all variables are converted.'))
  def __init__ (self, *args, **kwargs):
    vars = kwargs.pop('vars',None)
    if vars is not None:
      self._selected_vars = vars.split(',')
      info (_('Will look for variables: ') + ' '.join(self._selected_vars))
    else:
      self._selected_vars = None
    super(_SelectVars,self).__init__(*args,**kwargs)
  def _iter (self):
    found = set()
    for var in super(_SelectVars,self)._iter():
      if self._selected_vars is not None:
        if var.name not in self._selected_vars:
          continue
      found.add(var.name)
      yield var
    if self._selected_vars is None: return
    missing = set(self._selected_vars) - found
    if len(missing) > 0:
      warn(_('Unable to find variables: ') + ' '.join(missing))

