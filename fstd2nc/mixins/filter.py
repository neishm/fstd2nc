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
# Mixin for pre-filtering the RPN file records.

class _FilterRecords (_Buffer_Base):
  @classmethod
  def _cmdline_args (cls, parser):
    super(_FilterRecords,cls)._cmdline_args(parser)
    parser.add_argument('--filter', metavar='CONDITION', action='append', help=_("Subset RPN file records using the given criteria.  For example, to convert only 24-hour forecasts you could use --filter ip2==24"))
  def __init__ (self, *args, **kwargs):
    import numpy as np
    filter = kwargs.pop('filter',None)
    if filter is None:
      filter = []
    self._filters = tuple(filter)
    super(_FilterRecords,self).__init__(*args,**kwargs)
    if len(self._filters) == 0: return
    flags = np.ones(len(self._headers),dtype='bool')
    records = dict([(n,self._headers[n]) for n in self._headers.dtype.names])
    for cmd in self._filters:
      try:
        flags &= self._do_filter(records, cmd)
      except TypeError:
        error (_("unable to apply the filter: %s")%cmd)
    nrecs = sum(flags)
    self._headers = np.ma.empty(nrecs, dtype=self._headers_dtype)
    for k,v in list(records.items()):
        self._headers[k] = v[flags]
  @staticmethod
  def _do_filter (p, cmd):
    try:
      return eval(cmd, None, p)
    except SyntaxError:
      error (_("unable to parse the filter: %s")%cmd)
    except NameError as e:
      error (e.message)


