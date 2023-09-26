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
# Mixin for pre-filtering the RPN file records.

class FilterRecords (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    super(FilterRecords,cls)._cmdline_args(parser)
    parser.add_argument('--filter', metavar='CONDITION', action='append', help=_("Subset %s records using the given criteria.  For example, to convert only 24-hour forecasts you could use --filter ip2==24")%cls._format)
  def __init__ (self, *args, **kwargs):
    """
    filter : str or list, optional
        Subset RPN file records using the given criteria.  For example, to
        convert only 24-hour forecasts you could use filter="ip2==24"
    """
    import numpy as np
    import re
    filter = kwargs.pop('filter',None)
    if filter is None:
      filter = []
    if isinstance(filter,str):
      filter = [filter]
    self._filters = tuple(filter)
    super(FilterRecords,self).__init__(*args,**kwargs)
    if len(self._filters) == 0: return
    flags = True
    # Get copy of records with bytestrings and with unicode strings.
    records_with_bytestrings = self._headers
    records_with_strings = self._headers.copy()
    for key, value in records_with_strings.items():
      if value.dtype.char == 'S' and any(key in cmd for cmd in self._filters):
        new_dtype = 'U'+str(value.dtype.itemsize)
        value = value[:,None].view('B').astype('i4')
        # rstrip
        mask = np.cumsum(value[:,::-1]-32, axis=1)[:,::-1] == 0
        value[mask] = 0
        value = value.flatten()
        records_with_strings[key] = value.view(new_dtype)
    bytestrings_pattern = re.compile("b['\"]")
    # Loop over each filter and apply.
    for cmd in self._filters:
      try:
        # Check if the filter uses byte strings or standard (unicode) strings.
        if re.search(bytestrings_pattern, cmd):
          flags &= self._do_filter(records_with_bytestrings, cmd)
        else:
          flags &= self._do_filter(records_with_strings, cmd)
      except TypeError:
        error (_("unable to apply the filter: %s")%cmd)
    # To filter out unwanted records, unselected in the list.
    self._headers['selected'] = self._headers['selected'] & flags
  @staticmethod
  def _do_filter (p, cmd):
    # Allow numpy functions to work with the filters.
    import numpy
    p = dict(p, numpy=numpy, np=numpy)
    try:
      return eval(cmd, None, p)
    except SyntaxError:
      error (_("unable to parse the filter: %s")%cmd)
    except NameError as e:
      error (str(e))


