###############################################################################
# Copyright 2017-2024 - Climate Research Division
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
# Mixin for adding extended metadata for fst24 records.
#
class XMeta (BufferBase):

  @classmethod
  def _cmdline_args (cls, parser):
    import argparse
    super(XMeta,cls)._cmdline_args(parser)
    parser.add_argument('--extended-metadata', action='store_true', help=_('Add extended metadata to the variables (fst24 data only).'))

  def __init__ (self, *args, **kwargs):
    """
    extended_metadata : bool, optional
        Add extended metadata to the variables (fst24 data only).
    """
    self._extended_metadata = kwargs.pop('extended_metadata',False)
    super(XMeta,self).__init__(*args, **kwargs)

  def _makevars (self):
    from collections import OrderedDict
    import json

    super(XMeta,self)._makevars()
    if not self._extended_metadata: return
    if 'file_id' not in self._headers or 'address' not in self._headers or 'meta_length' not in self._headers:
      return

    # Add extended metadata to the variables.
    # Because of the way it's stored in the data segments (not headers), need
    # to do some extra file reading here to get it.
    last_filename = None
    last_f = None
    for var in self._varlist:

      # Use the first metadata instance for the whole variable.
      rec = var.record_id.flatten()[0]
      file_id = self._headers['file_id'][rec]
      if file_id < 0: continue
      address = self._headers['address'][rec]
      if address < 0: continue
      meta_length = self._headers['meta_length'][rec]
      if meta_length <= 18: continue
      filename = self._files[file_id]
      if filename == last_filename:
        f = last_f
      else:
        f = open(filename, 'rb')
      last_filename = filename
      last_f = f

      f.seek(address+72,0)
      meta = f.read(meta_length*4-72).decode().rstrip('\0')
      # This line would decode the JSON string, if it's needed in the future.
      # meta = json.loads(meta,object_pairs_hook=OrderedDict)
      var.atts['extended_metadata'] = meta
