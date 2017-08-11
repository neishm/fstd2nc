#!/usr/bin/env python

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


"""
Serve RPN standard files through a pydap server.
"""

# Handler for the FST data.
from pydap.handlers.lib import BaseHandler
class FST_Handler(BaseHandler):
  from fstd2nc import __version__
  extensions = r"^.*$"
  def __init__ (self, filepath):
    from pydap.model import DatasetType, GridType, BaseType
    from fstd2nc import Buffer
    from os.path import basename
    import numpy as np
    from collections import OrderedDict
    BaseHandler.__init__(self)
    self.filepath = filepath
    buf = Buffer(filepath)
    self.dataset = DatasetType(name=basename(filepath), attributes=dict(NC_GLOBAL=buf._metadata.get('global',{})))

    # First, scan into a dictionary
    vars = OrderedDict((var.name,var) for var in buf)

    # Based loosely on pydap's builtin netcdf handler.
    for var in vars.values():
      # Add grids.
      self.dataset[var.name] = GridType(var.name, var.atts)
      # Add array.
      self.dataset[var.name][var.name] = BaseType(var.name, var.array, tuple(var.axes.keys()), var.atts)
      # Add maps.
      for dim in var.axes.keys():
        if dim in vars:
          atts = vars[dim].atts
          array = vars[dim].array
        else:
          atts = {}
          array = np.array(var.axes[dim])
        self.dataset[var.name][dim] = BaseType(dim, array, None, atts)

# Hack this handler into the pydap interface.
# Can't register this in the normal way, because we can't provide specific
# extensions for the handler.
from pydap.handlers.lib import get_handler as pydap_get_handler
def get_handler(filepath, handlers=None):
  from rpnpy.librmn.fstd98 import isFST
  if isFST(str(filepath)):
    return FST_Handler(str(filepath))
  return pydap_get_handler(filepath, handlers)
from pydap.handlers import lib as pydap_lib
pydap_lib.get_handler = get_handler


# After hacking in this handler, invoke the pydap server.
if __name__ == '__main__':
  from pydap.wsgi.app import main
  main()
