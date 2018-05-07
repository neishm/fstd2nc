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
from fstd2nc.mixins import BufferBase


#################################################
# Mixin for removing variables and degenerate axes from the stream.

class RemoveStuff (BufferBase):

  @classmethod
  def _cmdline_args (cls, parser):
    super(RemoveStuff,cls)._cmdline_args(parser)
    parser.add_argument('--omit', metavar='NAME1,NAME2,...', help=_("Omit some axes or derived variables from the output.  Note that axes can only be omitted if they have a length of 1."))

  def __init__ (self, *args, **kwargs):
    import numpy as np
    omit = kwargs.pop('omit',None)
    if omit is None:
      omit = []
    if isinstance(omit,str):
      omit = omit.split(',')
    self._omit = tuple(omit)
    super(RemoveStuff,self).__init__(*args,**kwargs)

  def _iter (self):
    from fstd2nc.mixins import _var_type, _iter_type
    for var in super(RemoveStuff,self)._iter():
      if var.name in self._omit:
        continue

      # Remove omitted axes.
      for omit in self._omit:
        axisnames = list(var.axes.keys())
        if omit in axisnames and len(var.axes[omit]) == 1:
          i = axisnames.index(omit)
          var.axes.pop(omit)
          if isinstance(var,_var_type):
            shape = var.array.shape
            shape = shape[:i] + shape[i+1:]
            var.array = var.array.reshape(shape)
          if isinstance(var,_iter_type):
            shape = var.record_id.shape
            shape = shape[:i] + shape[i+1:]
            var.record_id = var.record_id.reshape(shape)

      # Remove all references to omitted vars from the metadata.
      for key,val in list(var.atts.items()):
        if isinstance(val,str):
          val = val.split()
          for omit in self._omit:
            if omit in val:
              val[val.index(omit)] = None
          val = filter(None,val)
          var.atts[key] = ' '.join(val)
        # Special case - list of object references
        elif isinstance(val,list):
          val = [v for v in val if v.name not in self._omit]
          var.atts[key] = val

      yield var
