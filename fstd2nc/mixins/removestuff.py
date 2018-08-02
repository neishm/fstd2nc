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
    parser.add_argument('--exclude', metavar='NAME,NAME,...', help=_("Exclude some axes or derived variables from the output.  Note that axes will only be excluded if they have a length of 1."))

  def __init__ (self, *args, **kwargs):
    import numpy as np
    exclude = kwargs.pop('exclude',None)
    if exclude is None:
      exclude = []
    if isinstance(exclude,str):
      exclude = exclude.replace(',', ' ')
      exclude = exclude.split()
    self._exclude = tuple(exclude)
    super(RemoveStuff,self).__init__(*args,**kwargs)

  def _iter (self):
    from fstd2nc.mixins import _var_type, _iter_type
    for var in super(RemoveStuff,self)._iter():

      # Omit the excluded variables from the data stream.
      # Note: will not remove the var associated with an axis if the axis is
      # non-degenerate.
      if var.name in self._exclude:
        if list(var.axes.keys()) != [var.name]:
          continue
        elif len(var.axes[var.name]) == 1:
          continue

      # Remove the excluded axes.
      for exclude in self._exclude:
        axisnames = list(var.axes.keys())
        if exclude in axisnames and len(var.axes[exclude]) == 1:
          i = axisnames.index(exclude)
          var.axes.pop(exclude)
          if isinstance(var,_var_type):
            shape = var.array.shape
            shape = shape[:i] + shape[i+1:]
            var.array = var.array.reshape(shape)
          if isinstance(var,_iter_type):
            shape = var.record_id.shape
            shape = shape[:i] + shape[i+1:]
            var.record_id = var.record_id.reshape(shape)

      # Remove all references to excluded vars from the metadata.
      for key,val in list(var.atts.items()):
        # Special case - list of object references
        if isinstance(val,list):
          val = [v for v in val if not hasattr(v,'name') or v.name not in self._exclude]
          var.atts[key] = val

      yield var
