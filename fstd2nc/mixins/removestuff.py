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
# Mixin for removing variables and degenerate axes from the stream.

class RemoveStuff (BufferBase):

  @classmethod
  def _cmdline_args (cls, parser):
    super(RemoveStuff,cls)._cmdline_args(parser)
    parser.add_argument('--exclude', metavar=_('NAME,NAME,...'), help=_("Exclude some axes, attributes,  or derived variables from the output.  For instance, excluding 'leadtime,reftime' can help for netCDF tools that don't recognize leadtime and reftime as valid coordinates.  Note that axes will only be excluded if they have a length of 1."))

  def __init__ (self, *args, **kwargs):
    """
    exclude : str or list, optional
        Exclude some axes, attributes, or derived variables from the output.
        For instance, excluding 'leadtime,reftime' can help for netCDF tools
        that don't recognize leadtime and reftime as valid coordinates.
        Note that axes will only be excluded if they have a length of 1.
    """
    import numpy as np
    exclude = kwargs.pop('exclude',None)
    if exclude is None:
      exclude = []
    if isinstance(exclude,str):
      exclude = exclude.replace(',', ' ')
      exclude = exclude.split()
    self._exclude = tuple(exclude)
    super(RemoveStuff,self).__init__(*args,**kwargs)

  def _do_exclude (self, var):
    from collections import OrderedDict
    from fstd2nc.mixins import _iter_type, _var_type

    # Remove excluded axes.
    for axisname in self._exclude:
      axis = var.getaxis(axisname)
      if axis is not None and len(axis) == 1:
        i = var.dims.index(axisname)
        var.axes.pop(i)
        if isinstance(var,_var_type):
          shape = var.array.shape
          shape = shape[:i] + shape[i+1:]
          var.array = var.array.reshape(shape)
        if isinstance(var,_iter_type):
          shape = var.record_id.shape
          shape = shape[:i] + shape[i+1:]
          var.record_id = var.record_id.reshape(shape)
    # Remove other dependencies.
    var.deps = [d for d in var.deps if d.name not in self._exclude]

    # Remove all references to excluded vars from the metadata.
    for key,val in list(var.atts.items()):
      # Handle list of meta variables.
      if isinstance(val,list):
        # Remove any meta variables that are to be excluded.
        newval = [v for v in val if v.name not in self._exclude]
        # Recursively process any meta variables.
        newval = list(map(self._do_exclude,newval))
        var.atts[key] = newval
      # Handle key/value entries for meta variables.
      # (E.g. for formula_terms).
      elif isinstance(val,OrderedDict):
        newval = []
        for (k,v) in newval:
          # Remove any meta variables that are to be excluded.
          if isinstance(k,_var_type) and k.name in self._exclude:
            continue
          if isinstance(v,_var_type) and v.name in self._exclude:
            continue
          # Recursively process any meta variables.
          if isinstance(k,_var_type):
            k = self._do_exclude(k)
          if isinstance(v,_var_type):
            v = self._do_exclude(v)
          newval.append((k,v))
        var.atts[key] = OrderedDict(newval)
      # Remove excluded attributes.
      if key in self._exclude:
        var.atts.pop(key)
    # Remove from silent dependencies.
    var.deps = [self._do_exclude(v) for v in var.deps]

    return var

  def _makevars (self):
    from fstd2nc.mixins import _var_type, _iter_type

    super(RemoveStuff,self)._makevars()

    # Filter out excluded variables.
    self._varlist = [v for v in self._varlist if v.name not in self._exclude]

    # Apply exclusions to the axes and metadata.
    self._varlist = list(map(self._do_exclude,self._varlist))

