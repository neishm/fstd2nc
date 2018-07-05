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

#####################################################################
# Mixin for an ensemble axis.

class Ensembles (BufferBase):

  # Define any command-line arguments for reading FSTD files.
  @classmethod
  def _cmdline_args (cls, parser):
    super(Ensembles,cls)._cmdline_args(parser)
    parser.add_argument('--ensembles', action='store_true', help=_('Collect different etikets for the same variable together into an "ensemble" axis.'))

  def __init__ (self, *args, **kwargs):
    ensembles = kwargs.pop('ensembles',False)
    if ensembles:
      self._outer_axes = ('etiket',) + self._outer_axes
      # Don't split the variable on different etikets.
      kwargs['ignore_etiket'] = True
    super(Ensembles,self).__init__(*args,**kwargs)

  def _iter (self):
    from fstd2nc.mixins import _var_type, _modify_axes
    from collections import OrderedDict
    import numpy as np
    handled_etikets = dict()
    for var in super(Ensembles,self)._iter():
      if 'etiket' not in var.axes:
        yield var
        continue
      etikets = var.axes['etiket']
      coordinates = var.atts.get('coordinates',[])
      if etikets not in handled_etikets:
        array = np.array(var.axes['etiket'])
        # Strip out trailing whitespace.
        array[:] = map(str.rstrip,array)
        #array = array.view('|S1').reshape(n,strlen)
        # Encode it as 2D character array for netCDF file output.
        n = len(array)
        array = array.view('|S1').reshape(n,-1)
        n, strlen = array.shape
        axes = OrderedDict([('ensemble_id',etikets),('ensemble_strlen',tuple(range(strlen)))])
        etiket_var = _var_type('ensemble',{},axes,array)
        yield etiket_var
        handled_etikets[etikets] = etiket_var
      var.axes = _modify_axes(var.axes, etiket='ensemble_id')
      coordinates.append(handled_etikets[etikets])
      if len(coordinates) > 0:
        var.atts['coordinates'] = coordinates
      yield var
