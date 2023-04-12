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

#####################################################################
# Mixin for an ensemble axis.

class Ensembles (BufferBase):

  # Define any command-line arguments for reading FSTD files.
  @classmethod
  def _cmdline_args (cls, parser):
    super(Ensembles,cls)._cmdline_args(parser)
    parser.add_argument('--ensembles', action='store_true', help=_('Collect different etikets for the same variable together into an "ensemble" axis.'))

  def __init__ (self, *args, **kwargs):
    """
    ensembles : bool, optional
        Collect different etikets for the same variable together into an
        "ensemble" axis.
    """
    ensembles = kwargs.pop('ensembles',False)
    if ensembles:
      self._outer_axes = ('etiket',) + self._outer_axes
      # Don't split the variable on different etikets.
      kwargs['ignore_etiket'] = True
    super(Ensembles,self).__init__(*args,**kwargs)

  def _makevars (self):
    from fstd2nc.mixins import _var_type, _dim_type
    import numpy as np
    super(Ensembles,self)._makevars()
    for etikets, varlist in self._iter_axes(name='etiket',varlist=True):
      # Python3: convert bytes to str.
      array = [str(arr.decode()) for arr in etikets.array]
      array = np.array(array,dtype=np.char.string_)
      # Strip out trailing whitespace.
      array[:] = [arr.rstrip() for arr in array]
      # Encode it as 2D character array for netCDF file output.
      n = len(array)
      array = array.view('|S1').reshape(n,-1)
      n, strlen = array.shape
      ensemble_id = _dim_type('ensemble_id', n)
      ensemble_strlen = _dim_type('ensemble_strlen',strlen)
      etiket_var = _var_type (name='ensemble', atts={}, axes=[ensemble_id, ensemble_strlen], array=array)
      for var in varlist:
        # Replace 'etiket' dimension of variable with ensemble_id dimension.
        var.axes[var.dims.index('etiket')] = ensemble_id
        # Add the etiket strings as a coordinate variable.
        coordinates = var.atts.get('coordinates',[])
        coordinates.append(etiket_var)
        var.atts['coordinates'] = coordinates
