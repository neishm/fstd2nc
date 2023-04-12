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
# Mixin for compute diagnostic fields on-the-fly for xarray interface.
#

class DiagHacks (BufferBase):
  # Override to_xarray() to allow extra fields to be generated.
  def to_xarray (self, **kwargs):
    add_px = kwargs.pop('add_px',False)
    add_gz = kwargs.pop('add_gz',False)
    dataset = super(DiagHacks,self).to_xarray(**kwargs)
    if add_px and 'PX' not in dataset.variables:
      dataset['PX'] = _compute (dataset, 'p')
      dataset['PX'].attrs.update(standard_name='air_pressure',units='hPa')
    if add_gz and 'GZ' not in dataset.variables:
      dataset['GZ'] = _compute (dataset, 'z')
      dataset['GZ'].attrs.update(standard_name='geopotential_height',units='m')
    return dataset

def _compute (x, field, samplefield=None):
  from xarray.ufuncs import exp, log

  # Find an appropriate set of levels.
  level = None
  for v in x.variables:
    if samplefield is not None and v != samplefield:
      continue
    for coord in x[v].coords:
      if coord.startswith('level'):
        l = x[coord]
        if hasattr(l,'formula') and hasattr(l,'formula_terms'):
          level = l
  if level is None:
    raise ValueError("No appropriate levels found.  You must include at least one 3D reference field in your dataset.")

  # Apply the pressure formula for these levels.
  terms = level.formula_terms.replace(': ', ':')
  terms = dict([t.split(':') for t in terms.split()])
  for v in terms.values():
    if v not in x:
      raise ValueError ("Missing field '%s', cannot calculate the diagnostic field."%v)
  terms = dict((k,x[v]) for k,v in terms.items())
  locals().update(terms)
  exec (level.formula)
  try:
    var = locals()[field]
  except KeyError:
     raise ValueError ("Cannot calculate the diagnostic field (wrong model coordinate type).")
  dims = var.dims
  # Make time the first dimension.
  if any (d.startswith('time') for d in dims):
    dims = [d for d in dims if d.startswith('time')] + \
           [d for d in dims if not d.startswith('time')]
    var = var.transpose(*dims)
  return var
