###############################################################################
# Copyright 2017-2021 - Climate Research Division
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
# Mixin for on-the-fly grid interpolation.
#

class Interp (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    from argparse import SUPPRESS
    super(Interp,cls)._cmdline_args(parser)
    #parser.add_argument('--interp', metavar="GRTYP,PARAM=VAL,PARAM=VAL,...", help=_("Interpolate to the specified grid."))
    parser.add_argument('--interp', metavar="GRTYP,PARAM=VAL,PARAM=VAL,...", help=SUPPRESS)

  def __init__ (self, *args, **kwargs):
    """
    interp : str or list, optional
        Interpolate to the specified grid.
    """
    import rpnpy.librmn.all as rmn
    import tempfile
    from os import path
    import numpy as np
    interp = kwargs.pop('interp',None)
    # Extract interpolation grid.
    def number(x):
      try: return int(x)
      except ValueError: return float(x)
    if interp is not None:
      interp = interp.split(',')
      grtyp = interp[0]
      grid_args = [number(v) for v in interp[1:] if '=' not in v]
      grid_kwargs = dict(v.split('=') for v in interp[1:] if '=' in v)
      grid_kwargs = dict((k,number(v)) for k,v in grid_kwargs.items())
      if not hasattr(rmn,'defGrid_'+grtyp):
        error(_("Unknown grid '%s'")%grtyp)
      self._interp_grid = getattr(rmn,'defGrid_'+grtyp)(*grid_args,**grid_kwargs)
      # Generate temporary file with target grid info.
      try: # Python 3
        self._interp_tmpdir = tempfile.TemporaryDirectory()
        gridfile = path.join(self._interp_tmpdir.name,"grid.fst")
      except AttributeError: # Python 2 (no auto cleanup)
        self._interp_tmpdir = tempfile.mkdtemp()
        gridfile = path.join(self._interp_tmpdir,"grid.fst")
      iun = rmn.fstopenall(gridfile, rmn.FST_RW)
      rmn.writeGrid (iun, self._interp_grid)
      rmn.fstcloseall (iun)
      # Handle case where a single file was passed from the Python interface.
      if isinstance(args[0],str):
        args = list(args)
        args[0] = [args[0]]
        args = tuple(args)
      # Add this grid file to the inputs.
      args[0].append(gridfile)
      super(Interp,self).__init__(*args,**kwargs)
      # Update records to use this grid.
      ismeta = self._headers['ismeta']
      for key in ('grtyp','ni','nj','ig1','ig2','ig3','ig4'):
        self._headers[key][:] = np.where(ismeta, self._headers[key], self._interp_grid[key])
      # Set up some options for ezsint.
      rmn.ezsetopt (rmn.EZ_OPT_EXTRAP_DEGREE, rmn.EZ_EXTRAP_VALUE)
      rmn.ezsetopt (rmn.EZ_OPT_EXTRAP_VALUE, self._fill_value)
    else:
      super(Interp,self).__init__(*args,**kwargs)

  # Add fill value to the data.
  # Modified from code in from mask mixin.
  # Set this fill value for any interpolated data, since it may contain
  # points outside the original boundary.
  def _makevars (self):
    from fstd2nc.mixins import _iter_type
    super(Interp,self)._makevars()
    if hasattr(self,'_interp_grid'):
      for var in self._varlist:
        if isinstance(var,_iter_type):
          var.atts['_FillValue'] = var.dtype.type(self._fill_value)

  # Handle grid interpolation
  def _fstluk (self, rec_id, dtype=None, rank=None, dataArray=None, _grid_cache={}):
    import rpnpy.librmn.all as rmn
    import numpy as np
    grid_params = ('ni','nj','grtyp','ig1','ig2','ig3','ig4')
    if not hasattr(self,'_interp_grid') or self._headers['ismeta'][rec_id] == 1:
      return super(Interp,self)._fstluk (rec_id, dtype, rank, dataArray)
    with self._lock:
      prm = super(Interp,self)._fstluk(rec_id, dtype, rank, dataArray)
      cache_key = tuple(prm[k] for k in grid_params)
      # Check if we've already defined the input grid in a previous call.
      if cache_key in _grid_cache:
        ingrid = _grid_cache[cache_key]
      else:
        ingrid = dict((k,prm[k]) for k in grid_params)
        ingrid['iunit'] = self._meta_funit
        ingrid = rmn.ezqkdef (**ingrid)
        _grid_cache[cache_key] = ingrid
      # Propogate any fill values to the interpolated grid.
      in_mask = np.zeros(prm['d'].shape, order='F', dtype='float32')
      in_mask[prm['d']==self._fill_value] = 1.0
      d = rmn.ezsint (self._interp_grid, ingrid, prm['d'])
      out_mask = rmn.ezsint (self._interp_grid, ingrid, in_mask)
    d[out_mask!=0] = self._fill_value
    # Return the data and metadata for the interpolated field.
    prm = dict((k,self._interp_grid[k]) for k in prm.keys() if k in self._interp_grid)
    prm['d'] = d
    return prm

  # Handle grid interpolations from raw binary array.
  def _decode (self, data, _grid_cache={}):
    import rpnpy.librmn.all as rmn
    import numpy as np
    from fstd2nc.extra import decode_headers, decode
    if not hasattr(self,'_interp_grid'):
      return super(Interp,self)._decode (data)
    grid_params = ('ni','nj','grtyp','ig1','ig2','ig3','ig4')
    prm = decode_headers(data[:72])
    prm = dict((k,v[0]) for k,v in prm.items())
    prm['d'] = super(Interp,self)._decode (data)
    with self._lock:
      cache_key = tuple(prm[k] for k in grid_params)
      # Check if we've already defined the input grid in a previous call.
      if cache_key in _grid_cache:
        ingrid = _grid_cache[cache_key]
      else:
        ingrid = dict((k,(str(prm[k].decode()) if k=='grtyp' else int(prm[k]))) for k in grid_params)
        ingrid['iunit'] = self._meta_funit
        ingrid = rmn.ezqkdef (**ingrid)
        _grid_cache[cache_key] = ingrid

    # Propogate any fill values to the interpolated grid.
    in_mask = np.zeros(prm['d'].shape, order='F', dtype='float32')
    in_mask[prm['d']==self._fill_value] = 1.0
    d = rmn.ezsint (self._interp_grid, ingrid, prm['d'])
    out_mask = rmn.ezsint (self._interp_grid, ingrid, in_mask)
    d[out_mask!=0] = self._fill_value
    # Return the data for the interpolated field.
    return d

