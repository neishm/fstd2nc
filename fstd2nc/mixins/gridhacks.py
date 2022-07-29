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
  def _decode (self, data, unused, _grid_cache={}):
    import rpnpy.librmn.all as rmn
    import numpy as np
    from fstd2nc.extra import decode_headers, decode
    if not hasattr(self,'_interp_grid'):
      return super(Interp,self)._decode (data, unused)
    grid_params = ('ni','nj','grtyp','ig1','ig2','ig3','ig4')
    prm = decode_headers(data[:72])
    prm = dict((k,v[0]) for k,v in prm.items())
    prm['d'] = super(Interp,self)._decode (data, unused)
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


#################################################
# Mixin for yin/yang grid subsetting.
#

class YinYang (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    from argparse import SUPPRESS
    super(YinYang,cls)._cmdline_args(parser)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--yin', action='store_true', help=_('Select first subgrid from a supergrid.'))
    group.add_argument('--yang', action='store_true', help=_('Select second subgrid from a supergrid.'))

  def __init__ (self, *args, **kwargs):
    """
    yin : bool, optional
        Select first subgrid from a supergrid.
    yang : bool, optional
        Select second subgrid from a supergrid.
    """
    import rpnpy.librmn.all as rmn
    import tempfile
    from os import path
    import numpy as np
    from fstd2nc.extra import structured_array
    from collections import OrderedDict
    self._yin = kwargs.pop('yin',False)
    self._yang = kwargs.pop('yang',False)

    # Load the metadata from the file(s).
    super(YinYang,self).__init__(*args,**kwargs)

    # Update yin-yang records to appear as regular rotated grids.
    if not self._yin and not self._yang:
      return

    # Find all unique YY input grids.
    mask = (self._headers['ismeta']==0) & (self._headers['grtyp']==b'U')
    grids = OrderedDict((key,self._headers[key][mask]) for key in ('ig1','ig2','ig3','ig4'))
    grids, rec_ids = np.unique(structured_array(grids), return_index=True)
    rec_ids = np.where(mask)[0][rec_ids]
    nomvars = self._headers['nomvar'][rec_ids]
    source_grids = dict()
    for nomvar, grid in zip(nomvars,grids):
      ig1 = grid['ig1']
      ig2 = grid['ig2']
      ig3 = grid['ig3']
      ig4 = grid['ig4']
      for handle in rmn.fstinl(self._meta_funit, nomvar=nomvar.decode()):
        prm = rmn.fstprm(handle)
        if prm['grtyp'] == b'U' and prm['ig1'] == ig1 and prm['ig2'] == ig2 and prm['ig3'] == ig3 and prm['ig4'] == ig4:
          break
      source_grids[(ig1,ig2,ig3,ig4)] = rmn.readGrid (self._meta_funit, prm)

    # Generate temporary file with target grid info.
    try: # Python 3
      self._yy_tmpdir = tempfile.TemporaryDirectory()
      gridfile = path.join(self._yy_tmpdir.name,"grid.fst")
    except AttributeError: # Python 2 (no auto cleanup)
      self._yy_tmpdir = tempfile.mkdtemp()
      gridfile = path.join(self._yy_tmpdir,"grid.fst")
    iun = rmn.fstopenall(gridfile, rmn.FST_RW)

    if self._yin: yy_ind = 0
    else: yy_ind = 1

    # Write all target grids and set the headers appropriately.
    for (ig1,ig2,ig3,ig4), source_grid in source_grids.items():
      dest_grid = source_grid['subgrid'][yy_ind]
      rmn.writeGrid (iun, dest_grid)
      submask = mask & (self._headers['ig1'] == ig1) & (self._headers['ig2'] == ig2) & (self._headers['ig3'] == ig3) & (self._headers['ig4'] == ig4)
      for key in ('grtyp','ni','nj','ig1','ig2','ig3','ig4'):
        self._headers[key][submask] = dest_grid[key]
    rmn.fstcloseall (iun)

    # Hack in these extra coordinate records, and refresh metadata file unit.
    rmn.fstcloseall (self._meta_funit)
    self._meta_filenames.append(gridfile)
    self._meta_funit = rmn.fstopenall(self._meta_filenames,rmn.FST_RO)

  # Handle grid interpolation
  def _fstluk (self, rec_id, dtype=None, rank=None, dataArray=None):
    import rpnpy.librmn.all as rmn
    import numpy as np
    out = super(YinYang,self)._fstluk (rec_id, dtype, rank, dataArray)
    # Check original grtyp to see if need to split.
    if out['grtyp'] == 'U':
      if self._yin:
        out['d'] = out['d'][:,:out['nj']//2]
      elif self._yang:
        out['d'] = out['d'][:,out['nj']//2:]
    return out

  # Handle grid interpolations from raw binary array.
  def _decode (self, data, unused):
    from fstd2nc.extra import decode_headers
    if not self._yin and not self._yang:
      return super(YinYang,self)._decode (data, unused)
    prm = decode_headers(data[:72])
    prm = dict((k,v[0]) for k,v in prm.items())
    d = super(YinYang,self)._decode (data, unused)
    if prm['grtyp'] == b'U' and self._yin:
      d = d[:,:prm['nj']//2]
    elif prm['grtyp'] == b'U' and self._yang:
      d = d[:,prm['nj']//2:]
    return d

#################################################
# Mixin for grid cropping.
#

class Crop (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    from argparse import SUPPRESS
    super(Crop,cls)._cmdline_args(parser)
    parser.add_argument('--crop-to-smallest-grid', action='store_true', help=_('Crop grids to the smaller (inner core) domain for LAM outputs.'))

  def __init__ (self, *args, **kwargs):
    """
    crop_to_smallest_room : bool, optional
        Crop grids to the smaller (inner core) domain for LAM outputs.
    """
    import rpnpy.librmn.all as rmn
    import numpy as np
    from fstd2nc.extra import structured_array
    from collections import OrderedDict
    self._crop_to_smallest_grid = kwargs.pop('crop_to_smallest_grid',False)

    # Load the metadata from the file(s).
    super(Crop,self).__init__(*args,**kwargs)

    if not self._crop_to_smallest_grid:
      return

    # Keep track of cropping regions for the data.
    self._headers['i0'] = np.zeros(self._nrecs,'uint16')
    self._headers['iN'] = np.array(self._headers['ni'],'uint16')
    self._headers['j0'] = np.zeros(self._nrecs,'uint16')
    self._headers['jN'] = np.array(self._headers['nj'],'uint16')

    # Find all unique grids.
    ###
    # Skip supergrids.
    mask = (self._headers['ismeta']==0) & (self._headers['grtyp'] != b'U')
    grids = OrderedDict((key,self._headers[key][mask]) for key in ('ig1','ig2','ig3','ig4'))
    grids, rec_ids = np.unique(structured_array(grids), return_index=True)
    rec_ids = np.where(mask)[0][rec_ids]
    source_grids = dict()
    for rec_id in rec_ids:
      prm = dict()
      prm['grtyp'] = self._headers['grtyp'][rec_id].decode()
      prm['ni'] = int(self._headers['ni'][rec_id])
      prm['nj'] = int(self._headers['nj'][rec_id])
      prm['ig1'] = int(self._headers['ig1'][rec_id])
      prm['ig2'] = int(self._headers['ig2'][rec_id])
      prm['ig3'] = int(self._headers['ig3'][rec_id])
      prm['ig4'] = int(self._headers['ig4'][rec_id])
      key = tuple(prm[k] for k in ('ig1','ig2','ig3','ig4'))
      source_grids[key] = rmn.readGrid (self._meta_funit, prm)

    ###

    # Group by ig1ref, etc.
    matches = dict()
    for grid in source_grids.values():
      key = (grid['ig1ref'],grid['ig2ref'],grid['ig3ref'],grid['ig4ref'])
      matches.setdefault(key,[]).append(grid)

    # Find the smallest grid for each projection.
    smallest = dict()
    for key, grids in matches.items():
      shapes = [(len(grid['ax'].flatten()),len(grid['ay'].flatten())) for grid in grids]
      ind = shapes.index(min(shapes))
      smallest[key] = ind

    # Ok, now look at all the source grids, and figure out which ones are
    # candidates for cropping.
    for (ig1,ig2,ig3,ig4), grid in source_grids.items():
      key = (grid['ig1ref'],grid['ig2ref'],grid['ig3ref'],grid['ig4ref'])
      smallest_grid = matches[key][smallest[key]]
      # Nothing to do if already the smallest grid.
      if grid['id'] == smallest_grid['id']: continue
      # Make sure the coordinates are actually compatible.
      ind = np.searchsorted(grid['ax'].flatten(), smallest_grid['ax'].flatten())
      i0, iN = ind[0], ind[-1]+1
      ind = np.searchsorted(grid['ay'].flatten(), smallest_grid['ay'].flatten())
      j0, jN = ind[0], ind[-1]+1
      if iN-i0 != smallest_grid['ni']: continue
      if jN-j0 != smallest_grid['nj']: continue
      if np.any(grid['ax'].flatten()[i0:iN] != smallest_grid['ax'].flatten()): continue
      if np.any(grid['ay'].flatten()[j0:jN] != smallest_grid['ay'].flatten()): continue
      # Able to crop, so update the headers to point to the cropped coordinates.
      submask = mask & (self._headers['ig1'] == ig1) & (self._headers['ig2'] == ig2) & (self._headers['ig3'] == ig3) & (self._headers['ig4'] == ig4)
      for key in ('grtyp','ni','nj','ig1','ig2','ig3','ig4'):
        self._headers[key][submask] = smallest_grid[key]
      self._headers['i0'][submask] = i0
      self._headers['iN'][submask] = iN
      self._headers['j0'][submask] = j0
      self._headers['jN'][submask] = jN

  # Handle cropping
  def _fstluk (self, rec_id, dtype=None, rank=None, dataArray=None):
    out = super(Crop,self)._fstluk (rec_id, dtype, rank, dataArray)
    if not self._crop_to_smallest_grid: return out
    if self._headers['ismeta'][rec_id] == 1: return out
    # Check if cropping necessary.
    ni = self._headers['ni'][rec_id]
    nj = self._headers['nj'][rec_id]
    if out['d'].shape == (ni,nj): return out
    i0 = self._headers['i0'][rec_id]
    iN = self._headers['iN'][rec_id]
    j0 = self._headers['j0'][rec_id]
    jN = self._headers['jN'][rec_id]
    out['d'] = out['d'][i0:iN,j0:jN]
    return out

  # Handle cropping from raw binary array.
  def _decode (self, data, rec_id):
    import numpy as np
    from fstd2nc.extra import decode_headers
    d = super(Crop,self)._decode (data, rec_id)
    if not self._crop_to_smallest_grid: return d
    # Check if cropping necessary.
    ni = self._headers['ni'][rec_id]
    nj = self._headers['nj'][rec_id]
    if d.shape == (ni,nj): return d
    i0 = self._headers['i0'][rec_id]
    iN = self._headers['iN'][rec_id]
    j0 = self._headers['j0'][rec_id]
    jN = self._headers['jN'][rec_id]
    d = d[i0:iN,j0:jN]
    return d

