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

############################################################
# Mixin for questionable modifications to the headers table.
#

class GridHacks (BufferBase):
  # Switch on hacks.
  def _enable_hacks (self):
    if not hasattr(self,'_hacks'):
      self._hacks = dict()
      self._read_record = self.__read_record
      self._fstluk = self.__fstluk

  # Hack in record data that doesn't actually exist in the source file.
  def __read_record (self, rec_id):
    if rec_id not in getattr(self,'_hacks',{}):
      return super(GridHacks,self)._read_record(rec_id)
    return self._hacks[rec_id]['d'].T

  def __fstluk (self, rec_id):
    if rec_id not in getattr(self,'_hacks',{}):
      return super(GridHacks,self)._fstluk(rec_id)
    return self._hacks[rec_id]

  # Helper method - "write" a grid descriptor into the headers table.
  # Based loosely on writeGrid from rpnpy.
  def _writeGrid (self, gid):
    import numpy as np
    if gid['grtyp'] not in ('Z','#','Y','U'): return None
    self._enable_hacks()
    # Define some parameters for the grid.
    prm = dict()
    prm['typvar'] = 'X '
    prm['ip1'] = gid['tag1']
    prm['ip2'] = gid['tag2']
    prm['ip3'] = gid.get('tag3',0)
    prm['grtyp'] = gid['grref']
    prm['ig1'] = gid['ig1ref']
    prm['ig2'] = gid['ig2ref']
    prm['ig3'] = gid['ig3ref']
    prm['ig4'] = gid['ig4ref']
    prm['datyp'] = gid.get('datyp',5)
    prm['nbits'] = gid.get('nbits',32)
    prm['etiket'] = gid.get('etiket','ETIKET').ljust(12)
    prm['ni'] = 1; prm['nj'] = 1; prm['nk'] = 1
    prm['ismeta'] = True
    nrecs = len(self._headers['name'])
    new_recs = []
    if 'ax' in gid:
      new_recs.append(dict(prm, nomvar='>>  ', ni=gid['ni'], d=gid['ax']))
    if 'ay' in gid:
      new_recs.append(dict(prm, nomvar='^^  ', nj=gid['nj'], d=gid['ay']))
    if 'axy' in gid:
      new_recs.append(dict(prm, nomvar='^>  ', ni=gid['ni'], nj=gid['nj'], d=gid['axy']))
    for k,v in self._headers.items():
      self._headers[k] = np.zeros_like(v, shape=nrecs + len(new_recs))
      self._headers[k][:nrecs] = v[:]
    for i in range(len(new_recs)):
      for k,v in new_recs[i].items():
        if k in self._headers:
          if isinstance(v,str): v = v.encode()
          self._headers[k][nrecs+i] = v
      self._hacks[nrecs+i] = new_recs[i]
    self._nrecs += len(new_recs)


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
    import numpy as np
    interp = kwargs.pop('interp',None)
    super(Interp,self).__init__(*args,**kwargs)
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
      # Hack the new grid descriptors into the headers.
      self._writeGrid (self._interp_grid)

      # Store original and modified versions of the grid descriptors.
      ismeta = self._headers['ismeta']
      self._original_grid = dict()
      self._modified_grid = dict()
      for key in ('grtyp','ni','nj','ig1','ig2','ig3','ig4'):
        self._original_grid[key] = self._headers[key]
        self._modified_grid[key] = np.array(self._headers[key])
        self._modified_grid[key][:] = np.where(ismeta, self._headers[key], self._interp_grid[key])

      # Set up some options for ezsint.
      rmn.ezsetopt (rmn.EZ_OPT_EXTRAP_DEGREE, rmn.EZ_EXTRAP_VALUE)
      rmn.ezsetopt (rmn.EZ_OPT_EXTRAP_VALUE, self._fill_value)


  def _makevars (self):
    from fstd2nc.mixins import _iter_type
    import numpy as np

    if not hasattr(self,'_interp_grid'):
      return super(Interp,self)._makevars()

    # Run _makevars chain with original grid descriptors to allow
    # xycoords to give us source grid ids,
    # switch out the grid descriptors in the table, then let xycoords
    # construct the target grid axes for us.
    super(Interp,self)._makevars()
    self._source_gids = np.array(self._gids)
    self._headers.update(self._modified_grid)
    super(Interp,self)._makevars()

    # Add fill value to the data.
    # Modified from code in from mask mixin.
    # Set this fill value for any interpolated data, since it may contain
    # points outside the original boundary.
    for var in self._varlist:
      if isinstance(var,_iter_type):
        var.atts['_FillValue'] = var.dtype.type(self._fill_value)

  # Handle grid interpolations from raw binary array.
  def _decode (self, data, rec_id, _grid_cache={}):
    import rpnpy.librmn.all as rmn
    import numpy as np
    if self._headers['ismeta'][rec_id] or not hasattr(self,'_interp_grid'):
      return super(Interp,self)._decode (data, rec_id)
    # Retrieve an active librmn grid id associated with this grid.
    # (must be supplied by xycoords mixin).
    ingrid = int(self._source_gids[rec_id])
    if ingrid < 0:
      raise ValueError("Source data is not on a recognized grid.  Unable to interpolate.")
    d = super(Interp,self)._decode (data, rec_id).T
    with self._lock:
      # Propogate any fill values to the interpolated grid.
      in_mask = np.zeros(d.shape, order='F', dtype='float32')
      in_mask[d==self._fill_value] = 1.0
      d = rmn.ezsint (self._interp_grid, ingrid, d)
      out_mask = rmn.ezsint (self._interp_grid, ingrid, in_mask)
      d[out_mask!=0] = self._fill_value
      # Return the data for the interpolated field.
      return d.T


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
    # Do an early call to _makevars to get this info from xycoords.
    self._makevars()
    gids = self._gids
    # Get the Z grid for one half of the YY input.
    if self._yin: yy_ind = 0
    else: yy_ind = 1
    for gid in np.unique(gids):
      if gid<0: continue
      source_grid = rmn.decodeGrid(int(gid))
      dest_grid = source_grid['subgrid'][yy_ind]
      self._writeGrid (dest_grid)
      mask = (self._headers['ismeta']==0) & (self._headers['grtyp']==b'U')
      # Need to get ig1/ig2/ig3/ig4 from a sample record, because apparently
      # these values get mangles in librmn?
      rec_id = np.where(gids==gid)[0][0]
      ig1 = int(self._headers['ig1'][rec_id])
      ig2 = int(self._headers['ig2'][rec_id])
      ig3 = int(self._headers['ig3'][rec_id])
      ig4 = int(self._headers['ig4'][rec_id])
      submask = mask & (self._headers['ig1'] == ig1) & (self._headers['ig2'] == ig2) & (self._headers['ig3'] == ig3) & (self._headers['ig4'] == ig4)
      for key in ('grtyp','ni','nj','ig1','ig2','ig3','ig4'):
        self._headers[key][submask] = dest_grid[key]


  # Handle grid interpolations from raw binary array.
  def _decode (self, data, unused):
    if not self._yin and not self._yang:
      return super(YinYang,self)._decode (data, unused)
    prm = self._decode_headers(data[:72])
    prm = dict((k,v[0]) for k,v in prm.items())
    d = super(YinYang,self)._decode (data, unused).T
    if prm['grtyp'] == b'U' and self._yin:
      d = d[:,:prm['nj']//2]
    elif prm['grtyp'] == b'U' and self._yang:
      d = d[:,prm['nj']//2:]
    return d.T

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
    crop_to_smallest_grid : bool, optional
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

  # Handle cropping from raw binary array.
  def _decode (self, data, rec_id):
    import numpy as np
    d = super(Crop,self)._decode (data, rec_id)
    if not self._crop_to_smallest_grid: return d
    # Check if cropping necessary.
    ni = self._headers['ni'][rec_id]
    nj = self._headers['nj'][rec_id]
    if d.shape == (nj,ni): return d
    i0 = self._headers['i0'][rec_id]
    iN = self._headers['iN'][rec_id]
    j0 = self._headers['j0'][rec_id]
    jN = self._headers['jN'][rec_id]
    d = d[j0:jN,i0:iN]
    return d

