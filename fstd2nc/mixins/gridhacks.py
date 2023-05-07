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
      if hasattr(v,'mask'):
        self._headers[k] = np.ma.zeros(shape=nrecs + len(new_recs), dtype=v.dtype)
      else:
        self._headers[k] = np.zeros(shape=nrecs + len(new_recs), dtype=v.dtype)
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


# Helper method - given an interpolation grid string, return a grid id.
def _get_interp_grid (interp):
  import rpnpy.librmn.all as rmn
  from fstd2nc.mixins.fstd import _lock
  # Extract interpolation grid.
  def number(x):
    try: return int(x)
    except ValueError: return float(x)
  with _lock:
    interp = interp.split(',')
    grtyp = interp[0]
    grid_args = [number(v) for v in interp[1:] if '=' not in v]
    grid_kwargs = dict(v.split('=') for v in interp[1:] if '=' in v)
    grid_kwargs = dict((k,number(v)) for k,v in grid_kwargs.items())
    if not hasattr(rmn,'defGrid_'+grtyp):
      error(_("Unknown grid '%s'")%grtyp)
    return getattr(rmn,'defGrid_'+grtyp)(*grid_args,**grid_kwargs)

# Keep track of valid grid ids (to detect if we have a problem with grid ids)
_valid_gids = set()

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
    if interp is not None:
      interp_grid = _get_interp_grid(interp)
      self._interp_grid = interp_grid
      # Hack the new grid descriptors into the headers.
      self._writeGrid (interp_grid)
      _valid_gids.add(interp_grid['id'])

      # Store original and modified versions of the grid descriptors.
      self._decoder_extra_args = self._decoder_extra_args + ('source_gid','dest_gid')
      self._ignore_atts = self._ignore_atts + ('source_gid','dest_gid')
      self._headers['source_gid'] = np.empty(self._nrecs,dtype=object)
      self._headers['dest_gid'] = np.empty(self._nrecs,dtype=object)
      ismeta = self._headers['ismeta']
      self._headers['dest_gid'][~ismeta] = interp_grid['id']

      # Set up some options for ezsint.
      import rpnpy.librmn.all as rmn
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
    self._headers['source_gid'][:] = np.array(self._gids)
    _valid_gids.update(g for g in self._gids if g >= 0)
    # Now, use interpolated grid descriptors.
    ismeta = self._headers['ismeta']
    for key in ('grtyp','ni','nj','ig1','ig2','ig3','ig4'):
      self._headers[key][:] = np.where(ismeta, self._headers[key], self._interp_grid[key])
    super(Interp,self)._makevars()

    # Add fill value to the data.
    # Modified from code in from mask mixin.
    # Set this fill value for any interpolated data, since it may contain
    # points outside the original boundary.
    for var in self._varlist:
      if isinstance(var,_iter_type):
        var.atts['_FillValue'] = var.dtype.type(self._fill_value)

  # Handle grid interpolations from raw binary array.
  def _decode (self, data, source_gid=None, dest_gid=None, **kwargs):
    import rpnpy.librmn.all as rmn
    import numpy as np
    from fstd2nc.mixins.fstd import _lock
    if source_gid is None or dest_gid is None:
      return super(Interp,self)._decode (data, **kwargs)
    if source_gid not in _valid_gids or dest_gid not in _valid_gids:
      error(_("Problem finding grid id.  It's possible that you're running this in a multi-processing environment, which does not support the 'interp' option."))
    # Retrieve an active librmn grid id associated with this grid.
    if source_gid < 0:
      raise ValueError("Source data is not on a recognized grid.  Unable to interpolate.")
    d = super(Interp,self)._decode (data, **kwargs).T
    with _lock:
      # Propogate any fill values to the interpolated grid.
      in_mask = np.zeros(d.shape, order='F', dtype='float32')
      in_mask[d==self._fill_value] = 1.0
      d = rmn.ezsint (dest_gid, source_gid, d)
      out_mask = rmn.ezsint (dest_gid, source_gid, in_mask)
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
    self._yin = kwargs.pop('yin',False)
    self._yang = kwargs.pop('yang',False)

    # Load the metadata from the file(s).
    super(YinYang,self).__init__(*args,**kwargs)

    # Update yin-yang records to appear as regular rotated grids.
    if not self._yin and not self._yang:
      return

    # Run _makevars early to generate grid ids with xycoords mixin.
    # Silence warnings from makevars, which might not be relevant to the final
    # state after we make our grid modifications.
    import fstd2nc
    streams = fstd2nc.stdout.streams
    fstd2nc.stdout.streams = ()
    self._makevars()
    fstd2nc.stdout.streams = streams
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
      # these values get mangles in librmn for U grids?
      rec_id = np.where(gids==gid)[0][0]
      ig1 = int(self._headers['ig1'][rec_id])
      ig2 = int(self._headers['ig2'][rec_id])
      ig3 = int(self._headers['ig3'][rec_id])
      ig4 = int(self._headers['ig4'][rec_id])
      submask = mask & (self._headers['ig1'] == ig1) & (self._headers['ig2'] == ig2) & (self._headers['ig3'] == ig3) & (self._headers['ig4'] == ig4)
      for key in ('grtyp','ni','nj','ig1','ig2','ig3','ig4'):
        self._headers[key][submask] = dest_grid[key]

  def _decoder_scalar_args (self):
    kwargs = super(YinYang,self)._decoder_scalar_args()
    if self._yin: kwargs['yin'] = True
    if self._yang: kwargs['yang'] = True
    return kwargs

  # Handle grid interpolations from raw binary array.
  def _decode (self, data, yin=False, yang=False, **kwargs):
    if not yin and not yang:
      return super(YinYang,self)._decode (data, **kwargs)
    prm = self._decode_headers(data[:72])
    prm = dict((k,v[0]) for k,v in prm.items())
    d = super(YinYang,self)._decode (data, **kwargs).T
    if prm['grtyp'] == b'U' and yin:
      d = d[:,:prm['nj']//2]
    elif prm['grtyp'] == b'U' and yang:
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
    self._crop_to_smallest_grid = kwargs.pop('crop_to_smallest_grid',False)

    # Load the metadata from the file(s).
    super(Crop,self).__init__(*args,**kwargs)

    if not self._crop_to_smallest_grid:
      return

    # Keep track of cropping regions for the data.
    self._decoder_extra_args = self._decoder_extra_args + ('crop_j','crop_i')
    self._headers['crop_j'] = np.empty(self._nrecs,object)
    self._headers['crop_i'] = np.empty(self._nrecs,object)
    self._ignore_atts = self._ignore_atts + ('crop_j','crop_i')

    # Run _makevars early to generate grid ids with xycoords mixin.
    # Silence warnings from makevars, which might not be relevant to the final
    # state after we make our grid modifications.
    import fstd2nc
    streams = fstd2nc.stdout.streams
    fstd2nc.stdout.streams = ()
    self._makevars()
    fstd2nc.stdout.streams = streams
    gids = self._gids

    # Find all available grids, group them by their projection parameters.
    # I.e., lump together all grids that only differ in their x/y extent.
    matches = dict()
    for gid in np.unique(gids):
      if gid<0: continue
      grid = rmn.decodeGrid(int(gid))
      # Skip supergrids.
      if grid['grtyp'] == 'U': continue
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
    for key in smallest.keys():
      smallest_grid = matches[key][smallest[key]]
      for grid in matches[key]:
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
        submask = (self._headers['ig1'] == grid['tag1']) & (self._headers['ig2'] == grid['tag2']) & (self._headers['ig3'] == grid['tag3'])
        for key in ('grtyp','ni','nj','ig1','ig2','ig3','ig4'):
          self._headers[key][submask] = smallest_grid[key]
        self._headers['crop_j'][submask] = slice(j0,jN)
        self._headers['crop_i'][submask] = slice(i0,iN)

  # Handle cropping from raw binary array.
  def _decode (self, data, crop_j=None, crop_i=None, **kwargs):
    d = super(Crop,self)._decode (data, **kwargs)
    if crop_j is not None:
      d = d[crop_j,:]
    if crop_i is not None:
      d = d[:,crop_i]
    return d


