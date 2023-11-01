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


# Define a lock for controlling threaded access to ezscint, etc.
from threading import RLock
_lock = RLock()
del RLock


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
      new_recs.append(dict(prm, nomvar='>>  ', ni=gid['ni'], d=gid['ax'], dtype=gid['ax'].dtype))
    if 'ay' in gid:
      new_recs.append(dict(prm, nomvar='^^  ', nj=gid['nj'], d=gid['ay'], dtype=gid['ay'].dtype))
    if 'axy' in gid:
      new_recs.append(dict(prm, nomvar='^>  ', ni=gid['ni'], nj=gid['nj'], d=gid['axy'], dtype=gid['axy'].dtype))
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
  # Extract interpolation grid.
  def number(x):
    try: return int(x)
    except ValueError: return float(x)
  # Degenerate case: already have the grid info in a dictionary.
  if hasattr(interp,'keys'):
    return interp
  # Degenerate case: have a grid id, just need to look up the info.
  if isinstance(interp,int):
    with _lock:
      return rmn.decodeGrid(interp)
  # Main case: Have grid info encoded as a string.
  with _lock:
    interp = interp.split(',')
    grtyp = interp[0]
    grid_args = [number(v) for v in interp[1:] if '=' not in v]
    grid_kwargs = dict(v.split('=') for v in interp[1:] if '=' in v)
    grid_kwargs = dict((k,number(v)) for k,v in grid_kwargs.items())
    if not hasattr(rmn,'defGrid_'+grtyp):
      error(_("Unknown grid '%s'")%grtyp)
    return getattr(rmn,'defGrid_'+grtyp)(*grid_args,**grid_kwargs)

# Helper method - pack up grid id for transport to other processes.
# Converts to a dictionary, then to a tuple of key/value pairs.
from fstd2nc.mixins import vectorize
@vectorize
def _pack_grid (gid):
  import rpnpy.librmn.all as rmn
  def _pack_item (item):
    k, v = item
    if k == 'subgrid':
      return (k,tuple(map(_pack_dict,v)))
    if hasattr(v,'dumps'):  # Numpy array
      return (k,v.dumps())
    return (k,v)
  def _pack_dict (d):
    d = d.copy()
    d.pop('id',None)
    d.pop('subgridid',None)
    return tuple(map(_pack_item,d.items()))
  if gid < 0: return None
  with _lock:
    grid = rmn.decodeGrid(int(gid))
    return _pack_dict(grid)
del vectorize

# Helper method - unpack a grid after transport to another processor.
# Returns an active grid id that's valid for this processor.
_grid_lookup = {}
def _unpack_grid (grid):
  import rpnpy.librmn.all as rmn
  from pickle import loads
  def _unpack_item (item):
    k, v = item
    if k == 'subgrid':
      return (k,list(map(_unpack_dict,v)))
    if isinstance(v,bytes):  # Numpy array
      return (k,loads(v))
    return (k,v)
  def _unpack_dict (d):
    return dict(map(_unpack_item,d))
  if grid is None: return -1
  # Check if grid already handled.
  for gid, g in _grid_lookup.items():
    if g is grid: return gid
  # Otherwise, start unpacking it.
  packed_grid = grid
  grid = _unpack_dict(grid)
  with _lock:
    # Special case: super grid
    if 'subgrid' in grid:
      subgrids = [rmn.encodeGrid(subgrid)['id'] for subgrid in grid['subgrid']]
      gid = rmn.ezgdef_supergrid(grid['ni'], grid['nj'], grid['grtyp'], grid['grref'], grid['version'], subgrids)
    else:
      gid = rmn.encodeGrid(grid)['id']
  _grid_lookup[gid] = packed_grid
  return gid


class Interp (BufferBase):
  # List of known fields that should be treated as vectors for the purpose of
  # interpolation.
  _vector_fields=("UU,VV","UT1,VT1","URT1,VRT1","UUW,VVW","UU2W,VV2W","UUI,VVI","UDST,VDST","UUN,VVN","UUX,VVX")

  @classmethod
  def _cmdline_args (cls, parser):
    from argparse import SUPPRESS
    super(Interp,cls)._cmdline_args(parser)
    #parser.add_argument('--interp', metavar="GRTYP,PARAM=VAL,PARAM=VAL,...", help=_("Interpolate to the specified grid."))
    parser.add_argument('--interp', metavar="GRTYP,PARAM=VAL,PARAM=VAL,...", help=SUPPRESS)
    #parser.add_argument('--vector-fields', metavar="UFIELD,VFIELD", help=_("Specify a pair of fields to be treated as vectors.  Only needed if --interp is used on vector fields."))
    parser.add_argument('--vector-fields', metavar="UFIELD,VFIELD", help=SUPPRESS)
    #parser.add_argument('--scalar-fields', metavar="VAR1,VAR2,...", help=_("Disable vector interpolation for the specified fields."))
    parser.add_argument('--scalar-fields', metavar="VAR1,VAR2,...", help=SUPPRESS)

  def __init__ (self, *args, **kwargs):
    """
    interp : str or int, optional
        Interpolate to the specified grid.
    vector_fields : str or list, optional
        Specify a pair of fields to be treated as vectors.  Only needed if
        interp is used on vector fields.
    scalar_fields : str or list, optional
        Disable vector interpolation for the specified fields.
    """
    import rpnpy.librmn.all as rmn
    import numpy as np
    from fstd2nc.extra import structured_array
    from collections import OrderedDict

    self._decoder_data = self._decoder_data + (('u_data',('u_address','u_length','u_d')),('v_data',('v_address','v_length','v_d')))

    interp = kwargs.pop('interp',None)

    # Determine which fields should be interpolated as vectors.
    vector_fields = kwargs.pop('vector_fields',None)
    if vector_fields is None: vector_fields = []
    if isinstance(vector_fields,str):
      vector_fields = [vector_fields]
    vector_fields = tuple(vector_fields) + self._vector_fields
    vector_fields = [v.split(',') if isinstance(v,str) else v for v in vector_fields]
    scalar_fields = kwargs.pop('scalar_fields',None)
    if scalar_fields is None: scalar_fields = []
    if isinstance(scalar_fields,str):
      scalar_fields = scalar_fields.split(',')
    vector_fields = [(u,v) for u,v in vector_fields if u not in scalar_fields and v not in scalar_fields]
    b = lambda s: s.ljust(4).encode()
    vector_fields = [(b(u),b(v)) for u,v in vector_fields]

    super(Interp,self).__init__(*args,**kwargs)
    # Extract interpolation grid.
    if interp is not None:
      interp_grid = _get_interp_grid(interp)
      self._interp_grid = interp_grid
      # Hack the new grid descriptors into the headers.
      self._writeGrid (interp_grid)

      # Store original and modified versions of the grid descriptors.
      self._decoder_extra_args = self._decoder_extra_args + ('source_grid','dest_grid')
      self._ignore_atts = self._ignore_atts + ('source_grid','dest_grid')
      self._headers['source_grid'] = np.empty(self._nrecs,dtype=object)

      # Set up some options for ezsint.
      import rpnpy.librmn.all as rmn
      rmn.ezsetopt (rmn.EZ_OPT_EXTRAP_DEGREE, rmn.EZ_EXTRAP_VALUE)
      rmn.ezsetopt (rmn.EZ_OPT_EXTRAP_VALUE, self._fill_value)

      # Locate wind components.
      is_wind = np.isin(self._headers['nomvar'],vector_fields)
      if not np.any(is_wind): return
      # Set U,V record info for the interpolator.
      nrecs = len(self._headers['name'])
      self._headers['u_address'] = np.empty(nrecs,'int32')
      self._headers['u_address'][:] = -1
      self._headers['u_length'] = np.empty(nrecs,'int32')
      self._headers['u_length'][:] = -1
      self._headers['v_address'] = np.empty(nrecs,'int32')
      self._headers['v_address'][:] = -1
      self._headers['v_length'] = np.empty(nrecs,'int32')
      self._headers['v_length'][:] = -1
      for u,v in vector_fields:
        this_uv = np.isin(self._headers['nomvar'],(u,v))
        if np.sum(this_uv) == 0: continue
        # The following is adapted from masks mixin
        nomvar = self._headers['nomvar'][this_uv]
        typvar = self._headers['typvar'][this_uv]
        etiket = self._headers['etiket'][this_uv]
        datev = self._headers['datev'][this_uv]
        ip1 = self._headers['ip1'][this_uv]
        ip2 = self._headers['ip2'][this_uv]
        ip3 = self._headers['ip3'][this_uv]
        dltf = self._headers['dltf'][this_uv]
        # Figure out how to pair up the wind components
        # Requires O(n log n) time, which is better than O(n^2) for naive lookup
        # on each record.
        keys = OrderedDict([('dltf',dltf),('etiket',etiket),('datev',datev),('ip1',ip1),('ip2',ip2),('ip3',ip3),('typvar',typvar),('nomvar',nomvar)])
        ind = np.argsort(structured_array(keys))
        nomvar = nomvar[ind]
        typvar = typvar[ind]
        etiket = etiket[ind]
        datev = datev[ind]
        ip1 = ip1[ind]
        ip2 = ip2[ind]
        ip3 = ip3[ind]
        dltf = dltf[ind]
        is_paired = (etiket[:-1] == etiket[1:]) & (datev[:-1] == datev[1:]) & (ip1[:-1] == ip1[1:]) & (ip2[:-1] == ip2[1:]) & (ip3[:-1] == ip3[1:])
        is_paired_u = is_paired & (nomvar[:-1] == u) & (nomvar[1:] == v)
        is_paired_u = np.where(is_paired_u)[0]
        is_paired_v = is_paired_u + 1  # Assume V component follows U
        u_recid = np.arange(nrecs)[this_uv][ind][is_paired_u]
        v_recid = np.arange(nrecs)[this_uv][ind][is_paired_v]
        # If this is U record, then V record goes into auxiliary data, and
        # vice versa.
        self._headers['v_address'][u_recid] = self._headers['address'][v_recid]
        self._headers['v_length'][u_recid] = self._headers['length'][v_recid]
        self._headers['u_address'][v_recid] = self._headers['address'][u_recid]
        self._headers['u_length'][v_recid] = self._headers['length'][u_recid]

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
    self._headers['source_grid'][:] = _pack_grid(self._gids)
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

  def _decoder_scalar_args (self):
    args = super(Interp,self)._decoder_scalar_args()
    if hasattr(self,'_interp_grid'):
      args['dest_grid'] = _pack_grid(self._interp_grid['id'])
    return args

  # Handle grid interpolations from raw binary array.
  @classmethod
  def _postproc (cls, data, u_data=None, v_data=None, source_grid=None, dest_grid=None, **kwargs):
    import rpnpy.librmn.all as rmn
    import numpy as np
    # If no interpolation requested, nothing to do.
    if source_grid is None or dest_grid is None:
      return super(Interp,cls)._postproc (data, **kwargs)
    source_grid = _unpack_grid(source_grid)
    dest_grid = _unpack_grid(dest_grid)
    if source_grid < 0:
      raise ValueError("Source data is not on a recognized grid.  Unable to interpolate.")
    # Handle other post-processing of the input fields that occur before this
    # mixin.
    d = super(Interp,cls)._postproc (data, **kwargs).T
    # Decode auxiliary field(s), if they're needed.
    if u_data is not None:
      u = super(Interp,cls)._postproc (u_data, **kwargs).T
    else:
      u = None
    if v_data is not None:
      v = super(Interp,cls)._postproc (v_data, **kwargs).T
    else:
      v = None
    with _lock:
      # Propogate any fill values to the interpolated grid.
      fill_value = kwargs.get('fill_value')
      in_mask = np.zeros(d.shape, order='F', dtype='float32')
      in_mask[d==fill_value] = 1.0
      # Is this a U component? (i.e. use this + extra V component?)
      if v is not None:
        d, _ = rmn.ezuvint (dest_grid, source_grid, d, v)
      # Is this a V component? (i.e. use this + extra U component?)
      elif u is not None:
        _, d = rmn.ezuvint (dest_grid, source_grid, u, d)
      else:
        d = rmn.ezsint (dest_grid, source_grid, d)
      out_mask = rmn.ezsint (dest_grid, source_grid, in_mask)
      d[out_mask!=0] = fill_value
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
    yin = kwargs.pop('yin',False)
    yang = kwargs.pop('yang',False)

    # Load the metadata from the file(s).
    super(YinYang,self).__init__(*args,**kwargs)

    # Update yin-yang records to appear as regular rotated grids.
    if not yin and not yang:
      return
    # Add yin/yang selection boolean flags as columns.
    self._decoder_extra_args = self._decoder_extra_args + ('yin','yang')
    if yin:
      self._headers['yin'] = (self._headers['grtyp'] == b'U')
    elif yang:
      self._headers['yang'] = (self._headers['grtyp'] == b'U')

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
    if yin: yy_ind = 0
    else: yy_ind = 1
    for gid in np.unique(gids):
      if gid<0: continue
      source_grid = rmn.decodeGrid(int(gid))
      if 'subgrid' not in source_grid: continue   # Not a YY grid?
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

  # Handle grid interpolations from raw binary array.
  @classmethod
  def _postproc (cls, data, yin=False, yang=False, **kwargs):
    if not yin and not yang:
      return super(YinYang,cls)._postproc (data, **kwargs)
    d = super(YinYang,cls)._postproc (data, **kwargs)
    nj, ni = d.shape[-2:]
    if yin:
      return d[:nj//2,:]
    elif yang:
      return d[nj//2:,:]

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
    self._decoder_extra_args = self._decoder_extra_args + ('crop_j0','crop_jN','crop_i0','crop_iN')
    self._headers['crop_j0'] = np.empty(self._nrecs,'int32')
    self._headers['crop_jN'] = np.empty(self._nrecs,'int32')
    self._headers['crop_i0'] = np.empty(self._nrecs,'int32')
    self._headers['crop_iN'] = np.empty(self._nrecs,'int32')
    # Initialize with default values, to avoid processing garbage indices.
    self._headers['crop_j0'][:] = 0
    self._headers['crop_jN'][:] = self._headers['nj']
    self._headers['crop_i0'][:] = 0
    self._headers['crop_iN'][:] = self._headers['ni']
    self._ignore_atts = self._ignore_atts + ('crop_j0','crop_jN','crop_i0','crop_iN')

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
        self._headers['crop_j0'][submask] = j0
        self._headers['crop_jN'][submask] = jN
        self._headers['crop_i0'][submask] = i0
        self._headers['crop_iN'][submask] = iN

  # Handle cropping of the data.
  @classmethod
  def _postproc (cls, data, crop_j0=None, crop_jN=None, crop_i0=None, crop_iN=None, **kwargs):
    d = super(Crop,cls)._postproc (data, **kwargs)
    if crop_j0 is not None and crop_jN is not None:
      d = d[crop_j0:crop_jN,:]
    if crop_i0 is not None and crop_iN is not None:
      d = d[:,crop_i0:crop_iN]
    return d


