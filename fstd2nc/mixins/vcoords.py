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

# Decode ip1 information
from fstd2nc.mixins import vectorize
@vectorize
def decode_ip1 (ip1):
  from rpnpy.librmn.fstd98 import DecodeIp
  import numpy as np
  dtype = np.dtype([('kind','int32'),('level','float32')])
  out = np.empty(1,dtype=dtype)
  r1, r2, r3 = DecodeIp(ip1,0,0)
  out['kind'] = r1.kind
  if r1.v1 == 0 and r1.v2 != 0 :
    out['level'] = r1.v2
  else:
    out['level'] = r1.v1
  return out



#################################################
# Mixin for handling vertical coordinates.

class VCoords (BufferBase):

  @classmethod
  def _cmdline_args (cls, parser):
    super(VCoords,cls)._cmdline_args(parser)
    parser.add_argument('--strict-vcoord-match', action='store_true', help=_("Require the IP1/IP2/IP3 parameters of the vertical coordinate to match the IG1/IG2/IG3 paramters of the field in order to be used.  The default behaviour is to use the vertical record anyway if it's the only one in the file."))
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--diag-as-model-level', action='store_true', help=_("Treat diagnostic (near-surface) data as model level '1.0'.  This is the default behaviour."))
    group.add_argument('--split-diag-level', action='store_true', help=_("Put the diagnostic (near-surface) data in a separate variable, away from the 3D model output.  Suffices will be added to distinguish the different types of levels (i.e. _diag_level and _model_levels for diagnostic height and hybrid levels respectively)."))
    group.add_argument('--ignore-diag-level', action='store_true', help=_("Ignore data on diagnostic (near-surface) height."))
    parser.add_argument('--thermodynamic-levels', '--tlev', action='store_true', help=_("Only convert data that's on 'thermodynamic' vertical levels."))
    parser.add_argument('--momentum-levels', '--mlev', action='store_true', help=_("Only convert data that's on 'momentum' vertical levels."))
    parser.add_argument('--vertical-velocity-levels', '--wlev', action='store_true', help=_("Only convert data that's on 'vertical velocity' levels."))

  def __init__ (self, *args, **kwargs):
    """
    strict_vcoord_match : bool, optional
        Require the IP1/IP2/IP3 parameters of the vertical coordinate to
        match the IG1/IG2/IG3 paramters of the field in order to be used.
        The default behaviour is to use the vertical record anyway if it's
        the only one in the file.
    diag_as_model_level : bool, optional
        Treat diagnostic (near-surface) data as model level '1.0'.
        This is the default behaviour.
    split_diag_level : bool, optional
        Put the diagnostic (near-surface) data in a separate variable, away
        from the 3D model output.  Suffices will be added to distinguish
        the different types of levels (i.e. _diag_level and _model_levels for
        diagnostic height and hybrid levels respectively).
    ignore_diag_level : bool, optional
        Ignore data on diagnostic (near-surface) height.
    thermodynamic_levels : bool, optional
        Only convert data that's on 'thermodynamic' vertical levels.
    momentum_levels : bool, optional
        Only convert data that's on 'momentum' vertical levels.
    vertical_velocity_levels : bool, optional
        Only convert data that's on 'vertical velocity' levels.
    """
    from collections import OrderedDict
    import numpy as np
    from rpnpy.vgd.base import vgd_fromlist, vgd_get, vgd_free
    from rpnpy.vgd import VGDError
    from rpnpy.librmn.fstd98 import fstinl, fstprm, fstluk
    self._strict_vcoord_match = kwargs.pop('strict_vcoord_match',False)
    self._diag_as_model_level = kwargs.pop('diag_as_model_level',False)
    self._split_diag_level = kwargs.pop('split_diag_level',False)
    self._ignore_diag_level = kwargs.pop('ignore_diag_level',False)
    self._thermodynamic_levels = kwargs.pop('thermodynamic_levels',False) or kwargs.pop('tlev',False)
    self._momentum_levels = kwargs.pop('momentum_levels',False) or kwargs.pop('mlev',False)
    self._vertical_velocity_levels = kwargs.pop('vertical_velocity_levels',False) or kwargs.pop('wlev',False)

    # Use decoded IP1 values as the vertical axis.
    self._outer_axes = ('level',) + self._outer_axes
    # Tell the decoder not to process vertical records as variables.
    self._meta_records = self._meta_records + (b'!!',b'!!SF')
    self._maybe_meta_records = self._maybe_meta_records + (b'HY',)
    super(VCoords,self).__init__(*args,**kwargs)
    # Don't group records across different level 'kind'.
    # (otherwise can't create a coherent vertical axis).
    self._var_id = self._var_id + ('kind',)
    self._human_var_id = self._human_var_id + ('%(special_level)s','%(level_descr)s','vgrid%(kind)s')

    # Pre-scan the raw headers for special vertical records.
    # (these aren't available in the data stream, because we told the decoder
    # to ignore them).
    vrecs = OrderedDict()
    # Keep track of any diagnostic levels, e.g. from GEM.
    diag_ip1 = []
    for vcoord_nomvar in ('HY','!!'):
      for handle in fstinl(self._meta_funit, nomvar=vcoord_nomvar):
        header = fstprm(handle)
        key = (header['ip1'],header['ip2'])
        # For old HY records, there's no matching ipX/igX codes.
        if header['nomvar'].strip() == 'HY': key = 'HY'
        if key in vrecs: continue
        vrecs[key] = header
        # Extract diagnostic-level value from vertical coordinate.
        if key == 'HY': continue  # not applicable for HY coordinates.
        prm = fstluk(header,rank=3)
        try:
          vgd_id = vgd_fromlist(prm['d'])
        except VGDError:
          # Problems opening this vertical coordinate?
          continue
        try:
          diag_ip1.append(vgd_get(vgd_id,'DIPT'))
          diag_ip1.append(vgd_get(vgd_id,'DIPM'))
          diag_ip1.append(vgd_get(vgd_id,'DIPW'))
        except (KeyError,VGDError):
          #warn(_("Unable to parse diagnostic levels from the vertical coordinate"))
          pass
        vgd_free (vgd_id)
        # Done exracting diagnostic levels.

    self._vrecs = vrecs

    # Check if there are multiple '!!' records that all contain
    # identical coordinates.
    # In that case, only need to keep one copy in the lookup table.
    # This enables us to use that one unique coordinate for unstrict matching.
    if len(self._vrecs) > 1 and not self._strict_vcoord_match:
      unique_vrecs = OrderedDict()
      handled = set()
      for key, header in self._vrecs.items():
        coords = tuple(fstluk(header)['d'].flatten())
        if coords in handled: continue
        handled.add(coords)
        unique_vrecs[key] = header
      if len(unique_vrecs) == 1:
        self._vrecs = unique_vrecs
      del unique_vrecs, handled, coords

    fields = self._headers
    # Provide 'level' and 'kind' information to the decoder.
    decoded = np.concatenate(decode_ip1(fields['ip1']))
    # Only use first set of levels (can't handle ranges yet).
    fields['level'] = decoded['level']
    fields['kind'] = decoded['kind']

    # Pre-filter the diagnostic level data?
    # Start by treating diagnostic level as model level, then
    # revert this later if it ends up not working.
    if not self._split_diag_level:
      for ip1_val in diag_ip1:
          mask = (fields['ip1'] == ip1_val)
          if self._ignore_diag_level:
            fields['dltf'] |= mask
          else:
            #fields['ip1'][mask] = 93423264
            fields['level'][mask] = 1.0
            fields['kind'][mask] = 5

    # Keep only "thermodynamic" levels?
    if self._thermodynamic_levels:
      for key, header in self._vrecs.items():
        if key == 'HY': continue  # not applicable for HY coordinates.
        prm = fstluk(header,rank=3)
        vgd_id = vgd_fromlist(prm['d'])
        try:
          # Don't touch non-model levels.
          mask = (fields['kind'] != 5) & (fields['kind'] != 21)
          for ip1_t in vgd_get(vgd_id,'VIPT'):
            # Add this level to the mask (unless it's right at the surface).
            mask |= (fields['ip1'] == ip1_t) & (fields['level'] != 1.0) & (fields['level'] != 0.0)
          fields['dltf'] |= ~mask
        except (KeyError,VGDError):
          warn(_("Unable to parse thermodynamic levels from the vertical coordinate"))
        vgd_free (vgd_id)
    # Keep only "momentum" levels?
    if self._momentum_levels:
      for key, header in self._vrecs.items():
        if key == 'HY': continue  # not applicable for HY coordinates.
        prm = fstluk(header,rank=3)
        vgd_id = vgd_fromlist(prm['d'])
        try:
          # Don't touch non-model levels.
          mask = (fields['kind'] != 5) & (fields['kind'] != 21)
          for ip1_m in vgd_get(vgd_id,'VIPM'):
            # Add this level to the mask (unless it's right at the surface).
            mask |= (fields['ip1'] == ip1_m) & (fields['level'] != 1.0) & (fields['level'] != 0.0)
          fields['dltf'] |= ~mask
        except (KeyError,VGDError):
          warn(_("Unable to parse momentum levels from the vertical coordinate"))
        vgd_free (vgd_id)
    # Keep only "vertical velocity" levels?
    if self._vertical_velocity_levels:
      for key, header in self._vrecs.items():
        if key == 'HY': continue  # not applicable for HY coordinates.
        prm = fstluk(header,rank=3)
        vgd_id = vgd_fromlist(prm['d'])
        try:
          # Don't touch non-model levels.
          mask = (fields['kind'] != 5) & (fields['kind'] != 21)
          for ip1_w in vgd_get(vgd_id,'VIPW'):
            # Add this level to the mask (unless it's right at the surface).
            mask |= (fields['ip1'] == ip1_w) & (fields['level'] != 1.0) & (fields['level'] != 0.0)
          fields['dltf'] |= ~mask
        except (KeyError,VGDError):
          warn(_("Unable to parse vertical velocity levels from the vertical coordinate"))
        vgd_free (vgd_id)

    # Apply special labels for some levels.
    # Used to describe a variable when it's split into multiple coordinate types.
    fields['level_descr'] = np.empty(len(fields['kind']), dtype='|S20')
    fields['level_descr'][:] = 'unknown_levels'
    fields['level_descr'][fields['kind'] == 0] = 'depth_levels'
    fields['level_descr'][fields['kind'] == 1] = 'model_levels'
    fields['level_descr'][fields['kind'] == 2] = 'pressure_levels'
    fields['level_descr'][fields['kind'] == 3] = 'generic_levels'
    fields['level_descr'][fields['kind'] == 4] = 'height_levels'
    fields['level_descr'][fields['kind'] == 5] = 'model_levels'
    fields['level_descr'][fields['kind'] == 21] = 'model_levels'
    fields['special_level'] = np.empty(len(fields['kind']), dtype='|S20')
    fields['special_level'][:] = fields['level_descr'][:]
    # Discern ocean depths from ocean bottom.
    fields['special_level'][fields['kind'] == 0] = 'depth_levels'
    fields['special_level'][(fields['kind'] == 1) & (fields['level'] == 1.0)] = 'bottom_level'
    # Discern atmospheric levels from diagnostic levels.
    for ip1_val in diag_ip1:
      fields['special_level'][fields['ip1'] == ip1_val] = 'diag_level'

  def _makevars (self):
    super(VCoords,self)._makevars()
    # Call this makevars routine from a separate wrapper, to make it easier
    # to re-run just that one part without redoing the whole pipeline.
    # I.e., for case where our initial attempt at fitting the diagnostic level
    # into the model levels fails.
    try:
      return self._vcoords_makevars()
    except ValueError:
      super(VCoords,self)._makevars()
      return self._vcoords_makevars()

  def _vcoords_makevars (self):
    from fstd2nc.mixins import _var_type, _axis_type, _dim_type
    from collections import OrderedDict
    import numpy as np
    from rpnpy.vgd.base import vgd_fromlist, vgd_get, vgd_free
    from rpnpy.vgd.const import VGD_KEYS, VGD_OPR_KEYS
    from rpnpy.vgd import VGDError
    from rpnpy.librmn.fstd98 import fstinl, fstprm, fstluk

    vrecs = self._vrecs

    # If this becomes True by the end, then we need to rerun this routine
    # again with updated information.
    rerun = False

    # Scan through the data, and look for any use of vertical coordinates.
    vaxes = OrderedDict()
    prefs = dict()  # Reference pressure scalar variables.
    for var in self._varlist:
      level_axis = var.getaxis('level')
      if level_axis is None: continue
      # Degenerate vertical axis (0mb or 0hy?)
      if 'ip1' in var.atts and var.atts['ip1'] in (0, 99614720):
        if len(level_axis) == 1:
          i = var.dims.index('level')
          var.axes.pop(i)
          # Remove the axis if it's an "outer" axis (not along ni,nj,nk).
          if i < var.record_id.ndim:
            var.record_id = var.record_id.squeeze(axis=i)
          continue
      # No vertical axis?
      if 'kind' not in var.atts:
        continue
      # Dummy vertical dimension (no coordinate values?)
      if isinstance(level_axis,_dim_type):
        continue
      # Decode the vertical coordinate.
      levels = tuple(level_axis.array)
      kind = var.atts['kind']

      # Get version info for coordinate.
      key = (var.atts['ig1'],var.atts['ig2'])
      # If we're dealing with the old 'HY' records, then we don't match on
      # ig1/ig2.
      if key not in vrecs and 'HY' in vrecs:
        key = 'HY'
      # If there's only a single vertical record, then match that regardless of keys.
      if len(vrecs) == 1 and not self._strict_vcoord_match:
        key = list(vrecs.keys())[0]
      version = None
      if key in vrecs:
        code = vrecs[key]['ig1']
        if code // 1000 == kind:
          version = code % 1000

      # Special case - detect code 1002 (eta instead of sigma).
      if kind == 1 and 'HY' in vrecs:
        version = 2

      # Only need to provide one copy of the vertical axis.
      if (id(level_axis),kind,version) not in vaxes:

        # Keep track of any extra arrays needed for this axis.
        coordinates = []
        internal_atts = OrderedDict()
        atts = OrderedDict()
        name = 'level'
        new_axis = _axis_type(name, atts, level_axis.array)

        # Check if we have a vertical coordinate record to use.
        if key in vrecs:
          header = vrecs[key]
          # Add in metadata from the coordinate.
          atts.update(self._get_header_atts(header))
          # Add type-specific metadata.
          if header['nomvar'].strip() == '!!':
            # Get A and B info.
            try:
              vgd_id = vgd_fromlist(fstluk(header,rank=3)['d'])
              # Check if SLEVE coordinate.
              try:
                vgd_get(vgd_id, 'RFLS')
                sleve = True
              except VGDError:
                sleve = False
            except VGDError:
              warn (_("Problem decoding !! record."))
              continue
            # Partial definition of a/b coordinates for formula_terms reference.
            coordA = _var_type('a', {}, [new_axis], None)
            coordB = _var_type('b', {}, [new_axis], None)
            coordC = _var_type('c', {}, [new_axis], None)

            # First, patch up VGD_OPR_KEYS to include SLEVE coordinates?
            #TODO: remove this once rpnpy has full support for this.
            original_opr_double_keys = list(VGD_OPR_KEYS['get_double_1d'])
            original_opr_float_keys = list(VGD_OPR_KEYS['get_float_1d'])
            if sleve:
              VGD_OPR_KEYS['get_double_1d'].extend(['CC_M','CC_T'])
            # Also, patch in support for vertical velocity levels?
            #TODO: remove this once rpnpy has full support for this.
            if (kind,version) == (21,2):
              VGD_OPR_KEYS['get_float_1d'].extend(['VCDW'])
              VGD_OPR_KEYS['get_double_1d'].extend(['CA_W','CB_W'])
              if sleve:
                VGD_OPR_KEYS['get_double_1d'].extend(['CC_W'])

            # Add all parameters for this coordinate.
            for vgdkey in VGD_KEYS:
              try:
                val = vgd_get(vgd_id,vgdkey)
                # Skip multidimensional arrays (can't be encoded as metadata).
                if getattr(val,'ndim',1) > 1: continue
                internal_atts[vgdkey] = val
              except (KeyError,VGDError):
                pass  # Some keys not available in some vgrids?
            # Restore the rpnpy tables after SLEVE "patch", or may get
            # segmentation fault if a non-sleve coordinate is subsequently
            # read.
            VGD_OPR_KEYS['get_double_1d'] = original_opr_double_keys
            VGD_OPR_KEYS['get_float_1d'] = original_opr_float_keys

            # Put this information in the final output file?
            for k,v in internal_atts.items():
              if self._rpnstd_metadata_list is None or k in self._rpnstd_metadata_list:
                atts[k] = v

            # Attempt to fill in A/B coefficients (if available).
            try:
              # Special case for "tlift" levels (code 5003)?
              # Why do I need to do this???
              if (kind,version) == (5,3):
                internal_atts['CA_M'][-2] = (internal_atts['CA_T'][-1]+internal_atts['CA_T'][-2])/2
                internal_atts['CB_M'][-2] = (internal_atts['CB_T'][-1]+internal_atts['CB_T'][-2])/2

              all_z = list(internal_atts['VCDM'])+list(internal_atts['VCDT'])
              all_a = list(internal_atts['CA_M'])+list(internal_atts['CA_T'])
              all_b = list(internal_atts['CB_M'])+list(internal_atts['CB_T'])
              if sleve:
                all_c = list(internal_atts['CC_M'])+list(internal_atts['CC_T'])
              if 'VCDW' in internal_atts:
                all_z.extend(list(internal_atts['VCDW']))
                all_a.extend(list(internal_atts['CA_W']))
                all_b.extend(list(internal_atts['CB_W']))
                if sleve:
                  all_c.extend(list(internal_atts['CC_W']))
              A = []
              B = []
              C = []
              for z in levels:
                ind = all_z.index(z)
                A.append(all_a[ind])
                B.append(all_b[ind])
                if sleve: C.append(all_c[ind])
              coordA.array = np.asarray(A)
              coordB.array = np.asarray(B)
              coordC.array = np.asarray(C)
            except (KeyError,ValueError,VGDError):
              warn (_("Unable to get A/B coefficients for %s.")%var.name)
              coordA = None
              coordB = None
              coordC = None
            vgd_free (vgd_id)
        ###

        # Get metadata that's specific to this axis.
        atts['axis'] = 'Z'
        # Reference: http://web-mrb.cmc.ec.gc.ca/science//si/eng/si/libraries/rmnlib/fstd/main.html#RTFToC11
        # Also: https://wiki.cmc.ec.gc.ca/wiki/Vgrid/vcode
        if kind == 0:
          if var.atts.get('typvar',None) == 'P@' :   # masked ocean variable
            name = 'depth'
            atts['standard_name'] = 'depth'
            atts['units'] = 'm'
            atts['positive'] = 'down'
          else:
            # height [m] (metres)
            name = 'height'
            atts['standard_name'] = 'height'
            atts['units'] = 'm'
            atts['positive'] = 'up'
        elif kind == 1 and version in (1,None):
          # sigma [sg] (0.0->1.0)
          atts['standard_name'] = 'atmosphere_sigma_coordinate'
          atts['positive'] = 'down'
          atts['formula'] = 'p = sigma * ps'
          atts['formula_terms'] = OrderedDict([('sigma',new_axis),('ps','P0')])
        elif kind == 1 and version == 2:
          # eta levels
          atts['standard_name'] = 'atmosphere_sigma_coordinate'
          atts['positive'] = 'down'
          atts['formula'] = 'p = sigma * (ps-ptop) + ptop'
          atts['formula_terms'] = OrderedDict([('sigma',new_axis),('ps','P0'),('ptop','PT')])
        elif kind == 2:
          # pressure [mb] (millibars)
          name = 'pres'
          atts['standard_name'] = 'air_pressure'
          atts['units'] = 'hPa'
          atts['positive'] = 'down'
        elif kind == 3:
          # arbitrary code
          name = 'sfclevel'
          atts['standard_name'] = 'model_level_number'
          atts['positive'] = 'down'
        elif kind == 4:
          # height [M] (metres) with respect to ground level
          name = 'height'
          atts['standard_name'] = 'height'
          atts['units'] = 'm'
          atts['positive'] = 'up'
        elif kind == 5:
          # hybrid coordinates [hy] (0.0->1.0)
          atts['positive'] = 'down'
          if key in vrecs:
            if header['nomvar'].strip() == '!!':
              if internal_atts['LOGP']:
                # Not really a "standard" name, but there's nothing in the
                # CF conventions document on how to encode this.
                # I just merged the atmosphere_ln_pressure_coordinate and
                # atmosphere_hybrid_sigma_pressure_coordinate together.
                # http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#dimensionless-v-coord
                atts['standard_name'] = 'atmosphere_hybrid_sigma_ln_pressure_coordinate'
                # Document the formula to follow, since it's not in the conventions.
                #TODO: update this once there's an actual convention to follow!
                if version in (1,2,3,4,5) and None not in (coordA,coordB):
                  atts['formula'] = "p = exp(a+b*log(ps/pref)) / 100.0"
                  atts['formula_terms'] = OrderedDict([('a',coordA),('b',coordB),('ps','P0'),('pref','pref')])
                  coordinates.extend([coordA,coordB])
                elif sleve and None not in (coordA,coordB,coordC):
                  atts['formula'] = "p = exp(a+b*log(ps/pref)+c*log(psl/pref)) / 100.0"
                  atts['formula_terms'] = OrderedDict([('a',coordA),('b',coordB),('c',coordC),('ps','P0'),('psl','P0LS'),('pref','pref')])
                  coordinates.extend([coordA,coordB,coordC])
                # Try getting reference pressure as a scalar.
                try:
                  pref = internal_atts['PREF'] / 100.0
                  if pref not in prefs:
                    prefs[pref] = _var_type('pref',{'units':'hPa'},[],np.array(pref))
                  atts['formula_terms']['pref'] = prefs[pref]
                except (KeyError,VGDError):
                  pass # Don't have PREF available for some reason?
              else:
                atts['standard_name'] = 'atmosphere_hybrid_sigma_pressure_coordinate'
                if None not in (coordA,coordB):
                  coordA.array /= 100.0  # Use hPa for final units.
                  atts['formula'] = 'p = ap + b*ps'
                  coordA.name = 'ap'
                  atts['formula_terms'] = OrderedDict([('ap',coordA),('b',coordB),('ps','P0')])
                  coordinates.extend([coordA,coordB])

            # Not a '!!' coordinate, so must be 'HY'?
            else:
              atts['standard_name'] = 'atmosphere_hybrid_sigma_pressure_coordinate'
              # Get A and B info.
              eta = np.asarray(levels)
              ptop = decode_ip1(header['ip1'])['level']
              # Conversion taken from 'ig_to_hybref' function in librmn:
              pref = float(header['ig1'])
              rcoef = header['ig2']/1000.0
              # Apply the formula to compue A & B (from old fstd_core.c code):
              etatop = ptop/pref
              B = ((eta - etatop) / (1 - etatop)) ** rcoef
              A = pref * (eta - B)
              coordA = _var_type('ap', {}, [new_axis], np.asarray(A))
              coordB = _var_type('b', {}, [new_axis], np.asarray(B))
              coordinates.extend([coordA,coordB])
              atts['formula'] = 'p = ap + b*ps'
              atts['formula_terms'] = OrderedDict([('ap',coordA),('b',coordB),('ps','P0')])
              # Add extra HY record metadata.
              atts.update(ptop=ptop, rcoef=rcoef, pref=pref)
        elif kind == 6:
          # theta [th]
          name = 'theta'
          atts['standard_name'] = 'air_potential_temperature'
          atts['units'] = 'K'
          atts['positive'] = 'up'
        elif kind == 7:
          name = 'depth'
          atts['standard_name'] = 'depth'
          atts['units'] = 'm'
          atts['positive'] = 'down'
        elif kind == 21:
          if sleve:
            atts['standard_name'] = 'atmosphere_sleve_coordinate'
            atts['positive'] = 'up'
            if None not in (coordA,coordB,coordC):
              atts['formula'] = 'z = az + b1*zsurf1 + b2*zsurf2'
              coordA.name = 'az'
              coordB.name = 'b2'
              coordC.name = 'b1'
              atts['formula_terms'] = OrderedDict([('az',coordA),('b1',coordC),('zsurf1','MELS'),('b2',coordB),('zsurf2','ME')])
              coordinates.extend([coordA,coordB,coordC])

          else:
            atts['standard_name'] = 'atmosphere_hybrid_height_coordinate'
            atts['positive'] = 'up'
            if None not in (coordA,coordB):
              atts['formula'] = 'z = a + b*orog'
              atts['formula_terms'] = OrderedDict([('a',coordA),('b',coordB),('orog','ME')])
              coordinates.extend([coordA,coordB])

        # Add this vertical axis.
        if len(coordinates) > 0:
          atts['coordinates'] = coordinates
        # Update axis name.
        new_axis.name = name
        # Now have a fully defined axis to use.
        vaxes[(id(level_axis),kind,version)] = new_axis

      # Set the vertical axis for this variable.
      vaxis = vaxes[(id(level_axis),kind,version)]
      var.axes[var.dims.index('level')] = vaxis

    # Detect mixture of diagnostic / model levels, and print a warning.
    model_vars = set()
    diag_vars = set()
    for var in self._varlist:
      dims = var.dims
      if 'level' in var.dims:
        model_vars.add(var.name)
      if 'height' in var.dims:
        diag_vars.add(var.name)
    mixed_vars = model_vars & diag_vars
    if len(mixed_vars) > 0:
      warn (_("Mixture of model / height levels found.  This will cause multiple definitions of variables in the output.  If this is undesired, you could try using --ignore-diag-level or --diag-as-model-level."))

    # Detect different timesteps for diagnostic / model levels, and split
    # those fields.
    # Also detect where *only* diagnostic level was found, and revert that
    # back to the correct type.
    if not self._diag_as_model_level:
      for var in self._varlist:
        if not hasattr(var,'record_id'): continue
        valid_records = var.record_id[var.record_id>=0]
        # Ugly hack to find all related records.
        # There may be duplicate records that aren't in this set, but would
        # need to be transformed identically.
        # E.g., if have multiple copies of the diagnostic level, then we want
        # them all to stay equivalent or suddenly get multiple versions of the
        # field showing up on the 'rerun' iteration.
        nomvar = self._headers['nomvar'][valid_records][0]
        valid_records = np.where(self._headers['nomvar'] == nomvar)[0]
        ip1 = self._headers['ip1'][valid_records]
        decoded = np.concatenate(decode_ip1(ip1))
        if not any(decoded['kind']==4): continue
        if np.all(decoded['kind'] == self._headers['kind'][valid_records]): continue
        if np.any(var.record_id<0):
          warn (_("Having trouble treating %s diagnostic level as model level - splitting into a separate field.")%var.name)
          self._headers['kind'][valid_records] = decoded['kind']
          self._headers['level'][valid_records] = decoded['level']
          rerun = True
        # For single height level (no vertical structure), reset the level
        # type.
        if var.record_id.size == 1:
          self._headers['kind'][valid_records] = decoded['kind']
          self._headers['level'][valid_records] = decoded['level']
          rerun = True

    if rerun:
      raise ValueError
