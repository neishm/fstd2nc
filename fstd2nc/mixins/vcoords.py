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
  _vcoord_nomvars = ('HY','!!')

  @classmethod
  def _cmdline_args (cls, parser):
    super(VCoords,cls)._cmdline_args(parser)
    parser.add_argument('--strict-vcoord-match', action='store_true', help=_("Require the IP1/IP2/IP3 parameters of the vertical coordinate to match the IG1/IG2/IG3 paramters of the field in order to be used.  The default behaviour is to use the vertical record anyway if it's the only one in the file."))
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--diag-as-model-level', action='store_true', help=_("Treat diagnostic (near-surface) data as model level '1.0'.  Normally, this data goes in a separate variable because it has incompatible units for the vertical coordinate.  Use this option if your variables are getting split with suffixes '_vgrid4' and '_vgrid5', and you'd rather keep both sets of levels together in one variable."))
    group.add_argument('--ignore-diag-level', action='store_true', help=_("Ignore data on diagnostic (near-surface) height."))

  # Need to extend _headers_dtype before __init__.
  def __new__ (cls, *args, **kwargs):
    obj = super(VCoords,cls).__new__(cls, *args, **kwargs)
    obj._headers_dtype = obj._headers_dtype + [('level','float32'),('kind','int32')]
    return obj

  def __init__ (self, *args, **kwargs):
    """
    strict_vcoord_match : bool, optional
        Require the IP1/IP2/IP3 parameters of the vertical coordinate to
        match the IG1/IG2/IG3 paramters of the field in order to be used.
        The default behaviour is to use the vertical record anyway if it's
        the only one in the file.
    diag_as_model_level : bool, optional
        Treat diagnostic (near-surface) data as model level '1.0'.
        Normally, this data goes in a separate variable because it has
        incompatible units for the vertical coordinate.  Use this option if
        your variables are getting split with suffixes '_vgrid4' and
        '_vgrid5', and you'd rather keep both sets of levels together in one
        variable.
    ignore_diag_level : bool, optional
        Ignore data on diagnostic (near-surface) height.
    """
    from collections import OrderedDict
    import numpy as np
    from rpnpy.vgd.base import vgd_fromlist, vgd_get, vgd_free
    from rpnpy.vgd import VGDError
    from rpnpy.librmn.fstd98 import fstinl, fstprm, fstluk
    self._strict_vcoord_match = kwargs.pop('strict_vcoord_match',False)
    self._diag_as_model_level = kwargs.pop('diag_as_model_level',False)
    self._ignore_diag_level = kwargs.pop('ignore_diag_level',False)

    # Use decoded IP1 values as the vertical axis.
    self._outer_axes = ('level',) + self._outer_axes
    # Tell the decoder not to process vertical records as variables.
    self._meta_records = self._meta_records + self._vcoord_nomvars
    super(VCoords,self).__init__(*args,**kwargs)
    # Don't group records across different level 'kind'.
    # (otherwise can't create a coherent vertical axis).
    self._var_id = self._var_id + ('kind',)
    self._human_var_id = self._human_var_id + ('vgrid%(kind)s',)

    # Pre-scan the raw headers for special vertical records.
    # (these aren't available in the data stream, because we told the decoder
    # to ignore them).
    vrecs = OrderedDict()
    for vcoord_nomvar in self._vcoord_nomvars:
      for handle in fstinl(self._meta_funit, nomvar=vcoord_nomvar):
        header = fstprm(handle)
        key = (header['ip1'],header['ip2'])
        # For old HY records, there's no matching ipX/igX codes.
        if header['nomvar'].strip() == 'HY': key = 'HY'
        if key in vrecs: continue
        vrecs[key] = header
    self._vrecs = vrecs

    # Pre-filter the diagnostic level data?
    if self._diag_as_model_level or self._ignore_diag_level:
      for key, header in self._vrecs.items():
        prm = fstluk(header,rank=3)
        vgd_id = vgd_fromlist(prm['d'])
        try:
          ip1_t = vgd_get(vgd_id,'DIPT')
          ip1_m = vgd_get(vgd_id,'DIPM')
          mask = (self._headers['ip1'] == ip1_t) | (self._headers['ip1'] == ip1_m)
          if len(self._vrecs) > 1 or self._strict_vcoord_match:
            mask &= (self._headers['ig1'] == prm['ip1']) & (self._headers['ig2'] == prm['ip2'])
          if self._diag_as_model_level:
            self._headers['ip1'][mask] = 93423264
          elif self._ignore_diag_level:
            self._headers['dltf'] |= mask
        except (KeyError,VGDError):
          warn(_("Unable to parse diagnostic levels from the vertical coordinate"))
        vgd_free (vgd_id)

    fields = self._headers
    # Provide 'level' and 'kind' information to the decoder.
    decoded = np.concatenate(decode_ip1(fields['ip1']))
    # Only use first set of levels (can't handle ranges yet).
    fields['level'] = decoded['level']
    fields['kind'] = decoded['kind']

  def _makevars (self):
    from fstd2nc.mixins import _var_type, _axis_type, _dim_type
    from collections import OrderedDict
    import numpy as np
    from rpnpy.vgd.base import vgd_fromlist, vgd_get, vgd_free
    from rpnpy.vgd.const import VGD_KEYS
    from rpnpy.vgd import VGDError
    from rpnpy.librmn.fstd98 import fstinl, fstprm, fstluk

    vrecs = self._vrecs

    super(VCoords,self)._makevars()

    # Scan through the data, and look for any use of vertical coordinates.
    vaxes = OrderedDict()
    prefs = dict()  # Reference pressure scalar variables.
    for var in self._varlist:
      level_axis = var.getaxis('level')
      if level_axis is None: continue
      # Degenerate vertical axis?
      if 'ip1' in var.atts and var.atts['ip1'] == 0:
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
      # Only need to provide one copy of the vertical axis.
      if (id(level_axis),kind) not in vaxes:
        # Keep track of any extra arrays needed for this axis.
        coordinates = []
        # Get metadata that's specific to this axis.
        name = 'level'
        atts = OrderedDict()
        atts['axis'] = 'Z'
        new_axis = _axis_type(name, atts, level_axis.array)
        # Reference: http://web-mrb.cmc.ec.gc.ca/science//si/eng/si/libraries/rmnlib/fstd/main.html#RTFToC11
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
        elif kind == 1:
          # sigma [sg] (0.0->1.0)
          atts['standard_name'] = 'atmosphere_sigma_coordinate'
          atts['units'] = 'sigma_level'   # units defined for compliancy with COARDS
          atts['positive'] = 'down'
          atts['formula_terms'] = OrderedDict([('sigma',new_axis),('ps','P0')])
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
          atts['units'] = 'level'  # units defined for compliancy with COARDS
          atts['positive'] = 'down'
        elif kind == 4:
          # height [M] (metres) with respect to ground level
          name = 'height'
          atts['standard_name'] = 'height'
          atts['units'] = 'm'
          atts['positive'] = 'up'
        elif kind == 5:
          # hybrid coordinates [hy] (0.0->1.0)
          atts['units'] = 'level'  # units defined for compliancy with COARDS
          atts['positive'] = 'down'
          key = (var.atts['ig1'],var.atts['ig2'])
          # If we're dealing with the old 'HY' records, then we don't match on
          # ig1/ig2.
          if key not in vrecs and 'HY' in vrecs:
            key = 'HY'
          # If there's only a single vertical record, then match that regardless of keys.
          if len(vrecs) == 1 and not self._strict_vcoord_match:
            key = list(vrecs.keys())[0]
          # Check if we have a vertical coordinate record to use.
          if key in vrecs:
            header = vrecs[key]
            # Add in metadata from the coordinate.
            atts.update(self._get_header_atts(header))
            # Add type-specific metadata.
            if header['nomvar'].strip() == '!!':
              # Get A and B info.
              vgd_id = vgd_fromlist(fstluk(header,rank=3)['d'])
              # Partial definition of a/b coordinates for formula_terms reference.
              coordA = _var_type('a', {}, [new_axis], None)
              coordB = _var_type('b', {}, [new_axis], None)

              if vgd_get (vgd_id,'LOGP'):
                # Not really a "standard" name, but there's nothing in the
                # CF convensions document on how to encode this.
                # I just merged the atmosphere_ln_pressure_coordinate and
                # atmosphere_hybrid_sigma_pressure_coordinate together.
                # http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#dimensionless-v-coord
                atts['standard_name'] = 'atmosphere_hybrid_sigma_ln_pressure_coordinate'
                # Document the formula to follow, since it's not in the conventions.
                #TODO: update this once there's an actual convention to follow!
                atts['formula'] = "p = exp(a+b*log(ps/pref))"
                atts['formula_terms'] = OrderedDict([('a',coordA),('b',coordB),('ps','P0'),('pref','pref')])
                # Try getting reference pressure as a scalar.
                try:
                  pref = vgd_get(vgd_id,'PREF')
                  if pref not in prefs:
                    prefs[pref] = _var_type('pref',{'units':'Pa'},[],np.array(pref))
                  atts['formula_terms']['pref'] = prefs[pref]
                except (KeyError,VGDError):
                  pass # Don't have PREF available for some reason?
              else:
                atts['standard_name'] = 'atmosphere_hybrid_sigma_pressure_coordinate'
                atts['formula_terms'] = OrderedDict([('ap',coordA),('b',coordB),('ps','P0')])
              # Add all parameters for this coordinate.
              internal_atts = OrderedDict()
              for key in VGD_KEYS:
                try:
                  val = vgd_get(vgd_id,key)
                  # Skip multidimensional arrays (can't be encoded as metadata).
                  if getattr(val,'ndim',1) > 1: continue
                  internal_atts[key] = val
                except (KeyError,VGDError):
                  pass  # Some keys not available in some vgrids?
              # Put this information in the final output file?
              for k,v in internal_atts.items():
                if self._rpnstd_metadata_list is None or k in self._rpnstd_metadata_list:
                  atts[k] = v
              # Attempt to fill in A/B coefficients (if available).
              try:
                all_z = list(internal_atts['VCDM'])+list(internal_atts['VCDT'])
                all_a = list(internal_atts['CA_M'])+list(internal_atts['CA_T'])
                all_b = list(internal_atts['CB_M'])+list(internal_atts['CB_T'])
                A = []
                B = []
                for z in levels:
                  ind = all_z.index(z)
                  A.append(all_a[ind])
                  B.append(all_b[ind])
                coordA.array = np.asarray(A)
                coordB.array = np.asarray(B)
                coordinates.extend([coordA,coordB])
              except (KeyError,ValueError,VGDError):
                warn (_("Unable to get A/B coefficients for %s.")%var.name)
                atts.pop('formula',None)
                atts.pop('formula_terms',None)
              vgd_free (vgd_id)
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
              A = pref * 100. * (eta - B)
              coordA = _var_type('a', {}, [new_axis], np.asarray(A))
              coordB = _var_type('b', {}, [new_axis], np.asarray(B))
              coordinates.extend([coordA,coordB])
              atts['formula_terms'] = OrderedDict([('ap',coordA),('b',coordB),('ps','P0')])
              # Add extra HY record metadata.
              atts.update(ptop=ptop, rcoef=rcoef, pref=pref)
        elif kind == 6:
          # theta [th]
          name = 'theta'
          atts['standard_name'] = 'air_potential_temperature'
          atts['units'] = 'K'
          atts['positive'] = 'up'

        # Add this vertical axis.
        if len(coordinates) > 0:
          atts['coordinates'] = coordinates
        # Update axis name.
        new_axis.name = name
        # Now have a fully defined axis to use.
        vaxes[(id(level_axis),kind)] = new_axis
      # Set the vertical axis for this variable.
      vaxis = vaxes[(id(level_axis),kind)]
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
      warn (_("Mixture of model / height levels found.  This will cause multiple definitions of variables in the output.  You may be able to resolve this by using --ignore-diag-level or --diag-as-model-level."))

