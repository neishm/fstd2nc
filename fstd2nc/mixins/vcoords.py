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
  out['level'] = r1.v1
  return out



#################################################
# Mixin for handling vertical coordinates.

class VCoords (BufferBase):
  _vcoord_nomvars = ('HY','!!')

  # Need to extend _headers_dtype before __init__.
  def __new__ (cls, *args, **kwargs):
    obj = super(VCoords,cls).__new__(cls, *args, **kwargs)
    obj._headers_dtype = obj._headers_dtype + [('level','float32'),('kind','int32')]
    return obj

  def __init__ (self, *args, **kwargs):
    import numpy as np
    # Use decoded IP1 values as the vertical axis.
    self._outer_axes = ('level',) + self._outer_axes
    # Tell the decoder not to process vertical records as variables.
    self._meta_records = self._meta_records + self._vcoord_nomvars
    super(VCoords,self).__init__(*args,**kwargs)
    # Don't group records across different level 'kind'.
    # (otherwise can't create a coherent vertical axis).
    self._var_id = self._var_id + ('kind',)
    self._human_var_id = self._human_var_id + ('vgrid%(kind)s',)
    fields = self._headers
    # Provide 'level' and 'kind' information to the decoder.
    decoded = np.concatenate(decode_ip1(fields['ip1']))
    # Only use first set of levels (can't handle ranges yet).
    fields['level'] = decoded['level']
    fields['kind'] = decoded['kind']
  # Add vertical axis as another variable.
  def _iter (self):
    from fstd2nc.mixins import _var_type, _modify_axes
    from collections import OrderedDict
    import numpy as np
    from rpnpy.vgd.base import vgd_fromlist, vgd_get, vgd_free
    from rpnpy.vgd.const import VGD_KEYS
    from rpnpy.vgd import VGDError
    from rpnpy.librmn.fstd98 import fstinl, fstprm, fstluk
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

    # Scan through the data, and look for any use of vertical coordinates.
    vaxes = OrderedDict()
    for var in super(VCoords,self)._iter():
      # Degenerate vertical axis?
      if 'ip1' in var.atts and var.atts['ip1'] == 0:
        if 'level' in var.axes and len(var.axes['level']) == 1:
          i = list(var.axes).index('level')
          del var.axes['level']
          # Remove the axis if it's an "outer" axis (not along ni,nj,nk).
          if i < var.record_id.ndim:
            var.record_id = var.record_id.squeeze(axis=i)
          yield var
          continue
      # No vertical axis?
      if 'level' not in var.axes or 'kind' not in var.atts:
        yield var
        continue
      # Decode the vertical coordinate.
      levels = var.axes['level']
      kind = var.atts['kind']
      # Only need to provide one copy of the vertical axis.
      if (levels,kind) not in vaxes:
        # Keep track of any extra arrays needed for this axis.
        ancillary_variables = []
        coordinates = []
        # Get metadata that's specific to this axis.
        name = 'zaxis'
        atts = OrderedDict()
        atts['axis'] = 'Z'
        # Reference: http://web-mrb.cmc.ec.gc.ca/science//si/eng/si/libraries/rmnlib/fstd/main.html#RTFToC11
        if kind == 0:
          # height [m] (metres)
          name = 'height'
          atts['standard_name'] = 'height'
          atts['units'] = 'm'
          atts['positive'] = 'up'
        elif kind == 1:
          # sigma [sg] (0.0->1.0)
          name = 'sigma'
          atts['standard_name'] = 'atmosphere_sigma_coordinate'
          atts['units'] = 'sigma_level'   # units defined for compliancy with COARDS
          atts['positive'] = 'down'
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
          # Check if we have a vertical coordinate record to use.
          if key in vrecs:
            header = vrecs[key]
            # Add in metadata from the coordinate.
            atts.update(self._get_header_atts(header))
            # Add type-specific metadata.
            if header['nomvar'].strip() == '!!':
              # Get A and B info.
              vgd_id = vgd_fromlist(fstluk(header,rank=3)['d'])
              if vgd_get (vgd_id,'LOGP'):
                name = 'zeta'
                # Not really a "standard" name, but there's nothing in the
                # CF convensions document on how to encode this.
                # I just merged the atmosphere_ln_pressure_coordinate and
                # atmosphere_hybrid_sigma_pressure_coordinate together.
                # http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#dimensionless-v-coord
                atts['standard_name'] = 'atmosphere_hybrid_sigma_ln_pressure_coordinate'
              else:
                name = 'eta'
                atts['standard_name'] = 'atmosphere_hybrid_sigma_pressure_coordinate'
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
              if not self._minimal_metadata:
                atts.update(internal_atts)
              # Attempt to fill in A/B ancillary data (if available).
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
                ancA = _var_type(name+'_A', {}, {name:levels}, np.asarray(A))
                ancB = _var_type(name+'_B', {}, {name:levels}, np.asarray(B))
                ancillary_variables.extend([ancA,ancB])
                coordA = _var_type('a', {}, {name:levels}, np.asarray(A))
                coordB = _var_type('b', {}, {name:levels}, np.asarray(B))
                coordinates.extend([coordA,coordB])
              except (KeyError,ValueError,VGDError):
                warn (_("Unable to get A/B coefficients."))
              vgd_free (vgd_id)
            # Not a '!!' coordinate, so must be 'HY'?
            else:
              name = 'eta'
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
              ancB = _var_type(name+'_B', {}, {name:levels}, B)
              ancA = _var_type(name+'_A', {}, {name:levels}, A)
              ancillary_variables.extend([ancA,ancB])
              coordA = _var_type('a', {}, {name:levels}, np.asarray(A))
              coordB = _var_type('b', {}, {name:levels}, np.asarray(B))
              coordinates.extend([coordA,coordB])
              # Add extra HY record metadata.
              atts.update(ptop=ptop, rcoef=rcoef, pref=pref)
        elif kind == 6:
          # theta [th]
          name = 'theta'
          atts['standard_name'] = 'air_potential_temperature'
          atts['units'] = 'K'
          atts['positive'] = 'up'

        # Add this vertical axis.
        axes = OrderedDict([(name,levels)])
        if len(ancillary_variables) > 0:
          atts['ancillary_variables'] = ' '.join(v.name for v in ancillary_variables)
          atts['coordinates'] = ' '.join(v.name for v in coordinates)
        array = np.asarray(levels)
        vaxes[(levels,kind)] = _var_type(name,atts,axes,array)
        yield vaxes[(levels,kind)]
        # Add any ancillary data needed for the axis.
        for anc in ancillary_variables:
          yield anc
        # Add coordinates.
        for coord in coordinates:
          yield coord
      # Get the vertical axis.
      vaxis = vaxes[(levels,kind)]
      # Modify the variable's dimension name to match the axis name.
      var.axes = _modify_axes(var.axes, level=vaxis.name)
      yield var

