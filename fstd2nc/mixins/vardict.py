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
from fstd2nc.mixins.vcoords import decode_ip1

#################################################
# Mixin for utilizing variable information from CMC-RPN variable dictionary.
#
class VarDict (BufferBase):

  @classmethod
  def _cmdline_args (cls, parser):
    import argparse
    super(VarDict,cls)._cmdline_args(parser)
    parser.add_argument('--vardict', action='append', help=_('Use metadata from the specified variable dictionary (XML format).'))
    parser.add_argument('--opdict', action='store_true', help=_('Similar to above, but use the standard CMC-RPN operational dictionary.'))

  @classmethod
  def _check_args (cls, parser, args):
    from os import environ
    from os.path import exists

    super(VarDict,cls)._check_args(parser,args)

    if args.opdict:
      if 'CMCCONST' in environ:
        f = environ['CMCCONST']+'/opdict/ops.variable_dictionary.xml'
        if not exists(f):
          parser.error(_("Unable to find $CMCCONST/opdict/ops.variable_dictionary.xml"))
      elif 'AFSISIO' in environ:
        f = environ['AFSISIO']+'/datafiles/constants/opdict/ops.variable_dictionary.xml'
        if not exists(f):
          parser.error(_("Unable to find $AFSISIO/datafiles/constants/opdict/ops.variable_dictionary.xml"))
      else:
        parser.error(_("Neither $CMCCONST nor $AFSISIO defined.  Can't find operational dictionary."))

    if args.vardict is not None:
      for f in args.vardict:
        if not exists(f):
          parser.error(_("Unable to find '%s'")%f)

  def __init__ (self, *args, **kwargs):
    """
    vardict : str or list, optional
        Use metadata from the specified CMC-RPN variable dictionary (XML format).
    opdict : bool, optional
        Similar to above, but use the standard CMC-RPN operational dictionary.
    """
    from os import environ
    from fstd2nc.mixins import _var_type, _axis_type, _dim_type
    from collections import OrderedDict
    from xml.etree import ElementTree as ET

    vardicts = kwargs.pop('vardict',None) or []
    if kwargs.pop('opdict',False):
      if 'CMCCONST' in environ:
        f = environ['CMCCONST']+'/opdict/ops.variable_dictionary.xml'
      elif 'AFSISIO' in environ:
        f = environ['AFSISIO']+'/datafiles/constants/opdict/ops.variable_dictionary.xml'
      vardicts.append(f)

    super(VarDict,self).__init__(*args, **kwargs)

    metadata = OrderedDict()
    ip1_axis = OrderedDict()
    ip3_axis = OrderedDict()

    for vardict in vardicts:
      try:
        metvars = list(ET.parse(vardict).getroot())
      except ET.ParseError:
        error (_("Invalid dictionary file '%s'"%vardict.name))

      for metvar in metvars:
        if metvar.attrib.get('usage','current') != 'current': continue
        nomvar = metvar.findtext('nomvar')
        if nomvar is None: continue
        var = metadata.setdefault(nomvar,OrderedDict())
        atts = metvar.find('nomvar').attrib
        coords = (atts.get('ip1',None), atts.get('ip3',None))
        d = var.setdefault(coords,OrderedDict())
        for desc in metvar.iterfind("description/short"):
          if desc.attrib.get('lang','') == 'en':
            d['long_name'] = desc.text
        for desc in metvar.iterfind("description/long"):
          if desc.attrib.get('lang','') == 'en':
            if desc.text is not None:
              d['definition_opdict'] = desc.text
        units = metvar.find("measure/real/units")
        if units is not None:
          d['units'] = units.text

      for nomvar in metadata.keys():
        descs = [d['long_name'] for d in metadata[nomvar].values() if 'long_name' in d]
        n = max(len(desc) for desc in descs)
        while len(set(desc[:n] for desc in descs)) > 1: n = n - 1
        long_name = descs[0][:n].rstrip().rstrip('(').rstrip()
        descs = [d['definition_opdict'] for d in metadata[nomvar].values() if 'definition_opdict' in d]
        if len(descs) > 0:
          n = max(len(desc) for desc in descs)
          while len(set(desc[:n] for desc in descs)) > 1: n = n - 1
          definition_opdict = descs[0][:n].rstrip().rstrip('(').rstrip()
        else:
          definition_opdict = None
        for (ip1,ip3), d in metadata[nomvar].items():
          if ip1 is not None:
            ip1 = int(ip1)
            ip1_axis.setdefault(nomvar,OrderedDict())
            if ip1 not in ip1_axis[nomvar]:
              ip1_axis[nomvar][ip1] = d['long_name'][len(long_name):].strip().lstrip('(').rstrip(')')
          if ip3 is not None:
            ip3 = int(ip3)
            ip3_axis.setdefault(nomvar,OrderedDict())
            if ip3 not in ip3_axis[nomvar]:
              ip3_axis[nomvar][ip3] = d['long_name'][len(long_name):].strip().lstrip('(').rstrip(')')
        metadata[nomvar] = list(metadata[nomvar].values())[0]
        metadata[nomvar]['long_name'] = long_name
        if definition_opdict is not None:
          metadata[nomvar]['definition_opdict'] = definition_opdict

    self._vardict = metadata
    self._vardict_ip1_axis = ip1_axis
    self._vardict_ip3_axis = ip3_axis

  def _makevars (self):
    from fstd2nc.mixins import _var_type, _axis_type, _dim_type
    from collections import OrderedDict
    import numpy as np

    handled_agg_codes = dict()

    super(VarDict,self)._makevars()

    for var in self._varlist:

      # Apply metadata from the dictionary files.
      if var.name in self._vardict:
        var.atts.update(self._vardict[var.name])

      levels = var.getaxis('level')

      # Look for variables with ip1 codes.
      if var.atts.get('kind',None) != 3 or levels is None:
        continue

      codes = tuple(levels.array)
      coordinates = []

      if var.name in self._vardict_ip1_axis:
        # Generate the list of surface types.
        lookup = OrderedDict((decode_ip1(ip1)['level'][0], name) for ip1,name in self._vardict_ip1_axis[var.name].items())
        codenames = tuple(lookup.get(code,"unknown") for code in codes)
        if codenames not in handled_agg_codes:
          array = np.array(codenames,dtype=np.char.string_).view('|S1').reshape(len(codes),-1)
          sfctype = _dim_type('sfctype',array.shape[0])
          sfctype_strlen = _dim_type('sfctype_strlen',array.shape[1])
          surface_type = _var_type("surface_type",{},[sfctype,sfctype_strlen],array)
          handled_agg_codes[codenames] = surface_type
        surface_type = handled_agg_codes[codenames]
        # "Levels" are actually surface type ids.
        var.axes[var.dims.index('level')] = surface_type.getaxis('sfctype')
        # Add the area type to the list of auxiliary coordinates.
        coordinates.append(surface_type)

      if len(coordinates) > 0:
        var.deps.extend(coordinates)


