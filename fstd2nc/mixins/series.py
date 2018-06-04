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


#################################################
# Mixin for handling timeseries data.
#
# This is a first attempt at handling these 'series' files as output from
# the GEM model, and may be incomplete / inaccurate.  Please correct this
# section if you notice anything wrong.
#
# There are two types of timeseries records I've seen:
#
# - typvar='T', grtyp='Y'.
#   Here, ni corresponds to horizontal points (like the usual Y-grid).
#   There should be corresponding '^^' and '>>' fields in this case.
#
# - Vertical profiles, which have typvar='T', grtype='+'.
#   This data uses a different meaning for ni and nj.
#   Here, 'ni' is actually the # of vertical levels, and 'nj' is the number of
#   forecast times.  The horizontal points are split per record, and enumerated
#   by the 'ip3' parameter.
#   ig1/ig2 is set to zero in my sample - coincidentally matches ip1/ip2 of
#   !! record.
#   ig3/ig4 give some kind of horizontal coordinate info (?).
#   ip1/ip2 seem to match ip1/ip2 of '^^', '>>' records.
#   'HH' record gives forecast hours corresponding to nj.
#   'SH' and 'SV' give some kind of vertical info corresponding to ni, but
#   with one extra level?
#   'STNS' gives the names of the stations (corresponding to ip3 numbers?)

class Series (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    super(Series,cls)._cmdline_args(parser)
    #group = parser.add_argument_group(_('Options for profile data'))
    group = parser
    group.add_argument('--profile-momentum-vars', metavar='VAR1,VAR2,...', help=_('Comma-separated list of variables that use momentum levels.'))
    group.add_argument('--profile-thermodynamic-vars', metavar='VAR1,VAR2,...', help=_('Comma-separated list of variables that use thermodynamic levels.'))
    group.add_argument('--missing-bottom-profile-level', action='store_true', help=_('Assume the bottom level of the profile data is missing.'))
    group.add_argument('--missing-top-profile-level', action='store_true', help=_('ASsume the top level of the profile data is missing.'))

  # Need to extend _headers_dtype before __init__.
  def __new__ (cls, *args, **kwargs):
    obj = super(Series,cls).__new__(cls, *args, **kwargs)
    obj._headers_dtype = obj._headers_dtype + [('station_id','int32')]
    return obj

  def __init__ (self, *args, **kwargs):
    import numpy as np
    momentum_vars = kwargs.pop('profile_momentum_vars',None)
    if momentum_vars is None:
      momentum_vars = []
    if isinstance(momentum_vars,str):
      momentum_vars = momentum_vars.split(',')
    thermo_vars = kwargs.pop('profile_thermodynamic_vars',None)
    if thermo_vars is None:
      thermo_vars = []
    if isinstance(thermo_vars,str):
      thermo_vars = thermo_vars.split(',')
    self._momentum_vars = momentum_vars
    self._thermo_vars = thermo_vars
    self._missing_bottom_profile_level = kwargs.pop('missing_bottom_profile_level',False)
    self._missing_top_profile_level = kwargs.pop('missing_top_profile_level',False)

    # Don't process series time/station/height records as variables.
    self._meta_records = self._meta_records + ('HH','STNS','SV','SH')
    # Add station # as another axis.
    self._outer_axes = ('station_id',) + self._outer_axes
    super(Series,self).__init__(*args,**kwargs)

    fields = self._headers
    nrecs = len(fields)
    # Identify timeseries records for further processing.
    is_series = (fields['typvar'] == 'T ') & ((fields['grtyp'] == '+') | (fields['grtyp'] == 'Y') | (fields['grtyp'] == 'T'))
    # More particular, data that has one station per record.
    is_split_series = (fields['typvar'] == 'T ') & (fields['grtyp'] == '+')

    # For timeseries data, station # is provided by 'ip3'.
    station_id = np.ma.array(np.array(fields['ip3']))
    # For non-timeseries data, ignore this info.
    station_id.mask = ~is_split_series
    fields['station_id'] = station_id
    # For timeseries data, the usual leadtime (from deet*npas) is not
    # used.  Instead, we will get forecast info from nj coordinate.
    if 'leadtime' in fields.dtype.names:
      fields['leadtime'] = np.ma.asarray(fields['leadtime'])
      fields['leadtime'].mask = np.ma.getmaskarray(fields['leadtime']) | is_series
    # Similarly, the 'reftime' is not used either.
    if 'reftime' in fields.dtype.names:
      fields['reftime'] = np.ma.asarray(fields['reftime'])
      fields['reftime'].mask = np.ma.getmaskarray(fields['reftime']) | is_series

    # True grid identifier is in ip1/ip2?
    # Overwrite the original ig1,ig2,ig3,ig4 values, which aren't actually grid
    # identifiers in this case (they're just the lat/lon coordinates of each
    # station?)
    fields['ig1'][is_split_series] = fields['ip1'][is_split_series]
    fields['ig2'][is_split_series] = fields['ip2'][is_split_series]
    fields['ig3'][is_split_series] = 0
    fields['ig4'][is_split_series] = 0
    # Do not treat the ip1 value any further - it's not really vertical level.
    # Set it to 0 to indicate a degenerate vertical axis.
    fields['ip1'][is_series] = 0

  def _iter (self):
    from fstd2nc.mixins import _iter_type, _var_type, _modify_axes
    from collections import OrderedDict
    import numpy as np
    from datetime import timedelta
    from rpnpy.librmn.fstd98 import fstlir

    forecast_hours = None
    created_time_axis = False  # To only create squashed time axis once.
                               # (for --squashed-forecasts option).
    station = None             # To attach the station names as coordinates.
    momentum = thermo = None   # To attach the vertical axes.

    # Get station and forecast info.
    # Need to read from original records, because this into isn't in the
    # data stream.
    station_header = fstlir(self._meta_funit, nomvar='STNS')
    # Create forecast axis.
    forecast_header = fstlir (self._meta_funit, nomvar='HH')
    if forecast_header is not None:
        atts = OrderedDict(units='hours')
        array = forecast_header['d'].flatten()
        axes = OrderedDict(forecast=tuple(array))
        forecast = _var_type('forecast',atts,axes,array)
        forecast_hours = list(array)
        if getattr(self,'_squash_forecasts',False) is False:
          yield forecast
    # Extract vertical coordinates.
    for vertvar in ('SH','SV'):
      header = fstlir (self._meta_funit, nomvar=vertvar)
      if header is None: continue
      array = header['d'].squeeze()
      # Drop the top or bottom levels to match the profile data?
      if self._missing_bottom_profile_level:
        array = array[:-1]
      if self._missing_top_profile_level:
        array = array[1:]
      if array.ndim != 1: continue
      atts = OrderedDict(self._get_header_atts(header))
      atts['kind'] = 5
      var = _var_type(vertvar,atts,{'level':tuple(array)},array)
      if vertvar == 'SH': thermo = tuple(array)
      if vertvar == 'SV': momentum = tuple(array)
      yield var

    for var in super(Series,self)._iter():

      # Hook in the station names as coordinate information.
      if 'station_id' in var.axes and station_header is not None:
        if station is None:
          atts = OrderedDict()
          array = station_header['d'].transpose()
          # Subset the stations to match the IP3 values found in the file
          # (in case we don't have records for all the stations).
          station_id = var.axes['station_id']
          indices = np.array(var.axes['station_id'],dtype=int) - 1
          array = array[indices,:]
          # Re-cast array as string.
          # I don't know why I have to subtract 128 - maybe something to do with
          # how the characters are encoded in the file?
          # This isn't always needed.  Have test files for both cases.
          # Need help making this more robust!
          if array.flatten()[0] >= 128:
            array -= 128
          array = array.view('|S1')
          nstations, strlen = array.shape
          # Strip out trailing whitespace.
          array = array.flatten().view('|S%d'%strlen)
          array[:] = map(str.rstrip,array)
          array = array.view('|S1').reshape(nstations,strlen)
          # Encode it as 2D character array for netCDF file output.
          axes = OrderedDict([('station_id',station_id),('station_strlen',tuple(range(strlen)))])
          station = _var_type('station',atts,axes,array)
          yield station
        var.atts['coordinates'] = [station]

      if not isinstance(var,_iter_type) or var.atts.get('typvar') != 'T':
        yield var
        continue

      # 'Y' data should be handled fine by _XYCoords - just give a more
      # specific name to the ni axis for clarity.
      if var.atts.get('grtyp') == 'Y':
        var.axes = _modify_axes(var.axes, i=('station_id',tuple(range(1,len(var.axes['i'])+1))))

      # '+' data has different meanings for the axes.
      # ni is actually vertical level.
      # nj is actually forecast time.
      if var.atts.get('grtyp') == '+' and forecast_hours is not None:
        # Remove degenerate vertical axis.
        if 'level' in var.axes:
          var.record_id = var.record_id.squeeze(axis=list(var.axes.keys()).index('level'))
          var.axes.pop('level')
        # Try to map to thermodynamic or momentum levels.
        level_def = 'level'
        if var.name in self._momentum_vars and momentum is not None:
          if len(var.axes['i']) == len(momentum):
            level_def = ('level',momentum)
          else:
            warn (_("Wrong number of momentum levels found in the data."))
        if var.name in self._thermo_vars and thermo is not None:
          if len(var.axes['i']) == len(thermo):
            level_def = ('level',thermo)
          else:
            warn (_("Wrong number of thermodynamic levels found in the data."))
        if level_def == 'level':
          warn (_("Unable to find the vertical coordinates for %s."%var.name))
        var.axes = _modify_axes(var.axes, i=level_def, j=('forecast',forecast_hours))

        # Some support for squashing forecasts.
        if getattr(self,'_squash_forecasts',False) is True:
          # Can only do this for a single date of origin, because the time
          # axis and forecast axis are not adjacent for this type of data.
          if len(var.axes['time']) == 1:
            var.record_id = var.record_id.squeeze(axis=list(var.axes.keys()).index('time'))
            time = var.axes.pop('time')[0]
            # Convert pandas times (if using pandas for processing the headers)
            time = np.datetime64(time,'s')
            forecast = var.axes['forecast']
            var.axes = _modify_axes(var.axes, forecast=('time',tuple(time+np.timedelta64(int(h*3600),'s') for h in forecast)))
            if not created_time_axis:
              yield _var_type('time',OrderedDict([('standard_name','time'),('long_name','Validity time'),('axis','T')]),{'time':var.axes['time']},np.array(var.axes['time']))
              # Include forecast and reftime auxiliary coordinates (emulate
              # what's done in the dates mixin)
              leadtime = _var_type('leadtime',OrderedDict([('standard_name','forecast_period'),('long_name','Lead time (since forecast_reference_time)'),('units','hours')]),{'time':var.axes['time']},np.array(forecast))
              yield leadtime
              reftime = _var_type('reftime',OrderedDict([('standard_name','forecast_reference_time')]),{},np.array(time))
              yield reftime
              created_time_axis = True
            # Add leadtime and reftime as auxiliary coordinates.
            coords = var.atts.get('coordinates',[])
            coords.extend([leadtime,reftime])
            var.atts['coordinates'] = coords
          else:
            warn(_("Can't use datev for timeseries data with multiple dates of origin.  Try re-running with the --dateo option."))
      # Remove 'kind' information for now - still need to figure out vertical
      # coordinates (i.e. how to map SV/SH here).
      var.atts.pop('kind',None)
      yield var

