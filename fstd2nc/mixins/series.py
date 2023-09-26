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
#   'HH' record gives forecast hours corresponding to nj.
#   'SH' and 'SV' give some kind of vertical info corresponding to ni, but
#   with one extra level sometimes (due to a bug in rpnphy).
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

  # Don't process series time/station/height records as variables.
  @classmethod
  def _meta_records(cls):
    return super(Series,cls)._meta_records() + (b'STNS',)
  @classmethod
  def _maybe_meta_records(cls):
    return super(Series,cls)._maybe_meta_records() + (b'HH',b'SV',b'SH')

  def __init__ (self, *args, **kwargs):
    """
    profile_momentum_vars : str or list, optional
        List of variables that use momentum levels.
    profile_thermodynamic_vars : str or list, optional
        List of variables that use thermodynamic levels.
    missing_bottom_profile_level : bool, optional
        Assume the bottom level of the profile data is missing.
    """
    import numpy as np
    momentum_vars = kwargs.pop('profile_momentum_vars',None)
    if momentum_vars is None:
      momentum_vars = []
    if isinstance(momentum_vars,str):
      momentum_vars = momentum_vars.replace(',',' ')
      momentum_vars = momentum_vars.split()
    thermo_vars = kwargs.pop('profile_thermodynamic_vars',None)
    if thermo_vars is None:
      thermo_vars = []
    if isinstance(thermo_vars,str):
      thermo_vars = thermo_vars.replace(',',' ')
      thermo_vars = thermo_vars.split()
    self._momentum_vars = momentum_vars
    self._thermo_vars = thermo_vars
    self._missing_bottom_profile_level = kwargs.pop('missing_bottom_profile_level',False)

    super(Series,self).__init__(*args,**kwargs)

    fields = self._headers
    nrecs = self._nrecs
    # Identify timeseries records for further processing.
    is_series = (fields['typvar'] == b'T ') & ((fields['grtyp'] == b'+') | (fields['grtyp'] == b'Y') | (fields['grtyp'] == b'T'))
    if not np.any(is_series):
      return
    # More particular, data that has one station per record.
    is_split_series = (fields['typvar'] == b'T ') & (fields['grtyp'] == b'+')

    # Add station # as another axis.
    self._outer_axes = ('station_id',) + self._outer_axes

    # For timeseries data, station # is provided by 'ip3'.
    station_id = np.ma.array(np.array(fields['ip3']), dtype='int32')
    # For non-timeseries data, ignore this info.
    station_id.mask = ~is_split_series
    fields['station_id'] = station_id
    # For timeseries data, the usual leadtime (from deet*npas) is not
    # used.  Instead, we will get forecast info from nj coordinate.
    if 'leadtime' in fields:
      fields['leadtime'] = np.ma.asarray(fields['leadtime'])
      fields['leadtime'].mask = np.ma.getmaskarray(fields['leadtime']) | is_series
    # Similarly, the 'reftime' is not used either.
    if 'reftime' in fields:
      fields['reftime'] = np.ma.asarray(fields['reftime'])
      fields['reftime'].mask = np.ma.getmaskarray(fields['reftime']) | is_series

    # Overwrite the original ig1,ig2,ig3,ig4 values, which aren't actually grid
    # identifiers in this case (they're just the lat/lon coordinates of each
    # station?)
    fields['ig1'][is_series] = 0
    fields['ig2'][is_series] = 0
    fields['ig3'][is_series] = 0
    fields['ig4'][is_series] = 0
    # Do not treat the ip1 value any further - it's not really vertical level.
    # Set it to 0 to indicate a degenerate vertical axis.
    fields['ip1'][is_series] = 0

  def _makevars (self):
    from fstd2nc.mixins import _var_type, _axis_type, _dim_type
    from fstd2nc.mixins.dates import stamp2datetime_scalar
    from collections import OrderedDict
    import numpy as np

    forecast_axis = None       # To attach the forecast axis.
    station = None             # To attach the station names as coordinates.
    momentum = thermo = None   # To attach the vertical axes.

    super(Series,self)._makevars()
    if 'station_id' not in self._headers: return

    # Get station and forecast info.
    # Need to read from original records, because this into isn't in the
    # data stream.
    station_header = self._fstlir(nomvar=b'STNS', grtyp=b'T')
    if station_header is not None:
      array = station_header['d'].transpose()
      # Re-cast array as string.
      # I don't know why I have to subtract 128 - maybe something to do with
      # how the characters are encoded in the file?
      # This isn't always needed.  Have test files for both cases.
      # Need help making this more robust!
      if array.flatten()[0] >= 128:
        array -= 128
      array = array.view('|S1')
      nstations, strlen = array.shape
      array = array.flatten().view('|S%d'%strlen)
      # Strip out trailing whitespace.
      # Python3: convert bytes to str
      array[:] = [str(arr.decode()).rstrip() for arr in array]
      array = array.view('|S1').reshape(nstations,strlen)
      station_id = _dim_type('station_id',nstations)
      station_strlen = _dim_type('station_strlen',strlen)
      # Encode it as 2D character array for netCDF file output.
      station = _var_type('station',{},[station_id,station_strlen],array)
    # Create forecast axis.
    forecast_header = self._fstlir (nomvar=b'HH  ', grtyp=b'T')
    if forecast_header is not None:
      atts = OrderedDict(units='hours')
      # Note: the information in 'HH' is actually the hour of validity.
      # Need to subtract the hour from the date of origin in order to get
      # the leadtime.
      starting_hour = stamp2datetime_scalar(forecast_header['datev']).hour - forecast_header['npas']*forecast_header['deet']
      array = forecast_header['d'].flatten() - starting_hour
      forecast_timedelta = np.array(array*3600,'timedelta64[s]')
      forecast_axis = _axis_type('forecast',atts,array)
    # Extract vertical coordinates.
    for vertvar in (b'SH  ',b'SV  '):
      header = self._fstlir (nomvar=vertvar, grtyp=b'T')
      if header is None: continue
      array = header['d'].squeeze()
      # Drop the top or bottom levels to match the profile data?
      if self._missing_bottom_profile_level:
        array = array[:-1]
      if array.ndim != 1: continue
      atts = OrderedDict(self._get_header_atts(header))
      if vertvar == b'SH  ': thermo = _axis_type('level',atts,array)
      if vertvar == b'SV  ': momentum = _axis_type('level',atts,array)


    # 'Y' data should be handled fine by _XYCoords - just give a more
    # specific name to the ni axis for clarity.
    for var in self._varlist:
      if var.atts.get('typvar') == 'T' and var.atts.get('grtyp') == 'Y':
        dims = var.dims
        iaxis = var.getaxis('i')
        if iaxis is not None and station is not None and len(iaxis) == station.shape[0]:
          var.axes[dims.index('i')] = station.axes[0]

    # Remove degenerate vertical axis for '+' data.
    # (The one coming from IP1, which is not used.)
    for var in self._varlist:
      if var.atts.get('typvar') == 'T' and var.atts.get('grtyp') == '+':
        dims = var.dims
        if 'level' in dims:
          var.record_id = var.record_id.squeeze(axis=dims.index('level'))
          var.axes.pop(dims.index('level'))


    # For '+' data, ni is the vertical level, and nj is the forecast.
    known_levels = dict()
    for var in self._varlist:

      if var.atts.get('typvar') != 'T': continue
      if var.atts.get('grtyp') != '+': continue

      dims = var.dims

      # The j dimension is actually the forecast time.
      jaxis = var.getaxis('j')
      if jaxis is not None and forecast_axis is not None and len(jaxis) == len(forecast_axis):
        var.axes[dims.index('j')] = forecast_axis

      # The i dimension is actually the vertical coordinate for this type of
      # data.
      iaxis = var.getaxis('i')
      if iaxis is not None:
        # If there's only 1 level (degenerate), then remove that dimension.
        if len(iaxis) == 1:
          var.axes.pop(dims.index('i'))
          continue
        # Try to map to thermodynamic or momentum levels.
        level = iaxis
        level.name = 'level'
        if var.name in self._momentum_vars and momentum is not None:
          if len(level) == len(momentum):
            level = momentum
          else:
            warn (_("Wrong number of momentum levels found in the data."))
        if var.name in self._thermo_vars and thermo is not None:
          if len(level) == len(thermo):
            level = thermo
          else:
            warn (_("Wrong number of thermodynamic levels found in the data."))
        if level is iaxis:
          warn (_("Unable to find the vertical coordinates for %s."%var.name))
          # Attach a generic level dimension.
          nlev = len(level)
          if nlev not in known_levels:
            known_levels[nlev] = _dim_type('level',nlev)
          level = known_levels[nlev]
        else:
          # Found vertical levels, now define the level kind so VCoords
          # mixin can add more metadata.
          var.atts['kind'] = 5
        var.axes[dims.index('i')] = level

    # Some support for squashing forecasts.
    if getattr(self,'_squash_forecasts',False) is True:
      known_squashed_forecasts = dict()
      known_leadtimes = dict()
      known_reftimes = dict()
      for var in self._varlist:
        # Can only do this for a single date of origin, because the time
        # axis and forecast axis are not adjacent for this type of data.
        time = var.getaxis('time')
        forecast = var.getaxis('forecast')
        if time is None or forecast is None: continue
        if len(time) != 1:
          warn(_("Can't use datev for timeseries data with multiple dates of origin.  Try re-running with the --dateo option."))
          continue
        var.record_id = var.record_id.squeeze(axis=var.dims.index('time'))
        var.axes.pop(var.dims.index('time'))
        key = (id(time),id(forecast))
        if key not in known_squashed_forecasts:
          time0 = time.array[0]
          # Convert pandas times (if using pandas for processing the headers)
          time0 = np.datetime64(time0,'s')
          # Calculate the date of validity
          forecast_timedelta = np.array(forecast.array*3600,'timedelta64[s]')
          squashed_times_array = time0+forecast_timedelta
          time = _axis_type('time',OrderedDict([('standard_name','time'),('long_name','Validity time'),('axis','T')]),squashed_times_array)
          known_squashed_forecasts[key] = time
          # Include forecast and reftime auxiliary coordinates (emulate
          # what's done in the dates mixin)
          leadtime = _var_type('leadtime',OrderedDict([('standard_name','forecast_period'),('long_name','Lead time (since forecast_reference_time)'),('units','hours')]),[time],forecast.array)
          reftime = _var_type('reftime',OrderedDict([('standard_name','forecast_reference_time')]),{},np.array(time0))
          known_leadtimes[key] = leadtime
          known_reftimes[key] = reftime
        var.axes[var.dims.index('forecast')] = known_squashed_forecasts[key]
        # Add leadtime and reftime as auxiliary coordinates.
        var.deps.extend([known_leadtimes[key],known_reftimes[key]])

    # Hook in the station names as coordinate information.
    if station is not None:
      for station_id, varlist in self._iter_axes('station_id', varlist=True):
        # Try to use the provided station coordinate, if it has a consistent
        # length.
        if len(station_id) == station.shape[0]:
          for var in varlist:
            var.axes[var.dims.index('station_id')] = station.axes[0]
          station_id = station.axes[0]
          station_coord = station
        # Otherwise, need to construct a new coordinate with the subset of
        # stations used.
        # Assume station_ids start at 1 (not 0).
        else:
          indices = station_id.array - 1
          array = station.array[indices,:]
          # Use _axis_type instead of _dim_type to retain the station_id values.
          station_id = _axis_type('station_id',{},station_id.array)
          axes = [station_id,station.axes[1]]
          station_coord = _var_type('station',{},axes,array)
        # Attach the station as a coordinate.
        for var in varlist:
          var.deps.append(station_coord)

  # Decode station data back into table.
  def _unmakevars (self):
    from fstd2nc.mixins import _iter_type, _chunk_type, _var_type, _axis_type, _dim_type
    import numpy as np

    # If no series data defined in the dataset, then nothing to do here.
    if not any(var.name == 'station' for var in self._varlist):
      return super(Series,self)._unmakevars()

    # Helper - add a coordinate record
    # Copied from xycoords mixin.
    def add_coord (name,nj,ni,values,**atts):
      import dask.array as da
      dims = []
      if nj > 1: dims.append(_dim_type('j',nj))
      if ni > 1: dims.append(_dim_type('i',ni))
      # Need to carefully encode the arrays so they are scalar object containing the data.
      array = np.empty(1,dtype=object)
      array[0] = da.from_array(values, chunks=-1)
      dtype = array[0].dtype
      array = array.squeeze()
      self._varlist.append(_iter_type(name,atts,dims,dtype,array))

    # Encode stations.
    for var in self._varlist:
      if var.name == 'station' and len(var.dims) == 1 and hasattr(var,'array'):
        array = var.array
        nstations = len(array)
        array = array.view('|S1').reshape(nstations,-1)
        # Pad out to 128 characters?  (Is this needed?)
        station_len = array.shape[-1]
        if station_len < 128:
          new_array = np.zeros((nstations,128),dtype='|S1')
          new_array[:,:station_len] = array
          array = new_array
        station_len = array.shape[-1]
        # Pad with spaces.
        array[array==b''] = b' '
        add_coord('STNS',nstations,station_len,array,typvar='T',datyp=7,nbits=8,grtyp='T')
        break
    self._varlist = [var for var in self._varlist if var.name != 'station']

    # Encode lat/lon axes.
    for var in self._varlist:
      if 'station_id' not in var.dims: continue
      nstations = len(var.getaxis('station_id'))
      for coord in var.atts.get('coordinates',[]):
        if coord.name == 'lat':
          add_coord('^^',1,nstations,coord.array,typvar='T',datyp=5,nbits=32,grtyp='T')
        if coord.name == 'lon':
          add_coord('>>',1,nstations,coord.array,typvar='T',datyp=5,nbits=32,grtyp='T')
      break

    # Encode vertical coordinates.
    # Note: can only encode one set of vertical coordinates in a series
    # dataset.
    for var in self._varlist:
      if 'station_id' not in var.dims: continue
      if 'level' in var.dims:
        levels = var.getaxis('level')
        if hasattr(levels,'atts'):
          if 'VCDT' in levels.atts and 'VCDM' in levels.atts:
            # Add thermodynamic coordinate variable.
            array = levels.atts['VCDT'][:-2]  # Skip 1.0 level and diag level?
            add_coord('SH',1,len(array),array,typvar='T',datyp=5,nbits=32,grtyp='T')
            # Add momentum coordinate variable.
            array = levels.atts['VCDM'][:-1]  # Skip 1.0 level
            add_coord('SV',1,len(array),array,typvar='T',datyp=5,nbits=32,grtyp='T')
            break

    # If the data does not contain a forecast axis, then reconstruct one
    # based on leadtime.
    # Only do this if more than one forecast time is available.
    #TODO: Move some of this logic to dates mixin?
    forecast_axes = {}
    for var in self._varlist:
      if 'time' not in var.dims: continue
      time = var.getaxis('time')
      if var.name == 'leadtime' and hasattr(var,'array'):
        forecasts = np.unique(var.array)
        # Convert forecast to timedelta64 if it's in hours.
        if forecasts.dtype != 'timedelta64[ns]':
          forecasts = forecasts * np.timedelta64(3600,'s')
        forecast_axes[id(time)] = _axis_type('forecast',var.atts,forecasts)
        break
    time0 = None
    for var in self._varlist:
      if var.name == 'reftime' and hasattr(var,'array'):
        time0 = np.atleast_1d(var.array)
        time0 = _axis_type('time',{},time0)
        break
    for var in self._varlist:
      if 'station_id' not in var.dims: continue
      if 'time' not in var.dims: continue
      # Only do this for '+' grid, where time axis is fastest varying.
      itime = var.dims.index('time')
      if itime < var.dims.index('station_id'): continue
      time = var.getaxis('time')
      if id(time) not in forecast_axes: continue
      if time0 is None: continue
      forecast = forecast_axes[id(time)]
      if len(time0) * len(forecast) != len(time): continue
      var.axes = var.axes[:itime] + [time0,forecast] + var.axes[itime+1:]
      var.array = var.array.reshape(var.shape)

    # Encode forecast info.
    # Note: similar to vertical info, can only encode one forecast variable
    # here.
    for var in self._varlist:
      if 'station_id' not in var.dims: continue
      if 'forecast' in var.dims:
        forecasts = var.getaxis('forecast')
        if hasattr(forecasts,'array'):
          forecasts = forecasts.array
          # Convert forecast to units of hours if it's a timedelta64.
          if forecasts.dtype == 'timedelta64[ns]':
            forecasts = forecasts / np.timedelta64(3600,'s')
          add_coord('HH',1,len(forecasts),forecasts,typvar='T',datyp=5,nbits=32,grtyp='T')
          break

    # Check if dimensions need to be transposed.
    # Last dims should be station_id or station_id,forecast.
    for var in self._varlist:
      if 'station_id' not in var.dims: continue
      if not hasattr(var,'array'): continue
      dims = var.dims
      order = [i for i,d in enumerate(dims) if d not in ('level','station_id','forecast')]
      order.append(dims.index('station_id'))
      if 'level' in dims:
        order.append(dims.index('level'))
      if 'forecast' in dims:
        order.append(dims.index('forecast'))
      if order != sorted(order):
        var.axes = [var.axes[o] for o in order]
        var.array = var.array.transpose(order)

    # Final preparation of station-based variables.
    for var_ind, var in enumerate(self._varlist):
      axes = list(var.axes)
      if 'station_id' not in var.dims: continue
      # Process the variable.
      # Case 1: Have a forecast axis.
      # Each station is in a separate record, which contains all forecasts.
      if 'forecast' in var.dims:
        if 'level' in var.dims:
          ni = len(var.getaxis('level'))
          axes[var.dims.index('level')] = _dim_type('i',ni)
        nj = len(var.getaxis('forecast'))
        axes[var.dims.index('forecast')] = _dim_type('j',nj)
        # Enumerate the stations (starting at 1) and store in ip3.
        axis = var.getaxis('station_id')
        axis_ind = var.dims.index('station_id')
        newaxis = _axis_type('ip3',{},np.arange(len(axis))+1)
        axes[axis_ind] = newaxis
        num_outer_axes = axis_ind + 1
        var.atts['grtyp'] = '+'
        # Patch up lat/lon coordinates.
        for coord_ind, coord in enumerate(var.atts.get('coordinates',[])):
          if coord.name == 'lat':
            encoded_lat = _var_type(coord.name, coord.atts, [newaxis], coord.array)
            var.atts['coordinates'][coord_ind] = encoded_lat
          if coord.name == 'lon':
            encoded_lon = _var_type(coord.name, coord.atts, [newaxis], coord.array)
            var.atts['coordinates'][coord_ind] = encoded_lon
      else:
        ni = len(var.getaxis('station_id'))
        axes[var.dims.index('station_id')] = _dim_type('i',ni)
        num_outer_axes = var.dims.index('station_id')
        var.atts['grtyp'] = 'Y'
      outer_shape = var.shape[:num_outer_axes]
      record_id = np.zeros(outer_shape, dtype=object)
      for ind in np.ndindex(tuple(outer_shape)):
        record_id[ind] = var.array[ind]
      var.atts['typvar'] = 'T'
      var = _iter_type (var.name, var.atts, axes, var.array.dtype, record_id)
      self._varlist[var_ind] = var

    super(Series,self)._unmakevars()

    # Encode lat/lon in ig3/ig4 values.
    encode_latlon = (self._headers['grtyp'] == '+')
    if 'ig3' not in self._headers.keys():
      self._headers['ig3'] = np.ma.masked_all(self._nrecs,dtype='int32')
    self._headers['ig3'][encode_latlon] = np.round((self._headers['lat']+100)*100)[encode_latlon]
    if 'ig4' not in self._headers.keys():
      self._headers['ig4'] = np.ma.masked_all(self._nrecs,dtype='int32')
    self._headers['ig4'][encode_latlon] = np.round(self._headers['lon']*100)[encode_latlon]

