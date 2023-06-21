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


# Convert an RPN date stamp to datetime object.
# Returns None for invalid stamps.
# Scalar version, using simple datetime objects.
def stamp2datetime_scalar (date):
  from rpnpy.rpndate import RPNDate
  dummy_stamps = (0, 10101011, 101010101)
  if date not in dummy_stamps:
    return RPNDate(int(date)).toDateTime().replace(tzinfo=None)
  else:
    return None
# Vectorized version, using datetime64 objects.
from fstd2nc.mixins import vectorize
@vectorize
def stamp2datetime64 (date):
  import numpy as np
  date = stamp2datetime_scalar (date)
  if date is not None:
    return np.asarray(date, dtype='datetime64[s]')[()]
  else:
    return None
# Convert datetime64 objects back to an RPN date stamp.
@vectorize
def datetime2stamp (date):
  from rpnpy.rpndate import RPNDate
  import numpy as np
  from datetime import datetime
  if date is None: return 0
  if not isinstance(date,datetime):
    stamp = date - np.datetime64('1970-01-01T00:00:00')
    stamp /= np.timedelta64(1,'s')
    date = datetime.utcfromtimestamp(stamp)
  return RPNDate(date).datev

#################################################
# Mixin for handling dates/times.

class Dates (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    super(Dates,cls)._cmdline_args(parser)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--datev', '--squash-forecasts', action='store_true', default=True, dest='squash_forecasts', help=_('Use the date of validity for the "time" axis.  This is the default.'))
    group.add_argument('--dateo', '--forecast-axis', action='store_false', dest='squash_forecasts', help=_('Use the date of original analysis for the time axis, and put the forecast times into a separate "forecast" axis.'))

  def __init__ (self, *args, **kwargs):
    """
    forecast_axis : bool, optional
        Use the date of original analysis for the time axis, and put the
        forecast times into a separate "forecast" axis.
    """
    import numpy as np
    squash_forecasts = kwargs.pop('squash_forecasts',None)
    if squash_forecasts is None:
      squash_forecasts = kwargs.pop('datev',None)
    forecast_axis = kwargs.pop('forecast_axis',None)
    if forecast_axis is None:
      forecast_axis = kwargs.pop('dateo',None)
    if squash_forecasts is None and forecast_axis is not None:
      squash_forecasts = not forecast_axis
    if squash_forecasts is None:
      squash_forecasts = True
    self._squash_forecasts = squash_forecasts

    if self._squash_forecasts:
      self._outer_axes = ('time',) + self._outer_axes
      # Define some auxiliary coordinates along the time axis.
      self._outer_coords['leadtime'] = ('time',)
      self._outer_coords['reftime'] = ('time',)
    else:
      self._outer_axes = ('time','leadtime') + self._outer_axes
    super(Dates,self).__init__(*args,**kwargs)

    # Get any extra (derived) fields needed for doing the decoding.
    # Calculate the forecast (in hours).
    fields = self._headers
    fields['leadtime']=fields['deet']*fields['npas']/3600.
    dateo = stamp2datetime64(fields['dateo'])
    datev = stamp2datetime64(fields['datev'])
    # Convert date stamps to datetime objects, filtering out dummy values.
    dateo = np.ma.asarray(dateo, dtype='datetime64[s]')
    datev = np.ma.asarray(datev, dtype='datetime64[s]')
    dateo.mask = np.isnat(dateo)
    datev.mask = np.isnat(datev)
    # Where there are dummy dates, ignore the forecast information too.
    forecast = np.ma.asarray(fields['leadtime'], dtype='float32')
    forecast.mask = np.ma.getmaskarray(forecast) | (np.ma.getmaskarray(dateo) & np.ma.getmaskarray(datev) & (fields['deet'] == 0))
    fields['leadtime'] = forecast
    fields['reftime'] = dateo
    fields['reftime'].mask = forecast.mask
    # Time axis
    if self._squash_forecasts:
      fields['time'] = datev
    else:
      fields['time'] = dateo

  # Add time and forecast axes to the data stream.
  def _makevars (self):
    import numpy as np
    super(Dates,self)._makevars()
    for coord in self._iter_objects():
      # Add metadata to auxiliary coordinates.
      if coord.name == 'leadtime':
        coord.atts['standard_name'] = 'forecast_period'
        coord.atts['long_name'] = "Lead time (since forecast_reference_time)"
        coord.atts['units'] = 'hours'
      if coord.name == 'reftime':
        coord.atts['standard_name'] = 'forecast_reference_time'
        coord.array = np.asarray(coord.array,dtype='datetime64[s]')
        # Special case: reftimes are all identical.
        # Convert to a scalar.
        reftimes = set(coord.array)
        if len(reftimes) == 1:
          coord.axes.pop(0)
          coord.array = coord.array[0]
    for axis in self._iter_axes():
      # Add metadata to time axis.
      if axis.name == 'time':
        axis.atts['standard_name'] = 'time'
        axis.atts['long_name'] = 'Validity time'
        axis.atts['axis'] = 'T'
      # When used as an axis, rename 'leadtime' to 'forecast' for backwards
      # compatibility with previous versions of the converter.
      if axis.name == 'leadtime':
        axis.name = 'forecast'
        axis.atts['units'] = 'hours'

  def _unmakevars (self):
    import numpy as np
    # Re-attach leadtime, reftime as coordinates for variables.
    # First, find all leadtime/reftime coordinates.
    leadtimes = dict()
    reftimes = dict()
    # Helper method - find the time dimension of a variable.
    # Returned as an id.
    def get_time (var):
      time = [axis for axis in var.axes if axis.name == 'time']
      if len(time) == 0:
        return None
      else:
        return id(time[0])

    for var in self._varlist:
      time = get_time(var)
      if var.name == 'leadtime':
        leadtimes[time] = var
      if var.name == 'reftime':
        reftimes[time] = var

    # Remove these from the varlist.
    self._varlist = [var for var in self._varlist if var.name not in ('leadtime','reftime')]
    # Now, add these coordinates into the vars.
    for var in self._varlist:
      time = get_time(var)
      if time in leadtimes:
        var.atts.setdefault('coordinates',[]).append(leadtimes[time])
      if time in reftimes:
        var.atts.setdefault('coordinates',[]).append(reftimes[time])

    # Continue processing the variables.
    super (Dates,self)._unmakevars()

    # Get leadtime column (may be coming from forecast axis).
    if 'forecast' in self._headers.keys():
      self._headers['leadtime'] = self._headers['forecast']
    # Convert leadtime to units of hours if it's a timedelta64.
    if 'leadtime' in self._headers.keys():
      if self._headers['leadtime'].dtype == 'timedelta64[ns]':
        self._headers['leadtime'] = self._headers['leadtime'] / np.timedelta64(3600,'s')

    # Get datev, dateo, npas, using time, leadtime, deet.
    if 'deet' not in self._headers.keys():
      self._headers['deet'] = np.empty(self._nrecs, dtype='int32')
      self._headers['deet'][:] = 60  # Set default value of 60s
    self._headers['npas'] = np.ma.array(self._headers['leadtime']*3600 / self._headers['deet'], dtype='int32')
    if 'time' in self._headers.keys():
      if 'forecast' in self._headers.keys():
        self._headers['dateo'] = self._headers['time']
        self._headers['datev'] = self._headers['time'] + self._headers['npas'] * self._headers['deet'] * np.timedelta64(1,'s')
      else:
        self._headers['dateo'] = self._headers['time']
        self._headers['datev'] = self._headers['time']

    # Convert dateo and datev into RPN date stamps.
    if 'datev' in self._headers.keys():
      self._headers['datev'] = datetime2stamp(self._headers['datev'])
    if 'dateo' in self._headers.keys():
      self._headers['dateo'] = datetime2stamp(self._headers['dateo'])

    # Set ip2 values.
    self._headers['ip2'] = self._headers['npas'] * self._headers['deet'] // 3600
