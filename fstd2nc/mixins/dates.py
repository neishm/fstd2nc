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


# Convert an RPN date stamp to datetime object.
# Returns None for invalid stamps.
from fstd2nc.mixins import vectorize
@vectorize
def stamp2datetime (date):
  from rpnpy.rpndate import RPNDate
  dummy_stamps = (0, 10101011)
  if date not in dummy_stamps:
    return RPNDate(int(date)).toDateTime().replace(tzinfo=None)
  else:
    return None


#################################################
# Mixin for handling dates/times.

class Dates (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    super(Dates,cls)._cmdline_args(parser)
    parser.add_argument('--datev', '--squash-forecasts', action='store_true', dest='squash_forecasts', help=_('Use the date of validity for the "time" axis.  Otherwise, the default is to use the date of original analysis, and the forecast length goes in a "forecast" axis.'))

  # Need to extend _headers_dtype before __init__.
  def __new__ (cls, *args, **kwargs):
    obj = super(Dates,cls).__new__(cls, *args, **kwargs)
    obj._headers_dtype = obj._headers_dtype + [('time','datetime64[s]'),('forecast','float32')]
    return obj

  def __init__ (self, *args, **kwargs):
    import numpy as np
    self._squash_forecasts = kwargs.pop('squash_forecasts',False) or kwargs.pop('datev',False)
    if self._squash_forecasts:
      self._outer_axes = ('time',) + self._outer_axes
    else:
      self._outer_axes = ('time','forecast') + self._outer_axes
    super(Dates,self).__init__(*args,**kwargs)

    # Get any extra (derived) fields needed for doing the decoding.
    # Calculate the forecast (in hours).
    fields = self._headers
    fields['forecast']=fields['deet']*fields['npas']/3600.
    # Time axis
    if self._squash_forecasts:
      dates = fields['datev']
    else:
      dates = fields['dateo']
    # Convert date stamps to datetime objects, filtering out dummy values.
    dates = stamp2datetime(dates)
    dates = np.ma.asarray(dates, dtype='datetime64[s]')
    dates.mask = np.isnat(dates)
    # Where there are dummy dates, ignore the forecast information too.
    forecast = np.ma.asarray(fields['forecast'])
    forecast.mask = np.ma.getmaskarray(forecast) | (np.ma.getmaskarray(dates) & (fields['deet'] == 0))
    fields['forecast'] = forecast
    fields['time'] = dates

  # Add time and forecast axes to the data stream.
  def _iter (self):
    from fstd2nc.mixins import _iter_type, _var_type
    from collections import OrderedDict
    import numpy as np
    # Keep track of all time and forecast axes found in the data.
    time_axes = set()
    forecast_axes = set()
    for var in super(Dates,self)._iter():
      if not isinstance(var,_iter_type):
        yield var
        continue
      if 'time' in var.axes:
        times = var.axes['time']
        if times not in time_axes:
          time_axes.add(times)
          atts = OrderedDict([('axis','T')])
          axes = OrderedDict([('time',var.axes['time'])])
          # Add the time axis to the data stream.
          yield _var_type('time',atts,axes,np.asarray(times,dtype='datetime64[s]'))
      if 'forecast' in var.axes:
        forecasts = var.axes['forecast']
        if forecasts not in forecast_axes:
          forecast_axes.add(forecasts)
          atts = OrderedDict(units='hours')
          axes = OrderedDict([('forecast',var.axes['forecast'])])
          # Add the forecast axis to the data stream.
          yield _var_type('forecast',atts,axes,np.asarray(forecasts))
      yield var

