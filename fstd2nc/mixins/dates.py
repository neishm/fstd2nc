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
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--datev', '--squash-forecasts', action='store_true', default=True, dest='squash_forecasts', help=_('Use the date of validity for the "time" axis.  This is the default.'))
    group.add_argument('--dateo', '--forecast-axis', action='store_false', dest='squash_forecasts', help=_('Use the date of original analysis for the time axis, and put the forecast times into a separate "forecast" axis.'))

  # Need to extend _headers_dtype before __init__.
  def __new__ (cls, *args, **kwargs):
    obj = super(Dates,cls).__new__(cls, *args, **kwargs)
    obj._headers_dtype = obj._headers_dtype + [('time','datetime64[s]'),('leadtime','float32'),('reftime','datetime64[s]')]
    return obj

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
    dateo = stamp2datetime(fields['dateo'])
    datev = stamp2datetime(fields['datev'])
    # Convert date stamps to datetime objects, filtering out dummy values.
    dateo = np.ma.asarray(dateo, dtype='datetime64[s]')
    datev = np.ma.asarray(datev, dtype='datetime64[s]')
    dateo.mask = np.isnat(dateo)
    datev.mask = np.isnat(datev)
    # Where there are dummy dates, ignore the forecast information too.
    forecast = np.ma.asarray(fields['leadtime'])
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
    for coord in self._iter_coords():
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

