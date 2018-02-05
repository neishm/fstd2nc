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
  # Need to extend _headers_dtype before __init__.
  def __new__ (cls, *args, **kwargs):
    obj = super(Series,cls).__new__(cls, *args, **kwargs)
    obj._headers_dtype = obj._headers_dtype + [('station_id','int32')]
    return obj

  def __init__ (self, *args, **kwargs):
    import numpy as np
    # Don't process series time/station/height records as variables.
    self._meta_records = self._meta_records + ('HH','STNS')
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
    # For timeseries data, the usual 'forecast' axis (from deet*npas) is not
    # used.  Instead, we will get forecast info from nj coordinate.
    if 'forecast' in fields.dtype.names:
      fields['forecast'] = np.ma.asarray(fields['forecast'])
      fields['forecast'].mask = np.ma.getmaskarray(fields['forecast']) | is_series
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

    # Get station and forecast info.
    # Need to read from original records, because this into isn't in the
    # data stream.
    header = fstlir(self._meta_funit, nomvar='STNS')
    if header is not None:
        atts = OrderedDict()
        array = header['d'].transpose()
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
        axes = OrderedDict([('station_id',tuple(range(1,nstations+1))),('station_strlen',tuple(range(strlen)))])
        station = _var_type('station',atts,axes,array)
        yield station
    # Create forecast axis.
    header = fstlir (self._meta_funit, nomvar='HH')
    if header is not None:
        atts = OrderedDict(units='hours')
        array = header['d'].flatten()
        axes = OrderedDict(forecast=tuple(array))
        forecast = _var_type('forecast',atts,axes,array)
        forecast_hours = list(array)
        if getattr(self,'_squash_forecasts',False) is False:
          yield forecast

    for var in super(Series,self)._iter():

      if not isinstance(var,_iter_type) or var.atts.get('typvar') != 'T':
        yield var
        continue

      # Vertical coordinates for series data.
      if var.name in ('SH','SV'):
        ind = int(var.record_id.flatten()[0])
        array = self._fstluk(ind)['d'].squeeze()
        if array.ndim != 1: continue
        var.atts['kind'] = 5
        yield _var_type(var.name,var.atts,{'level':tuple(array)},array)
        continue

      # 'Y' data should be handled fine by _XYCoords - just give a more
      # specific name to the ni axis for clarity.
      if var.atts.get('grtyp') == 'Y':
        var.axes = _modify_axes(var.axes, i=('station_id',tuple(range(1,len(var.axes['i'])+1))))

      # '+' data has different meanings for the axes.
      if var.atts.get('grtyp') == '+' and forecast_hours is not None:
        # Remove degenerate vertical axis.
        if 'level' in var.axes:
          var.record_id = var.record_id.squeeze(axis=list(var.axes.keys()).index('level'))
          var.axes.pop('level')
        # ni is actually vertical level.
        # nj is actually forecast time.
        var.axes = _modify_axes(var.axes, i='level', j=('forecast',forecast_hours))

        # Some support for squashing forecasts.
        if getattr(self,'_squash_forecasts',False) is True:
          # Can only do this for a single date of origin, because the time
          # axis and forecast axis are not adjacent for this type of data.
          if len(var.axes['time']) == 1:
            var.record_id = var.record_id.squeeze(axis=list(var.axes.keys()).index('time'))
            time = var.axes.pop('time')[0]
            # Convert pandas times (if using pandas for processing the headers)
            time = np.datetime64(time,'s')
            var.axes = _modify_axes(var.axes, forecast=('time',tuple(time+np.timedelta64(int(h*3600),'s') for h in var.axes['forecast'])))
            if not created_time_axis:
              yield _var_type('time',{},{'time':var.axes['time']},np.array(var.axes['time']))
              created_time_axis = True
          else:
            warn(_("Can't squash forecast axis for timeseries data with multiple dates of origin."))
      # Remove 'kind' information for now - still need to figure out vertical
      # coordinates (i.e. how to map SV/SH here).
      var.atts.pop('kind',None)
      yield var

