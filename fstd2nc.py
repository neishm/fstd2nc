#!/usr/bin/env python

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


"""
Functionality for converting between FSTD and netCDF files.
"""

__version__ = "0.20170728.1"

# Enable multi-language support.
from gettext import gettext as _
import gettext
from os import path, environ
try:
  import fstd2nc_locale
  gettext.bindtextdomain('fstd2nc', path.dirname(fstd2nc_locale.__file__))
  del fstd2nc_locale
except ImportError: pass
gettext.textdomain('fstd2nc')
# Check for special CMCLNG environment variable
if environ.get('CMCLNG') == 'francais':
  environ['LANGUAGE'] = 'fr_CA'
del gettext, path, environ

# Check for bundled rpnpy package.
# Fall back to this one if no standard rpnpy package available.
try:
  # Importing the module will set up the appropriate search paths.
  import fstd2nc_deps
  # Don't need a reference to the module after the paths are set.
  del fstd2nc_deps
except ImportError:
  pass

# How to handle warning messages.
# E.g., can either pass them through warnings.warn, or simply print them.
def warn (msg, _printed=set()):
  if msg not in _printed:
    print (_("Warning: %s")%msg)
    _printed.add(msg)

# Override default dtype for "binary" data.
# The one example I've seen has "float32" data encoded in it.
# https://wiki.cmc.ec.gc.ca/wiki/Talk:Python-RPN/2.0/examples#Plot_GIOPS_Forecast_Data_with_Basemap
def dtype_fst2numpy (datyp, nbits=None):
  from rpnpy.librmn.fstd98 import dtype_fst2numpy
  if datyp == 0:
    warn (_("Raw binary records detected.  The values may not be properly decoded if you're opening on a different platform."))
    datyp = 5
  return dtype_fst2numpy(datyp,nbits)

# Modify gdll to handle supergrids.
# Simply loops over each subgrid, then re-stacks them together.
def gdll (gdid):
  from rpnpy.librmn.interp import ezgprm, ezget_subgridids, gdll
  import numpy as np
  grtyp = ezgprm(gdid)['grtyp']
  if grtyp != 'U':
    return gdll(gdid)
  lat = []
  lon = []
  for subgrid in ezget_subgridids(gdid):
    subll = gdll(subgrid)
    lat.append(subll['lat'])
    lon.append(subll['lon'])
  # Stack them together again, to match the shape of the variable.
  lat = np.concatenate(lat, axis=1)
  lon = np.concatenate(lon, axis=1)
  return {'lat':lat, 'lon':lon}

# Modify gdgaxes to handle supergrids.
# Simply loops over each subgrid, then re-stacks them together.
def gdgaxes (gdid):
  from rpnpy.librmn.interp import ezgprm, ezget_subgridids, gdgaxes
  import numpy as np
  grtyp = ezgprm(gdid)['grtyp']
  if grtyp != 'U':
    return gdgaxes(gdid)
  ax = []
  ay = []
  for subgrid in ezget_subgridids(gdid):
    subll = gdgaxes(subgrid)
    ax.append(subll['ax'])
    ay.append(subll['ay'])
  # Stack them together again, to match the shape of the variable.
  ax = np.concatenate(ax, axis=1)
  ay = np.concatenate(ay, axis=1)
  return {'ax':ax, 'ay':ay}

# Convert an RPN date stamp to datetime object.
# Returns None for invalid stamps.
def stamp2datetime (date, cache={}):
  from rpnpy.rpndate import RPNDate
  dummy_stamps = (0, 10101011)
  if date not in cache:
    if date not in dummy_stamps:
      cache[date] = RPNDate(d).toDateTime().replace(tzinfo=None)
    else:
      cache[date] = None
  return cache[date]



# The type of data returned by the Buffer iterator.
class _var_type (object):
  __slots__ = ('name','atts','axes','array')
  def __init__ (self, name, atts, axes, array):
    self.name = name
    self.atts = atts
    self.axes = axes
    self.array = array
  def __iter__ (self):
    return (getattr(self,n) for n in self.__slots__)

# An internal type used in _iter methods.
class _iter_type (object):
  __slots__ = ('name','atts','axes','dtype','record_keys')
  def __init__ (self, name, atts, axes, dtype, record_keys):
    self.name = name
    self.atts = atts
    self.axes = axes
    self.dtype = dtype
    self.record_keys = record_keys
  def __iter__ (self):
    return (getattr(self,n) for n in self.__slots__)

# Helper method - modify the axes of an array.
def _modify_axes (axes, **kwargs):
  from collections import OrderedDict
  new_axes = OrderedDict()
  for name,values in axes.items():
    # Check if modifications requested.
    if name in kwargs:
      # Can either modify just the axis name, or the name and values too.
      if isinstance(kwargs[name],str):
        name = kwargs[name]
      else:
        name,values = kwargs[name]
    new_axes[name] = values
  return new_axes


# Define a class for encoding / decoding FSTD data.
# Each step is placed in its own "mixin" class, to make it easier to patch in 
# new behaviour if more exotic FSTD files are encountered in the future.
class _Buffer_Base (object):

  # Names of records that should be kept separate (never grouped into
  # multidimensional arrays).
  _meta_records = ()

  # Attributes which could potentially be used as axes.
  _outer_axes = ()

  # Attributes which uniquely identify a variable.
  _var_id = ('nomvar','ni','nj','nk')

  # Similar to above, but a human-readable version of the id.
  # Could be used for appending suffixes to variable names to make them unique.
  # Uses string formatting operations on the variable metadata.
  _human_var_id = ('%(nomvar)s', '%(ni)sx%(nj)s', '%(nk)sL')

  # Record parameters which should not be used as nc variable attributes.
  # (They're either internal to the file, or part of the data, not metadata).
  _ignore_atts = ('swa','lng','dltf','ubc','xtra1','xtra2','xtra3','key','shape','d')

  # Define any command-line arguments for reading FSTD files.
  @classmethod
  def _cmdline_args (cls, parser):
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('--minimal-metadata', action='store_true', help=_("Don't include RPN record attributes and other internal information in the output metadata."))
    parser.add_argument('--ignore-typvar', action='store_true', help=_('Tells the converter to ignore the typvar when deciding if two records are part of the same field.  Default is to split the variable on different typvars.'))
    parser.add_argument('--ignore-etiket', action='store_true', help=_('Tells the converter to ignore the etiket when deciding if two records are part of the same field.  Default is to split the variable on different etikets.'))

  # Do some checks on the command-line arguments after parsing them.
  @classmethod
  def _check_args (cls, parser, args):
    return  # Nothing to check here.

  # Clean up a buffer (close any attached files, etc.)
  def __del__ (self):
    from rpnpy.librmn.fstd98 import fstcloseall
    if hasattr(self,'_funit'):
      istat = fstcloseall(self._funit)

  # Extract metadata from a particular header.
  def _get_header_atts (self, header):
    for n,v in header.items():
      if n in self._ignore_atts: continue
      if isinstance(v,str):
        v = v.strip()
      yield (n,v)

  # Determine the dtype for an array.
  # Input: array of record indices.
  # Output: dtype that covers all datyps from the records.
  def _get_dtype (self, record_keys):
    import numpy as np
    dtype_list = [dtype_fst2numpy(self._params[k]['datyp'],self._params[k]['nbits']) for k in record_keys.flatten() if k >= 0]
    return np.result_type(*dtype_list)


  ###############################################
  # Basic flow for reading data

  def __init__ (self, filename, minimal_metadata=False, ignore_typvar=False, ignore_etiket=False):
    """
    Read raw records from FSTD files, into the buffer.
    Multiple files can be read simultaneously.
    """
    from rpnpy.librmn.fstd98 import fstopenall, fstinl, fstprm
    from rpnpy.librmn.const import FST_RO
    from collections import OrderedDict
    self._minimal_metadata = minimal_metadata
    if not ignore_typvar:
      # Insert typvar value just after nomvar.
      self._var_id = self._var_id[0:1] + ('typvar',) + self._var_id[1:]
      self._human_var_id = self._human_var_id[0:1] + ('%(typvar)s',) + self._human_var_id[1:]
    if not ignore_etiket:
      # Insert etiket value just after nomvar.
      self._var_id = self._var_id[0:1] + ('etiket',) + self._var_id[1:]
      self._human_var_id = self._human_var_id[0:1] + ('%(etiket)s',) + self._human_var_id[1:]

    # Scan the record headers.
    if isinstance(filename,str):
      filelist = [filename]
    else:
      filelist = list(filename)
    funit = fstopenall(filelist, FST_RO)
    keys = fstinl(funit)
    self._params = OrderedDict()
    for i,key in enumerate(keys):
      prm = fstprm(key)
      self._params[int(key)] = prm
    self._funit = funit

  # Collect the list of params from all the FSTD records, and concatenate them
  # into arrays.  Result is a single dictionary containing the vectorized
  # parameters of all records.
  # The sole purpose of this routine is to put the metadata in a structure
  # that's more efficient for doing bulk (vectorized) manipulations, instead
  # of manipulating each record param dictionary one at a time.
  # Subclasses may also insert extra (non-FSTD) parameters in here which are
  # needed for doing the variable decoding.
  def _vectorize_params (self):
    from collections import OrderedDict
    import numpy as np
    # Make sure the parameter names are consistent for all records.
    if len(set(map(frozenset,self._params.values()))) != 1:
      raise ValueError(("Inconsistent parameter names for the records."))
    fields = OrderedDict()
    for k,prm in self._params.items():
      for n,v in prm.items():
        fields.setdefault(n,[]).append(v)
    for n,v in list(fields.items()):
      # Treat array-like objects specially (don't want numpy to try to
      # read the data).
      if hasattr(v[0],'__array__'):
        v2 = v
        v = np.empty(len(self._params),dtype='O')
        v[:] = v2
      fields[n] = np.asarray(v)
    return fields


  # Internal iterator.
  def _iter (self):
    from collections import OrderedDict, namedtuple
    import numpy as np

    # Degenerate case: no data in buffer
    if len(self._params) == 0: return

    records = self._vectorize_params()

    # Group the records by variable.
    var_records = OrderedDict()
    var_id_type = namedtuple('var_id', self._var_id)
    for i,k in enumerate(records['key']):
      # Get the unique variable identifiers.
      var_id = var_id_type(*[records[n][i] for n in self._var_id])
      # Ignore coordinate records.
      if var_id.nomvar.strip() in self._meta_records: continue
      var_records.setdefault(var_id,[]).append(i)

    # Loop over each variable and construct the data & metadata.
    for var_id, rec_ids in var_records.items():
      # Get the metadata for each record.
      atts = OrderedDict()
      for n in records.keys():
        if n in self._outer_axes or n in self._ignore_atts: continue
        v = records[n][rec_ids]
        # Remove missing values before continuing.
        v = np.ma.compressed(v)
        if len(v) == 0: continue
        # Only use attributes that are consistent across all variable records.
        if len(set(v)) > 1: continue
        v = v[0]
        # Trim string attributes (remove whitespace padding).
        if isinstance(v,str): v = v.strip()
        atts[n] = v

      # Get the axis coordinates.
      axes = OrderedDict()
      for n in self._outer_axes:
        values = records[n][rec_ids]
        # Remove missing values before continuing.
        values = np.ma.compressed(values)
        # Ignore axes that have no actual coordinate values.
        if len(values) == 0: continue
        # Get all unique values (sorted).
        values = tuple(sorted(set(values)))
        axes[n] = values

      # Construct a multidimensional array to hold the record keys.
      record_keys = np.empty(map(len,axes.values()), dtype='int32')

      # Assume missing data (nan) unless filled in later.
      record_keys[()] = -1
      
      # Arrange the record keys in the appropriate locations.
      for rec_id in rec_ids:
        index = tuple(axes[n].index(records[n][rec_id]) for n in axes.keys())
        record_keys[index] = records['key'][rec_id]

      # Check if we have full coverage along all axes.
      have_data = [k >= 0 for k in record_keys.flatten()]
      if not np.all(have_data):
        warn (_("Missing some records for %s.")%var_id.nomvar)

      # Add dummy axes for the ni,nj,nk record dimensions.
      axes['k'] = tuple(range(var_id.nk))
      axes['j'] = tuple(range(var_id.nj))
      axes['i'] = tuple(range(var_id.ni))

      var = _iter_type( name = var_id.nomvar.strip(), atts = atts,
                        axes = axes,
                        dtype = self._get_dtype(record_keys),
                        record_keys = record_keys )
      yield var

      #TODO: Find a minimum set of partial coverages for the data.
      # (e.g., if we have surface-level output for some times, and 3D output
      # for other times).

  # How to read the data.
  @staticmethod
  def _fstluk (key, dtype=None, rank=None, dataArray=None):
    from rpnpy.librmn.fstd98 import fstluk
    # Use local overrides for dtype.
    if dtype is None and isinstance(key,dict):
      dtype = dtype_fst2numpy(key['datyp'],key['nbits'])
    return fstluk (key, dtype, rank, dataArray)

  #
  ###############################################




#####################################################################
# Mixins for different features / behaviour for the conversions.


#################################################
# Selecting for particular fields.
class _SelectVars (_Buffer_Base):
  @classmethod
  def _cmdline_args (cls, parser):
    super(_SelectVars,cls)._cmdline_args(parser)
    parser.add_argument('--vars', metavar='VAR1,VAR2,...', help=_('Comma-separated list of variables to convert.  By default, all variables are converted.'))
  def __init__ (self, *args, **kwargs):
    vars = kwargs.pop('vars',None)
    if vars is not None:
      self._selected_vars = vars.split(',')
      print (_('Will look for variables: ') + ' '.join(self._selected_vars))
    else:
      self._selected_vars = None
    super(_SelectVars,self).__init__(*args,**kwargs)
  def _iter (self):
    found = set()
    for var in super(_SelectVars,self)._iter():
      if self._selected_vars is not None:
        if var.name not in self._selected_vars:
          continue
      found.add(var.name)
      yield var
    if self._selected_vars is None: return
    missing = set(self._selected_vars) - found
    if len(missing) > 0:
      warn(_('Unable to find variables: ') + ' '.join(missing))


#################################################
# Pre-filtering the RPN file records.
class _FilterRecords (_Buffer_Base):
  @classmethod
  def _cmdline_args (cls, parser):
    super(_FilterRecords,cls)._cmdline_args(parser)
    parser.add_argument('--filter', metavar='CONDITION', action='append', help=_("Subset RPN file records using the given criteria.  For example, to convert only 24-hour forecasts you could use --filter ip2==24"))
  def __init__ (self, *args, **kwargs):
    filter = kwargs.pop('filter',None)
    if filter is None:
      filter = []
    self._filters = tuple(filter)
    super(_FilterRecords,self).__init__(*args,**kwargs)
  @staticmethod
  def _do_filter (p, cmd):
    try:
      return eval(cmd, None, p)
    except SyntaxError:
      print (_("Error: unable to parse the filter: %s")%cmd)
      exit(1)
    except NameError as e:
      print (_("Error: %s")%e.message)
      exit(1)
  def _vectorize_params (self):
    import numpy as np
    records = super(_FilterRecords,self)._vectorize_params()
    if len(self._filters) == 0: return records
    for cmd in self._filters:
      flags = []
      # Loop over each record, and apply the filter.
      for i in range(len(records['nomvar'])):
        p = dict((k,v[i]) for k,v in records.items())
        flags.append(self._do_filter(p, cmd))
      flags = np.array(flags)
      for k,v in list(records.items()):
        try:
          records[k] = v[flags]
        except IndexError:
          print (_("Error: unable to apply the filter: %s")%cmd)
          exit(1)
    return records


#################################################
# Logic for handling masks.

class _Masks (_Buffer_Base):
  @classmethod
  def _cmdline_args (cls, parser):
    super(_Masks,cls)._cmdline_args(parser)
    parser.add_argument('--fill-value', type=float, default=1e30, help=_("The fill value to use for masked (missing) data.  Gets stored as '_FillValue' attribute in the metadata.  Default is '%(default)s'."))
  def __init__ (self, *args, **kwargs):
    self._fill_value = kwargs.pop('fill_value',1e30)
    super(_Masks,self).__init__(*args,**kwargs)
    # Look for any mask records, and attach them to the associated field.
    masks = dict()
    # Look for masks, pull them out of the list of records.
    for key,prm in list(self._params.items()):
      if prm['typvar'] == '@@':
        masks[(prm['datev'],prm['etiket'],prm['ip1'],prm['ip2'],prm['nomvar'])] = key
        del self._params[key]
    # Store the mask records.
    self._masks = dict()
    for key,prm in list(self._params.items()):
      mask_key = masks.get((prm['datev'],prm['etiket'],prm['ip1'],prm['ip2'],prm['nomvar']),None)
      if mask_key is not None:
        self._masks[key] = mask_key

  # Apply the fill value to the data.
  def _iter (self):
    for var in super(_Masks,self)._iter():
      if not isinstance(var,_iter_type):
        yield var
        continue
      if any(k in self._masks for k in var.record_keys.flatten() if k >= 0):
        var.atts['_FillValue'] = var.dtype.type(self._fill_value)
      yield var

  # Apply the mask data
  def _fstluk (self, key, dtype=None, rank=None, dataArray=None):
    from rpnpy.librmn.fstd98 import fstluk
    import numpy as np
    prm = super(_Masks,self)._fstluk(key, dtype, rank, dataArray)
    if isinstance(key,dict):
      key = int(key['key'])
    mask_key = self._masks.get(key,None)
    if mask_key is not None:
      mask = fstluk(mask_key, rank=rank)['d']
      prm['d'] = np.ma.asarray(prm['d'])
      prm['d'].mask = (mask == 0)
    return prm


#################################################
# Logic for handling date field.

class _Dates (_Buffer_Base):
  @classmethod
  def _cmdline_args (cls, parser):
    super(_Dates,cls)._cmdline_args(parser)
    parser.add_argument('--squash-forecasts', action='store_true', help=_('Use the date of validity for the "time" axis.  Otherwise, the default is to use the date of original analysis, and the forecast length goes in a "forecast" axis.'))

  def __init__ (self, *args, **kwargs):
    self._squash_forecasts = kwargs.pop('squash_forecasts',False)
    if self._squash_forecasts:
      self._outer_axes = ('time',) + self._outer_axes
    else:
      self._outer_axes = ('time','forecast') + self._outer_axes
    super(_Dates,self).__init__(*args,**kwargs)

  # Get any extra (derived) fields needed for doing the decoding.
  def _vectorize_params (self):
    import numpy as np
    fields = super(_Dates,self)._vectorize_params()
    # Calculate the forecast (in hours).
    fields['forecast']=fields['deet']*fields['npas']/3600.
    # Time axis
    if self._squash_forecasts:
      dates = map(int,fields['datev'])
    else:
      dates = map(int,fields['dateo'])
    # Convert date stamps to datetime objects, filtering out dummy values.
    dates = map(stamp2datetime,dates)
    dates = np.ma.masked_equal(dates,None)
    # Where there are dummy dates, ignore the forecast information too.
    forecast = np.ma.asarray(fields['forecast'])
    forecast.mask = np.ma.getmaskarray(forecast) | (np.ma.getmaskarray(dates) & (fields['deet'] == 0))
    fields['forecast'] = forecast
    fields['time'] = dates
    return fields

  # Add time and forecast axes to the data stream.
  def _iter (self):
    from collections import OrderedDict
    import numpy as np
    # Keep track of all time and forecast axes found in the data.
    time_axes = set()
    forecast_axes = set()
    for var in super(_Dates,self)._iter():
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
          yield _var_type('time',atts,axes,np.asarray(times))
      if 'forecast' in var.axes:
        forecasts = var.axes['forecast']
        if forecasts not in forecast_axes:
          forecast_axes.add(forecasts)
          atts = OrderedDict(units='hours')
          axes = OrderedDict([('forecast',var.axes['forecast'])])
          # Add the forecast axis to the data stream.
          yield _var_type('forecast',atts,axes,np.asarray(forecasts))
      yield var



#################################################
# Logic for handling timeseries data.
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

class _Series (_Buffer_Base):
  def __init__ (self, *args, **kwargs):
    # Don't process series time/station/height records as variables.
    self._meta_records = self._meta_records + ('HH','STNS')
    # Add station # as another axis.
    self._outer_axes = ('station_id',) + self._outer_axes
    super(_Series,self).__init__(*args,**kwargs)

  def _vectorize_params (self):
    import numpy as np
    fields = super(_Series,self)._vectorize_params()
    nrecs = len(fields['nomvar'])
    # Identify timeseries records for further processing.
    is_series = (fields['typvar'] == 'T ') & ((fields['grtyp'] == '+') | (fields['grtyp'] == 'Y') | (fields['grtyp'] == 'T'))
    # More particular, data that has one station per record.
    is_split_series = (fields['typvar'] == 'T ') & (fields['grtyp'] == '+')

    # For timeseries data, station # is provided by 'ip3'.
    station_id = np.ma.array(fields['ip3'])
    # For non-timeseries data, ignore this info.
    station_id.mask = ~is_split_series
    fields['station_id'] = station_id
    # For timeseries data, the usual 'forecast' axis (from deet*npas) is not
    # used.  Instead, we will get forecast info from nj coordinate.
    if 'forecast' in fields:
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
    return fields
  def _iter (self):
    from rpnpy.librmn.fstd98 import fstluk
    from collections import OrderedDict
    import numpy as np
    from datetime import timedelta

    forecast_hours = None
    created_time_axis = False  # To only create squashed time axis once.
                               # (for --squashed-forecasts option).

    # Get station and forecast info.
    # Need to read from original records, because this into isn't in the
    # data stream.
    for header in self._params.values():
      if header['typvar'].strip() != 'T': continue
      nomvar = header['nomvar'].strip()
      # Create station axis.
      if nomvar == 'STNS':
        atts = OrderedDict()
        array = fstluk(header)['d'].transpose()
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
      if nomvar == 'HH':
        atts = OrderedDict(units='hours')
        array = fstluk(header)['d'].flatten()
        axes = OrderedDict(forecast=tuple(array))
        forecast = _var_type('forecast',atts,axes,array)
        forecast_hours = list(array)
        if getattr(self,'_squash_forecasts',False) is False:
          yield forecast

    for var in super(_Series,self)._iter():

      if not isinstance(var,_iter_type) or var.atts.get('typvar') != 'T':
        yield var
        continue

      # Vertical coordinates for series data.
      if var.name in ('SH','SV'):
        key = int(var.record_keys.flatten()[0])
        array = fstluk(key)['d'].squeeze()
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
          var.record_keys = var.record_keys.squeeze(axis=list(var.axes.keys()).index('level'))
          var.axes.pop('level')
        # ni is actually vertical level.
        # nj is actually forecast time.
        var.axes = _modify_axes(var.axes, i='level', j=('forecast',forecast_hours))

        # Some support for squashing forecasts.
        if getattr(self,'_squash_forecasts',False) is True:
          # Can only do this for a single date of origin, because the time
          # axis and forecast axis are not adjacent for this type of data.
          if len(var.axes['time']) == 1:
            var.record_keys = var.record_keys.squeeze(axis=list(var.axes.keys()).index('time'))
            time = var.axes.pop('time')[0]
            var.axes = _modify_axes(var.axes, forecast=('time',tuple(time+timedelta(hours=float(h)) for h in var.axes['forecast'])))
            if not created_time_axis:
              yield _var_type('time',{},{'time':var.axes['time']},np.array(var.axes['time']))
              created_time_axis = True
          else:
            warn(_("Can't squash forecast axis for timeseries data with multiple dates of origin."))
      # Remove 'kind' information for now - still need to figure out vertical
      # coordinates (i.e. how to map SV/SH here).
      var.atts.pop('kind',None)
      yield var

#################################################
# Logic for handling vertical coordinates.

class _VCoords (_Buffer_Base):
  _vcoord_nomvars = ('HY','!!')
  def __init__ (self, *args, **kwargs):
    # Use decoded IP1 values as the vertical axis.
    self._outer_axes = ('level',) + self._outer_axes
    # Tell the decoder not to process vertical records as variables.
    self._meta_records = self._meta_records + self._vcoord_nomvars
    super(_VCoords,self).__init__(*args,**kwargs)
    # Don't group records across different level 'kind'.
    # (otherwise can't create a coherent vertical axis).
    self._var_id = self._var_id + ('kind',)
    self._human_var_id = self._human_var_id + ('vgrid%(kind)s',)
  def _vectorize_params (self):
    from rpnpy.librmn.fstd98 import DecodeIp
    import numpy as np
    fields = super(_VCoords,self)._vectorize_params()
    # Provide 'level' and 'kind' information to the decoder.
    decoded = map(DecodeIp,fields['ip1'],fields['ip2'],fields['ip3'])
    rp1 = zip(*decoded)[0]
    levels = np.array([r.v1 for r in rp1])
    kind = np.array([r.kind for r in rp1])
    # Only use first set of levels (can't handle ranges yet).
    fields['level'] = levels
    fields['kind'] = kind
    return fields
  # Add vertical axis as another variable.
  def _iter (self):
    from collections import OrderedDict
    import numpy as np
    from rpnpy.vgd.base import vgd_fromlist, vgd_get, vgd_free
    from rpnpy.vgd.const import VGD_KEYS
    from rpnpy.vgd import VGDError
    from rpnpy.librmn.fstd98 import DecodeIp, fstluk
    # Pre-scan the raw headers for special vertical records.
    # (these aren't available in the data stream, because we told the decoder
    # to ignore them).
    vrecs = OrderedDict()
    for header in self._params.values():
      if header['nomvar'].strip() not in self._vcoord_nomvars: continue
      key = (header['ip1'],header['ip2'])
      # For old HY records, there's no matching ipX/igX codes.
      if header['nomvar'].strip() == 'HY': key = 'HY'
      if key in vrecs: continue
      vrecs[key] = header

    # Scan through the data, and look for any use of vertical coordinates.
    vaxes = OrderedDict()
    for var in super(_VCoords,self)._iter():
      # Degenerate vertical axis?
      if 'ip1' in var.atts and var.atts['ip1'] == 0:
        if 'level' in var.axes and len(var.axes['level']) == 1:
          i = list(var.axes).index('level')
          del var.axes['level']
          # Remove the axis if it's an "outer" axis (not along ni,nj,nk).
          if i < var.record_keys.ndim:
            var.record_keys = var.record_keys.squeeze(axis=i)
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
          atts['positive'] = 'down'
        elif kind == 2:
          # pressure [mb] (millibars)
          name = 'pres'
          atts['standard_name'] = 'air_pressure'
          atts['units'] = 'hPa'
          atts['positive'] = 'down'
        elif kind == 3:
          # arbitrary code
          name = 'code'
          atts.pop('axis',None)  # Not really a vertical axis?
        elif kind == 4:
          # height [M] (metres) with respect to ground level
          name = 'height'
          atts['standard_name'] = 'height'
          atts['units'] = 'm'
          atts['positive'] = 'up'
        elif kind == 5:
          # hybrid coordinates [hy] (0.0->1.0)
          atts['positive'] = 'down'
          key = (var.atts['ig1'],var.atts['ig2'])
          if header['nomvar'].strip() == 'HY': key = 'HY'
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
                A = _var_type(name+'_A', {}, {name:levels}, np.asarray(A))
                B = _var_type(name+'_B', {}, {name:levels}, np.asarray(B))
                ancillary_variables.extend([A,B])
              except (KeyError,ValueError,VGDError):
                warn (_("Unable to get A/B coefficients."))
              vgd_free (vgd_id)
            # Not a '!!' coordinate, so must be 'HY'?
            else:
              name = 'eta'
              atts['standard_name'] = 'atmosphere_hybrid_sigma_pressure_coordinate'
              # Get A and B info.
              eta = np.asarray(levels)
              ptop = DecodeIp(header['ip1'],0,0)[0].v1
              # Conversion taken from 'ig_to_hybref' function in librmn:
              pref = float(header['ig1'])
              rcoef = header['ig2']/1000.0
              # Apply the formula to compue A & B (from old fstd_core.c code):
              etatop = ptop/pref
              B = ((eta - etatop) / (1 - etatop)) ** rcoef
              A = pref * 100. * (eta - B)
              B = _var_type(name+'_B', {}, {name:levels}, B)
              A = _var_type(name+'_A', {}, {name:levels}, A)
              ancillary_variables.extend([A,B])
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
        array = np.asarray(levels)
        vaxes[(levels,kind)] = _var_type(name,atts,axes,array)
        yield vaxes[(levels,kind)]
        # Add any ancillary data needed for the axis.
        for anc in ancillary_variables:
          yield anc
      # Get the vertical axis.
      vaxis = vaxes[(levels,kind)]
      # Modify the variable's dimension name to match the axis name.
      var.axes = _modify_axes(var.axes, level=vaxis.name)
      yield var


#################################################
# Logic for handling lat/lon coordinates.

class _XYCoords (_Buffer_Base):
  # Special records that contain coordinate info.
  # We don't want to output these directly as variables, need to decode first.
  _xycoord_nomvars = ('^^','>>','^>')
  # These are extra copies of the lat/lon information.
  _xycoord_nomvars = _xycoord_nomvars + ('LA','LO')
  # Grids that can be read directly from '^^','>>' records, instead of going
  # through ezqkdef (in fact, may crash ezqkdef if you try decoding them).
  _direct_grids = ('X','Y','T','+')

  @classmethod
  def _cmdline_args (cls, parser):
    super(_XYCoords,cls)._cmdline_args(parser)
    parser.add_argument('--subgrid-axis', action='store_true', help=_('For data on supergrids, split the subgrids along a "subgrid" axis.  The default is to leave the subgrids stacked together as they are in the RPN file.'))

  def __init__ (self, *args, **kwargs):
    self._subgrid_axis = kwargs.pop('subgrid_axis',False)
    # Tell the decoder not to process horizontal records as variables.
    self._meta_records = self._meta_records + self._xycoord_nomvars
    super(_XYCoords,self).__init__(*args,**kwargs)
    # Variables must have an internally consistent horizontal grid.
    self._var_id = self._var_id + ('grtyp',)
    self._human_var_id = self._human_var_id + ('%(grtyp)s',)
    # Also, must have consistent igX records for a variable.
    if 'ig1' not in self._var_id:
      self._var_id = self._var_id + ('ig1','ig2','ig3','ig4')
      self._human_var_id = self._human_var_id + ('grid_%(ig1)s_%(ig2)s_%(ig3)s_%(ig4)s',)

  # Helper method - look up a coordinate record for the given variable.
  # Need this for manual lookup of 'X' grids, since ezqkdef doesn't support
  # them?
  def _find_coord (self, var, coordname):
    for header in self._params.values():
      if header['nomvar'].strip() != coordname: continue
      if header['ip1'] != var.atts['ig1']: continue
      if header['ip2'] != var.atts['ig2']: continue
      if header['ip3'] != var.atts['ig3']: continue
      return header
    raise KeyError("Unable to find matching '%s' for '%s'"%(coordname,var.name))


  # Add horizontal coordinate info to the data stream.
  def _iter (self):
    from collections import OrderedDict
    from rpnpy.librmn.interp import ezqkdef, EzscintError, ezget_nsubgrids
    from rpnpy.librmn.fstd98 import fstluk
    import numpy as np

    # Scan through the data, and look for any use of horizontal coordinates.
    grids = OrderedDict()
    for var in super(_XYCoords,self)._iter():
      # Don't touch derived variables.
      if not isinstance(var,_iter_type):
        yield var
        continue
      # Get grid parameters.
      ni = int(var.atts['ni'])
      nj = int(var.atts['nj'])
      grtyp = var.atts['grtyp']
      ig1 = int(var.atts['ig1'])
      ig2 = int(var.atts['ig2'])
      ig3 = int(var.atts['ig3'])
      ig4 = int(var.atts['ig4'])
      # Uniquely identify the grid for this variable.
      #
      # Use a looser identifier for timeseries data (ni/nj have different
      # meanings here (not grid-related), and could have multiple grtyp
      # values ('+','Y') that should share the same lat/lon info.
      if var.atts['typvar'].strip() == 'T':
        key = ('T',ig1,ig2)
      else:
        key = (grtyp,ni,nj,ig1,ig2,ig3,ig4)
      if grtyp in ('Y','+'): key = key[1:]
      # Check if we already defined this grid.
      if key not in grids:

        try:
          # Get basic information about this grid.
          # First, handle non-ezqkdef grids.
          if grtyp in self._direct_grids:
            lat = self._fstluk(self._find_coord(var,'^^'),rank=3)['d'].squeeze(axis=2)
            lon = self._fstluk(self._find_coord(var,'>>'),rank=3)['d'].squeeze(axis=2)
            ll = {'lat':lat, 'lon':lon}
          # Everything else should be handled by ezqkdef.
          else:
            gdid = ezqkdef (ni, nj, grtyp, ig1, ig2, ig3, ig4, self._funit)
            # For supergrids, need to loop over each subgrid to get the lat/lon
            ll = gdll(gdid)
        except (TypeError,EzscintError,KeyError):
          warn(_("Unable to get grid info for '%s'")%var.name)
          yield var
          continue


        # Find X/Y coordinates (if applicable).
        try:
          # Can't do this for direct grids (don't have a gdid defined).
          if grtyp in self._direct_grids: raise TypeError
          xycoords = gdgaxes(gdid)
          ax = xycoords['ax'].transpose()
          ay = xycoords['ay'].transpose()
          # Convert from degenerate 2D arrays to 1D arrays.
          ax = ax[0,:]
          ay = ay[:,0]
          xaxis = _var_type('x',{'axis':'X'},{'x':tuple(ax)},ax)
          yaxis = _var_type('y',{'axis':'Y'},{'y':tuple(ay)},ay)
        except (TypeError,EzscintError):
          # Can't get X/Y coords for this grid?
          xaxis = yaxis = None

        # Construct lat/lon fields.
        latarray = ll['lat'].transpose() # Switch from Fortran to C order.
        latatts = OrderedDict()
        latatts['long_name'] = 'latitude'
        latatts['standard_name'] = 'latitude'
        latatts['units'] = 'degrees_north'
        lonarray = ll['lon'].transpose() # Switch from Fortran to C order.
        lonatts = OrderedDict()
        lonatts['long_name'] = 'longitude'
        lonatts['standard_name'] = 'longitude'
        lonatts['units'] = 'degrees_east'

        # Case 1: lat/lon can be resolved into 1D Cartesian coordinates.
        # Calculate the mean lat/lon arrays in double precision.
        meanlat = np.mean(np.array(latarray,dtype=float),axis=1,keepdims=True)
        meanlon = np.mean(np.array(lonarray,dtype=float),axis=0,keepdims=True)
        if latarray.shape[1] > 1 and lonarray.shape[1] > 1 and np.allclose(latarray,meanlat) and np.allclose(lonarray,meanlon):
          # Reduce back to single precision for writing out.
          meanlat = np.array(meanlat,dtype=latarray.dtype).squeeze()
          meanlon = np.array(meanlon,dtype=lonarray.dtype).squeeze()
          # Ensure monotonicity of longitude field.
          # (gdll may sometimes wrap last longitude to zero).
          # Taken from old fstd_core.c code.
          if meanlon[-2] > meanlon[-3] and meanlon[-1] < meanlon[-2]:
            meanlon[-1] += 360.
          latarray = meanlat
          lonarray = meanlon
          lat = _var_type('lat',latatts,{'lat':tuple(latarray)},latarray)
          lon = _var_type('lon',lonatts,{'lon':tuple(lonarray)},lonarray)
          gridaxes = [('lat',tuple(latarray)),('lon',tuple(lonarray))]

        # Case 2: lat/lon are series of points.
        elif latarray.shape[0] == 1 and lonarray.shape[0] == 1 and ('i' in var.axes or 'station_id' in var.axes):
          latarray = latarray[0,:]
          lonarray = lonarray[0,:]
          # Special case for station data
          if 'station_id' in var.axes:
            gridaxes = [('station_id',var.axes['station_id'])]
          else:
            gridaxes = [('i',var.axes['i'])]
          lat = _var_type('lat',latatts,OrderedDict(gridaxes),latarray)
          lon = _var_type('lon',lonatts,OrderedDict(gridaxes),lonarray)

        # Case 3: General 2D lat/lon fields on X/Y coordinate system.
        elif xaxis is not None and yaxis is not None:
          gridaxes = [('y',tuple(yaxis.array)),('x',tuple(xaxis.array))]
          # Special case: have supergrid data, and the user wants to split it?
          if grtyp == 'U' and self._subgrid_axis:
            ngrids = ezget_nsubgrids(gdid)
            ny = len(yaxis.array)//ngrids
            yaxis.array = yaxis.array[:ny]
            yaxis.axes['y'] = tuple(yaxis.array)
            gridaxes = [('subgrid',tuple(range(ngrids))), ('y',tuple(yaxis.array)), ('x',tuple(xaxis.array))]
            latarray = latarray.reshape(ngrids,ny,-1)
            lonarray = lonarray.reshape(ngrids,ny,-1)
          yield yaxis
          yield xaxis
          lat = _var_type('lat',latatts,OrderedDict(gridaxes),latarray)
          lon = _var_type('lon',lonatts,OrderedDict(gridaxes),lonarray)

        # Case 4: General 2D lat/lon fields with no coordinate system.
        elif 'i' in var.axes and 'j' in var.axes:
          gridaxes = [('j',var.axes['j']),('i',var.axes['i'])]
          lat = _var_type('lat',latatts,OrderedDict(gridaxes),latarray)
          lon = _var_type('lon',lonatts,OrderedDict(gridaxes),lonarray)

        else:
          warn(_("Unhandled lat/lon coords for '%s'")%var.name)
          yield var
          continue

        yield lat
        yield lon
        grids[key] = gridaxes

      gridaxes = grids[key]

      # Update the var's horizontal coordinates.
      if len(gridaxes) == 1:
        var.axes = _modify_axes(var.axes, i=gridaxes[0])
      elif len(gridaxes) == 2:
        var.axes = _modify_axes(var.axes, j=gridaxes[0], i=gridaxes[1])
      elif len(gridaxes) == 3:
        var.axes = _modify_axes(var.axes, k=gridaxes[0], j=gridaxes[1], i=gridaxes[2])
      else:
        warn(_("Unusual grid axes for '%s' - ignoring.")%var.name)

      # For 2D lat/lon, need to reference them as coordinates in order for
      # netCDF viewers to display the field properly.
      if 'lat' not in var.axes or 'lon' not in var.axes:
        var.atts['coordinates'] = 'lon lat'

      yield var



#################################################
# Remove extraneous dimensions from the output.

class _NoNK (_Buffer_Base):
  def _iter (self):
    for var in super(_NoNK,self)._iter():
      if 'k' in var.axes and len(var.axes['k']) == 1:
        del var.axes['k']
      if 'j' in var.axes and len(var.axes['j']) == 1:
        del var.axes['j']
      yield var


#################################################
# Add netCDF metadata to the variables

class _netCDF_Atts (_Buffer_Base):
  @classmethod
  def _cmdline_args (cls, parser):
    import argparse
    super(_netCDF_Atts,cls)._cmdline_args(parser)
    parser.add_argument('--metadata-file', type=argparse.FileType('r'), action='append', help=_('Use metadata from the specified file.  You can repeat this option multiple times to build metadata from different sources.'))
  def __init__ (self, *args, **kwargs):
    import ConfigParser
    from collections import OrderedDict
    metadata_file = kwargs.pop('metadata_file',None)
    if metadata_file is None:
      metafiles = []
    else:
      metafiles = metadata_file
    metadata = OrderedDict()
    configparser = ConfigParser.SafeConfigParser()
    for metafile in metafiles:
      configparser.readfp(metafile)
    for varname in configparser.sections():
      metadata[varname] = OrderedDict(configparser.items(varname))
    self._metadata = metadata
    super(_netCDF_Atts,self).__init__(*args,**kwargs)
  def _iter (self):
    from collections import OrderedDict

    axis_renames = {}

    for var in super(_netCDF_Atts,self)._iter():
      orig_name = var.name
      # Add extra metadata provided by the user?
      if var.name in self._metadata:
        var.atts.update(self._metadata[var.name])
        # Rename the field? (Using special 'rename' key in the metadata file).
        if 'rename' in var.atts:
          var.name = var.atts.pop('rename')
          # Also rename any axis with this name.
          axis_renames[orig_name] = var.name

      # Check if any of the axes in this variable need to be renamed.
      axis_names, axis_values = zip(*var.axes.items())
      axis_names = [axis_renames.get(n,n) for n in axis_names]
      var.axes = OrderedDict(zip(axis_names,axis_values))

      yield var


#################################################
# Logic for reading/writing FSTD data from/to netCDF files.

class _netCDF_IO (_netCDF_Atts):
  @classmethod
  def _cmdline_args (cls, parser):
    super(_netCDF_IO,cls)._cmdline_args(parser)
    parser.add_argument('--time-units', choices=['seconds','minutes','hours','days'], default='hours', help=_('The units for the output time axis.  Default is %(default)s.'))
    parser.add_argument('--reference-date', metavar=_('YYYY-MM-DD'), help=_('The reference date for the output time axis.  The default is the starting date in the RPN file.'))

  @classmethod
  def _check_args (cls, parser, args):
    from datetime import datetime
    super(_netCDF_IO,cls)._check_args(parser,args)
    # Parse the reference date into a datetime object.
    if args.reference_date is not None:
      try:
        datetime.strptime(args.reference_date,'%Y-%m-%d')
      except ValueError:
        parser.error(_("Unable to to parse the reference date '%s'.  Expected format is '%s'")%(args.reference_date,_('YYYY-MM-DD')))

  def __init__ (self, *args, **kwargs):
    self._time_units = kwargs.pop('time_units','hours')
    self._reference_date = kwargs.pop('reference_date',None)
    self._unique_names = kwargs.pop('unique_names',True)
    super(_netCDF_IO,self).__init__(*args,**kwargs)

  def _iter (self):
    from datetime import datetime
    import numpy as np
    from netCDF4 import date2num

    if self._reference_date is None:
      reference_date = None
    else:
      reference_date = datetime.strptime(self._reference_date,'%Y-%m-%d')

    varlist = super(_netCDF_IO,self)._iter()
    if self._unique_names:
      varlist = list(varlist)
      self._fix_names(varlist)

    for var in varlist:

      # Modify time axes to be relative units instead of datetime objects.
      if var.name in var.axes and isinstance(var,_var_type) and isinstance(var.array[0],datetime):
        units = '%s since %s'%(self._time_units, reference_date or var.array[0])
        var.atts.update(units=units)
        var.array = np.asarray(date2num(var.array,units=units))

      yield var


  def _fix_names (self, varlist):
    from collections import OrderedDict

    # List of metadata keys that are internal to the FSTD file.
    internal_meta = list(self._vectorize_params().keys())

    # Generate unique axis names.
    axis_table = dict()
    for var in varlist:
      for axisname, axisvalues in var.axes.items():
        axisvalues = tuple(axisvalues)
        if axisname not in axis_table:
          axis_table[axisname] = []
        if axisvalues not in axis_table[axisname]:
          axis_table[axisname].append(axisvalues)
    axis_renames = dict()
    for axisname, axisvalues_list in axis_table.items():
      if len(axisvalues_list) == 1: continue
      warn (_("Multiple %s axes.  Appending integer suffixes to their names.")%axisname)
      for i,axisvalues in enumerate(axisvalues_list):
        axis_renames[(axisname,axisvalues)] = axisname+str(i+1)

    # Apply axis renames.
    def rename_axis ((axisname,axisvalues)):
      key = (axisname,tuple(axisvalues))
      if key in axis_renames:
        return (axis_renames[key],axisvalues)
      return (axisname,axisvalues)
    for var in varlist:
      # If this is a coordinate variable, use same renaming rules as the
      # dimension name.
      if isinstance(var,_var_type) and var.name in var.axes:
        var.name, axisvalues = rename_axis((var.name,var.axes[var.name]))
      var.axes = OrderedDict(map(rename_axis,var.axes.items()))

    # Generate a string-based variable id.
    # Only works for true variables from the FSTD source
    # (needs metadata like etiket, etc.)
    def get_var_id (var):
      out = []
      for fmt in self._human_var_id:
        out.append(fmt%var.atts)
      return tuple(out)

    # Generate unique variable names.
    var_table = dict()
    for i, var in enumerate(varlist):
      if var.name not in var_table:
        var_table[var.name] = []
      # Identify the variables by their index in the master list.
      var_table[var.name].append(i)
    for varname, var_indices in var_table.items():
      # Only need to rename variables that are non-unique.
      if len(var_indices) == 1: continue
      try:
        var_ids = [get_var_id(varlist[i]) for i in var_indices]
      except KeyError:
        # Some derived axes may not have enough metadata to generate an id,
        # so the best we can do is append an integer suffix.
        var_ids = [(str(r),) for r in range(1,len(var_indices)+1)]

      var_ids = zip(*var_ids)

      # Omit parts of the var_id that are invariant over all the variables.
      var_ids = [var_id for var_id in var_ids if len(set(var_id)) > 1]
      # Starting from the rightmost key, remove as many keys as possible while
      # maintaining uniqueness.
      for j in reversed(range(len(var_ids))):
        test = var_ids[:j] + var_ids[j+1:]
        if len(set(zip(*test))) == len(var_indices):
          var_ids = test

      var_ids = zip(*var_ids)

      var_ids = ['_'.join(var_id) for var_id in var_ids]

      warn (_("Multiple definitions of %s.  Adding unique suffixes %s.")%(varname, ', '.join(var_ids)))

      # Apply the name changes.
      for i, var_id in zip(var_indices, var_ids):
        var = varlist[i]
        orig_varname = var.name
        var.name = var.name + '_' + var_id
        # Apply the name changes to any metadata that references this variable.
        for othervar in varlist:
          # Must match axes.
          if not set(var.axes.keys()) <= set(othervar.axes.keys()): continue
          for key,val in list(othervar.atts.items()):
            if not isinstance(val,str): continue
            val = val.split()
            if orig_varname in val:
              val[val.index(orig_varname)] = var.name
            othervar.atts[key] = ' '.join(val)

    for var in varlist:
      # Names must start with a letter or underscore.
      if not var.name[0].isalpha():
        warn(_("Renaming '%s' to '_%s'.")%(var.name,var.name))
        var.name = '_'+var.name

      # Strip out FSTD-specific metadata?
      if self._minimal_metadata:
        for n in internal_meta:
          var.atts.pop(n,None)


  def write_nc_file (self, filename, nc_format='NETCDF4', global_metadata=None):
    """
    Write the records to a netCDF file.
    Requires the netCDF4 package.
    """
    from netCDF4 import Dataset
    import numpy as np
    f = Dataset(filename, "w", format=nc_format)

    # Apply global metadata (from config files and global_metadata argument).
    if 'global' in self._metadata:
      f.setncatts(self._metadata['global'])
    if global_metadata is not None:
      f.setncatts(global_metadata)

    for var in self._iter():

      for axisname, axisvalues in var.axes.items():
        # Only need to create each dimension once (even if it's in multiple
        # variables).
        if axisname not in f.dimensions:
          # Special case: make the time dimension unlimited.
          if axisname == 'time':
            f.createDimension(axisname, None)
          else:
            f.createDimension(axisname, len(axisvalues))

      dimensions = list(var.axes.keys())

      # Write the variable.
      # Easy case: already have the data.
      if isinstance(var,_var_type):
        v = f.createVariable(var.name, datatype=var.array.dtype, dimensions=dimensions)
        # Write the metadata.
        v.setncatts(var.atts)
        v[()] = var.array
        continue
      # Hard case: only have the record indices, need to loop over the records.
      v = f.createVariable(var.name, datatype=self._get_dtype(var.record_keys), dimensions=dimensions)
      # Write the metadata.
      v.setncatts(var.atts)
      # Write the data.
      data_shape = tuple(map(len,var.axes.values()))[var.record_keys.ndim:]
      for ind in np.ndindex(var.record_keys.shape):
        k = int(var.record_keys[ind])
        if k < 0: continue
        try:
          data = self._fstluk(k)['d'].transpose().reshape(data_shape)
          v[ind] = data
        except (IndexError,ValueError):
          warn(_("Internal problem with the script - unable to get data for '%s'")%var.name)
          continue
    # We need to explicitly state that we're using CF conventions in our
    # output files, or some utilities (like IDV) won't accept the data.
    f.Conventions = "CF-1.6"

    f.close()


#################################################
# A user-friendly iterator for using the Buffer in other Python scripts.

# Lightweight wrapper for the data.
# Allows the user to load the data through np.asarray or by slicing it.
class _Array (object):
  # Set some common attributes for the object.
  def __init__ (self, buffer, var):
    from functools import reduce
    self._buffer = buffer
    self._record_keys = var.record_keys
    # Expected shape and type of the array.
    self.shape = tuple(map(len,var.axes.values()))
    self.ndim = len(self.shape)
    self.size = reduce(int.__mul__, self.shape, 1)
    self.dtype = var.dtype
  def __getitem__ (self, key):
    import numpy as np
    # Coerce key into a tuple of slice objects.
    if not isinstance(key,tuple):
      if hasattr(key,'__len__'): key = tuple(key)
      else: key = (key,)
    if len(key) == 1 and hasattr(key[0],'__len__'):
      key = tuple(key[0])
    if Ellipsis in key:
      i = key.index(Ellipsis)
      key = key[:i] + (slice(None),)*(self.ndim-len(key)+1) + key[i+1:]
    key = key + (slice(None),)*(self.ndim-len(key))
    if len(key) > self.ndim:
      raise ValueError(("Too many dimensions for slicing."))
    final_shape = tuple(len(np.arange(n)[k]) for n,k in zip(self.shape,key) if not isinstance(k,int))
    data = np.ma.empty(final_shape, dtype=self.dtype)
    outer_ndim = self._record_keys.ndim
    record_keys = self._record_keys.__getitem__(key[:outer_ndim])
    # Final shape of each record (in case we did any reshaping along ni,nj,nk
    # dimensions).
    record_shape = self.shape[self._record_keys.ndim:]
    # Iterate of each record.
    for ind in np.ndindex(record_keys.shape):
      k = int(record_keys[ind])
      if k >= 0:
        data[ind] = self._buffer._fstluk(k)['d'].transpose().reshape(record_shape)[(Ellipsis,)+key[outer_ndim:]]
      else:
        data.mask = np.ma.getmaskarray(data)
        data.mask[ind] = True
    return data
  def __array__ (self):
    return self.__getitem__(())

class _Iter (_Buffer_Base):
  def __iter__ (self):
    """
    Processes the records into multidimensional variables.
    Iterates over (name, atts, axes, array) tuples.
    Note that array may not be a true numpy array (values are not yet loaded
    in memory).  To load the array, pass it to numpy.asarray().
    """
    for var in self._iter():
      if isinstance(var, _iter_type):
        array = _Array(self, var)
        var = _var_type (var.name, var.atts, var.axes, array)
      yield var


# Default interface for I/O.
class Buffer (_Iter,_netCDF_IO,_FilterRecords,_NoNK,_XYCoords,_VCoords,_Series,_Dates,_Masks,_SelectVars):
  """
  High-level interface for FSTD data, to treat it as multi-dimensional arrays.
  Contains logic for dealing with most of the common FSTD file conventions.
  """


# Command-line invocation:
def _fstd2nc_cmdline (buffer_type=Buffer):
  from argparse import ArgumentParser
  from sys import stdout, exit, argv
  from os.path import exists, isdir
  from glob import glob
  from rpnpy.librmn.fstd98 import isFST, FSTDError
  parser = ArgumentParser(description=_("Converts an RPN standard file (FSTD) to netCDF format."))
  parser.add_argument('infile', nargs='+', metavar='<fstd_file>', help=_('The RPN standard file(s) to convert.'))
  parser.add_argument('outfile', metavar='<netcdf_file>', help=_('The name of the netCDF file to create.'))
  buffer_type._cmdline_args(parser)
  parser.add_argument('--nc-format', choices=['NETCDF4','NETCDF4_CLASSIC','NETCDF3_CLASSIC','NETCDF3_64BIT_OFFSET','NETCDF3_64BIT_DATA'], default='NETCDF4', help=_('Which variant of netCDF to write.  Default is %(default)s.'))
  parser.add_argument('-f', '--force', action='store_true', help=_("Overwrite the output file if it already exists."))
  parser.add_argument('--no-history', action='store_true', help=_("Don't put the command-line invocation in the netCDF metadata."))
  args = parser.parse_args()
  buffer_type._check_args(parser, args)
  args = vars(args)
  infiles = args.pop('infile')
  outfile = args.pop('outfile')
  nc_format = args.pop('nc_format')
  force = args.pop('force')
  no_history = args.pop('no_history')

  # Apply wildcard expansion to filenames.
  infiles = [f for filepattern in infiles for f in (glob(filepattern) or filepattern)]

  # Make sure input file exists
  for infile in infiles:
    if not exists(infile):
      print (_("Error: '%s' does not exist!")%(infile))
      exit(1)
    if not isFST(infile) and not isdir(infile):
      print (_("Error: '%s' is not an RPN standard file!")%(infile))
      exit(1)

  try:
    buf = buffer_type(infiles, **args)
  except FSTDError:
    print (_("Error: problem opening one or more input directories."))
    exit(1)

  # Check if output file already exists
  if exists(outfile) and not force:
    overwrite = False
    if stdout.isatty():
      while True:
        print (_("Warning: '%s' already exists!  Overwrite? (y/n):")%(outfile)),
        try: ans = raw_input()
        except NameError: ans = input()
        if ans.lower() in ('y','yes','o','oui'):
          overwrite = True
          break
        if ans.lower() in ('n','no','non'):
          overwrite = False
          break
        print (_("Sorry, invalid response."))
    if overwrite is False:
      print (_("Refusing to overwrite existing file '%s'.")%(outfile))
      exit(1)

  # Append the command invocation to the netCDF metadata?
  if no_history:
    global_metadata = None
  else:
    from datetime import datetime
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    command = list(argv)
    # Any arguments with spaces should be surrounded by quotes.
    for i,c in enumerate(command):
      if " " in c:
        command[i] = "'"+c+"'"
    command = " ".join(command)
    history = timestamp + ": " + command
    global_metadata = {"history":history}

  buf.write_nc_file(outfile, nc_format, global_metadata)

# Command-line invocation with error trapping.
# Hides the Python stack trace when the user aborts the command.
def _fstd2nc_cmdline_trapped (*args, **kwargs):
  try:
    _fstd2nc_cmdline (*args, **kwargs)
  except KeyboardInterrupt:
    print (_("Aborted by user."))
    exit(1)


if __name__ == '__main__':
  _fstd2nc_cmdline_trapped()

