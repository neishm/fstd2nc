
from pygeode.axis import Axis, XAxis, YAxis, ZAxis, Lat, Lon, Hybrid, Pres
from pygeode.timeaxis import StandardTime
class Dateo(Axis): pass   # Used for writing to FSTD files, not needed for reading
class Forecast(Axis): pass
class NPASAxis(Axis): pass  # Used for writing to FSTD files, not needed for reading
class IAxis(Axis): name = 'i'
class JAxis(Axis): name = 'j'
class KAxis(Axis): name = 'k'
class IP1Axis(Axis): name = 'ip1'

# Vertical coordinates
class Height_wrt_SeaLevel(ZAxis):
  name = 'height'
  atts = dict(ZAxis.atts, units='m', standard_name='height_above_sea_level')

class Height_wrt_Ground(ZAxis):
  name = 'height'
  atts = dict(ZAxis.atts, units='m', standard_name='height')

class Sigma(ZAxis):
  name = 'sigma'
  atts = dict(ZAxis.atts, standard_name='atmosphere_sigma_coordinate')
  plotatts = dict(ZAxis.plotatts, plotorder=-1, plotscale='log')

class LogHybrid(ZAxis):
  name = 'zeta'
  atts = dict(ZAxis.atts, standard_name='atmosphere_hybrid_sigma_log_pressure_coordinate')  # Not really a standard

class Theta(ZAxis):
  name = 'theta'
  atts = dict(ZAxis.atts, units='K', standard_name='air_potential_temperature')


# Create a multi-dimensional variable from a set of records
from pygeode.var import Var
class FSTD_Var (Var):
  def __init__ (self, records, squash_forecasts=False):
    from pygeode.var import Var
    from pygeode.formats import fstd_core
    import numpy as np

    name = str(records[0]['nomvar']).rstrip()

    # Get the dates, forecast hours, and levels.
    dates = fstd_core.stamp2date(records['dateo'])
    forecasts = records['deet']*records['npas']
    if squash_forecasts:
      dates += forecasts
      forecasts[:] = 0

    levels, kind = fstd_core.decode_levels(records['ip1'])

    # Get unique values of these arrays
    dates, idate = np.unique(dates, return_inverse=True)
    forecasts, iforecast = np.unique(forecasts, return_inverse=True)
    # Decode the numeric values of the levels, just to get to proper order
    levels, ilevel = np.unique(levels, return_index=True)
    ip1 = records['ip1'][ilevel]

    # Construct a multidimensional array of data functions.
    # One function per date,forecast,level.
    data_funcs = np.empty([len(dates),len(forecasts),len(levels)], dtype='O')
    data_funcs[:] = None
    data_funcs[idate,iforecast,ilevel] = records['data_func']

    if np.any(data_funcs == None):
      print "Unable to construct a full field for %s - missing some expected records:"%name
      from datetime import datetime, timedelta
      for i,date in enumerate(dates):
        # Convert to printable date
        date = int(date)
        date = datetime(year=1980,month=1,day=1) + timedelta(seconds=date)
        for j,forecast in enumerate(forecasts):
          forecast = int(forecast)/3600.
          for k,level in enumerate(levels):
            header = "%s - %3gh - level = %8g"%(date,forecast,level)
            if data_funcs[i,j,k] is None:
              print "%s: MISSING"%header
            else:
              print "%s: found"%header
      raise ValueError ("Missing some records that are needed to fill out the (time,forecast,level) dimensions.")

    self.data_funcs = data_funcs

    # Construct the time/forecast axes
    taxis = StandardTime (startdate={'year':1980,'month':1}, values=dates, units='seconds')
    faxis = Forecast (values=forecasts/3600.)
    faxis.atts['deet'] = int(records[0]['deet'])

    nomvar = str(records[0]['nomvar']).rstrip()
    typvar = str(records[0]['typvar']).rstrip()
    etiket = str(records[0]['etiket']).rstrip()
    grtyp = str(records[0]['grtyp']).rstrip()
    ni = int(records[0]['ni'])
    nj = int(records[0]['nj'])
    nk = int(records[0]['nk'])
    ig1 = int(records[0]['ig1'])
    ig2 = int(records[0]['ig2'])
    ig3 = int(records[0]['ig3'])
    ig4 = int(records[0]['ig4'])

    # Construst the i,j,k,z axes

    iaxis = IAxis(ni)
    jaxis = JAxis(nj)
    kaxis = KAxis(nk)

    atts = dict(nomvar=nomvar, typvar=typvar, etiket=etiket, grtyp=grtyp, ig1=ig1, ig2=ig2, ig3=ig3, ig4=ig4)

    # Vertical axis
    # (to be decoded later)
    zaxis = IP1Axis(values=ip1)

    # Determine the dtype to use
    # Use the first dtype found.
    datyp = int(records[0]['datyp'])
    nbits = int(records[0]['nbits'])
    dtype = {1:'float', 2:'uint', 3:'a', 4:'int', 5:'float', 134:'float', 130:'uint', 132:'int', 133:'float'}[datyp]
    if dtype == 'a':
      dtype += str(nbits/8)
    else:
      dtype += '64' if nbits > 32 else '32'

    # Finish initializing
    Var.__init__(self, [taxis,faxis,zaxis,kaxis,jaxis,iaxis], dtype=dtype, name=name, atts=atts)

  def getview (self, view, pbar):

    import numpy as np
    out = np.empty(view.shape, dtype=self.dtype)

    itimes = view.integer_indices[0]
    iforecasts = view.integer_indices[1]
    ilevels = view.integer_indices[2]

    sl_k = view.slices[3]
    sl_j = view.slices[4]
    sl_i = view.slices[5]

    for out_t, in_t in enumerate(itimes):
      for out_f, in_f in enumerate(iforecasts):
        for out_l, in_l in enumerate(ilevels):
          data = self.data_funcs[in_t,in_f,in_l]()
          data = data[sl_k,:,:]
          data = data[:,sl_j,:] # Avoid triggering "advanced" numpy indexing
          data = data[:,:,sl_i]
          out[out_t,out_f,out_l,:,:,:] = data

    return out

del Var

# Helper function - preload a record.
def preload(record):
  data = record['data_func']()
  def data_func(): return data
  record['data_func'] = data_func

# Attach lat/lon arrays to a list of FSTD variables
# For 1D lat/lon, use these arrays as axes.
# For 2D lat/lon, create 2D coordinate variables.
def attach_latlon (varlist, latlon_arrays):
  from pygeode.var import Var

  ydim, xdim = 4, 5

  handled_latlon_vars = {}
  extra_coord_vars = []

  for var in varlist:
    # Generate key into lat/lon lookup table
    key = tuple(var.atts[a] for a in ['grtyp','ig1','ig2','ig3','ig4'])
    key += (var.shape[xdim], var.shape[ydim])   # ni, nj
    if key not in latlon_arrays: continue

    # Attach 1D lat/lon arrays as axes
    axes = list(var.axes)
    ax, ay, lat, lon = latlon_arrays[key]
    if lat.ndim == 1: axes[ydim] = Lat(lat)
    else: axes[ydim] = YAxis(ay,name='y')
    if lon.ndim == 1: axes[xdim] = Lon(lon)
    else: axes[xdim] = XAxis(ax,name='x')
    var.axes = tuple(axes)

    # Convert 2D lat/lon arrays to variables
    if key in handled_latlon_vars: continue
    if lat.ndim == 2:
      lat_var = Var([axes[ydim],axes[xdim]], values=lat, name="lat")
      extra_coord_vars.append(lat_var)
    if lon.ndim == 2:
      lon_var = Var([axes[ydim],axes[xdim]], values=lon, name="lon")
      extra_coord_vars.append(lon_var)
    handled_latlon_vars[key] = True

  varlist.extend (extra_coord_vars)

# Attach vertical axes to the variables
def attach_vertical_axes (varlist, vertical_records):
  from pygeode.formats import fstd_core
  import numpy as np

  zdim = 2

  bangbang_cache = {}

  for var in varlist:
    # Skip derived fields
    if not isinstance(var,FSTD_Var): continue

    axes = list(var.axes)

    # Decode the values
    ip1 = var.axes[zdim].values
    levels, kind = fstd_core.decode_levels(ip1)

    if kind == 0:
      axes[zdim] = Height_wrt_SeaLevel(levels)
    elif kind == 1:
      axes[zdim] = Sigma(levels)
    elif kind == 2:
      axes[zdim] = Pres(levels)
    elif kind == 3:
      axes[zdim] = ZAxis(levels)
    elif kind == 4:
      axes[zdim] = Height_wrt_Ground(levels)
    elif kind == 5:
      # Find a vertical record that matches
      # First, look for a !! record
      match = (vertical_records['ip1'] == var.atts['ig1']) & (vertical_records['ip2'] == var.atts['ig2']) & (vertical_records['ip3'] == var.atts['ig3'])
      if any(match):
        bangbang_record = vertical_records[match]
        key = int(bangbang_record['ip1'][0]), int(bangbang_record['ip2'][0]), int(bangbang_record['ip3'][0])
        if key not in bangbang_cache:
          bangbang_cache[key] = fstd_core.get_loghybrid_table(bangbang_record)
        table = bangbang_cache[key]
        # Determine the specific A & B for this axis
        A, B = fstd_core.get_loghybrid_a_b(ip1, *table)
        axes[zdim] = LogHybrid(values=levels, A=A, B=B)
        # Store some auxiliary vertical information
        # (needed for going back to FSTD format)
        for i,a in enumerate(['ip1_m','a_m','b_m','ip1_t','a_t','b_t']):
          axes[zdim].atts[a] = table[i]

      else:
        # Otherwise, look for a HY record (looser search criteria)
        match = (vertical_records['nomvar'] == 'HY  ')
        if any(match):
          hy_record = vertical_records[match]
          ptop, rcoef, pref, A, B = fstd_core.get_hybrid_a_b(hy_record, levels)
          axes[zdim] = Hybrid(values=levels, A=A, B=B)
          axes[zdim].atts = dict(ptop=ptop, rcoef=rcoef, pref=pref)
    elif kind == 6:
      axes[zdim] = Theta(levels)

    # Save the new axis information back into the variable
    var.axes = tuple(axes)

# Reduce the dimensionality of the given FSTD variable
def reduce_dimensionality (var, squash_forecasts=False):
  # Skip derived fields
  if not isinstance(var, FSTD_Var): return var

  remove_axes = []
  # Forecast (axis 1)
  if var.shape[1] == 1:
    if squash_forecasts:
      remove_axes += [1]
  # Vertical (axis 2)
  # Surface fields have a 'pressure' coordinate with a value of 0hPa
  if var.shape[2] == 1:
    if isinstance(var.axes[2], Pres) and var.axes[2].values == [0.]:
      remove_axes += [2]
  # K axis (axis 3)
  if var.shape[3] == 1:
    remove_axes += [3]

  return var.squeeze(*remove_axes)

#####################################################################
#
# Open a file for read access.  Returns a generic 'Dataset' object.
#
#####################################################################
def open (filename, squash_forecasts=False):

  from pygeode.formats import fstd_core
  import numpy as np

  # What header attributes define a unique variable
  unique_var_atts = ['nomvar', 'typvar', 'etiket', 'ni', 'nj', 'nk', 'grtyp', 'ip3', 'ig1', 'ig2', 'ig3', 'ig4']

  # Read the records
  records = fstd_core.read_records(filename)

  # Construct all possible lat/lon arrays from info in the records
  latlon_arrays = fstd_core.get_latlon(records)

  # Remove the coordinate records (not needed anymore).
  # Keep vertical descriptors in a separate array, though, because we'll need
  # to evaluate them once we know the vertical dimension of the variables.
  nomvar = records['nomvar']
  is_coord = (nomvar == '>>  ') | (nomvar == '^^  ') | (nomvar == 'HY  ') | (nomvar == '!!  ')
  vertical_records = records[ (nomvar == 'HY  ') | (nomvar == '!!  ') ]
  records = records[-is_coord]
  del nomvar, is_coord

  # Group the records together
  all_keys = records[unique_var_atts]
  unique_keys, var_indices = np.unique(all_keys, return_inverse=True)

  var_bins = [ records[var_indices==i] for i in range(len(unique_keys)) ]

  del all_keys, unique_keys, var_indices

  # Create the variables
  varlist = []
  for var_records in var_bins:
    var = FSTD_Var (var_records, squash_forecasts)
    varlist.append(var)

  # Attach any lat/lon coordinates
  attach_latlon (varlist, latlon_arrays)

  # Attach vertical coordinates
  attach_vertical_axes (varlist, vertical_records)

  # Dimensionality reduction
  varlist = [reduce_dimensionality(var,squash_forecasts) for var in varlist]

  # Return the variables as a dataset
  from pygeode.dataset import Dataset
  return Dataset(varlist)

#####################################################################


# Encode the time axis
def encode_time_axis (varlist):
  from pygeode.timeaxis import StandardTime
  from pygeode.timeutils import reltime
  from pygeode.formats import fstd_core
  import numpy as np
  for i,var in enumerate(varlist):
    if not var.hasaxis(StandardTime): continue
    time = var.getaxis(StandardTime)
    seconds = reltime (time, startdate=dict(year=1980,month=1,day=1), units='seconds')
    seconds = np.asarray(seconds,dtype=int)
    values = fstd_core.date2stamp(seconds)
    taxis = Dateo(values=values)
    varlist[i] = var.replace_axes(time=taxis)

# Convert to FSTD-compatible axes
# (e.g., detect hybrid / log-hybrid axes, forecast axis)
def detect_fstd_axes (varlist):
  for varnum,var in enumerate(varlist):
    axes = list(var.axes)
    # Look for serialized FSTD axes (represented by metadata)
    for i,axis in enumerate(axes):
      if isinstance(axis,(Height_wrt_SeaLevel,Height_wrt_Ground,Sigma,Hybrid,LogHybrid,Theta,Forecast,XAxis,YAxis,KAxis)): continue
      standard_name = axis.atts.get('standard_name','')
      if standard_name == 'height_above_sea_level':
        axes[i] = Height_wrt_SeaLevel(values=axis.values)
      elif standard_name == 'height':
        axes[i] = Height_wrt_Ground(values=axis.values)
      elif standard_name == 'atmosphere_sigma_coordinate':
        axes[i] = Sigma(values=axis.values)
      elif standard_name == 'atmosphere_hybrid_sigma_pressure_coordinate':
        axes[i] = Hybrid(values=axis.values, A=axis.auxarrays['A'], B=axis.auxarrays['B'])
        axes[i].atts = dict(**axis.atts)
      elif standard_name == 'atmosphere_hybrid_sigma_log_pressure_coordinate':
        axes[i] = LogHybrid(values=axis.values, A=axis.auxarrays['A'], B=axis.auxarrays['B'])
      elif standard_name == 'air_potential_temperature':
        axes[i] = Theta(values=axis.values)

      elif axis.name == 'forecast':
        axes[i] = Forecast(values=axis.values, atts=axis.atts)

      replacements = {}
      for original_axis, new_axis in zip(var.axes,axes):
        if new_axis is not original_axis:
          replacements[original_axis.name] = new_axis
      if len(replacements) > 0:
        varlist[varnum] = var.replace_axes(**replacements)

# Encode the forecast axis
def encode_forecast_axis (varlist):
  for i,var in enumerate(varlist):
    if not var.hasaxis(Forecast): continue
    forecast = var.getaxis(Forecast)
    if 'deet' not in forecast.atts: continue
    deet = forecast.atts['deet']
    npas = forecast.values * 3600 / deet
    npas_axis = NPASAxis(values=npas, atts={'deet':deet})
    varlist[i] = var.replace_axes(forecast=npas_axis)

# Encode vertical information into FSTD records
def encode_vertical (varlist):
  vertical_records = []
  from pygeode.formats import fstd_core
  import numpy as np
  from warnings import warn

  # First, check for any hybrid records
  hy_records = {}
  for varnum,var in enumerate(varlist):
    if var.hasaxis(Hybrid):
      eta = var.getaxis(Hybrid)
      if 'ptop' not in eta.atts or 'rcoef' not in eta.atts or 'pref' not in eta.atts:
        warn ("Not enough information to construct an HY record");
        continue
      ptop = eta.atts['ptop']
      rcoef = eta.atts['rcoef']
      pref = eta.atts['pref']
      key = (ptop,rcoef,pref)
      if key not in hy_records:
        hy_record = fstd_core.make_hy_record(ptop,rcoef,pref)
        hy_record['etiket'] = var.atts.get('etiket','        ')
        # Use the same date as the field
        if var.hasaxis(Dateo):
          hy_record['dateo'] = var.getaxis(Dateo).values[0]
        # Use the same forecast time as the field
        if var.hasaxis(NPASAxis):
          npas_axis = var.getaxis(NPASAxis)
          deet = npas_axis.atts['deet']
          npas = npas_axis.values[0]
          hy_record['npas'] = npas
          hy_record['deet'] = deet
        # Dummy data function - just give a zero.
        hy_record['data_func'] = lambda: np.array([[[0]]],dtype='float32')
        hy_records[key] = hy_record

  if len(hy_records) > 1:
    warn ("Multiple Hybrid axes detected.  The resulting file may not work the way you expect.")

  vertical_records.extend(hy_records.values())

  #TODO: !! record encoding
  #TODO: convert to IP1Axis

  # Convert from list to array
  if len(vertical_records) > 0:
    vertical_records = np.concatenate(vertical_records)
  else:
    vertical_records = np.empty([0], dtype=fstd_core.record_descr)

  return vertical_records

# Encode latitude/longitude information into FSTD records
def encode_latlon (varlist):
  latlon_records = {}

  for var in varlist:
    pass
#TODO

# Check for incompatible axes
def check_fstd_axes (varlist):
  compatible_axes = StandardTime, Forecast, IP1Axis, KAxis, JAxis, IAxis
  for var in varlist:
    incompatible_axes = [a for a in var.axes if not isinstance(a,compatible_axes)]
    if len(incompatible_axes) > 0:
      raise TypeError, "Cannot fit the following axes from var '%s' into an FSTD structure: %s"%(var.name,incompatible_axes)

#TODO: check for repeated axes


#####################################################################
#
# Save a dataset into an FSTD file.
#
#####################################################################
def save (filename, varlist):
  from pygeode.dataset import Dataset
  if isinstance(varlist,Dataset):
    varlist = varlist.vars

  varlist = list(varlist)

  # Encode the time axes
  encode_time_axis (varlist)

  # Convert to FSTD-compatible axes
  # (e.g., detect hybrid / log-hybrid axes)
  detect_fstd_axes (varlist)

  # Encode the forecast axis
  encode_forecast_axis (varlist)

  # Extract vertical information
  vertical_records = encode_vertical(varlist)

  #TODO

  # Extract horizontal information
  #TODO

  # We should now have a subset of (StandardTime,Forecast,IP1Axis,KAxis,JAxis,IAxis)

  # Fill in missing degenerate axes, and put them in the expected order
  #TODO

  # Convert variables to record arrays
  #TODO

  # Merge coordinate records and variable records
  #TODO

  # Save to the file
  #TODO
  return varlist

