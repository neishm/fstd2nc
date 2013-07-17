
from pygeode.axis import Axis, XAxis, YAxis, ZAxis, Lat, Lon, Hybrid, Pres, Height
from pygeode.timeaxis import StandardTime
class Forecast(Axis): pass
class IAxis(Axis): pass
class JAxis(Axis): pass
class KAxis(Axis): pass

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
    levels, ilevel = np.unique(levels, return_inverse=True)

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
    zaxis = ZAxis(values=levels) #TODO

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

  for i,var in enumerate(varlist):
    # Generate key into lat/lon lookup table
    key = tuple(var.atts[a] for a in ['grtyp','ig1','ig2','ig3','ig4'])
    key += (var.shape[xdim], var.shape[ydim])   # ni, nj
    if key not in latlon_arrays: continue

    # Attach 1D lat/lon arrays as axes
    axes = list(var.axes)
    ax, ay, lat, lon = latlon_arrays[key]
    if lat.ndim == 1: axes[ydim] = Lat(lat)
    else: axes[ydim] = YAxis(ay)
    if lon.ndim == 1: axes[xdim] = Lon(lon)
    else: axes[xdim] = XAxis(ax)
    var.axes = tuple(axes)

    # Convert 2D lat/lon arrays to variables
    if key in handled_latlon_vars: continue
    if lat.ndim == 2:
      lat_var = Var([axes[ydim],axes[xdim]], values=lat, name="latitudes")
      extra_coord_vars.append(lat_var)
    if lon.ndim == 2:
      lon_var = Var([axes[ydim],axes[xdim]], values=lon, name="longitudes")
      extra_coord_vars.append(lon_var)
    handled_latlon_vars[key] = True

  varlist.extend (extra_coord_vars)

# Reduce the dimensionality of the given FSTD variable
def reduce_dimensionality (var, squash_forecasts=False):
  remove_axes = []
  # Forecast (axis 1)
  if var.shape[1] == 1:
    if squash_forecasts:
      remove_axes += [1]
  # Vertical (axis 2)
  if var.shape[2] == 1:
    if isinstance(var.axes[2], Height) and var.axes[2].values == [0]:
      remove_axes += [2]
  # K axis (axis 3)
  if var.shape[3] == 1:
    remove_axes += [3]

  return var.squeeze(*remove_axes)

# Open a file for read access.  Returns a generic 'Dataset' object.
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
  nomvar = records['nomvar']
  is_coord = (nomvar == '>>  ') | (nomvar == '^^  ') | (nomvar == 'HY  ') | (nomvar == '!!  ')
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

  # Dimensionality reduction
  varlist = [reduce_dimensionality(var,squash_forecasts) for var in varlist]

  # Return the variables as a dataset
  from pygeode.dataset import Dataset
  return Dataset(varlist)

