
from pygeode.axis import Axis, XAxis, YAxis, ZAxis, Lat, Lon, Hybrid, Pres, Height
from pygeode.timeaxis import StandardTime
class Forecast(Axis): pass
class IAxis(Axis): pass
class JAxis(Axis): pass
class KAxis(Axis): pass

# Create a multi-dimensional variable from a set of records
from pygeode.var import Var
class FSTD_Var (Var):
  def __init__ (self, records, coords, squash_forecasts=False):
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

    # Construst the i,j,k,z axes

    iaxis, jaxis, kaxis, zaxis = None, None, None, None

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

    atts = dict(nomvar=nomvar, typvar=typvar, etiket=etiket)

    kaxis = KAxis(nk)

    # Lat/Lon grid?
    if grtyp == 'A':
      x0 = 0
      dx = 360./ni
      x = x0 + np.arange(ni) * dx
      if ig1 == 0:  # Global
        dy = 180./nj
        y_south = -90 + dy/2
        y_north = 90 - dy/2
      elif ig1 == 1:# Northern Hemisphere
        dy = 90./nj
        y_south = 0 + dy/2
        y_north = 90 - dy/2
      else:         # Southern Hemisphere
        dy = 90./nj
        y_south = -90 + dy/2
        y_north = 0 - dy/2

      if ig2 == 0:  # South -> North
        y = y_south + np.arange(nj) * dy
      else:         # North -> South
        y = y_north - np.arange(nj) * dy

      iaxis = Lon(values=x)
      jaxis = Lat(values=y)

      del x0, dx, y_south, y_north, dy, x, y

    elif grtyp == 'B':
      x0 = 0
      dx = 360./(ni-1)
      x = x0 + np.arange(ni) * dx
      if ig1 == 0:  # Global
        dy = 180./(nj-1)
        y_south = -90
        y_north = 90
      elif ig1 == 1:# Northern Hemisphere
        dy = 90./(nj-1)
        y_south = 0
        y_north = 90
      else:         # Southern Hemisphere
        dy = 90./(nj-1)
        y_south = -90
        y_north = 0

      if ig2 == 0:  # South -> North
        y = y_south + np.arange(nj) * dy
      else:         # North -> South
        y = y_north - np.arange(nj) * dy

      iaxis = Lon(values=x)
      jaxis = Lat(values=y)

      del x0, dx, y_south, y_north, dy, x, y

    elif grtyp == 'G':
      from pygeode.quadrulepy import legendre_compute
      from math import pi
      x0 = 0
      dx = 360./ni
      x = x0 + np.arange(ni) * dx

      y, w = legendre_compute(nj)
      y = np.arcsin(y) / pi + 0.5  # Range (0,1)

      if ig1 == 0:  # Global
        y = -90 + y*180
      elif ig1 == 1:# Northern Hemisphere
        y = 0 + y*90
      else:         # Southern Hemisphere
        y = -90 + y*90

      if ig2 == 0:  # South -> North
        pass
      else:         # North -> South
        y = y[::-1]

      iaxis = Lon(values=x)
      jaxis = Lat(values=y)

      del x0, dx, x, y, w

    elif grtyp == 'L': pass #TODO

    else:
      from warnings import warn
      warn ("Unable to attach meaningful horizontal axes to %s"%name)
      iaxis = IAxis(ni)
      jaxis = JAxis(nj)

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
    from pygeode.var import Var
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


# Open a file for read access.  Returns a generic 'Dataset' object.
def open (filename):

  from pygeode.formats import fstd_core
  import numpy as np

  # What header attributes define a unique variable
  unique_var_atts = ['nomvar', 'typvar', 'etiket', 'ni', 'nj', 'nk', 'grtyp', 'ip3', 'ig1', 'ig2', 'ig3', 'ig4']

  # Read the records
  records = fstd_core.read_records(filename)

  # Pull out the coordinate records
  nomvar = records['nomvar']
  is_coord = (nomvar == '>>  ') | (nomvar == '^^  ') | (nomvar == 'HY  ') | (nomvar == '!!  ')
  coords = records[is_coord]
  records = records[-is_coord]
  del nomvar, is_coord

  # Preload the coords
  map(preload, coords)

  # Group the records together
  all_keys = records[unique_var_atts]
  unique_keys, var_indices = np.unique(all_keys, return_inverse=True)

  var_bins = [ records[var_indices==i] for i in range(len(unique_keys)) ]

  # Create the variables
  varlist = []
  for var_records in var_bins:
    is_var_coord = coords[['ip1','ip2','ip3']] == var_records[['ig1','ig2','ig3']][0]
    var_coords = coords[is_var_coord]
    var = FSTD_Var (var_records, var_coords)
    varlist.append(var)
    del var_records, is_var_coord, var_coords

  #TODO

