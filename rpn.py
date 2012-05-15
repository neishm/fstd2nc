from pygeode.tools import load_lib

lib = load_lib("plugins/rpn/librpn.so")

######################################################################
# This first section defines a basic interface to get the raw data out.
######################################################################

from ctypes import Structure, c_int, c_uint, c_char, c_longlong, c_ulonglong, c_float, POINTER, byref, c_void_p
Nomvar = c_char * 5
Etiket = c_char * 13
Typvar = c_char * 3


"""
typedef struct {
  int status;
  int size;
  unsigned long long int data;
  int deet;
  int npak;
  int ni;
  char grtyp;
  int nj;
  int datyp;
  int nk;
  int npas;
  int ig4;
  int ig2;
  int ig1;
  int ig3;
  Etiket etiket;
  Typvar typvar;
  Nomvar nomvar;
  int ip1;
  int ip2;
  int ip3;
  long long int dateo;
  unsigned int checksum;
} RecordHeader;
"""

class RecordHeader (Structure):
  _fields_ = [('status', c_int), ('size', c_int), ('data', c_ulonglong),
              ('deet', c_int), ('npak', c_int), ('ni', c_int),
              ('grtyp', c_char), ('nj', c_int), ('datyp', c_int),
              ('nk', c_int), ('npas', c_int), ('ig4', c_int),
              ('ig2', c_int), ('ig1', c_int), ('ig3', c_int),
              ('etiket', Etiket), ('typvar', Typvar), ('nomvar', Nomvar),
              ('ip1', c_int), ('ip2', c_int), ('ip3', c_int),
              ('dateo', c_longlong), ('checksum', c_uint)]


# What header attributes define a unique variable
unique_var_atts = 'nomvar', 'typvar', 'etiket', 'ni', 'nj', 'nk', 'grtyp', 'ip3', 'ig1', 'ig2', 'ig3', 'ig4'


from pygeode.axis import Axis

# Raw axes from the file
class T(Axis): pass  # time axis
class F(Axis): name = 'forecast'  # forecast axis
class Z(Axis): pass  # levels (over ip1)
class K(Axis): pass  # ?? (nk)
class J(Axis): pass  # latitudes?
class I(Axis): pass  # longitudes?

del Axis



# Collect a list of headers into a PyGeode Var
from pygeode.var import Var
class RPN_Var (Var):
  def __init__ (self, filename, headers):
    from pygeode.var import Var
    import numpy as np

    name = headers[0].nomvar.strip()

    # Get unique time steps, forecasts, and levels
    # dateo, ip2, ip1?
    dateo_list = []
    ip2_list = []
    ip1_list = []

    for h in headers:
      if h.dateo not in dateo_list:
        dateo_list.append(h.dateo)
      if h.ip2 not in ip2_list:
        ip2_list.append(h.ip2)
      if h.ip1 not in ip1_list:
        ip1_list.append(h.ip1)

    # Encode as integers
    dateo_list = np.array(dateo_list, dtype='int64')
    ip2_list = np.array(ip2_list, dtype='int64')
    ip1_list = np.array(ip1_list, dtype='int64')

    # Define the axes
    t = T(dateo_list)
    f = F(ip2_list)
    z = Z(ip1_list)
    k = K(headers[0].nk)
    j = J(headers[0].nj)
    i = I(headers[0].ni)

    # Set some metadata (things that uniquely define the variable)
    atts = {}
    for att in unique_var_atts:
      atts[att] = getattr(headers[0], att)

    # Generate a matrix for the order of headers
    header_order = np.empty((len(t),len(f),len(z)), dtype=int)
    header_order[:] = -1
    for hi, h in enumerate(headers):
      ti = np.where(dateo_list == h.dateo)[0][0]
      fi = np.where(ip2_list == h.ip2)[0][0]
      zi = np.where(ip1_list == h.ip1)[0][0]
      header_order[ti,fi,zi] = hi

    # Make sure we have all the records we need
    assert np.all(header_order >= 0), "some records are missing from variable '%s'"%name

    self.filename = filename
    self.headers = headers
    self.header_order = header_order

    Var.__init__ (self, [t,f,z,k,j,i], name=name, atts=atts, dtype='float32')

  from pygeode.tools import need_full_axes
  @need_full_axes(I,J,K)
  def getview (self, view, pbar):
    from pygeode.tools import point
    import numpy as np
    from ctypes import c_int32, byref

    out = np.empty(view.shape, self.dtype)

    # Time, level slicing relative to *whole* variable
    TI = view.integer_indices[0]
    FI = view.integer_indices[1]
    ZI = view.integer_indices[2]

    # Generate a list of all the headers we need
    headers = []
    for ti in TI:
      for fi in FI:
        for zi in ZI:
          hi = self.header_order[ti,fi,zi]
          headers.append(self.headers[hi])
    headers = (RecordHeader*len(headers))(*headers)

    recsize = self.shape[3] * self.shape[4] * self.shape[5]
    lib.read_data (self.filename, len(headers), headers, recsize, point(out))

    pbar.update(100)

    return out
  del need_full_axes

del Var


# Helper method - group headers based on what variable the represent
def collect_headers (headers):
  var_headers = {}
  for h in headers:
    # Generate a unique key for this variable
    key = tuple(getattr(h,att) for att in unique_var_atts)
    # Include part of the ip1 parameter (must have consistent vertical coordinate type)
    key += (h.ip1>>24,)

    if key not in var_headers:
      var_headers[key] = []
    # Add this header to the variable
    var_headers[key].append(h)

  return var_headers


# Open an RPN file, return a list of raw variables (minimal processing)
def rawopen (filename):
  from os.path import exists
  assert exists(filename), "Can't find %s"%filename
  nrecs = lib.get_num_records(filename)
  headers = (RecordHeader*nrecs)()
  lib.get_record_headers (filename, headers)
  vardict = collect_headers(headers)
  # Generate PyGeode vars based on the header information
  varlist = [RPN_Var(filename, headers) for headers in vardict.values()]
  return varlist


######################################################################
#  The stuff below wraps a basic RPN variable to have a more useful
#  geophysical representation.
######################################################################

# More specific sub-axes:
from pygeode.axis import ZAxis
class Height (ZAxis): pass
class Sigma (ZAxis):
  plotorder = -1
  plotscale = 'log'

from pygeode.axis import Pres

class Height_wrt_Ground (ZAxis): pass

#TODO: sub-class this to allow a Hybrid.fromvar()-type thing
from pygeode.axis import Hybrid

class Theta (ZAxis): pass

# How the ip1 types get mapped to the above classes:
#
# kind    : level type
# 0       : height [m] (metres)
# 1       : sigma  [sg] (0.0->1.0)
# 2       : pressure [mb] (millibars)
# 3       : arbitrary code
# 4       : height [M] (metres) with respect to ground level
# 5       : hybrid coordinates [hy] (0.0->1.0)
# 6       : theta [th]
#
vertical_type = {0:Height, 1:Sigma, 2:Pres, 3:ZAxis, 4:Height_wrt_Ground, 5:Hybrid, 6:Theta}

# Axis with an associated type
# Slightly more specific than I or J axis, as we now have coordinate values,
#  a grid type, and associated ig1, ig2, ig3, ig4.
# We still don't interpret the values in any geophysical sense, though.
from pygeode.axis import Axis
class Horizontal(Axis):
  @classmethod
  def fromvar (cls, coordvar):
    # Get the coordinate values
    values = coordvar.get().squeeze()

    # Get the metadata
    atts = dict(**coordvar.atts)
    zaxis = coordvar.getaxis('Z')
    assert len(zaxis) == 1
    atts['ip1'] = int(zaxis.values[0])
    faxis = coordvar.getaxis('F')
    assert len(faxis) == 1
    atts['ip2'] = int(faxis.values[0])

    # Instantiate
    axis = cls(values, atts=atts)
    return axis
del Axis

class XCoord(Horizontal): pass  # '>>'
class YCoord(Horizontal): pass  # '^^'

from pygeode.axis import TAxis
from pygeode.timeaxis import StandardTime



  
# Helper method - check if a coordinate is compatible with a variable
# (on raw variables)
def uses_coord (v, coord):
  if v.atts['ig1'] != coord.atts['ip1']: return False
  if v.atts['ig2'] != coord.atts['ip2']: return False
  if v.atts['ig3'] != coord.atts['ip3']: return False
  return True

# Helper method - apply the coordinate variables ('>>','^^', 'HY') to the variables
def attach_xycoords (varlist):

  # Get the various coords
  xcoords = [XCoord.fromvar(v) for v in varlist if v.name == '>>']
  ycoords = [YCoord.fromvar(v) for v in varlist if v.name == '^^']

  varlist = [v for v in varlist if v.name not in ('>>','^^')]

  # Loop over the other (non-coordinate) variables, and apply whatever
  # coords we can.
  for i,v in enumerate(varlist):
    icoord = v.getaxis(I)
    jcoord = v.getaxis(J)
    for xcoord in xcoords:
      if uses_coord (v, xcoord): icoord = xcoord
    for ycoord in ycoords:
      if uses_coord (v, ycoord): jcoord = ycoord

    varlist[i] = v.replace_axes(I = icoord, J = jcoord)

  return varlist

# Helper method - decode an ip1 value
def decode_ip1 (ip1):
  import numpy as np
  ip1 = np.array(ip1, dtype=np.int)
  ip1 &= 0x0FFFFFF;
  exp = ip1>>20
  ip1 &= 0x00FFFFF
  ip1 = np.asarray(ip1,'float32')
  ip1 *= 10000;
  ip1 /= 10.**exp
  return ip1


# Helper method - decode the vertical coordinate (ip1)
def decode_zcoord (varlist):
  import numpy as np

  # The HY record has no dimensionality of its own, so we can't create
  # particular axes until we match it to particular variables.
  hycoords= [v for v in varlist if v.name == 'HY']

  varlist = [v for v in varlist if v.name not in ('HY',)]

  # Loop over the other (non-coordinate) variables, and apply whatever
  # coords we can.
  for i,v in enumerate(varlist):
    zcoord = v.getaxis(Z)

    # First, decode the ip1 values
    ip1 = v.getaxis(Z).values
    coordtype = ip1[0]>>24
    ip1 = decode_ip1 (ip1)

    # Set a default axis - just in case we don't get anything proper
    zcoord = ZAxis(ip1)

    # Special case - hybrid axis
    # we need the special 'HY' record to instantiate this
    if vertical_type[coordtype] is Hybrid:
     for hycoord in hycoords:
      # Assume we only need to match the 'etiket' parameter to use HY
      if v.atts['etiket'] == hycoord.atts['etiket']:
        p0 = hycoord.atts['ig1'] * 100.
        hy_ip1 = hycoord.getaxis(Z).values[0]
        plid = decode_ip1(hy_ip1) * 100.
        exp = hycoord.atts['ig2'] / 1000.
        etatop = plid / p0

        zvalues = np.array(zcoord.values)
        zvalues[zvalues<etatop] = etatop  # fix bug where z value is slightly less than etatop
        B = ((zvalues - etatop)/(1-etatop))**exp
        A = p0 * (zvalues - B)
        zcoord = Hybrid(zvalues,A,B)
        del p0, hy_ip1, plid, exp, etatop, zvalues, B, A

    # Other vertical coordinates (easy case)
    else:
      zcoord = vertical_type[coordtype](zcoord.values)

    varlist[i] = v.replace_axes(Z = zcoord)

  return varlist


# Helper method - decode time axis
def decode_timeaxis (varlist):
  import numpy as np
  from warnings import warn

  # Make a copy of the list
  varlist = list(varlist)

  # Iterate over each variable
  for i, v in enumerate(varlist):

    dateo = np.array(v.getaxis(T).values)
    dateo = dateo/4 * 5
    # Case 0: degenerate time axis
    if np.all(dateo == 0):
      warn ("degenerate time axis detected", stacklevel=3)
      taxis = TAxis(list(dateo))
    # Case 1: old style
    elif np.all(dateo < 123200000):
      dateo /= 10;  # ignore operation run digit
      hour = dateo % 100; dateo /= 100
      year = 1900 + (dateo % 100); dateo /= 100
      day = dateo % 100; dateo /= 100
      month = dateo
      badmonths = (month < 1) + (month > 12)
      if np.any(badmonths):
        warn("Invalid months detected.  Resetting to 1.", stacklevel=3)
        month[badmonths] = 1
      baddays = (day < 1) + (day > 31)
      if np.any(baddays):
        warn("Invalid days detected.  Resetting to 1.", stacklevel=3)
        day[baddays] = 1
      taxis = StandardTime (year=year, month=month, day=day, units='hours')
    # Case 2: new style
    else:
      dateo -= 123200000;
      dateo *= 4;  # now have # seconds since Jan 1, 1980
      # Convert to hours
      dateo /= 3600.
      taxis = StandardTime (values=dateo, units='hours', startdate={'year':1980})

    varlist[i] = v.replace_axes(T = taxis)

  return varlist

# Helper method - decode the latitudes / longitudes for A/B/G grids
def decode_latlon (varlist):
  import numpy as np
  from pygeode.axis import gausslat, Lat, Lon
  from warnings import warn

  varlist = list(varlist)

  for i, v in enumerate(varlist):

    grtyp = v.atts['grtyp']

    ni = v.atts['ni']
    nj = v.atts['nj']

    if grtyp == 'A':  # lat-lon
      lon = Lon(360./ni * np.arange(ni))
      lat = Lat(-90 + 180./nj * (np.arange(nj)+0.5)) #1/2 gridpoint offset

    elif grtyp == 'B':  # lat-lon, with poles
      lon = Lon(360./(ni-1) * np.arange(ni))
      lat = Lat(-90 + 180./(nj-1) * np.arange(nj) - 90)

    elif grtyp == 'G':  # gaussian grid
      if v.atts['ig1'] != 0:
        warn("can only handle global gaussian grid right now",stacklevel=2)
        continue
      lon = Lon(360./ni * np.arange(ni))
      lat = gausslat(nj)

    # Cheap hack - Z grid on regular lat/lon?
    elif grtyp == 'Z' and v.hasaxis(XCoord) and v.hasaxis(YCoord):
      xcoord = v.getaxis(XCoord)
      ycoord = v.getaxis(YCoord)
      ig1, ig2, ig3, ig4 = map(xcoord.atts.get, ['ig1','ig2','ig3','ig4'])
      # This trick only works for 'E' coordinates
      if xcoord.atts['grtyp'] != 'E': continue
      # Also only works for a particular kind of 'E' coordinates?
      if ig1 != 900 or ig2 != 0 or ig3 != 43200 or ig4 != 43200: continue
      lon = Lon(xcoord.values)
      lat = Lat(ycoord.values)

    # Unhandled grid type
    else: continue

    # Check for hemispheric data?
    if grtyp in ('A','B','G') and v.atts['ig1'] != 0:
      warn("can't handle northern/southern hemispheric grids yet.",stacklevel=2)
      continue

    # North to south?
    if grtyp in ('A','B','G') and v.atts['ig2'] == 1:
      lat = Lat(lat.values[::-1])

    varlist[i] = v.replace_axes (I = lon, xcoord = lon, J = lat, ycoord = lat, ignore_mismatch=True)

  return varlist

# Helper method - remove degenerate axes
def remove_degenerate_axes (varlist):
  varlist = list(varlist)

  for i, v in enumerate(varlist):
    remove_axes = []

    # Remove k dimension?
    if v.hasaxis(K):
      if len(v.getaxis(K)) == 1:
        remove_axes.append(v.whichaxis(K))

    # Remove levels (if only level is 0m above ground)
    if v.hasaxis(Height):
      height = v.getaxis(Height).values
      if len(height) == 1 and height[0] == 0:
        remove_axes.append(v.whichaxis(Height))
    # Remove degenerate time axis?
    if not v.hasaxis(StandardTime):
      taxis = v.getaxis(TAxis)
      if len(taxis) == 1 and taxis[0] == 0:
        remove_axes.append(v.whichaxis(TAxis))

    if len(remove_axes) > 0:
      varlist[i] = v.squeeze(*remove_axes)

  return varlist


# Filter function
# Check for any XCoord/YCoord axes, generate latitude and longitude fields
def add_latlon (varlist):
  from pygeode.var import Var
  from warnings import warn
  import numpy as np

  # Pull out the variables
  varlist = list(varlist)

  # Computed lat/lon fields
  lats = []
  lons = []

  # Keys of computed fields
  computed_keys = []

  for v in varlist:
    if not v.hasaxis(XCoord) or not v.hasaxis(YCoord): continue
    xaxis = v.getaxis(XCoord)
    yaxis = v.getaxis(YCoord)

    grtyp = xaxis.atts['grtyp']
    ig1 = xaxis.atts['ig1']
    ig2 = xaxis.atts['ig2']
    ig3 = xaxis.atts['ig3']
    ig4 = xaxis.atts['ig4']

    # Determine if we already have these fields
    key = (grtyp, ig1, ig2, ig3, ig4, len(xaxis), len(yaxis))
    if key in computed_keys: continue

    # Convert the X/Y coordinates to double precision
    # (for more accurate intermediate calculations)
    xvalues = np.asarray(xaxis.values, dtype='float64')
    yvalues = np.asarray(yaxis.values, dtype='float64')

    # Polar stereographic?
    #TODO: handle new-style encoding of d60/dgrw
    if grtyp in ('N', 'S'):
      d60, dgrw = ig3, ig4
      nhem = 1 if grtyp == 'N' else 2
      lat, lon = llfxy (xvalues, yvalues, d60, dgrw, nhem)
    elif grtyp == 'E':
      # Based on igaxg.f
      lg1 = (ig1<<2) | (ig3&3)
      lg2 = (ig2<<2) | (ig4&3)
      if lg2 > 3600: lg2 -= 7201
      lg3 = (ig3>>2)
      if lg3 < 3559: lg3 += 16384
      lg4 = (ig4>>2)
      xg1 = (lg1 - 3600.) / 40.
      xg2 = (lg3 - 3600.) / 40.
      xg3 = lg2 / 40.
      xg4 = lg4 / 40.
      lat, lon = rotated_ll (xvalues, yvalues, xg1,xg2,xg3,xg4)
    else:
      warn ("can't handle grtyp '%s' yet"%grtyp)
      continue

    # Convert the raw lat/lon values to a Var
    lat = Var([yaxis,xaxis], values=lat, name='latitudes')
    lon = Var([yaxis,xaxis], values=lon, name='longitudes')

    lons.append(lon)
    lats.append(lat)
    computed_keys.append(key)

  # Append these to the list of vars
  varlist.extend(lats)
  varlist.extend(lons)

  return varlist



######################################################################
# You're probably looking for this
######################################################################
def open (filename):
  from pygeode.dataset import Dataset

  # Get the raw (unprocessed) variables
  vars = rawopen(filename)

  # Do some filtering
  vars = attach_xycoords(vars)
  vars = decode_zcoord(vars)
  vars = decode_timeaxis(vars)
  vars = decode_latlon(vars)
  vars = remove_degenerate_axes(vars)
  vars = add_latlon(vars)

  # Convert to a dataset
  dataset = Dataset(vars)

  return dataset




# Miscellaneous functions
# (could be moved to a utility library)

# Polar stereographic conversion.
# Create a 2D array of latitudes and longitudes from the given 1D arrays of coordinates

# based on some random code found at
# http://www.windatlas.ca/scripts/xy_ll.js
# which was probably based on the LLFXY function in the RPN libraries?
def llfxy (x,y,d60,dgrw,nhem):
  import numpy as np

  # Reshape x and y into 2D arrays
  X = x.reshape(1,-1)
  Y = y.reshape(-1,1)

  rdtodg=180. / np.pi
  re=1.866025*6371000/d60
  re2=re*re

  dlon = np.arctan2(Y, X) * rdtodg

  dlon[:,x<0] += 180. * np.sign(Y)

  dlon -= dgrw
  dlon[dlon>180] -= 360.
  dlon[dlon<-180] += 360.

  r2 = X*X + Y*Y
  dlat = (re2-r2) / (re2+r2)
  dlat = np.arcsin(dlat) * rdtodg

  if nhem == 2:
    dlat = -(dlat)
    dlon = -(dlon)


  return dlat, dlon

# Convert points on a map to 3D cartesian points
def ll2p (lat, lon):
  import numpy as np
  # Convert to radians
  lat = lat * (np.pi/180)
  lon = lon * (np.pi/180)
  # Convert from spherical coordinates to cartesian coordinates
  Px = - np.cos(lat) * np.cos(lon)
  Py = - np.cos(lat) * np.sin(lon)
  Pz =   np.sin(lat)

  return Px, Py, Pz

# Convert 3D cartesian points to points on a map
def p2ll (Px, Py, Pz):
  import numpy as np
  lat = np.arcsin(Pz)
  lon = np.arctan2(-Py, -Px)
  # Convert from radians to degrees
  lat *= (180/np.pi)
  lon *= (180/np.pi)

# Compute the cross-product of 2 vectors
def cross_product ((Px, Py, Pz), (Qx, Qy, Qz)):
  # Take the cross-product of the 2 vectors
  Rx = Py*Qz - Qy*Pz
  Ry = Pz*Qx - Qz*Px
  Rz = Px*Qy - Qx*Py
  return Rx, Ry, Rz

# Normalize a vector
def normalize (Px, Py, Pz):
  import numpy as np
  M = np.sqrt(Px**2 + Py**2 + Pz**2)
  Px = Px / M
  Py = Py / M
  Pz = Pz / M
  return (Px, Py, Pz)


# Calculate the lat/lon coordinates for a rotated grid
def rotated_ll (x, y, lat1, lon1, lat2, lon2):

  # Reshape the grid x and y points into 2D arrays
  X = x.reshape(1,-1)
  Y = y.reshape(-1,1)

  # 3D 'X' coordinate
  Px, Py, Pz = ll2p (lat1, lon1)

  # Another point, in the 3D X/Y plane
  Qx, Qy, Qz = ll2p (lat2, lon2)

  # Compute the 3D 'Z' coordinate
  # (take the cross-product of the 2 vectors in the X/Y plane)
  Rx, Ry, Rz = cross_product ((Px, Py, Pz), (Qx, Qy, Qz))
  # Normalize
  Rx, Ry, Rz = normalize (Rx, Ry, Rz)
  

  # Compute the 3D 'Y' coordinate
  # (take the cross-product of the Z and X axes)
  Qx, Qy, Qz = cross_product ((Rx, Ry, Rz), (Px, Py, Pz))

  # Compute the model grid points in 3D space
  Gx, Gy, Gz = ll2p (Y, X)

  # Apply the rotation
  Hx = Px * Gx + Qx * Gy + Rx * Gz
  Hy = Py * Gx + Qy * Gy + Ry * Gz
  Hz = Pz * Gx + Qz * Gy + Rz * Gz

  # Re-normalize (there may be some numerical drift)
  Hx, Hy, Hz = normalize (Hx, Hy, Hz)

  # Convert back to lat/lon coordinates
  return p2ll (Hx, Hy, Hz)

