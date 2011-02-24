from pygeode.tools import load_lib

lib = load_lib("rpn/librpn.so")


from ctypes import Structure, c_int, c_char, c_longlong, c_ulonglong, c_float, POINTER, byref, c_void_p
Nomvar = c_char * 5
Etiket = c_char * 13
Typvar = c_char * 3

"""
typedef struct {
  Nomvar nomvar;
  Etiket etiket;
  Typvar typvar;
  char grtyp;
  int ip1;
  int ip2;
  int ip3;
  int ig1;
  int ig2;
  int ig3;
  int ig4;
  int nt;
  int nz;
  int nk;
  int nj;
  int ni;
  int deet;
  int npas;
  long long *t;
  int ip1kind;
  float *z;
  unsigned long long *offsets;
} Varinfo_entry;

"""
class Varinfo_entry (Structure):
  _fields_ = [('nomvar', Nomvar), ('etiket', Etiket),
              ('typvar', Typvar), ('grtyp', c_char),
              ('ip2', c_int), ('ip3', c_int),
              ('ig1', c_int), ('ig2', c_int), ('ig3', c_int), ('ig4', c_int),
              ('nt', c_int), ('nz', c_int),
              ('nk', c_int), ('nj', c_int), ('ni', c_int),
              ('deet', c_int), ('npas', c_int), ('t', POINTER(c_longlong)),
              ('ip1kind', c_int), ('z', POINTER(c_float)),
              ('offsets', POINTER(c_ulonglong))]
  def __str__ (self):
    return "%4s  %12s  %2s %1s %4d %4d %4d %4d %4d"%(
           self.nomvar, self.etiket, self.typvar, self.grtyp,
           self.nt, self.nz, self.nk, self.nj, self.ni
    )

MAX_NVARS = 1000
class Varinfo (Structure):
  _fields_ = [('nvars', c_int), ('var', Varinfo_entry*MAX_NVARS)]

Varinfo_p = POINTER(Varinfo)
# hook in a delete function?
def Varinfo_p_del(self):
  self.rpnlib.free_varinfo(self)
Varinfo_p.__del__ = Varinfo_p_del
Varinfo_p.rpnlib = lib

lib.get_varinfo.restype = Varinfo_p
lib.fopen.restype = c_void_p

from pygeode.var import Var
from pygeode.axis import Axis, ZAxis

# Raw axes from the file
class T(Axis): pass  # time axis
class Z(Axis): pass  # levels (over ip1)
class K(Axis): pass  # ?? (nk)
class J(Axis): pass  # latitudes?
class I(Axis): pass  # longitudes?

# More specific sub-axes:
class Height (ZAxis): pass
class Sigma (ZAxis):
  plotorder = -1
  plotscale = 'log'

class Height_wrt_Ground (ZAxis): pass

class Theta (ZAxis): pass

del Axis, ZAxis


#kind    : level type
# 0       : height [m] (metres)
# 1       : sigma  [sg] (0.0->1.0)
# 2       : pressure [mb] (millibars)
# 3       : arbitrary code
# 4       : height [M] (metres) with respect to ground level
# 5       : hybrid coordinates [hy] (0.0->1.0)
# 6       : theta [th]


# Wrap a Varinfo_entry into a pygeode Var
class RPN_Var (Var):
  def __init__ (self, var_, filename):
    from pygeode.var import Var
    t = T(var_.nt)
    z = Z(var_.nz)
    k = K(var_.nk)
    j = J(var_.nj)
    i = I(var_.ni)
    self.var_ = var_
    self.filename = filename
    Var.__init__ (self, [t,z,k,j,i], name=var_.nomvar.strip(), dtype='float32')

  from pygeode.tools import need_full_axes
  @need_full_axes(I,J,K)
  def getview (self, view, pbar):
    from pygeode.tools import point
    import numpy as np
    from ctypes import c_int32, byref

    file = lib.fopen(self.filename, "rb")

    out = np.empty(view.shape, self.dtype)

    # Time, level slicing relative to *whole* variable
    TI = view.integer_indices[0]
    TI = (c_int32*len(TI))(*TI)
    ZI = view.integer_indices[1]
    ZI = (c_int32*len(ZI))(*ZI)

    lib.read_data (file, byref(self.var_), len(TI), TI, len(ZI), ZI, point(out))

    lib.fclose (file)

    pbar.update(100)

    return out
  del need_full_axes

  def getview_old_but_working (self, view, pbar):
#  def getview (self, view, pbar):
    from pygeode.tools import point
    import numpy as np
    from ctypes import c_int32, byref

    file = lib.fopen(self.filename, "rb")
    ksl = view.slices[2]
    jsl = view.slices[3]
    isl = view.slices[4]

    out = np.empty(view.shape, self.dtype)

    # Don't slice over the horizontal axes
    # (we need to read in the whole record)
    bigview = view.only_slice(T, Z)
    viewloop = list(bigview.loop_mem())
    for i,smallview in enumerate(viewloop):
      assert smallview.shape[2:] == bigview.shape[2:], "horizontal axes are too large?!"
      tmp = np.empty(smallview.shape, dtype=self.dtype)
      # Time, level slicing relative to *whole* variable
      TI = smallview.integer_indices[0]
      TI = (c_int32*len(TI))(*TI)
      ZI = smallview.integer_indices[1]
      ZI = (c_int32*len(ZI))(*ZI)
      # Time, level slicing relative to the input view
      tsl, zsl = smallview.map_to(bigview.clip()).slices[0:2]
      lib.read_data (file, byref(self.var_), len(TI), TI, len(ZI), ZI, point(tmp))
      pbar.update(100.*(i+1)/len(viewloop))
      out[tsl,zsl,:,:,:] = tmp[:,:,ksl,jsl,isl]


    lib.fclose (file)
    return out

del Var

# Higher-level wrapper for vars
def wrap (var, stuff):
  from pygeode.axis import Lat, Lon, gausslat, Pres, ZAxis, Hybrid
  from pygeode.timeaxis import StandardTime
  import numpy as np

  # Time axis
  taxis = var.axes[0]
  nt = len(taxis)
  date = np.array(var.var_.t[:nt])
  year, date = divmod(date, 10000000000)
  month, date = divmod(date, 100000000)
  day, date = divmod(date, 1000000)
  hour, date = divmod(date, 10000)
  minute, second = divmod(date, 100)
  taxis = StandardTime(year=year, month=month, day=day, hour=hour, minute=minute, second=second, units = 'days')

  # Horizontal axes
  grtyp = var.var_.grtyp
  ig1 = var.var_.ig1
  ig2 = var.var_.ig2
  ig3 = var.var_.ig3
  ig4 = var.var_.ig4
  xaxis = var.axes[4]
  yaxis = var.axes[3]
  ni = len(xaxis)
  nj = len(yaxis)
  if grtyp == 'A':  # lat-lon
    xaxis = 360./ni * np.arange(ni)
    yaxis = 180./nj * (np.arange(nj)+0.5) #1/2 gridpoint offset

  elif grtyp == 'B':  # lat-lon, with poles
    xaxis = 360./(ni-1) * np.arange(ni)
    yaxis = 180./(nj-1) * np.arange(nj) - 90


  elif grtyp == 'G':  # gaussian grid
    assert ig1 == 0, "can only handle global gaussian grid right now"
    xaxis = 360./ni * np.arange(ni)
    yaxis = gausslat(nj).values

  elif grtyp == 'Z':  # other lat-lon mesh?
    xaxis = [v for v in stuff if v.name == '>>']
    if len(xaxis) != 1: return None
    xaxis = xaxis[0].squeeze().get()
    yaxis = [v for v in stuff if v.name == '^^']
    if len(yaxis) != 1: return None
    yaxis = yaxis[0].squeeze().get()

  else:
    return None  # ignore this variable?

  # Adjust the order of the latitudes?
  if grtyp in ('A','B','G'):
    # Global?
    if ig1 == 0: pass
    # Northern Hemisphere?
    elif ig1 == 1:
      yaxis = (yaxis+90)/2
    # Southern Hemisphere?
    elif ig1 == 2:
      yaxis = (yaxis-90)/2

    # South to north?
    if ig2 == 0: pass
    # North to south?
    elif ig2 == 1: yaxis = yaxis[::-1]

  xaxis = Lon(xaxis)
  yaxis = Lat(yaxis)

  #TODO

  # Vertical axis
  #kind    : level type
  # 0       : height [m] (metres)
  # 1       : sigma  [sg] (0.0->1.0)
  # 2       : pressure [mb] (millibars)
  # 3       : arbitrary code
  # 4       : height [M] (metres) with respect to ground level
  # 5       : hybrid coordinates [hy] (0.0->1.0)
  # 6       : theta [th]
  ip1kind = var.var_.ip1kind
  zaxis = var.var_.z[:var.var_.nz]
  zclass = {0:Height, 1:Sigma, 2:Pres, 3:ZAxis, 4:Height_wrt_Ground, 5:Hybrid, 6:Theta}[ip1kind]
  if zclass is Hybrid:
    hy = [v for v in stuff if v.name == 'HY']
    if len(hy) != 1: return None
    hy = hy[0]
    p0 = hy.var_.ig1 * 100.
    assert hy.var_.nz == 1
    plid = hy.var_.z[0] * 100.
    exp = hy.var_.ig2 / 1000.
    etatop = plid / p0
    zaxis = np.array(zaxis)
    zaxis[np.where(zaxis<etatop)] = etatop  # fix bug where zaxis is slightly less than etatop
    B = ((zaxis - etatop)/(1-etatop))**exp
    A = p0 * (zaxis - B)
    zaxis = Hybrid(zaxis,A,B)
  else:
    zaxis = zclass(zaxis)


  #TODO

  newvar = var.replace_axes (t=taxis, z=zaxis, i=xaxis, j=yaxis)

  remove_axes = []
  # Remove k dimension?
  if var.var_.nk == 1: remove_axes.append(2)
  # Remove levels (if only level is 0m above ground)
  if type(zaxis) == Height and len(zaxis) == 1 and zaxis.values == [0]:
    remove_axes.append(1)

  if len(remove_axes) > 0: newvar = newvar.squeeze(*remove_axes)

  return newvar




def open (filename):
  from pygeode.dataset import Dataset
  from os.path import exists

  assert exists(filename), "Can't find %s"%filename

  rawvars = []
  vinf_p = lib.get_varinfo(filename)
  vinf = vinf_p[0]
  for var_ in vinf.var[:vinf.nvars]:
    rawvars.append(RPN_Var(var_,filename))
  vars = [wrap(v,rawvars) for v in rawvars]
  vars = filter(None, vars)  # remove the unhandled vars
  return Dataset(vars, print_warnings=False)


