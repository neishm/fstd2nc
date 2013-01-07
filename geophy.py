from ctypes import Structure, c_int, c_longlong, c_char, c_void_p, byref
from pygeode.tools import load_lib
lib = load_lib("plugins/rpn/libgeophy.so")
lib.fopen.restype = c_void_p

class Geophy_File:
  def __init__ (self, filename):
    self.lib = lib  # So we can close the file upon cleanup
    f = lib.fopen(filename, "rb")
    assert f, "unable to open file '%s'"%filename
    self.file = f
  def __del__ (self):
    if hasattr(self,'file'):
      self.lib.fclose(self.file)

from pygeode.rpn.bmf import IAxis, JAxis, KAxis
from pygeode.var import Var
class Geophy_Var(Var):
  def __init__ (self, file, geonm1, geonm2, geonm5, geopar1, geopar2, geopar3, ni, nj, base_offset):
    from pygeode.var import Var
    self.file = file
    iaxis = IAxis(ni)
    jaxis = JAxis(nj)
    kaxis = KAxis(int(geopar3))
    name = geonm2.value.strip()
    etiket = geonm1.value.strip()
    interp_type = geonm5.value.strip()
    if name == '00': name = etiket
    self.offset = base_offset + (int(geopar1)-1)*4
    Var.__init__ (self, axes=[kaxis,jaxis,iaxis], name=name, atts={'long_name':etiket, 'interpolation_type':interp_type}, dtype='float32')

  from pygeode.tools import need_full_axes
  @need_full_axes(IAxis,JAxis,KAxis)
  def getview (self, view, pbar):
    import numpy as np
    from pygeode.tools import point
    out = np.empty(view.shape, dtype=self.dtype)
    lib.read_data(self.file.file, self.offset, self.size, point(out))
    pbar.update(100)
    return out
  del need_full_axes

del Var

def open (filename, ni, nj):
  from pygeode.dataset import Dataset

  f = Geophy_File(filename)
  p_bgeo_top = c_int()
  p_bgeo_siz = c_int()
  ret = lib.get_file_info (f.file, byref(p_bgeo_top), byref(p_bgeo_siz))
  assert ret > 0, "problem opening %s"%filename
  p_bgeo_top = p_bgeo_top.value
  p_bgeo_siz = p_bgeo_siz.value
  base_offset = 16  # account for p_bgeo_top, p_bgeo_siz header
  base_offset += 4 + p_bgeo_top*16*3 + 4  # account for geonm header
  base_offset += 4 + p_bgeo_top*4*3 + 4  # account for geopar header
  base_offset += 4  # account for leading size of geobus
#  lib.test_base_offset(f.file, base_offset-4, p_bgeo_siz)

  geonm1 = ((c_char*16)*p_bgeo_top)()
  geonm2 = ((c_char*16)*p_bgeo_top)()
  geonm5 = ((c_char*16)*p_bgeo_top)()
  geopar1 = (c_int*p_bgeo_top)()
  geopar2 = (c_int*p_bgeo_top)()
  geopar3 = (c_int*p_bgeo_top)()
  ret = lib.read_metadata (f.file, p_bgeo_top, geonm1, geonm2, geonm5, geopar1, geopar2, geopar3)
  assert ret > 0, "problem reading metadata"

  varlist = []
  for n1, n2, n5, p1, p2, p3 in zip(geonm1, geonm2, geonm5, geopar1, geopar2, geopar3):
    varlist.append(Geophy_Var(f, n1, n2, n5, p1, p2, p3, ni, nj, base_offset))

  return Dataset(varlist)
