# Interface for reading "BMF" files
# (binary files used internally in the GEM model)

from ctypes import Structure, c_int, c_longlong, c_char, c_void_p, byref
from pygeode.tools import load_lib
lib = load_lib("plugins/rpn/libbmf.so")
lib.fopen.restype = c_void_p

"""
typedef struct {
  char nom[5];
  int ni, istart, iend;
  int nj, jstart, jend;
  int nk, kstart, kend;
  int time1, time2;
  int hgrid, vgrid;
  int dtyp, scat;
  int ndata;
  long long data_start;
} BMF_Header;
"""

class BMF_Header (Structure):
  _fields_ = [('nom',c_char*5),
              ('ni',c_int),('istart',c_int),('iend',c_int),
              ('nj',c_int),('jstart',c_int),('jend',c_int),
              ('nk',c_int),('kstart',c_int),('kend',c_int),
              ('time1',c_int),('time2',c_int),
              ('hgrid',c_int),('vgrid',c_int),
              ('dtyp',c_int),('scat',c_int),
              ('ndata',c_int),
              ('data_start', c_longlong)
             ]

class BMF_File:
  def __init__ (self, filename):
    self.lib = lib  # So we can close the file upon cleanup
    f = lib.fopen(filename, "rb")
    assert f, "unable to open file '%s'"%filename
    self.file = f
    self.nrec = lib.num_records(self.file)
  def __del__ (self):
    if hasattr(self,'file'):
      self.lib.fclose(self.file)

from pygeode.axis import XAxis, YAxis, ZAxis
class IAxis(XAxis): name = 'i'
class JAxis(YAxis): name = 'j'
class KAxis(ZAxis): name = 'k'
del XAxis, YAxis, ZAxis

from pygeode.var import Var
# A var from a single segment
class BMF_Var(Var):
  def __init__ (self, file, header):
    from pygeode.var import Var
    import numpy as np
    self.file = file
    self.header = header
    self.time = "%08d%06d"%(self.header.time1,self.header.time2)
    name = str(header.nom).strip()
    ivals = np.arange(header.istart, header.iend+1)
    jvals = np.arange(header.jstart, header.jend+1)
    kvals = np.arange(header.kstart, header.kend+1)

    axes = [KAxis(kvals), JAxis(jvals), IAxis(ivals)]
    dtyp = header.dtyp
    dtyp = ('int','float')[dtyp%10]+str(dtyp//10*8)
    Var.__init__ (self, axes=axes, name=name, dtype=dtyp)

  from pygeode.tools import need_full_axes
  @need_full_axes(IAxis,JAxis,KAxis)
  def getview (self, view, pbar):
    import numpy as np
    out = np.empty(view.shape, dtype=self.dtype)
    lib.read_data(self.file.file, byref(self.header), out.ctypes.data)
    pbar.update(100)
    return out
  del need_full_axes
del Var


# Low-level interface - no axis manipulation done
def rawopen (filename):
  # Load test file
  f = BMF_File(filename)

  # Gather headers
  headers = (BMF_Header * f.nrec)()
  lib.read_all_headers(f.file, headers)

  # Convert each header to a variable
  varpieces = [BMF_Var(f,h) for h in headers]

  # Get all variable names (in original file order)
  varnames = []
  for var in varpieces:
    if var.name not in varnames:
      varnames.append(var.name)

  # Organize the pieces by name
  old_varpieces = varpieces
  varpieces = {}
  for var in old_varpieces:
    if var.name not in varpieces: varpieces[var.name] = []
    varpieces[var.name].append(var)
  del old_varpieces

  # Concatenate variable pieces together
  from pygeode.concat import concat
  import numpy as np
  varlist = []
  for name in varnames:
    pieces = varpieces[name]
    i_sections = sorted(set(var.header.istart for var in pieces))
    j_sections = sorted(set(var.header.jstart for var in pieces))
    k_sections = sorted(set(var.header.kstart for var in pieces))
    ni = len(i_sections)
    nj = len(j_sections)
    nk = len(k_sections)

    # Collect the pieces into a grid
    grid = np.empty([nk,nj,ni], dtype=object)
    for var in pieces:
      i = i_sections.index(var.header.istart)
      j = j_sections.index(var.header.jstart)
      k = k_sections.index(var.header.kstart)
      grid[k,j,i] = var
    # Concatenate
    imerge = np.empty([nk,nj],dtype=object)
    for k in range(nk):
      for j in range(nj):
        imerge[k,j] = concat(list(grid[k,j,:]), iaxis='iaxis')
    jmerge = np.empty([nk],dtype=object)
    for k in range(nk):
      jmerge[k] = concat(list(imerge[k,:]), iaxis='jaxis')
    kmerge = concat(list(jmerge), iaxis='kaxis')
    varlist.append(kmerge)

  return varlist


def open (filename):
  from pygeode.dataset import Dataset
  #TODO: further annotations?
  return Dataset(rawopen(filename))


# Open all tiles with the given prefix
def open_tiles (prefix):
  from glob import glob
  from re import match
  from pygeode.concat import concat
  from pygeode.dataset import Dataset
  import numpy as np
  files = glob(prefix+"-??-??")
  assert len(files) > 0, "no matches found for prefix '%s'"%prefix
  xtiles = [match(prefix+"-(..)-..$", file).group(1) for file in files]
  xtiles = sorted(set(int(xtile) for xtile in xtiles))
  ytiles = [match(prefix+"-..-(..)$", file).group(1) for file in files]
  ytiles = sorted(set(int(ytile) for ytile in ytiles))

  grid = np.empty([len(xtiles),len(ytiles)], dtype=object)
  grid[()] = None

  for xtile in xtiles:
    for ytile in ytiles:
      print (xtile, ytile)
      data = open("%s-%02d-%02d"%(prefix,xtile,ytile))
      grid[xtile,ytile] = data

  varlist = []
  for var in list(grid[0,0].vars):

    columns = []
    for xtile in xtiles:
      rows = []
      for ytile in ytiles:
        data = grid[xtile,ytile][var.name]
        ivals = data.iaxis.values
        jvals = data.jaxis.values
        # Adjust the coords?
        for previous_xtile in range(xtile):
          ivals = ivals + len(grid[previous_xtile,ytile][var.name].iaxis)
        for previous_ytile in range(ytile):
          jvals = jvals + len(grid[xtile,previous_ytile][var.name].jaxis)
        data = data.replace_axes(iaxis=IAxis(ivals), jaxis=JAxis(jvals))
        rows.append(data)
      columns.append(concat(rows))
    data = concat(columns)
    varlist.append(data)

  return Dataset(varlist)


