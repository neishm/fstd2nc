###############################################################################
# Copyright 2017-2021 - Climate Research Division
#                       Environment and Climate Change Canada
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
Optional helper functions.
These functions provide information about the FST files using alternative
approaches (not the standard librmn/rpnpy functions).
This is solely for performance considerations.
"""


def nrecs (f):
  '''
  Get the total number of records for the specified file.
  Same value as librmn function fstnbr.

  Parameters
  ----------
  f : file object
      The Python file object to scan for headers.
  '''
  import numpy as np
  f.seek(0x34, 0)
  return int(np.fromfile(f, '>i4', 1)[0].astype('uint32'))


def all_params (f, out=None):
  '''
  Extract record headers from the specified file.
  Returns a dictionary similar to fstprm, only the entries are
  vectorized over all records instead of 1 record at a time.
  NOTE: This includes deleted records as well.  You can filter them out using
        the 'dltf' flag.

  Parameters
  ----------
  f : file object
      The Python file object to scan for headers.
  out : dict or pandas DataFrame, optional
      Object to store the headers.  By default a dictionary will be created.
  '''
  import numpy as np
  # Get the raw (packed) parameters.
  pageaddr = 27; pageno = 0
  raw = []
  pageno_list = []
  recno_list = []
  while pageaddr > 0:
    f.seek(pageaddr*8-8, 0)
    page = np.fromfile(f, '>i4', 8+256*18).astype('uint32')
    params = page[8:].reshape(256,9,2)
    nent = page[5]
    raw.append(params[:nent])
    recno_list.extend(list(range(nent)))
    pageno_list.extend([pageno]*nent)
    pageaddr = page[4]; pageno += 1
  raw = np.concatenate(raw)
  # Start unpacking the pieces.
  # Reference structure (from qstdir.h):
  # 0      word deleted:1, select:7, lng:24, addr:32;
  # 1      word deet:24, nbits: 8, ni:   24, gtyp:  8;
  # 2      word nj:24,  datyp: 8, nk:   20, ubc:  12;
  # 3      word npas: 26, pad7: 6, ig4: 24, ig2a:  8;
  # 4      word ig1:  24, ig2b:  8, ig3:  24, ig2c:  8;
  # 5      word etik15:30, pad1:2, etik6a:30, pad2:2;
  # 6      word etikbc:12, typvar:12, pad3:8, nomvar:24, pad4:8;
  # 7      word ip1:28, levtyp:4, ip2:28, pad5:4;
  # 8      word ip3:28, pad6:4, date_stamp:32;
  nrecs = raw.shape[0]
  if out is None:
    out = {}
    out['lng'] = np.empty(nrecs, dtype='int32')
    out['dltf'] = np.empty(nrecs, dtype='ubyte')
    out['swa'] =  np.empty(nrecs, dtype='uint32')
    out['deet'] = np.empty(nrecs, dtype='int32')
    out['nbits'] = np.empty(nrecs, dtype='byte')
    out['grtyp'] = np.empty(nrecs, dtype='|S1')
    out['ni'] = np.empty(nrecs, dtype='int32')
    out['nj'] = np.empty(nrecs, dtype='int32')
    out['datyp'] = np.empty(nrecs, dtype='ubyte')
    out['nk'] = np.empty(nrecs, dtype='int32')
    out['ubc'] = np.empty(nrecs, dtype='uint16')
    out['npas'] = np.empty(nrecs, dtype='int32')
    out['ig1'] = np.empty(nrecs, dtype='int32')
    out['ig2'] = np.empty(nrecs, dtype='int32')
    out['ig3'] = np.empty(nrecs, dtype='int32')
    out['ig4'] = np.empty(nrecs, dtype='int32')
    out['etiket'] = np.empty(nrecs,dtype='|S12')
    out['typvar'] = np.empty(nrecs,dtype='|S2')
    out['nomvar'] = np.empty(nrecs,dtype='|S4')
    out['ip1'] = np.empty(nrecs, dtype='int32')
    out['ip2'] = np.empty(nrecs, dtype='int32')
    out['ip3'] = np.empty(nrecs, dtype='int32')
    out['datev'] = np.empty(nrecs, dtype='int32')
    out['dateo'] = np.empty(nrecs, dtype='int32')
    out['xtra1'] = np.empty(nrecs, dtype='uint32')
    out['xtra2'] = np.empty(nrecs, dtype='uint32')
    out['xtra3'] = np.empty(nrecs, dtype='uint32')
    out['key'] = np.empty(nrecs, dtype='int32')

  temp8 = np.empty(nrecs, dtype='ubyte')
  temp32 = np.empty(nrecs, dtype='int32')

  np.divmod(raw[:,0,0],2**24, temp8, out['lng'])
  np.divmod(temp8,128, out['dltf'], temp8)
  out['swa'][:] = raw[:,0,1]
  np.divmod(raw[:,1,0],256, out['deet'], out['nbits'])
  np.divmod(raw[:,1,1],256, out['ni'], out['grtyp'].view('ubyte'))
  np.divmod(raw[:,2,0],256, out['nj'], out['datyp'])
  np.divmod(raw[:,2,1],4096, out['nk'], out['ubc'])
  out['npas'][:] = raw[:,3,0]//64
  np.divmod(raw[:,3,1],256, out['ig4'], temp32)
  out['ig2'][:] = (temp32 << 16) # ig2a
  np.divmod(raw[:,4,0],256, out['ig1'], temp32)
  out['ig2'] |= (temp32 << 8) # ig2b
  np.divmod(raw[:,4,1],256, out['ig3'], temp32)
  out['ig2'] |= temp32 # ig2c
  etik15 = raw[:,5,0]//4
  etik6a = raw[:,5,1]//4
  et = raw[:,6,0]//256
  etikbc, _typvar = divmod(et, 4096)
  _nomvar = raw[:,6,1]//256
  np.divmod(raw[:,7,0],16, out['ip1'], temp8)
  out['ip2'][:] = raw[:,7,1]//16
  out['ip3'][:] = raw[:,8,0]//16
  date_stamp = raw[:,8,1]
  # Reassemble and decode.
  # (Based on fstd98.c)
  etiket_bytes = np.empty((nrecs,12),dtype='ubyte')
  for i in range(5):
    etiket_bytes[:,i] = ((etik15 >> ((4-i)*6)) & 0x3f) + 32
  for i in range(5,10):
    etiket_bytes[:,i] = ((etik6a >> ((9-i)*6)) & 0x3f) + 32
  etiket_bytes[:,10] = ((etikbc >> 6) & 0x3f) + 32
  etiket_bytes[:,11] = (etikbc & 0x3f) + 32
  out['etiket'][:] = etiket_bytes.flatten().view('|S12')
  nomvar_bytes = np.empty((nrecs,4),dtype='ubyte')
  for i in range(4):
    nomvar_bytes[:,i] = ((_nomvar >> ((3-i)*6)) & 0x3f) + 32
  out['nomvar'][:] = nomvar_bytes.flatten().view('|S4')
  typvar_bytes = np.empty((nrecs,2),dtype='ubyte')
  typvar_bytes[:,0] = ((_typvar >> 6) & 0x3f) + 32
  typvar_bytes[:,1] = ((_typvar & 0x3f)) + 32
  out['typvar'][:] = typvar_bytes.flatten().view('|S2')
  out['datev'][:] = (date_stamp >> 3) * 10 + (date_stamp & 0x7)
  # Note: this dateo calculation is based on my assumption that
  # the raw stamps increase in 5-second intervals.
  # Doing it this way to avoid a gazillion calls to incdat.
  date_stamp = date_stamp - (out['deet']*out['npas'])//5
  out['dateo'][:] = (date_stamp >> 3) * 10 + (date_stamp & 0x7)
  out['xtra1'][:] = out['datev']
  out['xtra2'][:] = 0
  out['xtra3'][:] = 0
  # Calculate the handles (keys)
  # Based on "MAKE_RND_HANDLE" macro in qstdir.h.
  out['key'][:] = ((np.array(recno_list)&0x1FF)<<10) | ((np.array(pageno_list)&0xFFF)<<19)

  return out

# Lightweight test for FST files.
# Uses the same test for fstd98 random files from wkoffit.c (librmn 16.2).
#
# The 'isFST' test from rpnpy calls c_wkoffit, which has a bug when testing
# many small (non-FST) files.  Under certain conditions the file handles are
# not closed properly, which causes the application to run out of file handles
# after testing ~1020 small non-FST files.
def maybeFST(filename):
  with open(filename, 'rb') as f:
    buf = f.read(16)
    if len(buf) < 16: return False
    # Same check as c_wkoffit in librmn
    return buf[12:] == b'STDR'
