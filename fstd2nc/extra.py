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

from rpnpy.librmn import librmn
import ctypes as ct
import numpy as np
import numpy.ctypeslib as npc
librmn.compact_float.argtypes = (npc.ndpointer(dtype='int32'), npc.ndpointer(dtype='int32'), npc.ndpointer(dtype='int32'), ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.POINTER(ct.c_double))
librmn.compact_double.argtypes = (npc.ndpointer(dtype='int32'), npc.ndpointer(dtype='int32'), npc.ndpointer(dtype='int32'), ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.POINTER(ct.c_double))
librmn.compact_integer.argtypes = (npc.ndpointer(dtype='int32'), ct.c_void_p, npc.ndpointer(dtype='int32'), ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int)
librmn.ieeepak_.argtypes = (npc.ndpointer(dtype='int32'), npc.ndpointer(dtype='int32'), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int))
librmn.compact_char.argtypes = (npc.ndpointer(dtype='int32'), ct.c_void_p, npc.ndpointer(dtype='int32'), ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int)
librmn.c_armn_uncompress32.argtypes = (npc.ndpointer(dtype='int32'), npc.ndpointer(dtype='int32'), ct.c_int, ct.c_int, ct.c_int, ct.c_int)
librmn.c_armn_compress_setswap.argtypes = (ct.c_int,)
librmn.armn_compress.argtypes = (npc.ndpointer(dtype='int32'),ct.c_int,ct.c_int,ct.c_int,ct.c_int,ct.c_int)
librmn.c_float_unpacker.argtypes = (npc.ndpointer(dtype='int32'),npc.ndpointer(dtype='int32'),npc.ndpointer(dtype='int32'),ct.c_int,ct.POINTER(ct.c_int))

def decode (data):
  '''
  Decodes the raw FSTD data into the final 2D array of values.
  Similar to fstluk, but here the data is already loaded in memory.
  The data should also include the header information at the
  beginning of the array.

  Parameters
  ----------
  data : array
      The encoded data
  '''
  import rpnpy.librmn.all as rmn
  import numpy as np
  from fstd2nc.mixins import dtype_fst2numpy
  data = data.view('>i4').astype('i4')
  ni, nj, nk = data[3]>>8, data[4]>>8, data[5]>>12
  nelm = ni*nj*nk
  datyp = int(data[4]%256) & 191  # Ignore +64 mask.
  nbits = int(data[2]%256)
  dtype = dtype_fst2numpy (datyp, nbits)
  if nbits <= 32:
    work = np.empty(nelm,'int32')
  else:
    work = np.empty(nelm,'int64').view('int32')
  # Strip header
  data = data[20:]
  # Extend data buffer for in-place decompression.
  if datyp in (129,130,134):
    d = np.empty(nelm + 100, dtype='int32')
    d[:len(data)] = data
    data = d
  shape = (nj,ni)
  ni = ct.c_int(ni)
  nj = ct.c_int(nj)
  nk = ct.c_int(nk)
  nelm = ct.c_int(nelm)
  npak = ct.c_int(-nbits)
  nbits = ct.c_int(nbits)
  zero = ct.c_int(0)
  one = ct.c_int(1)
  two = ct.c_int(2)
  tempfloat = ct.c_double(99999.0)

  #print (ni, nj, nk, nbits, datyp, dtype)
  if datyp == 0:
    work = data
  elif datyp == 1:
    if nbits.value <= 32:
      librmn.compact_float(work, data, data[3:], nelm, nbits, 24, 1, 2, 0, ct.byref(tempfloat))
    else:
      raise Exception
      librmn.compact_double(work, data, data[3:], nelm, nbits, 24, 1, 2, 0, ct.byref(tempfloat))
  elif datyp == 2:
    librmn.compact_integer(work, None, data, nelm, nbits, 0, 1, 2)
  elif datyp == 3:
    raise Exception
  elif datyp == 4:
    librmn.compact_integer(work, None, data, nelm, nbits, 0, 1, 4)
  elif datyp == 5:
    librmn.ieeepak_(work, data, ct.byref(nelm), ct.byref(one), ct.byref(npak), ct.byref(zero), ct.byref(two))
  elif datyp == 6:
    librmn.c_float_unpacker(work, data, data[3:], nelm, ct.byref(nbits));
  elif datyp == 7:
    ier = librmn.compact_char(work, None, data, nelm, 8, 0, 1, 10)
    work = work.view('B')[:len(work)] #& 127
  elif datyp == 8:
    raise Exception
  elif datyp == 129:
    librmn.armn_compress(data[5:],ni,nj,nk,nbits,2)
    librmn.compact_float(work,data[1:],data[5:],nelm,nbits.value+64*max(16,nbits.value),0,1,2,0,ct.byref(tempfloat))
  elif datyp == 130:
    #librmn.c_armn_compress_setswap(0)
    librmn.armn_compress(data[1:],ni,nj,nk,nbits,2)
    #librmn.c_armn_compress_setswap(1)
    work[:] = data[1:].astype('>i4').view('>H')[:nelm.value]
  elif datyp == 133:
    librmn.c_armn_uncompress32(work, data[1:], ni, nj, nk, nbits)
  elif datyp == 134:
    librmn.armn_compress(data[4:],ni,nj,nk,nbits,2);
    librmn.c_float_unpacker(work,data[1:],data[4:],nelm,ct.byref(nbits))
  else:
    raise Exception(datyp)
  return work.view(dtype)[:nelm.value].reshape(shape).T


def decode_headers (raw):
  '''
  Decode record headers from a raw byte array.
  Returns a dictionary similar to fstprm, only the entries are
  vectorized over all records instead of 1 record at a time.
  NOTE: This includes deleted records as well.  You can filter them out using
        the 'dltf' flag.

  Parameters
  ----------
  raw : numpy array (dtype='B')
      The raw array of headers to decode.
  '''
  import numpy as np
  raw = raw.view('>i4').astype('uint32').reshape(-1,9,2)
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

  temp8 = np.empty(nrecs, dtype='ubyte')
  temp32 = np.empty(nrecs, dtype='int32')

  np.divmod(raw[:,0,0],2**24, temp8, out['lng'])
  out['lng'] *= 2 # Convert from 8-byte to 4-byte units.
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
  return out


def raw_headers (filename):
  '''
  Extract record headers from the specified file.
  Returns a dictionary similar to fstprm, only the entries are
  vectorized over all records instead of 1 record at a time.
  NOTE: This includes deleted records as well.  You can filter them out using
        the 'dltf' flag.

  Parameters
  ----------
  filename : string
      The file to scan for headers.
  '''
  import numpy as np
  import os
  if not os.path.exists(filename):
    return None
  f = open(filename,'rb')
  # Use same check as maybeFST
  magic = f.read(16)
  if len(magic) < 16 or magic[12:] != b'STDR':
    f.close()
    return None
  # Get the raw (packed) parameters.
  pageaddr = 27; pageno = 0
  raw = []
  pageno_list = []
  recno_list = []
  while pageaddr > 0:
    f.seek(pageaddr*8-8, 0)
    page = np.fromfile(f, '>i4', 8+256*18)
    params = page[8:].reshape(256,9,2)
    nent = page[5]
    raw.append(params[:nent].view('B').flatten())
    recno_list.extend(list(range(nent)))
    pageno_list.extend([pageno]*nent)
    pageaddr = page[4]; pageno += 1
  raw = np.concatenate(raw)
  f.close()
  return raw

# Get the block size for the filesystem where the given file resides.
# May be useful for determining the best way to read chunks of data from
# FSTD files.
def blocksize (filename):
  import subprocess
  try:
    # Use the command-line "stat" program instead of os.stat, because the
    # latter is giving incorrect results for some filesystems.
    # E.g., on a GPFS filesystem with 16MB block size, os.stat is giving a size
    # of 256KB.
    return int(subprocess.check_output(['stat', '-c', '%s', '-f', filename]))
  except OSError:
    return 4096  # Give some default value if 'stat' not available.


# Return the given arrays as a single structured array.
# Makes certain operations easier, such as finding unique combinations.
# Takes a dictionary of arrays, returns a single structured array.
def structured_array (data):
  import numpy as np
  dtype = [(key,value.dtype) for key,value in data.items()]
  n = len(list(data.values())[0])
  out = np.ma.empty(n, dtype=dtype)
  for key in data.keys():
    out[key] = data[key]
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

