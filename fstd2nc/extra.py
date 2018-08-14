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
Optional helper functions.
Note: These rely on some assumptions about the internal structures of librmn,
      and may fail for future libray verisons.
      These have been tested for librmn 15.2 and 16.2.
"""

from ctypes import Structure, POINTER, c_void_p, c_uint32, c_int32, c_int, c_uint, c_byte, c_char_p
from rpnpy.librmn import librmn


# From fnom.h
MAXFILES = 1024
class attributs(Structure):
  _fields_ = [
     ('stream',c_uint,1), ('std',c_uint,1), ('burp',c_uint,1), ('rnd',c_uint,1), ('wa',c_uint,1), ('ftn',c_uint,1),
     ('unf',c_uint,1), ('read_only',c_uint,1), ('old',c_uint,1), ('scratch',c_uint,1), ('notpaged',c_uint,1),
     ('pipe',c_uint,1), ('write_mode',c_uint,1), ('remote',c_uint,1), ('padding',c_uint,18),
  ]

class general_file_info(Structure):
  _fields_ = [
    ('file_name',c_char_p),            # complete file name
    ('subname',c_char_p),              # sub file name for cmcarc files
    ('file_type',c_char_p),            # file type and options
    ('iun',c_int),                    # fnom unit number
    ('fd',c_int),                     # file descriptor
    ('file_size',c_int),              # file size in words
    ('eff_file_size',c_int),          # effective file size in words
    ('lrec',c_int),                   # record length when appliable
    ('open_flag',c_int),              # open/close flag
    ('attr',attributs),
  ]


Fnom_General_File_Desc_Table = (general_file_info*MAXFILES).in_dll(librmn,'Fnom_General_File_Desc_Table')

# From rpnmacros.h
word = c_uint32

# From qstdir.h
MAX_XDF_FILES = 1024
ENTRIES_PER_PAGE = 256
MAX_DIR_PAGES = 1024
MAX_PRIMARY_LNG = 16
MAX_SECONDARY_LNG = 8
max_dir_keys = word*MAX_PRIMARY_LNG
max_info_keys = word*MAX_SECONDARY_LNG

class stdf_dir_keys(Structure): pass  #defined further below

class xdf_dir_page(Structure):
  _fields_ = [
        ('lng',word,24), ('idtyp',word,8), ('addr',word,32), # XDF record header
        ('reserved1',word,32), ('reserved2',word,32),
        ('nxt_addr',word,32), ('nent',word,32),
        ('chksum',word,32), ('reserved3',word,32),
        ('entry',stdf_dir_keys*ENTRIES_PER_PAGE),
  ]
# idtyp:     id type (usualy 0)
# lng:       header length (in 64 bit units)
# addr:      address of directory page (origin 1, 64 bit units)
# reserved1: idrep (4 ascii char 'DIR0')
# reserved2: reserved (0)
# nxt_addr:  address of next directory page (origin 1, 64 bit units)
# nent:      number of entries in page
# chksum:    checksum (not valid when in core)
# page_no, record_no, file_index: handle templage
# entry:     (real allocated dimension will be ENTRIES_PER_PAGE * primary_len)

class full_dir_page(Structure):
  pass
page_ptr = POINTER(full_dir_page)
full_dir_page._fields_ = [
        ('next_page',page_ptr),
        ('prev_page',page_ptr),
        ('modified',c_int),
        ('true_file_index',c_int),
        ('dir',xdf_dir_page),
]

class file_record(Structure):
  _fields_ = [
        ('lng',word,24), ('idtyp',word,8), ('addr',word,32),        #XDF record header
        ('data',word*2),                        # primary keys, info keys, data
  ]


stdf_dir_keys._fields_ = [
        ('lng',word,24), ('select',word,7), ('deleted',word,1), ('addr',word,32),
        ('nbits',word,8), ('deet',word,24), ('gtyp',word,8), ('ni',word,24),
        ('datyp',word,8), ('nj',word,24), ('ubc',word,12), ('nk',word,20),
        ('pad7',word,6), ('npas',word,26), ('ig2a',word,8), ('ig4',word,24),
        ('ig2b',word,8), ('ig1',word,24), ('ig2c',word,8), ('ig3',word,24),
        ('pad1',word,2), ('etik15',word,30), ('pad2',word,2), ('etik6a',word,30),
        ('pad3',word,8), ('typvar',word,12), ('etikbc',word,12), ('pad4',word,8), ('nomvar',word,24),
        ('levtyp',word,4), ('ip1',word,28), ('pad5',word,4), ('ip2',word,28),
        ('pad6',word,4), ('ip3',word,28), ('date_stamp',word,32),
  ]


class key_descriptor(Structure):
  _fields_ = [
        ('ncle',word,32), ('reserved',word,8), ('tcle',word,6), ('lcle',word,5), ('bit1',word,13),
  ]


class file_header(Structure):
  _fields_ = [
        ('lng',word,24), ('idtyp',word,8), ('addr',word,32),  # standard XDF record header
        ('vrsn',word),     ('sign',word),               #char[4]
        ('fsiz',word,32), ('nrwr',word,32),
        ('nxtn',word,32), ('nbd',word,32),
        ('plst',word,32), ('nbig',word,32),
        ('lprm',word,16), ('nprm',word,16), ('laux',word,16), ('naux',word,16),
        ('neff',word,32), ('nrec',word,32),
        ('rwflg',word,32), ('reserved',word,32),
        ('keys',key_descriptor*1024),
  ]
# idtyp:     id type (usualy 0)
# lng:       header length (in 64 bit units)
# addr:      address (exception: 0 for a file header)
# vrsn:      XDF version
# sign:      application signature
# fsiz:      file size (in 64 bit units)
# nrwr:      number of rewrites
# nxtn:      number of extensions
# nbd:       number of directory pages
# plst:      address of last directory page (origin 1, 64 bit units)
# nbig:      size of biggest record
# nprm:      number of primary keys
# lprm:      length of primary keys (in 64 bit units)
# naux:      number of auxiliary keys
# laux:      length of auxiliary keys
# neff:      number of erasures
# nrec:      number of valid records
# rwflg:     read/write flag
# reserved:  reserved
# keys:      key descriptor table


class file_table_entry(Structure):
  _fields_ = [
        ('dir_page',page_ptr*MAX_DIR_PAGES), # pointer to directory pages
        ('cur_dir_page',page_ptr),           # pointer to current directory page
        ('build_primary',c_void_p),       # pointer to primary key building function
        ('build_info',c_void_p),           # pointer to info building function
        ('scan_file',c_void_p),            # pointer to file scan function
        ('file_filter',c_void_p),          # pointer to record filter function
        ('cur_entry',POINTER(word)),              # pointer to current directory entry
        ('header',POINTER(file_header)),          # pointer to file header
        ('nxtadr',c_int32),                # next write address (in word units)
        ('primary_len',c_int),
        # length in 64 bit units of primary keys (including 64 bit header)
        ('info_len',c_int),                 #length in 64 bit units of info keys
        ('link',c_int),                     # file index to next linked file,-1 if none
        ('cur_info',POINTER(general_file_info)),
                                      # pointer to current general file desc entry
        ('iun',c_int),                      # FORTRAN unit number, -1 if not open, 0 if C file
        ('file_index',c_int),               # index into file table, -1 if not open
        ('modified',c_int),                 # modified flag
        ('npages',c_int),                   # number of allocated directory pages
        ('nrecords',c_int),                 # number of records in file
        ('cur_pageno',c_int),               # current page number
        ('page_record',c_int),              # record number within current page
        ('page_nrecords',c_int),            # number of records in current page
        ('file_version',c_int),             # version number
        ('valid_target',c_int),             # last search target valid flag
        ('xdf_seq',c_int),                  # file is sequential xdf
        ('valid_pos',c_int),                # last position valid flag (seq only)
        ('cur_addr',c_int),                 # current address (WA, sequential xdf)
        ('seq_bof',c_int),                  # address (WA) of first record (seq xdf)
        ('fstd_vintage_89',c_int),          # old standard file flag
        ('head_keys',max_dir_keys),       # header & primary keys for last record
        ('info_keys',max_info_keys),      # info for last read/written record
        ('cur_keys',max_dir_keys),        # keys for current operation
        ('target',max_dir_keys),          # current search target
        ('srch_mask',max_dir_keys),       # permanent search mask for this file
        ('cur_mask',max_dir_keys),        # current search mask for this file
  ]
file_table_entry_ptr = POINTER(file_table_entry)

librmn.file_index.argtypes = (c_int,)
librmn.file_index.restype = c_int
file_table = (file_table_entry_ptr*MAX_XDF_FILES).in_dll(librmn,'file_table')


def all_params (funit, out=None):
  '''
  Extract parameters for *all* records.
  Returns a dictionary similar to fstprm, only the entries are
  vectorized over all records instead of 1 record at a time.
  NOTE: This includes deleted records as well.  You can filter them out using
        the 'dltf' flag.
  '''
  from ctypes import cast
  import numpy as np
  # Get the raw (packed) parameters.
  index = librmn.file_index(funit)
  raw = []
  file_index_list = []
  pageno_list = []
  recno_list = []
  while index >= 0:
    f = file_table[index].contents
    for pageno in range(f.npages):
      page = f.dir_page[pageno].contents
      params = cast(page.dir.entry,POINTER(c_uint32))
      params = np.ctypeslib.as_array(params,shape=(ENTRIES_PER_PAGE,9,2))
      nent = page.dir.nent
      raw.append(params[:nent])
      recno_list.extend(list(range(nent)))
      pageno_list.extend([pageno]*nent)
      file_index_list.extend([index]*nent)
    index = f.link
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
  out['key'][:] = (np.array(file_index_list)&0x3FF) | ((np.array(recno_list)&0x1FF)<<10) | ((np.array(pageno_list)&0xFFF)<<19)

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
