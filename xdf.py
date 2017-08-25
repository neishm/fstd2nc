import numpy as np
from ctypes import Structure, POINTER, addressof, cast, sizeof, c_void_p, c_uint32, c_int32, c_int, c_byte

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
        ('cur_info',c_void_p),
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

import fstd2nc_deps
from rpnpy.librmn.fstd98 import fstopenall, fstcloseall
from rpnpy.librmn.const import FST_RO
from rpnpy.librmn import librmn
from glob import glob

librmn.file_index.argtypes = (c_int,)
librmn.file_index.restype = c_int

file_table = (file_table_entry_ptr*MAX_XDF_FILES).in_dll(librmn,'file_table')

#iun = fstopenall(glob("/wrk6/neish/mn129/model/2010010?00_024"),FST_RO)
from rpnpy.librmn.base import fnom, fclos
from rpnpy.librmn.fstd98 import fstouv, fstfrm
iun = fnom("/wrk6/neish/mn129/model/2010010100_024",FST_RO)
istat = fstouv(iun,FST_RO)
print '??', iun

from rpnpy.librmn.fstd98 import fstinl
#keys = fstinl(iun)

index = librmn.file_index(iun)

def print_structure(s):
  for f in s._fields_:
    n = f[0]
    print '%s = %s'%(n,getattr(s,n))

f = file_table[index].contents
print_structure(f)
print '---'
print_structure(f.header.contents)
p = f.dir_page[0].contents
print '---'
print_structure(p)
print '---'
print_structure(p.dir)
print '---'
print_structure(p.dir.entry[0])
print p.dir.entry

#np.ctypeslib.as_array(p.dir.entry)
#data = (c_byte*(sizeof(stdf_dir_keys)*ENTRIES_PER_PAGE)).from_address(addressof(p.dir.entry[0]))

# Convert an array of entries to parameters.
# Reference for structure:
# 0      word deleted:1, select:7, lng:24, addr:32;
# 1      word deet:24, nbits: 8, ni:   24, gtyp:  8;
# 2      word nj:24,  datyp: 8, nk:   20, ubc:  12;
# 3      word npas: 26, pad7: 6, ig4: 24, ig2a:  8;
# 4      word ig1:  24, ig2b:  8, ig3:  24, ig2c:  8;
# 5      word etik15:30, pad1:2, etik6a:30, pad2:2;
# 6      word etikbc:12, typvar:12, pad3:8, nomvar:24, pad4:8;
# 7      word ip1:28, levtyp:4, ip2:28, pad5:4;
# 8      word ip3:28, pad6:4, date_stamp:32;
#TODO
data = cast(p.dir.entry,POINTER(c_uint32))
data = np.ctypeslib.as_array(data,shape=(ENTRIES_PER_PAGE,9,2))[:p.dir.nent]


#fstcloseall(iun)
fstfrm(iun)
fclos(iun)
