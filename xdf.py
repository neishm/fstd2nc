from ctypes import Structure, POINTER, c_void_p, c_uint32, c_int32, c_int

# From rpnmacros.h
word = c_uint32

# From qstdir.h
MAX_DIR_PAGES = 1024
MAX_PRIMARY_LNG = 16
MAX_SECONDARY_LNG = 8
max_dir_keys = word*MAX_PRIMARY_LNG
max_info_keys = word*MAX_SECONDARY_LNG

class xdf_dir_page(Structure):
  _fields_ = [
        ('idtyp',word,8), ('lng',word,24), ('addr',word,32), # XDF record header
        ('reserved1',word,32), ('reserved2',word,32),
        ('nxt_addr',word,32), ('nent',word,32),
        ('chksum',word,32), ('reserved3',word,32),
        ('entry',word*2),
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
full_dir_page._fields__ = [
        ('next_page',page_ptr),
        ('prev_page',page_ptr),
        ('modified',c_int),
        ('true_file_index',c_int),
        ('dir',xdf_dir_page),
]

class file_table_entry(Structure):
  _fields_ = [
        ('dir_page',page_ptr*MAX_DIR_PAGES), # pointer to directory pages
        ('cur_dir_page',page_ptr),           # pointer to current directory page
        ('build_primary',c_void_p),       # pointer to primary key building function
        ('build_info',c_void_p),           # pointer to info building function
        ('scan_file',c_void_p),            # pointer to file scan function
        ('file_filter',c_void_p),          # pointer to record filter function
        ('cur_entry',word),              # pointer to current directory entry
        ('header',c_void_p),          # pointer to file header
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

import fstd2nc_deps
from rpnpy.librmn.fstd98 import fstopenall, fstcloseall
from rpnpy.librmn.const import FST_RO
from rpnpy.librmn import librmn
from glob import glob

iun = fstopenall(glob("/wrk6/neish/mn129/model/2010010?00_024"),FST_RO)

p = POINTER(file_table_entry).in_dll(librmn,'file_table')
print p
print dir(p)
print bool(p)
print p[1].file_version

fstcloseall(iun)
