###############################################################################
# Copyright 2017-2023 - Climate Research Division
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

from fstd2nc.stdout import _, info, warn, error
from fstd2nc.mixins import BufferBase

# Override dtype_numpy2fst to work around issue with numpy >= 1.20.
# TODO: remove this once public release of Pyton-RPN is updated.
def dtype_numpy2fst (dtype, compress=False, missing=False):
  import rpnpy.librmn.all as rmn
  from rpnpy.librmn.const import FST_DATYP2NUMPY_LIST, FST_DATYP2NUMPY_LIST64
  assert compress is False and missing is False
  for (i, d) in FST_DATYP2NUMPY_LIST.items():
    if dtype == d: return i
  for (i, d) in FST_DATYP2NUMPY_LIST64.items():
    if dtype == d: return i
  return 0

#################################################
# Mixin for adding a compatibility layer to the netCDF output file, so it can
# also function as a valid FSTD file.


class FSTD_Compat (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    import argparse
    super(FSTD_Compat,cls)._cmdline_args(parser)
    parser.add_argument('--fstd-compat', action='store_true', help=_('Adds a compatibility layer to the netCDF output file, so it can also function as a valid FSTD file.  EXPERIMENTAL.'))

  def __init__ (self, *args, **kwargs):
    """
    fstd_compat : bool, optional
        Adds a compatibility layer to the netCDF output file, so it can
        also function as a valid FSTD file.  EXPERIMENTAL.
    """

    # Check if compatibility interface should be activated.
    fstd_compat = kwargs.pop('fstd_compat', False)
    if fstd_compat:
      self.to_netcdf = self._to_netcdf_compat
      self._fstluk = self._fstluk_compat
    super(FSTD_Compat,self).__init__(*args,**kwargs)

  # Simulate 'fstluk' call (returning data + parameters).
  # This will also keep track of which FST records were used in the
  # conversion.  The header information for these records will be added to the
  # output file.
  def _fstluk_compat (self, rec_id, dtype=None, rank=None, dataArray=None):

    # Keep track of which records were used.
    if not hasattr(self,'_used_rec_ids'):
      self._used_rec_ids = []
    self._used_rec_ids.append(rec_id)

    prm = super(FSTD_Compat,self)._fstluk(rec_id)
    if dtype is not None:
      prm['d'] = prm['d'].astype(dtype)

    return prm

  def _to_netcdf_compat (self, filename, nc_format='NETCDF4', global_metadata=None, zlib=False, compression=4, progress=False, turbo=False):
    """
    Write the records to a netCDF file.
    Requires the netCDF4 package.
    """
    from fstd2nc.mixins import _var_type, _ProgressBar, _FakeBar
    from netCDF4 import Dataset
    import numpy as np
    import os
    import rpnpy.librmn.all as rmn

    # This only works with an uncompressed netCDF4 file.
    nc_format = 'NETCDF4'
    zlib = False

    # Get a minimal netCDF4 header, to be used later.
    Dataset(filename, "w", format='NETCDF4').close() # Get a netCDF4 header.
    with open(filename,'rb') as f:
      nc_header = np.fromfile(f,'B')
    os.remove(filename)  # Got the netCDF4, now start over.
    # Start the file with an FSTD file header followed by a netCDF4 header.
    # The FSTD header must be at the very beginning of the file in order for
    # librmn to recognize the file, but the netCDF4 header is more flexible.
    iun = rmn.fstopenall(filename, rmn.FST_RW)
    rmn.fstcloseall(iun)  # Write FSTD header.
    # Write the netCDF4 header in a very particular location so it can be
    # detected.
    # The netCDF4 library will check the beginning of the file, but also check
    # byte positions at powers of two starting at 512, 1024, 2048, etc.
    # Here we put it at offset 32768 (hex 0x8000), which is just after the FSTD
    # header and the first FSTD directory "page".
    with open(filename,'r+b') as f:
      f.seek(0x8000,0)
      nc_header.tofile(f)

    # Now, open the file as netCDF4 format.
    f = Dataset(filename, "r+")

    # Apply global metadata (from config files and global_metadata argument).
    if 'global' in getattr(self,'_metadata',{}):
      f.setncatts(self._metadata['global'])
    if global_metadata is not None:
      f.setncatts(global_metadata)

    # Collect all the records that will be read/written.
    # List of (key,recshape,ncvar,ncind).
    # Note: derived variables (with values stored in memory) will be written
    # immediately, bypassing this list.
    io = []

    self._makevars()

    # Define the dimensions.
    for axis in self._iter_axes():
      # Special case: make the time dimension unlimited.
      if axis.name == 'time' and self._time_unlimited:
        f.createDimension(axis.name, None)
      else:
        f.createDimension(axis.name, len(axis))

    # When writing data to netCDF4, need to make sure it's aligned properly.
    # Also need to get the file offset when the data is written.
    def write_data (v, ind, array, npad=[0]):
      from os.path import getsize
      import rpnpy.librmn.all as rmn
      # Align the data on 8-byte boundaries, which is what FSTD expects.
      f.sync()
      if getsize(filename) % 8 != 0:
        # Add dummy global attributes as a way to pad out the file.
        # The chosen size for the attributes is from trial-and-error. It seems
        # to get the right padding after a few iterations.
        # Ideally it would be nice to close the file, add a few bytes of
        # padding, update the netCDF4 header and re-open with the required
        # alignment, but there's a bug when re-opening a netCDF4 file with
        # big-endian variables (https://github.com/Unidata/netcdf-c/issues/1802).
        # The current approach of writing attributes wastes more space (a few
        # dozen kilobytes?) but it seems to do the job.
        #f.setncattr('_pad%05d'%npad[0], np.zeros(4096+(8-alignment),dtype='B'))
        #f.setncattr('_pad%05d'%npad[0], np.zeros(0+(npad[0]%7),dtype='B'))
        padding = 1
        while True:
          f.setncattr('_pad%05d'%npad[0], np.zeros(padding,dtype='B'))
          f.sync()
          if getsize(filename) % 8 == 0: break
          f.delncattr('_pad%05d'%npad[0])
          padding += 127
        npad[0] += 1
      # Next, write the data
      address = getsize(filename)
      v[ind] = array
      f.sync()
      # Finally, locate *where* the data was written to file.
      array = array.flatten()
      with open(filename,'rb') as f2:
        f2.seek(address,0)
        test = np.fromfile(f2,v.dtype,array.size)
        if not np.all(test==array):
          address = getsize(filename) - array.size*array.dtype.itemsize
          f2.seek(address,0)
          test = np.fromfile(f2,v.dtype,array.size)
          if not np.all(test==array):
            return None
      # Determine which FST datyp would decode this data.
      datyp = dtype_numpy2fst(v.dtype.newbyteorder('='), compress=False, missing=False)
      # Floating-point must be IEEE, which is what netCDF4 would have written.
      if datyp == 1:
        datyp = 5
      nbits = v.dtype.itemsize*8
      # <Address of chunk>, <length>, <RPN datyp>, <nbits>
      return address, array.size*v.dtype.itemsize, datyp, nbits

    # Addresses of arrays not associated with FSTD records.
    # Keep track of these in case they *can* be associated with records later.
    # E.g., lat/lon coordinate arrays might exactly match data from ^^ or >>.
    direct_addresses = {}

    # Generate the variable structures.
    for var in self._iter_objects():

      # Write the variable.
      # Easy case: already have the data.
      if hasattr(var,'array'):
        # Tell netCDF4 to use big-endian encoding of data, where applicable.
        # FSTD requires big-endian encoding (and netCDF4 can work with either).
        if var.array.dtype.itemsize > 1:
          v = f.createVariable(var.name, datatype=var.array.dtype.newbyteorder('>'), endian='big', dimensions=var.dims, zlib=zlib, complevel=compression)
        else:
          v = f.createVariable(var.name, datatype=var.array.dtype, dimensions=var.dims, zlib=zlib, complevel=compression)
        # Write the metadata.
        v.setncatts(var.atts)
        direct_addresses[var.array.tobytes()] = write_data (v, (), var.array)
        continue
      # Hard case: only have the record indices, need to loop over the records.
      # Get the shape of a single record for the variable.
      if hasattr(var, 'chunks'):
        continue  #TODO
      elif hasattr(var,'record_id'):
        record_shape = var.shape[var.record_id.ndim:]
      else:
        continue
      # Use this as the "chunk size" for the netCDF file, to improve I/O
      # performance.
      chunksizes = (1,)*(len(var.axes)-len(record_shape)) + record_shape
      if hasattr(self,'_fill_value') and var.dtype.name.startswith('float32'):
        fill_value = self._fill_value
      else:
        fill_value = None
      # netCDF3 can't handle unsigned ints, so cast to signed.
      dtype = var.dtype
      if dtype.name.startswith('uint') and nc_format.startswith('NETCDF3'):
        warn (_("netCDF3 does not support unsigned ints.  Converting %s to signed int.")%var.name)
        dtype = np.dtype(dtype.name[1:])
      # Tell netCDF4 to use big-endian encoding of data, where applicable.
      # FSTD requires big-endian encoding (and netCDF4 can work with either).
      v = f.createVariable(var.name, datatype=dtype.newbyteorder('>'), endian='big', dimensions=var.dims, zlib=zlib, complevel=compression, chunksizes=chunksizes, fill_value=fill_value)
      # Turn off auto scaling of variables - want to encode the values as-is.
      # 'scale_factor' and 'add_offset' will only be applied when *reading* the
      # the file after it's created.
      v.set_auto_scale(False)
      # Write the metadata.
      v.setncatts(var.atts)
      # Write the data.
      if hasattr(var,'record_id'):
        indices = list(np.ndindex(var.record_id.shape))
        keys = map(int,var.record_id.flatten())
      else:
        indices = list(var.keys())
        keys = list(var.chunks.values())
        record_shape = None  # Reshaping with chunked data not supported.
      for r, ind in zip(keys,indices):
        if r >= 0:
          io.append((r,record_shape,v,ind))

    # Check if no data records exist and no coordinates were converted.
    if len(io) == 0 and len(f.variables) == 0:
      warn(_("No relevant FST records were found."))

    # Now, do the actual transcribing of the data.
    # Read/write the data in the same order of records in the RPN file(s) to
    # improve performance.
    Bar = _ProgressBar if (progress is True and len(io) > 0) else _FakeBar
    bar = Bar(_("Saving netCDF file"), suffix="%(percent)d%% [%(myeta)s]")
    chunk_addresses = {}
    chunk_sizes = {}
    for r,shape,v,ind in bar.iter(sorted(io)):
      try:
        data = self._fstluk(r,dtype=v.dtype.newbyteorder('='))['d'].transpose().reshape(shape)
        addr = write_data (v, ind, data)
        if addr is not None:
          chunk_addresses[r] = addr
        else:
          warn(_("Problem writing compatible record for %s:%s.  Writing separate netCDF / FSTD versions instead.")%(v.name,ind))
      except (IndexError,ValueError):
        warn(_("Internal problem with the script - unable to get data for '%s'")%v.name)
        continue

    # Clean up attributes used for padding.
    for att in list(f.ncattrs()):
      if att.startswith('_pad'):
        f.delncattr(att)

    f.close()

    # Get list of all records that were relevant to the conversion.
    used_rec_ids = set(self._used_rec_ids)
    # Include all metadata records.  Difficult to know for certain which ones
    # were used, since librmn may have read them internally in some routines
    # like horizontal / vertical grid extraction.
    ismeta = self._headers['ismeta']
    used_rec_ids.update(np.where(ismeta==1)[0])

    # Prepare the file for writing FSTD structures.
    with open(filename,'r+b') as f:
      # Disable "aux keys" in FSTD header.
      # Normally, FSTD will expect two 32-bit integers of zeros preceding the
      # actual data (in two unused "aux keys").  If those zeros aren't there,
      # then librmn will abort if you try to read the data.
      # We can't control the bytes just before the data in this scenario,
      # since the file layout is being dictated by netCDF.  As a workaround,
      # we use a non-standard (but still valid) FSTD header that sets the
      # number of "aux keys" to zero, effectively disabling that check.
      f.seek(0x2f,0)
      f.write(b'\0')
      # Pad the end of file to align with 8-byte boundary.
      f.seek(0,2)
      while (f.tell()%8!=0):
        f.write(b'0')
      # Update file size info.
      filesize = f.tell()
      f.seek(0x10,0)
      np.array([filesize//8],'>i4').tofile(f)

    # Find the FST records that weren't directly written to netCDF variables.
    # Examples include coordinates (>>,^^,!!), and mask (typvar=@@) records.
    unwritten_rec_ids = used_rec_ids - set(chunk_addresses.keys())

    # Write anything which isn't already in the file.
    iun = rmn.fstopenall(filename,rmn.FST_RW)
    for rec_id in unwritten_rec_ids:
      d = self._fstluk(rec_id)
      # This data might have already been written as a netCDF coordinate array.
      # E.g., could have got out "lat" coordinate from "^^", but the extraction
      # would have been opaque to us (handled by gdll routine in librmn).
      try:
        chunk_addresses[rec_id] = direct_addresses[d['d'].transpose().tobytes()]
        continue
      except KeyError:
        pass
      # If not, then write it using the original FSTD parameters.
      # This data wasn't used for the netCDF4 interface, so don't have to worry
      # about making it netCDF-compatible.

      #TODO: clean up this hack once public release of Python-RPN is updated.
      from rpnpy.librmn import fstd98
      backup = fstd98.dtype_numpy2fst
      fstd98.dtype_numpy2fst = dtype_numpy2fst
      rmn.fstecr(iun,d)
      fstd98.dtype_numpy2fst = backup
      # Keep track of where the data was written.
      prm = rmn.fstprm(rmn.fstinl(iun)[-1])
      chunk_addresses[rec_id] = (prm['swa']-1)*8+72, (prm['lng']*4)-96, d['datyp'], d['nbits']
    rmn.fstcloseall(iun)

    # By this point, all data is in the file.  This includes netCDF-only data,
    # FSTD-only data, and shared netCDF/FSTD data arrays.
    # The next step is to construct the FSTD record headers to annotate the
    # data so it's accessible from the FSTD interface.

    # First, generate the record header data in the low-level format that will
    # be used in the file.

    rec_ids, addresses = zip(*sorted(chunk_addresses.items()))
    addresses, sizes, datyps, nbits = zip(*addresses)
    # Transform addresses to 64-bit units, rewound for "header", origin at 1.
    addresses = np.asarray(addresses) // 8 - 9 + 1
    # Transform sizes to 64-bit units, using some padding.
    sizes = np.asarray(sizes) // 8 + 12
    datyps = np.asarray(datyps)
    nbits = np.asarray(nbits)
    # Encode the record headers.
    from fstd2nc.extra import structured_array
    headers = structured_array(self._headers)
    headers = headers[list(rec_ids)]
    nrecs = len(headers)
    buf = np.zeros((nrecs,18),'>i4')
    # deleted, select, size
    buf[:,0] = 0x01000000 + sizes
    # address
    buf[:,1] = addresses
    # deet, nbits
    buf[:,2] = (headers['deet']<<8) + np.asarray(nbits,'b').view('B')
    # ni, grtyp
    buf[:,3] = (headers['ni']<<8) + np.asarray(headers['grtyp'],'|S1').view('B')
    # nj, datyp
    buf[:,4] = (headers['nj']<<8) + datyps
    # nk, ubc
    # (nk and ubc no longer guaranteed to be kept in the headers).
    #buf[:,5] = (headers['nk']<<12) + headers['ubc']
    buf[:,5] = (1<<12) + 0
    # npas, pad7
    buf[:,6] = (headers['npas']<<6)
    # ig4, ig2a
    buf[:,7] = (headers['ig4']<<8) + (headers['ig2']>>16)
    # ig1, ig2b
    buf[:,8] = (headers['ig1']<<8) + ((headers['ig2']>>8)%256)
    # ig3, ig2c
    buf[:,9] = (headers['ig3']<<8) + (headers['ig2']%256)
    # etik15, pad1
    etiket = np.asarray(np.array(headers['etiket']).reshape(-1,1).view('B'),'int32')
    etiket -= 32
    buf[:,10] = (etiket[:,0]<<26) + (etiket[:,1]<<20) + (etiket[:,2]<<14) + (etiket[:,3]<<8) + (etiket[:,4]<<2)
    # etik6a, pad2
    buf[:,11] = (etiket[:,5]<<26) + (etiket[:,6]<<20) + (etiket[:,7]<<14) + (etiket[:,8]<<8) + (etiket[:,9]<<2)
    # etikbc, typvar, pad3
    typvar = np.asarray(np.array(headers['typvar']).reshape(-1,1).view('B'),'int32')
    typvar -= 32
    buf[:,12] = (etiket[:,10]<<26) + (etiket[:,11]<<20) + (typvar[:,0]<<14) + (typvar[:,1]<<8)
    # nomvar, pad4
    nomvar = np.asarray(np.array(headers['nomvar']).reshape(-1,1).view('B'),'int32')
    nomvar -= 32
    buf[:,13] = (nomvar[:,0]<<26) + (nomvar[:,1]<<20) + (nomvar[:,2]<<14) + (nomvar[:,3]<<8)
    # ip1, levtyp
    buf[:,14] = (headers['ip1']<<4)
    # ip2, pad5
    buf[:,15] = (headers['ip2']<<4)
    # ip3, pad6
    buf[:,16] = (headers['ip3']<<4)
    # date_stamp
    buf[:,17] = ((headers['datev']//10)<<3) + (headers['datev']%10)

    # Next, write these raw record headers to the file, dividing them into
    # FSTD directory pages.

    with open(filename,'r+b') as f:
      f.seek(0x18,0)
      np.array([nrecs],'>i4').tofile(f)  # number of extensions (?)
      f.seek(0x1c,0)
      np.array([1],'>i4').tofile(f)  # Reset number of pages to one to start
      f.seek(0x24,0)
      np.array(np.max(sizes),'>i4').tofile(f)  # maximum data length
      f.seek(0x34,0)
      np.array([nrecs],'>i4').tofile(f)  # total num records

      for rec0 in range(0,nrecs+1,256):

        nrecs = len(addresses[rec0:rec0+256])

        # Move to page location.
        if rec0 == 0:
          f.seek(0xd0)
        else:
          f.seek(0,2)
          while (f.tell()%8!=0):
            f.write(b'0')
        page = f.tell()//8+1

        # Write page header.
        np.array([2308],'>i4').tofile(f)  # idtyp, header length
        np.array([page],'>i4').tofile(f)  # address of this page
        np.array([0,0],'>i4').tofile(f)  # reserved
        np.array([0,nrecs],'>i4').tofile(f)  # next page, num records in page
        checksum = nrecs ^ 0
        for b in buf[rec0:rec0+256].flatten():
          checksum ^= b
        np.array([checksum,0],'>i4').tofile(f)  # checksum, reserved
        buf[rec0:rec0+256].tofile(f)
        # Pad last page to 256 entries.
        nrecs_written = buf[rec0:rec0+256].shape[0]
        if nrecs_written < 256:
          np.zeros((256-nrecs_written,18),'>i4').tofile(f)

        if rec0 != 0:
          # Update total number of pages.
          f.seek(0x1c,0)
          npages = np.fromfile(f,'>i4',1)[0]
          npages = npages + 1
          f.seek(-4,1)
          np.array([npages],'>i4').tofile(f)

          # Update pointer to last page.
          f.seek(0x20,0)
          prev_page = np.fromfile(f,'>i4',1)[0]
          f.seek(-4,1)
          np.array([page],'>i4').tofile(f)

          # Link this page to the next.
          f.seek((prev_page-1)*8+16,0)
          np.array([page],'>i4').tofile(f)

          # Update checksum.
          f.seek(4,1)
          checksum = np.fromfile(f,'>i4',1)[0]
          checksum ^= page
          f.seek(-4,1)
          np.array([checksum],'>i4').tofile(f)

      # Update file size info (FST header)
      f.seek(0,2)
      filesize = f.tell()
      f.seek(0x10,0)
      np.array([filesize//8],'>i4').tofile(f)
      # Update file size info (netCDF header)
      # (disabled, since it seems to corrupt the netCDF file under some circumstances?)
      # (could be related to deletion of pad attributes)
      #f.seek(0x8028,0)
      #np.array([filesize],'<i4').tofile(f)
