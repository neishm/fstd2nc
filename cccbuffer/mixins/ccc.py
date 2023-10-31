
# Create a Buffer mixin to work with CCC data.
from fstd2nc.mixins import BufferBase
from fstd2nc.stdout import _, info, warn, error
class CCCMixin(BufferBase):
  _format = _("CCC file")
  _format_singular = _("a CCC file")
  _format_plural = _("CCC file(s)")

  _outer_axes = ('time','level')
  _inner_axes = ('ilat','ilg')
  _var_id = ('name','ilat','ilg','khem','superlabel')
  _human_var_id = ('%(name)s','%(ilg)sx%(ilat)s','khem%(khem)s')

  # This method reads header information from a file, and returns it in a raw
  # format close to how it appears in the file.
  # The idea is to get in the file, grab the metadata, and return it as
  # quickly as possible (for the case where many files are being scanned).
  # Decoding will be done in a later step in bulk over all the metadata for
  # all the input files as once.
  # For FSTD files, this speeds things up, but for CCC files (where there is
  # no fast way to scan all the headers) this doesn't really make a difference.
  @staticmethod
  def _raw_headers (filename):
    import numpy as np
    import os
    if not os.path.exists(filename):
      return None
    f = open(filename,'rb')
    # Check for some common markers for CCC files.
    magic = np.fromfile(f,'>i4',19)
    if magic[0] != 64 or magic[-2] != 64:
      f.close()
      return None
    # Read raw headers, and keep track of indices.
    # First record already read from above.
    raw = [magic]
    while len(raw[-1]) > 0:
      f.seek(raw[-1][-1]+4,1)
      raw.append(np.fromfile(f,'>i4',19))
    raw.pop()
    raw = np.array(raw,'>i4')
    # Hijack part of the raw buffer to include the addresses.
    # First, fill with the lengths of each record (plus header)
    raw[0,-2] = 0
    raw[1:,-2] = raw[:-1,-1] + 80
    # Then, addresses are the cumulative sum of these lengths.
    raw[:,-2] = np.cumsum(raw[:,-2])
    return raw.view('B')

  # This method will take the collection of raw headers and decode them
  # into their final fields.  Result is stored internally in the _headers
  # attribute, which is a dictionary of 1D arrays (of size nrecs).
  @staticmethod
  def _decode_headers (raw):
    import numpy as np
    raw = raw.view('>i4')
    headers_int = np.array(raw[:,1:17],'>i4').view('>i8')
    headers_str = headers_int.view('|S8')
    out = {}
    # Standard CCC record paramters
    out['kind'] = np.array(headers_str[:,0])
    out['time'] = np.array(headers_int[:,1],'uint64')
    out['name'] = np.array(headers_str[:,2])
    out['level'] = np.array(headers_int[:,3],'uint64').view('int64')
    out['ilg'] = np.array(headers_int[:,4],'uint64')
    out['ilat'] = np.array(headers_int[:,5],'uint64')
    out['khem'] = np.array(headers_int[:,6],'uint64')
    out['pack'] = np.array(headers_int[:,7],'uint64')
    # Extra info (raw data length and position in file)
    length = raw[:,-1] + 80
    out['length'] = np.array(length,'uint32')
    out['address'] = np.array(raw[:,-2],'uint64')
    out['dtype'] = np.array([np.float32]*len(raw))
    return out

  # This method takes a raw *data* record, and decodes it into the final
  # array of values.  Uses the 'address' and 'length' entries from _headers to
  # determine the location and length of each raw data record.
  # This method is called from the xarray interface, and to_netcdf().
  @staticmethod
  def _decode (data):
    import numpy as np
    header = data[4:68].view('>i8')
    nlon = header[4]
    nlat = header[5]
    pack = np.maximum(header[7],1)
    nbits = 64 // pack
    lng = data[72:76].view('>i4')[0]
    # 64-bit floats?
    if pack == 1:
      return data[76:-4].view('>f8').reshape(nlat,nlon)
    # Packed data?
    xmin = data[76:84].view('>f8')[0]
    xmax = data[84:92].view('>f8')[0]
    dtype = {1:'>i8', 2:'>i4', 4:'>H', 8:'B'}[pack]
    data = data[92:-4].view(dtype).reshape(nlat,nlon)
    data = data * (xmax-xmin) / (2**(nbits-1)) + xmin
    return np.array(data,'float32')

