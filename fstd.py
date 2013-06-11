
# Wrapper for FSTD file.  Allows the file to be automatically closed when
# there are no more references to it.
class FSTD_File (object):
  def __init__ (self, filename, mode='r'):
    from pygeode.formats import fstd_core
    if mode == 'r':
      self.iun = fstd_core.open_readonly(filename)
    elif mode == 'w':
      self.iun = fstd_core.open_write(filename)
    else: raise ValueError ("Unknown file mode '%s'"%mode)
    self.mode = mode
  def __del__ (self):
    from pygeode.formats import fstd_core
    fstd_core.close (self.iun)
  def get_records (self):
    from pygeode.formats import fstd_core
    if self.mode != 'r':
      raise IOError ("File not set for read access")
    # Get the record headers
    records = fstd_core.get_record_headers(self.iun)
    # Fill in the data functions.
    # Done here, since it's not easy to implement on the C side.
    map(self.make_data_func, records)
    return records

  def write_records (self, records):
    from pygeode.formats import fstd_core
    if self.mode != 'w':
      raise IOError ("File not set for write access")
    fstd_core.write_records (self.iun, records)

  # Helper function - create a data function for the given handle.
  @staticmethod
  def make_data_func (header):
    handle = header['handle']
    datyp = header['datyp']
    nbits = header['nbits']
    dtypes = {1:'float', 2:'uint', 3:'a', 4:'int', 5:'float', 134:'float', 130:'uint', 132:'int', 133:'float'}
    if datyp not in dtypes:
      raise ValueError ("Can't handle datyp %d"%datyp)
    dtype = dtypes[datyp]
    if dtype in ('float','uint','int'):
      dtype += str(nbits)
    elif dtype == 'a':
      dtype += str(nbits/8)
    ni, nj, nk = header['ni'], header['nj'], header['nk']

    def data_func ():
      from pygeode.formats import fstd_core
      import numpy as np
      out = np.empty([nk,nj,ni], dtype=dtype)
      fstd_core.read_record(handle,out)
      return out

    header['data_func'] = data_func



# Open a file for read access.  Returns a generic 'Dataset' object.
def open (filename):
  f = FSTD_File (filename, mode='r')

  headers = f.get_records()

  # Test the access
  h = headers[headers['nomvar']=='>>  '][0]
  print h
  data = h['data_func']()
  print data

  # Try writing back out
  test = FSTD_File ("writeback.fst", mode='w')
  test.write_records (headers)

