# Inferface for reading/writing FSTD files within PyGeode.

from pygeode.axis import Axis

# Low-level axes that are directly mapped to FSTD files.
class DateoAxis(Axis): name = 'dateo'
class IAxis(Axis): name = 'i'
class JAxis(Axis): name = 'j'
class KAxis(Axis): name = 'k'
class IP1Axis(Axis): name = 'ip1'
class IP2Axis(Axis): name = 'ip2'
class IP3Axis(Axis): name = 'ip3'
class NPASAxis(Axis): name = 'npas'

# High-level axes that are determined from the specific context of the data.
# Only define things not already available in the standard PyGeode module.

# Forecast axis
class Forecast(Axis):
  name = 'forecast'
  atts = dict(units='hours')

# Vertical coordinates
from pygeode.axis import ZAxis
class Height_wrt_SeaLevel(ZAxis):
  name = 'height'
  atts = dict(ZAxis.atts, units='m', standard_name='height_above_sea_level')

class Height_wrt_Ground(ZAxis):
  name = 'height'
  atts = dict(ZAxis.atts, units='m', standard_name='height')

class Sigma(ZAxis):
  name = 'sigma'
  atts = dict(ZAxis.atts, standard_name='atmosphere_sigma_coordinate')
  plotatts = dict(ZAxis.plotatts, plotorder=-1, plotscale='log')

class LogHybrid(ZAxis):
  name = 'zeta'
  atts = dict(ZAxis.atts, standard_name='atmosphere_hybrid_sigma_log_pressure_coordinate')  # Not really a standard
  plotatts = dict(ZAxis.plotatts, plotorder=-1, plotscale='log')

class Theta(ZAxis):
  name = 'theta'
  atts = dict(ZAxis.atts, units='K', standard_name='air_potential_temperature')

del Axis, ZAxis


# Helper class - collects data references into a numpy-like object.
# Makes it more convenient for accessing the underlying data.
class _Array (object):
  def __init__ (self, shape, dtype, inner_dims, data_funcs):
    # Expected shape and type of the array.
    self.shape = tuple(map(int,shape))
    self.ndim = len(self.shape)
    self.size = reduce(int.__mul__, self.shape, 1)
    self.dtype = dtype
    # Boolean list, indicating which dimensions are part of the inner data_funcs.
    self._inner_dims = inner_dims
    # Array of references to the data from the FSTD records.
    self._data_funcs = data_funcs
  def __getitem__ (self, key):
    from itertools import product
    # Coerce key into a tuple of slice objects.
    if not isinstance(key,tuple):
      key = (key,)
    if Ellipsis in key:
      i = key.index(Ellipsis)
      key = key[:i] + (slice(None),)*(self.ndim-len(key)+1) + key[i+1:]
    key = key + (slice(None),)*(self.ndim-len(key))
    # Collect keys into inner / outer parts
    outer_keys = []
    inner_keys = []
    for k,is_inner in zip(key,self._inner_dims):
      if is_inner: inner_keys.append(k)
      else: outer_keys.append(k)
    outer_keys = tuple(outer_keys)
    inner_keys = tuple(inner_keys)
    data_funcs = self._data_funcs
    # Apply the outer keys.
    if len(outer_keys) > 0:
      data_funcs = data_funcs[outer_keys]
    # Apply the inner keys.
    if inner_keys != (slice(None),)*len(inner_keys):
      #print "Inner slicing triggered."
      if hasattr(data_funcs,'shape'):
        for ind in product(*map(range,data_funcs.shape)):
          f = data_funcs[ind]
          data_funcs[ind] = lambda old_f=f: old_f()[inner_keys]
      # data_funcs may be completely collapsed, in which case there is no
      # 'shape' attribute (no longer an ndarray).
      else:
        f = data_funcs
        data_funcs = lambda old_f=f: old_f()[inner_keys]

    new_shape = [range(s)[k] for s,k in zip(self.shape,key)]
    # Check for dimensions that have been reduced out.
    inner_dims = [i for i,s in zip(self._inner_dims,new_shape) if not isinstance(s,int)]
    new_shape = [len(s) for s in new_shape if not isinstance(s,int)]
    return _Array (new_shape, self.dtype, inner_dims, data_funcs)
  def __array__ (self):
    import numpy as np
    from itertools import product
    data = np.empty(self.shape, dtype=object)
    # Get indices of all dimensions, in preparation for iterating.
    indices = [[slice(None)] if is_inner else range(s) for is_inner,s in zip(self._inner_dims, self.shape)]
    if hasattr(self._data_funcs,'flatten'):
      for i,ind in enumerate(product(*indices)):
        data[ind] = self._data_funcs.flatten()[i]()
    # data_funcs is degenerate (single element)?
    else:
      data[:] = self._data_funcs()

    return data


# Define a class for encoding / decoding FSTD data.
# Each step is placed in its own method, to make it easier to patch in new
# behaviour if more exotic FSTD files are encountered in the future.
class _Buffer_Base (object):
  """
  Interface for reading / writing FSTD files.
  """

  # Names of records that should be kept separate (never grouped into
  # multidimensional arrays).
  _meta_records = ('>>', '^^', 'HY', '!!', 'HH', 'STNS', 'SH')

  # Attributes which could potentially be used as axes.
  _outer_axes = ()

  # Attributes which uniquely identify a variable.
  _var_id = ('nomvar','ni','nj','nk')

  # Attributes which should be completely ignored when decoding.
  # They're either not implemented, or are internal info for the file.
  _ignore_atts = ('pad','swa','lng','dltf','ubc','extra1','extra2','extra3','data_func')

  def __init__ (self):
    import numpy as np
    dtype = [("dateo", "i4"), ("deet", "i4"), ("npas", "i4"), ("ni", "i4"), ("nj", "i4"), ("nk", "i4"), ("nbits", "i4"), ("datyp", "i4"), ("ip1", "i4"), ("ip2", "i4"), ("ip3", "i4"), ("typvar", "a2"), ("nomvar", "a4"), ("etiket", "a12"), ("grtyp", "a2"), ("ig1", "i4"), ("ig2", "i4"), ("ig3", "i4"), ("ig4", "i4"), ("swa", "i4"), ("lng", "i4"), ("dltf", "i4"), ("ubc", "i4"), ("extra1", "i4"), ("extra2", "i4"), ("extra3", "i4")]
    self.headers = np.array([],dtype=dtype)
    self._data_funcs = []

  def clear (self):
    """
    Removes all existing data from the buffer.
    """
    import numpy as np
    self.headers = np.array([],dtype=self.headers.dtype)
    self._data_funcs = []


  ###############################################
  # Basic flow for reading data

  def read_file (self, filename):
    """
    Read raw records from an FSTD file, into the buffer.
    Multiple files can be read sequentially.
    """
    from pygeode.formats import fstd_core
    import numpy as np
    # Read the data
    records = fstd_core.read_records(filename)
    # Store the headers and data interfaces.
    headers = np.zeros(len(records),dtype=self.headers.dtype)
    for n in headers.dtype.names:
      headers[n] = records[n]
    self.headers = np.concatenate((self.headers,headers))
    self._data_funcs.extend(records['data_func'])


  # Decode the record headers into a dictionary, and
  # add extra (derived) attributes that are useful for
  # processing the data into multidimensional arrays.
  def _get_fields (self):
    from collections import OrderedDict
    import numpy as np
    fields = OrderedDict()
    for n in self.headers.dtype.names:
      fields[n] = self.headers[n]
    fields['data_func'] = np.asarray(self._data_funcs)
    return fields


  def __iter__ (self):
    """
    Processes the records into multidimensional variables.
    Iterates over (name, atts, axes, array) tuples.
    Note that array may not be a true numpy array (values are not yet loaded
    in memory).  To load the array, pass it to numpy.array().
    """
    from pygeode.formats import fstd_core
    from collections import OrderedDict, namedtuple
    import numpy as np

    records = self._get_fields()

    # Group the records by variable.
    var_records = OrderedDict()
    var_id_type = namedtuple('var_id', self._var_id)
    for i in range(len(records['nomvar'])):
      # Get the unique variable identifiers.
      var_id = var_id_type(*[records[n][i] for n in self._var_id])
      # Ignore coordinate records.
      if var_id.nomvar.strip() in self._meta_records: continue
      var_records.setdefault(var_id,[]).append(i)

    var_type = namedtuple('var', ('name','atts','axes','array'))

    # Loop over each variable and construct the data & metadata.
    for var_id, rec_ids in var_records.iteritems():
      # Get the metadata for each record.
      atts = OrderedDict()
      for n in records.iterkeys():
        if n in self._outer_axes or n in self._ignore_atts: continue
        v = records[n][rec_ids]
        # Only use attributes that are consistent across all variable records.
        if len(set(v)) > 1: continue
        v = v[0]
        # Trim string attributes (remove whitespace padding).
        if isinstance(v,str): v = v.strip()
        atts[n] = v

      # Get the axis coordinates.
      axes = OrderedDict((n,tuple(sorted(set(records[n][rec_ids])))) for n in self._outer_axes)
      axes['k'] = range(var_id.nk)
      axes['j'] = range(var_id.nj)
      axes['i'] = range(var_id.ni)

      # Construct a multidimensional array to hold the data functions.
      data_funcs = np.empty(map(len,axes.values()[:-3]), dtype='O')

      # Assume missing data (nan) unless filled in later.
      def missing_data(ni=var_id.ni, nj=var_id.nj, nk=var_id.nk):
        import numpy as np
        data = np.empty((ni,nj,nk), dtype=float)
        data[()] = np.float('nan')
        return data
      data_funcs[()] = missing_data
      
      # Arrange the data funcs in the appropriate locations.
      for rec_id in rec_ids:
        index = tuple(axes[n].index(records[n][rec_id]) for n in self._outer_axes)
        data_funcs[index] = records['data_func'][rec_id]

      # Check if we have full coverage along all axes.
      have_data = np.not_equal(data_funcs,None)
      if not np.all(have_data):
        from warnings import warn
        warn ("Missing some records for %s.")

      data = _Array (shape = data_funcs.shape+(var_id.nk,var_id.nj,var_id.ni),
                     dtype = 'float32',
                     inner_dims = [False]*data_funcs.ndim+[True,True,True],
                     data_funcs = data_funcs
             )

      yield var_type(var_id.nomvar.strip(), atts, axes, data)

      #TODO: Find a minimum set of partial coverages for the data.
      # (e.g., if we have surface-level output for some times, and 3D output
      # for other times).

  #
  ###############################################



  ###############################################
  # Basic flow for writing data

  # Add some derived data to this object.  The data must follow the same format
  # as the output of __iter__().
  def _add_data (self, data):
    import numpy as np
    from itertools import product
    from collections import OrderedDict

    # Get the type for each attribute (so we can use appropriate fill values
    # when certain attributes are missing.
    att_type = OrderedDict()
    for var_id, atts, axes, data_funcs in data:
      for n,v in atts.iteritems():
        att_type[n] = type(v)
      for n,v in axes.iteritems():
        if n in ('k','j','i'): continue
        att_type[n] = type(v)
    att_type['data_func'] = object

    # Construct the records.
    # The values will be stored in chunks, to be concatenated at the end.
    records = OrderedDict((n,[]) for n in att_type.iterkeys())

    # Loop over each variable.
    for varname, atts, axes, array in data:

      # Make sure we have nk,nj,ni dimensions, and in the right order.
      for n in 'k','j','i':
        if n not in axes:
          raise KeyError("'%s' axis not found in the data.")
      if axes.keys()[-3:] != ['k','j','i']:
        raise ValueError("Wrong dimension order - expected (nk,nj,ni) dimensions at the end.")

      # Make sure we have a nomvar.
      if 'nomvar' not in atts:
        atts['nomvar'] = varname[:4]

      # Wrap the data array into callable functions.
      data_funcs = []
      for ind in product(*map(range,array.shape[:-3])):
        data_funcs.append(lambda a=array[ind]: np.asarray(a))

      current_records = OrderedDict()

      # Add coordinate info.
      for coords in product(*axes.values()[:-3]):
        for n,c in zip(axes.keys(),coords):
          current_records.setdefault(n,[]).append(c)

      # Add metadata.
      nrecs = len(data_funcs)
      for n,v in atts.iteritems():
        current_records[n] = [v]*nrecs

      # Add in the data references.
      current_records['data_func'] = data_funcs

      # Check for anything that's not defined (apply a mask).
      for n,t in att_type.iteritems():
        if n not in current_records:
          bad = np.zeros(nrecs, dtype=t)
          current_records[n] = np.ma(bad,mask=True)

      # Add these records to the master list.
      for n,v in records.iteritems():
        v.append(current_records[n])

    # Convert to numpy arrays.
    records = OrderedDict((n,np.ma.concatenate(v)) for n,v in records.iteritems())

    # Process this into headers.
    self._put_fields(nrecs, records)

  # Final preparation of records, into proper headers.
  def _put_fields (self, nrecs, fields):
    import numpy as np

    # Pad out the string records with spaces
    fields['nomvar'] = [s.upper().ljust(4) for s in fields['nomvar']]
    fields['etiket'] = [s.upper().ljust(12) for s in fields['etiket']]
    fields['typvar'] = [s.upper().ljust(2) for s in fields['typvar']]
    fields['grtyp'] = map(str.upper, fields['grtyp'])

    headers = np.zeros(nrecs,dtype=self.headers.dtype)
    for n in headers.dtype.names:
      if n in fields:
        headers[n] = fields[n]
    self.headers = np.concatenate((self.headers,headers))
    self._data_funcs.extend(fields['data_func'])



  def write_file (self, filename):
    """
    Write the records from this object into the specified FSTD file.
    """
    from pygeode.formats.fstd_core import record_descr
    from pygeode.formats import fstd_core
    import numpy as np

    # Create a numpy structured array to hold the data.
    records = np.zeros(len(self.headers),dtype=record_descr)

    # Fill in what we have from our existing records.
    for n in record_descr.names:
      if n in self.headers.dtype.names:
        records[n] = self.headers[n]
    records['data_func'] = self._data_funcs

    # Write out the data.
    fstd_core.write_records (filename, records)


  #
  ###############################################




#####################################################################
# Mixins for different behaviour / assumptions about FSTD data.


#################################################
# Logic for handling date field.

class _Dates (_Buffer_Base):

  def __init__ (self, squash_forecasts=False, *args, **kwargs):
    super(_Dates,self).__init__(*args,**kwargs)
    self.squash_forecasts = squash_forecasts
    if squash_forecasts:
      self._outer_axes = ('time',) + self._outer_axes
    else:
      self._outer_axes = ('time','forecast') + self._outer_axes

  # Get any extra (derived) fields needed for doing the decoding.
  def _get_fields (self):
    from pygeode.formats import fstd_core
    from collections import OrderedDict
    from datetime import datetime, timedelta
    import numpy as np
    fields = super(_Dates,self)._get_fields()
    # Calculate the forecast (in hours) and date of validity.
    fields['forecast']=fields['deet']*fields['npas']/3600.
    dateo = fstd_core.stamp2date(fields['dateo'])
    datev = dateo + fields['deet']*fields['npas']
    # This isn't really needed by the decoder, but it provided for
    # convenience to the user.
    fields['datev'] = np.array(fstd_core.date2stamp(datev))
    # Time axis
    if self.squash_forecasts:
      dates = datev
    else:
      dates = dateo
    date0 = datetime(year=1980,month=1,day=1)
    dt = np.array([timedelta(seconds=int(s)) for s in dates])
    fields['time'] = date0 + dt
    return fields

  # Prepare date info for output.
  def _put_fields (self, nrecs, fields):
    from pygeode.formats import fstd_core
    import numpy as np
    from datetime import datetime

    # Check for missing record attributes
    for att in 'deet','npas','dateo','forecast':
      if att not in fields:
        fields[att] = np.ma.masked_all(nrecs)
    if 'time' not in fields:
      fields['time'] = np.ma.masked_all(nrecs,dtype='O')

    # Mask out invalid deet values.
    fields['deet'] = np.ma.masked_equal(fields['deet'],0)

    # Set default values.
    date0 = datetime(year=1980,month=1,day=1)
    fields['forecast'] = np.ma.filled(fields['forecast'],0)
    fields['time'] = np.ma.filled(fields['time'],date0)

    # Get dateo from 'time' field, wherever it's missing.
    data = fstd_core.date2stamp([t.seconds for t in (fields['time']-date0)])
    data = np.array(data)
    bad_dateo = np.ma.getmaskarray(fields['dateo'])
    fields['dateo'][bad_dateo] = data[bad_dateo]

    # Calculate npas from forecast, wherever npas is missing.
    bad_npas = np.ma.getmaskarray(fields['npas'])
    data = fields['forecast']*3600./fields['deet']
    fields['npas'][bad_npas] = np.asarray(data[bad_npas],dtype='int32')

    # Fill in missing values
    for att in 'deet','npas','dateo','forecast':
      fields[att] = np.ma.filled(fields[att],0)

    super(_Dates,self)._put_fields(nrecs, fields)


# Logic for handling vertical coordinates.
class _VCoords (_Buffer_Base):
  # Don't group records across different level 'kind'.
  # (otherwise can't create a coherent vertical axis).
  def __init__ (self, *args, **kwargs):
    super(_VCoords,self).__init__(*args,**kwargs)
    # Split data by level kind.
    self._var_id = self._var_id + ('kind',)
    # Use decoded IP1 values as the vertical axis.
    self._outer_axes += ('level',)
  def _get_fields (self):
    from pygeode.formats import fstd_core
    fields = super(_VCoords,self)._get_fields()
    # Provide 'level' and 'kind' information to the decoder.
    levels, kind = fstd_core.decode_levels(fields['ip1'])
    fields['level'] = levels
    fields['kind'] = kind
    return fields


#################################################
# Logic for reading/writing FSTD data from/to netCDF files.
# Note: this is not strictly an FSTD thing, but it's
# provided here for convenience.
class _netCDF_IO (_Buffer_Base):
  # Add a 'time' variable to the processed data.
  # NOTE: this should ideally go in the '_Dates' mixin, but it's put here
  # due to a dependence on the netCDF4 module.
  def __iter__ (self):
    from datetime import datetime
    from collections import OrderedDict
    import numpy as np
    from netCDF4 import date2num
    # Keep track of all time axes found in the data.
    time_axes = OrderedDict()
    for var in super(_netCDF_IO,self).__iter__():
      yield var
      if 'time' not in var.axes: continue
      taxis = var.axes['time']
      if not isinstance(taxis[0],datetime): continue
      if taxis in time_axes: continue
      units = 'hours since %s'%(taxis[0])
      atts = OrderedDict(units=units)
      axes = OrderedDict(time=taxis)
      array = np.asarray(date2num(taxis,units=units))
      time_axes[taxis] = type(var)('time',atts,axes,array)
    # Append the time axes to the data.
    for times in time_axes.values():
      yield times

  def write_nc_file (self, filename):
    """
    Write the records to a netCDF file.
    Requires the netCDF4 package.
    """
    from netCDF4 import Dataset
    import numpy as np
    f = Dataset(filename, "w", format="NETCDF4")

    # Only write one copy of each unique dimension.
    dims = dict()

    for varname, atts, axes, array in iter(self):
      for axisname, axisvalues in axes.items():
        # Only need to create each dimension once (even if it's in multiple
        # variables).
        if axisname not in f.dimensions:
          f.createDimension(axisname, len(axisvalues))
      # Write the variable.
      v = f.createVariable(varname, datatype=array.dtype, dimensions=axes.keys())
      v.setncatts(atts)
      v[:] = np.asarray(array[:])
    f.close()

# Default interface for I/O.
class Buffer (_Dates,_VCoords,_netCDF_IO):
  """
  High-level interface for FSTD data, to treat it as multi-dimensional arrays.
  Contains logic for dealing with most of the common FSTD file conventions.
  """




#####################################################################
# Everything below this point is old code, which needs to be
# reorganized or removed.

def old():
    from pygeode.var import Var
    from pygeode.formats import fstd_core
    import numpy as np

    name = str(records[0]['nomvar']).rstrip()

    # Get the dates, forecast hours, and levels.
    dates = fstd_core.stamp2date(records['dateo'])
    forecasts = records['deet']*records['npas']
    if squash_forecasts:
      dates += forecasts
      forecasts[:] = 0

    levels, kind = fstd_core.decode_levels(records['ip1'])

    # Get unique values of these arrays
    dates, idate = np.unique(dates, return_inverse=True)
    forecasts, iforecast = np.unique(forecasts, return_inverse=True)
    # Decode the numeric values of the levels, just to get to proper order
    levels, indlevel, ilevel = np.unique(levels, return_index=True, return_inverse=True)
    ip1 = records['ip1'][indlevel]

    # Construct a multidimensional array of data functions.
    # One function per date,forecast,level.
    data_funcs = np.empty([len(dates),len(forecasts),len(levels)], dtype='O')
    data_funcs[:] = None
    data_funcs[idate,iforecast,ilevel] = records['data_func']

    if any(f is None for f in data_funcs.flatten()):
      print "Unable to construct a full field for %s - missing some expected records:"%name
      from datetime import datetime, timedelta
      for i,date in enumerate(dates):
        # Convert to printable date
        date = int(date)
        date = datetime(year=1980,month=1,day=1) + timedelta(seconds=date)
        for j,forecast in enumerate(forecasts):
          forecast = int(forecast)/3600.
          for k,level in enumerate(levels):
            header = "%s - %3gh - level = %8g"%(date,forecast,level)
            if data_funcs[i,j,k] is None:
              print "%s: MISSING"%header
            else:
              print "%s: found"%header
      raise ValueError ("Missing some records that are needed to fill out the (time,forecast,level) dimensions.")

    self.data_funcs = data_funcs

    # Construct the time/forecast axes
    taxis = StandardTime (startdate={'year':1980,'month':1}, values=dates, units='seconds')
    faxis = Forecast (values=forecasts/3600.)
    faxis.atts['deet'] = int(records[0]['deet'])

    nomvar = str(records[0]['nomvar']).rstrip()
    typvar = str(records[0]['typvar']).rstrip()
    etiket = str(records[0]['etiket']).rstrip()
    grtyp = str(records[0]['grtyp']).rstrip()
    ni = int(records[0]['ni'])
    nj = int(records[0]['nj'])
    nk = int(records[0]['nk'])
    ig1 = int(records[0]['ig1'])
    ig2 = int(records[0]['ig2'])
    ig3 = int(records[0]['ig3'])
    ig4 = int(records[0]['ig4'])

    # Construct the i,j,k,z axes

    iaxis = IAxis(ni)
    jaxis = JAxis(nj)
    kaxis = KAxis(nk)

    # Vertical axis
    # (to be decoded later)
    zaxis = IP1Axis(values=ip1)

    # Determine the dtype to use
    # Use the first dtype found.
    datyp = int(records[0]['datyp'])
    nbits = int(records[0]['nbits'])
    dtype = {0:'float', 1:'float', 2:'uint', 3:'a', 4:'int', 5:'float', 134:'float', 130:'uint', 132:'int', 133:'float'}[datyp]
    if dtype == 'a':
      dtype += str(nbits/8)
    else:
      dtype += '64' if nbits > 32 else '32'

    atts = dict(nomvar=nomvar, typvar=typvar, etiket=etiket, grtyp=grtyp, ig1=ig1, ig2=ig2, ig3=ig3, ig4=ig4, datyp=datyp, nbits=nbits)

    # Finish initializing
    Var.__init__(self, [taxis,faxis,zaxis,kaxis,jaxis,iaxis], dtype=dtype, name=name, atts=atts)


# Helper function - preload a record.
def preload(record):
  data = record['data_func']()
  def data_func(): return data
  record['data_func'] = data_func

# Locate and apply mask (0/1) records to corresponding data records.
def apply_masks (records, fill_value):
  import numpy as np
  from warnings import warn
  # Exit early if no masks found.
  is_mask = (records['datyp']==2) & (records['nbits']==1)
  if not np.any(is_mask): return records
  # Group the data and masks together.
  matchers = ('dateo','ni','nj','nk','ip1','ip2','ip3','nomvar','etiket','grtyp','ig1','ig2','ig3','ig4')
  groups = {}
  for i in range(len(records)):
    key = tuple(records[m][i] for m in matchers)
    groups.setdefault(key,[]).append(i)
  masked = {}
  for key, indices in groups.iteritems():
    the_mask = [i for i in indices if records['datyp'][i] == 2 and records['nbits'][i] == 1]
    non_mask = sorted(set(indices)-set(the_mask))
    # If there's no mask for this record, then nothing to do
    if len(the_mask) == 0: continue
    # If there's no non-mask data that fits, then there's nothing to apply the
    # mask to (so, leave it alone).
    if len(non_mask) == 0: continue
    # If there's multiple masks, then don't know which one to use
    if len(the_mask) > 1:
      warn("Multiple masks found for '%s' - don't know what to do."%records['nomvar'][the_mask[0]])
      continue
    # Okay, now we should have exactly one mask to apply.
    def do_mask(x,y):
      x[y==0] = fill_value
      return x
    the_mask = the_mask[0]
    for i in non_mask:
      records['data_func'][i] = lambda f=records['data_func'][i], m=records['data_func'][the_mask]: do_mask(f(),m())
    # Applied the mask, so we don't need the mask record anymore.
    indices.remove(the_mask)
  # Collect all the indices that we will keep.
  return records[sorted(i for indices in groups.itervalues() for i in indices)]


# Attach lat/lon arrays to a list of FSTD variables
# For 1D lat/lon, use these arrays as axes.
# For 2D lat/lon, create 2D coordinate variables.
def attach_latlon (varlist, latlon_arrays):
  from pygeode.var import Var

  ydim, xdim = 4, 5

  handled_latlon_vars = {}
  extra_coord_vars = []

  for var in varlist:
    # Generate key into lat/lon lookup table
    key = tuple(var.atts[a] for a in ['grtyp','ig1','ig2','ig3','ig4'])
    key += (var.shape[xdim], var.shape[ydim])   # ni, nj
    if key not in latlon_arrays: continue

    # Attach 1D lat/lon arrays as axes
    axes = list(var.axes)
    ax, ay, lat, lon = latlon_arrays[key]
    if lat.ndim == 1: axes[ydim] = Lat(lat)
    elif ay is not None: axes[ydim] = YAxis(ay,name='y')
    else: axes[ydim] = YAxis(lat.shape[0],name='y')
    if lon.ndim == 1: axes[xdim] = Lon(lon)
    elif ax is not None: axes[xdim] = XAxis(ax,name='x')
    else: axes[xdim] = XAxis(lat.shape[1],name='x')
    var.axes = tuple(axes)

    # Convert 2D lat/lon arrays to variables
    if key in handled_latlon_vars: continue
    if lat.ndim == 2:
      lat_var = Var([axes[ydim],axes[xdim]], values=lat, name="lat")
      extra_coord_vars.append(lat_var)
    if lon.ndim == 2:
      lon_var = Var([axes[ydim],axes[xdim]], values=lon, name="lon")
      extra_coord_vars.append(lon_var)
    handled_latlon_vars[key] = True

  varlist.extend (extra_coord_vars)

# Attach vertical axes to the variables
def attach_vertical_axes (varlist, vertical_records):
  from pygeode.formats import fstd_core
  import numpy as np

  zdim = 2

  bangbang_cache = {}

  for var in varlist:
    # Skip derived fields
    if not isinstance(var,FSTD_Var): continue

    axes = list(var.axes)

    # Decode the values
    ip1 = var.axes[zdim].values
    levels, kind = fstd_core.decode_levels(ip1)

    if kind == 0:
      axes[zdim] = Height_wrt_SeaLevel(levels)
    elif kind == 1:
      axes[zdim] = Sigma(levels)
    elif kind == 2:
      axes[zdim] = Pres(levels)
    elif kind == 3:
      axes[zdim] = ZAxis(levels)
    elif kind == 4:
      axes[zdim] = Height_wrt_Ground(levels)
    elif kind == 5:
      # Find a vertical record that matches
      # First, look for a !! record
      match = (vertical_records['ip1'] == var.atts['ig1']) & (vertical_records['ip2'] == var.atts['ig2']) & (vertical_records['ip3'] == var.atts['ig3'])
      if any(match):
        bangbang_record = vertical_records[match]
        key = int(bangbang_record['ip1'][0]), int(bangbang_record['ip2'][0]), int(bangbang_record['ip3'][0])
        if key not in bangbang_cache:
          data = bangbang_record[0]['data_func']()
          bangbang_cache[key] = fstd_core.decode_loghybrid_table(data)
        table = bangbang_cache[key]
        # Determine the specific A & B for this axis
        A, B = fstd_core.get_loghybrid_a_b(ip1, table)
        axes[zdim] = LogHybrid(values=levels, A=A, B=B)
        # Store some auxiliary vertical information
        # (needed for going back to FSTD format)
        axes[zdim].atts.update(table)

      else:
        # Otherwise, look for a HY record (looser search criteria)
        match = (vertical_records['nomvar'] == 'HY  ')
        if any(match):
          hy_record = vertical_records[match]
          ptop, rcoef, pref, A, B = fstd_core.get_hybrid_a_b(hy_record, levels)
          axes[zdim] = Hybrid(values=levels, A=A, B=B)
          axes[zdim].atts = dict(ptop=ptop, rcoef=rcoef, pref=pref)
    elif kind == 6:
      axes[zdim] = Theta(levels)

    # Save the new axis information back into the variable
    var.axes = tuple(axes)

# Reduce the dimensionality of the given FSTD variable
def reduce_dimensionality (var, squash_forecasts=False):
  # Skip derived fields
  if not isinstance(var, FSTD_Var): return var

  remove_axes = []
  # Forecast (axis 1)
  if var.shape[1] == 1:
    if squash_forecasts:
      remove_axes += [1]
  # Vertical (axis 2)
  # Surface fields have a 'pressure' coordinate with a value of 0hPa
  if var.shape[2] == 1:
    if isinstance(var.axes[2], Pres) and var.axes[2].values == [0.]:
      remove_axes += [2]
  # K axis (axis 3)
  if var.shape[3] == 1:
    remove_axes += [3]

  return var.squeeze(*remove_axes)

#####################################################################
#
# Open a file for read access.  Returns a generic 'Dataset' object.
#
#####################################################################
def open (filename, squash_forecasts=False, print_warnings=True, raw_list=False, fill_value=1e30):

  from pygeode.formats import fstd_core
  import numpy as np

  # What header attributes define a unique variable
  unique_var_atts = ['nomvar', 'typvar', 'etiket', 'ni', 'nj', 'nk', 'grtyp', 'ip3', 'ig1', 'ig2', 'ig3', 'ig4']

  # Read the records
  records = fstd_core.read_records(filename)
  # Filter out any missing records
  records = records[records['nomvar']!='']

  # Warn about any raw binary records.
  raw_binary_records = records[records['datyp'] == 0]
  if len(raw_binary_records) > 0:
    from warnings import warn
    raw_nomvars = list(set(raw_binary_records['nomvar']))
    warn ("Raw binary records detected for %s.  The values may not be properly decoded if you're opening on a different platform."%raw_nomvars)

  # Locate and apply mask (0/1) records to corresponding data records.
  # (e.g., for RIOPS data).
  records = apply_masks (records, fill_value)

  # Construct all possible lat/lon arrays from info in the records
  latlon_arrays = fstd_core.get_latlon(records)

  # Remove the coordinate records (not needed anymore).
  # Keep vertical descriptors in a separate array, though, because we'll need
  # to evaluate them once we know the vertical dimension of the variables.
  nomvar = records['nomvar']
  is_coord = (nomvar == '>>  ') | (nomvar == '^^  ') | (nomvar == 'HY  ') | (nomvar == '!!  ')
  vertical_records = records[ (nomvar == 'HY  ') | (nomvar == '!!  ') ]
  records = records[-is_coord]
  del nomvar, is_coord

  # Group the records together
#  all_keys = records[unique_var_atts]
  # Above doesn't work in numpy >= 1.8, due to a bug when using record arrays
  # that contain object fields (http://github.com/numpy/numpy/issues/3256)
  all_keys = np.array(zip(*[records[a] for a in unique_var_atts]),dtype=[(a,records.dtype[a]) for a in unique_var_atts])
  unique_keys, var_indices = np.unique(all_keys, return_inverse=True)

  var_bins = [ records[var_indices==i] for i in range(len(unique_keys)) ]

  del all_keys, unique_keys, var_indices

  # Create the variables
  varlist = []
  for var_records in var_bins:
    var = FSTD_Var (var_records, squash_forecasts)
    # Add fill value, in case we applied any masks.
    var.atts['_FillValue'] = fill_value
    varlist.append(var)

  # Attach any lat/lon coordinates
  attach_latlon (varlist, latlon_arrays)

  # Attach vertical coordinates
  attach_vertical_axes (varlist, vertical_records)

  # Dimensionality reduction
  varlist = [reduce_dimensionality(var,squash_forecasts) for var in varlist]

  if raw_list is True:
    return varlist

  # Return the variables as a dataset
  from pygeode.dataset import Dataset
  return Dataset(varlist, print_warnings=print_warnings)

#####################################################################


# Encode the time axis
def encode_time_axis (varlist):
  from pygeode.timeaxis import StandardTime
  from pygeode.timeutils import reltime
  from pygeode.formats import fstd_core
  import numpy as np
  for i,var in enumerate(varlist):
    if not var.hasaxis(StandardTime): continue
    time = var.getaxis(StandardTime)
    seconds = reltime (time, startdate=dict(year=1980,month=1,day=1), units='seconds')
    seconds = np.asarray(seconds,dtype=int)
    values = fstd_core.date2stamp(seconds)
    taxis = Dateo(values=values)
    varlist[i] = var.replace_axes(time=taxis)

# Convert to FSTD-compatible axes
# (e.g., detect hybrid / log-hybrid axes, forecast axis)
def detect_fstd_axes (varlist):
  for varnum,var in enumerate(varlist):
    axes = list(var.axes)
    # Look for serialized FSTD axes (represented by metadata)
    for i,axis in enumerate(axes):
      if isinstance(axis,(Height_wrt_SeaLevel,Height_wrt_Ground,Sigma,Hybrid,LogHybrid,Theta,Forecast,XAxis,YAxis,KAxis)): continue
      standard_name = axis.atts.get('standard_name','')
      if standard_name == 'height_above_sea_level':
        axes[i] = Height_wrt_SeaLevel(values=axis.values)
      elif standard_name == 'height':
        axes[i] = Height_wrt_Ground(values=axis.values)
      elif standard_name == 'atmosphere_sigma_coordinate':
        axes[i] = Sigma(values=axis.values)
      elif standard_name == 'atmosphere_hybrid_sigma_pressure_coordinate':
        axes[i] = Hybrid(values=axis.values, A=axis.auxarrays['A'], B=axis.auxarrays['B'], atts=axis.atts)
        axes[i].atts = dict(**axis.atts)
      elif standard_name == 'atmosphere_hybrid_sigma_log_pressure_coordinate':
        axes[i] = LogHybrid(values=axis.values, A=axis.auxarrays['A'], B=axis.auxarrays['B'], atts=axis.atts)
      elif standard_name == 'air_potential_temperature':
        axes[i] = Theta(values=axis.values)

      elif axis.name == 'forecast':
        axes[i] = Forecast(values=axis.values, atts=axis.atts)

      replacements = {}
      for original_axis, new_axis in zip(var.axes,axes):
        if new_axis is not original_axis:
          replacements[original_axis.name] = new_axis
      if len(replacements) > 0:
        varlist[varnum] = var.replace_axes(**replacements)

# Encode the forecast axis
def encode_forecast_axis (varlist):
  for i,var in enumerate(varlist):
    if not var.hasaxis(Forecast): continue
    forecast = var.getaxis(Forecast)
    if 'deet' not in forecast.atts: continue
    deet = forecast.atts['deet']
    # Special case: analysis (deet could be zero)
    if len(forecast) == 1 and forecast.values[0] == 0:
      npas = [0]
    else:
      npas = forecast.values * 3600 / deet
    npas_axis = NPASAxis(values=npas)
    varlist[i] = var.replace_axes(forecast=npas_axis)
    varlist[i].atts['deet'] = deet

# Encode vertical information into FSTD records
def encode_vertical (varlist):
  vertical_records = []
  from pygeode.formats import fstd_core
  import numpy as np
  from warnings import warn

  # First, check for any hybrid records
  hy_records = {}
  for varnum,var in enumerate(varlist):
    if not var.hasaxis(Hybrid): continue
    eta = var.getaxis(Hybrid)
    if 'ptop' not in eta.atts or 'rcoef' not in eta.atts or 'pref' not in eta.atts:
      warn ("Not enough information to construct an HY record");
      continue
    ptop = eta.atts['ptop']
    rcoef = eta.atts['rcoef']
    pref = eta.atts['pref']
    key = (ptop,rcoef,pref)
    if key not in hy_records:
      hy_record = fstd_core.make_hy_record(ptop,rcoef,pref)
      hy_record['etiket'] = var.atts.get('etiket','        ')
      # Use the same date as the field
      if var.hasaxis(Dateo):
        hy_record['dateo'] = var.getaxis(Dateo).values[0]
      # Use the same forecast time as the field
      if var.hasaxis(NPASAxis):
        npas_axis = var.getaxis(NPASAxis)
        deet = var.atts['deet']
        npas = npas_axis.values[0]
        hy_record['npas'] = npas
        hy_record['deet'] = deet
      hy_records[key] = hy_record

  if len(hy_records) > 1:
    warn ("Multiple Hybrid axes detected.  The resulting file may not work the way you expect.")

  vertical_records.extend(hy_records.values())

  # Check for log-hybrid levels
  bangbang_records = {}
  for varnum,var in enumerate(varlist):
    if not var.hasaxis(LogHybrid): continue
    zeta = var.getaxis(LogHybrid)
    required_atts = 'kind', 'version', 'ptop', 'pref', 'rcoef1', 'rcoef2', 'ref_name', 'ip1_m', 'a_m', 'b_m', 'ip1_t', 'a_t', 'b_t'
    if any(att not in zeta.atts for att in required_atts):
      warn ("Not enough information to construct a !! record");
      continue
    kind = zeta.atts['kind']
    version = zeta.atts['version']
    if (kind != 5 or version != 2):
      warn ("Only support vgrid kind=5, version=2.  Found: kind=%d, version=%d.  Not encoding !! record"%(kind,version))
      continue
    # Define a unique key for this vertical coordinate
    key = zeta.atts.copy()
    for a in 'ip1_m', 'a_m', 'b_m', 'ip1_t', 'a_t', 'b_t':
      key[a] = tuple(key[a])
    key = tuple(sorted(key.items()))
    if key not in bangbang_records:
      bangbang_records[key] = fstd_core.make_bangbang_record (zeta.atts)
      # Link to variable through IP1,IP2
      bangbang_records[key]['ip1'] = var.atts.get('ig1',0)
      bangbang_records[key]['ip2'] = var.atts.get('ig2',0)

  vertical_records.extend(bangbang_records.values())

  # Convert vertical axes to IP1Axis
  for varnum,var in enumerate(varlist):
    if var.hasaxis(ZAxis):
      zaxis = var.getaxis(ZAxis)
      if isinstance(zaxis,Height_wrt_SeaLevel): kind = 0
      elif isinstance(zaxis,Sigma): kind = 1
      elif isinstance(zaxis,Pres): kind = 2
      elif isinstance(zaxis,Height_wrt_Ground): kind = 4
      elif isinstance(zaxis,(Hybrid,LogHybrid)): kind = 5
      elif isinstance(zaxis,Theta): kind = 6
      else:
        kind = 3
        if zaxis.__class__ != ZAxis:
          from warnings import warn
          warn ("Vertical coordinate not recognized;  Encoding a generic Z axis")
      ip1 = fstd_core.encode_levels(zaxis.values,kind)
      ip1_axis = IP1Axis(values=ip1)
      varlist[varnum] = var.replace_axes({zaxis.name:ip1_axis})

  # Convert from list to array
  if len(vertical_records) > 0:
    vertical_records = np.concatenate(vertical_records)
  else:
    vertical_records = np.empty([0], dtype=fstd_core.record_descr)

  return vertical_records

# Encode latitude/longitude information into FSTD records
# Set up unique IG1,IG2,IG3,IG4 values for the grid
#TODO: support non-cartesian lat/lon projections
def encode_latlon (varlist):
  from pygeode.formats import fstd_core
  import numpy as np

  latlon_records = {}

  for i,var in enumerate(varlist):

    # Get the 1D lat/lon arrays
    lat = var.getaxis(Lat).values
    lon = var.getaxis(Lon).values

    xcoord = np.array(lon, dtype='float32')
    # Fix last longitude?
    if xcoord[-1] < xcoord[-2]:
      xcoord[-1] += 360
    ycoord = np.array(lat, dtype='float32')

    # Set grid IG1,IG2,IG3,IG4
    # (Hard-coded as non-rotated Z grid)
    grid_ig1 = 900
    grid_ig2 = 0
    grid_ig3 = 43200
    grid_ig4 = 43200

    # Create pseudo-unique grid identifiers
    ix = hash((tuple(xcoord.flatten()), tuple(ycoord.flatten()), grid_ig1, grid_ig2, grid_ig3, grid_ig4))
    ix1 = (ix%65536) + 32768
    ix2 = ((ix/65536)%65536) + 32768

    var.atts['ig1'] = ix1
    var.atts['ig2'] = ix2
    var.atts['grtyp'] = 'Z'

    # Add the horizontal records (if they haven't been handled yet)
    key = ('E', ix1, ix2)
    if key not in latlon_records:

      lat_record = np.zeros([1], dtype=fstd_core.record_descr)
      lat_record['nomvar'] = '^^  '
      lat_record['typvar'] = 'X '
      lat_record['etiket'] = var.atts.get('etiket','ETIKET      ')
      lat_record['grtyp'] = 'E'
      lat_record['ni'] = 1
      lat_record['nj'] = len(ycoord)
      lat_record['nk'] = 1
      if var.hasaxis(Dateo):
        lat_record['dateo'] = var.getaxis(Dateo).values[0]
      if var.hasaxis(NPASAxis):
        npas_axis = var.getaxis(NPASAxis)
        deet = var.atts['deet']
        npas = npas_axis.values[0]
        lat_record['npas'] = npas
        lat_record['deet'] = deet
      lat_record['datyp'] = 5
      lat_record['nbits'] = 32
      lat_record['ip1'] = ix1
      lat_record['ip2'] = ix2
      lat_record['ig1'] = grid_ig1
      lat_record['ig2'] = grid_ig2
      lat_record['ig3'] = grid_ig3
      lat_record['ig4'] = grid_ig4
      lat_record['data_func'] = lambda ycoord=ycoord: ycoord

      lon_record = np.zeros([1], dtype=fstd_core.record_descr)
      lon_record['nomvar'] = '>>  '
      lon_record['typvar'] = 'X '
      lon_record['etiket'] = var.atts.get('etiket','ETIKET      ')
      lon_record['grtyp'] = 'E'
      lon_record['ni'] = len(xcoord)
      lon_record['nj'] = 1
      lon_record['nk'] = 1
      if var.hasaxis(Dateo):
        lon_record['dateo'] = var.getaxis(Dateo).values[0]
      if var.hasaxis(NPASAxis):
        npas_axis = var.getaxis(NPASAxis)
        deet = var.atts['deet']
        npas = npas_axis.values[0]
        lon_record['npas'] = npas
        lon_record['deet'] = deet
      lon_record['datyp'] = 5
      lon_record['nbits'] = 32
      lon_record['ip1'] = ix1
      lon_record['ip2'] = ix2
      lon_record['ig1'] = grid_ig1
      lon_record['ig2'] = grid_ig2
      lon_record['ig3'] = grid_ig3
      lon_record['ig4'] = grid_ig4
      lon_record['data_func'] = lambda xcoord=xcoord: xcoord

      latlon_records[key] = (lon_record,lat_record)

  if len(latlon_records) == 0:
    return np.empty([0], dtype=fstd_core.record_descr)
  return np.concatenate(sum(latlon_records.values(),()))


# Coerce the variables into the expected FSTD dimensions
# Fill in any missing axes (where possible)
def normalize_fstd_axes (varlist):
  for i,var in enumerate(varlist):

    # Add any missing FSTD axes
    newaxes = []
    if not var.hasaxis(Dateo):
      newaxes.append(Dateo(values=[0]))
    if not var.hasaxis(NPASAxis):
      newaxes.append(NPASAxis(values=[0]))
    var.atts.setdefault('deet',1800) #TODO: better default
    if not var.hasaxis(IP1Axis):
      newaxes.append(IP1Axis(values=[0]))
    if not var.hasaxis(KAxis):
      newaxes.append(KAxis(values=[0]))
    if not var.hasaxis(YAxis):
      raise ValueError("missing y axis for '%s'"%var.name)
    if not var.hasaxis(XAxis):
      raise ValueError("missing x axis for '%s'"%var.name)
    if len(newaxes) > 0:
      var = var.extend(0,*newaxes)

    # Put the axes in the expected order
    var = var.transpose(Dateo,NPASAxis,IP1Axis,KAxis,YAxis,XAxis)
    assert var.naxes == 6

    varlist[i] = var

# Convert part of a variable to record data
def make_data_func (var,i_dateo,i_npas,i_ip1):
  ni = var.shape[-1]
  nj = var.shape[-2]
  nk = var.shape[-3]
  return lambda: var[i_dateo,i_npas,i_ip1,:,:,:].reshape(ni,nj,nk)

# Convert variables to FSTD records
def create_records (varlist):
  from pygeode.formats import fstd_core
  import numpy as np

  records = []
  for var in varlist:

    nomvar = (var.name.upper()+'    ')[:4]
    typvar = (var.atts.get('typvar','P')+'  ')[:2]
    etiket = (var.atts.get('etiket','ETIKET')+'            ')[:12]
    ni = var.shape[-1]
    nj = var.shape[-2]
    nk = var.shape[-3]
    deet = var.atts.get('deet',1800) #TODO: better default
    datyp = var.atts.get('datyp',5)
    if 'nbits' in var.atts:
      nbits = var.atts['nbits']
    else:
      nbits = {'float32':32,'float64':64}[var.dtype.name]
      datyp = 5
    grtyp = var.atts.get('grtyp','X')
    ig1 = var.atts.get('ig1',0)
    ig2 = var.atts.get('ig2',0)
    ig3 = var.atts.get('ig3',0)
    ig4 = var.atts.get('ig4',0)

    for i_dateo,dateo in enumerate(var.axes[0].values):
      for i_npas,npas in enumerate(var.axes[1].values):
        for i_ip1,ip1 in enumerate(var.axes[2].values):
          record = np.zeros([1], dtype=fstd_core.record_descr)
          record['nomvar'] = nomvar
          record['typvar'] = typvar
          record['etiket'] = etiket
          record['ni'] = ni
          record['nj'] = nj
          record['nk'] = nk
          record['dateo'] = dateo
          record['ip1'] = ip1
          record['ip2'] = npas*deet/3600
          record['deet'] = deet
          record['npas'] = npas
          record['datyp'] = datyp
          record['nbits'] = nbits
          record['grtyp'] = grtyp
          record['ig1'] = ig1
          record['ig2'] = ig2
          record['ig3'] = ig3
          record['ig4'] = ig4
          record['data_func'] = make_data_func (var, i_dateo, i_npas, i_ip1)

          records.append(record)

  if len(records) == 0: return np.empty([0], dtype=fstd_core.record_descr)
  return np.concatenate(records)

# Check for incompatible axes
def check_fstd_axes (varlist):
  compatible_axes = StandardTime, Forecast, IP1Axis, KAxis, JAxis, IAxis
  for var in varlist:
    incompatible_axes = [a for a in var.axes if not isinstance(a,compatible_axes)]
    if len(incompatible_axes) > 0:
      raise TypeError, "Cannot fit the following axes from var '%s' into an FSTD structure: %s"%(var.name,incompatible_axes)

#TODO: check for repeated axes


#####################################################################
#
# Save a dataset into an FSTD file.
#
#####################################################################
def prepare_records (varlist):
  from pygeode.dataset import Dataset
  from pygeode.formats import fstd_core
  import numpy as np

  if isinstance(varlist,Dataset):
    varlist = varlist.vars

  varlist = list(varlist)

  # Encode the time axes
  encode_time_axis (varlist)

  # Convert to FSTD-compatible axes
  # (e.g., detect hybrid / log-hybrid axes)
  detect_fstd_axes (varlist)

  # Encode the forecast axis
  encode_forecast_axis (varlist)

  # Extract horizontal information
  latlon_records = encode_latlon (varlist)

  # Extract vertical information
  vertical_records = encode_vertical(varlist)

  # We should now have a subset of (StandardTime,Forecast,IP1Axis,KAxis,JAxis,IAxis)

  # Fill in missing degenerate axes, and put them in the expected order
  normalize_fstd_axes (varlist)

  # Convert variables to record arrays
  var_records = create_records (varlist)

  # Merge coordinate records and variable records
  records = np.concatenate([latlon_records, vertical_records, var_records])

  return records

#####################################################################
#
# Save a dataset into an FSTD file.
#
#####################################################################
def save (filename, varlist):
  from os.path import exists
  from os import remove
  from pygeode.formats import fstd_core

  records = prepare_records(varlist)

  # Save to the file
  if exists(filename): remove(filename)
  fstd_core.write_records (filename, records)
