#!/usr/bin/env python

"""
Functionality for converting between FSTD and netCDF files.
"""

# Helper classes for lazy-array evaluation.

# Common base class for all the numpy-like arrays defined below.
class _Array_Base (object):
  # Set some common attributes for the object.
  def __init__ (self, shape, dtype):
    # Expected shape and type of the array.
    self.shape = tuple(map(int,shape))
    self.ndim = len(self.shape)
    self.size = reduce(int.__mul__, self.shape, 1)
    self.dtype = dtype

# Data from a record.
# Has a 'shape' like a numpy array, but values aren't loaded from file until
# the array is sliced, or passed through np.asarray().
class _Record_Array (_Array_Base):
  def __init__ (self, fstdfile, params):
    from rpnpy.librmn.fstd98 import dtype_fst2numpy
    shape = (params['ni'],params['nj'],params['nk'])
    dtype = dtype_fst2numpy(params['datyp'], params['nbits'])
    _Array_Base.__init__(self, shape, dtype)
    # Keep a reference to the file so it doesn't get closed until after this
    # object is destroyed.
    # Otherwise, the key will become invalid.
    self._fstdfile = fstdfile
    self._key = params['key']
  def __getitem__ (self, key):
    return self.__array__().__getitem__(key)
  def __array__ (self):
    # Adapted from rpnpy.librmn.fstd98.fstluk
    import ctypes as _ct
    import numpy as _np
    from rpnpy.librmn import proto as _rp
    from rpnpy.librmn.fstd98 import FSTDError
    import numpy.ctypeslib as _npc
    (cni, cnj, cnk) = (_ct.c_int(), _ct.c_int(), _ct.c_int())
    data = _np.empty(self.shape, dtype=self.dtype, order='FORTRAN')
    _rp.c_fstluk.argtypes = (_npc.ndpointer(dtype=self.dtype), _ct.c_int,
                             _ct.POINTER(_ct.c_int), _ct.POINTER(_ct.c_int),
                             _ct.POINTER(_ct.c_int))
    istat = _rp.c_fstluk(data, self._key, _ct.byref(cni), _ct.byref(cnj),
                             _ct.byref(cnk))
    if istat < 0:
      raise FSTDError()
    return data

# Missing data.
class _NaN_Array (_Array_Base):
  def __init__ (self,*shape):
    _Array_Base.__init__(self,shape,'float32')
  def __getitem__ (self,key):
    return self.__array__().__getitem__(key)
  def __array__ (self):
    data = np.empty(self.shape,dtype=self.dtype,order='FORTRAN')
    data[()] = float('nan')
    return data

# Combine multiple array-like objects into a higher-dimensional array-like
# object.
# E.g., If you have a numpy array of dimensions (time,level), of type 'object',
# and each element is a _Record_Array object of dimenions (ni,nj,nk), then
# this will give you a new array-like object of dimenions (time,level,ni,nj,nk)
# that can be manipulated like a numpy array (e.g. sliced, tranposed).
# The values won't be loaded until the object is passed to numpy.asarray().
class _Array (_Array_Base):
  @classmethod
  # Create the object from an array of array-like objects.
  def create (cls, data):
    import numpy as np
    # Special case: have a single array (not an array of arrays).
    if data.dtype.name != 'object':
      outer_shape = ()
      inner_shape = data.shape
      dtype = data.dtype
    # Usual case (array of array-like objects).
    else:
      outer_shape = data.shape
      # Inner array objects must all have the same shape.
      inner_shape = set(map(np.shape,data.flatten()))
      if len(inner_shape) > 1:
        raise ValueError ("Different shapes for inner array objects.  Found shapes: %s"%list(inner_shape))
      inner_shape = inner_shape.pop()
      dtype = np.result_type(*data.flatten())
    shape = tuple(outer_shape) + tuple(inner_shape)
    # Define a map from outer axes to inner axes.
    inner_dimids = [None]*len(outer_shape) + range(len(inner_shape))
    inner_slices = [slice(None)]*len(inner_shape)
    return cls(shape, dtype, inner_dimids, data, inner_slices)
  def __init__ (self, shape, dtype, inner_dimids, data, inner_slices):
    assert len(shape) == len(inner_dimids)
    _Array_Base.__init__(self, shape, dtype)
    # Check if we have a dengenerate outer array.
    self._degenerate = not any(i is None for i in inner_dimids)
    # Check outer dimensions
    if not self._degenerate:
      assert sum(i is None for i in inner_dimids) == data.ndim
    # Map the full axes to the inner (record) axes.
    # Set to 'None' for outer axes not part of the record.
    self._inner_dimids = tuple(inner_dimids)
    # Array of references to the data from the FSTD records.
    self._data = data
    # Extra slicing to be done to the record data after reading it in.
    self._inner_slices = tuple(inner_slices)
  def __getitem__ (self, key):
    from itertools import product
    # Coerce key into a tuple of slice objects.
    if not isinstance(key,tuple):
      if hasattr(key,'__len__'): key = tuple(key)
      else: key = (key,)
    if len(key) == 1 and hasattr(key[0],'__len__'):
      key = tuple(key[0])
    if Ellipsis in key:
      i = key.index(Ellipsis)
      key = key[:i] + (slice(None),)*(self.ndim-len(key)+1) + key[i+1:]
    key = key + (slice(None),)*(self.ndim-len(key))
    if len(key) > self.ndim:
      raise ValueError("Too many dimensions for slicing.")
    shape = []
    inner_dimids = []
    outer_slices = []
    inner_slices = list(self._inner_slices)
    for sl,n,i in zip(key,self.shape,self._inner_dimids):
      # Only retain dimenions that aren't reduced out.
      if not isinstance(sl,int):
        shape.append(len(range(n)[sl]))
        inner_dimids.append(i)
      # Outer slice?
      if i is None:
        outer_slices.append(sl)
      # Inner slice?
      else:
        S = inner_slices[i]
        start = S.start or 0
        stop = S.stop
        step = S.step or 1
        if isinstance(sl,int):
          inner_slices[i] = start + step*sl
        else:
          if sl.stop is not None:
            stop = start + step*sl.stop
          start = start + step*(sl.start or 0)
          step = step * (sl.step or 1)
          inner_slices[i] = slice(start,stop,step)

    data = self._data.__getitem__(tuple(outer_slices))
    return _Array(shape, self.dtype, inner_dimids, data, inner_slices)

  # Remove degenerate dimensions from the array.
  def squeeze (self, axis=None):
    if isinstance(axis,int):
      axis = (axis,)
    if axis is None:
      axis = [i for i,s in enumerate(self.shape) if s == 1]
    key = [slice(None)]*self.ndim
    for a in axis:
      if self.shape[a] > 1:
        raise ValueError("Can only squeeze axes of length 1.")
      key[a] = 0
    return self.__getitem__(tuple(key))

  # Change the order of the dimensions
  def transpose (self, *axes):
    if len(axes) == 0:
      axes = tuple(range(self.ndim-1,-1,-1))
    if len(axes) == 1 and hasattr(axes[0],'__len__'):
      axes = tuple(axes[0])
    if len(axes) != self.ndim:
      raise ValueError("Wrong number of dimenions for transpose.")
    if sorted(axes) != range(self.ndim):
      raise ValueError("Bad axis arguments.")
    shape = [self.shape[a] for a in axes]
    inner_dimids = [self._inner_dimids[a] for a in axes]
    return _Array(shape, self.dtype, inner_dimids, self._data, self._inner_slices)

  # Allow this object to be loaded into a numpy array.
  def __array__ (self):
    import numpy as np
    from itertools import product
    # Final order of inner axes
    trans = [i for i in self._inner_dimids if i is not None]
    # Adjust the axis numbers in transpose operation
    # to account for any slicing that was done.
    for i,s in list(enumerate(self._inner_slices))[::-1]:
      if isinstance(s,int):
        trans = [t if t<i else t-1 for t in trans if t != i]
    # Outer array is degenerate?
    if self._degenerate:
      return np.asarray(self._data[self._inner_slices]).transpose(*trans)
    # Get indices of all dimensions, in preparation for iterating.
    indices = [range(s) if i is None else [slice(None)] for i,s in zip(self._inner_dimids, self.shape)]
    data = np.empty(self.shape, dtype=self.dtype)
    for i,ind in enumerate(product(*indices)):
      data[ind] = np.asarray(self._data.flatten()[i][self._inner_slices]).transpose(*trans)
    return data

# Helper class - keep an FSTD file open until all references are destroyed.
class _FSTD_File (object):
  def __init__ (self, filename, mode):
    import rpnpy.librmn.all as rmn
    from rpnpy.librmn.base import fnom
    from rpnpy.librmn.fstd98 import fstouv
    self.filename = filename
    self.funit = fnom(filename, mode)
    fstouv(self.funit, mode)
  def __del__ (self):
    from rpnpy.librmn.fstd98 import fstfrm
    from rpnpy.librmn.base import fclos
    istat = fstfrm(self.funit)
    istat = fclos(self.funit)


# Define a class for encoding / decoding FSTD data.
# Each step is placed in its own "mixin" class, to make it easier to patch in 
# new behaviour if more exotic FSTD files are encountered in the future.
class _Buffer_Base (object):

  # Names of records that should be kept separate (never grouped into
  # multidimensional arrays).
  _meta_records = ('HH', 'STNS', 'SH')

  # Attributes which could potentially be used as axes.
  _outer_axes = ()

  # Attributes which uniquely identify a variable.
  _var_id = ('nomvar','ni','nj','nk')

  # Record parameters which should not be used as nc variable attributes.
  # (They're either internal to the file, or part of the data, not metadata).
  _ignore_atts = ('swa','lng','dltf','ubc','xtra1','xtra2','xtra3','key','shape','d')

  # Define any command-line arguments for reading FSTD files.
  @classmethod
  def _input_args (cls, parser):
    pass

  def __init__ (self):
    """
    Create a new, empty buffer.
    """
    self._params = []

  # Extract metadata from a particular header.
  def _get_header_atts (self, header):
    for n,v in header.iteritems():
      if n in self._ignore_atts: continue
      if isinstance(v,str):
        v = v.strip()
      yield (n,v)

  def clear (self):
    """
    Removes all existing data from the buffer.
    """
    self._params = []


  ###############################################
  # Basic flow for reading data

  def read_file (self, filename):
    """
    Read raw records from an FSTD file, into the buffer.
    Multiple files can be read sequentially.
    """
    import numpy as np
    from rpnpy.librmn.fstd98 import fstinl, fstprm, fstluk
    from rpnpy.librmn.const import FST_RO
    # Read the data
    f = _FSTD_File(filename,FST_RO)
    keys = fstinl(f.funit)
    for i,key in enumerate(keys):
      prm = fstprm(key)
      prm['d'] = _Record_Array(f,prm)
      self._params.append(prm)

  # Collect the list of params from all the FSTD records, and concatenate them
  # into arrays.  Result is a single dictionary containing the vectorized
  # parameters of all records.
  # The sole purpose of this routine is to put the metiadata in a structure
  # that's more efficient for doing bulk (vectorized) manipulations, instead
  # of manipulating each record param dictionary one at a time.
  # Subclasses may also insert extra (non-FSTD) parameters in here which are
  # needed for doing the variable decoding.
  def _vectorize_params (self):
    from collections import OrderedDict
    import numpy as np
    # Make sure the parameter names are consistent for all records.
    if len(set(map(frozenset,self._params))) != 1:
      raise ValueError("Inconsistent parameter names for the records.")
    fields = OrderedDict()
    for prm in self._params:
      for n,v in prm.iteritems():
        fields.setdefault(n,[]).append(v)
    for n,v in fields.items():
      # Treat array-like objects specially (don't want numpy to try to
      # read the data).
      if hasattr(v[0],'__array__'):
        v2 = v
        v = np.empty(len(self._params),dtype='O')
        v[:] = v2
      fields[n] = np.asarray(v)
    return fields


  def __iter__ (self):
    """
    Processes the records into multidimensional variables.
    Iterates over (name, atts, axes, array) tuples.
    Note that array may not be a true numpy array (values are not yet loaded
    in memory).  To load the array, pass it to numpy.asarray().
    """
    from collections import OrderedDict, namedtuple
    import numpy as np

    # Degenerate case: no data in buffer
    if len(self._params) == 0: return

    records = self._vectorize_params()

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
      axes['i'] = range(var_id.ni)
      axes['j'] = range(var_id.nj)
      axes['k'] = range(var_id.nk)

      # Construct a multidimensional array to hold the data functions.
      data = np.empty(map(len,axes.values()[:-3]), dtype='O')

      # Assume missing data (nan) unless filled in later.
      data.fill(_NaN_Array(var_id.ni, var_id.nj, var_id.nk))
      
      # Arrange the data funcs in the appropriate locations.
      for rec_id in rec_ids:
        index = tuple(axes[n].index(records[n][rec_id]) for n in self._outer_axes)
        data[index] = records['d'][rec_id]

      # Check if we have full coverage along all axes.
      have_data = np.not_equal(data,None)
      if not np.all(have_data):
        from warnings import warn
        warn ("Missing some records for %s.")

      data = _Array.create (data)

      # Put the i,j,k dimensions in C order.
      data = data.transpose(range(data.ndim-3)+[data.ndim-1,data.ndim-2,data.ndim-3])
      axes = axes.items()
      axes = axes[:-3] + axes[-3:][::-1]
      axes = OrderedDict(axes)

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
    for var_id, atts, axes, array in data:
      for n,v in atts.iteritems():
        att_type[n] = type(v)
      for n,v in axes.iteritems():
        if n in ('k','j','i'): continue
        att_type[n] = type(v)
    att_type['d'] = object

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

      # Wrap the data array so we defer reading the values.
      data = []
      for ind in product(*map(range,array.shape[:-3])):
        data.append(_Array(array)[ind])

      current_records = OrderedDict()

      # Add coordinate info.
      for coords in product(*axes.values()[:-3]):
        for n,c in zip(axes.keys(),coords):
          current_records.setdefault(n,[]).append(c)

      # Add metadata.
      nrecs = len(data)
      for n,v in atts.iteritems():
        current_records[n] = [v]*nrecs

      # Add in the data references.
      current_records['d'] = data

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
    self._params.extend(list(self._unvectorize_params(nrecs, records)))


  # Take vectorized version of record param dictionary, and convert it back
  # to a list of individual record params (ready for writing to file).
  # Sub-classes may also do extra work here to re-encode FSTD params from
  # their own metadata.
  def _unvectorize_params (self, nrecs, fields):
    raise NotImplementedError #TODO
    import numpy as np
    # Pad out the string records with spaces
    fields['nomvar'] = [s.upper().ljust(4) for s in fields['nomvar']]
    fields['etiket'] = [s.upper().ljust(12) for s in fields['etiket']]
    fields['typvar'] = [s.upper().ljust(2) for s in fields['typvar']]
    fields['grtyp'] = map(str.upper, fields['grtyp'])

    params = []
    for i in range(nrecs):
      prm = dict()
      #TODO


  def write_file (self, filename):
    """
    Write the records from this object into the specified FSTD file.
    """
    raise NotImplementedError #TODO: update this code.
    from rpnpy.librmn.fstd98 import fstecr
    from rpnpy.librmn.const import FST_RW

    f = _FSTD_File(filename, FST_RW)
    for header, data_func in zip(self._headers, self._data_funcs):
      meta = dict([(n,header[n]) for n in header.dtype.names])
      data = data_func().transpose(2,1,0)  # Put back in Fortran order.
      fstecr(f.funit, data, meta)

  #
  ###############################################




#####################################################################
# Mixins for different features / behaviour for the conversions.


#################################################
# Selecting for particular fields.
class _SelectVars (_Buffer_Base):
  @classmethod
  def _input_args (cls, parser):
    super(_SelectVars,cls)._input_args(parser)
    parser.add_argument('--vars', metavar='VAR1,VAR2,...', help='Comma-seperated list of variables to convert.  By default, all variables are converted.')
  def __init__ (self, vars=None, *args, **kwargs):
    if vars is not None:
      self._selected_vars = vars.split(',')
      print 'Looking for variables:', ' '.join(self._selected_vars)
    else:
      self._selected_vars = None
    super(_SelectVars,self).__init__(*args,**kwargs)
  def __iter__ (self):
    for var in super(_SelectVars,self).__iter__():
      if self._selected_vars is not None:
        if var.name not in self._selected_vars:
          continue
      yield var


#################################################
# TODO: Logic for handling masks.
#
# Need some sample files first.
# In the old PyGeode-RPN code, I detected a mask by
# looking for fields with datyp=2 and nbits=1, and
# then pairing that with a data field with the same
# nomvar/etiket/etc.
# However, I can't find the original test file (RIOPS 2?)
# that I based this code on, and I don't know how
# general my assumptions were in the first place.


#################################################
# Logic for handling date field.

class _Dates (_Buffer_Base):
  @classmethod
  def _input_args (cls, parser):
    super(_Dates,cls)._input_args(parser)
    parser.add_argument('--squash-forecasts', action='store_true', help='Use the date of validity for the "time" axis.  Otherwise, the default is to use the date of original analysis, and the forecast length goes in a "forecast" axis.')

  def __init__ (self, squash_forecasts=False, *args, **kwargs):
    self._squash_forecasts = squash_forecasts
    if squash_forecasts:
      self._outer_axes = ('time',) + self._outer_axes
    else:
      self._outer_axes = ('time','forecast') + self._outer_axes
    super(_Dates,self).__init__(*args,**kwargs)

  # Get any extra (derived) fields needed for doing the decoding.
  def _vectorize_params (self):
    from rpnpy.rpndate import  RPNDate
    import numpy as np
    fields = super(_Dates,self)._vectorize_params()
    # Calculate the forecast (in hours).
    fields['forecast']=fields['deet']*fields['npas']/3600.
    # Time axis
    if self._squash_forecasts:
      dates = map(int,fields['datev'])
    else:
      dates = map(int,fields['dateo'])
    dates = [RPNDate(d).toDateTime().replace(tzinfo=None) if d > 0 else None for d in dates]
    dates = np.ma.masked_equal(dates,None)
    fields['time'] = dates
    return fields

  # Add time and forecast axes to the data stream.
  def __iter__ (self):
    from collections import OrderedDict
    import numpy as np
    # Keep track of all time and forecast axes found in the data.
    time_axes = set()
    forecast_axes = set()
    for var in super(_Dates,self).__iter__():
      if 'time' in var.axes:
        times = var.axes['time']
        if times not in time_axes:
          time_axes.add(times)
          atts = OrderedDict()
          axes = OrderedDict([('time',var.axes['time'])])
          # Add the time axis to the data stream.
          yield type(var)('time',atts,axes,np.asarray(times))
      if 'forecast' in var.axes:
        forecasts = var.axes['forecast']
        if forecasts not in forecast_axes:
          forecast_axes.add(forecasts)
          atts = OrderedDict(units='hours')
          axes = OrderedDict([('forecast',var.axes['forecast'])])
          # Add the forecast axis to the data stream.
          yield type(var)('forecast',atts,axes,np.asarray(forecasts))
      yield var

  # Prepare date info for output.
  def _unvectorize_params (self, nrecs, fields):
    import numpy as np
    from rpnpy.rpndate import RPNDate

    # Check for missing record attributes
    for att in 'deet','npas','dateo','forecast':
      if att not in fields:
        fields[att] = np.ma.masked_all(nrecs)
    if 'time' not in fields:
      fields['time'] = np.ma.masked_all(nrecs,dtype='O')

    # Mask out invalid deet values.
    fields['deet'] = np.ma.masked_equal(fields['deet'],0)

    # Set default values.
    fields['forecast'] = np.ma.filled(fields['forecast'],0)

    # Get dateo from 'time' field, wherever it's missing.
    data = [RPNDate(d).dateo if d is not None else 0 for d in fields['time']]
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

    return super(_Dates,self)._unvectorize_params(nrecs, fields)


#################################################
# Logic for handling vertical coordinates.

class _VCoords (_Buffer_Base):
  _vcoord_nomvars = ('HY','!!')
  def __init__ (self, *args, **kwargs):
    # Don't group records across different level 'kind'.
    # (otherwise can't create a coherent vertical axis).
    self._var_id = self._var_id + ('kind',)
    # Also, must have consistent igX records for a variable.
    if 'ig1' not in self._var_id:
      self._var_id = self._var_id + ('ig1','ig2','ig3','ig4')
    # Use decoded IP1 values as the vertical axis.
    self._outer_axes = ('level',) + self._outer_axes
    # Tell the decoder not to process vertical records as variables.
    self._meta_records = self._meta_records + self._vcoord_nomvars
    super(_VCoords,self).__init__(*args,**kwargs)
  def _vectorize_params (self):
    from rpnpy.librmn.fstd98 import DecodeIp
    import numpy as np
    fields = super(_VCoords,self)._vectorize_params()
    # Provide 'level' and 'kind' information to the decoder.
    decoded = map(DecodeIp,fields['ip1'],fields['ip2'],fields['ip3'])
    rp1 = zip(*decoded)[0]
    levels = np.array([r.v1 for r in rp1])
    kind = np.array([r.kind for r in rp1])
    # Only use first set of levels (can't handle ranges yet).
    fields['level'] = levels
    fields['kind'] = kind
    return fields
  def _unvectorize_params (self,nrecs,fields):
    import numpy as np
    from rpnpy.librmn.fstd98 import EncodeIp
    from rpnpy.librmn.proto import FLOAT_IP
    # Check for missing fields.
    for att in 'kind','level','ip1':
      if att not in fields:
        fields[att] = np.ma.masked_all(nrecs)
    # Re-encode ip1 values from kind,level information.
    have_data = ~np.ma.getmaskarray(fields['level']) & ~np.ma.getmaskarray(fields['kind'])
    level = fields['level'][have_data]
    kind = fields['kind'][have_data]
    # Doesn't handle level ranges.
    rp1 = [FLOAT_IP(z,z,k) for z,k in zip(level,kind)]
    rp2 = [FLOAT_IP(0.0,0.0,0)] * len(rp1)
    rp3 = rp2
    encoded = map(EncodeIp, rp1, rp2, rp3)
    ip1 = zip(*encoded)[0]
    fields['ip1'][have_data] = np.array(ip1)
    # Set default ip1 values.
    fields['ip1'] = np.ma.filled(fields['ip1'],0)
    # Continue re-encoding other parts of the field.
    return super(_VCoords,self)._unvectorize_params(nrecs,fields)

  # Add vertical axis as another variable.
  def __iter__ (self):
    from collections import OrderedDict
    import numpy as np
    from rpnpy.vgd.base import vgd_fromlist, vgd_get, vgd_free
    from rpnpy.vgd.const import VGD_KEYS
    from rpnpy.vgd import VGDError
    from rpnpy.librmn.fstd98 import DecodeIp
    # Pre-scan the raw headers for special vertical records.
    # (these aren't available in the data stream, because we told the decoder
    # to ignore them).
    vrecs = OrderedDict()
    for header in self._params:
      if header['nomvar'].strip() not in self._vcoord_nomvars: continue
      key = (header['ip1'],header['ip2'])
      # For old HY records, there's no matching ipX/igX codes.
      if header['nomvar'].strip() == 'HY': key = 'HY'
      if key in vrecs: continue
      vrecs[key] = header

    # Scan through the data, and look for any use of vertical coordinates.
    vaxes = OrderedDict()
    for var in super(_VCoords,self).__iter__():
      # Degenerate vertical axis?
      if 'ip1' in var.atts and var.atts['ip1'] == 0:
        del var.axes['level']
      if 'level' not in var.axes or 'kind' not in var.atts:
        yield var
        continue
      # Decode the vertical coordinate.
      levels = var.axes['level']
      kind = var.atts['kind']
      # Only need to provide one copy of the vertical axis.
      if (levels,kind) not in vaxes:
        # Keep track of any extra arrays needed for this axis.
        ancillary_variables = []
        # Get metadata that's specific to this axis.
        # Adapted from pygeode.formats.fstd, pygeode.axis, and
        # pygeode.format.cfmeta modules.
        name = 'zaxis'
        atts = OrderedDict()
        atts['axis'] = 'Z'
        if kind == 0:
          name = 'height'
          atts['standard_name'] = 'height_above_sea_level'
          atts['units'] = 'm'
          atts['positive'] = 'up'
        elif kind == 1:
          name = 'sigma'
          atts['standard_name'] = 'atmosphere_sigma_coordinate'
          atts['positive'] = 'down'
        elif kind == 2:
          name = 'pres'
          atts['standard_name'] = 'air_pressure'
          atts['units'] = 'hPa'
          atts['positive'] = 'down'
        elif kind == 3:
          name = 'lev'  # For backwards compatibility with PyGeode.
        elif kind == 4:
          name = 'height'
          atts['standard_name'] = 'height'
          atts['units'] = 'm'
          atts['positive'] = 'up'
        elif kind == 5:
          atts['positive'] = 'down'
          key = (var.atts['ig1'],var.atts['ig2'])
          if header['nomvar'].strip() == 'HY': key = 'HY'
          if key in vrecs:
            header = vrecs[key]
            # Add in metadata from the coordinate.
            atts.update(self._get_header_atts(header))
            # Add type-specific metadata.
            if header['nomvar'].strip() == '!!':
              # Get A and B info.
              vgd_id = vgd_fromlist(header['d'][...])
              if vgd_get (vgd_id,'LOGP'):
                name = 'zeta'
                # Not really a "standard" name, but used for backwards
                # compatibility with PyGeode.
                atts['standard_name'] = 'atmosphere_hybrid_sigma_log_pressure_coordinate'
              else:
                name = 'eta'
                atts['standard_name'] = 'atmosphere_hybrid_sigma_pressure_coordinate'
              # Add all parameters for this coordinate.
              for key in VGD_KEYS:
                try:
                  val = vgd_get(vgd_id,key)
                  # Skip multidimensional arrays (can't be encoded as metadata).
                  if getattr(val,'ndim',1) > 1: continue
                  atts[key] = val
                except (KeyError,VGDError):
                  pass  # Some keys not available in some vgrids?
              # Some attribute aliases that are needed for reverse-compatibility
              # with old PyGeode.formats.fstd
              aliases = OrderedDict([('CA_M','a_m'),('CA_T','a_t'),('CB_M','b_m'),('CB_T','b_t'),('VIPM','ip1_m'),('VIPT','ip1_t'),('VERS','version'),('RC_1','rcoef1'),('RC_2','rcoef2'),('RFLD','ref_name')])
              for oldname,newname in aliases.iteritems():
                if oldname in atts: atts[newname] = atts[oldname]
              # Attempt to fill in A/B ancillary data (if available).
              try:
                all_z = list(atts['VCDM'])+list(atts['VCDT'])
                all_a = list(atts['CA_M'])+list(atts['CA_T'])
                all_b = list(atts['CB_M'])+list(atts['CB_T'])
                A = []
                B = []
                for z in levels:
                  ind = all_z.index(z)
                  A.append(all_a[ind])
                  B.append(all_b[ind])
                A = type(var)(name+'_A', {}, {name:levels}, np.asarray(A))
                B = type(var)(name+'_B', {}, {name:levels}, np.asarray(B))
                ancillary_variables.extend([A,B])
              except (KeyError,ValueError,VGDError):
                from warnings import warn
                warn ("Unable to get A/B coefficients.")
              vgd_free (vgd_id)
            # Not a '!!' coordinate, so must be 'HY'?
            else:
              name = 'eta'
              atts['standard_name'] = 'atmosphere_hybrid_sigma_pressure_coordinate'
              # Get A and B info.
              eta = np.asarray(levels)
              ptop = DecodeIp(header['ip1'],0,0)[0].v1
              # Conversion taken from 'ig_to_hybref' function in librmn:
              pref = float(header['ig1'])
              rcoef = header['ig2']/1000.0
              # Apply the formula to compue A & B (from old fstd_core.c code):
              etatop = ptop/pref
              B = ((eta - etatop) / (1 - etatop)) ** rcoef
              A = pref * 100. * (eta - B)
              B = type(var)(name+'_B', {}, {name:levels}, B)
              A = type(var)(name+'_A', {}, {name:levels}, A)
              ancillary_variables.extend([A,B])
              # Add extra HY record metadata.
              atts.update(ptop=ptop, rcoef=rcoef, pref=pref)
        elif kind == 6:
          name = 'theta'
          atts['standard_name'] = 'air_potential_temperature'
          atts['units'] = 'K'
          atts['positive'] = 'up'

        # Add this vertical axis.
        axes = OrderedDict([(name,levels)])
        if len(ancillary_variables) > 0:
          atts['ancillary_variables'] = ' '.join(v.name for v in ancillary_variables)
        array = np.asarray(levels)
        vaxes[(levels,kind)] = type(var)(name,atts,axes,array)
        yield vaxes[(levels,kind)]
        # Add any ancillary data needed for the axis.
        for anc in ancillary_variables:
          yield anc
      # Get the vertical axis.
      vaxis = vaxes[(levels,kind)]
      # Modify the variable's dimension name to match the axis name.
      axes = OrderedDict(((vaxis.name if n == 'level' else n),v) for n,v in var.axes.iteritems())
      var = type(var)(var.name,var.atts,axes,var.array)
      yield var


#################################################
# Logic for handling lat/lon coordinates.

class _XYCoords (_Buffer_Base):
  _xycoord_nomvars = ('^^','>>')
  def __init__ (self, *args, **kwargs):
    # Variables must have an internally consistent horizontal grid.
    self._var_id = self._var_id + ('grtyp',)
    # Also, must have consistent igX records for a variable.
    if 'ig1' not in self._var_id:
      self._var_id = self._var_id + ('ig1','ig2','ig3','ig4')
    # Tell the decoder not to process horizontal records as variables.
    self._meta_records = self._meta_records + self._xycoord_nomvars
    super(_XYCoords,self).__init__(*args,**kwargs)
  # Add horizontal coordinate info to the data stream.
  def __iter__ (self):
    from collections import OrderedDict
    from rpnpy.librmn.interp import ezgdef_fmem, gdll, EzscintError
    import numpy as np

    # Scan through the data, and look for any use of horizontal coordinates.
    latlon = OrderedDict()
    for var in super(_XYCoords,self).__iter__():
      # Don't touch variables that have no horizontal extent.
      if 'i' not in var.axes or 'j' not in var.axes:
        yield var
        continue
      # Get basic information about this grid.
      gridinfo = OrderedDict()
      for n in ('ni','nj','grtyp','ig1','ig2','ig3','ig4'):
        v = var.atts[n]
        if isinstance(v,str):
          gridinfo[n] = v
        else:
          gridinfo[n] = int(v)  # ezgdef_fmem is very picky about types.
      # Check if we already defined this grid.
      key = tuple(gridinfo.values())
      if key not in latlon:
        # Remember the associated '>>','^^' metadata for later.
        xatts = OrderedDict()
        yatts = OrderedDict()
        # Check for reference grid data.
        for header in self._params:
          nomvar = header['nomvar'].strip()
          if nomvar not in self._xycoord_nomvars: continue
          k1 = (header['ip1'],header['ip2'],header['ip3'])
          k2 = (var.atts['ig1'],var.atts['ig2'],var.atts['ig3'])
          if k1 == k2:
            gridinfo['grref'] = header['grtyp'].strip()
            gridinfo['ig1'] = int(header['ig1'])
            gridinfo['ig2'] = int(header['ig2'])
            gridinfo['ig3'] = int(header['ig3'])
            gridinfo['ig4'] = int(header['ig4'])
            if nomvar == '>>':
              gridinfo['ax'] = header['d'][...]
              xatts = OrderedDict(self._get_header_atts(header))
            if nomvar == '^^':
              gridinfo['ay'] = header['d'][...]
              yatts = OrderedDict(self._get_header_atts(header))
        try:
          # Get the lat & lon data.
          gdid = ezgdef_fmem (**gridinfo)
          ll = gdll(gdid)
        except (TypeError,EzscintError):
          from warnings import warn
          warn("Unable to get grid info for '%s'"%var.name)
          yield var
          continue
        latarray = ll['lat'].transpose() # Switch from Fortran to C order.
        latatts = OrderedDict()
        latatts['long_name'] = 'latitude'
        latatts['standard_name'] = 'latitude'
        latatts['units'] = 'degrees_north'
        lonarray = ll['lon'].transpose() # Switch from Fortran to C order.
        lonatts = OrderedDict()
        lonatts['long_name'] = 'longitude'
        lonatts['standard_name'] = 'longitude'
        lonatts['units'] = 'degrees_east'
        axes = OrderedDict([('j',var.axes['j']),('i',var.axes['i'])])
        if 'ax' in gridinfo and 'ay' in gridinfo:
          axes = OrderedDict([('y',tuple(gridinfo['ay'].flatten())),('x',tuple(gridinfo['ax'].flatten()))])
        lat = type(var)('lat',latatts,axes,latarray)
        lon = type(var)('lon',lonatts,axes,lonarray)
        # Try to resolve lat/lon to 1D Cartesian coordinates, if possible.
        # Calculate the mean lat/lon arrays in double precision.
        meanlat = np.mean(np.array(latarray,dtype=float),axis=1,keepdims=True)
        meanlon = np.mean(np.array(lonarray,dtype=float),axis=0,keepdims=True)
        if np.allclose(latarray,meanlat) and np.allclose(lonarray,meanlon):
          # Reduce back to single precision for writing out.
          meanlat = np.array(meanlat,dtype=latarray.dtype).squeeze()
          meanlon = np.array(meanlon,dtype=lonarray.dtype).squeeze()
          # Ensure monotonicity of longitude field.
          # (gdll may sometimes wrap last longitude to zero).
          # Taken from old fstd_core.c code.
          if meanlon[-2] > meanlon[-3] and meanlon[-1] < meanlon[-2]:
            meanlon[-1] += 360.
          latatts.update(yatts)
          lonatts.update(xatts)
          lat = type(var)('lat',latatts,{'lat':tuple(meanlat)},meanlat)
          lon = type(var)('lon',lonatts,{'lon':tuple(meanlon)},meanlon)
        # Add x and y as variables, with the coord values and header metadata.
        if 'x' in lat.axes:
          yield type(var)('x',xatts,{'x':axes['x']},np.array(axes['x']))
          yield type(var)('y',yatts,{'y':axes['y']},np.array(axes['y']))
        yield lat
        yield lon
        latlon[key] = (lat,lon)
      lat,lon = latlon[key]
      # Update the var's horizontal coordinates.
      axes = OrderedDict()
      for axisname,axisvalues in var.axes.iteritems():
        if axisname == 'j':
          axisname,axisvalues = lat.axes.items()[0]
        elif axisname == 'i':
          axisname,axisvalues = lon.axes.items()[-1]
        axes[axisname] = axisvalues
      # For 2D lat/lon, add these to the metadata.
      if 'lat' not in axes or 'lon' not in axes:
        var.atts['coordinates'] = 'lon lat'

      yield type(var)(var.name,var.atts,axes,var.array)


#################################################
# Remove extraneous dimensions from the output.

class _NoNK (_Buffer_Base):
  def __iter__ (self):
    for var in super(_NoNK,self).__iter__():
      axes = var.axes
      array = var.array
      if 'k' in axes and len(axes['k']) == 1:
        array = array.squeeze(axis=axes.keys().index('k'))
        del axes['k']
      yield type(var)(var.name,var.atts,axes,array)


#################################################
# Logic for reading/writing FSTD data from/to netCDF files.

class _netCDF_IO (_Buffer_Base):
  @classmethod
  def _input_args (cls, parser):
    super(_netCDF_IO,cls)._input_args(parser)
    parser.add_argument('--time-units', choices=['seconds','minutes','hours','days'], default='hours', help='The units of time for the netCDF file.  Default is %(default)s.')
    parser.add_argument('--buffer_size', type=int, default=1, help='How much data to write at a time (in MBytes).  Default is %(default)s.')

  def __init__ (self, time_units='hours', buffer_size=1, *args, **kwargs):
    self._time_units = time_units
    self._buffer_size = int(buffer_size)
    super(_netCDF_IO,self).__init__(*args,**kwargs)

  def __iter__ (self):
    from datetime import datetime
    from collections import Counter, OrderedDict
    import numpy as np
    from netCDF4 import date2num
    # Keep track of all used variable / dimension names
    varcount = Counter()
    dimcount = Counter()
    # Mappings from original dimension names to unique names
    dimmap = dict()

    for var in super(_netCDF_IO,self).__iter__():

      # Modify time axes to be relative units instead of datetime objects.
      if var.name in var.axes and isinstance(var.array[0],datetime):
        units = '%s since %s'%(self._time_units,var.array[0])
        var.atts.update(units=units)
        array = np.asarray(date2num(var.array,units=units))
        var = type(var)(var.name,var.atts,var.axes,array)

      # Make variable / dimension names unique, by appending integer suffices
      # when necessary.
      # Also check for illegal characters in names.
      varname = var.name
      # Name must start with alphanumeric, or underscore.
      if not varname[0].isalnum():
        if not varname.startswith('_'):
          varname = '_'+varname
      varcount.update([varname])
      # Do we need to add an integer suffix?
      num = varcount.get(varname)
      if num > 1:
        varname = var.name+str(num)
      # Do we need to clobber the dimension names?
      axes = []
      for name,values in var.axes.iteritems():
        values = tuple(values)  # Values need to be hashable for lookup table.
        # Name must start with alphanumeric, or underscore.
        if not name[0].isalnum():
          if not name.startswith('_'):
            name = '_'+name
        # Append integer suffix to make unique dimension name?
        if (name,values) not in dimmap:
          dimcount.update([name])
          num = dimcount.get(name)
          if num > 1:
            dimmap[(name,values)] = name+str(num)
          else:
            dimmap[(name,values)] = name
        name = dimmap[(name,values)]
        axes.append((name,values))
      axes = OrderedDict(axes)
      var = type(var)(varname,var.atts,axes,var.array)
      # (end of name checks)

      yield var


  def write_nc_file (self, filename):
    """
    Write the records to a netCDF file.
    Requires the netCDF4 package.
    """
    from netCDF4 import Dataset
    import numpy as np
    from itertools import product
    f = Dataset(filename, "w", format="NETCDF4")

    for varname, atts, axes, array in iter(self):
      for axisname, axisvalues in axes.items():
        # Only need to create each dimension once (even if it's in multiple
        # variables).
        if axisname not in f.dimensions:
          f.createDimension(axisname, len(axisvalues))
      # Write the variable.
      v = f.createVariable(varname, datatype=array.dtype, dimensions=axes.keys())
      v.setncatts(atts)
      # Don't write too much at a time.
      a = 0
      check = array.size
      while check > self._buffer_size*1E6:
        check /= array.shape[a]
        a = a + 1
      for ind in product(*map(range,array.shape[:a])):
        v[ind] = np.asarray(array[ind])
    f.close()


# Default interface for I/O.
class Buffer (_netCDF_IO,_NoNK,_XYCoords,_VCoords,_Dates,_SelectVars):
  """
  High-level interface for FSTD data, to treat it as multi-dimensional arrays.
  Contains logic for dealing with most of the common FSTD file conventions.
  """


# Command-line invocation:
def _fstd2nc (buffer_type):
  from argparse import ArgumentParser
  parser = ArgumentParser(description="Converts an RPN standard file (FSTD) to netCDF format.")
  parser.add_argument('infile', metavar='<fstd_file>', help='The FSTD file to convert.')
  parser.add_argument('outfile', metavar='<netcdf_file>', help='The name of the netCDF file to create.')
  buffer_type._input_args(parser)
  args = parser.parse_args()
  args = vars(args)
  infile = args.pop('infile')
  outfile = args.pop('outfile')
  buf = buffer_type(**args)
  buf.read_file(infile)
  buf.write_nc_file(outfile)

if __name__ == '__main__':
  _fstd2nc (buffer_type=Buffer)

