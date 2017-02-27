"""
High-level interface for FSTD data.
"""

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

  # Names of records that should be kept separate (never grouped into
  # multidimensional arrays).
  _meta_records = ('>>', '^^', 'HH', 'STNS', 'SH')

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
    self._headers = np.array([],dtype=dtype)
    self._data_funcs = []

  # Extract metadata from a particular header.
  def _get_header_atts (self, header):
    for n in header.dtype.names:
      if n in self._ignore_atts: continue
      v = header[n]
      if isinstance(v,str):
        v = v.strip()
      yield (n,v)

  def clear (self):
    """
    Removes all existing data from the buffer.
    """
    import numpy as np
    self._headers = np.array([],dtype=self._headers.dtype)
    self._data_funcs = []


  ###############################################
  # Basic flow for reading data

  def read_file (self, filename):
    """
    Read raw records from an FSTD file, into the buffer.
    Multiple files can be read sequentially.
    """
    import fstd_core
    import numpy as np
    # Read the data
    records = fstd_core.read_records(filename)
    # Store the headers and data interfaces.
    headers = np.zeros(len(records),dtype=self._headers.dtype)
    for n in headers.dtype.names:
      headers[n] = records[n]
    self._headers = np.concatenate((self._headers,headers))
    self._data_funcs.extend(records['data_func'])


  # Decode the record headers into a dictionary, and
  # add extra (derived) attributes that are useful for
  # processing the data into multidimensional arrays.
  def _get_fields (self):
    from collections import OrderedDict
    import numpy as np
    fields = OrderedDict()
    for n in self._headers.dtype.names:
      fields[n] = self._headers[n]
    fields['data_func'] = np.asarray(self._data_funcs)
    return fields


  def __iter__ (self):
    """
    Processes the records into multidimensional variables.
    Iterates over (name, atts, axes, array) tuples.
    Note that array may not be a true numpy array (values are not yet loaded
    in memory).  To load the array, pass it to numpy.array().
    """
    import fstd_core
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

    headers = np.zeros(nrecs,dtype=self._headers.dtype)
    for n in headers.dtype.names:
      if n in fields:
        headers[n] = fields[n]
    self._headers = np.concatenate((self._headers,headers))
    self._data_funcs.extend(fields['data_func'])



  def write_file (self, filename):
    """
    Write the records from this object into the specified FSTD file.
    """
    import fstd_core
    import numpy as np

    # Create a numpy structured array to hold the data.
    records = np.zeros(len(self._headers),dtype=fstd_core.record_descr)

    # Fill in what we have from our existing records.
    for n in fstd_core.record_descr.names:
      if n in self._headers.dtype.names:
        records[n] = self._headers[n]
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
    self.squash_forecasts = squash_forecasts
    if squash_forecasts:
      self._outer_axes = ('time',) + self._outer_axes
    else:
      self._outer_axes = ('time','forecast') + self._outer_axes
    super(_Dates,self).__init__(*args,**kwargs)

  # Get any extra (derived) fields needed for doing the decoding.
  def _get_fields (self):
    import fstd_core
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

  # Add a time axis to the data stream.
  def __iter__ (self):
    from collections import OrderedDict
    import numpy as np
    # Keep track of all time axes found in the data.
    time_axes = set()
    for var in super(_Dates,self).__iter__():
      if 'time' in var.axes:
        times = var.axes['time']
        if times not in time_axes:
          time_axes.add(times)
          atts = OrderedDict()
          axes = OrderedDict([('time',var.axes['time'])])
          # Add the time axis to the data stream.
          yield type(var)('time',atts,axes,np.asarray(times))
      yield var

  # Prepare date info for output.
  def _put_fields (self, nrecs, fields):
    import fstd_core
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
    self._outer_axes += ('level',)
    # Tell the decoder not to process vertical records as variables.
    self._meta_records = self._meta_records + self._vcoord_nomvars
    super(_VCoords,self).__init__(*args,**kwargs)
  def _get_fields (self):
    import fstd_core
    fields = super(_VCoords,self)._get_fields()
    # Provide 'level' and 'kind' information to the decoder.
    levels, kind = fstd_core.decode_levels(fields['ip1'])
    fields['level'] = levels
    fields['kind'] = kind
    return fields
  def _put_fields (self,nrecs,fields):
    import numpy as np
    import fstd_core
    # Check for missing fields.
    for att in 'kind','level','ip1':
      if att not in fields:
        fields[att] = np.ma.masked_all(nrecs)
    # Re-encode ip1 values from kind,level information.
    have_data = ~np.ma.getmaskarray(fields['level']) & ~np.ma.getmaskarray(fields['kind'])
    level = fields['level'][have_data]
    kind = fields['kind'][have_data]
    fields['ip1'][have_data] = fstd_core.encode_levels(level,kind)
    # Set default ip1 values.
    fields['ip1'] = np.ma.filled(fields['ip1'],0)
    # Continue re-encoding other parts of the field.
    super(_VCoords,self)._put_fields(nrecs,fields)

  # Add vertical axis as another variable.
  def __iter__ (self):
    from collections import OrderedDict
    import numpy as np
    # Pre-scan the raw headers for special vertical records.
    # (these aren't available in the data stream, because we told the decoder
    # to ignore them).
    vrecs = OrderedDict()
    for header in self._headers:
      if header['nomvar'].strip() not in self._vcoord_nomvars: continue
      key = (header['ip1'],header['ip2'])
      if key in vrecs: continue
      vrecs[key] = header

    # Scan through the data, and look for any use of vertical coordinates.
    vaxes = OrderedDict()
    for var in super(_VCoords,self).__iter__():
      # Check if this variable uses a vertical coordinate, and determine
      # what that coordinate is.
      if 'level' not in var.axes or 'kind' not in var.atts:
        yield var
        continue
      levels = var.axes['level']
      kind = var.atts['kind']
      # Only need to provide one copy of the vertical axis.
      if (levels,kind) not in vaxes:
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
          if key in vrecs:
            header = vrecs[key]
            # Add in metadata from the coordinate.
            atts.update(self._get_header_atts(header))
            # Add type-specific metadata.
            #TODO: more robust check.
            if header['nomvar'].strip() == '!!':
              name = 'zeta'
              atts['standard_name'] = 'atmosphere_hybrid_sigma_log_pressure_coordinate'
            else:
              name = 'eta'
              atts['standard_name'] = 'atmosphere_hybrid_sigma_pressure_coordinate'
        # Add this vertical axis.
        axes = OrderedDict([(name,levels)])
        array = np.asarray(levels)
        vaxes[(levels,kind)] = type(var)(name,atts,axes,array)
        yield vaxes[(levels,kind)]
      # Get the vertical axis
      vaxis = vaxes[(levels,kind)]
      # Modify the dimension name to match the axis name.
      axes = OrderedDict(((vaxis.name if n == 'level' else n),v) for n,v in var.axes.iteritems())
      var = type(var)(var.name,var.atts,axes,var.array)
      yield var



#################################################
# Logic for reading/writing FSTD data from/to netCDF files.
# Note: this is not strictly an FSTD thing, but it's
# provided here for convenience.
class _netCDF_IO (_Buffer_Base):
  def __iter__ (self):
    from datetime import datetime
    from collections import OrderedDict
    import numpy as np
    from netCDF4 import date2num
    # Keep track of all time axes found in the data.
    for var in super(_netCDF_IO,self).__iter__():
      # Modify time axes to be relative units instead of datetime objects.
      if var.name == 'time':
        units = 'hours since %s'%(var.array[0])
        var.atts.update(units=units)
        array = np.asarray(date2num(var.array,units=units))
        var = type(var)('time',var.atts,var.axes,array)
      yield var


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
class Buffer (_netCDF_IO,_VCoords,_Dates):
  """
  High-level interface for FSTD data, to treat it as multi-dimensional arrays.
  Contains logic for dealing with most of the common FSTD file conventions.
  """


