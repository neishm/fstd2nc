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

  # Helper method - iterate over all data functions.
  def _iter_data_funcs (self):
    from itertools import product
    # Get indices of all dimensions, in preparation for iterating.
    indices = [[slice(None)] if is_inner else range(s) for is_inner,s in zip(self._inner_dims, self.shape)]
    if hasattr(self._data_funcs,'flatten'):
      for i,ind in enumerate(product(*indices)):
        yield ind, self._data_funcs.flatten()[i]
    # data_funcs is degenerate (single element)?
    else:
      yield [slice(None)]*self.ndim, self._data_funcs()

  # Allow this object to be loaded into a numpy array.
  def __array__ (self):
    import numpy as np
    from itertools import product
    data = np.empty(self.shape, dtype=object)
    for ind,data_func in self._iter_data_funcs():
      data[ind] = data_func()
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
  _ignore_atts = ('swa','lng','dltf','ubc','xtra1','xtra2','xtra3')

  def __init__ (self):
    import numpy as np
    dtype = [("dateo", "i4"), ("deet", "i4"), ("npas", "i4"), ("ni", "i4"), ("nj", "i4"), ("nk", "i4"), ("nbits", "i4"), ("datyp", "i4"), ("ip1", "i4"), ("ip2", "i4"), ("ip3", "i4"), ("typvar", "a2"), ("nomvar", "a4"), ("etiket", "a12"), ("grtyp", "a2"), ("ig1", "i4"), ("ig2", "i4"), ("ig3", "i4"), ("ig4", "i4"), ("swa", "i4"), ("lng", "i4"), ("dltf", "i4"), ("ubc", "i4"), ("xtra1", "i4"), ("xtra2", "i4"), ("xtra3", "i4"), ("datev", "i4")]
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
    import numpy as np
    from rpnpy.librmn.fstd98 import fstinl, fstprm, fstluk
    from rpnpy.librmn.const import FST_RO
    # Read the data
    f = _FSTD_File(filename,FST_RO)
    keys = fstinl(f.funit)
    headers = np.zeros(len(keys),dtype=self._headers.dtype)
    prms = dict()
    data_funcs = []
    for i,key in enumerate(keys):
      prm = fstprm(key)
      for n,v in prm.iteritems():
        prms.setdefault(n,[]).append(v)
      def data_func (_f=f, _key=key):
        atts = fstluk(_key,rank=3)
        data = atts['d']
        # Transpose the axes to C-order.
        data = data.transpose(2,1,0)
        return data
      data_funcs.append(data_func)
    for n in headers.dtype.names:
      headers[n] = prms[n]
    self._headers = np.concatenate((self._headers,headers))
    self._data_funcs.extend(data_funcs)

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
    from rpnpy.rpndate import  RPNDate
    import numpy as np
    fields = super(_Dates,self)._get_fields()
    # Calculate the forecast (in hours).
    fields['forecast']=fields['deet']*fields['npas']/3600.
    # Time axis
    if self.squash_forecasts:
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
  def _put_fields (self, nrecs, fields):
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
    from rpnpy.librmn.fstd98 import DecodeIp
    import numpy as np
    fields = super(_VCoords,self)._get_fields()
    # Provide 'level' and 'kind' information to the decoder.
    decoded = map(DecodeIp,fields['ip1'],fields['ip2'],fields['ip3'])
    rp1 = zip(*decoded)[0]
    levels = np.array([r.v1 for r in rp1])
    kind = np.array([r.kind for r in rp1])
    # Only use first set of levels (can't handle ranges yet).
    fields['level'] = levels
    fields['kind'] = kind
    return fields
  def _put_fields (self,nrecs,fields):
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
    super(_VCoords,self)._put_fields(nrecs,fields)

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
    for header, data_func in zip(self._headers,self._data_funcs):
      if header['nomvar'].strip() not in self._vcoord_nomvars: continue
      key = (header['ip1'],header['ip2'])
      # For old HY records, there's no matching ipX/igX codes.
      if header['nomvar'].strip() == 'HY': key = 'HY'
      if key in vrecs: continue
      vrecs[key] = (header,data_func)

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
            header, data_func = vrecs[key]
            # Add in metadata from the coordinate.
            atts.update(self._get_header_atts(header))
            # Add type-specific metadata.
            if header['nomvar'].strip() == '!!':
              # Get A and B info.
              vgd_id = vgd_fromlist(data_func().transpose(2,1,0))
              if vgd_get (vgd_id,'LOGP'):
                name = 'zeta'
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
# Logic for reading/writing FSTD data from/to netCDF files.
# Note: this is not strictly an FSTD thing, but it's
# provided here for convenience.
class _netCDF_IO (_Buffer_Base):
  def __iter__ (self):
    from datetime import datetime
    from collections import OrderedDict
    import numpy as np
    from netCDF4 import date2num
    for var in super(_netCDF_IO,self).__iter__():
      # Modify time axes to be relative units instead of datetime objects.
      if var.name in var.axes and isinstance(var.array[0],datetime):
        units = 'hours since %s'%(var.array[0])
        var.atts.update(units=units)
        array = np.asarray(date2num(var.array,units=units))
        var = type(var)(var.name,var.atts,var.axes,array)
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
      # Write one FSTD record at a time to avoid out-of-memory errors.
      if isinstance(array,_Array):
        for ind, data_func in array._iter_data_funcs():
          v[ind] = data_func()
      # Otherwise, data is already in memory so no problems?
      else:
        v[:] = np.asarray(array[:])
    f.close()

# Default interface for I/O.
class Buffer (_netCDF_IO,_VCoords,_Dates):
  """
  High-level interface for FSTD data, to treat it as multi-dimensional arrays.
  Contains logic for dealing with most of the common FSTD file conventions.
  """


