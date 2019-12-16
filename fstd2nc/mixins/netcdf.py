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

from fstd2nc.stdout import _, info, warn, error
from fstd2nc.mixins import BufferBase

#################################################
# Mixin for adding netCDF metadata to the variables

class netCDF_Atts (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    import argparse
    super(netCDF_Atts,cls)._cmdline_args(parser)
    parser.add_argument('--metadata-file', type=argparse.FileType('r'), action='append', help=_('Use metadata from the specified file.  You can repeat this option multiple times to build metadata from different sources.'))
    parser.add_argument('--rename', metavar="OLDNAME=NEWNAME,...", help=_('Apply the specified name changes to the variables.'))
  @classmethod
  def _check_args (cls, parser, args):
    super(netCDF_Atts,cls)._check_args(parser,args)
    if args.rename is not None:
      try:
        dict(r.split('=') for r in args.rename.split(','))
      except ValueError:
        parser.error(_("Unable to parse the rename arguments."))

  def __init__ (self, *args, **kwargs):
    """
    metadata_file : str or list, optional
        Use metadata from the specified file(s).
    rename : str or dict, optional
        Apply the specified name changes to the variables.
    """
    try:
      import ConfigParser
    except ImportError:
      import configparser as ConfigParser
    from collections import OrderedDict
    metadata_file = kwargs.pop('metadata_file',None)
    if metadata_file is None:
      metafiles = []
    elif isinstance(metadata_file,str):
      metafiles = [metadata_file]
    else:
      metafiles = metadata_file
    # Open the files if only the filename is provided.
    metafiles = [open(m,'r') if isinstance(m,str) else m for m in metafiles]
    metadata = OrderedDict()
    # Set some global defaults.
    # We need to explicitly state that we're using CF conventions in our
    # output files, or some utilities (like IDV) won't accept the data.
    metadata['global'] = OrderedDict(Conventions = "CF-1.6")

    configparser = ConfigParser.SafeConfigParser()
    for metafile in metafiles:
      configparser.readfp(metafile)
    for varname in configparser.sections():
      metadata.setdefault(varname,OrderedDict()).update(configparser.items(varname))
      # Detect numerical values
      for k,v in list(metadata[varname].items()):
        try:
          metadata[varname][k] = float(v) # First, try converting to float.
          metadata[varname][k] = int(v) # Try further conversion to int.
        except ValueError: pass
    self._metadata = metadata
    # Check for renames.
    # Will override any renames specified in the metadata file.
    rename = kwargs.pop('rename',None)
    if rename is None:
      rename = {}
    if isinstance(rename,str):
      rename = [r.split('=') for r in rename.split(',')]
      rename = [(k.strip(),v.strip()) for k,v in rename]
      rename = dict(rename)
    for oldname, newname in rename.items():
      self._metadata.setdefault(oldname,OrderedDict())['rename'] = newname
    super(netCDF_Atts,self).__init__(*args,**kwargs)

  def _makevars (self):

    super(netCDF_Atts,self)._makevars()

    # Extract variable rename requests from the user-supplied metadata.
    renames = {}
    for varname, atts in self._metadata.items():
      if 'rename' in atts:
        renames[varname] = atts.pop('rename')

    # Apply the user-supplied metadata.
    for var in self._varlist:
      # Add extra metadata provided by the user?
      if var.name in self._metadata:
        var.atts.update(self._metadata[var.name])

    # Apply renames.
    for obj in self._iter_objects():
      if obj.name in renames:
        obj.name = renames[obj.name]



#################################################
# Mixin for reading/writing FSTD data from/to netCDF files.

class netCDF_IO (BufferBase):
  @classmethod
  def _cmdline_args (cls, parser):
    super(netCDF_IO,cls)._cmdline_args(parser)
    parser.add_argument('--time-units', choices=['seconds','minutes','hours','days'], default='hours', help=_('The units for the output time axis.  Default is %(default)s.'))
    parser.add_argument('--reference-date', metavar=_('YYYY-MM-DD'), help=_('The reference date for the output time axis.  The default is the starting date in the RPN file.'))

  @classmethod
  def _check_args (cls, parser, args):
    from datetime import datetime
    super(netCDF_IO,cls)._check_args(parser,args)
    # Parse the reference date into a datetime object.
    if args.reference_date is not None:
      try:
        datetime.strptime(args.reference_date,'%Y-%m-%d')
      except ValueError:
        parser.error(_("Unable to to parse the reference date '%s'.  Expected format is '%s'")%(args.reference_date,_('YYYY-MM-DD')))

  def __init__ (self, *args, **kwargs):
    self._time_units = kwargs.pop('time_units','hours')
    self._reference_date = kwargs.pop('reference_date',None)
    self._unique_names = kwargs.pop('unique_names',True)
    super(netCDF_IO,self).__init__(*args,**kwargs)

  def _makevars (self):
    from fstd2nc.mixins import _var_type
    from datetime import datetime
    import numpy as np
    from netCDF4 import date2num

    super(netCDF_IO,self)._makevars()

    if self._reference_date is None:
      reference_date = None
    else:
      reference_date = datetime.strptime(self._reference_date,'%Y-%m-%d')

    # Check if time axis can be made an unlimited dimension.
    # Can only work if it only appears as the outermost dimension, otherwise
    # the netCDF4 module will crash (maybe a bug with netCDF4?)
    self._time_unlimited = True
    for var in self._varlist:
      if 'time' not in var.dims: continue
      if var.dims.index('time') > 0:
        self._time_unlimited = False

    if self._unique_names:
      self._fix_names()

    for var in self._iter_objects():
      # Modify time axes to be relative units instead of datetime objects.
      # Also attach relevant metadata.
      if hasattr(var,'array') and isinstance(var.array.reshape(-1)[0],np.datetime64):
        # Convert from np.datetime64 to datetime.datetime
        # .tolist() only returns a datetime object for datetime64[s], not
        # for datetime64[ns].
        var.array = np.asarray(var.array, dtype='datetime64[s]')
        # https://stackoverflow.com/a/13703930/9947646
        var.array = np.array(var.array.tolist())
        units = '%s since %s'%(self._time_units, reference_date or var.array.reshape(-1)[0])
        var.atts.update(units=units, calendar='gregorian')
        var.array = np.asarray(date2num(var.array,units=units))

    for obj in self._iter_objects():
      # Encode the attributes so they're ready for writing to netCDF.
      # Handles things like encoding coordinate objects to a string.
      if hasattr(obj,'atts'):
        self._encode_atts(obj)


  # Helper method - prepare attributes for writing to netCDF.
  @staticmethod
  def _encode_atts (obj):
    from collections import OrderedDict
    atts = obj.atts
    for attname, attval in list(atts.items()):
      # Detect veriable reference, convert to string.
      if hasattr(attval,'name'):
        # Store the dependency for future use.
        obj.deps.append(attval)
        atts[attname] = attval.name
      # Detect list of objects, convert to space-separated string.
      elif isinstance(attval,list):
        if any(hasattr(v,'name') for v in attval):
          # Store the dependencies for future use.
          obj.deps.extend(attval)
          # Convert attribute to string.
          atts[attname] = ' '.join(v.name for v in attval)
        # Remove the attribute if the list of objects is empty
        elif len(attval) == 0:
          atts.pop(attname)
      # Detect dictionaries, convert to "key1: value1 key2: value2"
      elif isinstance(attval,OrderedDict):
        if len(attval) > 0:
          # Store the dependencies for future use.
          obj.deps.extend([k for k in attval.keys() if hasattr(k,'name')])
          obj.deps.extend([v for v in attval.values() if hasattr(v,'name')])
          # Convert attribute to string.
          attval = [getattr(k,'name',k)+': '+getattr(v,'name',v) for k,v in attval.items()]
          atts[attname] = ' '.join(attval)
        # Remove the attribute if the list of objects is empty
        else:
          atts.pop(attname)


  def _fix_names (self):

    # List of metadata keys that are internal to the FSTD file.
    internal_meta = self._headers.dtype.names

    # Generate unique axis names.
    axis_table = dict()
    for axis in self._iter_axes():
      if axis.name not in axis_table:
        axis_table[axis.name] = []
      axis_table[axis.name].append(axis)
    for axisname, axis_list in axis_table.items():
      if len(axis_list) == 1: continue
      warn (_("Multiple %s axes.  Appending integer suffixes to their names.")%axisname)
      for i,axis in enumerate(axis_list):
        axis.name = axis.name+str(i+1)

    # Generate a string-based variable id.
    # Only works for true variables from the FSTD source
    # (needs metadata like etiket, etc.)
    def get_var_id (var):
      out = []
      for fmt in self._human_var_id:
        out.append(fmt%var.atts)
      return tuple(out)

    # Generate unique variable names.
    var_table = dict()
    for var in self._iter_objects():
      if var.name not in var_table:
        var_table[var.name] = []
      # Identify the variables by their index in the master list.
      var_table[var.name].append(var)

    for varname, var_list in var_table.items():
      # Only need to rename variables that are non-unique.
      if len(var_list) == 1: continue
      try:
        var_ids = [get_var_id(v) for v in var_list]
      except KeyError:
        # Some derived axes may not have enough metadata to generate an id,
        # so the best we can do is append an integer suffix.
        var_ids = [(str(r),) for r in range(1,len(var_list)+1)]

      var_ids = zip(*var_ids)

      # Omit parts of the var_id that are invariant over all the variables.
      var_ids = [var_id for var_id in var_ids if len(set(var_id)) > 1]
      # Starting from the rightmost key, remove as many keys as possible while
      # maintaining uniqueness.
      for j in reversed(range(len(var_ids))):
        test = var_ids[:j] + var_ids[j+1:]
        if len(set(zip(*test))) == len(var_list):
          var_ids = test

      var_ids = zip(*var_ids)

      var_ids = ['_'.join(var_id) for var_id in var_ids]

      warn (_("Multiple definitions of %s.  Adding unique suffixes %s.")%(varname, ', '.join(var_ids)))

      # Apply the name changes.
      for var, var_id in zip(var_list, var_ids):
        var.name = var.name + '_' + var_id

    for var in self._iter_objects():
      # Names must start with a letter or underscore.
      if not var.name[0].isalpha():
        warn(_("Renaming '%s' to '_%s'.")%(var.name,var.name))
        var.name = '_'+var.name

      # Strip out FSTD-specific metadata?
      if not hasattr(var,'atts'): continue
      if self._rpnstd_metadata_list is not None:
        for n in internal_meta:
          if n not in self._rpnstd_metadata_list:
            var.atts.pop(n,None)


  def to_netcdf (self, filename, nc_format='NETCDF4', global_metadata=None, zlib=False, compression=4, progress=False):
    """
    Write the records to a netCDF file.
    Requires the netCDF4 package.
    """
    from fstd2nc.mixins import _var_type, _ProgressBar, _FakeBar
    from netCDF4 import Dataset
    import numpy as np
    f = Dataset(filename, "w", format=nc_format)

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

    # Generate the variable structures.
    for var in self._iter_objects():

      # Write the variable.
      # Easy case: already have the data.
      if hasattr(var,'array'):
        v = f.createVariable(var.name, datatype=var.array.dtype, dimensions=var.dims, zlib=zlib, complevel=compression)
        # Write the metadata.
        v.setncatts(var.atts)
        v[()] = var.array
        continue
      # Hard case: only have the record indices, need to loop over the records.
      # Get the shape of a single record for the variable.
      if not hasattr(var,'record_id'): continue
      record_shape = var.shape[var.record_id.ndim:]
      # Use this as the "chunk size" for the netCDF file, to improve I/O
      # performance.
      chunksizes = (1,)*var.record_id.ndim + record_shape
      v = f.createVariable(var.name, datatype=var.dtype, dimensions=var.dims, zlib=zlib, complevel=compression, chunksizes=chunksizes, fill_value=getattr(self,'_fill_value',None))
      # Turn off auto scaling of variables - want to encode the values as-is.
      # 'scale_factor' and 'add_offset' will only be applied when *reading* the
      # the file after it's created.
      v.set_auto_scale(False)
      # Write the metadata.
      v.setncatts(var.atts)
      # Write the data.
      indices = list(np.ndindex(var.record_id.shape))
      # Sort the indices by FSTD key, so we're reading the records in the same
      # order as they're found on disk.
      keys = map(int,var.record_id.flatten())
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
    for r,shape,v,ind in bar.iter(sorted(io)):
      try:
        data = self._fstluk(r,dtype=v.dtype)['d'].transpose().reshape(shape)
        v[ind] = data
      except (IndexError,ValueError):
        warn(_("Internal problem with the script - unable to get data for '%s'")%v.name)
        continue

    f.close()

  # Alias "to_netcdf" as "write_nc_file" for backwards compatibility.
  write_nc_file = to_netcdf
