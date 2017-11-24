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
from fstd2nc.mixins import Buffer_Base

#################################################
# Mixin for adding netCDF metadata to the variables

class netCDF_Atts (Buffer_Base):
  @classmethod
  def _cmdline_args (cls, parser):
    import argparse
    super(netCDF_Atts,cls)._cmdline_args(parser)
    parser.add_argument('--metadata-file', type=argparse.FileType('r'), action='append', help=_('Use metadata from the specified file.  You can repeat this option multiple times to build metadata from different sources.'))
  def __init__ (self, *args, **kwargs):
    import ConfigParser
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
    configparser = ConfigParser.SafeConfigParser()
    for metafile in metafiles:
      configparser.readfp(metafile)
    for varname in configparser.sections():
      metadata[varname] = OrderedDict(configparser.items(varname))
      # Detect numerical values
      for k,v in list(metadata[varname].items()):
        try:
          metadata[varname][k] = float(v) # First, try converting to float.
          metadata[varname][k] = int(v) # Try further conversion to int.
        except ValueError: pass
    self._metadata = metadata
    super(netCDF_Atts,self).__init__(*args,**kwargs)
  def _iter (self):
    from collections import OrderedDict

    axis_renames = {}

    for var in super(netCDF_Atts,self)._iter():
      orig_name = var.name
      # Add extra metadata provided by the user?
      if var.name in self._metadata:
        var.atts.update(self._metadata[var.name])
        # Rename the field? (Using special 'rename' key in the metadata file).
        if 'rename' in var.atts:
          var.name = var.atts.pop('rename')
          # Also rename any axis with this name.
          axis_renames[orig_name] = var.name

      # Check if any of the axes in this variable need to be renamed.
      axis_names, axis_values = zip(*var.axes.items())
      axis_names = [axis_renames.get(n,n) for n in axis_names]
      var.axes = OrderedDict(zip(axis_names,axis_values))

      yield var


#################################################
# Mixin for reading/writing FSTD data from/to netCDF files.

class netCDF_IO (Buffer_Base):
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

  def _iter (self):
    from fstd2nc.mixins import _var_type
    from datetime import datetime
    import numpy as np
    from netCDF4 import date2num

    if self._reference_date is None:
      reference_date = None
    else:
      reference_date = datetime.strptime(self._reference_date,'%Y-%m-%d')

    varlist = super(netCDF_IO,self)._iter()
    if self._unique_names:
      varlist = list(varlist)
      self._fix_names(varlist)

    for var in varlist:

      # Modify time axes to be relative units instead of datetime objects.
      if var.name in var.axes and isinstance(var,_var_type) and isinstance(var.array[0],np.datetime64):
         # Convert from np.datetime64 to datetime.datetime
        var.array = var.array.tolist()
        units = '%s since %s'%(self._time_units, reference_date or var.array[0])
        var.atts.update(units=units)
        var.array = np.asarray(date2num(var.array,units=units))

      yield var


  def _fix_names (self, varlist):
    from fstd2nc.mixins import _var_type
    from collections import OrderedDict

    # List of metadata keys that are internal to the FSTD file.
    internal_meta = self._headers.dtype.names

    # Generate unique axis names.
    axis_table = dict()
    for var in varlist:
      for axisname, axisvalues in var.axes.items():
        axisvalues = tuple(axisvalues)
        if axisname not in axis_table:
          axis_table[axisname] = []
        if axisvalues not in axis_table[axisname]:
          axis_table[axisname].append(axisvalues)
    axis_renames = dict()
    for axisname, axisvalues_list in axis_table.items():
      if len(axisvalues_list) == 1: continue
      warn (_("Multiple %s axes.  Appending integer suffixes to their names.")%axisname)
      for i,axisvalues in enumerate(axisvalues_list):
        axis_renames[(axisname,axisvalues)] = axisname+str(i+1)

    # Apply axis renames.
    def rename_axis ((axisname,axisvalues)):
      key = (axisname,tuple(axisvalues))
      if key in axis_renames:
        return (axis_renames[key],axisvalues)
      return (axisname,axisvalues)
    for var in varlist:
      # If this is a coordinate variable, use same renaming rules as the
      # dimension name.
      if isinstance(var,_var_type) and var.name in var.axes:
        var.name, axisvalues = rename_axis((var.name,var.axes[var.name]))
      var.axes = OrderedDict(map(rename_axis,var.axes.items()))

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
    for i, var in enumerate(varlist):
      if var.name not in var_table:
        var_table[var.name] = []
      # Identify the variables by their index in the master list.
      var_table[var.name].append(i)
    for varname, var_indices in var_table.items():
      # Only need to rename variables that are non-unique.
      if len(var_indices) == 1: continue
      try:
        var_ids = [get_var_id(varlist[i]) for i in var_indices]
      except KeyError:
        # Some derived axes may not have enough metadata to generate an id,
        # so the best we can do is append an integer suffix.
        var_ids = [(str(r),) for r in range(1,len(var_indices)+1)]

      var_ids = zip(*var_ids)

      # Omit parts of the var_id that are invariant over all the variables.
      var_ids = [var_id for var_id in var_ids if len(set(var_id)) > 1]
      # Starting from the rightmost key, remove as many keys as possible while
      # maintaining uniqueness.
      for j in reversed(range(len(var_ids))):
        test = var_ids[:j] + var_ids[j+1:]
        if len(set(zip(*test))) == len(var_indices):
          var_ids = test

      var_ids = zip(*var_ids)

      var_ids = ['_'.join(var_id) for var_id in var_ids]

      warn (_("Multiple definitions of %s.  Adding unique suffixes %s.")%(varname, ', '.join(var_ids)))

      # Apply the name changes.
      for i, var_id in zip(var_indices, var_ids):
        var = varlist[i]
        orig_varname = var.name
        var.name = var.name + '_' + var_id
        # Apply the name changes to any metadata that references this variable.
        for othervar in varlist:
          # Must match axes.
          if not set(var.axes.keys()) <= set(othervar.axes.keys()): continue
          for key,val in list(othervar.atts.items()):
            # Don't touch FSTD metadata.
            if key in internal_meta: continue
            # Can only modify string attributes.
            if not isinstance(val,str): continue
            val = val.split()
            if orig_varname in val:
              val[val.index(orig_varname)] = var.name
            othervar.atts[key] = ' '.join(val)

    for var in varlist:
      # Names must start with a letter or underscore.
      if not var.name[0].isalpha():
        warn(_("Renaming '%s' to '_%s'.")%(var.name,var.name))
        var.name = '_'+var.name

      # Strip out FSTD-specific metadata?
      if self._minimal_metadata:
        for n in internal_meta:
          var.atts.pop(n,None)


  def write_nc_file (self, filename, nc_format='NETCDF4', global_metadata=None, zlib=False):
    """
    Write the records to a netCDF file.
    Requires the netCDF4 package.
    """
    from fstd2nc.mixins import _var_type
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

    for var in self._iter():

      for axisname, axisvalues in var.axes.items():
        # Only need to create each dimension once (even if it's in multiple
        # variables).
        if axisname not in f.dimensions:
          # Special case: make the time dimension unlimited.
          if axisname == 'time':
            f.createDimension(axisname, None)
          else:
            f.createDimension(axisname, len(axisvalues))

      dimensions = list(var.axes.keys())

      # Write the variable.
      # Easy case: already have the data.
      if isinstance(var,_var_type):
        v = f.createVariable(var.name, datatype=var.array.dtype, dimensions=dimensions, zlib=zlib)
        # Write the metadata.
        v.setncatts(var.atts)
        v[()] = var.array
        continue
      # Hard case: only have the record indices, need to loop over the records.
      # Get the shape of a single record for the variable.
      record_shape = tuple(map(len,var.axes.values()))[var.record_id.ndim:]
      # Use this as the "chunk size" for the netCDF file, to improve I/O
      # performance.
      chunksizes = (1,)*var.record_id.ndim + record_shape
      v = f.createVariable(var.name, datatype=var.dtype, dimensions=dimensions, zlib=zlib, chunksizes=chunksizes)
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

    # Now, do the actual transcribing of the data.
    # Read/write the data in the same order of records in the RPN file(s) to
    # improve performance.
    bar = self._Bar(_("Saving netCDF file"), suffix="%(percent)d%% [%(myeta)s]")
    for r,shape,v,ind in bar.iter(sorted(io)):
      try:
        data = self._fstluk(r)['d'].transpose().reshape(shape)
        v[ind] = data.view(v.dtype)
      except (IndexError,ValueError):
        warn(_("Internal problem with the script - unable to get data for '%s'")%v.name)
        continue
    # We need to explicitly state that we're using CF conventions in our
    # output files, or some utilities (like IDV) won't accept the data.
    f.Conventions = "CF-1.6"

    f.close()


