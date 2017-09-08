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

# Grab necessary dependencies for running this stuff.
# (e.g., may need to hook in our own copy of rpnpy).
import fstd2nc
# Hook in translation support
from fstd2nc import _

"""
Serve RPN standard files through a pydap server.
"""

# Helper method - construct a Dataset object from the file path.
def make_dataset (filepath, buffer_cache={}, dataset_cache={}, mtimes={}, known_infiles={}):
  from fstd2nc import Buffer, _var_type
  from pydap.model import DatasetType, GridType, BaseType
  from os.path import basename, getmtime
  from glob import glob
  import numpy as np
  from collections import OrderedDict

  infiles = filepath
  buffer_args = dict()

  if filepath.endswith('.fstall'):
    # Read the extra arguments and parse.
    from argparse import ArgumentParser, Namespace
    from os import chdir, path
    import shlex
    # Change the working directory to where this file is, so that relative
    # paths work properly.
    chdir(path.dirname(filepath))
    # Parse the arguments from the file.
    parser = ArgumentParser()
    parser.add_argument('infile', nargs='+')
    Buffer._cmdline_args(parser)
    buffer_args = shlex.split(open(filepath).readline())
    buffer_args = parser.parse_args(buffer_args)
    buffer_args = vars(buffer_args)
    infiles = buffer_args.pop('infile')
    # Apply wildcard expansion to filenames.
    infiles = [f for filepattern in infiles for f in sorted(glob(filepattern)) or [filepattern]]
    # Make sure the filenames are strings (not unicode).
    infiles = map(str,infiles)

  if isinstance(infiles,str):
    mtime = getmtime(infiles)
  else:
    mtime = max(map(getmtime,infiles))
    # Look at modification time of control file.
    mtime = max(mtime,getmtime(filepath))

  # Return a cached version of the dataset if nothing about the file(s) have
  # changed since last time.
  if filepath in dataset_cache and mtime <= mtimes[filepath] and known_infiles[filepath] == infiles:
    return dataset_cache[filepath]

  mtimes[filepath] = mtime
  known_infiles[filepath] = infiles

  # Use the quick scan feature, and a private table for the Buffer.
  buffer_args['quick_scan'] = True
  buffer_args['private_table'] = True

  # Construct an fstd2nc Buffer object with the decoded FST data.
  buf = Buffer(infiles, **buffer_args)
  # Save a reference to the Buffer so the file reference(s) remain valid.
  buffer_cache[filepath] = buf

  # Construct a pydap Dataset object.
  global_metadata = buf._metadata.get('global',{})
  global_metadata['Conventions'] = "CF-1.6"
  dataset = DatasetType(name=basename(filepath), attributes=dict(NC_GLOBAL=global_metadata))
  # Save this so it can be immediately returned next time it's requested.
  dataset_cache[filepath] = dataset

  # Split into vars / dims.
  buf = list(buf)
  variables = OrderedDict((var.name,var) for var in buf if var.name not in var.axes)
  dims = OrderedDict((var.name,var) for var in buf if var.name in var.axes)

  # Add "dummy" dimensions (or they're not interpreted properly by some
  # clients like Panoply).
  for var in variables.values():
    for axisname, axisvalues in var.axes.items():
      if axisname not in dims:
        dims[axisname] = _var_type(axisname,{},axisvalues,np.array(axisvalues))

  # Based loosely on pydap's builtin netcdf handler.
  for var in variables.values():
    # Add grids.
    dataset[var.name] = GridType(var.name, var.atts)
    # Add array.
    dataset[var.name][var.name] = BaseType(var.name, var.array, tuple(var.axes.keys()), var.atts)
    # Add maps.
    for dim in var.axes.keys():
      if dim in dims:
        atts = dims[dim].atts
        array = dims[dim].array
        dataset[var.name][dim] = BaseType(dim, array, None, atts)

  for dim in dims.values():
    dataset[dim.name] = BaseType(dim.name, dim.array, None, dim.atts)
    # Handle unlimited dimension.
    if dim.name == 'time':
      dataset.attributes['DODS_EXTRA'] = {
        'Unlimited_Dimension': dim,
      }

  return dataset


# Handler for the FST data.
from pydap.handlers.lib import BaseHandler
class FST_Handler(BaseHandler):
  extensions=r'^.*\.fstall$'
  def __init__ (self, filepath):
    self.filepath = filepath
  # Only create the dataset object if needed.
  def __getattr__ (self, name):
    if name != 'dataset': raise AttributeError
    BaseHandler.__init__(self)
    self.dataset = make_dataset(self.filepath)
    return self.dataset

# Hack this handler into the pydap interface.
# Can't register this in the normal way, because we can't provide specific
# extensions for the handler.
from pydap.handlers.lib import get_handler as pydap_get_handler
def get_handler(filepath, handlers=None):
  from rpnpy.librmn.fstd98 import isFST
  if isFST(str(filepath)):
    return FST_Handler(str(filepath))
  return pydap_get_handler(filepath, handlers)
from pydap.wsgi import app as pydap_app
# Override the handler in the wsgi app (for running as pydap).
pydap_app.get_handler = get_handler


