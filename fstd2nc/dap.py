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

"""
Serve RPN standard files through a pydap server.
"""

# Helper method - construct a Dataset object from a buffer argument string.
def dataset_from_str (name, buffer_str, mtime, directory='.', buffer_cache={}, dataset_cache={}, mtimes={}, known_infiles={}):
  from fstd2nc import Buffer
  from fstd2nc.mixins import _var_type, _axis_type, _dim_type
  from pydap.model import DatasetType, GridType, BaseType
  from os.path import basename, getmtime
  import numpy as np
  from collections import OrderedDict
  from datetime import datetime
  from argparse import ArgumentParser
  from os import chdir, path
  import shlex
  from glob import glob

  # Set the directory (to properly evaluate relative paths).
  chdir(directory)
  # Parse the arguments from the string.
  parser = ArgumentParser()
  parser.add_argument('infile', nargs='+')
  Buffer._cmdline_args(parser)
  buffer_args = shlex.split(buffer_str)
  buffer_args = parser.parse_args(buffer_args)
  buffer_args = vars(buffer_args)
  infiles = buffer_args.pop('infile')
  # Apply wildcard expansion to filenames.
  infiles = [f for filepattern in infiles for f in sorted(glob(filepattern)) or [filepattern]]
  # Make sure the filenames are strings (not unicode).
  infiles = list(map(str,infiles))

  # Look at modification times of individual files.
  mtime = max(map(getmtime,infiles))

  # Return a cached version of the dataset if nothing about the file(s) have
  # changed since last time.
  if name in dataset_cache and mtime <= mtimes[name] and known_infiles[name] == infiles:
    return dataset_cache[name]

  mtimes[name] = mtime
  known_infiles[name] = infiles

  # Construct an fstd2nc Buffer object with the decoded FST data.
  buf = Buffer(infiles, **buffer_args)
  # Save a reference to the Buffer so the file reference(s) remain valid.
  buffer_cache[name] = buf

  # Get global metadata.
  global_metadata = buf._metadata.get('global',{})
  # Add history to global metadata.
  timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  history = timestamp + ": %s via Pydap+fstd2dap"%path.basename(name)
  global_metadata = {"history":history}

  # Construct a pydap Dataset object.
  dataset = DatasetType(name=path.basename(name), attributes=dict(NC_GLOBAL=global_metadata))
  # Save this so it can be immediately returned next time it's requested.
  dataset_cache[name] = dataset

  # Split into vars / dims.
  buf = list(buf)
  variables = OrderedDict((var.name,var) for var in buf if not isinstance(var,(_axis_type,_dim_type)))
  dims = OrderedDict((var.name,var) for var in buf if isinstance(var,(_axis_type,_dim_type)))

  # Based loosely on pydap's builtin netcdf handler.
  for var in variables.values():
    # Add grids.
    dataset[var.name] = GridType(var.name, var.atts)
    # Add array.
    dataset[var.name][var.name] = BaseType(var.name, var.array, var.dims, var.atts)
    # Add maps.
    for dim in var.dims:
      if dim in dims:
        if hasattr(dims[dim],'array'):
          array = dims[dim].array
          atts = dims[dim].atts
        else:
          # Add "dummy" dimensions (or they're not interpreted properly by some
          # clients like Panoply).
          array = np.arange(len(dims[dim]))
        if hasattr(dims[dim],'atts'):
          atts = dims[dim].atts
        else:
          atts = {}
        dataset[var.name][dim] = BaseType(dim, array, None, atts)

  for dim in dims.values():
    if hasattr(dim,'array'):
      array = dim.array
    else:
      # Add "dummy" dimensions (or they're not interpreted properly by some
      # clients like Panoply).
      array = np.arange(len(dim))
    if hasattr(dim,'atts'):
      atts = dim.atts
    else:
      atts = {}
    dataset[dim.name] = BaseType(dim.name, array, None, atts)
    # Handle unlimited dimension.
    if dim.name == 'time':
      dataset.attributes['DODS_EXTRA'] = {
        'Unlimited_Dimension': dim.name,
      }

  return dataset

# Helper method - construct a Dataset object from the file path.
def make_dataset (filepath):
  from os.path import getmtime, dirname

  # Input filenames and arguments stored in a '.combo' file?
  if filepath.endswith('.combo'):
    buffer_str = open(filepath).readline()
  # Otherwise, assume the filepath is for a single FSTD file.
  # The only argument is the filename itself.
  else:
    buffer_str = filepath

  name = filepath
  mtime = getmtime(filepath)
  return dataset_from_str(name, buffer_str, mtime, directory=dirname(filepath))


# Handler for the FST data.
from pydap.handlers.lib import BaseHandler
class FST_Handler(BaseHandler):
  extensions=r'^.*\.combo$'
  def __init__ (self, filepath):
    self.filepath = filepath
    self.additional_headers = []
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
  from fstd2nc.extra import maybeFST as isFST
  if isFST(str(filepath)):
    return FST_Handler(str(filepath))
  return pydap_get_handler(filepath, handlers)
from pydap.wsgi import app as pydap_app
# Override the handler in the wsgi app (for running as pydap).
pydap_app.get_handler = get_handler

# Turn off warning / information messages (can clog up the logs).
import fstd2nc
fstd2nc.stdout.streams=('error',)
from rpnpy.librmn.fstd98 import fstopt
fstopt('MSGLVL','ERRORS')
del fstopt
import warnings
warnings.simplefilter("ignore")
del warnings

