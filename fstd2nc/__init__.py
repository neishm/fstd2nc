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
Functionality for converting between FSTD and netCDF files.
"""

__version__ = "0.20171103.0"

# Enable multi-language support.
from gettext import gettext as _
import gettext
from os import path, environ
gettext.bindtextdomain('fstd2nc', path.join(path.dirname(__file__),'locale'))
gettext.textdomain('fstd2nc')
# Check for special CMCLNG environment variable
if environ.get('CMCLNG') == 'francais':
  environ['LANGUAGE'] = 'fr_CA'
del gettext, path, environ

# Check for bundled rpnpy package.
# Fall back to this one if no standard rpnpy package available.
try:
  # Importing the module will set up the appropriate search paths.
  import fstd2nc_deps
  # Don't need a reference to the module after the paths are set.
  del fstd2nc_deps
except ImportError:
  pass

# Module-level variable to control the amount of information printed.
stdout_streams=('info','warn','error')


# Information messages
def info (msg):
  if 'info' in stdout_streams:
    print (msg)

# How to handle warning messages.
# E.g., can either pass them through warnings.warn, or simply print them.
def warn (msg, _printed=set()):
  if msg not in _printed and 'warn' in stdout_streams:
    print (_("Warning: %s")%msg)
    _printed.add(msg)

# Error messages
def error (msg):
  from sys import exit
  if 'error' in stdout_streams:
    print (_("Error: %s")%msg)
  exit(1)


# Combin all the mixins to create a final interface for I/O.
from fstd2nc.mixins.select import SelectVars
from fstd2nc.mixins.masks import Masks
from fstd2nc.mixins.dates import Dates
from fstd2nc.mixins.series import Series
from fstd2nc.mixins.vcoords import VCoords
from fstd2nc.mixins.xycoords import XYCoords
from fstd2nc.mixins.misc import NoNK
from fstd2nc.mixins.filter import FilterRecords
from fstd2nc.mixins.netcdf import netCDF_IO
from fstd2nc.mixins.array import Iter

class Buffer (Iter,netCDF_IO,FilterRecords,NoNK,XYCoords,VCoords,Series,Dates,Masks,SelectVars):
  """
  High-level interface for FSTD data, to treat it as multi-dimensional arrays.
  Contains logic for dealing with most of the common FSTD file conventions.
  """


