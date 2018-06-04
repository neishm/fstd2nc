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

__version__ = "0.20180604.0"


# Check for bundled rpnpy package.
# Fall back to this one if no standard rpnpy package available.
try:
  # Importing the module will set up the appropriate search paths.
  import fstd2nc_deps
  # Don't need a reference to the module after the paths are set.
  del fstd2nc_deps
except ImportError:
  pass

# Combine all the mixins to create a final interface for I/O.
from fstd2nc.mixins.select import SelectVars
from fstd2nc.mixins.masks import Masks
from fstd2nc.mixins.dates import Dates
from fstd2nc.mixins.series import Series
from fstd2nc.mixins.sfc_codes import Sfc_Codes
from fstd2nc.mixins.vcoords import VCoords
from fstd2nc.mixins.xycoords import XYCoords
from fstd2nc.mixins.misc import NoNK
from fstd2nc.mixins.filter import FilterRecords
from fstd2nc.mixins.removestuff import RemoveStuff
from fstd2nc.mixins.pruneaxes import PruneAxes
from fstd2nc.mixins.netcdf import netCDF_Atts, netCDF_IO
from fstd2nc.mixins.array import XArray
from fstd2nc.mixins.iter import Iter

class Buffer (Iter,XArray,netCDF_IO,netCDF_Atts,PruneAxes,RemoveStuff,FilterRecords,NoNK,XYCoords,VCoords,Sfc_Codes,Series,Dates,Masks,SelectVars):
  """
  High-level interface for FSTD data, to treat it as multi-dimensional arrays.
  Contains logic for dealing with most of the common FSTD file conventions.
  """


