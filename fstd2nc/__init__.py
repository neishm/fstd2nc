###############################################################################
# Copyright 2017-2023 - Climate Research Division
#                       Environment and Climate Change Canada
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

__version__ = "0.20231105.1"


# Check for bundled rpnpy package.
# Fall back to this one if no standard rpnpy package available.
try:
  # Importing the module will set up the appropriate search paths.
  import fstd2nc_deps
  # Don't need a reference to the module after the paths are set.
  del fstd2nc_deps
except ImportError:
  pass

# Turn off annoying messages from librmn.
# Usually unwanted when using this code from Python.
# Can be turned back on from the command-line with --msglvl option.
try:
  import rpnpy.librmn.all as rmn
  rmn.fstopt("MSGLVL","ERRORS")
  del rmn
except ImportError:
  pass

# Combine all the mixins to create a final interface for I/O.
from fstd2nc.mixins.fstd import FSTD
from fstd2nc.mixins.select import SelectVars
from fstd2nc.mixins.ascii import ASCII
from fstd2nc.mixins.masks import Masks
from fstd2nc.mixins.dates import Dates
from fstd2nc.mixins.accum import Accum
from fstd2nc.mixins.ensembles import Ensembles
from fstd2nc.mixins.series import Series
from fstd2nc.mixins.sfc_codes import Sfc_Codes
from fstd2nc.mixins.vardict import VarDict
from fstd2nc.mixins.vcoords import VCoords
from fstd2nc.mixins.xycoords import XYCoords
from fstd2nc.mixins.mesh import Mesh
from fstd2nc.mixins.misc import NoNK
from fstd2nc.mixins.filter import FilterRecords
from fstd2nc.mixins.removestuff import RemoveStuff
from fstd2nc.mixins.gridhacks import GridHacks, Interp, YinYang, Crop
from fstd2nc.mixins.pruneaxes import PruneAxes
from fstd2nc.mixins.netcdf import netCDF_Atts, netCDF_IO
from fstd2nc.mixins.compat import FSTD_Compat
from fstd2nc.mixins.extern import ExternInput, ExternOutput
from fstd2nc.mixins.diaghacks import DiagHacks

class Buffer (DiagHacks,ExternOutput,FSTD_Compat,netCDF_IO,netCDF_Atts,PruneAxes,Crop,YinYang,Interp,GridHacks,RemoveStuff,FilterRecords,NoNK,Mesh,XYCoords,VCoords,Sfc_Codes,VarDict,Series,Ensembles,Accum,Dates,Masks,ASCII,SelectVars,ExternInput,FSTD):
  """
  High-level interface for FSTD data, to treat it as multi-dimensional arrays.
  Contains logic for dealing with most of the common FSTD file conventions.
  """
  def __init__ (self, filename, *args, **kwargs):
    super(Buffer,self).__init__(filename, *args,**kwargs)

# Dynamically generate final init docstring from the mixins.
def _docstring ():
  from fstd2nc.mixins import BufferBase
  base_doc = BufferBase.__init__.__doc__
  docstring = [base_doc.rstrip().strip('\n')]
  for cls in Buffer.__bases__[::-1]:
    doc = cls.__init__.__doc__
    if doc is None or doc == base_doc: continue
    docstring.append(doc.rstrip().strip('\n'))
  return '\n'.join(docstring)
try:
  # Python 2
  Buffer.__init__.__func__.__doc__ = _docstring()
except AttributeError:
  # Python 3
  Buffer.__init__.__doc__ = _docstring()
