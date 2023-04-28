
# Combine all the mixins to create a final interface for I/O.
from fstd2nc.mixins.select import SelectVars
from cccbuffer.mixins.ccc import CCCMixin
from cccbuffer.mixins.char import Char
from cccbuffer.mixins.superlabels import Superlabels
from cccbuffer.mixins.times import Times
from cccbuffer.mixins.levels import Levels
from cccbuffer.mixins.grid import Grid
from fstd2nc.mixins.filter import FilterRecords
from fstd2nc.mixins.removestuff import RemoveStuff
from fstd2nc.mixins.pruneaxes import PruneAxes
from fstd2nc.mixins.netcdf import netCDF_Atts, netCDF_IO
from fstd2nc.mixins.extern import ExternInput, ExternOutput
class Buffer (ExternOutput,netCDF_IO,netCDF_Atts,PruneAxes,RemoveStuff,FilterRecords,Grid,Levels,Times,Superlabels,Char,CCCMixin,SelectVars): pass


