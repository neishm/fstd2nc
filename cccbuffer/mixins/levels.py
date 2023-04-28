# Mixin for decoding levels from CCC files.

from fstd2nc.mixins import BufferBase
from fstd2nc.stdout import _, info, warn, error
class Levels(BufferBase):
  # This method adds additional command-line arguments needed for the CCC
  # file interpretation.
  @classmethod
  def _cmdline_args (cls, parser):
    super (Levels,cls)._cmdline_args(parser)
    parser.add_argument ('--level-type', choices=('pres','eta','code','auto'), help=_('How to interpret the level values.'), default='auto')

  # This method is responsible for setting up the initial table of information
  # for the records in the file.
  # It adds any specific metadata columns needed for the mixin.
  # In this case, we are adding/modifying fields related to the CCC file.
  def __init__ (self, *args, **kwargs):
    """
    level_type : str, optional
        How to interpret the level values.
        ('pres','eta','code','auto')
    """
    self._level_type = kwargs.pop('level_type','auto')
    super(Levels,self).__init__(*args,**kwargs)

  # This routine handles metadata and coordinates for the data
  # (after it has been structured into multi-dimensional arrays).
  # Here we annotate the level axis.
  def _makevars (self):
    from fstd2nc.mixins import _axis_type, _var_type
    import numpy as np
    from collections import OrderedDict
    super(Levels,self)._makevars()

    # Remove degenerate levels (single level, encoded as "0" or "1"?)
    for var in self._varlist:
      if 'level' in var.dims:
        ind = var.dims.index('level')
        axis = var.axes[ind]
        if len(axis) == 1 and axis.array[0] in (0,1):
          var.axes = var.axes[:ind] + var.axes[ind+1:]
          var.record_id = var.record_id.squeeze(axis=ind)

    for axis in self._iter_axes('level'):
      # Determine the type of vertical coordinate to use.
      level_type = self._level_type
      # Some heuristics for automatically matching to a coordinate type.
      if level_type == 'auto':
        level_type = 'code'
        if np.all(axis.array%5==0):
          level_type = 'pres'
        elif np.all((axis.array<=1000)|((axis.array>=99101)&(axis.array<=99999))):
          if np.any(axis.array>100):
            level_type = 'eta'
        if np.all(axis.array == np.arange(len(axis.array))+1):
          level_type = 'code'

      # Decode the vertical levels.
      if level_type != 'code':
        array = axis.array.astype('float32')
        array = np.where(axis.array<0, (-axis.array)%1000 * 10.0**(-(-axis.array//1000)-2), array)
        array = np.where(axis.array>=99101, (axis.array-99000)/10.0, array)
        axis.array = array
      # Keep track of the type of level chosen (needed later).
      axis.atts['level_type'] = level_type

    # Put the levels in the correct order.
    reorder = {}  # Keep track of reordering for each axis.
    for axis in self._iter_axes('level'):
      reorder[id(axis)] = np.argsort(axis.array)
    for var in self._varlist:
      if 'level' in var.dims:
        ind = var.dims.index('level')
        axis = var.axes[ind]
        var.record_id = var.record_id[(slice(None),)*ind+(reorder[id(axis)],)]
    for axis in self._iter_axes('level'):
      axis.array = axis.array[reorder[id(axis)]]

    # Annotate the vertical axis.
    for axis in self._iter_axes('level'):
      level_type = axis.atts.pop('level_type')
      if level_type == 'eta':
        # Scale the values.
        axis.array = axis.array / 1000.
        # Hard code plid / pref values for now.
        pref = 101320.0
        ptop = 50.0
        eta_top = ptop / pref
        power = 1.5
        b = _var_type ('b', {}, [axis], ((axis.array - eta_top) / (1.0 - eta_top)) ** power)
        ap = _var_type('ap', {}, [axis], pref*(axis.array - b.array))
        axis.atts['standard_name'] = 'atmosphere_hybrid_sigma_pressure_coordinate'
        axis.atts['axis'] = 'Z'
        axis.atts['positive'] = 'down'
        axis.atts['coordinates'] = [ap, b]
        axis.atts['formula_terms'] = OrderedDict([('ap',ap),('b',b),('ps','PS')])
      elif level_type == 'pres':
        array = axis.array.astype('float32')
        axis.array = array
        axis.name = 'pres'
        axis.atts['standard_name'] = 'air_pressure'
        axis.atts['units'] = 'mb'
        axis.atts['axis'] = 'Z'
        axis.atts['positive'] = 'down'

    #TODO: add level_descr attribute to variables.
