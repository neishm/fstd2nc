# Mixin for decoding superlabels from CCC files.
from fstd2nc.mixins import BufferBase
from fstd2nc.stdout import _, info, warn, error
class Superlabels(BufferBase):

  # This method is responsible for setting up the initial table of information
  # for the records in the file.
  # It adds any specific metadata columns needed for the mixin.
  # In this case, we are adding/modifying fields related to the CCC file.
  def __init__ (self, *args, **kwargs):
    import numpy as np
    super(Superlabels,self).__init__(*args,**kwargs)
    self._outer_axes = self._outer_axes+('type',)
    # Assign superlabels to records.
    # Superlabels are "sticky", they apply to all records until new superlabel
    # is found.
    is_superlabel = (self._headers['kind'] == b'LABL    ')
    closest = np.zeros(self._nrecs, 'int32')
    closest[is_superlabel] = 1
    closest = np.cumsum(closest) - 1
    superlabels = []
    # Load the superlabels
    for superlabel_ind in np.where(is_superlabel)[0]:
      with open(self._files[self._headers['file_id'][superlabel_ind]]) as f:
        f.seek(self._headers['address'][superlabel_ind]+76,0)
        superlabels.append(np.fromfile(f,'|S80',1)[0].decode().strip())
    superlabels.append('')  # Handle index -1 which indicates no superlabel.
    superlabels = np.ma.array(superlabels,object)
    self._headers['superlabel'] = superlabels[closest]
    # Now ignore the superlabel entries in the table (they are not data).
    self._headers['selected'][is_superlabel] = False
    # For TIME records, treat the superlabels as an outer axis.
    self._headers['type'] = np.ma.empty(self._nrecs,superlabels.dtype)
    is_time = (self._headers['kind'] == b'TIME    ')
    self._headers['type'][is_time] = self._headers['superlabel'][is_time]
    self._headers['type'][~is_time] = np.ma.masked
    self._headers['superlabel'][is_time] = ''
    self._headers['superlabel'][is_time] = np.ma.masked

  # This routine handles metadata and coordinates for the data
  # (after it has been structured into multi-dimensional arrays).
  # Here we assign superlabels as an attribute for variables.
  def _makevars (self):
    import numpy as np
    from fstd2nc.mixins import _axis_type, _var_type
    from collections import OrderedDict
    super(Superlabels,self)._makevars()
    for var in self._varlist:
      # Add superlabel as an attribute for the variable.
      if 'superlabel' in var.atts and var.atts['superlabel'] != '':
        var.atts['label'] = var.atts['superlabel']
    # Ensure we have unique names for variables.
    # Sometimes have the same variable name multiple times, with
    # different superlabels.
    var_table = {}
    for var in self._varlist:
      if var.name not in var_table:
        var_table[var.name] = []
      var_table[var.name].append(var)
    for var_name, var_list in var_table.items():
      if len(var_list) == 1: continue
      # If the variable name is used multiple times, add an integer suffix.
      for i, var in enumerate(var_list):
        var.name = var.name + '_%d'%(i+1)
