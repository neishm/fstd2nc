# Mixin for detecting use of tile categories encoded in 'level' value.

from fstd2nc.mixins import BufferBase
from fstd2nc.stdout import _, info, warn, error
class Categories(BufferBase):

  # This routine handles metadata and coordinates for the data
  # (after it has been structured into multi-dimensional arrays).
  # Here we annotate the level axis.
  def _makevars (self):
    from fstd2nc.mixins import _iter_type
    import numpy as np
    super(Categories,self)._makevars()

    # Remove 'level' axis for cases where it is used in conjunction with
    # 'type' axis, and result is block-diagonal.
    # In that case, 'level' is more of a tile category (water, ice, etc.?)
    # and is redundant information if the specific types are included.
    for var in self._varlist:
      if isinstance(var,_iter_type) and 'level' in var.dims and 'type' in var.dims:
        ilevel = var.dims.index('level')
        if np.isin(-1,var.record_id):
          # Check if only one 'level' is valid for each 'type'
          # (the rest having record_id of -1 indicating missing data).
          if np.all(np.sum(var.record_id>=0,axis=ilevel) == 1):
            # Remove 'level' in that case, it's not actually an independent
            # dimension.
            var.record_id = var.record_id.max(axis=ilevel)
            var.axes = var.axes[:ilevel] + var.axes[ilevel+1:]
