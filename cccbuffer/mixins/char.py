# Mixin for decoding CHAR records.
from fstd2nc.mixins import BufferBase
from fstd2nc.stdout import _, info, warn, error
class Char(BufferBase):

  # This routine handles metadata and coordinates for the data
  # (after it has been structured into multi-dimensional arrays).
  # Here we modify the dimensions and data type for character data.
  def _makevars (self):
    from fstd2nc.mixins import _dim_type
    import numpy as np
    super(Char,self)._makevars()
    for var in self._varlist:
      if var.atts['kind'] == 'CHAR':
        if 'ilat' not in var.dims: continue
        if 'ilg' not in var.dims: continue
        rows = var.getaxis('ilat')
        rows = _dim_type ('row',len(rows))
        cols = var.getaxis('ilg')
        cols = _dim_type ('col',len(cols)*8)
        var.axes[var.dims.index('ilat')] = rows
        var.axes[var.dims.index('ilg')] = cols
        var.dtype = np.dtype('|S1')

  # Modify the decoder to return character arrays for CHAR records.
  @classmethod
  def _decode (cls, data):
    import numpy as np
    kind = data[4:68].view('|S8')[0]
    data = super(Char,cls)._decode(data)
    if kind == b'CHAR    ':
      data = data.view('|S1')
    return data

