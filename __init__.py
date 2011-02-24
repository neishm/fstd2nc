
# Hook into pygeode.formats

from pygeode import formats
import rpn

setattr(formats,'rpn',rpn)
formats.__all__.append('rpn')

del formats, rpn
