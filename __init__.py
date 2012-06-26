
# Hook into pygeode.formats

from pygeode import formats
import rpn, bmf

setattr(formats,'rpn',rpn)
formats.__all__.append('rpn')
setattr(formats,'bmf',bmf)
formats.__all__.append('bmf')

del formats, rpn, bmf

