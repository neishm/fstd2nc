
# Hook into pygeode.formats

from pygeode import formats
import rpn, bmf, geophy

setattr(formats,'rpn',rpn)
formats.__all__.append('rpn')
setattr(formats,'bmf',bmf)
formats.__all__.append('bmf')
setattr(formats,'geophy',geophy)
formats.__all__.append('geophy')

del formats, rpn, bmf, geophy

