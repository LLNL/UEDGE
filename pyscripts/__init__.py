from .uedge import *
from . import *

try:
    from uedge.__version__ import __version__
    from uedge.__version__ import GitHash
except:
    try:
        from __version__ import __version__
        from __version__ import GitHash
    except:
        __version__ = 'unknown'
        GitHash=None
print('>>> UEDGE version:',__version__) 
print('>>> UEDGE Git SHA:,',GitHash)    
