from .uedge import *
from . import *

try:
    from uedge.__version__ import __version__
except:
    try:
        from __version__ import __version__
    except:
        __version__ = 'unknown'
print('# UEDGE version:',__version__)        

