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
from .UEDGESettings import *
from .UEDGESimulation import *    
from .UEDGEDoc import *
from .UEDGEToolBox import *       
from .UEDGEBas2Py import *
#__all__=["UEDGEToolBox","UEDGESettings"]

