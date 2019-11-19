from .uedge import *
try:
    from uedge.__version__ import __version__
except:
    try:
        from __version__ import __version__
    except:
        __version__ = 'unknown'

