from .uedge import *
try:
    from uedge.__version__ import __version__
    from uedge.__src__ import __src__
    import uedge.checkver
except:
    try:
        from __version__ import __version__
        from __src__ import __src__
        import checkver
    except:
        __version__ = 'unknown'
        __src__ = 'unknown'

