from .uedge import *
from os import path
import uedge.checkver
from uedge import __path__

with open(path.join(__path__[0],"VERSION")) as f:
    __version__ = f.read().replace('\n', '').strip()

# Load the startup file .uedgerc.py from cwd or home.
#
_homefile = path.join(path.expanduser('~'), '.uedgerc.py')
_localfile = path.join(path.expanduser('.'), '.uedgerc.py')

if path.exists(_localfile):
   with open(_localfile) as f:
      exec(open(_localfile).read())
elif path.exists(_homefile):
   with open(_homefile) as f:
      exec(open(_homefile).read())
      
#
# Ensure that the UEDGE internal version matches the Python release
#

uedge.bbb.uedge_ver = __version__
