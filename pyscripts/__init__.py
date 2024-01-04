from .uedge import *
from os import path
from pathlib import Path

print('''             __  ____________  ____________
            / / / / ____/ __ \/ ____/ ____/
           / / / / __/ / / / / / __/ __/   
          / /_/ / /___/ /_/ / /_/ / /___   
          \____/_____/_____/\____/_____/
''')


print('''
======================= NOTE =======================
If you publish or present UEDGE results, you are
required to properly cite the use of UEDGE according
to the instructions at:
        http://software.llnl.gov/UEDGE/citing
''')

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


#
# Load the startup file .uedgerc.py from cwd or home.
#
_homepath = path.expanduser('~')
_homefile = Path('{}/.uedgerc.py'.format(_homepath))
_localpath = path.expanduser('.')
_localfile = Path('{}/.uedgerc.py'.format(_localpath))

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

#
# Display text instructing on use of citations
#

