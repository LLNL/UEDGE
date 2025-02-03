from .uedge import *
from os import path
from uedge import __path__
pkg = "uedge"
try:
    import json
    import urllib.request
except:
    pass
try:
    import importlib.metadata
    thisver = importlib.metadata.version(pkg)
except:
    import pkg_resources
    thisver = pkg_resources.get_distribution(pkg).version

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

#
# Set up crmpath to point to the included 
#
#aph.crmdir[0] = '/'.join([__path__[0],'rates/aph'])
aph.crmdir = __path__[0]
#
# Check version: included from checkver
# 
try:
    contents = urllib.request.urlopen("https://pypi.org/pypi/" + pkg + "/json").read()
    data = json.loads(contents.decode())
    thatver = data["info"]["version"]

    if thisver < thatver:
        print()
        print("Uedge version " + thisver + ", an update is available to " + thatver)
        print()

except Exception as err:
    print()
    print("Error checking pypi version: {}".format(err))
    print()
