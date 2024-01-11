
import json

pkg = 'uedge'

try:
   import importlib.metadata
   thisver = importlib.metadata.version(pkg)
except:
   import pkg_resources
   thisver = pkg_resources.get_distribution(pkg).version

try:
   import urllib.request
   contents = urllib.request.urlopen('https://pypi.org/pypi/'+pkg+'/json').read()
   data = json.loads(contents.decode())
   thatver = data['info']['version']
except:
   import urllib
   contents = urllib.urlopen('https://pypi.org/pypi/'+pkg+'/json').read()
   data = json.loads(contents)
   thatver = str(data['info']['version'])


if thisver > thatver:
   #print('Uedge version '+thisver+' is newer than available with pip ('+thatver+')')
   pass
elif thisver == thatver:
   #print('Uedge version '+thisver+' is up-to-date')
   pass
elif thisver < thatver:
   print()
   print('Uedge version '+thisver+', an update is available to '+thatver)
   print()
else:
   print()
   print('Some error checking pypi version')
   print()
