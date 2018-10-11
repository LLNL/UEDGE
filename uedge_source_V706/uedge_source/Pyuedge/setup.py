##########################################################################
#
# setup.py file for uedge
#   $Id: setup.py,v 7.0 2018/02/28 18:32:42 meyer8 Exp $
#   gmake runs $(PYTHON) setup.py install
#
##########################################################################

import os
import string
from distutils.core import setup,Extension

cwd              = os.getcwd()
pywrapper        = os.environ['PYUTILS'] + '/pywrapper'
fortlib          = os.environ['FORTLIB'] # -L... -l...
arch             = os.environ['ARCH']
packages         = os.environ['Packages'] 

pybasis_dir      = pywrapper + '/../' + arch
fort_libs        = string.split(fortlib) # otherwise gcc gets 1 string with
                                    # embedded blanks, and does not see -l's
fort_libs        = fort_libs + ['-lblas']
fortpath         = string.replace(fort_libs[0],'-L','')
packages         = string.split(packages)
package_objects  = []
package_files    = []

for pkg in packages:
   package_objects.append(pkg+'/'+arch+'/'+pkg+'.sor.o')
   package_objects.append(pkg+'/'+arch+'/'+pkg+'.p.o')
   package_files.append(pkg+'/'+arch+'/'+pkg+'pymodule.c')

aux_objects = ['Csor/'+arch+'/handlers.o', 'Fsor/'+arch+'/aux_tot.o',
                arch+'/uedgepymod.o' ]

setup(
      # ext_package = "uedge2py", # creates subdir under build & install dirs
      ext_modules = [
        Extension('uedge',  package_files,
          extra_objects = package_objects + aux_objects,
          include_dirs = [pywrapper],
          extra_link_args = [pybasis_dir+'/pybasisC.so'] + fort_libs +
                  [ '-Wl,-rpath='+fortpath]
                 )
                    ]
     )
