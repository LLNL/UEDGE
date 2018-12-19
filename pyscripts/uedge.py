warp_version = "$Id: uedge.py,v 7.1 2018/10/18 18:13:14 meyer8 Exp $"
# import all of the neccesary packages
import __main__
import sys

#test for numpy
try:
  import numpy
  __main__.__dict__['with_numpy'] = 1
except:
  try:
    import Numeric
  except:
    raise SystemExit, "Neither numeric nor numpy present"

try:
  with_numpy = __main__.__dict__['with_numpy']
except KeyError:
  with_numpy = ('--with-numpy' in sys.argv)
  __main__.__dict__['with_numpy'] = with_numpy

if with_numpy:
  from numpy import *
  ArrayType = ndarray
  def gettypecode(x):
    return x.dtype.char
  def oldnonzero(a):
    return a.nonzero()[0]
else:
  from Numeric import *
  import MA
  def gettypecode(x):
    return x.typecode()
  oldnonzero = nonzero
import os.path
import time

if not with_numpy:
   #-- Set this to a more reasonable value
   MA.set_print_limit(10000)

# --- Since gist is only loaded on PE0 in the parallel
# --- environment, the parallel module must also be loaded in to check
# --- whether this is in fact running in parallel.

# --- This creates a logical lparallel which is true when running in parallel
# --- and false when running serial.
# --- This also creates a number of routines which are needed for parallel
# --- data handling. These routines all have sensible default results when
# --- not running in parallel.
from ueparallel import *

try:
  if me == 0 and sys.platform != 'mac':
    from gist import *
  else:
    from gistdummy import *
except ImportError:
  me = 0

# Print a greeting
# mmiah: Actually, lets not.
# print "Welcome to PyUedge"
# Import the uedgeC shared object which contains all of UEDGE
# from Forthon import *
import Forthon
try:
  import PyPDB
  from PyPDB import PW,PR
  from PyPDB.pypdb import *
except:
  # print "Unable to import PyPDB or * from PyPDB.pypdb."
  # print "Will proceed to try to import pypdb in case of old installation."
  try:
    from pypdb import *
  except:
    # print "pypdb not found."
    pass
#    raise
#  mmiah: Why are we reraising the exception?  We don't
#         necessarily need PyPDB
from uedgeC import *
from uedgeutils import *

# --- The UEDGE modules must be imported in the order below because of
# --- linking dependencies.
from compy import *
from aphpy import *
from apipy import *
from bbbpy import *
from flxpy import *
from grdpy import *
from svrpy import *
from wdfpy import *

# print "Packages loaded: ", Forthon.package()
# --- Add stuff to the path
import sys

# --- Set default runid to first filename in the command line, stripping off
# --- the .py suffix.
if sys.argv[0]:
  if sys.argv[0][-3:] == '.py':
    h,t = os.path.split(sys.argv[0][:-3])
    runid = t
    del h,t
  else:
    h,t = os.path.split(sys.argv[0])
    runid = t
    del h,t

# --- Check if the compiler was ifort - if so, set the stacksize unlimited
# --- The fcompname is not yet be available yet if Forthon is not up to date
try:
  if fcompname == 'ifort':
    import resource
    resource.setrlimit(resource.RLIMIT_STACK,(-1,-1))
except:
  pass
##############################################################################
######  Don't put anything below this line!!! ################################
##############################################################################

