try:
   import IPython
   from IPython.terminal.prompts import Prompts,Token
   from IPython.terminal.embed import InteractiveShellEmbed
except:
   pass

try:
   from traitlets.config.loader import Config
except:
   pass

import sys,os,__main__

import numpy as np
ArrayType = np.ndarray

if sys.hexversion >= 0x03000000:
    # --- With Python3, the so files of each Fortran package are imported
    # --- separately. The dlopen flag needs to be set so that cross references
    # --- among the packages can be satisfied.
    sys.setdlopenflags(os.RTLD_LAZY | os.RTLD_GLOBAL)

from . import uedgeC
from .uedgeC import *

from Forthon import *

if sys.hexversion >= 0x03000000:

    from .compy import *
    from .svrpy import *
    from .grdpy import *
    from .flxpy import *
    from .bbbpy import *
    from .wdfpy import *
    from .apipy import *
    from .aphpy import *
else:
    from wdfpy import *
    from svrpy import *
    from grdpy import *
    from flxpy import *
    from bbbpy import *
    from apipy import *
    from aphpy import *
    from compy import *

import time
import os.path
import __main__
warp_version = "$Id: uedge.py,v 7.4 2019/11/12 21:10:56 meyer8 Exp $"
# import all of the neccesary packages



def gettypecode(x):
    return x.dtype.char

def oldnonzero(a):
    return a.nonzero()[0]

# Print a greeting
# mmiah: Actually, lets not.
# print "Welcome to PyUedge"
# Import the uedgeC shared object which contains all of UEDGE
try:
    import PyPDB
    from PyPDB import PW, PR
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

# --- The UEDGE modules must be imported in the order below because of
# --- linking dependencies.

# --- Set default runid to first filename in the command line, stripping off
# --- the .py suffix.
if sys.argv[0]:
    if sys.argv[0][-3:] == '.py':
        h, t = os.path.split(sys.argv[0][:-3])
        runid = t
        del h, t
    else:
        h, t = os.path.split(sys.argv[0])
        runid = t
        del h, t

# --- Check if the compiler was ifort - if so, set the stacksize unlimited
# --- The fcompname is not yet be available yet if Forthon is not up to date
try:
    if fcompname == 'ifort':
        import resource
        resource.setrlimit(resource.RLIMIT_STACK, (-1, -1))
except:
    pass


try:

   class MyPrompt(Prompts):
     def in_prompt_tokens(self, cli=None):
         return [(Token.Prompt, 'UEDGE>>> ')]
     def out_prompt_tokens(self, cli=None):
         return [(Token.Prompt, 'UEDGE>>> ')]

   get_ipython
except:
   sys.ps1='UEDGE>>> '
else:
   ip = get_ipython()
   ip.prompts = MyPrompt(ip)


##############################################################################
######  Don't put anything below this line!!! ################################
##############################################################################
