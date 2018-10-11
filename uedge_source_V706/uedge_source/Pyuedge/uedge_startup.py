import os
cpu = os.uname()[0] # os.environ['CPU']

import sys
sys.path.append('../../pyUtils/'+cpu)
sys.path.append('../../pyUtils/scripts')
sys.path.append(cpu)

from pyBasis import *

from uedge import bbb, com, api, aph, grd, flx
