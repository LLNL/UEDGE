#import uedge as ue
from uedge import *

#import uefacets
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os

#-import some utilities for using OS
execfile(os.path.join(os.environ['HOME'], 'utils/python/osfun.py'))


#-in .bashrc: "export PYLIB=/home/umansky1/PyUEDGE/uedge/pylib"
execfile(os.environ['PYLIB']+"/plotmesh.py")
execfile(os.environ['PYLIB']+"/plotcontour.py")
execfile(os.environ['PYLIB']+"/plotvar.py")
execfile(os.environ['PYLIB']+"/paws.py")

plt.ion()


#-read UEDGE settings
execfile("rd_d3dHsm_in.py")


#-do a quick preliminary run to set all internals
bbb.restart=0; bbb.ftol=1e10; bbb.dtreal = 1e-6; bbb.exmain()


#-show grid
plotmesh(iso=1)
wait = raw_input("PAUSING, PRESS ENTER TO CONTINUE...")


#-run to steady state
bbb.restart=1; bbb.ftol=1e-8; 
bbb.isbcwdt=1
bbb.dtreal = 1e-14; bbb.itermx=30; bbb.exmain()
bbb.t_stop=1e0
bbb.rundt()
bbb.dtreal=1e20; bbb.isbcwdt=0; bbb.exmain()


#-show some results
plotvar(bbb.te/ev)
