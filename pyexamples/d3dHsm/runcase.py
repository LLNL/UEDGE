#-import uedge
from uedge import *

#-import hdf5 routines
from uedge.hdf5 import *

#-import graphics, math, etc.
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os

#-import some utilities for using OS
###execfile(os.path.join(os.environ['HOME'], 'utils/python/osfun.py'))


#-in .bashrc: "export PYLIB=/home/umansky1/PyUEDGE/uedge/pylib"
exec(open(os.environ['PYLIB']+"/plotmesh.py").read())
exec(open(os.environ['PYLIB']+"/plotcontour.py").read())
exec(open(os.environ['PYLIB']+"/plotvar.py").read())
exec(open(os.environ['PYLIB']+"/paws.py").read())
exec(open(os.environ['PYLIB']+"/osfun.py").read())

plt.ion()


#-read UEDGE settings
exec(open("rd_d3dHsm_in.py").read())


#-do a quick preliminary run to set all internals
bbb.restart=0; bbb.ftol=1e10; bbb.dtreal = 1e-6; bbb.exmain()


#-show grid
plotmesh(iso=1)
wait = input("PAUSING, PRESS ENTER TO CONTINUE...")


#-run to steady state
bbb.restart=1; bbb.ftol=1e-8; 
bbb.isbcwdt=1
bbb.dtreal = 1e-14; bbb.itermx=30; bbb.exmain()
bbb.t_stop=1e0
bbb.rundt()
bbb.dtreal=1e20; bbb.isbcwdt=0; bbb.exmain()


#-show some results
plotvar(bbb.te/bbb.ev)

#-export the solution in hdf5 file
hdf5_save('mycase.h5')

#-can be imported with this command
#hdf5_restore('mycase.h5')


###-refine the grid, interpolate to new grid, and restart:
#com.nxleg[0,0]=20; bbb.newgeo=1; bbb.icntnunk=0
#bbb.dtreal = 1e-14; bbb.isbcwdt=1; bbb.itermx=30; bbb.exmain()

###-time advance by another second
#bbb.t_stop=2e0; bbb.rundt()

###-now to steady state (infinite time)
#bbb.dtreal=1e20; bbb.isbcwdt=0; bbb.exmain()

###-show some results
#plotvar(bbb.te/bbb.ev)
