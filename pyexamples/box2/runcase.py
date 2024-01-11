#import uedge
from uedge import *

#-import hdf5 routines
from uedge.hdf5 import *

#-import graphics, math, etc.
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
from uedge.uedgeplots import *
from uedge.rundt import rundt

#-import some utilities for using OS
###execfile(os.path.join(os.environ['HOME'], 'utils/python/osfun.py'))

def paws():
    programPause = input("Press the <ENTER> key to continue...")

##-how to do this better?
#-in .bashrc: "export PYLIB=/home/umansky1/PyUEDGE/uedge/pylib"
#execfile(os.environ['PYLIB']+"/plotmesh.py")
#execfile(os.environ['PYLIB']+"/plotcontour.py")
#execfile(os.environ['PYLIB']+"/plotvar.py")
#execfile(os.environ['PYLIB']+"/osfun.py")
#execfile("../../plotmesh.py")
#execfile("../../pylib/plotvar.py")
#execfile("../../pylib/plotr.py")
#execfile("../../pylib/showrange.py")
#execfile("../../pylib/paws.py")


plt.ion()


#-read UEDGE settings
from box2_in import *


#-do a quick preliminary run to set all internals
bbb.restart=0; bbb.ftol=1e10; bbb.dtreal = 1e-6; bbb.exmain()


#-show grid
plotmesh()
wait = input("PAUSING, PRESS ENTER TO CONTINUE...")


#-this should be done in uefacets
#ev=1.6022e-19


if (0):
    hdf5_restore('mycase.h5')
    bbb.dtreal = 1e20; bbb.exmain()
else:
    #-set up some initial state
    bbb.ngs=1e14; bbb.ng=1e14
    bbb.nis=1e20; bbb.ni=1e20 
    bbb.ups=0.0;  bbb.up=0.0
    bbb.tes=bbb.ev;   bbb.te=bbb.ev
    bbb.tis=bbb.ev;   bbb.ti=bbb.ev
    #
    #-Note: if you make a gap here then it will change the logic of if-else!
    #
    #-run to steady state
    bbb.restart=1; bbb.ftol=1e-8; 
    bbb.isbcwdt=1
    bbb.dtreal = 1e-14; bbb.itermx=30; bbb.exmain()
    bbb.t_stop=1e0
    rundt()
    bbb.dtreal=1e20; bbb.isbcwdt=0; bbb.exmain()
    hdf5_save('mycase.h5')


    ###execfile('plotcontour.py')
    from plotcontour import *
    paws()


    ##-now refine the solution on a larger grid
    #com.nycore[0]=2
    #com.nysol[0]=6
    #com.nxleg[0,1]=8
    #com.nxcore[0,1]=8
    #bbb.restart=1; bbb.newgeo=1; bbb.icntnunk=0
    #bbb.dtreal = 1e-14; bbb.ftol=1e-10;
    #bbb.isbcwdt=1; bbb.itermx=30; bbb.exmain()
    #plotmesh()
    #paws()
    #bbb.t_stop=2e0; bbb.ftol=1e-8; bbb.rundt()
    #bbb.dtreal=1e20; bbb.isbcwdt=0; bbb.exmain()
    #execfile('plotcontour.py')
    #paws()

#==========================================================================#
