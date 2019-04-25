#import uedge as ue
from uedge import *

#import uefacets
import matplotlib.pyplot as plt
import numpy as np
import os

##-how to do this better?
execfile("../../pylib/plotmesh.py")
execfile("../../pylib/plotvar.py")
execfile("../../pylib/plotr.py")
execfile("../../pylib/showrange.py")
execfile("../../pylib/paws.py")


plt.ion()


#-read UEDGE settings
execfile("box2_in.py")


#-do a quick preliminary run to set all internals
bbb.restart=0; bbb.ftol=1e10; bbb.dtreal = 1e-6; bbb.exmain()


#-show grid
plotmesh()
wait = raw_input("PAUSING, PRESS ENTER TO CONTINUE...")


#-this should be done in uefacets
ev=1.6022e-19


if (0):
    ue.restore("box2.h5")
    bbb.dtreal = 1e20; bbb.exmain()
else:
    #-set up some initial state
    bbb.ngs=1e14; bbb.ng=1e14
    bbb.nis=1e20; bbb.ni=1e20 
    bbb.ups=0.0;  bbb.up=0.0
    bbb.tes=ev;   bbb.te=ev
    bbb.tis=ev;   bbb.ti=ev
    #
    #-Note: if you make a gap here then it will change the logic of if-else!
    #
    #-run to steady state
    bbb.restart=1; bbb.ftol=1e-8; 
    bbb.isbcwdt=1
    bbb.dtreal = 1e-14; bbb.itermx=30; bbb.exmain()
    bbb.t_stop=1e0
    bbb.rundt()
    bbb.dtreal=1e20; bbb.isbcwdt=0; bbb.exmain()


    execfile('plotcontour.py')
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
