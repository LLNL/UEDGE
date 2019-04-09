# Plotting UEDGE computational mesh
#
#-Expects imported modules:
#import matplotlib.pyplot as plt; import numpy as np
#
#
#-Usage:
# execfile("../pylib/plotmesh.py") 
# plotmesh()
#
#-Optional arguments:
#   iso (True/False) - True for equal aspect ratio
#
#
# First coding: MVU, 17-jul-17
#=========================================================#

def plotmesh(iso=False):

    fig,ax = plt.subplots(1)

    if (iso):
        plt.axes().set_aspect('equal', 'datalim')
    else:
        plt.axes().set_aspect('auto', 'datalim')


    for iy in np.arange(0,com.ny+2):
        for ix in np.arange(0,com.nx+2):
            plt.plot(com.rm[ix,iy,[1,2,4,3,1]],
                     com.zm[ix,iy,[1,2,4,3,1]], 
                     color="b")


    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
    fig.suptitle('UEDGE mesh')
    plt.grid(True)

    plt.show()

#=========================================================#
