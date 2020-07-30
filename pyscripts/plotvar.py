#=================================================================#
#================= Plotting UEDGE 2D data on mesh  ===============#
#
#-Usage:
# execfile("../pylib/plotmesh.py") 
# plotvar(var)
# plotvar(var, iso=True)
#
# Input arguments:
#   var(0:com.ny+2, 0:com.ny+2)
#
#-Optional arguments:
#   iso (True/False) - True for equal aspect ratio
#
#-Expects imported modules:
# import matplotlib.pyplot as plt; import numpy as np
#
# Usage example:
# fig1 = plt.figure(); plotvar(com.rm[:,:,0])
#
# First coding: MVU, 22-jul-17
#=================================================================#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def plotvar(var, iso=False, title="UEDGE data"):
    
    patches = []

    for iy in np.arange(0,com.ny+2):
        for ix in np.arange(0,com.nx+2):
            rcol=com.rm[ix,iy,[1,2,4,3]]
            zcol=com.zm[ix,iy,[1,2,4,3]]
            rcol.shape=(4,1)
            zcol.shape=(4,1)
            polygon = Polygon(np.column_stack((rcol,zcol)), True)
            patches.append(polygon)

    #-is there a better way to cast input data into 2D array?
    vals=np.zeros((com.nx+2)*(com.ny+2))

    for iy in np.arange(0,com.ny+2):
        for ix in np.arange(0,com.nx+2):
            k=ix+(com.nx+2)*iy
            vals[k] = var[ix,iy]

    p = PatchCollection(patches)
    p.set_array(np.array(vals))




    fig,ax = plt.subplots(1)

    ax.add_collection(p)
    ax.autoscale_view()
    plt.colorbar(p)
    
    fig.suptitle(title)
    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
    plt.grid(True)

    #if (iso):
    #    plt.axes().set_aspect('equal', 'datalim')
    #else:
    #    plt.axes().set_aspect('auto', 'datalim')

    plt.show()

#=================================================================#
