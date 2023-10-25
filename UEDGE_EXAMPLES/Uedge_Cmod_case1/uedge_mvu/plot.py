import os
import datetime
import copy
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.backends.backend_pdf import PdfPages
import h5py
from uedge import bbb, com, api
from uedge import __version__ as uedgeVersion



def plotmesh(iso=True, zshift=0.0, xlim=None, ylim=None, yinv=False, title="UEDGE grid", subtitle=None, show=True):

    #fig,ax = plt.subplots(1)

    if (iso):
        plt.axes().set_aspect('equal', 'datalim')
    else:
        plt.axes().set_aspect('auto', 'datalim')

    #plt.plot([np.min(com.rm),np.max(com.rm)], [np.min(com.zm),np.max(com.zm)])
        
    for iy in np.arange(0,com.ny+2):
        for ix in np.arange(0,com.nx+2):
            plt.plot(com.rm[ix,iy,[1,2,4,3,1]],
                     com.zm[ix,iy,[1,2,4,3,1]]+zshift, 
                     color="black", linewidth=0.5)
            
    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
    #fig.suptitle('UEDGE grid')
    #plt.title('UEDGE grid')
    plt.suptitle(title)
    plt.title(subtitle, loc="left")
    plt.grid(False)

    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)

    if yinv:
        plt.gca().invert_yaxis()
        
    if show:
        plt.show()



def plotvar(var, zshift=0.0, iso=True, grid=False, label=None, vmin=None, vmax=None, yinv=False, title="UEDGE data", subtitle=None, show=True):
    
    patches = []

    for iy in np.arange(0,com.ny+2):
        for ix in np.arange(0,com.nx+2):
            rcol=com.rm[ix,iy,[1,2,4,3]]
            zcol=com.zm[ix,iy,[1,2,4,3]]+zshift
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


     # Set vmin and vmax disregarding guard cells
    if not vmax:
        vmax = np.max(var)
    if not vmin:
        vmin = np.min(var)


    
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    ###p = PatchCollection(patches, cmap=cmap, norm=norm)
    p = PatchCollection(patches, norm=norm)
    p.set_array(np.array(vals))




    fig,ax = plt.subplots(1)

    ax.add_collection(p)
    ax.autoscale_view()
    plt.colorbar(p, label=label)

    if iso:
        plt.axis('equal')  # regular aspect-ratio
    
    fig.suptitle(title)
    ax.set_title(subtitle, loc="left")
    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')

    if grid:
        plt.grid(True)

    if yinv:
        plt.gca().invert_yaxis()
        
    #if (iso):
    #    plt.axes().set_aspect('equal', 'datalim')
    #else:
    #    plt.axes().set_aspect('auto', 'datalim')

    if show:
        plt.show()




def plotrprof(var, ixcut=-1, title="UEDGE data", subtitle=None, lines=True, dots=False, xlim=None, ylim=None, xlog=False, ylog=False, show=True):
    # Plotting radial profiles of UEDGE data
    #
    # Usage example:
    # plotrprof(bbb.te/ev, title="Te [eV]")
    # plotrprof(ni[:,:,1], title="Nn [m-3]")
    #==================================#
    
    fig,ax = plt.subplots(1)

    if (ixcut<0):
        ix0=bbb.ixmp
    else:
        ix0=ixcut

        

    if (lines):
        plt.plot(com.rm[ix0,:,0]-com.rm[ix0,com.iysptrx,0],var[ix0,:])
    
    if (dots):
        plt.plot(com.rm[ix0,:,0]-com.rm[ix0,com.iysptrx,0],var[ix0,:],"o")
        
    if xlim:
        plt.xlim(xlim)

    if ylim:
        plt.ylim(ylim)

    if ylog:
        plt.yscale('log')

    if xlog:
        plt.xscale('log')
        
    plt.xlabel('R-Rsep [m]')
    fig.suptitle(title)
    ax.set_title(subtitle, loc="right")
    plt.grid(True)

    if show:
        plt.show()





        
def show_flow(vx, vy, scale=1.0, ispec=0, color="red", xlim=None, ylim=None, title=None, subtitle=None):
#-show a vector field given by its UEDGE x,y components

    #-RZ components of the vector field                                                                                
    vr=np.zeros([com.nx+2,com.ny+2])
    vz=np.zeros([com.nx+2,com.ny+2])


    for ix in range(com.nx+1):
        for jy in range(com.ny+1):

	    #-unit vector in poloidal direction (given by R,Z components)                                             
            ex=np.array([com.rm[ix,jy,2]-com.rm[ix,jy,1],com.zm[ix,jy,2]-com.zm[ix,jy,1]])
            ex=ex/np.sqrt(np.dot(ex,ex))

            #-unit vectors in the radial direction (given by R,Z components)                                          
            ey=np.array([-ex[1],ex[0]])
            #ey=ey/np.sqrt(np.dot(ey,ey)) #-redundant


            #-change the direction in the private regions
            sign=1
            
            ##-inner lower leg PF
            if (ix>=0 and ix<=com.ixpt1[0] and jy<=com.iysptrx1[0]):
                sign=-1
                
            ##-inner upper leg PF
            if (ix>=com.ixpt2[0]+1 and ix<=com.ixlb[1]-1 and jy<=com.iysptrx1[0]+1):
                sign=-1

            ##-outer lower leg PF
            if (ix>=com.ixpt2[1]+1 and ix<=com.nx+1 and jy<=com.iysptrx1[0]):
                sign=-1
            
            ##-outer upper leg PF
            if (ix>=com.ixlb[1] and ix<=com.ixpt1[1] and jy<=com.iysptrx1[0]+1):
                sign=-1

            ey=ey*sign
                
                            
            #-flux vector (in RZ components)                                                                          
            vr[ix,jy]=vx[ix,jy]*ex[0] + vy[ix,jy]*ey[0]
            vz[ix,jy]=vx[ix,jy]*ex[1] + vy[ix,jy]*ey[1]


            #-set zero in guard cells                                                                                 
            if (ix==0 or ix==com.nx+1 or jy==0 or jy==com.ny+1):
                vr[ix,jy]=0.0
                vz[ix,jy]=0.0

                
    plotmesh(show=False, xlim=xlim, ylim=ylim, title=title, subtitle=subtitle)
    plt.quiver(com.rm[:,:,0], com.zm[:,:,0], vr*scale, vz*scale, color=color)
    plt.show()



def show_flow_ni(ispec=0, flux=False, scale=1.0, color="red", xlim=None, ylim=None, title=None, subtitle=None):
    #-show ion or neutral density flow
    
    ##-X,Y components of the the flux vector
    vx=bbb.fnix[:,:,ispec]/com.dy
    vy=bbb.fniy[:,:,ispec]/com.dx

    if not(flux):
        #-velocity vector
        vx=vx/bbb.ni[:,:,ispec]
        vy=vy/bbb.ni[:,:,ispec]
    
    show_flow(vx, vy, scale=scale, color=color, xlim=xlim, ylim=ylim, title=title, subtitle=subtitle)



    

def show_flow_test(option=2, scale=1.0, color="red", xlim=None, ylim=None, title=None, subtitle=None):
    #-show test flow
    
    ##-X,Y components of the the flux vector
    
    if (option==0):
        print("Using unit poloidal vector flow")
        vx=1.0+np.zeros([com.nx+2,com.ny+2])
        vy=0.0+np.zeros([com.nx+2,com.ny+2])
    elif (option==1):
        print("Using unit radial vector flow")
        vx=0.0+np.zeros([com.nx+2,com.ny+2])
        vy=1.0+np.zeros([com.nx+2,com.ny+2])    
    else:
        print("Using unit radial=poloidal vector flow")
        vx=1.0+np.zeros([com.nx+2,com.ny+2])
        vy=1.0+np.zeros([com.nx+2,com.ny+2])

    show_flow(vx, vy, scale=scale, color=color, xlim=xlim, ylim=ylim, title=title, subtitle=subtitle)
