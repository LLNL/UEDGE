import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.axes as ax
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rcdefaults
from matplotlib import interactive
from uedge import *

# This file defines a function to plot the UEDGE mesh, and
# then calls the function to plot the entire mesh.
# To use this file in a simple way, give the following ands:
#    read plotmesh
#    nf
# The function could then be used in a more sophisticated way
# to plot portions of the mesh, possibly with customized plot limits
# (by resetting acom.ny of r_min, r_max, z_min, and z_max):
#    call plotmesh(ixmn,ixmx,iymn,iymx)
#    nf
# where ixmn, ixmx, iymn, and iymx are integer variables or
# expressions.  Always give an "nf" and after reading the file
# plotmesh or calling the function plotmesh.

# DEFINE THE PLOT FUNCTION --

def plotmesh(ixmin=None,ixmax=None,iymin=None,iymax=None,
          r_min=None,r_max=None,z_min=None,z_max=None,title=None,
          block=True):
   """
   plotmesh(ixmin=<int>,ixmax=<int>,iymin=<int>,iymax=<int>
            title=<string>,r_min=<val>,r_max=<val>,z_min=<val>,z_max=<val>,
            block=<True|False>)
      where ixmin, ixmax, iymin, and iymax are integer variables or
      expressions used to plot a portion of the grid. title is used as
      both the title and the figure name. Block default is True.
      
      The plot axis limits may be specified with r_rmin,r_max,z_min,z_max.
   """

   if ixmin == None: ixmin = com.nxomit
   if ixmax == None: ixmax = (com.nxm-1)
   if iymin == None: iymin = 0
   if iymax == None: iymax = (com.ny-1)
   if r_min == None: r_min = com.rm.min()
   if r_max == None: r_max = com.rm.max()
   if z_min == None: z_min = com.zm.min()
   if z_max == None: z_max = com.zm.max()

   rcdefaults()
   if title == None:
    title='Uedge Grid'
   plt.figure(title,figsize=(4.0,8.0))
   plt.title(title)

   try:
      plt.plot(com.xlim,com.ylim,'k-',label='Limiter',linewidth=3)
      plt.plot(com.xlim,com.ylim,'y-',label='Limiter',linewidth=1)
      plt.plot(com.rbdry,com.zbdry,'b-',label='Last Closed')
   except:
      pass
   plt.ylim(z_min,z_max)
   plt.xlim(r_min,r_max)
   for ix in range(ixmax-ixmin+1):
       for iy in range(iymax-iymin+1):
           r0 = [com.rm[ix,iy,1],com.rm[ix,iy,2],com.rm[ix,iy,4],com.rm[ix,iy,3],com.rm[ix,iy,1]]
           z0 = [com.zm[ix,iy,1],com.zm[ix,iy,2],com.zm[ix,iy,4],com.zm[ix,iy,3],com.zm[ix,iy,1]]
           plt.plot(r0,z0,'k-',label='Grid',linewidth=1)
   plt.ylabel('Z (m)')
   plt.xlabel('R (m)')
   plt.axes().set_aspect('equal')



   
   #plt.show(block=block)
   plt.ion()
   plt.show()
   plt.pause(0.001)


def plotmeshval(val,ixmin=None,ixmax=None,iymin=None,iymax=None,
          r_min=None,r_max=None,z_min=None,z_max=None,title=None,units=None,
          block=True):
   """
   plotmeshval(val,ixmin=<int>,ixmax=<int>,iymin=<int>,iymax=<int>
            title=<string>,units=<string>,block=<True|False>)
      Display 2-D quantity using polyfill.
      where ixmin, ixmax, iymin, and iymax are integer variables or
      expressions used to plot a portion of the grid. title is used as
      both the title and the figure name. Units are displayed in the
      side colorbar. Block default is True.

      The plot axis limits may be specified with r_rmin,r_max,z_min,z_max.
   """
   if ixmin == None: ixmin = com.nxomit
   if ixmax == None: ixmax = (com.nxm-1)
   if iymin == None: iymin = 0
   if iymax == None: iymax = (com.ny-1)
   if r_min == None: r_min = com.rm.min()
   if r_max == None: r_max = com.rm.max()
   if z_min == None: z_min = com.zm.min()
   if z_max == None: z_max = com.zm.max()
   rcdefaults()
   if title == None:
    title='Uedge'
   fig, ax = plt.subplots()

   verts = np.array([])
   z = np.array([])
   for ix in range(ixmax-ixmin+1):
       for iy in range(iymax-iymin+1):
           v = []
           v.append([com.rm[ix,iy,1],com.zm[ix,iy,1]])
           v.append([com.rm[ix,iy,2],com.zm[ix,iy,2]])
           v.append([com.rm[ix,iy,4],com.zm[ix,iy,4]])
           v.append([com.rm[ix,iy,3],com.zm[ix,iy,3]])
           verts = np.append(verts,v)
           z = np.append(z,val[ix,iy])
   verts = verts.reshape(len(z),4,2)
   ax.set_title(title)
   ax.set_ylabel('Z (m)')
   ax.set_xlabel('R (m)')
   ax.set_aspect('equal')

   coll = PolyCollection(verts,array=z,cmap=cm.jet,edgecolors='face')
   ax.add_collection(coll)
   ax.autoscale_view()
   cbar = fig.colorbar(coll,ax=ax)
   #if units != None: cbar.ax.set_ylabel(units,rotation=-90,va='bottom')
   if units != None: cbar.ax.set_ylabel(units,va='bottom')
   plt.ylim(z_min,z_max)
   plt.xlim(r_min,r_max)

   #plt.show(block=block)
   plt.ion()
   plt.show()
   plt.pause(0.001)


def profile(rval,zval,title=None,style=None,linewidth=None,xlabel=None,ylabel=None):
   """
   profile(rval,zval,title=<None>)
      title is used as both the title and the figure name. 
      Interactive is turned on so subsequent calls go to the same plot
      Style encoded color, line, and marker.  See matplotlib documention.
      examples: black solid line  - style='k-'
                red circle marks  - style='ro'
                green x marks and dotted line - style='gx--'
      
   """

   rcdefaults()
   interactive(True)
   if title == None:
      title='Uedge Profile'
   if style == None:
      style = 'k-'
   if linewidth == None:
      lw = 1
     
   plt.figure(title,figsize=(4.0,8.0))
   plt.title(title)

   try:
      plt.plot(rval,zval,style,linewidth=lw)
   except:
      pass
   if ylabel != None:
      plt.ylabel(ylabel)
   if xlabel != None:
      plt.xlabel(xlabel)
   #plt.axes().set_aspect('equal')
   plt.ion()
   plt.show()
   plt.pause(0.001)



