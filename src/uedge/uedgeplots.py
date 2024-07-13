import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.axes as ax
import sys
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rcdefaults
from matplotlib import interactive
from uedge import *
from numpy import sin, cos
from scipy import spatial
from skimage.util import img_as_ubyte as bytescale


# This file defines a function to plot the UEDGE mesh, and # then calls the function to plot the entire mesh.
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


def plotmesh(ixmin=None, ixmax=None, iymin=None, iymax=None,
             r_min=None, r_max=None, z_min=None, z_max=None, title=None,
             block=False, figsize=(4.0, 8.0),xlabel=None,ylabel=None):
    """
    plotmesh(ixmin=<int>,ixmax=<int>,iymin=<int>,iymax=<int>
             title=<string>,r_min=<val>,r_max=<val>,z_min=<val>,z_max=<val>,
             block=<True|False>,xlabel=None,ylabel=None,zlabel=None)

       Plot the uedge grid.
       where ixmin, ixmax, iymin, and iymax are integer variables or
       expressions used to plot a portion of the grid. title is used as
       both the title and the figure name. Block default is True.

       The plot axis limits may be specified with r_rmin,r_max,z_min,z_max.
    """
    try:
        zrefl = com.zm
        zlim = com.ylim
        zreflbdry = com.zbdry
        if str(com.geometry) == str([b'uppersn         ']):
           zrefl = 2.0 * com.zmid - com.zm
           zlim = 2.0 * com.zmid - com.ylim
           zreflbdry = 2.0 * com.zmid - com.zbdry
    except:
        pass

    if ixmin == None:
        ixmin = com.nxomit
    if ixmax == None:
        ixmax = (com.nxm-1)
    if iymin == None:
        iymin = 0
    if iymax == None:
        iymax = (com.ny-1)
    if r_min == None:
        r_min = com.rm.min()
    if r_max == None:
        r_max = com.rm.max()
    if z_min == None:
        z_min = zrefl.min()
    if z_max == None:
        z_max = zrefl.max()

    rcdefaults()
    if title == None:
        title = 'Uedge Grid'
    fig,ax = plt.subplots(figsize=figsize)
    ax.set_title(title)

    try:
        ax.plot(com.xlim, zlim, 'k-', label='Limiter', linewidth=3)
        ax.plot(com.xlim, zlim, 'y-', label='Limiter', linewidth=1)
        ax.plot(com.rbdry, zreflbdry, 'b-', label='Last Closed')
    except:
        pass
    for ix in range(ixmax-ixmin+1):
        for iy in range(iymax-iymin+1):
            r0 = [com.rm[ix, iy, 1], com.rm[ix, iy, 2],
                  com.rm[ix, iy, 4], com.rm[ix, iy, 3], com.rm[ix, iy, 1]]
            z0 = [zrefl[ix, iy, 1], zrefl[ix, iy, 2],
                  zrefl[ix, iy, 4], zrefl[ix, iy, 3], zrefl[ix, iy, 1]]
            ax.plot(r0, z0, 'k-', label='Grid', linewidth=1)
    if ylabel == None: ax.set_ylabel('Z (m)')
    else: ax.set_ylabel(ylabel)
    if xlabel == None: ax.set_xlabel('R (m)')
    else: ax.set_xlabel(xlabel)
    ax.set_ylim(z_min, z_max)
    ax.set_xlim(r_min, r_max)
    ax.set_aspect('equal')

    plt.ion()
    plt.show(block=block)
    plt.pause(0.001)

def plotanymesh(verts, r_min=None, r_max=None, z_min=None, z_max=None, title=None,
             block=False, figsize=(4.0, 8.0),xlabel=None,ylabel=None):
    """
    plotanymesh(verts, title=<string>,r_min=<val>,r_max=<val>,z_min=<val>,z_max=<val>,
             block=<True|False>,xlabel=None,ylabel=None)

       Plot any polynomial NxM grid. verts dimensions are [0:N,0:M,0:nverts,0:2].
       Last dim is [:,:,:,0] is R array, [:,:,:,1] is Z array
       title is used as both the title and the figure name. Block default is True.

       The plot axis limits may be specified with r_rmin,r_max,z_min,z_max.
    """
    try:
        zrefl = com.zm
        zlim = com.ylim
        zreflbdry = com.zbdry
        if str(com.geometry) == str([b'uppersn         ']):
           zrefl = 2.0 * com.zmid - com.zm
           zlim = 2.0 * com.zmid - com.ylim
           zreflbdry = 2.0 * com.zmid - com.zbdry
    except:
        pass

    if r_min == None:
        r_min = np.min(verts[:,:,:,0])
    if r_max == None:
        r_max = np.max(verts[:,:,:,0])
    if z_min == None:
        z_min = np.min(verts[:,:,:,1])
    if z_max == None:
        z_max = np.max(verts[:,:,:,1])

    rcdefaults()
    if title == None:
        title = 'Grid'
    fig,ax = plt.subplots(figsize=figsize)
    ax.set_title(title)

    try:
        ax.plot(com.xlim, zlim, 'k-', label='Limiter', linewidth=3)
        ax.plot(com.xlim, zlim, 'y-', label='Limiter', linewidth=1)
        ax.plot(com.rbdry, zreflbdry, 'b-', label='Last Closed')
    except:
        pass
    s = verts.shape
    xlen = s[0]
    ylen = s[1]
    for ix in range(xlen):
        for iy in range(ylen):
            r0 = [verts[ix, iy, 0, 0], verts[ix, iy, 1, 0],
                  verts[ix, iy, 2, 0], verts[ix, iy, 3, 0], 
                  verts[ix, iy, 0, 0]]
            z0 = [verts[ix, iy, 0, 1], verts[ix, iy, 1, 1],
                  verts[ix, iy, 2, 1], verts[ix, iy, 3, 1], 
                  verts[ix, iy, 0, 1]]
            ax.plot(r0, z0, 'k-', label='Grid', linewidth=1)
    if ylabel == None: ax.set_ylabel('Z (m)')
    else: ax.set_ylabel(ylabel)
    if xlabel == None: ax.set_xlabel('R (m)')
    else: ax.set_xlabel(xlabel)
    ax.set_ylim(z_min, z_max)
    ax.set_xlim(r_min, r_max)
    ax.set_aspect('equal')

    plt.ion()
    plt.show(block=block)
    plt.pause(0.001)

def plotmeshval(val, ixmin=None, ixmax=None, iymin=None, iymax=None,
                r_min=None, r_max=None, z_min=None, z_max=None, title=None, units=None,
                block=False,xlabel=None,ylabel=None,zlabel=None,figsize=(5.0,8.0)):
    """
    plotmeshval(val,ixmin=<int>,ixmax=<int>,iymin=<int>,iymax=<int>
             title=<string>,units=<string>,block=<True|False>
             xlabel=None,ylabel=None,zlabel=None)

       Display Uedge 2-D quantity using polyfill.
       where ixmin, ixmax, iymin, and iymax are integer variables or
       expressions used to plot a portion of the grid. title is used as
       both the title and the figure name. Units are displayed in the
       side colorbar. Block default is True.

       The plot axis limits may be specified with r_rmin,r_max,z_min,z_max.
    """
    try:
        zrefl = com.zm
        zlim = com.ylim
        zreflbdry = com.zbdry
        if str(com.geometry) == str([b'uppersn         ']):
           zrefl = 2.0 * com.zmid - com.zm
           zlim = 2.0 * com.zmid - com.ylim
           zreflbdry = 2.0 * com.zmid - com.zbdry
    except:
        pass

    if ixmin == None:
        ixmin = com.nxomit
    if ixmax == None:
        ixmax = (com.nxm-1)
    if iymin == None:
        iymin = 0
    if iymax == None:
        iymax = (com.ny-1)
    if r_min == None:
        r_min = com.rm.min()
    if r_max == None:
        r_max = com.rm.max()
    if z_min == None:
        z_min = zrefl.min()
    if z_max == None:
        z_max = zrefl.max()
    rcdefaults()
    if title == None:
        title = 'Uedge'
    fig, ax = plt.subplots(figsize=figsize)

    verts = np.array([])
    z = np.array([])
    for ix in range(ixmax-ixmin+1):
        for iy in range(iymax-iymin+1):
            v = []
            v.append([com.rm[ix, iy, 1], zrefl[ix, iy, 1]])
            v.append([com.rm[ix, iy, 2], zrefl[ix, iy, 2]])
            v.append([com.rm[ix, iy, 4], zrefl[ix, iy, 4]])
            v.append([com.rm[ix, iy, 3], zrefl[ix, iy, 3]])
            verts = np.append(verts, v)
            z = np.append(z, val[ix, iy])
    verts = verts.reshape(len(z), 4, 2)
    ax.set_title(title)
    if ylabel == None: ax.set_ylabel('Z (m)')
    else: ax.set_ylabel(ylabel)
    if xlabel == None: ax.set_xlabel('R (m)')
    else: ax.set_xlabel(xlabel)
    try:
        ax.plot(com.xlim, zlim, 'k-', label='Limiter', linewidth=3)
        ax.plot(com.xlim, zlim, 'y-', label='Limiter', linewidth=1)
        ax.plot(com.rbdry, zreflbdry, 'b-', label='Last Closed')
    except:
        pass

    coll = PolyCollection(verts, array=z, cmap=cm.jet, edgecolors='face')
    ax.add_collection(coll)
    ax.autoscale_view()
    cbar = fig.colorbar(coll, ax=ax,label=zlabel)
    # if units != None: cbar.ax.set_ylabel(units,rotation=-90,va='bottom')
    if units != None:
        cbar.ax.set_ylabel(units, va='bottom')
    ax.set_ylim(z_min, z_max)
    ax.set_xlim(r_min, r_max)
    ax.set_aspect('equal')

    plt.ion()
    plt.show(block=block)
    plt.pause(0.001)

def plotanymeshval(verts,z, r_min=None, r_max=None, z_min=None, z_max=None, title=None, units=None,
                block=False,xlabel=None,ylabel=None,zlabel=None):
    """
    plotanymeshval(verts, val, title=<string>,units=<string>,block=<True|False>,
                   xlabel=None,ylabel=None,zlabel=None)

       Display 2-D (NxM) quantity, val, using polyfill of NxM polynomial grid verts 
       verts dimensions are [0:N,0:M,0:nverts,0:2].
       Last dim is [:,:,:,0] is R array, [:,:,:,1] is Z array
       title is used as both the title and the figure name. Units are displayed in the
       side colorbar. Block default is True.

       The plot axis limits may be specified with r_rmin,r_max,z_min,z_max.
    """
    try:
        zrefl = com.zm
        zlim = com.ylim
        zreflbdry = com.zbdry
        if str(com.geometry) == str([b'uppersn         ']):
           zrefl = 2.0 * com.zmid - com.zm
           zlim = 2.0 * com.zmid - com.ylim
           zreflbdry = 2.0 * com.zmid - com.zbdry
    except:
        pass

    if r_min == None:
        r_min = com.rm.min()
    if r_max == None:
        r_max = com.rm.max()
    if z_min == None:
        z_min = zrefl.min()
    if z_max == None:
        z_max = zrefl.max()
    rcdefaults()
    if title == None:
        title = 'Uedge'
    fig, ax = plt.subplots()

    ax.set_title(title)
    if ylabel == None: ax.set_ylabel('Z (m)')
    else: ax.set_ylabel(ylabel)
    if xlabel == None: ax.set_xlabel('R (m)')
    else: ax.set_xlabel(xlabel)
    try:
        ax.plot(com.xlim, zlim, 'k-', label='Limiter', linewidth=3)
        ax.plot(com.xlim, zlim, 'y-', label='Limiter', linewidth=1)
        ax.plot(com.rbdry, zreflbdry, 'b-', label='Last Closed')
    except:
        pass

    coll = PolyCollection(verts, array=z, cmap=cm.jet, edgecolors='face')
    ax.add_collection(coll)
    ax.autoscale_view()
    cbar = fig.colorbar(coll, ax=ax,label=zlabel)
    # if units != None: cbar.ax.set_ylabel(units,rotation=-90,va='bottom')
    if units != None:
        cbar.ax.set_ylabel(units, va='bottom')
    ax.set_ylim(z_min, z_max)
    ax.set_xlim(r_min, r_max)
    ax.set_aspect('equal')

    plt.ion()
    plt.show(block=block)
    plt.pause(0.001)

def mkdensityfile(filename, ival, renmin=None, renmax=None, samples=[500, 500, 500],
                  xrange=[-2.4, 2.4], yrange=[-2.4, 2.4], zrange=[0, 3.2], tree=None):
    """
    mkdensityfile(filename, ival,renmin=<float>,renmax=<float>, 
                  samples=[<xsamps>,<ysamps>,<zsamps>],
                  xrange=[xmin,xmax],yrange=[ymin,ymax],zrange=[zmin,zmax],
                  tree=<cKDTree object> )

       Output Povray include and density field file (df3) for rendering.
       where: renmin,renmax are the values scaled to 0,255 in the final 
                bytescaling
              samples is an array of three values giving the volume sampling 
                for the density file (def [500,500,500])
              xrange, yrange, zrange are the vessel dimensions of the sampled 
                volume (m) (def xrange[-2.4,2.4], yrange[-2.4,2.4], 
                zrange[0,3.2])
              tree is returned and may be reused for another call for efficiency

       The defaults are set for DIII-D and will sample the full torus at 
          about 1cm r resolution and .6cm in z.

    """
    try:
        zrefl = com.zm
        zlim = com.ylim
        zreflbdry = com.zbdry
        if str(com.geometry) == str([b'uppersn         ']):
           zrefl = 2.0 * com.zmid - com.zm
           zlim = 2.0 * com.zmid - com.ylim
           zreflbdry = 2.0 * com.zmid - com.zbdry
    except:
        pass

    if renmin == None:
        renmin = np.min(ival)
    if renmax == None:
        renmax = np.max(ival)
    nx, ny, nz = samples
    dims = np.array([nx, ny, nz], dtype=np.int16)

    file = open(filename, 'wb')
    if sys.byteorder == 'little':
       file.write(dims.byteswap(True))
    else:
       file.write(dims)
    rrm = com.rm[:, :, 0].ravel()
    rzm = zrefl[:, :, 0].ravel()
    if tree == None:
        tree = spatial.cKDTree(list(zip(rrm, rzm)))
    treelen = tree.data.shape[0]
    rpts = np.array([])
    zpts = np.array([])

    z, x, y = np.mgrid[
        zrange[0]:zrange[1]:complex(0, nz),
        xrange[0]:xrange[1]:complex(0, nx),
        yrange[0]:yrange[1]:complex(0, ny)
    ]
    r = (x*x + y*y)**0.5

    pts = list(zip(r.ravel(), z.ravel()))

    #d,i = tree.query(pts,k=1,distance_upper_bound=0.1)
    d, i = tree.query(pts, k=1)
    val = ival
    val[0, :] = 0
    val[-1, :] = 0
    val[:, 0] = 0
    val[:, -1] = 0
    vf = np.append(val.ravel(), [renmin])
    #dens = bytescale(np.average(vf[i],axis=1,weights=1./d),cmin=renmin,cmax=renmax)
    dens = bytescale((vf[i] - renmin)/(renmax - renmin))
    file.write(dens)
    file.close()

    return tree


def profile(rval, zval, title=None, style=None, linewidth=None, xlabel=None, ylabel=None, figsize=(4.0, 8.0), block=False,marker=None):
    """
    profile(xval,yval,title=<None>,style=<None>,linewidth=<None>,xlabel=<None>,ylabel=<None>,block=<True|False>,marker=<none>)
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
        title = 'Uedge Profile'
    if style == None:
        style = 'k-'
    if linewidth == None:
        lw = 1

    fig,ax = plt.subplots(figsize=figsize)
    ax.set_title(title)

    try:
        ax.plot(rval, zval, style, linewidth=lw,marker=marker)
    except:
        pass
    if ylabel != None:
        ax.set_ylabel(ylabel)
    if xlabel != None:
        ax.set_xlabel(xlabel)
    plt.ion()
    plt.show(block=block)
    plt.pause(0.001)
