# Holm10 Nov 5 2019, created from scratch
# 191203    - added showfig for displaying while in debug mode
#           - added projection options for all radial plots


from matplotlib.pyplot import ion
ion()



def display():
    ''' Shows the figure output while in debug mode '''
    from matplotlib.pyplot import show
    show(block=True)

def calcang(geo,nodeorder=range(5), ccw=True):
    """ Function calculating poloidal and radial angles at CC

    Parameters:
    geo       -   Dict geo must be in the global scope. Should contain the following data
                    rm - [nx,ny,5] array with radial position data
                    zm - as rm for the vertical data
    nodeorder   -   Indices in geo for [CC, SW, SE, NW, NE] nodes. rang(5) is default
    ccw         -   Switch whether to use counterclockwise (=True, default) or counterclockwise
                    (=False) directionality for poloidal direction 
                
    """
    from numpy import pi,zeros,shape,arctan
    
    # Unpack node order
    CC=nodeorder[0]
    SW=nodeorder[1]
    SE=nodeorder[2]
    NW=nodeorder[3]
    NE=nodeorder[4]


    # Angles to be stores
    angles=["angpol", "angrad", "angpolraw", "angradraw"]

    # Initialize arrays
    for dim in angles:
        geo[dim]=zeros(shape(geo["rm"][:,:,0]))

    # Calculate angles cell-wise
    for i in range(shape(geo["rm"][:,:,0])[0]):                  # Poloidal index
        for j in range(shape(geo["rm"][:,:,0])[1]):              # Radial index
            # Angle arrays: filled to contain N,E,S,W side center coordinates
            r,z=[],[]
            # Loop over each of the sidesi: N,E,S,W
            for side in [(NW,NE),(NE,SE),(SE,SW),(SW,NW)]: # Tuples for node coordinates
                r.append(0.5*(geo["rm"][i,j,side[0]]+geo["rm"][i,j,side[1]]))
                z.append(0.5*(geo["zm"][i,j,side[0]]+geo["zm"][i,j,side[1]]))

            # Side length arrays: filled to contain horizontal and vertical lengths
            dr,dz=[],[] 
            # Calculate side projections onto r & z coordinates over E-W and N-S nodes
            for side in [(1,3),(0,2)]:
                dr.append(r[side[0]]-r[side[1]])
                dz.append(z[side[0]]-z[side[1]])
            
            # Calculate angles in r-z coordinate system to x axis and store to datarr
            for k in range(2):
                # Pi-correction for when cell angles are against tilted "backwards"
                geo[angles[k]][i,j]=arctan(dz[k]/dr[k])+(dr[k]<0)*pi+(ccw==False)*(k==0 or k==2)*pi
                # TODO: What is radial?! Perp to poloidal direction, or radial in cell?!
                geo[angles[k+2]][i,j]=arctan(dz[k]/dr[k])+(ccw==False)*(k==0 or k==2)*pi






def ue_interpolate(val,interp=1):
    """ Interpolation routine from CC values to corner nodes of geo parameters

    Paramters:
    val     -       Array to be interpolated (including GC:s)
    interp  -       Choose model for interpolation to cell corner nodes.
                    =0 : Arithmetic mean
                    =1 (default) : Weighted L1 arithmetic mean
    """
    from numpy import sqrt,zeros
    from uedge import com,bbb

    ret=zeros((com.nx,com.ny,5))
    ret[:,:,0]=val[1:-1,1:-1]

    # Perform interpolation
    x1,x2,y1=com.ixpt1[0],com.ixpt2[0],com.iysptrx         # Helper indices
    xptind=[    [x1,y1,4],[x1,y1+1,2],[x1+1,y1,3],[x1+1,y1+1,1],    # LH X-pt
                [x2,y1,4],[x2,y1+1,2],[x2+1,y1,3],[x2+1,y1+1,1] ]   # RH X-pt
   
    for i in range(1,com.nx+1):   # Loop over poloidal cells
        for j in range(1,com.ny+1):   # Loop over radial cells
                                                    # Loop over [node, [neigh indices]]
            for n in [  [1,[    [i,j],  [bbb.ixm1[i,j],j],  [bbb.ixm1[i,j],j-1],    [i,j-1]     ]],
                        [2,[    [i,j],  [bbb.ixp1[i,j],j],  [bbb.ixp1[i,j],j-1],    [i,j-1]     ]],
                        [3,[    [i,j],  [bbb.ixm1[i,j],j],  [bbb.ixm1[i,j],j+1],    [i,j+1]     ]],
                        [4,[    [i,j],  [bbb.ixp1[i,j],j],  [bbb.ixp1[i,j],j+1],    [i,j+1]     ]]  ]:
                tot,w=0,0
                for neig in n[1]:
                    if interp==0: # Arithmetic interpolation
                        l=1
                    elif interp==1: # Weighted L1 arithmetic mean
                        l=sqrt( (com.rm[i,j,n[0]]-com.rm[neig[0],neig[1],0])**2+
                                (com.zm[i,j,n[0]]-com.rm[neig[0],neig[1],0])**2 ) #Distance between neigh CC and node
                    elif interp==2: # Harmonic mean
                        print("To be implemented")
                        exit(0)
                    tot+=val[neig[0],neig[1]]*l
                    w+=l
                ret[i-1,j-1]=tot/w
        
    return ret

def heatmap(Z,s=None,ax=False,zrange=False,cbar=True,cmap="magma",zoom="div",plotvessel=["sep","outline"],maintainaspect=True,
    xlabel=None,ylabel=None,title=None,units=None,zaxis="lin",showgrid=True):
    """Creates a heatmap of requested variable using polygons.
    heatmap(var,**keys)

    Variables:
    var:                2D or 3D array containing the cell-center values to be plotted (3D array controlled by s)

    Optional arguments:
    s[=0]:              Species index of 3D array to be plotted
    ax[=False]:         axis object on which to plot: if false, creates new figure and returns axis
    zrange:             Tuple defining the lower and upper Z-bounds. Defaults to min/max of Z
                        First tuple entry can be set to "min"/"max" to only limit either Z-boundary
    zaxis:              Z-axis type; "lin"/"log" (TODO log)
    cbar:               Boolean defining whether to plot a vertical colorbar on the same axis
    cmap[='magma']:     Colormap object to use with C, as defined by maplotlib.cm. 
    plotvessel:         variable defining whether to plot vessel outlines, grid boundaries, and separatrix
                        All are plotted is plotvessel=True. If only part of the boundaries should be plotted
                        plotvessel should be a list containing any or all of the following strings:
                            -"outline":    plots grid outline
                            -"sep":     plots separatrix
                            -"vessel"   plots vessel
    zoom[="div"]:   zoom area of the plot, one of the following:
                            -"ot":      outer target
                            -"it":      inner target
                            -"div":     divertor region
                            -"device":  whole device
    xlabel:         X-axis label string
    ylabel:         Y-axis label string
    title:          Plot title string
    units:          Colorbar units
    maintainaspect[=True]: Boolean whether to keep the aspect ration constant
    


    Returns: Figure object if ax=False [default], Void otherwise.
    """
    # TODO: implement levels-argument. If levels = N, plot N levels equidistant in Z. If levels=list, plot the requested levels
#   from uedge import *
    from uedge import com,bbb
    from matplotlib.patches import Polygon
    from matplotlib.colors import Normalize,LogNorm
    from matplotlib.pyplot import get_cmap,colorbar,figure
#    from plot_vector import pltves
    from numpy import shape,zeros,log10,floor,ceil
    from matplotlib.cm import ScalarMappable
   # from matplotlib.ticker import LogFormatter,LogLocator

    # Check if we have mutispecies array or not
    if len(Z.shape)==3:
        if s is None:   # Warn user if s is unset, and return sepcies index 0
            Z=Z[:,:,0]
            print('WARNING! Multi-species array requested. Plotting species index 0.')
        else:
            Z=Z[:,:,s]

    if ax is False:
        ret=figure()
        ax=ret.add_subplot(111)
    else:
        try:
            ax=ax.get_axes()[0]
        except:
            pass

    # Set zoom location
    if zoom=="ot":
        xlim = [com.rm[com.ixpt1[0]-2,0,0]*0.99,com.rm[-1,-1,0]*1.01]
        ylim = [com.zm.min()*0.99,com.zm[com.ixpt1[0],com.iysptrx+1,0]*1.01]
    elif zoom=="it":
        xlim = [com.rm[0,-1,0]*0.99,com.rm[com.ixpt2[0]+4,0,0]*1.01]
        ylim = [com.zm.min()*0.99,com.zm[0,-1,0]*1.01]
    elif zoom=="div":
        xlim = [com.rm[0,-1,0].min()*0.99,com.rm[-1,-1,0].max()*1.01]
        ylim = [com.zm.min()*0.99,com.zm[com.ixpt1[0],com.iysptrx+1,0]*1.01]
    elif zoom=="device":
        xlim,ylim= [com.rm.min(),com.rm.max()],[com.zm.min(),com.zm.max()]


    # Slab exception
    if bbb.mhdgeo==-1:
        plotvessel=False
        if zoom=="ot" or zoom=="div":
            ylim = [com.rm.min(),com.rm.max()]
            xlim = [com.zm[com.ixpt2[0],0,2],com.zm.max()]
        elif zoom=="device":
            ylim,xlim= [com.rm.min(),com.rm.max()],[com.zm.min(),com.zm.max()]
        else:
            print("Slab geometry is only compatible with 'ot' and 'device' zooming!")
            return


    # Set heatmap limits if requested
    if zrange is False:
        Zmax=Z[1:-1,1:-1].max()
        Zmin=Z[1:-1,1:-1].min()
    else:
        if isinstance(zrange[0],str): # Choosing to limit only one boundary
            if zrange[0]=="min":
                Zmin,Zmax=zrange[1],Z[1:-1,1:-1].max()
            elif zrange[0]=="max":
                Zmin,Zmax=Z[1:-1,1:-1].min(),zrange[1]
            else:
                print("zrange can only be defined as 'min'/'max'! Terminating...")
                exit(0)
        else:
            Zmin=zrange[0]
            Zmax=zrange[1]

    if zaxis=="lin":
        Zcol=(Z-Zmin)/(Zmax-Zmin)
    elif zaxis=="log":
        Zcol=((log10(Z)-floor(log10(Zmin)))/(floor(log10(Zmax))-floor(log10(Zmin))))
    else:
        print("Only valid zaxis options are 'lin' and 'log'!")
        return
        
    # Set colormap
    cmap=get_cmap(cmap)
    
    # Plot heatmap using polygons
    for i in range(shape(Z)[0]):    
        for j in range(shape(Z)[1]):    
            # Create polygon for each grid cell
            xy=zeros((4,2))
            if bbb.mhdgeo==-1:
                xy[:,1]=[com.rm[i,j,4],com.rm[i,j,2],com.rm[i,j,1],com.rm[i,j,3]]
                xy[:,0]=[com.zm[i,j,4],com.zm[i,j,2],com.zm[i,j,1],com.zm[i,j,3]]
            else:
                xy[:,0]=[com.rm[i,j,4],com.rm[i,j,2],com.rm[i,j,1],com.rm[i,j,3]]
                xy[:,1]=[com.zm[i,j,4],com.zm[i,j,2],com.zm[i,j,1],com.zm[i,j,3]]
            # Set color based on Z-value
#            col=cmap((Z[i,j]-Zmin)/Zmax)
            col=cmap(Zcol[i,j])
            # Plot the patch
            ax.add_patch(Polygon(xy,closed=True,facecolor=col,edgecolor=col))
    
    # Plot vessel if requested
    vesselparams={"plotoutline":True,"plotsep":True,"plotves":True}
    if plotvessel not in [True,False]: # Plot outlines if requested
        if "outline" not in plotvessel:
            vesselparams["plotoutline"]=False
        if "sep" not in plotvessel:
            vesselparams["plotsep"]=False
        if "vessel" not in plotvessel:
            vesselparams["plotves"]=False
    if plotvessel is not False:
        pltves(ax,**vesselparams)

    
    if bbb.mhdgeo==-1:    
        # Highlight sep and x-point
        ax.plot(com.zm[:-1,com.iysptrx+1,2],com.rm[:-1,com.iysptrx+1,2],"k-",linewidth=3) 
        ax.plot(com.zm[com.ixpt2[0]+1,:com.iysptrx+2,1][0],com.rm[com.ixpt2[0]+1,:com.iysptrx+2,1][0],"k-",linewidth=3) 
        ax.plot(com.zm[com.ixpt2[0],com.iysptrx+1:,2][0],com.rm[com.ixpt2[0],com.iysptrx+1:,2][0],"k--",linewidth=3) 

    # Set colorbar if requested
    if cbar is True:
        if zaxis=="lin":
            norm = Normalize(vmin=Zmin,vmax=Zmax)
        elif zaxis=="log":
            norm = Normalize(vmin=floor(log10(Zmin)),vmax=floor(log10(Zmax)))
            norm = LogNorm(vmin=Zmin,vmax=Zmax)
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
    
        # See if colorbar needs to be extended
        if zrange is False:
            extend="neither"
        elif zrange[0]>Z.min() and zrange[1]<Z.max():
            extend="both"
        elif zrange[0]>Z.min():
            extend="min"
        elif zrange[1]<Z.max():
            extend="max"
        else:
            extend="neither"
        cbar=colorbar(sm,ax=ax,extend=extend)
    
    # Plot grid if requested
    if showgrid is True:
        grid(ax=ax,color='grey',alpha=0.5,linewidth=0.1)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    # Check if aspect ratio should be preserved
    if maintainaspect is True:
        if bbb.mhdgeo!=-1:
            ax.set_aspect('equal', 'datalim')


    # Set additional labels as requested
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_xlabel(ylabel)
    if title is not None:
        ax.set_title(title)
    if units is not None:
        cbar.set_label(units)


    if "ret" in vars():
        ret.show()
        return ret 


    

def pltves(ax=False,plotoutline=True,plotves=True,plotsep=True):
    """ Creates the UEDGE grid and vessel on the supplied axis.

    Keyword arguments:
    ax[=False]:         axis object on which to plot: if false creates and returns figure
    plotoutline[=True]: bool defining whether to plot grid outlines or not
    plotsep[=True]:     bool defining whether to plot separatrix or not
    plotves[=True]:     bool defining whether to plot vessel outlines or not

    Returns: Figure object if ax=False [default], Void otherwise.
    """ 
   # from uedge import *    
    from uedge import com
    from matplotlib.pyplot import figure

    if not ax:
        ret=figure()
        ax=ret.add_subplot(111)

    if plotsep:
        sepx,sepy=[],[]
        for i in range(com.nx):  # Add each poloidal sep position
            sepx.append(com.rm[i,com.iysptrx,3])
            sepy.append(com.zm[i,com.iysptrx,3])
        # Complete RH Bound
        sepx.append(com.rm[-1,com.iysptrx,4])
        sepy.append(com.zm[-1,com.iysptrx,4])
                
        ax.plot(sepx,sepy,"k-",linewidth=0.5)
    if plotoutline:
        lplatex,lplatey,rplatex,rplatey,oboundx,oboundy,pfrboundx,pfrboundy,coreboundx,coreboundy=[],[],[],[],[],[],[],[],[],[]
        for i in range(com.nx):  # Add each poloidal sep position
            oboundx.append(com.rm[i,-1,3])
            oboundy.append(com.zm[i,-1,3])
        # Complete RH Bound
        oboundx.append(com.rm[-1,-1,4])
        oboundy.append(com.zm[-1,-1,4])
        ax.plot(oboundx,oboundy,"k-",linewidth=0.5)

        for i in range(com.ixpt1[0]+1):  # Add each poloidal sep position
            pfrboundx.append(com.rm[i,0,1])
            pfrboundy.append(com.zm[i,0,1])
        pfrboundx.append(com.rm[com.ixpt1[0],0,2])
        pfrboundy.append(com.zm[com.ixpt1[0],0,2])
        for i in range(com.ixpt2[0]+1,com.nx):  # Add each poloidal sep position
            pfrboundx.append(com.rm[i,0,1])
            pfrboundy.append(com.zm[i,0,1])
        # Complete RH Bound
        pfrboundx.append(com.rm[-1,0,2])
        pfrboundy.append(com.zm[-1,0,2])
        ax.plot(pfrboundx,pfrboundy,"k-",linewidth=0.5)
        
        for i in range(com.ixpt1[0]+1,com.ixpt2[0]+1):  # Add each poloidal sep position
            coreboundx.append(com.rm[i,0,1])
            coreboundy.append(com.zm[i,0,1])
        # Complete RH Bound
        coreboundx.append(com.rm[com.ixpt2[0],0,2])
        coreboundy.append(com.zm[com.ixpt2[0],0,2])
        ax.plot(coreboundx,coreboundy,"k-",linewidth=0.5)
        
        for j in range(com.ny+1):  # Add each poloidal sep position
            lplatex.append(com.rm[0,j,1])
            lplatey.append(com.zm[0,j,1])
        # Complete RH Bound
        lplatex.append(com.rm[0,-1,3])
        lplatey.append(com.zm[0,-1,3])
        ax.plot(lplatex,lplatey,"k-",linewidth=0.5)

        for j in range(com.ny+1):  # Add each poloidal sep position
            rplatex.append(com.rm[-1,j,1])
            rplatey.append(com.zm[-1,j,1])
        # Complete RH Bound
        rplatex.append(com.rm[-1,-1,3])
        rplatey.append(com.zm[-1,-1,3])
        ax.plot(rplatex,rplatey,"k-",linewidth=0.5)


    if plotves:
        ax.plot(com.xlim,com.ylim,"k-",linewidth=2)
    

    if "ret" in vars():
        ret.show()
        return ret 

def vector(poldata,raddata,ax=False,C=False,datascale=1, arrow_scale=10.,
        plotpol=True,plotrad=True,zoom="div",xlabel=None,ylabel=None,title=None,units=None,quiverunits=False,color=(0,0,0),cmap=False,plotvessel=["sep","outline"],
        maintainaspect=True,unitlength=False, norm=None,s=None):
    """Creates a quiver vector diagram on the supplied grid to the supplied axis.
    vecto(poldata,raddata,**keys)

    Variables:
    poldata:    2D array containing poloidal components
    raddata:    2D array containing radial components

    Optional arguments:
    s[=0]           Species index to be plotted for multi-species arrays
    ax[=False]:     axis object on which to plot: if False creates and returns fig 
    C[=0 -> "k"]:   2D array containing each arrow color in colormap ranges [0,1]
    cmap[="bwr"]:   Colormap object to use with C, as defined by maplotlib.cm.
    color[=(0,0,0)]:    (r,g,b) tubple defining arrow color. Overridden if C is defined
    datascale[=1]:      scaling value of data
    arrow_scale[=10]:   scaling value of arrow lengths 
    plotpol[=True]:     bool defining whether to plot poloidal contributions
    plotrad[=True]:     bool defining whether to plot radial contributions
    plotvessel:         variable defining whether to plot vessel outlines, grid boundaries, and separatrix
                            All are plotted is plotvessel=True. If only part of the boundaries should be plotted
                            plotvessel should be a list containing any or all of the following strings:
                            -"grid":    plots grid outline [default]
                            -"sep":     plots separatrix [default]
                            -"vessel"   plots vessel
    zoom[=div]:         zoom area of the plot, one of the following:
                            -"ot":      outer target
                            -"it":      inner target
                            -"div":     divertor region
                            -"device":  whole device
    maintainaspect[=True]:  Boolean whether to keep the aspect ration constant
    xlabel              xlabel string
    ylabel              ylabel string
    title               title string
    units               Units string
    quiverunits[=quiver default]: definition of quiver keyword argument unit. [default=quiver default]
    unitlength[=False]:     Boolean defining wheter to plot unit length arrows or not
    norm:                   Set the arrow scaling. Default scales to max in domain

    Returns: Figure object if ax=False [default], Void otherwise.
    """
    import matplotlib.cm as cm
    from matplotlib.pyplot import figure
    from matplotlib.colors import LinearSegmentedColormap
    from numpy import full,shape,cos,sin,floor,log10,ones
    from numpy.ma import array
    #plot_heatmap
#    from uedge import *
    from uedge import com,bbb

    if ax is False:
        ret=figure()
        ax=ret.add_subplot(111)
    else:
        try:
            ax=ax.get_axes()[0]
        except:
            pass


    # Check that shapes are correct
    if poldata.shape!=raddata.shape:
        return 'ERROR! Poloidal and radial arrays differ in shape. Aborting...'

    # Check if we have mutispecies array or not
    if len(poldata.shape)==3:
        if s is None:   # Warn user if s is unset, and return sepcies index 0
            poldata=poldata[:,:,0]
            raddata=raddata[:,:,0]
            print('WARNING! Multi-species array requested. Plotting species index 0.')
        else:
            poldata=poldata[:,:,s]
            raddata=raddata[:,:,s]
            

    # Calculate CC pol and rad angle from geom. data
    geo={"rm" : com.rm, "zm" : com.zm}
    calcang(geo)
    

    # Find which are to zoom to, mask out the remaining areas in order to get the arrow scaling right, and place legend
    mask=full(shape(geo["angrad"][1:-1,1:-1]),True)
    qx,qy=0.8,0.05
    if bbb.mhdgeo==-1:
        maintainaspect=False
        if zoom in ["it","div",'ot']:
            xlim,ylim= [com.rm.min()-0.02*com.rm.max(),com.rm.max()*1.02],[com.zm[com.ixpt2[0]+1,:,:].min()*0.98,com.zm.max()*1.02]
            mask[:com.ixpt2[0],:]=False
    elif zoom=="it":
        xlim = [com.rm[com.ixpt1[0]-2,0,0]*0.99,com.rm[-1,-1,0]*1.01]
        ylim = [com.zm.min()*0.99,com.zm[com.ixpt1[0],com.iysptrx+1,0]*1.01]
        mask[com.ixpt2[0]:,:]=False
    elif zoom=="ot":
        xlim = [com.rm[0,-1,0]*0.99,com.rm[com.ixpt2[0]+4,0,0]*1.01]
        ylim = [com.zm.min()*0.99,com.zm[0,-1,0]*1.01]
        mask[:com.ixpt1[0],:]=False
    elif zoom=="div":
        xlim = [com.rm[0,-1,0].min()*0.99,com.rm[-1,-1,0].max()*1.01]
        ylim = [com.zm.min()*0.99,com.zm[com.ixpt1[0],com.iysptrx+1,0]*1.01]
        mask[:com.ixpt1[0],:]=False
        mask[com.ixpt2[0]:,:]=False
    elif zoom=="device":
        xlim,ylim= [com.rm.min(),com.rm.max()],[com.zm.min(),com.zm.max()]
        mask=False
        qx,qy=0.5,0.5
    if bbb.mhdgeo==-1:
        xlim,ylim=ylim,xlim

    arroptions={"mask":mask} # Set mask as parameter
    #arroptions={} # Set mask as parameter

    # Calc carthesian contirbutions of poloidal and radial components
    # Also check whether to plot all contributions, and if data scaling is requested
    vix=datascale*(plotpol*poldata[1:-1,1:-1]*cos(geo["angpol"][1:-1,1:-1])+plotrad*raddata[1:-1,1:-1]*cos(geo["angrad"][1:-1,1:-1]))
    viy=datascale*(plotpol*poldata[1:-1,1:-1]*sin(geo["angpol"][1:-1,1:-1])+plotrad*raddata[1:-1,1:-1]*sin(geo["angrad"][1:-1,1:-1]))
            
    Xp, Yp = array(vix,**arroptions), array(viy,**arroptions) #Mask arrays
    
    vesselparams={"plotoutline":True,"plotsep":True,"plotves":True}
    if bbb.mhdgeo==-1:
        plotvessel=False
    if plotvessel not in [ True , False]: # Plot outlines if requested
        if "outline" not in plotvessel:
            vesselparams["plotioutline"]=False
        if "sep" not in plotvessel:
            vesselparams["plotsep"]=False
        if "vessel" not in plotvessel:
            vesselparams["plotves"]=False
    if plotvessel is not False:
        pltves(ax,**vesselparams)

    X, Y = com.rm[1:-1,1:-1,0],com.zm[1:-1,1:-1,0] # Parse geometric coords to X and Y

    # Calculate norm to determine arrow scaling 
    normarr=(Xp**2+Yp**2)**0.5

    if unitlength is False:
        normarr=round(normarr.max()/10**floor(log10(normarr.max())))*10**floor(log10(normarr.max())) # Calc scaling for arrows
    else:
        Xp,Yp=Xp/normarr,Yp/normarr
        norm=5

    if norm is None:
        norm=normarr


    # Set colormaps to default, requested single color, or colormap and color values
    if C is False:
        cmap=LinearSegmentedColormap("test",{   "red":((0,color[0],color[0]),(1,color[0],color[0])), 
                            "green": ((0,color[1],color[1]),(1,color[1],color[1])), 
                            "blue":((0,color[2],color[2]),(1,color[2],color[2]))        })
        C=ones(shape(X))
    else:
        if cmap is False:
            cmap="bwr"


    # Create quiver parameter dictionary
    quiverparams={"scale":arrow_scale*norm,"cmap":cmap} 

    if quiverunits is not False: # Set quiver units if requested
        quiverparams["units"]=quiverunits
    
    quiverparams["zorder"]=10

    if bbb.mhdgeo==-1:
        X,Y,Xp,Yp=Y,X,Yp,Xp

    Q=ax.quiver(X,Y,Xp,Yp,C,**quiverparams) # Plot quiver
    if bbb.mhdgeo==-1:    
        # Highlight sep and x-point
        ax.plot(com.zm[:-1,com.iysptrx+1,2],com.rm[:-1,com.iysptrx+1,2],"k-",linewidth=3) 
        ax.plot(com.zm[com.ixpt2[0]+1,:com.iysptrx+2,1][0],com.rm[com.ixpt2[0]+1,:com.iysptrx+2,1][0],"k-",linewidth=3) 
        ax.plot(com.zm[com.ixpt2[0],com.iysptrx+1:,2][0],com.rm[com.ixpt2[0],com.iysptrx+1:,2][0],"k--",linewidth=3) 
        # Plot vessel outline
        ax.plot(com.zm[:,0,2],com.rm[:,0,2],'k-')
        ax.plot(com.zm[:,-1,1],com.rm[:,-1,1],'k-')
        ax.plot(com.zm[0,:,1],com.rm[0,:,1],'k-')
        ax.plot(com.zm[-1,:,4],com.rm[-1,:,4],'k-')
    
    legtext="%.1E" % (norm) # Display arrow showing magnitudes


    # Set additional labels as requested
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_xlabel(ylabel)
    if title is not None:
        ax.set_title(title)
    if units is not None:
        legtext=legtext+' '+units

    if unitlength is not True:
        qk=ax.quiverkey(Q, qx, qy, norm, legtext, coordinates='axes') # Show arrow legend
    
    # Set axes
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if maintainaspect is True:
        ax.set_aspect('equal', 'datalim')

    if "ret" in vars():
        ret.show()
        return ret 





def grid(ax=False,showcc=False,showcut=False,showves=False,showplates=False,color='k',alpha=1,linewidth=0.2):
    """ Plots grid from UE execution 
    grid(**keys)

    Optional arguments:
    ax[=False]:         axis object on which to plot: if False creates and returns fig 
    showcc[=False]:     Show red dots at cell center coordinates
    showcut[=False]:    Mark the eight points around X-point (Left/Right Cuts in core, PFR and SOL)
    showves[=False]:    Plots the vessel outlines
    showplates[=False]: Shows the user-defined plates

    Returns: Figure object if ax=False [default], Void otherwise.
    """
    from matplotlib.pyplot import figure
    from uedge import com,grd,bbb


    # See if we are to add to axis or plot new fig
    if ax is False:
        ret=figure()
        ax=ret.add_subplot(111)
    else:
        try:
            ax=ax.get_axes()[0]
        except:
            pass

    for i in range(com.nx+2):   # Do poloidal dir
        for j in range(com.ny+2):   # Do radial dir
            pltx,plty=[],[]
            for k in [1,2,4,3,1]:
                if bbb.mhdgeo==-1:
                    plty.append(com.rm[i,j,k])
                    pltx.append(com.zm[i,j,k])
                else:
                    pltx.append(com.rm[i,j,k])
                    plty.append(com.zm[i,j,k])
                    
            ax.plot(pltx,plty,"-",linewidth=linewidth,color=color,alpha=alpha)
            if showcc:
                ax.plot(com.rm[i,j,0],com.zm[i,j,0],"r.")
    if showcut:    # Mark most important cells around X-point cut
        for xcell in [com.ixpt1[0],com.ixpt2[0]]: # Do each side of cut
            for dx in [0,1]: # Do both sides of cut
                for dy in [0,1]: # Do both sides of sep
                    if bbb.mhdgeo==-1:
                        ax.plot(com.zm[xcell+dx,com.iysptrx+dy,0],com.rm[xcell+dx,com.iysptrx+dy,0],"k.")
                    else:
                        ax.plot(com.rm[xcell+dx,com.iysptrx+dy,0],com.zm[xcell+dx,com.iysptrx+dy,0],"k.")
    if bbb.mhdgeo==-1:    
        # Highlight sep and x-point
        ax.plot(com.zm[:-1,com.iysptrx+1,2],com.rm[:-1,com.iysptrx+1,2],"-",linewidth=3,color=color,alpha=alpha) 
        ax.plot(com.zm[com.ixpt2[0]+1,:com.iysptrx+2,1][0],com.rm[com.ixpt2[0]+1,:com.iysptrx+2,1][0],"-",linewidth=3,color=color,alpha=alpha) 
        ax.plot(com.zm[com.ixpt2[0],com.iysptrx+1:,2][0],com.rm[com.ixpt2[0],com.iysptrx+1:,2][0],"--",linewidth=3,color=color,alpha=alpha) 
    

    if showves:
        ax.plot(com.xlim,com.ylim,"k-")
    if showplates:
        ax.plot(grd.rplate1,grd.zplate1,"ro-")
        ax.plot(grd.rplate2,grd.zplate2,"ro-")

            
    if bbb.mhdgeo==-1:
        ax.set_xlabel("Poloidal/Axial [m]")
        ax.set_ylabel("Radial  [m]")
    else:
        ax.set_aspect("equal","datalim")
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
       

    
    if "ret" in vars():
        ret.tight_layout() 
        ret.show()
        return ret 



def contours(Z,ax=False,cbar=True,cmap="YlOrRd",zoom="div",plotvessel=['sep','outline'],maintainaspect=True,
    labels=False,zaxis="lin"):
    """Creates a "quick and dirty" heatmap of requested variable. For serious plotting, use heatmap which is more rigorous but 2x slower
    cotours(var,**keys)

    Variables:
    var:      2D array containing the cell-center values to be plotted

    Optional arguments:
    ax[=False]:         axis object on which to plot: if false, creates new figure and returns axis
    cbar:               Boolean defining whether to plot a vertical colorbar on the same axis
    cmap[=YlOrRd]:      Colormap object to use with C, as defined by maplotlib.cm. 
    plotvessel[=True]:  variable defining whether to plot vessel outlines, grid boundaries, and separatrix
                        All are plotted is plotvessel=True. If only part of the boundaries should be plotted
                        plotvessel should be a list containing any or all of the following strings:
                            -"outline":    plots grid outline
                            -"sep":     plots separatrix
                            -"vessel"   plots vessel
    zoom[="div"]:   zoom area of the plot, one of the following:
                            -"ot":      outer target
                            -"it":      inner target
                            -"div":     divertor region
                            -"device":  whole device
    labels[=None]:     dictionary object for plot labels [default=None]. May contain any of the following keys:
                            -x: X-axis label string
                            -y: Y-axis label string
                            -title: Axis title string
                            -units: Arrow legend units
    maintainaspect[=True]: Boolean whether to keep the aspect ration constant
    


    Returns: Figure object if ax=False [default], Void otherwise.
    """
    from uedge import com
    from matplotlib.patches import Polygon
    from matplotlib.colors import Normalize
    from matplotlib.pyplot import get_cmap,colorbar,figure
    from numpy import zeros,where,concatenate,reshape, linspace
    from matplotlib.cm import ScalarMappable

    
    if not ax:
        ret=figure()
        ax=ret.add_subplot(111)




    # TODO :  Rewrite so that zlim can use same indices to find local mins and maxes for more efficient plotting
    # Set zoom location
    if zoom=="ot":
        xlim = [com.rm[com.ixpt1[0]-2,0,0]*0.99,com.rm[-1,-1,0]*1.01]
        ylim = [com.zm.min()*0.99,com.zm[com.ixpt1[0],com.iysptrx+1,0]*1.01]
    elif zoom=="it":
        xlim = [com.rm[0,-1,0]*0.99,com.rm[com.ixpt2[0]+4,0,0]*1.01]
        ylim = [com.zm.min()*0.99,com.zm[0,-1,0]*1.01]
    elif zoom=="div":
        xlim = [com.rm[0,-1,0].min()*0.99,com.rm[-1,-1,0].max()*1.01]
        ylim = [com.zm.min()*0.99,com.zm[com.ixpt1[0],com.iysptrx+1,0]*1.01]
    elif zoom=="device":
        xlim,ylim= [com.rm.min(),com.rm.max()],[com.zm.min(),com.zm.max()]

    # Set heatmap limits if requested
    Zmax=Z[1:-1,1:-1].max()
    Zmin=Z[1:-1,1:-1].min()
    if abs(Zmin-Zmax)<1e-3:
       Zmin,Zmax=Zmin*0.95,Zmax*1.05 
       levels=10
    # Set levels
    levels=linspace(Zmin,Zmax,3*com.nx)
    # Set colormap
    cmap=get_cmap(cmap)
    



    Zt=ue_interpolate(Z) # Interpolate values to nodes and remove GC:s
    Z=zeros((com.nx+2,com.ny+2,5))  # Create padded array
    Z[1:-1,1:-1,:]=Zt    # Array same shape as R & Z
    pl=zeros((com.nx+1,com.ny*2+1,3))

    

    # ADD CC AND LOWER NODES FOR ONE DIM, APPEND UPPER ONTOP
    for var in range(3):
    # Coose variable to parse
        if var==0:
            X=com.rm
        elif var==1:
            X=com.zm
        else:
            X=Z
    
        for i in range(1,com.nx+1): 
            for j in range(1,com.ny+1): 
                for k in range(2): 
                    pl[i-1,2*(j-1)+k,var]=X[i,j,k]
                    pl[i-1,2*(j-1)+k,var]=X[i,j,k] 

        pl[:,-1,var]=X[1:,-2,3]
        
        for j in range(com.ny):
            for k in range(2):
                pl[-1,2*j+k,var]=X[-2,j,k+2]

        pl[-1,-1,var]=X[-2,-2,4]


    pfr=concatenate((pl[:com.ixpt1[0],:(com.iysptrx+1)*2+1,:],pl[com.ixpt2[0]:,:(com.iysptrx+1)*2+1,:]),axis=0)
    sol=pl[:,2*(com.iysptrx):,:]
    core=concatenate((pl[com.ixpt1[0]:com.ixpt2[0],:2*(com.iysptrx+1),:],reshape(pl[com.ixpt1[0],:2*(com.iysptrx+1),:],(1,len(pl[0,:2*(com.iysptrx+1),0]),3))),axis=0)


    ax.contourf(pfr[:,:,0],pfr[:,:,1],pfr[:,:,2],levels,cmap=cmap)#,extend="both")
    ax.contourf(core[:,:,0],core[:,:,1],core[:,:,2],levels,cmap=cmap)#,extend="both")
    ax.contourf(sol[:,:,0],sol[:,:,1],sol[:,:,2],levels,cmap=cmap)#,extend="both")




    # Plot vessel if requested
    vesselparams={"plotoutline":True,"plotsep":True,"plotves":True}
    if plotvessel not in [True,False]: # Plot outlines if requested
        if "outline" not in plotvessel:
            vesselparams["plotoutline"]=False
        if "sep" not in plotvessel:
            vesselparams["plotsep"]=False
        if "vessel" not in plotvessel:
            vesselparams["plotves"]=False
    if plotvessel is not False:
        pltves(ax,**vesselparams)

    # Set colorbar if requested
    if cbar is True:
        norm = Normalize(vmin=Zmin,vmax=Zmax)
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

    cbar=colorbar(sm,ax=ax,extend="neither")

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    # Check if aspect ratio should be preserved
    if maintainaspect is True:
        ax.set_aspect('equal', 'datalim')


    # Set additional labels as requested
    if labels is not False:
        if "x" in labels:
            ax.set_xlabel(labels["x"])
        if "y" in labels:
            ax.set_ylabel(labels["y"])
        if "title" in labels:
            ax.set_title(labels["title"])
        if "units" in labels:
            cbar.set_label(labels["units"])


    if "ret" in vars():
        ret.show()
        return ret 



def species(Y,s):
    ''' Routine checking for multispecies arrays, returning the requested or default species, with warning if species not requested '''
    if len(Y.shape)==3:
        if s is None:
            print('WARNING! Multi-species array requested. Plotting species index 0.')
            return Y[:,:,0]
        else:
            return Y[:,:,s]
    else:
        return Y




def radialplotter(X,Y,ax,line,marker,color,linewidth,markersize,xlabel,ylabel,title,xlim,ylim,yaxis,legend):
    """ Routine containing main radial plotting options """
    from matplotlib.pyplot import figure

    # See if we are to add to axis or plot new fig
    if ax is False:
        ret=figure()
        ax=ret.add_subplot(111)
    else:
        try:
            ax=ax.get_axes()[0]
        except:
            pass

    # Plot
    if yaxis=="lin":
        ax.plot(X[1:-1],Y[1:-1],color=color,marker=marker,linestyle=line,linewidth=linewidth,markersize=markersize)
    elif yaxis=="log":
        ax.semilogy(X[1:-1],Y[1:-1],color=color,marker=marker,linestyle=line,linewidth=linewidth,markersize=markersize)
    else:
        print("Specify yaxis to be either 'lin' or 'log'!")
        return

    # Mark SEP
    if X.min()<0 and X.max()>0:
        ax.axvline(0,color="k")
    else:
        ax.axvline(1,color="k")

    # Set limits if requested
    if xlim is not False:
        ax.set_xlim(xlim)
    if ylim is not False:
        ax.set_ylim(ylim)
    

    # Plot the requested labels
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.get_figure().suptitle(title)

    # Set legend if requested
    if legend is not False:
        ax.legend(legend,loc='best',frameon=False,handlelength=3,facecolor=None)

    # Rerturn the figure if drawn
    if "ret" in vars():
        return ret 




def it(Y,s=None,ax=False,line='-',marker='o',color='k',linewidth=2,markersize=8,xlabel=None,ylabel=False,title="IT profile",xlim=False,ylim=False,yaxis="lin",legend=False,xaxis='it'):
    """ Plots the specified parameter along the inner target
    it(var,**keys)

    Variables:
    var:            Variable to be plotted at inner target

    Keyword parameters
    ax[=False]:     Axis where to plot. If False, new figure is created. If ax is figure object, plots on first axes object. 
                    If ax is axes object, plots on the requested axes.
    line[='-']:     Line style
    marker[='o']    Marker style
    color[='k']     Marker and line color
    linewidth[=2]   Plotted line width
    markersize[=8]  Plotted marker size
    xlabel:         String shown as x-axis label
    ylabel:         String shown as y-axis label
    title:          Plot title
    xlim:           Tuple containing X-axis limits
    ylim:           Tuple containing Y-axis limits
    yaxis:          Y-axis scale: 'lin' or 'log'
    xaxis[='it']    X-axis selector. Options are inner target ('it'), outer target ('ot'), midplane ('mp'),
                    or normalized Psi ('psi')
    s[=0]           Species index if multi-species array set as variable
    legend:         List of label titles, in order that lines were plotted


    Returns: Figure object if ax=False [default], Void otherwise.
    """
    from uedge import com
    if xaxis=='it':
        x=com.yylb[:,0]
        if xlabel is None:
            xlabel='Distance along IT [m]'
    elif xaxis=='ot':
        x=com.yyrb[:,0]
        if xlabel is None:
            xlabel='Distance along OT [m]'
    elif xaxis=='mp':
        x=com.yyc
        if xlabel is None:
            xlabel='Distance projected onto LFS-MP [m]'
    elif xaxis=='psi':
        x=com.psinormc
        if xlabel is None:
            xlabel=r'$\rm{\Psi_{norm}}$ [-]'
    else:
        return "Projection must be set to 'mp', 'it', 'ot', or 'psi'!"
        

    Y=species(Y,s)[1,:]

    ret=radialplotter(x,Y,ax,line,marker,color,linewidth,markersize,xlabel,ylabel,title,xlim,ylim,yaxis,legend)

    if ax is False:
        ret.show()
        return ret

def ot(Y,s=None,ax=False,line='-',marker='o',color='k',linewidth=2,markersize=8,xlabel=None,ylabel=False,title="OT profile",xlim=False,ylim=False,yaxis="lin",legend=False,mp=False,psi=False,xaxis='ot'):
    """ Plots the specified parameter along the outer target
    ot(var,**keys)

    Variables:
    var:            Variable to be plotted at outer target

    Keyword parameters
    ax[=False]:     Axis where to plot. If False, new figure is created. If ax is figure object, plots on first axes object. 
                    If ax is axes object, plots on the requested axes.
    line[='-']:     Line style
    marker[='o']     Marker style
    mp[=False]:     Project the target profiles to the midplane
    color[='k']     Marker and line color
    linewidth[=2]   Plotted line width
    markersize[=8]  Plotted marker size
    xlabel:         String shown as x-axis label
    ylabel:         String shown as y-axis label
    title:          Plot title
    xlim:           Tuple containing X-axis limits
    ylim:           Tuple containing Y-axis limits
    yaxis:          Y-axis scale: 'lin' or 'log'
    xaxis[='ot']    X-axis selector. Options are inner target ('it'), outer target ('ot'), midplane ('mp'),
                    or normalized Psi ('psi')
    s[=0]           Species index if multi-species array set as variable


    Returns: Figure object if ax=False [default], Void otherwise.
    """
    from uedge import com
    if xaxis=='it':
        x=com.yylb[:,0]
        if xlabel is None:
            xlabel='Distance along IT [m]'
    elif xaxis=='ot':
        x=com.yyrb[:,0]
        if xlabel is None:
            xlabel='Distance along OT [m]'
    elif xaxis=='mp':
        x=com.yyc
        if xlabel is None:
            xlabel='Distance projected onto LFS-MP [m]'
    elif xaxis=='psi':
        x=com.psinormc
        if xlabel is None:
            xlabel=r'$\rm{\Psi_{norm}}$ [-]'
    else:
        return "Projection must be set to 'mp', 'it', 'ot', or 'psi'!"


    Y=species(Y,s)[-2,:]

    ret=radialplotter(x,Y,ax,line,marker,color,linewidth,markersize,xlabel,ylabel,title,xlim,ylim,yaxis,legend)

    if ax is False:
        ret.show()
        return ret


def mp(Y,s=None,ax=False,line='-',marker='o',color='k',linewidth=2,markersize=8,xlabel=None,ylabel=False,title="Midplane profile",xlim=False,ylim=False,yaxis="lin",legend=False,xaxis='mp'):
    """ Plots the specified parameter along the midplane
    mp(var,**keys)

    Variables:
    var:            Variable to be plotted at the midplane

    Keyword parameters
    ax[=False]:     Axis where to plot. If False, new figure is created. If ax is figure object, plots on first axes object. 
                    If ax is axes object, plots on the requested axes.
    line[='-']:     Line style
    marker[='o']     Marker style
    color[='k']     Marker and line color
    linewidth[=2]   Plotted line width
    markersize[=8]  Plotted marker size
    xlabel:         String shown as x-axis label
    ylabel:         String shown as y-axis label
    title:          Plot title
    xlim:           Tuple containing X-axis limits
    ylim:           Tuple containing Y-axis limits
    yaxis:          Y-axis scale: 'lin' or 'log'
    xaxis[='mp']    X-axis selector. Options are inner target ('it'), outer target ('ot'), midplane ('mp'),
                    or normalized Psi ('psi')
    s[=0]           Species index if multi-species array set as variable
    legend:         List of label titles, in order that lines were plotted


    Returns: Figure object if ax=False [default], Void otherwise.
    """
    from uedge import bbb,com
    if xaxis=='it':
        x=com.yylb[:,0]
        if xlabel is None:
            xlabel='Distance along IT [m]'
    elif xaxis=='ot':
        x=com.yyrb[:,0]
        if xlabel is None:
            xlabel='Distance along OT [m]'
    elif xaxis=='mp':
        x=com.yyc
        if xlabel is None:
            xlabel='Distance along LFS-MP [m]'
    elif xaxis=='psi':
        x=com.psinormc
        if xlabel is None:
            xlabel=r'$\rm{\Psi_{norm}}$ [-]'
    else:
        return "Projection must be set to 'mp', 'it', 'ot', or 'psi'!"

    Y=species(Y,s)[bbb.ixmp,:]

    ret=radialplotter(x,Y,ax,line,marker,color,linewidth,markersize,xlabel,ylabel,title,xlim,ylim,yaxis,legend)

    if ax is False:
        ret.show()
        return ret

def row(Y,row,s=None,ax=False,line='-',marker='o',color='k',linewidth=2,markersize=8,xlabel=None,ylabel=None,title=None,xlim=False,ylim=False,yaxis="lin",xaxis="mp",legend=False):
    """ Plots the specified parameter along requested row projected onto the requested location
    mp(var,row,**keys)

    Variables:
    var:                Variable to be plotted at the requested row
    row:                Index of row to be plotted

    Keyword parameters
    ax[=False]:         Axis where to plot. If False, new figure is created. If ax is figure object, plots on first axes object. 
                        If ax is axes object, plots on the requested axes.
    line[='-']:         Line style
    marker[='o']        Marker style
    color[='k']         Marker and line color
    linewidth[=2]       Plotted line width
    markersize[=8]      Plotted marker size
    xlabel:             String shown as x-axis label
    ylabel:             String shown as y-axis label
    title:              Plot title
    xlim:               Tuple containing X-axis limits
    ylim:               Tuple containing Y-axis limits
    yaxis:              Y-axis scale: 'lin' or 'log'
    xaxis[='mp']    X-axis selector. Options are inner target ('it'), outer target ('ot'), midplane ('mp'),
                    or normalized Psi ('psi')
    s[=0]               Species index if multi-species array set as variable
    legend:         List of label titles, in order that lines were plotted


    Returns: Figure object if ax=False [default], Void otherwise.
    """
    """ Plots the specified parameter Y along the outer target


    Returns: Figure object if ax=False [default], Void otherwise.
    """
    from uedge import bbb,com
    if xaxis=='it':
        x=com.yylb[:,0]
        if xlabel is None:
            xlabel='Distance along IT [m]'
    elif xaxis=='ot':
        x=com.yyrb[:,0]
        if xlabel is None:
            xlabel='Distance along OT [m]'
    elif xaxis=='mp':
        x=com.yyc
        if xlabel is None:
            xlabel='Distance along LFS-MP [m]'
    elif xaxis=='psi':
        x=com.psinormc
        if xlabel is None:
            xlabel=r'$\rm{\Psi_{norm}}$ [-]'
    else:
        return "Projection must be set to 'mp', 'it', 'ot', or 'psi'!"

    if title is None:
        title='Radial profile for row={}'.format(row)

    Y=species(Y,s)

    ret=radialplotter(x,Y[row,:],ax,line,marker,color,linewidth,markersize,xlabel,ylabel,title,xlim,ylim,yaxis,legend)

    if ax is False:
        ret.show()
        return ret


def ftplotter(X,Y,ax,line,marker,color,linewidth,markersize,xlabel,ylabel,title,xlim,ylim,yaxis,legend):
    from matplotlib.pyplot import figure
    from uedge import com,bbb
    
    # Plot figure unless plotting to axis
    if not ax:
        ret=figure()
        ax=ret.add_subplot(111)
    else:
        try:
            ax=ax.get_axes()[0]
        except:
            pass

    # Plot with option for yaxis
    if yaxis=="lin":
        ax.plot(X[1:-1],Y[1:-1],color=color,marker=marker,linestyle=line,linewidth=linewidth,markersize=markersize)
    elif yaxis=="log":
        ax.semilogy(X[1:-1],Y[1:-1],color=color,marker=marker,linestyle=line,linewidth=linewidth,markersize=markersize)
    else:
        print("Specify yaxis to be either 'lin' or 'log'!")
        return

    # Mark xpt and midplane
    
    # TODO: Weight the cell faces by gx!
    if bbb.mhdgeo>0:
        ax.axvline((X[bbb.ixmp]),color="k")
        ax.axvline(0.5*(X[com.ixpt1[0]]+X[com.ixpt1[0]+1]),color="k")
    ax.axvline(0.5*(X[com.ixpt2[0]]+X[com.ixpt2[0]+1]),color="k")

    # Set axis limits if requested
    if xlim is not False:
        ax.set_xlim(xlim)
    if ylim is not False:
        ax.set_ylim(ylim)
    

    # Show requested labels
    if xlabel is not False:
        ax.set_xlabel(xlabel)
    if ylabel is not False:
        ax.set_ylabel(ylabel)
    if title is not False:
        ax.get_figure().suptitle(title)

    # Set legend if requested
    if legend is not False:
        ax.legend(legend,loc='best',frameon=False,handlelength=3,facecolor=None)


    # Return figure if created
    if "ret" in vars():
        ret.show()
        return ret 

def calcdist(para):
    from uedge import com,bbb
    from numpy import cumsum,ones,sqrt
    
    # Set scaling
    scale=ones((com.nx+2,com.ny+2))
    if para is True:
            scale=(com.bphi[:,:,0]**2+com.bpol[:,:,0]**2)**0.5

    # Calculate the distance in the poloidal direction
    x=ones((com.nx+2,com.ny+2))
    for i in range(com.ny+2):
        x[:,i]=cumsum(scale[:,i]/com.gxf[:,i])
    x[-1,:]=x[-2,:] # Fix GC
    x=x-x[bbb.ixmp,:] # Normalize to OMP
    
    return x


def sep(Y,s=None,ax=False,line='-',marker='o',color='k',linewidth=2,markersize=8,xlabel=False,ylabel=False,title=False,xlim=False,ylim=False,para=False,yaxis="lin",legend=False):
    """ Plots the specified parameter along (First flux tube outside) the separatrix
    sep(var,**keys)

    Variables:
    var:            Variable to be plotted along (just outside) the separatrix

    Keyword parameters
    ax[=False]:     Axis where to plot. If False, new figure is created. If ax is figure object, plots on first axes object. 
                    If ax is axes object, plots on the requested axes.
    line[='-']:     Line style
    marker[='o']     Marker style
    color[='k']     Marker and line color
    linewidth[=2]   Plotted line width
    markersize[=8]  Plotted marker size
    xlabel:         String shown as x-axis label
    ylabel:         String shown as y-axis label
    title:          Plot title
    xlim:           Tuple containing X-axis limits
    ylim:           Tuple containing Y-axis limits
    yaxis:          Y-axis scale: 'lin' or 'log'
    s[=0]           Species index if multi-species array set as variable
    legend:         List of label titles, in order that lines were plotted


    Returns: Figure object if ax=False [default], Void otherwise.
    """
    from uedge import com
    
    Y=species(Y,s)[:,com.iysptrx+1]

    if xlabel is False:
        if para is True:
            xlabel="Parallel distance along separatrix"
        else:
            xlabel="Poloidal distance along separatrix"
    if title is False:
        if para is True:
            title="Parallel distance along separatrix"
        else:
            title="Poloidal distance along separatrix"


    ret=ftplotter(calcdist(para)[:,com.iysptrx+1],Y,ax,line,marker,color,linewidth,markersize,xlabel,ylabel,title,xlim,ylim,yaxis,legend)


    if ax is False:
        ret.show()
        return ret 



def ft(Y,ft,s=None,ax=False,line='-',color='k',marker='o',linewidth=2,markersize=8,xlabel=False,ylabel=False,title=False,xlim=False,ylim=False,para=False,yaxis="lin",legend=False):
    """ Plots the specified parameter along the specified flux tube
    ft(var,ft,**keys)

    Variables:
    var:            Variable to be plotted along the requested flux tube
    ft:             Index of flux tube for which variable is to be plotted

    Keyword parameters
    ax[=False]:     Axis where to plot. If False, new figure is created. If ax is figure object, plots on first axes object. 
                    If ax is axes object, plots on the requested axes.
    line[='-']:     Line style
    marker[='o']     Marker style
    color[='k']     Marker and line color
    linewidth[=2]   Plotted line width
    markersize[=8]  Plotted marker size
    xlabel:         String shown as x-axis label
    ylabel:         String shown as y-axis label
    title:          Plot title
    xlim:           Tuple containing X-axis limits
    ylim:           Tuple containing Y-axis limits
    yaxis:          Y-axis scale: 'lin' or 'log'
    s[=0]           Species index if multi-species array set as variable
    legend:         List of label titles, in order that lines were plotted


    Returns: Figure object if ax=False [default], Void otherwise.
    """
    from uedge import com

    Y=species(Y,s)

    if xlabel is False:
        if para is True:
            xlabel="Parallel distance along FT "+str(ft)
        else:
            xlabel="Poloidal distance along FT "+str(ft)
    if title is False:
        if para is True:
            title="Parallel distance along FT "+str(ft)
        else:
            title="Poloidal distance along FT "+str(ft)


    ret=ftplotter(calcdist(para)[:,ft],Y[:,ft],ax,line,marker,color,linewidth,markersize,xlabel,ylabel,title,xlim,ylim,yaxis,legend)



    if ax is False:
        ret.show()
        return ret 



def plot(X,Y,ax=False,linewidth=2,markersize=8,xlabel=False,ylabel=False,title=False,xlim=False,ylim=False,yaxis="lin",marker='',line='-',color='k',legend=False,xaxis='lin'):
    """ Plots the specified x and y axes, including guard cells
    plot(x,y,**keys)

    Variables:
    x:              Variable to be plotted along the X-axis
    y:              Variable to be plotted along the Y-axis

    Keyword parameters
    ax[=False]:     Axis where to plot. If False, new figure is created. If ax is figure object, plots on first axes object. 
                    If ax is axes object, plots on the requested axes.
    line[='-']:     Line style
    marker[='']     Marker style
    color[='k']     Marker and line color
    linewidth[=2]   Plotted line width
    markersize[=8]  Plotted marker size
    xlabel:         String shown as x-axis label
    ylabel:         String shown as y-axis label
    title:          Plot title
    xlim:           Tuple containing X-axis limits
    ylim:           Tuple containing Y-axis limits
    yaxis:          Y-axis scale: 'lin' or 'log'
    legend:         List of label titles, in order that lines were plotted


    Returns: Figure object if ax=False [default], Void otherwise.
    """
    from matplotlib.pyplot import figure
    from uedge import com,bbb
    
    # Plot figure unless plotting to axis
    if not ax:
        ret=figure()
        ax=ret.add_subplot(111)
    else:
        try:
            ax=ax.get_axes()[0]
        except:
            pass


    # Plot with option for yaxis
    if yaxis=="lin":
        if xaxis=="lin":
            ax.plot(X[1:-1],Y[1:-1],color=color,marker=marker,linestyle=line,linewidth=linewidth,markersize=markersize)
        elif xaxis=="log":
            ax.semilogx(X[1:-1],Y[1:-1],color=color,marker=marker,linestyle=line,linewidth=linewidth,markersize=markersize)
        else:
            print("Specify xaxis to be either 'lin' or 'log'!")
            return
    elif yaxis=="log":
        if xaxis=="lin":
            ax.semilogy(X[1:-1],Y[1:-1],color=color,marker=marker,linestyle=line,linewidth=linewidth,markersize=markersize)
        elif xaxis=="log":
            ax.loglog(X[1:-1],Y[1:-1],color=color,marker=marker,linestyle=line,linewidth=linewidth,markersize=markersize)
        else:
            print("Specify xaxis to be either 'lin' or 'log'!")
            return
    else:
        print("Specify yaxis to be either 'lin' or 'log'!")
        return




    # Set axis limits if requested
    if xlim is not False:
        ax.set_xlim(xlim)
    if ylim is not False:
        ax.set_ylim(ylim)

    # Set legend if requested
    if legend is not False:
        ax.legend(legend,loc='best',frameon=False,handlelength=3,facecolor=None)

    # Show requested labels
    if xlabel is not False:
        ax.set_xlabel(xlabel)
    if ylabel is not False:
        ax.set_ylabel(ylabel)
    if title is not False:
        ax.get_figure().suptitle(title)

    # Return figure if created
    if "ret" in vars():
        ret.show()
        return ret 
   

