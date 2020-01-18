# Holm10 Nov 5 2019, created from scratch

def natsort(l): 
    from re import split
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


def create_im(var,plotter,figsize,subplot=None):
    ''' Creates a placeholder image, converts it into a PIL_image which is 
        returned as an image object
        Parameters:
            var         Dictionary with data to plot:
                            z:      List of parameters to be plotted
                                    len(z)<= product of dims
                            x:      List of parameters to be plotted for x-directional
                                    quantities if z not present
                            y:      List of parameters to be plotted for y-directional
                                    quantities if z not present
                            keys:   Keys to pass to plotting function
                            dims:   Tuple of subplot grid dimensions
            plotter     Function used for plotting
            figsize     Figure size tuple
    '''
    from PIL import Image
    from io import BytesIO
    from matplotlib.pyplot import figure,imshow,get_current_fig_manager,close
    # Create a placeholder figure and axis
    tmpfig=figure(figsize=figsize)

    # Check if subplots are requested
    ax=[]   # 
    if subplot is None:    # No
        ax.append(tmpfig.add_subplot(111))
    else:                           # Yes
        (i,j)=subplot   # Exctract subplot grid dimensions
        k=1
        for a in range(i):
            for b in range(j):
                ax.append(tmpfig.add_subplot(i,j,k))    # Create subplot grid 
                k+=1                            # Increase index of axis

    tmpfig.patch.set_facecolor('w') # Set the surrounding colors to be white for prettier animation
    try:
        z=var.pop('z')  # Extract data to be plotted from keys
        onearg=True
    except:
        try:
            x=var.pop('x')
        except:
            return 'Neither "z" nor "x" plotting objects found! Aborting...'
        try:
            y=var.pop('y')
        except:
            return 'No "y" plotting object found! Aborting...'
        if len(x)!=len(y):
            print('X and Y have different dimensions! Aborting...')
            retur -1
        onearg=False

    if onearg is True:
        for axis in range(len(z)):  # Loop through list of data
            plotter(z[axis],ax=ax[axis],**var) # Add a z-plot
    else:
        for axis in range(len(x)):
            plotter(x[axis],y[axis],ax=ax[axis],**var)  # Add a 2-directional plot

    canvas = get_current_fig_manager().canvas # Create a figure from canvas
    canvas.draw()
    pil_image = Image.frombytes('RGB', canvas.get_width_height(), canvas.tostring_rgb())    # Turn image into PIL file
    pil_image.save(BytesIO(), 'PNG') # Store PIL in buffer
    close() # Close pil
    return imshow(pil_image,animated=True) # Create image from PIL



def save_animation(ani,savename,output,fps,bitrate=1800):
    ''' Function saving the movie to the requested file type.
        Paramters
            ani         Animation object to be saved
            savename    Name (incl. path) of movie to be saved w/o file extension
            output      Output type: 'mp4' and 'gif' currently supported
            fps         Framerate requested
        Optional
            bitrate     Bitrate (1800 default)
    '''
    print(savename)
    from getpass import getuser
    from matplotlib.animation import writers
    if output=='mp4':
        Writer = writers['ffmpeg']
        writer = Writer(fps=fps, bitrate=bitrate,  metadata=dict(artist=getuser()))
        ani.save(savename+".mp4",writer=writer)

    elif output=='gif':
        Writer = writers['pillow']
        writer = Writer(fps=fps, metadata=dict(artist=getuser()), bitrate=1800)
        ani.save(savename+".gif",writer=writer)
    else:
        print( "Unknown format requested! Movie is not saved!")
                    



def animate_heatmapdt(casedir,varible,show=True,savename=False,interval=1000,figsize=(1.618*6,6),keys={},commands=[],steps="t="):
    """ Function creating a time-dependent heatmap animation based on file

    animate_heatmapdt(casedir,var,keys**)   
        Arguments:
        casedir:            Path (absolute or relative) to the directory containing the data and solutions folders. 
                            The data folder should contain the input file 'input.py'. The solutions folder should contain
                            an ascii-file named label_dtout.txt, where label is the variable label[0] defined in 'input.py'.
                            The dtout-file should be a csv file with the name of the savefile and the corresponding
                            time of the save.
        variable:                String of the variable that should be plotted (e.g 'bbb.ne', 'bbb.tg[:,:,1]/1.602e-19')
        
        Keyword arguments
        show[=True]:        Boolean determining whether to show plot or not
        savename[=False]:   If animation is to be saved, save must be set as path and save file name
        interval[=1000]:    Time in ms each frame is shown
        figsize             Tuple containing the figure width and height
        keys:               A dictionary containing the keyword arguments of heatmap
        commands:           List of strings containing UEDGE commands to be executed before evaluation,
                            e.g. turning coefficients on or off
        steps:              String with text to displayed in front of second col values in animation title

    """
    from pandas import read_csv
    from os import chdir,getcwd
    from uedge.contrib.ue_plot import heatmap
    from uedge.hdf5 import hdf5_restore
    from PIL import Image
    from StringIO import StringIO
    from matplotlib.pyplot import figure,imshow,show,get_current_fig_manager,close,axes,subplots_adjust
    from matplotlib.animation import ArtistAnimation,writers
    from getpass import getuser

    cwd=getcwd()    # Get original directory
    chdir(casedir+"/data")  # Go to data directory of specified case: required as required files might be stored here
    # Execute and store variables globally 
    import input as i
    i.restore_input()

    # Loop through list of commands in case coefficients need to be changed
    for i in commands:
        exec(i) in globals(),locals()
    
    
    svs=read_csv("../solutions/"+bbb.label[0].decode('UTF-8')+"_dtout.txt") # Get list of saves to animate
    # Output list
    print("Save repository and save times:")
    print(svs)
    svs=svs.values  # Create array

    # Extract original title for title timestamp
    if "labels" in keys.keys():
        if "title" in keys["labels"]:
            origtitle=keys["labels"]["title"]
    else:
        origtitle=""
    

    # Fig for dislpay
    fig=figure(figsize=figsize)
    fig.patch.set_visible(False)

    ims=[]  # Empty array to store frames 
    for frame in range(len(svs)):

        # Add timestamp to plot title
        if "labels" in keys.keys():
            if "title" in keys["labels"]:
                keys["labels"]["title"]=origtitle+" "+steps+"%.2E" % svs[frame,1]
            else:
                keys["labels"]["title"]=steps+"%.2E" % svs[frame,1]
        else:
            keys["labels"]={"title":steps+"%.2E" % svs[frame,1]}

        # Recover all parameters from case
        print("Frame "+str(frame+1)+" of "+str(len(svs)))
        print( "===========================")
        hdf5_restore("../solutions/"+svs[frame,0])
        bbb.ftol=1e20
        bbb.issfon=0
        bbb.exmain()

        # This is a dirty solution, but as the Forthon objects are only pointed to by python the values stored cannot (to my knowledge) be directly accessed by python (list of vars, etc). Instead one must execute a string statement to get the parameter passed onto the function

        exec('ims.append([create_im('+variable+'fwd,keys,figsize)])') # Append frame to list

    chdir(cwd) # Move back to original directory
    ani=ArtistAnimation(fig, ims, interval=interval, blit=True, repeat_delay=0) # Create animation

    # Turn off extra pair of axes, make the animation "fullscreen"
    ax=axes()
    ax.axis('off') 
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

    # Show and store as requested
    if savename is not False:
        
        Writer = writers['ffmpeg']
        writer = Writer(fps=round(1/(interval*1e-3)), metadata=dict(artist=getuser()), bitrate=1800)
        ani.save(savename+".mp4",writer=writer)
    if show:
        show()



def animate_vectordt(casedir,varx,vary,show=True,savename=False,interval=1000,figsize=(1.618*6,6),keys={},commands=[],heatmap=False,steps="t=",unitlength=False):
    """ Function creating a vector animation 
    
    animate_vectordt(casedir,var,keys**)   
        Arguments:
        casedir:            Path (absolute or relative) to the directory containing the data and solutions folders. 
                            The data folder should contain the input file 'input.py'. The solutions folder should contain
                            an ascii-file named label_dtout.txt, where label is the variable label[0] defined in 'input.py'.
                            The dtout-file should be a csv file with the name of the savefile and the corresponding
                            time of the save.
        var:                String of the variable that should be plotted (e.g 'bbb.ne', 'bbb.tg[:,:,1]/1.602e-19')
        
        Keyword arguments
        show[=True]:        Boolean determining whether to show plot or not
        savename[=False]:       If animation is to be saved, save must be set as path and save file name
        interval[=1000]:    Time in ms each frame is shown
        figsize             Tuple containing the figure width and height
        keys:        A dictionary containing the keyword arguments of heatmap
                            Must contain entry "var" with string of parameter to be plotted
        commands:           List of strings containing UEDGE commands to be executed before evaluation,
                            e.g. turning coefficients on or off
        steps:              String with text to displayed in front of second col values in animation title

    """
    from pandas import read_csv
    from os import chdir,getcwd
    from uedge.contrib.ue_plot import vector
    from uedge.hdf5 import hdf5_restore
    from PIL import Image
    from StringIO import StringIO
    from matplotlib.pyplot import figure,imshow,show,get_current_fig_manager,close,axes,subplots_adjust
    from matplotlib.animation import ArtistAnimation,writers
    from getpass import getuser


    cwd=getcwd()    # Get original directory
    chdir(casedir+"/data")  # Go to data directory of specified case: required as required files might be stored here
    # Execute and store variables globally 
    import input as i
    i.restore_input()    

    # Loop through list of commands in case coefficients need to be changed
    for i in commands:
        exec(i) in globals(),locals()
    
    
    svs=read_csv("../solutions/"+bbb.label[0].decode('UTF-8')+"_dtout.txt") # Get list of saves to animate
    # Output list
    print( "Save repository and save times:")
    print(svs)
    svs=svs.values  # Create array

    # Extract original title for title timestamp
    if "labels" in keys.keys():
        if "title" in keys["labels"]:
            origtitle=keys["labels"]["title"]
    else:
        origtitle=""
    
         
    if heatmap is not False:
        ex=heatmap["var"]

    # Fig for dislpay
    fig=figure(figsize=figsize)
    fig.patch.set_visible(False)

    ims=[]  # Empty array to store frames 
    for frame in range(len(svs)):

        # Add timestamp to plot title
        if "labels" in keys.keys():
            if "title" in keys["labels"]:
                keys["labels"]["title"]=origtitle+" "+steps+"%.2E" % svs[frame,1]
            else:
                keys["labels"]["title"]=steps+"%.2E" % svs[frame,1]
        else:
            keys["labels"]={"title":steps+"%.2E" % svs[frame,1]}


        # Recover all parameters from case
        print("Frame "+str(frame+1)+" of "+str(len(svs)))
        print( "===========================")
        hdf5_restore("../solutions/"+svs[frame,0])
        bbb.ftol=1e20
        bbb.issfon=0
        bbb.exmain()

        if heatmap is not False:
            
            exec('keys["heatmap"]=heatmap["var"]')

        # This is a dirty solution, but as the Forthon objects are only pointed to by python the values stored cannot (to my knowledge) be directly accessed by python (list of vars, etc). Instead one must execute a string statement to get the parameter passed onto the function

        exec('ims.append([create_im('+varx+','+vary+',keys,figsize)])') # Append frame to list

    chdir(cwd) # Move back to original directory
    ani=ArtistAnimation(fig, ims, interval=interval, blit=True, repeat_delay=0) # Create animation

    # Turn off extra pair of axes, make the animation "fullscreen"
    ax=axes()
    ax.axis('off') 
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

    # Show and store as requested
    if savename is not False:
        Writer = writers['ffmpeg']
        writer = Writer(fps=round(1/(interval*1e-3)), metadata=dict(artist=getuser()), bitrate=1800)
        ani.save(savename+".mp4",writer=writer)
    if show:
        show()



def ani_vec(varx,vary,path='.',norm=None,zoom='div',show_animation=True,savename=False,fps=10,figsize=(1.618*6,6),keys={},steps="",framerate=1,output='gif',database=None,subplot=None,commands=[],s=None):
    """ Function creating an animated vector plot  
        ani_vec(poldata,raddata,keys**)
    
        Arguments:
        poldata:            Variable name of x-directional quantities to be plotted as string
        raddata:            Variable name of y-directional quantities to be plotted as string

        Keyword arguments
        zoom[='div']:       Zoom of animation, options are: divertor ('div'), device ('device'), inner target ('it') or outer target ('ot')
        norm[=None]:        Normalization to an arrow length 
        s[=0]               Index to be plotted in case multispecies arrays are supplied
        path:               Path
        show[=True]:        Boolean determining whether to show plot or not
        save[=False]:       If animation is to be saved, save must be set as path and save file name
        interval[=1000]:    Time in ms each frame is shown
        figsize             Tuple containing the figure width and height
        keys:             A dictionary containing the keyword arguments of ft plot
        commands:           List of strings containing UEDGE commands to be executed before evaluation,
                            e.g. turning coefficients on or off
        steps:              String with text to displayed in front of second col values in animation title

    """
  #  from pandas import read_csv
    from uedge import bbb,com
    from os import chdir,getcwd,walk
    from matplotlib.animation import ArtistAnimation
    from matplotlib.pyplot import axes,subplots_adjust,figure,show,close
    from uedge.contrib.ue_plot import vector
    from importlib import reload
    from uedge.contrib.utils import readcase

    func=vector


    # Common keys
    keys['norm']=norm
    keys['s']=s
    keys['zoom']=zoom

    chdir(path)
    path=getcwd()
    # Get list of subdirectories in path
    dirs=natsort(next(walk(path))[1])
    dirs=dirs[::framerate]
    # Omit supporting file directories
    try:
        dirs.remove('grid')
    except:
        pass
    try:    
        dirs.remove('rates')
    except:
        pass
    try:
        dirs.remove('ignore')
    except:
        pass

    if len(dirs)==0:
        return 'No directories found! Aborting...'
    

    # Extract original title for repeat use
    try:
        origtitle=keys['title']
    except:
        origtitle=''

    # Fig for dislpay
    fig=figure(figsize=figsize)
    fig.patch.set_visible(False)

    ims=[]  # Empty array to store frames 
    
    if database is None:

        # Pass plotfunc-specific parameters to the function

        for frame in range(len(dirs)):
            # Verbose
            print( "Frame ", frame+1, " of ", len(dirs))
            print( "===========================")
        
            readcase('{}/{}'.format(path,dirs[frame]))

            # Execute any requested commands
            for cmd in commands:
                exec(cmd)

            bbb.ftol=1e20;bbb.issfon=0;bbb.exmain()

            # Add timestamp to plot title
            keys['title']=origtitle+' '+steps+str(bbb.ncore[0])

            # This is a dirty solution, but as the Forthon objects are only pointed to by python the values stored cannot (to my knowledge) be directly accessed by python (list of vars, etc). Instead one must execute a string statement to get the parameter passed onto the function

            exec('keys["x"]=['+varx+']')
            exec('keys["y"]=['+vary+']')
            
            ims.append([create_im(keys,func,figsize,subplot=subplot)]) # Append frame to list
    else:
        for frame in range(len(database)):
            # Verbose
            print( "Frame ", frame+1, " of ", len(dirs))
            print( "===========================")
 
            keys['x']=database[frame]['varx']
            keys['y']=database[frame]['vary']

            ims.append([create_im(keys,func,figsize,subplot=subplot)]) # Append frame to list

    # Turn off extra pair of axes, make the animation "fullscreen"
    ax=axes()
    ax.axis('off') 
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    ani=ArtistAnimation(fig, ims, interval=int(1000.0/fps), blit=True, repeat_delay=1000) # Create animation in fig

    chdir(path) # Go back to dir where function called

    if savename is not False:   # Show and store if requested
        save_animation(ani,savename,output,fps)
    if show_animation is True: 
        return ani
    else:
        close()

def ani_hm(variable,path='.',zoom='div',zrange=False,zaxis='lin',norm=None,show_animation=True,savename=False,fps=10,figsize=(1.618*6,6),keys={},steps="",framerate=1,output='gif',database=None,subplot=None,commands=[]):
    """ Function creating a flux-tube plot  animation 
        ani_hm(var,**keys)
        
        Arguments:
        var:                Variable name to be plotted as string
        
        Keyword arguments
        zoom[='div']:       Zoom of animation, options are: divertor ('div'), device ('device'), inner target ('it') or outer target ('ot')
        zrange:             Tuple specifying heatmap range
        zaxis[='lin']:      Tuple specifying heatmpap axis: 'lin' or 'log'
        show[=True]:        Boolean determining whether to show plot or not
        savename[=False]:   If animation is to be saved, save must be set as path and save file name
        interval[=1000]:    Time in ms each frame is shown
        figsize             Tuple containing the figure width and height
        keys:               A dictionary containing the keyword arguments of ft plot
        commands:           List of strings containing UEDGE commands to be executed before evaluation,
                            e.g. turning coefficients on or off
        steps:              String with text to displayed in front of second col values in animation title

    """
  #  from pandas import read_csv
    from uedge import bbb,com
    from os import chdir,getcwd,walk
    from matplotlib.animation import ArtistAnimation
    from matplotlib.pyplot import axes,subplots_adjust,figure,show,close
    from uedge.contrib.ue_plot import heatmap
    from uedge.contrib.utils import readcase

    # Set common keys
    keys['zaxis']=zaxis
    keys['zoom']=zoom
    keys['zrange']=zrange

    func=heatmap

    chdir(path)
    path=getcwd()
    # Get list of subdirectories in path
    dirs=natsort(next(walk(path))[1])
    dirs=dirs[::framerate]
    # Omit supporting file directories
    try:
        dirs.remove('grid')
    except:
        pass
    try:
        dirs.remove('rates')
    except:
        pass
    try:
        dirs.remove('ignore')
    except:
        pass

    if len(dirs)==0:
        return 'No directories found! Aborting...'
    # Extract original title for repeat use
    try:
        origtitle=keys['title']
    except:
        origtitle=''

    # Fig for dislpay
    fig=figure(figsize=figsize)
    fig.patch.set_visible(False)

    ims=[]  # Empty array to store frames 
    
    if database is None:

        # Pass plotfunc-specific parameters to the function

        for frame in range(len(dirs)):
            # Verbose
            print( "Frame ", frame+1, " of ", len(dirs))
            print( "===========================")

            readcase('{}/{}'.format(path,dirs[frame]))

            # Execute any requested commands
            for cmd in commands:
                exec(cmd)

            bbb.ftol=1e20;bbb.issfon=0;bbb.exmain()

            # Add timestamp to plot title
            keys['title']=origtitle+' '+steps+str(bbb.ncore[0])

            # This is a dirty solution, but as the Forthon objects are only pointed to by python the values stored cannot (to my knowledge) be directly accessed by python (list of vars, etc). Instead one must execute a string statement to get the parameter passed onto the function
            exec('keys["z"]=['+variable+']')
            
            ims.append([create_im(keys,func,figsize,subplot=subplot)]) # Append frame to list

    else:
        for frame in range(len(database)):
            # Verbose
            print("Frame ", frame+1, " of ", len(dirs))
            print( "===========================")
 
            keys['z']=database[frame]['var']

            ims.append([create_im(keys,func,figsize,subplot=subplot)]) # Append frame to list

    # Turn off extra pair of axes, make the animation "fullscreen"
    ax=axes()
    ax.axis('off') 
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    ani=ArtistAnimation(fig, ims, interval=int(1000.0/fps), blit=True, repeat_delay=1000) # Create animation in fig

    chdir(path) # Go back to dir where function called

    if savename is not False:   # Show and store if requested
        save_animation(ani,savename,output,fps)
    if show_animation is True: 
        return ani
    else:
        close()

def ani_ft(variable,plotft,path='.',yaxis='lin',ylim=False,show_animation=True,save=False,fps=10,figsize=(1.618*6,6),keys={},steps="",framerate=1,output='gif',database=None,subplot=None,commands=[]):
    """ Function creating a flux-tube plot  animation 
        ani_ft(var,ft,**keys)        

        Arguments:
        var:                Variable name to be plotted as string
        ft:                 Index of flux tube to be animated        

        Keyword arguments
        show[=True]:        Boolean determining whether to show plot or not
        save[=False]:       If animation is to be saved, save must be set as path and save file name
        interval[=1000]:    Time in ms each frame is shown
        figsize             Tuple containing the figure width and height
        keys:               A dictionary containing the keyword arguments of ft plot
        commands:           List of strings containing UEDGE commands to be executed before evaluation,
                            e.g. turning coefficients on or off
        steps:              String with text to displayed in front of second col values in animation title

    """
  #  from pandas import read_csv
    from uedge import bbb,com
    from os import chdir,getcwd,walk
    from matplotlib.animation import ArtistAnimation
    from matplotlib.pyplot import axes,subplots_adjust,figure,show,close
    from uedge.contrib.ue_plot import ft
    from uedge.contrib.utils import readcase

    func=ft

    # Set common key params
    keys['ft']=plotft
    keys['yaxis']=yaxis
    keys['ylim']=ylim

    chdir(path)
    path=getcwd()
    # Get list of subdirectories in path
    dirs=natsort(next(walk(path))[1])
    dirs=dirs[::framerate]
    # Omit supporting file directories
    try:
        dirs.remove('grid')
    except:    
        pass
    try:
        dirs.remove('rates')
    except:
        pass
    try:
        dirs.remove('ignore')
    except:
        pass

    if len(dirs)==0:
        return 'No directories found! Aborting...'
    # Extract original title for repeat use
    try:
        origtitle=keys['title']
    except:
        origtitle=''
    
    


    # Fig for dislpay
    fig=figure(figsize=figsize)
    fig.patch.set_visible(False)

    ims=[]  # Empty array to store frames 
    

    if 'ft' not in keys.keys():
        return 'Must pass the flux-tube index to be plotted as key "ft"! Aborting...'


    if database is None:
        for frame in range(len(dirs)):
            # Verbose
            print( "Frame ", frame+1, " of ", len(dirs))
            print( "===========================")

            readcase('{}/{}'.format(path,dirs[frame]))

            # Execute any requested commands
            for cmd in commands:
                exec(cmd)

            bbb.ftol=1e20;bbb.issfon=0;bbb.exmain()

            # Add timestamp to plot title
            keys['title']=origtitle+' '+steps+str(bbb.ncore[0])

            # This is a dirty solution, but as the Forthon objects are only pointed to by python the values stored cannot (to my knowledge) be directly accessed by python (list of vars, etc). Instead one must execute a string statement to get the parameter passed onto the function

            exec('keys["z"]=['+variable+']')
            
            ims.append([create_im(keys,func,figsize,subplot=subplot)]) # Append frame to list
    else:
        for frame in range(len(database)):
            # Verbose
            print( "Frame ", frame+1, " of ", len(dirs))
            print( "===========================")
 
            keys['z']=database[frame]['var']

            ims.append([create_im(keys,func,figsize,subplot=subplot)]) # Append frame to list

    # Turn off extra pair of axes, make the animation "fullscreen"
    ax=axes()
    ax.axis('off') 
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    ani=ArtistAnimation(fig, ims, interval=int(1000.0/fps), blit=True, repeat_delay=1000) # Create animation in fig

    chdir(path) # Go back to dir where function called

    if save is not False:   # Show and store if requested
        save_animation(ani,save,output,fps)
    if show_animation is True: 
        return ani
    else:
        close()


def ani_row(variable,plotrow,ylim=False,yaxis='lin',path='.',show_animation=True,savename=False,fps=10,figsize=(1.618*6,6),keys={},steps="",framerate=1,output='gif',database=None,subplot=None,commands=[]):
    """ Function creating a flux-tube plot  animation 
        ani_row(var,row,**keys)
        
        Arguments:
        var:                Variable name to be plotted as string
        row:            Index of row to be plotted        

        Keyword arguments
        show[=True]:        Boolean determining whether to show plot or not
        save[=False]:       If animation is to be saved, save must be set as path and save file name
        interval[=1000]:    Time in ms each frame is shown
        figsize             Tuple containing the figure width and height
        keys:             A dictionary containing the keyword arguments of ft plot
        commands:           List of strings containing UEDGE commands to be executed before evaluation,
                            e.g. turning coefficients on or off
        steps:              String with text to displayed in front of second col values in animation title

    """
  #  from pandas import read_csv
    from uedge import bbb,com
    from os import chdir,getcwd,walk
    from matplotlib.animation import ArtistAnimation
    from matplotlib.pyplot import axes,subplots_adjust,figure,show,close
    from uedge.contrib.ue_plot import row
    from uedge.contrib.utils import readcase

    func=row


    chdir(path)
    path=getcwd()
    # Get list of subdirectories in path
    dirs=natsort(next(walk(path))[1])
    dirs=dirs[::framerate]
    # Omit supporting file directories
    try:
        dirs.remove('grid')
    except:
        pass
    try:
        dirs.remove('rates')
    except:
        pass
    try:
        dirs.remove('ignore')
    except:
        pass

    if len(dirs)==0:
        return 'No directories found! Aborting...'

    # Extract original title for repeat use
    try:
        origtitle=keys['title']
    except:
        origtitle=''

    # Fig for dislpay
    fig=figure(figsize=figsize)
    fig.patch.set_visible(False)

    ims=[]  # Empty array to store frames 
    
    # Set common key params
    keys['row']=plotrow
    keys['yaxis']=yaxis
    keys['ylim']=ylim
    

    if database is None:

        for frame in range(len(dirs)):
            # Verbose
            print("Frame ", frame+1, " of ", len(dirs))
            print("===========================")

            readcase('{}/{}'.format(path,dirs[frame]))

            # Execute any requested commands
            for cmd in commands:
                exec(cmd)


            bbb.ftol=1e20;bbb.issfon=0;bbb.exmain()

            # Add timestamp to plot title
            keys['title']=origtitle+' '+steps+str(bbb.ncore[0])

            # This is a dirty solution, but as the Forthon objects are only pointed to by python the values stored cannot (to my knowledge) be directly accessed by python (list of vars, etc). Instead one must execute a string statement to get the parameter passed onto the function
            exec('keys["z"]=['+variable+']')
            
            ims.append([create_im(keys,func,figsize,subplot=subplot)]) # Append frame to list
    else:
        for frame in range(len(database)):
            # Verbose
            print("Frame ", frame+1, " of ", len(dirs))
            print("===========================")
 
            keys['z']=database[frame]['var']

            ims.append([create_im(keys,func,figsize,subplot=subplot)]) # Append frame to list
    # Turn off extra pair of axes, make the animation "fullscreen"
    ax=axes()
    ax.axis('off') 
    subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    ani=ArtistAnimation(fig, ims, interval=int(1000.0/fps), blit=True, repeat_delay=1000) # Create animation in fig

    chdir(path) # Go back to dir where function called

    if savename is not False:   # Show and store if requested
        save_animation(ani,savename,output,fps)
    if show_animation is True: 
        return ani
    else:
        close()


