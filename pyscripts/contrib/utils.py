# Utilities created to help emulate BASIS capabilities
# Created from scratch by holm10
#   191121 - Created 
from uedge import *

def readcase(casedir,dirpath='data'):
    '''
    Simple script reading the requested input file from folder casedir

    readcase(casedir,**keys)

    Required parameters:
    casedir     -   The name of the case folder containin dirpath/input.py,
                    or alternatively the path to this folder
    
    Keyword parameters:
    dirpath     -   Path to input.py within casedir, default is "data"
    '''
    from os import chdir,getcwd
    from os.path import exists, abspath, relpath
    from importlib import reload
    import sys

    # Check that path is set up correctly
    if exists(casedir) is False:
        print('Case folder {} not found!'.format(casedir))
        return
    elif exists('{}/{}'.format(casedir,dirpath)) is False:
        print('Input folder{}/{} not found!'.format(casedir,dirpath))
        return
    
    cwd=getcwd()                                        # Store original location
    chdir('{}/{}'.format(casedir,dirpath))  # Move to the requested location
    import input as i                              # This is ugly coding, but it is straightforward and works
    reload(i)                                           # Reload case: should not strictly be necessary
    i.restore_input()                                              # Read the new setup 

    # Fix the path dependencies to be absolute paths for reading case out of data
    newaphdir='{}/{}'.format(abspath(getcwd()),aph.aphdir[0].decode('utf-8'))
    newapidir='{}/{}'.format(abspath(getcwd()),api.apidir[0].decode('utf-8'))
    # Set grid path file in V7.08 and later
    try:
        newaeqdskfname='{}/{}'.format(abspath(getcwd()),com.aeqdskfname[0].decode('utf-8'))
        newgeqdskfname='{}/{}'.format(abspath(getcwd()),com.geqdskfname[0].decode('utf-8'))
    except:
        pass

    chdir(cwd)                                          # Go back to the original directory to enable recursive loads 

    # Shorten the concatenated path to fit into arrays
    aph.aphdir[0]=abspath(newaphdir)
    api.apidir[0]=abspath(newapidir)
    try:
        com.aeqdskfname[0]=abspath(newaeqdskfname)
        com.geqdskfname[0]=abspath(newgeqdskfname)
    except:
        pass


def write_csv(varlist,fname,header=None,probname='',descriptor=''):
    ''' Writes a list of 1D profiles to a csv file '''
    from numpy import savetxt,zeros 
    from time import ctime    

    # Set problem name if not given
    if probname is None:
        probname=fname


    # Find longest array to be written
    maxlen=0
    for col in varlist:
        maxlen=max(maxlen,len(col))
    
    # Create an array to accommodate all variables to be written
    vals=zeros((maxlen,len(varlist)))

    # Set up header
    headerp='Problem name: {}\n{}\n{}\n'.format(probname,ctime(),descriptor)    

    # Populate the array and header
    for ii in range(len(varlist)):
        vals[:len(varlist[ii]),ii]=varlist[ii]
        if header is not None:
            headerp+=(header[ii]+', '*(ii!=len(varlist)-1)).ljust(14)
    
    # Save to CSV file
    savetxt('{}.csv'.format(fname),vals,fmt='%.5E',header=headerp,comments='',delimiter=' , ')


def write_radial_csv(varlist, fname,probname=None,  location='mp',gc=False):
    ''' Writes a CSV file with radial data along columns 

    write_radial_csv(varlist, **keys)

    Required parameters:
    varlist     -   2D list containing the data to be written in the following form:
                    [   [title1,    data1],
                        [title2,    data2],
                        ...
                        [titleN,    dataN] ]
                    Here, title is a string and dataX is the data to be written.
                    The data can be a 1D or 2D array:
                        -1D arrays are directly written to the column
                        -The radial profile specified by location is written to the column

    Keyword parameters:
    probname    -   Problem name to be printed at first row, defaults to bbb.label
    fname       -   File name for printed file (script adds .csv)
    location    -   What location to pic from 2D arrays. The options are:
                        'mp' - midplane [default]
                        'ot' - outer target
                        'it' - inner target
                        int  - integer specifying index
    gc          -   Boolean whether to include guard cells or not
                        False - omits guard cells [default]
                        True  - includes guard cells

    ************************************
    WARNING! Fluxes are taken at the RHS boundary of the evaluated cell, e.g. out
    of the cell in front of the inner target for location='it'. To get the right
    fluxes, specify the flux as a 1D array at the appropriate cell index.
    ************************************
    '''
    from uedge import bbb,com
    from numpy import zeros,savetxt

    if probname==None:
        probname=bbb.label[0].decode('UTF-8')
    
    # Set the poloidal index
    if location=='mp':
        ix=bbb.ixmp
        descriptor='Midplane values'
    elif location=='ot':
        ix=com.nx
        descriptor='Outer target values'
    elif location=='it':
        ix=1
        descriptor='Inner target values'
    elif isinstance(location,int):
        ix=location
        descriptor='Row index {} values'.format(ix)
    else:
        print('Location "{}" not recognized!'.format(location))
        print('Use "mp", "ot", "it" or the row integer.')
        return

    # Initialize array for values
    vals=zeros((len(varlist),com.ny+2))
    header='Problem name: {}\n{}\n'.format(probname,descriptor)    
    

    # Loop throug list of titles and values
    for i in range(len(varlist)):
        if len(varlist[i][1].shape)>2:
            print('Variable {} has more than two dimensions!'.format(varlist[i][0]))
            print('Only supply 1D or 2D arrays. Aborting...')
            return
        # Append header
        header+=(varlist[i][0]+', '*(i!=len(varlist)-1)).ljust(14)
        # Store values
        if len(varlist[i][1].shape)==2:   # Store index from 2D profile
            vals[i,:]=varlist[i][1][ix,:]
        elif len(varlist[i][1].shape)==1: # 1D array, store full array
            vals[i,:]=varlist[i][1]

    # Drop guard cells unless they are requested
    if gc is False:
        vals=vals[:,1:-1]
    
    # Save to CSV file
    savetxt('{}.csv'.format(fname),vals.transpose(),fmt='%.5E',header=header,comments='',delimiter=' , ')


def default_csv_parameters():
    '''
    Default parameters to be written, based on MGroth file
    '''
    from uedge import bbb,com
    return  [   ['psin',    com.psinormc],
                ['dyc',     1/(com.gy+1e-30)],
                ['dyf',     1/(com.gyf+1e-30)],
                ['ne',      bbb.ne],
                ['ni',      bbb.ni[:,:,0]],
                ['te',      bbb.te],
                ['ti',      bbb.ti]    ]


