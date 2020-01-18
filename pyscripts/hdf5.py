
import numpy as np
import h5py
import uedge
from .uedge import bbb
from .uedge import com
from .uedge_lists import *
import time


def hdf5_restore(file):
    """
        Read a hdf5 file previously written from Uedge. This reads the file and puts
        the 6 standard variables into the correct format.
    """
    try:
        hf = h5py.File(file, 'r')
    except ValueError as error:
        print("Couldn't open hdf5 file ", file)
        print(error)
    except:
        print("Couldn't open hdf5 file ", file)
        return

    try:
        dummy = hf['bbb']   # force an exception if the group not there
        hfb = hf.get('bbb')
        try:
           print('File attributes:')
           print('     written on: ', hfb.attrs['ctime'])
           print('        by code: ', hfb.attrs['code'])
           print('    physics tag: ', np.char.strip(hfb.attrs['ver']))
           print(' python version: ', np.char.strip(hfb.attrs['pyver']))
        except:
            pass
        try:
            bbb.ngs[...] = np.array(hfb.get('ngs'))
        except ValueError as error:
            print("Couldn't read ngs from  ", file)
            print(error)
        except:
            print("Couldn't read ngs from  ", file)
        try:
            bbb.nis[...] = np.array(hfb.get('nis'))
        except ValueError as error:
            print("Couldn't read nis from  ", file)
            print(error)
        except:
            print("Couldn't read nis from  ", file)
        try:
            bbb.phis[...] = np.array(hfb.get('phis'))
        except ValueError as error:
            print("Couldn't read phis from  ", file)
            print(error)
        except:
            print("Couldn't read phis from  ", file)
        try:
            bbb.tes[...] = np.array(hfb.get('tes'))
        except ValueError as error:
            print("Couldn't read tes from  ", file)
            print(error)
        except:
            print("Couldn't read tes from  ", file)
        try:
            bbb.tis[...] = np.array(hfb.get('tis'))
        except ValueError as error:
            print("Couldn't read tis from  ", file)
            print(error)
        except:
            print("Couldn't read tis from  ", file)
        try:
            bbb.ups[...] = np.array(hfb.get('ups'))
        except ValueError as error:
            print("Couldn't read ups from  ", file)
            print(error)
        except:
            print("Couldn't read ups from  ", file)
        try:
            bbb.tgs[...] = np.array(hfb.get('tgs'))
        except ValueError as error:
            print("Couldn't read tgs from  ", file)
            print(error)
        except:
            print("Couldn't read tgs from  ", file)

    except:
        print("Old style hdf5 file")
        try:
            bbb.ngs[...] = np.array(hf.get('ngs@bbb'))
        except ValueError as error:
            print("Couldn't read ngs from  ", file)
            print(error)
        except:
            print("Couldn't read ngs from  ", file)
        try:
            bbb.nis[...] = np.array(hf.get('nis@bbb'))
        except ValueError as error:
            print("Couldn't read nis from  ", file)
            print(error)
        except:
            print("Couldn't read nis from  ", file)
        try:
            bbb.phis[...] = np.array(hf.get('phis@bbb'))
        except ValueError as error:
            print("Couldn't read phis from  ", file)
            print(error)
        except:
            print("Couldn't read phis from  ", file)
        try:
            bbb.tes[...] = np.array(hf.get('tes@bbb'))
        except ValueError as error:
            print("Couldn't read tes from  ", file)
            print(error)
        except:
            print("Couldn't read tes from  ", file)
        try:
            bbb.tis[...] = np.array(hf.get('tis@bbb'))
        except ValueError as error:
            print("Couldn't read tis from  ", file)
            print(error)
        except:
            print("Couldn't read tis from  ", file)
        try:
            bbb.ups[...] = np.array(hf.get('ups@bbb'))
        except ValueError as error:
            print("Couldn't read ups from  ", file)
            print(error)
        except:
            print("Couldn't read ups from  ", file)
        try:
            bbb.tgs[...] = np.array(hf.get('tgs@bbb'))
        except ValueError as error:
            print("Couldn't read tgs from  ", file)
            print(error)
        except:
            print("Couldn't read tgs from  ", file)

    hf.close()


def hdf5_save(file):
    """
    Save HDF5 output for restarting and plotting.
    """
    try:
        hf = h5py.File(file, 'w')
        hfb = hf.create_group('bbb')
        hfb.attrs['time'] = time.time()
        hfb.attrs['ctime'] = time.ctime()
        hfb.attrs['code'] = 'UEDGE'
        hfb.attrs['ver'] = bbb.uedge_ver
        try:
           hfb.attrs['pyver'] = uedge.__version__
        except:
           pass
    except ValueError as error:
        print("HDF5 file open failed to ", file)
        print(error)
    except:
        print("HDF5 file open failed to ", file)
        raise
    try:
        d = hfb.create_dataset('ngs', data=bbb.ngs)
        d.attrs['units'] = bbb.getvarunit('ngs')
        d.attrs['comment'] = bbb.getvardoc('ngs')
        d = hfb.create_dataset('ng', data=bbb.ng)
        d.attrs['units'] = bbb.getvarunit('ng')
        d.attrs['comment'] = bbb.getvardoc('ng')
    except ValueError as error:
        print("ng or ngs HDF5 write failed to ", file)
        print(error)
    except:
        print("ng or ngs HDF5 write failed to ", file)
    try:
        d = hfb.create_dataset('nis', data=bbb.nis)
        d.attrs['units'] = bbb.getvarunit('nis')
        d.attrs['comment'] = bbb.getvardoc('nis')
        d = hfb.create_dataset('ni', data=bbb.ni)
        d.attrs['units'] = bbb.getvarunit('ni')
        d.attrs['comment'] = bbb.getvardoc('ni')
    except ValueError as error:
        print("ni or nis HDF5 write failed to ", file)
        print(error)
    except:
        print("ni or nis HDF5 write failed to ", file)
    try:
        d = hfb.create_dataset('phis', data=bbb.phis)
        d.attrs['units'] = bbb.getvarunit('phis')
        d.attrs['comment'] = bbb.getvardoc('phis')
        d = hfb.create_dataset('phi', data=bbb.phi)
        d.attrs['units'] = bbb.getvarunit('phi')
        d.attrs['comment'] = bbb.getvardoc('phi')
    except ValueError as error:
        print("phi or phis HDF5 write failed to ", file)
        print(error)
    except:
        print("phi or phis HDF5 write failed to ", file)
    try:
        d = hfb.create_dataset('tes', data=bbb.tes)
        d.attrs['units'] = bbb.getvarunit('tes')
        d.attrs['comment'] = bbb.getvardoc('tes')
        d = hfb.create_dataset('te', data=bbb.te)
        d.attrs['units'] = bbb.getvarunit('te')
        d.attrs['comment'] = bbb.getvardoc('te')
    except ValueError as error:
        print("te or tes HDF5 write failed to ", file)
        print(error)
    except:
        print("te or tes HDF5 write failed to ", file)
    try:
        d = hfb.create_dataset('tis', data=bbb.tis)
        d.attrs['units'] = bbb.getvarunit('tis')
        d.attrs['comment'] = bbb.getvardoc('tis')
        d = hfb.create_dataset('ti', data=bbb.ti)
        d.attrs['units'] = bbb.getvarunit('ti')
        d.attrs['comment'] = bbb.getvardoc('ti')
    except ValueError as error:
        print("ti or tis HDF5 write failed to ", file)
        print(error)
    except:
        print("Couldn't write tis to  ", file)
    try:
        d = hfb.create_dataset('ups', data=bbb.ups)
        d.attrs['units'] = bbb.getvarunit('ups')
        d.attrs['comment'] = bbb.getvardoc('ups')
        d = hfb.create_dataset('up', data=bbb.up)
        d.attrs['units'] = bbb.getvarunit('up')
        d.attrs['comment'] = bbb.getvardoc('up')
    except ValueError as error:
        print("up or ups HDF5 write failed to ", file)

    except:
        print("up or ups HDF5 write failed to ", file)
    try:
        d = hfb.create_dataset('tgs', data=bbb.tgs)
        d.attrs['units'] = bbb.getvarunit('tgs')
        d.attrs['comment'] = bbb.getvardoc('tgs')
        d = hfb.create_dataset('tg', data=bbb.tg)
        d.attrs['units'] = bbb.getvarunit('tg')
        d.attrs['comment'] = bbb.getvardoc('tg')
    except ValueError as error:
        print("tg or tgs HDF5 write failed to ", file)
        print(error)
    except:
        print("tg or tgs HDF5 write failed to ", file)
    try:
        d = hfb.create_dataset('ev', data=bbb.ev)
        d.attrs['units'] = bbb.getvarunit('ev')
        d.attrs['comment'] = bbb.getvardoc('ev')
    except ValueError as error:
        print("ev HDF5 write failed to ", file)
        print(error)
    except:
        print("ev HDF5 write failed to ", file)

    try:
        hfc = hf.create_group('com')
        hfc.create_dataset('nx', data=com.nx)
        hfc.create_dataset('ny', data=com.ny)
        hfc.create_dataset('rm', data=com.rm)
        hfc.create_dataset('zm', data=com.zm)
    except ValueError as error:
        print("com HDF5 write failed to ", file)
        print(error)
    except:
        print("com HDF5 write failed to ", file)

    hf.close()


def hdf5_dump(file, packages=list_packages(objects=1), vars=None, globals=None):
    """
       Dump all variables from a list of package objects into a file.
       Default packages are output of uedge.uedge_lists.list_packages() 
       vars=[varlist] dump limited to intersection of varlist and packages
    """
    try:
        hf = h5py.File(file, 'w')
    except ValueError as error:
        print("Couldn't open hdf5 file ", file)
        print(error)
    except:
        print("Couldn't open hdf5 file ", file)
    for p in packages:
        hfg = hf.create_group(p.name())
        hfg.attrs['time'] = time.time()
        hfg.attrs['ctime'] = time.ctime()
        hfg.attrs['code'] = 'UEDGE'
        hfg.attrs['ver'] = bbb.uedge_ver
        try:
           hfg.attrs['pyver'] = uedge.__version__
        except:
           pass
        for v in list_variables(p, vars=vars):
            if p.allocated(v):
                try:
                    d = hfg.create_dataset(v, data=p.getpyobject(v))
                    d.attrs['units'] = bbb.getvarunit(v)
                    d.attrs['comment'] = bbb.getvardoc(v)
                except ValueError as error:
                    print("Couldn't write out: "+p.name()+'.'+v)
                    print(error)
                except:
                    print("Couldn't write out: "+p.name()+'.'+v)
            else:
                print(p.name()+'.'+v+" is not allocated")
    if globals != None:
        hfg = hf.create_group('globals')
        hfg.attrs['time'] = time.time()
        hfg.attrs['ctime'] = time.ctime()
        hfg.attrs['code'] = 'UEDGE'
        hfg.attrs['ver'] = bbb.uedge_ver
        try:
           hfg.attrs['pyver'] = uedge.__version__
        except:
           pass
        for v in list(set(globals.keys()) & set(vars)):
            try:
                d = hfg.create_dataset(v, data=globals[v])
                d.attrs['units'] = 'none'
                d.attrs['comment'] = 'Global Variable'
            except ValueError as error:
                print("Couldn't write out: "+p.name()+'.'+v)
                print(error)
            except:
                print("Couldn't write out: "+p.name()+'.'+v)

    hf.close()


def hdf5_restore_dump(file, vars=None, scope=globals()):
    """
       Dump all variables from a list of package objects into a file.
       Default packages are output of uedge.uedge_lists.list_packages() 
       vars=[varlist] dump limited to intersection of varlist and packages
    """
    prfileattrs = True
    try:
        hf = h5py.File(file, 'r')
    except:
        print("Couldn't open hdf5 file ", file)
    for gn in hf.keys():
        g = hf[gn]
        try:
            if prfileattrs:
                prfileattrs = False
                print('File attributes:')
                print('     written on: ', g.attrs['ctime'])
                print('        by code: ', g.attrs['code'])
                print( '       version: ', np.char.strip(g.attrs['ver']))
                print('    physics tag: ', np.char.strip(g.attrs['ver']))
                print(' python version: ', np.char.strip(g.attrs['pyver']))
        except:
            print('No file attributes, trying to restore')
        if vars == None:
            varlist = g.keys()
        else:
            varlist = list(set(g.keys()) & set(vars))
        for v in varlist:
            print(v)
    hf.close()
