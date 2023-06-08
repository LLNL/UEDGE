import numpy as np
import h5py
import uedge
try:
    import __version__ as pyv
    pyver = pyv.__version__
except:
    pyver = uedge.__version__
from .uedge import bbb
from .uedge import com
from .uedge_lists import *
import time
from Forthon import packageobject


def hdf5_restore(file):
    """
        Read a hdf5 file previously written from pyUedge. This reads the file recursively 
        and will attempt to restore all datasets. This will restore a file saved by either
        hdf5_save or hdf5_dump.
    """
    try:
        hf = h5py.File(file, 'r')
    except:
        print("Couldn't open hdf5 file ", file)
        raise

    try:
        dummy = hf['bbb']   # force an exception if the group not there
        hfb = hf.get('bbb')
        try:           
            for var in ['tes', 'tis', 'ups', 'nis', 'phis', 'ngs', 'tgs']:
                packageobject('bbb').__setattr__(var, hf['bbb'][var][()]) 
        except:
            raise

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
        try:
            bbb.tipers[...] = np.array(hf.get('tipers@bbb'))
        except ValueError as error:
            print("Couldn't read tipers from  ", file)
            print(error)
        except:
            print("Couldn't read tipers from  ", file)

    hf.close()
    return True


def hdf5_save(file, varlist=['bbb.ngs', 'bbb.ng',
                             'bbb.ni', 'bbb.nis',
                             'bbb.phi', 'bbb.phis',
                             'bbb.te', 'bbb.tes',
                             'bbb.ti', 'bbb.tis',
                             'bbb.up', 'bbb.ups',
                             'bbb.tg', 'bbb.tgs',
                             'bbb.ev', 'bbb.prad',
                             'bbb.pradhyd', 'bbb.tipers','com.nx',
                             'com.ny', 'com.rm', 'com.zm'],
              addvarlist=[]):
    """
    Save HDF5 output for restarting and plotting.
    varlist=[] a list of variables to save specified as strings.
               package prefix required. Default list is usual
               variable list. Example use: 
               varlist=['bbb.ni','bbb.te']
    addvarlist=[] a list of variables to save in addition to the ones
               in varlist. Syntax is the same as varlist parameter.
               Envisioned use is to add output in addition
               to the default list in varlist. 

    """
    grps = {}
    vars = {}
    try:
        hf = h5py.File(file, 'w')
        hfb = hf.create_group('bbb')
        grps['bbb'] = {'h5': hfb}
        grps['bbb']['vars'] = ['uedge_ver']
        grps['bbb']
        hfb.attrs['time'] = time.time()
        hfb.attrs['ctime'] = time.ctime()
        hfb.attrs['code'] = 'UEDGE'
        hfb.attrs['ver'] = bbb.uedge_ver
        try:
            hfb.attrs['pyver'] = pyver
            grps['bbb']['vars'].append('pyver')
        except:
            print("couldn\'t write pyver to header")
    except ValueError as error:
        print("HDF5 file open failed to ", file)
        print(error)
        raise
    except:
        print("HDF5 file open failed to ", file)
        raise

    for lvt in varlist:
        try:
            vt = lvt.split('.')
            if not vt[0] in grps.keys():
                hfb = hf.create_group(vt[0])
                grps[vt[0]] = {'h5': hfb}
                grps[vt[0]]['vars'] = []
            else:
                hfb = grps[vt[0]]['h5']
            pck = packagename2object(vt[0])
            po = pck.getpyobject(vt[1])
            if vt[1] in grps[vt[0]]['vars']:
                print(vt[1], " already written, skipping")
            else:
                grps[vt[0]]['vars'].append(vt[1])
                d = hfb.create_dataset(vt[1], data=po)
                d.attrs['units'] = pck.getvarunit(vt[1])
                d.attrs['comment'] = pck.getvardoc(vt[1])
        except ValueError as error:
            print("HDF5 write failed to ", file, ' var ', vt[1])
            print(error)
        except:
            print("HDF5 write failed to ", file, ' var ', vt[1])

    for lvt in addvarlist:
        try:
            vt = lvt.split('.')
            if not vt[0] in grps.keys():
                hfb = hf.create_group(vt[0])
                grps[vt[0]] = {'h5': hfb}
                grps[vt[0]]['vars'] = []
            else:
                hfb = grps[vt[0]]['h5']
            pck = packagename2object(vt[0])
            po = pck.getpyobject(vt[1])
            if vt[1] in grps[vt[0]]['vars']:
                print(vt[1], " already written, skipping")
            else:
                grps[vt[0]]['vars'].append(vt[1])
                d = hfb.create_dataset(vt[1], data=po)
                d.attrs['units'] = pck.getvarunit(vt[1])
                d.attrs['comment'] = pck.getvardoc(vt[1])
        except ValueError as error:
            print("HDF5 write failed to ", file, ' var ', vt[1])
            print(error)
        except:
            print("HDF5 write failed to ", file, ' var ', vt[1])

    hf.close()
    return True


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
        raise
    except:
        print("Couldn't open hdf5 file ", file)
        raise
    for p in packages:
        hfg = hf.create_group(p.name())
        hfg.attrs['time'] = time.time()
        hfg.attrs['ctime'] = time.ctime()
        hfg.attrs['code'] = 'UEDGE'
        hfg.attrs['ver'] = bbb.uedge_ver
        try:
            hfg.attrs['pyver'] = pyver
        except:
            pass
        for v in list_package_variables(p, vars=vars):
            if p.allocated(v):
                try:
                    d = hfg.create_dataset(v, data=p.getpyobject(v))
                    d.attrs['units'] = p.getvarunit(v)
                    d.attrs['comment'] = p.getvardoc(v)
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
            hfg.attrs['pyver'] = pyver
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
    return True


def h5py_dataset_iterator(g, prefix=''):
    for key in g.keys():
        item = g[key]
        path = '{}/{}'.format(prefix, key)
        if isinstance(item, h5py.Dataset):  # test for dataset
            yield (path, item)
        elif isinstance(item, h5py.Group):  # test for group (go down)
            # following yield is not python 2.7 compatible
            # yield from h5py_dataset_iterator(item, path)
            for (path, item) in h5py_dataset_iterator(item, path):
                yield (path, item)


def hdf5_restore_dump(file, scope=None, hdffile=None,quiet=False):
    """
       Restore all variables from a previously saved HDF5 file.
       This is called by hdf5_restore and the recommended way
       to restore.
    """
    if hdffile == None:
        try:
            hf = h5py.File(file, 'r')
        except:
            print("Couldn't open hdf5 file ", file)
            raise
    else:
        hf = hdffile
    try:
        try:
            g = hf['bbb']
            if not quiet:
                prfileattrs = False
                print('File attributes:')
                print('     written on: ', g.attrs['ctime'])
                print('        by code: ', g.attrs['code'])
                print('       version: ', np.char.strip(g.attrs['ver']))
                print('    physics tag: ', np.char.strip(g.attrs['ver']))
                print(' python version: ', np.char.strip(g.attrs['pyver']))
        except:
            if not quiet:
                print('No file attributes, trying to restore')

        for (path, dset) in h5py_dataset_iterator(hf):
            vt = path.split('/')
            if vt[1] == 'globals':
                if scope != None:
                    if dset.size > 1:
                        scope[vt[2]] = np.array(dset[()])
                    else:
                        scope[vt[2]] = dset[()]

            else:
                try:
                    pck = packagename2object(vt[1])
                    po = pck.getpyobject(vt[2])
                    if dset.size > 1:
                        po[...] = np.array(dset[()])
                    else:
                        setattr(pck, vt[2], dset[()])
                except ValueError as error:
                    print("Couldn't read tes from  ", file)
                    print(error)
                except:
                    print('Couldn\'t read dataset ', path)
                    raise
    except:
        print("Couldn't read hdf5 file ", file)
        raise

    if hdffile == None:
        hf.close()

    return True
