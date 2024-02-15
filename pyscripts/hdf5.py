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
    from os.path import exists
    from numpy import array
    from h5py import File 

    iprint = packageobject("com").__getattribute__("iprint")
    # Check that the requested file exists
    if not exists(file):
        raise Exception("{} not found: Aborting!".format(file))

    # The standard restore variables
    varlist = ['tes', 'tis', 'ups', 'nis', 'phis', 'ngs']
    # Newer versions of UEDGE introduce additional, newer variables,
    # which are to be restored if available:
    # - tgs: gaseous temperture for solving gaseous energy eqs
    # - tipers: perpendicular ion temperature
    for newvar in ['tgs', 'tipers']:
        # Ensure any available new variables also are restored
        if newvar in packageobject("bbb").varlist():
            varlist.append(newvar)

    with File(file, "r") as hf:
        # Check whether an (very) old version of the save file is used
        if "nis@bbb" in list(hf.keys()):
            if iprint > 0:
                print("\nReading old style, flat-structured, save-file '{}'.".format(file))
            for var in varlist:
                # Restore variable if available from save
                try:
                    packageobject('bbb').__setattr__(
                        var, 
                        array(hf.get('{}@bbb'.format(var)))
                    )
                except ValueError as error:
                    print("Couldn't read {} from {}".format(var, file))
                    print(error)
                except:
                    print("Couldn't read {} from {}".format(var, file))
            return True
        # Check whether file is UETOOLS-style file
        if "restore" in list(hf.keys()):
            if iprint > 0:
                print("\nReading UETOOLS style save-file '{}'.".format(file))
            restoregroup = hf['restore/bbb']
        # Check whether file is uedge.hdf5-style file
        elif "bbb" in list(hf.keys()):
            if iprint > 0:
                print("\nReading uedge.hdf5-style save-file '{}'.".format(file))
            restoregroup = hf['bbb']
        # Unknown file structure, abort
        else:
            raise Exception("\nUnable to read data from {}".file(file))
        # Read all the required restore variables
        for var in varlist:
            try:           
                packageobject('bbb').__setattr__(var, restoregroup[var][()]) 
            except:
                raise
    if iprint > 0:
        print(">>> Save read successfully")
    return True


def hdf5_save(file, varlist=['bbb.ngs', 'bbb.ng',
                             'bbb.ni', 'bbb.nis',
                             'bbb.phi', 'bbb.phis',
                             'bbb.te', 'bbb.tes',
                             'bbb.ti', 'bbb.tis',
                             'bbb.up', 'bbb.ups',
                             'bbb.tg', 'bbb.tgs',
                             'bbb.ev', 'bbb.prad',
                             'bbb.pradhyd', 'com.nx',
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
    try:
        bbb.tipers
        varlist.append('bbb.tipers')
    except:
        pass

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
    from numpy import char, array
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
                print('       version: ', char.strip(g.attrs['ver']))
                print('    physics tag: ', char.strip(g.attrs['ver']))
                print(' python version: ', char.strip(g.attrs['pyver']))
        except:
            if not quiet:
                print('No file attributes, trying to restore')

        for (path, dset) in h5py_dataset_iterator(hf):
            vt = path.split('/')
            if vt[1] == 'globals':
                if scope != None:
                    if dset.size > 1:
                        scope[vt[2]] = array(dset[()])
                    else:
                        scope[vt[2]] = dset[()]

            else:
                try:
                    pck = packagename2object(vt[1])
                    po = pck.getpyobject(vt[2])
                    if dset.size > 1:
                        po[...] = array(dset[()])
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
