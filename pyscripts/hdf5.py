
import numpy as np
import h5py
from uedge import bbb
from uedge import com



def hdf5_restore(file):
    """
        Read a hdf5 file previously written from Uedge. This reads the file and puts
        the 6 standard variables into the correct format.
    """
    try:
        hf = h5py.File(file,'r')
    except:
        print "Couldn't open hdf5 file ",file
        return


    try:
       hfg = hf.setdefault('bbb')
       print "New style hdf5 file"
       try:
           bbb.ngs[...]  = np.array(hfg.get('ngs'))
       except:
           print "Couldn't read ngs from  ",file
       try:
           bbb.nis[...]  = np.array(hfg.get('nis'))
       except:
           print "Couldn't read nis from  ",file
       try:
           bbb.phis[...]  = np.array(hfg.get('phis'))
       except:
           print "Couldn't read phis from  ",file
       try:
           bbb.tes[...]  = np.array(hfg.get('tes'))
       except:
           print "Couldn't read tes from  ",file
       try:
           bbb.tis[...]  = np.array(hfg.get('tis'))
       except:
           print "Couldn't read tis from  ",file
       try:
           bbb.ups[...]  = np.array(hfg.get('ups'))
       except:
           print "Couldn't read ups from  ",file
       try:
           bbb.tgs[...]  = np.array(hfg.get('tgs'))
       except:
           print "Couldn't read tgs from  ",file

    except:
       print "Old style hdf5 file"
       try:
           bbb.ngs[...]  = np.array(hf.get('ngs@bbb'))
       except:
           print "Couldn't read ngs from  ",file
       try:
           bbb.nis[...]  = np.array(hf.get('nis@bbb'))
       except:
           print "Couldn't read nis from  ",file
       try:
           bbb.phis[...]  = np.array(hf.get('phis@bbb'))
       except:
           print "Couldn't read phis from  ",file
       try:
           bbb.tes[...]  = np.array(hf.get('tes@bbb'))
       except:
           print "Couldn't read tes from  ",file
       try:
           bbb.tis[...]  = np.array(hf.get('tis@bbb'))
       except:
           print "Couldn't read tis from  ",file
       try:
           bbb.ups[...]  = np.array(hf.get('ups@bbb'))
       except:
           print "Couldn't read ups from  ",file
       try:
           bbb.tgs[...]  = np.array(hf.get('tgs@bbb'))
       except:
           print "Couldn't read tgs from  ",file

    hf.close()




# --------------------------------------------------------------------------------

def hdf5_save(file):
    """
    Save HDF5 output for restarting and plotting.
    """
    
    try:
        hf = h5py.File(file,'w')
        hfb = hf.create_group('bbb')
    except:
        print "HDF5 file open failed to ", file
    try:
        hfb.create_dataset('ng',data=bbb.ng)
        hfb.create_dataset('ngs',data=bbb.ngs)
    except:
        print "ng or ngs HDF5 write failed to ", file
    try:
        hfb.create_dataset('ni',data=bbb.ni)
        hfb.create_dataset('nis',data=bbb.nis)
    except:
        print "ni or nis HDF5 write failed to ", file
    try:
        hfb.create_dataset('phi',data=bbb.phi)
        hfb.create_dataset('phis',data=bbb.phis)
    except:
        print "phi or phis HDF5 write failed to ", file
    try:
        hfb.create_dataset('te',data=bbb.te)
        hfb.create_dataset('tes',data=bbb.tes)
    except:
        print "te or tes HDF5 write failed to ", file
    try:
        hfb.create_dataset('ti',data=bbb.ti)
        hfb.create_dataset('tis',data=bbb.tis)
    except:
        print "ti or tis HDF5 write failed to ", file
    try:
        hfb.create_dataset('up',data=bbb.up)
        hfb.create_dataset('ups',data=bbb.ups)
    except:
        print "up or ups HDF5 write failed to ", file
    try:
        hfb.create_dataset('tg',data=bbb.tg)
        hfb.create_dataset('tgs',data=bbb.tgs)
    except:
        print "tg or tgs HDF5 write failed to ", file

    try:
        hfb.create_dataset('ev',data=bbb.ev)
    except:
        print "bbb HDF5 write failed to ", file

    try:
        hfc = hf.create_group('com')
        hfc.create_dataset('nx',data=com.nx)
        hfc.create_dataset('ny',data=com.ny)
        hfc.create_dataset('rm',data=com.rm)
        hfc.create_dataset('zm',data=com.zm)
    except:
        print "com HDF5 write failed to ", file        

    hf.close()


