
import numpy as np
import h5py
from uedge import bbb



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






def hdf5_save(file):
    """
        Write the 6 standard variables into an hdf5 file.
    """
    try:
        hf = h5py.File(file,'w')
        hfg = hf.create_group('bbb')
    except:
        print "Couldn't open hdf5 file ",file
    try:
        hfg.create_dataset('ngs',data=bbb.ngs)
    except:
        print "Couldn't write ngs to  ",file
    try:
        hfg.create_dataset('nis',data=bbb.nis)
    except:
        print "Couldn't write nis to  ",file
    try:
        hfg.create_dataset('phis',data=bbb.phis)
    except:
        print "Couldn't write phis to  ",file
    try:
        hfg.create_dataset('tes',data=bbb.tes)
    except:
        print "Couldn't write tes to  ",file
    try:
        hfg.create_dataset('tis',data=bbb.tis)
    except:
        print "Couldn't write tis to  ",file
    try:
        hfg.create_dataset('ups',data=bbb.ups)
    except:
        print "Couldn't write ups to  ",file
    try:
        hfg.create_dataset('tgs',data=bbb.tgs)
    except:
        print "Couldn't write ups to  ",file

    hf.close()


