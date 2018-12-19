
import numpy as np
import netCDF4 as nc
from uedge import bbb



def cdf4_restore(file):
    """
        Read a cdf4 file previously written from Uedge. This reads the file and puts
        the 6 standard variables into the correct format.
    """
    try:
        cf = nc.Dataset(file)
    except:
        print "Couldn't open cdf4 file ",file
        return


    try:
        bbb.ngs[...]  = np.array(cf.variables['ngs'])
    except:
        print "Couldn't read ngs from  ",file
    try:
        bbb.nis[...]  = np.array(cf.variables['nis'])
    except:
        print "Couldn't read nis from  ",file
    try:
        bbb.phis[...]  = np.array(cf.variables['phis'])
    except:
        print "Couldn't read phis from  ",file
    try:
        bbb.tes[...]  = np.array(cf.variables['tes'])
    except:
        print "Couldn't read tes from  ",file
    try:
     bbb.tis[...]  = np.array(cf.variables['tis'])
    except:
        print "Couldn't read tis from  ",file
    try:
        bbb.ups[...]  = np.array(cf.variables['ups'])
    except:
        print "Couldn't read ups from  ",file
    try:
        bbb.tgs[...]  = np.array(cf.variables['tgs'])
    except:
        print "Couldn't read tgs from  ",file

    cf.close()






def cdf4_save(file):
    """
        Write the 6 standard variables into an cdf4 file.
    """
    try:
        hf = h5py.File(file,'w')
        hfg = hf.create_group('bbb')
    except:
        print "Couldn't open cdf4 file ",file
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


