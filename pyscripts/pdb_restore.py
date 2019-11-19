
import numpy as np
import pact.pdb as pdb
from .uedge import bbb



def pdb_restore(file):
    """
        Read a pdb file previously written from Uedge. This reads the file and puts
        the 6 standard variables into the correct format.
    """
    try:
        fp = pdb.open(file,'r')
    except:
        print("Couldn't open pdb file ",file)
    try:
        bbb.ngs[...]  = np.array(fp.read('ngs@bbb')).reshape(bbb.ngs.shape[::-1]).transpose()
    except:
        print("Couldn't read ngs from  ",file)
    try:
        bbb.nis[...]  = np.array(fp.read('nis@bbb')).reshape(bbb.nis.shape[::-1]).transpose()
    except:
        print("Couldn't read nis from  ",file)
    try:
        bbb.phis[...]  = np.array(fp.read('phis@bbb')).reshape(bbb.phis.shape[::-1]).transpose()
    except:
        print("Couldn't read phis from  ",file)
    try:
        bbb.tes[...]  = np.array(fp.read('tes@bbb')).reshape(bbb.tes.shape[::-1]).transpose()
    except:
        print("Couldn't read tes from  ",file)
    try:
        bbb.tis[...]  = np.array(fp.read('tis@bbb')).reshape(bbb.tis.shape[::-1]).transpose()
    except:
        print("Couldn't read tis from  ",file)
    try:
        bbb.ups[...]  = np.array(fp.read('ups@bbb')).reshape(bbb.ups.shape[::-1]).transpose()
    except:
        print("Couldn't read ups from  ",file)
    try:
        bbb.tgs[...]  = np.array(fp.read('tgs@bbb')).reshape(bbb.tgs.shape[::-1]).transpose()
    except:
        print("Couldn't read tgs from  ",file)

    fp.close()


