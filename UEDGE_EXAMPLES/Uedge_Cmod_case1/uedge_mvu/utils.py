'''
Methods of very general utility.
'''
import sys
#from distutils.util import strtobool
import numpy as np
#import shapely.geometry
from uedge import bbb, com, flx, grd, svr

def paws(message=""):
    print(message)
    input("Press the <ENTER> key to continue...")


def uelist(varname):
    """
    Show info on a variable
    """
    resbbb=bbb.listvar(varname)
    if (resbbb != None):
        print(resbbb)

    rescom=com.listvar(varname)
    if (rescom != None):
        print(rescom)

    ressvr=svr.listvar(varname)
    if (ressvr != None):
        print(ressvr)

    resflx=flx.listvar(varname)
    if (resflx != None):
        print(resflx)

    resgrd=grd.listvar(varname)
    if (resgrd != None):
        print(resgrd)
