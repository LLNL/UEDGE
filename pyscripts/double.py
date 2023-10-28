from uedge import *

def uedouble():
    com.nxleg=2*com.nxleg
    com.nxcore=2*com.nxcore
    com.nycore=2*com.nycore
    com.nysol=2*com.nysol
    if com.nxomit > 0:
        if com.geometry=="dnbot":
            com.nxomit = com.nxleg[0,0]+com.nxcore[0,0] + 2*com.nxxptxx + 1
        else:
            com.nxomit=2*(com.nxomit-2*com.nxxptxx) + 2*com.nxxpt  # assumes com.nxomit removes 1/2 SOL
    if com.nyomitmx == 1: com.nysol = 1
    if grd.kxmesh == 4:
        dxgasold=grd.dxgas
        alfxold=grd.alfx
        grd.alfx=alfxold/2.
        grd.dxgas=dxgasold*(exp(grd.alfx)-1)/(exp(alfxold)-1)
        grd.nxgas=2*grd.nxgas
    bbb.restart=1
    bbb.newgeo=1
    bbb.gengrid=1
    bbb.isnintp=1
    grd.ixdstar = com.nxcore[0,1]+1
