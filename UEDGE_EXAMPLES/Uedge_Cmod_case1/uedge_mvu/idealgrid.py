from fortranformat import FortranRecordWriter
import re
import numpy as np
from uedge import bbb, com, grd


def import_field(fp,nxm,nym):

    nvals=(nxm+2)*(nym+2)*5
    buffer=np.zeros(nvals)

    #-how many lines to read
    nlines=int(np.ceil(float(nvals)/3.))

    #-skip a blank line
    line=fp.readline()

    #-read data
    for i in range(0,nlines):
        line=fp.readline()
        columns=line.split()
        str=columns[0]; buffer[i*3]=float(re.sub('[dD]', 'e', str))
        str=columns[1]; buffer[i*3+1]=float(re.sub('[dD]', 'e', str))
        str=columns[2]; buffer[i*3+2]=float(re.sub('[dD]', 'e', str))

    #-convert data to 3D array
    fdata=np.zeros((nxm+2,nym+2,5))

    i=0
    for k in range(5):
        for jy in range(0,nym+2):
            for ix in range(0,nxm+2):
                fdata[ix,jy,k]=buffer[i]
                i=i+1

    return fdata


def export_field(fp, nxm, nym, fdata):

    #-skip a line                                                                                             
    fp.write("\n")

    nvals=(nxm+2)*(nym+2)*5
    buffer=np.zeros(nvals)
    nlines=int(np.ceil(float(nvals)/3.))


    i=0
    for k in range(5):
        for jy in range(0,nym+2):
            for ix in range(0,nxm+2):
                buffer[i]=fdata[ix,jy,k]
                i=i+1

    frw = FortranRecordWriter('(3e23.15)')


    for m in range(nlines):
        line = frw.write(buffer[3*m:3*m+3])+"\n"
        fp.write(line)




def idlg_Read(fname="gridue"):

    print("In idealgrid: idlg_read()")
    print("Importing data from ", fname)

    f = open(fname, 'r')
     
    #-read headerline
    line=f.readline()
    columns=line.split()
    nxm=int(columns[0])
    nym=int(columns[1])
    ixpt1=int(columns[2])
    ixpt2=int(columns[3])
    iysptrx1=int(columns[4])
    

    #import data, one field at a time
    rm=import_field(f,nxm,nym)
    zm=import_field(f,nxm,nym)
    psi=import_field(f,nxm,nym)
    br=import_field(f,nxm,nym)
    bz=import_field(f,nxm,nym)
    bpol=import_field(f,nxm,nym)
    bphi=import_field(f,nxm,nym)
    b=import_field(f,nxm,nym)


    #import runid
    line=f.readline()
    print("runid:",line)

    f.close()

    return rm,zm,psi,br,bz,bpol,bphi,b,nxm,nym,ixpt1,ixpt2,iysptrx1



def idlg_Write(fname="gridue", runid="runid"):
    ##-saving UEDGE geometry and magnetic data in gridue file-##

    print("In idealgrid: idlg_write()")
    print("Exporting data to ", fname)


    f = open(fname, 'w')

    frw = FortranRecordWriter('(I4,I4,I4,I4,I4)')
    headerline = frw.write([com.nxm,com.nym,com.ixpt1[0],com.ixpt2[0],com.iysptrx1[0]])+"\n"

    f.write(headerline)

    export_field(f, com.nxm, com.nym, com.rm)
    export_field(f, com.nxm, com.nym, com.zm)
    export_field(f, com.nxm, com.nym, com.psi)
    export_field(f, com.nxm, com.nym, com.br)
    export_field(f, com.nxm, com.nym, com.bz)
    export_field(f, com.nxm, com.nym, com.bpol)
    export_field(f, com.nxm, com.nym, com.bphi)
    export_field(f, com.nxm, com.nym, com.b)

    f.write("\n"+runid)

    f.close()



def idlg_ModRad(jymin=0, jymax=1, alfyt=0.0):
    ##-modify radial grid in a selected jy range-##

    print("Modifying radial grid...")

    #-range of affected radial cells (need at least two cells)
    ncell=jymax-jymin+1

    #-radial range for selected jy
    rmin=com.rm[0,jymin,1]
    rmax=com.rm[0,jymax,3]
    ##deltar=rmax-rmin

    #-calculate radial locations of vertices
    if (np.abs(alfyt) <= 1e-3):
        #-uniform
        rr=rmin + (rmax-rmin)*(np.arange(ncell+1))/(ncell)
    else:
        #-exponential
        rr=rmin + (rmax-rmin)*(np.exp(alfyt*(np.arange(ncell+1)))-1.0)/(np.exp(alfyt*(ncell))-1.0)


    for ix in range(com.nxm+2):
        com.rm[ix,jymin:jymax+1,1]=rr[0:ncell]
        com.rm[ix,jymin:jymax+1,2]=rr[0:ncell]
        com.rm[ix,jymin:jymax+1,3]=rr[1:ncell+1]
        com.rm[ix,jymin:jymax+1,4]=rr[1:ncell+1]
        com.rm[ix,jymin:jymax+1,0]=0.5*(rr[0:ncell]+rr[1:ncell+1])

    #-update B field on the grid
    idlg_SetB()


def idlg_ModPol(ixmin=0, ixmax=0, alfxt=0.0, smoothe=False, debug=False):
    ##-modify poloidal grid in selected ix range [ixmin,ixmax]-##

    print("Modifying poloidal grid...")

    if (np.abs(alfxt)>1.99):
        print("must have abs(alfxt)<2!")

    #-range of affected poloidal cells
    ncell=ixmax-ixmin+1

    #-range of poloidal coordinate (implies ixmax>=ixmin?)
    zmin=com.zm[ixmin,0,1]
    zmax=com.zm[ixmax,0,2]

    if debug:
        print("Poloidal index range", ixmin, ixmax)
        print("Poloidal coordinate range", zmin, zmax)


    if (np.abs(alfxt)<1e-3):
        #-uniform spacing
        zz=zmin + (zmax-zmin)*(np.arange(ncell+1))/(ncell)
    else:
        zz=zmin + 0.5*(zmax-zmin)*(1.0+np.arctanh(alfxt*((1.0*np.arange(ncell+1)/ncell)-0.5))/np.arctanh(alfxt*0.5))


    if (smoothe):
        zz_lin = zmin + (zmax-zmin)*(np.arange(ncell+1))/(ncell)
        zz[1]=0.5*(zz[1]+zz_lin[1])
        zz[ncell-1]=0.5*(zz[ncell-1]+zz_lin[ncell-1])

    if debug:
        print("Poloidal vertices: ", zz)    


    for jy in range(com.nym+2):
        com.zm[ixmin:ixmax+1,jy,1]=zz[0:ncell]
        com.zm[ixmin:ixmax+1,jy,3]=zz[0:ncell]
        com.zm[ixmin:ixmax+1,jy,2]=zz[1:ncell+1]
        com.zm[ixmin:ixmax+1,jy,4]=zz[1:ncell+1]
        com.zm[ixmin:ixmax+1,jy,0]=0.5*(zz[0:ncell]+zz[1:ncell+1])

    #-update B field on the grid
    idlg_SetB()




def idlg_SetB():
    #-fix the magnetic field for modified grid

    print("Setting magnetic field on the grid")

    btorfix=np.sqrt(grd.btfix**2-grd.bpolfix**2)

    for k in range(5):
        for ix in range(com.nxm+2):
            for jy in range(com.nym+2):                
                com.bphi[ix,jy,k] = btorfix*(com.rm[ix,jy,k]/grd.rmajfix)**grd.sigma_btor
                com.bpol[ix,jy,k] = grd.bpolfix*(com.rm[ix,jy,k]/grd.rmajfix)**grd.sigma_bpol
                com.b[ix,jy,k]    = np.sqrt(com.bphi[ix,jy,k]**2 + com.bpol[ix,jy,k]**2)

                ##-Note: poloidal flux per radian, correct only for cylindrical case
                com.psi[ix,jy,k]  = ((grd.bpolfix*grd.rmajfix**2)/(grd.sigma_bpol+2))*(com.rm[ix,jy,k]/grd.rmajfix)**(grd.sigma_bpol+2)

                com.br[ix,jy,k] = 0.
                com.bz[ix,jy,k] = - com.bpol[ix,jy,k]

