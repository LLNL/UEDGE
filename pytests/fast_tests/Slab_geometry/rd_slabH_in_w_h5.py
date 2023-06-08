#
# 
###########################################################################
# DESCRIPTION OF PROBLEM (slabH):
# Slab model with 4 hydrogen variables (ni,upi,te,ti) and a (6+2)*(10+2) 
# = 8x12 [poloidal*radial] mesh yielding 384 variables.  The +2 in the mesh 
# size description arises from one guard-cell at each end of the domain used 
# to set boundary conditions. This case starts from generic internal initial
# profiles (bbb.restart=0), which is generally the "hard way" to start if a
# similar solution has been saved previous as an HDF5 (or PDB) file; see
# ../tokgeo_H case. Also, this case has the neutrals frozen (bbb.isngon=0)
# and uses a simple internal hydrogen ionization function for any residual
# ionization from the frozen neutrals (small).
# Solver is the time-dependent, Newton-Krylov ODE routine VODPK 
# (svrpkg="vodpk") and the preconditioner approximate LU decomposition is
# provided by the ILUT sparse solver (premeth="ilut").
###########################################################################
import uedge
from uedge import *

bbb.mhdgeo=-1
bbb.isnion=1
bbb.isupon=1
bbb.isteon=1
bbb.istion=1
bbb.isngon=0
bbb.svrpkg="nksol"
bbb.premeth="ilut"
bbb.epscon1=1.e-2
bbb.ireorder=0
bbb.ncore=2.e19
bbb.tcoree=100
bbb.tcorei=100
bbb.isfixlb=2
bbb.recycp=.9
bbb.trange=4.e3
bbb.nsteps=1
bbb.n0g=1.e16
bbb.difni=1.
bbb.kye=1.
bbb.flalfe=0.21
bbb.flalfi=0.21
bbb.flalfgx=1.e10
bbb.flalfgy=1.e10
com.nycore=0
com.nysol=10
com.nxleg[0,0]=0
com.nxleg[0,1]=2
com.nxcore[0,0]=0
com.nxcore[0,1]=4
grd.zax=1.
grd.zaxpt=.75
grd.alfyt=-1.e-5
bbb.restart=0
bbb.ftol = 1.e-5
bbb.dtreal = 1.e2

bbb.allocate()      #allocates storage for arrays
bbb.restart = 1     #use initial solution in slabH.h5
from uedge.hdf5 import *
hdf5_restore("slabH.h5")

# perform on time-step
bbb.exmain()

#now should have starting profile at dtreal=1.e-9; restart from that soln
bbb.restart = 1

# Increase bbb.dtreal for more time-steps or read bbb.rdinitdt and bbb.rdcontdt

