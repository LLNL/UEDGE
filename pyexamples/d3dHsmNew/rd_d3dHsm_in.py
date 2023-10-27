#
#
###########################################################################
# DESCRIPTION OF PROBLEM (d3dHsm) from FACETS test suite:
# DIII-D single-null geometry with 5 variables (ni,upi,te,ti,ng) and a
# (16+2)*(8+2)=18x10 [poloidal*radial] mesh yielding 900 variables.
# Solver used is Newton Krylov (svrpkg="nksol") and preconditioner uses a
# direct banded solver for the LU decomposition (premeth="banded"). Iterates
# to steady-state solution from an initial profile file (HF5).
###########################################################################
###import uedge

from uedge import *
from uedge.hdf5 import *
from uedge.rundt import *


# Set the geometry
bbb.mhdgeo = 1 		        #=1 use MHD equilibrium files
#flx.aeqdskfname = "aeqdskd3d"   #name of EFIT 'a' file for flux-surface mesh
#flx.geqdskfname = "neqdskd3d"   #name of EFIT 'g' or 'n' file for flux-sur mesh
flx.psi0min1 = 0.98		#normalized flux on core bndry
flx.psi0min2 = 0.98		#normalized flux on pf bndry
flx.psi0sep = 1.00001	        #normalized flux at separatrix
flx.psi0max = 1.07		#normalized flux on outer wall bndry
bbb.ngrid = 1		        #number of mesh sequenc. (always set to 1)
com.nxleg[0,0] = 4		#pol. mesh pts from inner plate to x-point
com.nxcore[0,0] = 4		#pol. mesh pts from x-point to top on inside
com.nxcore[0,1] = 4		#pol. mesh pts from top to x-point on outside
com.nxleg[0,1] = 4		#pol. mesh pts from x-point to outer plate
com.nysol[0] = 6		#rad. mesh pts in SOL
com.nycore[0] = 2		#rad. mesh pts in core

# Finite-difference algorithms (upwind, central diff, etc.)
bbb.methn = 33		#ion continuty eqn
bbb.methu = 33		#ion parallel momentum eqn
bbb.methe = 33		#electron energy eqn
bbb.methi = 33		#ion energy eqn
bbb.methg = 33		#neutral gas continuity eqn

# Boundary conditions
bbb.ncore[0] = 2.5e19	#hydrogen ion density on core
##	iflcore = 0	 #flag; =0, fixed Te,i; =1, fixed power on core
bbb.tcoree = 100.	#core Te
bbb.tcorei = 100.	#core Ti
bbb.tedge = 2.		#fixed wall,pf Te,i if istewcon=1, etc
bbb.recycp[0] = 0.8	#hydrogen recycling coeff at plates

# Transport coefficients (m**2/s)
bbb.difni[0] = 1.	#D for radial hydrogen diffusion
bbb.kye = 1.		#chi_e for radial elec energy diffusion
bbb.kyi = 1.		#chi_i for radial ion energy diffusion
bbb.travis[0] = 1.	#eta_a for radial ion momentum diffusion

# Flux limits
bbb.flalfe = 0.21       #electron parallel thermal conduct. coeff
bbb.flalfi = 0.21	#ion parallel thermal conduct. coeff
bbb.flalfv = 1.		#ion parallel viscosity coeff
bbb.flalfgx = 1.e20	#neut. gas in poloidal direction
bbb.flalfgy = 1.e20	#neut. gas in radial direction

# Solver package
bbb.svrpkg = "nksol"	#Newton solver using Krylov method
bbb.premeth = "banded"	#Solution method for precond. Jacobian matrix

# Restart from a HDF5 or PDB savefile
bbb.restart = 1	    #Begin from savefile, not estimated profiles
bbb.allocate()      #allocates storage for arrays


# Atomic data switches
com.istabon = 0            #-analytic rates
###com.istabon = 10        #=10 hydrogen data file ehr2.dat



if (1):
    hdf5_restore("d3dHsm.h5")
    bbb.dtreal = 1e20; bbb.exmain()
else:
    #-set up some simple initial state
    bbb.ngs=1e14; bbb.ng=1e14
    bbb.nis=1e20; bbb.ni=1e20 
    bbb.ups=0.0;  bbb.up=0.0
    bbb.tes=bbb.ev;   bbb.te=bbb.ev
    bbb.tis=bbb.ev;   bbb.ti=bbb.ev

    #-short run to initialize everything
    bbb.dtreal = 1e-12; bbb.isbcwdt=1; bbb.exmain()

    ##-run to steady state
    rundt(dtreal=1e-10)
