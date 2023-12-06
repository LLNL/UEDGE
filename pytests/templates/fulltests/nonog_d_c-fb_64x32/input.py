# 
###########################################################################
# DESCRIPTION OF PROBLEM (d3dHmsCnog) from FACETS test suite:
# DIII-D single-null geometry with 6 hydrogen variables (ni,ng,upi,upg,te,ti)
# and 7 carbon variables (six charge-state densities ni and one ng).
# Solver used is Newton Krylov (svrpkg="nksol") and preconditioner uses an 
# iterative solver ILUT for Jacobian LU decomposition. Also includes tilted
# divertor plates wrt flux-surface normal, thus testing the nonorthogonal
# finite-volume difference stencil.  Iterates to steady-state solution from 
# an initial profile file (HDF5).
###########################################################################

# Import uedge into python and make variables active
from uedge import *
from uedge.hdf5 import *

# Set seec switch to new model, if available
try:
    bbb.oldseec = 0
except:
    pass

# Read grid from file
bbb.mhdgeo = 1                  #=1 use MHD equilibrium files
bbb.gengrid = 0         # reads mesh from gridue file
com.geometry = "snull"  # Defines a LSN geometry for the problem
bbb.GridFileName = "gridue_64x32" # Defines the gridue file name
com.isnonog = 1         # non orthogonal differencing

# Finite-difference algorithms (upwind, central diff, etc.)
bbb.methn = 33          #ion continuty eqn
bbb.methu = 33          #ion parallel momentum eqn
bbb.methe = 33          #electron energy eqn
bbb.methi = 33          #ion energy eqn
bbb.methg = 66

# Boundary conditions
bbb.ncore[0] = 2.e19    #hydrogen ion density on core
bbb.iflcore = 0         #flag; =0, fixed Te,i; =1, fixed power on core
bbb.tcoree = 250.       #core Te
bbb.tcorei = 250.       #core Ti
bbb.tedge = 2.          #fixed wall,pf Te,i if istewcon=1, etc
bbb.recycp[0] = 0.95    #hydrogen recycling coeff at plates
bbb.recycm = 0.1

# Transport coefficients (m**2/s)
bbb.difni[0] = 1.       #D for radial hydrogen diffusion
bbb.kye = 1.            #chi_e for radial elec energy diffusion
bbb.kyi = 1.            #chi_i for radial ion energy diffusion
bbb.travis[0] = 1.      #eta_a for radial ion momentum diffusion

# Flux limits
bbb.flalfe = 0.21       #electron parallel thermal conduct. coeff
bbb.flalfi = 0.21       #ion parallel thermal conduct. coeff
bbb.flalfv = 0.5        #ion parallel viscosity coeff
bbb.flalfgx = 1.        #neut. gas part. flux in poloidal direction
bbb.flalfgy = 1.        #neut. gas part. flux in radial direction
bbb.flalfgxy = 1.       #neut. gas part. flux in  mixed derivatives
bbb.flalftgx = 1.       #neut. gas thermal flux, poloidal direction
bbb.flalftgy = 1.       #neut. gas thermal flux, radial direction
bbb.lgmax = 0.1		#max scale length for flalfgx,y
bbb.lgtmax = 0.1	#max scale length for flalftgx,y
bbb.lgvmax = 0.1	#max scale length for flalfvgx,y


# Solver package
bbb.svrpkg = "nksol"    #Newton solver using Krylov method
bbb.premeth = "ilut"    #Solution method for precond. Jacobian matrix

# Parallel neutral momentum equation
bbb.isupgon[0] = 1
bbb.ineudif = 2             #=2 for evolving pg=ng*tg variable
bbb.isngon[0] = 0
com.ngsp = 1
com.nhsp = 2

bbb.ziin[1] = 0

## Impurity gas basics
com.ngsp = 2                #total number of gas species
bbb.isngon[1] = 1           #turns on impurity gas
bbb.ngbackg[1] = 1.e9	    #neutral impurity background for added source
bbb.ingb = 2		    #exponent for strength of ngbackg turn-on
bbb.istgcon[1] = 1	    #=1 for constant tg(2) at tgas(2)
bbb.tgas[1] = 1.	    #value for tg when istgcon=1
bbb.rcxighg = 0.            # best value; ratio of imp cx to hyd cx
bbb.kelighi[1] = 5.e-16     #elastic sig_v for imp_gas/h_ion
bbb.kelighg[1] = 5.e-16     #elastic sig_v for imp_gas/h_gas
bbb.n0g[1] = 1.e16          #imp. gas density normalization

# Impurity gas boundary conditions
bbb.recycp[1] = 0.01        #plate recycling of impurities
bbb.recycw[1] = 1e-4        #wall recycling; matwsi,o set above for hyd
bbb.isch_sput[1]=7          # Haasz/Davis chemical sputtering model
bbb.isph_sput[1]=3 	        # physical sputtering model

# Reduce C isputtering owing to boronization
bbb.fphysylb[1,0] = 0.5
bbb.fphysyrb[1,0] = 0.5
bbb.fchemylb[1,0] = 0.1
bbb.fchemyrb[1,0] = 0.1

bbb.t_wall = 300.
bbb.t_plat = 500.

## Impurity ions
bbb.isimpon = 6             #Use force-balance only
com.nzsp[0] = 6             #number chrg states impurity isotope #1

bbb.csfaclb[2:8,0] = 2.191
bbb.csfacrb[2:8,0] = 2.191
bbb.minu[2:8] = 12.
bbb.ziin[:2] = [1, 0]
bbb.ziin[2:8] = list(range(1,7))
bbb.znuclin[:2] = 1
bbb.znuclin[2:8] = 6
bbb.n0[2:8] = 1.e17

bbb.nzbackg = 1.e9          #background density for impurities
bbb.inzb = 2		    #exponent for switching on nzbackg
bbb.ismctab = 2             # use Braams" rate tables
com.mcfilename[0] = "C_rates.strahl"     # Imp rate file name

bbb.isnicore[7] = 3
bbb.curcore[7] = 0.

bbb.isnwcono[2:8] = 3
bbb.isnwconi[2:8] = 3
bbb.nwomin[2:8] = 1.e7
bbb.nwimin[2:8] = 1.e7

bbb.restart = 1             #Begin from savefile, not estimated profiles

# Filling newly allocated arrays as desired
bbb.ftol = 1.e-8

# Atomic data switches
com.istabon = 10            #=10 specifics hydrogen data file ehr2.dat

# Scale factor converting (upi-upg)**2 energy to thermal energy
bbb.cfnidh = 0.2

bbb.allocate()      #allocates storage for arrays
hdf5_restore("solution.h5")
