#
# 
###########################################################################
# DESCRIPTION OF PROBLEM (d3dHmsCnog) using Strahl imp data; from FACETS :
# DIII-D single-null geometry with 6 hydrogen variables (ni,ng,upi,upg,te,ti)
# and 7 carbon variables (six charge-state densities ni and one ng) on a
# (16+2)*(8+2)=18x10 [poloidal*radial] mesh yielding 2340 variables.
# Solver used is Newton Krylov (svrpkg="nksol") and preconditioner uses an 
# iterative solver ILUT for Jacobian LU decomposition. Also includes tilted
# divertor plates wrt flux-surface normal, thus testing the nonorthogonal
# finite-volume difference stencil.  Iterates to steady-state solution from 
# an initial profile file (HF5).
###########################################################################

# Import uedge into python and make variables active
import uedge
from uedge import *

# Begin uedge parameter input
# Set the geometry
bbb.mhdgeo = 1                  #=1 use MHD equilibrium files
##flx.aeqdskfname = "a110465.03500" #EFIT "a" file for flux-surface mesh
##flx.geqdskfname = "g110465.03500" #EFIT "g" or "n" file for flux-sur mesh
flx.psi0min1 = 0.96             #normalized flux on core bndry
flx.psi0min2 = 0.98             #normalized flux on pf bndry
flx.psi0sep = 1.00001           #normalized flux at separatrix
flx.psi0max = 1.07              #normalized flux on outer wall bndry
bbb.ngrid = 1                   #number of mesh sequenc. (always set to 1)
com.nxleg[0,0] = 4              #pol. mesh pts from inner plate to x-point
com.nxcore[0,0] = 4             #pol. mesh pts from x-point to top on inside
com.nxcore[0,1] = 4             #pol. mesh pts from top to x-point on outside
com.nxleg[0,1] = 4              #pol. mesh pts from x-point to outer plate
com.nysol[0] = 4                #rad. mesh pts in SOL
com.nycore[0] = 4               #rad. mesh pts in core
flx.alfcy = 2                   #>0 concentrates y-mesh near separatrix

# Mesh construction--non orthogonal mesh
##com.ismmon = 3          #controls transition from nonorthog to orthog mesh
##com.isnonog = 1         # non orthogonal differencing
##grd.istream = 0
##grd.iplate = 1
##grd.nsmooth = 3
##grd.wtmesh1 = 0.75
##grd.dmix0 = 1.0
# Set params for line-segments defining inner(1) and outer(2) plots
##grd.nplate1 = 2
##grd.nplate2 = 2

# Finite-difference algorithms (upwind, central diff, etc.)
bbb.methn = 33          #ion continuty eqn
bbb.methu = 33          #ion parallel momentum eqn
bbb.methe = 33          #electron energy eqn
bbb.methi = 33          #ion energy eqn
##bbb.methg = 66
bbb.methg = 33          

# Boundary conditions
bbb.ncore[0] = 2.e19    #hydrogen ion density on core
bbb.iflcore = 0         #flag; =0, fixed Te,i; =1, fixed power on core
bbb.tcoree = 600.       #core Te
bbb.tcorei = 600.       #core Ti
bbb.tedge = 2.          #fixed wall,pf Te,i if istewcon=1, etc
bbb.recycp[0] = 0.95    #hydrogen recycling coeff at plates

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
bbb.isph_sput[1]=3 	    # physical sputtering model

bbb.t_wall = 300.
bbb.t_plat = 500.

bbb.crmb = 2.


## Impurity ions
bbb.isimpon = 6             #Use force-balance only
com.nzsp[0] = 6             #number chrg states impurity isotope #1

bbb.csfaclb[2,0] = 2.191
bbb.csfaclb[3,0] = 2.191
bbb.csfaclb[4,0] = 2.191
bbb.csfaclb[5,0] = 2.191
bbb.csfaclb[6,0] = 2.191
bbb.csfaclb[7,0] = 2.191
bbb.csfacrb[2,0] = 2.191
bbb.csfacrb[3,0] = 2.191
bbb.csfacrb[4,0] = 2.191
bbb.csfacrb[5,0] = 2.191
bbb.csfacrb[6,0] = 2.191
bbb.csfacrb[7,0] = 2.191
bbb.csfaclb[2,1] = 2.191
bbb.csfaclb[3,1] = 2.191
bbb.csfaclb[4,1] = 2.191
bbb.csfaclb[5,1] = 2.191
bbb.csfaclb[6,1] = 2.191
bbb.csfaclb[7,1] = 2.191
bbb.csfacrb[2,1] = 2.191
bbb.csfacrb[3,1] = 2.191
bbb.csfacrb[4,1] = 2.191
bbb.csfacrb[5,1] = 2.191
bbb.csfacrb[6,1] = 2.191
bbb.csfacrb[7,1] = 2.191
bbb.minu[2] = 12.
bbb.minu[3] = 12.
bbb.minu[4] = 12.
bbb.minu[5] = 12.
bbb.minu[6] = 12.
bbb.minu[7] = 12.
bbb.ziin[2] = 1
bbb.ziin[3] = 2
bbb.ziin[4] = 3
bbb.ziin[5] = 4
bbb.ziin[6] = 5
bbb.ziin[7] = 6
bbb.znuclin[0] = 1
bbb.znuclin[1] = 1
bbb.znuclin[2] = 6
bbb.znuclin[3] = 6
bbb.znuclin[4] = 6
bbb.znuclin[5] = 6
bbb.znuclin[6] = 6
bbb.znuclin[7] = 6
bbb.n0[2] = 1.e17
bbb.n0[3] = 1.e17
bbb.n0[4] = 1.e17
bbb.n0[5] = 1.e17
bbb.n0[6] = 1.e17
bbb.n0[7] = 1.e17

bbb.nzbackg = 1.e9          #background density for impurities
bbb.inzb = 2		    #exponent for switching on nzbackg
bbb.ismctab = 2             # use Braams" rate tables
com.mcfilename[0] = "C_rates.strahl"     # Imp rate file name

bbb.isnicore[7] = 3
bbb.curcore[7] = 0.

bbb.isnwcono[2] = 3
bbb.isnwcono[3] = 3
bbb.isnwcono[4] = 3
bbb.isnwcono[5] = 3
bbb.isnwcono[6] = 3
bbb.isnwcono[7] = 3
bbb.isnwconi[2] = 3
bbb.isnwconi[3] = 3
bbb.isnwconi[4] = 3
bbb.isnwconi[5] = 3
bbb.isnwconi[6] = 3
bbb.isnwconi[7] = 3
bbb.nwomin[2] = 1.e7
bbb.nwomin[3] = 1.e7
bbb.nwomin[4] = 1.e7
bbb.nwomin[5] = 1.e7
bbb.nwomin[6] = 1.e7
bbb.nwomin[7] = 1.e7
bbb.nwimin[2] = 1.e7
bbb.nwimin[3] = 1.e7
bbb.nwimin[4] = 1.e7
bbb.nwimin[5] = 1.e7
bbb.nwimin[6] = 1.e7
bbb.nwimin[7] = 1.e7

bbb.restart = 1             #Begin from savefile, not estimated profiles

# Filling newly allocated arrays as desired
bbb.ftol = 1.e-8
bbb.allocate()      #allocates storage for arrays
from uedge.hdf5 import *
hdf5_restore("d3dHm_Cog_strahl_Rp95.h5")

# Atomic data switches
com.istabon = 10            #=10 specifics hydrogen data file ehr2.dat

# Scale factor converting (upi-upg)**2 energy to thermal energy
bbb.cfnidh = 0.2
