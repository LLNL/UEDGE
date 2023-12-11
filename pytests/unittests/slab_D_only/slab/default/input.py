# This is the test header file
# It is copied into the start of every test case
# created from this subdirectory. It is the 
# best place to include lines common to
# all test cases, such as rate paths etc.

from uedge import *
from uedge.hdf5 import hdf5_restore

try:
    bbb.oldseec = 0
    bbb.jhswitch = 0
except:
    pass

bbb.isfixlb[0]=0    # =1 fixes values on left boundary (nib, upb, teb, tib, yylb0)
                    # =2 for symmetry point
grd.radx = 4.5e-2       #outer "radial" wall
grd.rad0 = 0.
grd.radm = -5.0e-4      #minimum "radial" position
grd.alfyt=-2.5          #radial nonuniformity factor <0 => expanding
grd.za0 = 0.            #poloidal symmetry plane location
grd.zaxpt = 1.0         #poloidal location of flx.x-point
grd.zax = 1.5           #poloidal location of divertor plate
grd.alfxt=5.5

grd.btfix = 2.          #constant total B-field
grd.bpolfix = .2        #constant poloidal B-field
bbb.ngrid = 1
com.nxleg[0,0]=0
com.nxcore[0,0]=0
com.nxcore[0,1]=8
com.nxleg[0,1]=24
com.nycore[0]=1
com.nysol=16
com.isgrdsym = 1        #symmeterize the grid about z=0
bbb.isnfmiy = 1         #diff bbb.fmiy to preserve symmetry about X-pt

# Finite-difference algorithms [upwind, central diff, etc.]
bbb.methn = 33		#ion continuty eqn
bbb.methu = 33		#ion parallel momentum eqn
bbb.methe = 33		#electron energy eqn
bbb.methi = 33		#ion energy eqn
bbb.methg = 33		#neutral gas continuity eqn

# templates/slab_D_only/inputs/atoms/default.py
# Neutral gas propeties - Inertial neutrals
com.nhsp=2
bbb.ziin[1]=0.
bbb.ineudif = 2		    # use neutral pressure to calc gas vel
bbb.cfnidh = 0		    #=0 for no ion heating from cx bbb.up bbb.friction
# Atomic physics packages
com.istabon=10		    #Stotler tables


# templates/slab_D_only/inputs/boundary/core/density/default.py
bbb.ncore = 2.e19		#value of bbb.ncore at single point
bbb.isnicore = 1		#switches off bbb.ni=bbb.ncore for all separatrix pts.

# templates/slab_D_only/inputs/boundary/core/energy/default.py
bbb.iflcore = 1		    #specify core power 
bbb.pcoree = 1.e5
bbb.pcorei = 1.e5

# templates/slab_D_only/inputs/boundary/core/momentum/default.py
bbb.isupcore = 0		#parallel vel has bbb.up=0 on core bdry

# templates/slab_D_only/inputs/boundary/plates/energy/default.py
bbb.bcee = 5.
bbb.bcei = 2.5	        #energy transmission coeffs.

# templates/slab_D_only/inputs/boundary/plates/momentum/default.py
bbb.recycm = -0.5		# if close to -1, upg~bbb.upi on plates
bbb.isupss = 0		    #parallel vel sonic

# templates/slab_D_only/inputs/boundary/plates/recycling/default.py

# Add asymmetry for recycling at plates after bbb.allocate()
bbb.recycp = 0.99		#base plate recycling grd.coeff.
bbb.recycflb = 0.98		#additional recycling factor, left plate
bbb.recycfrb = 1.00		#additional recycling factor, right plate

# templates/slab_D_only/inputs/boundary/walls/energy/default.py
bbb.istepfc=0
bbb.istipfc=0	        #bbb.priv. flux has api.zero temp. deriv.
bbb.istewc=0
bbb.istiwc=0	        #wall has api.zero temp. deriv.

# templates/slab_D_only/inputs/curpot/default.py
# Currents and potential parameters
bbb.b0 = 1.		        # =1 for normal direction of B-field
bbb.rsigpl=1.e-9		#anomalous cross-field conductivity
bbb.cfjhf=1.	   	    #turn-on heat flow from current [bbb.fqp]
bbb.cfjve = 1.		    #makes bbb.vex = vix - bbb.cfjve*bbb.fqx
bbb.cfjpy = 0.		    #diamag. cur. in flx.y-direction
bbb.cfjp2 = 0.		    #diamag. cur. in 2-direction

bbb.newbcl = 1		    #Sheath BC [bee,i] from current equation
bbb.newbcr = 1		    #Sheath BC [bee,i] from current equation
bbb.isfdiax =0.		    #Factor to turn on diamag. contrib. to sheath

bbb.cfyef = 0.6	    	#ExB drift in flx.y-dir.
bbb.cf2ef = 0.6		    #ExB drift in 2-dir.
bbb.cfybf = 1.0e-20     #turns on bbb.vycb - radial grad_B drift
bbb.cf2bf = 1.0e-20     #turns on bbb.v2cb - perp grad_B drift [nearly pol]
bbb.cfqybf = 1.0e-20    #turns on bbb.vycb contrib to radial current
bbb.cfq2bf = 1.0e-20    #turns on bbb.v2cb contrib to perp["2"] current
bbb.cftef = 1.0e-20		#turns on bbb.v2ce for toroidal velocity
bbb.cftdd = 1.0e-20		#turns on bbb.v2dd [diagm vel] for toroidal vel
bbb.cfqym = 1.0e-20		#turns on inertial correction to bbb.fqy current
bbb.cfydd = 0.		    #Diamag. drift in flx.y-dir. [always=0]
bbb.cf2dd = 0.		    #Diamag. drift in 2-dir. [always=0]
bbb.cfrd = 0.		    #Resistive drift in flx.y and 2 dirs.
bbb.cfbgt = 0.	    	#Diamag. energy drift [always=0]
bbb.isnewpot = 0
bbb.rnewpot = 0.

# templates/slab_D_only/inputs/equations/default.py
bbb.isnion = 1
bbb.isupon = 1
bbb.isteon = 1
bbb.istion = 1
bbb.isphion = 1
bbb.isphiofft = 0
bbb.istgon = 0
bbb.isngon = 0
bbb.isngon[0] = 0
bbb.isupgon = 0
bbb.isupgon[0] = 1

# templates/slab_D_only/inputs/fluxlimits/default.py

# Flux limits
bbb.flalfe = 0.21       #electron parallel thermal conduct. grd.coeff
bbb.flalfi = 0.21       #ion parallel thermal conduct. grd.coeff
bbb.flalfv = 1.         #ion parallel viscosity grd.coeff
bbb.flalfgx = 1.        #neut. gas in poloidal direction
bbb.flalfgy = 1.        #neut. gas in radial direction
bbb.flalfgxy = 1.   	#nonorthogonal neut. flux limit
bbb.flalftgx = 1.		#neut. energy pol flux limit factor
bbb.flalftgy = 1.		#neut. energy rad flux limit factor
bbb.flalfvgx = 1.e20	#neut. viscosity pol flux limit factor
bbb.flalfvgy = 1.e20	#neut. viscosity rad flux limit factor

# templates/slab_D_only/inputs/radialtransport/default.py
# Transport coefficients
bbb.difni[0] = 0.333333
bbb.kye = 0.5
bbb.kyi = 0.5
bbb.travis[0]=1.
bbb.parvis[0]=1.
bbb.difutm = 0.4

# templates/slab_D_only/inputs/solver/default.py
bbb.svrpkg="nksol"
bbb.mfnksol = -3
bbb.epscon1 = .005
bbb.ftol = 1.e-8
bbb.premeth="ilut"
bbb.runtim=1.e-7
bbb.rlx=.4

# templates/restore/default.py
bbb.restart = 1			#Begin from savefile
bbb.allocate()			#allocate plasma/neutral arrays
hdf5_restore("solution.h5")    #read in the solution from pfb file