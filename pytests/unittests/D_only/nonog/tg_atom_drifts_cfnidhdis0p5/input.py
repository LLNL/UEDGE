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

bbb.mhdgeo = 1 		#use MHD equilibrium
bbb.gengrid=0		#read mesh from gridue file
com.geometry = "snull"
com.isnonog = 1		#nonorthogal differencing used

bbb.methn = 33		#ion continuty eqn
bbb.methu = 33		#ion parallel momentum eqn
bbb.methe = 33		#electron energy eqn
bbb.methi = 33		#ion energy eqn
bbb.methg = 66		#neutral gas continuity eqn

# templates/D_only/inputs/atoms/default.py
bbb.isupwo[1] = 0
bbb.ineudif = 2
com.ngsp=1
com.nhsp=2
bbb.ziin[1]=0
bbb.travis[1] = 0.	#shouldn't be used for neutrals - but to be sure

# templates/D_only/inputs/boundary/core/density/default.py
bbb.isnicore[0] = 1 	#use fixed-density BC on core
bbb.ncore[0] = 2.e19	#hydrogen ion density on core

# templates/D_only/inputs/boundary/core/energy/default.py
bbb.iflcore = 0		#flag=0, fixed Te,i;=1, fixed power on core
bbb.tcoree = 100.	#core Te 
bbb.tcorei = 100.	#core Ti

# templates/D_only/inputs/boundary/plates/particle/default.py
bbb.recycp[0] = 0.9		#hydrogen recycling grd.coeff at plates

# templates/D_only/inputs/boundary/walls/energy/default.py
bbb.istewc = 1 		#=1 sets constant wall te
bbb.tedge = 2.		#fixed wall,pf Te,i if istewcon=1, etc
bbb.istepfc = 3		#=3 sets scale length bbb.lyte
bbb.lyte = 0.03		#radial scale-length for Te on PF boundary


# templates/D_only/inputs/boundary/walls/particle/default.py
bbb.matwso[0] = 1		#recycle on main-chamber wall
bbb.isnwcono = 1		#if=1, set bbb.ni[,,com.ny+1]=bbb.nwallo
bbb.isnwconi = 1		#if=1, set PF wall bbb.ni=bbb.nwalli
bbb.allocate()		#bbb.allocate() space of bbb.nwallo,i
bbb.nwallo = 1.e18
bbb.nwalli = 1.e18

# templates/D_only/inputs/boundary/walls/recycling/default.py
bbb.recycw[0] = 0.9		#wall recycling if bbb.matwso,i=1

# templates/D_only/inputs/curpot/tg_atom_drifts_cfnidhdis0p5.py
bbb.istgon[0] = 1   # Turn on D0 temperature equation
bbb.cftiexclg = 0.  # Remove the Tg part in the Ti equation
bbb.cfdiss = 1.0
bbb.istgcore[0] = 2 #..albedo
bbb.istgpfc[0] = 4#2  #..zml flow change to 4 later
bbb.istgwc[0] = 4#2   #..flow
bbb.istglb[0] = 4
bbb.istgrb[0] = 4

bbb.recyce = 0.3
bbb.recycm = -0.9

bbb.cfnidhgy = 1.0
bbb.cfnidhg2 = 1.0
bbb.cfnidhdis = 0.5 #..consider drift heating for dissociation
bbb.cfnidhmol = 1.0 #..ignore molecular v in dissociation drift heating


bbb.ncore[0] = 2e18
# templates/D_only/inputs/boundary/walls/recycling/default.py
bbb.recycw[0:2] = [0.25, 0.95]		#wall recycling if bbb.matwso,i=1

# gas pumping
bbb.nwsor = 1  #.. number of sources (sinks).  defalut: =1
  #.. PRF
bbb.igasi[0] = 0.   #.. pumping. default =0.
bbb.xgasi[0] = 0.   #.. position from inner. default =0.
bbb.wgasi[0] = 100. #.. width. default =100.(m)
bbb.igspsori[0] = 1 #.. species. default =1
bbb.albdsi[0] = 0.25 #.. pumping flux = (1.-albdsi)*1/4*n*sqrt(8T/pi/me)
  #.. wall
bbb.igaso[0] = 0.   #.. default =0.
bbb.xgaso[0] = 0.   #.. position from inner. default =0.
bbb.wgaso[0] = 100. #.. width. default =100.(m)
bbb.igspsoro[0] = 1 #.. species. default =1
bbb.albdso[0] = 0.25 #.. pumping flux = (1.-albdsi)*1/4*n*sqrt(8T/pi/me)
  #.. target pumping for atoms
bbb.albedorb[0,0] = 0.25
bbb.albedolb[0,0] = 0.25


# Drifts and currents
bbb.isphion = 1
bbb.b0 = 1.0                    # =1 for normal direction of B-field
bbb.rsigpl=1.e-8                #anomalous cross-field conductivity
bbb.cfjhf=1.            #turn-on heat flow from current [bbb.fqp]
bbb.cfjve = 1.          #makes bbb.vex = vix - bbb.cfjve*bbb.fqx
bbb.cfjpy = 0.          #diamag. cur. in flx.y-direction
bbb.cfjp2 = 0.          #diamag. cur. in 2-direction
bbb.cfsigm = 1.0
#
bbb.newbcl=1
bbb.newbcr=1
bbb.isfdiax =1.         #Factor to turn on (ExB+diamag.) contrib. to sheath
#
bbb.cfyef = 1.0         #ExB drift in flx.y-dir.
bbb.cf2ef = 1.0         #ExB drift in 2-dir.
bbb.cfydd = 0.          #Diamag. drift in flx.y-dir. [always=0]
bbb.cf2dd = 0.          #Diamag. drift in 2-dir. [always=0]
bbb.cfrd = 0.           #Resistive drift in flx.y and 2 dirs.
bbb.cfbgt = 0.          #Diamag. energy drift [always=0]
bbb.cfybf = 1.0              #turns on bbb.vycb - radial grad_B drift
bbb.cf2bf = 1.0              #turns on bbb.v2cb - perp grad_B drift [nearly pol]
bbb.cfqybf = 1.0             #turns on bbb.vycb contrib to radial current
bbb.cfq2bf = 1.0             #turns on bbb.v2cb contrib to perp["2"] current
#..modified
bbb.cfqybbo=1.     #turn off Grad B current on boundary
bbb.cfqydbo=0.     #use full diagmagetic current on boundary to force j_r=0

bbb.cfniybbo=1.   # use to avoid artificial source at core boundary
bbb.cfniydbo=0.   # use to avoid artificial source at core boundary
bbb.cfeeybbo=1.   # ditto
bbb.cfeeydbo=0.   # ditto
#bbb.cfeixdbo=1.   # turn on BXgrad(T) drift in plate BC
#bbb.cfeexdbo=1.   # turn on diamagnetic drift in plate BC
bbb.cfeixdbo=0.
bbb.cfeexdbo=0.
bbb.cftef=1.0     #turns on v2ce for toroidal velocity
bbb.cftdd=1.0     #turns on v2dd (diamag vel) for toloidal velocity
bbb.cfqym=1.0     #turns on inertial correction to fqy current
#
bbb.isnewpot = 1
bbb.rnewpot = 1.
#bbb.iphibcc = 0         #set bbb.phi[,1] uniform poloidally
bbb.iphibcc=3     # =3 gives ey=eycore on core bdry
bbb.iphibcwi=0    # set ey=0 on inner wall if =0
                  # phi(PF)=phintewi*te(ix,0) on PF wall if =1
bbb.iphibcwo=0    # same for outer wall
bbb.isutcore=2    # =1, set dut/dy=0 on iy=0 (if iphibcc=0)
                # =0, toroidal angular momentum=lzcore on iy=0 (iphibcc=0)


# templates/D_only/inputs/equations/default.py
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


# Catch-all for turning off potential equaiton in slab geometry
if bbb.mhdgeo == -1:
	bbb.isphion = 0
	bbb.isphiofft = 1

# templates/D_only/inputs/fluxlimits/default.py
bbb.flalfe = 0.21	#electron parallel thermal conduct. grd.coeff
bbb.flalfi = 0.21	#ion parallel thermal conduct. grd.coeff
bbb.flalfv = 1.		#ion parallel viscosity grd.coeff
bbb.flalfgx = 1.	#neut. density poloidal diffusion
bbb.flalfgy = 1.	#neut. density radial diffusion
bbb.flalfvgx = 1.	#neut. momentum poloidal viscosity
bbb.flalfvgy = 1.	#neut. momentum radial viscosity
bbb.flalftgx = 1.	#neut. particle poloidal heat diffusion 
bbb.flalftgy = 1.	#neut. particle radial heat diffusion

# templates/D_only/inputs/radialtransport/default.py
bbb.difni[1] = 1.	#D for radial hydrogen diffusion
bbb.kye = 1.		#chi_e for radial elec energy diffusion
bbb.kyi = 1.		#chi_i for radial ion energy diffusion
bbb.travis[1] = 1.	#eta_a for radial ion momentum diffusion

# templates/D_only/inputs/solver/default.py
bbb.svrpkg = "nksol"	#Newton solver using Krylov method
bbb.premeth = "ilut"	#Solution method for precond. Jacobian matrix
bbb.mfnksol = 3

# templates/restore/default.py
bbb.restart = 1			#Begin from savefile
bbb.allocate()			#allocate plasma/neutral arrays
hdf5_restore("solution.h5")    #read in the solution from pfb file