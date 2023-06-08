#Coarse mesh [com.nx=8, com.ny=4] for DIII-D MHD equilibrium
#Uses diffusive neutrals, so five variables [bbb.ni,bbb.upi,Te,Ti,bbb.ng]
#
import sys,os
from uedge import *
# Set the com.geometry
bbb.mhdgeo = 1 		#use MHD equilibrium
#com.rm -flx.f aeqdsk neqdsk		#change names of MHD eqil. files 

flx.psi0min1 = 0.98		#normalized flux on core bndry
flx.psi0min2 = 0.98		#normalized flux on pf bndry
flx.psi0sep = 1.00001	#normalized flux at separatrix
flx.psi0max = 1.07		#normalized flux on outer wall bndry
bbb.ngrid = 1		#number of meshes [always set to 1]
com.nxleg[0,0] = 4		#pol. mesh pts from inner plate to flx.x-point
com.nxcore[0,0] = 4		#pol. mesh pts from flx.x-point to top on inside
com.nxcore[0,1] = 4		#pol. mesh pts from top to flx.x-point on outside
com.nxleg[0,1] = 4		#pol. mesh pts from flx.x-point to outer plate
com.nysol[0] = 6		#rad. mesh pts in SOL
com.nycore[0] =2		#rad. mesh pts in core

# Finite-difference algorithms [upwind, central diff, etc.]
bbb.methn = 33		#ion continuty eqn
bbb.methu = 33		#ion parallel momentum eqn
bbb.methe = 33		#electron energy eqn
bbb.methi = 33		#ion energy eqn
bbb.methg = 33		#neutral gas continuity eqn

# Boundary conditions
bbb.ncore[0] = 2.e19	#hydrogen ion density on core
bbb.iflcore = 0		#flag =0, fixed Te,i =1, fixed power on core
bbb.tcoree = 100.		#core Te 
bbb.tcorei = 100.		#core Ti
bbb.tedge = 2.		#fixed wall,pf Te,i if istewcon=1, etc
bbb.istepfc = 3		#=3 sets scale length bbb.lyte
bbb.lyte = 0.03		#radial scale-length for Te on PF boundary
bbb.recycp[0] = 0.98		#hydrogen recycling grd.coeff at plates
bbb.recycw[0] = 0.9		#wall recycling if bbb.matwso,i=1
bbb.matwso[0] = 1		#recycle on main-chamber wall
bbb.isnwcono = 1		#if=1, set bbb.ni[,,com.ny+1]=bbb.nwallo
bbb.isnwconi = 1		#if=1, set PF wall bbb.ni=bbb.nwalli
bbb.allocate()		#bbb.allocate() space of bbb.nwallo,i
bbb.nwallo = 1.e18
bbb.nwalli = 1.e18

# Transport coefficients
bbb.difni[0] = 1.		#D for radial hydrogen diffusion
bbb.kye = 1.		#chi_e for radial elec energy diffusion
bbb.kyi = 1.		#chi_i for radial ion energy diffusion
bbb.travis[0] = 1.		#eta_a for radial ion momentum diffusion

# Flux limits
bbb.flalfe = 0.21		#electron parallel thermal conduct. grd.coeff
bbb.flalfi = 0.21		#ion parallel thermal conduct. grd.coeff
bbb.flalfv = 1.		#ion parallel viscosity grd.coeff
bbb.flalfgx = 1.		#neut. density poloidal diffusion
bbb.flalfgy = 1.		#neut. density radial diffusion
bbb.flalfvgx = 1.		#neut. momentum poloidal viscosity
bbb.flalfvgy = 1.		#neut. momentum radial viscosity
bbb.flalftgx = 1.		#neut. particle poloidal heat diffusion 
bbb.flalftgy = 1.		#neut. particle radial heat diffusion

# Solver package
bbb.svrpkg = "nksol"	#Newton solver using Krylov method
bbb.premeth = "ilut"	#Solution method for precond. Jacobian matrix
bbb.mfnksol = 3

# Parallel neutral momentum equation
bbb.ineudif = 2
bbb.isupgon[0]=1
if bbb.isupgon[0] == 1:
   bbb.isngon[0]=0
   com.ngsp=1
   com.nhsp=2
   bbb.ziin[com.nhsp-1]=0
   bbb.travis[1] = 0.	#shouldn't be used for neutrals - to be sure
   # The following are probably default, set them anyway to be sure
   bbb.cngmom=0
   bbb.cmwall=0
   bbb.cngtgx=0
   bbb.cngtgy=0
   bbb.kxn=0
   bbb.kyn=0
##	bbb.ingb = 0

# Currents and potenial parameters
bbb.isphion = 1
bbb.b0 = 1.			# =1 for normal direction of B-field
bbb.rsigpl=1.e-8		#anomalous cross-field conductivity
bbb.cfjhf=1.		#turn-on heat flow from current [bbb.fqp]
bbb.cfjve = 1.		#makes bbb.vex = vix - bbb.cfjve*bbb.fqx
bbb.jhswitch=1		#Joule Heating switch
bbb.cfjpy = 0.		#diamag. cur. in flx.y-direction
bbb.cfjp2 = 0.		#diamag. cur. in 2-direction

bbb.newbcl=1
bbb.newbcr=1	#Sheath BC [bee,i] from current equation
bbb.isfdiax =1.		#Factor to turn on diamag. contrib. to sheath

bbb.cfyef = 1.0		#ExB drift in flx.y-dir.
bbb.cf2ef = 1.0		#ExB drift in 2-dir.
bbb.cfydd = 0.		#Diamag. drift in flx.y-dir. [always=0]
bbb.cf2dd = 0.		#Diamag. drift in 2-dir. [always=0]
bbb.cfrd = 0.		#Resistive drift in flx.y and 2 dirs.
bbb.cfbgt = 0.		#Diamag. energy drift [always=0]
bbb.cfybf = 1.              #turns on bbb.vycb - radial grad_B drift
bbb.cf2bf = 1.              #turns on bbb.v2cb - perp grad_B drift [nearly pol]
bbb.cfqybf = 1.             #turns on bbb.vycb contrib to radial current
bbb.cfq2bf = 1.             #turns on bbb.v2cb contrib to perp["2"] current
bbb.isnewpot = 1
bbb.rnewpot = 1.		
bbb.iphibcc = 3		#set bbb.phi[,1] uniform poloidally		

# Restart from bbb.a pfb savefile
bbb.restart = 1		#Begin from savefile, not estimated profiles
bbb.allocate()		#bbb.allocate() space for savevariables
