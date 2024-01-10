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

# templates/D_mols/inputs/atoms/default.py
bbb.isupwo[1] = 0
bbb.ineudif = 2
com.ngsp=1
com.nhsp=2
bbb.ziin[1]=0
bbb.travis[1] = 0.	#shouldn't be used for neutrals - but to be sure

# templates/D_mols/inputs/boundary/core/density/default.py
bbb.isnicore[0] = 1 	#use fixed-density BC on core
bbb.ncore[0] = 2.e19	#hydrogen ion density on core

# templates/D_mols/inputs/boundary/core/energy/default.py
bbb.iflcore = 0		#flag=0, fixed Te,i;=1, fixed power on core
bbb.tcoree = 100.	#core Te 
bbb.tcorei = 100.	#core Ti

# templates/D_mols/inputs/boundary/plates/particle/default.py
bbb.recycp[0] = 0.9		#hydrogen recycling grd.coeff at plates
#.. target pumping for atoms
bbb.albedorb[0,0] = 0.25
bbb.albedolb[0,0] = 0.25

# templates/D_mols/inputs/boundary/walls/energy/default.py
bbb.istewc = 1 		#=1 sets constant wall te
bbb.tedge = 2.		#fixed wall,pf Te,i if istewcon=1, etc
bbb.istepfc = 3		#=3 sets scale length bbb.lyte
bbb.lyte = 0.03		#radial scale-length for Te on PF boundary


# templates/D_mols/inputs/boundary/walls/particle/default.py
bbb.matwso[0] = 1		#recycle on main-chamber wall
bbb.isnwcono = 1		#if=1, set bbb.ni[,,com.ny+1]=bbb.nwallo
bbb.isnwconi = 1		#if=1, set PF wall bbb.ni=bbb.nwalli
bbb.allocate()		#bbb.allocate() space of bbb.nwallo,i
bbb.nwallo = 1.e18
bbb.nwalli = 1.e18


# gas pumping
bbb.nwsor = 2  #.. number of sources (sinks).  defalut: =1
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

#..new pumping for atoms
  #.. PRF
bbb.igasi[1] = 0.   #.. pumping. default =0.
bbb.xgasi[1] = 0.   #.. position from inner. default =0.
bbb.wgasi[1] = 100. #.. width. default =100.(m)
bbb.igspsori[1] = 2 #.. species. default =1
bbb.albdsi[1] = 0.9999 #.. pumping flux = (1.-albdsi)*1/4*n*sqrt(8T/pi/me)
  #.. wall
bbb.igaso[1] = 0.   #.. default =0.
bbb.xgaso[1] = 0.   #.. position from inner. default =0.
bbb.wgaso[1] = 100. #.. width. default =100.(m)
bbb.igspsoro[1] = 2 #.. species. default =1
bbb.albdso[1] = 0.9999 #.. pumping flux = (1.-albdsi)*1/4*n*sqrt(8T/pi/me)
  #.. recycling
bbb.matwsi[1] = 1   #.. default =0
#bbb.recycw[0] = 0.9    #.. wall recycling coef. for igsp=1 (for PRF and wall)
bbb.matwso[1] = 1   #.. coefficient is recycw

# templates/D_mols/inputs/boundary/walls/recycling/default.py
bbb.recycw[0:2] = [0.25, 0.95]

# templates/D_mols/inputs/equations/default.py
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

# templates/D_mols/inputs/fluxlimits/default.py
bbb.flalfe = 0.21	#electron parallel thermal conduct. grd.coeff
bbb.flalfi = 0.21	#ion parallel thermal conduct. grd.coeff
bbb.flalfv = 1.		#ion parallel viscosity grd.coeff
bbb.flalfgx = 1.	#neut. density poloidal diffusion
bbb.flalfgy = 1.	#neut. density radial diffusion
bbb.flalfvgx = 1.	#neut. momentum poloidal viscosity
bbb.flalfvgy = 1.	#neut. momentum radial viscosity
bbb.flalftgx = 1.	#neut. particle poloidal heat diffusion 
bbb.flalftgy = 1.	#neut. particle radial heat diffusion

# templates/D_mols/inputs/mols/tg_mol_istglb1.py
bbb.ishymol = 1
com.ngsp = 2
com.nhgsp = 2
bbb.isngon[1] = 1
bbb.istgcon[1] = -1
bbb.istgcore[1] = 2
bbb.istgpfc[1] = 0
bbb.istgwc[1] = 0
bbb.istglb[1] = 1
bbb.istgrb[1] = 0
bbb.tgwall[1] = 0.2

# templates/D_mols/inputs/radialtransport/default.py
bbb.difni[1] = 1.	#D for radial hydrogen diffusion
bbb.kye = 1.		#chi_e for radial elec energy diffusion
bbb.kyi = 1.		#chi_i for radial ion energy diffusion
bbb.travis[1] = 1.	#eta_a for radial ion momentum diffusion

# templates/D_mols/inputs/solver/default.py
bbb.svrpkg = "nksol"	#Newton solver using Krylov method
bbb.premeth = "ilut"	#Solution method for precond. Jacobian matrix
bbb.mfnksol = 3

# templates/restore/default.py
bbb.restart = 1			#Begin from savefile
bbb.allocate()			#allocate plasma/neutral arrays
hdf5_restore("solution.h5")    #read in the solution from pfb file