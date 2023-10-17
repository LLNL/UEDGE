###########################################################################
# DESCRIPTION OF PROBLEM (box2):
#
# This is a Python version of the box case from Andreas Holm
###########################################################################

from uedge import bbb, com, grd


#-Geometry
bbb.mhdgeo=-1 #-set cartesian geometry
bbb.isfixlb=2 #left boundary as sym. plane; no flux at cut

grd.radx= 4.e-2    #-outer "radial" wall
grd.rad0=0.0       #-location of 'radial' separ'x for cylinder or slab
grd.radm=-1.e-2    #-minimum "radial" position

grd.za0 = 0.    #-poloidal symmetry plane location
grd.zax=3.0     #-poloidal location of divertor plate   
grd.zaxpt=2.25  #-poloidal location of x-point
grd.alfyt=-2.0  #radial nonuniformity factor; < 0 => expanding
grd.alfxt=2.76  #poliodal nonuniformity factor; to make smooth
                #transition to exp. grid, alfxt should satisfy 
                #the eqn dzun = (zax-zaxpt+dzun)
                #               (1-exp(-alfxt/(nx-ixpt2+1))) /
                #               (1-exp(-alfxt))
                #where dzun = (zaxpt-za0)/ixpt2 and 
                #ixpt2 = ncore(1,2).

grd.btfix = 2.              #constant total B-field
grd.bpolfix = .2            #constant poloidal B-field


#-Grid
bbb.gengrid=1; #-Note: for slab the grid is not saved in gridue
bbb.ngrid=1
com.nycore[0]=2
com.nysol[0]=4
com.nxleg[0,1]=3
com.nxcore[0,1]=3



#-Boundary conditions
bbb.isnicore[0]=1 #-same density at all core points
bbb.ncore=1.1e19  #-density on core boundary
bbb.iflcore=1     #if=1, specify core power
bbb.tcoree=25.0   #-used if iflcore=0
bbb.tcorei=25.0   #-used if iflcore=0
bbb.pcoree = 2.5e4  #-used if iflcore=1        
bbb.pcorei = 2.5e4  #-used if iflcore=1

bbb.recycp=0.98  #-recycling coef at plates if ndatlb,rb=0
bbb.albdsi=0.99  #-albedos at inner gas source locations
bbb.albdso=0.99  #-albedos at inner gas source locations

bbb.istepfc=0; bbb.istipfc=0 #-priv. flux has zero temp. deriv.
bbb.istewc=0;  bbb.istiwc=0  #-wall has zero temp. deriv.
bbb.bcee = 4.; bbb.bcei = 2.5  #-energy transmission coeffs.        
bbb.bcen = 0.      #-energy transmission coefficint for neutrals 
bbb.isupss = 0     #-parallel vel sonic
bbb.isupcore = 0   #-parallel vel =0 on core bndry



#-Transport coefficients
bbb.difni=0.5
bbb.kye=0.7
bbb.kyi=0.7
bbb.travis=1.0
bbb.parvis=1.0



#-Flux limits
bbb.flalfe=0.2
bbb.flalfi=0.2
bbb.flalfgx=1.e0
bbb.flalfgy=1.e0
bbb.flalfgxy=1.e0
bbb.flalfv=0.5


# Finite difference algorithms
bbb.methe=33;bbb.methu=33;bbb.methg=33
bbb.methn=33;bbb.methi=33


#-Solver package
bbb.svrpkg = "nksol"    #Newton solver using Krylov method
bbb.mfnksol=-3
bbb.epscon1=0.005
bbb.ftol=1e-10
bbb.premeth = "ilut"    #Solution method for precond. Jacobian matrix
bbb.runtim=1e-07
bbb.rlx=0.9
###bbb.del=1.e-8 #-this one causes syntax error!


#-Neutral gas propeties
#bbb.tfcx=5.;bbb.tfcy=5.         #Franck-Condon temperatures
bbb.cngfx=1.;bbb.cngfy=1.       #turn-on grad(T_g) flux if =1
bbb.cngflox=1.;bbb.cngfloy=0.   #turn-on drift with ions if =1
bbb.cngmom = 1.             #ion-gas momentum transfer
bbb.eion = 5.               #birth energy of ions
bbb.ediss = 10.             #dissoc. energy lost from elecs (eion=2*ediss)
bbb.isrecmon = 1            #=1 turns on recombination
bbb.cfupcx=1.0              # factor multipling momentum cx
bbb.cfticx=1.0              # factor multipling cx terms in ion energy eqn



#-Parallel neutral momentum equation
bbb.isupgon[0]=1

if (bbb.isupgon[0] == 1):
    bbb.isngon=0               
    com.ngsp=1                
    com.nhsp=2                
    ###bbb.ziin[com.nhsp-1]=1
    bbb.ziin[0]=1
    bbb.ziin[1]=0

    #-the following are probably default, set them anyway to be sure
    bbb.cngmom=0
    bbb.cmwall=0
    bbb.cngtgx=0
    bbb.cngtgy=0
    bbb.kxn=0
    bbb.kyn=0


#-Currents and potential parameters
bbb.isphion=0
bbb.rsigpl=1.e-8            #anomalous cross-field conductivity
bbb.cfjhf=0.                #turn-on heat flow from current (fqp)
bbb.jhswitch=0              #Joule Heating switch


# Atomic physics packages
#com.istabon=10              #DEGAS rates
com.istabon=0              #-analytic rates


#-Misc
bbb.restart=0
