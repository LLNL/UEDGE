import os
import sys
import numpy as np
import matplotlib.pyplot as plt
###import uedge
from uedge import *
from uedge.hdf5 import *
import uedge_mvu.plot as mp
import uedge_mvu.utils as mu



def setGrid(aeqdsk='V1E_aeqdsk', geqdsk='V1E_geqdsk'):
    #execfile(os.environ['PYLIB']+"/plotmesh.py")
    #exec(open(os.environ['PYLIB']+"/plotmesh.py").read())

    #def setGrid(aeqdsk='V1E_aeqdsk', geqdsk='V1E_geqdsk'):
    # Set the geometry
    bbb.mhdgeo = 1 # =1 use MHD equilibrium files
    com.geometry = 'snull'
    bbb.gengrid = 0
    #com.isudsym = 0 # up-down symmetric

    #aeqdsk='aeqdsk'
    #geqdsk='geqdsk'

    #if aeqdsk:
     #   com.aeqdskfname = aeqdsk 
    #if geqdsk:
     #   com.geqdskfname = geqdsk

    #flx.xcutoff1 = 0.25 # set inner cutoff on the eqdsk data
    #flx.ycutoff1 = 0.0 # set lower cutoff on the eqdsk data
    #flx.psi0min1 = 0.86 # normalized flux on core bndry
    #flx.psi0min2 = 0.985 # normalized flux on pf bndry
    #flx.psi0sep = 1.00001 # normalized flux at separatrix
    #flx.psi0max = 1.065 # normalized flux on outer wall bndry
    #flx.psi0max_inner = 1.02

    #bbb.ngrid = 1 # number of mesh sequenc. (always set to 1)
    #com.nxleg[0,0] = 24 # pol. mesh pts from inner plate to x-point
    #com.nxcore[0,0] = 24 # pol. mesh pts from x-point to top on inside
    #com.nxcore[0,1] = 24 # pol. mesh pts from top to x-point on outside
    #com.nxleg[0,1] = 24 # pol. mesh pts from x-point to outer plate
    #com.nysol[0] = 18 # rad. mesh pts in SOL
    #com.nycore[0] = 18 # rad. mesh pts in core

    # Fix to avoid unnecessary interpolation
    bbb.isnintp = 1

    # Set correct strike points
    #grd.isspnew = 1 # =1 for user-defined strike point
    #grd.isztest = [0,2] # =0 for testing on rstrike only and 2 elements for inner/outer strike points
    #grd.rstrike[0] = 0.48 # inner R strike-point
    #grd.zstrike[1] = 0.22 

    # Set options for grid generator
    #flx.alfcy = 1.5
    #grd.slpxt = 1.2
    #flx.altsearch = 0
    #flx.istchkon = 1
    
    #flx.flxrun() # required for next line
    #grd.grdrun() # generate gridue file 
    #bbb.gengrid = 0 # use existing gridue file

    #com.nx=com.nxm
    #com.ny=com.nym


    #non-orthog grid
    com.ismmon=2; grd.istream=0; grd.iplate=1; grd.nsmooth=2; grd.wtmesh1=0.75; grd.dmix0=1.0
    #execfile('plate.DTT.py')
    #exec(open('plate.CMod.py').read())

    #flx.flxrun() # required for next line
    #grd.grdrun() # generate gridue file 
    #bbb.gengrid = 0 # use existing gridue file

    #com.nx=com.nxm
    #com.ny=com.nym

    com.isnonog=1; bbb.methg=66
    #execfile('plotmesh.py')
    #pm.plotmesh(iso=1)

def setBoundaryConditions():
    global bbb
    # Ion Mass
    bbb.minu[0]=2.0		#D-D plasma

    # Boundary conditions
    bbb.isnicore[0] = 1		#=1 gives fixed uniform density ncore
    #ncore(1) = 1.0e20	#ion density on core if isnicore=1
    bbb.ncore[0] = 2.0e20

    bbb.istewc = 3		#fixed Te=tewallo on outer wall
    bbb.istiwc = 3		#fixed Ti=tiwallo on outer wall
    bbb.isnwcono = 3		#fixed ni=nwallo on outer wall
    
    # Specify decay lengths at private flux wall
    bbb.isnwconi = 3
    bbb.lyni = 1e-2 
    bbb.nwimin = 1e10
        
    bbb.istepfc = 3
    bbb.lyte = 1e-2
        
    bbb.istipfc = 3 
    bbb.lyti = 1e-2
	
    bbb.isupcore = 1		#=1 sets d(up)/dy=0

    bbb.iflcore = 1		#flag; =0, fixed Te,i;
    bbb.pcoree = 1.8e6		#core elec power if iflcore=1
    bbb.pcorei = 1.8e6		#core ion power if iflcore=1

    bbb.recycp[0] = 1.0		#hydrogen recycling coeff at plates


def setSolver(experimental=False):
    global bbb
    # Solver package
    bbb.svrpkg = 'nksol'	#Newton solver using Krylov method
    bbb.premeth = 'ilut'	#Solution method for precond Jac'n matrix
    bbb.mfnksol = -3		# =3 for restricted increase in fnrm
    bbb.lfililut = 200
    bbb.lenplufac = 300	
    bbb.lenpfac = 200
    bbb.rlx = 0.1
    #del = 1.e-7
    bbb.epscon1 = 3e-3		#linear solver convergence
    bbb.n0g[0] = 1e18

    # Finite-difference algorithms (upwind, central diff, etc.)
    bbb.methn = 33 # ion continuity eqn
    bbb.methu = 33 # ion parallel momentum eqn
    bbb.methe = 33 # electron energy eqn
    bbb.methi = 33 # ion energy eqn
    bbb.methg = 66 # neutral gas continuity eqn


def setPhysics():
    global bbb, com, aph

    # Flux limits
    bbb.isflxlde = 1
    bbb.isflxldi = 2
    bbb.flalfe = 0.21		#electron parallel thermal conduct. coeff
    bbb.flalfi = 0.21		#ion parallel thermal conduct. coeff
    bbb.flalfv = 1.		#ion parallel viscosity coeff
    bbb.flalfgx = 1.		#neut. gas in poloidal direction
    bbb.flalfgy = 1.		#neut. gas in radial direction
    bbb.flalftgx = 1.		#neut power in poloidal direction
    bbb.flalftgy = 1.		#neut power in radial direction
    bbb.lgmax = 0.05		#max scale for gas particle diffusion
    bbb.lgtmax = 0.05		#max scale for gas thermal diffusion

    # Atomic physics
    #com.istabon = 0 # analytic rates
    # Atomic Physics
    com.istabon = 10		#Stotler tables
    bbb.isrecmon = 1		#=1 for recombination
    aph.aphdir='aph/'

    ### source number 1:
    bbb.jxsoro[0] = 1		# this wall source starts in mesh region 1
    bbb.jxsori[0] = 1		# this pf source starts in mesh region 1
    bbb.matwso[0] = 1.
    bbb.matwsi[0] = 1.
    bbb.wgaso[0] = 1e3
    bbb.wgasi[0] = 1e3
    
    bbb.recycw[0]=1.0

    bbb.istgwc=1
    bbb.istgpfc=2

    # Neutrals
    bbb.ineudif=1		
    bbb.cngfx=0.		#only effective when ineudif=1
    bbb.cngfy=0.		#only effective when ineudif=1
    bbb.cngflox=0.	#default#flux from convection
    bbb.cngfloy=0.	#default#flux from convection

    # Inertial neutrals
    #com.nhsp=2
    #bbb.ziin[1]=0.
    #bbb.isngon[0]=0
    #bbb.isupgon[0]=1
    bbb.cngmom=0
    bbb.cmwall=0
    bbb.cngtgx=0
    bbb.cngtgy=0
    bbb.cfbgt=0
    bbb.kxn=0
    bbb.kyn=0
    #bbb.recycm = -0.1	#latest recommmendation???

    bbb.allocate()
    bbb.nphygeo()
    #nwalli=1.e18




    bbb.kxe=1.0 #-parallel conduction factor                                                    
    bbb.islnlamcon=1 	#-Coulomb log                                                              
    bbb.lnlam=12
    bbb.isplflxlv=1 #=0, flalfv not active at ix=0 & nx;=1 active all ix                         
    bbb.isplflxlgx=1
    bbb.isplflxlgxy=1
    bbb.isplflxlvgx=1
    bbb.isplflxlvgxy=1
    bbb.iswflxlvgy=1
    bbb.isplflxltgx=1
    bbb.isplflxltgxy=1
    bbb.iswflxltgy=1 #-the main one, the only one that really matters                            
    bbb.isplflxl=1   # =0, flalfe,i not active at ix=0 & nx;=1 active all ix                     
    bbb.iswflxlvgy=1 # integer /1/ #=0, flalfvgy not active at iy=0 & ny;=1 active all iy        
    bbb.cngflox=1    # real /ngspmx*1./ #fac for x-flux from convection in ng-eqn.               
    bbb.isupwi=1     # integer /nispmx*1/ #=2 sets dup/dy=0 on inner wall                        
    bbb.isupwo=1     # integer /nispmx*1/ #=2 sets dup/dy=0 on outer wall

    

def setDChi(DOutFactor=1, DInFactor=1, chiInb=0.007, convec=False):
    """
    chiInb: constant chi on inboard side, increasing it decreases P_inner_xpt
    """
    global bbb
    # Transport coefficients (m**2/s) 
    bbb.difni[0] = 0.2		#D for radial hydrogen diffusion
    bbb.kye = 0.10		#chi_e for radial elec energy diffusion
    bbb.kyi = 0.10		#chi_i for radial ion energy diffusion
    bbb.travis[0] = 0.50	#eta_a for radial ion momentum diffusion
    #bbb.difutm = 0.50 		#toroidal diffusivity for potential
    #bbb.isbohmcalc = 0
    #bbb.facbee = 1.0 # factor for Bohm Te diffusion coeff 
    #bbb.facbei = 1.0 # factor for Bohm Ti diffusion coeff
    #bbb.facbni = 1.0 # factor for Bohm ni diffusion coeff
    #bbb.nphygeo() # calculate radial coords of grid centers, yyc
    
                   
        
    
def setup(impurity=None, experimental=False, convec=False):
    setGrid()
    setBoundaryConditions()
    setSolver(experimental=experimental)
    setPhysics()
    setDChi(convec=convec)
    
    
def restore(name='uerun'):
    global bbb
    h5file = name + '.h5'
    if os.path.isfile(h5file):
        bbb.label[0] = name 
        bbb.restart = 1
        hdf5_restore(name + '.h5')
        bbb.dtreal = 1e20
        bbb.exmain() 
    else:
        raise OSError(2, 'No such file', h5file) 
    
    
def postprocess(name='uerun'):
    if yesno('Write uerun.h5 and plots.pdf?'):
        hdf5_save(name + '.h5')
        plt.close() # for all_plots to have correct page size
        plotall(name + '.pdf')
            
            
    
            
def plotD():
    fac = 2
    Ds = [(1,1), (fac, 1), (1, fac), (fac, fac), (1/fac, 1), (1, 1/fac), (1/fac, 1/fac), (1/fac, fac), (fac, 1/fac)]
    nsepi = [1.8766e+20, 1.9678e+20, 1.8622e+20, 1.9460e+20, 1.8414e+20, 1.9151e+20, 1.8719e+20, 1.8340e+20, 2.0142e+20]
    nsepo = [1.6306e+20, 1.6842e+20, 1.6362e+20, 1.6818e+20, 1.5927e+20, 1.6253e+20, 1.5807e+20, 1.6084e+20, 1.6837e+20]
    plt.figure()
    Dos = list(zip(*Ds))[0]
    Dis = list(zip(*Ds))[1]
    plt.scatter(Dos, Dis, c=nsepo, edgecolors='black', s=200)
    plt.xlabel('Dof')
    plt.ylabel('Dif')
    plt.colorbar()
    ax = plt.gca()
    ax.set_facecolor('lightgray');
    plt.show()
    
  
    

# Code below will not run if runcase is imported
if __name__ == '__main__':
    ## Restore default base case
    # setup(pcore=5e6, impurity=None, impFrac=None, experimental=False)
    # restore('base_pcore5MW')
    # plotall()
    
    ## Adjust D such that n_LCFS = 1e20
    # adjustD()
    # rundt(1e-10)
    
    ## Start adding impurities
    # setup(pcore=1e6)
    # restore('base')
    # setup(pcore=1e6, impFrac=.0001)
    # bbb.dtreal=1e-10; bbb.ftol=1e10; bbb.exmain()
    # rundt(1e-10, ii1max=5000)
    
    ## Restore case with impurities
    setup()
    restore('cmod-nonorthog')
    # for i in np.arange(2.51, 2.7, .05):
    #     setup(pcore=i*1e6, impFrac=.001)
    #     rundt()
    #     hdf5_save(f'{i}MW_impFrac.001.h5')
    
    # Restore convective case
    #setup(pcore=5e6)
    #restore('base_pcore5MW')
    #setup(pcore=5e6, convec=True)
    #restore('5MW_vconv')
    
    # postprocess()
