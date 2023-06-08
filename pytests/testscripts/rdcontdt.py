# This file runs a time-dependent case using dtreal.  First, obtain a converged
# solution for a (usually small) dtreal; xuedge must report iterm=1 at the end.
# Then adjust control parameters in rdinitdt; read this file, which reads rdinitdt.
# If a mistake is made, to restart this file without a Jacobian evaluation,
# be sure to reset iterm=1 (=> last step was successful)

# IMPORT UEDGE (assuming starting from ipython before any imports)
from uedge import *
from numpy import zeros,sqrt
        
# IMPORT HDF5 routines for saving solutions below
from uedge.hdf5 import *

i_stor = 0
nfe_tot = 0
savefn = "savedt.hdf5"  # name of hdf5 savefile written every timestep

no = 0;yes = 1
echo = no

# Set precisions of floating point output
###import print_options
###print_options.set_float_precision(4)

nx=com.nx;ny=com.ny;nisp=com.nisp;ngsp=com.ngsp;numvar=bbb.numvar
isteon=bbb.isteon
if (i_stor==0):
   ni_stor = zeros((bbb.n_stor,nx+1+1,ny+1+1,nisp),"d")	# set time storage arrays
   up_stor = zeros((bbb.n_stor,nx+1+1,ny+1+1,nisp),"d")
   te_stor = zeros((bbb.n_stor,nx+1+1,ny+1+1),"d")
   ti_stor = zeros((bbb.n_stor,nx+1+1,ny+1+1),"d")
   ng_stor = zeros((bbb.n_stor,nx+1+1,ny+1+1,ngsp),"d")
   phi_stor = zeros((bbb.n_stor,nx+1+1,ny+1+1),"d")
   tim_stor = zeros((bbb.n_stor),"d")
   dtreal_stor = zeros((bbb.n_stor),"d")
   nfe_stor = zeros((bbb.n_stor),"l")
   dt_stor = (bbb.tstor_e - bbb.tstor_s)/(bbb.n_stor - 1)

i_stor = max(i_stor,1)			# set counter for storage arrays
bbb.dt_tot = max(bbb.dt_tot,0.)
nfe_tot = max(nfe_tot,0)
deldt_0 = bbb.deldt
isdtsf_sav = bbb.isdtsfscal

if (bbb.ipt==1 and bbb.isteon==1): 	# set ipt to te(nx,iysptrx+1) if no user value
   ipt = bbb.idxte[nx-1,com.iysptrx]  #note: ipt is local, bbb.ipt global

bbb.irev = -1         # forces second branch of irev in ii1 loop below
if (bbb.iterm == 1):  # successful initial run with dtreal
   bbb.dtreal = bbb.dtreal/bbb.mult_dt     # gives same dtreal after irev loop
else:                 # unsuccessful initial run; reduce dtreal
   bbb.dtreal = bbb.dtreal/(3*bbb.mult_dt) # causes dt=dt/mult_dt after irev loop
   
if (bbb.initjac == 0): bbb.newgeo=0
dtreal_sav = bbb.dtreal
bbb.itermx = bbb.itermxrdc
bbb.dtreal = bbb.dtreal/bbb.mult_dt	#adjust for mult. to follow; mult_dt in rdinitdt
bbb.dtphi = bbb.rdtphidtr*bbb.dtreal
neq=bbb.neq
svrpkg=bbb.svrpkg.tostring().strip()
#
bbb.ylodt = bbb.yl
bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
fnrm_old = sqrt(sum((bbb.yldot[0:neq]*bbb.sfscal[0:neq])**2))
if (bbb.initjac == 1): fnrm_old=1.e20
print(( "initial fnrm =",fnrm_old))

for ii1 in range( 1, bbb.ii1max+1):
   if (bbb.ismfnkauto==1): bbb.mfnksol = 3
   # adjust the time-step
   if (bbb.irev == 0):
      # Only used after a dt reduc. success. completes loop ii2 for fixed dt
      bbb.dtreal = min(3*bbb.dtreal,bbb.t_stop)	#first move forward after reduction
      bbb.dtphi = bbb.rdtphidtr*bbb.dtreal
      if (bbb.ismfnkauto==1 and bbb.dtreal > bbb.dtmfnk3): bbb.mfnksol = -3
      bbb.deldt =  3*bbb.deldt
   else:
      # either increase or decrease dtreal; depends on mult_dt
      bbb.dtreal = min(bbb.mult_dt*bbb.dtreal,bbb.t_stop)
      bbb.dtphi = bbb.rdtphidtr*bbb.dtreal
      if (bbb.ismfnkauto==1 and bbb.dtreal > bbb.dtmfnk3): bbb.mfnksol = -3
      bbb.deldt =  bbb.mult_dt*bbb.deldt
      
   bbb.dtreal = min(bbb.dtreal,bbb.dt_max)
   bbb.dtphi = bbb.rdtphidtr*bbb.dtreal
   if (bbb.ismfnkauto==1 and bbb.dtreal > bbb.dtmfnk3): bbb.mfnksol = -3
   bbb.deldt = min(bbb.deldt,deldt_0)
   bbb.deldt = max(bbb.deldt,bbb.deldt_min)
   nsteps_nk=1
   print('--------------------------------------------------------------------')
   print('--------------------------------------------------------------------')
   print(' ')
   print(('*** Number time-step changes = ',ii1,' New time-step = ', bbb.dtreal))
   print('--------------------------------------------------------------------')

   bbb.itermx = bbb.itermxrdc
   if (ii1>1  or  bbb.initjac==1):	# first time calc Jac if initjac=1
      if (bbb.irev == 1):      # decrease in bbb.dtreal
         if (bbb.numrev < bbb.numrevjmax and \
            bbb.numrfcum < bbb.numrevjmax+bbb.numfwdjmax): #dont recom bbb.jac
            bbb.icntnunk = 1	
            bbb.numrfcum = bbb.numrfcum + 1
         else:                          # force bbb.jac calc, reset numrev
            bbb.icntnunk = 0
            bbb.numrev = -1		      # yields api.zero in next statement
            bbb.numrfcum = 0
         bbb.numrev = bbb.numrev + 1
         bbb.numfwd = 0
      else:  # increase in bbb.dtreal
         if (bbb.numfwd < bbb.numfwdjmax and \
            bbb.numrfcum < bbb.numrevjmax+bbb.numfwdjmax): 	#dont recomp bbb.jac
            bbb.icntnunk = 1
            bbb.numrfcum = bbb.numrfcum + 1
         else:
            bbb.icntnunk = 0			#recompute jacobian for increase dt
            bbb.numfwd = -1
            bbb.numrfcum = 0
         bbb.numfwd = bbb.numfwd + 1
         bbb.numrev = 0			#bbb.restart counter for dt reversals
      bbb.isdtsfscal = isdtsf_sav
      bbb.ftol = min(bbb.ftol_dt, 0.01*fnrm_old)
      bbb.exmain() # take a single step at the present bbb.dtreal
      if (bbb.iterm == 1):
         bbb.dt_tot = bbb.dt_tot + bbb.dtreal
         nfe_tot = nfe_tot + bbb.nfe[0,0]
         bbb.ylodt = bbb.yl
         bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
         fnrm_old = sqrt(sum((bbb.yldot[0:neq-1]*bbb.sfscal[0:neq-1])**2))
         if (bbb.dt_tot>=0.9999999*bbb.t_stop  or  fnrm_old<bbb.ftol_min):
            print(' ')
            print('*****************************************************')
            print('**  SUCCESS: frnm < bbb.ftol; or dt_tot >= t_stop  **')
            print('*****************************************************')
            break

   bbb.icntnunk = 1
   bbb.isdtsfscal = 0
   for ii2 in range( 1, bbb.ii2max+1): #take ii2max steps at the present time-step
      if (bbb.iterm == 1):
         bbb.itermx = bbb.itermxrdc
         bbb.ftol = min(bbb.ftol_dt, 0.01*fnrm_old)
         bbb.exmain()
         if (bbb.iterm == 1):
            bbb.ylodt = bbb.yl
            bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
            fnrm_old = sqrt(sum((bbb.yldot[0:neq-1]*bbb.sfscal[0:neq-1])**2))
            print("Total time = ",bbb.dt_tot,"; Timestep = ",bbb.dtreal)
            print("variable index ipt = ",ipt, " bbb.yl[ipt] = ",bbb.yl[ipt])
            dtreal_sav = bbb.dtreal
            bbb.dt_tot = bbb.dt_tot + bbb.dtreal
            nfe_tot = nfe_tot + bbb.nfe[0,0]
            if (bbb.dt_tot>=0.999999999999*bbb.t_stop  or  fnrm_old<bbb.ftol_min):
               break
##       Store variables if a storage time has been crossed
            if (bbb.dt_tot >= dt_stor*i_stor and i_stor<=bbb.n_stor):
               i_stor1 = i_stor-1
               ni_stor[i_stor1,:,:,:] = ni
               up_stor[i_stor1,:,:,:] = up
               te_stor[i_stor1,:,:] = te
               ti_stor1[i_stor1,:,:] = ti
               ng_stor[i_stor1,:,:,:] = ng
               phi_stor1[i_stor1,:,:] = phi
               tim_stor[i_stor1] = bbb.dt_tot
               nfe_stor[i_stor1] = nfe_tot
               dtreal_stor[i_stor1] = bbb.dtreal
               i_stor = i_stor + 1
   ##          End of storage section
      
   if (bbb.dt_tot>=bbb.t_stop  or  fnrm_old<bbb.ftol_min): break   # need for both loops
   bbb.irev = bbb.irev-1
   if (bbb.iterm != 1):	#print bad eqn, cut dtreal by 3, set irev flag
            ####### a copy of idtroub script ########################
      oldecho=echo
      echo=no
      # integer ii
      # real8 ydmax 
      scalfac = bbb.sfscal
      if (svrpkg != "nksol"): scalfac = 1/(bbb.yl + 1.e-30)  # for time-dep calc.
      ydmax = 0.999999999*max(abs(bbb.yldot*scalfac))
      itrouble = 0
      for ii in range(neq):
         if (abs(bbb.yldot[ii]*scalfac[ii]) > ydmax):
            itrouble=ii
            print("** Fortran index of trouble making equation is:")
            print(itrouble+1)
            break
      print("** Number of variables is:")
      print("numvar = ", numvar)
      print(" ")
      iv_t = (itrouble).__mod__(numvar) + 1
      print("** Troublemaker equation is:")
      print("iv_t = ",iv_t)
      print(" ")
      print("** Troublemaker cell (ix,iy) is:")
      print(bbb.igyl[itrouble,])
      print(" ")
      print("** Timestep for troublemaker equation:")
      print(bbb.dtuse[itrouble])
      print(" ")
      print("** yl for troublemaker equation:")
      print(bbb.yl[itrouble])
      print(" ")
      echo=oldecho
      ######## end of idtroub script ##############################

      if (bbb.dtreal < bbb.dt_kill):
          print(' ')
          print('*************************************')
          print('**  FAILURE: time-step < dt_kill   **')
          print('*************************************')
          break
      bbb.irev = 1
      print('*** Converg. fails for bbb.dtreal; reduce time-step by 3, try again')
      print('----------------------------------------------------------------- ')
      bbb.dtreal = bbb.dtreal/(3*bbb.mult_dt)
      bbb.dtphi = bbb.rdtphidtr*bbb.dtreal
      if (bbb.ismfnkauto==1 and bbb.dtreal > bbb.dtmfnk3): bbb.mfnksol = -3
      bbb.deldt =  bbb.deldt/(3*bbb.mult_dt) 
      bbb.iterm = 1
echo = yes
