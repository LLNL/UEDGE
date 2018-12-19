# This version hand-processed following bas2py conversion
# This file runs a time-dependent case using dtreal.  First, read rdinitdt
# to set variables, and: read this file.
# If a mistake is made, to restart this file without a Jacobian evaluation,
# be sure to reset iterm=1 (=> last step was successful)

from uedge import *
import PW
          
# INITIALIZE PARAMS -- SHOULD BE DONE IN MASTER SCRIPT OR TERMINAL SESSION
# BEFORE INVOKING THIS SCRIPT
#execfile("rdinitdt.py")
no = 0;yes=1
echo = no

# allocate space -- SHOULD BE DONE IN MASTER SCRIPT OR TERMINAL SSESSION
# BEFORE INVOKING THIS SCRIPT
# bbb.allocate()

nx=com.nx;ny=com.ny;nisp=com.nisp;ngsp=com.ngsp;numvar=bbb.numvar
isteon=bbb.isteon
if (i_stor==0):
   ni_stor = zeros((n_stor,nx+1+1,ny+1+1,nisp),"d")	#set storage arrays
   up_stor = zeros((n_stor,nx+1+1,ny+1+1,nisp),"d")
   te_stor = zeros((n_stor,nx+1+1,ny+1+1),"d")
   ti_stor = zeros((n_stor,nx+1+1,ny+1+1),"d")
   ng_stor = zeros((n_stor,nx+1+1,ny+1+1,ngsp),"d")
   phi_stor = zeros((n_stor,nx+1+1,ny+1+1),"d")
   tim_stor = zeros((n_stor),"d")
   dtreal_stor = zeros((n_stor),"d")
   nfe_stor = zeros((n_stor),"l")
   dt_stor = (tstor_e - tstor_s)/(n_stor - 1)
   # character*80 s(1)
   sfn ="fs=PW.PW('pfdt_"+savefname+"' ,'w')"
   # character*80 s2(1)
   if (n_stor>0): s2 ="fs2=PW.PW(pftstor_"+savefname+" ,'w')"
i_stor = max(i_stor,1)			# set counter for storage arrays
dt_tot = max(dt_tot,0.)
nfe_tot = max(nfe_tot,0)
deldt_0 = bbb.deldt
# real8 fnrm_old
isdtsf_sav = bbb.isdtsfscal
if (ipt==1 and isteon==1): 	#reset ipt to te(nx,iysptrx+1) if no user value
   ipt = bbb.idxte[nx-1,com.iysptrx]
#
##exmain
irev = -1         # forces second branch of irev loop below
if (bbb.iterm == 1):  # successful initial run with dtreal
   bbb.dtreal = bbb.dtreal/mult_dt     # causes same dt after irev loop
else:             # unsuccessful initial run with dtreal; reduce dtreal
   bbb.dtreal = bbb.dtreal/(3*mult_dt) # causes dt=dt/mult_dt after irev loop
   
if (initjac == 0): bbb.newgeo=0
dtreal_sav = bbb.dtreal
bbb.itermx = 10
##ftol = ftol_dt 
##bbb.dtreal = bbb.dtreal/mult_dt		#adjust dtreal for multiplic. to follow
bbb.dtphi = rdtphidtr*bbb.dtreal
irev = -1
# allocate space -- is this really where we should do it?
#bbb.allocate()
#introduce some shortcut na mes
##yl = bbb.yl
##yldot=bbb.yldot
neq=bbb.neq
svrpkg=bbb.svrpkg.tostring().strip()
#
bbb.ylodt = bbb.yl
bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
fnrm_old = sqrt(sum((bbb.yldot[0:neq]*bbb.sfscal[0:neq])**2))
if (initjac == 1): fnrm_old=1.e20
print "initial fnrm =",fnrm_old

for ii1 in range( 1, ii1max+1):
   if (ismfnkauto==1): bbb.mfnksol = 3
   #set the time-step
   if (irev == 0):
      # Only used after a dt reduc. success. completes loop ii2 for fixed dt
      bbb.dtreal = min(3*bbb.dtreal,t_stop)	#first move forward after reduction
      bbb.dtphi = rdtphidtr*bbb.dtreal
      if (ismfnkauto==1 and bbb.dtreal > dtmfnk3): bbb.mfnksol = -3
      bbb.deldt =  3*bbb.deldt
      print ' increase dtreal to ',bbb.dtreal
   else:
      # Most common branch for either increasing to decreasing dtreal
      bbb.dtreal = min(mult_dt*bbb.dtreal,t_stop)
      bbb.dtphi = rdtphidtr*bbb.dtreal
      if (ismfnkauto==1 and bbb.dtreal > dtmfnk3): bbb.mfnksol = -3
      bbb.deldt =  mult_dt*bbb.deldt
      print ' decrease dtreal to ',bbb.dtreal,' by 1/mult_dt = ',1./mult_dt
   bbb.dtreal = min(bbb.dtreal,dt_max)
   bbb.dtphi = rdtphidtr*bbb.dtreal
   if (ismfnkauto==1 and bbb.dtreal > dtmfnk3): bbb.mfnksol = -3
   bbb.deldt = min(bbb.deldt,deldt_0)
   bbb.deldt = max(bbb.deldt,deldt_min)
   nsteps_nk=1
   print ' '
   print '*** Number time-step changes and dt: idtc; bbb.dtreal = ', ii1,';', bbb.dtreal
   print '------------------------------------------------------------------------------'
   ##   ftol = ftol_dt/(1 + 100*bbb.dtreal/dt_max) # gives better converg at large dt
   bbb.itermx = 10
   if (ii1>1  or  initjac==1):	# first time calc Jac if initjac=1
      bbb.icntnunk = 0
      bbb.isdtsfscal = isdtsf_sav
      bbb.ftol = min(ftol_dt, 0.01*fnrm_old)
      exmain() # take a single step at the present bbb.dtreal
      if (bbb.iterm == 1):
         dt_tot = dt_tot + bbb.dtreal
         nfe_tot = nfe_tot + bbb.nfe[0,0]
         bbb.ylodt = bbb.yl
         bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
         fnrm_old = sqrt(sum((bbb.yldot[0:neq-1]*bbb.sfscal[0:neq-1])**2))
         if (dt_tot>=0.999999999999*t_stop  or  fnrm_old<ftol_min):
            print ' total time this advance: dt_tot = ', dt_tot
            print ' '
            print '*****************************************************'
            print '**  SUCCESS: frnm < bbb.ftol; or dt_tot >= t_stop  **'
            print '*****************************************************'
            break

   #bbb.icntnunk=1 #turn off the automatic recalculation of initial Jacobians for subsequent steps
   bbb.isdtsfscal = 0
   for ii2 in range( 1, ii2max+1): #take ii2max steps at the present time-step
      if (bbb.iterm == 1):
         bbb.itermx = 10
         bbb.ftol = min(ftol_dt, 0.01*fnrm_old)
         bbb.exmain()
         if (bbb.iterm == 1):
            bbb.ylodt = bbb.yl
            bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
            fnrm_old = sqrt(sum((bbb.yldot[0:neq-1]*bbb.sfscal[0:neq-1])**2))
            print "* Step number at this dtreal:",ii2,"; initial fnrm =",fnrm_old
            print 'Number of func eval',bbb.nfe[0,0],'; nonlinear iters',bbb.nni[0,0],'; preJac eval = ',bbb.npe[0,0]
            dtreal_sav = bbb.dtreal
            dt_tot = dt_tot + bbb.dtreal
            nfe_tot = nfe_tot + bbb.nfe[0,0]
            if (dt_tot>=0.999999999999*t_stop  or  fnrm_old<ftol_min):
               print ' total time this advance: dt_tot = ', dt_tot
               print ' '
               print '*****************************************************'
               print '**  SUCCESS: frnm < bbb.ftol; or dt_tot >= t_stop  **'
               print '*****************************************************'
               break
            print "Total time = ",dt_tot,"; Timestep = ",bbb.dtreal,"; bbb.yl[ipt], ipt = ",bbb.yl[ipt],", ",ipt
            print " "
##       Store variables if a storage time has been crossed
            if (dt_tot >= dt_stor*i_stor and i_stor<=n_stor):
               i_stor1 = i_stor-1
               ni_stor[i_stor1,:,:,:] = ni
               up_stor[i_stor1,:,:,:] = up
               te_stor[i_stor1,:,:] = te
               ti_stor1[i_stor1,:,:] = ti
               ng_stor[i_stor1,:,:,:] = ng
               phi_stor1[i_stor1,:,:] = phi
               tim_stor[i_stor1] = dt_tot
               nfe_stor[i_stor1] = nfe_tot
               dtreal_stor[i_stor1] = bbb.dtreal
               i_stor = i_stor + 1
   ##          End of storage section
## Write out the time-dependent arrays in outer loop if i_stor > 1
                  
   exec(sfn) #see file rdinitdt - this creates pfdt_"+savefname+"
   fs.nis=bbb.nis
   fs.ups=bbb.ups
   fs.tes=bbb.tes
   fs.tis=bbb.tis
   fs.ngs=bbb.ngs
   fs.phis=bbb.phis
   fs.dtreal_sav=dtreal_sav
   fs.dt_tot=dt_tot
   fs.ii1=ii1
   fs.close()
   if (i_stor > 1):
      exec(s2)
      fs2.ni_stor=ni_stor;fs2.up_stor=up_stor;fs2.te_stor=te_stor
      fs2.ti_stor=ti_stor;fs2.ng_stor=ng_stor;fs2.phi_stor=phi_stor
      fs2.tim_stor=tim_stor
      fs2.close()
      
   if (dt_tot>=t_stop  or  fnrm_old<ftol_min): break   # need for both loops
   irev = irev-1
   if (bbb.iterm != 1):	#print bad eqn, cut dtreal by 3, set irev flag
            ####### just a copy of idtroub script ########################
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
            print "** Fortran index of trouble making equation is:"
            print itrouble+1
            break
      print "** Number of variables is:"
      print "numvar = ", numvar
      print " "
      iv_t = (itrouble).__mod__(numvar) + 1
      print "** Troublemaker equation is:"
      print "iv_t = ",iv_t
      print " "
      print "** Troublemaker cell (ix,iy) is:"
      print bbb.igyl[itrouble,]
      print " "
      print "** Timestep for troublemaker equation:"
      print bbb.dtuse[itrouble]
      print " "
      print "** yl for troublemaker equation:"
      print bbb.yl[itrouble]
      print " "
      echo=oldecho
      ######## end of idtroub script ##############################

      irev = 1
      print '*** Converg. fails for bbb.dtreal; reduce time-step by 3, try again'
      print '----------------------------------------------------------------- '
      bbb.dtreal = bbb.dtreal/(3*mult_dt)
      bbb.dtphi = rdtphidtr*bbb.dtreal
      if (ismfnkauto==1 and bbb.dtreal > dtmfnk3): bbb.mfnksol = -3
      bbb.deldt =  bbb.deldt/(3*mult_dt) 
      bbb.iterm = 1
echo = yes
