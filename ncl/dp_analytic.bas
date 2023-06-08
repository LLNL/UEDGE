#**************************************************************
function dp_analytic(rmaj,rmin,qsafe,bfil,ni,te)
#
#-Analytic formula for neoclassical particle diff. coeff.
#
# Inputs: 
#         Te, Ti in eV
#         Ni in 1/cm**3
#         R,a in cm
#         B in Gauss
#
#-MVU, 21-Oct-2004
#------------------------------------------------------------

  real coulog=10.
  real eps=rmin/rmaj #inverse aspect ratio
  real vte=4.19e7*sqrt(te) #cm/sec
  real wce=1.76e7*bfil #rad/sec

  real nue=2.91e-6*ni*coulog*(te**(-1.5)) # 1/sec, coll. freq.
  real nub=vte/abs(qsafe*rmaj) #1/s, bounce freq.

  real dperp=(vte**2/wce**2)*nue
  real dnc, nieff, nueff

   
   if (nue > nub) then #PS regime
     << "PS"
     dnc=dperp*(qsafe**2)    
   endif  

   if (nue <= nub .AND. nue >= nub*(eps**1.5)) then #plateau
     << "plateau"
     nieff=te**2*(4.19e7/2.91e-6)/(coulog*qsafe*rmaj)
     nueff=2.91e-6*nieff*coulog*(te**(-1.5))    
     dperp=(vte**2/wce**2)*nueff
     dnc=dperp*(qsafe**2)    
   endif  

   if (nue < nub*(eps**1.5)) then #banana
     << "banana"
     dnc=dperp*(qsafe**2)/(eps**1.5)    
   endif  

   ##<< "nue=" << nue << ", nub=" << nub << ", nustar=" << nue/nub

#
#
return dnc
endf
