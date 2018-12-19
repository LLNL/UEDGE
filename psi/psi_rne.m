c/*
c
c=======================================================
c       Plasma Surface Interaction (PSI) package
c  for simulation of plasma-wall sputtering processes.
c
c                 PSI Version 1.2   
c              Date: March 6, 2007
c
c
c Copyright 2007 by A.Yu. Pigarov and R.D. Smirnov
c
c
c                   Developers:
c          A.Yu. Pigarov and R.D. Smirnov
c
c
c         E-mails: apigarov@ucsd.edu
c                  rsmirnov@ucsd.edu
c
c    PSI package has been developed as part of
c   the Dust Transport in Tokamaks (DUSTT) code:
c         
c [1] A.Yu.Pigarov, S.I. Krasheninnikov, T.K. Soboleva,
c     and T.D. Rognlien "Dust-particle transport in
c     tokamak edge plasmas", Phys. Plasmas 12 (2005), 122508.
c [2] R.D. Smirnov, A.Yu. Pigarov, M. Rosenberg,
c     S.I. Krasheninnikov, and D.A. Mendis "Modelling of
c     dynamics and transport of carbon dust particles in 
c     tokamaks", Plasma Phys. Control. Fusion 49 (2007), 347-371.
c
c
c    The PSI package is distributed in the hope 
c that it will be useful to computational scientists in field
c of fusion plasma and other plasma related studies.
c
c    The PSI package is FREE SOFTWARE. 
c You can use, copy, and modify this software for any purpose
c and without fee provided that the above copyright
c notice appear in all copies.
c
c    The PSI package is distributed "AS IS", 
c i.e. without any warranty including all implied warranties
c of merchantability and fitness.
c   
c=======================================================
c
c*/
      subroutine psirne(e,amu, rn,re)
      implicit none
      
      double precision e,amu,rn,re
      double precision rmimi,rmama
      data rmimi/1.d-6/
      data rmama/0.99999999d0/


      call psijrne(e,amu, rn,re)
      
      if(rn > rmama) then
      rn=rmama
      endif
      if(rn < rmimi) then
      rn=rmimi
      endif

      if(re > rmama) then
      re=rmama
      endif
      if(re < rmimi) then
      re=rmimi
      endif

      end


      double precision function psi_rn(e,amu)
      implicit none
      double precision e,amu

      double precision rn,re

      call psirne(e,amu,rn,re)

      psi_rn=rn
      return

      end


      double precision function psi_re(e,amu)
      implicit none
      double precision e,amu

      double precision rn,re

      call psirne(e,amu,rn,re)

      psi_re=re
      return

      end


      subroutine psirnep(izpt,ampt,izsf,amsf)
      implicit none
      integer izpt,izsf
      double precision ampt,amsf
c---------------------------------------------------------
c  prepare the constants for ippj reflection coefficients
c----------------------------------------------------------
      common/psirnpar/epsC,ca,cs,ct,cbt,cg,csb,cal,cf
      double precision epsC,ca,cs,ct,cbt,cg,csb,cal,cf


      double precision psicel
      double precision fal,fmu,t,z1s1,z1s,z2s1,z2s,am1,s,gcor,r
      double precision zpt,zsf


      zpt = dble(izpt) 
      zsf = dble(izsf)
      fmu = amsf/ampt 
      fal = 1.d0 + 1.d0 / fmu 
      t   = dlog(zpt) 
      z1s1= dexp( t/3.d0 ) 
      z1s = z1s1*z1s1 
      t   = dlog(zsf) 
      z2s1= dexp( t/3.d0 ) 
      z2s = z2s1*z2s1 
      t   = dsqrt(1.d0+z1s/z2s) 
      epsC=0.032534/(zpt*zsf*z2s*fal*t) 
      am1 = dsqrt(ampt) 
      cf  = z1s/am1 
      cal = dlog(fal) 
      ca  = 0.0793d0 * z1s*fmu/am1 
      s   = (1.d0 + z2s) * fal/t 
      cs  = 0.2617d0 * s*dsqrt(s)
      gcor= psicel(izsf) 
      cs  = cs * ampt/zsf * gcor 
      csb = 61.47d0 * zpt*fal*fmu*t/z2s * 4.d0 
      if(izpt < 3) then
      cg = 100.d0 * zpt/zsf 
      else
      cg = 5.d0
      endif 
      if(izsf < 13) then
       r =  12.d0 +  7.d0 / zsf 
      else
       r = dlog(zsf) 
       r = r * (-1.19d0) 
       r =  9.76d0 +  58.5d0 * dexp(r)
      endif
      ct  =  1.d-9 / (ampt*z2s1) 
      cbt =  1.022d6 /(zsf*r) 
      csb = csb/r

      end



      double precision function psirnjf5(e)
      implicit none
      double precision e
 
      common/psirnpar/epsC,ca,cs,ct,cbt,cg,csb,cal,cf
      double precision epsC,ca,cs,ct,cbt,cg,csb,cal,cf


      double precision bet,t,eb,sb


      t   = ct * e 
      bet =  1.d0 / (  1.d0 + t ) 
      bet = (  2.d0 + t) * t * bet * bet 
      eb  = bet*cbt


      sb =  1.d0 + cg/eb + eb/(  1.d0 - bet) 
      sb = (dlog(sb) - bet ) / eb 
      sb = sb*csb


      psirnjf5=sb 
      return
 
      end





      double precision function psirnjft(eps)
      implicit none
      double  precision eps
c
c------------------------------------------
c   reverced reduced range of ions
c   for nuclear stopping only
c    w.d.wilson et al  pr b15 2458(1977)
c------------------------------------------
c
      double precision psiegamm
      double precision xx,x,x1,x2,s1,s2,s
      double precision c1,c2,c3,c4
 
      data c1 /-0.3732d0/ 
      data c2 /1.1776d0/ 
      data c3 /-2.d0/
      data c4 /0.56258d0/
 
      xx = eps 
      x  = dlog(xx*c2) 
      x1 = c1*x 
      x2 = c3*x 
      s1 = psiegamm(x1) 
      s2 = psiegamm(x2) 
       s = (s1-s2) / (c4*c2) 
       s =  2.d0 * xx / s


      psirnjft=s 
      return

      end




      double precision function psirnjfs(e)
      implicit none
      double precision e
c
c--------------------------------
c   correction factor of ippj
c--------------------------------
c
      common/psirnpar/epsC,ca,cs,ct,cbt,cg,csb,cal,cf
      double precision epsC,ca,cs,ct,cbt,cg,csb,cal,cf


      double precision psirnjf5,psirnjft 
      double precision eps,se,sa,sl,sb,ss,pe,f,su


      eps = epsC*e 
      se = dsqrt(eps) 
      sa = ca*se 
      sl = cs*sa 
      sb = psirnjf5(e) 
      ss = sb * sl/(sl+sb) 
      pe = psirnjft(eps) 
      f  = (pe+ss)/sa 
      su = (2.d0 * eps - 3.d0)/(eps + 1.d0) 
      f = f * dexp(su*cal) 
      f = f * cf


      psirnjfs=f
      return 


      end




      double precision function psirnjun(ep,am,qtn)
      implicit none
      double precision ep,am,qtn
c
c-------------------------------------------------
c   angular-dependence of particle reflection
c	   coefficients f(e,am)
c	 usage: r=r+(1-r)*f(e,mu)
c-------------------------------------------------
c
      common/psirnpar/epsC,ca,cs,ct,cbt,cg,csb,cal,cf
      double precision epsC,ca,cs,ct,cbt,cg,csb,cal,cf


      double precision ae,cot,amu,eps,act,sa1,sa2,s
      double precision amui,amua,eps00,cn1,cn2,cn3,cn4
 
      data amui  /1.d-8/ 
      data amua  /0.99999999d0/ 
      data eps00 /1.d-18/
      data cn1   /7.38d0/
      data cn2   /0.359d0/
      data cn3   /0.836d0/
      data cn4   /-0.087d0/


      eps = ep 
      amu = am


      if(eps < eps00) then
	psirnjun=0.d0
	return 
	endif


      if(amu < amui) then
	psirnjun=1.d0
	return 
	endif


      if(amu > amua) then
	psirnjun=0.d0
      return 
	endif


 	cot =  1.d0 - amu*amu 
      if(cot < 0.d0) cot = -cot 
	cot = amu / dsqrt(cot) 
      if(cot < 0.d0) cot = -cot 
      ae  = dlog(eps) 
      act = dlog(cot) 
      sa1 = cn1 * qtn * dexp(cn2*ae) 
      sa2 = cn3 * dexp(cn4*ae) 
      s   = sa1 * dexp(2.d0*sa2*act)
      s   =  1.0 / (1.0 + s) 
      if(s > 1.0)  s = 1.d0


      psirnjun=s
      return 


      end





      double precision function psirnjue(ep,am,qte)
      implicit none
      double precision ep,am,qte
c
c-------------------------------------------------
c   angular-dependence of energy reflection
c	   coefficients f(e,am)
c	 usage: r=r+(1-r)*f(e,mu)
c-------------------------------------------------
c
      common/psirnpar/epsC,ca,cs,ct,cbt,cg,csb,cal,cf
      double precision epsC,ca,cs,ct,cbt,cg,csb,cal,cf


      double precision ae,cot,amu,eps,act,sa1,sa2,s
      double precision amui,amua,eps00,ce1,ce2,ce3,ce4


      data amui  /1.d-4/ 
      data amua  /0.9999d0/ 
      data eps00 /1.d-12/ 
      data ce1   /17.9d0/ 
      data ce2   /0.453d0/ 
      data ce3   /0.771d0/ 
      data ce4   /-0.014d0/


      eps = ep 
      amu = am


      if(eps < eps00) then
	psirnjue=0.d0
	return 
	endif


      if(amu < amui) then
	psirnjue=1.d0
	return 
	endif


      if(amu > amua) then
	psirnjue=0.d0
      return 
	endif


 	cot =  1.d0 - amu*amu 
      if(cot < 0.d0) cot = -cot 
	cot = amu / dsqrt(cot) 
      if(cot < 0.d0) cot = -cot 
      ae  = dlog(eps) 
      act = dlog(cot) 
      sa1 = ce1 * qte * dexp(ce2*ae) 
      sa2 = ce3 * dexp(ce4*ae) 
      s   = sa1 * dexp(2.d0*sa2*act)
      s   =  1.d0 / (1.d0 + s) 
      if(s > 1.d0 )  s=1.d0


      psirnjue=s
      return
 
      end




      double precision function psirnjfe(x)
      implicit none
      double precision x
c
c-------------------------------------------
c   scalling for reduced energy-reflection
c		coefficient
c-------------------------------------------
c
      double precision s1,s2,s
      double precision c1,c2,c3,c4,c5
 
      data c1 /0.047d0/ 
      data c2 /0.597d0/ 
      data c3 /0.619d0/ 
      data c4 /1.5d0/ 
      data c5 /0.705d0/
 
      s1 = dlog(x/c1) 
      s2 = dlog(x/c3) 
      s  = dexp(c2*s1) 
      s  = s + dexp(c4*s2) 
      s  = c5/(1.d0 + s)


      psirnjfe=s 
      return


      end




      double precision function psirnjfd(x)
      implicit none
      double precision x
c
c--------------------------------------------
c   average energy of reflected particles
c--------------------------------------------
c
      double precision s1,s2,s
      double precision c1,c2,c3,c4,c5,cm
 
      data c1   /0.133d0/ 
      data c2   /0.285d0/
      data c3   /0.53d0/
      data c4   /85.0d0/
      data c5   /-1.46d0/
      data cm   /0.99999999999d0/


      s  = dlog(x/c1) 
      s1 =  1.d0 + dexp(c2*s) 
      s  = dlog(x/c4) 
      s2 =  1.d0 + dexp(c5*s) 
      s  =  1.d0/s1 + c3/s2 
      if( s > cm ) s=cm


      psirnjfd=s 
      return 


      end




      subroutine psijrne(e,amu, rn,re)
      implicit none
      double precision e,amu,rn,re
c
c------------------------------------------------
c   particle & energy reflection coefficients
c  r.ito et. al., ippj-am-41 nagoya japan 1985
c	( call psi_ripp_prp before )
c------------------------------------------------
c
      common/psirnpar/epsC,ca,cs,ct,cbt,cg,csb,cal,cf
      double precision epsC,ca,cs,ct,cbt,cg,csb,cal,cf


      double precision psirnjfe,psirnjfd,psirnjfs,
     ,                 psirnjun,psirnjue
      double precision ep,eps,fun,fue,rfe,f,rnn,ree,sn,se
      double precision fcor,emax,emin
c      double precision feps

      data fcor /0.999999d0/
c      data feps /1.d-15/
      data emax /1.0d6/
      data emin /1.0d-6/


 	    ep = e 
	if( ep > emax ) ep = emax
      if( ep < emin ) ep = emin


      eps = epsC * ep 
      ree = psirnjfe(eps) 
      rfe = psirnjfd(eps) 
      f   = psirnjfs(ep)
 
      ree = ree/f 
      rnn = ree/rfe
 
      sn =  1.d0 - rnn 
      se =  1.d0 - ree 
      fun = psirnjun(eps,amu, 1.d0) 
      fue = psirnjue(eps,amu, 1.d0) 
 	sn = sn * fun 
	se = se * fue


      f = sn*fcor 
      if(se > f) then
      se = f
      endif


      rn = rnn + sn 
      re = ree + se


      end

