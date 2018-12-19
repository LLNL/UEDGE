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
      integer function psi_yippu_prp(izp,mp,izt,mt, uth,CLight)
      implicit none
      integer izp,izt
      double precision mp,mt
      double precision uth     ! threshold energy 
      double precision CLight  ! critical ratio (mp/mt)

      common/psiyippui/Light
      integer    Light
      common/psiyippud/
     ,      pi0,
     ,      L_fb,L_gb1,L_gb2,
     ,      H_fb,H_gb_,H_gb1,H_gb2,H_gb3,
     ,      CY_,CKe_,
     ,      El_,Eth_t,
     ,      Z1_,Z2_,M1_,M2_
      double precision pi0
      double precision L_fb,L_gb1,L_gb2
      double precision H_fb,H_gb_,H_gb1,H_gb2,H_gb3
      double precision CY_,CKe_
      double precision El_,Eth_t
      double precision Z1_,Z2_,M1_,M2_

      double precision psi_zig_fs
      double precision psiens,psia0,psiwsub

      double precision z1,z2,m,m1,m2
      double precision sp,s,ss,f,ff,fp,g,gg
      double precision Gam,Al
      double precision A0,ro,us


      psi_yippu_prp=0

      z1 = dble(izp)
      m1 = mp
      z2 = dble(izt)
      m2 = mt

      Z1_ = z1
      Z2_ = z2
      M1_ = m1
      M2_ = m2

        m = m1/m2
      if(m < CLight) then
       Light=1
      else 
       Light=0
      endif

      pi0=dacos(-1.d0)

      s = m1+m2
      f = s*s 
      Gam = 4.d0*m1*m2/f

      sp = 1.d0/3.d0
      fp = z1**sp
       f = fp*fp
      sp = 2.d0*sp
      ff = z2**sp
      ss = f+ff
       s = dsqrt(ss)
      ss = ss*s

      Al  = 0.468503d0/s
      El_ = 30.74d0 * (m1+m2)/m2 * z1*z2 * s

      ro  = psiens(izt)
      A0  = psia0(izt)
      us  = psiwsub(izt)

      Eth_t = uth
 
      if(Light .ne. 0) then
       L_fb  = dsqrt(us) * (0.94 - 0.00133d0 * m2/m1)
       L_gb1 = Al / A0
       L_gb2 = 4.d0*us/(Gam*El_*El_)
      else
       H_fb  = psi_zig_fs(m2/m1)
       H_gb1 = pi0/180.d0 * 286.d0
       H_gb2 = (-0.45) / 2.d0
       g  = Al/A0
       gg = g * g * g
       g  = z1*z2/s
       gg = g*gg
       g  = dsqrt(gg)
       H_gb_ = g
       gg = 0.45d0 * dlog(g)
       H_gb3 = gg
      endif

      return
      end


      subroutine psi_yippu_fget(e, fopt,copt)
      implicit none
      double precision e
      double precision fopt,copt

      common/psiyippui/Light
      integer    Light
      common/psiyippud/
     ,      pi0,
     ,      L_fb,L_gb1,L_gb2,
     ,      H_fb,H_gb_,H_gb1,H_gb2,H_gb3,
     ,      CY_,CKe_,
     ,      El_,Eth_t,
     ,      Z1_,Z2_,M1_,M2_
      double precision pi0
      double precision L_fb,L_gb1,L_gb2
      double precision H_fb,H_gb_,H_gb1,H_gb2,H_gb3
      double precision CY_,CKe_
      double precision El_,Eth_t
      double precision Z1_,Z2_,M1_,M2_

      double precision emin 
      data emin /1.d-4/
      double precision smin 
      data smin /0.1d0/

      double precision em,s,f

      if(e < emin) then
       em=emin
      else
       em=e
      endif

      if(Light.gt.0) then
       fopt = L_fb
       s = dsqrt(L_gb2 * em)
       f = L_gb1 / dsqrt(s)
      else
       s = dsqrt(em / Eth_t) - 1.d0
       if(s < smin) then
         s=smin
       endif
       f = 1.d0 + 2.5 / s
       fopt = H_fb * f
       s = H_gb2 * dlog(em) + H_gb3
       f = H_gb1 * dexp(s)
      endif

      if(f*2.0 > pi0) then
       copt = 1.d0
      elseif(f < 0.d0) then
       copt = 0.d0
      else
       copt = dsin(f)
      endif

      end


      double precision function psi_yippu_fteta(e,amu)
      implicit none
      double precision e,amu
c
c---------------------------------------------
c  Angular dependence of sputtering yields
c	   of monoatomic solids
c
c    Y Yamamura, Y Itikawa, N. Itoh
c    IPPJ-AM-26 (1983)
c---------------------------------------------
c
      double precision seps
      data seps /1.d-18/

      double precision x,f,ff
      double precision fb,fbopt
   
      x = amu
      if(x < 0.d0) then
      x=-x
      endif

      if(x > 1.d0-seps) then
       psi_yippu_fteta=1.d0
      elseif(x < seps) then
       psi_yippu_fteta=0.d0
      else
         call psi_yippu_fget(e, fb,fbopt) 
         x = 1.d0/x
         f = dlog(x)
        ff = fbopt * (1.d0 - x) + f
         f = dexp(fb*ff)
       psi_yippu_fteta=f
      endif

      return
      end


