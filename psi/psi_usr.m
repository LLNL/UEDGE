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
      double precision function psi_usr_ft_f(t)
      implicit none
      double precision t

      Use (psiusr)

      common/psiusr/Temp,Delta,Mach,Epsilon,Ht,Acur
      double precision Temp,Delta,Mach,Epsilon,Ht,Acur

      double precision tt, x,g,s,ss


	 tt = t
      if(tt < Teps) then
         tt = Teps
      else
           x=1.d0-Teps
        if(tt > x) then
         tt = x
        endif
      endif


      x = Epsilon


       s = tt * tt
      ss = 1.d0 - s
      if(ss < 0.d0) ss = 0.d0
       s = x*s
      ss = x*ss


       x = s - Delta
       g=0.d0
      if(x .gt. 0.d0) then

       s = dsqrt(x) - Mach
      if(s < 0.d0) s = -s
      ss = ss + s*s
       g = dexp( -ss )
       g = g * tt
      endif

      psi_usr_ft_f=g

      return

      end




      double precision function psi_usr_ft(eps)
      implicit none
      double precision eps

      Use (psiusr)

      common/psiusr/Temp,Delta,Mach,Epsilon,Ht,Acur
      double precision Temp,Delta,Mach,Epsilon,Ht,Acur
      common/psiusri/Nt
      integer Nt

      double precision psi_usr_ft_f
      external psi_usr_ft_f

      integer i
      double precision s,t, f,ff


       Epsilon = eps


       t = 0.d0
       f = 0.d0
       s = 0.d0
      do 1 i=1,Nt
       t = t + Ht

       ff = psi_usr_ft_f(t)

       s = s + (f+ff)*0.5d0
       f = ff
 1    continue


      f =  s * Epsilon*Ht

      psi_usr_ft=f
      return

      end




      double precision function psi_usr_fyt_f(t,Fusr)
      double precision t
      external Fusr
      double precision Fusr

      Use (psiusr)

      common/psiusr/Temp,Delta,Mach,Epsilon,Ht,Acur
      double precision Temp,Delta,Mach,Epsilon,Ht,Acur
 
      double precision e,y,tt, x,g,s,ss


         tt = t
      if(tt < Teps) then
         tt = Teps
      else
          x=1.d0-Teps
       if(tt > x)	then
        tt = x
       endif
      endif


      x = Epsilon
      e = x*Temp

      y = Fusr(e,tt)

      g = dexp( -x )
      g = g * tt * y

       s = tt * tt
      ss = 1.d0 - s
      if(ss < 0.d0) ss = 0.d0
       s = x*s
      ss = x*ss


      x = s - Delta
      g = 0.d0
      if(x .gt. 0.d0) then

       s = dsqrt(x) - Mach
      if(s < 0.d0) s = -s
      ss = ss + s*s
       g = dexp( -ss )
       g = g * tt * y
      endif

      psi_usr_fyt_f=g

      return

      end




      double precision function psi_usr_fyt(eps,Fusr)
      double precision eps
      external Fusr
      double precision Fusr

      Use (psiusr)

      common/psiusr/Temp,Delta,Mach,Epsilon,Ht,Acur
      double precision Temp,Delta,Mach,Epsilon,Ht,Acur
      common/psiusri/Nt
      integer Nt

      external psi_usr_fyt_f
      double precision psi_usr_fyt_f

      integer i
      double precision s,t, f,ff


      Epsilon = eps


      t = 0.d0
      f = 0.d0
      s = 0.d0
      do 1 i=1,Nt
       t = t + Ht

      ff = psi_usr_fyt_f(t,Fusr)
       s = s + (f+ff)*0.5d0
       f = ff
1     continue


      f =  s * Epsilon*Ht

      psi_usr_fyt=f
      return

      end




      integer function psi_usr_fet(FET_fun, Ti,Fi,M0, Eth, g,gf,gs)
      implicit none
      external FET_fun
      double precision FET_fun
      double precision Ti  !eV
      double precision Fi  !Z e Phi/ k/Ti
      double precision M0  !v/vTi
      double precision Eth !eV
      double precision g,gf,gs

      Use (psiusr)

      common/psiusr/Temp,Delta,Mach,Epsilon,Ht,Acur
      double precision Temp,Delta,Mach,Epsilon,Ht,Acur
      common/psiusri/Nt
      integer Nt

      external psigrl,psigrlf
      external psi_usr_ft,psi_usr_fyt
      double precision psigrl,psigrlf,psi_usr_ft,psi_usr_fyt

      integer, parameter :: N_lim=16 
      double precision Lim(16)
      data Lim/
     , 1.e-3,3.e-3,1.e-2,3.e-2,
     , 1.e-1,3.e-1,1.e0 ,3.e0 ,
     , 1.e1 ,3.e1 ,1.e2 ,3.e2 ,
     , 1.e3 ,3.e3 ,1.e4 ,3.e4/

       double precision f,ff,fs,gt, x,eps0,eps1
       integer i

       psi_usr_fet=0


       if(Ti < Ti_eps) then
       psi_usr_fet=1
       return
       endif

       Delta= Fi
       Mach = M0

       Temp = Ti


       Nt = N_teta
       ACur = E_eps
       x = dble(Nt)
       Ht = 1.d0/x


       f = 0.d0
       eps0 = Fi
      do 1 i=1,N_lim
       eps1 = Lim(i)
       if(eps0 < eps1) then
        gt=psigrl(psi_usr_ft, eps0,eps1, ACur)
        f=f+gt
       endif
       eps0 = eps1
 1    continue



         x = Eth
      if(x < 0.d0) x = 0.d0
      if(x < Fi*Ti) then
      eps0=Fi
      else

      eps0 = x/Ti
      endif

       ff = 0.d0
      do 2 i=1,N_lim
       eps1 = Lim(i)
       if(eps0 < eps1) then
        gt=psigrlf(psi_usr_fyt,FET_fun, eps0,eps1, ACur)
        ff=ff+gt
       endif
       eps0 = eps1
 2    continue



      if(f .le. 0.d0) then
       fs=0.d0
      else
       fs = ff/f
      endif


      g  = f 
      gf = ff
      gs = fs


      return
      end

