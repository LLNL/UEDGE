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
c

      double precision function psiPbgf(g)
      implicit none
      double precision g
      double precision f

      if(g .gt. 0.5d0) then
       f=0.9d0
      else if(g .gt. 0.09d0) then
       f=0.5d0
      else
       f=0.15d0
      endif

      psiPbgf=f
      return

      end



      double precision function psiPbtf(t,tm)
      implicit none
      double precision t,tm
      double precision f,tmm
      double precision ctm
      double precision gc 
      double precision tm0
      ctm=0.5d0
      gc =0.2d0
      tm0=320.0d0

      tmm=tm*ctm
     
      if(t .lt. tm0+0.1d0) then
       f=0.d0
      else if(t .lt. tmm) then
       f=gc/(tmm-tm0)*(t-tm0)
      else if(t .lt. tm) then
       f=gc+(1.d0-gc)/(tm-tmm)*(t-tmm)
      else
       f=1.0d0
      endif

      psiPbtf=f
      return

      end



      double precision function psiPbt(iz,t)
      implicit none
      integer iz
      double precision t
c
c Pb(T), T=Kelvin 
c
      integer idt
      parameter (idt=1)
      integer ierr,psimatype
      double precision psi_Pb,psiPbgf,psiPbtf,psi_Tsub
      double precision tabs
      double precision f,g,gg,tm
      tabs=273.15d0
         ierr=psimatype(iz)
      if(ierr .lt. 0) then
       f=0.d0
      else
       f=psi_Pb(iz)
       if(idt .gt. 0) then
        gg=1.d0-f
        g=psiPbgf(gg)
        gg=gg*g
        tm=psi_Tsub(iz)+tabs
        g=psiPbtf(t,tm)
        f=f+gg*g
       endif
      endif

      psiPbt=f            
      return

      end
