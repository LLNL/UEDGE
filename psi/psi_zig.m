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
      double precision function psi_zig_fs(x)
      implicit none
      double precision x

      double precision Ft(14)
      data Ft/
     ,   1.67d0,1.75d0,1.80d0,1.83d0,1.85d0,
     ,   1.81d0,1.65d0,1.45d0,1.36d0,1.26d0, 1.0d0,
     ,   0.75d0,0.25d0, 0.d0/
      double precision Xt(14)
      data Xt/
     ,   0.01d0, 0.1d0,  0.2d0, 0.4d0, 0.7d0,
     ,    1.0d0, 2.0d0,  3.0d0, 4.0d0, 5.0d0, 7.0d0,
     ,   10.0d0,20.0d0, 50.0d0/

      double precision f
      double precision psilintd

      f = psilintd(Ft,Xt,14, x)

      psi_zig_fs=f
      return
      end



      double precision function psi_zig_lm(x)
      implicit none
      double precision x

      double precision Ft(5)
      data Ft/
     ,  24.000000000000000d0,
     ,   3.800000000000000d0,
     ,   1.309000000000000d0,
     ,   0.327000000000000d0,
     ,   0.500000000000000d0/

      double precision Xt(5)
      data Xt/
     ,   0.0550000000000000d0,
     ,   0.1666666666666667d0,
     ,   0.3333333333333333d0,
     ,   0.5000000000000000d0,
     ,   1.0000000000000000d0/

      double precision f
      double precision psilintd

      f = psilintd(Ft,Xt,5, x)

      psi_zig_lm=f
      return
      end
