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
      double precision function psi_lss_sn_TF(epsx)
      implicit none
      double precision epsx

c
c----------------------------------
c      Thomas-Fermi potential
c  stopping power cross section
c   for elastic scattering
c
c    n.matsunami et al
c    rad.eff.lett 50 39 1980
c-----------------------------------
c
      double precision, parameter :: LSS_SN_MIN=1.d-16

      double precision f,se,eps,ff

      if(epsx < LSS_SN_MIN) then
        eps = LSS_SN_MIN
      else 
        eps = epsx
      endif

      se = dsqrt(eps)
      f  = eps + 2.718
      ff = dlog(f)
      f  = ff * se * 3.441
      ff = -1.708 + 6.882 * se
      ff = 1.0 + 6.355 * se + eps*ff
      f  = f/ff

      psi_lss_sn_TF=f
      return

      end


      double precision function psi_lss_sn_KrC(epsx)
      implicit none
      double precision epsx

c
c----------------------------------
c	Kr-C potential
c   stopping power cross section
c    for elastic scattering
c
c    n.matsunami et al
c    rad.eff.lett 50 39 1980
c-----------------------------------
c
      double precision,  parameter :: LSS_EPS_SN_MIN=1.d-16

      double precision g,f,se,eps,ff
      double precision acc

      data acc/1.d-40/

      if(epsx < LSS_EPS_SN_MIN) then
       eps = LSS_EPS_SN_MIN
      else                      
       eps = epsx
      endif

      f  = 0.5 * dlog( 1.0 + 1.2288 * eps )

      se = dsqrt(eps)
      ff = eps + se * 0.1728
      ff = ff  + 0.008 * eps**0.1504

      if(ff .lt. 0.d0) then
       g = -ff
      else 
       g = ff
      endif

      if(g .lt. acc) then
       f = 0.d0
      else 
       f = f/ff
      endif

      psi_lss_sn_KrC=f
      return

      end


      double precision function psi_lss_sn_AKB(epsx)
      implicit none
      double precision epsx

c
c----------------------------------
c	Andersen potential
c   stopping power cross section
c    for elastic scattering
c
c       s kalbitzer et al.
c     z physica a278 223 1976
c  'the stopping power and ranges
c       of ions in matter'
c	 pp n4 1979 v4
c-----------------------------------
c
      double precision,  parameter :: LSS_EPS_SN_MIN=1.d-16
      double precision f,se,eps

      if(epsx < LSS_EPS_SN_MIN) then
       eps = LSS_EPS_SN_MIN
      else 
       eps = epsx
      endif

      if( eps .lt. 0.01d0 ) then
	f = dsqrt(eps) * 1.593
      else
       if (eps .gt. 10.d0) then
	f = dlog(0.47 * eps)/(eps * 2.d0)  
       else
        se=2.d0/3.d0
	f = 3.4*eps**se
	f = f+1.d0 + 6.8 * eps
	f = ( eps + 2.718284090d0 ) / f
	se = dsqrt(eps)
	f = f*1.7d0 * se
       endif
      endif

      psi_lss_sn_AKB=f
      return

      end



      double precision function psi_lss_sn(eps,j)
      implicit none
      double precision eps
      integer j

      double precision psi_lss_sn_KrC,psi_lss_sn_TF,psi_lss_sn_AKB
      double precision f

      f=0.d0
      if(j .eq. 1) then
       f = psi_lss_sn_KrC(eps)
      else
       if(j .eq. 2) then
        f = psi_lss_sn_TF (eps)
       else
        if(j .eq. 3) then 
         f = psi_lss_sn_AKB(eps)
        else  
         f = psi_lss_sn_TF (eps)
        endif
       endif
      endif

      psi_lss_sn=f
      return

      end

