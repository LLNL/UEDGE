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
      subroutine psidotab
      implicit none

      Use (psitab)
      Use (psitab_s)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     *            (yres_mx,yres_mx_s),(epsp_mx,epsp_mx_s)

      integer ierr
      integer psitabw

      ierr=0

      call psi_tmesh_ini

      call psimatinp(ierr)
      if(ierr > 0) goto 1

      call psitabicro(ierr)
      if(ierr > 0) goto 1

      call psi_tmesh_par
c
c--calculate tables of thermochemical parameters for each matter
c
      call psimatpar0
      call psimatpar1
c
c---calculate PSI data tables
c
      call psitabdig
      call psitabsee
      call psitabtheps
c------
      ierr=psitabw(Nmat_ptm)

1     print *, '***Error: cannot create PSI tables!'
c      stop
      return
      
      end
