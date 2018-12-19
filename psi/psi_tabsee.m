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
      subroutine psitabsee
      implicit none

      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      Use (psitabmxx)
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     ,            (yres_mx,yres_mx_s),(epsp_mx,epsp_mx_s)
     
      integer          psi_usr_mx_fet
      external         psi_usr_mx_fet
      integer          psi_usr_fet
      external         psi_usr_fet
      double precision psisee
      external psisee

      integer ierr
      integer i,j,k,imat
      double precision Eth,m0
      double precision g,gf,gs
      double precision t,at,phi
      double precision val

      do 10 imat=1,Nmat_pt
c
c---secondary e-e emission yield
c
       call psiseep(IZmatter(imat))
       Eth=0.d0
       m0=1.d0

       k=0
       at=atmi
       do 7 i=1,Nt_pt
        t=dexp(at)
        do 8 j=1,Nphi_pt 
         phi=Phimin_pt+(j-1)*(Phimax_pt-Phimin_pt)/(Nphi_pt-1)
         if (phi .gt. 0.d0) then
          ierr=psi_usr_fet(psisee, t,phi,m0,Eth, g,gf,gs)
         else
          ierr=psi_usr_mx_fet(psisee, t,Eth, g,gf,gs)
         endif
         if(ierr .gt. 1) then
          val=asml
         else
          if(gs .lt. sml) then
           val=asml
          elseif(gs .gt. big) then
           val=abig
          else
           val=dlog(gs)
          endif
         endif
         k=k+1
         gdse_mx(i,j)=val
8       continue
        at=at+datm
7      continue

       do 2 j=1,Nphi_pt
        do 1 i=1,Nt_pt
         dse_mx(i,j,imat)=gdse_mx(i,j)
1       continue
2      continue
10    continue

      end
