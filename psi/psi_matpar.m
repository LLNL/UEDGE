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
      subroutine psimatpar0

      implicit none

      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)
      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     ,            (yres_mx,yres_mx_s),(epsp_mx,epsp_mx_s)

      double precision psimassa,psiwsub,psia0
      double precision psi_Hmelt_ev,psi_Tmelt,psi_Tsub,psi_Aew

      integer i
      integer iz ! material atomic number
      double precision val

      do 1 i=1,Nmat_pt
       iz=IZmatter(i)
       val=psimassa(iz)
       mat_mas(i)=val
	
       val=psiwsub(iz)
       mat_wsub(i)=val
	
       val=psia0(iz)
       mat_a0(i)=val
	
       val=psi_Tmelt(iz)+273.15d0
       mat_Tmelt(i)=val
	
       val=psi_Hmelt_ev(iz)
       mat_Hmelt(i)=val
	
       val=psi_Tsub(iz)+273.15d0
       mat_Tsub(i)=val
	
       val=psi_Aew(iz)
       mat_Wf(i)=val
1     continue

      end


      subroutine psimatpar1

      implicit none

      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     ,            (yres_mx,yres_mx_s),(epsp_mx,epsp_mx_s)

      double precision psi_ro_t,psi_eBB_t,psi_cp_t
      double precision psi_kp_t,psi_gsub_t,psi_hsub_t
      double precision psiebbt,psiPbt

      integer iz
      integer i,j
      double precision val,td,gs,at

c
c--ro
c
      do 10 j=1,Nmat_pt
       iz=IZmatter(j)
       at=atdmi
       do 20 i=1,Nt_pt
        td=dexp(at)
        gs=psi_ro_t(iz,td)
        if(gs .lt. sml) then
         val=sml
        else if(gs .gt. big) then
         val=big
        else
         val=gs
        endif
        psi_ro(i,j)=val
        at=at+datdm
20     continue
10    continue

c
c--cp
c
      do 11 j=1,Nmat_pt
       iz=IZmatter(j)
       at=atdmi
       do 21 i=1,Nt_pt
        td=dexp(at)
        gs=psi_cp_t(iz,td)
        if(gs .lt. sml) then
         val=sml
        else if(gs .gt. big) then
         val=big
        else
         val=gs
        endif
        psi_cp(i,j)=val
        at=at+datdm
21     continue
11    continue

c
c--kp
c
      do 12 j=1,Nmat_pt
       iz=IZmatter(j)
       at=atdmi
       do 22 i=1,Nt_pt
        td=dexp(at)
        gs=psi_kp_t(iz,td)
        if(gs .lt. sml) then
         val=sml
        else if(gs .gt. big) then
         val=big
        else
         val=gs
        endif
        psi_kp(i,j)=val
        at=at+datdm
22     continue
12    continue

c
c--ep
c
      do 13 j=1,Nmat_pt
       iz=IZmatter(j)
       at=atdmi
       do 23 i=1,Nt_pt
        td=dexp(at)
cccc        gs=psi_eBB_t(iz,td)
c 2008.1.30
        gs=psiebbt(iz,td)

        if(gs .lt. sml) then
         val=sml
        else if(gs .gt. big) then
         val=big
        else
         val=gs
        endif
        psi_ep(i,j)=val
        at=at+datdm
23     continue
13    continue

c
c--hs
c
      do 14 j=1,Nmat_pt
       iz=IZmatter(j)
       at=atdmi
       do 24 i=1,Nt_pt
        td=dexp(at)
        gs=psi_hsub_t(iz,td)
        if(gs .lt. sml) then
         val=sml
        else if(gs .gt. big) then
         val=big
        else
         val=gs
        endif
        psi_hs(i,j)=val
        at=at+datdm
24     continue
14    continue

c
c--gs
c
      do 15 j=1,Nmat_pt
       iz=IZmatter(j)
       at=atdmi
       do 25 i=1,Nt_pt
        td=dexp(at)
        gs=psi_gsub_t(iz,td)
        if(gs .lt. sml) then
         val=asml
        else if(gs .gt. big) then
         val=abig
        else
         val=dlog(gs)
        endif
        psi_gs(i,j)=val
        at=at+datdm
25     continue
15    continue


c--Pb, coefficient of thermal emission 
c 2008.1.30

      do 16 j=1,Nmat_pt
       iz=IZmatter(j)
       at=atdmi
       do 26 i=1,Nt_pt
        td=dexp(at)
        gs=psiPbt(iz,td)

        if(gs .lt. sml) then
         val=sml
        else if(gs .gt. big) then
         val=big
        else
         val=gs
        endif
        psi_pbtd(i,j)=val
        at=at+datdm
26     continue
16    continue






      end


      double precision function psi_matpar_tf(t,tm)
      double precision t,tm

      double precision f,sincr1,sincr2,tm0
      data sincr1/0.9d0/, sincr2/0.8d0/, tm0/320.d0/
      
      if(t .le. tm0+0.1d0) then
       f=1.d0
      else if(t .lt. tm+0.1d0) then
       f=1.d0+(sincr1-1.d0)/(tm-tm0)*(t-tm0)
       if(f .lt. sincr1) f=sincr1
      else
c       f=sincr1+(sincr2-sincr1)/(tm*0.5d0)*(t-tm)
c       if(f .lt. 0.5d0) f=0.5d0
       f=sincr1
      endif

      psi_matpar_tf=f
      return

      end


      double precision function psi_ro_t(iz,t)
      implicit none
      integer iz
      double precision t
c
c---ro(T)
c   ro(g/cm^3) , T=Kelvin
c
      double precision psiens,psi_Tsub,psi_dens,psi_Tmelt
      double precision psimassa,psi_matpar_tf
      integer psimatype

      integer ierr,idt
      double precision f,g,mol1,tm

      data idt  /1/
      data mol1 /1.66053873d0/ !amu*10^24 or 10^24/mole 

      ierr=psimatype(iz)
      if(ierr .lt. 0) then
       f=0.d0
      else
       if(idt .gt. 0) then
        f=psi_dens(iz) !Angstrom^-3
        f=f*psimassa(iz)*mol1
       else
        f=psiens(iz)   !g/cm^3
       endif
       tm=psi_Tsub(iz)+273.15d0 ! boiling(sublimation) T, Kelvin
       g=psi_matpar_tf(t,tm)
       f=f*g
      endif
      
      psi_ro_t=f      
      return

      end


      double precision function psi_eBB_t(iz,t)
      implicit none
      integer iz
      double precision t
c
c semigeneral by Tanaka
c
c---eps_BB(T), T=Kelvin
c
      double precision psi_eBB,psi_matpar_tf,psi_Tsub
      integer psimatype

      integer ierr
      double precision f,g,tm

      ierr=psimatype(iz)
      if(ierr .lt. 0) then
       f=0.d0
      else
c       tm=psi_Tsub(iz)+273.15d0 ! sublimation T, Kelvin
       f=psi_eBB(iz)
c       g=psi_matpar_tf(t,tm)
c       f=f*g
      endif
      
      psi_eBB_t=f      
      return

      end


      double precision function psi_cp_t(iz,t)
      implicit none
      integer iz
      double precision t
c
c---cp(T): cp=J/ g grad K, T=Kelvin
c
      double precision psicpc
      integer psimatype

      integer ierr
      double precision f

      ierr=psimatype(iz)
      if(ierr .lt. 0) then
       f=0.d0
      else
       f=psicpc(iz,t)
      endif
      
      psi_cp_t=f      
      return

      end


      double precision function psi_kp_t(iz,t)
      implicit none
      integer iz
      double precision t
c
c---kappa(T): kappa=W/ m grad K, T=Kelvin
c
      double precision psikapc
      integer psimatype

      integer ierr
      double precision f

      ierr=psimatype(iz)
      if(ierr .lt. 0) then
       f=0.d0
      else
       f=psikapc(iz,t)
      endif
      
      psi_kp_t=f      
      return

      end


      double precision function psi_gsub_t(iz,t)
      implicit none
      integer iz
      double precision t
c
c---gamma_sublimation(T) in 10^20 atoms/cm^2/sec, T=Kelvin
c
      double precision psigsub
      integer psimatype,psicsub

      integer ierr,i
      double precision f

      ierr=psimatype(iz)
      if(ierr .lt. 0) then
       f=0.d0
      else
       i=psicsub(iz)
       f=psigsub(t)
      endif
      
      psi_gsub_t=f      
      return

      end


      double precision function psi_hsub_t(iz,t)
      implicit none
      integer iz
      double precision t
c
c---h_sublimation(T) in eV, T=Kelvin
c
      double precision psiesub
      integer psimatype,psicsub

      integer ierr,i
      double precision f,g

      ierr=psimatype(iz)
      if(ierr .lt. 0) then
       f=0.d0
      else
       i=psicsub(iz)
       g=psiesub(t)
       f=g
      endif
      
      psi_hsub_t=f      
      return

      end
