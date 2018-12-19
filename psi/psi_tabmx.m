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
      subroutine psitabdig
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

      integer k

      do 1 k=1,Npsi_pt
       call psitabdo(k)
       call psitabsave(k)
1     continue
      end

c ----------------------------------------------------------------------
      subroutine psitabdo(kpm)
      implicit none
      integer kpm

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
      integer          psi_bgh_prp_OLD
      external         psi_bgh_prp_OLD
      double precision psi_rn,psi_re
      external         psi_rn,psi_re
      double precision psi_bgh_yn,psi_bgh_ye,psi_bgh_get_wth
      external         psi_bgh_yn,psi_bgh_ye,psi_bgh_get_wth
      double precision psicspy,psiresyy,psiwres
      external         psicspy,psiresyy,psiwres
      double precision psimassa

      integer izp,izt
      double precision mp,mt

      integer ierr
      double precision Eth
      double precision g,gf,gs
      double precision t,at
      double precision td,atd
      double precision val
      integer i,j,k

        i=I_wal(kpm)
      izt=IZmatter(i)
       mt=psimassa(izt)

        j=I_prj(kpm)
      izp=IZprjctl(j)
       mp=AZprjctl(j)
      
c
c---reflection Rn
c
      call psirnep(izp,mp,izt,mt)
      Eth=0.d0

      at=atmi
      do 1 i=1,Nt_pt
       t=dexp(at)
       ierr=psi_usr_mx_fet(psi_rn, t,Eth, g,gf,gs)
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
       grn_mx(i)=val
       at=at+datm
1     continue

c
c---reflection Re
c
      call psirnep(izp,mp,izt,mt)
      Eth=0.d0

      at=atmi
      do 2 i=1,Nt_pt
       t=dexp(at)
       ierr=psi_usr_mx_fet(psi_re, t,Eth, g,gf,gs)
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
       gre_mx(i)=val
       at=at+datm
2     continue

c
c---physical sputtering yield
c
      ierr=psi_bgh_prp_OLD(izp,mp,izt,mt)
      if(ierr .gt. 0) then
       do 4 i=1,Nt_pt
        gypsp_mx(i)=asml
4      continue
      else
       Eth=psi_bgh_get_wth(1)
       at=atmi
       do 3 i=1,Nt_pt
        t=dexp(at)
        ierr=psi_usr_mx_fet(psi_bgh_yn, t,Eth, g,gf,gs)
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
        gypsp_mx(i)=val
        at=at+datm
3      continue

       at=atmi
       do 5 i=1,Nt_pt
        t=dexp(at)
        ierr=psi_usr_mx_fet(psi_bgh_ye, t,Eth, g,gf,gs)
        if(ierr .gt. 1) then
         val=0.d0
        else
         if(gs .lt. sml) then
          val=0.d0
         elseif(gs .gt. big) then
          val=0.d0
         else
          val=gs*dexp(-gypsp_mx(i))
         endif
        endif
        gepsp_mx(i)=val
        at=at+datm
5      continue
      endif

c
c---chemical sputtering yield
c
      call psiyspc(izp,mp,izt,mt)
      Eth=0.d0

      k=0
      atd=atdmi
      do 8 j=1,Ntd_pt
       td=dexp(atd)
       call psicsptd(td)
       at=atmi
       do 7 i=1,Nt_pt
        t=dexp(at)
        ierr=psi_usr_mx_fet(psicspy, t,Eth, g,gf,gs)
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
        gycsp_mx(i,j)=val
        at=at+datm
7      continue
       atd=atd+datdm
8     continue

c
c---RES yield
c
      call psiresc(izp,mp,izt,mt)
      Eth=psiwres(1)

      k=0
      atd=atdmi
      do 9 j=1,Ntd_pt
       td=dexp(atd)
       call psirestd(td)
       at=atmi
       do 10 i=1,Nt_pt
        t=dexp(at)
        ierr=psi_usr_mx_fet(psiresyy, t,Eth, g,gf,gs)
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
        gyres_mx(i,j)=val
        at=at+datm
10     continue
       atd=atd+datdm
9     continue

      end


      subroutine psitabsave(kpm)
      implicit none
      integer kpm

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

      integer i,j

      do 10 i=1,Nt_pt
       rn_mx(i,kpm)=grn_mx(i)
       re_mx(i,kpm)=gre_mx(i)
       ypsp_mx(i,kpm)=gypsp_mx(i)
       epsp_mx(i,kpm)=gepsp_mx(i)
10    continue

      do 2 j=1,Ntd_pt
       do 1 i=1,Nt_pt
        ycsp_mx(i,j,kpm)=gycsp_mx(i,j)
        yres_mx(i,j,kpm)=gyres_mx(i,j)
1      continue
2     continue

      end
