c/*
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
c*/
      subroutine psi_tmesh_ini
      implicit none

      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)

      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)
c
c---set maximal dimensions for meshes
c
      Nt_pt=Nt_ptm
      Ntd_pt=Ntd_ptm
      Nttm=Nt_ptm*Nt_ptm
      Nttdm=Nt_ptm*Ntd_ptm
      Ntphm=Nt_ptm*Nphi_ptm
      Nttdmc=-Ntd_ptm*(Nt_ptm+1)
      Ntm=Nt_ptm
      Nphi_pt=Nphi_ptm
c------------
      Nrd_pt=Nrd_ptm
      Ntde_pt=Ntde_ptm
      Nrtdem=Nrd_ptm*Ntde_ptm
c------------
      end


      subroutine psi_tmesh_par
      implicit none

      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)
c
c---log mesh on [Tmin,Tmax]
c
      atmi=dlog(Tmin_pt)
      atma=dlog(Tmax_pt)
      datm=(atma-atmi)/dble(Nt_pt-1)
c
c---log mesh on [Tdmin,Tdmax]
c
      atdmi=dlog(Tdmin_pt)
      atdma=dlog(Tdmax_pt)
      datdm=(atdma-atdmi)/dble(Ntd_pt-1)
c
c---linear mesh on [Phimin,Phimax]
c
      aphmi=Phimin_pt
      aphma=Phimax_pt
      daphm=(aphma-aphmi)/dble(Nphi_pt-1)
c
c---log mesh on [Tdmine,Tdmaxe]
c
      atdmie=dlog(Tdmine_pt)
      atdmae=dlog(Tdmaxe_pt)
      datdme=(atdmae-atdmie)/dble(Ntde_pt-1)
c
c---log mesh on [rdmin,rdmax]
c
      armi=dlog(rdmin_pt)
      arma=dlog(rdmax_pt)
      darm=(arma-armi)/dble(Nrd_pt-1)

      end
      
      
      subroutine psi_tdemesh_int(tde)
      implicit none
c
c---get the node for interpolation over Twall mesh
c
      
      real*8 tde
      
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      real*8 t

      t=tde

      if (t .le. Tdmine_pt) then
       iatdme=-1
       patdme=0.d0
      elseif (t .ge. Tdmaxe_pt) then
       iatdme=0
       patdme=1.d0
      else
ccc atdmie,atdmae,datdme,patdme
       patdme=(dlog(t)-atdmie)/datdme
       iatdme=int(patdme)
       patdme=patdme-dble(iatdme)
       iatdme=iatdme+1
      endif

      end

c-----

      subroutine psi_rdmesh_int(rdd)
      implicit none

      real*8 rdd
c
c---get the node for interpolation over rd mesh
c
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      real*8 r

      r=rdd

      if(r .le. rdmin_pt) then
       iarm=-1
       parm=0.d0
      elseif (r .ge. rdmax_pt) then
       iarm=0
       parm=1.d0
      else
ccc armi,arma,darm,parm,
       parm=(dlog(r)-armi)/darm
       iarm=int(parm)
       parm=parm-dble(iarm)
       iarm=iarm+1
      endif

      end



      subroutine psi_mesh_theps
      implicit none
c
c---get the Emissivity (rd, Td) via interpolation over2D
c   [rd,Td] mesh
c   Note,  psi_rdmesh_int(rd)
c          psi_tdemesh_int(td)
c   and    psi_tmesh_imat(imat)
c   must be called first
c

      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      Use (psiparloc)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      integer jarm,jatdme
      integer i
      integer i1,i2,i3,i4
      integer ii1,ii2,ii3,ii4
      real*8 g1,g2,g3,g4
      real*8 f,s,g,gtde
      
      if (iarm .ge. 1) then
       jarm=iarm
      elseif (iarm .eq. 0) then
       jarm=Nrd_pt-1
      else
       jarm=1
      endif


      if (iatdme .ge. 1) then
       jatdme=iatdme
      elseif (iatdme .eq. 0) then
       jatdme=Ntde_pt-1
      else
       jatdme=1
      endif

      g=1.d0-parm
      gtde=1.d0-patdme

      g1=patdme*g     !ix,iy+1
      g2=patdme*parm  !ix+1,iy+1
      g3=gtde*parm    !ix+1,iy
      g4=gtde*g       !ix,iy

      i1=jarm  + jatdme   *Nrd_pt   !ix,iy+1
      i2=jarm+1+ jatdme   *Nrd_pt   !ix+1,iy+1
      i3=jarm+1+(jatdme-1)*Nrd_pt   !ix+1,iy
      i4=jarm  +(jatdme-1)*Nrd_pt   !ix,iy

      i=(imat-1)*Nrtdem

      ii1=i1+i
      ii2=i2+i
      ii3=i3+i
      ii4=i4+i

      s=  dtheps_mx_s(ii1)*g1
      s=s+dtheps_mx_s(ii2)*g2
      s=s+dtheps_mx_s(ii3)*g3
      s=s+dtheps_mx_s(ii4)*g4

      if (s .lt. asml) then
       f=sml
      elseif (s .gt. abig) then
       f=big
      else
       f=s
      endif
      theps=f

      end



c===================================
      integer function psi_tmesh_imat(jmat)
      implicit none
      
      integer jmat
c
c---store current matter index for futher interpolation over Twall mesh
c
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      Use (psiparloc)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      integer ierr

      if (jmat .lt. 1) then
       ierr=1
      elseif (jmat .gt. Nmat_pt) then
       ierr=2
      else
       ierr=0
       imat=jmat
      endif

      psi_tmesh_imat=ierr
      return

      end


      subroutine psi_tmesh_int(tplasma)
      implicit none
      
      real*8 tplasma
c
c---get the node for interpolation over Tplasma mesh
c
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      real*8 t

      t=tplasma

      if (t .le. Tmin_pt) then
       iatm=-1
       patm=0.d0
      elseif (t .ge. Tmax_pt) then
       iatm=0
       patm=1.d0
      else
       patm=(dlog(t)-atmi)/datm
       iatm=int(patm)
       patm=patm-dble(iatm)
       iatm=iatm+1
      endif

      end



      subroutine psi_tmesh_yt
      implicit none
c
c---get the PSI-yields(Tplasma) via interpolation over Tplasma mesh
c   for all projectiles
c   Note, psi_tmesh_int(tplasma) and psi_tmesh_imat(imat) 
c   must be called first
c
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      Use (psiparloc)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      integer i,j,k,m,n
      real*8 f,s,g,p

      if(iatm.ge.1) then
       g=patm
       p=1.d0-g
       do 1 i=1,Nprj_pt
        j=I_pwi(i,imat)
        k=iatm+1
        n=(j-1)*Nt_pt
        m=iatm+n
        n=k+n
c        n=j+(k-1)*Npsi_pt
c        m=j+(iatm-1)*Npsi_pt

        s=rn_mx_s(m)*p+rn_mx_s(n)*g
        if(s .lt. asml) then
         f=sml
        else if(s .gt. abig) then
         f=big
        else
         f=dexp(s)
        endif
        frn  (i)=f

        s=re_mx_s(m)*p+re_mx_s(n)*g
        if(s .lt. asml) then
         f=sml
        else if(s .gt. abig) then
         f=big
        else
         f=dexp(s)
        endif
        fre  (i)=f

        s=ypsp_mx_s(m)*p+ypsp_mx_s(n)*g
        if(s .lt. asml) then
         f=sml
        else if(s .gt. abig) then
         f=big
        else
         f=dexp(s)
        endif
        fypsp(i)=f

        s=epsp_mx_s(m)*p+epsp_mx_s(n)*g
        if(s .lt. sml) then
         f=sml
        else if(s .gt. big) then
         f=big
        else
         f=s
        endif
        fepsp(i)=f
1      continue
      else
       if(iatm .eq. 0) then
        do 2 i=1,Nprj_pt
         j=I_pwi(i,imat)
         
         n=Nt_pt+(j-1)*Nt_pt
         
         s=rn_mx_s(n)
         if(s .lt. asml) then
          f=sml
         else if(s .gt. abig) then
          f=big
         else
          f=dexp(s)
         endif
         frn  (i)=f

         s=re_mx_s(n)
         if(s .lt. asml) then
          f=sml
         else if(s .gt. abig) then
          f=big
         else
          f=dexp(s)
         endif
         fre  (i)=f

         s=ypsp_mx_s(n)
         if(s .lt. asml) then
          f=sml
         else if(s .gt. abig) then
          f=big
         else
          f=dexp(s)
         endif
         fypsp(i)=f

         s=epsp_mx_s(n)
         if(s .lt. sml) then
          f=sml
         else if(s .gt. big) then
          f=big
         else
          f=s
         endif
         fepsp(i)=f
2       continue
       else
        do 3, i=1,Nprj_pt
         j=I_pwi(i,imat)

         n=1+(j-1)*Nt_pt

         s=rn_mx_s(n)
         if (s .lt. asml) then
          f=sml
         elseif (s .gt. abig) then
          f=big
         else
          f=dexp(s)
         endif
         frn  (i)=f

         s=re_mx_s(n)
         if (s .lt. asml) then
          f=sml
         elseif (s .gt. abig) then
          f=big
         else
          f=dexp(s)
         endif
         fre  (i)=f

         s=ypsp_mx_s(n)
         if (s .lt. asml) then
          f=sml
         elseif (s .gt. abig) then
          f=big
         else
          f=dexp(s)
         endif
         fypsp(i)=f

         s=epsp_mx_s(n)
         if (s .lt. sml) then
          f=sml
         elseif (s .gt. big) then
          f=big
         else
          f=s
         endif
         fepsp(i)=f
3       continue
       endif
      endif

      end


      subroutine psi_tmesh_ytt
      implicit none
c
c---get the PSI-yields(Tplasma,Twall) via interpolation over2D
c   [Tplasma,Twall] mesh
c   for all projectiles
c   Note, psi_tmesh_int(tplasma)
c         psi_tdmesh_int(twall)
c   and   psi_tmesh_imat(imat)
c   must be called first
c
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      Use (psiparloc)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      integer jatm,jatdm
      integer jatdmnt,jatdm1nt,jatm1
      integer i,j,k
      integer i1,i2,i3,i4
      integer ii1,ii2,ii3,ii4
      real*8 g1,g2,g3,g4
      real*8 f,s,g,gd
      
      if(iatm .ge. 1) then
       jatm=iatm
      elseif(iatm .eq. 0) then
       jatm=Nt_pt-1
      else
       jatm=1
      endif

      if(iatdm .ge. 1) then
       jatdm=iatdm
      elseif(iatdm .eq. 0) then
       jatdm=Ntd_pt-1
      else
       jatdm=1
      endif

      g =1.d0-patm
      gd=1.d0-patdm

      g1=patdm*g    !ix,iy+1
      g2=patdm*patm !ix+1,iy+1
      g3=gd*patm    !ix+1,iy
      g4=gd*g       !ix,iy

      jatdmnt=jatdm*Nt_pt
      jatdm1nt=jatdmnt-Nt_pt
      jatm1=jatm+1

      i1=jatm  + jatdmnt    !ix,iy+1
      i2=jatm1 + jatdmnt    !ix+1,iy+1
      i3=jatm1 + jatdm1nt   !ix+1,iy
      i4=jatm  + jatdm1nt   !ix,iy

      do 10 k=1,Nprj_pt
       j=I_pwi(k,imat)
       i=(j-1)*Nttdm
       
       ii1=i1+i
       ii2=i2+i
       ii3=i3+i
       ii4=i4+i
       
       s=  ycsp_mx_s(ii1)*g1
       s=s+ycsp_mx_s(ii2)*g2
       s=s+ycsp_mx_s(ii3)*g3
       s=s+ycsp_mx_s(ii4)*g4

       if(s .lt. asml) then
        f=sml
       else if(s .gt. abig) then
        f=big
       else
        f=dexp(s)
       endif
       fycsp(k)=f

       s=  yres_mx_s(ii1)*g1
       s=s+yres_mx_s(ii2)*g2 
       s=s+yres_mx_s(ii3)*g3 
       s=s+yres_mx_s(ii4)*g4
 
       if(s .lt. asml) then
        f=sml
       else if(s .gt. abig) then
        f=big
       else
        f=dexp(s)
       endif
       fyres(k)=f

10    continue

      end


      subroutine psi_mesh_yse
      implicit none
c
c---get the SEE coeff.(Tplasma,Phi) via interpolation over2D
c   [Tplasma,Phi] mesh
c   Note, psi_tmesh_int(tplasma)
c         psi_phmesh_int(phi)
c   and   psi_tmesh_imat(imat) 
c   must be called first
c
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      Use (psiparloc)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      integer jatm,japhm
      integer i
      integer i1,i2,i3,i4
      integer ii1,ii2,ii3,ii4
      real*8 g1,g2,g3,g4
      real*8 f,s,g,gph
      
      if(iatm .ge. 1) then
       jatm=iatm
      elseif(iatm .eq. 0) then
       jatm=Nt_pt-1
      else
       jatm=1
      endif

      if(iaphm .ge. 1) then
       japhm=iaphm
      elseif(iaphm .eq. 0) then
       japhm=Nphi_pt-1
      else
       japhm=1
      endif

      g  =1.d0-patm
      gph=1.d0-paphm

      g1=paphm*g     !ix,iy+1
      g2=paphm*patm  !ix+1,iy+1
      g3=gph*patm    !ix+1,iy
      g4=gph*g       !ix,iy

      i1=jatm  + japhm   *Nt_pt   !ix,iy+1
      i2=jatm+1+ japhm   *Nt_pt   !ix+1,iy+1
      i3=jatm+1+(japhm-1)*Nt_pt   !ix+1,iy
      i4=jatm  +(japhm-1)*Nt_pt   !ix,iy

      i=(imat-1)*Ntphm

      ii1=i1+i
      ii2=i2+i
      ii3=i3+i
      ii4=i4+i

      s=  dse_mx_s(ii1)*g1
      s=s+dse_mx_s(ii2)*g2
      s=s+dse_mx_s(ii3)*g3
      s=s+dse_mx_s(ii4)*g4
 
      if(s .lt. asml) then
       f=sml
      else if(s .gt. abig) then
       f=big
      else
       f=dexp(s)
      endif
      fysee=f

      end


      subroutine psi_phmesh_int(phi)
      implicit none
      
      real*8 phi
c
c---get the node for interpolation over Phi mesh
c
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      real*8 t

      t=phi
      if (t .le. Phimin_pt) then
       iaphm=-1
       paphm=0.d0
      elseif (t .ge. Phimax_pt) then
       iaphm=0
       paphm=1.d0
      else
       paphm=(t-aphmi)/daphm
       iaphm=int(paphm)
       paphm=paphm-dble(iaphm)
       iaphm=iaphm+1
      endif

      end



      subroutine psi_tdmesh_int(twall)
      implicit none
      
      real*8 twall
c
c---get the node for interpolation over Twall mesh
c
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      real*8 t

      t=twall

      if (t .le. Tdmin_pt) then
       iatdm=-1
       patdm=0.d0
      elseif (t .ge. Tdmax_pt) then
       iatdm=0
       patdm=1.d0
      else
       patdm=(dlog(t)-atdmi)/datdm
       iatdm=int(patdm)
       patdm=patdm-dble(iatdm)
       iatdm=iatdm+1
      endif

      end


      subroutine psi_tdmesh_rocp
      implicit none
c
c---get the matter(Twall) via interpolation over Twall mesh
c   Note, psi_tdmesh_int(twall) and psi_tdmesh_imat(imat) 
c   must be called first
c
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      Use (psiparloc)
      
      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      integer j,k,m
      real*8 p,g

      if(iatdm .ge. 1) then
       j=iatdm+1
       g=patdm
       p=1.d0-g
       k=(imat-1)*Ntd_pt
       m=j+k
       k=iatdm+k
       rotd=psi_ro_s(k)*p+psi_ro_s(m)*g
       cptd=psi_cp_s(k)*p+psi_cp_s(m)*g
      else
       if(iatdm .eq. 0) then
        k=Ntd_pt+(imat-1)*Ntd_pt
        rotd=psi_ro_s(k)
        cptd=psi_cp_s(k)
       else
        k=1+(imat-1)*Ntd_pt
        rotd=psi_ro_s(k)
        cptd=psi_cp_s(k)
       endif
      endif
      
      end


      subroutine psi_tdmesh_vall
      implicit none
c
c---get the matter(Twall) via interpolation over Twall mesh
c   Note, psi_tdmesh_int(twall) and psi_tdmesh_imat(imat) 
c   must be called first
c
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      Use (psiparloc)
      
      equivalence (psi_ro_s,psi_ro),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
      equivalence (dtheps_mx,dtheps_mx_s)

      integer j,k,m
      real*8 g,s,p

      if(iatdm .ge. 1) then
       j=iatdm+1
       g=patdm
       p=1.d0-g
       k=(imat-1)*Ntd_pt
       m=j+k
       k=iatdm+k
       rotd=psi_ro_s(k)*p+psi_ro_s(m)*g
       cptd=psi_cp_s(k)*p+psi_cp_s(m)*g
       kptd=psi_kp_s(k)*p+psi_kp_s(m)*g
       eptd=psi_ep_s(k)*p+psi_ep_s(m)*g
       hstd=psi_hs_s(k)*p+psi_hs_s(m)*g
       pbtd=psi_pbtd_s(k)*p+psi_pbtd_s(m)*g

       s=psi_gs_s(k)*p+psi_gs_s(m)*g
       if (s .lt. asml) then
        gstd=sml
       elseif (s .gt. abig) then
        gstd=big
       else
        gstd=dexp(s)
       endif
      else
       if (iatdm .eq. 0) then
        k=Ntd_pt+(imat-1)*Ntd_pt
        rotd=psi_ro_s(k)
        cptd=psi_cp_s(k)
        kptd=psi_kp_s(k)
        eptd=psi_ep_s(k)
        hstd=psi_hs_s(k)
        pbtd=psi_pbtd_s(k)
        s=psi_gs_s(k)
        if(s .lt. asml) then
         gstd=sml
        elseif(s .gt. abig) then
         gstd=big
        else
         gstd=dexp(s)
        endif
       else
        k=1+(imat-1)*Ntd_pt
        rotd=psi_ro_s(k)
        cptd=psi_cp_s(k)
        kptd=psi_kp_s(k)
        eptd=psi_ep_s(k)
        hstd=psi_hs_s(k)
        pbtd=psi_pbtd_s(k)
        s=psi_gs_s(k)
        if (s .lt. asml) then
         gstd=sml
        elseif (s .gt. abig) then
         gstd=big
        else
         gstd=dexp(s)
        endif
       endif
      endif

      end
