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
      integer function psitabw(Nmatmax)
      implicit none
      integer Nmatmax
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


      integer ierr,i
c     
      ierr=0
      if(Nmat_pt .gt. Nmatmax) then
      ierr=100
      print *,'psitabw:error in logic'
      goto 1
      endif
c     
      do 100 i=1,Nmat_pt
       call psitab_w(ierr,i)
       if (ierr .gt. 0) goto 1
100   continue
c     
1     psitabw=ierr
      return
c     
      end
c
      integer function psitabr(Nmatmax)
      implicit none
      integer Nmatmax
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


c
      integer ierr,i
c
      ierr=0
      if(Nmat_pt .gt. Nmatmax) then
       ierr=100
       print *,'psitabr:error in logic'
       goto 1
      endif
c
      do 100 i=1,Nmat_pt
       call psitab_r(ierr,i)
       if (ierr .gt. 0) goto 1
100   continue
c
1     psitabr=ierr
      return
      end
c      
      subroutine psitab_w(ierr,imatter)
      implicit none
      integer ierr,imatter

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
      integer Iprjj(Nprj_ptm),Iprjz(Nprj_ptm)
      real*8 tk(Ntd_ptm),tpl(Nt_ptm),tphi(Nphi_ptm)
      real*8 at
      integer i,j,k
      integer ican
      data ican/1/

      ierr=0

      if(imatter .lt. 1) goto 2
      if(imatter .gt. Nmat_pt) goto 2

      open(ican,file=Namfpsi(imatter), err=3)

      write(ican,*) 'PSI data file created for UEDGE'

      write(ican,100) 'matter=',IZmatter(imatter)
100   format(A10,i4)

      write(ican,101) 'number of projectiles=',Nprj_pt
101   format(A24,i4)
      write(ican,102) 'projectiles=',
     +                (IZprjctl(i),i=1,Nprj_pt)
102   format(A24,10i4)
      write(ican,103) 'projectile mass=',
     +                (AZprjctl(i),i=1,Nprj_pt)
103   format(A24,10F10.3)

      write(ican,104) 'number of nodes in Tplasma mesh=',
     +                Nt_pt
104   format(A36,I5)
      write(ican,105) 'Tplasma mesh min=',Tmin_pt
105   format(A24,d25.15)
      write(ican,105) 'Tplasma mesh max=',Tmax_pt

      write(ican,104) 'number of nodes in Twall mesh=',
     +                Ntd_pt
      write(ican,105) 'Twall mesh min=',Tdmin_pt
      write(ican,105) 'Twall mesh max=',Tdmax_pt

      write(ican,106) 'matter atom mass, amu=',
     +                mat_mas(imatter)
      write(ican,106) 'sublimation heat, eV=',
     +                mat_wsub(imatter)
      write(ican,106) 'lattice constant, Angstroms=', 
     +                mat_a0(imatter)
      write(ican,106) 'melting temperature, Kelvin=', 
     +                mat_Tmelt(imatter)
      write(ican,106) 'melting enthalpy, eV=', 
     +                mat_Hmelt(imatter) 
      write(ican,106) 'sublimation temperature, Kelvin=',
     +                mat_Tsub(imatter)
      write(ican,106) 'work function, eV=',
     +                mat_Wf(imatter)
106   format(A36,d25.15)

      at=atdmi
      do 5 i=1,Ntd_pt
      tk(i)=dexp(at)
      at=at+datdm
5     continue

      at=atmi
      do 4 i=1,Nt_pt
       tpl(i)=dexp(at)
       at=at+datm
4     continue

      at=aphmi
      do 203 i=1,Nphi_pt
      tphi(i)=at
      at=at+daphm
203   continue

      write(ican,107) 'ro, matter mass density,g/cm3'
      write(ican,107) 'cp, heat capacity,J/g grad K'
      write(ican,107) 'kp, heat conductivity,W/m grad K'
      write(ican,107) 'ep, Black-Body emissivity coefficient'
      write(ican,107) 'hs, thermal emission heat,eV'
      write(ican,107) 'gs,(LOG)thermal emission flux,10^20/cm2/sec'
      write(ican,107) 'pb, Coefficient of thermal emission'
      write(ican,107) 'T(K)   ro   cp   kp   ep   hs   gs  pbtd:'
107   format(A52)

      do 10 i=1,Ntd_pt
      write(ican,110) tk(i),
     +                psi_ro(i,imatter),
     +                psi_cp(i,imatter),
     +                psi_kp(i,imatter),
     +                psi_ep(i,imatter),
     +                psi_hs(i,imatter),
     +                psi_gs(i,imatter),
     +                psi_pbtd(i,imatter)
110   format(E14.3,7(d25.15))
10    continue

      do 12 i=1,Nprj_pt
       j=I_pwi(i,imatter)
       Iprjj(i)=j
       Iprjz(i)=IZprjctl(i)
12    continue

c
c--Rn
c
      write(ican,107) '<Rn>Mx, particle reflection coefficient'
      write(ican,108) 'Tplasma(eV), <Rn>Mx for projectiles:',
     +                (Iprjz(i),i=1,Nprj_pt)
108   format(A52,10I4)
      do 11 i=1,Nt_pt
       write(ican,111) tpl(i),(rn_mx(i,Iprjj(j)),j=1,Nprj_pt)      
111   format(d15.3,10(d25.15))
11    continue

c
c--Re
c
      write(ican,107) '<Re>Mx, energy reflection coefficient'
      write(ican,108) 'Tplasma(eV), <Re>Mx for projectiles:',
     +                (Iprjz(i),i=1,Nprj_pt)
      do 15 i=1,Nt_pt
      write(ican,111) tpl(i),(re_mx(i,Iprjj(j)),j=1,Nprj_pt)      
15    continue

c
c--Yphys
c
      write(ican,107) '<Yph>Mx, physical sputtering yield'
      write(ican,108) 'Tplasma(eV), <Yph>Mx, <Eph>Mx for projectiles:',
     +                (Iprjz(i),i=1,Nprj_pt)
      do 16 i=1,Nt_pt
      write(ican,113) tpl(i),(ypsp_mx(i,Iprjj(j)),epsp_mx(i,Iprjj(j)),
     +                        j=1,Nprj_pt)      
113   format(d15.3,20(d25.15))
16    continue

c
c--Ychem
c
      write(ican,107) '<Ych>Mx, chemical sputtering yield'
      write(ican,108) 'Tplasma,Twall(eV),<Ych>Mx for projectiles:',
     +                (Iprjz(i),i=1,Nprj_pt)
      do 18 k=1,Ntd_pt
       do 17 i=1,Nt_pt
       write(ican,112) tpl(i),tk(k),
     +                 (ycsp_mx(i,k,Iprjj(j)),j=1,Nprj_pt)
112    format(d15.3,d15.3,10(d25.15))      
17     continue
18    continue

c
c--Yres
c
      write(ican,107) '<Yres>Mx, RES,radiation enhanced sublimation'
      write(ican,108) 'Tplasma,Twall(eV),<Yres>Mx for projectiles:',
     +                (Iprjz(i),i=1,Nprj_pt)
      do 19 k=1,Ntd_pt
       do 20 i=1,Nt_pt
       write(ican,112) tpl(i),tk(k),
     +                 (yres_mx(i,k,Iprjj(j)),j=1,Nprj_pt)      
20     continue
19    continue

c
c--Ysee
c
      write(ican,207) '<Ysee>Mx, SEE,secondary electron-electron
     + emission'
      write(ican,207) 'Tplasma(eV),Phi/T,<Ysee>Mx'
      do 219 k=1,Nphi_pt
       do 220 i=1,Nt_pt
       write(ican,212) tpl(i),tphi(k),dse_mx(i,k,imatter)
220    continue
219   continue
207   format (A52)
212   format (d15.3,d15.3,d25.15)


c====== add on 2008.1.20 ============================================================
c
c-Emissivity
c
      write(ican,307) '<Eps>Mx, EMS, thermal emissivity from Mie'
      write(ican,307) 'Radius(micron),Tdust(K),<Ytheps>Mx'
      do 319 k=1,Ntde_pt
       do 320 i=1,Nrd_pt
        write(ican,312) rdk(i),tkkk(k),dtheps_mx(i,k,imatter)
320    continue
319   continue
307   format (A52)
312   format (d15.3,d15.3,d25.15)

c=============================================================================

      write(ican,107) 'That''s all folks'

      goto 1

2     ierr=1
      print *, 'psitab_w:error in imatter=',imatter
      goto 22
      
3     ierr=2
      print *, 'psitab_w:error in doing file=', 
     +         Namfpsi(imatter)

1     close (ican)
22    continue

      end


      subroutine psitab_r(ierr,imatter)
      implicit none
      
      integer ierr,imatter

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

      integer Ipjj(Nprj_ptm),Iprjj(Nprj_ptm),Iprjz(Nprj_ptm)
      real*8 Cjj(Nprj_ptm),C2jj(Nprj_ptm)
      real*8 Apij(Nprj_ptm)
      real*8 at,g,gg
      real*8 ad,d,dd,c
      integer i,ii,j,k,kk
      integer n
      character*1   a
      character*256 atxt

      integer ican
      data ican/1/

      ierr=0

      if (imatter .lt. 1) goto 2
      if (imatter .gt. Nmat_pt) goto 2

      open(ican, file=Namfpsi(imatter), action='read',
     +status='old', form='formatted', err=3)
      rewind(ican)

      read(ican,99) a !'PSI data file created for UEDGE'
99    format(A1)

      n=0
      read(ican,100) atxt,n !'matter=',IZmatter(imatter)
100   format(A10,I4)
      k=0
      if(n .lt. 1 ) k=k+1
      if(n .gt. 92) k=k+1
      if(k .gt. 0) then
       ierr=1
       print *,'improper matter=',n
       goto 3
      endif
      if(n .ne. IZmatter(imatter)) then
       ierr=2
       print *,'inconsistent matter=',n,
     +         'instead of =',IZmatter(imatter)
       goto 3
      endif

      n=0
      read(ican,101) atxt,n !'number of projectiles=',Nprj_pt
101   format(A24,i4)
      if(n .ne. Nprj_pt) then
       ierr=3
       print *,'inconsistent number of projectiles=',n,
     +         'instead of =',Nprj_pt
       goto 3
      endif

      do 30 i=1,Nprj_pt
       Ipjj(i)=0
       Iprjj(i)=0
30    continue      
      read(ican,102) atxt,(Ipjj(i),i=1,Nprj_pt) !'projectiles=',IZprjctl
102   format(A24,10i4)
      do 23 i=1,Nprj_pt
       j=Ipjj(i)
       if(j .lt. 0) then
        ierr=4
       else 
        if(j .gt. 92) then
         ierr=4
        endif
       endif
       if(ierr .gt. 0) then
        print *,'incorrect projectile=',j,' index=',i
       endif
23    continue
      if(ierr .gt. 0) goto 3
      do 21 i=1,Nprj_pt
       k=IZprjctl(i)
       kk=0
       do 22 j=1,Nprj_pt
        if(k .eq. Ipjj(j)) then
         kk=kk+1
         Iprjj(i)=j
        endif
22     continue
       if(kk .lt. 1) then
        ierr=4
        print *,'missing projectile=',k      
       else
        if(kk .gt. 1) then
         ierr=4
         print *,'repeated projectile=',k
        endif
       endif
21    continue
      if(ierr .gt. 0) goto 3

      do 32 i=1,Nprj_pt
       Apij(i)=0.0
32    continue
      read(ican,103) atxt,(Apij(i),i=1,Nprj_pt) !'projectile mass='
103   format(A24,10F10.3)
      do 33 i=1,Nprj_pt
      j=Iprjj(i)
      at=Apij(i)
      if(at .lt. 0.3) then
      ierr=5
      else if(at .gt. 1.e3) then
      ierr=5
      else
      g=AZprjctl(i)
      gg=(g+at)*0.25
      at=at-g
      if(at .lt. 0.0) at=-at
      if(at .gt. gg) then
      ierr=4
      endif
      endif
      if(ierr .gt. 0) then
      print *,'inconsistent list of projectile mass=',Apij(i),
     +        'instead of =',AZprjctl(i),ierr
      endif
33    continue
      if(ierr .gt. 0) goto 3

      n=0
      read(ican,104) atxt,n !'number of nodes in Tplasma mesh=',Nt_pt
104   format(A36,I5)
      if(n .ne. Nt_pt) then
      ierr=5
      print *,'inconsistent number of nodes in Tplasma mesh=',n,
     +        'instead of =',Nt_pt
      goto 3
      endif

      d=0.d0
      read(ican,105) atxt,d !'Tplasma mesh min=',Tmin_pt
105   format(A24,d25.15)
      if(d .lt. 0.d0) then
      ierr=6
      else
      dd=(d+Tmin_pt)*0.25
      ad=d-Tmin_pt
      if(ad .lt. 0.0) ad=-ad
      if(ad .gt. dd) then
      ierr=6
      endif
      endif
      if(ierr .gt. 0) then
      print *,'inconsistent Tplasma mesh min=',d,
     +        'instead of =',Tmin_pt
      goto 3
      endif

      d=0.d0
      read(ican,105) atxt,d !'Tplasma mesh max=',Tmax_pt
      if(d .lt. 0.d0) then
      ierr=7
      else
      dd=(d+Tmax_pt)*0.25
      ad=d-Tmax_pt
      if(ad .lt. 0.0) ad=-ad
      if(ad .gt. dd) then
      ierr=7
      endif
      endif
      if(ierr .gt. 0) then
      print *,'inconsistent Tplasma mesh max=',d,
     +        'instead of =',Tmax_pt
      goto 3
      endif

      n=0
      read(ican,104) atxt,n !'number of nodes in Twall mesh=',Ntd_pt
      if(n .ne. Ntd_pt) then
      ierr=8
      print *,'inconsistent number of nodes in Twall mesh=',n,
     +        'instead of =',Ntd_pt
      goto 3
      endif

      d=0.d0
      read(ican,105) atxt,d !'Twall mesh min=',Tdmin_pt
      if(d .lt. 0.d0) then
       ierr=9
      else
       dd=(d+Tdmin_pt)*0.25
       ad=d-Tdmin_pt
      if(ad .lt. 0.0) ad=-ad
       if(ad .gt. dd) then
        ierr=9
       endif
      endif
      if(ierr .gt. 0) then
       print *,'inconsistent Twall mesh min=',d,
     +         'instead of =',Tdmin_pt
      goto 3
      endif

      d=0.d0
      read(ican,105) atxt,d !'Twall mesh max=',Tdmax_pt
      if(d .lt. 0.d0) then
       ierr=10
      else
       dd=(d+Tdmax_pt)*0.25
       ad=d-Tdmax_pt
      if(ad .lt. 0.0) ad=-ad
       if(ad .gt. dd) then
        ierr=10
       endif
      endif
      if(ierr .gt. 0) then
       print *,'inconsistent Twall mesh max=',d,
     +         'instead of =',Tdmax_pt
       goto 3
      endif

      d=0.d0
      read(ican,106) atxt,d !'matter atom mass, amu=',mat_mas(imatter)
106   format(A36,d25.15)
      if(d .lt. 0.d0) then
       ierr=11
      else 
       if(d .gt. 1.e3) then
       ierr=11
       endif
      endif
      if(ierr .gt. 0) then
       print *,'improper matter mass:',d
      else
       mat_mas(imatter)=d
      endif

      d=0.d0
      read(ican,106) atxt,d !'sublimation heat, eV=',mat_wsub(imatter)
      if(d .lt. 0.d0) then
       ierr=11
      else 
       if(d .gt. 30.d0) then
       ierr=11
       endif
      endif
      if(ierr .gt. 0) then
       print *,'improper sublimation heat:',d
      else
       mat_wsub(imatter)=d
      endif

      d=0.d0
      read(ican,106) atxt,d !'lattice constant, Angstroms=', mat_a0(imatter)
      if(d .lt. 0.d0) then
       ierr=11
      else 
       if(d .gt. 30.d0) then
       ierr=11
       endif
      endif
      if(ierr .gt. 0) then
       print *,'improper lattice constant:',d
      else
       mat_a0(imatter)=d
      endif
         
      d=0.d0
      read(ican,106) atxt,d !'melting temperature, Kelvin=', mat_Tmelt(imatter)
      if(d .lt. 0.d0) then
       ierr=11
      else 
       if(d .gt. 1.d4) then
       ierr=11
       endif
      endif
      if(ierr .gt. 0) then
       print *,'improper melting temperature:',d
      else
       mat_Tmelt(imatter)=d
      endif

      d=0.d0
      read(ican,106) atxt,d !'melting enthalpy, eV=', mat_Hmelt(imatter)
      if(d .lt. 0.d0) then
      ierr=11
      else 
       if(d .gt. 30.d0) then
       ierr=11
       endif
      endif
      if(ierr .gt. 0) then
      print *,'improper melting enthalpy:',d
      else
      mat_Hmelt(imatter)=d
      endif

      d=0.d0
      read(ican,106) atxt,d !'sublimation temperature,Kelvin=',mat_Tsub(imatter)
      if(d .lt. 0.d0) then
       ierr=11
      else 
       if(d .gt. 1.d4) then
       ierr=11
       endif
      endif
      if(ierr .gt. 0) then
       print *,'improper sublimation temperature:',d
      else
       mat_Tsub(imatter)=d
      endif
      
      read(ican,106) atxt,d !'work function, eV=',mat_Wf(imatter)
      if(d .lt. 0.d0) then
       ierr=11
      else 
       if(d .gt. 30.d0) then
       ierr=11
       endif
      endif
      if(ierr .gt. 0) then
       print *,'improper work function:',d
      else
       mat_Wf(imatter)=d
      endif


      read(ican,107) atxt !'ro, matter mass density,g/cm3'
      read(ican,107) atxt !'cp, heat capacity,J/ g grad K'
      read(ican,107) atxt !'kp, heat conductivity,W/ m grad K'
      read(ican,107) atxt !'ep, Black-Body emissivity coefficient'
      read(ican,107) atxt !'hs, thermal emission heat,eV'
      read(ican,107) atxt !'gs, (LOG)thermal emission flux,10^20 atoms/cm^2/sec'
      read(ican,107) atxt !'pb, Coefficient of thermal emission'
      read(ican,107) atxt !'T(K)   ro   cp   kp   ep   hs   gs  pbtd:'
107   format(A52)

      do 10 i=1,Ntd_pt
      read(ican,110)  d,
     +                psi_ro(i,imatter),
     +                psi_cp(i,imatter),
     +                psi_kp(i,imatter),
     +                psi_ep(i,imatter),
     +                psi_hs(i,imatter),
     +                psi_gs(i,imatter),
     +                psi_pbtd(i,imatter)
110   format(E14.3,7(d25.15))
10    continue

      do 12 i=1,Nprj_pt
      j=I_pwi(i,imatter)
      Iprjj(i)=j
12    continue

c
c--Rn
c
      do 40 i=1,Nprj_pt
       Iprjz(i)=0
       Ipjj(i)=0
       Iprjj(i)=0
40    continue
      read(ican,107) atxt !'<Rn>Mx, particle reflection coefficient'
      read(ican,108) atxt,(Iprjz(i),i=1,Nprj_pt) 
                          !'Tplasma(eV),<Rn>Mx for projectiles:',(Iprjz(i),i=1,Nprj_pt)
108   format(A52,10I4)
      do 43 i=1,Nprj_pt
       j=Iprjz(i)
       if(j .lt. 0) then
        ierr=12
       else 
        if(j .gt. 92) then
         ierr=12
        endif
       endif
       if(ierr .gt. 0) then
        print *,'<Rn>Mx incorrect projectile=',j,' index=',i
       endif
43    continue
      if(ierr .gt. 0) goto 3
      do 41 i=1,Nprj_pt
       k=IZprjctl(i)
       kk=0
       do 42 j=1,Nprj_pt
        if(k .eq. Iprjz(j)) then
         kk=kk+1
         Iprjj(i)=j
        endif
42     continue
       if(kk .lt. 1) then
        ierr=4
        print *,'<Rn>Mx missing projectile=',k      
       else
        if(kk .gt. 1) then
         ierr=4
         print *,'<Rn>Mx repeated projectile=',k
        endif
       endif
41    continue
      if(ierr .gt. 0) goto 3
      do 44 i=1,Nprj_pt
      Ipjj(i)=I_pwi(i,imatter)
44    continue

      do 45 i=1,Nt_pt

      read(ican,111) d,(Cjj(j),j=1,Nprj_pt)       
111   format(d15.3,10(d25.15))

      do 46 j=1,Nprj_pt
      k=Ipjj(j)
      kk=Iprjj(j)
      rn_mx(i,k)=Cjj(kk)
46    continue

45    continue

c
c--Re
c
      do 50 i=1,Nprj_pt
50    continue
      read(ican,107) atxt !'<Re>Mx, particle reflection coefficient'
      read(ican,108) atxt,(Iprjz(i),i=1,Nprj_pt) 
                       !'Tplasma(eV),<Re>Mx for projectiles:',(Iprjz(i),i=1,Nprj_pt)
      do 53 i=1,Nprj_pt
      j=Iprjz(i)
      if(j .lt. 0) then
      ierr=12
      else 
       if(j .gt. 92) then
       ierr=12
       endif
      endif
      if(ierr .gt. 0) then
      print *,'<Re>Mx incorrect projectile=',j,' index=',i
      endif
53    continue
      if(ierr .gt. 0) goto 3
      do 51 i=1,Nprj_pt
      k=IZprjctl(i)
      kk=0
      do 52 j=1,Nprj_pt
      if(k .eq. Iprjz(j)) then
      kk=kk+1
      Iprjj(i)=j
      endif
52    continue
      if(kk .lt. 1) then
      ierr=4
      print *,'<Re>Mx missing projectile=',k      
      else
       if(kk .gt. 1) then
       ierr=4
       print *,'<Re>Mx repeated projectile=',k
       endif
      endif
51    continue
      if(ierr .gt. 0) goto 3
      do 54 i=1,Nprj_pt
      Ipjj(i)=I_pwi(i,imatter)
54    continue

      do 55 i=1,Nt_pt

      read(ican,111) d,(Cjj(j),j=1,Nprj_pt) 

      do 56 j=1,Nprj_pt
      k=Ipjj(j)
      kk=Iprjj(j)
      re_mx(i,k)=Cjj(kk)
56    continue

55    continue

c
c--Yphys
c
      do 60 i=1,Nprj_pt
       Iprjz(i)=0
       Ipjj(i)=0
       Iprjj(i)=0
60    continue
      read(ican,107) atxt !'<Yph>Mx, physical sputtering yield'
      read(ican,108) atxt,(Iprjz(i),i=1,Nprj_pt) !'Tplasma(eV),<Yph>Mx for projectiles:',(Iprjz(i),i=1,Nprj_pt)
      do 63 i=1,Nprj_pt
       j=Iprjz(i)
       if(j .lt. 0) then
        ierr=12
       else 
        if(j .gt. 92) then
         ierr=12
        endif
       endif
       if(ierr .gt. 0) then
        print *,'<Yph>Mx incorrect projectile=',j,' index=',i
       endif
63    continue
      if(ierr .gt. 0) goto 3
      do 61 i=1,Nprj_pt
       k=IZprjctl(i)
       kk=0
       do 62 j=1,Nprj_pt
        if(k .eq. Iprjz(j)) then
        kk=kk+1
        Iprjj(i)=j
       endif
62    continue
      if(kk .lt. 1) then
       ierr=4
       print *,'<Yph>Mx missing projectile=',k      
      else
       if(kk .gt. 1) then
        ierr=4
        print *,'<Yph>Mx repeated projectile=',k
       endif
      endif
61    continue
      if(ierr .gt. 0) goto 3
      do 64 i=1,Nprj_pt
       Ipjj(i)=I_pwi(i,imatter)
64    continue

      do 65 i=1,Nt_pt
       read(ican,113) d,(Cjj(j),C2jj(j),j=1,Nprj_pt) 
113    format(d15.3,20(d25.15))
       do 66 j=1,Nprj_pt
        k=Ipjj(j)
        kk=Iprjj(j)
        ypsp_mx(i,k)=Cjj(kk)
        epsp_mx(i,k)=C2jj(kk)
66     continue
65    continue

c
c--Ychem
c
      do 70 i=1,Nprj_pt
      Iprjz(i)=0
      Ipjj(i)=0
      Iprjj(i)=0
70    continue
      read(ican,107) atxt !'<Ych>Mx, chemical sputtering yield'
      read(ican,108) atxt,(Iprjz(i),i=1,Nprj_pt) 
               !'Tplasma(eV),Twall(eV),<Ych>Mx for projectiles:',(Iprjz(i),i=1,Nprj_pt)
      do 73 i=1,Nprj_pt
      j=Iprjz(i)
      if(j .lt. 0) then
      ierr=12
      else 
       if(j .gt. 92) then
       ierr=12
       endif
      endif
      if(ierr .gt. 0) then
      print *,'<Ych>Mx incorrect projectile=',j,' index=',i
      endif
73    continue
      if(ierr .gt. 0) goto 3
      do 71 i=1,Nprj_pt
      k=IZprjctl(i)
      kk=0
      do 72 j=1,Nprj_pt
      if(k .eq. Iprjz(j)) then
      kk=kk+1
      Iprjj(i)=j
      endif
72    continue
      if(kk .lt. 1) then
      ierr=4
      print *,'<Ych>Mx missing projectile=',k      
      else
       if(kk .gt. 1) then
       ierr=4
       print *,'<Ych>Mx repeated projectile=',k
       endif
      endif
71    continue
      if(ierr .gt. 0) goto 3
      do 74 i=1,Nprj_pt
      Ipjj(i)=I_pwi(i,imatter)
74    continue

      do 77 ii=1,Ntd_pt
      do 75 i=1,Nt_pt

      read(ican,112) d,d,(Cjj(j),j=1,Nprj_pt) 
112   format(d15.3,d15.3,10(d25.15))

      do 76 j=1,Nprj_pt
      k=Ipjj(j)
      kk=Iprjj(j)
      ycsp_mx(i,ii,k)=Cjj(kk)
76    continue

75    continue
77    continue

c
c--Yres
c
      do 80 i=1,Nprj_pt
       Iprjz(i)=0
       Ipjj(i)=0
       Iprjj(i)=0
80    continue
      read(ican,107) atxt !'<Yres>Mx, RES, radiation enhanced sublimation'
      read(ican,108) atxt,(Iprjz(i),i=1,Nprj_pt) 
            !'Tplasma(eV),Twall(eV),<Yres>Mx for projectiles:',(Iprjz(i),i=1,Nprj_pt)
      do 83 i=1,Nprj_pt
       j=Iprjz(i)
       if(j .lt. 0) then
        ierr=12
       else 
        if(j .gt. 92) then
         ierr=12
        endif
       endif
       if(ierr .gt. 0) then
        print *,'<Yres>Mx incorrect projectile=',j,' index=',i
       endif
83    continue
      if(ierr .gt. 0) goto 3
      do 81 i=1,Nprj_pt
       k=IZprjctl(i)
       kk=0
       do 82 j=1,Nprj_pt
        if(k .eq. Iprjz(j)) then
         kk=kk+1
         Iprjj(i)=j
        endif
82     continue
       if(kk .lt. 1) then
        ierr=4
        print *,'<Yres>Mx missing projectile=',k      
       else
        if(kk .gt. 1) then
         ierr=4
         print *,'<Yres>Mx repeated projectile=',k
        endif
       endif
81    continue
      if(ierr .gt. 0) goto 3
      do 84 i=1,Nprj_pt
       Ipjj(i)=I_pwi(i,imatter)
84    continue

      do 87 ii=1,Ntd_pt
      do 85 i=1,Nt_pt

      read(ican,112) d,d,(Cjj(j),j=1,Nprj_pt) 

      do 86 j=1,Nprj_pt
      k=Ipjj(j)
      kk=Iprjj(j)
      yres_mx(i,ii,k)=Cjj(kk)
86    continue

85    continue
87    continue

c
c--Ysee
c
      read(ican,207) atxt !'<Ysee>Mx, SEE, secondary electron-electron emission'
      read(ican,207) atxt !'Tplasma(eV),Twall(eV),<Ysee>Mx'

      do 287 ii=1,Nphi_pt
       do 285 i=1,Nt_pt
        read(ican,212) d,d,c
        dse_mx(i,ii,imatter)=c
285    continue
287   continue
207   format (A52)
212   format (d15.3,d15.3,d25.15)

c
c--Emissivity
c
      read(ican,307) atxt !'<Eps>Mx, EMS, thermal emissivity from Mie'
      read(ican,307) atxt !'Radius,Tdust(K),<Ytheps>Mx'

      do 319 ii=1,Ntde_pt
       do 320 i=1,Nrd_pt
        read(ican,312) d,d,c
        dtheps_mx(i,ii,imatter)=c
320    continue
319   continue
307   format (A52)
312   format (d15.3,d15.3,d25.15)


      read(ican,107) atxt !'That's all folks'

      goto 1

2     ierr=1
      print *, 'psitab_r:error in imatter=',imatter
      goto 220
      
3     ierr=ierr+1000
      print *, 'psitab_r:error=',ierr,' in reading file=',
     +          Namfpsi(imatter)

1     close (ican)
220   continue

      end
