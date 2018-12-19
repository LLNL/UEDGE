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
c-------------------------------------------
c  Formular for Physical sputtering yield
c       given by
c
c    J. Bohdansky
c    Nuclear Instrumental Methods,
c    Sec B, 2 (1984) 587
c--------------------------------------------
c
      integer function psi_bgh_OLD(izp,mp,izt,mt)
      implicit none
      integer izp,izt
      double precision mp,mt
c
c-------------------------------------
c      table from:
c       W Eckstein J.Bohdansky J Roth
c       Physical sputtering
c       Atomic and plasma-material interaction data for fusion,
c       Supplement to the journal Nucl.Fus
c       V1  1991 p51
c-------------------------------------
c
      common/psibghc/Q_,El_,Eth_,
     *               Z1_,Z2_,M1_,M2_
      double precision Q_,El_,Eth_
      double precision Z1_,Z2_,M1_,M2_

c -- mask

c         h      d     t     he     c     o    self
c--- Be  
c--- B    
c--- C   
c--- Al  
c--- Ti  
c--- Fe  
c--- Ni  
c--- Cu    
c--- Mo 
c--- W 

      real  BGH_wth(60)
      data BGH_wth/ 
     ,    20.0,  9.0,  21.0, 30.0, 40.0, 70.0,  25.0, 
     ,    20.0, 20.0,  21.0, 30.0, 40.0,        35.0,
     ,    35.0, 30.0,  30.0, 29.0,              42.0,
     ,    53.0, 24.0,  22.0, 20.0,              42.0,
     ,    80.0, 50.0,  40.0, 25.0, 35.0,        40.0,
     ,    64.0, 44.0,  40.0, 33.0, 35.0,        40.0,
     ,    50.0, 31.0,  25.0, 22.0, 28.0, 48.0,  34.4,
     ,    55.0, 34.0,  30.0, 19.0,              36.0,
     ,   199.0, 90.0,  70.0, 46.0, 55.0,        64.0,
     ,   443.0,220.0, 140.0,110.0, 80.0, 40.0,  65.0/


      real  BGH_q(60)
      data BGH_q/
     ,     0.10,  0.3,  0.24, 0.59,  1.6,  1.3,   1.4,
     ,     0.10, 0.15,  0.24, 0.59,  1.6,        0.94,
     ,    0.035, 0.10,  0.20, 0.32,               1.5,
     ,    0.043,0.093,   0.2, 0.34,               5.4,
     ,    0.017,0.055,   0.1,0.125,  1.0,         3.7,
     ,    0.042, 0.13,  0.21, 0.44,  3.2,        13.0,
     ,    0.041,0.115,  0.22, 0.45,  2.5, 1.55,  11.7,
     ,      0.1, 0.23,  0.35, 0.81,              17.0,
     ,    0.007,0.023, 0.045, 0.12, 0.93,        18.0,
     ,    0.007,0.019, 0.038,0.106, 0.93,  2.2,  20.0/

      integer, parameter :: n=60
      integer BGH_cod(60)
      data BGH_cod/
     ,    104, 9304,  9404,  204,  604,  804,   404,
     ,    105, 9305,  9405,  205,  605,         505,
     ,    106, 9306,  9406,  206,               606,
     ,    113, 9313,  9413,  213,              1313,
     ,    122, 9322,  9422,  222,  622,        2222,
     ,    126, 9326,  9426,  226,  626,        2626,
     ,    128, 9328,  9428,  228,  628,  828,  2828,
     ,    129, 9329,  9429,  229,              2929,
     ,    142, 9342,  9442,  242,  642,        4242,
     ,    174, 9374,  9474,  274,  674,  874,  7474/

      integer, parameter :: BGH_Nt = 10
      integer BGH_tcod(10)
      data BGH_tcod/4, 5, 6, 13, 22, 26, 28, 29, 42, 74/

      integer, parameter :: BGH_Np = 12
      integer BGH_pcod(12)
c       he  be  b  c  o  al  ti  fe  ni  cu  mo  w
      data BGH_pcod/
     ,	 2,  4, 5, 6, 8, 13, 22, 26, 28, 29, 42, 74/

      double precision z1,z2,m1,m2
      double precision s,f,ff
      integer iz,l


      do 10 iz=1,BGH_Nt
       if(BGH_tcod(iz) .eq. izt) then
       l=iz
       goto 2
       endif
  10  continue
      psi_bgh_OLD=1
      return
   2  continue

      if (izp .eq. 1) then
       l = 0
       if (mp .lt. 0.9) then
        l = 0
       elseif (mp .lt. 1.1) then
        l = 1 ! h 
       elseif (mp .lt. 2.1) then
        l = 93 ! d
       elseif (mp .lt. 3.1) then
        if (mp .gt. 2.9) then
         l = 94 ! t
        else
         l = 0
        endif
       else
       l = 0
       endif
       if(l .eq. 0) then
         psi_bgh_OLD=3
         return
       endif
      else
       do 20 iz=1,BGH_Nt
        if(BGH_pcod(iz) .eq. izp) then
         l=iz
         goto 3
        endif
20     continue
       psi_bgh_OLD=2
       return
3      continue
      endif

      l = l*100 + izt
      do 30 iz=1,n
       if(BGH_cod(iz) .eq. l) goto 4
30    continue
      psi_bgh_OLD=4
      return

4     continue

      z1 = dble(izp)
      m1 = mp
      z2 = dble(izt)
      m2 = mt

      Z1_ = z1
      Z2_ = z2
      M1_ = m1
      M2_ = m2

      s = 2.d0/3.d0
      f = z1**s
      ff= z2**s
      s = dsqrt(f+ff)
      f = 30.74 * (m1+m2)/m2 * z1*z2 * s

      El_ = f
      Eth_= dble(BGH_wth(iz))
      Q_  = dble(BGH_q(iz))

      psi_bgh_OLD=0
      return

      end


      double precision function psi_bgh_as(m1,m2)
      implicit none
      double precision m1,m2
c
c------------------------------------------
c surface correction factor for sputtering
c
c   J. Bohdansky
c   Nuclear Instrumental Methods,
c   Sec B, 2 (1984) 587
c-------------------------------------------
c
      double precision, parameter :: Ka_surf=0.2d0
      double precision, parameter :: Ka_surf_min=0.2d0

      double precision as,al,sp

      as = 1.d0 / (1.d0 + Ka_surf * m2/m1)

      sp = 2.d0/3.d0;
      al = 0.3d0 * (m2/m1)**sp
      if(al .lt. Ka_surf_min) then
      al = Ka_surf_min
      endif
      as = as * al

      psi_bgh_as=as
      return
      end


      double precision function psi_bgh_wth(us,m1,m2)
      implicit none
      double precision us,m1,m2

c
c----------------------------------
c sputtering threshold energy
c
c    J. Bohdansky
c    Nuclear Instrumental Methods,
c    Sec B, 2 (1984) 587
c----------------------------------
c
      double precision, parameter :: GLIGHT=0.2d0

      double precision seps 
      double precision sp,Gam
      double precision f,ff

      data seps/1.d-8/

      if(m1 .le. GLIGHT*m2) then  
       ff = m1+m2
        f = ff*ff
       Gam= 4.d0*m1*m2/f
       ff = Gam * (1.d0-Gam)
       if(ff .lt. seps) ff=seps
        f = 1.d0/ff   
      else   
       sp = 0.4d0
       ff = m1/m2
        f = 8.d0 * ff**sp
      endif

      psi_bgh_wth=us*f
      return

      end


      double precision function psi_bgh_fy(eps)
      implicit none
      double precision eps
c
c----------------------------------
c   J. Bohdansky
c   Nuclear Instrumental Methods,
c   Sec B, 2 (1984) 587
c----------------------------------
c
      double precision seps
      double precision x,f,ae,ff

      data seps/1.d-18/

      psi_bgh_fy=0.d0

      if(eps < 1.d0+seps) then
      return
      endif

       x = 1.d0/eps
      ff = 1.d0 - x
       f = ff*ff
      if(f < seps) then
      return
      endif

      ae = 2.d0 * dlog(x) / 3.d0;
      ff = 1.d0 - dexp(ae)
      if(ff < seps)  then
      return
      endif

      psi_bgh_fy=f*ff
      return
 
      end


      double precision function psi_bgh_Q(izp,mp,izt,mt,us)
      implicit none
      integer izp,izt
      double precision mp,mt,us

      double precision psifunsig
      double precision psi_zig_lm

      double precision z1,z2,m1,m2
      double precision sp,s,f,ff
      double precision gam,el,q

      z1 = dble(izp)
      m1 = mp;
      z2 = dble(izt)
      m2 = mt

      s = m1+m2
      f = s*s
      gam = 4.d0*m1*m2/f

      sp = 2.d0/3.d0
       f = z1**sp
      ff = z2**sp
       s = dsqrt(f+ff)

      el = 30.74 * (m1+m2)/m2 * z1*z2 * s

      sp = 0.055

      ff = el/us
      f  = gam * (ff*ff)
      s  = f**(-sp)
      q  = s * gam * ff

      f  = psifunsig(sp)
      ff = psi_zig_lm(sp)
      q  = q * f * 0.25 / ( ff * (1.d0 - 2.d0*sp) )

      psi_bgh_Q=q
      return

      end



      integer function psi_bgh_prp_OLD(izp,mp,izt,mt)
      implicit none
      integer izp,izt
      double precision mp,mt

      common/psibghc/Q_,El_,Eth_,
     *               Z1_,Z2_,M1_,M2_
      double precision Q_,El_,Eth_
      double precision Z1_,Z2_,M1_,M2_

      double precision, parameter :: SLIGHT=0.2d0

      integer psi_bgh_OLD
      integer psi_bgh_prp_,psi_yippu_prp

      integer j,ierr
      double precision s

      psi_bgh_prp_OLD=0

      j = psi_bgh_OLD(izp,mp,izt,mt)
      if (j .gt. 0) then
       ierr = psi_bgh_prp_(izp,mp,izt,mt)
       if (ierr .gt. 0) then 
       psi_bgh_prp_OLD=100+ierr
       return
       endif
      endif

      s=SLIGHT
      ierr = psi_yippu_prp(izp,mp,izt,mt, Eth_,s)
      if (ierr .gt. 0) then
      psi_bgh_prp_OLD=200+ierr
      return
      endif

      return
      end



      integer function psi_bgh_prp(izp,mp,izt,mt)
      implicit none
      integer izp,izt
      double precision mp,mt

      common/psibghc/Q_,El_,Eth_,
     *               Z1_,Z2_,M1_,M2_
      double precision Q_,El_,Eth_
      double precision Z1_,Z2_,M1_,M2_

      double precision, parameter :: SLIGHT=0.2d0

      integer psi_bgh_prp_,psi_yippu_prp

      integer ierr
      double precision s

      psi_bgh_prp=0

      ierr = psi_bgh_prp_(izp,mp,izt,mt)
      if(ierr.gt.0) then
      psi_bgh_prp=100+ierr
      return
      endif

      s=SLIGHT
      ierr = psi_yippu_prp(izp,mp,izt,mt, Eth_,s)
      if(ierr .gt. 0) then
      psi_bgh_prp=200+ierr
      return
      endif

      return
      end


      integer function psi_bgh_prp_(izp,mp,izt,mt)
      implicit none
      integer izp,izt
      double precision mp,mt

      common/psibghc/Q_,El_,Eth_,
     *               Z1_,Z2_,M1_,M2_
      double precision Q_,El_,Eth_
      double precision Z1_,Z2_,M1_,M2_

      double precision z1,z2,m1,m2
      double precision sp,s,f,ff
      double precision q,Gam,us,as

      double precision psiwsub
      double precision psi_bgh_wth,psi_bgh_as,psi_bgh_Q

      z1 = dble(izp)
      m1 = mp
      z2 = dble(izt)
      m2 = mt

      Z1_ = z1
      M1_ = m1
      Z2_ = z2
      M2_ = m2

        s = m1+m2
        f = s*s
      Gam = 4.d0*m1*m2/f

      sp = 2.d0/3.d0
       f = z1**sp
      ff = z2**sp
       s = dsqrt(f+ff)

      El_ = 30.74 * (m1+m2)/m2 * z1*z2 * s

      us  = psiwsub(izt)

      Eth_ = psi_bgh_wth(us,m1,m2)

      as = psi_bgh_as(m1,m2)

      q  = psi_bgh_Q(izp,mp,izt,mt,us)

c     Q_ = 0.042*q*as/us
      Q_ = q*as

      psi_bgh_prp_=0
      return

      end


      double precision function psi_bgh_ye(e,amu)
      implicit none
      double precision e,amu
      double precision f,ff
      
      double precision psi_bgh_yn,psi_bgh_en
      external         psi_bgh_yn,psi_bgh_en
      
      f=psi_bgh_yn(e,amu)
      ff=psi_bgh_en(e)
      psi_bgh_ye=f*ff
      return
      
      end
      


      double precision function psi_bgh_yn(e,amu)
      implicit none
      double precision e,amu

      common/psibghc/Q_,El_,Eth_,
     *               Z1_,Z2_,M1_,M2_
      double precision Q_,El_,Eth_
      double precision Z1_,Z2_,M1_,M2_

      double precision psi_bgh_fy
      double precision psi_lss_sn_TF,psi_yippu_fteta

      double precision f,ff,ft

      psi_bgh_yn=0.d0

      if(e .le. Eth_) then
       return
      endif

      ff = psi_bgh_fy (e/Eth_)
       f = psi_lss_sn_TF(e/El_)
      ft = psi_yippu_fteta(e,amu)

      psi_bgh_yn=f*ff*Q_*ft
      return

      end


      double precision function psi_bgh_get_wth(dummy)
      implicit none
      integer dummy

      common/psibghc/Q_,El_,Eth_,
     *               Z1_,Z2_,M1_,M2_
      double precision Q_,El_,Eth_
      double precision Z1_,Z2_,M1_,M2_

      psi_bgh_get_wth=Eth_
      return

      end


      double precision function psi_bgh_en(ec)
      implicit none
      double precision ec
c
c--------------------------------------------
c   averaged energy of sputtered particles
c    for Thompson spectrum
c--------------------------------------------
c
      common/psibghc/Q_,El_,Eth_,
     *               Z1_,Z2_,M1_,M2_
      double precision Q_,El_,Eth_
      double precision Z1_,Z2_,M1_,M2_

      double precision ebd,sd,d,s,redf
      double precision psiwsub
      external psiwsub
      

      data redf /0.3d0/

      ebd = psiwsub(int(Z2_))
      sd = ebd/Eth_*ec*redf
      d  = ebd/(ebd+sd)
      sd = 1.d0-d
      s  = 2.d0*(dlog(d)+sd)
      sd = sd*sd
      s  = s/sd + 1.d0
      s  = (-s) * ebd

      psi_bgh_en=s
      return

      end
