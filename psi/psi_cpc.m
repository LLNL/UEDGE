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
      double precision function psicpc(iz,tk)
      implicit none
c -------------------------
c    termal capacity
c	  T [grad K]
c	 Cp [J/ g grad K]
c -------------------------
      
      Use (psitab)
      Use (psitab_s)
      Use (psitabint)
      
      integer iz
      double precision tk
      double precision psicpcbe,psicpcb,psicpcc
      double precision psicpcfe,psicpcmo,psicpcw
      double precision psicpcbef,psicpcbf,psicpccf
      double precision psicpcfef,psicpcmof,psicpcwf
      double precision psicpcli
      
c --- Li      
      if(iz .eq. 3)then
       psicpc=psicpcli(tk)

c --- Be
      elseif(iz .eq. 4) then
       psicpc=psicpcbe(tk)

c --- B      
      elseif(iz .eq. 5) then
       psicpc=psicpcb(tk)

c --- C      
      elseif(iz .eq. 6) then
       psicpc=psicpcc(tk)

c --- Fe      
      elseif(iz .eq. 26) then
       psicpc=psicpcfe(tk)

c --- Mo      
      elseif(iz .eq. 42) then
       psicpc=psicpcmo(tk)

c --- W      
      elseif(iz .eq. 74) then
       psicpc=psicpcw(tk)

c --- DUMMY
      else
       psicpc=0.d0
      endif 

      return
      end
c-----------------------------------------------------------------------

      double precision function psikapc(iz,tk)
      implicit none
      
      integer iz
      double precision tk
c  ------------------------------
c     parsing elements (currently C and W)
c	  T [grad K]
c     kappa [W/ m grad K]
c  ------------------------------
      double precision psikapcc,psikapcw

      if    (iz .eq. 3) then  ! Li not correct
       psikapc=0.d0
      elseif(iz .eq. 4) then  ! Be not correct
       psikapc=0.d0
      elseif(iz .eq. 5) then  ! B not correct
       psikapc=0.d0
      elseif(iz .eq. 6) then  ! C
       psikapc=psikapcc(tk)
      elseif(iz .eq. 26) then ! Fe not correct
       psikapc=0.d0
      elseif(iz .eq. 42) then ! Mo not correct
       psikapc=0.d0
      elseif(iz .eq. 74) then ! W 
       psikapc=psikapcw(tk)
      else                    ! DUMMY
       psikapc=0.d0
      endif 

      return
      end
c-----------------------------------------------------------------------

      double precision function psikapcc(tk)
      implicit none
c  ------------------------------
c    thermal conductivity for C 
c         T [grad K]
c     kappa [W/ m grad K]
c  ------------------------------
      integer nt
      real tkabs
      real st(8),sl(8)
      double precision tk
      double precision t
      double precision psilint

      data st
     ,    /   20.0 ,  200.0 ,  600.0 ,  900.0 , 1200.0 ,
     ,      1500.0 , 1800.0 , 2500.0 /
      data sl
     ,    /  114.0 ,   86.0 ,   61.0 ,   47.0 ,   40.0 ,
     ,        34.0 ,   25.0 ,   10.0 /

      data tkabs /20.0/
      data nt    /8   /

      t = tk - tkabs;
      psikapcc = psilint(sl,st, nt,t)
      return
      end

c---------------------------------------------
      double precision function psikapcw(tk)
      implicit none
      double precision tk
c  ------------------------------
c    thermal conductivity for  W
c     T     [grad K]
c     kappa [W/ m grad K]
c  ------------------------------
      integer nt
      real tkabs
      real st(12),sl(12)
      double precision t
      double precision psilint

      data st
     ,    /  273.0 ,  573.0 ,  773.0 , 1273.0 , 1773.0 ,
     ,      2273.0 , 2773.0 , 3273.0 , 3673.0 , 3723.0 ,
     ,      3773.0 , 4273.0 /
      data sl
     ,    /  174.0 ,  139.0 ,  127.0 ,  111.0 ,  103.0 ,
     ,        96.3 ,   92.1 ,   91.1 ,   88.9 ,   70.0 ,
     ,        71.1 ,   75.3 /
      data tkabs /0.0/
      data nt    /12 /

      t = tk - tkabs;
      psikapcw = psilint(sl,st, nt, t)
      return
      
      end
c----------------------------------------------------------------

      double precision function psicpcli(tk)
      implicit none
      double precision tk
cc  ------------------------------
cc    thermal capacity for Li(ref)
cc         T [grad K] in JANAF
cc        data Cp [J/ mol grad K]-->J/g/K
cc  ------------------------------
      integer nt
      real tkabs
      real st(22),sl(22)
      double precision t
      double precision psilint
      double precision psimassa
      data st
     ,    /  298.15 ,  300.0  ,  350.0  ,  400.0  ,  450.0  ,
     ,       453.690,  500.0  ,  600.0  ,  700.0  ,  800.0  ,
     ,       900.0  , 1000.0  , 1100.0  , 1200.0  , 1300.0  ,
     ,      1400.0  , 1500.0  , 1600.0  , 1620.120, 1700.0  ,
     ,      1800.0  , 6000.0  /
      data sl
     ,    /   31.300,   31.284,   31.000,   30.711,   30.414,
     ,        30.392,   30.125,   29.539,   28.987,   28.937,
     ,        28.886,   28.836,   28.786,   28.736,   28.702,
     ,        28.619,   28.535,   28.451,   28.431,   28.368,
     ,        28.326,   28.326/
      data tkabs /0.0/
      data nt    /22 /

      t = tk - tkabs
      psicpcli = psilint(sl,st, nt,t)/psimassa(3) ! J/(K mol)->J/(g K)
      return

      end

c----------------------------------------------------------------

      double precision function psicpcbe(tk)
      implicit none
c  ---------------------------------------------------------
c    thermal capacity for Be
c         T [grad K] in JANAF Table 
c        data Cp [J/ mol grad K]-->J/g/K
c  ---------------------------------------------------------
      integer nt
      real tkabs
      real st(40),sl(40)
      double precision tk
      double precision t
      double precision psimassa
      double precision psilint

      data st
     ,    /  298.15 ,  300.0  ,  350.0  ,  400.0  ,  450.0  ,
     ,       500.0  ,  600.0  ,  700.0  ,  800.0  ,  900.0  ,
     ,      1000.0  , 1100.0  , 1150.000, 1150.001, 1200.0  ,
     ,      1300.0  , 1400.0  , 1500.0  , 1560.0  , 1600.0  ,
     ,      1700.0  , 1800.0  , 1900.0  , 2000.0  , 2100.0  ,
     ,      2200.0  , 2300.0  , 2400.0  , 2500.0  , 2600.0  ,
     ,      2700.0  , 2741.437, 2800.0  , 2900.0  , 3000.0  ,
     ,      3100.0  , 3200.0  , 3300.0  , 3400.0  , 3500.0  /
      data sl
     ,    /  16.380 ,  16.472 ,  18.520 ,  19.965 ,  21.061 ,
     ,       21.943 ,  23.336 ,  24.463 ,  25.458 ,  26.384 ,
     ,       27.274 ,  28.147 ,  28.582 ,  27.907 ,  28.014 ,
     ,       28.229 ,  28.444 ,  28.659 ,  28.788 ,  28.874 ,
     ,       29.089 ,  29.304 ,  29.519 ,  29.734 ,  29.949 ,
     ,       30.164 ,  30.379 ,  30.594 ,  30.809 ,  31.024 ,
     ,       31.239 ,  31.328 ,  31.454 ,  31.669 ,  31.884 ,
     ,       32.099 ,  32.314 ,  32.529 ,  32.744 ,  32.959 /

      data tkabs /0.0/
      data nt    /40 /

      t = tk - tkabs

      psicpcbe = psilint(sl,st, nt, t)/psimassa(4) ! J/(K mol)->J/(g K)
      return
		
      end
					      
c----------------------------------------------------------------
							      
      double precision function psicpcb(tk)
      implicit none
      double precision tk
c
c  ------------------------------
c    thermal capacity for B
c         T [grad K]
c        data Cp [J/ mol grad K]-->J/g/K
c  ------------------------------
c
      double precision psilint
      real tkabs
      integer nt
      real st(21),sl(21)
      double precision t
      double precision psimassa
      external psimassa
											      
      data st
     ,    /  298.15 ,  300.0  ,  350.0  ,  400.0  ,  450.0  ,
     ,       500.0  ,  600.0  ,  700.0  ,  800.0  ,  900.0  ,
     ,      1000.0  , 1100.0  , 1200.0  , 1300.0  , 1400.0  ,
     ,      1500.0  , 1600.0  , 1700.0  , 1750.000, 1750.001,
     ,      6000.0  /
      data sl
     ,    /  11.315 ,  11.405 ,  13.654 ,  15.693 ,  17.361 ,
     ,       18.722 ,  20.778 ,  22.249 ,  23.361 ,  24.245 ,
     ,       24.978 ,  25.606 ,  26.161 ,  26.663 ,  27.125 ,
     ,       27.557 ,  27.966 ,  28.356 ,  28.546 ,  31.750 ,
     ,       31.750 /
      data tkabs /0.0/
      data nt    /21 /
      
      t = tk - tkabs
      psicpcb = psilint(sl,st, nt, t)/psimassa(5) ! J/(K mol)->J/(g K)
      return
      end
c----------------------------------------------------------------


      double precision function psicpcc(tk)
      implicit none
c  ------------------------------
c    termal capacity for C(ref) JANAF
c	  T [grad K]
c    data Cp [J/ mol grad K]--> J/g/K
c  ------------------------------
      integer nt
      real tkabs
      real st(64),sl(64)
      double precision tk
      double precision t
      double precision psimassa
      double precision psilint

      data st
     ,    /  100.0 ,  200.0 ,  250.0 ,  298.15,  300.0 ,
     ,       350.0 ,  400.0 ,  450.0 ,  500.0 ,  600.0 ,
     ,       700.0 ,  800.0 ,  900.0 , 1000.0 , 1100.0 ,
     ,      1200.0 , 1300.0 , 1400.0 , 1500.0 , 1600.0 ,
     ,      1700.0 , 1800.0 , 1900.0 , 2000.0 , 2100.0 ,
     ,      2200.0 , 2300.0 , 2400.0 , 2500.0 , 2600.0 ,
     ,      2700.0 , 2800.0 , 2900.0 , 3000.0 , 3100.0 ,
     ,      3200.0 , 3300.0 , 3400.0 , 3500.0 , 3600.0 ,
     ,      3700.0 , 3800.0 , 3900.0 , 4000.0 , 4100.0 ,
     ,      4200.0 , 4300.0 , 4400.0 , 4500.0 , 4600.0 ,
     ,      4700.0 , 4800.0 , 4900.0 , 5000.0 , 5100.0 ,
     ,      5200.0 , 5300.0 , 5400.0 , 5500.0 , 5600.0 ,
     ,      5700.0 , 5800.0 , 5900.0 , 6000.0 /

      data sl
     ,    /  1.674 ,  5.006 ,  6.816 ,  8.517 ,  8.581 ,
     ,      10.241 , 11.817 , 13.289 , 14.623 , 16.844 ,
     ,      18.537 , 19.827 , 20.824 , 21.610 , 22.244 ,
     ,      22.766 , 23.204 , 23.578 , 23.904 , 24.191 ,
     ,      24.448 , 24.681 , 24.895 , 25.094 , 25.278 ,
     ,      25.453 , 25.618 , 25.775 , 25.926 , 26.071 ,
     ,      26.212 , 26.348 , 26.481 , 26.611 , 26.738 ,
     ,      26.863 , 26.986 , 27.106 , 27.225 , 27.342 ,
     ,      27.459 , 27.574 , 27.688 , 27.801 , 27.913 ,
     ,      28.024 , 28.134 , 28.245 , 28.354 , 28.462 ,
     ,      28.570 , 28.678 , 28.785 , 28.893 , 28.999 ,
     ,      29.106 , 29.211 , 29.317 , 29.422 , 29.528 ,
     ,      29.632 , 29.737 , 29.842 , 29.946 /
      data tkabs /0.0/
      data nt    /64 /

      t = tk - tkabs
      psicpcc = psilint(sl,st, nt, t)/psimassa(6)
      return
      end
c----------------------------------------------------------------

      double precision function psicpcfe(tk)
      implicit none
      double precision tk
c  ------------------------------
c    thermal capacity for Fe(ref)
c         T [grad K] in JANAF Table 
c        data Cp [J/ mol grad K]-->J/g/K
c  ------------------------------
      double precision psilint
      real tkabs
      integer nt
      real st(14),sl(14)
      double precision t
      double precision psimassa
      external psimassa

      data st
     ,    /  298.15 ,  300.0  ,  350.0  ,  400.0  ,  450.0  ,
     ,       500.0  ,  600.0  ,  700.0  ,  800.0  ,  900.0  ,
     ,      1000.0  , 1100.0  , 1200.000,  6000.0 /
      data sl
     ,    /  25.096 ,  25.142 ,  26.254 ,  27.384 ,  28.542 ,
     ,       29.702 ,  32.033 ,  34.363 ,  36.694 ,  39.024 ,
     ,       41.355 ,  43.689 ,  46.024 ,  46.024 /

      data tkabs /0.0/
      data nt    /14 /

      t = tk - tkabs
      psicpcfe = psilint(sl,st, nt, t)/psimassa(26) ! J/(K mol)->J/(g K)

      return
      end
c----------------------------------------------------------------


      double precision function psicpcmo(tk)
      implicit none
c  ------------------------------
c    thermal capacity for Mo(cr, 100-2150K), W(l,2150-4952K)
c         T [grad K] in JANAF Table
c        data Cp [J/ mol grad K]-->J/g/K
c  ------------------------------
      integer nt
      real tkabs
      real st(25),sl(25)
      double precision tk
      double precision t
      double precision psilint
      double precision psimassa

      data st
     ,    /  298.15 ,  300.0  ,  350.0  ,  400.0  ,  450.0  ,
     ,       500.0  ,  600.0  ,  700.0  ,  800.0  ,  900.0  ,
     ,      1000.0  , 1100.0  , 1200.0  , 1300.0  , 1400.0  ,
     ,      1500.0  , 1600.0  , 1700.0  , 1800.0  , 1900.0  ,
     ,      2000.0  , 2100.0  , 2150.000, 2150.001, 6000.0  / 
      data sl
     ,    /  23.932 ,  23.958 ,  24.580 ,  25.080 ,  25.500 ,
     ,       25.850 ,  26.460 ,  26.980 ,  27.440 ,  27.890 ,
     ,       28.370 ,  28.900 ,  29.490 ,  30.140 ,  30.860 ,
     ,       31.650 ,  32.500 ,  33.420 ,  34.420 ,  35.490 ,
     ,       36.650 ,  37.900 ,  38.550 ,  37.656 ,  37.656 /

      data tkabs /0.0/
      data nt    /25 /

      t = tk - tkabs
      psicpcmo = psilint(sl,st, nt, t)/psimassa(42) ! J/(K mol)->J/(g K)

      return
      end
c-----------------------------------------------------------------------


      double precision function psicpcw(tk)
      implicit none
c  ------------------------------
c    thermal capacity for W(cr,100-2450K), W(l,2450-6000K) 
c         T [grad K] JANAF Table
c        data Cp [J/ mol grad K]-->J/g/K
c  ------------------------------
      integer nt
      real tkabs
      real st(28),sl(28)
      double precision tk
      double precision t
      double precision psimassa
      double precision psilint

      data st
     ,    /  298.15 ,  300.0  ,  350.0  ,  400.0  ,  450.0  ,
     ,       500.0  ,  600.0  ,  700.0  ,  800.0  ,  900.0  ,
     ,      1000.0  , 1100.0  , 1200.0  , 1300.0  , 1400.0  ,
     ,      1500.0  , 1600.0  , 1700.0  , 1800.0  , 1900.0  ,
     ,      2000.0  , 2100.0  , 2200.0  , 2300.0  , 2400.0  ,
     ,      2450.000, 2450.001, 6000.0  /
      data sl
     ,    /  24.296 ,  24.310 ,  24.642 ,  24.928 ,  25.144 ,
     ,       25.359 ,  25.790 ,  26.229 ,  26.669 ,  27.112 ,
     ,       27.564 ,  28.017 ,  28.472 ,  28.930 ,  29.393 ,
     ,       29.862 ,  30.334 ,  30.806 ,  31.284 ,  31.768 ,
     ,       32.254 ,  32.745 ,  33.238 ,  33.735 ,  34.233 ,
     ,       34.483 ,  35.564 ,  35.564 /


      data tkabs /0.0/
      data nt    /28 /

      t = tk - tkabs
      psicpcw = psilint(sl,st, nt, t)/psimassa(74) ! J/(K mol)->J/(g K)
      return
      end
