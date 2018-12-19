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
      double precision function psigsub(tk)
      implicit none
      double precision tk
c ---------------------------------
c Flux from thermal sublimation
c [Flux]=10^20 atoms/cm^2/sec
c [tk]=kelvin
c ---------------------------------
      common/psigsu/aln10,         !ln(10)
     ,              msub,          ! atom mass, amu
     ,              asub,bsub,     !coefs for sublimation equation
     ,              znum
      integer znum
      real*8 aln10
      real*8 msub
      real*8 asub,bsub
      real*8 g,f

      f=bsub-asub/tk
      g=f*aln10
      f=dexp(g)

      g=2.634776d0*f/dsqrt(msub*tk)

      psigsub=g
      return
      end


      double precision function psiesub(tk)
      implicit none
      
      double precision tk
c -------------------------------------------
c specific heat of thermal sublimation, eV
c [tk]=kelvin
c -------------------------------------------
      common/psigsu/aln10,         !ln(10)
     ,              msub,          ! atom mass, amu
     ,              asub,bsub,     !coefs for sublimation equation
     ,              znum
      integer znum
      real*8 aln10
      real*8 msub,dels
      real*8 asub,bsub
      real*8 tktoev
      data tktoev /11604.5d0/
      real*8 psi_ef_li,psi_ef_be,psi_ef_b,psi_ef_c,psi_ef_fe,psi_ef_mo,
     ,       psi_ef_w

c---Li
      if (znum .eq. 3) then
       dels=psi_ef_li(tk)
c---Be
      elseif (znum .eq. 4) then
       dels=psi_ef_be(tk)
c---B
      elseif (znum .eq. 5) then
       dels=psi_ef_b(tk)
c---C
      elseif (znum .eq. 6) then
       dels=psi_ef_c(tk)
c---Fe
      elseif (znum .eq. 26) then
       dels=psi_ef_fe(tk)
c---Mo
      elseif (znum .eq. 42) then
       dels=psi_ef_mo(tk)
c---W
      elseif (znum .eq. 74) then
       dels=psi_ef_w(tk)
c---Dummy      
      else
       dels=0.d0
      endif

      psiesub=dels+2.d0*tk/tktoev
      return
      end


      integer function psicsub(iz)
      implicit none
c --------------------------------------------------------------
c -- prepare constants for saturated vapor pressure
c --------------------------------------------------------------
      integer iz

      common/psigsu/aln10,         !ln(10)
     ,              msub,          ! atom mass, amu
     ,              asub,bsub,     !coefs for sublimation equation
     ,              znum
      integer znum
      real*8 aln10
      real*8 msub
      real*8 asub,bsub
      real*8 psimassa
      external psimassa
      integer k

      k=0
      aln10=dlog(10.d0)
      
c---Li
      if(iz .eq. 3) then
       asub=8194.4     !7676.
       bsub=10.217     !9.78
       msub=psimassa(3)
       znum=3
       k=k+1
       
c---BE
      elseif(iz .eq. 4) then
       asub=16795.     !16720.
       bsub=11.319     !11.61
       msub=psimassa(4)
       znum=4
       k=k+1

c---B
      elseif(iz .eq. 5) then
       asub=29013.     !25576.
       bsub=12.345     !9.4817
       msub=psimassa(5)
       znum=5
       k=k+1

c---C
      elseif(iz .eq. 6) then
       asub=38286.     !40181.
       bsub=14.522     !14.8
       msub=psimassa(6)
       znum=6
       k=k+1

c---Fe
      elseif(iz .eq. 26) then
       asub=21415.     !18259.
       bsub=12.187     !10.947
       msub=psimassa(26)
       znum=26
       k=k+1

c---Mo
      elseif(iz .eq. 42) then
       asub=34087.     !31232.
       bsub=12.204     !10.951
       msub=psimassa(42)
       znum=42
       k=k+1

c---W
      elseif(iz .eq. 74) then
       asub=44522.     !44485.
       bsub=12.626     !12.74
       msub=psimassa(74)
       znum=74
       k=k+1

c---DUMMY
      else
       asub=0.
       bsub=0.
       msub=1.
       znum=0
       k=k+1
      endif

      psicsub=k
      return

      end

c-----------------------------------------------------------------------
      real*8 function psi_ef_c(tk)
      implicit none
c  -----------------------------------------------------------
c    enthalpy of formation C(ref->gas) JANAF (averaged over C1-C5)
c	  T [grad K]
c    data dH [eV]
c  -----------------------------------------------------------
      integer nt
      data nt /60/
      real st(60),sl(60)
      real*8 tk
      real*8 psilint

      data st
     ,    /  200.0 ,  298.15,  300.0 ,  400.0 ,  500.0 ,
     ,       600.0 ,  700.0 ,  800.0 ,  900.0 , 1000.0 ,
     ,      1100.0 , 1200.0 , 1300.0 , 1400.0 , 1500.0 ,
     ,      1600.0 , 1700.0 , 1800.0 , 1900.0 , 2000.0 ,
     ,      2100.0 , 2200.0 , 2300.0 , 2400.0 , 2500.0 ,
     ,      2600.0 , 2700.0 , 2800.0 , 2900.0 , 3000.0 ,
     ,      3100.0 , 3200.0 , 3300.0 , 3400.0 , 3500.0 ,
     ,      3600.0 , 3700.0 , 3800.0 , 3900.0 , 4000.0 ,
     ,      4100.0 , 4200.0 , 4300.0 , 4400.0 , 4500.0 ,
     ,      4600.0 , 4700.0 , 4800.0 , 4900.0 , 5000.0 ,
     ,      5100.0 , 5200.0 , 5300.0 , 5400.0 , 5500.0 ,
     ,      5600.0 , 5700.0 , 5800.0 , 5900.0 , 6000.0 /

      data sl
     ,    / 7.413  , 7.428  , 7.428  , 7.439  , 7.447  ,
     ,      7.452  , 7.455  , 7.453  , 7.440  , 7.396  ,
     ,      7.291  , 7.092  , 6.783  , 6.378  , 5.921  ,
     ,      5.461  , 5.042  , 4.681  , 4.378  , 4.135  ,
     ,      3.936  , 3.777  , 3.649  , 3.543  , 3.458  ,
     ,      3.388  , 3.327  , 3.278  , 3.236  , 3.199  ,
     ,      3.168  , 3.140  , 3.116  , 3.095  , 3.075  ,
     ,      3.057  , 3.041  , 3.025  , 3.012  , 2.998  ,
     ,      2.985  , 2.973  , 2.962  , 2.952  , 2.941  ,
     ,      2.930  , 2.920  , 2.910  , 2.910  , 2.890  ,
     ,      2.880  , 2.870  , 2.860  , 2.851  , 2.841  ,
     ,      2.832  , 2.823  , 2.814  , 2.804  , 2.795  /

      psi_ef_c = psilint(sl,st, nt,tk)
      return
      end


c-----------------------------------------------------------------------
      real*8 function psi_ef_b(tk)
      implicit none
c  -----------------------------------------------------------
c    enthalpy of formation B(glass,l->gas) JANAF (averaged over B1-B2)
c	  T [grad K]
c    data dH [eV]
c  -----------------------------------------------------------
      integer nt
      data nt /47/
      real st(47),sl(47)
      real*8 tk
      real*8 psilint

      data st
     ,    /  200.0  ,  250.0  ,  298.15 ,  300.0  ,  350.0  ,
     ,       400.0  ,  450.0  ,  500.0  ,  600.0  ,  700.0  ,
     ,       800.0  ,  900.0  , 1000.0  , 1100.0  , 1200.0  ,
     ,      1300.0  , 1400.0  , 1500.0  , 1600.0  , 1700.0  ,
     ,      1750.000, 1750.001, 1800.0  , 1900.0  , 2000.0  ,
     ,      2100.0  , 2200.0  , 2300.0  , 2400.0  , 2500.0  ,
     ,      2600.0  , 2700.0  , 2800.0  , 2900.0  , 3000.0  ,
     ,      3100.0  , 3200.0  , 3300.0  , 3400.0  , 3500.0  ,
     ,      3600.0  , 3700.0  , 3800.0  , 3900.0  , 4000.0  ,
     ,      4100.0  , 6000.0  /

      data sl
     ,    / 5.792  , 5.799  , 5.804  , 5.804  , 5.808  ,
     ,      5.812  , 5.814  , 5.815  , 5.816  , 5.815  ,
     ,      5.813  , 5.811  , 5.806  , 5.802  , 5.796  ,
     ,      5.790  , 5.784  , 5.777  , 5.770  , 5.762  ,
     ,      5.759  , 5.251  , 5.246  , 5.234  , 5.223  ,
     ,      5.212  , 5.200  , 5.189  , 5.177  , 5.166  ,
     ,      5.154  , 5.143  , 5.131  , 5.119  , 5.107  ,
     ,      5.095  , 5.083  , 5.071  , 5.059  , 5.046  ,
     ,      5.034  , 5.021  , 5.008  , 4.995  , 4.982  ,
     ,      4.969  , 4.717  /

      psi_ef_b = psilint(sl,st, nt,tk)
      return
      end

c-----------------------------------------------------------------------
      real*8 function psi_ef_li(tk)
      implicit none
c  -----------------------------------------------------------
c    enthalpy of formation Li(cr,l->gas) JANAF (averaged over Li1-Li2)
c	  T [grad K]
c    data dH [eV]
c  -----------------------------------------------------------
      integer nt
      data nt /23/
      real st(23),sl(23)
      real*8 tk
      real*8 psilint

      data st
     ,    /  100.0  ,  200.0  ,  250.0  ,  298.15 ,  300.0  ,
     ,       350.0  ,  400.0  ,  450.0  ,  453.690,  453.691,
     ,       500.0  ,  600.0  ,  700.0  ,  800.0  ,  900.0  ,
     ,      1000.0  , 1100.0  , 1200.0  , 1300.0  , 1400.0  ,
     ,      1500.0  , 1600.0  , 6000.0  /

      data sl
     ,    / 1.651  , 1.654  , 1.653  , 1.651  , 1.651  ,
     ,      1.649  , 1.646  , 1.642  , 1.641  , 1.615  ,
     ,      1.606  , 1.596  , 1.587  , 1.577  , 1.566  ,
     ,      1.553  , 1.539  , 1.524  , 1.508  , 1.491  ,
     ,      1.474  , 1.456  , 0.685  /

      psi_ef_li = psilint(sl,st, nt,tk)
      return
      end

c-----------------------------------------------------------------------
      real*8 function psi_ef_be(tk)
      implicit none
c  -----------------------------------------------------------
c    enthalpy of formation Be(gl,l->g) JANAF (averaged over Be1-Be2)
c	  T [grad K]
c    data dH [eV]
c  -----------------------------------------------------------
      integer nt
      data nt /31/
      real st(31),sl(31)
      real*8 tk
      real*8 psilint

      data st
     ,    /  100.0  ,  200.0  ,  298.15 ,  300.0  ,  400.0  ,
     ,       500.0  ,  600.0  ,  700.0  ,  800.0  ,  900.0  ,
     ,      1000.0  , 1100.0  , 1150.000, 1150.001, 1200.0  ,
     ,      1300.0  , 1400.0  , 1500.0  , 1600.0  , 1700.0  ,
     ,      1800.0  , 1900.0  , 2000.0  , 2100.0  , 2200.0  ,
     ,      2300.0  , 2400.0  , 2500.0  , 2600.0  , 2700.0  ,
     ,      6000.0  /

      data sl
     ,    / 3.335  , 3.351  , 3.358  , 3.358  , 3.361  ,
     ,      3.360  , 3.358  , 3.355  , 3.351  , 3.345  ,
     ,      3.339  , 3.332  , 3.328  , 3.167  , 3.164  ,
     ,      3.156  , 3.148  , 3.140  , 3.132  , 3.123  ,
     ,      3.115  , 3.106  , 3.097  , 3.088  , 3.078  ,
     ,      3.068  , 3.058  , 3.048  , 3.037  , 3.026  ,
     ,      2.675  /

      psi_ef_be = psilint(sl,st, nt,tk)
      return
      end

c-----------------------------------------------------------------------
      real*8 function psi_ef_fe(tk)
      implicit none
c  -----------------------------------------------------------
c    enthalpy of formation Fe(gl,l->gas) JANAF
c	  T [grad K]
c    data dH [eV]
c  -----------------------------------------------------------
      integer nt
      data nt /37/
      real st(37),sl(37)
      real*8 tk
      real*8 psilint

      data st
     ,    /  100.0  ,  200.0  ,  250.0  ,  298.15 ,  300.0  ,
     ,       350.0  ,  400.0  ,  450.0  ,  500.0  ,  600.0  ,
     ,       700.0  ,  800.0  ,  900.0  , 1000.0  , 1100.0  ,
     ,      1200.000, 1200.001, 1300.0  , 1400.0  , 1500.0  ,
     ,      1600.0  , 1700.0  , 1800.0  , 1900.0  , 2000.0  ,
     ,      2100.0  , 2200.0  , 2300.0  , 2400.0  , 2500.0  ,
     ,      2600.0  , 2700.0  , 2800.0  , 2900.0  , 3000.0  ,
     ,      3100.0  , 6000.0  /
      data sl
     ,    / 4.299  , 4.304  , 4.306  , 4.306  , 4.306  ,
     ,      4.306  , 4.305  , 4.304  , 4.302  , 4.295  ,
     ,      4.286  , 4.273  , 4.255  , 4.229  , 4.194  ,
     ,      4.164  , 4.067  , 4.042  , 4.017  , 3.993  ,
     ,      3.968  , 3.944  , 3.919  , 3.895  , 3.871  ,
     ,      3.848  , 3.825  , 3.802  , 3.779  , 3.757  ,
     ,      3.735  , 3.713  , 3.692  , 3.671  , 3.650  ,
     ,      3.630  , 3.040  /

      psi_ef_fe = psilint(sl,st, nt,tk)
      return
      end

c-----------------------------------------------------------------------
      real*8 function psi_ef_mo(tk)
      implicit none
c  -----------------------------------------------------------
c    enthalpy of formation Mo(gl,l->gas) JANAF
c	  T [grad K]
c    data dH [eV]
c  -----------------------------------------------------------
      integer nt
      data nt /56/
      real st(56),sl(56)
      real*8 tk
      real*8 psilint

      data st
     ,    /  100.0  ,  200.0  ,  250.0  ,  298.15 ,  300.0  ,
     ,       350.0  ,  400.0  ,  450.0  ,  500.0  ,  600.0  ,
     ,       700.0  ,  800.0  ,  900.0  , 1000.0  , 1100.0  ,
     ,      1200.0  , 1300.0  , 1400.0  , 1500.0  , 1600.0  ,
     ,      1700.0  , 1800.0  , 1900.0  , 2000.0  , 2100.0  ,
     ,      2150.000, 2150.001, 2200.0  , 2300.0  , 2400.0  ,
     ,      2500.0  , 2600.0  , 2700.0  , 2800.0  , 2900.0  ,
     ,      3000.0  , 3100.0  , 3200.0  , 3300.0  , 3400.0  ,
     ,      3500.0  , 3600.0  , 3700.0  , 3800.0  , 3900.0  ,
     ,      4000.0  , 4100.0  , 4200.0  , 4300.0  , 4400.0  ,
     ,      4500.0  , 4600.0  , 4700.0  , 4800.0  , 4900.0  ,
     ,      6000.0  /

      data sl
     ,    / 6.830  , 6.832  , 6.831  , 6.830  , 6.830  ,
     ,      6.828  , 6.826  , 6.824  , 6.821  , 6.815  ,
     ,      6.809  , 6.803  , 6.795  , 6.788  , 6.780  ,
     ,      6.771  , 6.762  , 6.752  , 6.741  , 6.729  ,
     ,      6.717  , 6.703  , 6.689  , 6.674  , 6.658  ,
     ,      6.650  , 6.218  , 6.210  , 6.194  , 6.178  ,
     ,      6.163  , 6.149  , 6.135  , 6.121  , 6.109  ,
     ,      6.097  , 6.086  , 6.076  , 6.067  , 6.060  ,
     ,      6.053  , 6.048  , 6.044  , 6.042  , 6.040  ,
     ,      6.041  , 6.043  , 6.046  , 6.051  , 6.057  ,
     ,      6.066  , 6.075  , 6.086  , 6.099  , 6.113  ,
     ,      6.271  /

      psi_ef_mo = psilint(sl,st, nt,tk)
      return
      end


c-----------------------------------------------------------------------
      real*8 function psi_ef_w(tk)
      implicit none
c  -----------------------------------------------------------
c    enthalpy of formation W(gl,l->gas) JANAF
c	  T [grad K]
c    data dH [eV]
c  -----------------------------------------------------------
      integer nt
      data nt /66/
      real st(66),sl(66)
      real*8 tk
      real*8 psilint

      data st
     ,    /  100.0  ,  200.0  ,  250.0  ,  298.15 ,  300.0  ,
     ,       350.0  ,  400.0  ,  450.0  ,  500.0  ,  600.0  ,
     ,       700.0  ,  800.0  ,  900.0  , 1000.0  , 1100.0  ,
     ,      1200.0  , 1300.0  , 1400.0  , 1500.0  , 1600.0  ,
     ,      1700.0  , 1800.0  , 1900.0  , 2000.0  , 2100.0  ,
     ,      2200.0  , 2300.0  , 2400.0  , 2450.000, 2450.001,
     ,      2500.0  , 2600.0  , 2700.0  , 2800.0  , 2900.0  ,
     ,      3000.0  , 3100.0  , 3200.0  , 3300.0  , 3400.0  ,
     ,      3500.0  , 3600.0  , 3700.0  , 3800.0  , 3900.0  ,
     ,      4000.0  , 4100.0  , 4200.0  , 4300.0  , 4400.0  ,
     ,      4500.0  , 4600.0  , 4700.0  , 4800.0  , 4900.0  ,
     ,      5000.0  , 5100.0  , 5200.0  , 5300.0  , 5400.0  ,
     ,      5500.0  , 5600.0  , 5700.0  , 5800.0  , 5900.0  ,
     ,      6000.0  /

      data sl
     ,    / 8.822  , 8.823  , 8.822  , 8.820  , 8.820  ,
     ,      8.819  , 8.818  , 8.817  , 8.817  , 8.820  ,
     ,      8.827  , 8.837  , 8.849  , 8.863  , 8.877  ,
     ,      8.891  , 8.903  , 8.914  , 8.923  , 8.930  ,
     ,      8.935  , 8.939  , 8.942  , 8.943  , 8.942  ,
     ,      8.941  , 8.939  , 8.935  , 8.934  , 8.447  ,
     ,      8.445  , 8.439  , 8.433  , 8.428  , 8.422  ,
     ,      8.416  , 8.411  , 8.405  , 8.400  , 8.396  ,
     ,      8.391  , 8.387  , 8.383  , 8.380  , 8.377  ,
     ,      8.375  , 8.373  , 8.372  , 8.372  , 8.372  ,
     ,      8.372  , 8.374  , 8.376  , 8.379  , 8.383  ,
     ,      8.388  , 8.394  , 8.400  , 8.408  , 8.417  ,
     ,      8.427  , 8.439  , 8.451  , 8.465  , 8.481  ,
     ,      8.497  /

      psi_ef_w = psilint(sl,st, nt,tk)
      return
      end

