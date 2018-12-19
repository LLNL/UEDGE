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
C
C RES from :
C   J. Roth, E.Vietzke, A.A.Haaz, Atomic and Plasma-Material Interaction
C   Data for Fusion, Nucl. Fusion Supplement, Vol-1 (1991) 63-78.
C-----------------------------------------------------------------------

      subroutine psiresc(izpt,ampt,izsf,amsf)
      implicit none
      
      integer izpt,izsf
      real*8 ampt,amsf
c
c  constants for RES on carbon
c
      real*8 c1yres,wres,ethres,etfres,pires,ethy,fcres,tdres
      common/psiresy/c1yres,wres,ethres,etfres,
     ,               pires,ethy,fcres,tdres

      real*8 pwr2by3,zsum,am2byam1
      real*8 z1,z2,am1,am2,es,s
      real*8 psiwsub,psiresef,psirescf
      external psiwsub,psiresef,psirescf

      z1=dble(izpt)
      z2=dble(izsf)
      am1=ampt
      am2=amsf
      es=psiwsub(izsf)
      
      wres=psiresef(izsf)*11604.5d0
      c1yres=psirescf(izsf)*(am1**1.18d0)

      pwr2by3 = 2.0d0 / 3.0d0
      zsum = z1**pwr2by3 + z2**pwr2by3
      am2byam1 = am2 / am1

      ethres = (7.0d0/am2byam1**0.54d0 + 0.15d0*am2byam1**1.12d0) * es
      etfres = 30.74d0 * (am1+am2)/am2 * z1*z2*dsqrt(zsum)


          s = -1.d0
      pires = dacos(s)

        s  = 4.d0*am1+am2
        s  = s*s
        s  = s/(4.d0*am2*am1)
      ethy = es*s

      fcres = dsqrt(am2*ethy/am1)
      fcres = dsqrt(fcres/(20.d0*z1))

      end


      real*8 function psiwres(idummy)
      implicit none
      integer idummy

      common/psiresy/c1yres,wres,ethres,etfres,
     ,               pires,ethy,fcres,tdres
      real*8 c1yres,wres,ethres,etfres,
     ,       pires,ethy,fcres,tdres

      psiwres=ethres
      return

      end



      real*8 function psiyresgc(eo,amu,td)
      implicit none
      real*8 eo,amu,td

      common/psiresy/c1yres,wres,ethres,etfres,
     ,               pires,ethy,fcres,tdres
      real*8 c1yres,wres,ethres,etfres,
     ,       pires,ethy,fcres,tdres
      real*8 psiresugl
      real*8 yldpar,yldres
      real*8 stoppwr,stoppwr1,stoppwr2,ethbyeo,eobyetf

      yldres=0.d0

      if(eo .gt. ethres) then
       yldpar = c1yres * dexp(-wres/(td+0.2d0))
       ethbyeo = ethres / eo
       eobyetf = eo / etfres
       stoppwr1 = (0.5d0*dlog(1.0d0+1.2288d0*eobyetf))
       stoppwr2 = (eobyetf + 0.1728d0*dsqrt(eobyetf) +
     @             0.008d0*eobyetf**0.1504d0)
       stoppwr = stoppwr1 / stoppwr2
       yldres = yldpar * stoppwr * (1.0d0-ethbyeo**0.6666666666d0)
     @        * ( 1.0d0 - ethbyeo ) ** 2
       yldres=yldres*psiresugl(eo,amu)
       if (yldres .le. 0.0d0) yldres = 0.d0
      endif
      
      psiyresgc=yldres
      return

      end


      real*8 function psiyresfc(gamma)
      implicit none
      real*8 gamma

      real*8 f

      f=(gamma*1.0d4)**(-0.1)

      psiyresfc=f
      return

      end


      real*8 function psiresugl(e,amu)
      implicit none
      real*8 e,amu

      common/psiresy/c1yres,wres,ethres,etfres,
     ,               pires,ethy,fcres,tdres
      real*8 c1yres,wres,ethres,etfres,
     ,       pires,ethy,fcres,tdres
      real*8 gmax
      data gmax/2.5d0/

      real*8 g,ss,sss,t,x

         x=e/ethres
         ss = amu
      if(ss .lt. 0.d0) ss = -ss
      if(ss .gt. 0.999) then
       g=1.d0
      else if(ss .lt. 1.d-6) then
         g=0.d0
      else if(x .lt. 4.00000001d0) then
         g=1.d0
      else
         sss= ss*ss
         t  = 1.d0-sss
         if(t .lt. 0.d0) then
         sss= -t 
         else 
         sss=t
       endif
         t  = dsqrt(sss)/ss
         sss= 2.d0/pires
         t  = datan(t)*sss
         if(t .lt. 0.d0 ) t = -t
         if(t .gt. 1.d0 ) then
         t = 0.d0
         else
         t = 1.d0 - t
       endif
         t  = t*dsqrt(x - 4.d0)
         sss= fcres * dsqrt(t)
         g  = dexp( (-sss) * dlog(ss) )
       endif

      if(g .gt. gmax) g=gmax

      psiresugl=g
      return

      end


      subroutine psirestd(td)
      implicit none
      real*8 td

      common/psiresy/c1yres,wres,ethres,etfres,
     ,               pires,ethy,fcres,tdres
      real*8 c1yres,wres,ethres,etfres,
     ,                 pires,ethy,fcres,tdres

      tdres=td
      end



      real*8 function psiresyy(eo,amu)
      implicit none
      real*8 eo,amu

      common/psiresy/c1yres,wres,ethres,etfres,
     ,               pires,ethy,fcres,tdres
      real*8 c1yres,wres,ethres,etfres,
     ,                 pires,ethy,fcres,tdres
      real*8 psiyresgc
      real*8 y

      y=psiyresgc(eo,amu,tdres)

      psiresyy=y
      return

      end


      real*8 function psiyresc(eo,amu,td,gamma)
c
c  calculates RES yield on carbon
c  eo   -> incident ion energy (in eV).
c  td   -> temperature of the target plate (K).
c  gamma -> incident ion flux (10^20/cm^2/sec).
c
      implicit none
      real*8 eo,amu,td,gamma

      real*8 psiyresgc,psiyresfc
      real*8 y,yf

      y=psiyresgc(eo,amu,td)
      yf=psiyresfc(gamma)

      psiyresc=y*yf
      return

      end

c
c      ===============================
c                Table Mask
c      ===============================
c       "h " ,"he" ,"li" ,"be" ,"b ",
c       "c " ,"n " ,"o " ,"f " ,"ne",
c       "na" ,"mg" ,"al" ,"si" ,"p ",
c       "s " ,"cl" ,"ar" ,"k " ,"ca",
c       "sc" ,"ti" ,"v " ,"cr" ,"mn",
c       "fe" ,"co" ,"ni" ,"cu" ,"zn",
c       "ga" ,"ge" ,"as" ,"se" ,"br",
c       "kr" ,"rb" ,"sr" ,"y " ,"zr",
c       "nb" ,"mo" ,"tc" ,"ru" ,"rh",
c       "pd" ,"ag" ,"cd" ,"in" ,"sn",
c       "sb" ,"te" ,"i " ,"xe" ,"cs",
c       "ba" ,"la" ,"ce" ,"pr" ,"nd",
c       "pm" ,"sm" ,"eu" ,"gd" ,"tb",
c       "dy" ,"ho" ,"er" ,"tm" ,"yb",
c       "lu" ,"hf" ,"ta" ,"w " ,"re",
c       "os" ,"ir" ,"pt" ,"au" ,"hg",
c       "tl" ,"pb" ,"bi" ,"po" ,"at",
c       "rn" ,"fr" ,"ra" ,"ac" ,"th",
c       "pa" ,"u "
c
      double precision function psiresef(iz)
      implicit none
      integer iz
c-----------------------------------------
c   effective binding energy for RES, eV
c   (needs update)
c-----------------------------------------
      real ef(92)
      data ef /
     , 0.01, 0.01, 1.1 , 1.8 , 1.0 ,
     , 0.78, 1.0 , 1.0 , 0.84, 0.02,
     , 1.0 , 1.0 , 1.0 , 0.84, 1.0 ,
     , 1.0 , 1.0 , 0.08, 0.93, 1.0 ,
     , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
     , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
     , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
     , 0.12, 0.85, 1.0 , 1.0 , 1.0 ,
     , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
     , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
     , 1.0 , 1.0 , 1.0 , 0.16, 0.80,
     , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
     , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
     , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
     , 1.0 , 1.0 , 1.0 , 6.0 , 1.0 ,
     , 1.0 , 1.0 , 1.0 , 1.0 , 0.67,
     , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 ,
     , 0.20, 1.0 , 1.0 , 1.0 , 1.0 ,
     , 1.0 , 1.0 /

      psiresef=dble(ef(iz))
      return
      end

      double precision function psirescf(iz)
      implicit none
      integer iz
c-----------------------------------------
c   effective binding energy for RES, eV
c   (needs update)
c-----------------------------------------
      real cf(92)
      data cf /
     , 1.0  , 1.0  , 1.e7 , 2.e7 , 1.e2 ,
     , 54.0 , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.e2 , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.e2 , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.e2 , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  ,
     , 1.0  , 1.0  /

      psirescf=dble(cf(iz))
      return
      end
