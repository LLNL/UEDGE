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
      subroutine psiseep(iz)
      implicit none
      
      integer iz
      double precision psiseedm,psiseeem,psiseegc
      external psiseedm,psiseeem,psiseegc
      double precision seedm,seeem,seegc
      common /seemat/ seedm,seeem,seegc
      
      seedm=psiseedm(iz)
      seeem=psiseeem(iz)
      seegc=psiseegc(iz)
      
      end

      double precision function psieseecw(e)
      implicit none
      
      double precision e
      double precision erat
      double precision seedm,seeem,seegc
      common /seemat/ seedm,seeem,seegc

      erat=e/seeem
      psieseecw=7.3984d0*seedm*erat*dexp(-2.d0*dsqrt(erat))
      return
      end


      double precision function psiphsee (amu)
      implicit none

      double precision amu
      double precision seedm,seeem,seegc
      common /seemat/ seedm,seeem,seegc

      psiphsee=amu**(-seegc)
      return
      end


      double precision function psisee(e,amu)
      implicit none

      double precision e,amu
      double precision psieseecw,psiphsee
      external psieseecw,psiphsee

      psisee=psieseecw(e)*psiphsee(amu)
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
      double precision function psiseedm(iz)
      implicit none
      integer iz
c-----------------------------------------
c   max SEE yield of chemical element iz
c-----------------------------------------
      real dm(92)
      data dm /
     , 0.0    , 0.0     , 0.59   , 0.5     , 1.2    ,
     , 1.0    , 0.0     , 0.0    , 0.0     , 0.0    ,
     , 0.82   , 0.95    , 1.0    , 1.1     , 0.0    ,
     , 0.0    , 0.0     , 0.0    , 0.0     , 0.0    ,
     , 0.0    , 0.9     , 0.0    , 0.0     , 0.0    ,
     , 1.3    , 0.0     , 1.35   , 0.0     , 0.0    ,
     , 0.0    , 0.0     , 0.0    , 0.0     , 0.0    ,
     , 0.0    , 0.9     , 0.0    , 0.0     , 1.1    ,
     , 1.2    , 1.25    , 0.0    , 0.0     , 0.0    ,
     , 1.3    , 0.0     , 0.0    , 0.0     , 1.35   ,
     , 1.3    , 0.0     , 0.0    , 0.0     , 0.0    ,
     , 0.0    , 0.0     , 0.0    , 0.0     , 0.0    ,
     , 0.0    , 0.0     , 0.0    , 0.0     , 0.0    ,
     , 0.0    , 0.0     , 0.0    , 0.0     , 0.0    ,
     , 0.0    , 0.0     , 1.3    , 1.4     , 0.0    ,
     , 0.0    , 0.0     , 1.8    , 0.0     , 0.0    ,
     , 1.7    , 1.1     , 0.0    , 0.0     , 0.0    ,
     , 0.0    , 0.0     , 0.0    , 0.0     , 1.1    ,
     , 0.0    , 0.0     /

      psiseedm=dble(dm(iz))
      return
      end


      double precision function psiseeem(iz)
      implicit none
      integer iz
c----------------------------------------------------
c   energy at max SEE yoeld of chemical element iz
c----------------------------------------------------
      real em(92)
      data em /
     , 1.0    , 1.0     , 85.    , 200.    , 150.   ,
     , 300.   , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 300.   , 300.    , 300.   , 250.    , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 280.    , 1.0    , 1.0     , 1.0    ,
     , 400.   , 1.0     , 550.   , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 350.    , 1.0    , 1.0     , 350.   ,
     , 375.   , 375.    , 1.0    , 1.0     , 1.0    ,
     , 250.   , 1.0     , 1.0    , 1.0     , 500.   ,
     , 600.   , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 600.   , 650.    , 1.0    ,
     , 1.0    , 1.0     , 700.   , 1.0     , 1.0    ,
     , 650.   , 500.    , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 800.   ,
     , 1.0    , 1.0     /

      psiseeem=dble(em(iz))
      return
      end


      double precision function psiseegc(iz)
      implicit none
      integer iz
c----------------------------------------------------
c   energy at max SEE yoeld of chemical element iz
c----------------------------------------------------
      real gc(92)
      data gc /
     , 1.5    , 1.5     , 1.5    , 1.3     , 1.3    ,
     , 1.3    , 1.3     , 1.3    , 1.3     , 1.3    ,
     , 1.2    , 1.2     , 1.2    , 1.2     , 1.2    ,
     , 1.1    , 1.1     , 1.1    , 1.1     , 1.1    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     , 1.0    , 1.0     , 1.0    ,
     , 1.0    , 1.0     /

      psiseegc=dble(gc(iz))
      return
      end
