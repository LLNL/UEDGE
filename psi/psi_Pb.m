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
c
c      ===============================
c		Table Mask
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


      double precision function psi_Pb(iz)
      implicit none
      integer iz
c
c---------------------------------------------------------
c  coefficient in R-D formular for barrier penetrability
c---------------------------------------------------------
c
      real am(92)
      data am /
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.25, 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 ,
     , 0.5 , 0.5 /

      psi_Pb=dble(am(iz))
      return

      end

