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
      double precision function psipspty(e0)
      implicit none
      double precision e0     

      double precision psiypspt
      double precision amu
      data amu /0.75d0/

      psipspty=psiypspt(e0,amu)

      return
      end



      subroutine psipsptc(izpt,ampt,izsf,amsf)
      implicit none
      integer izpt,izsf
      double precision ampt,amsf

      common/psipspr/q,w,fc,ebd,eth,gamm,pii0
      double precision q,w,fc,ebd,eth,gamm,pii0

      double precision psiwsub

      integer i
      double precision s,ss,c
      double precision zpt,zsf

      zpt = dble(izpt)
      zsf = dble(izsf)
      ebd = psiwsub(izsf)

         s = ampt + amsf
        ss = s*s
      gamm = 4.d0*ampt*amsf/ss
        s  = 4.d0*ampt+amsf
        s  = s*s
        s  = s/(4.d0*amsf*ampt)
       eth = ebd*s

      i=0
      if(izpt .eq. 1 ) i=i+1
      if(ampt .lt. 1.5d0) i=i+1
      if(i .eq. 2) then
      c = 2000.d0
      else
      c = 400.d0
      endif

      s  = dsqrt(zpt)
      ss = dsqrt(s)
      s  = zpt / ss
      ss = (ampt-0.8d0) / amsf
      ss = ss*dsqrt(ss)
      q  = zsf + (-1.8d0)
      q  = q*q
      q  = q* c * ss*s / (ebd*eth)
      w  = 50.d0 * s*zsf/eth - 1.d0

         s = -1.d0
      pii0 = dacos( s )

      fc = dsqrt( amsf*eth/ampt )
      fc = dsqrt( fc/(20.d0*zpt))

      end


      double precision function psiypspt(e,amu)
      implicit none
      double precision e,amu
c
c---------------------------------------------
c      dsput-code for phisical sputtering 
c---------------------------------------------
c
      common/psipspr/q,w,fc,ebd,eth,gamm,pii0
      double precision q,w,fc,ebd,eth,gamm,pii0

      integer noamu
      double precision s,ss,sss,t,x,g,seps1

      data noamu /1/
      data seps1 / 1.000000001d0/

      x = e/eth
      if(x < seps1 ) then
      psiypspt=0.d0
      return
      endif

      s  = x + w
      t  = s * s
      s  = (x - 1.d0) / t
      s  = s*q

      if(noamu > 0) then

      g=1.d0

      else

 	   ss = amu
      if(ss < 0.d0) ss = -ss
      if(ss > 0.999) then
       g=1.d0
      else if(ss < 1.d-6) then
	 g=0.d0
      else if( x < 4.00000001d0) then
	 g=1.d0
      else
	 sss= ss*ss
	 t  = 1.d0-sss
	 if(t < 0.d0) then
         sss= -t 
         else 
         sss=t
       endif
	 t  = dsqrt(sss)/ss
	 sss= 2.d0/pii0
	 t  = datan(t)*sss
	 if(t < 0.d0 ) t = -t
	 if(t > 1.d0 ) then
         t = 0.d0
         else
         t = 1.d0 - t
       endif
	 t  = t*dsqrt(x - 4.d0)
	 sss= fc * dsqrt(t)
	 g  = dexp( (-sss) * dlog(ss) )
       endif

      endif !noamu

      psiypspt=s*g
      return 
      end



      double precision function psiepspt(ec)
      implicit none
      double precision ec
c
c--------------------------------------------
c   averaged energy of sputtered particles
c	    for Thompson spectrum
c--------------------------------------------
c
      common/psipspr/q,w,fc,ebd,eth,gamm,pii0
      double precision q,w,fc,ebd,eth,gamm,pii0

      double precision sd,d,s,redf

      data redf /0.3d0/

      sd = ebd/eth*ec*redf
      d  = ebd/(ebd+sd)
      sd = 1.d0-d
      s  = dlog(d*d)+2.d0*sd
      sd = sd*sd
      s  = s/sd + 1.d0
      s  = (-s) * ebd

      psiepspt=s
      return

      end
