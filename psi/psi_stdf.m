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
      double precision function psilint(ss,ee,n, te)

      implicit none

      real ss(*),ee(*)
      integer n
      double precision te
c-----------------------------------------------
c Linear interpolation over 1D non-uniform mesh
c-----------------------------------------------
      integer i,k
      double precision s1,s2,x1,x2,x,d,f

      if(n < 2) then
       psilint=dble(ss(1))
       return
      endif

      x=te
      k=0
      do 10 i=1,n
       d=ee(i)
       if(x < d) then
        k=i
        goto 1
       endif
  10  continue

   1  if(k .eq. 0) then
      psilint=dble(ss(n))
      return
      endif

      if(k .eq. 1) then
       psilint=dble(ss(1))
       return
      endif

        i = k-1
       s1 = dble(ss(i))
       x1 = dble(ee(i))
       s2 = dble(ss(k))
       x2 = dble(ee(k))
        f = (x-x1)/(x2-x1)
      if(f < 0.d0) then
       f=0.d0
      endif
      if(f > 1.d0) then
       f=1.d0
      endif

      psilint = s1+(s2-s1)*f
      return

      end


      double precision function psilintd(ss,ee,n, te)
      implicit none

      double precision ss(*),ee(*)
      integer n
      double precision te
c-----------------------------------------------
c Linear interpolation over 1D non-uniform mesh
c-----------------------------------------------
      integer i,k
      double precision s1,s2,x1,x2,x,d,f

      if(n < 2) then
       psilintd=ss(1)
       return
      endif

      x=te
      k=0
      do 10 i=1,n
       d=ee(i)
       if(x < d) then
        k=i
        goto 1
       endif
  10  continue

   1  if(k .eq. 0) then
       psilintd=ss(n)
       return
      endif

      if(k .eq. 1) then
       psilintd=ss(1)
       return
      endif

        i = k-1
       s1 = ss(i)
       x1 = ee(i)
       s2 = ss(k)
       x2 = ee(k)
        f = (x-x1)/(x2-x1)
      if(f < 0.d0) then
       f=0.d0
      endif
      if(f > 1.d0) then
       f=1.d0
      endif

      psilintd = s1+(s2-s1)*f
      return

      end


      double precision function psiE1(x)
      implicit none
      double precision x
c--------------------------------
c   EXPONENTIAL INTEGRAL E1(X)
c--------------------------------
      double precision xhi, t, y
      data xhi /174.673d0/

      if( x .le.  0.d0 ) then
       psiE1 = 0.d0
       return 
      endif

      if( x > 4.d0 ) then     
       if(x .ge. xhi) then
        psiE1 = 0.d0
        return 
       endif
       goto 1 
      endif

      t = x *  0.5d0 +  (-1.d0) 
       y =    7.420476912680064d-1 + t * ( 8.646647167633864d-1 +
     & t * ( -2.969970751450805d-1 + t * ( 1.077745279390218d-1 +
     & t * ( -3.571913487537507d-2 + t * ( 1.053060346814121d-2 +
     & t * ( -2.760601413314358d-3 + t * ( 6.476865075230550d-4 +
     & t * ( -1.370898716058566d-4 + t * ( 2.638302391770205d-5 +
     & t * ( -4.649805703397782d-6 + t * ( 7.553168050042681d-7 +
     & t * ( -1.137209741906139d-7 + t * ( 1.592437332203632d-8 +
     & t * ( -2.090265558347686d-9 + t * ( 2.724676194855373d-10+
     & t * ( -3.158414049758465d-11)))))))))))))))) 
       y = y - dlog(x) 
       psiE1 = y
       return
 
1      t =  14.5d0 / ( x +  3.25d0 ) +  (-1.0d0) 
       y =    9.237553078077841d-1 + t * ( -8.506991549845731d-2 +
     & t * ( -1.063748751165807d-2 + t * ( -2.142130549996267d-3 +
     & t * ( -4.058921304226172d-4 + t * ( -9.164504841773118d-5 +
     & t * ( -1.950763781198026d-5 + t * ( -4.720903675164850d-6 +
     & t * ( -1.055648511722045d-6 + t * ( -2.703104833900864d-7 +
     & t * ( -6.118194012663569d-8 + t * ( -1.668705862724518d-8 +
     & t * ( -3.663650270142246d-9 + t * ( -1.119940934924999d-9 +
     & t * ( -2.211936213417643d-10+ t * ( -6.804206766157201d-11+
     & t * ( -1.048399621792147d-11+
     & t * ( -1.181683914463340d-11))))))))))))))))) 
       y = y * dexp( -x ) / x 
      psiE1 = y
      return
 
      end



      double precision function psiegamm(x)
      implicit none
      double precision x 
c------------------------------
c  ENCOMPLETE GAMMA-FUNCTION
c------------------------------
      double precision  psiE1
      double precision  f, p, q, t, z, fq 
      double precision cei,weps,sam,peps

      data cei  /0.57721566490153d0/ 
      data weps /83.5d0/ 
      data sam  /1.d75/ 
      data peps /1.d-8/

      if(x .ge. 1.d0) then
       if(x > weps) then
        psiegamm = sam
        return 
       else
        psiegamm = psiE1(x)
        return
       endif 
      endif

      if(x .eq. 0.d0) then
       psiegamm = -cei
       return 
      endif

      if(x < -weps) then
       psiegamm = -sam
       return 
      endif
c
c     ANALITICAL EXTENSION AT X<0
c
      p = cei 
      if(x <  0.d0) then
      z = -x  
      else 
      z = x
      endif 
      p  = p + dlog(z) 
      f  = 1.d0 
      t  = -z


1     q  = t + p 
      fq = p - q 
      if ( fq <  0.d0 ) fq = -fq 
      if ( fq < peps ) goto 2 
      p  = q 
      t  = (-z) * f * t 
      f  = f + 1.d0 
      fq = f * f 
      t  = t/fq 
      goto 1
2     psiegamm = -p 
      return

      end


      double precision function psierf(xx)
c
c erf(x)
c
      implicit none
      double precision  xx

      double precision t,f,ff

      double precision  a1 /0.0705230784/
      double precision  a2 /0.0422820123/
      double precision  a3 /0.0092705272/
      double precision  a4 /0.0001510143/
      double precision  a5 /0.0002765672/
      double precision  a6 /0.0000430638/

      if(xx < 0.d0) then
       t = -xx
      else
       t = xx
      endif

      ff = 1.d0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6)))))
      ff=ff*ff
      ff=ff*ff
      ff=ff*ff
      ff=ff*ff
      f =1.d0/ff
      ff = 1.d0 - f

      if(xx < 0.d0) then
       f = -ff
      else
       f=ff
      endif

      psierf=f
      return

      end




      double precision function psilngamma(xx)
      implicit none
      double precision xx
c
c---------------------------------------
c natural logarithm of gamma function
c	    ln[gamma(x)]
c---------------------------------------
c
      double precision xsmall /1.0d-15/
      double precision xbig   /6.8d6/
      double precision xvbig  /4.29d73/
      double precision gbig   /7.231d75/
      double precision aln2pi /9.189385332046727d-1/


      double precision f,g,x,t
      integer i
      integer m


      x = xx


      if(x .le. 0.d0) then
        f = 0.d0
      else
      if(x .le. xsmall) then
        f = -dlog(x)
      else
       if(x .gt. 15.d0) then

        if(x .gt. xbig) then
           if(x .gt. xvbig) then
	    f = gbig
           else
	    g = dlog(x)
	    f = (x - 0.5)*g - x + aln2pi
           endif 
        else
           f=x*x
           t=450.d0/f - 1.d0
           f=((((2.002019273379824d-14)*t-6.451144077929628d-12)*t+
     &           3.899788998764847d-9 )*t-6.165020494506090d-6 )*t+
     &           8.332716440657866d-2


           g=dlog(x)
           t=f/x
           f=(x - 0.5)*g - x + aln2pi + t
        endif


      else
        m = int(x)
        t = x - dble(m)
        m = m - 1
        g = 1.d0


        if(m < 0) then
         g = g/x
        endif

        if(m > 0) then
         do 5 i=1,m
          g = g * ( x - dble(i) )
 5       continue
       endif
      endif
      endif


        t = 2.0*t - 1.d0


        f = ((((((((((((((((((((   (-1.243191705600000d-10)*t+
     +    3.622882508800000d-10)*t - 4.030909644800000d-10)*t+
     +    1.265236705280000d-9 )*t - 5.419466096640000d-9 )*t+
     +    1.613133578240000d-8 )*t - 4.620920340480000d-8 )*t+
     +    1.387603440435200d-7 )*t - 4.179652784537600d-7 )*t+
     +    1.253148247777280d-6 )*t - 3.754930502328320d-6 )*t+
     +    1.125234962812416d-5 )*t - 3.363759801664768d-5 )*t+
     +    1.009281733953869d-4 )*t - 2.968901194293069d-4 )*t+
     +    9.157859942174304d-4 )*t - 2.422595384546340d-3 )*t+
     +    9.040334940477911d-3 )*t - 1.341185057058971d-2 )*t+
     +    1.037033634220705d-1 )*t + 1.616919872444243d-2 )*t+
     +    8.862269254527580d-1

c       f = ((((((((((((((((       (-8.965837291520000d-9 )*t+
c    +    2.612707393536000d-8 )*t - 3.802866827264000d-8 )*t+
c    +    1.173294768947000d-7 )*t - 4.275076254106000d-7 )*t+
c    +    1.276176602829000d-6 )*t - 3.748495971011000d-6 )*t+
c    +    1.123829871408000d-5 )*t - 3.364018663166000d-5 )*t+
c    +    1.009331480887000d-4 )*t - 2.968895120407000d-4 )*t+
c    +    9.157850115110000d-4 )*t - 2.422595461409000d-3 )*t+
c    +    9.040335037321000d-3 )*t - 1.341185056618000d-2 )*t+
c    +    1.037033634184000d-1 )*t - 1.616919872437000d-2 )*t+
c    +    8.862269254528000d-1

        t = f*g
        f = dlog(t)
      endif


      psilngamma=f
      return

      end




      double precision function psifunsig(xx)
      implicit none
      double precision xx
c
c---------------------
c Sigmund function
c---------------------
c
      double precision psidigamma


      double precision seps /1.d-8/
      double precision x,f,f1,fx

      x=xx

      if(x < seps) then
       f = 1.d0


      elseif(x > 1.d0) then
       f  = 0.d0


      else
       f  = 1.d0
       f1 = psidigamma(f)


       f  = 1.d0 - x
       fx = psidigamma(f)


       f = x / (f1 - fx)
      endif

      psifunSIG=f
      return
      end





      double precision function psidigamma(xx)
      implicit none
      double precision xx
c
c---------------------------------------
c	  digamma function
c  derivative of ln[gamma(x)] function
c	but only for  0 < x <= 1
c---------------------------------------
c


      double precision xsmall /1.0d-15/
      double precision xone   /0.9999999d0/
      double precision egam   /0.5772156649015329d0/


      double precision f,g,x,t


      if(xx .gt. xone) then
       x = xone
      else 
       x = xx
      endif


      if(x .le. 0.d0) then
       f = 0.d0


      elseif(x .le. xsmall) then
       f = (-1.d0)/x


      else
       t = 2.0*x - 1.d0


       f = ((((((((((((((((((((  (-1.243191705600000d-10)*t+
     +  3.622882508800000d-10)*t - 4.030909644800000d-10)*t+
     +  1.265236705280000d-9 )*t - 5.419466096640000d-9 )*t+
     +  1.613133578240000d-8 )*t - 4.620920340480000d-8 )*t+
     +  1.387603440435200d-7 )*t - 4.179652784537600d-7 )*t+
     +  1.253148247777280d-6 )*t - 3.754930502328320d-6 )*t+
     +  1.125234962812416d-5 )*t - 3.363759801664768d-5 )*t+
     +  1.009281733953869d-4 )*t - 2.968901194293069d-4 )*t+
     +  9.157859942174304d-4 )*t - 2.422595384546340d-3 )*t+
     +  9.040334940477911d-3 )*t - 1.341185057058971d-2 )*t+
     +  1.037033634220705d-1 )*t + 1.616919872444243d-2 )*t+
     +  8.862269254527580d-1
c      f = ((((((((((((((((      (-8.965837291520000d-9 )*t+
c    +  2.612707393536000d-8 )*t - 3.802866827264000d-8 )*t+
c    +  1.173294768947000d-7 )*t - 4.275076254106000d-7 )*t+
c    +  1.276176602829000d-6 )*t - 3.748495971011000d-6 )*t+
c    +  1.123829871408000d-5 )*t - 3.364018663166000d-5 )*t+
c    +  1.009331480887000d-4 )*t - 2.968895120407000d-4 )*t+
c    +  9.157850115110000d-4 )*t - 2.422595461409000d-3 )*t+
c    +  9.040335037321000d-3 )*t - 1.341185056618000d-2 )*t+
c    +  1.037033634184000d-1 )*t - 1.616919872437000d-2 )*t+
c    +  8.862269254528000d-1

       g = (((((((((((((((((((       (-21.0*1.243191705600000d-10)*t+
     + 20.0*3.622882508800000d-10)*t - 19.0*4.030909644800000d-10)*t+
     + 18.0*1.265236705280000d-9 )*t - 17.0*5.419466096640000d-9 )*t+
     + 16.0*1.613133578240000d-8 )*t - 15.0*4.620920340480000d-8 )*t+
     + 14.0*1.387603440435200d-7 )*t - 13.0*4.179652784537600d-7 )*t+
     + 12.0*1.253148247777280d-6 )*t - 11.0*3.754930502328320d-6 )*t+
     + 10.0*1.125234962812416d-5 )*t -  9.0*3.363759801664768d-5 )*t+
     +  8.0*1.009281733953869d-4 )*t -  7.0*2.968901194293069d-4 )*t+
     +  6.0*9.157859942174304d-4 )*t -  5.0*2.422595384546340d-3 )*t+
     +  4.0*9.040334940477911d-3 )*t -  3.0*1.341185057058971d-2 )*t+
     +  2.0*1.037033634220705d-1 )*t +  1.0*1.616919872444243d-2
c      g = (((((((((((((((         (-17.0*8.965837291520000d-9 )*t+
c    +16.0*2.612707393536000d-8 )*t -15.0*3.802866827264000d-8 )*t+
c    +14.0*1.173294768947000d-7 )*t -13.0*4.275076254106000d-7 )*t+
c    +12.0*1.276176602829000d-6 )*t -11.0*3.748495971011000d-6 )*t+
c    +10.0*1.123829871408000d-5 )*t - 9.0*3.364018663166000d-5 )*t+
c    + 8.0*1.009331480887000d-4 )*t - 7.0*2.968895120407000d-4 )*t+
c    + 6.0*9.157850115110000d-4 )*t - 5.0*2.422595461409000d-3 )*t+
c    + 4.0*9.040335037321000d-3 )*t - 3.0*1.341185056618000d-2 )*t+
c    + 2.0*1.037033634184000d-1 )*t - 1.0*1.616919872437000d-2


       t = (2.d0*x*g)/f
       f = -1.d0/x + t
      endif

      psidigamma=f
      return

      end




      double precision function psigrl(ext,a,b,eps)
      implicit none
      double precision ext
      external ext
      double precision a,b,eps
c----------------------------------------
c    integration of function ext(x)
c  from x=a to x=b using gauss method
c----------------------------------------
      double precision w(12)
      data w/
     , 1.01228536290376d-1, 2.22381034453374d-1,
     , 3.13706645877887d-1, 3.62683783378362d-1,
     , 2.71524594117540d-2, 6.22535239386480d-2,
     , 9.51585116824930d-2, 1.24628971255534d-1,
     , 1.49595988816577d-1, 1.69156519395003d-1,
     , 1.82603415044924d-1, 1.89450610455069d-1/

      double precision x(12)
      data x/
     , 9.60289856497536d-1, 7.96666477413627d-1,
     , 5.25532409916329d-1, 1.83434642495650d-1,
     , 9.89400934991650d-1, 9.44575023073233d-1,
     , 8.65631202387832d-1, 7.55404408355003d-1,
     , 6.17876244402644d-1, 4.58016777657227d-1,
     , 2.81603550779259d-1, 9.50125098376370d-2/

      double precision, parameter :: PMI_DBL_SMALL=1.d-111 

      common/psigrlcom/k
      integer k

      integer i
      double precision constan,delta
      double precision u,c1,c2,aa,bb,s8,s16,sum,y,ay,epss

      k=0
      constan=PMI_DBL_SMALL

      delta = a-b
      if(delta .lt. 0.d0) then
      delta = -delta
      endif
      delta = delta*constan
      if(delta .lt. constan) then
      delta = constan
      endif

      sum = 0.d0
      aa  = a

 1    continue
      if(k .eq. 0) then
      y = b-aa
       if(y < 0.d0 ) then
       ay = -y
       else
       ay = y
       endif
      if(ay .lt. delta) goto 2
      endif

      bb = y+aa
      c1 = (aa+bb)*0.5d0
      c2 = c1-aa

      s8 = 0.d0
      do 10 i=1,4
       u = x(i)*c2
      s8 = s8 + w(i)*(ext(c1+u)+ext(c1-u))
  10  continue
      s8 = s8*c2

      s16 = 0.d0
      do 20 i=5,12
       u  = x(i)*c2
      s16 = s16 + w(i)*(ext(c1+u)+ext(c1-u))
  20  continue
      s16=s16*c2

      epss = s16
      if(epss .lt. 0.d0 ) then
      epss = -epss
      endif
      epss = eps * (1.d0 + epss)

      u = s16-s8
      if(u .lt. 0.d0 ) then
       u = -u
      endif

      if(u .le. epss) then
      sum = sum+s16
       aa = bb
	k = 0
      else
       y = y*0.5
       if(y .lt. 0.d0 ) then
       ay = -y
       else
       ay = y
       endif
       if(ay .lt. delta) then
       print *, "psigrl:too high accuracy required"
       goto 2
       endif
	k = 1
      endif

      goto 1

2     psigrl=sum
      return

      end


      double precision function psigrlf(ext,fun,a,b,eps)
      implicit none
      double precision ext,fun
      external ext,fun
      double precision a,b,eps
c----------------------------------------
c    integration of function ext(x,fun)
c  from x=a to x=b using gauss method
c----------------------------------------
      double precision w(12)
      data w/
     , 1.01228536290376d-1, 2.22381034453374d-1,
     , 3.13706645877887d-1, 3.62683783378362d-1,
     , 2.71524594117540d-2, 6.22535239386480d-2,
     , 9.51585116824930d-2, 1.24628971255534d-1,
     , 1.49595988816577d-1, 1.69156519395003d-1,
     , 1.82603415044924d-1, 1.89450610455069d-1/

      double precision x(12)
      data x/
     , 9.60289856497536d-1, 7.96666477413627d-1,
     , 5.25532409916329d-1, 1.83434642495650d-1,
     , 9.89400934991650d-1, 9.44575023073233d-1,
     , 8.65631202387832d-1, 7.55404408355003d-1,
     , 6.17876244402644d-1, 4.58016777657227d-1,
     , 2.81603550779259d-1, 9.50125098376370d-2/

      double precision, parameter :: PMI_DBL_SMALL=1.d-111 

      common/psigrlcom/k
      integer k

      integer i
      double precision constan,delta
      double precision u,c1,c2,aa,bb,s8,s16,sum,y,ay,epss

      k=0
      constan=PMI_DBL_SMALL


      delta = a-b
      if(delta .lt. 0.d0) then
      delta = -delta
      endif
      delta = delta*constan
      if(delta .lt. constan) then
      delta = constan
      endif

      sum = 0.d0
      aa  = a


 1    continue
      if(k .eq. 0) then
      y = b-aa
       if(y .lt. 0.d0 ) then
       ay = -y
       else
       ay = y
       endif
      if(ay .lt. delta) goto 2
      endif


      bb = y+aa
      c1 = (aa+bb)*0.5d0
      c2 = c1-aa


      s8 = 0.d0
      do 10 i=1,4
       u = x(i)*c2
      s8 = s8 + w(i)*(ext(c1+u,fun)+ext(c1-u,fun))
  10  continue
      s8 = s8*c2


      s16 = 0.d0
      do 20 i=5,12
       u  = x(i)*c2
      s16 = s16 + w(i)*(ext(c1+u,fun)+ext(c1-u,fun))
  20  continue
      s16=s16*c2


      epss = s16
      if(epss .lt. 0.d0 ) then
      epss = -epss
      endif
      epss = eps * (1.d0 + epss)


 	 u = s16-s8
      if(u .lt. 0.d0 ) then
       u = -u
      endif

      if(u .le. epss) then
      sum = sum+s16
       aa = bb
	k = 0
      else
       y = y*0.5
       if(y .lt. 0.d0 ) then
       ay = -y
       else
       ay = y
       endif
       if(ay .lt. delta) then
       print *, "psigrlf:too high accuracy required"
       goto 2
       endif
	k = 1
      endif


      goto 1

2     psigrlf=sum
      return

      end


      double precision function psi_root(ext,a,b,eps,jnet)
      implicit none
      double precision ext
      external ext
      double precision a,b,eps
      integer jnet
c--------------------------------
c   let's find a root: ext(x)=0
c--------------------------------
      double precision epsr
      data epsr/1.d-100/

      integer kmx
      data kmx/100000/

      integer k,kc
      integer ia,ib,ix
      double precision a1,b1,x,f,d,eps1

      k=0
      kc=0

      a1 = a
      b1 = b
      d  = b1-a1
      if( d .lt. 0.d0 ) d = -d

      if(d .lt. epsr) then
       jnet = 0
       psi_root=a1
       return
      endif

      eps1 = eps * d
      if(eps1 .lt. epsr) eps1 = epsr

      f = ext(a1)
      if(f .lt. 0.d0) then
       ia = 1
      else
       ia = 2
      endif

      f = ext(b1)
      if(f .lt. 0.d0) then
       ib = 1
      else
       ib = 2
      endif

      if(ia .eq. ib) then
       jnet = 2
       psi_root=a1
       return
      endif

      if(b1 .lt. a1) then
       a1 = b
       b1 = a
       ix = ia
       ia = ib
       ib = ix
      endif

 99   continue
      x = 0.5d0 * (a1+b1)

      k=k+1
      if(k .gt. kmx) then
       kc = 1
       goto 999
      endif

      f = ext(x)
      if(f .lt. 0.d0) then
       ix = 1
      else
       ix = 2
      endif

      if(ia .ne. ix) then
       if(ib .ne. ix) then
        kc = 3
        goto 999
       endif
       b1 = x
      elseif(ib .ne. ix) then
       a1 = x
      else
       kc = 3
       goto 999
      endif

      d = d * 0.5d0

      if(d .gt. eps1) goto 99

 999  jnet = kc
      psi_root=x
      return

      end
