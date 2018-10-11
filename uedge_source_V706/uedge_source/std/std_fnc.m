c-----------------------------------------------------------------------
c Linear interpolation over 1D non-uniform mesh
c-----------------------------------------------------------------------
      real*8 function std_lint(ss,ee,n, te)

      implicit none

      integer n
      integer i,k
      real*8 ss(*),ee(*)
      real*8 te
      real*8 s1,s2,x1,x2,x,d,f

      if(n .lt. 2) then
       std_lint=dble(ss(1))
       return
      endif

      x=te
      k=0
      do 10 i=1,n
       d=ee(i)
       if(x .lt. d) then
        k=i
        goto 1
       endif
  10  continue

   1  if(k .eq. 0) then
       std_lint=dble(ss(n))
       return
      endif

      if(k .eq. 1) then
       std_lint=dble(ss(1))
       return
      endif

        i = k-1
       s1 = dble(ss(i))
       x1 = dble(ee(i))
       s2 = dble(ss(k))
       x2 = dble(ee(k))
        f = (x-x1)/(x2-x1)

      if(f .lt. 0.d0) then
       f=0.d0
      endif
      if(f .gt. 1.d0) then
       f=1.d0
      endif

      std_lint = s1+(s2-s1)*f
      return

      end



c-----------------------------------------------------------------------
c   EXPONENTIAL INTEGRAL E1(X)
c-----------------------------------------------------------------------
      real*8 function std_E1(x)

      implicit none

      real*8 x
      real*8 xhi, t, y
      data xhi /174.673d0/

      if( x .le.  0.d0 ) then
       std_E1 = 0.d0
       return 
      endif

      if( x .gt. 4.d0 ) then     
       if(x .ge. xhi) then
         std_E1 = 0.d0
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
       std_E1 = y
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
      std_E1 = y
      return
 
      end



      real*8 function std_egamm(x)

      implicit none

      real*8 x 
c
c------------------------------
c  ENCOMPLETE GAMMA-FUNCTION
c------------------------------
c
      real*8  std_E1
      real*8  f, p, q, t, z, fq 
      real*8 cei,weps,sam,peps

      data cei  /0.57721566490153d0/ 
      data weps /83.5d0/ 
      data sam  /1.d75/ 
      data peps /1.d-8/

      if(x .ge. 1.d0) then
       if(x .gt. weps) then
        std_egamm = sam
        return 
       else
        std_egamm = std_E1(x)
        return
       endif 
      endif

      if(x .eq. 0.d0) then
       std_egamm = -cei
       return 
      endif

      if(x .lt. -weps) then
       std_egamm = -sam
       return 
      endif
c
c     ANALITICAL EXTENSION AT X<0
c
      p = cei 
      if(x .lt. 0.d0) then
       z = -x  
      else 
       z = x
      endif 
      p  = p + dlog(z) 
      f  = 1.d0 
      t  = -z

1     q  = t + p 
      fq = p - q 
      if ( fq .lt.  0.d0 ) fq = -fq 
      if ( fq .lt. peps ) goto 2 
      p  = q 
      t  = (-z) * f * t 
      f  = f + 1.d0 
      fq = f * f 
      t  = t/fq 
      goto 1

2     std_egamm = -p 
      return

      end


c-----------------------------------------------------------------------
c   Error function: erf(x)
c-----------------------------------------------------------------------
      real*8 function std_erf(xx)
      implicit none

      real*8  xx,t,f,ff

      real*8  a1 /0.0705230784d0/
      real*8  a2 /0.0422820123d0/
      real*8  a3 /0.0092705272d0/
      real*8  a4 /0.0001510143d0/
      real*8  a5 /0.0002765672d0/
      real*8  a6 /0.0000430638d0/


      if(xx .lt. 0.d0) then
      t = -xx
      else
      t = xx
      endif

      ff = 1.d0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6)))))
      ff = ff*ff
      ff = ff*ff
      ff = ff*ff
      ff = ff*ff
      f  = 1.d0/ff
      ff = 1.d0 - f

      if(xx .lt. 0.d0) then
       f = -ff
      else
       f =  ff
      endif

      std_erf = f

      return
      end
