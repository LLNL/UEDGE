c============================================
      real*8 function radiation(u)
      implicit none
      
      Use (const)
      
      real*8 u,lambda,q,tmp
      
      real*8 rad,temp
      integer iz
      common /dust1/ rad,temp
      common /dust2/ iz

      complex*16 e
      
      real*8 planck
      real*8 qabs_w,qabs_cu,qabs_mo,qabs_li,qabs_be,qabs_ss,qabs_c
      complex*16 eps_w,eps_cu,eps_mo,eps_li,eps_be,eps_ss,eps_c
      external planck
      external qabs_w,qabs_cu,qabs_mo,qabs_li,qabs_be,qabs_ss,qabs_c
      external eps_w,eps_cu,eps_mo,eps_li,eps_be,eps_ss,eps_c
      
      integer iappro

c ----------------------------------------------------------------------
c -----------  I/O SPECIFICATIONS FOR SUBROUTINE  MIEV0  ---------------
c ----------------------------------------------------------------------
      integer   maxang,momdim
      parameter ( maxang=3, momdim=200)
      logical   anyang,perfct,prnt(2)
      integer   ipolzn,numang,nmom
      real*8    gqsc,mimcut,pmom(0:momdim,4),qext,qsca,spike,
     *          xmu(maxang),xx
      complex*16 crefin,sforw,sback,s1(maxang),s2(maxang),
     *           tforw(2),tback(2)
c ----------------------------------------------------------------------

c --------------- LOCAL VARIABLES --------------------------------------

c     .. Local Scalars ..
      logical   nopmom
      integer   i,j,ncas,npquan
      real*8    n,k
      real*8    degpol,fnorm,i1,i2,inten,t1,t2,testin
c     ..
c     .. Local Arrays ..
      real*8    angle(maxang)
c     ..
c     .. External Subroutines ..
      external  miev0
c     ..
c ----------------------------------------------------------------------
      nopmom = .false.
c                            ** Set MIEV0 input variables
      mimcut = 0.d0
      numang = maxang

      do 10, i=1,numang
       angle(i)=(i-1)*180.d0/(numang-1)
       xmu  (i)=dcos(pi/180.d0*angle(i))
10    continue

      perfct = .false.
      anyang = .true.
      ipolzn = +1234
      nmom   = 200

      if (nopmom) nmom=0

      prnt(1) = .false.
      prnt(2) = .false.

      
      lambda=hplanck*c/(u*kboltz*temp)
      xx=2.d0*pi*rad/lambda

c-------------------
c    e=eps_material
c------------------
      
c --- Li      
      if(iz .eq. 3)then
       e=eps_li(lambda,rad,temp)

c --- Be
      elseif(iz .eq. 4) then
       e=eps_be(lambda,rad,temp)

c --- B      
c      elseif(iz .eq. 5) then
c       e=eps_b(lambda,rad,temp)

c --- C      
      elseif(iz .eq. 6) then
       e=eps_c(lambda,rad,temp)

c --- Fe      
      elseif(iz .eq. 26) then
       e=eps_ss(lambda,rad,temp)

c --- Cu      
      elseif(iz .eq. 26) then
       e=eps_cu(lambda,rad,temp)

c --- Mo      
      elseif(iz .eq. 42) then
       e=eps_mo(lambda,rad,temp)

c --- W      
      elseif(iz .eq. 74) then
       e=eps_w(lambda,rad,temp)

c --- DUMMY
      else
       e=1.d0
      endif

c--------------------

      crefin=cdsqrt(e)
      
c  --- Mie absorbtion ---
      call miev0( xx,crefin,perfct,mimcut,anyang,numang,xmu,nmom,
     *            ipolzn,momdim,prnt,qext,qsca,gqsc,pmom,sforw,sback,
     *            s1,s2,tforw,tback,spike)
      q=qext-qsca
      
c-------------------
c    e=qabs_"material" approximation
c------------------
      iappro=0
      if (iappro .eq. 1) then

c --- Li      
       if (iz .eq. 3) then
        q=qabs_li(lambda,rad,temp)

c --- Be
       elseif(iz .eq. 4) then
        q=qabs_be(lambda,rad,temp)

c --- B      
c       elseif(iz .eq. 5) then
c        q=qabs_b(lambda,rad,temp)

c --- C      
       elseif(iz .eq. 6) then
        q=qabs_c(lambda,rad,temp)

c --- Fe      
       elseif(iz .eq. 26) then
        q=qabs_ss(lambda,rad,temp)

c --- Cu      
       elseif(iz .eq. 26) then
        q=qabs_cu(lambda,rad,temp)

c --- Mo      
       elseif(iz .eq. 42) then
        q=qabs_mo(lambda,rad,temp)

c --- W      
       elseif(iz .eq. 74) then
        q=qabs_w(lambda,rad,temp)

c --- DUMMY
       else
        q=1.d0
       endif
      endif

      tmp=planck(u,temp)*q
      radiation=tmp
      
      end
c========================================================================
      REAL*8 FUNCTION PLANCK(U,TEMP)
      IMPLICIT NONE
      
      Use (const)
      
      REAL*8 U,TEMP
      REAL*8 TMP
      
      TMP=KBOLTZ*TEMP
      TMP=TMP*TMP
      TMP=TMP*TMP
      TMP=TMP*U*U*U
      TMP=2.D0*PI*TMP/(C*C*HPLANCK*HPLANCK*HPLANCK*(DEXP(U)-1.D0))
      
      PLANCK=TMP
      RETURN
      END
    
c==========================================================================
