c==================================================================
      real*8 function qabs_li(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 x,rho1,rho2,rhoabs
      real*8 tmp,tmp1,tmp2
      real*8 n1,n2,e1,e2
      real*8 c0
      complex*16 n,e
      complex*16 alphae,alpham,alpha,beta,s,a,b,c1,c2
      complex*16 sina,cosa,cota,cotb,expia
      
      complex*16 eps_li
      external eps_li
      
      e=eps_li(lambda,rad,temp)
      n=cdsqrt(e)
      
      n1=dreal(n)
      n2=dimag(n)
      e1=dreal(e)
      e2=dimag(e)
      
      x=2.d0*pi*rad/lambda
      rho1=x*n1
      rho2=x*n2
      rhoabs=rho1*rho1+rho2*rho2
c --- separate small/large x model ---
      if (x.le.1.d0) then
       beta=n*x
       cotb=cdcos(beta)/cdsin(beta)
       s=(beta*beta-1.d0+beta*cotb)/(beta*(1.d0-beta*cotb))
       alphae=3.d0/(4.d0*pi)*(1.d0-x*s/(2.d0*n))/(1.d0+x*s/n)
       tmp1=16.d0*pi/3.d0*x*dimag(alphae)
c       tmp1=12.d0*x*e2/((e1+2.d0)**2+e2*e2)
       tmp1=tmp1+6.d0*x*
     &     ((rho1*dsinh(2.d0*rho2)+rho2*dsin(2.d0*rho1))/
     &     (rhoabs*(dcosh(2.d0*rho2)-dcos(2.d0*rho1)))
     &     -2.d0*rho1*rho2/(rhoabs*rhoabs))
       tmp=tmp1
      else
       if (n2.ge.1.d0) then
        tmp2=16.d0/3.d0*n1/(n1*n1+n2*n2)
c        tmp2=4.d0*n1/(n1*n1+n2*n2+1.d0+2.d0*n1)
        tmp=tmp2
       else
c        tmp=8.d0/3.d0*x*n2
        c0=4.d0*n2*x
        if (c0.le.1.d-4) then
         tmp=c0*(2.d0/3.d0-0.25d0*c0)
        else
         tmp=1.d0-2.d0*(1.d0-(1.d0+c0)*dexp(-c0))/(c0*c0)
        endif
       endif
      endif
c --- dipol model for any x ---
c      alpha=x
c      beta=x*n
c      sina=cdsin(alpha)
c      cosa=cdcos(alpha)
c      cota=cosa/sina
c      cotb=one/cdtan(beta)
c      expia=cdexp(ione*alpha)
c      s=(beta*beta-one+beta*cotb)/(beta*(one-beta*cotb))
c      a=(alpha*alpha-one+alpha*cota)/(alpha*(one-alpha*cota))
c      b=-(one-alpha*alpha+ione*alpha)/(alpha*(one+ione*alpha))
c      c1=ione*expia*(sina-alpha*cosa)/(one+alpha)
c      c2=-ione*expia*((alpha*alpha-one)*sina+alpha*cosa)/
c     &   (one-alpha*alpha+ione*alpha)
c      alphae=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c2*(one-s/(a*n))/(one-s/(b*n))
c      alpham=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c1*(one-a/(s*n))/(one-b/(s*n))
c      tmp=16.d0*pi/3.d0*x*dimag(alphae+alpham)
      qabs_li=tmp
      return
      end
c      
c
c
      complex*16 function eps_li(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      real*8 omega,omega2,omegap,omegap2,gamma,gamma2
      real*8 e1,e2
      complex*16 e
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 fomegap_li,feta_li
      external fomegap_li,feta_li
      
      omega=2.d0*pi*c/lambda
      omegap=fomegap_li(temp)
      omegap2=omegap*omegap
      gamma=omegap2*feta_li(temp)*eps0
      omega2=omega*omega
      gamma2=gamma*gamma
      e1=1.d0-omegap2/(omega2+gamma2)
      e2=omegap2*gamma/(omega*(omega2+gamma2))
      e=dcmplx(e1,e2)
      
      eps_li=e
      
      return
      end
c     
c
c
c
      real*8 function feta_li(temp)
      implicit none
      
      real*8 temp
      
      integer n,i
      real*8 tmp
      parameter (n=4)
      real*8 coeff1(n),coeff2(n),coeff3(n),coeff4(n),coeff5(n)
      data coeff1 /-1.173d0, 3.193d0, 7.549d0,-17.43d0 /
      data coeff2 /0.0139d0, 2.904d0,-8.494d0, 38.64d0 /
      data coeff3 /0.1575d0, 2.314d0,-1.962d0, 1.127d0 /
      data coeff4 / 1.395d0, 0.622d0,-0.228d0, 0.430d0 /
      data coeff5 / 1.620d0, 0.634d0, 0.258d0, 0.314d0 /
      
      real*8 copt
      data copt / 1.00d0 /
      
      if     (temp.le.40.0d0   ) then
       tmp=dlog10(0.074d0)
      elseif (temp.le.81.06d0  ) then
       tmp=dlog10(temp/40.0d0)
       tmp=((coeff1(4)*tmp+coeff1(3))*tmp+coeff1(2))*tmp+coeff1(1)
      elseif (temp.le.92.295d0 ) then
       tmp=dlog10(temp/81.06d0)
       tmp=((coeff2(4)*tmp+coeff2(3))*tmp+coeff2(2))*tmp+coeff2(1)
      elseif (temp.le.453.7d0  ) then
       tmp=dlog10(temp/92.295d0)
       tmp=((coeff3(4)*tmp+coeff3(3))*tmp+coeff3(2))*tmp+coeff3(1)
      elseif (temp.le.1080.5d0 ) then
       tmp=dlog10(temp/453.7d0)
       tmp=((coeff4(4)*tmp+coeff4(3))*tmp+coeff4(2))*tmp+coeff4(1)
      elseif (temp.le.2200.0d0 ) then
       tmp=dlog10(temp/1080.5d0)
       tmp=((coeff5(4)*tmp+coeff5(3))*tmp+coeff5(2))*tmp+coeff5(1)
      else
       tmp=dlog10(85.59d0)
      endif
      tmp=10.d0**tmp
      
      feta_li=copt*tmp*1.d-08
      return
      end
c


c      real*8 function feta_li(temp)
c      integer i, m
c      parameter (m=24)
c      real*8 temp, tmp
c      real*8 t(m), r(m)
c      data t/
c     *        300.0d0 ,  350.0d0 ,  400.0d0 ,  450.0d0 , 453.70d0 ,
c     *       453.71d0 ,  500.0d0 ,  600.0d0 ,  700.0d0 ,  800.0d0 ,
c     *        900.0d0 , 1000.0d0 , 1100.0d0 , 1200.0d0 , 1300.0d0 ,
c     *       1400.0d0 , 1500.0d0 , 1600.0d0 , 1700.0d0 , 1800.0d0 ,
c     *       1900.0d0 , 2000.0d0 , 2100.0d0 , 2200.0d0 /
c      data r/
c     *       9.540d0  , 11.44d0  , 13.39d0  , 15.43d0  , 15.58d0  ,
c     *       24.80d0  , 26.33d0  , 29.34d0  , 32.10d0  , 34.71d0  ,
c     *       37.22d0  , 39.69d0  , 42.13d0  , 44.61d0  , 47.41d0  ,
c     *       49.97d0  , 53.00d0  , 56.34d0  , 60.03d0  , 64.12d0  ,
c     *       68.67d0  , 73.73d0  , 79.44d0  , 85.59d0  /
c      
c      real*8 copt
c      data copt / 1.00d0 /
c      
c      if (temp.lt.t(1)) then
c       tmp=r(1)+(r(2)-r(1))*(temp-t(1))/(t(2)-t(1))
c      elseif (temp.gt.t(m)) then 
c       tmp=r(m)+(r(m)-r(m-1))*(temp-t(m))/(t(m)-t(m-1))
c      else
c       i=2
c10     if (temp.gt.t(i)) then
c        i=i+1
c        goto 10
c       endif
c       tmp=r(i-1)+(r(i)-r(i-1))*(temp-t(i-1))/(t(i)-t(i-1))
c      endif
c      
c      feta_li=copt*tmp*1.d-6   ! [Ohm*m]
c      return
c      end
c
c

c
c
      real*8 function fomegap_li(temp)
      implicit none
      
      real*8 temp
      real*8 tmp
      
      tmp=1.0d+16
      fomegap_li=tmp
      return
      end
      
c======================================================================
      real*8 function qabs_be(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 x,rho1,rho2,rhoabs
      real*8 tmp,tmp1,tmp2
      real*8 n1,n2,e1,e2
      real*8 c0
      complex*16 n,e
      complex*16 alphae,alpham,alpha,beta,s,a,b,c1,c2
      complex*16 sina,cosa,cota,cotb,expia
      
      complex*16 eps_be
      external eps_be
      
      e=eps_be(lambda,rad,temp)
      n=cdsqrt(e)
      
      n1=dreal(n)
      n2=dimag(n)
      e1=dreal(e)
      e2=dimag(e)
      
      x=2.d0*pi*rad/lambda
      rho1=x*n1
      rho2=x*n2
      rhoabs=rho1*rho1+rho2*rho2
c --- separate small/large x model ---
      if (x.le.1.d0) then
       beta=n*x
       cotb=cdcos(beta)/cdsin(beta)
       s=(beta*beta-1.d0+beta*cotb)/(beta*(1.d0-beta*cotb))
       alphae=3.d0/(4.d0*pi)*(1.d0-x*s/(2.d0*n))/(1.d0+x*s/n)
       tmp1=16.d0*pi/3.d0*x*dimag(alphae)
c       tmp1=12.d0*x*e2/((e1+2.d0)**2+e2*e2)
       tmp1=tmp1+6.d0*x*
     &     ((rho1*dsinh(2.d0*rho2)+rho2*dsin(2.d0*rho1))/
     &     (rhoabs*(dcosh(2.d0*rho2)-dcos(2.d0*rho1)))
     &     -2.d0*rho1*rho2/(rhoabs*rhoabs))
       tmp=tmp1
      else
       if (n2.ge.1.d0) then
        tmp2=16.d0/3.d0*n1/(n1*n1+n2*n2)
c        tmp2=4.d0*n1/(n1*n1+n2*n2+1.d0+2.d0*n1)
        tmp=tmp2
       else
c        tmp=8.d0/3.d0*x*n2
        c0=4.d0*n2*x
        if (c0.le.1.d-4) then
         tmp=c0*(2.d0/3.d0-0.25d0*c0)
        else
         tmp=1.d0-2.d0*(1.d0-(1.d0+c0)*dexp(-c0))/(c0*c0)
        endif
       endif
      endif
c --- dipol model for any x ---
c      alpha=x
c      beta=x*n
c      sina=cdsin(alpha)
c      cosa=cdcos(alpha)
c      cota=cosa/sina
c      cotb=one/cdtan(beta)
c      expia=cdexp(ione*alpha)
c      s=(beta*beta-one+beta*cotb)/(beta*(one-beta*cotb))
c      a=(alpha*alpha-one+alpha*cota)/(alpha*(one-alpha*cota))
c      b=-(one-alpha*alpha+ione*alpha)/(alpha*(one+ione*alpha))
c      c1=ione*expia*(sina-alpha*cosa)/(one+alpha)
c      c2=-ione*expia*((alpha*alpha-one)*sina+alpha*cosa)/
c     &   (one-alpha*alpha+ione*alpha)
c      alphae=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c2*(one-s/(a*n))/(one-s/(b*n))
c      alpham=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c1*(one-a/(s*n))/(one-b/(s*n))
c      tmp=16.d0*pi/3.d0*x*dimag(alphae+alpham)
      qabs_be=tmp
      return
      end
c      
c
c
      complex*16 function eps_be(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      real*8 omega,omega2,omegap,omegap2,gamma,gamma2
      real*8 e1,e2
      complex*16 e,ed
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 fomegap_be,feta_be
      complex*16 epsd_be
      external fomegap_be,feta_be,epsd_be
      
      omega=2.d0*pi*c/lambda
      omegap=fomegap_be(temp)
      omegap2=omegap*omegap
      gamma=omegap2*feta_be(temp)*eps0
      omega2=omega*omega
      gamma2=gamma*gamma
      e1=1.d0-omegap2/(omega2+gamma2)
      e2=omegap2*gamma/(omega*(omega2+gamma2))
      e=dcmplx(e1,e2)
      ed=epsd_be(omega)
c      e=e+ed
      
      eps_be=e
      
      return
      end
c     
c
c
      real*8 function feta_be(temp)
      integer i, m
      parameter (m=14)
      real*8 temp, tmp
      real*8 t(m), r(m)
c      data t/
c     *        300.0d0 ,  350.0d0 ,  400.0d0 ,  450.0d0 ,  500.0d0 ,
c     *        550.0d0 ,  600.0d0 ,  650.0d0 ,  700.0d0 ,  750.0d0 ,
c     *        800.0d0 ,  850.0d0 ,  900.0d0 ,  950.0d0 , 1000.0d0 ,
c     *       1100.0d0 , 1200.0d0 , 1300.0d0 , 1400.0d0 , 1500.0d0 /
c      data r/
c     *       3.730d0  , 5.190d0  , 6.730d0  , 8.300d0  , 9.910d0  ,
c     *       11.50d0  , 13.20d0  , 14.80d0  , 16.50d0  , 18.30d0  ,
c     *       20.00d0  , 21.80d0  , 23.70d0  , 25.60d0  , 27.50d0  ,
c     *       31.50d0  , 35.70d0  , 40.10d0  , 44.80d0  , 49.90d0  /
      data t/
     *         293.0d0 ,  525.0d0 ,  718.0d0 ,  896.0d0 , 1060.0d0 ,
     *        1217.0d0 , 1364.0d0 , 1505.0d0 , 1564.0d0 , 1564.1d0 ,
     *        1607.0d0 , 1762.0d0 , 1908.0d0 , 2062.0d0 /

      data r/
     *           4.6d0 ,   11.2d0 ,   18.5d0 ,   26.5d0 ,   37.4d0 ,
     *          50.4d0 ,   66.3d0 ,   88.9d0 ,  110.4d0 ,  161.9d0 ,
     *         164.2d0 ,  165.4d0 ,  165.5d0 ,  165.7d0 /
      
      real*8 copt
c      data copt / 2.57d0 /
      data copt / 1.97d0 /
            if (temp.lt.t(1)) then
       tmp=r(1)+(r(2)-r(1))*(temp-t(1))/(t(2)-t(1))
      elseif (temp.gt.t(m)) then 
       tmp=r(m)+(r(m)-r(m-1))*(temp-t(m))/(t(m)-t(m-1))
      else
       i=2
10     if (temp.gt.t(i)) then
        i=i+1
        goto 10
       endif
       tmp=r(i-1)+(r(i)-r(i-1))*(temp-t(i-1))/(t(i)-t(i-1))
      endif
      
      feta_be=copt*tmp*1.d-8   ! [Ohm*m]
      return
      end
c
c
c
      real*8 function fomegap_be(temp)
      implicit none
      
      real*8 temp
      real*8 tmp
      
      tmp=8.15d+15
      fomegap_be=tmp
      return
      end
c
c
c
      complex*16 function epsd_be(omega)
      implicit none
      Use (const)
      
      real*8 omega
      complex*16 tmp
      integer i,n
      parameter (n=4)
      real*8 f(n),gama(n),omegai(n),omegap,omegap2
      data f      /    0.031,     0.140,     0.530,     0.130/
      data gama   /2.528d+15, 5.158d+15, 6.767d+15, 2.738d+15/
      data omegai /1.519d+14, 1.568d+15, 4.836d+15, 6.995d+15/
      data omegap /2.812d+16/
      
      omegap2=omegap*omegap
      tmp=dcmplx(0.d0,0.d0)
      do 10, i=1,n
       tmp=tmp+f(i)*omegap2/((omegai(i)*omegai(i)-omega*omega)
     &                       +ione*omega*gama(i))
10    continue
      tmp=dconjg(tmp)
      
      epsd_be=tmp
      return
      end

c==================================================================
      real*8 function qabs_c(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 x,rho1,rho2,rhoabs
      real*8 tmp,tmp1,tmp2
      real*8 n1,n2,e1,e2
      real*8 c0
      complex*16 n,e
      complex*16 alphae,alpham,alpha,beta,s,a,b,c1,c2
      complex*16 sina,cosa,cota,cotb,expia
      
      complex*16 eps_c
      external   eps_c
      
      e=eps_c(lambda,rad,temp)
      n=cdsqrt(e)
      
      n1=dreal(n)
      n2=dimag(n)
      e1=dreal(e)
      e2=dimag(e)
      
      x=2.d0*pi*rad/lambda
      rho1=x*n1
      rho2=x*n2
      rhoabs=rho1*rho1+rho2*rho2
c --- separate small/large x model ---
      if (x.le.1.d0) then
       beta=n*x
       cotb=cdcos(beta)/cdsin(beta)
       s=(beta*beta-1.d0+beta*cotb)/(beta*(1.d0-beta*cotb))
       alphae=3.d0/(4.d0*pi)*(1.d0-x*s/(2.d0*n))/(1.d0+x*s/n)
       tmp1=16.d0*pi/3.d0*x*dimag(alphae)
c       tmp1=12.d0*x*e2/((e1+2.d0)**2+e2*e2)
       tmp1=tmp1+6.d0*x*
     &     ((rho1*dsinh(2.d0*rho2)+rho2*dsin(2.d0*rho1))/
     &     (rhoabs*(dcosh(2.d0*rho2)-dcos(2.d0*rho1)))
     &     -2.d0*rho1*rho2/(rhoabs*rhoabs))
       tmp=tmp1
      else
       if (n2.ge.1.d0) then
        tmp2=16.d0/3.d0*n1/(n1*n1+n2*n2)
c        tmp2=4.d0*n1/(n1*n1+n2*n2+1.d0+2.d0*n1)
        tmp=tmp2
       else
c        tmp=8.d0/3.d0*x*n2
        c0=4.d0*n2*x
        if (c0.le.1.d-4) then
         tmp=c0*(2.d0/3.d0-0.25d0*c0)
        else
         tmp=1.d0-2.d0*(1.d0-(1.d0+c0)*dexp(-c0))/(c0*c0)
        endif
       endif
      endif
c --- dipol model for any x ---
c      alpha=x
c      beta=x*n
c      sina=cdsin(alpha)
c      cosa=cdcos(alpha)
c      cota=cosa/sina
c      cotb=one/cdtan(beta)
c      expia=cdexp(ione*alpha)
c      s=(beta*beta-one+beta*cotb)/(beta*(one-beta*cotb))
c      a=(alpha*alpha-one+alpha*cota)/(alpha*(one-alpha*cota))
c      b=-(one-alpha*alpha+ione*alpha)/(alpha*(one+ione*alpha))
c      c1=ione*expia*(sina-alpha*cosa)/(one+alpha)
c      c2=-ione*expia*((alpha*alpha-one)*sina+alpha*cosa)/
c     &   (one-alpha*alpha+ione*alpha)
c      alphae=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c2*(one-s/(a*n))/(one-s/(b*n))
c      alpham=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c1*(one-a/(s*n))/(one-b/(s*n))
c      tmp=16.d0*pi/3.d0*x*dimag(alphae+alpham)

      qabs_c=tmp
      
      return
      end
c      
c
c
      complex*16 function eps_c(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      real*8 n,k
      complex*16 e
      
      real*8 rad,temp
cc      common /dust/ rad,temp

      real*8   ngraph,kgraph
      external ngraph,kgraph
      
      
      n=ngraph(lambda,rad,temp)
      k=kgraph(lambda,rad,temp)
      e=dcmplx(n,k)
      e=e*e
      
      eps_c=e
      
      return
      end

      
c=========================================================================
      real*8 function qabs_cu(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 x,rho1,rho2,rhoabs
      real*8 tmp,tmp1,tmp2
      real*8 n1,n2,e1,e2
      real*8 c0
      complex*16 n,e
      complex*16 alphae,alpham,alpha,beta,s,a,b,c1,c2
      complex*16 sina,cosa,cota,cotb,expia
      
      complex*16 eps_cu
      external eps_cu
      
      e=eps_cu(lambda,rad,temp)
      n=cdsqrt(e)
      
      n1=dreal(n)
      n2=dimag(n)
      e1=dreal(e)
      e2=dimag(e)
      
      x=2.d0*pi*rad/lambda
      rho1=x*n1
      rho2=x*n2
      rhoabs=rho1*rho1+rho2*rho2
c --- separate small/large x model ---
      if (x.le.1.d0) then
       beta=n*x
       cotb=cdcos(beta)/cdsin(beta)
       s=(beta*beta-1.d0+beta*cotb)/(beta*(1.d0-beta*cotb))
       alphae=3.d0/(4.d0*pi)*(1.d0-x*s/(2.d0*n))/(1.d0+x*s/n)
       tmp1=16.d0*pi/3.d0*x*dimag(alphae)
c       tmp1=12.d0*x*e2/((e1+2.d0)**2+e2*e2)
       tmp1=tmp1+6.d0*x*
     &     ((rho1*dsinh(2.d0*rho2)+rho2*dsin(2.d0*rho1))/
     &     (rhoabs*(dcosh(2.d0*rho2)-dcos(2.d0*rho1)))
     &     -2.d0*rho1*rho2/(rhoabs*rhoabs))
       tmp=tmp1
      else
       if (n2.ge.1.d0) then
        tmp2=16.d0/3.d0*n1/(n1*n1+n2*n2)
c        tmp2=4.d0*n1/(n1*n1+n2*n2+1.d0+2.d0*n1)
        tmp=tmp2
       else
c        tmp=8.d0/3.d0*x*n2
        c0=4.d0*n2*x
        if (c0.le.1.d-4) then
         tmp=c0*(2.d0/3.d0-0.25d0*c0)
        else
 	 tmp=1.d0-2.d0*(1.d0-(1.d0+c0)*dexp(-c0))/(c0*c0)
        endif
       endif
      endif
c --- dipol model for any x ---
c      alpha=x
c      beta=x*n
c      sina=cdsin(alpha)
c      cosa=cdcos(alpha)
c      cota=cosa/sina
c      cotb=one/cdtan(beta)
c      expia=cdexp(ione*alpha)
c      s=(beta*beta-one+beta*cotb)/(beta*(one-beta*cotb))
c      a=(alpha*alpha-one+alpha*cota)/(alpha*(one-alpha*cota))
c      b=-(one-alpha*alpha+ione*alpha)/(alpha*(one+ione*alpha))
c      c1=ione*expia*(sina-alpha*cosa)/(one+alpha)
c      c2=-ione*expia*((alpha*alpha-one)*sina+alpha*cosa)/
c     &   (one-alpha*alpha+ione*alpha)
c      alphae=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c2*(one-s/(a*n))/(one-s/(b*n))
c      alpham=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c1*(one-a/(s*n))/(one-b/(s*n))
c      tmp=16.d0*pi/3.d0*x*dimag(alphae+alpham)
      qabs_cu=tmp
      return
      end
c      
c
c
      complex*16 function eps_cu(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      real*8 omega,omega2,omegap,omegap2,gamma,gamma2
      real*8 e1,e2
      complex*16 e,ed
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 fomegap_cu,feta_cu
      complex*16 epsd_cu
      external fomegap_cu,feta_cu,epsd_cu
      
      omega=2.d0*pi*c/lambda
      omegap=fomegap_cu(temp)
      omegap2=omegap*omegap
      gamma=omegap2*feta_cu(temp)*eps0
      omega2=omega*omega
      gamma2=gamma*gamma
      e1=1.d0-omegap2/(omega2+gamma2)
      e2=omegap2*gamma/(omega*(omega2+gamma2))
      e=dcmplx(e1,e2)
      ed=epsd_cu(omega)
c      e=e+ed
      
      eps_cu=e
      
      return
      end
c     
c
c
      real*8 function feta_cu(temp)
      integer i, m
      parameter (m=18)
      real*8 temp, tmp
      real*8 t(m), r(m)
      data t/
     *        300.0d0 ,  350.0d0 ,  400.0d0 ,  500.0d0 ,  600.0d0 ,
     *        700.0d0 ,  800.0d0 ,  900.0d0 , 1000.0d0 , 1100.0d0 ,
     *       1200.0d0 , 1300.0d0 ,1357.60d0 ,1357.61d0 , 1400.0d0 ,
     *       1500.0d0 , 1600.0d0 , 1700.0d0 /
      data r/
     *       1.725d0  , 2.063d0  , 2.402d0  , 3.090d0  , 3.792d0  ,
     *       4.514d0  , 5.262d0  , 6.041d0  , 6.858d0  , 7.717d0  ,
     *       8.626d0  , 9.592d0  ,10.171d0  , 21.01d0  , 21.43d0  ,
     *       22.42d0  , 23.42d0  , 24.41d0  /
      
      real*8 copt
      data copt / 2.38d0 /
      
      if (temp.lt.t(1)) then
       tmp=r(1)+(r(2)-r(1))*(temp-t(1))/(t(2)-t(1))
      elseif (temp.gt.t(m)) then 
       tmp=r(m)+(r(m)-r(m-1))*(temp-t(m))/(t(m)-t(m-1))
      else
       i=2
10     if (temp.gt.t(i)) then
        i=i+1
        goto 10
       endif
       tmp=r(i-1)+(r(i)-r(i-1))*(temp-t(i-1))/(t(i)-t(i-1))
      endif
      
      feta_cu=copt*tmp*1.d-8   ! [Ohm*m]
      return
      end
c
c
c
      real*8 function fomegap_cu(temp)
      implicit none
      
      real*8 temp
      real*8 tmp
      
      tmp=1.2d+16
      fomegap_cu=tmp
      return
      end
c
c
c
      complex*16 function epsd_cu(omega)
      implicit none
      Use (const)
      
      real*8 omega
      complex*16 tmp
      integer i,n
      parameter (n=4)
      real*8 f(n),gama(n),omegai(n),omegap,omegap2
      data f      /    0.061,     0.104,     0.723,     0.638/
      data gama   /5.743d+14, 1.604d+15, 4.881d+15, 6.540d+15/
      data omegai /4.421d+14, 4.492d+15, 8.052d+15, 1.699d+16/
      data omegap /1.645d+16/
      
      omegap2=omegap*omegap
      tmp=dcmplx(0.d0,0.d0)
      do 10, i=1,n
       tmp=tmp+f(i)*omegap2/((omegai(i)*omegai(i)-omega*omega)
     &                       +ione*omega*gama(i))
10    continue
      tmp=dconjg(tmp)
      
      epsd_cu=tmp
      return
      end
      
c=================================================================
      real*8 function qabs_ss(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 x,rho1,rho2,rhoabs
      real*8 tmp,tmp1,tmp2
      real*8 n1,n2,e1,e2
      real*8 c0
      complex*16 n,e
      complex*16 alphae,alpham,alpha,beta,s,a,b,c1,c2
      complex*16 sina,cosa,cota,cotb,expia
      
      complex*16 eps_ss
      external eps_ss
      
      e=eps_ss(lambda,rad,temp)
      n=cdsqrt(e)
      
      n1=dreal(n)
      n2=dimag(n)
      e1=dreal(e)
      e2=dimag(e)
      
      x=2.d0*pi*rad/lambda
      rho1=x*n1
      rho2=x*n2
      rhoabs=rho1*rho1+rho2*rho2
c --- separate small/large x model ---
      if (x.le.1.d0) then
       beta=n*x
       cotb=cdcos(beta)/cdsin(beta)
       s=(beta*beta-1.d0+beta*cotb)/(beta*(1.d0-beta*cotb))
       alphae=3.d0/(4.d0*pi)*(1.d0-x*s/(2.d0*n))/(1.d0+x*s/n)
       tmp1=16.d0*pi/3.d0*x*dimag(alphae)
c       tmp1=12.d0*x*e2/((e1+2.d0)**2+e2*e2)
       tmp1=tmp1+6.d0*x*
     &     ((rho1*dsinh(2.d0*rho2)+rho2*dsin(2.d0*rho1))/
     &     (rhoabs*(dcosh(2.d0*rho2)-dcos(2.d0*rho1)))
     &     -2.d0*rho1*rho2/(rhoabs*rhoabs))
       tmp=tmp1
      else
       if (n2.ge.1.d0) then
        tmp2=16.d0/3.d0*n1/(n1*n1+n2*n2)
c        tmp2=4.d0*n1/(n1*n1+n2*n2+1.d0+2.d0*n1)
        tmp=tmp2
       else
c        tmp=8.d0/3.d0*x*n2
        c0=4.d0*n2*x
        if (c0.le.1.d-4) then
         tmp=c0*(2.d0/3.d0-0.25d0*c0)
        else
         tmp=1.d0-2.d0*(1.d0-(1.d0+c0)*dexp(-c0))/(c0*c0)
        endif
       endif
      endif
c --- dipol model for any x ---
c      alpha=x
c      beta=x*n
c      sina=cdsin(alpha)
c      cosa=cdcos(alpha)
c      cota=cosa/sina
c      cotb=one/cdtan(beta)
c      expia=cdexp(ione*alpha)
c      s=(beta*beta-one+beta*cotb)/(beta*(one-beta*cotb))
c      a=(alpha*alpha-one+alpha*cota)/(alpha*(one-alpha*cota))
c      b=-(one-alpha*alpha+ione*alpha)/(alpha*(one+ione*alpha))
c      c1=ione*expia*(sina-alpha*cosa)/(one+alpha)
c      c2=-ione*expia*((alpha*alpha-one)*sina+alpha*cosa)/
c     &   (one-alpha*alpha+ione*alpha)
c      alphae=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c2*(one-s/(a*n))/(one-s/(b*n))
c      alpham=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c1*(one-a/(s*n))/(one-b/(s*n))
c      tmp=16.d0*pi/3.d0*x*dimag(alphae+alpham)
      qabs_ss=tmp
      return
      end
c      
c
c
      complex*16 function eps_ss(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      real*8 omega,omega2,omegap,omegap2,gamma,gamma2
      real*8 e1,e2
      complex*16 e
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 fomegap_ss,feta_ss
      external fomegap_ss,feta_ss
      
      omega=2.d0*pi*c/lambda
      omegap=fomegap_ss(temp)
      omegap2=omegap*omegap
      gamma=omegap2*feta_ss(temp)*eps0
      omega2=omega*omega
      gamma2=gamma*gamma
      e1=1.d0-omegap2/(omega2+gamma2)
      e2=omegap2*gamma/(omega*(omega2+gamma2))
      e=dcmplx(e1,e2)
      
      eps_ss=e
      
      return
      end
c     
c
c
      real*8 function feta_ss(temp)
      integer i, m
      parameter (m=7)
      real*8 temp, tmp
      real*8 t(m), r(m)
      data t/
     *        298.0d0 ,  373.0d0 ,  473.0d0 ,  673.0d0 ,  873.0d0 ,
     *       1073.0d0 , 1273.0d0 /
      data r/
     *       72.00d0  , 77.60d0  , 85.00d0  , 97.60d0  , 107.2d0  ,
     *       114.1d0  , 119.6d0  /
      
      real*8 copt
      data copt / 1.00d0 /
            if (temp.lt.t(1)) then
       tmp=r(1)+(r(2)-r(1))*(temp-t(1))/(t(2)-t(1))
      elseif (temp.gt.t(m)) then 
       tmp=r(m)+(r(m)-r(m-1))*(temp-t(m))/(t(m)-t(m-1))
      else
       i=2
10     if (temp.gt.t(i)) then
        i=i+1
        goto 10
       endif
       tmp=r(i-1)+(r(i)-r(i-1))*(temp-t(i-1))/(t(i)-t(i-1))
      endif
      
      feta_ss=copt*tmp*1.d-8   ! [Ohm*m]
      return
      end
c
c
c
      real*8 function fomegap_ss(temp)
      implicit none
      
      real*8 temp
      real*8 tmp
      
      tmp=2.266d+16
      fomegap_ss=tmp
      return
      end
      
c================================================================
      real*8 function qabs_mo(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 x,rho1,rho2,rhoabs
      real*8 tmp,tmp1,tmp2
      real*8 n1,n2,e1,e2
      real*8 c0
      complex*16 n,e
      complex*16 alphae,alpham,alpha,beta,s,a,b,c1,c2
      complex*16 sina,cosa,cota,cotb,expia
      
      complex*16 eps_mo
      external eps_mo
      
      e=eps_mo(lambda,rad,temp)
      n=cdsqrt(e)
      
      n1=dreal(n)
      n2=dimag(n)
      e1=dreal(e)
      e2=dimag(e)
      
      x=2.d0*pi*rad/lambda
      rho1=x*n1
      rho2=x*n2
      rhoabs=rho1*rho1+rho2*rho2
c --- separate small/large x model ---
      if (x.le.1.d0) then
       beta=n*x
       cotb=cdcos(beta)/cdsin(beta)
       s=(beta*beta-1.d0+beta*cotb)/(beta*(1.d0-beta*cotb))
       alphae=3.d0/(4.d0*pi)*(1.d0-x*s/(2.d0*n))/(1.d0+x*s/n)
       tmp1=16.d0*pi/3.d0*x*dimag(alphae)
c       tmp1=12.d0*x*e2/((e1+2.d0)**2+e2*e2)
       tmp1=tmp1+6.d0*x*
     &     ((rho1*dsinh(2.d0*rho2)+rho2*dsin(2.d0*rho1))/
     &     (rhoabs*(dcosh(2.d0*rho2)-dcos(2.d0*rho1)))
     &     -2.d0*rho1*rho2/(rhoabs*rhoabs))
       tmp=tmp1
      else
       if (n2.ge.1.d0) then
        tmp2=16.d0/3.d0*n1/(n1*n1+n2*n2)
c        tmp2=4.d0*n1/(n1*n1+n2*n2+1.d0+2.d0*n1)
        tmp=tmp2
       else
c        tmp=8.d0/3.d0*x*n2
        c0=4.d0*n2*x
        if (c0.le.1.d-4) then
         tmp=c0*(2.d0/3.d0-0.25d0*c0)
        else
  	 tmp=1.d0-2.d0*(1.d0-(1.d0+c0)*dexp(-c0))/(c0*c0)
        endif
       endif
      endif
c --- dipol model for any x ---
c      alpha=x
c      beta=x*n
c      sina=cdsin(alpha)
c      cosa=cdcos(alpha)
c      cota=cosa/sina
c      cotb=one/cdtan(beta)
c      expia=cdexp(ione*alpha)
c      s=(beta*beta-one+beta*cotb)/(beta*(one-beta*cotb))
c      a=(alpha*alpha-one+alpha*cota)/(alpha*(one-alpha*cota))
c      b=-(one-alpha*alpha+ione*alpha)/(alpha*(one+ione*alpha))
c      c1=ione*expia*(sina-alpha*cosa)/(one+alpha)
c      c2=-ione*expia*((alpha*alpha-one)*sina+alpha*cosa)/
c     &   (one-alpha*alpha+ione*alpha)
c      alphae=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c2*(one-s/(a*n))/(one-s/(b*n))
c      alpham=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c1*(one-a/(s*n))/(one-b/(s*n))
c      tmp=16.d0*pi/3.d0*x*dimag(alphae+alpham)
      qabs_mo=tmp
      return
      end
c      
c
c
      complex*16 function eps_mo(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      real*8 omega,omega2,omegap,omegap2,gamma,gamma2
      real*8 e1,e2
      complex*16 e
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 fomegap_mo,feta_mo
      external fomegap_mo,feta_mo
      
      omega=2.d0*pi*c/lambda
      omegap=fomegap_mo(temp)
      omegap2=omegap*omegap
      gamma=omegap2*feta_mo(temp)*eps0
      omega2=omega*omega
      gamma2=gamma*gamma
      e1=1.d0-omegap2/(omega2+gamma2)
      e2=omegap2*gamma/(omega*(omega2+gamma2))
      e=dcmplx(e1,e2)
      
      eps_mo=e
      
      return
      end
c     
c
c
      real*8 function feta_mo(temp)
      implicit none
      
      real*8 temp
      
      integer n,i
      real*8 tmp
      parameter (n=6)
      real*8 coeff1(n),coeff2(n),coeff3(n),coeff4(n)
      data coeff1 / 7.00000d-04, 0.00000d-00, 1.48000d-06, 0.00000d-00,
     &              0.00000d-00, 4.12000d-10/
      data coeff2 / 1.07700d-01,-8.19320d-03, 1.67780d-04, 0.00000d-00,
     &             -1.07940d-09, 0.00000d-00/
      data coeff3 / 1.18700d-01,-1.61590d-02, 3.54610d-04,-1.38840d-06,
     &              1.98030d-09, 0.00000d-00/
      data coeff4 /-1.70210d-00, 2.33190d-02, 2.55070d-06,-2.59300d-10,
     &              0.00000d-00, 0.00000d-00/
      
      real*8 copt
      data copt / 1.28d0 /
      
      tmp=0.d0
      do 10, i=1,n
       if     (temp.le.1.d0) then
        tmp=tmp+coeff1(i)*temp**(i-1)
       elseif (temp.le.30.d0) then
        tmp=tmp+coeff1(i)*temp**(i-1)
       elseif (temp.le.100.d0) then
        tmp=tmp+coeff2(i)*temp**(i-1)
       elseif (temp.le.250.d0) then
        tmp=tmp+coeff3(i)*temp**(i-1)
       elseif (temp.le.2894.d0) then
        tmp=tmp+coeff4(i)*temp**(i-1)
       else
        tmp=97.d0
       endif
10    continue
      
      feta_mo=copt*tmp*1.d-08
      return
      end
c
c
c
      real*8 function fomegap_mo(temp)
      implicit none
      
      real*8 temp
      real*8 tmp
      
      tmp=1.13d+16
      fomegap_mo=tmp
      return
      end
      
      
c=====================================================================
      real*8 function qabs_w(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 x,rho1,rho2,rhoabs
      real*8 tmp,tmp1,tmp2
      real*8 n1,n2,e1,e2
      real*8 c0
      complex*16 n,e
      complex*16 alphae,alpham,alpha,beta,s,a,b,c1,c2
      complex*16 sina,cosa,cota,cotb,expia
      
      complex*16 eps_w
      external eps_w
      
      e=eps_w(lambda,rad,temp)
      n=cdsqrt(e)
      
      n1=dreal(n)
      n2=dimag(n)
      e1=dreal(e)
      e2=dimag(e)
      
      x=2.d0*pi*rad/lambda
      rho1=x*n1
      rho2=x*n2
      rhoabs=rho1*rho1+rho2*rho2
c --- separate small/large x model ---
      if (x.le.1.d0) then
       beta=n*x
       cotb=cdcos(beta)/cdsin(beta)
       s=(beta*beta-1.d0+beta*cotb)/(beta*(1.d0-beta*cotb))
       alphae=3.d0/(4.d0*pi)*(1.d0-x*s/(2.d0*n))/(1.d0+x*s/n)
       tmp1=16.d0*pi/3.d0*x*dimag(alphae)
c       tmp1=12.d0*x*e2/((e1+2.d0)**2+e2*e2)
       tmp1=tmp1+6.d0*x*
     &     ((rho1*dsinh(2.d0*rho2)+rho2*dsin(2.d0*rho1))/
     &     (rhoabs*(dcosh(2.d0*rho2)-dcos(2.d0*rho1)))
     &     -2.d0*rho1*rho2/(rhoabs*rhoabs))
       tmp=tmp1
      else
       if (n2.ge.1.d0) then
        tmp2=16.d0/3.d0*n1/(n1*n1+n2*n2)
c        tmp2=4.d0*n1/(n1*n1+n2*n2+1.d0+2.d0*n1)
        tmp=tmp2
       else
c        tmp=8.d0/3.d0*x*n2
        c0=4.d0*n2*x
        if (c0.le.1.d-4) then
         tmp=c0*(2.d0/3.d0-0.25d0*c0)
        else
	 tmp=1.d0-2.d0*(1.d0-(1.d0+c0)*dexp(-c0))/(c0*c0)
        endif
       endif
      endif
c --- dipol model for any x ---
c      alpha=x
c      beta=x*n
c      sina=cdsin(alpha)
c      cosa=cdcos(alpha)
c      cota=cosa/sina
c      cotb=one/cdtan(beta)
c      expia=cdexp(ione*alpha)
c      s=(beta*beta-one+beta*cotb)/(beta*(one-beta*cotb))
c      a=(alpha*alpha-one+alpha*cota)/(alpha*(one-alpha*cota))
c      b=-(one-alpha*alpha+ione*alpha)/(alpha*(one+ione*alpha))
c      c1=ione*expia*(sina-alpha*cosa)/(one+alpha)
c      c2=-ione*expia*((alpha*alpha-one)*sina+alpha*cosa)/
c     &   (one-alpha*alpha+ione*alpha)
c      alphae=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c2*(one-s/(a*n))/(one-s/(b*n))
c      alpham=9.d0*ione/(8.d0*pi*alpha**3.d0)*
c     &       c1*(one-a/(s*n))/(one-b/(s*n))
c      tmp=16.d0*pi/3.d0*x*dimag(alphae+alpham)
      qabs_w=tmp
      return
      end
c      
c
c
      complex*16 function eps_w(lambda,rad,temp)
      implicit none
      
      Use (const)
      
      real*8 lambda
      real*8 omega,omega2,omegap,omegap2,gamma,gamma2
      real*8 e1,e2
      complex*16 e,ed
      
      real*8 rad,temp
cc      common /dust/ rad,temp
      
      real*8 fomegap_w,feta_w
      complex*16 epsd_w
      external fomegap_w,feta_w,epsd_w
      
      omega=2.d0*pi*c/lambda
      omegap=fomegap_w(temp)
      omegap2=omegap*omegap
      gamma=omegap2*feta_w(temp)*eps0
      omega2=omega*omega
      gamma2=gamma*gamma
      e1=1.d0-omegap2/(omega2+gamma2)
      e2=omegap2*gamma/(omega*(omega2+gamma2))
      e=dcmplx(e1,e2)
      ed=epsd_w(omega)
c      e=e+ed
      
      eps_w=e
      
      return
      end
c     
c
c
      real*8 function feta_w(temp)
      implicit none
      
      real*8 temp
      
      integer n,k,i,j
      real*8 tmp
      parameter (n=6, k=4)
      real*8 coeff1(n),coeff2(n),coeff3(n),coeff4(n)
      data coeff1 / 1.50000d-05, 0.00000d-00, 7.00000d-07, 0.00000d-00,
     &              0.00000d-00, 5.20000d-10/
      data coeff2 / 1.44070d-01,-1.16651d-02, 2.41437d-04, 0.00000d-00,
     &             -3.66335d-09, 0.00000d-00/
      data coeff3 /-1.06871d-00, 2.06884d-02, 1.27971d-06, 8.53101d-09,
     &             -5.14195d-12, 0.00000d-00/
      data coeff4 /-1.72573d-00, 2.14350d-02, 5.74811d-06,-1.13698d-09,
     &              1.11670d-13, 0.00000d-00/
      
      real*8 t(k), r(k)
      data t /3660.d0, 4000.d0, 4500.d0, 5000.d0/
      data r / 131.d0,  135.d0,  151.d0,  160.d0/
     
      real*8 copt
      data copt / 2.10d0 /
      
      tmp=0.d0
      do 10, i=1,n
       if     (temp.le.1.d0) then
        tmp=tmp+coeff1(i)*temp**(i-1)
       elseif (temp.le.40.d0) then
        tmp=tmp+coeff1(i)*temp**(i-1)
       elseif (temp.le.90.d0) then
        tmp=tmp+coeff2(i)*temp**(i-1)
       elseif (temp.le.750.d0) then
        tmp=tmp+coeff3(i)*temp**(i-1)
       elseif (temp.le.3660.d0) then
        tmp=tmp+coeff4(i)*temp**(i-1)
       elseif (temp.le.t(k)) then
        j=1
15      j=j+1
        if (temp.gt.t(j)) goto 15
        tmp=r(j-1)+(temp-t(j-1))*(r(j)-r(j-1))/(t(j)-t(j-1))
       else
        tmp=r(k)+(r(k)-r(k-1))*(temp-t(k))/(t(k)-t(k-1))
       endif
10    continue
      
      feta_w=copt*tmp*1.d-08
      return
      end
c
c
c
      real*8 function fomegap_w(temp)
      implicit none
      
      real*8 temp
      real*8 tmp
      
      tmp=9.12d+15
      fomegap_w=tmp
      return
      end
c
c
c
      complex*16 function epsd_w(omega)
      implicit none
      
      Use (const)
      
      real*8 omega
      complex*16 tmp
      integer i,n
      parameter (n=4)
      real*8 f(n),gama(n),omegai(n),omegap,omegap2
      data f      /    0.054,     0.166,     0.706,     2.590/
      data gama   /8.052d+14, 1.946d+15, 5.062d+15, 8.866d+15/
      data omegai /1.525d+15, 2.912d+15, 5.439d+15, 1.139d+16/
      data omegap /2.008d+16/
      
      omegap2=omegap*omegap
      tmp=dcmplx(0.d0,0.d0)
      do 10, i=1,n
       tmp=tmp+f(i)*omegap2/((omegai(i)*omegai(i)-omega*omega)
     &                       +ione*omega*gama(i))
10    continue
      tmp=dconjg(tmp)
      
      epsd_w=tmp
      return
      end
      
      
c================================================================
      real*8 function fomegap_mat(iz)
      implicit none
      integer iz
      real*8 tmp
      integer n 
      parameter (n=92)
      real*8 f(n)
      data f / ! dummy data except for Li Be Cu SS Mo W
     ,    1.d16,   1.d16,   1.d16,  8.15d15,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    2.266d16,1.d16,   1.d16,   1.2d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16, 1.13d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,  9.12d15,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16,   1.d16,    1.d16,   1.d16,
     ,    1.d16,   1.d16/

      tmp=f(iz)
      fomegap_mat=tmp

      return
      end
c===============================================================



