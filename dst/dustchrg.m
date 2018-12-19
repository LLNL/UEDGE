c Calculation of Fc and Fo 
c 
c
c=======================================================================
      subroutine dstdrag

      implicit none

c      Use (dustdim)
      Use (dustcgm)
      Use (dustinp)
      Use (dusttrj)
      Use (dustcur)
      Use (dustcom)
      Use (dustcntrl)
      Use (std_cns)

      integer m
      real*8 vn,xi,ff,tmp
      real*8 tetoti,titote
      real*8 aui2,aui2exp,xerf
      real*8 titotz(nimpzm)
      real*8 vniz(nimpzm),auz2(nimpzm),auz2exp(nimpzm),zerf(nimpzm),xiz
      real*8 std_erf,onedrag,zdrag

      common /nrmcns/ xerf,aui2,aui2exp
      common /txtoty/ tetoti,titote
      common /zimpur/ zerf,auz2,auz2exp,titotz

      tetoti=te/ti
      titote=ti/te

      vn=2.d0*tierg*invmpnrm
      aui=(uzp-vz)*(uzp-vz)
      aui=aui+(urp-vr)*(urp-vr)
      aui=aui+(uthp-vt)*(uthp-vt)
      aui2=aui/vn
      aui=dsqrt(aui2)
      aui2exp=dexp(-aui2)
      vn=dsqrt(vn)
      xerf=std_erf(aui)

      if (jimp.gt.0) then
       do 10, m=1,nimpz
        vniz(m)=2.d0*tiz(m)*ev2erg/(mz(m)*mol1nrm)
        auz(m)=(uziz(m)-vz)*(uziz(m)-vz)
        auz(m)=auz(m)+(uriz(m)-vr)*(uriz(m)-vr)
        auz(m)=auz(m)+(uthiz(m)-vt)*(uthiz(m)-vt)
        auz2(m)=auz(m)/vniz(m)
        auz(m)=dsqrt(auz2(m))
        auz2exp(m)=dexp(-auz2(m))
        vniz(m)=dsqrt(vniz(m))
        zerf(m)=std_erf(auz(m))
        titotz(m)=ti/tiz(m)
10     continue
      endif

      xi=lambda*tetoti

      call calxi(xi)

      lambda=xi*titote

      ff=onedrag(xi)

      if (aui. gt. 0.d0) then
       tmp=1.d0/(aui*vn)
       ffz=ff*(uzp -vz)*tmp
       ffr=ff*(urp -vr)*tmp
       fft=ff*(uthp-vt)*tmp
      else
       ffz=0.d0
       ffr=0.d0
       fft=0.d0
      endif

      if (jimp.gt.0) then
       do 20, m=1,nimpz
        xiz=xi*titotz(m)*zimp(m)
        ff=zdrag(m,xiz)
        if (auz(m) .gt. 0.d0) then
         tmp=1.d0/(auz(m)*vniz(m))
         ffzz(m)=ff*(uziz(m) -vz)*tmp
         ffrz(m)=ff*(uriz(m) -vr)*tmp
         fftz(m)=ff*(uthiz(m)-vt)*tmp
        else
         ffzz(m)=0.d0
         ffrz(m)=0.d0
         fftz(m)=0.d0
        endif
20     continue
      endif

      end
c==================================================================
       real*8 function onedrag(xi)

       implicit none

c       Use (dustdim)
       Use (dustcgm)
       Use (dusttrj)
       Use (dustcur)
       Use (std_cns)

c------
       real*8 xi
       real*8 fC,fO,ff
       real*8 dnume1,dnume2,dnume3
       real*8 c1,c2,c3,c4
       real*8 gu,xlng,coeff
       real*8 aui2,aui2exp,xerf,xisqr
       real*8 coulog,std_erf
        
       common /nrmcns/ xerf,aui2,aui2exp

       coeff=2.d-8*np*rd*rd*tierg

c---------------------
c Coefficient of Fc 
c --------------------
       if (xi .ge. 0.d0) then 
        if(aui .lt. 5.d-3)then
         fC=pisqrt*onethird*(8.d0+4.d0*xi)*aui
        else
         dnume1=aui*(aui2+.5d0+xi)*aui2exp
         dnume2=aui2*(aui2+1.d0)
     +         -.25d0-(.5d0-aui2)*xi
         dnume2=dnume2*pisqrt*xerf
         dnume3=(dnume1+dnume2)*pisqrt
         fC=dnume3/aui2
        endif
       else
        if(aui.lt.5.d-3)then
         fC=pisqrt*onethird*(8.d0-4.d0*xi)*aui*dexp(xi)
        else
         xisqr=dsqrt(-xi)
         c3=2.d0*aui2
         c1=1.d0+c3
         c2=1.d0-c3
         c3=c1
         c4=xisqr*c2/aui
         c1=c3+c4
         c2=c3-c4
         c3=aui+xisqr
         c4=aui-xisqr
         dnume1=dexp(-c3*c3)*c1+dexp(-c4*c4)*c2
         dnume1=dnume1*pisqrt/aui
         c1=1.d0+2.d0*(xi+aui2)
         c2=1.d0+2.d0*(xi-aui2)
         dnume2=c1-.5d0*c2/aui2
         dnume2=pi0*dnume2*(std_erf(c3)+std_erf(c4))
         fC=.25d0*(dnume1+dnume2)
        endif
       endif
c-----------------------------------------------------------------------
c Coefficient of Fo
c-----------------------------------------------------------------------

        xlng=coulog(xi)

        if(aui .lt. 5.d-3) then 
         gu=2.d0*onethird*aui*invpisqrt
        else 
         dnume1=.5d0*xerf-aui*aui2exp*invpisqrt
         gu=dnume1/aui2
        endif

        fO=xi*xi*(2.d0*pi0)*xlng*gu
        ff=coeff*(fC+fO)
        onedrag=ff
       return
       end 


c==================================================================
      real*8 function zdrag(m,xiz)

      implicit none

c      Use (dustdim)
      Use (dustcgm)
      Use (dusttrj)
      Use (dustcur)
      Use (std_cns)

c------
      integer m
      real*8 xiz

      real*8 fC,fO,ff
      real*8 dnume1,dnume2,dnume3
      real*8 c1,c2,c3,c4,ct
      real*8 gu,xlng,coeff,xizsqr
      real*8 titotz(nimpzm)
      real*8 auz2(nimpzm),auz2exp(nimpzm),zerf(nimpzm)
      real*8 zcoulog,std_erf
      real*8 auzt,auz2t,auz2expt,zerft,invauz,invauz2

      common /zimpur/ zerf,auz2,auz2exp,titotz

      coeff=2.d-8*niz(m)*rd*rd*tizerg(m)
      
      auzt=auz(m)
      auz2t=auz2(m)
      auz2expt=auz2exp(m)
      zerft=zerf(m)
      invauz=1.d0/auzt
      invauz2=invauz*invauz

c---------------------
c Coefficient of Fc 
c --------------------
      if (xiz .ge. 0.d0) then 
       if(auzt .lt. 5.d-3)then
        fC=pisqrt*onethird*(8.d0+4.d0*xiz)*auzt
       else
        dnume1=auzt*(auz2t+.5d0+xiz)*auz2expt
        dnume2=auz2t*(auz2t+1.d0)
     +        -.25d0-(.5d0-auz2t)*xiz
        dnume2=dnume2*pisqrt*zerft
        dnume3=(dnume1+dnume2)*pisqrt
        fC=dnume3*invauz2
       endif
      else
       if(auzt .lt. 5.d-3)then
        fC=pisqrt*onethird*(8.d0-4.d0*xiz)*auzt*dexp(xiz)
       else
        xizsqr=dsqrt(-xiz)
        c3=2.d0*auz2t
        c1=1.d0+c3
        c2=1.d0-c3
        c3=c1
        c4=xizsqr*c2*invauz
        c1=c3+c4
        c2=c3-c4
        c3=auzt+xizsqr
        c4=auzt-xizsqr
        dnume1=dexp(-c3*c3)*c1+dexp(-c4*c4)*c2
        dnume1=dnume1*pisqrt*invauz
        c1=1.d0+xiz+xiz
        c2=auz2t+auz2t
        ct=c1-c2
        c1=c1+c2
        c2=ct
        dnume2=c1-.5d0*c2*invauz2
        dnume2=pi0*dnume2*(std_erf(c3)+std_erf(c4))
        fC=.25d0*(dnume1+dnume2)
       endif
      endif
c--------------------
c Coefficient of Fo
c--------------------

      xlng=zcoulog(m,xiz)

      if (auzt .lt. 5.d-3) then 
       gu=twothirds*auzt*invpisqrt
      else 
       dnume1=.5d0*zerft-auzt*auz2expt*invpisqrt
       gu=dnume1*invauz2
      endif

      fO=xiz*xiz*twopi*xlng*gu
      ff=coeff*(fC+fO)

      zdrag=ff

      return
      end 


C========================================================================
c      calculation of coulomb logarithm
c------------------------------------------------------------------
      real*8 function coulog(xi)

      implicit none

c      Use (dustdim)
      Use (dustcgm)
      Use (dusttrj)
      Use (std_cns)

      real*8 xi
c--
      real*8 debyee,debyes,invdbys,b90,b902
      real*8 eta,a,g,rdnrm
      real*8 tetoti,titote
      real*8 aui2,aui2exp,xerf
      real*8 xlng

      common /nrmcns/ xerf,aui2,aui2exp
      common /txtoty/ tetoti,titote

      rdnrm=rd*1.d-4
c--------
      g=1.d0/(3.d0+2.d0*aui2)
      debyee=5.873d8*dsqrt(teerg/ne)
      invdbys=1.d0+3.d0*tetoti*g
      invdbys=debyee*debyee/invdbys
      debyes=dsqrt(invdbys)
      a=rdnrm/debyes
      if (xi.ge.0.d0) then
       eta=(1.d0+a*(1.d0+.40824829d0*dsqrt(tetoti)))*debyes
       b90=rdnrm*xi*g
       b902=b90*b90
       xlng=(b902+eta*eta)/(b902+rdnrm*rdnrm)
       xlng=.5d0*dlog(xlng)
      else
       xlng=-dlog(a)+xi*(0.306852819d0-xi*9.657359d-2)   ! (1-ln2)=0.3068..
      endif
      xlng=max(0.d0,xlng)
        
      coulog=xlng

      return
      end
c------------------------------------------------------------------
      real*8 function zcoulog(m,xiz)

      implicit none

c      Use (dustdim)
      Use (dustcgm)
      Use (dusttrj)
      Use (std_cns)

      integer m
      real*8 xiz
c--
      real*8 debyee,debyes,invdbys,b90,b902
      real*8 eta,a,g,rdnrm
      real*8 tetoti,titote
      real*8 aui2,aui2exp,xerf
      real*8 xlng
      real*8 titotz(nimpzm)
      real*8 auz2(nimpzm),auz2exp(nimpzm),zerf(nimpzm)

      common /nrmcns/ xerf,aui2,aui2exp
      common /zimpur/   zerf,auz2,auz2exp,titotz
      common /txtoty/ tetoti,titote

      rdnrm=rd*1.d-4
c--------
      g=1.d0/(3.d0+2.d0*auz2(m))
      debyee=3.4492d17*teerg/ne
      invdbys=1.d0+3.d0*tetoti*titotz(m)*g
      invdbys=debyee/invdbys
      debyes=dsqrt(invdbys)
      a=rdnrm/debyes
      if (xiz.ge.0.d0) then
       eta=(1.d0+a*(1.d0+.40824829d0*dsqrt(tetoti*titotz(m))))*debyes
       b90=rdnrm*xiz*g
       b902=b90*b90
       xlng=(b902+eta*eta)/(b902+rdnrm*rdnrm)
       xlng=.5d0*dlog(xlng)
      else
       xlng=-dlog(a)+xiz*(0.306852819d0-xiz*9.657359d-2)
      endif
      xlng=max(0.d0,xlng)
        
      zcoulog=xlng

      return
      end

c========================================================================
c      calculation for xi=-phi/Ti
c------------------------------------------------------------------------
      subroutine calxi(xi)

      implicit none
        
c      Use (dustdim)
      Use (dustinp)
      Use (dustcgm)
      Use (dusttrj)
      Use (dustcur)
      Use (dustcom)
      Use (dustcntrl)
      Use (std_cns)
      Use (psitab)
c      Use (psitab_s)
      Use (psiparloc)

      real*8 xi

      integer m
      integer itr,maxitr
      real*8 abserr,err,deltaxi,xirng,ximin,ximax,xiold
      real*8 curr,game0,gamp0,gamthe0,gamz0(nimpzm)
      real*8 planck3,boltz
      real*8 c1,wfunc,wfuncti,boltd,psb0,psb
      real*8 titotd
      real*8 nenrm4,npnrm4
      real*8 a1,a2,u2
      real*8 a1z(nimpzm),a2z(nimpzm),u2z(nimpzm)
      real*8 aui2,aui2exp,xerf
      real*8 titotz(nimpzm)
      real*8 auz2(nimpzm),auz2exp(nimpzm),zerf(nimpzm)
      real*8 dstcurr
      external dstcurr

      common /nrmcns/   xerf,aui2,aui2exp
      common /currents/ gamz0,game0,gamp0,gamthe0,wfuncti,
     +                  titotd,a1,a2,u2,a1z,a2z,u2z
      common /zimpur/   zerf,auz2,auz2exp,titotz


      data abserr    /+1.d-2/         ! absolute error
      data deltaxi   /+1.d-1/         ! +- range of xi
      data xirng     /+1.d+2/         ! absolute range of xi 
      data maxitr    /30/             ! max number of itterations
      data boltz     /1.3806503d-16/  ! Boltzmann const
      data planck3   /2.90916139d-79/ ! Planck const ^3
      data psb0      /.25d0/          ! barrier coefficient

      nenrm4=.25d-20*ne
      npnrm4=.25d-20*np
      boltd=boltz*td
      titotd=tierg/boltd
      wfunc=cwfunc*ev2erg
      wfuncti=tierg/(wfunc+wfunc)
      
c      if (td .lt. .5d0*tsub) then
c       psb=psb0
c      else
c       if (td .lt. tsub) then
c        psb=psb0+(2.d0*td/tsub-1.d0)*(1.d0-psb0)
c       else
c        psb=1.d0
c       endif
c      endif

       call psi_tdmesh_int(td)
       call psi_tdmesh_vall
       psb=pbtd


      c1=4.d-20*psb*pi0*menrm/planck3
      gamthe0=c1*boltd*boltd*dexp(-wfunc/boltd)
      game0=nenrm4*dsqrt(8.d0*teerg/(pi0*menrm))*ksip
      gamp0=npnrm4*dsqrt(2.d0*tierg*invmpnrm)*ksip

      if (aui .lt. 5.d-3) then
       a1=invpisqrt*(1.d0-aui2)
       a2=2.d0*invpisqrt*(1.d0+onethird*aui2)
      else
       a1=aui*xerf+aui2exp*invpisqrt
       a2=xerf/aui
       u2=1.d0/aui2
      endif

      if (jimp .gt. 0) then
       do 10, m=1,nimpz
        gamz0(m)=.25d-20*niz(m)*ksip
     *          *dsqrt(2.d0*tiz(m)*ev2erg/(mz(m)*mol1nrm))
        if (auz(m) .lt. 5.d-3) then
         a1z(m)=invpisqrt*(1.d0-auz2(m))
         a2z(m)=2.d0*invpisqrt*(1.d0+onethird*auz2(m))
        else
         a1z(m)=auz(m)*zerf(m)+auz2exp(m)*invpisqrt
         a2z(m)=zerf(m)/auz(m)
         u2z(m)=1.d0/auz2(m)
        endif
10     continue
      endif

      ximin=xi-deltaxi
      ximax=xi+deltaxi

      curr=dstcurr(ximax)*dstcurr(ximin)
      if (curr .gt. 0.d0) then
       ximax=xirng
       ximin=-xirng
      endif

      itr=0
100   continue
      xiold=xi
      curr=dstcurr(xi)

      if (curr .gt. 0.d0) then
       ximax=xi
      elseif (curr .lt. 0.d0) then
       ximin=xi
      else
       return
      endif

      xi=.5d0*(ximax+ximin)

      itr=itr+1
      err=abs(xi-xiold)

      if (err .lt. abserr) then
       return
      elseif (itr .ge. maxitr) then
       print *, '***Warning: desired precision of xi was not reached!'
       print *, 'xi=',xi,'xiold=',xiold,'I=',curr
       return
      else
       goto 100
      endif

      end
cc===================================================================

c-----------------------------------------------------------------------
c      Calculation of currents to dust particle
c-----------------------------------------------------------------------
      real*8 function dstcurr (xi)

      implicit none

c      Use (dustdim)
      Use (dustcntrl)
      Use (dustinp)
      Use (dustcur)
      Use (dustcom)
      Use (dustcgm)
      Use (std_cns)
      Use (psitab)
      Use (psitab_s)
      Use (psiparloc)

      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     +            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)

        
      integer m
      real*8 xi,xisq,xitd,u2
      real*8 f1,f2
      real*8 upx,umx
      real*8 tetoti,titote,titotd
      real*8 game0,gamp0,gamthe0,gamz0(nimpzm)
      real*8 a1,a2
      real*8 a1z(nimpzm),a2z(nimpzm),u2z(nimpzm)
      real*8 xite,wfuncti,xiwti,xizz
      real*8 titotz(nimpzm)
      real*8 auz2(nimpzm),auz2exp(nimpzm),zerf(nimpzm)
      real*8 std_erf

      common /currents/ gamz0,game0,gamp0,gamthe0,wfuncti,
     +                  titotd,a1,a2,u2,a1z,a2z,u2z
      common /txtoty/ tetoti,titote
      common /zimpur/ zerf,auz2,auz2exp,titotz

      xite=xi*titote
      call psi_phmesh_int(-xite)
      call psi_mesh_yse

      if(xi .ge. 0.d0) then
       game=game0*dexp(-xite)
       gamp=gamp0*(a1+a2*(.5d0+xi))
       if (jthe .gt. 0) then
        gamthe=gamthe0
       else
        gamthe=0.d0
       endif
       if (jsee .gt. 0) then
        gamsee=game*fysee
       else
        gamsee=0.d0
       endif
       if (jimp .gt. 0) then
        do 10, m=1,nimpz
         gamz(m)=(a1z(m)+a2z(m)*(.5d0+zimp(m)*xi*titotz(m)))
     *          *gamz0(m)
10      continue
       endif
      else
       if (aui .lt. 5.d-5) then
        gamp=2.d0*invpisqrt*dexp(xi)
        gamp=gamp*(1.d0+aui*aui*onethird*(1.d0-2.d0*xi))
       else
        xisq=dsqrt(-xi)
        upx=aui+xisq
        umx=aui-xisq
        f1=1.d0+(.5d0+xi)*u2
        f1=f1*(std_erf(umx)+std_erf(upx))
        f2=upx*dexp(-umx*umx)+umx*dexp(-upx*upx)
        f2=f2*u2*invpisqrt
        gamp=.5d0*aui*(f1+f2)
       endif

       game=game0*(1.d0-titote*xi)
       gamp=gamp0*gamp
       if (jthe .gt. 0) then
        xitd=xi*titotd
        gamthe=gamthe0*(1.d0-xitd)*dexp(xitd)
       else
        gamthe=0.d0
       endif
       if (jsee .gt. 0) then
        xiwti=xi*wfuncti
        gamsee=game*fysee*(1.d0-xiwti)*dexp(xiwti)
       else
        gamsee=0.d0
       endif

       if (jimp .gt. 0) then
        do 40, m=1,nimpz
         xizz=zimp(m)*xi*titotz(m)
         if (auz(m) .lt. 5.d-5) then
          gamz(m)=2.d0*invpisqrt*dexp(xizz)
          gamz(m)=gamz(m)*(1.d0+auz2(m)*onethird*(1.d0-2.d0*xizz))
         else
          xisq=dsqrt(-xizz)
          upx=auz(m)+xisq
          umx=auz(m)-xisq
          f1=1.d0+(.5d0+xizz)*u2z(m)
          f1=f1*(std_erf(umx)+std_erf(upx))
          f2=upx*dexp(-umx*umx)+umx*dexp(-upx*upx)
          f2=f2*u2z(m)*invpisqrt
          gamz(m)=.5d0*auz(m)*(f1+f2)
         endif
         gamz(m)=gamz(m)*gamz0(m)
40      continue
       endif
      endif
c--
      dstcurr=gamp+gamthe+gamsee-game

      if (jimp .gt. 0) then
       do 60, m=1,nimpz
        dstcurr=dstcurr+gamz(m)*zimp(m)
60     continue
      endif

      return
      end
