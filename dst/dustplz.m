      integer function idstplz(namf)
      implicit none
      character*256 namf

c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)

      integer idstplf
      integer ker

      ker=idstplf(namf)
      if(ker .gt. 0) goto 1

      call dstpll

1     idstplz=ker
      return

      end


      subroutine dstpll

      implicit none

      Use (std_cns)
c      Use (dustdim)
      Use (dustcom)
      Use (dustinp)
      Use (dustcgm)
      Use (dustcntrl)

      integer ix,iy
      integer kdof,kdop 
      real*8  ex,ey,efr,efz,eft
      real*8  vps,vys,vcr,vcz,vct
      real*8  pa2dyne
      real*8  nublb

      data kdof/1/
      data kdop/1/
      data pa2dyne/10.d0/

      if(kdof .eq. 1) then
       eft=0.d0
       do  10 ix=1,nx
        do 110 iy=1,ny
         frm(ix,iy)=erm(ix,iy)
         fzm(ix,iy)=ezm(ix,iy)
         ftm(ix,iy)=eft
110     continue
10     continue
      else
       do  1 ix=1,nx
        do 11 iy=1,ny
         ex=exm(ix,iy)
         ey=eym(ix,iy)
         call dstvcp(ix,iy, ex,ey, efr,efz,eft)
         frm(ix,iy)=efr
         fzm(ix,iy)=efz
         ftm(ix,iy)=eft
11      continue
1      continue
      endif

      do  2 ix=1,nx
       do 22 iy=1,ny
        vps=vpm(ix,iy)
        vys=vym(ix,iy)
        call dstvcp(ix,iy,vps,vys,vcr,vcz,vct)
        crm(ix,iy)=vcr
        czm(ix,iy)=vcz
        ctm(ix,iy)=vct
22     continue
2     continue

      do  3 ix=1,nx
       do 33 iy=1,ny
        vps=vpam(ix,iy)
        vys=vyam(ix,iy)
        call dstvcp(ix,iy,vps,vys,vcr,vcz,vct)
        cram(ix,iy)=vcr
        czam(ix,iy)=vcz
        ctam(ix,iy)=vct
33     continue
3     continue

      if(kdop .eq.1) then
       eft=0.d0
       do  4 ix=1,nx
        do 44 iy=1,ny
         pprm(ix,iy)=prm(ix,iy)*pa2dyne
         ppzm(ix,iy)=pzm(ix,iy)*pa2dyne
         pptm(ix,iy)=eft
44      continue
4      continue
      endif

c --- calculate blob collision frequency in cells
      nblb=npm(ixpt1(1)+1,1)
      tblb=tim(ixpt1(1)+1,1)
      do 5, ix=ixpt1(1)+2,ixpt2(1)
       if (nblb .gt. npm(ix,1)) then
        nblb=npm(ix,1)
        tblb=tim(ix,1)
       endif
5     continue
      if (kmesh .eq. 'dnull') then
       do 52, ix=ixpt1(2)+1,ixpt2(2)
        if (nblb .gt. npm(ix,1)) then
         nblb=npm(ix,1)
         tblb=tim(ix,1)
        endif
52     continue
      endif
      tblb=tblb*ev2erg

      tmfp=0.d0
      do 6, ix=1,nx
       do 66, iy=1,ny
        nublb=0.5d0*ksib*nem(ix,iy)*vcvm(ix,iy)/(nblb*rblb)
        if (tmfp .lt. nublb) then
         tmfp=tmfp+nublb
        endif
        blfm(ix,iy)=nublb
66     continue
6     continue
      tmfp=1.d0/tmfp

      end


      subroutine dstvcp(i,j,up,vy,vcr,vcz,vct)
      implicit none
      integer i,j
      real*8 up,vy,vcr,vcz,vct

c      Use (dustdim)
      Use (dustcom)

      real*8 rs1,zs1,rs2,zs2
      real*8 rs3,zs3,rs4,zs4
      real*8 r12,z12,r34,z34
      real*8 r13,z13,r24,z24
      real*8 rper,zper,rpar,zpar
      real*8 cbr,cbz,cbt,cbp,s

      rs1=rrm1(i,j)
      zs1=zzm1(i,j)
      rs2=rrm2(i,j)
      zs2=zzm2(i,j)
      rs3=rrm3(i,j)
      zs3=zzm3(i,j)
      rs4=rrm4(i,j)
      zs4=zzm4(i,j)

      r12=rs2-rs1
      z12=zs2-zs1
      s=dsqrt(r12*r12+z12*z12)
      r12=r12/s
      z12=z12/s
      r34=rs4-rs3
      z34=zs4-zs3
      s=dsqrt(r34*r34+z34*z34)
      r34=r34/s
      z34=z34/s
      rpar=r12+r34
      zpar=z12+z34
      s=dsqrt(rpar*rpar+zpar*zpar)
      zpar=zpar/s
      rpar=rpar/s
      
      r13=rs3-rs1
      z13=zs3-zs1
      s=dsqrt(r13*r13+z13*z13)
      r13=r13/s
      z13=z13/s
      r24=rs4-rs2
      z24=zs4-zs2
      s=dsqrt(r24*r24+z24*z24)
      r24=r24/s
      z24=z24/s
      rper=r13+r24
      zper=z13+z24
      s=dsqrt(rper*rper+zper*zper)
      zper=zper/s
      rper=rper/s

      cbr=(brm1(i,j)+brm2(i,j)+brm3(i,j)+brm4(i,j))*0.25d0
      cbz=(bzm1(i,j)+bzm2(i,j)+bzm3(i,j)+bzm4(i,j))*0.25d0
      cbt=(btm1(i,j)+btm2(i,j)+btm3(i,j)+btm4(i,j))*0.25d0
      cbp=dsqrt(cbr*cbr+cbz*cbz)
      s=dsqrt(cbp*cbp+cbt*cbt)
      cbp=cbp/s
      cbt=cbt/s
      if(cbt .lt. 0.d0) cbt=-cbt

      vcr=rpar*cbp*up+rper*vy
      vcz=zpar*cbp*up+zper*vy
      vct=cbt*up

      end


      subroutine dstcgz(i,j)
      implicit none
      
      integer i,j
      integer l,m
c
c---this function prepares
c   the plasma parameters of a cell (i,j)
c   and stores them in common dstcgp
c
      Use (std_cns)
c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)

        ne = nem (i,j)   !cm-3
        te = tem (i,j)   !eV
        np = npm (i,j)   !cm-3
        ti = tim (i,j)   !eV
       uzp = czm (i,j)   !cm/s
       urp = crm (i,j)   !cm/s
      uthp = ctm (i,j)   !cm/s
        er = frm (i,j)   !V/m
        ez = fzm (i,j)   !V/m
       eth = ftm (i,j)   !V/m
        na = nam (i,j)   !cm-3
        ta = tam (i,j)   !eV
       uza = czam(i,j)   !cm/s
       ura = cram(i,j)   !cm/s
      utha = ctam(i,j)   !cm/s
        nm = nmm (i,j)   !cm-3
        tm = tmm (i,j)   !eV
       uzm = czmm(i,j)   !cm/s
       urm = crmm(i,j)   !cm/s
      uthm = ctmm(i,j)   !cm/s
       prp = pprm(i,j)   !dyne/cm2
       pzp = ppzm(i,j)   !dyne/cm2
      pthp = pptm(i,j)   !dyne/cm2
       vcv = vcvm(i,j)   !cm/s
       blf = blfm(i,j)   !1/s
      teerg=te*ev2erg
      tierg=ti*ev2erg
      taerg=ta*ev2erg

      if (jimp .gt. 0) then
       do 1, l=1,nimp
        nia   (l) = niam (i,j,l)
        tia   (l) = tiam (i,j,l)
        uria  (l) = criam(i,j,l)
        uzia  (l) = cziam(i,j,l)
        uthia (l) = ctiam(i,j,l)
        tiaerg(l) = tia(l)*ev2erg
1      continue
       do 2, m=1,nimpz
        niz   (m) = nizm (i,j,m)
        tiz   (m) = tizm (i,j,m)
        uriz  (m) = crizm(i,j,m)
        uziz  (m) = czizm(i,j,m)
        uthiz (m) = ctizm(i,j,m)
        tizerg(m) = tiz(m)*ev2erg
2      continue
      endif

      end
