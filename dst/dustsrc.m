      subroutine dstsrcini

      implicit none

c      Use (dustdim)
      Use (dustcom)
      Use (dustinp)
      Use (dustcntrl)
      Use (std_cns)
      Use (psitab)
      Use (psitab_s)
      Use (psiparloc)
      
      integer flg
      integer i,j,k,l,m,ix,iy
      integer tp1,tp2,tp3
      real*8 rdistm(frdim),dndistm(frdim),dgdistm(frdim)
      real*8 frfunc(ffrdim+1),fgdist(ffrdim+1)
      real*8 dnsum,dgsum,ro,rwfunc
      real*8 fluxtdc,fluxtsc,fluxtrc,fluxtp,fluxta
      real*8 flxc,flxp,flxa,flxsc,flxrc,dstic
      real*8 ran
      real*8 dl,ds
      real*8 logrd,logradst,fr,ravr3,wrd2
      real*8 vdust

      equivalence (psi_ro,psi_ro_s),(psi_cp,psi_cp_s),(psi_kp,psi_kp_s)
      equivalence (psi_ep,psi_ep_s),(psi_gs,psi_gs_s),(psi_hs,psi_hs_s)
      equivalence (psi_pbtd,psi_pbtd_s)

      equivalence (rn_mx,rn_mx_s),(re_mx,re_mx_s)
      equivalence (ypsp_mx,ypsp_mx_s),(ycsp_mx,ycsp_mx_s),
     *            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)

      srcwi  =0.d0
      srcwo  =0.d0
      srcdib =0.d0
      srcdob =0.d0
      srcpib =0.d0
      srcpob =0.d0
      srcdit =0.d0
      srcdot =0.d0
      srcpit =0.d0
      srcpot =0.d0
      fluxtdc=0.d0
      fluxtsc=0.d0
      fluxtrc=0.d0
      fluxtp =0.d0
      fluxta =0.d0
      dnsum  =0.d0
      dgsum  =0.d0
      nisp   =0
      
      if (jrds .gt. 0) then
c --- calculation of the integral dust initial size distribution
       logradst=dlog(radst)/aln10
       wrd2=wrd*wrd
       frdist(1)=0.d0
       fgdist(1)=0.d0
       frfunc(1)=0.d0
       do 70, i=2,ffrdim+1
        logrd=minlogrd+(i-1)*fdlogrd
        rwfunc=dexp(1.5d0*(logrd-0.5d0*fdlogrd)*aln10)
        logrd=logrd-logradst
        if (jrdf .eq. 1) then
         fr=dexp(-0.5d0*logrd*logrd/wrd2)
        elseif (jrdf .eq. 2) then
         fr=dexp((rpow+1.0d0)*(logrd+logradst-0.5d0*fdlogrd)*aln10)
        else
         fr=dexp(-0.5d0*logrd*logrd/wrd2)
        endif
        frdist(i)=frdist(i-1)+fr
        fgdist(i)=fgdist(i-1)
     +           +fr*dexp(3.d0*(logrd+logradst-0.5d0*fdlogrd)*aln10)
c        frfunc(i)=frfunc(i-1)+fr*rwfunc  ! weighted Rd distribution
        frfunc(i)=frfunc(i-1)+1.d0        ! uniform log(Rd) distribution
70     continue
c --- construction of the dust size distribution
       k=2
       do 50, i=1,frdim
        ran=frfunc(ffrdim+1)*dble(i)/dble(frdim)
54      if ((ran .gt. frfunc(k)) .and. (k .lt. ffrdim+1)) then
         k=k+1
         goto 54
        endif
        logrd=minlogrd+fdlogrd*(k-2+(ran-frfunc(k-1))
     +                          /(frfunc(k)-frfunc(k-1)))
        dndist(i)=frdist(k-1)+(frdist(k)-frdist(k-1))
     +                       *((logrd-minlogrd)/fdlogrd-k+2)
        dgdist(i)=fgdist(k-1)+(fgdist(k)-fgdist(k-1))
     +                       *((logrd-minlogrd)/fdlogrd-k+2)
        if (i .gt. 1) then
         dndist(i)=dndist(i)-dnsum
         dgdist(i)=dgdist(i)-dgsum
        endif
        rdist(i)=dexp(onethird*dlog(dgdist(i)/dndist(i)))
        dnsum=dnsum+dndist(i)
        dgsum=dgsum+dgdist(i)
50     continue
       ravr3=dgsum/dnsum
c --- renormalizing
       do 53, i=1,frdim
        dndist(i)=dndist(i)/dnsum
        dgdist(i)=dgdist(i)/dgsum
53     continue
c --- random mixing
       k=frdim
       do 51, i=1,frdim
        call std_rndd2x(ran)
        j=idint(k*ran)+1
        rdistm (i)=rdist (j)
        dndistm(i)=dndist(j)
        dgdistm(i)=dgdist(j)
        rdist (j)=rdist (k)
        dndist(j)=dndist(k)
        dgdist(j)=dgdist(k)
        k=k-1
51     continue
       do 52, i=1,frdim
        rdist (i)=rdistm (i)
        dndist(i)=dndistm(i)
        dgdist(i)=dgdistm(i)
52     continue
      else
       do 55, i=1,frdim
        rdist(i)=rd0
        dndist(i)=1.d0/frdim
        dgdist(i)=1.d0/frdim
55     continue
       ravr3=rd0*rd0*rd0
      endif

c --- calculation of the integral distribution of initial dust speed
      fvdist(1)=0.d0
      do 60, i=2,ffvdim+1
       vdust=(i-1)*fdvd
       fvdist(i)=fvdist(i-1)
     +          +vdust*vdust*vdust*dexp(-vdust*vdust/(vdt*vdt))
60    continue
      fvdist(ffvdim+1)=1.d0/fvdist(ffvdim+1)
      do 65, i=1,ffvdim
       fvdist(i)=fvdist(i)*fvdist(ffvdim+1)
65    continue
      fvdist(ffvdim+1)=1.d0

c --- calculation of the sputtered fluxes from surfaces
      do 10, i=1,nw
       call dstscll(i,ix,iy,tp1,tp2,tp3)
c --- set wall temperature
       twall(ix,iy)=tw0
       
       if (jsrc .gt. 0) then
        flxsc=0.d0
        flxrc=0.d0
        do 30, l=1,nimp
         flxsc=flxsc+dabs(flia(ix,iy,l))
30      continue
        do 40, m=1,nimpz
         flxrc=flxrc+dabs(fliz(ix,iy,m))
40      continue
        dstic=1.d0
        flxp=dabs(flp(ix,iy))
        flxa=dabs(fla(ix,iy)+flm(ix,iy))
       else
        flxc=0.d0
        flxp=0.d0
        flxa=0.d0
        flg=0
        if (tp1 .eq. 2) then
         dl=dsqrt((rrm3(ix,iy)-rrm4(ix,iy))*(rrm3(ix,iy)-rrm4(ix,iy))+
     +            (zzm3(ix,iy)-zzm4(ix,iy))*(zzm3(ix,iy)-zzm4(ix,iy)))
         ds=pi0*(rrm3(ix,iy)+rrm4(ix,iy))*dl
         flg=1
        elseif (tp1 .eq. 3) then
         dl=dsqrt((rrm1(ix,iy)-rrm2(ix,iy))*(rrm1(ix,iy)-rrm2(ix,iy))+
     +            (zzm1(ix,iy)-zzm2(ix,iy))*(zzm1(ix,iy)-zzm2(ix,iy)))
         ds=pi0*(rrm1(ix,iy)+rrm2(ix,iy))*dl
         flg=1
        elseif (tp1 .eq. 1) then
         if (((tp2 .eq. 1) .and. (tp3 .eq. 1)) .or.
     +       ((tp2 .eq. 2) .and. (tp3 .eq. 2))) then
          dl=dsqrt((rrm1(ix,iy)-rrm3(ix,iy))*(rrm1(ix,iy)-rrm3(ix,iy))+
     +             (zzm1(ix,iy)-zzm3(ix,iy))*(zzm1(ix,iy)-zzm3(ix,iy)))
          ds=pi0*(rrm1(ix,iy)+rrm3(ix,iy))*dl
          flg=1
         elseif (((tp2 .eq. 1) .and. (tp3 .eq. 2)) .or.
     +       ((tp2 .eq. 2) .and. (tp3 .eq. 1))) then
          dl=dsqrt((rrm2(ix,iy)-rrm4(ix,iy))*(rrm2(ix,iy)-rrm4(ix,iy))+
     +             (zzm2(ix,iy)-zzm4(ix,iy))*(zzm2(ix,iy)-zzm4(ix,iy)))
          ds=pi0*(rrm2(ix,iy)+rrm4(ix,iy))*dl
          flg=1
         endif
        endif
        if (flg .eq. 0) then
         ds=0.d0
         print *, '*** Error - dust source is out of the domain!'
c         stop
         return
        endif
	
        sbnd (ix,iy)=ds		! cm2

        call psi_tdmesh_int(twall(ix,iy))
        call dstflxsp(ix,iy,tp1,tp2,tp3,flxc,flxp,flxa)
        flxsc=flxc*sbnd(ix,iy)	! dust_carbon_atoms/s
        flxrc=0.d0
        dstic=1.d0
        flxp =flxp*sbnd(ix,iy)
        flxa =flxa*sbnd(ix,iy)

        flp(ix,iy)=flxp
        fnp(ix,iy)=npm(ix,iy)
       endif

       flg=0
       if (tp2 .eq. 1) then         ! inner
        if (tp1 .eq. 2) then        ! wall
         if (winisp .eq. 0) flg=1
        elseif (tp1 .eq. 3) then    ! PF
         if (tp3 .eq. 1) then       ! bottom
          if (pibnisp .eq. 0) flg=1
         elseif (tp3 .eq. 2) then   ! top
          if (pitnisp .eq. 0) flg=1
         endif
        elseif (tp1 .eq. 1) then    ! plate
         if (tp3 .eq. 1) then       ! bottom
          if (dibnisp .eq. 0) flg=1
         elseif (tp3 .eq. 2) then   ! top
          if (ditnisp .eq. 0) flg=1
         endif
        endif
       elseif (tp2 .eq. 2) then     ! outer
        if     (tp1 .eq. 2) then    ! wall
         if (wonisp .eq. 0) flg=1
        elseif (tp1 .eq. 3) then    ! PF
         if (tp3 .eq. 1) then       ! bottom
          if (pobnisp .eq. 0) flg=1
         elseif (tp3 .eq. 2) then   ! top
          if (potnisp .eq. 0) flg=1
         endif
        elseif (tp1 .eq. 1) then    ! plate
         if (tp3 .eq. 1) then       ! bottom
          if (dobnisp .eq. 0) flg=1
         elseif (tp3 .eq. 2) then   ! top
          if (dotnisp .eq. 0) flg=1
         endif
        endif
       endif
       if (cells .eq. 1) then
        if ((ix .ne. ixs) .or. (iy .ne. iys)) then
         flg=1
        endif
       endif
       if (flg .eq. 1) then
        flxsc=0.d0
        flxrc=0.d0
        flxp =0.d0
        flxa =0.d0
        dstic=0.d0
       endif
       
       fluxc(ix,iy)=dstsc*flxsc+dstrc*flxrc+dstic*flxic	! dust_carbon_atoms/s
       fluxtdc=fluxtdc+fluxc(ix,iy)
       fluxtsc=fluxtsc+flxsc
       fluxtrc=fluxtrc+flxrc
       fluxtp =fluxtp +flxp
       fluxta =fluxta +flxa

       flg=0
       if (tp2 .eq. 1) then         ! inner
        if (tp1 .eq. 2) then        ! wall
         nnispcl(ix,iy)=winisp/frdim
         srcwi=srcwi+fluxc(ix,iy)
         flg=1
        elseif (tp1 .eq. 3) then    ! PF
         if (tp3 .eq. 1) then       ! bottom
          nnispcl(ix,iy)=pibnisp/frdim
          srcpib=srcpib+fluxc(ix,iy)
          flg=1
         elseif (tp3 .eq. 2) then   ! top
          nnispcl(ix,iy)=pitnisp/frdim
          srcpit=srcpit+fluxc(ix,iy)
          flg=1
         endif
        elseif (tp1 .eq. 1) then    ! plate
         if (tp3 .eq. 1) then       ! bottom
          nnispcl(ix,iy)=dibnisp/frdim
          srcdib=srcdib+fluxc(ix,iy)
          flg=1
         elseif (tp3 .eq. 2) then   ! top
          nnispcl(ix,iy)=ditnisp/frdim
          srcdit=srcdit+fluxc(ix,iy)
          flg=1
         endif
        endif
       elseif (tp2 .eq. 2) then     ! outer
        if (tp1 .eq. 2) then        ! wall
         nnispcl(ix,iy)=wonisp/frdim
         srcwo=srcwo+fluxc(ix,iy)
         flg=1
        elseif (tp1 .eq. 3) then    ! PF
         if (tp3 .eq. 1) then       ! bottom
          nnispcl(ix,iy)=pobnisp/frdim
          srcpob=srcpob+fluxc(ix,iy)
          flg=1
         elseif (tp3 .eq. 2) then   ! top
          nnispcl(ix,iy)=potnisp/frdim
          srcpot=srcpot+fluxc(ix,iy)
          flg=1
         endif
        elseif (tp1 .eq. 1) then    ! plate
         if (tp3 .eq. 1) then       ! bottom
          nnispcl(ix,iy)=dobnisp/frdim
          srcdob=srcdob+fluxc(ix,iy)
          flg=1
         elseif (tp3 .eq. 2) then   ! top
          nnispcl(ix,iy)=dotnisp/frdim
          srcdot=srcdot+fluxc(ix,iy)
          flg=1
         endif
        endif
       endif
       if (cells .eq. 1) then
        if ((ix .ne. ixs) .or. (iy .ne. iys)) then
         nnispcl(ix,iy)=0
        endif
       endif
       if (flg .eq. 1) then
        nispcl(ix,iy)=nnispcl(ix,iy)*frdim
        nisp=nisp+nispcl(ix,iy)
       else
        print *, '*Error: non-boundary source cell',ix,iy
c        stop
        return
       endif
        
10    continue
     
      print *, '  Gp_out= ',fluxtp ,'   Ga_inp= ',fluxta
      dstcflx=fluxtdc
      print *, '  Gs_inp= ',fluxtsc,'   Gr_out= ',fluxtrc
      print *, '  Gc_inp= ',dstcflx
      print *, ' Gwi_inp= ',srcwi  ,'  Gwo_inp= ',srcwo
      print *, 'Gdib_inp= ',srcdib ,' Gdob_inp= ',srcdob
      print *, 'Gpib_inp= ',srcpib ,' Gpob_inp= ',srcpob
      print *, 'Gdit_inp= ',srcdit ,' Gdot_inp= ',srcdot
      print *, 'Gpit_inp= ',srcpit ,' Gpot_inp= ',srcpot

      do 14, i=1,nxm
       score(i)=0.d0
       ssptr(i)=0.d0
14    continue

      do 15, i=ixpt1(1)+1,ixpt2(1)
       ix=i
       iy=1
       dl=dsqrt((rrm1(ix,iy)-rrm2(ix,iy))*(rrm1(ix,iy)-rrm2(ix,iy))+
     +          (zzm1(ix,iy)-zzm2(ix,iy))*(zzm1(ix,iy)-zzm2(ix,iy)))
       score(i)=pi0*(rrm1(ix,iy)+rrm2(ix,iy))*dl

       iy=iysptrx(1)
       dl=dsqrt((rrm3(ix,iy)-rrm4(ix,iy))*(rrm3(ix,iy)-rrm4(ix,iy))+
     +          (zzm3(ix,iy)-zzm4(ix,iy))*(zzm3(ix,iy)-zzm4(ix,iy)))
       ssptr(i)=pi0*(rrm3(ix,iy)+rrm4(ix,iy))*dl
15    continue
      if (kmesh .eq. 'dnull') then
       do 16, i=ixpt1(2)+1,ixpt2(2)
        ix=i
        iy=1
        dl=dsqrt((rrm1(ix,iy)-rrm2(ix,iy))*(rrm1(ix,iy)-rrm2(ix,iy))+
     +           (zzm1(ix,iy)-zzm2(ix,iy))*(zzm1(ix,iy)-zzm2(ix,iy)))
        score(i)=pi0*(rrm1(ix,iy)+rrm2(ix,iy))*dl

        iy=iysptrx(1)
        dl=dsqrt((rrm3(ix,iy)-rrm4(ix,iy))*(rrm3(ix,iy)-rrm4(ix,iy))+
     +           (zzm3(ix,iy)-zzm4(ix,iy))*(zzm3(ix,iy)-zzm4(ix,iy)))
        ssptr(i)=pi0*(rrm3(ix,iy)+rrm4(ix,iy))*dl
16     continue
      endif

      do 20, i=1,nw
       call dstscll(i,ix,iy,tp1,tp2,tp3)
       call psi_tdmesh_int(twall(ix,iy))
       call psi_tdmesh_vall
       call psi_tdmesh_rocp
       ro=rotd
       wes0(ix,iy)=fluxc(ix,iy)/dstcflx
       ndst(ix,iy)=3.d0*fluxc(ix,iy)*mol1nrm*md/(pi4*ravr3*1.d-12*ro)
20    continue

      end
c-----------------------------------------------------------------------

      subroutine dstscll(k, ix,iy,tp1,tp2,tp3)
      implicit none
c
c --- produce number and type of a boundary cell
c
c      Use (dustdim)
      Use (dustcom)

      integer k
      integer ix,iy,tp1,tp2,tp3

      ix=ixcl(k)
      iy=iycl(k)
      tp1=tp1cl(ix,iy)
      tp2=tp2cl(ix,iy)
      tp3=tp3cl(ix,iy)

      end
c-----------------------------------------------------------------------

      subroutine dstflxsp(ix,iy,tp1,tp2,tp3,cflux,pflux,aflux)
c
c --- calculates spatterd fluxes from walls
c
      implicit none

c      Use (dustdim)
      Use (dustcom)
      Use (dustinp)
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
     ,            (yres_mx,yres_mx_s),(dse_mx,dse_mx_s)
      equivalence (epsp_mx,epsp_mx_s)
	       
      integer ix,iy,tp1,tp2,tp3
      integer ker
      real*8  cflux,pflux,aflux
      real*8  eis,eas,fish,etai
      real*8  flxa,flxp,flxt,flxc,vta
      real*8  amu,bmag,br,bz,bt,brm,bzm,btm,vnp,vna
      real*8  yresrc,ycsprc
      real*8  psiyresfc,psiycsptf
      integer idstcgp

      ker=idstcgp(ix,iy)
      if(ker .gt. 0) then
       write(*,110)'***Error: geometry of boundary cell, code',ker
110    format(a35,I6)
       goto 101
      endif
      if ( typecl(ix,iy) .ne. (tp1*10+tp2)*10+tp3 ) then
       print *, '***Error: incorrect type of boundary cell',ix,iy
       goto 101
      endif

      brm=.25d0*(brm1(ix,iy)+brm2(ix,iy)+brm3(ix,iy)+brm4(ix,iy))
      bzm=.25d0*(bzm1(ix,iy)+bzm2(ix,iy)+bzm3(ix,iy)+bzm4(ix,iy))
      btm=.25d0*(btm1(ix,iy)+btm2(ix,iy)+btm3(ix,iy)+btm4(ix,iy))
      bmag=brm*brm+bzm*bzm+btm*btm
      bmag=dsqrt(bmag)
      br=brm/bmag
      bz=bzm/bmag
      bt=btm/bmag

c    - 1/4 atom thermal velocity,cm/s
      vta=0.25*dsqrt((8.d0/pi0)*(tam(ix,iy)*ev2erg)/(mol1*1.d-24*mp))

      if (tp1 .eq. 2) then
       vnp=dabs(vym(ix,iy))
       vna=vta
      elseif (tp1 .eq. 3) then
       vnp=dabs(vym(ix,iy))
       vna=vta
      elseif (tp1 .eq. 1) then
       if (((tp2 .eq. 1) .and. (tp3 .eq. 1)) .or.
     *     ((tp2 .eq. 2) .and. (tp3 .eq. 2))) then
        amu=br*r13+bz*z13
        vnp=dabs(vpm(ix,iy)*amu)
        vna=vta
       elseif (((tp2 .eq. 1) .and. (tp3 .eq. 2)) .or.
     *         ((tp2 .eq. 2) .and. (tp3 .eq. 1))) then
        amu=br*r42+bz*z42
        vnp=dabs(vpm(ix,iy)*amu)
        vna=vta
       endif
      endif

c    - plasma particle flux on surface
      flxp=vnp*npm(ix,iy)      ! cm-2*s-1
c    - flux of neutral atoms
      flxa=nam(ix,iy)*vna      !atom flux cm-2*s-1 when [na]=cm^-3
      flxt=(flxp+flxa)*1.d-20
      flxc=0.d0
c    - reduction coefficient for RES and chemical sruttering yeilds
      yresrc=psiyresfc(flxt)
      ycsprc=psiycsptf(flxt,twall(ix,iy))

      fish=-.5d0*dlog(2*pi0*me/mp*(1.d0+tim(ix,iy)/tem(ix,iy)))	         ! sheath potential drop, e*fi/Te
      etai=2.d0+0.5d0*(1.d0+tem(ix,iy)/tim(ix,iy))
      eis=etai*tim(ix,iy)+fish*tem(ix,iy)

      call psi_tmesh_int(.5d0*eis)
      call psi_tmesh_yt
      call psi_tmesh_ytt

c    - associated with chemical sputtering
      flxc=flxc+fycsp(ip)*flxp*ycsprc
c    - associated with RES
      flxc=flxc+fyres(ip)*flxp*yresrc
c    - associated with physical sputtering
      flxc=flxc+fypsp(ip)*flxp

c    -- atom projectiles
      eas=2.d0*tam(ix,iy) ! everaged energy of atom striking surface

      call psi_tmesh_int(.5d0*eas)
      call psi_tmesh_yt
      call psi_tmesh_ytt

c    - associated with chemical sputtering
      flxc=flxc+fycsp(ip)*flxa*ycsprc
c    - associated with RES
      flxc=flxc+fyres(ip)*flxa*yresrc
c    - associated with physical sputtering
      flxc=flxc+fypsp(ip)*flxa
      
      cflux=flxc
      pflux=flxp
      aflux=flxa

      goto 102
101   cflux=0.d0
      pflux=0.d0
      aflux=0.d0
102   end
