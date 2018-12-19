c-----------------------------------------------------------------------
c     trajectory mode
c-----------------------------------------------------------------------
      integer function idstsrs(ix,iy,r,z,vr,vz,vt)
      implicit none

      real*8 r,z,vr,vz,vt
      integer ix,iy

c      Use (dustdim)
      Use (dustinp)
      Use (dustcom)
      Use (dustcgm)

      integer std_irnd1,idstcgp
      integer n,ker,kcer,kplt
      integer tp1,tp2,tp3
      integer cellnum
      integer flg
      integer idstinc
      real*8  g,gg
      real*8  peps,delta

      data peps  /1.d-8/
      data delta /1.d-6/

      ker=0

c---make offset in random number sequence
      call std_rnd1sk(irndsk)

c---indices of source cell
      if (inipnt .lt. 0) then
10     cellnum=idstinc(rini,zini)
       iy=cellnum/nxm
       ix=cellnum-iy*nxm
       flg=0
       if (cellnum .lt. 0) then
        rini=rini+delta
        goto 10
       endif
       if (kmesh .eq. 'dnull') then
        if ((ix .gt. ixrb1) .and. (ix .le. ixlb2)) flg=1
       endif
       if ((cellnum .eq. 0) .or. (flg .eq. 1)) then
        ker=12
        print *, '***Error: initial particle position is out of mesh!',
     *           rini,zini
        goto 101
       endif
      elseif (inipnt .eq. 0) then
c---random choise of a boundary cell
       n=nw-1
       n=std_irnd1(n)+1
       ix=ixcl(n)
       iy=iycl(n)
       
      elseif (inipnt .gt. 0) then
c--- specified input for cell
       ix=inipnt/100
       iy=inipnt-ix*100
c--- check input cell number
       kcer=0
       if (ix .lt.  1) kcer=kcer+1
       if (ix .gt. nx) kcer=kcer+1
       if (iy .lt.  1) kcer=kcer+1
       if (iy .gt. ny) kcer=kcer+1
       if (kmesh .eq. 'dnull') then
        if ((ix .gt. ixrb1) .and. (ix .le. ixlb2)) then
         kcer=kcer+1
        endif
       endif

       if(kcer .gt. 0) then
        ker=10
        goto 101
       endif
      endif
      
      ker=idstcgp(ix,iy)         ! cell parameters
      if (ker .gt. 0) goto 101   ! geometry error

      if (inipnt .lt. 0) then
       vnt=0.d0
       vnr=0.d0
       vnz=1.d0
       r=rini
       z=zini
      else
       if (inipnt .eq. 0) then
c --- random point inside the cell
        call std_rndd1x(g)
        g=g*(1.d0-peps)
        gg=1.d0-g
       elseif (inipnt .gt. 0) then
c --- fixed point at the cell center
        g=0.5d0
        gg=1.d0-g
       endif
c --- type of cell
       tp1=tp1cl(ix,iy)
       tp2=tp2cl(ix,iy)
       tp3=tp3cl(ix,iy)
c --- 
       if (tp1 .eq. 1) then         ! plate
        if (tp2 .eq. 1) then        ! inner
         if (tp3 .eq. 1) then       ! bottom
          kplt=1
         elseif (tp3 .eq. 2) then   ! top
          kplt=2
         else                       ! unknown
          kplt=-1
         endif
        elseif (tp2 .eq. 2) then    ! outer
         if (tp3 .eq. 1) then       ! bottom
          kplt=2
         elseif (tp3 .eq. 2) then   ! top
          kplt=1
         else                       ! unknown
          kplt=-1
         endif
        endif
       elseif (tp1 .eq. 2) then     ! wall
        kplt=4
       elseif (tp1 .eq. 3) then     ! PF
        kplt=3
       elseif (tp1 .eq. 0) then     ! other
        kplt=0
       else                         ! unknown
        kplt=-1
       endif

c ---coordinates of birth point and normal
       vnt=0.d0
       if (kplt .eq. 1) then
        vnr=r13
        vnz=z13
        r=rc1*g+rc3*gg
        z=zc1*g+zc3*gg
       elseif (kplt .eq. 2) then
        vnr=r42
        vnz=z42
        r=rc2*g+rc4*gg
        z=zc2*g+zc4*gg
       elseif (kplt .eq. 3) then
        vnr=r12
        vnz=z12
        r=rc1*g+rc2*gg
        z=zc1*g+zc2*gg
       elseif (kplt .eq. 4) then
        vnr=r34
        vnz=z34
        r=rc3*g+rc4*gg
        z=zc3*g+zc4*gg
       elseif (kplt. eq. 0) then
        vnr=0.d0
        vnz=1.d0
        r=.25d0*(rc1+rc2+rc3+rc4)
        z=.25d0*(zc1+zc2+zc3+zc4)
       else
        ker=11
        goto 101
       endif

       r=r+vnr*peps
       z=z+vnz*peps
      endif
c
c---flight direction
c
      if (inipnt .lt. 0) then
       vp0=dsqrt(vrini*vrini+vzini*vzini+vtini*vtini)
       if (vp0 .gt. 0.d0) then
        vr=vrini/vp0
        vz=vzini/vp0
        vt=vtini/vp0
       else
        vr=0.d0
        vz=0.d0
        vt=0.d0
       endif
      else
       if(inivel .lt. 1) then
c ---random direction
        call rfldif(vr,vz,vt, vnr,vnz,vnt)
       else
c ---input angles vs normal
        call dirfix(vr,vz,vt, unorm0,uphi0, vnr,vnz,vnt)
       endif
      endif
      goto 102
      
101   print *, '***Error: incorrect dust start cell ',ix,iy

102   idstsrs=ker
      return
      end

c-----------------------------------------------------------------------
c     statistical regime
c-----------------------------------------------------------------------
      integer function idstsrst(m,ix,iy,r,z,vr,vz,vt)
      implicit none

      real*8 r,z,vr,vz,vt
      integer m,ix,iy

c      Use (dustdim)
      Use (dustinp)
      Use (dustcom)
      Use (dustcgm)
      Use (std_cns)

      integer std_irnd1,idstcgp
      integer tp1,tp2,tp3
      integer ker,kplt
      real*8  g,gg
      real*8  peps
      real*8  vrd,vzd,vtd,vrf,vzf,vtf
      real*8  vnrmr,vnrmz,vnrmt
      real*8  unorm,uphi,vabs,amu

      data peps /1.d-8/

      ker=0
      
c --- make offset in random number sequence
      call std_rnd1sk(irndsk)

c --- incremental points at the cell
      gg=(m-.5d0)/nispcl(ix,iy)
      g=1.d0-gg
      
      ker=idstcgp(ix,iy)       ! cell parameters
      if(ker .gt. 0) goto 101  ! geometry error

      tp1=tp1cl(ix,iy)
      tp2=tp2cl(ix,iy)
      tp3=tp3cl(ix,iy)
      
      if (tp1 .eq. 1) then         ! plate
       if (tp2 .eq. 1) then        ! inner
        if (tp3 .eq. 1) then       ! bottom
         kplt=1
        elseif (tp3 .eq. 2) then   ! top
         kplt=2
        else                       ! unknown
         kplt=-1
        endif
       elseif (tp2 .eq. 2) then    ! outer
        if (tp3 .eq. 1) then       ! bottom
         kplt=2
        elseif (tp3 .eq. 2) then   ! top
         kplt=1
        else                       ! unknown
         kplt=-1
        endif
       endif
      elseif (tp1 .eq. 2) then     ! wall
       kplt=4
      elseif (tp1 .eq. 3) then     ! PF
       kplt=3
      elseif (tp1 .eq. 0) then     ! other
       kplt=0
      else                         ! unknown
       kplt=-1
      endif

c --- coordinates of birth point and normal
      vnt=0.d0
      if (kplt .eq. 1) then
       vnr=r13
       vnz=z13
       r=rc1*g+rc3*gg
       z=zc1*g+zc3*gg
      elseif (kplt .eq. 2) then
       vnr=r42
       vnz=z42
       r=rc2*g+rc4*gg
       z=zc2*g+zc4*gg
      elseif (kplt .eq. 3) then
       vnr=r12
       vnz=z12
       r=rc1*g+rc2*gg
       z=zc1*g+zc2*gg
      elseif (kplt .eq. 4) then
       vnr=r34
       vnz=z34
       r=rc3*g+rc4*gg
       z=zc3*g+zc4*gg
      else
       ker=10
       goto 101
      endif

      r=r+vnr*peps
      z=z+vnz*peps

c---flight direction
      if (inivel .eq. 0) then
       call rfldif(vr,vz,vt, vnr,vnz,vnt)
      else
       unorm=unorm1*u2pi
       uphi=uphi1*u2pi
       vnrmr=dsin(unorm)*dsin(uphi)
       vnrmt=dsin(unorm)*dcos(uphi)
       vnrmz=dcos(unorm)
       call dirfix(vrf,vzf,vtf, unorm0,uphi0, vnrmr,vnrmz,vnrmt)
       call rfldif(vrd,vzd,vtd, vnrmr,vnrmz,vnrmt)
       vr=vrd*(1.d0-dirvc)+vrf*dirvc
       vz=vzd*(1.d0-dirvc)+vzf*dirvc
       vt=vtd*(1.d0-dirvc)+vtf*dirvc
       vabs=dsqrt(vr*vr+vz*vz+vt*vt)
       if (vabs .gt. 0.d0) then
        vr=vr/vabs
        vz=vz/vabs
        vt=vt/vabs
       else
        vr=0.d0
        vz=0.d0
        vt=0.d0
       endif
c    - check outward velocity direction
       amu=vr*vnr+vz*vnz+vt*vnt
       if (amu .lt. 0.d0) then
        call rfldif(vr,vz,vt, vnr,vnz,vnt)
       endif
      endif
      goto 102
      
101   print *, '***Error: incorrect source cell ',ix,iy

102   idstsrst=ker
      return
      end
