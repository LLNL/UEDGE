      integer function idstcgp(i,j)
c --------------------------------------------
c   this function prepares
c   the geometry parameters of a cell (i,j)
c   and stores them in common dstcgm.
c   It returns  >0 if error occurs
c --------------------------------------------
      implicit none

c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)

      integer i,j
      integer k
      real*8  s,s1,s2

      k=0
      icgx=i
      icgy=j

      if(i .lt.  1) k=k+1
      if(i .gt. nx) k=k+1
      if(j .lt.  1) k=k+1
      if(j .gt. ny) k=k+1
      if(k .gt.  0) goto 77

      rc1=rrm1(i,j)
      zc1=zzm1(i,j)
      rc2=rrm2(i,j)
      zc2=zzm2(i,j)
      rc3=rrm3(i,j)
      zc3=zzm3(i,j)
      rc4=rrm4(i,j)
      zc4=zzm4(i,j)

c---normal eq. coefs for 1-2
      r12=zc2-zc1
      z12=rc1-rc2
      d12=-(rc1*r12+zc1*z12)
      s=dsqrt(r12*r12+z12*z12)
      r12=r12/s
      z12=z12/s
      d12=d12/s
      s1=r12*rc3+z12*zc3+d12
      s2=r12*rc4+z12*zc4+d12
      s=s1*s2
      if (s .le. 0.d0) then
       print *, 'Cell ix=',i,' iy=',j,' is not convex'
       k=k+1
      endif
      if(s1 .lt. 0.d0) then
       r12=-r12
       z12=-z12
       d12=-d12
      endif

c---normal eq. coefs for 1-3
      r13=zc3-zc1
      z13=rc1-rc3
      d13=-(rc1*r13+zc1*z13)
      s=dsqrt(r13*r13+z13*z13)
      r13=r13/s
      z13=z13/s
      d13=d13/s
      s1=r13*rc2+z13*zc2+d13
      s2=r13*rc4+z13*zc4+d13
      s=s1*s2
      if (s .le. 0.d0) then
       print *, 'Cell ix=',i,' iy=',j,' is not convex'
       k=k+1
      endif
      if(s1 .lt. 0.d0) then
       r13=-r13
       z13=-z13
       d13=-d13
      endif

c---normal eq. coefs for 3-4
      r34=zc4-zc3
      z34=rc3-rc4
      d34=-(rc3*r34+zc3*z34)
      s=dsqrt(r34*r34+z34*z34)
      r34=r34/s
      z34=z34/s
      d34=d34/s
      s1=r34*rc1+z34*zc1+d34
      s2=r34*rc2+z34*zc2+d34
      s=s1*s2
      if (s .le. 0.d0) then
       print *, 'Cell ix=',i,' iy=',j,' is not convex'
       k=k+1
      endif
      if(s1 .lt. 0.d0) then
       r34=-r34
       z34=-z34
       d34=-d34
      endif

c---normal eq. coefs for 4-2
      r42=zc2-zc4
      z42=rc4-rc2
      d42=-(rc4*r42+zc4*z42)
      s=dsqrt(r42*r42+z42*z42)
      r42=r42/s
      z42=z42/s
      d42=d42/s
      s1=r42*rc1+z42*zc1+d42
      s2=r42*rc3+z42*zc3+d42
      s=s1*s2
      if (s .le. 0.d0) then
       print *, 'Cell ix=',i,' iy=',j,' is not convex'
       k=k+1
      endif
      if(s1 .lt. 0.d0) then
       r42=-r42
       z42=-z42
       d42=-d42
      endif

77    idstcgp=k
      return

      end


      integer function idstinc(r,z)
c ----------------------------------------------------------------------
c   this function determines to which cell point (r,z) belongs.
c   It returns 0 if (r,z) is out of modelling domain;
c   or returns the cell code as iy*nxm+ix; 
c   or -1 if the point is on the boundary
c ----------------------------------------------------------------------
      implicit none

c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)

      integer idstcgp,idstwrz
      integer i,j,k
      integer lkod,lcod
      integer knum(8)
      real*8  r,z

      k=0
      do 1 i=1,nx
       do 11 j=1,ny
        lkod=idstcgp(i,j)
        if(lkod .gt. 0) then
         lkod=-2
         goto 77
        endif
        lcod=idstwrz(r,z)
        if(lcod .gt. 0) then
         k=k+1
         if(k .le. 8) then
          knum(k)=j*nxm+i
         endif
        endif
11     continue
1     continue
      
      if(k .eq. 0) then
       lkod=0
      elseif(k .eq. 1) then
       lkod=knum(k)
      else
       lkod=-1
      endif

77    idstinc=lkod
      return

      end


      integer function idstwrz(r,z)
c ----------------------------------------------------------------------
c   if (r,z) is inside then all distances should be of the same sign
c ----------------------------------------------------------------------
      implicit none

c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)

      integer k,kk
      real*8  r,z
      real*8  s

      k=0
      kk=0

      s=r12*r+z12*z+d12
      if(s .ge. 0.d0) then
       k=k+1
      endif
      if(s .le. 0.d0) then
       kk=kk+1
      endif

      s=r13*r+z13*z+d13
      if(s .ge. 0.d0) then
       k=k+1
      endif
      if(s .le. 0.d0) then
       kk=kk+1
      endif

      s=r34*r+z34*z+d34
      if(s .ge. 0.d0) then
       k=k+1
      endif
      if(s .le. 0.d0) then
       kk=kk+1
      endif

      s=r42*r+z42*z+d42
      if(s .ge. 0.d0) then
       k=k+1
      endif
      if(s .le. 0.d0) then
       kk=kk+1
      endif

      if ((k .eq. 4) .or. (kk .eq. 4)) then
       idstwrz=1
      else
       idstwrz=0
      endif
      
      return

      end



      integer function idstpin(rt1,zt1,rt2,zt2)
c ----------------------------------------------------------------------
c   the particle moves straightforward in cell from point p1=(rt1,zt1)
c   to point p2=(rt2,zt2). Point p1 must be inside the cell.
c   The current cell parameters must be pre-calculated by idstcgp.
c   The function checks if p2 is inside (idstpin=0), or outside (idstpin>0).
c   In the case when p2 is outside, it calcualtes the distance to the
c   nearest boundary and assignes idstpin the index of this boundary:
c   idstpin=1 (1-2);=2 (1-3);=3 (3-4);=4 (4-2)
c   The common block dstcgm contains values of minimal distance ppt and
c   coordinates of intersection point (rrt,zzt). This common block
c   receives also the unity normal vector (vnr,vnz)
c   to the boudary surface (directed into the current cell).
c ----------------------------------------------------------------------
      implicit none

c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)

      integer k
      real*8  rt1,zt1,rt2,zt2
      real*8  seps,sepsm,pmx
      real*8  s,ss,p,pp,cr,cz,pf
      real*8  drrt,dzzt,dspl

      data seps  / 0.d-99 /
      data sepsm /-0.d-99 /
      data pmx   / 1.d+111/
      data dspl  / 1.d-4  /

      k=0
      cr=rt2-rt1
      cz=zt2-zt1
      pf=dsqrt(cr*cr+cz*cz)
      cr=cr/pf
      cz=cz/pf
      pp=pmx

c---intersection with 1-2
      s=r12*rt2+z12*zt2+d12
      if(s .lt. 0.d0) then
       if(s .gt. sepsm) then
        p=pf
       else
        ss=r12*cr+z12*cz
        p=-(r12*rt1+z12*zt1+d12)/ss
       endif
       if(p .le. pp) then
        k=1
        pp=p
        vnr=r12
        vnz=z12
        drrt=dspl*r12
        dzzt=dspl*z12
       endif
      else
       if(s .lt. seps) then
        p=pf
        if(p .le. pp) then
         k=1
         pp=p
         vnr=r12
         vnz=z12
         drrt=dspl*r12
         dzzt=dspl*z12
        endif
       endif
      endif

c---intersection with 1-3
      s=r13*rt2+z13*zt2+d13
      if(s .lt. 0.d0) then
       if(s .gt. sepsm) then
        p=pf
       else
        ss=r13*cr+z13*cz
        p=-(r13*rt1+z13*zt1+d13)/ss
       endif
       if(p .le. pp) then
        k=2
        pp=p
        vnr=r13
        vnz=z13
        drrt=dspl*r13
        dzzt=dspl*z13
       endif
      else
       if(s .lt. seps) then
        p=pf
        if(p .le. pp) then
         k=2
         pp=p
         vnr=r13
         vnz=z13
         drrt=dspl*r13
         dzzt=dspl*z13
        endif
       endif
      endif

c---intersection with 3-4
      s=r34*rt2+z34*zt2+d34
      if(s .lt. 0.d0) then
       if(s .gt. sepsm) then
       p=pf
       else
       ss=r34*cr+z34*cz
       p=-(r34*rt1+z34*zt1+d34)/ss
       endif
       if(p .le. pp) then
        k=3
        pp=p
        vnr=r34
        vnz=z34
        drrt=dspl*r34
        dzzt=dspl*z34
       endif
      else
       if(s .lt. seps) then
       p=pf
        if(p .le. pp) then
         k=3
         pp=p
         vnr=r34
         vnz=z34
         drrt=dspl*r34
         dzzt=dspl*z34
        endif
       endif
      endif

c---intersection with 4-2
      s=r42*rt2+z42*zt2+d42
      if(s .lt. 0.d0) then
       if(s .gt. sepsm) then
        p=pf
       else
        ss=r42*cr+z42*cz
        p=-(r42*rt1+z42*zt1+d42)/ss
       endif
       if(p .le. pp) then
        k=4
        pp=p
        vnr=r42
        vnz=z42
        drrt=dspl*r42
        dzzt=dspl*z42
       endif
      else
       if(s .lt. seps) then
        p=pf
        if(p .le. pp) then
         k=4
         pp=p
         vnr=r42
         vnz=z42
         drrt=dspl*r42
         dzzt=dspl*z42
        endif
       endif
      endif

      if(k .gt. 0) then
       ppt=pp
       rrt=rt1+cr*pp
       zzt=zt1+cz*pp
       rrtin=rrt+drrt
       zztin=zzt+dzzt
       rrtout=rrt-drrt
       zztout=zzt-dzzt
       vnt=0.d0
      endif
 
      idstpin=k
      return

      end


      integer function idstcrsprp()
      implicit none
c
c   Prepares arrays of neighbor cells and transitions
c
c     single null:
c      kcrs=  0 - plasma transition
c             1 - wall collision
c             2 - enter core
c             3 - cross separatrix toward core
c             4 - PF wall collision
c             5 - inner plate collision
c             6 - outer plate collision
c             7 - PF in->out transition
c             8 - PF out->in transition
c            -1 - cell is out of mesh
c     double null:
c      kcrs=  0 - plasma transition
c            11 - innerwall collision
c            12 - outer wall collision
c             2 - enter core
c             3 - cross first separatrix toward core
c            41 - bottom PF wall collision
c            42 - top PF wall collision
c            51 - inner bottom plate collision
c            51 - inner top plate collision
c            61 - outer bottom plate collision
c            62 - outer top plate collision
c            71 - PF in->out bottom transition
c            72 - PF & SOL in->out top transition
c            81 - PF out->in bottom transition
c            82 - PF & SOL out->in top transition
c            -1 - cell is out of mesh
      
c      Use (dustdim)
      Use (dustcom)
      
      integer i,j,k
      integer ker
      
      ker=0
      
      do 10, i=1,nxm
       do 20, j=1,nym
        do 30, k=1,4
         ixcrs(i,j,k)=-1
         iycrs(i,j,k)=-1
         kcrs (i,j,k)=-1
30      continue
20     continue
10    continue

      if (kmesh .eq. 'snull') then
       do 11, i=1,nx
        do 21, j=1,ny
         if (i .eq. 1) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=4
           kcrs (i,j,2)=5
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=5
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=5
           kcrs (i,j,3)=1
           kcrs (i,j,4)=0
          endif
         elseif (i .lt. ixpt1(1)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=4
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=1
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt1(1)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt2(1)+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=4
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=7
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt2(1)+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=7
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=1
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt1(1)+1) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt2(1)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=2
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt2(1)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. iysptrx(1)+1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=3
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=1
           kcrs (i,j,4)=0
          endif
         elseif (i .lt. ixpt2(1)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=2
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. iysptrx(1)+1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=3
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=1
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt2(1)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt1(1)+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=2
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt1(1)+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. iysptrx(1)+1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=3
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=1
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt2(1)+1) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt1(1)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=4
           kcrs (i,j,2)=8
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt1(1)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=8
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=1
           kcrs (i,j,4)=0
          endif
         elseif (i .lt. nx) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=4
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=1
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. nx) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=4
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=6
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=6
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=1
           kcrs (i,j,4)=6
          endif
         else
          print *, '***Error: incorrect mesh cell ',i,j
          goto 1001
         endif
21      continue
11     continue
      elseif (kmesh .eq. 'dnull') then
       do 12, i=1,nx
        do 22, j=1,ny
         if (i .eq. 1) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=41
           kcrs (i,j,2)=51
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=51
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=51
           kcrs (i,j,3)=11
           kcrs (i,j,4)=0
          endif
         elseif (i .lt. ixpt1(1)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=41
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=11
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt1(1)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt2(2)+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=41
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=71
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt2(2)+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=71
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=11
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt1(1)+1) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt2(2)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=2
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt2(2)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. iysptrx(1)+1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=3
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=11
           kcrs (i,j,4)=0
          endif
         elseif (i .lt. ixpt2(1)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=2
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. iysptrx(1)+1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=3
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=11
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt2(1)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt1(2)+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=2
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt1(2)+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. iysptrx(1)+1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt1(2)+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=3
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=72
          elseif (j .le. iysptrx(2)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt1(2)+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=72
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=11
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt2(1)+1) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt1(2)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=42
           kcrs (i,j,2)=72
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(2)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt1(2)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=72
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=11
           kcrs (i,j,4)=0
          endif
         elseif (i .lt. ixrb1) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=42
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=11
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixrb1) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=42
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=52
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=52
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=11
           kcrs (i,j,4)=52
          endif
         elseif (i .le. ixlb2) then
          continue
         elseif (i .eq. ixlb2+1) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=42
           kcrs (i,j,2)=62
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=62
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=62
           kcrs (i,j,3)=12
           kcrs (i,j,4)=0
          endif
         elseif (i .lt. ixpt1(2)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=42
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=12
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt1(2)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt2(1)+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=42
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=82
          elseif (j .le. iysptrx(2)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt2(1)+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=82
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=12
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt1(2)+1) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt2(1)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=2
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt2(1)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. iysptrx(1)+1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt2(1)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=3
           kcrs (i,j,2)=82
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(2)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt2(1)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=82
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=12
           kcrs (i,j,4)=0
          endif
         elseif (i .lt. ixpt2(2)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=2
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. iysptrx(1)+1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=3
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=12
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt2(2)) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt1(1)+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=2
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=ixpt1(1)+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. iysptrx(1)+1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=3
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=12
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. ixpt2(2)+1) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt1(1)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=41
           kcrs (i,j,2)=81
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .le. iysptrx(1)) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=ixpt1(1)
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=81
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=12
           kcrs (i,j,4)=0
          endif
         elseif (i .lt. nx) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=41
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=0
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i+1
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=12
           kcrs (i,j,4)=0
          endif
         elseif (i .eq. nx) then
          if (j .eq. 1) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i
           iycrs(i,j,1)=j
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=41
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=61
          elseif (j .lt. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j+1
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=0
           kcrs (i,j,4)=61
          elseif (j .eq. ny) then
           ixcrs(i,j,1)=i
           ixcrs(i,j,2)=i-1
           ixcrs(i,j,3)=i
           ixcrs(i,j,4)=i
           iycrs(i,j,1)=j-1
           iycrs(i,j,2)=j
           iycrs(i,j,3)=j
           iycrs(i,j,4)=j
           kcrs (i,j,1)=0
           kcrs (i,j,2)=0
           kcrs (i,j,3)=12
           kcrs (i,j,4)=61
          endif
         else
          print *, '***Error: incorrect mesh cell ',i,j
          goto 1001
         endif
22      continue
12     continue
      else
       print *, '***Error: unknown mesh geometry ',kmesh
       goto 1001
      endif
      goto 1002

1001  ker=1

1002  idstcrsprp=ker

      return
      end
      