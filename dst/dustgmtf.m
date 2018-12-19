      integer function idstgfi(namf)
c ----------------------------------------------------------------------
c   this function reads uedge geometry data from ascii file nam
c   and stores data in the common block dustgm.
c   nam is a character array of len
c ----------------------------------------------------------------------
      implicit none

c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)
      Use (std_cns)

      integer i,j,k,n,ix,iy
      integer ker,ierr
      integer tpclg,knw,knwt
      integer idstcgp,idstcrsprp
      real*8 r1,z1,t1,r2,z2,t2,r3,z3,t3,r4,z4,t4
      real*8 rv12,zv12,rv13,zv13,rv42,zv42,rv43,zv43
      real*8 rv14,zv14,rv23,zv23
      real*8 s,ravr,sz1,sz2,sz3,sz4,sz5,sz6
      real*8 mtocm
      character*256 namf
      character*7 a,b
      character*1 c

      data mtocm /100.0d0/

      ker=0
      open(3,file=namf)

      read(3,107) a
107   format(a8)
      read(3,108) a
108   format(a5)
      kmesh=a
      if (kmesh .eq. 'snull') then
       read(3,1021) a,a,a,a,a,a
1021   format(6(a7))
       read (3,1011) nx,ny,ixpt1(1),ixpt2(1),iysptrx(1),ixtop
1011   format(6I7)
      elseif (kmesh .eq. 'dnull') then
       read(3,1022) a,a,a,a,a,a,a,a,a,a
1022   format(10(a7))
       read (3,1012) nx,ny,ixpt1(1),ixpt2(1),ixpt1(2),ixpt2(2),
     *              iysptrx(1),iysptrx(2),ixrb1,ixlb2
1012   format(10I7)
      else
       write (*,*) '***Error: unknown mesh geometry ',kmesh
c       stop
       return
      endif

      k=0
      if(nx .lt.   1) k=k+1
      if(nx .gt. nxm) k=k+1
      if(ny .lt.   1) k=k+1
      if(ny .gt. nym) k=k+1
      if(k  .gt.   0) then
       ker=6
       write (*,*) '***Error: mesh size is out of limits: nx=',nx,
     *             ' ny=',ny
       goto 1003
      endif

      read (3,104) a,a,b,b,b,b,b,b,b,b
104   format(2(a7),8(a7))

      do  1 ix=1,nx
       do 11 iy=1,ny
        read (3,100) i,j,r1,z1,r2,z2
        read (3,105) r3,z3,r4,z4
100     format(2(I7),4(d24.15))
105     format(4(d24.15))
        k=0
        if(i .eq. ix) k=k+1
        if(j .eq. iy) k=k+1
        if(k .ne. 2) goto 1001
        rrm1(ix,iy)=r1*mtocm
        zzm1(ix,iy)=z1*mtocm
        rrm2(ix,iy)=r2*mtocm
        zzm2(ix,iy)=z2*mtocm
        rrm3(ix,iy)=r3*mtocm
        zzm3(ix,iy)=z3*mtocm
        rrm4(ix,iy)=r4*mtocm
        zzm4(ix,iy)=z4*mtocm
11     continue
1     continue
      goto 1002

1001  ker=1
      print *,'err=1 in idsgfi:','i=',i,'j=',j
      goto 1003

1002  read (3,106) a,a,b,b,b,b,b,b,b,b,b,b,b,b
106   format(2(a7),12(a7))

      do  2 ix=1,nx
       do 22 iy=1,ny
        read (3,200) i,j,r1,z1,t1,r2
        read (3,201) z2,t2,r3,z3
        read (3,201) t3,r4,z4,t4
200     format(2(I7),4(d24.15))
201     format(4(d24.15))
        k=0
        if(i .eq. ix) k=k+1
        if(j .eq. iy) k=k+1
        if(k .ne. 2) goto 1004
        brm1(ix,iy)=r1
        bzm1(ix,iy)=z1
        btm1(ix,iy)=t1
        brm2(ix,iy)=r2
        bzm2(ix,iy)=z2
        btm2(ix,iy)=t2
        brm3(ix,iy)=r3
        bzm3(ix,iy)=z3
        btm3(ix,iy)=t3
        brm4(ix,iy)=r4
        bzm4(ix,iy)=z4
        btm4(ix,iy)=t4
22     continue
2     continue

      do 113, ix=1,nxm
       do 112, iy=1,nym
        typecl(ix,iy)=0
        tp1cl (ix,iy)=0
        tp2cl (ix,iy)=0
        tp3cl (ix,iy)=0
112    continue
113   continue
      do 114, k=1,ncellm
       ixcl(k)=0
       iycl(k)=0
114   continue

      read(3,119) c
      read(3,119) c
      read(3,119) c
      read(3,119) c
      read(3,119) c
      read(3,119) c
119   format (A1)
      read(3,120) knw
120   format (I7)
      
      if (kmesh .eq. 'snull') then
       knwt=2*(nx+ny)+ixpt1(1)-ixpt2(1)-8
      elseif (kmesh .eq. 'dnull') then
       knwt=2*nx+4*ny+ixpt1(1)+ixpt1(2)-ixpt2(1)-ixpt2(2)-20
      else
       print *, '***Error: unknown mesh geometry ',kmesh
      endif
      
      if (knw .ne. knwt) then
       ker=8
       print *, '***Error: inconsistent number of boundary cells',knw
       goto 1003
      endif
      nw=knw
      read(3,119) c
      do 14, k=1,nw
       read (3,118) tpclg,i,j
118    format (3(I7))
       ixcl(k)=i
       iycl(k)=j
       typecl(i,j)=tpclg
       n=tpclg/100
       tp1cl(i,j)=n
       if ((n .lt. 1) .or. (n .gt. 3)) goto 1009
       tpclg=tpclg-tp1cl(i,j)*100
       n=tpclg/10
       tp2cl(i,j)=n
       if ((n .lt. 1) .or. (n .gt. 2)) goto 1009
       tpclg=tpclg-tp2cl(i,j)*10
       n=tpclg
       tp3cl(i,j)=n
       if ((n .lt. 1) .or. (n .gt. 2)) goto 1009
c     ---check the indexes
       if (kmesh .eq. 'snull') then
        if (tp1cl(i,j) .eq. 1) then
         if ((iycl(k) .le. 1) .or. (iycl(k) .ge. ny)) goto 1009
         if (tp2cl(i,j) .eq. 1) then
          if (ixcl(k) .ne. 1) goto 1009
         elseif (tp2cl(i,j) .eq. 2) then
          if (ixcl(k) .ne. nx) goto 1009
         endif
        elseif (tp1cl(i,j) .eq. 2) then
         if (iycl(k) .ne. ny) goto 1009
         if ((ixcl(k) .le. 1) .or. (ixcl(k) .ge. nx)) goto 1009
        elseif (tp1cl(i,j) .eq. 3) then
         if (iycl(k) .ne. 1) goto 1009
         if ((ixcl(k) .le. 1) .or. (ixcl(k) .ge. nx)) goto 1009
         if ((ixcl(k) .gt. ixpt1(1)) .and.
     *       (ixcl(k) .le. ixpt2(1))) goto 1009
        endif
       elseif (kmesh .eq. 'dnull') then
        if (tp1cl(i,j) .eq. 1) then
         if ((iycl(k) .le. 1) .or. (iycl(k) .ge. ny)) goto 1009
         if (tp2cl(i,j) .eq. 1) then
          if (tp3cl(i,j) .eq. 1) then
           if (ixcl(k) .ne. 1) goto 1009
          elseif(tp3cl(i,j) .eq. 2) then
           if (ixcl(k) .ne. ixrb1) goto 1009
          endif
         elseif (tp2cl(i,j) .eq. 2) then
          if (tp3cl(i,j) .eq. 1) then
           if (ixcl(k) .ne. nx) goto 1009
          elseif(tp3cl(i,j) .eq. 2) then
           if (ixcl(k) .ne. ixlb2+1) goto 1009
          endif
         endif
        elseif (tp1cl(i,j) .eq. 2) then
         if (iycl(k) .ne. ny) goto 1009
         if (tp2cl(i,j) .eq. 1) then
          if ((ixcl(k) .le. 1) .or.
     *        (ixcl(k) .ge. ixrb1)) goto 1009
         elseif (tp2cl(i,j) .eq. 2) then
          if ((ixcl(k) .le. ixlb2+1) .or.
     *        (ixcl(k) .ge. nx)) goto 1009
         endif
        elseif (tp1cl(i,j) .eq. 3) then
         if (iycl(k) .ne. 1) goto 1009
         if (tp2cl(i,j) .eq. 1) then
          if (tp3cl(i,j) .eq. 1) then
           if ((ixcl(k) .le. 1) .or.
     *         (ixcl(k) .gt. ixpt1(1))) goto 1009
          elseif (tp3cl(i,j) .eq. 2) then
           if ((ixcl(k) .le. ixpt2(1)) .or.
     *         (ixcl(k) .ge. ixrb1)) goto 1009
          endif
         elseif (tp2cl(i,j) .eq. 2) then
          if (tp3cl(i,j) .eq. 1) then
           if ((ixcl(k) .le. ixpt2(2)) .or.
     *         (ixcl(k) .ge. nx)) goto 1009
          elseif (tp3cl(i,j) .eq. 2) then
           if ((ixcl(k) .le. ixlb2+1) .or.
     *         (ixcl(k) .gt. ixpt1(2))) goto 1009
          endif
         endif
        endif
       endif
14    continue
      goto 1010
      
1009  ker=9
      write (*,140) '***Error: boundary cell type is incorrect',
     *              i,j,typecl(i,j)
140   format (A,3(I7))
      goto 1003
      
1010  do 3, ix=1,nx
       do 33, iy=1,ny
        rv12=rrm2(ix,iy)-rrm1(ix,iy)
        zv12=zzm2(ix,iy)-zzm1(ix,iy)
        rv13=rrm3(ix,iy)-rrm1(ix,iy)
        zv13=zzm3(ix,iy)-zzm1(ix,iy)
        rv42=rrm2(ix,iy)-rrm4(ix,iy)
        zv42=zzm2(ix,iy)-zzm4(ix,iy)
        rv43=rrm3(ix,iy)-rrm4(ix,iy)
        zv43=zzm3(ix,iy)-zzm4(ix,iy)
        rv14=rrm4(ix,iy)-rrm1(ix,iy)
        zv14=zzm4(ix,iy)-zzm1(ix,iy)
        rv23=rrm3(ix,iy)-rrm2(ix,iy)
        zv23=zzm3(ix,iy)-zzm2(ix,iy)
        
        s=.5d0*(dabs(rv12*zv13-rv13*zv12)+dabs(rv42*zv43-rv43*zv42))
        ravr=.25d0*(rrm1(ix,iy)+rrm2(ix,iy)+rrm3(ix,iy)+rrm4(ix,iy))
        vol(ix,iy)=2.d0*pi0*ravr*s
      
        sz1=dsqrt(rv12*rv12+zv12*zv12)
        sz2=dsqrt(rv13*rv13+zv13*zv13)
        sz3=dsqrt(rv42*rv42+zv42*zv42)
        sz4=dsqrt(rv43*rv43+zv43*zv43)
        sz5=dsqrt(rv14*rv14+zv14*zv14)
        sz6=dsqrt(rv23*rv23+zv23*zv23)
        sz1=min(sz1,sz4)
        sz2=min(sz2,sz3)
        sz1=min(sz1,sz2)
        sz5=min(sz5,sz6)
        msize(ix,iy)=0.1d0*min(sz1,sz5)
c --- check UEDGE mesh
        ierr=idstcgp(ix,iy)
        if (ierr .gt. 0) then
         ker=3
        endif
33     continue
3     continue

c --- prepare mesh transition arrays
      ierr=idstcrsprp()
      if (ierr .gt. 0) then
       ker=4
      endif
      goto 1003
      

1004  ker=2
      write (*,*) 'err=2 in idsgfi:','i=',i,'j=',j

1003  close(3)

      idstgfi=ker
      return

      end
