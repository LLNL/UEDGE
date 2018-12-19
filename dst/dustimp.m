      integer function idstimp(namf)
      implicit none

c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)

      character*256 namf
      integer idstimf    ! read impurity and surface file
      integer ker

      ker=idstimf(namf)
      if(ker .gt. 0) goto 1

      if (jimp .gt. 0) call dstiml    ! find velocity components of impurities

1     idstimp=ker

      return
      end


      subroutine dstiml
c ----------------------------------------------------------------------
      implicit none

c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)

      integer ix,iy,l,m
      real*8  vps,vys,vcr,vcz,vct

      do 6, l=1,nimp
       do  5 ix=1,nx
        do 55 iy=1,ny
         vps=vpiam(ix,iy,l)
         vys=vyiam(ix,iy,l)
         call dstvcp(ix,iy,vps,vys,vcr,vcz,vct)
         criam(ix,iy,l)=vcr
         cziam(ix,iy,l)=vcz
         ctiam(ix,iy,l)=vct
55      continue
5      continue
6     continue
      do 7, m=1,nimpz
       do  9 ix=1,nx
        do 99 iy=1,ny
         vps=vpizm(ix,iy,m)
         vys=vyizm(ix,iy,m)
         call dstvcp(ix,iy,vps,vys,vcr,vcz,vct)
         crizm(ix,iy,m)=vcr
         czizm(ix,iy,m)=vcz
         ctizm(ix,iy,m)=vct
99      continue
9      continue
7     continue
      end


      integer function idstimf(namf)
c ----------------------------------------------------------------------
c   this function reads uedge impurity and fluxes data from ascii file namf
c   and stores data in the common block dustimp and dustflx.
c   nam is a character array of len
c ----------------------------------------------------------------------
      implicit none
 
c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)
      Use (std_cns)
      Use (psitab)
c      Use (psitab_s)
      Use (psitabint)

      character*256 namf
      character*7 a,b
      integer i,j,k,l,m,n,ix,iy
      integer knx,kny,kixpt1(1:2),kixpt2(1:2),kiysptrx(1:2),kixtop
      integer knw,knwt
      integer kixrb1,kixlb2
      integer knz,knza,kiz,kzz,kznuc
      integer kia,kza
      integer in,izn
      integer ker
      real*8  kmz,kzpot,kmia
      real*8  mtocm,m2tocm2,m3tocm3
      real*8  ng,tg,vpg,vyg
      real*8  flg,fng,flmg,fnmg,sbndg

      data m3tocm3 /1.d6/
      data m2tocm2 /1.d4/
      data mtocm   /1.d2/
      
      real*8 tfact
      data tfact /1.0d0/

      ker=0
      if (jsrc .gt. 0) then
       open(3,file=namf)
       read(3,101) a,a
101    format(2(a7))
       read (3,102) knz,knza
       nimpz=knz
       nimp=knza
102    format(2I7)
       read(3,103) a,a,a
103    format(3(a7))
       do 18, l=1,nimp
        read(3,19) kia,kza,kmia
        ia    (l)=kia
        zatm  (l)=kza
        mia   (l)=kmia
        mianrm(l)=mia(l)*mol1nrm
        k=0
        do 30, i=1,Nprj_pt
         if (zatm(l) .eq. IZprjctl(i)) k=i
30      continue
        if (k .gt. 0) then
         iia(l)=k
        else
         print *, 'No PSI data for impurity neutrals number ',l
c         stop
         return
        endif
        if (mia(l) .ne. AZprjctl(k)) then
         print *, '*** Warning: mass of impurity atom ',l,' different fr
     *om PSI data'
        endif
18     continue
19     format (2I7,D24.15)
       read(3,130) a,a,a,a,a
130    format(5(a7))
       do 20, m=1,nimpz
        read(3,21) kiz,kzz,kznuc,kmz,kzpot
        iz   (m)=kiz
        zimp (m)=kzz
        znuc (m)=kznuc
        mz   (m)=kmz
        mznrm(m)=mz(m)*mol1nrm
        zpot (m)=kzpot
        k=0
        do 40, i=1,Nprj_pt
         if (znuc(m) .eq. IZprjctl(i)) k=i
40      continue
        if (k .gt. 0) then
         iiz(m)=k
        else
         print *, "No PSI data for impurity ion number ",m
c         stop
         return
        endif
        if (mz(m) .ne. AZprjctl(k)) then
         print *, "*** Warning: mass of impurity ion ",m," different fro
     *m PSI data"
        endif
20     continue
21     format(3I7,2D24.15)
       read(3,124) a
124    format(a8)
       read(3,123) a
123    format(a5)
       if (a .ne. kmesh) then
        print *, '***Error: impurity data is inconsistent with the mesh 
     *type!'
c        stop
        return
       endif
       if (kmesh .eq. 'snull') then
        read(3,1041) a,a,a,a,a,a
1041    format(6(a7))
        read (3,1051) knx,kny,kixpt1(1),kixpt2(1),kiysptrx(1),kixtop
1051    format(6I7)
       elseif (kmesh .eq. 'dnull') then
        read(3,1042) a,a,a,a,a,a,a,a,a,a
1042    format(10(a7))
        read (3,1052) knx,kny,kixpt1(1),kixpt2(1),kixpt1(2),kixpt2(2),
     *                kiysptrx(1),kiysptrx(2),kixrb1,kixlb2
1052    format(10I7)
       else
        print *, 'Error: unknown impurity geometry',kmesh
       endif
       read(3,106) a
106    format(a7)
       read (3,107) knw
       nw=knw
107    format(I7)

       k=0
       if(knx .lt.   1) k=k+1
       if(knx .gt. nxm) k=k+1
       if(kny .lt.   1) k=k+1
       if(kny .gt. nym) k=k+1
       if(k .gt. 0) then
        ker=6
        write (*,109) 'idstimf:',knx,kny
109     format(a8,6(I7))
        goto 1003
       endif

       k=0
       if (kmesh .eq. 'snull') then
        if(knx        .ne.        nx) k=k+1
        if(kny        .ne.        ny) k=k+1
        if(kixpt1(1)  .ne.  ixpt1(1)) k=k+1
        if(kixpt2(1)  .ne.  ixpt2(1)) k=k+1
        if(kiysptrx(1).ne.iysptrx(1)) k=k+1
        if(kixtop     .ne.     ixtop) k=k+1
        if(k .gt. 0) then
         ker=7
         write (*,1081) 'idstplf1',nx,ny,
     *                  ixpt1(1),ixpt2(1),iysptrx(1),ixtop
         write (*,1081) 'idstplf2',knx,kny,
     *                  kixpt1(1),kixpt2(1),kiysptrx(1),kixtop
1081     format(a8,6(I7))
         goto 1003
        endif
       elseif (kmesh .eq. 'dnull') then
        if(knx        .ne.        nx) k=k+1
        if(kny        .ne.        ny) k=k+1
        if(kixpt1(1)  .ne.  ixpt1(1)) k=k+1
        if(kixpt1(2)  .ne.  ixpt1(2)) k=k+1
        if(kixpt2(1)  .ne.  ixpt2(1)) k=k+1
        if(kixpt2(2)  .ne.  ixpt2(2)) k=k+1
        if(kiysptrx(1).ne.iysptrx(1)) k=k+1
        if(kiysptrx(2).ne.iysptrx(2)) k=k+1
        if(kixrb1     .ne.     ixrb1) k=k+1
        if(kixlb2     .ne.     ixlb2) k=k+1
        if(k .gt. 0) then
         ker=7
         write (*,1082) 'idstplf1',nx,ny,
     *                  ixpt1(1),ixpt2(1),ixpt2(1),ixpt2(2),
     *                  iysptrx(1),iysptrx(2),ixrb1,ixlb2
         write (*,1082) 'idstplf2',knx,kny,
     *                  kixpt1(1),kixpt2(1),kixpt2(1),kixpt2(2),
     *                  kiysptrx(1),kiysptrx(2),kixrb1,kixlb2
1082     format(a8,10(I7))
         goto 1003
        endif
       endif
        
       if (kmesh .eq. 'snull') then
        knwt=2*(nx+ny)+ixpt1(1)-ixpt2(1)-8
       elseif (kmesh .eq. 'dnull') then
        knwt=2*nx+4*ny+ixpt1(1)+ixpt1(2)-ixpt2(1)-ixpt2(2)-20
       endif
       if(nw .ne. knwt) then
        ker=8
        write (*,110) '***Error: number of boundary cells is incorrect i
     *n *.uedgm:',nw,knwt
110     format (A8,2(I7))
        goto 1003
       endif
        
       read (3,111) a,a
111    format(2A1)
       read (3,112) a,a,b,b,b
112    format(5(A7))

       do 2, k=1,nw
        read (3,113) i,j,flg,fng,sbndg
c        call dstidxcnv(i,j)
        flp (i,j)=flg
        fnp (i,j)=fng/m3tocm3
        sbnd(i,j)=sbndg*m2tocm2
2      continue
113    format (2(I7),3(D24.15))

       read (3,111) a,a
       read (3,114) a,a,b,b,b,b
114    format(6(a7))

       do 3, k=1,nw
        read (3,115) i,j,flg,fng,flmg,fnmg
c        call dstidxcnv(i,j)
        fla(i,j)=flg
        fna(i,j)=fng/m3tocm3
        flm(i,j)=flmg
        fnm(i,j)=fnmg/m3tocm3
3      continue
115    format (2(I7),4(D24.15))

       do 4, l=1,nimp
        read (3,111) a,a
        read (3,116) in
        read (3,117) a,a,b,b
        do 5, k=1,nw
         read (3,118) i,j,flg,fng
c         call dstidxcnv(i,j)
         flia(i,j,in)=flg
         fnia(i,j,in)=fng/m3tocm3
5       continue
4      continue
       do 6, m=1,nimpz
        read (3,111) a,a
        read (3,116) izn
        read (3,117) a,a,b,b
        do 7, k=1,nw
         read (3,118) i,j,flg,fng
c         call dstidxcnv(i,j)
         fliz(i,j,izn)=flg
         fniz(i,j,izn)=fng/m3tocm3
7       continue
6      continue
116    format(I7)
117    format(4(a7))
118    format(2(I7),2(D24.15))
      endif

      if (jimp .le. 0) goto 1003
      do 8, l=1,nimp
       read (3,119) a
       read (3,116) in
       read (3,120) a,a,b,b,b
       do  9 ix=1,nx
        do 10 iy=1,ny
         read (3,100) i,j,ng,tg,vyg
         k=0
         if(i .eq. ix) k=k+1
         if(j .eq. iy) k=k+1
         if(k .ne. 2) goto 1001
         niam(ix,iy,in)=ng/m3tocm3
         tiam(ix,iy,in)=tg*tfact
         vyiam(ix,iy,in)=vyg*dsqrt(tfact)*mtocm
10      continue
9      continue
8     continue
      
      do 11, m=1,nimpz
       read (3,119) a
       read (3,116) izn
       read (3,121) a,a,b,b,b,b
       do 12 ix=1,nx
        do 13 iy=1,ny
         read (3,122) i,j,ng,tg,vpg,vyg
         k=0
         if(i .eq. ix) k=k+1
         if(j .eq. iy) k=k+1
         if(k .ne. 2) goto 1001
         nizm(ix,iy,izn)=ng/m3tocm3
         tizm(ix,iy,izn)=tg*tfact
         vpizm(ix,iy,izn)=vpg*dsqrt(tfact)*mtocm
         vyizm(ix,iy,izn)=vyg*dsqrt(tfact)*mtocm
13      continue
12     continue
11    continue
119   format(a7)
120   format(5(a7))
121   format(6(a7))
122   format(2(I7),4(d24.15))
100   format(2(I7),3(d24.15))

      goto 1003

1001  ker=1
      print *,'err=1 in idstimf:','i=',i,'j=',j

1003  if (jsrc .gt. 0) close(3)

      idstimf=ker
      return

      end


      subroutine dstinul
      implicit none
      
c      Use (dustdim)
      Use (dustcom)
      
      integer ix,iy
      integer k,l,m

      do 2, l=1,nimpm
       do 4, ix=1,nxm
        do 5, iy=1,nym
         flia (ix,iy,l)=0.d0
         fnia (ix,iy,l)=0.d0
         niam (ix,iy,l)=0.d0
         tiam (ix,iy,l)=0.d0
         vpiam(ix,iy,l)=0.d0
         vyiam(ix,iy,l)=0.d0
         criam(ix,iy,l)=0.d0
         cziam(ix,iy,l)=0.d0
         ctiam(ix,iy,l)=0.d0
5       continue
4      continue
2     continue

      do 6, m=1,nimpzm
       do 8, ix=1,nxm
        do 9, iy=1,nym
         fliz (ix,iy,m)=0.d0
         fniz (ix,iy,m)=0.d0
         nizm (ix,iy,m)=0.d0
         tizm (ix,iy,m)=0.d0
         vpizm(ix,iy,m)=0.d0
         vyizm(ix,iy,m)=0.d0
         crizm(ix,iy,m)=0.d0
         czizm(ix,iy,m)=0.d0
         ctizm(ix,iy,m)=0.d0
9       continue
8      continue
6     continue

      do 10, ix=1,nxm
       do 11, iy=1,nym
        flp   (ix,iy)=0.d0
        fla   (ix,iy)=0.d0
        flm   (ix,iy)=0.d0
        fnp   (ix,iy)=0.d0
        fna   (ix,iy)=0.d0
        fnm   (ix,iy)=0.d0
        sbnd  (ix,iy)=0.d0
11     continue
10    continue

      end


      subroutine dstidxcnv(ix,iy)
      implicit none
      
c      Use (dustdim)
      Use (dustcom)
      
      integer ix,iy
      
      if (ix .eq. 1 ) ix=ix+1
      if (ix .eq. nx) ix=ix-1
      if (kmesh .eq. 'snull') then
       if ((iy .eq. 1 ) .and.
     *   ( (ix .le. ixpt1(1)) .or.
     *     (ix .gt. ixpt2(1)) )) iy=iy+1
      elseif (kmesh .eq. 'dnull') then
       if ((iy .eq. 1 ) .and.
     *     ( (ix .le. ixpt1(1)) .or.
     *      ((ix .gt. ixpt2(1)) .and. (ix .le. ixpt1(2))) .or.
     *       (ix .gt. ixpt2(2)) )) iy=iy+1
      endif
      if (iy .eq. ny) iy=iy-1
      
      end
      