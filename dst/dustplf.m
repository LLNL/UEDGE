      integer function idstplf(namf)
      implicit none
      character*256 namf
c
c---this function reads uedge plasma data from ascii file namf
c   and stores data in the common block dustgm.
c   nam is a character array of len
c
c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)

      character*7 a,b
      integer i,j,k,ix,iy
      integer knx,kny,kixpt1(1:2),kixpt2(1:2),kiysptrx(1:2),kixtop
      integer kixrb1,kixlb2
      integer ker
      real*8  mtocm,m3tocm3
      real*8  neg,teg,exg,eyg,npg,tig,vpg,vyg,
     ,        nag,tag,vapg,vayg,nmg,tmg,vpmg,vymg,
     ,        erg,ezg,pxg,pyg,prg,pzg,vrg,vzg,vcg

      data m3tocm3 /1.0d6/
      data mtocm   /1.0d2/
      
      real*8 tfact
      data   tfact /1.0d0/

      ker=0
      open(7,file=namf)

      read(7,109) a
109   format(a8)
      read(7,108) a
108   format(a5)
      if (a .ne. kmesh) then
       print *, '***Error: plasma data is inconsistent with the mesh typ
     *e!'
c       stop
       return
      endif
      if (kmesh .eq. 'snull') then
       read(7,1021) a,a,a,a,a,a
1021   format(6(a7))
       read (7,1011) knx,kny,kixpt1(1),kixpt2(1),kiysptrx(1),kixtop
1011   format(6I7)
      elseif (kmesh .eq. 'dnull') then
       read(7,1022) a,a,a,a,a,a,a,a,a,a
1022   format(10(a7))
       read (7,1012) knx,kny,kixpt1(1),kixpt2(1),kixpt1(2),kixpt2(2),
     *              kiysptrx(1),kiysptrx(2),kixrb1,kixlb2
1012   format(10I7)
      else
       print *, '***Error: unknown plasma geometry',kmesh
      endif

      k=0
      if(knx .lt.   1) k=k+1
      if(knx .gt. nxm) k=k+1
      if(kny .lt.   1) k=k+1
      if(kny .gt. nym) k=k+1
      if(k .gt. 0) then
       ker=6
       write (*,103) 'idstplf:',kmesh,knx,kny
103    format(a8,a5,2(I7))
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
        write (*,1061) 'idstplf1',nx,ny,
     *                 ixpt1(1),ixpt2(1),iysptrx(1),ixtop
        write (*,1061) 'idstplf2',knx,kny,
     *                 kixpt1(1),kixpt2(1),kiysptrx(1),kixtop
1061    format(a8,6(I7))
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
        write (*,1062) 'idstplf1',nx,ny,
     *                 ixpt1(1),ixpt2(1),ixpt2(1),ixpt2(2),
     *                 iysptrx(1),iysptrx(2),ixrb1,ixlb2
        write (*,1062) 'idstplf2',knx,kny,
     *                 kixpt1(1),kixpt2(1),kixpt2(1),kixpt2(2),
     *                 kiysptrx(1),kiysptrx(2),kixrb1,kixlb2
1062    format(a8,10(I7))
        goto 1003
       endif
      endif

      read (7,104) a,a,b,b,b,b,b,b,b,b
104   format(10(a7))
      read (7,107) a,a,b,b,b,b,b,b
      read (7,107) a,a,b,b,b,b,b,b
107   format(8(a7))

      do  1, ix=1,nx
       do 11, iy=1,ny
        read (7,100) i,j,neg,teg,exg
        read (7,105) eyg,npg,tig
        read (7,105) vpg,vyg,nag
        read (7,105) tag,vapg,vayg
        read (7,105) nmg,tmg,vpmg
        read (7,105) vymg,erg,ezg
        read (7,105) pxg,pyg,prg
        read (7,105) pzg,vrg,vzg
100     format(2(I7),3(d24.15))
105     format(3(d24.15))

        k=0
        if(i .eq. ix) k=k+1
        if(j .eq. iy) k=k+1
        if(k .ne. 2) goto 1001

         nem(ix,iy)=neg/m3tocm3
         tem(ix,iy)=teg*tfact
         exm(ix,iy)=exg
         eym(ix,iy)=eyg
         npm(ix,iy)=npg/m3tocm3
         tim(ix,iy)=tig*tfact
         vpm(ix,iy)=vpg*dsqrt(tfact)*mtocm
         vym(ix,iy)=vyg*dsqrt(tfact)*mtocm
         nam(ix,iy)=nag/m3tocm3
         tam(ix,iy)=tag*tfact
        vpam(ix,iy)=vapg*dsqrt(tfact)*mtocm
        vyam(ix,iy)=vayg*dsqrt(tfact)*mtocm
         nmm(ix,iy)=nmg/m3tocm3
         tmm(ix,iy)=tmg*tfact
        vpmm(ix,iy)=vpmg*dsqrt(tfact)*mtocm
        vymm(ix,iy)=vymg*dsqrt(tfact)*mtocm
         erm(ix,iy)=erg
         ezm(ix,iy)=ezg
         pxm(ix,iy)=pxg
         pym(ix,iy)=pyg
         prm(ix,iy)=prg
         pzm(ix,iy)=pzg
        vvrm(ix,iy)=vrg*dsqrt(tfact)*mtocm
        vvzm(ix,iy)=vzg*dsqrt(tfact)*mtocm
11     continue
1     continue

      read(7,110) a,a,b
110   format(3(A7))
      do  2, ix=1,nx
       do 22, iy=1,ny
        read (7,111) i,j,vcg
111     format(2(I7),d24.15)
        vcvm(ix,iy)=dabs(vcg)*dsqrt(tfact)*mtocm
22     continue
2     continue

      goto 1003

1001  ker=1
      print *,'err=1 in idstplf:','i=',i,'j=',j

1003  close(7)

      idstplf=ker
      return

      end


      subroutine dstpnul

      implicit none

c      Use (dustdim)
      Use (dustcom)
      Use (dustcgm)

      integer ix,iy

      do  1, ix=1,nxm
       do 11, iy=1,nym
        nem(ix,iy)=0.d0
        tem(ix,iy)=0.d0
        exm(ix,iy)=0.d0
        eym(ix,iy)=0.d0
        npm(ix,iy)=0.d0
        tim(ix,iy)=0.d0
        vpm(ix,iy)=0.d0
        vym(ix,iy)=0.d0
        nam(ix,iy)=0.d0
        tam(ix,iy)=0.d0
       vpam(ix,iy)=0.d0
       vyam(ix,iy)=0.d0
        nmm(ix,iy)=0.d0
        tmm(ix,iy)=0.d0
       vpmm(ix,iy)=0.d0
       vymm(ix,iy)=0.d0
        erm(ix,iy)=0.d0
        ezm(ix,iy)=0.d0
        pxm(ix,iy)=0.d0
        pym(ix,iy)=0.d0
        prm(ix,iy)=0.d0
        pzm(ix,iy)=0.d0
       vvrm(ix,iy)=0.d0
       vvzm(ix,iy)=0.d0
       vcvm(ix,iy)=0.d0
 
        crm(ix,iy)=0.d0
        czm(ix,iy)=0.d0
        ctm(ix,iy)=0.d0
        frm(ix,iy)=0.d0
        fzm(ix,iy)=0.d0
        ftm(ix,iy)=0.d0
       cram(ix,iy)=0.d0
       czam(ix,iy)=0.d0
       ctam(ix,iy)=0.d0
       crmm(ix,iy)=0.d0
       czmm(ix,iy)=0.d0
       ctmm(ix,iy)=0.d0
       pprm(ix,iy)=0.d0
       ppzm(ix,iy)=0.d0
       pptm(ix,iy)=0.d0
       
       twall(ix,iy)=0.d0
11     continue
1     continue

      end
