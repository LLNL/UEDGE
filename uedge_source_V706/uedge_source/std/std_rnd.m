      subroutine std_rnd1f
c*** generator of random-sequenty      rnd1
      implicit none
      integer iy,iy0,jy
      common /rnd1c/iy
      data iy0/2003/
c----start from the beginning
      iy=iy0
      return
      entry std_rnd1s(jy)
c----save point
       jy=iy
       return
      entry std_rnd1r(jy)
c---restore point
       iy=jy
       return
      end

      subroutine std_rnd2f
c*** generator of random-sequenty      rnd2
      implicit none
      integer iy,iy0,jy
      common /rnd2c/iy
      data iy0/2011/
c----start from the beginning
      iy=iy0
      return
      entry std_rnd2s(jy)
c----save point
       jy=iy
       return
      entry std_rnd2r(jy)
c---restore point
       iy=jy
       return
      end

      subroutine std_rnd3f
c*** generator of random-sequenty      rnd3
      implicit none
      integer iy,iy0,jy
      common /rnd3c/iy
      data iy0/2017/
c----start from the beginning
      iy=iy0
      return
      entry std_rnd3s(jy)
c----save point
       jy=iy
       return
      entry std_rnd3r(jy)
c---restore point
       iy=jy
       return
      end

      subroutine std_rnd4f
c*** generator of random-sequenty      rnd4
      implicit none
      integer iy,iy0,jy
      common /rnd4c/iy
      data iy0/2027/
c----start from the beginning
      iy=iy0
      return
      entry std_rnd4s(jy)
c----save point
       jy=iy
       return
      entry std_rnd4r(jy)
c---restore point
       iy=jy
       return
      end

      subroutine std_rnd5f
c*** generator of random-sequenty      rnd5
      implicit none
      integer iy,iy0,jy
      common /rnd5c/iy
      data iy0/2029/
c----start from the beginning
      iy=iy0
      return
      entry std_rnd5s(jy)
c----save point
       jy=iy
       return
      entry std_rnd5r(jy)
c---restore point
       iy=jy
       return
      end

      subroutine std_rnd6f
c*** generator of random-sequenty      rnd6
      implicit none
      integer iy,iy0,jy
      common /rnd6c/iy
      data iy0/2039/
c----start from the beginning
      iy=iy0
      return
      entry std_rnd6s(jy)
c----save point
       jy=iy
       return
      entry std_rnd6r(jy)
c---restore point
       iy=jy
       return
      end

      subroutine std_rnd7f
c*** generator of random-sequenty      rnd7
      implicit none
      integer iy,iy0,jy
      common /rnd7c/iy
      data iy0/2053/
c----start from the beginning
      iy=iy0
      return
      entry std_rnd7s(jy)
c----save point
       jy=iy
       return
      entry std_rnd7r(jy)
c---restore point
       iy=jy
       return
      end

      subroutine std_rnd8f
c*** generator of random-sequenty      rnd8
      implicit none
      integer iy,iy0,jy
      common /rnd8c/iy
      data iy0/2063/
c----start from the beginning
      iy=iy0
      return
      entry std_rnd8s(jy)
c----save point
       jy=iy
       return
      entry std_rnd8r(jy)
c---restore point
       iy=jy
       return
      end


      subroutine std_rnd1sk(idummy)
      implicit none
      integer idummy
c
c---skip idummy random points
c
      common /rnd1c/iy
      integer iy,i
c
      if(idummy .gt. 0) then
      do 111 i=1,idummy

      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif

111   continue
      endif

      end



      real*8 function std_rnd1(idummy)
      implicit none
      integer idummy
c
c---random point is x=rnd1(idummy)
c
      common /rnd1c/iy
      integer iy
      real*8 halfm
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      std_rnd1=dble(iy)*halfm
      return
c
      end



      real*8 function std_rndd1(idummy)
      implicit none
      integer idummy
c
c---random point is x=rnd1(idummy)
c
      common /rnd1c/iy
      integer iy
      real*8 halfm
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      std_rndd1=dble(iy)*halfm
      return
c
      end



      subroutine std_rnd1x(x)
      implicit none
      real*8 x
c
c---random point is x
c
      common /rnd1c/iy
      integer iy
      real*8 halfm
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
      return
c
      end



      subroutine std_rndd1x(x)
      implicit none
      real*8 x
c
c---random point is x
c
      common /rnd1c/iy
      integer iy
      real*8 halfm
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
      return
c
      end

      subroutine std_rndd2x(x)
      implicit none
      real*8 x
c
c---random point is x
c
      common /rnd2c/iy
      integer iy
      real*8 halfm
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
      return
c
      end

      subroutine std_rndd3x(x)
      implicit none
      real*8 x
c
c---random point is x
c
      common /rnd3c/iy
      integer iy
      real*8 halfm
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
      return
c
      end

      subroutine std_rndd4x(x)
      implicit none
      real*8 x
c
c---random point is x
c
      common /rnd4c/iy
      integer iy
      real*8 halfm
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
      return
c
      end

      subroutine std_rndd5x(x)
      implicit none
      real*8 x
c
c---random point is x
c
      common /rnd5c/iy
      integer iy
      real*8 halfm
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
      return
c
      end

      subroutine std_rndd6x(x)
      implicit none
      real*8 x
c
c---random point is x
c
      common /rnd6c/iy
      integer iy
      real*8 halfm
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
      return
c
      end

      subroutine std_rndd7x(x)
      implicit none
      real*8 x
c
c---random point is x
c
      common /rnd7c/iy
      integer iy
      real*8 halfm
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
      return
c
      end

      subroutine std_rndd8x(x)
      implicit none
      real*8 x
c
c---random point is x
c
      common /rnd8c/iy
      integer iy
      real*8 halfm
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
      return
c
      end



      integer function std_irnd1(n)
      implicit none
      integer n
c
c---random index within n-numbers is ix=irnd1(n)
c
      common /rnd1c/iy
      integer iy
      real*8 halfm,sn,x
c
      data halfm/4.656612873077393d-010/
c
      sn=n
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
c
      std_irnd1=int(x*sn)+1
      return
c
      end



      subroutine std_rnd1gt(g,*)
c
c---random indicator
c
      implicit none

      integer iy
      real*8 g
      real*8 halfm,x

      common /rnd1c/iy
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
c
      if(x.gt.g) return 1
      return
c
      end



      subroutine std_rnd1lt(g,*)
c
c---random indicator
c
      implicit none

      integer*4 iy
      real*8 g
      real*8 halfm,x

      common /rnd1c/iy
c
      data halfm/4.656612873077393d-010/
c
      iy=iy*843314861
      if(iy.gt.1693666955) then
       iy=(iy-1073741824)-1073741824
       iy=iy+453816693
      else
       iy=iy+453816693
       if(iy.lt.0) iy=(iy+1073741824)+1073741824
      endif
c
      x=dble(iy)*halfm
c
      if(x.lt.g) return 1
      return
c
      end
