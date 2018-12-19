      subroutine dstrjifo(ierr)
      implicit none
      integer ierr

      Use (dustinp)

      ierr=0

      if(jtrj .gt. 0) then

      open(3,file=namfo,err=1)
      write(3,*)'r z th vr vz vt rd td md w t fp fm qp qm v lp lt Qd xli
     *qud Pth'
      endif
      goto 2

1     ierr=1
2     return
      end


      subroutine dstrjifc

      implicit none

      Use (dustinp)
      Use (dusttrj)


      if(jtrj .gt. 0) then

      write (3,*) 'npt=',kjm
      close(3)

      endif

      end



      subroutine dstrjif1

      implicit none

      Use (dustinp)
      Use (dusttrj)
      Use (std_cns)

      real*8 vsum

      if(jtrj .gt. 0) then
       vsum=dsqrt(vrb*vrb+vzb*vzb+vtb*vtb)
       kjms=-1
       write(3,1020) 1.d-2*rb,1.d-2*zb,thb*pi2u,vrb,vzb,vtb,rdb,tdb,
     *               msdb,wtb,tb,fpb,fmb,qpb,qmb,vsum,spb,stb,fidb,
     *               xliqudb,lghtb
1020  format(1p21e25.15e3)
      endif

      end



      subroutine dstrjifr(k)

      implicit none

      integer k,i

      Use (dustinp)
      Use (dusttrj)
      Use (std_cns)

      real*8 vsum

      if(jtrj .gt. 0) then

       if(jrartrj .lt. 1) then
        kjms=kjms+1
         if(kjms .ge. 1) then
         kjm=kjm+1
         if(kjm .le. jtrjfm) then
          vsum=dsqrt(vvr*vvr+vvz*vvz+vvt*vvt)
          write(3,1030) 1.d-2*rr,1.d-2*zz,thth*pi2u,vvr,vvz,vvt,rdd,tdd,
     *                  msdd,wtt,tt,fpdd,fmdd,qpdd,qmdd,vsum,spp,stt,
     *                  fidd,xliqudd,lghtd
1030      format(1p21e25.15e3)
         else
          kjm=jtrjfm
         endif
         kjms=0
        endif
       endif

       if(jrartrj .ge. 1) then
        do 33 i=1,k
         kjms=kjms+1
         if(kjms .ge. jrartrj) then
          kjm=kjm+1
          if(kjm .le. jtrjfm) then
           vsum=dsqrt(vra(i)*vra(i)+vza(i)*vza(i)+vta(i)*vta(i))
           write(3,1030) 1.d-2*ra(i),1.d-2*za(i),tha(i)*pi2u,vra(i),
     *                   vza(i),vta(i),rda(i),tda(i),msda(i),wta(i),
     *                   tta(i),fpa(i),fma(i),qpa(i),qma(i),
     *                   vsum,spa(i),sta(i),fida(i),xliquda(i),lghta(i)
          else
           kjm=jtrjfm
          endif
          kjms=0
         endif
33      continue
       endif

      endif

      end
