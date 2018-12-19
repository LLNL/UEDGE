c
c	'3D modeling of dust dymanics in tokamak 
c	  edge plasmas in real geometry'
c	   San-Diego, version 8.0
c
      subroutine dustt

      implicit none

c      Use (dustdim)
      Use (dustinp)
      Use (dustcom)
      Use (dustcgm)
      Use (dustcntrl)
      Use (dustout)

      integer idstplz,idstimp,idstgfi
      integer ker

      ker=0
c---prepare constants
      call std_cnst
      
c---prepare random number generators
      call std_rnd1f
      call std_rnd2f
      call std_rnd3f
      call std_rnd4f
      call std_rnd5f
      call std_rnd6f
      call std_rnd7f
      call std_rnd8f

c---assign values to input parameters
      call dstinp
      call dstpnul
      call dstinul
      
c---prepare PSI parameters
      call dstprr

c---read geometry file into common block
      ker=idstgfi(fnamgmt) 
      if(ker .ne. 0) then
      write(*,105)'gmf error=',ker,fnamgmt
105   format(a7,I6,1x,a56)
      goto 101
      endif

c---read plasma data file into common block
      ker=idstplz(fnamplz) 
      if(ker .ne. 0) then
      write(*,106)'plz error=',ker,fnamplz
      goto 101
      endif

c---read impurity and sources data file into common block
      ker=idstimp(fnamimp) 
      if(ker .ne. 0) then
       write(*,106)'imp error=',ker,fnamimp
       goto 101
      endif

106   format(a7,I6,1x,a56)

c---zeroing the statistic/output arrays
      call dststn
      call dstrji

c---init output file
      call dstrjifo(ker)
      if(ker .gt. 0) then
      write(*,107)'ofile error=',ker,namfo
107   format(a7,I6,1x,a56)
      goto 101
      endif

c---dust source initialization for stat mode
      if (jstt .gt. 0) then
       call dstsrcini
      endif

c---generate M-C process
      call dsttrj

c---output data
      if (jstt .gt. 0) then
       call dststo
      endif
c---close output file
      call dstrjifc

101   print *,'done'
c      stop
      return

      end
      