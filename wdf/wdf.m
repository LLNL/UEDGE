
      subroutine write_degas
      implicit none
Use(Dim)                # nx,ny from com package
Use(Xpoint_indices)     # ixpt1,ixpt2 from com package
Use(Share)              # geometry from com package
Use(Dimwdf)             # nptsvb,nptshb,nptskb,nptsw,npsegxz
Use(Auxw)               # ixpt1b,ixpt2b,ixtop1b,ixtop2b,
                        # nohzsb,novzsb,nosegsxzb
Use(Options)            # iswdfon

c     main driver routine for wdf package

      if (mod((ixpt2(1)-ixpt1(1)),2) .ne. 0) then
         call remark(" ")
         call remark("*** Error:  write_degas procedure only valid when")
         call remark("               ixpt2-ixpt1 is an even number")
         call remark(" ")
         call xerrab("")
      endif

c     allocate arrays for DEGAS

      ixpt1b  = ixpt1(1)
      ixtop1b = ixpt1b  + (ixpt2(1)-ixpt1(1))/2
      ixtop2b = ixtop1b + 1
      ixpt2b  = ixpt2(1)
      if (geometry .eq. "dnbot") then
	 ixtop1b=ixtop1b-1	# we omit guard cells at the midplane
	 ixtop2b=ixtop2b+1	# of up/down symmetric double nulls
      endif

      novzsb  = max(ixtop1b+2,nx-ixtop2b+3)
      nohzsb  = 2*(ny+1)
      nosegsxzb = 2*novzsb+nohzsb+3

      nptsvb = nohzsb + 1
      nptshb = novzsb + 1
      nptsw  = 2*nptshb + nptsvb + 1
      npsegxz = nptsw - 1
      call gallot("Degas1",0)
      call gallot("Degas2",0)
      call remark("***** allocated DEGAS arrays *****")

c     read data from grd and bbb packages
      if (iswdfon==1) then
         call grd2wdf	# write grd-wdf file from grd package
         call readgrd	# read grid data from grd package
         call bbb2wdf	# write bbb-wdf file from bbb package
         call readbbb	# read plasma data from bbb package
      endif
c
      call degasgrid	# generate grid and wall arrays for DEGAS
      call defaultz	# generate kzone arrays for DEGAS
      call ueplasma	# write UEDGE plasma into DEGAS arrays
      call cgsunits	# convert from UEDGE to DEGAS units
      call write_namelist	# write DEGAS namelist input file

      return
      end

c----------------------------------------------------------------------c

      subroutine readgrd
      implicit none
Use(Dimwdf)
Use(Eqdsk)
Use(Linkgrd)
Use(Degas1)
Use(Degas2)
      integer iunit, ios
      external remark, xerrab, gallot, allot
      external rdgrd1, rdgrd2, rdgrd3

c **************** read the output of the grd package ******************

      data iunit /55/
      open (iunit, file='grd-wdf', form='unformatted', iostat=ios,
     &      status='old')
      if (ios .ne. 0) then
         call xerrab("**** grd-wdf file not found")
      endif

c     read dimensioning parameters and allot storage space for arrays --

      read (iunit) idimw,jdimw,nixw,noregsw
      call gallot("Linkgrd",0)
      call rdgrd1 (iunit)

      call rdgrd2 (iunit)

      read(iunit) nlimw
      call allot("xlimw",nlimw)
      call allot("ylimw",nlimw)

      call rdgrd3 (iunit)

      close(iunit)

      return
      end

c----------------------------------------------------------------------c

      subroutine rdgrd1 (iunit)
      implicit none
Use(Dimwdf)
Use(Linkgrd)
      integer iunit

      read (iunit) cmeshxw,cmeshyw,
     &             ilmaxw,ixpointw,jminw,jmaxw,jsptrxw,jaxisw

      return
      end

c----------------------------------------------------------------------c

      subroutine rdgrd2 (iunit)
      implicit none
Use(Dimwdf)
Use(Eqdsk)
      integer iunit

      read (iunit) bcentrw,rcentrw,rmagxw,zmagxw,simagxw,sibdryw,
     &             rgrid1w,xdimw,zdimw

      return
      end

c----------------------------------------------------------------------c

      subroutine rdgrd3 (iunit)
      implicit none
Use(Dimwdf)
Use(Eqdsk)
      integer iunit

      read (iunit) xlimw,ylimw
      read (iunit) eshotw,etimew,rsepsw,zsepsw,
     &              rvsinw,zvsinw,rvsoutw,zvsoutw

      return
      end

c----------------------------------------------------------------------c

      subroutine degasgrid
      implicit none
Use(Dim)           # nxm,nym from com package
Use(RZ_grid_info)  # rm,zm from com package
Use(Linkbbb)       # nxbbb,nybbb from com package
Use(Dimwdf)
Use(Auxw)          # ixpt1b,ixpt2b,ixtop1b,ixtop2b
Use(Linkgrd)
Use(Degas1)
Use(Degas2)
Use(Eqdsk)
      integer ix,iy,i,j,itot,n
c
c     *** copy UEDGE (rm,zm) onto DEGAS (gridx,gridz)
c

c First, the inboard half of the mesh: 
      do iy = 0, nybbb
         j = nybbb - iy + 1
         i = 1
         gridx(j,i,1) = rm(ixtop1b,iy,4) - rgrid1w
         gridz(j,i,1) = zm(ixtop1b,iy,4)
         do ix = ixtop1b, ixpt1b+1, -1
            i = i + 1
            gridx(j,i,1) = rm(ix,iy,3) - rgrid1w
            gridz(j,i,1) = zm(ix,iy,3)
         enddo
         i = i + 1
         gridx(j,i,1) = gridx(j,i-1,1)
         gridz(j,i,1) = gridz(j,i-1,1)
         i = i + 1
         gridx(j,i,1) = rm(ixpt1b,iy,4) - rgrid1w
         gridz(j,i,1) = zm(ixpt1b,iy,4)
         do ix = ixpt1b, 1, -1
            i = i + 1
            gridx(j,i,1) = rm(ix,iy,3) - rgrid1w
            gridz(j,i,1) = zm(ix,iy,3)
         enddo
      enddo
c fill out the arrays by duplicating plate data
      while (i .lt. nptshb)
         i = i + 1
         do j = 1, nybbb+1
            gridx(j,i,1) = gridx(j,i-1,1)
            gridz(j,i,1) = gridz(j,i-1,1)
         enddo
      endwhile
c Next, the outboard half of the mesh: 
      do iy = 0, nybbb
         j = nybbb + iy + 3
         i = 1
         gridx(j,i,1) = rm(ixtop2b,iy,3) - rgrid1w
         gridz(j,i,1) = zm(ixtop2b,iy,3)
         do ix = ixtop2b, ixpt2b
            i = i + 1
            gridx(j,i,1) = rm(ix,iy,4) - rgrid1w
            gridz(j,i,1) = zm(ix,iy,4)
         enddo
         i = i + 1
         gridx(j,i,1) = gridx(j,i-1,1)
         gridz(j,i,1) = gridz(j,i-1,1)
         i = i + 1
         gridx(j,i,1) = rm(ixpt2b+1,iy,3) - rgrid1w
         gridz(j,i,1) = zm(ixpt2b+1,iy,3)
         do ix = ixpt2b+1, nxbbb
            i = i + 1
            gridx(j,i,1) = rm(ix,iy,4) - rgrid1w
            gridz(j,i,1) = zm(ix,iy,4)
         enddo
      enddo
c fill out the arrays by duplicating plate data
      while (i .lt. nptshb)
         i = i + 1
         do j = nybbb+3, 2*nybbb+3
            gridx(j,i,1) = gridx(j,i-1,1)
            gridz(j,i,1) = gridz(j,i-1,1)
         enddo
      endwhile
c Finally, the "magnetic axis" for core and private flux regions
      j = nybbb + 2
      do i = 1, ixtop1b-ixpt1b+2
         gridx(j,i,1) = rmagxw - rgrid1w
         gridz(j,i,1) = zmagxw
      enddo
      do i = ixtop1b-ixpt1b+3, nptshb
         gridx(j,i,1) = 0.5*(gridx(j-1,nptshb,1)+gridx(j+1,nptshb,1))
         gridz(j,i,1) = 0.5*(gridz(j-1,nptshb,1)+gridz(j+1,nptshb,1))
      enddo

c Wall construction:
      nowals = 3

c First, clockwise around the outside of the mesh
      j=1
      itot=0
      do i=1,nptshb
         itot=itot+1
         xwall(itot,1)=gridx(j,i,1)
         zwall(itot,1)=gridz(j,i,1)
      enddo
c
      i=nptshb
      do j=1,nptsvb
         itot=itot+1
         xwall(itot,1)=gridx(j,i,1)
         zwall(itot,1)=gridz(j,i,1)
      enddo
c
      j=nptsvb
      do i=nptshb,1,-1
         itot=itot+1
         xwall(itot,1)=gridx(j,i,1)
         zwall(itot,1)=gridz(j,i,1)
      enddo
c re-connect with starting point
      itot=itot+1
      xwall(itot,1)=xwall(1,1)
      zwall(itot,1)=zwall(1,1)
      nosegsxz(1)=itot-1
c set default values for kwmat and twall -
      do n = 1, nosegsxz(1)
         twall(n,1,1) = 300.		# room temperature
         kwmat(n,1)   = "'c'"		# carbon
      enddo

c Second, clockwise around the private flux
      itot=1
      i=nptshb
      j=nybbb+2                         # start at "magnetic axis" point
      xwall(itot,2)=gridx(j,i,1)
      zwall(itot,2)=gridz(j,i,1)
      j=j-1                             # jump to inboard half of p.f.
      do i=ixtop1b+3,ixtop1b-ixpt1b+3,-1
         itot=itot+1
         xwall(itot,2)=gridx(j,i,1)
         zwall(itot,2)=gridz(j,i,1)
      enddo
      j=j+2                             # jump to outboard half of p.f.
      do i=ixpt2b-ixtop2b+4,nxbbb-ixtop2b+4
         itot=itot+1
         xwall(itot,2)=gridx(j,i,1)
         zwall(itot,2)=gridz(j,i,1)
      enddo
c re-connect with starting point
      itot=itot+1
      xwall(itot,2)=xwall(1,2)
      zwall(itot,2)=zwall(1,2)
      nosegsxz(2)=itot-1
c set default values for kwmat and twall -
      do n = 1, nosegsxz(2)
         twall(n,1,2) = 300.		# room temperature
         kwmat(n,2)   = "'c'"		# carbon
      enddo

c Third, clockwise around the core
      itot=0
      j=nybbb+1
      do i=ixtop1b-ixpt1b+2,1,-1
         itot=itot+1
         xwall(itot,3)=gridx(j,i,1)
         zwall(itot,3)=gridz(j,i,1)
      enddo
      j=j+2
      do i=1,ixtop1b-ixpt1b+2
         itot=itot+1
         xwall(itot,3)=gridx(j,i,1)
         zwall(itot,3)=gridz(j,i,1)
      enddo
c re-connect with starting point
      itot=itot+1
      xwall(itot,3)=xwall(1,3)
      zwall(itot,3)=zwall(1,3)
      nosegsxz(3)=itot-1
c set default values for kwmat and twall -
      do n = 1, nosegsxz(3)
         twall(n,1,3) = 300.		# room temperature
         kwmat(n,3)   = "'c'"		# carbon
      enddo

      return
      end


c----------------------------------------------------------------------c

      subroutine readbbb
      implicit none
Use(Dim)
Use(Linkbbb)
      integer iunit, ios
      external remark, xerrab, gallot, allot

c **************** read the output of the bbb package ******************

      data iunit /55/
      open (iunit, file='bbb-wdf', form='unformatted', iostat=ios,
     &      status='old')
      if (ios .ne. 0) then
         call xerrab("**** bbb-wdf file not found")
      endif

      read (iunit) nxbbb,nybbb,nycorebbb,nysolbbb,nxleg1bbb,nxcore1bbb,
     &            nxleg2bbb,nxcore2bbb
      read (iunit) nibbb,tibbb,nebbb,tebbb,vflowxbbb,vflowybbb,vflowzbbb,
     &             fnixbbb,fngysibbb,fngysobbb
      read (iunit) geometrybbb

      close (iunit)

      return
      end

c----------------------------------------------------------------------c

      subroutine ueplasma
      implicit none
Use(Dim)                # from com package
Use(Xpoint_indices)     # from com package
Use(Linkbbb)
Use(Dimwdf)
Use(Auxw)               # ixpt1b,ixpt2b,ixtop1b,ixtop2b,nosegsxzb,novzsb
Use(Linkgrd)
Use(Degas1)
Use(Degas2)
      integer ix,iy,jh,iv,iseg

c     Convert UEDGE plasma information to DEGAS namelist arrays

      do ix=1,nxbbb
	 do iy=1,nybbb
c translate UEDGE cell index (ix,iy) to DEGAS cell index (iv,jh)
	    if ((1 .le. ix) .and. (ix .le. ixpt1b)) then
	       jh = nybbb + 1 - iy
	       iv = ixtop1b + 3 - ix
	    elseif ((ixpt1b+1 .le. ix) .and. (ix .le. ixtop1b)) then
	       jh = nybbb + 1 - iy
	       iv = ixtop1b + 1 - ix
	    elseif ((ixtop2b .le. ix) .and. (ix .le. ixpt2b)) then
	       jh = nybbb + 2 + iy
	       iv = ix - ixtop2b + 1
	    elseif ((ixpt2b+1 .le. ix) .and. (ix .le. nxbbb)) then
	       jh = nybbb + 2 + iy
	       iv = ix - ixtop2b + 3
	    endif
	    denihvt(jh,iv,1,1) = nibbb(ix,iy)
	    tihvt(jh,iv,1,1)   = tibbb(ix,iy)
	    denehvt(jh,iv,1)   = nebbb(ix,iy)
	    tehvt(jh,iv,1)     = tebbb(ix,iy)
	    vflowx(jh,iv,1)    = vflowxbbb(ix,iy)
	    vflowy(jh,iv,1)    = vflowybbb(ix,iy)
	    vflowz(jh,iv,1)    = vflowzbbb(ix,iy)
         enddo
      enddo

      do iy=1,nybbb
c translate UEDGE plate index (0,iy) to DEGAS wall segment index (iseg)
	 iseg = novzsb + 1 + (nybbb + 1 - iy)
	 currxzt(iseg,1,1,1) = - fnixbbb(0,iy)
c translate UEDGE plate index (nxm,iy) to DEGAS wall segment index (iseg)
	 iseg = novzsb + 1 + nybbb + 2 + (iy)
	 currxzt(iseg,1,1,1) = fnixbbb(nxbbb,iy)
      enddo

c translate UEDGE private flux wall index (ix) to DEGAS private region
c wall segment index (iseg) on wall number 2
      do ix=1,ixpt1b
	 iseg = 1 + ix
	 currxzt(iseg,1,2,1) = - abs(fngysibbb(ix))# NOTE negative for puffing
      enddo
      do ix=ixpt2b+1,nxbbb
	 iseg = 1 + ixpt1b + 1 + (ix-ixpt2b)
	 currxzt(iseg,1,2,1) = - abs(fngysibbb(ix))
      enddo

c translate UEDGE outermost flux wall index (ix) to DEGAS external
c wall segment index (iseg) on wall number 1
      do ix=1,ixtop1b
	 iseg = ixtop1b-ix+1
	 currxzt(iseg,1,1,1) = - abs(fngysobbb(ix))
      enddo
      do ix=ixtop2b,nxbbb
	 iseg = nosegsxzb-ix+ixtop2b-1
	 currxzt(iseg,1,1,1) = - abs(fngysobbb(ix))
      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine defaultz
      integer jb,jj,ib,ii
Use(Dimwdf)
Use(Auxw)               # ixpt1b
Use(Degas1)
Use(Degas2)
Use(Eqdsk)
Use(Options)
Use(Linkgrd)

c.....Initialize some variables with default values.
      endmark = " $"
      fname   = "degas.in"
      idline  = "plasma information from UEDGE"

c.....Define default plasma zone boundaries.
c.....along horizontal direction
      jb=0
      do jj=1,nptsvb,1
         jb = jb +1
         lboun1(jb) = jj
      enddo
      nohbs = jb
      nohzs = nohbs -1
c.....along vertical direction
      ib=0
      do ii=1,nptshb,1
         ib = ib + 1
         lboun2(ib) = ii
      enddo
      novbs = ib
      novzs = novbs - 1

      ivnull = ixtop1b - ixpt1b + 2     # vertical index of non-physical zone
                                        # between core and private flux

c.....set kzone arrays
      call setkz

c.....Set default values for other DEGAS input parameters -
      xlen = xdimw
      zlen = zdimw
      rmajor = rgrid1w + xdimw/2.

      return
      end

c----------------------------------------------------------------------c

      subroutine setkz
      implicit none
Use(Dimwdf)
Use(Degas1)
Use(Degas2)
Use(Options)
      integer i,j,iv,jh,kzvalue

c.....Define kzone arrays for DEGAS
      do jh=1,nohzs
         do j=lboun1(jh),lboun1(jh+1)-1
            kzone1(j,1)=jh
         enddo
      enddo
c
      do i=2,nptshb-1
         do j=1,nptsvb-1
            kzone1(j,i)=kzone1(j,1)
         enddo
      enddo
c
      do iv=1,novzs
         do i=lboun2(iv),lboun2(iv+1)-1
            kzone2(1,i)=iv
         enddo
      enddo
c
      do i=1,nptshb-1
         kzvalue=kzone2(1,i)
         do j=1,nptsvb-1
            kzone2(j,i)=kzone2(1,i)
            if (kzvalue .eq. ivnull) then
               kzone1(j,i)=-1
               kzone2(j,i)=-1
            endif
         enddo
      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine frrate
c.....For neutral atoms,
c.....set fictitious reaction rate = frr01, for denihvt .eq. zero
c.....                             =     0, for denihvt .ne. zero
c.....and similarly for neutral molecules,
c.....set fictitious reaction rate = frr02, for denihvt .eq. zero
c.....                             =     0, for denihvt .ne. zero
      integer j,i
Use(Dimwdf)
Use(Options)
Use(Degas1)
Use(Degas2)

      do j=1,nohzs
         do i=1,novzs
            if (denihvt(j,i,1,1) .eq. 0.) then
               ficrrhvt(j,i,1,1) = frr01
               ficrrhvt(j,i,1,2) = frr02
            else
               ficrrhvt(j,i,1,1) = 0.
               ficrrhvt(j,i,1,2) = 0.
            endif
         enddo
      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine cgsunits
      implicit none
Use(Dimwdf)
Use(Degas1)
Use(Degas2)
      integer i,j,k,l,m,iv,jh
      real ev
      data ev/1.6022e-19/

c     Convert lengths from [m] to [cm] for DEGAS

      do k=1,nptskb
         do j=1,nptsvb
            do i=1,nptshb
               gridx(j,i,k)=100.*gridx(j,i,k)
               gridz(j,i,k)=100.*gridz(j,i,k)
            enddo
         enddo
      enddo

      do m=1,npw
         do l=1,nptsw
            xwall(l,m)=100.*xwall(l,m)
            zwall(l,m)=100.*zwall(l,m)
         enddo
      enddo

      xlen=xlen*100.
      zlen=zlen*100.
      rmajor=rmajor*100.

c convert plasma from UEDGE to DEGAS units :
      do jh=1,nohzs
         do iv=1,novzs
	    denihvt(jh,iv,1,1) = denihvt(jh,iv,1,1) * 1.0e-06
	    tihvt(jh,iv,1,1)   = tihvt(jh,iv,1,1)   / ev
	    denehvt(jh,iv,1)   = denehvt(jh,iv,1)   * 1.0e-06
	    tehvt(jh,iv,1)     = tehvt(jh,iv,1)     / ev
	    vflowx(jh,iv,1)    = vflowx(jh,iv,1)    * 100.
	    vflowy(jh,iv,1)    = vflowy(jh,iv,1)    * 100.
	    vflowz(jh,iv,1)    = vflowz(jh,iv,1)    * 100.
	 enddo
      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine write_namelist
      integer nunit,i,jwall,k,j,iv,jh
Use(Dimwdf)
Use(Options)
Use(Degas1)
Use(Degas2)
      external remark

c     Write output file to be read by DEGAS code; 
c     file name is given by the variable 'fname' in group Io

      data nunit /66/
      open (nunit, file=fname, form='formatted', status='unknown')

      write(nunit,1000)idline
 1000 format(1a80)
c    
c     first degas namelist /input/
c    
      write(nunit,998)
  998 format(1x,"$input")
      write(nunit, 999) mparts,maxerrs,ncount
  999 format(" mparts=  ",i5," maxerrs= ",i5," ncount=  ",i5)
      write(nunit,1001) lmesh,lprofh,lproft,lprofv,
     &                   loutput,lrflct1,lrflct2,
     &                   lrrulet,lsputter,lsymetry,lunit
 1001 format(" lmesh=   ",i5," lprofh=  ",i5," lproft=  ",i5," lprofv=  ",i5/,
     &       " loutput= ",i5," lrflct1= ",i5," lrflct2= ",i5/,
     &       " lrrulet= ",i5," lsputter=",i5," lsymetry=",i5," lunit=   ",i5)
      write(nunit,1002) nexit,nhyd(1),nocols,nohbs,nohzs,norows
 1002 format(" nexit=   ",i5," nhyd=    ",i5," nocols=  ",i5," nohbs=   ",i5/,
     &       " nohzs=   ",i5," norows=  ",i5)
      write(nunit,1003)(nosegsxz(i),i=1,nowals)
 1003 format(" nosegsxz=",5(i5,1x))
      write(nunit,1004) nosplits,novbs,novzs,nowals,nptshb,nptsvb
 1004 format(" nosplits=",i5," novbs=   ",i5," novzs=   ",i5," nowals=  ",i5/,
     &       " nptshb=  ",i5," nptsvb=  ",i5)
      write(nunit,1005) h2frac,plsang,rmajor,scrrmax,xlen,zlen
 1005 format(" h2frac= ",f6.4," plsang= ",f6.2," rmajor= ",f6.2/,
     &       " scrrmax=",f6.4," xlen=   ",f6.2," zlen=   ",f6.2)
      write(nunit,1006) lmaxwell,shethp
 1006 format(" lmaxwell=",i5," shethp= ",f6.2)
      write(nunit,1007)(t0puff(i),i=1,npns)
 1007 format(" t0puff= ",5(f6.3,9x))
c    
      write(nunit,1009)endmark
 1009 format(1a8)
c    
c     second degas namelist /arrays/
c    
      write(nunit,2000)
 2000 format(1x,"$arrays")

c.....recycling current at divertor plates and walls
      do jwall=1,nowals
        write(nunit,2001)jwall
        write(nunit,2002)(currxzt(i,1,jwall,1),i=1,nosegsxz(jwall))
      enddo
 2001 format(1x,"currxzt(1,1,",i1,",1)=")
 2002 format(18x,1p6e10.3)

c.....electron density
      do i=1,novzs
      write(nunit,2010)i,(denehvt(k,i,1),k=1,min0(6,nohzs))
         do j=7,nohzs,6
         write(nunit,2011)(denehvt(k,i,1),k=j,min0(j+5,nohzs))
         enddo
      enddo
 2010 format(" denehvt(1,",i2,",1)=",1p6e10.3)
 2011 format(17x,1p6e10.3)

c.....ion density
      do i=1,novzs
      write(nunit,2020)i,(denihvt(k,i,1,1),k=1,min0(6,nohzs))
         do j=7,nohzs,6
         write(nunit,2021)(denihvt(k,i,1,1),k=j,min0(j+5,nohzs))
         enddo
      enddo
 2020 format(" denihvt(1,",i2,",1,1)=",1p6e10.3)
 2021 format(19x,1p6e10.3)

c.....fictitious reaction rate for neutrals

      if (frr01 .ne. 0.) then
         do i=1,novzs
         write(nunit,2030)i,(ficrrhvt(k,i,1,1),k=1,min0(6,nohzs))
            do j=7,nohzs,6
            write(nunit,2031)(ficrrhvt(k,i,1,1),k=j,min0(j+5,nohzs))
            enddo
         enddo
 2030 format(" ficrrhvt(1,",i2,",1,1)=",1p6e9.2)
 2031 format(20x,1p6e9.2)
      endif

      if (frr02 .ne. 0.) then
         do i=1,novzs
         write(nunit,2032)i,(ficrrhvt(k,i,1,2),k=1,min0(6,nohzs))
            do j=7,nohzs,6
            write(nunit,2033)(ficrrhvt(k,i,1,2),k=j,min0(j+5,nohzs))
            enddo
         enddo
 2032 format(" ficrrhvt(1,",i2,",1,2)=",1p6e9.2)
 2033 format(20x,1p6e9.2)
      endif

c.....wall absorption fraction
      do jwall=1,nowals
        write(nunit,2040)jwall
        write(nunit,2041)(frabsorb(i,1,jwall,1),i=1,nosegsxz(jwall))
      enddo
 2040 format(1x,"frabsorb(1,1,",i1,",1)=")
 2041 format(19x,5f9.4)

c.....2-d grid coordinates
      do iv=1,nptshb
          write(nunit,2050)iv
          write(nunit,2051)(gridx(jh,iv,1),jh=1,nptsvb)
          write(nunit,2052)iv
          write(nunit,2051)(gridz(jh,iv,1),jh=1,nptsvb)
      enddo
 2050 format(1x,'gridx(1,',i3,',1)=')
 2051 format(5(1x,f10.2))
 2052 format(1x,'gridz(1,',i3,',1)=')

c.....flight splitting factors for plasma zones
      if (nosplits .gt. 1) then
        do i=1,novzs
          write(nunit,2055)i
          write(nunit,2056)(ksplzone(j,i,1),j=1,nohzs)
        enddo
 2055   format(1x,'ksplzone(1,',i3,',1)=')
 2056   format(20(1x,i2))
      endif

c.....wall materials
      do jwall=1,nowals
        write(nunit,2060)jwall
        write(nunit,2061)(kwmat(i,jwall),i=1,nosegsxz(jwall))
      enddo
 2060 format(1x,"kwmat(1,",i1,")= ")
 2061 format(7(1x,1a8,1x))

c.....kzone arrays
        do j=1,nptshb-1
          write(nunit,2070)j
          write(nunit,2071)(kzone1(i,j),i=1,nptsvb-1)
        enddo
        do j=1,nptshb-1
          write(nunit,2072)j
          write(nunit,2071)(kzone2(i,j),i=1,nptsvb-1)
        enddo
 2070   format(1x,'kzone1(1,',i3,')=')
 2071   format(20(1x,i2))
 2072   format(1x,'kzone2(1,',i3,')=')

c.....plasma zone boundaries
        write(nunit,2080)
        write(nunit,2081)(lboun1(i),i=1,nohbs)
        write(nunit,2082)
        write(nunit,2081)(lboun2(i),i=1,novbs)
 2080   format(1x,"lboun1=")
 2081   format(24(1x,i2))
 2082   format(1x,"lboun2=")

c.....electron temperature
      do i=1,novzs
      write(nunit,2090)i,(tehvt(k,i,1),k=1,min0(6,nohzs))
         do j=7,nohzs,6
         write(nunit,2091)(tehvt(k,i,1),k=j,min0(j+5,nohzs))
         enddo
      enddo
 2090 format(" tehvt(1,",i2,",1)=",1p6e10.3)
 2091 format(15x,1p6e10.3)

c.....ion temperature
      do i=1,novzs
      write(nunit,2100)i,(tihvt(k,i,1,1),k=1,min0(6,nohzs))
         do j=7,nohzs,6
         write(nunit,2101)(tihvt(k,i,1,1),k=j,min0(j+5,nohzs))
         enddo
      enddo
 2100 format(" tihvt(1,",i2,",1,1)=",1p6e10.3)
 2101 format(17x,1p6e10.3)

c.....ion flow velocity
      do i=1,novzs
      write(nunit,2102)i,(vflowx(k,i,1),k=1,min0(6,nohzs))
         do j=7,nohzs,6
         write(nunit,2103)(vflowx(k,i,1),k=j,min0(j+5,nohzs))
         enddo
      enddo
 2102 format(" vflowx(1,",i2,",1)= ",1p6e10.3)
 2103 format(17x,1p6e10.3)
      do i=1,novzs
      write(nunit,2104)i,(vflowy(k,i,1),k=1,min0(6,nohzs))
         do j=7,nohzs,6
         write(nunit,2105)(vflowy(k,i,1),k=j,min0(j+5,nohzs))
         enddo
      enddo
 2104 format(" vflowy(1,",i2,",1)= ",1p6e10.3)
 2105 format(17x,1p6e10.3)
      do i=1,novzs
      write(nunit,2106)i,(vflowz(k,i,1),k=1,min0(6,nohzs))
         do j=7,nohzs,6
         write(nunit,2107)(vflowz(k,i,1),k=j,min0(j+5,nohzs))
         enddo
      enddo
 2106 format(" vflowz(1,",i2,",1)= ",1p6e10.3)
 2107 format(17x,1p6e10.3)

c.....wall temperature
      do jwall=1,nowals
        write(nunit,2110)jwall
        write(nunit,2111)(twall(i,1,jwall),i=1,nosegsxz(jwall))
      enddo
 2110 format(1x,"twall(1,1,",i1,")=")
 2111 format(5(1x,f10.2))

c.....wall node coordinates
        do jwall=1,nowals
          write(nunit,2120)jwall
          write(nunit,2121)(xwall(i,jwall),i=1,nosegsxz(jwall)+1)
          write(nunit,2122)jwall
          write(nunit,2121)(zwall(i,jwall),i=1,nosegsxz(jwall)+1)
        enddo
 2120   format(1x,'xwall(1,',i1,')=')
 2121   format(5(1x,f10.2))
 2122   format(1x,'zwall(1,',i1,')=')

c    
      write(nunit,3000)endmark
 3000 format(1a8)

      close(nunit)

      call remark("***** Wrote DEGAS input file degas.in *****")

      return
      end

c----------------------------------------------------------------------c

