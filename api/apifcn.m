c-----------------------------------------------------------------------
      subroutine getatau(nx, ny, uu, gx, ixpt1, ixpt2, iysptrx,
     .                   atau, tau1, tau2)
      implicit none

c     This subroutine calculates the time for ions
c     to escape along open flux surfaces to either the left (inboard) or
c     right (outboard) divertor plate.

c ... Input arguments:
      integer nx   # poloidal dimension of mesh (excluding boundaries)
      integer ny   # radial   dimension of mesh (excluding boundaries)
      real uu(0:nx+1,0:ny+1)   # poloidal ion velocity
      real gx(0:nx+1,0:ny+1)   # 1/(x-width) of primary mesh cells
      integer ixpt1    # ix of last private-flux cell before cut on left
      integer ixpt2    # ix of last core-plasma cell before cut on right
      integer iysptrx  # iy of cell just below the separatrix

c ... Output arguments:
      real atau(0:nx+1,0:ny+1)   # lifetime of impurity
      real tau1(0:nx+1,0:ny+1)   # time to escape to inboard div. plate
      real tau2(0:nx+1,0:ny+1)   # time to escape to outboard div. plate

c     local variables --
      integer ix, iy
      real u, dxt, taumax

c     initialization --
c     Upper limit for particle confinement time
      taumax = 1.0e+00
c     Default values for core plasma confinement
      do ix = ixpt1+1, ixpt2
         do iy = 1, iysptrx
            tau1(ix,iy) = taumax
            tau2(ix,iy) = taumax
         enddo
      enddo

c     ------------------------------------------------------------------
c     For particles which flow toward the left plate --
c     ------------------------------------------------------------------

c     on private flux surfaces:
c
      do iy=1,iysptrx
	 tau1(0,iy) = 0.
	 do ix=1,ixpt1
	    u = uu(ix-1,iy)
	    dxt = 0.5/gx(ix-1,iy) + 0.5/gx(ix,iy)
	    if (u .lt. 0.) then
	       tau1(ix,iy) = min (taumax, tau1(ix-1,iy) + dxt / abs(u))
	    else
	       tau1(ix,iy) = taumax
	    endif
	 enddo
	 u = uu(ixpt2,iy)
	 dxt = 0.5/gx(ixpt1,iy) + 0.5/gx(ixpt2+1,iy)
	 if (u .lt. 0.) then
	    tau1(ixpt2+1,iy) = min (taumax, tau1(ixpt1,iy) + dxt/abs(u))
	 else
	    tau1(ixpt2+1,iy) = taumax
	 endif
	 do ix=ixpt2+2,nx
	    u = uu(ix-1,iy)
	    dxt = 0.5/gx(ix-1,iy) + 0.5/gx(ix,iy)
	    if (u .lt. 0.) then
	       tau1(ix,iy) = min (taumax, tau1(ix-1,iy) + dxt / abs(u))
	    else
	       tau1(ix,iy) = taumax
	    endif
	 enddo
	 tau1(nx+1,iy) = taumax
      enddo

c     on external flux surfaces :
c
      do iy=iysptrx+1,ny
	 tau1(0,iy) = 0.
	 do ix=1,nx
	    u = uu(ix-1,iy)
	    dxt = 0.5/gx(ix-1,iy) + 0.5/gx(ix,iy)
	    if (u .lt. 0.) then
	       tau1(ix,iy) = min (taumax, tau1(ix-1,iy) + dxt / abs(u))
	    else
	       tau1(ix,iy) = taumax
	    endif
	 enddo
	 tau1(nx+1,iy) = taumax
      enddo

c     ------------------------------------------------------------------
c     For particles which flow toward the right plate --
c     ------------------------------------------------------------------

c     on private flux surfaces:
c
      do iy=1,iysptrx
	 tau2(nx+1,iy) = 0.
	 do ix=nx,ixpt2+1,-1
	    u = uu(ix,iy)
	    dxt = 0.5/gx(ix,iy) + 0.5/gx(ix+1,iy)
	    if (u .gt. 0.) then
	       tau2(ix,iy) = min (taumax, tau2(ix+1,iy) + dxt / abs(u))
	    else
	       tau2(ix,iy) = taumax
	    endif
	 enddo
	 u = uu(ixpt2,iy)
	 dxt = 0.5/gx(ixpt1,iy) + 0.5/gx(ixpt2+1,iy)
	 if (u .gt. 0.) then
	    tau2(ixpt1,iy) = min (taumax, tau2(ixpt2+1,iy) + dxt/abs(u))
	 else
	    tau2(ixpt1,iy) = taumax
	 endif
	 do ix=ixpt1-1,1,-1
	    u = uu(ix,iy)
	    dxt = 0.5/gx(ix,iy) + 0.5/gx(ix+1,iy)
	    if (u .gt. 0.) then
	       tau2(ix,iy) = min (taumax, tau2(ix+1,iy) + dxt / abs(u))
	    else
	       tau2(ix,iy) = taumax
	    endif
	 enddo
	 tau2(0,iy) = taumax
      enddo

c     on external flux surfaces :
c
      do iy=iysptrx+1,ny
	 tau2(nx+1,iy) = 0.
	 do ix=nx,1,-1
	    u = uu(ix,iy)
	    dxt = 0.5/gx(ix,iy) + 0.5/gx(ix+1,iy)
	    if (u .gt. 0.) then
	       tau2(ix,iy) = min (taumax, tau2(ix+1,iy) + dxt / abs(u))
	    else
	       tau2(ix,iy) = taumax
	    endif
	 enddo
	 tau2(0,iy) = taumax
      enddo

c     ------------------------------------------------------------------
c     Combined loss to either plate --
c     ------------------------------------------------------------------

      do ix = 1, nx
         do iy = 1, ny
            atau(ix,iy) = min (tau1(ix,iy), tau2(ix,iy))
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine getprad(nx, ny, ngsp, te, ne, ng, afrac, atau,
     .                   prad, na, ntau, nratio)
      implicit none

c ... Input arguments:
      integer nx     # poloidal dimension of mesh (excluding boundaries)
      integer ny     # radial   dimension of mesh (excluding boundaries)
      integer ngsp   # number of gas species
      real te(0:nx+1,0:ny+1)      # electron temperature
      real ne(0:nx+1,0:ny+1)      # electron density
      real ng(0:nx+1,0:ny+1,1:ngsp)    # gas density
      real afrac(0:nx+1,0:ny+1)   # atomic concentration of impurity
      real atau(0:nx+1,0:ny+1)    # lifetime of impurity

c ... Output arguments:
      real prad(0:nx+1,0:ny+1)    # electron energy loss due to
                                  # impurity radiation
      real na(0:nx+1,0:ny+1)      # atomic density of impurity
      real ntau(0:nx+1,0:ny+1)    # confinement parameter for impurity
      real nratio(0:nx+1,0:ny+1)  # (neutral density) / (electron dens)

c     local variables --
      integer ix,iy,is

c     procedures --
      real emissbs
      external emissbs

c     compute energy loss for electrons --

      do ix=1,nx
         do iy=1,ny
            na(ix,iy) = afrac(ix,iy) * ne(ix,iy)
            ntau(ix,iy) = atau(ix,iy) * ne(ix,iy)
            nratio(ix,iy) = 0.
            do is=1,ngsp
               nratio(ix,iy) = nratio(ix,iy) + ng(ix,iy,is)
            enddo
            nratio(ix,iy) = nratio(ix,iy)/ne(ix,iy)
            prad(ix,iy) = na(ix,iy) * ne(ix,iy) *
     &            emissbs (te(ix,iy), nratio(ix,iy), ntau(ix,iy))
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine sapitim (timingsw)
      implicit none
      integer timingsw
      Use(Timing)

      istimingon = timingsw

      return
      end
c-----------------------------------------------------------------------
      subroutine wapitim
c ... Write timing data for the api package.
      implicit none
      Use(Timing)

      write(*,902) 'Impur.:  physics w/o bookkeeping = ', ttimpc
 902  format(a36,f10.4,20x,' sec')
      write(*,*) '(included in above f & Jac numbers)'
      ttimpc = 0.

      return
      end
c-----------------------------------------------------------------------

      subroutine readrates(apidir,impfname)
      implicit none
Use(Dim)
Use(Emissivities)
      character*(*) apidir
      character*(*) impfname
      integer impunit, ios

c     Reads Isler's excitation rate data
c     Based on Gary Porter's read_rates script

c ... Function:
      integer utgetcl   # defined in the Basis runtime library

c ... Local variable:
      integer stringlen, MAXSTRING
      parameter (MAXSTRING=500)
      character*(MAXSTRING) apidirx

c ... Get length of string containing directory name (omitting trailing
c     blanks). Basis function basfilex expands $,~.
      call basfilex(apidir,apidirx)
      stringlen = utgetcl(apidirx)

      call freeus(impunit)
      open (impunit, file=apidirx(1:stringlen) // '/' // impfname,
     .     form='formatted', iostat=ios, status='old')
      if (ios .ne. 0) then
         write(*,*) 
     .      '*** Input file ',impfname,' not found'
         call xerrab("")
      else
         write(*,*)
     .      '*** Reading from impurity excitation rate file: ',impfname
      endif

      read (impunit,*) ntemp
      read (impunit,*) nlam
      read (impunit,*) nden

      call gchange("Emissivities",0)

      call readrates1(impunit)

      return
      end
c-----------------------------------------------------------------------

      subroutine readrates1(impunit)
      implicit none
      integer impunit
Use(Dim)
Use(Emissivities)
      integer ii,jj,kk
      character*8 zdum

      read (impunit, '(6f12.2)') (lamb(ii),ii=1,nlam)
      write (*, '(6f12.2)') (lamb(ii),ii=1,nlam)
      do ii=1,ntemp
         read (impunit, '(a8)') zdum
         read (impunit, '(a8)') zdum
         read (impunit, '(f9.2)') etemp(ii)
         read (impunit, '(a8)') zdum
         do jj=1,nden
            read (impunit, '(a8)') zdum
            read (impunit, '(1x,e12.3,f8.2)') eden(jj),etemp(ii)
            read (impunit, '(6e12.3)') (rate(kk,ii,jj),kk=1,nlam)
         enddo
      enddo

      close(impunit)

c     Convert to m^-3 for UEDGE compatibility
      do ii=1,nden
         eden(ii)=eden(ii)*1e6
      enddo

      write (*,*)
     .  'Emission rate [ph/m^3/s] is rate*(appropriate density)'

      return
      end
c-----------------------------------------------------------------------

      subroutine calcrates(ne,te,density)
      implicit none
Use(Dim)
Use(Emissivities)
      real ne(0:nx+1,0:ny+1),te(0:nx+1,0:ny+1),density(0:nx+1,0:ny+1)
      real ev
      real newrate
      integer ii,ix,iy,ij,indxden,indxte

      data ev /1.6022e-19/

c     Calculates emission rates for each cell and wavelength:
c     Interpolate on (ne,te) in Isler data to get coeff for each line.
c     Multiply by appropriate input density to get emissivity.

      do ix=0,nx+1
        do iy=0,ny+1
          indxden=1
          indxte=1
          do ij=1,nden
            if (ne(ix,iy) > eden(ij)) indxden=ij
          enddo
          do ij=1,ntemp
            if (te(ix,iy)/ev > etemp(ij)) indxte=ij
          enddo
          do ii=1,nlam
	    if (indxte .eq. 1) then
	      newrate=0.
	    elseif ((indxden .eq. 1) .and. (indxte .lt. ntemp)) then
	      newrate=rate(ii,indxte,1) +
     .              (te(ix,iy)/ev-etemp(indxte))*
     .              (rate(ii,indxte+1,1)-rate(ii,indxte,1))/
     .              (etemp(indxte+1)-etemp(indxte))
	      newrate=newrate*ne(ix,iy)/eden(1)
	    elseif ((indxden .ge. nden) .and. (indxte .ge. ntemp)) then
	      newrate=rate(ii,ntemp,nden)*ne(ix,iy)/eden(nden)
	    elseif ((indxden .ge. nden) .and. (indxte .lt. ntemp)) then
	      newrate=rate(ii,indxte,nden) +
     .              (te(ix,iy)/ev-etemp(indxte))*
     .  	    (rate(ii,indxte+1,nden)-rate(ii,indxte,nden))/
     .   	    (etemp(indxte+1)-etemp(indxte))
	      newrate=newrate*ne(ix,iy)/eden(nden)
	    elseif ((indxden .lt. nden) .and. (indxte .ge. ntemp)) then
	      newrate=rate(ii,ntemp,indxden) +
     .  	    (ne(ix,iy)-eden(indxden))*
     .  	    (rate(ii,ntemp,indxden+1)-rate(ii,ntemp,indxden))/
     .  	    (eden(indxden+1)-eden(indxden))
	    else
	      newrate=rate(ii,indxte,indxden) +
     .  	    (te(ix,iy)/ev-etemp(indxte))*
     .  	    (rate(ii,indxte+1,indxden)-rate(ii,indxte,indxden))/
     .  	    (etemp(indxte+1)-etemp(indxte)) +
     .   	    (ne(ix,iy)-eden(indxden))*
     .  	    (rate(ii,indxte,indxden+1)-rate(ii,indxte,indxden))/
     .  	    (eden(indxden+1)-eden(indxden))
	    endif
	    emiss(ii,ix,iy)=newrate*density(ix,iy)
	  enddo
	enddo
      enddo

      return
      end
c-----------------------------------------------------------------------

      real function lineintegral(arg,rvertex,zvertex)
      implicit none
Use(Pixels)
      real arg(nrpix,nzpix)
      real rvertex(2),zvertex(2)
      real r1,r2,z1,z2,rp,zp
      integer iv1,iv2,jv1,jv2,k,ii,jj
      integer j,jbeg,jend,jstep
      integer i,ibeg,iend,istep

c MER 02 Jul 2003
c Compute the line integral of the 2-d array arg(nrpix,nzpix) along the
c path (rvertex(1:2),zvertex(1:2)).  We assume that arg(ii,jj) is the
c pixel representation of some UEDGE array arg_ue(ix,iy), obtained by mapping
c the UEDGE cells (rm,zm) onto a rectangular pixel domain via the DCE
c subroutine rzxform, e.g.,
c
c    arg = rzxform(arg_ue,rm,zm,nrpix,nzpix,rminpix,rmaxpix,zminpix,zmaxpix)
c
c The path integral is obtained via an algorithm which finds the
c intersection points of the line-of-sight (LOS) with horizontal
c and vertical pixel boundaries; these points then define the path
c length of the LOS in each cell.

c Initialize arrays:
      npd=0
      rp1=0.
      zp1=0.
      rp2=0.
      zp2=0.
      wt=0.

c Find the indices of the pixels containing the LOS vertices:
      r1=rvertex(1)
      z1=zvertex(1)
      r2=rvertex(2)
      z2=zvertex(2)
      iv1=int((r1-rminpix)/drpix)+1	# rpixel zone containing vertex 1
      jv1=int((z1-zminpix)/dzpix)+1	# zpixel zone containing vertex 1
      iv2=int((r2-rminpix)/drpix)+1	# rpixel zone containing vertex 2
      jv2=int((z2-zminpix)/dzpix)+1	# zpixel zone containing vertex 2

c Find the intersection of the LOS with the horizontal cell boundaries
c between the LOS vertices:
      if (jv2 .ge. jv1) then
         jbeg=max(0,jv1)       # first horiz. boundary above vertex 1
         jend=min(nzpix,jv2-1) # last  horiz. boundary below vertex 2
	 jstep=1
      else
         jbeg=min(nzpix,jv1-1) # last  horiz. boundary below vertex 1
         jend=max(0,jv2)       # first horiz. boundary above vertex 2
	 jstep=-1
      endif
      do j=jbeg,jend,jstep
	 zp=zminpix+j*dzpix
	 rp=r1+(r2-r1)*(zp-z1)/(z2-z1)
	 jj=j		#zp lies between zpixel zones jj and jj+1
	 if (rp .gt. rminpix) then
	    ii=int((rp-rminpix)/drpix)+1  #rp lies in rpixel zone ii
         else
            ii=0
         endif
c Assign this point to pixels in rows jj and jj+1
         if ( (1 .le. ii  ) .and. (ii   .le. nrpix) .and.
     .        (1 .le. jj+1) .and. (jj+1 .le. nzpix) ) then
           if (npd(ii,jj+1)==0) then
	     rp1(ii,jj+1)=rp
	     zp1(ii,jj+1)=zp
	     npd(ii,jj+1)=npd(ii,jj+1)+1
	   elseif (npd(ii,jj+1)==1) then
	     rp2(ii,jj+1)=rp
	     zp2(ii,jj+1)=zp
	     npd(ii,jj+1)=npd(ii,jj+1)+1
  	   else
             write (*,*) 'error:  ',
     .         'tried to assign more than 2 points to jj+1 pixel (',
     .          ii,',',jj+1,')'
         call xerrab("")
	   endif
         endif
	 if ( (1 .le. ii) .and. (ii .le. nrpix) .and.
     .        (1 .le. jj) .and. (jj .le. nzpix) ) then
	   if (npd(ii,jj)==0) then
	     rp1(ii,jj)=rp
	     zp1(ii,jj)=zp
	     npd(ii,jj)=npd(ii,jj)+1
	   elseif (npd(ii,jj)==1) then
	     rp2(ii,jj)=rp
	     zp2(ii,jj)=zp
	     npd(ii,jj)=npd(ii,jj)+1
	   else
	     write (*,*) 'error:  ',
     .         'tried to assign more than 2 points to jj pixel (',
     .          ii,',',jj,')'
         call xerrab("")
	   endif
         endif
      enddo

c Find the intersection of the LOS with the vertical cell boundaries
c between the LOS vertices:
      if (iv2 .ge. iv1) then
	 ibeg=max(0,iv1)       # first vertical boundary right of vertex 1
	 iend=min(nrpix,iv2-1) # last  vertical boundary left  of vertex 2
	 istep=1
      else
         ibeg=min(nrpix,iv1-1) # last  vertical boundary left  of vertex 1
         iend=max(0,iv2)       # first vertical boundary right of vertex 2
	 istep=-1
      endif
      do i=ibeg,iend,istep
	 rp=rminpix+i*drpix
	 zp=z1+(z2-z1)*(rp-r1)/(r2-r1)
	 ii=i		#rp lies between rpixel zones ii and ii+1
	 if (zp .gt. zminpix) then
	    jj=int((zp-zminpix)/dzpix)+1  #zp lies in zpixel zone jj
         else
            jj=0
         endif
c Assign this point to pixels in column ii and ii+1
         if ( (1 .le. ii+1) .and. (ii+1 .le. nrpix) .and.
     .        (1 .le. jj  ) .and. (jj   .le. nzpix) ) then
	   if (npd(ii+1,jj)==0) then
	     rp1(ii+1,jj)=rp
	     zp1(ii+1,jj)=zp
	     npd(ii+1,jj)=npd(ii+1,jj)+1
	   elseif (npd(ii+1,jj)==1) then
	     rp2(ii+1,jj)=rp
	     zp2(ii+1,jj)=zp
	     npd(ii+1,jj)=npd(ii+1,jj)+1
	   else
	     write (*,*) 'error:  ',
     .         'tried to assign more than 2 points to ii+1 pixel (',
     .          ii+1,',',jj,')'
         call xerrab("")
	   endif
         endif
         if ( (1 .le. ii) .and. (ii .le. nrpix) .and.
     .        (1 .le. jj) .and. (jj .le. nzpix) ) then
	   if (npd(ii,jj)==0) then
	     rp1(ii,jj)=rp
	     zp1(ii,jj)=zp
	     npd(ii,jj)=npd(ii,jj)+1
	   elseif (npd(ii,jj)==1) then
	     rp2(ii,jj)=rp
	     zp2(ii,jj)=zp
	     npd(ii,jj)=npd(ii,jj)+1
	   else
	     write (*,*) 'error:  ',
     .         'tried to assign more than 2 points to ii pixel (',
     .          ii,',',jj,')'
         call xerrab("")
	   endif
         endif
      enddo

c End-point contributions:
c Assign LOS vertices to pixels (iv1,jv1) and (iv2,jv2):
      if ( (1 .le. iv1) .and. (iv1 .le. nrpix) .and.
     .     (1 .le. jv1) .and. (jv1 .le. nzpix) ) then
        if (npd(iv1,jv1)==1) then
 	  rp2(iv1,jv1)=r1
          zp2(iv1,jv1)=z1
          npd(iv1,jv1)=npd(iv1,jv1)+1
        else
          write (*,*) 'error in end-point 1 pixel (',iv1,',',jv1,')'
         call xerrab("")
        endif
      endif
      if ( (1 .le. iv2) .and. (iv2 .le. nrpix) .and.
     .     (1 .le. jv2) .and. (jv2 .le. nzpix) ) then
        if (npd(iv2,jv2)==1) then
          rp2(iv2,jv2)=r2
          zp2(iv2,jv2)=z2
          npd(iv2,jv2)=npd(iv2,jv2)+1
        else
          write (*,*) 'error in end-point 2 pixel (',iv2,',',jv2,')'
         call xerrab("")
        endif
      endif

c Compute the weights (path lengths) for each pixel and sum contributions
c to lineintegral:
      lineintegral=0.
      do ii=max(1,min(nrpix,iv1)),max(1,min(nrpix,iv2)),istep
         do jj=max(1,min(nzpix,jv1)),max(1,min(nzpix,jv2)),jstep
            wt(ii,jj) = sqrt( (rp1(ii,jj)-rp2(ii,jj))**2 + 
     .                        (zp1(ii,jj)-zp2(ii,jj))**2 )
            lineintegral=lineintegral+wt(ii,jj)*arg(ii,jj)
         enddo
      enddo

      return
      end

