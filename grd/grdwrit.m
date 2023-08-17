c!include "../mppl.h"
c------------------------------------------------------------------

      subroutine grd2wdf

Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Aeqflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Curves)
Use(Linkco)
      integer nunit

# ********* Write the file to be read by the wdf package *************

      call freeus (nunit)
      open (nunit, file='grd-wdf', form='unformatted', status='unknown')
      write(nunit) idim,jdim,nix,noregs
      write(nunit) cmeshx,cmeshy,
     &              ilmax,ixpoint,jmin,jmax,jsptrx,jaxis
      write(nunit) bcentr,rcentr,rmagx,zmagx,simagx,sibdry,
     &              rgrid1,xdim,zdim
      write(nunit) nlim
      write(nunit) xlim,ylim
      write(nunit) eshot,etime,rseps,zseps,
     &              rvsin,zvsin,rvsout,zvsout
      close (nunit)

      return
      end

c-----------------------------------------------------------------------

      subroutine writeue
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)
Use(Share)   # geometry,islimon

      external writesn, writedn

c     Write the grid data file for the UEDGE code

      if ((geometry .eq. "snull").or.(geometry .eq. "uppersn")) then
         if (islimon .ne. 0) then
            call writelim ('gridue', runid)
         else
            call writesn ('gridue', runid)
         endif
      else
         call writedn ('gridue', runid)
      endif

      return
      end

c     ------------------------------------------------------------------

      subroutine writelim (fname, runidarg)
      implicit none
Use(Dim)          # nxm,nym
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)
Use(Share)        # nxxpt
      character*(*) fname, runidarg
      external gallot, wrdnbot, wrdndata

c     Write a grid data file for a limiter configuration

c     Map the PLANET/DEGAS grid into the B2/UEDGE grid --

c        write the full single null divertor configuration
c        excluding the PLANET guard cells near the x-points,
c        but including B2 guard cells at the limiter surface.
         nxm = (ilmax(1)-3) + (ilmax(2)-3) + 2 + 4*nxxpt
         nym = jmax(2) - jmin(2)
         call gallot("RZ_grid_info",0)
         call wrlim (fname, runidarg)

      return
      end

c     ------------------------------------------------------------------

      subroutine writedn (fname, runidarg)
      implicit none
Use(Dim)              # nxm,nym
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)
Use(Share)            # nxxpt
Use(UEgrid)           # ixtop
Use(Xpoint_indices)   # ixpt1,ixpt2,iysptrx1
Use(RZ_grid_info)     # rm,zm,rmt,zmt
      character*(*) fname, runidarg
      integer ixpt1t, ixtopt, ixpt2t, nxmt, iysptrxt, nymt
      external gallot

c     Write a grid data file

c     Map the PLANET/DEGAS grid into the B2/UEDGE grid --

      if (geometry .eq. "dnbot" .or. geometry == "isoleg") then
c        write the bottom half of the full double null geometry,
c        excluding the PLANET guard cells near the x-points,
c        but including B2 guard cells at the symmetry plane.
         if (ishalfm==0) then
            nxm = (ilmax(1)-3) + (ilmax(2)-3) + 2 + 4*nxxpt
         else
            nxm = (ilmax(1)-3) + (ilmax(2)-3) + 1 + 2*nxxpt
         endif
         nym = jmax(2) - jmin(2)
         call gallot("RZ_grid_info",0)
         call wrdnbot (fname, runidarg)
      elseif (geometry .eq. 'dnull') then
c        First, write half-space data in usual 'dnbot' format
         if (ishalfm==0) then
            nxm = (ilmax(1)-3) + (ilmax(2)-3) + 2 + 4*nxxpt
         else
            nxm = (ilmax(1)-3) + (ilmax(2)-3) + 1 + 2*nxxpt
         endif
         nym = jmax(2) - jmin(2)
         call gallot("RZ_grid_info",0)
         call wrdnbot (fname, runidarg)
ccc   NOTE: wrdnbot above contains un-necessary calls to magnetics
ccc         and to writedata
c        Save characteristic indices for half-space data:
         ixpt1t=ixpt1(1)
         ixtopt=ixtop
         ixpt2t=ixpt2(1)
         nxmt=nxm
         iysptrxt=iysptrx1(1)
         nymt=nym
c        Then, map data into full-space arrays, assuming symmetry:
         if (ishalfm==0) then
            nxm = 2*((ilmax(1)-3) + (ilmax(2)-3) + 2 + 4*nxxpt) - 2
         else
            nxm = 2*((ilmax(1)-3) + (ilmax(2)-3) + 1 + 2*nxxpt) - 2
         endif
         nym = jmax(2) - jmin(2)
         call gchange("RZ_grid_info",0)
         call mapdnbot2dnull(ixpt1t,ixtopt,ixpt2t,nxmt,iysptrxt,nymt)
         call add_guardc_tp
         call magnetics(0,nxm+1,1,nym)
         call symmetrize_magnetics
         if (isgriduehdf5 .eq. 1) then
            call parsestr('import uedge.gridue as gue;gue.write_gridue()')
         else
            call writednf (fname, runidarg)
         endif
      else
c        write the outboard half of the full double null geometry,
c        excluding the PLANET guard cells near the x-points.
         nxm = 2*(ilmax(2)-3) + 4*nxxpt
         nym = jmax(2) - jmin(2)
         call gallot("RZ_grid_info",0)
         call wrdndata (fname, runidarg)
      endif

      return
      end

c     ------------------------------------------------------------------

      subroutine wrdnbot (fname, runidarg)
      implicit none
Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm
Use(Xpoint_indices) # ixpt1,ixpt2,iysptrx1,iysptrx2
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(UEgrid)         # ixtop
Use(Inmesh)
Use(Linkco)
Use(Share)          # nxxpt,epslon,ishalfm
Use(Refinex)        # isrefxptn
      character*(*) fname, runidarg
      integer ix, iy, ii, jj
      external remark, xerrab, magnetics, writedata

#     Map the inboard half of the cmesh to the B2 mesh --
      iy = 0
      do jj = jmax(1),jmin(1)+1,-1
         iy = iy + 1
         ix = 0
         do ii = ilmax(1),ixpoint(3,1)+1,-1
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii-1,jj)
            zm(ix,iy,2) = cmeshy(ii-1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj-1)
            zm(ix,iy,3) = cmeshy(ii,jj-1)
            rm(ix,iy,4) = cmeshx(ii-1,jj-1)
            zm(ix,iy,4) = cmeshy(ii-1,jj-1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         ixpt1(1) = ix
         do ii=ixpoint(1,1),2,-1
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii-1,jj)
            zm(ix,iy,2) = cmeshy(ii-1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj-1)
            zm(ix,iy,3) = cmeshy(ii,jj-1)
            rm(ix,iy,4) = cmeshx(ii-1,jj-1)
            zm(ix,iy,4) = cmeshy(ii-1,jj-1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         if (ishalfm==0) then
c        insert a guard cell at the symmetry plane
         ix = ix + 1
         rm(ix,iy,1) = cmeshx(1,jj)
         zm(ix,iy,1) = cmeshy(1,jj)
         rm(ix,iy,2) = cmeshx(1,jj)
         zm(ix,iy,2) = cmeshy(1,jj) + epslon
         rm(ix,iy,3) = cmeshx(1,jj-1)
         zm(ix,iy,3) = cmeshy(1,jj-1)
         rm(ix,iy,4) = cmeshx(1,jj-1)
         zm(ix,iy,4) = cmeshy(1,jj-1) + epslon
         rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
         zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         endif      # end if-test on ishalfm
      enddo

      ixtop = ix


c     Map the outboard half of the cmesh to the B2 mesh --
      iy = 0
      do jj = jmin(2),jmax(2)-1
         iy = iy + 1
         ix = ixtop
c        insert a guard cell at the symmetry plane
         ix = ix + 1
         rm(ix,iy,1) = cmeshx(1,jj)
         zm(ix,iy,1) = cmeshy(1,jj) + epslon
         rm(ix,iy,2) = cmeshx(1,jj)
         zm(ix,iy,2) = cmeshy(1,jj)
         rm(ix,iy,3) = cmeshx(1,jj+1)
         zm(ix,iy,3) = cmeshy(1,jj+1) + epslon
         rm(ix,iy,4) = cmeshx(1,jj+1)
         zm(ix,iy,4) = cmeshy(1,jj+1)
         rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
         zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         do ii = 1,ixpoint(1,2)-1
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii+1,jj)
            zm(ix,iy,2) = cmeshy(ii+1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj+1)
            zm(ix,iy,3) = cmeshy(ii,jj+1)
            rm(ix,iy,4) = cmeshx(ii+1,jj+1)
            zm(ix,iy,4) = cmeshy(ii+1,jj+1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         ixpt2(1) = ix
         do ii = ixpoint(3,2),ilmax(2)-1
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii+1,jj)
            zm(ix,iy,2) = cmeshy(ii+1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj+1)
            zm(ix,iy,3) = cmeshy(ii,jj+1)
            rm(ix,iy,4) = cmeshx(ii+1,jj+1)
            zm(ix,iy,4) = cmeshy(ii+1,jj+1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
      enddo

      if (ix .ne. nxm-4*nxxpt) then
         call xerrab("*** ix indexing error in subroutine wrdnbot")
      endif
      if ( (jsptrx(2)-jmin(2)) .ne. (jmax(1)-jsptrx(1)) ) then
         call xerrab("*** iy indexing error in subroutine wrdnbot")
      endif
      iysptrx1(1) = jsptrx(2) - jmin(2)
      iysptrx2(1) = iysptrx1(1)

#     Add extra cells around the X-point if nxxpt > 0
      if (nxxpt .gt. 0) then
         if (isrefxptn==0) then
            call add_xptg       # uses interpolated flux surfaces
         elseif (isrefxptn==1) then
            call refine_xpt     # uses accurate flux surfaces
         endif
      endif

#     Next, evaluate the magnetic field and flux at each grid point --
      call magnetics(1,nxm,1,nym)

#     Finally, write out the data --
      call writedata (fname, runidarg)

      return
      end

c     ------------------------------------------------------------------

      subroutine wrdndata (fname, runidarg)
      implicit none
Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm
Use(Xpoint_indices) # ixpt1,ixpt2,iysptrx1,iysptrx2
Use(Dimflxgrd)	#nlim,nxefit,nyefit,jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)
Use(Share)          # nxxpt
Use(Refinex)        # isrefxptn
      character*(*) fname, runidarg
      integer ix, iy, ii, jj, ixlast
      external remark, xerrab, magnetics, writedata

#     Map the outboard half of the cmesh to the B2 mesh --
      iy = 0
      do jj = jmin(2),jmax(2)-1
         iy = iy + 1
         ix = nxm-4*nxxpt +1
         do ii = ilmax(2),ixpoint(3,2)+1,-1
            ix = ix - 1
            rm(ix,iy,1) = cmeshx(ii-1,jj)
            zm(ix,iy,1) = cmeshy(ii-1,jj)
            rm(ix,iy,2) = cmeshx(ii,jj)
            zm(ix,iy,2) = cmeshy(ii,jj)
            rm(ix,iy,3) = cmeshx(ii-1,jj+1)
            zm(ix,iy,3) = cmeshy(ii-1,jj+1)
            rm(ix,iy,4) = cmeshx(ii,jj+1)
            zm(ix,iy,4) = cmeshy(ii,jj+1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         ixpt2(1) = ix - 1
         do ii = ixpoint(1,2),2,-1
            ix = ix - 1
            rm(ix,iy,1) = cmeshx(ii-1,jj)
            zm(ix,iy,1) = cmeshy(ii-1,jj)
            rm(ix,iy,2) = cmeshx(ii,jj)
            zm(ix,iy,2) = cmeshy(ii,jj)
            rm(ix,iy,3) = cmeshx(ii-1,jj+1)
            zm(ix,iy,3) = cmeshy(ii-1,jj+1)
            rm(ix,iy,4) = cmeshx(ii,jj+1)
            zm(ix,iy,4) = cmeshy(ii,jj+1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo

#        Impose up/down symmetry -
         do ii = 2,ixpoint(1,2)
            ix = ix - 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = zdim - cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii-1,jj)
            zm(ix,iy,2) = zdim - cmeshy(ii-1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj+1)
            zm(ix,iy,3) = zdim - cmeshy(ii,jj+1)
            rm(ix,iy,4) = cmeshx(ii-1,jj+1)
            zm(ix,iy,4) = zdim - cmeshy(ii-1,jj+1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         ixpt1(1) = ix - 1
         do ii = ixpoint(3,2)+1,ilmax(2)
            ix = ix - 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = zdim - cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii-1,jj)
            zm(ix,iy,2) = zdim - cmeshy(ii-1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj+1)
            zm(ix,iy,3) = zdim - cmeshy(ii,jj+1)
            rm(ix,iy,4) = cmeshx(ii-1,jj+1)
            zm(ix,iy,4) = zdim - cmeshy(ii-1,jj+1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         ixlast = ix
      enddo
      if (ixlast .ne. 1) then
         call xerrab("*** ix indexing error in subroutine wrdndata ***")
      endif
      iysptrx1(1) = jsptrx(2) - jmin(2)
      iysptrx2(1) = iysptrx1(1)

#     Add extra cells around the X-point if nxxpt > 0
      if (nxxpt .gt. 0) then
         if (isrefxptn==0) then
            call add_xptg       # uses interpolated flux surfaces
         elseif (isrefxptn==1) then
            call refine_xpt     # uses accurate flux surfaces
         endif
      endif

#     Next, evaluate the magnetic field and flux at each grid point --
      call magnetics(1,nxm,1,nym)

#     Finally, write out the data --
      if (isgriduehdf5 .eq. 1) then
        call parsestr('import uedge.gridue as gue;gue.write_gridue()')
      else
        call writedata (fname, runidarg)
      endif

      return
      end

c     ------------------------------------------------------------------

      subroutine writesn (fname, runidarg)
      implicit none
Use(Dim)          # nxm,nym
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(Inmesh)
Use(Linkco)
Use(Share)        # nxxpt, isgriduehdf5
      character*(*) fname, runidarg
      external gallot, wrsndata

########################################################################
#     Write a grid data file for a single-null divertor

#     Dimension the B2 code grid arrays,
#     but exclude the PLANET guard cells near the x-points.
      nxm = (ilmax(1)-3) + (ilmax(2)-3) + 4*nxxpt
      nym = jmax(1) - jmin(1)

      call gallot("RZ_grid_info",0)

      call wrsndata (fname, runidarg)
      

      return
      end

c     ------------------------------------------------------------------

      subroutine wrsndata (fname, runidarg)
      implicit none
Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm
Use(Xpoint_indices) # ixpt1,ixpt2,iysptrx1,iysptrx2
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(UEgrid)         # ixtop
Use(Inmesh)
Use(Linkco)
Use(Share)          # nxxpt,ishalfm
Use(Refinex)        # isrefxptn
      character*(*) fname, runidarg
      integer ix, iy, ii, jj, ixlast, ixmin
      external remark, xerrab, magnetics, writedata

#     Map the inboard half of the cmesh to the B2 mesh --
      iy = 0
      do jj = jmax(1),jmin(1)+1,-1
         iy = iy + 1
         ix = 0
         do ii = ilmax(1),ixpoint(3,1)+1,-1
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii-1,jj)
            zm(ix,iy,2) = cmeshy(ii-1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj-1)
            zm(ix,iy,3) = cmeshy(ii,jj-1)
            rm(ix,iy,4) = cmeshx(ii-1,jj-1)
            zm(ix,iy,4) = cmeshy(ii-1,jj-1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         ixpt1(1) = ix
         do ii=ixpoint(1,1),2,-1
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii-1,jj)
            zm(ix,iy,2) = cmeshy(ii-1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj-1)
            zm(ix,iy,3) = cmeshy(ii,jj-1)
            rm(ix,iy,4) = cmeshx(ii-1,jj-1)
            zm(ix,iy,4) = cmeshy(ii-1,jj-1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         ixtop = ix
      enddo

#     Map the outboard half of the cmesh to the B2 mesh --

      iy = 0
      do jj = jmin(2),jmax(2)-1
         iy = iy + 1
         ix = ixtop
         do ii = 1, ixpoint(1,2)-1
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii+1,jj)
            zm(ix,iy,2) = cmeshy(ii+1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj+1)
            zm(ix,iy,3) = cmeshy(ii,jj+1)
            rm(ix,iy,4) = cmeshx(ii+1,jj+1)
            zm(ix,iy,4) = cmeshy(ii+1,jj+1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         ixpt2(1) = ix
         do ii = ixpoint(3,2),ilmax(2)-1
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii+1,jj)
            zm(ix,iy,2) = cmeshy(ii+1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj+1)
            zm(ix,iy,3) = cmeshy(ii,jj+1)
            rm(ix,iy,4) = cmeshx(ii+1,jj+1)
            zm(ix,iy,4) = cmeshy(ii+1,jj+1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         ixlast = ix
      enddo

      if (ixlast .ne. nxm-4*nxxpt) then
         call xerrab("*** ix indexing error in subroutine wrsndata ***")
      endif
      iysptrx1(1) = jmax(1) - jsptrx(1)
      iysptrx2(1) = iysptrx1(1)

#     Add extra cells around the X-point if nxxpt > 0
      if (nxxpt .gt. 0) then
         if (isrefxptn==0) then
            call add_xptg       # uses interpolated flux surfaces
         elseif (isrefxptn==1) then
            call refine_xpt     # uses accurate flux surfaces
         endif
      endif

#     Next, evaluate the magnetic field and flux at each grid point --
      ixmin = 1
      if (ishalfm==1) ixmin = ixtop + 1
      call magnetics(ixmin,nxm,1,nym)

#     Finally, write out the data --

      if (isgriduehdf5 .eq. 1) then
        call parsestr('import uedge.gridue as gue;gue.write_gridue()')
      else
        call writedata (fname, runidarg)
      endif

      return
      end

c     -------------------------------------------------------------------

      subroutine writedata (fname, runidarg)
      implicit none

Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm,psi,br,bz,bpol,bphi,b
Use(Xpoint_indices) # ixpt1,ixpt2,iysptrx1

      character*(*) fname, runidarg
      integer nunit, ix, iy, n

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

      write(nunit,1999) nxm,nym,ixpt1(1),ixpt2(1),iysptrx1(1)
 1999 format(5i4)
      write(nunit,2000)
      write(nunit,2001) (((rm(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((zm(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((psi(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((br(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((bz(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((bpol(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((bphi(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((b(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2002) runidarg

 2000 format()
ifelse([WORDSIZE],64,\
 2001 format(1p3e23.15)
,\
 2001 format(1p3d23.15)
)\
 2002 format(a60)
      close (nunit)

      write(*,*) 'Wrote file "', fname, '" with runid:  ', runidarg
      write(*,*)

      return
      end

c     ------------------------------------------------------------------

      subroutine magnetics(ixmin,ixmax,iymin,iymax)
      implicit none
      integer ixmin,ixmax,iymin,iymax
Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm,psi,br,bz,bpol,bphi,b
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)     # niwk,nwk
Use(Share)          # spheromak, isbphicon
      real psi00, psi10, psi01
      real B2VAhL, fpoloidal
      external B2INhT, B2VAhL, fpoloidal

      integer ix, iy, n

#     Evaluate the magnetic field and flux at B2 grid points
#
#     note: this mimics the procedures in the PLOTMED code

#     For tokamaks we assume the toroidal magnetic field is just
#     the vacuum field; for the spheromak configuration we use
#               B_toroidal = I_poloidal / R
#     where I_poloidal is the flux function, F, in the Grad-Shafronov
#     equation for the MHD equilibrium.

#     Set up the 2-d spline representation --
      call s2copy (nxefit,nyefit,fold,1,nxefit,bscoef,1,nxefit) 
      ldf = nxefit
      iflag = 1	# let B2INhT define the knots
      call B2INhT (xold, nxefit, yold, nyefit, kxord, kyord,
     .             xknot, yknot, bscoef, ldf, work, iflag)


#     Evaluate the spline and partial derivatives --
      do ix=ixmin,ixmax
         do iy=iymin,iymax
            do n=1,4
            psi00 =B2VAhL(rm(ix,iy,n),zm(ix,iy,n),0,0,xknot,yknot,
     .            nxefit,nyefit,kxord,kyord,bscoef,ldf,work,iflag)
            psi10 =B2VAhL(rm(ix,iy,n),zm(ix,iy,n),1,0,xknot,yknot,
     .            nxefit,nyefit,kxord,kyord,bscoef,ldf,work,iflag)
            psi01 =B2VAhL(rm(ix,iy,n),zm(ix,iy,n),0,1,xknot,yknot,
     .            nxefit,nyefit,kxord,kyord,bscoef,ldf,work,iflag)
            psi(ix,iy,n)  =   psi00
            bz(ix,iy,n)   = - psi10/rm(ix,iy,n)
            br(ix,iy,n)   =   psi01/rm(ix,iy,n)
            bpol(ix,iy,n) =   sqrt(bz(ix,iy,n)**2+br(ix,iy,n)**2)
            if (isbphicon == 0) then
               bphi(ix,iy,n) = bcentr*rcentr/rm(ix,iy,n)
            else
               bphi(ix,iy,n) = bcentr
            endif
            if (spheromak==1) then
               bphi(ix,iy,n) = fpoloidal(psi(ix,iy,n))/rm(ix,iy,n)
            endif
            b(ix,iy,n)    =   sqrt(bpol(ix,iy,n)**2+bphi(ix,iy,n)**2)
            enddo
            psi(ix,iy,0)  = 0.25*(psi(ix,iy,1)+psi(ix,iy,2)
     &                           +psi(ix,iy,3)+psi(ix,iy,4))  
            br(ix,iy,0)   = 0.25*(br(ix,iy,1)+br(ix,iy,2)
     &                           +br(ix,iy,3)+br(ix,iy,4))  
            bz(ix,iy,0)   = 0.25*(bz(ix,iy,1)+bz(ix,iy,2)
     &                           +bz(ix,iy,3)+bz(ix,iy,4))  
            bpol(ix,iy,0) = 0.25*(bpol(ix,iy,1)+bpol(ix,iy,2)
     &                           +bpol(ix,iy,3)+bpol(ix,iy,4))  
            bphi(ix,iy,0) = 0.25*(bphi(ix,iy,1)+bphi(ix,iy,2)
     &                           +bphi(ix,iy,3)+bphi(ix,iy,4))  
            b(ix,iy,0)    = 0.25*(b(ix,iy,1)+b(ix,iy,2)
     &                           +b(ix,iy,3)+b(ix,iy,4))  
         enddo
      enddo

      return
      end

c     ------------------------------------------------------------------

      subroutine wrlim (fname, runidarg)
      implicit none
Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm
Use(Xpoint_indices) # ixpt1,ixpt2,iysptrx1,iysptrx2
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)
Use(Dimensions)
Use(UEgrid)         # ixtop
Use(Inmesh)
Use(Linkco)
Use(Share)          # iy_lims, nxxpt
Use(Refinex)        # isrefxptn
      character*(*) fname, runidarg
      integer ix, iy, ii, jj, jjlim1
      real cmx, cmy
      external remark, xerrab, magnetics, writedata

ccc MER 97/12/29
c     Correction to ensure that inboard and outboard segments of flux 
c     surfaces iy < iy_lims meet at exactly the same point.
c     NOTE: iy_lims should already have been set by subroutine setlimindex
      jjlim1 = max( jmin(1), jmax(1)-iy_lims+1 )
      do jj = jjlim1,jmax(1)
         cmx = 0.5 * (cmeshx(1,jj)+cmeshx(1,jmax(2)-jj+jmin(1)))
         cmy = 0.5 * (cmeshy(1,jj)+cmeshy(1,jmax(2)-jj+jmin(1)))
         cmeshx(1,jj) = cmx
         cmeshx(1,jmax(2)-jj+jmin(1)) = cmx
         cmeshy(1,jj) = cmy
         cmeshy(1,jmax(2)-jj+jmin(1)) = cmy
      enddo
ccc MER 97/12/29

#     Map the inboard half of the cmesh to the B2 mesh --
      iy = 0
      do jj = jmax(1),jmin(1)+1,-1
         iy = iy + 1
         ix = 0
         do ii = ilmax(1),ixpoint(3,1)+1,-1
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii-1,jj)
            zm(ix,iy,2) = cmeshy(ii-1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj-1)
            zm(ix,iy,3) = cmeshy(ii,jj-1)
            rm(ix,iy,4) = cmeshx(ii-1,jj-1)
            zm(ix,iy,4) = cmeshy(ii-1,jj-1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         ixpt1(1) = ix
         do ii=ixpoint(1,1),2,-1
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii-1,jj)
            zm(ix,iy,2) = cmeshy(ii-1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj-1)
            zm(ix,iy,3) = cmeshy(ii,jj-1)
            rm(ix,iy,4) = cmeshx(ii-1,jj-1)
            zm(ix,iy,4) = cmeshy(ii-1,jj-1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
c        insert a guard cell at the top of the inboard mesh
         ix = ix + 1
ccc   MER 97/10/10  interpolate (extrapolation produces cell overlap)
         rm(ix,iy,1)=rm(ix-1,iy,2)
     .                    -epslon_lim*(rm(ix-1,iy,2)-rm(ix-1,iy,1))
         zm(ix,iy,1)=zm(ix-1,iy,2)
     .                    -epslon_lim*(zm(ix-1,iy,2)-zm(ix-1,iy,1))
         rm(ix,iy,2)=rm(ix-1,iy,2)
         zm(ix,iy,2)=zm(ix-1,iy,2)
         rm(ix,iy,3)=rm(ix-1,iy,4)
     .                    -epslon_lim*(rm(ix-1,iy,4)-rm(ix-1,iy,3))
         zm(ix,iy,3)=zm(ix-1,iy,4)
     .                    -epslon_lim*(zm(ix-1,iy,4)-zm(ix-1,iy,3))
         rm(ix,iy,4)=rm(ix-1,iy,4)
         zm(ix,iy,4)=zm(ix-1,iy,4)
         rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
         zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
c     Now go back and correct the previous cell to avoid overlap
         rm(ix-1,iy,2)=rm(ix,iy,1)
         zm(ix-1,iy,2)=zm(ix,iy,1)
         rm(ix-1,iy,4)=rm(ix,iy,3)
         zm(ix-1,iy,4)=zm(ix,iy,3)
         rm(ix-1,iy,0) = 0.25*(rm(ix-1,iy,1)+rm(ix-1,iy,2)
     &                                 +rm(ix-1,iy,3)+rm(ix-1,iy,4))
         zm(ix-1,iy,0) = 0.25*(zm(ix-1,iy,1)+zm(ix-1,iy,2)
     &                                 +zm(ix-1,iy,3)+zm(ix-1,iy,4))
      enddo

      ixtop = ix

c     Map the outboard half of the cmesh to the B2 mesh --
      iy = 0
      do jj = jmin(2),jmax(2)-1
         iy = iy + 1
         ix = ixtop
c        insert a guard cell at the top of the outboard mesh
         ix = ix + 1
ccc   MER 97/10/10  interpolate (extrapolation produces cell overlap)
         rm(ix,iy,1)=cmeshx(1,jj)
         zm(ix,iy,1)=cmeshy(1,jj)
         rm(ix,iy,2)=cmeshx(1,jj)
     .                     +epslon_lim*(cmeshx(2,jj)-cmeshx(1,jj))
         zm(ix,iy,2)=cmeshy(1,jj)
     .                     +epslon_lim*(cmeshy(2,jj)-cmeshy(1,jj))
         rm(ix,iy,3)=cmeshx(1,jj+1)
         zm(ix,iy,3)=cmeshy(1,jj+1)
         rm(ix,iy,4)=cmeshx(1,jj+1)
     .                 +epslon_lim*(cmeshx(2,jj+1)-cmeshx(1,jj+1))
         zm(ix,iy,4)=cmeshy(1,jj+1)
     .                 +epslon_lim *(cmeshy(2,jj+1)-cmeshy(1,jj+1))
         rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
         zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         ix = ix + 1
         rm(ix,iy,1)=rm(ix-1,iy,2)
         zm(ix,iy,1)=zm(ix-1,iy,2)
         rm(ix,iy,2)=cmeshx(2,jj)
         zm(ix,iy,2)=cmeshy(2,jj)
         rm(ix,iy,3)=rm(ix-1,iy,4)
         zm(ix,iy,3)=zm(ix-1,iy,4)
         rm(ix,iy,4)=cmeshx(2,jj+1)
         zm(ix,iy,4)=cmeshy(2,jj+1)
         rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
         zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         do ii = 2,ixpoint(1,2)-1    # MER 97/10/10 start at ii=2 now
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii+1,jj)
            zm(ix,iy,2) = cmeshy(ii+1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj+1)
            zm(ix,iy,3) = cmeshy(ii,jj+1)
            rm(ix,iy,4) = cmeshx(ii+1,jj+1)
            zm(ix,iy,4) = cmeshy(ii+1,jj+1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
         ixpt2(1) = ix
         do ii = ixpoint(3,2),ilmax(2)-1
            ix = ix + 1
            rm(ix,iy,1) = cmeshx(ii,jj)
            zm(ix,iy,1) = cmeshy(ii,jj)
            rm(ix,iy,2) = cmeshx(ii+1,jj)
            zm(ix,iy,2) = cmeshy(ii+1,jj)
            rm(ix,iy,3) = cmeshx(ii,jj+1)
            zm(ix,iy,3) = cmeshy(ii,jj+1)
            rm(ix,iy,4) = cmeshx(ii+1,jj+1)
            zm(ix,iy,4) = cmeshy(ii+1,jj+1)
            rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                     +rm(ix,iy,3)+rm(ix,iy,4))
            zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                     +zm(ix,iy,3)+zm(ix,iy,4))
         enddo
      enddo

      if (ix .ne. nxm-4*nxxpt) then
         call xerrab("*** ix indexing error in subroutine wrlim")
      endif
      if ( (jsptrx(2)-jmin(2)) .ne. (jmax(1)-jsptrx(1)) ) then
         call xerrab("*** iy indexing error in subroutine wrlim")
      endif
      iysptrx1(1) = jsptrx(2) - jmin(2)
      iysptrx2(1) = iysptrx1(1)

#     Add extra cells around the X-point if nxxpt > 0
      if (nxxpt .gt. 0) then
         if (isrefxptn==0) then
            call add_xptg       # uses interpolated flux surfaces
         elseif (isrefxptn==1) then
            call refine_xpt     # uses accurate flux surfaces
         endif
      endif

#     Next, evaluate the magnetic field and flux at each grid point --
      call magnetics(1,nxm,1,nym)

#     Finally, write out the data --
      call writedata (fname, runidarg)

      return
      end

c     ------------------------------------------------------------------

      real function fpoloidal (psi)
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)
      real psi

c     local variables:
      integer i1
      real dpsi, f1, f2, lambda_edge, psi1

#     Evaluate the poloidal current function, F(psi) = R*B_toroidal,
#     for EFIT equilibria.

#     The EFIT array fpol(nxefit) in the geqdsk file gives F(psi) on a
#     uniform psi mesh; the first element is the value on the magnetic
#     axis and the last element is the value on the separatrix.

      dpsi = (sibdry-simagx)/(nxefit-1)
      i1 = int( (psi-simagx)/dpsi ) + 1

      if (i1 .lt. nxefit) then    # interpolate
         psi1 = simagx + (i1-1)*dpsi
         f1 = fpol(i1)
         f2 = fpol(i1+1)
         fpoloidal = f1 + (f2-f1)*(psi-psi1)/dpsi
      else                         # extrapolate
         lambda_edge = (fpol(nxefit) - fpol(nxefit-1))/dpsi
         fpoloidal = fpol(nxefit) + lambda_edge*(psi-sibdry)
      endif

      return
      end

c     ------------------------------------------------------------------

      real function pressure (psi)
      implicit none
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)
      real psi

c     local variables:
      integer i1
      real dpsi, f1, f2, psi1

#     Evaluate the plasma pressure function, P(psi), for EFIT equilibria.

#     The EFIT array pres(nxefit) in the geqdsk file gives the plasma
#     pressure on a uniform psi mesh; the first element is the value
#     on the magnetic axis and the last element is the value on the
#     separatrix.

      dpsi = (sibdry-simagx)/(nxefit-1)
      i1 = int( (psi-simagx)/dpsi ) + 1

      if (i1 .lt. nxefit) then	# interpolate
         psi1 = simagx + (i1-1)*dpsi
         f1 = pres(i1)
         f2 = pres(i1+1)
         pressure = f1 + (f2-f1)*(psi-psi1)/dpsi
      else			# return separatrix pressure
         pressure = pres(nxefit)
      endif

      return
      end

c     ------------------------------------------------------------------

      real function psif(r,z)
      implicit none
      real r,z
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)
      real B2VAhL
      external B2VAhL

      psif = B2VAhL(r,z,0,0,xknot,yknot,nxefit,nyefit,
     .                      kxord,kyord,bscoef,ldf,work,iflag)

      return
      end

c     ------------------------------------------------------------------

      real function brf(r,z)
      implicit none
      real r,z
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)
      real B2VAhL, dpsidz
      external B2VAhL

      dpsidz = B2VAhL(r,z,0,1,xknot,yknot,nxefit,nyefit,
     .                      kxord,kyord,bscoef,ldf,work,iflag)
      brf = dpsidz / r

      return
      end

c     ------------------------------------------------------------------

      real function bzf(r,z)
      implicit none
      real r,z
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)
      real B2VAhL, dpsidr
      external B2VAhL

      dpsidr = B2VAhL(r,z,1,0,xknot,yknot,nxefit,nyefit,
     .                      kxord,kyord,bscoef,ldf,work,iflag)
      bzf = - dpsidr / r

      return
      end

c     -------------------------------------------------------------------

      subroutine writednf (fname, runidarg)
      implicit none
      character*(*) fname, runidarg

Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm,psi,br,bz,bpol,bphi,b
Use(Xpoint_indices) # ixlb,ixpt1,ixmdp,ixpt2,ixrb,iysptrx1,iysptrx2

      integer nunit, ix, iy, n

c     Write gridue data file for full double-null configuration

      call freeus (nunit)
      open (nunit, file=fname, form='formatted', status='unknown')

      write(nunit,1999) nxm, nym
      write(nunit,1999) iysptrx1(1), iysptrx2(1)
      write(nunit,1999) ixlb(1), ixpt1(1), ixmdp(1), ixpt2(1), ixrb(1)
      write(nunit,1999) iysptrx1(2), iysptrx2(2)
      write(nunit,1999) ixlb(2), ixpt1(2), ixmdp(2), ixpt2(2), ixrb(2)
 1999 format(5i4)
      write(nunit,2000)
      write(nunit,2001) (((rm(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((zm(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((psi(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((br(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((bz(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((bpol(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((bphi(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2000)
      write(nunit,2001) (((b(ix,iy,n),ix=0,nxm+1),iy=0,nym+1),n=0,4)
      write(nunit,2002) runidarg

 2000 format()
ifelse([WORDSIZE],64,\
 2001 format(1p3e23.15)
,\
 2001 format(1p3d23.15)
)\
 2002 format(a60)
      close (nunit)

      write(*,*) 'Wrote file "', fname, '" with runid:  ', runidarg
      write(*,*)

      return
      end

c     ------------------------------------------------------------------

      subroutine add_xptg

c  -- This subroutine adds nxxpt poloidal mesh points (per quadrant)
c  -- around the X-point; call from the various grid-writing routines

      implicit none

Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm
Use(Xpoint_indices) # ixpt1,ixpt2
Use(Dimflxgrd)	#jdim,npts,noregs      
Use(Dimensions)
Use(UEgrid)         # ixtop
Use(Inmesh)
Use(Linkco)
Use(Share)          # nxxpt
Use(Refinex)        # alfxpt

c  ** Local variables
      integer ix, iy, ii, jj, ixpt1n, ixpt2n, ixpt2t
      real frac, dr0, dz0, drs, dzs

c ... Store orignal mesh in rmt, zmt
      do jj = 0, 4
        do iy = 1, nym
          do ix = 1, nxm-4*nxxpt
            rmt(ix,iy,jj) = rm(ix,iy,jj)
            zmt(ix,iy,jj) = zm(ix,iy,jj)
          enddo
        enddo
      enddo

c ... First redo the mesh for ix=ixpt1 region

      do ii = 0, nxxpt
        frac = (float(ii)/float(nxxpt+1))**alfxpt
        do iy = 1, nym
          do jj = 1, 3, 2
            if (iy .eq. 1 .or. iy .eq. nym) then
              rm(ixpt1(1)+ii,iy,jj) = rmt(ixpt1(1),iy,jj) + frac*
     .                        (rmt(ixpt1(1),iy,jj+1) - rmt(ixpt1(1),iy,jj))
              zm(ixpt1(1)+ii,iy,jj) = zmt(ixpt1(1),iy,jj) + frac*
     .                        (zmt(ixpt1(1),iy,jj+1) - zmt(ixpt1(1),iy,jj))
            else  # smooth in iy with (1 2 1) filter
              dr0 = rmt(ixpt1(1),iy,jj+1) - rmt(ixpt1(1),iy,jj)
              dz0 = zmt(ixpt1(1),iy,jj+1) - zmt(ixpt1(1),iy,jj)
              if (abs(dr0) .gt. abs(dz0)) then  # avoid dr0 --> 0
                drs = (frac/4) * (
     .              2*(rmt(ixpt1(1),iy  ,jj+1) - rmt(ixpt1(1),iy  ,jj))+
     .                (rmt(ixpt1(1),iy+1,jj+1) - rmt(ixpt1(1),iy+1,jj))+
     .                (rmt(ixpt1(1),iy-1,jj+1) - rmt(ixpt1(1),iy-1,jj)) ) 
                rm(ixpt1(1)+ii,iy,jj) = rmt(ixpt1(1),iy,jj) + drs
                zm(ixpt1(1)+ii,iy,jj) = zmt(ixpt1(1),iy,jj) + dz0*drs/dr0
              else                              # avoid dz0 --> 0 
                dzs = (frac/4) * (
     .              2*(zmt(ixpt1(1),iy  ,jj+1) - zmt(ixpt1(1),iy  ,jj))+
     .                (zmt(ixpt1(1),iy+1,jj+1) - zmt(ixpt1(1),iy+1,jj))+
     .                (zmt(ixpt1(1),iy-1,jj+1) - zmt(ixpt1(1),iy-1,jj)) ) 
                zm(ixpt1(1)+ii,iy,jj) = zmt(ixpt1(1),iy,jj) + dzs
                rm(ixpt1(1)+ii,iy,jj) = rmt(ixpt1(1),iy,jj) + dr0*dzs/dz0
              endif
            endif
          enddo
        enddo
      enddo

      do ii = 0, nxxpt-1
        do iy = 1, nym
          do jj = 2, 4, 2
            rm(ixpt1(1)+ii,iy,jj) = rm(ixpt1(1)+1+ii,iy,jj-1)
            zm(ixpt1(1)+ii,iy,jj) = zm(ixpt1(1)+1+ii,iy,jj-1)
          enddo
        enddo
      enddo

      do iy = 1, nym
        do jj = 2, 4, 2
          rm(ixpt1(1)+nxxpt,iy,jj) = rmt(ixpt1(1),iy,jj)
          zm(ixpt1(1)+nxxpt,iy,jj) = zmt(ixpt1(1),iy,jj)
        enddo
      enddo

      # Correct inside vertices at iy=1 and iy=nym to fit interior cells
      do ii = 0, nxxpt
         do jj = 3, 4
            rm(ixpt1(1)+ii,1,jj) = rm(ixpt1(1)+ii,2,jj-2)
            zm(ixpt1(1)+ii,1,jj) = zm(ixpt1(1)+ii,2,jj-2)
            rm(ixpt1(1)+ii,nym,jj-2) = rm(ixpt1(1)+ii,nym-1,jj)
            zm(ixpt1(1)+ii,nym,jj-2) = zm(ixpt1(1)+ii,nym-1,jj)
         enddo
      enddo

c ... Redo the mesh for the ix=ixpt1+1 region (note use of ixpt1n)

      ixpt1n = ixpt1(1) + nxxpt
      do ii = nxxpt, 0, -1   # start at ixpt1+2 and fill below
        frac =(float(nxxpt-ii)/float(nxxpt+1))**alfxpt
        do iy = 1, nym
          do jj = 1, 3, 2
            if (iy .eq. 1 .or. iy .eq. nym) then
              rm(ixpt1n+1+ii,iy,jj+1) = rmt(ixpt1(1)+2,iy,jj) - frac*
     .                      (rmt(ixpt1(1)+1,iy,jj+1) - rmt(ixpt1(1)+1,iy,jj))
              zm(ixpt1n+1+ii,iy,jj+1) = zmt(ixpt1(1)+2,iy,jj) - frac*
     .                      (zmt(ixpt1(1)+1,iy,jj+1) - zmt(ixpt1(1)+1,iy,jj))
            else  # smooth in iy with (1 2 1) filter
              dr0 = rmt(ixpt1(1)+1,iy,jj+1) - rmt(ixpt1(1)+1,iy,jj)
              dz0 = zmt(ixpt1(1)+1,iy,jj+1) - zmt(ixpt1(1)+1,iy,jj)
              if (abs(dr0) .gt. abs(dz0)) then  # avoid dr0 --> 0
                drs = (frac/4) * (
     .              2*(rmt(ixpt1(1)+1,iy  ,jj+1) - rmt(ixpt1(1)+1,iy  ,jj))+
     .                (rmt(ixpt1(1)+1,iy+1,jj+1) - rmt(ixpt1(1)+1,iy+1,jj))+
     .                (rmt(ixpt1(1)+1,iy-1,jj+1) - rmt(ixpt1(1)+1,iy-1,jj)) ) 
                rm(ixpt1n+1+ii,iy,jj+1) = rmt(ixpt1(1)+2,iy,jj)- drs
                zm(ixpt1n+1+ii,iy,jj+1)= zmt(ixpt1(1)+2,iy,jj)- dz0*drs/dr0
              else                              # avoid dz0 --> 0 
                dzs = (frac/4) * (
     .              2*(zmt(ixpt1(1)+1,iy  ,jj+1) - zmt(ixpt1(1)+1,iy  ,jj))+
     .                (zmt(ixpt1(1)+1,iy+1,jj+1) - zmt(ixpt1(1)+1,iy+1,jj))+
     .                (zmt(ixpt1(1)+1,iy-1,jj+1) - zmt(ixpt1(1)+1,iy-1,jj)) ) 
                zm(ixpt1n+1+ii,iy,jj+1) = zmt(ixpt1(1)+2,iy,jj)- dzs
                rm(ixpt1n+1+ii,iy,jj+1)= rmt(ixpt1(1)+2,iy,jj)- dr0*dzs/dz0
              endif
            endif
          enddo
        enddo
      enddo

      do ii = 1, nxxpt
        do iy = 1, nym
          do jj = 1, 3, 2
            rm(ixpt1n+1+ii,iy,jj) = rm(ixpt1n+ii,iy,jj+1)
            zm(ixpt1n+1+ii,iy,jj) = zm(ixpt1n+ii,iy,jj+1)
          enddo
        enddo
      enddo

      do iy = 1, nym
        do jj = 1, 3, 2
          rm(ixpt1n+1,iy,jj) = rmt(ixpt1(1)+1,iy,jj)
          zm(ixpt1n+1,iy,jj) = zmt(ixpt1(1)+1,iy,jj)
        enddo
      enddo

      # Correct inside vertices at iy=1 and iy=nym to fit interior cells
      do ii = 0, nxxpt
         do jj = 3, 4
            rm(ixpt1n+1+ii,1,jj) = rm(ixpt1n+1+ii,2,jj-2)
            zm(ixpt1n+1+ii,1,jj) = zm(ixpt1n+1+ii,2,jj-2)
            rm(ixpt1n+1+ii,nym,jj-2) = rm(ixpt1n+1+ii,nym-1,jj)
            zm(ixpt1n+1+ii,nym,jj-2) = zm(ixpt1n+1+ii,nym-1,jj)
         enddo
      enddo

c ... Redo the mesh for the ix=ixpt2 region (note use of ixpt2t)

      ixpt2t = ixpt2(1) + 2*nxxpt
      do ii = 0, nxxpt
        frac =(float(ii)/float(nxxpt+1))**alfxpt
        do iy = 1, nym
          do jj = 1, 3, 2
            if (iy .eq. 1 .or. iy .eq. nym) then
              rm(ixpt2t+ii,iy,jj) = rmt(ixpt2(1),iy,jj) + frac*
     .                        (rmt(ixpt2(1),iy,jj+1) - rmt(ixpt2(1),iy,jj))
              zm(ixpt2t+ii,iy,jj) = zmt(ixpt2(1),iy,jj) + frac*
     .                        (zmt(ixpt2(1),iy,jj+1) - zmt(ixpt2(1),iy,jj))
            else  # smooth in iy with (1 2 1) filter
              dr0 = rmt(ixpt2(1),iy,jj+1) - rmt(ixpt2(1),iy,jj)
              dz0 = zmt(ixpt2(1),iy,jj+1) - zmt(ixpt2(1),iy,jj)
              if (abs(dr0) .gt. abs(dz0)) then  # avoid dr0 --> 0
                drs = (frac/4) * (
     .              2*(rmt(ixpt2(1),iy  ,jj+1) - rmt(ixpt2(1),iy  ,jj))+
     .                (rmt(ixpt2(1),iy+1,jj+1) - rmt(ixpt2(1),iy+1,jj))+
     .                (rmt(ixpt2(1),iy-1,jj+1) - rmt(ixpt2(1),iy-1,jj)) ) 
                rm(ixpt2t+ii,iy,jj) = rmt(ixpt2(1),iy,jj) + drs
                zm(ixpt2t+ii,iy,jj) = zmt(ixpt2(1),iy,jj) + dz0*drs/dr0
              else                              # avoid dz0 --> 0 
                dzs = (frac/4) * (
     .              2*(zmt(ixpt2(1),iy  ,jj+1) - zmt(ixpt2(1),iy  ,jj))+
     .                (zmt(ixpt2(1),iy+1,jj+1) - zmt(ixpt2(1),iy+1,jj))+
     .                (zmt(ixpt2(1),iy-1,jj+1) - zmt(ixpt2(1),iy-1,jj)) ) 
                zm(ixpt2t+ii,iy,jj) = zmt(ixpt2(1),iy,jj) + dzs
                rm(ixpt2t+ii,iy,jj) = rmt(ixpt2(1),iy,jj) + dr0*dzs/dz0
              endif
            endif
          enddo
        enddo
      enddo

      do ii = 0, nxxpt-1
        do iy = 1, nym
          do jj = 2, 4, 2
            rm(ixpt2t+ii,iy,jj) = rm(ixpt2t+1+ii,iy,jj-1)
            zm(ixpt2t+ii,iy,jj) = zm(ixpt2t+1+ii,iy,jj-1)
          enddo
        enddo
      enddo

      do iy = 1, nym
        do jj = 2, 4, 2
          rm(ixpt2t+nxxpt,iy,jj) = rmt(ixpt2(1),iy,jj)
          zm(ixpt2t+nxxpt,iy,jj) = zmt(ixpt2(1),iy,jj)
        enddo
      enddo

      # Correct inside vertices at iy=1 and iy=nym to fit interior cells
      do ii = 0, nxxpt
         do jj = 3, 4
            rm(ixpt2t+ii,1,jj) = rm(ixpt2t+ii,2,jj-2)
            zm(ixpt2t+ii,1,jj) = zm(ixpt2t+ii,2,jj-2)
            rm(ixpt2t+ii,nym,jj-2) = rm(ixpt2t+ii,nym-1,jj)
            zm(ixpt2t+ii,nym,jj-2) = zm(ixpt2t+ii,nym-1,jj)
         enddo
      enddo

c ... Redo the mesh for the ix=ixpt2+1 region (note use of ixpt2n)

      ixpt2n = ixpt2(1) + 3*nxxpt
      do ii = nxxpt, 0, -1   # start at ixpt2(1)+2 and fill up
        frac = (float(nxxpt-ii)/float(nxxpt+1))**alfxpt
        do iy = 1, nym
          do jj = 1, 3, 2
            if (iy .eq. 1 .or. iy .eq. nym) then
              rm(ixpt2n+1+ii,iy,jj+1) = rmt(ixpt2(1)+2,iy,jj) - frac*
     .                      (rmt(ixpt2(1)+1,iy,jj+1) - rmt(ixpt2(1)+1,iy,jj))
              zm(ixpt2n+1+ii,iy,jj+1) = zmt(ixpt2(1)+2,iy,jj) - frac*
     .                      (zmt(ixpt2(1)+1,iy,jj+1) - zmt(ixpt2(1)+1,iy,jj))
            else  # smooth in iy with (1 2 1) filter
              dr0 = rmt(ixpt2(1)+1,iy,jj+1) - rmt(ixpt2(1)+1,iy,jj)
              dz0 = zmt(ixpt2(1)+1,iy,jj+1) - zmt(ixpt2(1)+1,iy,jj)
              if (abs(dr0) .gt. abs(dz0)) then  # avoid dr0 --> 0
                drs = (frac/4) * (
     .              2*(rmt(ixpt2(1)+1,iy  ,jj+1) - rmt(ixpt2(1)+1,iy  ,jj))+
     .                (rmt(ixpt2(1)+1,iy+1,jj+1) - rmt(ixpt2(1)+1,iy+1,jj))+
     .                (rmt(ixpt2(1)+1,iy-1,jj+1) - rmt(ixpt2(1)+1,iy-1,jj)) ) 
                rm(ixpt2n+1+ii,iy,jj+1) = rmt(ixpt2(1)+2,iy,jj)- drs
                zm(ixpt2n+1+ii,iy,jj+1)=zmt(ixpt2(1)+2,iy,jj)- dz0*drs/dr0
              else                              # avoid dz0 --> 0 
                dzs = (frac/4) * (
     .              2*(zmt(ixpt2(1)+1,iy  ,jj+1) - zmt(ixpt2(1)+1,iy  ,jj))+
     .                (zmt(ixpt2(1)+1,iy+1,jj+1) - zmt(ixpt2(1)+1,iy+1,jj))+
     .                (zmt(ixpt2(1)+1,iy-1,jj+1) - zmt(ixpt2(1)+1,iy-1,jj)) ) 
                zm(ixpt2n+1+ii,iy,jj+1) = zmt(ixpt2(1)+2,iy,jj)- dzs
                rm(ixpt2n+1+ii,iy,jj+1)=rmt(ixpt2(1)+2,iy,jj)- dr0*dzs/dz0
              endif
            endif
          enddo
        enddo
      enddo

      do ii = 1, nxxpt  # 0, nxxpt-1
        do iy = 1, nym
          do jj = 1, 3, 2
            rm(ixpt2n+1+ii,iy,jj) = rm(ixpt2n+ii,iy,jj+1)
            zm(ixpt2n+1+ii,iy,jj) = zm(ixpt2n+ii,iy,jj+1)
          enddo
        enddo
      enddo

      do iy = 1, nym
        do jj = 1, 3, 2
          rm(ixpt2n+1,iy,jj) = rmt(ixpt2(1)+1,iy,jj)
          zm(ixpt2n+1,iy,jj) = zmt(ixpt2(1)+1,iy,jj)
        enddo
      enddo

      # Correct inside vertices at iy=1 and iy=nym to fit interior cells
      do ii = 0, nxxpt
         do jj = 3, 4
            rm(ixpt2n+1+ii,1,jj) = rm(ixpt2n+1+ii,2,jj-2)
            zm(ixpt2n+1+ii,1,jj) = zm(ixpt2n+1+ii,2,jj-2)
            rm(ixpt2n+1+ii,nym,jj-2) = rm(ixpt2n+1+ii,nym-1,jj)
            zm(ixpt2n+1+ii,nym,jj-2) = zm(ixpt2n+1+ii,nym-1,jj)
         enddo
      enddo

c ... Now reset the unchanged cells by shifting indices

      do jj = 1, 4
        do iy = 1, nym
          do ix = 1, ixpt1(1)-1
            rm(ix,iy,jj) = rmt(ix,iy,jj)
            zm(ix,iy,jj) = zmt(ix,iy,jj)
          enddo
          do ix = ixpt1(1)+2, ixpt2(1)-1 
            rm(ix+2*nxxpt,iy,jj) = rmt(ix,iy,jj)
            zm(ix+2*nxxpt,iy,jj) = zmt(ix,iy,jj)
          enddo
          do ix = ixpt2(1)+2, nxm-4*nxxpt
            rm(ix+4*nxxpt,iy,jj) = rmt(ix,iy,jj)
            zm(ix+4*nxxpt,iy,jj) = zmt(ix,iy,jj)
          enddo
        enddo
      enddo

c ... Finally, reset rm(,,0) and zm(,,0), and also ixpt1,2
               
      do iy = 1, nym
        do ix = 1, nxm
          rm(ix,iy,0) = 0.25*( rm(ix,iy,1)+rm(ix,iy,2)+rm(ix,iy,3)+
     .                                                 rm(ix,iy,4) )
          zm(ix,iy,0) = 0.25*( zm(ix,iy,1)+zm(ix,iy,2)+zm(ix,iy,3)+
     .                                                 zm(ix,iy,4) )
        enddo
      enddo

      ixpt1(1) = ixpt1n
      ixtop    = ixtop + 2*nxxpt
      ixpt2(1) = ixpt2n

      return
      end
c ***** End of subroutine add_xptg *******************************

c----------------------------------------------------------------------c

      subroutine refine_xpt
      implicit none
Use(Dim)		# nym
Use(Dimflxgrd)	#npts
Use(Dimensions)		#
Use(Refinex)		#
Use(RZ_grid_info)	#
Use(Share)              # nxxpt

c     This is the driver subroutine for refining the mesh around
c     the x-point, using a more accurate flux surface representation.

c     Allocate space for working arrays --
      call gallot("Refinex",0)

c     Do the refinement --
      call refinexm

      return
      end

c----------------------------------------------------------------------c

      subroutine refinexm
      implicit none
Use(Dim)		# nym
Use(Dimflxgrd)	#npts
Use(Dimensions)		#
Use(Refinex)		# rsu,zsu,rsx,zsx,rflux,zflux,dsflux,nflux,alfxpt,
                        # nsmoothx,nxmod
Use(RZ_grid_info)	#
Use(Share)              # nxxpt,ishalfm
Use(UEgrid)             # ixtop
Use(Xpoint_indices)     # ixpt1,ixpt2

c     Refine the mesh near the x-point

c     local variables --

      integer ix, iy, iv, ivx, ivu, quadrant, idum, ierr, ibeg, iend
      integer i, ii, n, ixptold, ixptnew, idir, nrcell, k, qbegin
      real frac, ds, ff, ds_iyo, ds_iym, ds_min, frac2

c     Store original mesh in rmt,zmt
      do iv = 0, 4
        do iy = 1, nym
          do ix = 1, nxm-4*nxxpt
            rmt(ix,iy,iv) = rm(ix,iy,iv)
            zmt(ix,iy,iv) = zm(ix,iy,iv)
          enddo
        enddo
      enddo

c     NOTE: Quadrants are numbered clockwise about the x-point,
c     starting with number 1 for the lower-left quadrant.

      qbegin = 1
      if (ishalfm==1) qbegin = 3
      do quadrant=qbegin,4

c     Set characteristic indices for the xpoint location and
c     the 'upstream' poloidal direction
c     NOTE: The ixpt1 and ixpt2 below refer to the ORIGINAL mesh
         if (quadrant==1) then
            ixptold = ixpt1(1)
            ixptnew = ixptold+nxxpt
            idir = -1
            ivx  =  2
            ivu  =  1
            alfxpt = alfxptl
            alfxpt2 = alfxpt2l
         elseif (quadrant==2) then
            ixptold = ixpt1(1)+1
            ixptnew = ixptold+nxxpt
            idir =  1
            ivx  =  1
            ivu  =  2
            alfxpt = alfxptu
            alfxpt2 = alfxpt2u
         elseif (quadrant==3) then
            ixptold = ixpt2(1)
            ixptnew = ixptold+3*nxxpt
            idir = -1
            ivx  =  2
            ivu  =  1
            alfxpt = alfxptu
            alfxpt2 = alfxpt2u
         elseif (quadrant==4) then
            ixptold = ixpt2(1)+1
            ixptnew =  ixptold+3*nxxpt
            idir =  1
            ivx  =  1
            ivu  =  2
            alfxpt = alfxptl
            alfxpt2 = alfxpt2l
         endif

c     get angle-like 'upstream' (rsu,zsu) and 'x-point' (rsx,zsx) ref surfaces 
         call rsurface(quadrant)

c     Calulate the distance between ref. surfaces at iy=0 and iy=nym to use as
c     an alternate distance measure when constructing rmm,zmm below
         call calc_dsflux(quadrant,0,ibeg,iend)
           ds_iyo = dsflux(iend)
         call calc_dsflux(quadrant,nym,ibeg,iend)
           ds_iym = dsflux(iend)
         ds_min = min(ds_iyo, ds_iym)


         do iy=0,nym

c     Calc. distance, dsflux, along a flux surface between ref. surfaces
c     Computes common-block variables dsflux, rflux, and zflux
            call calc_dsflux(quadrant,iy,ibeg,iend)

c     interpolate along rflux,zflux --
c     Start with the point nearest (but not on) the 'x-point' reference
c     surface and move poloidally 'upstream'
            nrcell = nxmod + nxxpt
            ii=nrcell        # index of modified surface(s) in this quadrant
            do n = 1, nrcell-1 
               ii = ii - 1
               frac = (float(n)/float(nrcell))**alfxpt
               frac2 = (float(n)/float(nrcell-1))**alfxpt2
               ds = frac * (frac2*dsflux(iend) + (1-frac2)*ds_min)
               do i = ibeg, iend 
                  if ( (dsflux(i) .lt. ds) .and.
     .                              (dsflux(i+1) .ge. ds) ) then
                     ff = (ds - dsflux(i))/(dsflux(i+1)-dsflux(i))
                     rmm(iy,ii) = rflux(i) + ff * (rflux(i+1)-rflux(i))
                     zmm(iy,ii) = zflux(i) + ff * (zflux(i+1)-zflux(i))
                     break
                  endif
               enddo	# end loop over flux surface data points
            enddo	# end loop over modified mesh points

         enddo		# end loop over flux surfaces

c     smooth the newly-constructed angle-like surfaces:
         do k=1, nsmoothx
            do ii=1, nrcell-1
               call smoothx(rmm(0:nym,ii),zmm(0:nym,ii),0,nym,0,nym,quadrant)
            enddo
         enddo

c     copy data for modified cells into rm,zm
c     First, the un-modified x-point reference surface:
         do iv=ivx,ivx+2,2
            do iy=1,nym
               ix=ixptnew
               rm(ix,iy,iv) = rmt(ixptold,iy,iv)
               zm(ix,iy,iv) = zmt(ixptold,iy,iv)
            enddo
         enddo

c     Next, the un-modified upstream reference surface:
         do iv=ivu,ivu+2,2
            do iy=1,nym
               ix=ixptnew+idir*(nrcell-1)
               rm(ix,iy,iv) = rmt(ixptold+idir*(nxmod-1),iy,iv)
               zm(ix,iy,iv) = zmt(ixptold+idir*(nxmod-1),iy,iv)
            enddo
         enddo

c     Finally, the modified surfaces:
         do ii=1,nrcell-1
            do iy=1,nym       # some redundant assignment here
               ix=ixptnew+ii*idir
               rm(ix,iy,ivx) = rmm(iy-1,ii)
               zm(ix,iy,ivx) = zmm(iy-1,ii)
               call setvadj(ix,iy,ivx,rmm(iy-1,ii),zmm(iy-1,ii))
               rm(ix,iy,ivx+2) = rmm(iy,ii)
               zm(ix,iy,ivx+2) = zmm(iy,ii)
               call setvadj(ix,iy,ivx+2,rmm(iy,ii),zmm(iy,ii))
            enddo

         enddo

      enddo		# end loop over quadrants

c     Now reset the unchanged cells by shifting indices

      do iv = 1, 4
        do iy = 1, nym
          do ix = 1, ixpt1(1)-nxmod
            rm(ix,iy,iv) = rmt(ix,iy,iv)
            zm(ix,iy,iv) = zmt(ix,iy,iv)
          enddo
          do ix = ixpt1(1)+1+nxmod, ixpt2(1)-nxmod 
            rm(ix+2*nxxpt,iy,iv) = rmt(ix,iy,iv)
            zm(ix+2*nxxpt,iy,iv) = zmt(ix,iy,iv)
          enddo
          do ix = ixpt2(1)+1+nxmod, nxm-4*nxxpt
            rm(ix+4*nxxpt,iy,iv) = rmt(ix,iy,iv)
            zm(ix+4*nxxpt,iy,iv) = zmt(ix,iy,iv)
          enddo
        enddo
      enddo

c     Re-compute cell-center coordinates
      do iy = 1, nym
        do ix = 1, nxm
          rm(ix,iy,0) = 0.25*( rm(ix,iy,1)+rm(ix,iy,2)+rm(ix,iy,3)+
     .                                                 rm(ix,iy,4) )
          zm(ix,iy,0) = 0.25*( zm(ix,iy,1)+zm(ix,iy,2)+zm(ix,iy,3)+
     .                                                 zm(ix,iy,4) )
        enddo
      enddo

c     Finally, reset x-point and top-of-mesh indices
      ixpt1(1) = ixpt1(1) +   nxxpt
      ixtop    = ixtop    + 2*nxxpt
      ixpt2(1) = ixpt2(1) + 3*nxxpt

      return
      end

c----------------------------------------------------------------------c

      subroutine calc_dsflux(quad,iy_idx,ibeg,iend)

      implicit none
Use(Dim)		# nym
Use(Dimflxgrd)	#npts
Use(Dimensions)		#
Use(Refinex)		# rsu,zsu,rsx,zsx,rflux,zflux,dsflux,nflux


c     Calulates distance along flux surfaces, dsflux, for refinex routine

c...  Input variables
      integer quad    # defines which mesh quadrant
      integer iy_idx  # flux surface index

c...  Output variables
      integer ibeg,iend   # beginning and ending indices of r,z,dsflux arrays

c     local variables --

      integer idum, ierr, i
      real rbeg, zbeg, rend, zend, fuzz

      data fuzz /1.0e-08/

c     get rflux and zflux data, with index increasing toward the x-point
      call fluxcurve(quad,iy_idx)

c     compute intersection with 'upstream' reference surface
      call intersect2(rsu,zsu,0,nym+2,rflux,zflux,1,nflux,
     .                                  rbeg,zbeg,idum,ibeg,fuzz,ierr)
c     use the intersection as the starting point on the flux surface:
      rflux(ibeg)=rbeg
      zflux(ibeg)=zbeg

c     compute intersection with 'x-point' reference surface
      call intersect2(rsx,zsx,0,nym+2,rflux,zflux,1,nflux,
     .                                 rend,zend,idum,iend,fuzz,ierr)
c     use the intersection as the end point on the flux surface:
      iend=iend+1
      rflux(iend)=rend
      zflux(iend)=zend

c     compute distance along r,zflux between r,zsu and r,zsx
      dsflux(ibeg) = 0.
      do i = ibeg+1, iend
         dsflux(i) = dsflux(i-1) +
     .         sqrt( (rflux(i)-rflux(i-1))**2+(zflux(i)-zflux(i-1))**2 )
      enddo

      return
      end

c----------------------------------------------------------------------c


      subroutine setvadj(ix,iy,iv,rs,zs)
      implicit none
      integer ix,iy,iv
      real rs,zs
Use(Dim)
Use(RZ_grid_info)

c     Set vertices of adjacent cells that share the same data as (ix,iy,iv).
c     NOTE:
c     This procedure is not valid along the cuts in the mesh or along
c     the outermost guard cell boundaries of the mesh (08 Mar 1999)

      if (iv==1) then
         rm(ix-1,iy,2) = rs
         zm(ix-1,iy,2) = zs
         rm(ix,iy-1,3) = rs
         zm(ix,iy-1,3) = zs
         rm(ix-1,iy-1,4) = rs
         zm(ix-1,iy-1,4) = zs
      elseif (iv==2) then
         rm(ix+1,iy,1) = rs
         zm(ix+1,iy,1) = zs
         rm(ix,iy-1,4) = rs
         zm(ix,iy-1,4) = zs
         rm(ix+1,iy-1,3) = rs
         zm(ix+1,iy-1,3) = zs
      elseif (iv==3) then
         rm(ix-1,iy,4) = rs
         zm(ix-1,iy,4) = zs
         rm(ix,iy+1,1) = rs
         zm(ix,iy+1,1) = zs
         rm(ix-1,iy+1,2) = rs
         zm(ix-1,iy+1,2) = zs
      elseif (iv==4) then
         rm(ix+1,iy,3) = rs
         zm(ix+1,iy,3) = zs
         rm(ix,iy+1,2) = rs
         zm(ix,iy+1,2) = zs
         rm(ix+1,iy+1,1) = rs
         zm(ix+1,iy+1,1) = zs
      endif

      return
      end

c----------------------------------------------------------------------c

      subroutine rsurface(quadrant)
      implicit none
      integer quadrant
Use(Dim)               # nym
Use(Dimflxgrd)	#npts
Use(Dimensions)        #
Use(Refinex)           # rsu,zsu
Use(Xpoint_indices)    # ixpt1,ixpt2
Use(RZ_grid_info)

c     local variables --
      integer iy

c     Define angle-like reference surface arrays (for each quadrant):

      if (quadrant==1) then

c     Upstream reference surface is WEST edge of IXPT1-nxmod+1:
         rsu(1) = rmt(ixpt1(1)-nxmod+1,1,1)
         zsu(1) = zmt(ixpt1(1)-nxmod+1,1,1)
         do iy=1, nym
            rsu(iy+1) = rmt(ixpt1(1)-nxmod+1,iy,3)
            zsu(iy+1) = zmt(ixpt1(1)-nxmod+1,iy,3)
         enddo
c     Get endpoints by extrapolation:
         rsu(0)    = 2*rsu(1)    - rsu(2)
         zsu(0)    = 2*zsu(1)    - zsu(2)
         rsu(nym+2) = 2*rsu(nym+1) - rsu(nym)
         zsu(nym+2) = 2*zsu(nym+1) - zsu(nym)

c     X-point reference surface is EAST edge of IXPT1:
         rsx(1) = rmt(ixpt1(1),1,2)
         zsx(1) = zmt(ixpt1(1),1,2)
         do iy=1, nym
            rsx(iy+1) = rmt(ixpt1(1),iy,4)
            zsx(iy+1) = zmt(ixpt1(1),iy,4)
         enddo
c     Get endpoints by extrapolation:
         rsx(0)    = 2*rsx(1)    - rsx(2)
         zsx(0)    = 2*zsx(1)    - zsx(2)
         rsx(nym+2) = 2*rsx(nym+1) - rsx(nym)
         zsx(nym+2) = 2*zsx(nym+1) - zsx(nym)

      elseif (quadrant==2) then

c     Upstream reference surface is EAST edge of IXPT1+nxmod:
         rsu(1) = rmt(ixpt1(1)+nxmod,1,2)
         zsu(1) = zmt(ixpt1(1)+nxmod,1,2)
         do iy=1, nym
            rsu(iy+1) = rmt(ixpt1(1)+nxmod,iy,4)
            zsu(iy+1) = zmt(ixpt1(1)+nxmod,iy,4)
         enddo
c     Get endpoints by extrapolation:
         rsu(0)    = 2*rsu(1)    - rsu(2)
         zsu(0)    = 2*zsu(1)    - zsu(2)
         rsu(nym+2) = 2*rsu(nym+1) - rsu(nym)
         zsu(nym+2) = 2*zsu(nym+1) - zsu(nym)


c     X-point reference surface is WEST edge of IXPT1+1:
         rsx(1) = rmt(ixpt1(1)+1,1,1)
         zsx(1) = zmt(ixpt1(1)+1,1,1)
         do iy=1, nym
            rsx(iy+1) = rmt(ixpt1(1)+1,iy,3)
            zsx(iy+1) = zmt(ixpt1(1)+1,iy,3)
         enddo
c     Get endpoints by extrapolation:
         rsx(0)    = 2*rsx(1)    - rsx(2)
         zsx(0)    = 2*zsx(1)    - zsx(2)
         rsx(nym+2) = 2*rsx(nym+1) - rsx(nym)
         zsx(nym+2) = 2*zsx(nym+1) - zsx(nym)

      elseif (quadrant==3) then

c     Upstream reference surface is WEST edge of IXPT2-nxmod+1:
         rsu(1) = rmt(ixpt2(1)-nxmod+1,1,1)
         zsu(1) = zmt(ixpt2(1)-nxmod+1,1,1)
         do iy=1, nym
            rsu(iy+1) = rmt(ixpt2(1)-nxmod+1,iy,3)
            zsu(iy+1) = zmt(ixpt2(1)-nxmod+1,iy,3)
         enddo
c     Get endpoints by extrapolation:
         rsu(0)    = 2*rsu(1)    - rsu(2)
         zsu(0)    = 2*zsu(1)    - zsu(2)
         rsu(nym+2) = 2*rsu(nym+1) - rsu(nym)
         zsu(nym+2) = 2*zsu(nym+1) - zsu(nym)

c     X-point reference surface is EAST edge of IXPT2:
         rsx(1) = rmt(ixpt2(1),1,2)
         zsx(1) = zmt(ixpt2(1),1,2)
         do iy=1, nym
            rsx(iy+1) = rmt(ixpt2(1),iy,4)
            zsx(iy+1) = zmt(ixpt2(1),iy,4)
         enddo
c     Get endpoints by extrapolation:
         rsx(0)    = 2*rsx(1)    - rsx(2)
         zsx(0)    = 2*zsx(1)    - zsx(2)
         rsx(nym+2) = 2*rsx(nym+1) - rsx(nym)
         zsx(nym+2) = 2*zsx(nym+1) - zsx(nym)

      elseif (quadrant==4) then

c     Upstream reference surface is EAST edge of IXPT2+nxmod:
         rsu(1) = rmt(ixpt2(1)+nxmod,1,2)
         zsu(1) = zmt(ixpt2(1)+nxmod,1,2)
         do iy=1, nym
            rsu(iy+1) = rmt(ixpt2(1)+nxmod,iy,4)
            zsu(iy+1) = zmt(ixpt2(1)+nxmod,iy,4)
         enddo
c     Get endpoints by extrapolation:
         rsu(0)    = 2*rsu(1)    - rsu(2)
         zsu(0)    = 2*zsu(1)    - zsu(2)
         rsu(nym+2) = 2*rsu(nym+1) - rsu(nym)
         zsu(nym+2) = 2*zsu(nym+1) - zsu(nym)

c     X-point reference surface is WEST edge of IXPT2+1:
         rsx(1) = rmt(ixpt2(1)+1,1,1)
         zsx(1) = zmt(ixpt2(1)+1,1,1)
         do iy=1, nym
            rsx(iy+1) = rmt(ixpt2(1)+1,iy,3)
            zsx(iy+1) = zmt(ixpt2(1)+1,iy,3)
         enddo
c     Get endpoints by extrapolation:
         rsx(0)    = 2*rsx(1)    - rsx(2)
         zsx(0)    = 2*zsx(1)    - zsx(2)
         rsx(nym+2) = 2*rsx(nym+1) - rsx(nym)
         zsx(nym+2) = 2*zsx(nym+1) - zsx(nym)

      endif

      return
      end

c----------------------------------------------------------------------c

      subroutine fluxcurve(quadrant,iy)
      implicit none
      integer iy, quadrant
Use(Dim)         # nym
Use(Dimflxgrd)	 #jdim,npts
Use(Comflxgrd)   #jsptrx
Use(Dimensions)  #
Use(Curves)      # xcurveg, ycurveg, npointg
Use(Refinex)
## Use(Inmesh)      # jsptrxg
Use(Transfm)     # ijump

c     Copy the data for a single flux surface to work arrays r,zflux.
c     Order the data with increasing index towards the x-point.

c     local variables --
      integer j, n, ibeg, iend

c     NOTE: Quadrants are numbered clockwise about the x-point,
c     starting with number 1 for the lower-left quadrant.

      if (quadrant == 1) then		# lower-left quadrant
         j=nym-iy+1
         if (j > jsptrx(1)) then        # pf region
            ibeg=npointg(j)
            iend=ijump(j)+1
         else                           # SOL region
            ibeg=npointg(j)
            iend=1
         endif
         nflux=ibeg-iend+1
         do n=1,nflux
            rflux(n)=xcurveg(ibeg-n+1,j)
            zflux(n)=ycurveg(ibeg-n+1,j)
         enddo
      elseif (quadrant == 2) then	# upper-left quadrant
         j=nym-iy+1
	 if (j > jsptrx(1)) then       # core region
            ibeg=1
            iend=ijump(j)
         else                           # SOL region
            ibeg=1
            iend=npointg(j)
         endif
         nflux=iend-ibeg+1
         do n=1,nflux
            rflux(n)=xcurveg(ibeg+n-1,j)
            zflux(n)=ycurveg(ibeg+n-1,j)
         enddo
      elseif (quadrant == 3) then	# upper-right quadrant
         j=nym+iy+3
	 if (j < jsptrx(2)) then       # core region
            ibeg=1
            iend=ijump(j)
         else                           # SOL region
            ibeg=1
            iend=npointg(j)
         endif
         nflux=iend-ibeg+1
         do n=1,nflux
            rflux(n)=xcurveg(ibeg+n-1,j)
            zflux(n)=ycurveg(ibeg+n-1,j)
         enddo
      elseif (quadrant == 4) then	# lower-right quadrant
         j=nym+iy+3
         if (j < jsptrx(2)) then        # pf region
            ibeg=npointg(j)
            iend=ijump(j)+1
         else                           # SOL region
            ibeg=npointg(j)
            iend=1
         endif
         nflux=ibeg-iend+1
         do n=1,nflux
            rflux(n)=xcurveg(ibeg-n+1,j)
            zflux(n)=ycurveg(ibeg-n+1,j)
         enddo
      endif

      return
      end

c----------------------------------------------------------------------c

      subroutine smoothx(rmm,zmm,nd1,nd2,iy1,iy2,quadrant)
      implicit none
      integer nd1,nd2,iy1,iy2,quadrant
      real rmm(nd1:nd2),zmm(nd1:nd2)
Use(Dimflxgrd)	#jdim,npts,noregs
Use(Comflxgrd)  # jmin,jmax,jsptrx
Use(Dimensions)
Use(Curves)	# npointg,xcurveg,ycurveg
##Use(Inmesh)	# jmin,jmax,jsptrx
Use(Mmod)	# nupstream,rupstream,zupstream as temporary arrays
                # ndnstream,rdnstream,zdnstream as temporary arrays
		# xcrv,ycrv,dsc,fuzzm,wtold
Use(Transfm)	# ijump
      external intersect2, remark, xerrab, twixt

c     Local variables --
      integer ii,j,i1crv,i2crv,icu,kcu,icp,kcp,ierr,istart,iend,iy
      integer region
      real xcu,ycu,xcp,ycp,xjm1,yjm1,distnew,frac

c     Smooth spatial irregularities in angle-like surfaces of the mesh
c     while ensuring that the smoothed points lie exactly on flux surfaces.

      xjm1=rmm(iy1)
      yjm1=zmm(iy1)

      do iy=iy1+1,iy2-1

c     Relation between iy and j flux surface indices and region is:
         if ((quadrant==1) .or. (quadrant==2)) then
            region=1
            j=jmax(1)-iy
         elseif ((quadrant==3) .or. (quadrant==4)) then
            region=2
            j=jmin(2)+iy
         endif

c     Temporary variables for (xcurveg,ycurveg) surface data
         if ((quadrant==1) .and. (j > jsptrx(region))) then      # pf surface
            i1crv=1+ijump(j)
            i2crv=npointg(j)
         elseif ((quadrant==2) .and. (j > jsptrx(region))) then  # core surface
            i1crv=1
            i2crv=ijump(j)
         elseif ((quadrant==3) .and. (j < jsptrx(region))) then  # core surface
            i1crv=1
            i2crv=ijump(j)
         elseif ((quadrant==4) .and. (j < jsptrx(region))) then  # pf surface
            i1crv=1+ijump(j)
            i2crv=npointg(j)
         else                                                     # SOL surface
            i1crv=1
            i2crv=npointg(j)
         endif
###      i1crv=1		# start index
###      i2crv=npointg(j)	# end index
         do ii=i1crv,i2crv
            xcrv(ii)=xcurveg(ii,j)
            ycrv(ii)=ycurveg(ii,j)
         enddo
   
c     Straight line between iy-1 and iy+1 points on angle-like surface
	 nupstream=2
	 rupstream(1)=rmm(iy-1)
	 zupstream(1)=zmm(iy-1)
	 rupstream(2)=rmm(iy+1)
	 zupstream(2)=zmm(iy+1)

c     Intersection of straight line with x,ycurveg is (xcu,ycu):
         call intersect2(rupstream,zupstream,1,nupstream,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xcu,ycu,kcu,icu,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,886) j,iy
 886        format 
     & ("smoothx: no intersec'n of j=",i3," with str line for iy=",i3)
            call xerrab("")
         endif
c     The point (x,ycu) lies between icu and icu+1 on x,ycrv.

c     Segmented line through iy-1, iy and iy+1 points on angle-like surface
         ndnstream=3
         rdnstream(1)=rmm(iy-1)
         zdnstream(1)=zmm(iy-1)
         rdnstream(2)=rmm(iy)
         zdnstream(2)=zmm(iy)
         rdnstream(3)=rmm(iy+1)
         zdnstream(3)=zmm(iy+1)

c     Intersection of segmented line with x,ycurveg is (xcp,ycp):
         call intersect2(rdnstream,zdnstream,1,ndnstream,
     &                   xcrv(i1crv:i2crv),ycrv(i1crv:i2crv),i1crv,i2crv,
     &                   xcp,ycp,kcp,icp,fuzzm,ierr)
         if (ierr .ne. 0) then
            write (STDOUT,887) j,iy
 887        format 
     & ("smoothx: no intersec'n of j=",i3," with seg line for iy=",i3)
            call xerrab("")
         endif
c     The point (x,ycp) lies between icp and icp+1 on x,ycrv.

c     Update meshpoint at iy-1, now that we are done using it:
         rmm(iy-1)=xjm1
         zmm(iy-1)=yjm1

         if (icp .lt. icu) then		# x,ycp is  upstream  from x,ycu
            istart=icp
            xcrv(istart)=xcp
            ycrv(istart)=ycp
            iend=icu+1
            xcrv(iend)=xcu
            ycrv(iend)=ycu
            dsc(istart)=0.
            do ii=istart+1,iend
               dsc(ii)=dsc(ii-1)
     &                     + sqrt( (xcrv(ii)-xcrv(ii-1))**2
     &                            +(ycrv(ii)-ycrv(ii-1))**2 )
            enddo
            distnew=(1-wtold)*dsc(iend)
            do ii=istart+1,iend
               if ( distnew .le. dsc(ii) ) then
c        new meshpoint lies between ii-1 and ii
                  frac = (distnew-dsc(ii-1))
     &                         /(dsc(ii)-dsc(ii-1))
                  xjm1=xcrv(ii-1) + frac*(xcrv(ii)-xcrv(ii-1))
                  yjm1=ycrv(ii-1) + frac*(ycrv(ii)-ycrv(ii-1))
                  break
               endif
            enddo
         elseif (icp .gt. icu) then	# x,ycp is downstream from x,ycu
            istart=icu
            xcrv(istart)=xcu
            ycrv(istart)=ycu
            iend=icp+1
            xcrv(iend)=xcp
            ycrv(iend)=ycp
            dsc(istart)=0.
            do ii=istart+1,iend
               dsc(ii)=dsc(ii-1)
     &                     + sqrt( (xcrv(ii)-xcrv(ii-1))**2
     &                            +(ycrv(ii)-ycrv(ii-1))**2 )
            enddo
            distnew=wtold*dsc(iend)
            do ii=istart+1,iend
               if (distnew .le. dsc(ii)) then
c        new meshpoint lies between ii-1 and ii
                  frac = (distnew-dsc(ii-1))
     &                         /(dsc(ii)-dsc(ii-1))
                  xjm1=xcrv(ii-1) + frac*(xcrv(ii)-xcrv(ii-1))
                  yjm1=ycrv(ii-1) + frac*(ycrv(ii)-ycrv(ii-1))
                  break
               endif
            enddo
         else	# both x,ycu and x,ycp lie on same segment of x,ycrv
            xjm1=(1-wtold)*xcu+wtold*xcp
            yjm1=(1-wtold)*ycu+wtold*ycp
         endif

      enddo

      rmm(iy2-1)=xjm1
      zmm(iy2-1)=yjm1

      return
      end

c----------------------------------------------------------------------c

      subroutine mapdnbot2dnull(ixpt1t,ixtopt,ixpt2t,nxmt,iysptrxt,nymt)
      implicit none
      integer ixpt1t, ixtopt, ixpt2t, nxmt, iysptrxt, nymt
Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm,rmt,zmt
Use(Xpoint_indices) # ixlb,ixpt1,ixmdp,ixpt2,ixrb,iysptrx1,iysptrx2
Use(UEgrid)         # ixtop
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)     #
      integer ix, iy, iv, ivt, nxshift

c     Store half-space data in temporary arrays:
      do iv=0,4
         call s2copy(nxm+2, nym+2, rm(0:nxm+1,0:nym+1,iv), 1, nxm+2, 
     .                              rmt(0:nxm+1,0:nym+1,iv), 1, nxm+2)
         call s2copy(nxm+2, nym+2, zm(0:nxm+1,0:nym+1,iv), 1, nxm+2, 
     .                              zmt(0:nxm+1,0:nym+1,iv), 1, nxm+2)
      enddo

c     The left half of dnbot maps directly to the lower-left quadrant
c     of the dnull configuration:
      ixlb(1)=0
      ixpt1(1)=ixpt1t
      ixmdp(1)=ixtopt-1
      do ix=ixlb(1),ixmdp(1)
         do iy=1,nym
            do iv=0,4
               rm(ix,iy,iv)=rmt(ix,iy,iv)
               zm(ix,iy,iv)=zmt(ix,iy,iv)
            enddo
         enddo
      enddo

c     The right half of dnbot maps to the lower-right quadrant
c     of the dnull configuration by shifting the ix index:
      nxshift  = nxm - nxmt
      ixrb(2)  = nxmt + nxshift
      ixpt2(2) = ixpt2t + nxshift
      ixmdp(2) = ixtopt+1 + nxshift
      do ix=ixmdp(2)+1,ixrb(2)+1
         do iy=1,nym
            do iv=0,4
               rm(ix,iy,iv)=rmt(ix-nxshift,iy,iv)
               zm(ix,iy,iv)=zmt(ix-nxshift,iy,iv)
            enddo
         enddo
      enddo

c     NOTE: the iv indices map as follows:
c        1(lower)->2(upper)
c        3(lower)->4(upper)
c        2(lower)->1(upper)
c        4(lower)->3(upper)

c     The left half of dnbot maps to the upper-left quadrant
c     of the dnull configuration using the relation:
c              ix_upper = ixmdp(1) + ixtopt - ix_lower
      ixpt2(1) = ixmdp(1) + ixtopt - ixpt1t - 1
      ixrb(1)  = ixpt2(1) + ixpt1t
      do ix=ixmdp(1)+1,ixrb(1)+1
         do iy=1,nym
            do iv=0,4
               if (iv==0) ivt=0
               if (iv==1) ivt=2
               if (iv==2) ivt=1
               if (iv==3) ivt=4
               if (iv==4) ivt=3
               rm(ix,iy,iv)=rmt(ixmdp(1)+ixtopt-ix,iy,ivt)
               zm(ix,iy,iv)=2*zmagx - zmt(ixmdp(1)+ixtopt-ix,iy,ivt)
            enddo
         enddo
      enddo

c     The right half of dnbot maps to the upper-right quadrant
c     of the dnull configuration using the relation:
c              ix_upper = ixmdp(2) + ixtopt - ix_lower + 2
      ixlb(2)  = ixmdp(2) + ixtopt - (nxmt+1) + 2
      ixpt1(2) = ixmdp(2) + ixtopt - (ixpt2t+1) + 2
      do ix=ixlb(2),ixmdp(2)
         do iy=1,nym
            do iv=0,4
               if (iv==0) ivt=0
               if (iv==1) ivt=2
               if (iv==2) ivt=1
               if (iv==3) ivt=4
               if (iv==4) ivt=3
               rm(ix,iy,iv)=rmt(ixmdp(2)+ixtopt-ix+2,iy,ivt)
               zm(ix,iy,iv)=2*zmagx - zmt(ixmdp(2)+ixtopt-ix+2,iy,ivt)
            enddo
         enddo
      enddo

c     Set x-point locations in iy:
      iysptrx1(1) = iysptrxt
      iysptrx2(1) = iysptrxt
      iysptrx1(2) = iysptrxt
      iysptrx2(2) = iysptrxt

      return
      end

c     -------------------------------------------------------------------

      subroutine add_guardc_tp
      implicit none
Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm
Use(Xpoint_indices) # ixlb,ixrb
Use(Share)          # epslon
      integer ix,iy

c     Create guard cells at target plates:
      do iy=1,nym
         ix=ixlb(1)      # at the lower left target plate
         rm(ix,iy,1)=rm(ix+1,iy,1)+epslon*(rm(ix+1,iy,1)-rm(ix+1,iy,2))
         zm(ix,iy,1)=zm(ix+1,iy,1)+epslon*(zm(ix+1,iy,1)-zm(ix+1,iy,2))
         rm(ix,iy,2)=rm(ix+1,iy,1)
         zm(ix,iy,2)=zm(ix+1,iy,1)
         rm(ix,iy,3)=rm(ix+1,iy,3)+epslon*(rm(ix+1,iy,3)-rm(ix+1,iy,4))
         zm(ix,iy,3)=zm(ix+1,iy,3)+epslon*(zm(ix+1,iy,3)-zm(ix+1,iy,4))
         rm(ix,iy,4)=rm(ix+1,iy,3)
         zm(ix,iy,4)=zm(ix+1,iy,3)
         rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                  +rm(ix,iy,3)+rm(ix,iy,4))
         zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                  +zm(ix,iy,3)+zm(ix,iy,4))
         ix=ixrb(1)+1    # at the upper left target plate
         rm(ix,iy,1)=rm(ix-1,iy,2)
         zm(ix,iy,1)=zm(ix-1,iy,2)
         rm(ix,iy,2)=rm(ix-1,iy,2)+epslon*(rm(ix-1,iy,2)-rm(ix-1,iy,1))
         zm(ix,iy,2)=zm(ix-1,iy,2)+epslon*(zm(ix-1,iy,2)-zm(ix-1,iy,1))
         rm(ix,iy,3)=rm(ix-1,iy,4)
         zm(ix,iy,3)=zm(ix-1,iy,4)
         rm(ix,iy,4)=rm(ix-1,iy,4)+epslon*(rm(ix-1,iy,4)-rm(ix-1,iy,3))
         zm(ix,iy,4)=zm(ix-1,iy,4)+epslon*(zm(ix-1,iy,4)-zm(ix-1,iy,3))
         rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                  +rm(ix,iy,3)+rm(ix,iy,4))
         zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                  +zm(ix,iy,3)+zm(ix,iy,4))
         ix=ixlb(2)      # at the upper right target plate
         rm(ix,iy,1)=rm(ix+1,iy,1)+epslon*(rm(ix+1,iy,1)-rm(ix+1,iy,2))
         zm(ix,iy,1)=zm(ix+1,iy,1)+epslon*(zm(ix+1,iy,1)-zm(ix+1,iy,2))
         rm(ix,iy,2)=rm(ix+1,iy,1)
         zm(ix,iy,2)=zm(ix+1,iy,1)
         rm(ix,iy,3)=rm(ix+1,iy,3)+epslon*(rm(ix+1,iy,3)-rm(ix+1,iy,4))
         zm(ix,iy,3)=zm(ix+1,iy,3)+epslon*(zm(ix+1,iy,3)-zm(ix+1,iy,4))
         rm(ix,iy,4)=rm(ix+1,iy,3)
         zm(ix,iy,4)=zm(ix+1,iy,3)
         rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                  +rm(ix,iy,3)+rm(ix,iy,4))
         zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                  +zm(ix,iy,3)+zm(ix,iy,4))
         ix=ixrb(2)+1    # at the lower right target plate
         rm(ix,iy,1)=rm(ix-1,iy,2)
         zm(ix,iy,1)=zm(ix-1,iy,2)
         rm(ix,iy,2)=rm(ix-1,iy,2)+epslon*(rm(ix-1,iy,2)-rm(ix-1,iy,1))
         zm(ix,iy,2)=zm(ix-1,iy,2)+epslon*(zm(ix-1,iy,2)-zm(ix-1,iy,1))
         rm(ix,iy,3)=rm(ix-1,iy,4)
         zm(ix,iy,3)=zm(ix-1,iy,4)
         rm(ix,iy,4)=rm(ix-1,iy,4)+epslon*(rm(ix-1,iy,4)-rm(ix-1,iy,3))
         zm(ix,iy,4)=zm(ix-1,iy,4)+epslon*(zm(ix-1,iy,4)-zm(ix-1,iy,3))
         rm(ix,iy,0) = 0.25*(rm(ix,iy,1)+rm(ix,iy,2)
     &                                  +rm(ix,iy,3)+rm(ix,iy,4))
         zm(ix,iy,0) = 0.25*(zm(ix,iy,1)+zm(ix,iy,2)
     &                                  +zm(ix,iy,3)+zm(ix,iy,4))
      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine symmetrize_magnetics
      implicit none
Use(Dim)            # nxm,nym
Use(Xpoint_indices) # ixlb,ixpt1,ixmdp,ixpt2,ixrb
Use(RZ_grid_info)   # psi,br,bz,bphi,bpol,b

      integer ix, iy, iv, ivt

c     Symmetrize the double-null magnetic configuration by mapping
c     lower-half values to the upper-half.

c     NOTE: the iv indices map as follows:
c        1(lower)->2(upper)
c        3(lower)->4(upper)
c        2(lower)->1(upper)
c        4(lower)->3(upper)

c     The upper-left quadrant maps from the lower-left quadrant
c     via the relation:
c              ix_upper = 2*ixmdp(1) - ix_lower + 1
      do ix=ixmdp(1)+1,ixrb(1)+1
         do iy=1,nym
            do iv=0,4
               if (iv==0) ivt=0
               if (iv==1) ivt=2
               if (iv==2) ivt=1
               if (iv==3) ivt=4
               if (iv==4) ivt=3
               psi(ix,iy,iv)  =  psi(2*ixmdp(1)-ix+1,iy,ivt)
               br(ix,iy,iv)   = -br(2*ixmdp(1)-ix+1,iy,ivt)
               bz(ix,iy,iv)   =  bz(2*ixmdp(1)-ix+1,iy,ivt)
               bphi(ix,iy,iv) =  bphi(2*ixmdp(1)-ix+1,iy,ivt)
               bpol(ix,iy,iv) =  bpol(2*ixmdp(1)-ix+1,iy,ivt)
               b(ix,iy,iv)    =  b(2*ixmdp(1)-ix+1,iy,ivt)
            enddo
         enddo
      enddo

c     The upper-right quadrant maps from the lower-right quadrant
c     via the relation:
c              ix_upper = 2*ixmdp(2) - ix_lower + 1
      do ix=ixlb(2),ixmdp(2)
         do iy=1,nym
            do iv=0,4
               if (iv==0) ivt=0
               if (iv==1) ivt=2
               if (iv==2) ivt=1
               if (iv==3) ivt=4
               if (iv==4) ivt=3
               psi(ix,iy,iv)  =  psi(2*ixmdp(2)-ix+1,iy,ivt)
               br(ix,iy,iv)   = -br(2*ixmdp(2)-ix+1,iy,ivt)
               bz(ix,iy,iv)   =  bz(2*ixmdp(2)-ix+1,iy,ivt)
               bphi(ix,iy,iv) =  bphi(2*ixmdp(2)-ix+1,iy,ivt)
               bpol(ix,iy,iv) =  bpol(2*ixmdp(2)-ix+1,iy,ivt)
               b(ix,iy,iv)    =  b(2*ixmdp(2)-ix+1,iy,ivt)
            enddo
         enddo
      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine mapdnbot
      implicit none
Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm
Use(Xpoint_indices) # ixlb,ixpt1,ixmdp,ixpt2,ixrb
Use(Dnull_temp)     # ixpt1b,ixtopb,ixpt2b,nxmb,nymb,rmb,zmb
      integer ix, iy, iv, nxshift

c     Map bottom-half mesh data (from dnbot) onto full double-null mesh

c     The left half of dnbot maps directly to the lower-left quadrant
c     of the dnull configuration:
      ixlb(1)=0
      ixpt1(1)=ixpt1b
      ixmdp(1)=ixtopb-1
      do ix=ixlb(1),ixmdp(1)
         do iy=1,nym
            do iv=0,4
               rm(ix,iy,iv)=rmb(ix,iy,iv)
               zm(ix,iy,iv)=zmb(ix,iy,iv)
            enddo
         enddo
      enddo

c     The right half of dnbot maps to the lower-right quadrant
c     of the dnull configuration by shifting the ix index:
      nxshift  = nxm - nxmb
      ixrb(2)  = nxmb + nxshift
      ixpt2(2) = ixpt2b + nxshift
      ixmdp(2) = ixtopb+1 + nxshift
      do ix=ixmdp(2)+1,ixrb(2)+1
         do iy=1,nym
            do iv=0,4
               rm(ix,iy,iv)=rmb(ix-nxshift,iy,iv)
               zm(ix,iy,iv)=zmb(ix-nxshift,iy,iv)
            enddo
         enddo
      enddo

      return
      end

c----------------------------------------------------------------------c

      subroutine mapdntop
      implicit none
Use(Dim)            # nxm,nym
Use(RZ_grid_info)   # rm,zm
Use(Xpoint_indices) # ixlb,ixpt1,ixmdp,ixpt2,ixrb
Use(Dnull_temp)     # ixpt1u,ixtopu,ixpt2u,nxmu,nymu,rmu,zmu
Use(Dimflxgrd)	#nlim,nxefit,nyefit
Use(Comflxgrd)
Use(Dimensions)     #
      integer ix, iy, iv, ivt, nxshift

c     Map upper-half mesh data (from dntop) onto full double-null mesh
c     NOTE: In transforming the zm coordinates we use zmid rather than
c           zmagx because the eqdsk "reflection" routine, convert_eqdsk,
c           invoked for the 'dntop' mesh uses z=zmid as the reflection plane.

c     NOTE: the iv indices map as follows:
c        1(lower)->2(upper)
c        3(lower)->4(upper)
c        2(lower)->1(upper)
c        4(lower)->3(upper)

c     The left half of dntop maps to the upper-left quadrant
c     of the dnull configuration by shifting & reflecting the ix index:
c              ix_upper = nxshift - ix_lower
      nxshift  = ixmdp(1) + ixtopu
      ixrb(1)  = nxshift  - 1
      ixpt2(1) = ixrb(1)  - ixpt1u
      do ix=ixmdp(1)+1,ixrb(1)+1
         do iy=1,nym
            do iv=0,4
               if (iv==0) ivt=0
               if (iv==1) ivt=2
               if (iv==2) ivt=1
               if (iv==3) ivt=4
               if (iv==4) ivt=3
               rm(ix,iy,iv)=rmu(nxshift-ix,iy,ivt)
               zm(ix,iy,iv)=2*zmid - zmu(nxshift-ix,iy,ivt)
            enddo
         enddo
      enddo

c     The right half of dntop maps to the upper-right quadrant
c     of the dnull configuration by shifting & reflecting the ix index:
c              ix_upper = nxshift - ix_lower
      nxshift  = ixmdp(2) + ixtopu + 2
      ixlb(2)  = nxshift  - (nxmu+1)
      ixpt1(2) = nxshift  - (ixpt2u+1)
      do ix=ixlb(2),ixmdp(2)
         do iy=1,nym
            do iv=0,4
               if (iv==0) ivt=0
               if (iv==1) ivt=2
               if (iv==2) ivt=1
               if (iv==3) ivt=4
               if (iv==4) ivt=3
               rm(ix,iy,iv)=rmu(nxshift-ix,iy,ivt)
               zm(ix,iy,iv)=2*zmid - zmu(nxshift-ix,iy,ivt)
            enddo
         enddo
      enddo

      return
      end

c----------------------------------------------------------------------c





