c-----------------------------------------------------------------------

      subroutine grdrun
      implicit none
Use(Dimflxgrd)	#npts
Use(Comflxgrd)
Use(Dimensions)	#ndata,nbkpt,nwdim,niwdim

      external grdinit,readflx,ingrd,grdgen,writeue,remark

c     main driver routine for grd package

      call grdinit

### No need to read flx-grd file (MER 08 Feb 2010)
### These variables are now in com package and accessible to both flx and grd
      if (isfw==1) then	# option to read flx-grd file 
         call readflx
      else		# option to transfer data via com package variables
         call gallot("Curves",0)	#uses npts,jdim
         call gallot("Transfm",0)	#uses npts,jdim,nalpha
         call gallot("Spline",0)	#uses npts,jdim,mseg
         call gallot("Transit",0)	#uses npts
         call gallot("System",0)	#uses mseg,jdim,npts,noregs
c     ...Set angle-like indices and dimensioning parameter idim
         call setidim
         call gchange("Inmesh",0)	#uses idim,noregs,npts
					#save seedxp,seedxpxl for kxmesh=0
         call gallot("Linkco",0)	#uses idim,jdim
         call gchange("Mmod",0)	#uses idim,jdim,npts
				#save nplate,rplate,zplate data when iplate=1
c     ...Set storage for SLATEC spline routine FC --
         ndata=npts
         nbkpt=npts
         nwdim=8*npts+11	# >=5*(nbkpt-1)+2*max(ndata,nbkpt)+nbkpt+16
         niwdim=2*npts-6	# >=2*(nbkpt-3)
         call gallot("Argfc",0)	#uses npts,nconst,nwdim,niwdim
c     ...Initialize x,ycurveg arrays with data from flx pkg in x,ycurve
         call copyflx
      endif		# end of option for flx-grd file

      call ingrd
      call grdgen
      call writeue

      write(*,*) '***** Grid generation has been completed'

      return
      end

