c     ------------------------------------------------------------------

      subroutine flxrun
      implicit none
Use(Dimflxgrd)
Use(Comflxgrd)

#     main driver routine for flx package
      call inflx
      call flxgen
### No need to write flx-grd file (MER 08 Feb 2010)
### These variables are now in com package and accessible to both flx and grd
      if (isfw==1) call flxfin

      return
      end

