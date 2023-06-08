c!include "api.h"
      subroutine imp_sorc_walls (nxd, nzsptd, xcwi, xcwo, syi, syo,
     .                                     ixp1i, ixp1o, fnzysi, fnzyso)

c ... Calculate impurity-source profiles along the inner and outer walls
c     for a number of impurity sources.

      implicit none

c ... Input arguments:
      integer nxd               # number of physical cells along wall
      integer nzsptd            # total number of impurity species
      integer ixp1i(0:nxd+1)    # poloidal index of cell physically adjacent
                                # to right face of cells along inner wall
      integer ixp1o(0:nxd+1)    # poloidal index of cell physically adjacent
                                # to right face of cells along outer wall
      real xcwi(0:nxd+1)        # x-coord. of cell center on inner wall
      real xcwo(0:nxd+1)        # x-coord. of cell center on outer wall
      real syi(0:nxd+1)         # surface areas of inner-wall cells
      real syo(0:nxd+1)         # surface areas of outer-wall cells

c ... Output arguments:
      real fnzysi(0:nxd+1,nzsptd)   # source profiles along inner wall
      real fnzyso(0:nxd+1,nzsptd)   # source profiles along outer wall

c ... Common block:
      Use(Dim)                # nxpt
      Use(Xpoint_indices)     # ixlb,ixrb
      Use(Sources_at_walls)   # nzsor,iszsorlb,ximpi,ximpo,wimpi,wimpo,
                              # impsori,impsoro
                              # jxzsori,jxzsoro
                              # ixzbegi,ixbego,ixzendi,ixzendo

c ... Local variables:
      integer iz, ix, isor, jxlbi, jxrbi, jxlbo, jxrbo

c ... Check that array limits are not exceeded.
      if (nzsor .gt. NZSORMX) call xerrab(
     .  '*** nzsor > NZSORMX; enlarge impurity-source arrays.')

c ... Zero output arrays.
      do iz = 1, nzsptd
         do ix = 0, nxd+1
            fnzysi(ix,iz) = 0.
            fnzyso(ix,iz) = 0.
         enddo
      enddo

c ... Find the beginning and ending ix-indices for each wall
      do isor = 1, nzsor
         do iz = 1, nzsptd
            jxlbo = jxzsoro(iz,isor)
c     On the outboard wall, sources begin and end in the same mesh region:
            jxrbo = jxlbo
            ixzbego(iz,isor) = ixlb(jxlbo)
            ixzendo(iz,isor) = ixrb(jxrbo)+1

            jxlbi = jxzsori(iz,isor)
c     On private flux walls, sources may begin and end in different mesh regions
c     if nxpt > 1 as in a full double-null configuration.
            if (jxlbi==1) then
               jxrbi = nxpt
            else
               jxrbi = jxlbi - 1
            endif
            ixzbegi(iz,isor) = ixlb(jxlbi)
            ixzendi(iz,isor) = ixrb(jxrbi)+1
         enddo
      enddo

c ... Loop over impurity sources, doing the outer and inner walls.
      do isor = 1, nzsor
         call imp_sorc (nxd, nzsptd, iszsorlb(1:,isor),
     .                  ixzbego(1:,isor), ixzendo(1:,isor), ixp1o,
     .                  ximpo(1:,isor), wimpo(1:,isor), impsoro(1:,isor),
     .                                               xcwo, syo, fnzyso)
         call imp_sorc (nxd, nzsptd, iszsorlb(1:,isor),
     .                  ixzbegi(1:,isor), ixzendi(1:,isor), ixp1i,
     .                  ximpi(1:,isor), wimpi(1:,isor), impsori(1:,isor),
     .                                               xcwi, syi, fnzysi)
      enddo
   
      return
      end
c----------------------------------------------------------------------c
      subroutine imp_sorc (nx, nzspt, iszsorlb, ixbeg, ixend, ixplus1,
     .                           ximp, wimp, impsor, xgrd, sy, fnzys)

c ... Calculate impurity-source profiles along a single wall for a set
c     of impurities.

      implicit none

c ... Input arguments:
      integer nx               # number of physical cells along wall
      integer nzspt            # number of impurity species
      integer iszsorlb(nzspt)  # =1 if origin for ximp is at left bndry
      integer ixbeg(nzspt)     # beginning ix-index for wall
      integer ixend(nzspt)     # ending ix-index for wall
      integer ixplus1(0:nx+1)  # poloidal index of cell physically adjacent
                               # to right face
      real ximp(nzspt)         # centers of impurity-source profiles
      real wimp(nzspt)         # widths of impurity-source profiles
      real impsor(nzspt)       # impurity source strengths (Amps)
      real xgrd(0:nx+1)        # cell locations
      real sy(0:nx+1)          # surface areas of cells

c ... In-out argument:
      real fnzys(0:nx+1,nzspt)

c ... Common block:
      Use(Physical_constants2) # qe2

c ... Local variables:
      integer iz, ix
      real pi, sycos, xnot, arg

      pi = Pi

c ... Begin loop over impurity species.
      do iz = 1, nzspt

c ... Initialize normalization factor and center of profile
         sycos = 0.
         xnot = iszsorlb(iz)*ximp(iz)
     .            + (1-iszsorlb(iz))*(xgrd(ixend(iz))-ximp(iz))

c ... Calculate normalization factor
         ix = ixbeg(iz)
         do
            arg = (xgrd(ix) - xnot) * pi / (wimp(iz) + 1.e-20)
            if (abs(arg) .lt. pi/2.) sycos = sycos + cos(arg) * sy(ix)
            if (ix==ixend(iz)) break
            ix = ixplus1(ix)
         enddo

c ... Add impurity-source profile to any previous contributions.
         ix = ixbeg(iz)
         do
            arg = (xgrd(ix) - xnot) * pi / (wimp(iz) + 1.e-20)
            if (abs(arg) .lt. pi/2. .and. sycos .gt. 0.) then
               fnzys(ix,iz) = fnzys(ix,iz) +
     .                 impsor(iz) * cos(arg) * sy(ix) / (sycos * qe2)
            endif
            if (ix==ixend(iz)) break
            ix = ixplus1(ix)
         enddo

c ... End loop over impurity species.
      enddo
         
      return
      end
