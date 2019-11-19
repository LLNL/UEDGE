c!include "../mppl.h"
c!include "../sptodp.h"
      function ssmin (n, sx, incx)
      integer           n, incx
      real              sx(*)
      real ssmin
*     ------------------------------------------------------------------
*     Compute minimum value in the s.p. real vector described by
*     (n, sx, incx). BLAS storage conventions.
*     n must be positive.
*     ------------------------------------------------------------------
      integer           i, ix, m
      real              t
ccc   intrinsic         min, max, mod
      integer           ibase
      ibase(n, incx)  = max(1, 1+(1-n)*incx)
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         stop 'ssmin'
*        - - - - - - -
      else if (incx .eq. 1) then
*        use unrolled loops for unit increments.
         m = mod(n, 6)
         t = sx(1)
         do 1, i = 2, m
            t = min(t, sx(i))
    1    continue
         do 2, i = m+1, n, 6
            t = min(t, sx(i), sx(i+1), sx(i+2),
     .            sx(i+3), sx(i+4), sx(i+5))
    2    continue
         ssmin = t
         return
*        - - - - - - -
      else
         ix = ibase(n, incx)
         t = sx(ix)
         ix = ix + incx
         do 3, i = 2, n
            t = min(t, sx(ix))
            ix = ix + incx
    3    continue
         ssmin = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/06/14---
      end
      function smax (n, sx, incx)
      integer           n, incx
      real              sx(*)
      real smax
*     ------------------------------------------------------------------
*     Compute maximum value in the s.p. real vector described by
*     (n, sx, incx). BLAS storage conventions.
*     n must be positive.
*     ------------------------------------------------------------------
      integer           i, ix, m
      real              t
ccc   intrinsic         max, mod
      integer           ibase
      ibase(n, incx)  = max(1, 1+(1-n)*incx)
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         stop 'smax'
*        - - - - - - -
      else if (incx .eq. 1) then
*        use unrolled loops for unit increments.
         m = mod(n, 6)
         t = sx(1)
         do 1, i = 2, m
            t = max(t, sx(i))
    1    continue
         do 2, i = m+1, n, 6
            t = max(t, sx(i), sx(i+1), sx(i+2),
     .            sx(i+3), sx(i+4), sx(i+5))
    2    continue
         smax = t
         return
*        - - - - - - -
      else
         ix = ibase(n, incx)
         t = sx(ix)
         ix = ix + incx
         do 3, i = 2, n
            t = max(t, sx(ix))
            ix = ix + incx
    3    continue
         smax = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/06/14---
      end
      function samax (n, sx, incx)
      integer           n, incx
      real              sx(*)
      real samax
*     ------------------------------------------------------------------
*     Compute maximum absolute value in the s.p. real vector described
*     by (n, sx, incx). BLAS storage conventions.
*     n must be positive.
*     ------------------------------------------------------------------
      integer           i, ix, m
      real              t
ccc   intrinsic         max, mod, abs
      integer           ibase
      ibase(n, incx)  = max(1, 1+(1-n)*incx)
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         stop 'samax'
*        - - - - - - -
      else if (incx .eq. 1) then
*        use unrolled loops for unit increments.
         m = mod(n, 5)
         t = abs(sx(1))
         do 1, i = 2, m
            t = max(t, abs(sx(i)))
    1    continue
         do 2, i = m+1, n, 5
            t = max(t, abs(sx(i)), abs(sx(i+1)),
     .            abs(sx(i+2)), abs(sx(i+3)), abs(sx(i+4)))
    2    continue
         samax = t
         return
*        - - - - - - -
      else
         ix = ibase(n, incx)
         t = abs(sx(ix))
         ix = ix + incx
         do 3, i = 2, n
            t = max(t, abs(sx(ix)))
            ix = ix + incx
    3    continue
         samax = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/06/14---
      end
      subroutine sfill (n, sa, sx, incx)
      integer           n, incx
      real              sa, sx(*)
*     ------------------------------------------------------------------
*     Set all elements of the s.p. real vector x equal to sa. x is
*     described by (n, sx, incx), with BLAS storage conventions.
*     Return without effects if n .lt. 1.
*     Undefined if incx .eq. 0.
*     ------------------------------------------------------------------
      integer           i, ix
ccc   intrinsic         max
      integer           ibase
      ibase(n, incx)  = max(1, 1+(1-n)*incx)
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         return
*        - - - - - - -
      else if (incx .eq. 1) then
         do 1, i = 1, n
            sx(i) = sa
    1    continue
         return
*        - - - - - - -
      else
         ix = ibase(n, incx)
         do 3, i = 1, n
            sx(ix) = sa
            ix = ix + incx
    3    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/09/30---
      end
      subroutine sbias (n, sa, sx, incx)
      integer           n, incx
      real              sa, sx(*)
*     ------------------------------------------------------------------
*     Replace the s.p. real vector x by sa + x. x is described by
*     (n, sx, incx), with BLAS storage conventions.
*     Return without effects if n .lt. 1 or if sa .eq. 0.
*     Undefined if incx .eq. 0.
*     ------------------------------------------------------------------
      integer           i, ix
ccc   intrinsic         max
      integer           ibase
      ibase(n, incx)  = max(1, 1+(1-n)*incx)
*     ------------------------------------------------------------------
      if (n .lt. 1 .or. sa .eq. 0) then
         return
*        - - - - - - -
      else if (incx .eq. 1) then
         do 1, i = 1, n
            sx(i) = sa + sx(i)
    1    continue
         return
*        - - - - - - -
      else
         ix = ibase(n, incx)
         do 3, i = 1, n
            sx(ix) = sa + sx(ix)
            ix = ix + incx
    3    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/09/30---
      end
      function sadot (n, sx, incx, sy, incy)
      integer           n, incx, incy
      real              sx(*), sy(*)
      real sadot
*     ------------------------------------------------------------------
*     Compute the inner product of absolute values in two s.p. real
*     vectors described by (n, sx, incx) and (n, sy, incy).
*     BLAS storage conventions.
*     Return with sadot = 0 when n .lt. 1.
*     ------------------------------------------------------------------
      integer           i, ix, iy
      real              t
ccc   intrinsic         max, abs
      integer           ibase
      ibase(n, incx)  = max(1, 1+(1-n)*incx)
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         sadot = 0
         return
*        - - - - - - -
      else if (incx .eq. 1 .and. incy .eq. 1) then
         t = 0
         do 1, i = 1, n
            t = t + abs(sx(i)*sy(i))
    1    continue
         sadot = t
         return
*        - - - - - - -
      else if (incx .eq. incy .and. incx .gt. 0) then
         t = 0
         do 3, i = 1, n*incx, incx
            t = t + abs(sx(i)*sy(i))
    3    continue
         sadot = t
         return
*        - - - - - - -
      else
         t = 0
         ix = ibase(n, incx)
         iy = ibase(n, incy)
         do 4, i = 1, n
            t = t + abs(sx(ix)*sy(iy))
            ix = ix + incx
            iy = iy + incy
    4    continue
         sadot = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-87/07/26---
      end
      subroutine sadd (n, sx, incx, sy, incy)
      integer           n, incx, incy
      real              sx(*), sy(*)
*     ------------------------------------------------------------------
*     Add two s.p. real vectors: y := y + x. x and y are described
*     by (n, sx, incx) and (n, sy, incy). BLAS storage conventions.
*     Return without effects if n .lt. 1.
*     Undefined if incy .eq. 0.
*     ------------------------------------------------------------------
      integer           i, ix, iy
ccc   intrinsic         max
      integer           ibase
      ibase(n, incx)  = max(1, 1+(1-n)*incx)
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         return
*        - - - - - - -
      else if (incx .eq. 1 .and. incy .eq. 1) then
         do 1, i = 1, n
            sy(i) = sy(i) + sx(i)
    1    continue
         return
*        - - - - - - -
      else if (incx .eq. incy .and. incx .gt. 0) then
         do 3, i = 1, n*incx, incx
            sy(i) = sy(i) + sx(i)
    3    continue
         return
*        - - - - - - -
      else
         ix = ibase(n, incx)
         iy = ibase(n, incy)
         do 4, i = 1, n
            sy(iy) = sy(iy) + sx(ix)
            ix = ix + incx
            iy = iy + incy
    4    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/09/30---
      end
      subroutine ssub (n, sx, incx, sy, incy)
      integer           n, incx, incy
      real              sx(*), sy(*)
*     ------------------------------------------------------------------
*     Subtract two s.p. real vectors: y := y - x. x and y are described
*     by (n, sx, incx) and (n, sy, incy). BLAS storage conventions.
*     Return without effects if n .lt. 1.
*     Undefined if incy .eq. 0.
*     ------------------------------------------------------------------
      integer           i, ix, iy
ccc   intrinsic         max
      integer           ibase
      ibase(n, incx)  = max(1, 1+(1-n)*incx)
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         return
*        - - - - - - -
      else if (incx .eq. 1 .and. incy .eq. 1) then
         do 1, i = 1, n
            sy(i) = sy(i) - sx(i)
    1    continue
         return
*        - - - - - - -
      else if (incx .eq. incy .and. incx .gt. 0) then
         do 3, i = 1, n*incx, incx
            sy(i) = sy(i) - sx(i)
    3    continue
         return
*        - - - - - - -
      else
         ix = ibase(n, incx)
         iy = ibase(n, incy)
         do 4, i = 1, n
            sy(iy) = sy(iy) - sx(ix)
            ix = ix + incx
            iy = iy + incy
    4    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/09/30---
      end
      subroutine scopym (n, sx, incx, sy, incy)
      integer           n, incx, incy
      real              sx(*), sy(*)
*     ------------------------------------------------------------------
*     Copy the negative of the s.p. real vector described by (n, sx,
*     incx) to the vector (n, sy, incy). BLAS storage conventions.
*     Return without effects if n .lt. 1.
*     Undefined if incy .eq. 0.
*     ------------------------------------------------------------------
      integer           i, ix, iy
ccc   intrinsic         max
      integer           ibase
      ibase(n, incx)  = max(1, 1+(1-n)*incx)
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         return
*        - - - - - - -
      else if (incx .eq. 1 .and. incy .eq. 1) then
         do 1, i = 1, n
            sy(i) = - sx(i)
    1    continue
         return
*        - - - - - - -
      else if (incx .eq. incy .and. incx .gt. 0) then
         do 3, i = 1, n*incx, incx
            sy(i) = - sx(i)
    3    continue
         return
*        - - - - - - -
      else
         ix = ibase(n, incx)
         iy = ibase(n, incy)
         do 4, i = 1, n
            sy(iy) = - sx(ix)
            ix = ix + incx
            iy = iy + incy
    4    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/09/30---
      end
      subroutine szpxy (n, sx, incx, sy, incy, sz, incz)
      integer           n, incx, incy, incz
      real              sx(*), sy(*), sz(*)
*     ------------------------------------------------------------------
*     Replace z by z + x * y. x, y, and z are s.p. real vectors,
*     described by (n, sx, incx), (n, sy, incy), and (n, sz, incz).
*     BLAS storage conventions.
*     The operations + and * are applied pointwise.
*     Return without effects if n .lt. 1.
*     Undefined if incz .eq. 0.
*     ------------------------------------------------------------------
      integer           i, ix, iy, iz
ccc   intrinsic         max
      integer           ibase
      ibase(n, incx)  = max(1, 1+(1-n)*incx)
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         return
*        - - - - - - -
      else if (incx .eq. 1 .and. incy .eq. 1 .and. incz .eq. 1) then
         do 1, i = 1, n
            sz(i) = sz(i) + sx(i) * sy(i)
    1    continue
         return
*        - - - - - - -
      else
         ix = ibase(n, incx)
         iy = ibase(n, incy)
         iz = ibase(n, incz)
         do 3, i = 1, n
            sz(iz) = sz(iz) + sx(ix) * sy(iy)
            ix = ix + incx
            iy = iy + incy
            iz = iz + incz
    3    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/09/30---
      end
      subroutine szmxy (n, sx, incx, sy, incy, sz, incz)
      integer           n, incx, incy, incz
      real              sx(*), sy(*), sz(*)
*     ------------------------------------------------------------------
*     Replace z by z - x * y. x, y, and z are s.p. real vectors,
*     described by (n, sx, incx), (n, sy, incy), and (n, sz, incz).
*     BLAS storage conventions.
*     The operations - and * are applied pointwise.
*     Return without effects if n .lt. 1.
*     Undefined if incz .eq. 0.
*     ------------------------------------------------------------------
      integer           i, ix, iy, iz
ccc   intrinsic         max
      integer           ibase
      ibase(n, incx)  = max(1, 1+(1-n)*incx)
*     ------------------------------------------------------------------
      if (n .lt. 1) then
         return
*        - - - - - - -
      else if (incx .eq. 1 .and. incy .eq. 1 .and. incz .eq. 1) then
         do 1, i = 1, n
            sz(i) = sz(i) - sx(i) * sy(i)
    1    continue
         return
*        - - - - - - -
      else
         ix = ibase(n, incx)
         iy = ibase(n, incy)
         iz = ibase(n, incz)
         do 3, i = 1, n
            sz(iz) = sz(iz) - sx(ix) * sy(iy)
            ix = ix + incx
            iy = iy + incy
            iz = iz + incz
    3    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/09/30---
      end
      function s2min (n1, n2, sx, inc1x, inc2x)
      integer           n1, n2, inc1x, inc2x
      real              sx(*)
      real s2min
*     ------------------------------------------------------------------
*     Compute minimum value in the s.p. real matrix described by
*     (n1, n2, sx, inc1x, inc2x). BLAS**2 storage conventions.
*     n1 and n2 must be positive.
*     ------------------------------------------------------------------
      integer           i, i2, ibase
      real              t, ssmin
ccc   intrinsic         max, min
      external          ssmin
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         stop 's2min'
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x) then
         s2min = ssmin (n1*n2, sx, inc1x)
         return
*        - - - - - - -
      else
         i2 = ibase(n2, inc2x)
         t = ssmin (n1, sx(i2), inc1x)
         i2 = i2 + inc2x
         do 1, i = 2, n2
            t = min (t, ssmin (n1, sx(i2), inc1x))
            i2 = i2 + inc2x
    1    continue
         s2min = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      function s2max (n1, n2, sx, inc1x, inc2x)
      integer           n1, n2, inc1x, inc2x
      real              sx(*)
      real s2max
*     ------------------------------------------------------------------
*     Compute maximum value in the s.p. real matrix described by
*     (n1, n2, sx, inc1x, inc2x). BLAS**2 storage conventions.
*     n1 and n2 must be positive.
*     ------------------------------------------------------------------
      integer           i, i2, ibase
      real              t, smax
ccc   intrinsic         max
      external          smax
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         stop 's2max'
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x) then
         s2max = smax (n1*n2, sx, inc1x)
         return
*        - - - - - - -
      else
         i2 = ibase(n2, inc2x)
         t = smax (n1, sx(i2), inc1x)
         i2 = i2 + inc2x
         do 1, i = 2, n2
            t = max (t, smax (n1, sx(i2), inc1x))
            i2 = i2 + inc2x
    1    continue
         s2max = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      function s2amax (n1, n2, sx, inc1x, inc2x)
      integer           n1, n2, inc1x, inc2x
      real              sx(*)
      real s2amax
*     ------------------------------------------------------------------
*     Compute maximum absolute value in the s.p. real matrix described
*     by (n1, n2, sx, inc1x, inc2x). BLAS**2 storage conventions.
*     n1 and n2 must be positive.
*     ------------------------------------------------------------------
      integer           i, i2, ibase
      real              t, samax
ccc   intrinsic         max
      external          samax
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         stop 's2amax'
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x) then
         s2amax = samax (n1*n2, sx, inc1x)
         return
*        - - - - - - -
      else
         i2 = ibase(n2, inc2x)
         t = samax (n1, sx(i2), inc1x)
         i2 = i2 + inc2x
         do 1, i = 2, n2
            t = max (t, samax (n1, sx(i2), inc1x))
            i2 = i2 + inc2x
    1    continue
         s2amax = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      function s2sum (n1, n2, sx, inc1x, inc2x)
      integer           n1, n2, inc1x, inc2x
      real              sx(*)
      real s2sum
*     ------------------------------------------------------------------
*     Compute sum of values in the s.p. real matrix described by
*     (n1, n2, sx, inc1x, inc2x). BLAS**2 storage conventions.
*     Return with s2sum = 0 when n1 .lt. 1 or n2 .lt. 1.
*     ------------------------------------------------------------------
      integer           i, i2, ibase
      real              t, ssum
ccc   intrinsic         max
      external          ssum
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         s2sum = 0
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x) then
         s2sum = ssum (n1*n2, sx, inc1x)
         return
*        - - - - - - -
      else
         t = 0
         i2 = ibase(n2, inc2x)
         do 1, i = 1, n2
            t = t + ssum (n1, sx(i2), inc1x)
            i2 = i2 + inc2x
    1    continue
         s2sum = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      function s2asum (n1, n2, sx, inc1x, inc2x)
      integer           n1, n2, inc1x, inc2x
      real              sx(*)
      real s2asum
*     ------------------------------------------------------------------
*     Compute sum of absolute values in the s.p. real matrix described
*     by (n1, n2, sx, inc1x, inc2x). BLAS**2 storage conventions.
*     Return with s2asum = 0 when n1 .lt. 1 or n2 .lt. 1.
*     ------------------------------------------------------------------
      integer           i, i2, ibase
      real              t, sasum
ccc   intrinsic         max
      external          sasum
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         s2asum = 0
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x) then
         s2asum = sasum (n1*n2, sx, inc1x)
         return
*        - - - - - - -
      else
         t = 0
         i2 = ibase(n2, inc2x)
         do 1, i = 1, n2
            t = t + sasum (n1, sx(i2), inc1x)
            i2 = i2 + inc2x
    1    continue
         s2asum = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      subroutine s2fill (n1, n2, sa, sx, inc1x, inc2x)
      integer           n1, n2, inc1x, inc2x
      real              sa, sx(*)
*     ------------------------------------------------------------------
*     Set all elements of the s.p. real matrix x equal to sa. x is
*     described by (n1, n2, sx, inc1x, inc2x), with BLAS**2 storage
*     conventions.
*     Return without effects if n1 .lt. 1 or n2 .lt. 1.
*     Undefined if inc1x .eq. 0 or inc2x .eq. 0.
*     ------------------------------------------------------------------
      integer           i, i2, ibase
ccc   intrinsic         max
      external          sfill
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x) then
         call sfill (n1*n2, sa, sx, inc1x)
         return
*        - - - - - - -
      else
         i2 = ibase(n2, inc2x)
         do 1, i = 1, n2
            call sfill (n1, sa, sx(i2), inc1x)
            i2 = i2 + inc2x
    1    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      subroutine s2scal (n1, n2, sa, sx, inc1x, inc2x)
      integer           n1, n2, inc1x, inc2x
      real              sa, sx(*)
*     ------------------------------------------------------------------
*     Replace the s.p. real matrix x by sa * x. x is described by
*     (n1, n2, sx, inc1x, inc2x), with BLAS**2 storage conventions.
*     Return without effects if min(n1, n2) .lt. 1 .or. sa .eq. 1.
*     Undefined if inc1x .eq. 0 or inc2x .eq. 0.
*     ------------------------------------------------------------------
      integer           i, i2, ibase
ccc   intrinsic         max
      external          sscal
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1 .or. sa .eq. 1) then
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x) then
         call sscal (n1*n2, sa, sx, inc1x)
         return
*        - - - - - - -
      else
         i2 = ibase(n2, inc2x)
         do 1, i = 1, n2
            call sscal (n1, sa, sx(i2), inc1x)
            i2 = i2 + inc2x
    1    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      function s2dot (n1, n2, sx, inc1x, inc2x, sy, inc1y, inc2y)
      integer           n1, n2, inc1x, inc2x, inc1y, inc2y
      real              sx(*), sy(*)
      real s2dot
*     ------------------------------------------------------------------
*     Compute the inner product of two s.p. real matrices described
*     by (n1, n2, sx, inc1x, inc2x) and (n1, n2, sy, inc1y, inc2y).
*     BLAS**2 storage conventions.
*     Return with s2dot = 0 when n1 .lt. 1 .or. n2 .lt. 1.
*     ------------------------------------------------------------------
      integer           i, i2x, i2y, ibase
      real              t
ccc   intrinsic         max
      real              sdot
      external          sdot
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         s2dot = 0
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x .and. n1*inc1y .eq. inc2y) then
         s2dot = sdot (n1*n2, sx, inc1x, sy, inc1y)
         return
*        - - - - - - -
      else
         t = 0
         i2x = ibase(n2, inc2x)
         i2y = ibase(n2, inc2y)
         do 1, i = 1, n2
            t = t + sdot (n1, sx(i2x), inc1x, sy(i2y), inc1y)
            i2x = i2x + inc2x
            i2y = i2y + inc2y
    1    continue
         s2dot = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-86/02/21---
      end
      function s2adot (n1, n2, sx, inc1x, inc2x, sy, inc1y, inc2y)
      integer           n1, n2, inc1x, inc2x, inc1y, inc2y
      real              sx(*), sy(*)
      real s2adot
*     ------------------------------------------------------------------
*     Compute the inner product of absolute values in two s.p. real
*     matrices described by (n1, n2, sx, inc1x, inc2x) and
*     (n1, n2, sy, inc1y, inc2y). BLAS**2 storage conventions.
*     Return with s2adot = 0 when n1 .lt. 1 .or. n2 .lt. 1.
*     ------------------------------------------------------------------
      integer           i, i2x, i2y, ibase
      real              t
ccc   intrinsic         max
      real              sadot
      external          sadot
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         s2adot = 0
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x .and. n1*inc1y .eq. inc2y) then
         s2adot = sadot (n1*n2, sx, inc1x, sy, inc1y)
         return
*        - - - - - - -
      else
         t = 0
         i2x = ibase(n2, inc2x)
         i2y = ibase(n2, inc2y)
         do 1, i = 1, n2
            t = t + sadot (n1, sx(i2x), inc1x, sy(i2y), inc1y)
            i2x = i2x + inc2x
            i2y = i2y + inc2y
    1    continue
         s2adot = t
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-87/07/26---
      end
      subroutine s2add (n1, n2, sx, inc1x, inc2x, sy, inc1y, inc2y)
      integer           n1, n2, inc1x, inc2x, inc1y, inc2y
      real              sx(*), sy(*)
*     ------------------------------------------------------------------
*     Add two s.p. real matrices: y := y + x. x and y are described
*     by (n1, n2, sx, inc1x, inc2x) and (n1, n2, sy, inc1y, inc2y),
*     with BLAS**2 storage conventions.
*     Return without effects if n1 .lt. 1 .or. n2 .lt. 1.
*     Undefined if inc1y .eq. 0 or inc2y .eq. 0.
*     ------------------------------------------------------------------
      integer           i, i2x, i2y, ibase
ccc   intrinsic         max
      external          sadd
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x .and. n1*inc1y .eq. inc2y) then
         call sadd (n1*n2, sx, inc1x, sy, inc1y)
         return
*        - - - - - - -
      else
         i2x = ibase(n2, inc2x)
         i2y = ibase(n2, inc2y)
         do 1, i = 1, n2
            call sadd (n1, sx(i2x), inc1x, sy(i2y), inc1y)
            i2x = i2x + inc2x
            i2y = i2y + inc2y
    1    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      subroutine s2sub (n1, n2, sx, inc1x, inc2x, sy, inc1y, inc2y)
      integer           n1, n2, inc1x, inc2x, inc1y, inc2y
      real              sx(*), sy(*)
*     ------------------------------------------------------------------
*     Subtract two s.p. real matrices: y := y - x. x, y are described
*     by (n1, n2, sx, inc1x, inc2x) and (n1, n2, sy, inc1y, inc2y),
*     with BLAS**2 storage conventions.
*     Return without effects if n1 .lt. 1 .or. n2 .lt. 1.
*     Undefined if inc1y .eq. 0 or inc2y .eq. 0.
*     ------------------------------------------------------------------
      integer           i, i2x, i2y, ibase
ccc   intrinsic         max
      external          ssub
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x .and. n1*inc1y .eq. inc2y) then
         call ssub (n1*n2, sx, inc1x, sy, inc1y)
         return
*        - - - - - - -
      else
         i2x = ibase(n2, inc2x)
         i2y = ibase(n2, inc2y)
         do 1, i = 1, n2
            call ssub (n1, sx(i2x), inc1x, sy(i2y), inc1y)
            i2x = i2x + inc2x
            i2y = i2y + inc2y
    1    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      subroutine s2copy (n1, n2, sx, inc1x, inc2x, sy, inc1y, inc2y)
      integer           n1, n2, inc1x, inc2x, inc1y, inc2y
      real              sx(*), sy(*)
*     ------------------------------------------------------------------
*     Copy the s.p. real matrix described by (n1, n2, sx, inc1x,
*     inc2x) to the matrix described by (n1, n2, sy, inc1y, inc2y),
*     with BLAS**2 storage conventions.
*     Return without effects if n1 .lt. 1 .or. n2 .lt. 1.
*     Undefined if inc1y .eq. 0 or inc2y .eq. 0.
*     ------------------------------------------------------------------
      integer           i, i2x, i2y, ibase
ccc   intrinsic         max
      external          scopy
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x .and. n1*inc1y .eq. inc2y) then
         call scopy (n1*n2, sx, inc1x, sy, inc1y)
         return
*        - - - - - - -
      else
         i2x = ibase(n2, inc2x)
         i2y = ibase(n2, inc2y)
         do 1, i = 1, n2
            call scopy (n1, sx(i2x), inc1x, sy(i2y), inc1y)
            i2x = i2x + inc2x
            i2y = i2y + inc2y
    1    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      subroutine s2copm (n1, n2, sx, inc1x, inc2x, sy, inc1y, inc2y)
      integer           n1, n2, inc1x, inc2x, inc1y, inc2y
      real              sx(*), sy(*)
*     ------------------------------------------------------------------
*     Copy the negative of the s.p. real matrix described by
*     (n1, n2, sx, inc1x, inc2x) to the matrix described by
*     (n1, n2, sy, inc1y, inc2y), with BLAS**2 storage conventions.
*     Return without effects if n1 .lt. 1 .or. n2 .lt. 1.
*     Undefined if inc1y .eq. 0 or inc2y .eq. 0.
*     ------------------------------------------------------------------
      integer           i, i2x, i2y, ibase
ccc   intrinsic         max
      external          scopym
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x .and. n1*inc1y .eq. inc2y) then
         call scopym (n1*n2, sx, inc1x, sy, inc1y)
         return
*        - - - - - - -
      else
         i2x = ibase(n2, inc2x)
         i2y = ibase(n2, inc2y)
         do 1, i = 1, n2
            call scopym (n1, sx(i2x), inc1x, sy(i2y), inc1y)
            i2x = i2x + inc2x
            i2y = i2y + inc2y
    1    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/09/30---
      end
      subroutine s2axpy (n1, n2, sa, sx, inc1x, inc2x,
     .                  sy, inc1y, inc2y)
      integer           n1, n2, inc1x, inc2x, inc1y, inc2y
      real              sa, sx(*), sy(*)
*     ------------------------------------------------------------------
*     Set y := y + sa * x. x and y are s.p. real matrices, described
*     by (n1, n2, sx, inc1x, inc2x) and (n1, n2, sy, inc1y, inc2y).
*     BLAS**2 storage conventions.
*     Return without effects if min(n1, n2) .lt. 1 or sa .eq. 0.
*     Undefined if inc1y .eq. 0 or inc2y .eq. 0.
*     ------------------------------------------------------------------
      integer           i, i2x, i2y, ibase
ccc   intrinsic         max
      external          saxpy
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x .and. n1*inc1y .eq. inc2y) then
         call saxpy (n1*n2, sa, sx, inc1x, sy, inc1y)
         return
*        - - - - - - -
      else
         i2x = ibase(n2, inc2x)
         i2y = ibase(n2, inc2y)
         do 1, i = 1, n2
            call saxpy (n1, sa, sx(i2x), inc1x, sy(i2y), inc1y)
            i2x = i2x + inc2x
            i2y = i2y + inc2y
    1    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      subroutine s2zpxy (n1, n2, sx, inc1x, inc2x, sy, inc1y, inc2y,
     .                  sz, inc1z, inc2z)
      integer           n1, n2, inc1x, inc2x, inc1y, inc2y,
     .                  inc1z, inc2z
      real              sx(*), sy(*), sz(*)
*     ------------------------------------------------------------------
*     Set z := z + x * y (pointwise). x, y, and z are s.p. real
*     matrices, described by (n1, n2, sx, inc1x, inc2x), (n1, n2,
*     sy, inc1y, in2y), and (n1, n2, sz, inc1z, inc2z). BLAS**2 storage
*     conventions.
*     Return without effects if n1 .lt. 1 or n2 .lt. 1.
*     Undefined if inc1z .eq. 0 or inc2z .eq. 0.
*     ------------------------------------------------------------------
      integer           i, i2x, i2y, i2z, ibase
ccc   intrinsic         max
      external          szpxy
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x .and. n1*inc1y .eq. inc2y
     .      .and. n1*inc1z .eq. inc2z) then
         call szpxy (n1*n2, sx, inc1x, sy, inc1y, sz, inc1z)
         return
*        - - - - - - -
      else
         i2x = ibase(n2, inc2x)
         i2y = ibase(n2, inc2y)
         i2z = ibase(n2, inc2z)
         do 1, i = 1, n2
            call szpxy (n1, sx(i2x), inc1x, sy(i2y), inc1y,
     .            sz(i2z), inc1z)
            i2x = i2x + inc2x
            i2y = i2y + inc2y
            i2z = i2z + inc2z
    1    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      subroutine s2zmxy (n1, n2, sx, inc1x, inc2x, sy, inc1y, inc2y,
     .                  sz, inc1z, inc2z)
      integer           n1, n2, inc1x, inc2x, inc1y, inc2y,
     .                  inc1z, inc2z
      real              sx(*), sy(*), sz(*)
*     ------------------------------------------------------------------
*     Set z := z - x * y (pointwise). x, y, and z are s.p. real
*     matrices, described by (n1, n2, sx, inc1x, inc2x), (n1, n2,
*     sy, inc1y, in2y), and (n1, n2, sz, inc1z, inc2z). BLAS**2 storage
*     conventions.
*     Return without effects if n1 .lt. 1 or n2 .lt. 1.
*     Undefined if inc1z .eq. 0 or inc2z .eq. 0.
*     ------------------------------------------------------------------
      integer           i, i2x, i2y, i2z, ibase
ccc   intrinsic         max
      external          szmxy
      ibase(n1, inc1x) = max(1, 1+(1-n1)*inc1x)
*     ------------------------------------------------------------------
      if (n1 .lt. 1 .or. n2 .lt. 1) then
         return
*        - - - - - - -
      else if (n1*inc1x .eq. inc2x .and. n1*inc1y .eq. inc2y
     .      .and. n1*inc1z .eq. inc2z) then
         call szmxy (n1*n2, sx, inc1x, sy, inc1y, sz, inc1z)
         return
*        - - - - - - -
      else
         i2x = ibase(n2, inc2x)
         i2y = ibase(n2, inc2y)
         i2z = ibase(n2, inc2z)
         do 1, i = 1, n2
            call szmxy (n1, sx(i2x), inc1x, sy(i2y), inc1y,
     .            sz(i2z), inc1z)
            i2x = i2x + inc2x
            i2y = i2y + inc2y
            i2z = i2z + inc2z
    1    continue
         return
*        - - - - - - -
      end if
*     ---------------------------------------------------BJB-83/04/06---
      end
      INTEGER FUNCTION IDAMAX_U(N,DX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
*
*     IDAMAX finds the index of element having max. absolute value.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS
*     ..
      IDAMAX_U = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX_U = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DMAX = DABS(DX(1))
         DO I = 2,N
            IF (DABS(DX(I)).GT.DMAX) THEN
               IDAMAX_U = I
               DMAX = DABS(DX(I))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         IX = 1
         DMAX = DABS(DX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DABS(DX(IX)).GT.DMAX) THEN
               IDAMAX_U = I
               DMAX = DABS(DX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
      DOUBLE PRECISION FUNCTION DASUM_U(N,DX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
*
*     DASUM takes the sum of the absolute values.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS,MOD
*     ..
      DASUM_U = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
*        code for increment equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,6)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DABS(DX(I))
            END DO
            IF (N.LT.6) THEN
               DASUM_U = DTEMP
               RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,6
            DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) +
     $              DABS(DX(I+2)) + DABS(DX(I+3)) +
     $              DABS(DX(I+4)) + DABS(DX(I+5))
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DTEMP = DTEMP + DABS(DX(I))
         END DO
      END IF
      DASUM_U = DTEMP
      RETURN
      END
      SUBROUTINE DAXPY_U(N,DA,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     DAXPY constant times a vector plus a vector.
*     uses unrolled loops for increments equal to one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DY(I) + DA*DX(I)
            END DO
         END IF
         IF (N.LT.4) RETURN
         MP1 = M + 1
         DO I = MP1,N,4
            DY(I) = DY(I) + DA*DX(I)
            DY(I+1) = DY(I+1) + DA*DX(I+1)
            DY(I+2) = DY(I+2) + DA*DX(I+2)
            DY(I+3) = DY(I+3) + DA*DX(I+3)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
         END DO
      END IF
      RETURN
      END
      SUBROUTINE DCOPY_U(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     DCOPY copies a vector, x, to a vector, y.
*     uses unrolled loops for increments equal to one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF   
         MP1 = M + 1
         DO I = MP1,N,7
            DY(I) = DX(I)
            DY(I+1) = DX(I+1)
            DY(I+2) = DX(I+2)
            DY(I+3) = DX(I+3)
            DY(I+4) = DX(I+4)
            DY(I+5) = DX(I+5)
            DY(I+6) = DX(I+6)
         END DO
      ELSE      
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DY(IY) = DX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT_U(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     DDOT forms the dot product of two vectors.
*     uses unrolled loops for increments equal to one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      DDOT_U = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DX(I)*DY(I)
            END DO
            IF (N.LT.5) THEN
               DDOT_U=DTEMP
            RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     $            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      DDOT_U = DTEMP
      RETURN
      END
      DOUBLE PRECISION FUNCTION DNRM2_U(N,X,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
*     ..
*
*  Purpose
*  =======
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*  Further Details
*  ===============
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
*
      DNRM2_U = NORM
      RETURN
*
*     End of DNRM2.
*
      END
      SUBROUTINE DSCAL_U(N,DA,DX,INCX)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  Purpose
*  =======
*
*     DSCAL scales a vector by a constant.
*     uses unrolled loops for increment equal to one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
      END IF
      RETURN
      END
      SUBROUTINE DSWAP_U(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     interchanges two vectors.
*     uses unrolled loops for increments equal one.
*
*  Further Details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*       code for both increments equal to 1
*
*
*       clean-up loop
*
         M = MOD(N,3)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DX(I)
               DX(I) = DY(I)
               DY(I) = DTEMP
            END DO
            IF (N.LT.3) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,3
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
            DTEMP = DX(I+1)
            DX(I+1) = DY(I+1)
            DY(I+1) = DTEMP
            DTEMP = DX(I+2)
            DX(I+2) = DY(I+2)
            DY(I+2) = DTEMP
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DX(IX)
            DX(IX) = DY(IY)
            DY(IY) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
c     subroutine checksizes
c 
c 
c     integer i
c     real x
c 
c     print *,"size of integer = ",sizeof(i)
c     print *,"size of real = ",sizeof(x)
c 
c     return
c     end


      
