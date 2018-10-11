      real function ssum(n, x, incx)
c ... Sums the n elements of a single precision vector x.  The increment
c     incx between elements of x must be non-zero but can be positive or
c     negative to cause elements to be added in either the forward or
c     reverse directions.

      implicit none

      integer n, incx
      real x(*)
      integer ib, ie, i

      ssum = 0.
      if (n .lt. 0) return
      ib = 1
      ie = 1 + (n-1) * abs(incx)
      if (incx .lt. 0) then
         ib = ie
         ie = 1
      endif
      do i = ib, ie, incx
         ssum = ssum + x(i)
      enddo

      return
      end
