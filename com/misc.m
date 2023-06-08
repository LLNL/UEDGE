      function scal10 (y)
*     ------------------------------------------------------------------
*     This service routine can be used to compute a proper scale for
*     graphical output. scal10(y) returns the least in magnitude
*     real value r with abs(r) .ge. abs(y) and which has a decimal
*     representation of the form
*           r  =  0,  or
*           r  =  sgn(y) * i1 * 10.0 ** i2
*     where i1 is 2 or 5 or 10, and i2 is an integer.
*     ------------------------------------------------------------------
      real y
      integer i1, i2
      real r1
      real scal10
      real one
      data one /1./
*     ------------------------------------------------------------------
      if (y .eq. 0) then
         scal10 = 0
      else
         i2 = nint (log10(abs(y)) - 0.5)
         r1 = abs(y) * 10.0e0 ** (-i2)
*        r1 is in the range (1 .. 10), to machine accuracy.
         if (r1 .le. 2) then
            i1 = 2
         else if (r1 .le. 5) then
            i1 = 5
         else
            i1 = 10
         end if
         scal10 = sign(one, y) * i1 * 10.0 ** (i2)
* xlf95 wants the above.  Cannot replace one by 1.0
      end if
*     ------------------------------------------------------------------
      return
*     ------------------------------------------------------------------
      end

c *********************************************************************
      subroutine quadsvr (neq,xs,xs1,xs2,xs3,ys1,ys2,ys3,yq1,yq2)

cc    This routine calculates two guesses to all of the plasma
cc    variables base on a quadratic for the control parameter xs.
cc    The solution closest to the present one is in yq1, and the
cc    second one is in yq2.

      implicit none

c **--input variables--**
      integer neq           # number of equations
      real xs,xs1,xs2,xs3   # new xs & previous values xs1,xs2, & xs3
      real ys1(neq),ys2(neq),ys3(neq)  # previous solutions
                            	       # correponding to xs1, xs2, xs3

c **--output variables--**    # Two solutions to quadratic, yq1 and yq2
      real yq1(neq),yq2(neq)  # yq1 is closest to solution ys1

c **--local variables--**
      integer iq
      real dn1,dn2,dn3,aq,bq,cq,discr,srt,yqtemp
      real epsi
      real zero
      data zero /0/
      data epsi/1.e-30/

      do iq = 1, neq
         dn1 = xs1 / ( (ys2(iq)-ys1(iq))*(ys3(iq)-ys1(iq)) + epsi)
         dn2 = xs2 / ( (ys1(iq)-ys2(iq))*(ys3(iq)-ys2(iq)) + epsi)
         dn3 = xs3 / ( (ys1(iq)-ys3(iq))*(ys2(iq)-ys3(iq)) + epsi)
         aq = dn1 + dn2 + dn3
         bq = - (ys2(iq)+ys3(iq))*dn1 - (ys1(iq)+ys3(iq))*dn2 - 
     .                                  (ys1(iq)+ys2(iq))*dn3
         cq = ys2(iq)*ys3(iq)*dn1 + ys1(iq)*ys3(iq)*dn2 + 
     .                              ys1(iq)*ys2(iq)*dn3 - xs
         discr = max( bq*bq - 4*aq*cq, zero)
* xlf95 wants the above.  Cannot replace zero by 0.
         srt = sqrt(discr)
         yq1(iq) = 0.5 * ( -bq + srt ) / aq
         yq2(iq) = 0.5 * ( -bq - srt ) / aq
         if ( abs(ys1(iq)-yq1(iq)) .gt. abs(ys1(iq)-yq2(iq)) ) then
            yqtemp = yq1(iq)
            yq1(iq) = yq2(iq)
            yq2(iq) = yqtemp
         endif
      enddo 

      return
      end
c ****end of subroutine quadsvr**********

c ****Function gettime ******************

      real(Size4) function gettime(clk)
      real(Size4) clk, rcount, rrate
      integer count, rate
      gettime = 0
      call system_clock(count, rate)
      if (rate .ne. 0) then
         rcount = count
         rrate = rate
         gettime = rcount / rrate
      endif
      return
      end
c ****end of function gettime**********
