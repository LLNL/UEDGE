      SUBROUTINE MNBRAK(iprint,maxmag,AX,BX,CX,FA,FB,FC,FUNC)
c ... Given a function FUNC, and given distinct initial points AX and BX,
c     this routine searches in the downhill direction (defined by the
c     function as evaluated at the initial points) and returns new points
c     AX, BX, CX that bracket a minimum of the function.  Also returned
c     are the function values at the three points, FA, FB, and FC.

      implicit none

c ... Input variables:
      integer iprint
      real maxmag

c ... In-out variables:
      real AX, BX

c ... Output variables:
      real CX, FA, FB, FC

c ... External function:
      external FUNC
      real FUNC

c ... Local variables:
      real ULIM, U, R, Q, FU, DUM

c ... GOLD is the default ratio by which successive intervals are
c     magnified; GLIMIT is the maximum magnification allowed for a
c     parabolic-fit step.
      real GOLD, GLIMIT, TINY
      PARAMETER (GOLD=1.618034, TINY=1.E-20)
      GLIMIT = maxmag

c ... Description of the algorithm is in "Numerical Recipes in C" by
c     W.H. Press, et al., Sec. 10.1.
      FA=FUNC(AX)
      FB=FUNC(BX)
      IF(FB.GT.FA)THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BX+GOLD*(BX-AX)
      FC=FUNC(CX)
1     if (iprint .gt. 1) write(*,90) 'mnbrak:  2 older points have f(',
     .                      AX, ') = ', FA,
     .                               '                         and f(',
     .                      BX, ') = ', FB
      if (iprint .gt. 1) write(*,90) 'mnbrak:  new point has f(', CX,
     .                      ') = ', FC
 90   format(2(a,f11.7))
      IF(FB.GE.FC)THEN
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            GO TO 1
          ELSE IF(FU.GT.FB)THEN
            CX=U
            FC=FU
            GO TO 1
          ENDIF
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
            FU=FUNC(U)
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
          U=ULIM
          FU=FUNC(U)
        ELSE
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
        GO TO 1
      ENDIF
      if (iprint .ge. 1) write(*,92)
     .   'mnbrak:  bracketing points have f(', AX, ') = ', FA, ',',
     .   '                                f(', BX, ') = ', FB, ',',
     .   '                            and f(', CX, ') = ', FC
      RETURN
 92   format(a,f11.7,a,f11.7,a)
      END
