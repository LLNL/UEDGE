      real FUNCTION BRENT(iprint,ITMAX,AX,BX,CX,F,TOL,XMIN)
c ... Given a function F, and given a bracketing triplet of abscissas
c     AX, BX, CX (such that BX is between AX and CX, and F(BX) is less
c     than both F(AX) and F(CX)), this routine isolates the minimum to
c     a fractional precision of about TOL using Brent's method.  The
c     abscissa of the minimum is returned as XMIN, and the minimum
c     function value is returned as BRENT, the returned function value.
c ... This version of BRENT uses an additional (first) argument ITMAX,
c     the maximum allowed number of iterations, and returns to the
c     Basis parser by calling kaboom if convergence is not achieved in
c     ITMAX iterations.

      implicit none

c ... Input variables:
      integer iprint
      integer ITMAX
      real AX, BX, CX, TOL

c ... Output variable:
      real XMIN

c ... External function:
      external F
      real F

c ... Local variables:
      integer ITER
      real A, B, D, ETEMP, FU, FV, FW, FX, P, Q, R, TOL1, TOL2
      real U, V, W, X, XM, E

c ... CGOLD is the golden ratio; ZEPS is a small number that protects
c     against trying to achieve fractional accuracy for a minimum that
c     happens to be exactly zero.
      real CGOLD, ZEPS
      PARAMETER (CGOLD=.3819660)
      ZEPS = tol

c ... Description of the algorithm is in "Numerical Recipes in C" by
c     W.H. Press, et al., Sec. 10.2.
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      FX=F(X)
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX
        if (iprint .gt. 1) write(*,90) 'brent:  iteration ', ITER,
     .                        ' -- test point has f(', X, ') = ', FX
 90     format(a,i2,2(a,f11.7))
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3
        IF(ABS(E).GT.TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.*(Q-R)
          IF(Q.GT.0.) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR. 
     *        P.GE.Q*(B-X)) GOTO 1
          D=P/Q
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
1       IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
2       IF(ABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=F(U)
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
11    CONTINUE
      call xerrab ('*** brent exceeded maximum iterations.')
3     XMIN=X
      BRENT=FX
      if (iprint .ge. 1) write(*,91) 'brent:  final point has f(',
     .                      X, ') = ', FX
 91   format(2(a,f11.7))
      RETURN
      END
