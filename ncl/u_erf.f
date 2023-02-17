      FUNCTION U_ERF(x)
!***********************************************************************
!U_ERF evaluates the error function
!  W.A.Houlberg 1/99
!Input:
!  x-argument of error function
!***********************************************************************
      IMPLICIT NONE
      REAL           U_ERF
!Declaration of input variables
      REAL           x
!Declaration of local variables
      INTEGER        i
      REAL           a,                       b,
     #               c,                       d,
     #               del,                     gln,
     #               h,                       x2
      gln=5.723649e-1
      x2=x**2
      U_ERF=0.0
      IF(x2.lt.1.5) THEN
        a=0.5
        b=2.0
        del=b
        DO i=1,100
          a=a+1.0
          del=del*x2/a
          b=b+del
          IF(ABS(del).lt.ABS(b)*3.0e-7) THEN
            U_ERF=b*EXP(-x2+0.5*LOG(x2)-gln)
            IF(x.lt.0.0) U_ERF=-U_ERF
            GOTO 1000
          ENDIF
        ENDDO
      ELSE
        b=x2+0.5
        c=1.0/1.0e-30
        d=1.0/b
        h=d
        DO i=1,100
          a=-i*(i-0.5)
          b=b+2.0
          d=a*d+b
          IF(ABS(d).lt.1.0e-30) d=1.0e-30
          c=b+a/c
          IF(ABS(c).lt.1.0e-30) c=1.0e-30
          d=1.0/d
          del=d*c
          h=h*del
          IF(ABS(del-1.0).lt.3.0e-7) THEN
            U_ERF=1.0-EXP(-x2+0.5*LOG(x2)-gln)*h
            IF(x.lt.0.0) U_ERF=-U_ERF
            GOTO 1000
          ENDIF
        ENDDO
      ENDIF
 1000 RETURN
      END
