      SUBROUTINE U_LU_DECOMP(a,n,ndim,indx,d,iflag)
!***********************************************************************
!U_LU_DECOMP performs an LU decomposition of the matrix a ans is called
!   prior to U_LU_BACKSUB to solve linear equations or to invert a
!   matrix
!References:
!  Flannery, Teukolsky, Vetterling, Numerical Recipes
!  W.A.Houlberg 1/99
!Input:
!  a(ndim,ndim)-square matrix, overwritten on return
!  n-number of equations to be solved
!  ndim-first dimension of a array
!Output:
!  indx-vector showing row permutations due to partial pivoting
!  d-flag for number of row exchanges
!   =1.0 even number
!   =-1.0 odd number
!  iflag-error flag.
!       =0 no errors
!       =1 singular matrix
!***********************************************************************
      IMPLICIT NONE
!Declaration of parameters
      INTEGER         nmax
      PARAMETER      (nmax=100)
!Declaration of input
      INTEGER         iflag,                indx(*),
     #                n,                    ndim
      REAL            a(ndim,*),            d
!Declaration of local variables
      INTEGER         i,                    imax,
     #                j,                    k
      REAL            aamax,                dum,
     #                sum,                  vv(nmax)
!Initialization
      iflag=0
      d=1.0
!Loop over rows to get the implicit scaling information
      DO i=1,n
        aamax=0.0
        DO j=1,n
          IF(ABS(a(i,j)).gt.aamax) aamax=ABS(a(i,j))
        ENDDO    
        IF(aamax.eq.0.0) THEN
          iflag=1
          GOTO 1000
        ENDIF
        vv(i)=1.0/aamax
      ENDDO   
!Loop over columns using Crout's method
      DO j=1,n
        DO i=1,j-1
          sum=a(i,j)
          DO k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          ENDDO   
          a(i,j)=sum
        ENDDO    
!  Search for largest pivot element using dum as a figure of merit
        aamax=0.0
        DO i=j,n
          sum=a(i,j)
          DO k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          ENDDO    
          a(i,j)=sum
          dum=vv(i)*ABS(sum)
          IF(dum.ge.aamax) THEN
            imax=i
            aamax=dum
          ENDIF
        ENDDO   
        IF(j.ne.imax) THEN
!         Interchange rows
          DO k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          ENDDO    
          d=-d
          vv(imax)=vv(j)
        ENDIF
        indx(j)=imax
        IF(a(j,j).eq.0.0) THEN
          iflag=1
          RETURN
        ENDIF
        IF(j.ne.n) THEN
!         Divide by pivot element
          dum=1.0/a(j,j)
          DO i=j+1,n
            a(i,j)=a(i,j)*dum
          ENDDO   
        ENDIF
      ENDDO   
 1000 RETURN
      END
