      SUBROUTINE RARRAY_COPY(n,x,incx,y,incy)
!***********************************************************************
!RARRAY_COPY copies elements of the array x into y
!  W.A.Houlberg 12/98
!Input:  
!  n-number of elements to be copied
!  x-array to be copied
!  incx-increment in x index
!  incy-increment in y index
!Output:
!  y-new array
!***********************************************************************
      IMPLICIT NONE
!Declaration of input variables
      INTEGER        incx,                    incy,
     #               n
      REAL           x(*)
!Declaration of output variables
      REAL           y(*)
!Declaration of local variables
      INTEGER        i,                       ix,
     #               iy
      ix=1
      iy=1
      DO i=1,n
        y(iy)=x(ix)
        ix=ix+incx
        iy=iy+incy
      ENDDO   
      RETURN
      END
