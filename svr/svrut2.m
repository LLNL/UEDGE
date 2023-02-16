      real FUNCTION R1MACH9 (IDUM)
      INTEGER IDUM
C-----------------------------------------------------------------------
C This routine computes the unit roundoff of the machine.
C This is defined as the smallest positive machine number
C u such that  1.0 + u .ne. 1.0
C
C Subroutines/functions called by R1MACH9.. None
C-----------------------------------------------------------------------
      real U, COMP
C
      U = 1.0E0
 10   U = U*0.5E0
      COMP = 1.0E0 + U
      IF (COMP .NE. 1.0E0) GO TO 10
      R1MACH9 = U*2.0E0
      RETURN
C----------------------- End of Function R1MACH9 ------------------------
      END
