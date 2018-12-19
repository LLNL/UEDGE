      subroutine xerrab (msg)
      character*(*) msg
c ... Output the message in character string "msg" to standard out,
c     and then return control to the user.
c     Note that "remark" and "kaboom" are part of the Basis runtime
c     library.
      call remark(msg)
      call kaboom(0)
      return
      end
*DECK XERMShG
      SUBROUTINE XERMShG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C-----------------------------------------------------------------------
C Subroutines XERMShG, XSETFP, XSETUNP, and parmsetget, as
C given here, constitute a simplified version of the SLATEC error
C handling package.  Written by A. C. Hindmarsh, 18 November 1992.
c Modified by Gary R. Smith, 5 January 1994.
C
C All arguments are input arguments.
C LIBRAR = Library name (character array).  Prefixed to message.
C SUBROU = Routine name (character array).  Prefixed to message.
C MESSG  = The message (character array).
C NERR   = Integer error number.  Prefixed to message.
C LEVEL  = The error level..
C          0 or 1 means recoverable (control returns to caller).
C          2 means fatal (run is aborted--see note below).
C
C Note..  This routine has been simplified in the following ways..
C 1. A single prefix line is printed with NERR, SUBROU, and LIBRAR.
C 2. The message in MESSG is printed, unaltered, on lines of up to 72
C    characters each using a format of (A).
C 3. If LEVEL = 2, control passes to the statement   STOP
C    to abort the run.  This statement may be machine-dependent.
C
C For a different default logical unit number, change the data
C statement in subroutine parmsetget.
C For a different run-abort command, change the statement following
C statement 100 at the end.
C-----------------------------------------------------------------------
C Subroutine called by XERMShG.. parmsetget
C Function routines called by XERMShG.. None
C Intrinsic function used by XERMShG.. LEN
C-----------------------------------------------------------------------
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      INTEGER NERR, LEVEL
      INTEGER I1, I2, IL, LENMSG, LLEN, LUNIT, MESFLG, NLINES
      PARAMETER (LLEN = 72)
C
C Get message print flag and logical unit number. ----------------------
      call parmsetget (1, LUNIT, .false.)
      call parmsetget (2, MESFLG, .false.)
      IF (MESFLG .EQ. 0) GO TO 100
C Write NERR, SUBROU, and LIBRAR. --------------------------------------
      I1 = LEN(SUBROU)
      I2 = LEN(LIBRAR)
      WRITE (LUNIT, 10) NERR, SUBROU(1:I1), LIBRAR(1:I2)
  10  FORMAT(/,'***Error number ',I6,' from ',A,' in library ',A,'***')
C Write the message. ---------------------------------------------------
      LENMSG = LEN(MESSG)
      NLINES = ( (LENMSG - 1)/LLEN ) + 1
      DO 20 IL = 1,NLINES
        I1 = 1 + (IL - 1)*LLEN
        I2 = MIN(IL*LLEN,LENMSG)
        WRITE (LUNIT,'(A)') MESSG(I1:I2)
  20    CONTINUE
C Abort the run if LEVEL = 2. ------------------------------------------
 100  IF (LEVEL .NE. 2) RETURN
      call kaboom(0)
C----------------------- End of Subroutine XERMShG ----------------------
      END
*DECK XSETUNP
      SUBROUTINE XSETUNP (LUN)
C-----------------------------------------------------------------------
C This routine resets the logical unit number for messages.
C
C Subroutine called by XSETUNP.. parmsetget
C Function routines called by XSETUNP.. None
C-----------------------------------------------------------------------
      INTEGER LUN
C
      IF (LUN .GT. 0) call parmsetget (1,LUN,.true.)
      RETURN
C----------------------- End of Subroutine XSETUNP ---------------------
      END
*DECK XSETFP
      SUBROUTINE XSETFP (MFLAG)
C-----------------------------------------------------------------------
C This routine resets the print control flag MFLAG.
C
C Subroutine called by XSETFP.. parmsetget
C Function routines called by XSETFP.. None
C-----------------------------------------------------------------------
      INTEGER MFLAG
C
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1)
     .   call parmsetget (2,MFLAG,.true.)
      RETURN
C----------------------- End of Subroutine XSETFP ----------------------
      END
*DECK parmsetget
      SUBROUTINE parmsetget (IPAR, IVALUE, ISET)
      LOGICAL ISET
      INTEGER IPAR, IVALUE
C-----------------------------------------------------------------------
C parmsetget saves and recalls one of two error message parameters:
C   LUNIT, the logical unit number to which messages are printed, and
C   MESFLG, the message print flag.
C This is a modification of the SLATEC library routine J4SAVE.
C
C Saved local variables..
C  LUNIT  = Logical unit number for messages.
C           The default is 6 (machine-dependent).
C  MESFLG = Print control flag..
C           1 means print all messages (the default).
C           0 means no printing.
C
C On input..
C   IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
C   IVALUE = The value to be set for the parameter, if ISET = .true.
c            or the existing value of the parameter, if ISET = .false.
C   ISET   = Logical flag to indicate whether to read or write.
C            If ISET = .true., the parameter will be given
C            the value IVALUE.  If ISET = .false., the parameter
C            will be unchanged, and IVALUE contains the parameter
c            value on exit.
C
C Subroutines/functions called by parmsetget.. None
C-----------------------------------------------------------------------
      INTEGER LUNIT, MESFLG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this routine.
C-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/6/, MESFLG/1/
C
      IF (IPAR .EQ. 1) THEN
        IF (ISET) THEN
           LUNIT = IVALUE
        ELSE
           IVALUE = LUNIT
        ENDIF
      ENDIF
C
      IF (IPAR .EQ. 2) THEN
        IF (ISET) THEN
           MESFLG = IVALUE
        ELSE
           IVALUE = MESFLG
        ENDIF
      ENDIF
C
      RETURN
C----------------------- End of Subroutine parmsetget -------------------------
      END
