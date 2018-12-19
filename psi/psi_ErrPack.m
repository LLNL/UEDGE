      SUBROUTINE ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error
c        after making symbolic dump (machine-specific)

c     .. Scalar Arguments ..

      CHARACTER MESSAG*(*)
      LOGICAL   FATAL
c     ..
c     .. Local Scalars ..

      LOGICAL   MSGLIM
      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

ccccc EXTERNAL  SYMDUMP
c     ..
      SAVE      MAXMSG, NUMMSG, MSGLIM
      DATA      NUMMSG / 0 /,  MAXMSG / 100 /,  MSGLIM /.FALSE./


      IF( FATAL ) THEN

         WRITE( *, '(//,2A,//)' ) ' ****** ERROR *****  ', MESSAG

c                                 ** Example symbolic dump call for Cray
ccccc    CALL SYMDUMP( '-B -c3' )

         STOP

      END IF


      NUMMSG = NUMMSG + 1

      IF( MSGLIM ) RETURN

      IF( NUMMSG.LE.MAXMSG ) THEN

         WRITE( *, '(/,2A,/)' ) ' ****** WARNING *****  ', MESSAG

      ELSE

         WRITE( *, 9000 )
         MSGLIM = .True.
      END IF


      RETURN

 9000 FORMAT( / , / , ' ****** TOO MANY WARNING MESSAGES --  ',
     &      'They will no longer be printed *******', / , / )
      END

      LOGICAL FUNCTION WrtBad( VarNam )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

c     .. Scalar Arguments ..

      CHARACTER VarNam*(*)
c     ..
c     .. Local Scalars ..

      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

      EXTERNAL  ErrMsg
c     ..
      SAVE      NUMMSG, MAXMSG
      DATA      NUMMSG / 0 /, MAXMSG / 50 /


      WrtBad = .TRUE.
      NUMMSG = NUMMSG + 1
      WRITE( *, '(3A)' ) ' ****  Input variable  ', VarNam,
     &                   '  in error  ****'

      IF( NUMMSG.EQ.MAXMSG )
     &    CALL ErrMsg( 'Too many input errors.  Aborting...',.TRUE.)

      RETURN
      END

      LOGICAL FUNCTION WrtDim( DimNam, Minval )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

c     .. Scalar Arguments ..

      CHARACTER DimNam*(*)
      INTEGER   Minval
c     ..

      WRITE( *, '(3A,I7)' ) ' ****  Symbolic dimension  ', DimNam,
     &                      '  should be increased to at least ', Minval
      WrtDim = .TRUE.
      RETURN
      END

      LOGICAL FUNCTION TstBad( VarNam, RelErr )

c       Write name (VarNam) of variable failing self-test and its
c       percent error from the correct value;  return  'FALSE'.

c     .. Scalar Arguments ..

      CHARACTER VarNam*(*)
      REAL      RelErr
c     ..

      TstBad = .FALSE.
      WRITE( *, '(/,3A,1P,E11.2,A)' ) ' *** Output variable ', VarNam,
     &   ' differed by ', 100.*RelErr,
     &   ' per cent from correct value.  Self-test failed.'
      RETURN
      END
