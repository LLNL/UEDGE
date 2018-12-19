      SUBROUTINE MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU,
     &                  NMOM, IPOLZN, MOMDIM, PRNT, QEXT, QSCA, GQSC,
     &                  PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK,
     &                  SPIKE )

c    Computes Mie scattering and extinction efficiencies; asymmetry
c    factor;  forward- and backscatter amplitude;  scattering
c    amplitudes vs. scattering angle for incident polarization parallel
c    and perpendicular to the plane of scattering;
c    coefficients in the Legendre polynomial expansions of either the
c    unpolarized phase function or the polarized phase matrix;
c    some quantities needed in polarized radiative transfer;  and
c    information about whether or not a resonance has been hit.

c    Input and output variables are described in file MIEV.DOC

c    CALLING TREE:

c        MIEV0
c            TESTMI
c                TSTBAD
c                MIPRNT
c                ERRMSG
c            CKINMI
c                WRTBAD
c                WRTDIM
c                ERRMSG
c            SMALL1
c            SMALL2
c            ERRMSG
c            BIGA
c                CONFRA
c                    ERRMSG
c            LPCOEF
c                LPCO1T
c                LPCO2T
c                ERRMSG
c            MIPRNT


c      I N T E R N A L   V A R I A B L E S
c      -----------------------------------

c  AN,BN           Mie coefficients a-sub-n, b-sub-n ( Ref. 1, Eq. 16 )
c  ANM1,BNM1       Mie coefficients  a-sub-(n-1),
c                     b-sub-(n-1);  used in GQSC sum
c  ANP             Coeffs. in S+ expansion ( Ref. 2, p. 1507 )
c  BNP             Coeffs. in S- expansion ( Ref. 2, p. 1507 )
c  ANPM            Coeffs. in S+ expansion ( Ref. 2, p. 1507 )
c                     when  MU  is replaced by  - MU
c  BNPM            Coeffs. in S- expansion ( Ref. 2, p. 1507 )
c                     when  MU  is replaced by  - MU
c  CALCMO(K)       TRUE, calculate moments for K-th phase quantity
c                     (derived from IPOLZN)
c  CBIGA(N)        Bessel function ratio A-sub-N (Ref. 2, Eq. 2)
c                     ( COMPLEX version )
c  CDENAN,         (COMPLEX) denominators of An,Bn
c   CDENBN
c  CIOR            Complex index of refraction with negative
c                     imaginary part (Van de Hulst convention)
c  CIORIV          1 / cIoR
c  COEFF           ( 2N + 1 ) / ( N ( N + 1 ) )
c  FN              Floating point version of index in loop performing
c                     Mie series summation
c  LITA,LITB(N)    Mie coefficients An, Bn, saved in arrays for
c                     use in calculating Legendre moments PMOM
c  MAXTRM          Max. possible no. of terms in Mie series
c  MM              + 1 and  - 1,  alternately.
c  MIM             Magnitude of imaginary refractive index
c  MRE             Real part of refractive index
c  MAXANG          Max. possible value of input variable NUMANG
c  NANGD2          (NUMANG+1)/2 ( no. of angles in 0-90 deg; ANYANG=F )
c  NOABS           TRUE, sphere non-absorbing (determined by MIMCUT)
c  NP1DN           ( N + 1 ) / N
c  NPQUAN          Highest-numbered phase quantity for which moments are
c                     to be calculated (the largest digit in IPOLZN
c                     if  IPOLZN .NE. 0)
c  NTRM            No. of terms in Mie series
c  PASS1           TRUE on first entry, FALSE thereafter; for self-test
c  PIN(J)          Angular function pi-sub-n ( Ref. 2, Eq. 3 )
c                     at J-th angle
c  PINM1(J)        pi-sub-(n-1) ( see PIn ) at J-th angle
c  PSINM1          Ricatti-Bessel function psi-sub-(n-1), argument XX
c  PSIN            Ricatti-Bessel function psi-sub-n of argument XX
c                     ( Ref. 1, p. 11 ff. )
c  RBIGA(N)        Bessel function ratio A-sub-N (Ref. 2, Eq. 2)
c                     ( REAL version, for when imag refrac index = 0 )
c  RIORIV          1 / Mre
c  RN              1 / N
c  RTMP            (REAL) temporary variable
c  SP(J)           S+  for J-th angle  ( Ref. 2, p. 1507 )
c  SM(J)           S-  for J-TH angle  ( Ref. 2, p. 1507 )
c  SPS(J)          S+  for (NUMANG+1-J)-th angle ( ANYANG=FALSE )
c  SMS(J)          S-  for (NUMANG+1-J)-th angle ( ANYANG=FALSE )
c  TAUN            Angular function tau-sub-n ( Ref. 2, Eq. 4 )
c                     at J-th angle
c  TCOEF           N ( N+1 ) ( 2N+1 ) (for summing TFORW,TBACK series)
c  TWONP1          2N + 1
c  YESANG          TRUE if scattering amplitudes are to be calculated
c  ZETNM1          Ricatti-Bessel function  zeta-sub-(n-1) of argument
c                     XX  ( Ref. 2, Eq. 17 )
c  ZETN            Ricatti-Bessel function  zeta-sub-n of argument XX


      IMPLICIT  NONE

c ----------------------------------------------------------------------
c --------  I / O SPECIFICATIONS FOR SUBROUTINE MIEV0  -----------------
c ----------------------------------------------------------------------
      LOGICAL  ANYANG, PERFCT, PRNT(*)
      INTEGER  IPOLZN, MOMDIM, NUMANG, NMOM
      REAL*8   GQSC, MIMCUT, PMOM( 0:MOMDIM, * ), QEXT, QSCA, SPIKE,
     &         XMU(*), XX
      COMPLEX*16  CREFIN, SFORW, SBACK, S1(*), S2(*), TFORW(*), TBACK(*)
c ----------------------------------------------------------------------

c                                  ** NOTE --  MAXTRM = 10100  is neces-
c                                  ** sary to do some of the test probs,
c                                  ** but 1100 is sufficient for most
c                                  ** conceivable applications
c     .. Parameters ..

      INTEGER   MAXANG, MXANG2
      PARAMETER ( MAXANG = 501, MXANG2 = MAXANG / 2 + 1 )
      INTEGER   MAXTRM
      PARAMETER ( MAXTRM = 100100 )
      REAL*8    ONETHR
      PARAMETER ( ONETHR = 1.D0 / 3.D0 )
c     ..
c     .. Local Scalars ..

      LOGICAL   NOABS, PASS1, YESANG
      INTEGER   I, J, N, NANGD2, NPQUAN, NTRM
      REAL*8    CHIN, CHINM1, COEFF, DENAN, DENBN, FN, MIM, MM, MRE,
     &          NP1DN, PSIN, PSINM1, RATIO, RIORIV, RN, RTMP, TAUN,
     &          TCOEF, TWONP1, XINV
      COMPLEX*16   AN, ANM1, ANP, ANPM, BN, BNM1, BNP, BNPM, CDENAN,
     &                 CDENBN, CIOR, CIORIV, CTMP, ZET, ZETN, ZETNM1
c     ..
c     .. Local Arrays ..

      LOGICAL   CALCMO( 4 )
      REAL*8    PIN( MAXANG ), PINM1( MAXANG ), RBIGA( MAXTRM )
      COMPLEX*16   CBIGA( MAXTRM ), LITA( MAXTRM ), LITB( MAXTRM ),
     &          SM( MAXANG ), SMS( MXANG2 ), SP( MAXANG ), SPS( MXANG2 )
c     ..
c     .. External Subroutines ..

      EXTERNAL  BIGA, CKINMI, ERRMSG, LPCOEF, MIPRNT, SMALL1, SMALL2,
     &          TESTMI
c     ..
      SAVE      PASS1

c     .. Statement Functions ..

      REAL*8      SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = DREAL( CTMP )**2 + DIMAG( CTMP )**2
c     ..
      DATA      PASS1 / .TRUE. /


c                    ** Save some input variables and replace them
c                    ** with values needed to do the self-test
      
      IF( PASS1 ) CALL TESTMI( .FALSE., XX, CREFIN, MIMCUT, PERFCT,
     &                         ANYANG, NMOM, IPOLZN, NUMANG, XMU, QEXT,
     &                         QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW,
     &                         TBACK, PMOM, MOMDIM )
      
   10 CONTINUE
c                                        ** Check input and calculate
c                                        ** certain variables from input

      CALL CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, MOMDIM, NMOM,
     &             IPOLZN, ANYANG, XMU, CALCMO, NPQUAN )


      IF( PERFCT .AND. XX.LE.0.1D0 ) THEN
c                                            ** Use totally-reflecting
c                                            ** small-particle limit

         CALL SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK,
     &                S1, S2, TFORW, TBACK, LITA, LITB )

         NTRM = 2
         GO TO  100

      END IF


c      NOABS = .TRUE.
      NOABS = .FALSE.

      IF( .NOT.PERFCT ) THEN
      
         CIOR = CREFIN

         IF( DIMAG(CIOR).GT.0.D0 ) CIOR = DCONJG( CIOR )

         MRE    = DREAL( CIOR )
         MIM    = -DIMAG( CIOR )
         NOABS  = MIM.LE.MIMCUT
         CIORIV = 1.D0 / CIOR
         RIORIV = 1.D0 / MRE
      
         IF( XX*MAX( 1.D0, CDABS(CIOR) ).LE.0.1D0 ) THEN

c                                    ** Use general-refractive-index
c                                    ** small-particle limit
c                                    ** ( Ref. 2, p. 1508 )
      
            CALL SMALL2( XX, CIOR, .NOT.NOABS, NUMANG, XMU, QEXT, QSCA,
     &                   GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK, LITA,
     &                   LITB )
      
            NTRM = 2
            GO TO  100
      
         END IF
      
      END IF


      NANGD2 = ( NUMANG + 1 ) / 2
      YESANG = NUMANG.GT.0

c                              ** Estimate number of terms in Mie series
c                              ** ( Ref. 2, p. 1508 )
      IF( XX.LE.8.D0 ) THEN

         NTRM = XX + 4.D0*XX**ONETHR + 1.D0

      ELSE IF( XX.LT.4200.D0 ) THEN

         NTRM = XX + 4.05D0*XX**ONETHR + 2.D0

      ELSE

         NTRM = XX + 4.D0*XX**ONETHR + 2.D0

      END IF

      IF( NTRM+1 .GT. MAXTRM )
     &    CALL ERRMSG('MIEV0--PARAMETER MaxTrm TOO SMALL',.TRUE.)

c                            ** Calculate logarithmic derivatives of
c                            ** J-Bessel-fcn., A-sub-(1 to NTrm)

      IF( .NOT.PERFCT ) CALL BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA,
     &                             CBIGA )

c                            ** Initialize Ricatti-Bessel functions
c                            ** (psi,chi,zeta)-sub-(0,1) for upward
c                            ** recurrence ( Ref. 1, Eq. 19 )
      XINV   = 1.D0 / XX
      PSINM1 = DSIN( XX )
      CHINM1 = DCOS( XX )
      PSIN   = PSINM1*XINV - CHINM1
      CHIN   = CHINM1*XINV + PSINM1
      ZETNM1 = DCMPLX( PSINM1, CHINM1 )
      ZETN   = DCMPLX( PSIN, CHIN )
c                                     ** Initialize previous coeffi-
c                                     ** cients for -GQSC- series
      ANM1 = ( 0.D0, 0.D0 )
      BNM1 = ( 0.D0, 0.D0 )
c                             ** Initialize angular function  pi
c                             ** and sums for S+, S- ( Ref. 2, p. 1507 )
      IF( ANYANG ) THEN

         DO 20 J = 1, NUMANG
            PINM1( J ) = 0.D0
            PIN  ( J ) = 1.D0
            SP   ( J ) = ( 0.D0, 0.D0 )
            SM   ( J ) = ( 0.D0, 0.D0 )
   20    CONTINUE

      ELSE

         DO 30 J = 1, NANGD2
            PINM1( J ) = 0.D0
            PIN  ( J ) = 1.D0
            SP   ( J ) = ( 0.D0, 0.D0 )
            SM   ( J ) = ( 0.D0, 0.D0 )
            SPS  ( J ) = ( 0.D0, 0.D0 )
            SMS  ( J ) = ( 0.D0, 0.D0 )
   30    CONTINUE

      END IF

c                         ** Initialize Mie sums for efficiencies, etc.
      QSCA   = 0.D0
      GQSC   = 0.D0
      SFORW  = ( 0.D0, 0.D0 )
      SBACK  = ( 0.D0, 0.D0 )
      TFORW( 1 ) = ( 0.D0, 0.D0 )
      TBACK( 1 ) = ( 0.D0, 0.D0 )
c ---------  LOOP TO SUM MIE SERIES  -----------------------------------
      
      MM     = 1.D0
      SPIKE  = 1.D0
      
      DO 60  N = 1, NTRM
c                           ** Compute various numerical coefficients
         FN     = N
         RN     = 1.D0 / FN
         NP1DN  = 1.D0 + RN
         TWONP1 = 2*N + 1
         COEFF  = TWONP1 / ( FN * ( N + 1 ) )
         TCOEF  = TWONP1 * ( FN * ( N + 1 ) )
c                              ** Calculate Mie series coefficients
         IF( PERFCT ) THEN
c                                   ** Totally-reflecting case
            AN = ( ( FN*XINV )*PSIN - PSINM1 ) /
     &           ( ( FN*XINV )*ZETN - ZETNM1 )
            BN = PSIN / ZETN

         ELSE IF( NOABS ) THEN
c                                      ** No-absorption case

            CDENAN = ( RIORIV*RBIGA(N) + ( FN*XINV ) ) * ZETN - ZETNM1
            AN   = ( ( RIORIV*RBIGA(N) + ( FN*XINV ) ) * PSIN - PSINM1 )
     &             / CDENAN
            CDENBN = ( MRE*RBIGA(N) + ( FN*XINV ) ) * ZETN - ZETNM1
            BN   = ( ( MRE*RBIGA(N) + ( FN*XINV ) ) * PSIN - PSINM1 )
     &             / CDENBN

         ELSE
c                                       ** Absorptive case

            CDENAN = ( CIORIV*CBIGA( N ) + ( FN*XINV ) )*ZETN - ZETNM1
            CDENBN =   ( CIOR*CBIGA( N ) + ( FN*XINV ) )*ZETN - ZETNM1
            AN   = ( ( CIORIV*CBIGA( N ) + ( FN*XINV ) )*PSIN - PSINM1 )
     &             / CDENAN
            BN     = ( ( CIOR*CBIGA( N ) + ( FN*XINV ) )*PSIN - PSINM1 )
     &             / CDENBN

            QSCA   = QSCA + TWONP1*( SQ( AN ) + SQ( BN ) )

         END IF
c                       ** Save Mie coefficients for *PMOM* calculation

         LITA( N ) = AN
         LITB( N ) = BN


         IF( .NOT.PERFCT .AND. N.GT.XX ) THEN
c                                               ** Flag resonance spikes
            DENAN  = CDABS( CDENAN )
            DENBN  = CDABS( CDENBN )
            RATIO  = DENAN / DENBN

            IF( RATIO.LE.0.2D0 .OR. RATIO.GE.5.D0 )
     &          SPIKE = DMIN1( SPIKE, DENAN, DENBN )

         END IF
c                                  ** Increment Mie sums for non-angle-
c                                  ** dependent quantities

         SFORW      = SFORW      + TWONP1*( AN + BN )
         TFORW( 1 ) = TFORW( 1 ) + TCOEF *( AN - BN )
         SBACK      = SBACK      + ( MM*TWONP1 )*( AN - BN )
         TBACK( 1 ) = TBACK( 1 ) + ( MM*TCOEF ) *( AN + BN )
         GQSC       = GQSC + ( FN - RN ) * DREAL( ANM1*DCONJG( AN )
     &              + BNM1*DCONJG( BN ) )
     &              + COEFF * DREAL( AN*DCONJG( BN ) )


         IF( YESANG ) THEN
c                                      ** Put Mie coefficients in form
c                                      ** needed for computing S+, S-
c                                      ** ( Ref. 2, p. 1507 )
            ANP= COEFF*( AN + BN )
            BNP= COEFF*( AN - BN )

c                                      ** Increment Mie sums for S+, S-
c                                      ** while upward recursing
c                                      ** angular functions pi and tau
            IF( ANYANG ) THEN
c                                         ** Arbitrary angles

c                                              ** vectorizable loop
               DO 40 J = 1, NUMANG
                  RTMP   = ( XMU( J ) * PIN( J ) ) - PINM1( J )
                  TAUN   = FN * RTMP - PINM1( J )
                  SP( J ) = SP( J ) + ANP * ( PIN( J ) + TAUN )
                  SM( J ) = SM( J ) + BNP * ( PIN( J ) - TAUN )
                  PINM1( J ) = PIN( J )
                  PIN( J ) = ( XMU( J ) * PIN( J ) ) + NP1DN * RTMP
   40          CONTINUE

            ELSE
c                                  ** Angles symmetric about 90 degrees
               ANPM = MM*ANP
               BNPM = MM*BNP
c                                          ** vectorizable loop
               DO 50 J = 1, NANGD2
                  RTMP   = ( XMU( J ) * PIN( J ) ) - PINM1( J )
                  TAUN   = FN * RTMP - PINM1( J )
                  SP ( J ) = SP ( J ) + ANP * ( PIN( J ) + TAUN )
                  SMS( J ) = SMS( J ) + BNPM *( PIN( J ) + TAUN )
                  SM ( J ) = SM ( J ) + BNP * ( PIN( J ) - TAUN )
                  SPS( J ) = SPS( J ) + ANPM *( PIN( J ) - TAUN )
                  PINM1( J ) = PIN( J )
                  PIN( J ) = ( XMU( J ) * PIN( J ) ) + NP1DN * RTMP
   50          CONTINUE

            END IF

         END IF
c                          ** Update relevant quantities for next
c                          ** pass through loop

         MM   = - MM
         ANM1 = AN
         BNM1 = BN
c                           ** Upward recurrence for Ricatti-Bessel
c                           ** functions ( Ref. 1, Eq. 17 )

         ZET    = ( TWONP1*XINV ) * ZETN - ZETNM1
         ZETNM1 = ZETN
         ZETN   = ZET
         PSINM1 = PSIN
         PSIN   = DREAL( ZETN )

   60 CONTINUE

c ---------- END LOOP TO SUM MIE SERIES --------------------------------

      QEXT = 2.D0 / XX**2 * DREAL( SFORW )

      IF( PERFCT .OR. NOABS ) THEN

         QSCA = QEXT

      ELSE

         QSCA = 2.D0/ XX**2 * QSCA

      END IF

      GQSC   = 4.D0/ XX**2 * GQSC
      SFORW  = 0.5D0*SFORW
      SBACK  = 0.5D0*SBACK
      TFORW( 2 ) = 0.5D0*(  SFORW + 0.25D0*TFORW( 1 ) )
      TFORW( 1 ) = 0.5D0*(  SFORW - 0.25D0*TFORW( 1 ) )
      TBACK( 2 ) = 0.5D0*(  SBACK + 0.25D0*TBACK( 1 ) )
      TBACK( 1 ) = 0.5D0*( -SBACK + 0.25D0*TBACK( 1 ) )


      IF( YESANG ) THEN
c                                ** Recover scattering amplitudes
c                                ** from S+, S- ( Ref. 1, Eq. 11 )

         IF( ANYANG ) THEN
c                                         ** vectorizable loop
            DO 70 J = 1, NUMANG
               S1( J ) = 0.5D0*( SP( J ) + SM( J ) )
               S2( J ) = 0.5D0*( SP( J ) - SM( J ) )
   70       CONTINUE

         ELSE
c                                         ** vectorizable loop
            DO 80 J = 1, NANGD2
               S1( J ) = 0.5D0*( SP( J ) + SM( J ) )
               S2( J ) = 0.5D0*( SP( J ) - SM( J ) )
   80       CONTINUE
c                                         ** vectorizable loop
            DO 90 J = 1, NANGD2
               S1( NUMANG + 1 - J ) = 0.5D0*( SPS( J ) + SMS( J ) )
               S2( NUMANG + 1 - J ) = 0.5D0*( SPS( J ) - SMS( J ) )
   90       CONTINUE

         END IF

      END IF
c                                         ** Calculate Legendre moments

  100 CONTINUE
      IF( NMOM.GT.0 ) CALL LPCOEF( NTRM, NMOM, IPOLZN, MOMDIM, CALCMO,
     &                             NPQUAN, LITA, LITB, PMOM )


      IF( DIMAG( CREFIN ).GT.0.D0 ) THEN
c                                         ** Take complex conjugates
c                                         ** of scattering amplitudes

         SFORW = DCONJG( SFORW )
         SBACK = DCONJG( SBACK )

         DO 110 I = 1, 2
            TFORW( I ) = DCONJG( TFORW( I ) )
            TBACK( I ) = DCONJG( TBACK( I ) )
  110    CONTINUE

         DO 120 J = 1, NUMANG
            S1( J ) = DCONJG( S1( J ) )
            S2( J ) = DCONJG( S2( J ) )
  120    CONTINUE

      END IF


      IF( PASS1 ) THEN
c                           ** Compare test case results with
c                           ** correct answers and abort if bad;
c                           ** otherwise restore user input and proceed

         CALL TESTMI( .TRUE., XX, CREFIN, MIMCUT, PERFCT, ANYANG, NMOM,
     &                IPOLZN, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW,
     &                SBACK, S1, S2, TFORW, TBACK, PMOM, MOMDIM )

         PASS1  = .FALSE.
         GO TO  10

      END IF


      IF( PRNT( 1 ) .OR. PRNT( 2 ) ) CALL MIPRNT( PRNT, XX, PERFCT,
     &    CREFIN, NUMANG, XMU, QEXT, QSCA, GQSC, NMOM, IPOLZN, MOMDIM,
     &    CALCMO, PMOM, SFORW, SBACK, TFORW, TBACK, S1, S2 )

      RETURN

      END

      SUBROUTINE CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, MOMDIM,
     &                   NMOM, IPOLZN, ANYANG, XMU, CALCMO, NPQUAN )

c        Check for bad input to MIEV0 and calculate CALCMO, NPQUAN

c     Routines called :  ERRMSG, WRTBAD, WRTDIM


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   ANYANG, PERFCT
      INTEGER   IPOLZN, MAXANG, MOMDIM, NMOM, NPQUAN, NUMANG
      REAL*8    XX
      COMPLEX*16   CREFIN
c     ..
c     .. Array Arguments ..

      LOGICAL   CALCMO( * )
      REAL*8    XMU( * )
c     ..
c     .. Local Scalars ..

      CHARACTER STRING*4
      LOGICAL   INPERR
      INTEGER   I, IP, J, L
c     ..
c     .. External Functions ..

      LOGICAL   WRTBAD, WRTDIM
      EXTERNAL  WRTBAD, WRTDIM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..

      INPERR = .FALSE.

      IF( NUMANG.GT.MAXANG ) INPERR = WRTDIM( 'MaxAng', NUMANG )
      IF( NUMANG.LT.0 ) INPERR = WRTBAD( 'NUMANG' )

      IF( XX.LT.0.D0 ) INPERR = WRTBAD( 'XX' )

      IF( .NOT.PERFCT .AND. DREAL( CREFIN ).LE.0.D0 )
     &    INPERR = WRTBAD( 'CREFIN' )

      IF( MOMDIM.LT.1 ) INPERR = WRTBAD( 'MOMDIM' )


      IF( NMOM.NE.0 ) THEN

         IF( NMOM.LT.0 .OR. NMOM.GT.MOMDIM ) INPERR = WRTBAD( 'NMOM' )

         IF( ABS( IPOLZN ).GT.4444 ) INPERR = WRTBAD( 'IPOLZN' )

         NPQUAN = 0

         DO 10 L = 1, 4
            CALCMO( L ) = .FALSE.
   10    CONTINUE

         IF( IPOLZN.NE.0 ) THEN
c                                 ** Parse out IPOLZN into its digits
c                                 ** to find which phase quantities are
c                                 ** to have their moments calculated

            WRITE( STRING, '(I4)' ) ABS( IPOLZN )

            DO 20 J = 1, 4
               IP = ICHAR( STRING( J:J ) ) - ICHAR( '0' )

               IF( IP.GE.1 .AND. IP.LE.4 ) CALCMO( IP ) = .TRUE.

               IF( IP.EQ.0 .OR. ( IP.GE.5 .AND. IP.LE.9 ) )
     &             INPERR = WRTBAD( 'IPOLZN' )

               NPQUAN = MAX( NPQUAN, IP )
   20       CONTINUE

         END IF

      END IF


      IF( ANYANG ) THEN
c                                ** Allow for slight imperfections in
c                                ** computation of cosine
         DO 30 I = 1, NUMANG

            IF( XMU( I ).LT.-1.00001D0 .OR. XMU( I ).GT.1.00001D0 )
     &          INPERR = WRTBAD( 'XMU' )

   30    CONTINUE

      ELSE

         DO 40 I = 1, ( NUMANG + 1 ) / 2

            IF( XMU( I ).LT.-0.00001D0 .OR. XMU( I ).GT.1.00001D0 )
     &          INPERR = WRTBAD( 'XMU' )

   40    CONTINUE

      END IF


      IF( INPERR ) CALL ERRMSG( 'MIEV0--Input error(S).  Aborting...',
     &                          .TRUE. )

c      IF( XX.GT.20000.D0 .OR. DREAL( CREFIN ).GT.10.D0 .OR.
c     &    DABS( DIMAG( CREFIN ) ).GT.10.D0 )
c     &    CALL ERRMSG( 'MIEV0--XX or CREFIN outside tested range',
c     &    .FALSE.)

      RETURN
      END

      SUBROUTINE LPCOEF( NTRM, NMOM, IPOLZN, MOMDIM, CALCMO, NPQUAN, A,
     &                   B, PMOM )

c         Calculate Legendre polynomial expansion coefficients (also
c         called moments) for phase quantities ( Ref. 5 formulation )

c     INPUT:  NTRM                    Number terms in Mie series
c             NMOM, IPOLZN, MOMDIM    MIEV0 arguments
c             CALCMO                  Flags calculated from IPOLZN
c             NPQUAN                  Defined in MIEV0
c             A, B                    Mie series coefficients

c     OUTPUT: PMOM                   Legendre moments (MIEV0 argument)

c     Routines called :  ERRMSG, LPCO1T, LPCO2T

c     *** NOTES ***

c         (1)  Eqs. 2-5 are in error in Dave, Appl. Opt. 9,
c         1888 (1970).  Eq. 2 refers to M1, not M2;  eq. 3 refers to
c         M2, not M1.  In eqs. 4 and 5, the subscripts on the second
c         term in square brackets should be interchanged.

c         (2)  The general-case logic in this subroutine works correctly
c         in the two-term Mie series case, but subroutine LPCO2T
c         is called instead, for speed.

c         (3)  Subroutine  LPCO1T, to do the one-term case, is never
c         called within the context of MIEV0, but is included for
c         complete generality.

c         (4)  Some improvement in speed is obtainable by combining the
c         310- and 410-loops, if moments for both the third and fourth
c         phase quantities are desired, because the third phase quantity
c         is the real part of a complex series, while the fourth phase
c         quantity is the imaginary part of that very same series.  But
c         most users are not interested in the fourth phase quantity,
c         which is related to circular polarization, so the present
c         scheme is usually more efficient.


c           ** Definitions of local variables ***

c      AM(M)       Numerical coefficients  a-sub-m-super-l
c                     in Dave, Eqs. 1-15, as simplified in Ref. 5.

c      BI(I)       Numerical coefficients  b-sub-i-super-l
c                     in Dave, Eqs. 1-15, as simplified in Ref. 5.

c      BIDEL(I)    1/2 Bi(I) times factor capital-del in Dave

c      CM,DM()     Arrays C and D in Dave, Eqs. 16-17 (Mueller form),
c                     calculated using recurrence derived in Ref. 5

c      CS,DS()     Arrays C and D in Ref. 4, Eqs. A5-A6 (Sekera form),
c                     calculated using recurrence derived in Ref. 5

c      C,D()       Either CM,DM or CS,DS, depending on IPOLZN

c      EVENL       True for even-numbered moments;  false otherwise

c      IDEL        1 + little-del  in Dave

c      MAXTRM      Max. no. of terms in Mie series

c      MAXMOM      Max. no. of non-zero moments

c      NUMMOM      Number of non-zero moments

c      RECIP(K)    1 / K


      IMPLICIT  NONE

c     .. Parameters ..

      INTEGER   MAXTRM, MAXMOM, MXMOM2, MAXRCP
      PARAMETER ( MAXTRM = 100102, MAXMOM = 2*MAXTRM, MXMOM2 = MAXMOM/2,
     &            MAXRCP = 4*MAXTRM + 2 )
c     ..
c     .. Scalar Arguments ..

      INTEGER   IPOLZN, MOMDIM, NMOM, NPQUAN, NTRM
c     ..
c     .. Array Arguments ..

      LOGICAL   CALCMO( * )
      REAL*8    PMOM( 0:MOMDIM, * )
      COMPLEX*16   A( * ), B( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   EVENL, PASS1
      INTEGER   I, IDEL, IMAX, J, K, L, LD2, M, MMAX, NUMMOM
      REAL*8    SUM
c     ..
c     .. Local Arrays ..

      REAL*8    AM( 0:MAXTRM ), BI( 0:MXMOM2 ), BIDEL( 0:MXMOM2 ),
     &          RECIP( MAXRCP )
      COMPLEX*16   C( MAXTRM ), CM( MAXTRM ), CS( MAXTRM ),
     &          D( MAXTRM ), DM( MAXTRM ), DS( MAXTRM )
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, LPCO1T, LPCO2T
c     ..
c     .. Equivalences ..

      EQUIVALENCE ( C, CM ), ( D, DM )
c     ..
      SAVE      PASS1, RECIP

      DATA      PASS1 / .TRUE. /


      IF( PASS1 ) THEN

         DO 10 K = 1, MAXRCP
            RECIP( K ) = 1.D0 / K
   10    CONTINUE

         PASS1  = .FALSE.

      END IF


      DO 30 J = 1, MAX( 1, NPQUAN )

         DO 20 L = 0, NMOM
            PMOM( L, J ) = 0.D0
   20    CONTINUE

   30 CONTINUE


      IF( NTRM.EQ.1 ) THEN

         CALL LPCO1T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

         RETURN

      ELSE IF( NTRM.EQ.2 ) THEN

         CALL LPCO2T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

         RETURN

      END IF


      IF( NTRM + 2.GT.MAXTRM )
     &    CALL ERRMSG('LPCoef--PARAMETER MaxTrm too small',.TRUE.)

c                                     ** Calculate Mueller C, D arrays
      CM( NTRM + 2 ) = ( 0.D0, 0.D0 )
      DM( NTRM + 2 ) = ( 0.D0, 0.D0 )
      CM( NTRM + 1 ) = ( 1.D0 - RECIP( NTRM+1 ) ) * B( NTRM )
      DM( NTRM + 1 ) = ( 1.D0 - RECIP( NTRM+1 ) ) * A( NTRM )
      CM( NTRM ) = ( RECIP( NTRM ) + RECIP( NTRM+1 ) ) * A( NTRM ) +
     &             ( 1.D0 - RECIP( NTRM ) )*B( NTRM-1 )
      DM( NTRM ) = ( RECIP( NTRM ) + RECIP( NTRM+1 ) ) * B( NTRM ) +
     &             ( 1.D0 - RECIP( NTRM ) )*A( NTRM-1 )

      DO 40 K = NTRM-1, 2, -1
         CM( K ) = CM( K+2 ) - ( 1.D0 + RECIP(K+1) ) * B( K+1 )
     $                       + ( RECIP(K) + RECIP(K+1) ) * A( K )
     $                       + ( 1.D0 - RECIP(K) ) * B( K-1 )
         DM( K ) = DM( K+2 ) - ( 1.D0 + RECIP(K+1) ) * A( K+1 )
     $                       + ( RECIP(K) + RECIP(K+1) ) * B( K )
     $                       + ( 1.D0 - RECIP(K) ) * A( K-1 )
   40 CONTINUE

      CM( 1 ) = CM( 3 ) + 1.5D0 * ( A( 1 ) - B( 2 ) )
      DM( 1 ) = DM( 3 ) + 1.5D0 * ( B( 1 ) - A( 2 ) )


      IF( IPOLZN.GE.0 ) THEN

         DO 50 K = 1, NTRM + 2
            C( K ) = ( 2*K - 1 ) * CM( K )
            D( K ) = ( 2*K - 1 ) * DM( K )
   50    CONTINUE

      ELSE
c                                    ** Compute Sekera C and D arrays
         CS( NTRM + 2 ) = ( 0.D0, 0.D0 )
         DS( NTRM + 2 ) = ( 0.D0, 0.D0 )
         CS( NTRM + 1 ) = ( 0.D0, 0.D0 )
         DS( NTRM + 1 ) = ( 0.D0, 0.D0 )

         DO 60 K = NTRM, 1, -1
            CS( K ) = CS( K+2 ) + ( 2*K + 1 ) * ( CM( K+1 ) - B( K ) )
            DS( K ) = DS( K+2 ) + ( 2*K + 1 ) * ( DM( K+1 ) - A( K ) )
   60    CONTINUE

         DO 70 K = 1, NTRM + 2
            C( K ) = ( 2*K - 1 ) * CS( K )
            D( K ) = ( 2*K - 1 ) * DS( K )
   70    CONTINUE

      END IF


      IF( IPOLZN.LT.0 ) NUMMOM = MIN( NMOM, 2*NTRM - 2 )
      IF( IPOLZN.GE.0 ) NUMMOM = MIN( NMOM, 2*NTRM )

      IF( NUMMOM.GT.MAXMOM )
     &    CALL ERRMSG('LPCoef--PARAMETER MaxTrm too small',.TRUE.)


c                          ** Loop over moments

      DO 240 L = 0, NUMMOM

         LD2 = L / 2
         EVENL  = MOD( L, 2 ).EQ.0
c                                    ** Calculate numerical coefficients
c                                    ** a-sub-m and b-sub-i in Dave
c                                    ** double-sums for moments
         IF( L.EQ.0 ) THEN

            IDEL = 1

            DO 80 M = 0, NTRM
               AM( M ) = 2.D0 * RECIP( 2*M + 1 )
   80       CONTINUE

            BI( 0 ) = 1.D0

         ELSE IF( EVENL ) THEN

            IDEL = 1

            DO 90 M = LD2, NTRM
               AM( M ) = ( 1.D0 + RECIP( 2*M - L + 1 ) ) * AM( M )
   90       CONTINUE

            DO 100 I = 0, LD2 - 1
               BI( I ) = ( 1.D0 - RECIP( L - 2*I ) ) * BI( I )
  100       CONTINUE

            BI( LD2 ) = ( 2.D0 - RECIP( L ) ) * BI( LD2 - 1 )

         ELSE

            IDEL = 2

            DO 110 M = LD2, NTRM
               AM( M ) = ( 1.D0 - RECIP( 2*M + L + 2 ) ) * AM( M )
  110       CONTINUE

            DO 120 I = 0, LD2
               BI( I ) = ( 1.D0 - RECIP( L + 2*I + 1 ) ) * BI( I )
  120       CONTINUE

         END IF
c                                     ** Establish upper limits for sums
c                                     ** and incorporate factor capital-
c                                     ** del into b-sub-i
         MMAX = NTRM - IDEL
         IF( IPOLZN.GE.0 ) MMAX = MMAX + 1
         IMAX = MIN( LD2, MMAX - LD2 )

         IF( IMAX.LT.0 ) GO TO  250

         DO 130 I = 0, IMAX
            BIDEL( I ) = BI( I )
  130    CONTINUE

         IF( EVENL ) BIDEL( 0 ) = 0.5D0*BIDEL( 0 )

c                                    ** Perform double sums just for
c                                    ** phase quantities desired by user
         IF( IPOLZN.EQ.0 ) THEN

            DO 150 I = 0, IMAX
c                                           ** vectorizable loop

               SUM = 0.D0

               DO 140 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     $                     ( DREAL( C(M-I+1) * DCONJG( C(M+I+IDEL) ) )
     $                     + DREAL( D(M-I+1) * DCONJG( D(M+I+IDEL) ) ) )
  140          CONTINUE

               PMOM( L, 1 ) = PMOM( L, 1 ) + BIDEL( I ) * SUM

  150       CONTINUE

            PMOM( L, 1 ) = 0.5D0*PMOM( L, 1 )
            GO TO  240

         END IF


         IF( CALCMO( 1 ) ) THEN

            DO 170 I = 0, IMAX

               SUM = 0.D0
c                                           ** vectorizable loop
               DO 160 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     $                        DREAL( C(M-I+1) * DCONJG( C(M+I+IDEL) ) )
  160          CONTINUE

               PMOM( L, 1 ) = PMOM( L, 1 ) + BIDEL( I ) * SUM

  170       CONTINUE

         END IF


         IF( CALCMO( 2 ) ) THEN

            DO 190 I = 0, IMAX

               SUM = 0.D0
c                                           ** vectorizable loop
               DO 180 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     $                        DREAL( D(M-I+1) * DCONJG( D(M+I+IDEL) ) )
  180          CONTINUE

               PMOM( L, 2 ) = PMOM( L, 2 ) + BIDEL( I ) * SUM

  190       CONTINUE

         END IF


         IF( CALCMO( 3 ) ) THEN

            DO 210 I = 0, IMAX

               SUM = 0.D0
c                                           ** vectorizable loop
               DO 200 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     $                     ( DREAL( C(M-I+1) * DCONJG( D(M+I+IDEL) ) )
     $                     + DREAL( C(M+I+IDEL) * DCONJG( D(M-I+1) ) ) )
  200          CONTINUE

               PMOM( L, 3 ) = PMOM( L, 3 ) + BIDEL( I ) * SUM

  210       CONTINUE

            PMOM( L, 3 ) = 0.5D0*PMOM( L, 3 )
         END IF


         IF( CALCMO( 4 ) ) THEN

            DO 230 I = 0, IMAX

               SUM= 0.D0
c                                           ** vectorizable loop
               DO 220 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     $                     ( DIMAG( C(M-I+1) * DCONJG( D(M+I+IDEL) ) )
     $                     + DIMAG( C(M+I+IDEL) * DCONJG( D(M-I+1) ) ) )
  220          CONTINUE

               PMOM( L, 4 ) = PMOM( L, 4 ) + BIDEL( I ) * SUM

  230       CONTINUE

            PMOM( L, 4 ) = - 0.5D0 * PMOM( L, 4 )

         END IF

  240 CONTINUE


  250 CONTINUE

      RETURN
      END

      SUBROUTINE LPCO1T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

c         Calculate Legendre polynomial expansion coefficients (also
c         called moments) for phase quantities in special case where
c         no. terms in Mie series = 1

c        INPUT:  NMOM, IPOLZN, MOMDIM     MIEV0 arguments
c                CALCMO                   Flags calculated from IPOLZN
c                A(1), B(1)               Mie series coefficients

c        OUTPUT: PMOM                     Legendre moments


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      INTEGER   IPOLZN, MOMDIM, NMOM
c     ..
c     .. Array Arguments ..

      LOGICAL   CALCMO( * )
      REAL*8    PMOM( 0:MOMDIM, * )
      COMPLEX*16   A( * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER   L, NUMMOM
      REAL*8    A1SQ, B1SQ
      COMPLEX*16   A1B1C, CTMP
c     ..
c     .. Statement Functions ..

      REAL*8    SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = DREAL( CTMP )**2 + DIMAG( CTMP )**2
c     ..


      A1SQ   = SQ( A( 1 ) )
      B1SQ   = SQ( B( 1 ) )
      A1B1C  = A( 1 ) * DCONJG( B( 1 ) )


      IF( IPOLZN.LT.0 ) THEN

         IF( CALCMO( 1 ) ) PMOM( 0, 1 ) = 2.25D0*B1SQ

         IF( CALCMO( 2 ) ) PMOM( 0, 2 ) = 2.25D0*A1SQ

         IF( CALCMO( 3 ) ) PMOM( 0, 3 ) = 2.25D0*DREAL( A1B1C )

         IF( CALCMO( 4 ) ) PMOM( 0, 4 ) = 2.25D0*DIMAG( A1B1C )

      ELSE

         NUMMOM = MIN( NMOM, 2 )

c                             ** Loop over moments
         DO 10  L = 0, NUMMOM

            IF( IPOLZN.EQ.0 ) THEN

               IF( L.EQ.0 ) PMOM( L, 1 ) = 1.5D0*( A1SQ + B1SQ )

               IF( L.EQ.1 ) PMOM( L, 1 ) = 1.5D0*DREAL( A1B1C )

               IF( L.EQ.2 ) PMOM( L, 1 ) = 0.15D0*( A1SQ + B1SQ )

               GO TO  10

            END IF


            IF( CALCMO( 1 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 1 ) = 2.25D0*( A1SQ + B1SQ / 3.D0)

               IF( L.EQ.1 ) PMOM( L, 1 ) = 1.5D0*DREAL( A1B1C )

               IF( L.EQ.2 ) PMOM( L, 1 ) = 0.3D0*B1SQ

            END IF


            IF( CALCMO( 2 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 2 ) = 2.25D0*( B1SQ + A1SQ / 3.D0 )

               IF( L.EQ.1 ) PMOM( L, 2 ) = 1.5D0*DREAL( A1B1C )

               IF( L.EQ.2 ) PMOM( L, 2 ) = 0.3D0*A1SQ

            END IF


            IF( CALCMO( 3 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 3 ) = 3.0D0*DREAL( A1B1C )

               IF( L.EQ.1 ) PMOM( L, 3 ) = 0.75D0*( A1SQ + B1SQ )

               IF( L.EQ.2 ) PMOM( L, 3 ) = 0.3D0*DREAL( A1B1C )

            END IF


            IF( CALCMO( 4 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 4 ) = -1.5D0*DIMAG( A1B1C )

               IF( L.EQ.1 ) PMOM( L, 4 ) = 0.D0

               IF( L.EQ.2 ) PMOM( L, 4 ) = 0.3D0*DIMAG( A1B1C )

            END IF


   10    CONTINUE

      END IF

      RETURN
      END

      SUBROUTINE LPCO2T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

c         Calculate Legendre polynomial expansion coefficients (also
c         called moments) for phase quantities in special case where
c         no. terms in Mie series = 2

c        INPUT:  NMOM, IPOLZN, MOMDIM     MIEV0 arguments
c                CALCMO                   Flags calculated from IPOLZN
c                A(1-2), B(1-2)           Mie series coefficients

c        OUTPUT: PMOM                     Legendre moments


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      INTEGER   IPOLZN, MOMDIM, NMOM
c     ..
c     .. Array Arguments ..

      LOGICAL   CALCMO( * )
      REAL*8    PMOM( 0:MOMDIM, * )
      COMPLEX*16   A( * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER   L, NUMMOM
      REAL*8    A2SQ, B2SQ, PM1, PM2
      COMPLEX*16  A2C, B2C, CA, CAC, CAT, CB, CBC, CBT, CG, CH, CTMP
c     ..
c     .. Statement Functions ..

      REAL*8    SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = DREAL( CTMP )**2 + DIMAG( CTMP )**2
c     ..


      CA   = 3.D0*A( 1 ) - 5.D0*B( 2 )
      CAT  = 3.D0*B( 1 ) - 5.D0*A( 2 )
      CAC  = DCONJG( CA )
      A2SQ = SQ( A( 2 ) )
      B2SQ = SQ( B( 2 ) )
      A2C  = DCONJG( A( 2 ) )
      B2C  = DCONJG( B( 2 ) )


      IF( IPOLZN.LT.0 ) THEN

c                                   ** Loop over Sekera moments
         NUMMOM = MIN( NMOM, 2 )

         DO 10 L = 0, NUMMOM

            IF( CALCMO( 1 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 1 ) = 0.25D0 * ( SQ( CAT )
     &                                   + (100.D0/3.D0)* B2SQ )

               IF( L.EQ.1 ) PMOM( L, 1 ) = (5.D0/3.D0)*DREAL( CAT*B2C )

               IF( L.EQ.2 ) PMOM( L, 1 ) = (10.D0/3.D0)*B2SQ

            END IF


            IF( CALCMO( 2 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 2 ) = 0.25D0 * ( SQ( CA )
     &                                   + (100.D0/3.D0) * A2SQ )

               IF( L.EQ.1 ) PMOM( L, 2 ) = (5.D0/3.D0)*DREAL( CA*A2C )

               IF( L.EQ.2 ) PMOM( L, 2 ) = (10.D0/3.D0)*A2SQ

            END IF


            IF( CALCMO( 3 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 3 ) = 0.25D0 * DREAL( CAT * CAC
     &                                   + (100.D0/3.D0) * B(2) * A2C )

               IF( L.EQ.1 ) PMOM( L, 3 ) = 5.D0/6.D0*
     &                                     DREAL( B(2)*CAC + CAT*A2C )

               IF( L.EQ.2 ) PMOM( L, 3 ) = 10.D0/3.D0* DREAL( B(2)*A2C )

            END IF


            IF( CALCMO( 4 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 4 ) = -0.25D0 * DIMAG( CAT * CAC
     &                                   + (100.D0/3.D0)* B(2) * A2C )

               IF( L.EQ.1 ) PMOM( L, 4 ) = -5.D0/ 6.D0*
     &                                     DIMAG( B(2)*CAC + CAT*A2C )

               IF( L.EQ.2 ) PMOM( L, 4 ) = -10.D0/3.D0* DIMAG(B(2)*A2C)

            END IF


   10    CONTINUE


      ELSE

         CB  = 3.D0*B( 1 ) + 5.D0*A( 2 )
         CBT = 3.D0*A( 1 ) + 5.D0*B( 2 )
         CBC = DCONJG( CB )
         CG  = ( CBC*CBT + 10.D0*( CAC*A( 2 ) + B2C*CAT ) ) / 3.D0
         CH  = 2.D0*( CBC*A( 2 ) + B2C*CBT )

c                               ** Loop over Mueller moments
         NUMMOM = MIN( NMOM, 4 )

         DO 20 L = 0, NUMMOM


            IF( IPOLZN.EQ.0 .OR. CALCMO( 1 ) ) THEN

               IF( L.EQ.0 ) PM1 = 0.25D0*SQ( CA ) + SQ( CB ) / 12.D0
     &                          + (5.D0/3.D0)*DREAL( CA*B2C ) +5.D0*B2SQ

               IF( L.EQ.1 ) PM1 = DREAL( CB * ( CAC / 6.D0+ B2C ) )

               IF( L.EQ.2 ) PM1 = SQ( CB ) / 30.D0 + (20.D0/7.D0)*B2SQ
     &                          + (2.D0/3.D0)*DREAL( CA*B2C )

               IF( L.EQ.3 ) PM1 = (2.D0/7.D0) * DREAL( CB*B2C )

               IF( L.EQ.4 ) PM1 = (40.D0/63.D0) * B2SQ

               IF( CALCMO( 1 ) ) PMOM( L, 1 ) = PM1

            END IF


            IF( IPOLZN.EQ.0 .OR. CALCMO( 2 ) ) THEN

               IF( L.EQ.0 ) PM2 = 0.25D0*SQ( CAT ) + SQ( CBT ) / 12.D0
     &                           + ( 5.D0/ 3.D0) * DREAL( CAT*A2C )
     &                           + 5.D0*A2SQ

               IF( L.EQ.1 ) PM2 = DREAL( CBT *
     &                                 ( DCONJG( CAT ) / 6.D0+ A2C ) )

               IF( L.EQ.2 ) PM2 = SQ( CBT ) / 30.D0
     &                          + ( 20.D0/7.D0) * A2SQ
     &                          + ( 2.D0/3.D0) * DREAL( CAT*A2C )

               IF( L.EQ.3 ) PM2 = (2.D0/7.D0) * DREAL( CBT*A2C )

               IF( L.EQ.4 ) PM2 = (40.D0/63.D0) * A2SQ

               IF( CALCMO( 2 ) ) PMOM( L, 2 ) = PM2

            END IF


            IF( IPOLZN.EQ.0 ) THEN

               PMOM( L, 1 ) = 0.5D0*( PM1 + PM2 )
               GO TO  20

            END IF


            IF( CALCMO( 3 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 3 ) = 0.25D0 * DREAL( CAC*CAT + CG
     &                                   + 20.D0* B2C * A(2) )

               IF( L.EQ.1 ) PMOM( L, 3 ) = DREAL( CAC*CBT + CBC*CAT
     &                                   + 3.D0*CH ) / 12.D0

               IF( L.EQ.2 ) PMOM( L, 3 ) = 0.1D0 * DREAL( CG
     &                                   + (200.D0/7.D0) * B2C * A(2) )

               IF( L.EQ.3 ) PMOM( L, 3 ) = DREAL( CH ) / 14.D0

               IF( L.EQ.4 ) PMOM( L, 3 ) = 40.D0/63.D0* DREAL( B2C*A(2))

            END IF


            IF( CALCMO( 4 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 4 ) = 0.25D0 * DIMAG( CAC*CAT + CG
     &                                   + 20.D0* B2C * A(2) )

               IF( L.EQ.1 ) PMOM( L, 4 ) = DIMAG( CAC*CBT + CBC*CAT
     &                                   + 3.D0*CH ) / 12.D0

               IF( L.EQ.2 ) PMOM( L, 4 ) = 0.1D0 * DIMAG( CG
     &                                   + (200.D0/7.D0) * B2C * A(2) )

               IF( L.EQ.3 ) PMOM( L, 4 ) = DIMAG( CH ) / 14.D0

               IF( L.EQ.4 ) PMOM( L, 4 ) = 40.D0/63.D0* DIMAG( B2C*A(2))

            END IF


   20    CONTINUE

      END IF

      RETURN
      END

      SUBROUTINE BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA, CBIGA )

c        Calculate logarithmic derivatives of J-Bessel-function

c     Input :  CIOR, XX, NTRM, NOABS, YESANG  (defined in MIEV0)

c    Output :  RBIGA or CBIGA  (defined in MIEV0)

c    Routines called :  CONFRA


c    INTERNAL VARIABLES :

c       CONFRA     Value of Lentz continued fraction for CBIGA(NTRM),
c                     used to initialize downward recurrence.
c       DOWN       = True, use down-recurrence.  False, do not.
c       F1,F2,F3   Arithmetic statement functions used in determining
c                     whether to use up-  or down-recurrence
c                     ( Ref. 2, Eqs. 6-8 )
c       MRE        Real refractive index
c       MIM        Imaginary refractive index
c       REZINV     1 / ( MRE * XX ); temporary variable for recurrence
c       ZINV       1 / ( CIOR * XX ); temporary variable for recurrence


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   NOABS, YESANG
      INTEGER   NTRM
      REAL*8    XX
      COMPLEX*16   CIOR
c     ..
c     .. Array Arguments ..

      REAL*8    RBIGA( * )
      COMPLEX*16   CBIGA( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   DOWN
      INTEGER   N
      REAL*8    MIM, MRE, REZINV, RTMP
      COMPLEX*16   CTMP, ZINV
c     ..
c     .. External Functions ..

      COMPLEX*16   CONFRA
      EXTERNAL  CONFRA
c     ..
c     .. Statement Functions ..

      REAL*8    F1, F2
ccc, F3
c     ..
c     .. Statement Function definitions ..

      F1( MRE ) = -8.D0 + MRE**2*( 26.22D0 +
     &            MRE*( -0.4474D0 + MRE**3*( 2.04D-3 - 1.75D-4*MRE ) ) )

      F2( MRE ) = 3.9D0 + MRE*( -10.8D0 + 13.78D0*MRE )

ccc      F3( MRE ) = -15.04D0 + MRE*( 8.42D0 + 16.35D0*MRE )
c     ..

c                                  ** Decide whether BigA can be
c                                  ** calculated by up-recurrence
      MRE = DREAL( CIOR )
      MIM = DABS( DIMAG( CIOR ) )

      IF( MRE.LT.1.0 .OR. MRE.GT.10.0 .OR. MIM.GT.10.0 ) THEN

         DOWN = .TRUE.

      ELSE IF( YESANG ) THEN

         DOWN = .TRUE.
         IF( MIM*XX .LT. F2( MRE ) ) DOWN = .FALSE.

      ELSE

         DOWN = .TRUE.
         IF( MIM*XX .LT. F1( MRE ) ) DOWN = .FALSE.

      END IF


      ZINV   = 1.D0 / ( CIOR*XX )
      REZINV = 1.D0 / ( MRE*XX )


      IF( DOWN ) THEN
c                          ** Compute initial high-order BigA using
c                          ** Lentz method ( Ref. 1, pp. 17-20 )

         CTMP = CONFRA( NTRM, ZINV )

c                                   *** Downward recurrence for BigA
c                                   *** ( Ref. 1, Eq. 22 )
         IF( NOABS ) THEN
c                                            ** No-absorption case
            RBIGA( NTRM ) = DREAL( CTMP )

            DO 10 N = NTRM, 2, -1
               RBIGA( N - 1 ) = ( N*REZINV ) -
     &                          1.D0 / ( ( N*REZINV ) + RBIGA( N ) )
   10       CONTINUE

         ELSE
c                                            ** Absorptive case
            CBIGA( NTRM ) = CTMP

            DO 20 N = NTRM, 2, -1
               CBIGA( N-1 ) = (N*ZINV) - 1.D0/ ( (N*ZINV) + CBIGA( N ) )
   20       CONTINUE

         END IF


      ELSE
c                              *** Upward recurrence for BigA
c                              *** ( Ref. 1, Eqs. 20-21 )
         IF( NOABS ) THEN
c                                            ** No-absorption case
            RTMP = SIN( MRE*XX )
            RBIGA( 1 ) = - REZINV + RTMP /
     &                   ( RTMP*REZINV - COS( MRE*XX ) )

            DO 30 N = 2, NTRM
               RBIGA( N ) = -( N*REZINV ) +
     &                      1.D0 / ( ( N*REZINV ) - RBIGA( N - 1 ) )
   30       CONTINUE

         ELSE
c                                                ** Absorptive case
            CTMP = CDEXP( - (0.D0,2.D0)*CIOR*XX )
            CBIGA( 1 ) = - ZINV + (1.D0-CTMP) /
     $                 ( ZINV * (1.D0-CTMP) - (0.D0,1.D0)*(1.D0+CTMP) )

            DO 40 N = 2, NTRM
              CBIGA( N ) = - (N*ZINV) + 1.D0 / ((N*ZINV) - CBIGA( N-1 ))
   40       CONTINUE

         END IF

      END IF

      RETURN
      END

      COMPLEX*16 FUNCTION CONFRA( N, ZINV )

c         Compute Bessel function ratio A-sub-N from its
c         continued fraction using Lentz method ( Ref. 1, pp. 17-20 )

c         ZINV = Reciprocal of argument of A


c    I N T E R N A L    V A R I A B L E S
c    ------------------------------------

c    CAK      Term in continued fraction expansion of A (Ref. 1, Eq. 25)
c    CAPT     Factor used in Lentz iteration for A (Ref. 1, Eq. 27)
c    CDENOM   Denominator in capT  ( Ref. 1, Eq. 28B )
c    CNUMER   Numerator   in capT  ( Ref. 1, Eq. 28A )
c    CDTD     Product of two successive denominators of capT factors
c                 ( Ref. 1, Eq. 34C )
c    CNTN     Product of two successive numerators of capT factors
c                 ( Ref. 1, Eq. 34B )
c    EPS1     Ill-conditioning criterion
c    EPS2     Convergence criterion
c    KK       Subscript k of cAk  ( Ref. 1, Eq. 25B )
c    KOUNT    Iteration counter ( used only to prevent runaway )
c    MAXIT    Max. allowed no. of iterations
c    MM       + 1  and - 1, alternately


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      INTEGER   N
      COMPLEX*16   ZINV
c     ..
c     .. Local Scalars ..

      INTEGER   KK, KOUNT, MAXIT, MM
      REAL*8    EPS1, EPS2
      COMPLEX*16   CAK, CAPT, CDENOM, CDTD, CNTN, CNUMER
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
      DATA      EPS1 / 1.D-2 / , EPS2 / 1.D-8 /
      DATA      MAXIT / 10000 /


c                                      *** Ref. 1, Eqs. 25a, 27
      CONFRA = ( N + 1 ) * ZINV
      MM     = - 1
      KK     = 2*N + 3
      CAK    = ( MM*KK ) * ZINV
      CDENOM = CAK
      CNUMER = CDENOM + 1.D0 / CONFRA
      KOUNT  = 1

   10 CONTINUE
      KOUNT = KOUNT + 1

      IF( KOUNT.GT.MAXIT )
     &    CALL ERRMSG('ConFra--Iteration failed to converge',.TRUE.)

c                                         *** Ref. 2, Eq. 25b
      MM  = - MM
      KK  = KK + 2
      CAK = ( MM*KK ) * ZINV
c                                         *** Ref. 2, Eq. 32
      IF( CDABS( CNUMER / CAK ).LE.EPS1 .OR.
     &    CDABS( CDENOM / CAK ).LE.EPS1 ) THEN

c                                  ** Ill-conditioned case -- stride
c                                  ** two terms instead of one

c                                         *** Ref. 2, Eqs. 34
         CNTN   = CAK * CNUMER + 1.D0
         CDTD   = CAK * CDENOM + 1.D0
         CONFRA = ( CNTN / CDTD ) * CONFRA
c                                             *** Ref. 2, Eq. 25b
         MM  = - MM
         KK  = KK + 2
         CAK = ( MM*KK ) * ZINV
c                                         *** Ref. 2, Eqs. 35
         CNUMER = CAK + CNUMER / CNTN
         CDENOM = CAK + CDENOM / CDTD
         KOUNT  = KOUNT + 1
         GO TO  10

      ELSE
c                                ** Well-conditioned case

c                                        *** Ref. 2, Eqs. 26, 27
         CAPT   = CNUMER / CDENOM
         CONFRA = CAPT * CONFRA
c                                    ** Check for convergence
c                                    ** ( Ref. 2, Eq. 31 )

         IF (      DABS( DREAL (CAPT) - 1.D0 ).GE.EPS2
     $        .OR. DABS( DIMAG (CAPT) )       .GE.EPS2 )  THEN

c                                        *** Ref. 2, Eqs. 30A-B
            CNUMER = CAK + 1.D0 / CNUMER
            CDENOM = CAK + 1.D0 / CDENOM
            GO TO  10

         END IF

      END IF


      RETURN

      END

      SUBROUTINE MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                   QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM,
     &                   SFORW, SBACK, TFORW, TBACK, S1, S2 )

c         Print scattering quantities of a single particle


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   PERFCT
      INTEGER   IPOLZN, MOMDIM, NMOM, NUMANG
      REAL*8    GQSC, QEXT, QSCA, XX
      COMPLEX*16   CREFIN, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      LOGICAL   CALCMO( * ), PRNT( * )
      REAL*8    PMOM( 0:MOMDIM, * ), XMU( * )
      COMPLEX*16   S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      CHARACTER FMAT*22
      INTEGER   I, J, M
      REAL*8    FNORM, I1, I2
c     ..


      IF( PERFCT ) WRITE( *, '(''1'',10X,A,1P,E11.4)' )
     &    'Perfectly Conducting Case, size parameter =', XX

      IF( .NOT.PERFCT ) WRITE( *, '(''1'',10X,3(A,1P,E11.4))' )
     &    'Refractive Index:  Real ', DREAL( CREFIN ), '  Imag ',
     &    DIMAG( CREFIN ), ',   Size Parameter =', XX


      IF( PRNT( 1 ) .AND. NUMANG.GT.0 ) THEN

         WRITE( *, '(/,A)' )
     &      '    cos(angle)  ------- S1 ---------  ------- S2 ---------'
     &      // '  --- S1*conjg(S2) ---   i1=S1**2   i2=S2**2  (i1+i2)/2'
     &      // '  DEG POLZN'

         DO 10 I = 1, NUMANG
            I1 = DREAL( S1( I ) )**2 + DIMAG( S1( I ) )**2
            I2 = DREAL( S2( I ) )**2 + DIMAG( S2( I ) )**2
            WRITE( *, '( I4, F10.6, 1P,10E11.3 )'   )
     $              I, XMU(I), S1(I), S2(I), S1(I)*DCONJG(S2(I)),
     $              I1, I2, 0.5D0*(I1+I2), (I2-I1)/(I2+I1)
   10    CONTINUE

      END IF


      IF( PRNT( 2 ) ) THEN

         WRITE ( *, '(/,A,9X,A,17X,A,17X,A,/,(0P,F7.2, 1P,6E12.3) )' )
     $           '  Angle', 'S-sub-1', 'T-sub-1', 'T-sub-2',
     $               0.D0,     SFORW,    TFORW(1),  TFORW(2),
     $              180.D0,    SBACK,    TBACK(1),  TBACK(2)
         WRITE ( *, '(/,4(A,1P,E11.4))' )
     $           ' Efficiency Factors,  extinction:', QEXT,
     $                              '   scattering:', QSCA,
     $                              '   absorption:', QEXT-QSCA,
     $                           '   rad. pressure:', QEXT-GQSC

         IF( NMOM.GT.0 ) THEN

            WRITE( *, '(/,A)' ) ' Normalized moments of :'

            IF( IPOLZN.EQ.0 ) WRITE( *, '(''+'',27X,A)' )
     &          'Phase Fcn'

            IF( IPOLZN.GT.0 ) WRITE( *, '(''+'',33X,A)' )
     &          'M1           M2          S21          D21'

            IF( IPOLZN.LT.0 ) WRITE( *, '(''+'',33X,A)' )
     &          'R1           R2           R3           R4'

            FNORM = 4.D0/ ( XX**2 * QSCA )

            DO 30  M = 0, NMOM

               WRITE( *, '(A,I4)' ) '      Moment no.', M

               DO 20 J = 1, 4

                  IF( CALCMO( J ) ) THEN
                     WRITE( FMAT, '(A,I2,A)' )
     &                      '( ''+'', T', 24+(J-1)*13, ', 1P,E13.4 )'
                     WRITE( *, FMAT ) FNORM * PMOM( M, J )
                  END IF

   20          CONTINUE
   30       CONTINUE

         END IF

      END IF


      RETURN

      END

      SUBROUTINE SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW,
     &                   SBACK, S1, S2, TFORW, TBACK, A, B )

c       Small-particle limit of Mie quantities in totally reflecting
c       limit ( Mie series truncated after 2 terms )

c        A,B       First two Mie coefficients, with numerator and
c                  denominator expanded in powers of XX ( a factor
c                  of XX**3 is missing but is restored before return
c                  to calling program )  ( Ref. 2, p. 1508 )

      IMPLICIT  NONE

c     .. Parameters ..

      REAL*8      TWOTHR, FIVTHR, FIVNIN
      PARAMETER ( TWOTHR = 2.D0/3.D0, FIVTHR = 5.D0/3.D0,
     &            FIVNIN = 5.D0/9.D0 )
c     ..
c     .. Scalar Arguments ..

      INTEGER   NUMANG
      REAL*8    GQSC, QEXT, QSCA, XX
      COMPLEX*16   SBACK, SFORW
c     ..
c     .. Array Arguments ..

      REAL*8    XMU( * )
      COMPLEX*16   A( * ), B( * ), S1( * ), S2( * ),
     &                 TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      INTEGER   J
      REAL*8    RTMP
      COMPLEX*16   CTMP
c     ..
c     .. Statement Functions ..

      REAL*8    SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = DREAL( CTMP )**2 + DIMAG( CTMP )**2
c     ..


      A( 1 ) = DCMPLX( 0.D0, TWOTHR*( 1.D0 - 0.2D0*XX**2 ) ) /
     &         DCMPLX( 1.D0 - 0.5D0*XX**2, TWOTHR*XX**3 )

      B( 1 ) = DCMPLX( 0.D0, - ( 1.D0 - 0.1D0*XX**2 ) / 3.D0) /
     &         DCMPLX( 1.D0 + 0.5D0*XX**2, - XX**3 / 3.D0)

      A( 2 ) = DCMPLX( 0.D0,   XX**2 / 30.D0)
      B( 2 ) = DCMPLX( 0.D0, - XX**2 / 45.D0)

      QSCA = 6.D0* XX**4 *( SQ( A(1) ) + SQ( B(1) ) +
     &             FIVTHR*( SQ( A(2) ) + SQ( B(2) ) ) )
      QEXT = QSCA
      GQSC = 6.D0* XX**4 *DREAL( A(1)*DCONJG( A(2) + B(1) ) +
     &         ( B(1) + FIVNIN*A(2) )*DCONJG( B(2) ) )

      RTMP   = 1.5D0 * XX**3
      SFORW  = RTMP*( A(1) + B(1) + FIVTHR*( A(2) + B(2) ) )
      SBACK  = RTMP*( A(1) - B(1) - FIVTHR*( A(2) - B(2) ) )
      TFORW( 1 ) = RTMP*( B(1) + FIVTHR*( 2.D0*B(2) - A(2) ) )
      TFORW( 2 ) = RTMP*( A(1) + FIVTHR*( 2.D0*A(2) - B(2) ) )
      TBACK( 1 ) = RTMP*( B(1) - FIVTHR*( 2.D0*B(2) + A(2) ) )
      TBACK( 2 ) = RTMP*( A(1) - FIVTHR*( 2.D0*A(2) + B(2) ) )


      DO 10 J = 1, NUMANG
         S1( J ) = RTMP * ( A(1) + B(1) * XMU(J) + FIVTHR *
     $           ( A(2) * XMU(J) + B(2) * ( 2.D0*XMU(J)**2 - 1.D0)) )
         S2( J ) = RTMP * ( B(1) + A(1) * XMU(J) + FIVTHR *
     $           ( B(2) * XMU(J) + A(2) * ( 2.D0*XMU(J)**2 - 1.D0)) )
   10 CONTINUE

c                                     ** Recover actual Mie coefficients
      A( 1 ) = XX**3 * A(1)
      A( 2 ) = XX**3 * A(2)
      B( 1 ) = XX**3 * B(1)
      B( 2 ) = XX**3 * B(2)

      RETURN
      END

      SUBROUTINE SMALL2( XX, CIOR, CALCQE, NUMANG, XMU, QEXT, QSCA,
     &                   GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK,
     &                   A, B )

c       Small-particle limit of Mie quantities for general refractive
c       index ( Mie series truncated after 2 terms )

c        A,B       First two Mie coefficients, with numerator and
c                  denominator expanded in powers of XX ( a factor
c                  of XX**3 is missing but is restored before return
c                  to calling program )  ( Ref. 2, p. 1508 )

c        CIORSQ    Square of refractive index


      IMPLICIT  NONE

c     .. Parameters ..

      REAL*8      TWOTHR, FIVTHR
      PARAMETER ( TWOTHR = 2.D0/3.D0, FIVTHR = 5.D0/3.D0)
c     ..
c     .. Scalar Arguments ..

      LOGICAL   CALCQE
      INTEGER   NUMANG
      REAL*8    GQSC, QEXT, QSCA, XX
      COMPLEX*16   CIOR, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      REAL*8    XMU( * )
      COMPLEX*16   A( * ), B( * ), S1( * ), S2( * ),
     &                 TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      INTEGER   J
      REAL*8    RTMP
      COMPLEX*16   CIORSQ, CTMP
c     ..
c     .. Statement Functions ..

      REAL*8    SQ
c     ..
c     .. Statement Function definitions ..

      SQ( CTMP ) = DREAL( CTMP )**2 + DIMAG( CTMP )**2
c     ..


      CIORSQ = CIOR**2
      CTMP   = DCMPLX( 0.D0, TWOTHR )*( CIORSQ - 1.D0)
      A( 1 ) = CTMP*( 1.D0- 0.1D0*XX**2 +
     &         ( CIORSQ / 350.D0 + 1.D0/280.D0)*XX**4 ) /
     &         ( CIORSQ + 2.D0+ ( 1.D0- 0.7D0*CIORSQ )*XX**2 -
     &         ( CIORSQ**2 / 175.D0- 0.275D0*CIORSQ + 0.25D0 )*XX**4 +
     &         XX**3 * CTMP * ( 1.D0- 0.1D0*XX**2 ) )

      B( 1 ) = ( XX**2 / 30.D0 )*CTMP*( 1.D0+
     &         ( CIORSQ / 35.D0 - 1.D0/ 14.D0)*XX**2 ) /
     &         ( 1.D0- ( CIORSQ / 15.D0 - 1.D0/6.D0)*XX**2 )

      A( 2 ) = ( 0.1D0*XX**2 )*CTMP*( 1.D0- XX**2 / 14.D0 ) /
     &         ( 2.D0*CIORSQ + 3.D0- ( CIORSQ / 7.D0- 0.5D0 ) * XX**2 )

      QSCA = 6.D0* XX**4 * ( SQ( A(1) ) + SQ( B(1) ) +
     &                     FIVTHR * SQ( A(2) ) )
      GQSC = 6.D0*XX**4*DREAL( A(1)*DCONJG( A(2) + B(1) ) )
      QEXT = QSCA
      IF( CALCQE ) QEXT = 6.D0*XX * DREAL( A(1) + B(1) + FIVTHR*A(2) )

      RTMP   = 1.5D0 * XX**3
      SFORW  = RTMP*( A(1) + B(1) + FIVTHR*A(2) )
      SBACK  = RTMP*( A(1) - B(1) - FIVTHR*A(2) )
      TFORW( 1 ) = RTMP*( B(1) - FIVTHR*A(2) )
      TFORW( 2 ) = RTMP*( A(1) + 2.D0*FIVTHR*A(2) )
      TBACK( 1 ) = TFORW(1)
      TBACK( 2 ) = RTMP*( A(1) - 2.D0*FIVTHR*A(2) )


      DO 10 J = 1, NUMANG
         S1( J ) = RTMP*( A(1) + ( B(1) + FIVTHR*A(2) )*XMU( J ) )
         S2( J ) = RTMP*( B(1) + A(1)*XMU( J ) +
     &                    FIVTHR*A(2)*( 2.D0*XMU( J )**2 - 1.D0) )
   10 CONTINUE

c                                     ** Recover actual Mie coefficients
      A( 1 ) = XX**3 * A(1)
      A( 2 ) = XX**3 * A(2)
      B( 1 ) = XX**3 * B(1)
      B( 2 ) = ( 0.D0, 0.D0)

      RETURN
      END

      SUBROUTINE TESTMI( COMPAR, XX, CREFIN, MIMCUT, PERFCT, ANYANG,
     &                   NMOM, IPOLZN, NUMANG, XMU, QEXT, QSCA, GQSC,
     &                   SFORW, SBACK, S1, S2, TFORW, TBACK, PMOM,
     &                   MOMDIM )

c         Set up to run test case when  COMPAR = False;  when  = True,
c         compare Mie code test case results with correct answers
c         and abort if even one result is inaccurate.

c         The test case is :  Mie size parameter = 10
c                             refractive index   = 1.5 - 0.1 i
c                             scattering angle = 140 degrees
c                             1 Sekera moment

c         Results for this case may be found among the test cases
c         at the end of reference (1).

c         *** NOTE *** When running on some computers, esp. in single
c         precision, the Accur criterion below may have to be relaxed.
c         However, if Accur must be set larger than 10**-3 for some
c         size parameters, your computer is probably not accurate
c         enough to do Mie computations for those size parameters.

c     Routines called :  ERRMSG, MIPRNT, TSTBAD


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL   ANYANG, COMPAR, PERFCT
      INTEGER   IPOLZN, MOMDIM, NMOM, NUMANG
      REAL*8    GQSC, MIMCUT, QEXT, QSCA, XX
      COMPLEX*16   CREFIN, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      REAL*8    PMOM( 0:MOMDIM, * ), XMU( * )
      COMPLEX*16   S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   ANYSAV, OKEY, PERSAV
      INTEGER   IPOSAV, M, N, NMOSAV, NUMSAV
      REAL*8    ACCUR, CALC, EXACT, MIMSAV, TESTGQ, TESTQE, TESTQS,
     &          XMUSAV, XXSAV
      COMPLEX*16   CRESAV, TESTS1, TESTS2, TESTSB, TESTSF
c     ..
c     .. Local Arrays ..

      LOGICAL   CALCMO( 4 ), PRNT( 2 )
      REAL*8    TESTPM( 0:1 )
      COMPLEX*16   TESTTB( 2 ), TESTTF( 2 )
c     ..
c     .. External Functions ..

      LOGICAL   TSTBAD
      EXTERNAL  TSTBAD
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, MIPRNT
c     ..
c     .. Statement Functions ..

      LOGICAL   WRONG
c     ..
      SAVE      XXSAV, CRESAV, MIMSAV, PERSAV, ANYSAV, NMOSAV, IPOSAV,
     &          NUMSAV, XMUSAV

      DATA      TESTQE / 2.459791D0 /,
     &          TESTQS / 1.235144D0 /,
     &          TESTGQ / 1.139235D0 /,
     &          TESTSF / ( 61.4947600D0, -3.1779940D0 ) /,
     &          TESTSB / ( 1.49343400D0,  0.2963657D0 ) /,
     &          TESTS1 / (-0.15483800D0, -1.1289720D0 ) /,
     &          TESTS2 / ( 0.05669755D0,  0.5425681D0 ) /,
     &          TESTTF / ( 12.9523800D0, -136.64360D0 ),
     &                   ( 48.5423800D0,  133.46560D0 ) /,
     &          TESTTB / ( 41.8841400D0, -15.578330D0 ),
     &                   ( 43.3775800D0, -15.281960D0 ) /,
     &          TESTPM / 227.1975D0, 183.6898D0 /

      DATA      ACCUR / 1.D-4 /
c     ..
c     .. Statement Function definitions ..

      WRONG( CALC, EXACT ) = DABS( ( CALC - EXACT ) / EXACT ).GT.ACCUR
c     ..


      IF( .NOT.COMPAR ) THEN
c                                   ** Save certain user input values
         XXSAV  = XX
         CRESAV = CREFIN
         MIMSAV = MIMCUT
         PERSAV = PERFCT
         ANYSAV = ANYANG
         NMOSAV = NMOM
         IPOSAV = IPOLZN
         NUMSAV = NUMANG
         XMUSAV = XMU( 1 )
c                                   ** Reset input values for test case
         XX     = 10.D0
         CREFIN = ( 1.5D0, -0.1D0 )
         MIMCUT = 0.D0
         PERFCT = .FALSE.
         ANYANG = .TRUE.
         NMOM   = 1
         IPOLZN = -1
         NUMANG = 1
         XMU( 1 ) = -0.7660444D0

      ELSE
c                                    ** Compare test case results with
c                                    ** correct answers and abort if bad
         OKEY = .TRUE.

         IF( WRONG( QEXT,TESTQE ) )
     &       OKEY = TSTBAD( 'QEXT', DABS( ( QEXT - TESTQE ) / TESTQE ) )

         IF( WRONG( QSCA,TESTQS ) )
     &       OKEY = TSTBAD( 'QSCA', DABS( ( QSCA - TESTQS ) / TESTQS ) )

         IF( WRONG( GQSC,TESTGQ ) )
     &       OKEY = TSTBAD( 'GQSC', DABS( ( GQSC - TESTGQ ) / TESTGQ ) )

         IF( WRONG( DREAL( SFORW ),DREAL( TESTSF ) ) .OR.
     &       WRONG( DIMAG( SFORW ),DIMAG( TESTSF ) ) )
     &       OKEY = TSTBAD( 'SFORW', CDABS( ( SFORW - TESTSF ) / TESTSF ))

         IF( WRONG( DREAL( SBACK ),DREAL( TESTSB ) ) .OR.
     &       WRONG( DIMAG( SBACK ),DIMAG( TESTSB ) ) )
     &       OKEY = TSTBAD( 'SBACK', CDABS( ( SBACK - TESTSB ) / TESTSB ))

         IF( WRONG( DREAL( S1(1) ),DREAL( TESTS1 ) ) .OR.
     &       WRONG( DIMAG( S1(1) ),DIMAG( TESTS1 ) ) )
     &       OKEY = TSTBAD( 'S1', CDABS( ( S1(1) - TESTS1 ) / TESTS1 ))

         IF( WRONG( DREAL( S2(1) ),DREAL( TESTS2 ) ) .OR.
     &       WRONG( DIMAG( S2(1) ),DIMAG( TESTS2 ) ) )
     &       OKEY = TSTBAD( 'S2', CDABS( ( S2(1) - TESTS2 ) / TESTS2 ))


         DO 10  N = 1, 2

            IF( WRONG( DREAL( TFORW(N) ),DREAL( TESTTF(N) ) ) .OR.
     &          WRONG( DIMAG( TFORW(N) ),
     &          DIMAG( TESTTF(N) ) ) ) OKEY = TSTBAD( 'TFORW',
     &          CDABS( ( TFORW(N) - TESTTF(N) ) / TESTTF(N) ) )
     
            IF( WRONG( DREAL( TBACK(N) ),DREAL( TESTTB(N) ) ) .OR.
     &          WRONG( DIMAG( TBACK(N) ),
     &          DIMAG( TESTTB(N) ) ) ) OKEY = TSTBAD( 'TBACK',
     &          CDABS( ( TBACK(N) - TESTTB(N) ) / TESTTB(N) ) )

   10    CONTINUE


         DO 20 M = 0, 1

            IF ( WRONG( PMOM(M,1), TESTPM(M) ) )
     $           OKEY =  TSTBAD( 'PMOM', DABS( (PMOM(M,1)-TESTPM(M)) /
     $                                      TESTPM(M) ) )

   20    CONTINUE


         IF( .NOT.OKEY ) THEN

            PRNT( 1 ) = .TRUE.
            PRNT( 2 ) = .TRUE.
            CALCMO( 1 ) = .TRUE.
            CALCMO( 2 ) = .FALSE.
            CALCMO( 3 ) = .FALSE.
            CALCMO( 4 ) = .FALSE.

            CALL MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                   QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM,
     &                   SFORW, SBACK, TFORW, TBACK, S1, S2 )
            CALL ERRMSG( 'MIEV0 -- Self-test failed',.TRUE.)

         END IF
c                                       ** Restore user input values
         XX     = XXSAV
         CREFIN = CRESAV
         MIMCUT = MIMSAV
         PERFCT = PERSAV
         ANYANG = ANYSAV
         NMOM   = NMOSAV
         IPOLZN = IPOSAV
         NUMANG = NUMSAV
         XMU( 1 ) = XMUSAV

      END IF

      RETURN
      END

