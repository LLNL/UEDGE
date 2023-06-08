c!include "../mppl.h"
c!include "../sptodp.h"
      SUBROUTINE VODPK (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, PSOL, MF, RPAR, IPAR,
     1  scalpk, efacnin)
      EXTERNAL F, JAC, PSOL
      real Y, T, TOUT, RTOL, ATOL, RWORK, RPAR, scalpk, efacnin
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF, IPAR
      DIMENSION Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     1          RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C This is the 30 November 1993 version of
C VODPK.. variable-coefficient ordinary differential equation solver
C         with the preconditioned Krylov method GMRES for the solution
C         of linear systems.
C
C Modified for use in LEDGE and UEDGE by T. Rognlien et al.
C
C Modified by A. C. Hindmarsh: added failure flag in F routine.
C
C This version is in single precision.
C
C VODPK solves the initial value problem for stiff or nonstiff
C systems of first order ode-s,
C     dy/dt = f(t,y) ,  or, in component form,
C     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
C VODPK is a package based on the VODE and LSODPK packages, and on
C the October 23, 1978 version of the ODEPACK user interface
C standard, with minor modifications.
C-----------------------------------------------------------------------
C Note to LLNL Cray users..
C For maximum efficiency, appropriate compiler optimization directives
C have been inserted for Cray compilers (but not CIVIC).
C-----------------------------------------------------------------------
C References..
C 1. P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, "VODE, A Variable-
C    Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10  (1989),
C    pp., 1038-1051.  Also LLNL report UCRL-98412, June 1988.
C 2. P. N. Brown and A. C. Hindmarsh, "Reduced Storage Matrix Methods
C    in Stiff ODE Systems," J. Appl. Math. & Comp., 31 (1989), pp.40-91.
C    Also LLNL report UCRL-95088, Rev. 1, June 1987.
C 3. G. D. Byrne, "Pragmatic Experiments with Krylov Methods in the
C    Stiff ODE Setting," Numerical ODEs, J. Cash and I. Duff, eds.,
C    Oxford Univ. Press, London, to appear.
C-----------------------------------------------------------------------
C Authors..
C               George D. Byrne
C               Technical Computing Division
C               Exxon Research and Engineering Co.
C               P. O. Box 998
C               Annandale, NJ 08801-0998
C and
C               Alan C. Hindmarsh and Peter N. Brown
C               Computing and Mathematics Research Division, L-316
C               Lawrence Livermore National Laboratory
C               Livermore, CA 94550
C-----------------------------------------------------------------------
C Introduction.
C
C This is a modification of the VODE package which incorporates
C the preconditioned Krylov subspace iterative method SPIGMR for the
C linear algebraic systems that arise in the case of stiff systems.
C SPIGMR denotes a scaled preconditioned incomplete version of the
C GMRES (Generalized Minimum Residual) method.
C
C The linear systems that are solved have the form
C   A * x  = b ,  where  A = I - hrl1 * (df/dy) .
C here hrl1 is a scalar, I is the identity matrix, and df/dy is the
C Jacobian matrix of partial derivatives of f with respect to y
C (an NEQ by NEQ matrix).
C
C The particular Krylov method is chosen by setting the second digit,
C MITER, in the method flag MF.
C Currently, the values of MITER have the following meanings..
C
C          1 means SPIGMR, a scaled, preconditioned, incomplete version
C            of GMRES, a generalized minimum residual method.
C            This is the best choice in general.
C
C          9 means that only a user-supplied matrix P (approximating A)
C            will be used, with no Krylov iteration done internally to
C            VODPK.  This option allows the user to provide the complete
C            linear system solution algorithm, if desired.
C
C The user can apply preconditioning to the linear system A*x = b,
C by means of arbitrary matrices (the preconditioners).
C
C     In the case of SPIGMR, one can apply left and right
C preconditioners P1 and P2, and the basic iterative method is then
C applied to the matrix (P1-inverse)*A*(P2-inverse) instead of to the
C matrix A.  The product P1*P2 should be an approximation to A
C such that linear systems with P1 or P2 are easier to solve than with
C A alone.  Preconditioning from the left only or right only means using
C P2 = I  or  P1 = I, respectively.
C
C     If the Jacobian  J = df/dy  splits in a natural way into a sum
C J = J1 + J2, then one possible choice of preconditioners is
C            P1 = I - hrl1 * J1  and  P2 = I - hrl1 * J2
C provided each of these is easy to solve (or to approximately solve).
C
C NOTE:  To achieve an efficient solution, the preconditioned Krylov
C methods in VODPK generally require a thoughtful choice of
C preconditioners.  If the ODE system produces linear systems that are
C not amenable to solution by such iterative methods, the cost can be
C higher than with a solver that uses sparse direct methods.  However,
C for many systems, careful use of VODPK can be highly effective.
C
C See Ref. 2 for more details on the methods and applications.
C-----------------------------------------------------------------------
C Summary of usage.
C
C Communication between the user and the VODPK package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See full description (below)
C for details, including optional communication, nonstandard options,
C and instructions for special situations.  See also the example
C program embedded in the comments below.
C
C A. First provide a subroutine of the form..
C     SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR, IFAIL)
C     DIMENSION Y(NEQ), YDOT(NEQ), RPAR(*), IPAR(*)
C which supplies the vector function f by loading YDOT(i) with f(i).
C It should set IFAIL .ne. 0 if unable to evaluate this function.
C
C B. Next determine (or guess) whether or not the problem is stiff.
C Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
C whose real part is negative and large in magnitude, compared to the
C reciprocal of the t span of interest.  If the problem is nonstiff,
C use method flag MF = 10.  If it is stiff, MF should be 21.
C
C The following four parameters must also be set.
C  IWORK(1) = LWP  = length of real array WP for preconditioning.
C  IWORK(2) = LIWP = length of integer array IWP for preconditioning.
C  IWORK(3) = JPRE = preconditioner type flag..
C                  = 0 for no preconditioning (P1 = P2 = I)
C                  = 1 for left-only preconditioning (P2 = I)
C                  = 2 for right-only preconditioning (P1 = I)
C                  = 3 for two-sided preconditioning
C  IWORK(4) = JACFLG = flag for whether JAC is called.
C                    = 0 if JAC is not to be called,
C                    = 1 if JAC is to be called.
C  Use JACFLG = 1 if JAC computes any nonconstant data for use in
C  preconditioning, such as Jacobian elements.  See next paragraph.
C  The arrays WP and IWP are work arrays under the user-s control,
C  for use in the routines that perform preconditioning operations.
C
C C. If the problem is stiff, you must supply two routines that deal
C with the preconditioning of the linear systems to be solved.
C These are as follows..
C
C     SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY, V, HRL1, WP, IWP,
C                     IER, RPAR, IPAR)
C     DIMENSION Y(NEQ), YSV(NEQ), REWT(NEQ), FTY(NEQ), V(NEQ),
C               WP(*), IWP(*), RPAR(*), IPAR(*)
C
C        This routine is optional, and is to evaluate and preprocess
C     any parts of the Jacobian matrix df/dy involved in the
C     preconditioners P1 and P2.
C     The Y and FTY arrays contain the current values of y and f(t,y),
C     respectively, and YSV also contains the current value of y.
C     The array V is work space of length NEQ.
C     JAC must multiply all computed Jacobian elements by the scalar
C     -hrl1, add the identity matrix I, and do any factorization
C     operations called for, in preparation for solving linear systems
C     with a coefficient matrix of P1 or P2.  The matrix P1*P2
C     should be an approximation to  I - hrl1 * (df/dy), where hrl1 is a
C     scalar stored in HRL1.
C     JAC should return IER = 0 if successful, and IER .ne. 0 if not.
C     (If IER .ne. 0, a smaller time step will be tried.)
C
C     SUBROUTINE PSOL (NEQ, T, Y, FTY, WK, HRL1, WP, IWP, B, LR,
C                      IER, RPAR, IPAR)
C     DIMENSION Y(NEQ), FTY(NEQ), WK(NEQ), WP(*), IWP(*), B(NEQ),
C               RPAR(*), IPAR(*)
C
C        This routine must solve a linear system with b (stored in B)
C     as right-hand side and one of the preconditioning matrices, P1 or
C     P2, as coefficient matrix, and return the solution vector in B.
C     LR is a flag concerning left vs. right preconditioning, input
C     to PSOL.  PSOL is to use P1 if LR = 1, and P2 if LR = 2.
C
C        PSOL can use data generated in the JAC routine and stored in
C     WP and IWP.  WK is a work array of length NEQ.
C     The argument HRL1 is the current value of the scalar appearing
C     in the linear system.  If the old value, at the time of the last
C     JAC call, is needed, it must have been saved by JAC in WP.
C     on return, PSOL should set the error flag  IER as follows..
C        IER = 0 if PSOL was successful,
C        IER .gt. 0 if a recoverable error occurred, meaning that the
C              time step will be retried,
C        IER .lt. 0 if an unrecoverable error occurred, meaning that the
C              solver is to stop immediately.
C
C D. Write a main program which calls subroutine VODPK once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages
C by VODPK.  on the first call to VODPK, supply arguments as follows..
C F      = name of subroutine for right-hand side vector f.
C          This name must be declared EXTERNAL in calling program.
C NEQ    = number of first order ODE-s.
C Y      = array of initial values, of length NEQ.
C T      = the initial value of the independent variable.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          The estimated local error in Y(i) will be controlled so as
C          to be roughly less (in magnitude) than
C             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution.. Actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of Y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1.
C IOPT   = 0 to indicate no optional input used.
C RWORK  = real work array of length at least..
C             20 + 16*NEQ           for MF = 10,
C             61 + 17*NEQ + LWP     for MF = 21.
C LRW    = declared length of RWORK (in user-s DIMENSION statement).
C IWORK  = integer work array of length at least..
C             30            for MF = 10,
C             35 + LIWP     for MF = 21.
C LIW    = declared length of IWORK (in user-s DIMENSION statement).
C JAC,PSOL = names of subroutines for preconditioning.  These names
C            must be declared EXTERNAL in the user-s calling program.
C MF     = method flag.  Standard values are..
C          10 for nonstiff (Adams) method.
C          21 for stiff (BDF) method, with SPIGMR.
C
C RPAR, IPAR  User-specified arrays used to communicate real and integer
C             parameters (respectively) to user-supplied subroutines.
C             to user-supplied subroutines.  If RPAR is a vector, then
C             it must be dimensioned in the user's main program.  If it
C             is unused or it is a scalar, then it need not be
C             dimensioned.
C
C IPAR     User-specified array used to communicate integer parameter
C          to user-supplied subroutines.  The comments on dimensioning
C          RPAR apply to IPAR.
C
C Note that the user-s main (calling) program must declare arrays
C Y, RWORK, IWORK, and possibly ATOL, RPAR, and IPAR.
C
C E. The output from the first call (or any call) is..
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable (normally TOUT).
C ISTATE = 2  if VODPK was successful, negative otherwise.
C          -1 means excess work done on this call (perhaps wrong MF).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all input).
C          -5 means repeated convergence failures (perhaps bad JAC
C             or PSOL routine supplied or wrong choice of MF or
C             tolerances, or this solver is inappropriate).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 means an unrecoverable error occurred in JAC or PSOL.
C          -8 means f returned IFAIL .ne. 0 and VODPK was unable to
C             avoid this condition.
C
C F. To continue the integration after a successful return, simply
C reset TOUT and call VODPK again.  No other parameters need be reset.
C
C-----------------------------------------------------------------------
C Example problem.
C An ODE system is generated from the following 2-species diurnal
C kinetics advection-diffusion PDE system in 2 space dimensions:
C
C dc(i)/dt = Kh*(d/dx)**2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)
C                 + Ri(c1,c2,t)      for i = 1,2,   where
C   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
C   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
C   Kv(z) = Kv0*exp(z/5) ,
C Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
C vary diurnally.   The problem is posed on the square
C   0 .le. x .le. 20,    30 .le. z .le. 50   (all in km),
C with homogeneous Neumann boundary conditions, and for time t in
C   0 .le. t .le. 86400 sec (1 day).
C The PDE system is treated by central differences on a uniform
C 10 x 10 mesh, with simple polynomial initial profiles.
C The problem is solved with VODPK, with the BDF/GMRES method and
C the block-diagonal part of the Jacobian as a left preconditioner.
C-----------------------------------------------------------------------
C      EXTERNAL FEX, JACBD, SOLBD
C      COMMON /PCOM/ MX,MZ,MM,Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO
C      DIMENSION Y(2, 10, 10), RWORK(3861), IWORK(230)
C      DATA DKH/4.0E-6/, VEL/0.001E0/, DKV0/1.0E-8/, HALFDA/4.32E4/,
C     1  PI/3.1415926535898E0/, TWOHR/7200.0E0/, RTOL/1.0E-5/,
C     2  FLOOR/100.0E0/, LRW/3861/, LIW/230/, MF/21/, JPRE/1/, JACFLG/1/
C
C Load Common block of problem parameters.
C      MX = 10
C      MZ = 10
C      MM = MX*MZ
C      Q1 = 1.63E-16
C      Q2 = 4.66E-16
C      A3 = 22.62E0
C      A4 = 7.601E0
C      OM = PI/HALFDA
C      C3 = 3.7E16
C      DX = 20.0E0/(MX - 1.0E0)
C      DZ = 20.0E0/(MZ - 1.0E0)
C      HDCO = DKH/DX**2
C      HACO = VEL/(2.0E0*DX)
C      VDCO = (1.0E0/DZ**2)*DKV0
C Set other input arguments.
C      ATOL = RTOL*FLOOR
C      NEQ = 2*MX*MZ
C      IWORK(1) = 4*MX*MZ
C      IWORK(2) = NEQ
C      IWORK(3) = JPRE
C      IWORK(4) = JACFLG
C      T = 0.0E0
C      TOUT = TWOHR
C      ISTATE = 1
C Set initial profiles.
C      DO 20 JZ = 1,MZ
C        Z = 30.0E0 + (JZ - 1.0E0)*DZ
C        CZ = (0.1E0*(Z - 40.0E0))**2
C        CZ = 1.0E0 - CZ + 0.5E0*CZ**2
C        DO 10 JX = 1,MX
C          X = (JX - 1.0E0)*DX
C          CX = (0.1E0*(X - 10.0E0))**2
C          CX = 1.0E0 - CX + 0.5E0*CX**2
C          Y(1,JX,JZ) = 1.0E6*CX*CZ
C          Y(2,JX,JZ) = 1.0E12*CX*CZ
C 10       CONTINUE
C 20     CONTINUE
C
C Loop over output points, call VODPK, print sample solution values.
C      DO 70 IOUT = 1,12
C        CALL VODPK (FEX, NEQ, Y, T, TOUT, 1, RTOL, ATOL, 1, ISTATE, 0,
C     1            RWORK, LRW, IWORK, LIW, JACBD, SOLBD, MF, RPAR, IPAR)
C        WRITE(6,50) T,IWORK(11),IWORK(14),RWORK(11)
C 50     FORMAT(/' t =',E10.2,5X,'no. steps =',I5,
C     1                      '   order =',I3,'   stepsize =',E10.2)
C        WRITE(6,60) Y(1,1,1), Y(1,5,5), Y(1,10,10),
C     1              Y(2,1,1), Y(2,5,5), Y(2,10,10)
C 60     FORMAT('  c1 (bot.left/middle/top rt.) =',3E12.3/
C     1         '  c2 (bot.left/middle/top rt.) =',3E12.3)
C        IF (ISTATE .NE. 2) STOP
C        TOUT = TOUT + TWOHR
C 70     CONTINUE
C
C Print final statististics.
C      LENRW = IWORK(17)
C      LENIW = IWORK(18)
C      NST = IWORK(11)
C      NFE = IWORK(12)
C      NPE = IWORK(13)
C      NPS = IWORK(24)
C      NNI = IWORK(20)
C      NLI = IWORK(23)
C      AVDIM = FLOAT(NLI)/FLOAT(NNI)
C      NCFN = IWORK(21)
C      NCFL = IWORK(25)
C      WRITE (6,80) LENRW,LENIW,NST,NFE,NPE,NPS,NNI,NLI,AVDIM,NCFN,NCFL
C 80   FORMAT(//' final statistics..'/
C     1 ' RWORK size =',I5,5X,' IWORK size =',I4/
C     2 ' number of steps        =',I5,5X,'number of f evals.     =',I5/
C     3 ' number of prec. evals. =',I5,5X,'number of prec. solves =',I5/
C     4 ' number of nonl. iters. =',I5,5X,'number of lin. iters.  =',I5/
C     5 ' average Krylov subspace dimension (NLI/NNI)  =',F8.4/
C     6 ' number of conv. failures..  nonlinear =',I3,'  linear =',I3)
C      STOP
C      END
C
C      SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR, IFAIL)
C      DIMENSION Y(2,*), YDOT(2,*)
C      COMMON /PCOM/ MX,MZ,MM,Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO
C
C Set diurnal rate coefficients.
C      S = SIN(OM*T)
C      IF (S .GT. 0.0E0) THEN
C        Q3 = EXP(-A3/S)
C        Q4 = EXP(-A4/S)
C      ELSE
C        Q3 = 0.0E0
C        Q4 = 0.0E0
C      ENDIF
C Loop over all grid points.
C      DO 20 JZ = 1,MZ
C        ZDN = 30.0E0 + (JZ - 1.5E0)*DZ
C        ZUP = ZDN + DZ
C        CZDN = VDCO*EXP(0.2E0*ZDN)
C        CZUP = VDCO*EXP(0.2E0*ZUP)
C        IBLOK0 = (JZ-1)*MX
C        IDN = -MX
C        IF (JZ .EQ. 1) IDN = MX
C        IUP = MX
C        IF (JZ .EQ. MZ) IUP = -MX
C        DO 10 JX = 1,MX
C          IBLOK = IBLOK0 + JX
C          C1 = Y(1,IBLOK)
C          C2 = Y(2,IBLOK)
C Set kinetic rate terms.
C          QQ1 = Q1*C1*C3
C          QQ2 = Q2*C1*C2
C          QQ3 = Q3*C3
C          QQ4 = Q4*C2
C          RKIN1 = -QQ1 - QQ2 + 2.0E0*QQ3 + QQ4
C          RKIN2 = QQ1 - QQ2 - QQ4
C Set vertical diffusion terms.
C          C1DN = Y(1,IBLOK+IDN)
C          C2DN = Y(2,IBLOK+IDN)
C          C1UP = Y(1,IBLOK+IUP)
C          C2UP = Y(2,IBLOK+IUP)
C          VERTD1 = CZUP*(C1UP - C1) - CZDN*(C1 - C1DN)
C          VERTD2 = CZUP*(C2UP - C2) - CZDN*(C2 - C2DN)
C Set horizontal diffusion and advection terms.
C          ILEFT = -1
C          IF (JX .EQ. 1) ILEFT = 1
C          IRIGHT = 1
C          IF (JX .EQ. MX) IRIGHT = -1
C          C1LT = Y(1,IBLOK+ILEFT)
C          C2LT = Y(2,IBLOK+ILEFT)
C          C1RT = Y(1,IBLOK+IRIGHT)
C          C2RT = Y(2,IBLOK+IRIGHT)
C          HORD1 = HDCO*(C1RT - 2.0E0*C1 + C1LT)
C          HORD2 = HDCO*(C2RT - 2.0E0*C2 + C2LT)
C          HORAD1 = HACO*(C1RT - C1LT)
C          HORAD2 = HACO*(C2RT - C2LT)
C Load all terms into YDOT.
C          YDOT(1,IBLOK) = VERTD1 + HORD1 + HORAD1 + RKIN1
C          YDOT(2,IBLOK) = VERTD2 + HORD2 + HORAD2 + RKIN2
C 10       CONTINUE
C 20     CONTINUE
C      RETURN
C      END
C
C      SUBROUTINE JACBD (F, NEQ, T, Y, YSV, REWT, F0, F1, HRL1,
C     1    BD, IPBD, IER, RPAR, IPAR)
C      EXTERNAL F
C      DIMENSION Y(2, *), YSV(*), REWT(*), F0(*), F1(*), BD(2, 2, *),
C     1          IPBD(2, *)
C      COMMON /PCOM/ MX,MZ,MM,Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO
C
C Compute diagonal Jacobian blocks, multiplied by -HRL1
C   (using q3 and q4 values computed on last F call).
C      DO 20 JZ = 1,MZ
C        ZDN = 30.0E0 + (JZ - 1.5E0)*DZ
C        ZUP = ZDN + DZ
C        CZDN = VDCO*EXP(0.2E0*ZDN)
C        CZUP = VDCO*EXP(0.2E0*ZUP)
C        DIAG = -(CZDN + CZUP + 2.0E0*HDCO)
C        IBLOK0 = (JZ-1)*MX
C        DO 10 JX = 1,MX
C          IBLOK = IBLOK0 + JX
C          C1 = Y(1,IBLOK)
C          C2 = Y(2,IBLOK)
C          BD(1,1,IBLOK) = -HRL1*( (-Q1*C3 - Q2*C2) + DIAG )
C          BD(1,2,IBLOK) = -HRL1*( -Q2*C1 + Q4 )
C          BD(2,1,IBLOK) = -HRL1*( Q1*C3 - Q2*C2 )
C          BD(2,2,IBLOK) = -HRL1*( (-Q2*C1 - Q4) + DIAG )
C 10       CONTINUE
C 20     CONTINUE
C Add identity matrix and do LU decompositions on blocks.
C      DO 40 IBLOK = 1,MM
C        BD(1,1,IBLOK) = BD(1,1,IBLOK) + 1.0E0
C        BD(2,2,IBLOK) = BD(2,2,IBLOK) + 1.0E0
C        CALL SGEFA (BD(1,1,IBLOK), 2, 2, IPBD(1,IBLOK), IER)
C        IF (IER .NE. 0) RETURN
C 40     CONTINUE
C      RETURN
C      END
C
C      SUBROUTINE SOLBD (NEQ, T, Y, F0, WK, HRL1, BD, IPBD, V, LR, IER,
C     1    RPAR, IPAR)
C      DIMENSION BD(2,2,*), IPBD(2,*), V(2,*)
C      COMMON /PCOM/ MX,MZ,MM,Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO
C Solve the block-diagonal system Px = v using LU factors stored in BD
C and pivot data in IPBD, and return the solution in V.
C      IER = 0
C      DO 10 I = 1,MM
C        CALL SGESL (BD(1,1,I), 2, 2, IPBD(1,I), V(1,I), 0)
C 10     CONTINUE
C      RETURN
C      END
C
C The output of this program, on a Cray-1 in single precision,
C is as follows:
C
C t =  7.20e+03     no. steps =  194   order =  5   stepsize =  1.17e+02
C  c1 (bot.left/middle/top rt.) =   1.047e+04   2.964e+04   1.119e+04
C  c2 (bot.left/middle/top rt.) =   2.527e+11   7.154e+11   2.700e+11
C
C t =  1.44e+04     no. steps =  227   order =  5   stepsize =  2.73e+02
C  c1 (bot.left/middle/top rt.) =   6.659e+06   5.316e+06   7.301e+06
C  c2 (bot.left/middle/top rt.) =   2.582e+11   2.057e+11   2.833e+11
C
C t =  2.16e+04     no. steps =  252   order =  5   stepsize =  4.21e+02
C  c1 (bot.left/middle/top rt.) =   2.665e+07   1.036e+07   2.931e+07
C  c2 (bot.left/middle/top rt.) =   2.993e+11   1.028e+11   3.313e+11
C
C t =  2.88e+04     no. steps =  291   order =  4   stepsize =  2.13e+02
C  c1 (bot.left/middle/top rt.) =   8.702e+06   1.292e+07   9.650e+06
C  c2 (bot.left/middle/top rt.) =   3.380e+11   5.029e+11   3.751e+11
C
C t =  3.60e+04     no. steps =  321   order =  5   stepsize =  9.90e+01
C  c1 (bot.left/middle/top rt.) =   1.404e+04   2.029e+04   1.561e+04
C  c2 (bot.left/middle/top rt.) =   3.387e+11   4.894e+11   3.765e+11
C
C t =  4.32e+04     no. steps =  374   order =  4   stepsize =  4.44e+02
C  c1 (bot.left/middle/top rt.) =  -5.457e-09  -4.365e-09  -6.182e-09
C  c2 (bot.left/middle/top rt.) =   3.382e+11   1.355e+11   3.804e+11
C
C t =  5.04e+04     no. steps =  393   order =  5   stepsize =  5.22e+02
C  c1 (bot.left/middle/top rt.) =   3.396e-12   2.798e-12   3.789e-12
C  c2 (bot.left/middle/top rt.) =   3.358e+11   4.930e+11   3.864e+11
C
C t =  5.76e+04     no. steps =  407   order =  5   stepsize =  3.54e+02
C  c1 (bot.left/middle/top rt.) =   7.738e-12   6.455e-12   8.598e-12
C  c2 (bot.left/middle/top rt.) =   3.320e+11   9.650e+11   3.909e+11
C
C t =  6.48e+04     no. steps =  419   order =  5   stepsize =  5.90e+02
C  c1 (bot.left/middle/top rt.) =  -2.018e-11  -1.680e-11  -2.243e-11
C  c2 (bot.left/middle/top rt.) =   3.313e+11   8.922e+11   3.963e+11
C
C t =  7.20e+04     no. steps =  432   order =  5   stepsize =  5.90e+02
C  c1 (bot.left/middle/top rt.) =  -2.837e-11  -2.345e-11  -3.166e-11
C  c2 (bot.left/middle/top rt.) =   3.330e+11   6.186e+11   4.039e+11
C
C t =  7.92e+04     no. steps =  444   order =  5   stepsize =  5.90e+02
C  c1 (bot.left/middle/top rt.) =  -4.861e-14  -4.433e-14  -5.162e-14
C  c2 (bot.left/middle/top rt.) =   3.334e+11   6.669e+11   4.120e+11
C
C t =  8.64e+04     no. steps =  456   order =  5   stepsize =  5.90e+02
C  c1 (bot.left/middle/top rt.) =   2.511e-15   2.071e-15   2.802e-15
C  c2 (bot.left/middle/top rt.) =   3.352e+11   9.107e+11   4.163e+11
C
C
C final statistics..
C RWORK size = 3861      IWORK size = 230
C number of steps        =  456     number of f evals.     = 1317
C number of prec. evals. =   82     number of prec. solves = 1226
C number of nonl. iters. =  571     number of lin. iters.  =  743
C average Krylov subspace dimension (NLI/NNI)  =  1.3012
C number of conv. failures..  nonlinear =  0  linear =  0
C-----------------------------------------------------------------------
C Full description of user interface to VODPK.
C
C The user interface to VODPK consists of the following parts.
C
C i.   The call sequence to subroutine VODPK, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      Following these descriptions are
C        * a description of optional input available through the
C          call sequence,
C        * a description of optional output (in the work arrays), and
C        * instructions for interrupting and restarting a solution.
C
C ii.  Descriptions of other routines in the VODPK package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      COMMON, and obtain specified derivatives of the solution y(t).
C
C iii. Descriptions of COMMON blocks to be declared in overlay
C      or similar environments.
C
C iv.  Description of two routines in the VODPK package, either of
C      which the user may replace with the user-s own version, if
C      desired.  These relate to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part i.  Call Sequence.
C
C The call sequence parameters used for input only are
C     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW,
C     JAC, PSOL, MF,
C and those used for both input and output are
C     Y, T, ISTATE.
C The work arrays RWORK and IWORK are also used for conditional and
C optional input and optional output.  (The term output here refers
C to the return from subroutine VODPK to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 in the input.
C
C The descriptions of the call arguments are as follows.
C
C F      = The name of the user-supplied subroutine defining the
C          ODE system.  The system must be put in the first-order
C          form dy/dt = f(t,y), where f is a vector-valued function
C          of the scalar t and the vector y.  Subroutine F is to
C          compute the function f.  It is to have the form
C               SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR, IFAIL)
C               DIMENSION Y(NEQ), YDOT(NEQ), RPAR(*), IPAR(*)
C          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
C          is output.  Y and YDOT are arrays of length NEQ.
C          (In the DIMENSION statement above, NEQ  can be replaced by
C          *  to make  Y  and  YDOT  assumed size arrays.)
C
C          The IFAIL argument is an output flag to indicate whether
C          f(t,y) could not be evaluated for the given arguments. 
C          The input value of IFAIL is 0 and a return value of 0 
C          indicates success.  The F routine should set IFAIL .ne. 0 
C          if it failed to evaluate this function.
C
C          Subroutine F should not alter Y or T.
C          F must be declared EXTERNAL in the calling program.
C
C          Subroutine F may access user-defined real and integer
C          work arrays RPAR and IPAR, which are to be dimensioned
C          in the user-s calling (main) program.
C
C          If quantities computed in the F routine are needed
C          externally to VODPK, an extra call to F should be made
C          for this purpose, for consistent and accurate results.
C          If only the derivative dy/dt is needed, use VINDY instead.
C
C NEQ    = The size of the ODE system (number of first order
C          ordinary differential equations).  Used only for input.
C          NEQ may not be increased during the problem, but
C          can be decreased (with ISTATE = 3 in the input).
C
C Y      = A real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 1), and only for output on other calls.
C          On the first call, Y must contain the vector of initial
C          values.  In the output, Y contains the computed solution
C          evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to
C          F, JAC, and PSOL.
C
C T      = The independent variable.  In the input, T is used only on
C          the first call, as the initial point of the integration.
C          In the output, after each call, T is the value at which a
C          computed solution Y is evaluated (usually the same as TOUT).
C          On an error return, T is the farthest point reached.
C
C TOUT   = The next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 1), TOUT may be equal
C          to T for one call, then should .ne. T for the next call.
C          For the initial T, an input value of TOUT .ne. T is used
C          in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal t interval, whose endpoints are
C          TCUR - HU and TCUR.  (See optional output, below, for
C          TCUR and HU.)
C
C ITOL   = An indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = A relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = An absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C          The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector e = (e(i)) of estimated local errors
C          in Y, according to an inequality of the form
C                      rms-norm of ( e(i)/EWT(i) )   .le.   1,
C          where       EWT(i) = RTOL(i)*abs(Y(i)) + ATOL(i),
C          and the rms-norm (root-mean-square norm) here is
C          rms-norm(v) = sqrt(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
C          is a vector of weights which must always be positive, and
C          the values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting
C          user-supplied routines for the setting of EWT and/or for
C          the norm calculation.  See Part iv below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = An index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note..  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at T = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          the state of the calculation.
C
C          In the input, the values of ISTATE are as follows.
C          1  means this is the first call for the problem
C             (initializations will be done).  See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
C             and any of the optional input except H0.
C
C          Note..  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful to include
C          the initial conditions in the output.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 1 in the input.
C
C          In the output, ISTATE has the following values and meanings.
C           1  means nothing was done, as TOUT was equal to T with
C              ISTATE = 1 in the input.
C           2  means the integration was performed successfully.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again.
C              (The excess work step counter will be reset to 0.)
C              In addition, the user may increase MXSTEP to avoid
C              this error return.  (See optional input below.)
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note.. If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note..  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by a poor preconditioner matrix.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i)=0.0)
C              was requested on a variable which has now vanished.
C              The integration was successful as far as T.
C          -7  means an unrecoverable error occurred in JAC or PSOL.
C              Either JAC returned IER .ne. 0, or PSOL returned
C              IER .lt. 0.
C          -8  means the F routine returned IFAIL .ne. 0, and the
C              integrator was unable to avoid this condition.
C              The integration was successful as far as T.
C
C          Note..  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other input, before
C          calling the solver again.
C
C IOPT   = An integer flag to specify whether or not any optional
C          input is being used on this call.  Input only.
C          The optional input is listed separately below.
C          IOPT = 0 means no optional input is being used.
C                   Default values will be used in all cases.
C          IOPT = 1 means optional input is being used.
C
C RWORK  = A real working array (single precision).
C          The length of RWORK must be at least
C             20 + NYH*(MAXORD + 1) + 3*NEQ + LENK + LWP   where
C          NYH    = the initial value of NEQ,
C          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
C                   smaller value is given as an optional input),
C          LENK = length of work space for Krylov-related data..
C          LENK = 0                                 if MITER = 0,
C          LENK = NEQ*(MAXL+3+MIN(1,MAXL-KMP))
C                  + (MAXL+3)*MAXL + 1              if MITER = 1,
C          LENK = 3*NEQ                             if MITER = 9.
C          LWP = length of real user work space for preconditioning.
C          (See JAC/PSOL.)
C          (See the MF description for METH and MITER.)
C          Thus if MAXORD has its default value and NEQ is constant,
C          this length is..
C             20 + 16*NEQ                    for MF = 10,
C             61 + 24*NEQ + LWP              for MF = 11,
C             20 + 19*NEQ + LWP              for MF = 19,
C             20 + 9*NEQ                     for MF = 20,
C             61 + 17*NEQ + LWP              for MF = 21,
C             20 + 12*NEQ + LWP              for MF = 29
C          The first 20 words of RWORK are reserved for conditional
C          and optional input and optional output.
C
C          The following word in RWORK is a conditional input..
C            RWORK(1) = TCRIT = critical value of t which the solver
C                       is not to overshoot.  Required if ITASK is
C                       4 or 5, and ignored otherwise.  (See ITASK.)
C
C LRW    = The length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = An integer work array.  The length of IWORK must be at least
C             30        if MITER = 0  (MF = 10, 20), or
C             30 + LIWP  otherwise (MF = 11, 21, 19, 29).
C          LIWP = length of integer user work space for preconditioning.
C          (See conditional input list following).
C
C          The first 30 words of IWORK are reserved for conditional and
C          optional input and optional output.
C
C          The following 4 words in IWORK are conditional input,
C          required if MITER .ge. 1..
C
C          IWORK(1) = LWP  = length of real array WP for use in
C                     preconditioning (part of RWORK array).
C          IWORK(2) = LIWP = length of integer array IWP for use in
C                     preconditioning (part of IWORK array).
C                     The arrays WP and IWP are work arrays under the
C                     user-s control, for use in the routines that
C                     perform preconditioning operations (JAC and PSOL).
C          IWORK(3) = JPRE = preconditioner type flag..
C                   = 0 for no preconditioning (P1 = P2 = I
C                   = 1 for left-only preconditioning (P2 = I)
C                   = 2 for right-only preconditioning (P1 = I)
C                   = 3 for two-sided preconditioning
C          IWORK(4) = JACFLG = flag for whether JAC is called.
C                   = 0 if JAC is not to be called,
C                   = 1 if JAC is to be called.
C                     Use JACFLG = 1 if JAC computes any nonconstant
C                     data needed in preconditioning operations,
C                     such as some of the Jacobian elements.
C
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note..  The work arrays must not be altered between calls to VODPK
C for the same problem, except possibly for the conditional and
C optional input, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside VODPK between calls, if
C desired (but not for use by F or JAC).
C
C JAC    = The name of the user-supplied routine (MITER = 1 or 4) to
C          compute the Jacobian matrix, df/dy, as a function of
C          the scalar t and the vector y.  It is to have the form
C
C               SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY, V, HRL1,
C                               WP, IWP, IER, RPAR, IPAR)
C               EXTERNAL F
C               DIMENSION Y(NEQ), YSV(NEQ), REWT(NEQ), FTY(NEQ),
C                         V(NEQ), WP(*), IWP(*), RPAR(*), IPAR(*)
C
C          This routine must evaluate and preprocess any parts of the
C          Jacobian matrix df/dy used in the preconditioners P1, P2 .
C          The Y and FTY arrays contain the current values of y and
C          f(t,y), respectively, and YSV also contains the current
C          value of y.  The array V is work space of length
C          NEQ for use by JAC.  REWT is the array of reciprocal error
C          weights (1/ewt).  JAC must multiply all computed Jacobian
C          elements by the scalar -hrl1, add the identity matrix I and
C          do any factorization operations called for, in preparation
C          for solving linear systems with a coefficient matrix of
C          P1 or P2.  The matrix P1*P2 should be an approximation to
C          I - hrl1 * (df/dy), where hrl1 is stored in HRL1.  JAC should
C          return IER = 0 if successful, and IER .ne. 0 if not.
C          (If IER .ne. 0, a smaller time step will be tried.)
C          The arrays WP (of length LWP) and IWP (of length LIWP)
C          are for use by JAC and PSOL for work space and for storage
C          of data needed for the solution of the preconditioner
C          linear systems.  Their lengths and contents are under the
C          user-s control.
C          The JAC routine may save relevant Jacobian elements (or
C          approximations) used in the preconditioners, along with the
C          value of hrl1, and use these to reconstruct preconditioner
C          matrices later without reevaluationg those elements.
C          This may be cost-effective if JAC is called with hrl1
C          considerably different from its earlier value, indicating
C          that a corrector convergence failure has occurred because
C          of the change in hrl1, not because of changes in the
C          value of the Jacobian.  In doing this, use the saved and
C          current values of hrl1 to decide whether to use saved
C          or reevaluated elements.
C          JAC may alter V, but not Y, YSV, REWT, FTY, or HRL1.
C          JAC must be declared external in the calling program.
C
C PSOL   = the name of the user-supplied routine for the
C          solution of preconditioner linear systems.
C          It is to have the form
C            SUBROUTINE PSOL (NEQ, T, Y, FTY, WK, HRL1, WP, IWP, B, LR,
C                             IER, RPAR, IPAR)
C            DIMENSION Y(NEQ), FTY(NEQ), WK(NEQ), WP(*), IWP(*),
C                      B(NEQ), RPAR(*), IPAR(*)
C          This routine must solve a linear system with b (stored in B)
C          as right-hand side and one of the preconditioning matrices,
C          P1 or P2, as coefficient matrix, and return the solution
C          vector in B.  LR is a flag concerning left vs. right
C          preconditioning, input to PSOL.  PSOL is to use P1 if LR = 1
C          and P2 if LR = 2.  In the case miter = 9 (no Krylov
C          iteration), LR will be 0, and PSOL is to return in B the
C          desired approximate solution to A * x = b, where
C          A = I - hrl1 * (df/dy).  (hrl1 is stored in HRL1.)  PSOL can
C          use data generated in the JAC routine and stored in WP and
C          IWP.  The Y and FTY arrays contain the current values of y
C          and f(t,y), respectively.  The array WK is work space of
C          length NEQ for use by PSOL.
C          The argument HRL1 is the current value of the scalar appear-
C          ing in the linear system.  If the old value, as of the last
C          JAC call, is needed, it must have been saved by JAC in WP.
C          On return, PSOL should set the error flag IER as follows..
C            IER = 0 if PSOL was successful,
C            IER .gt. 0 on a recoverable error, meaning that the
C                   time step will be retried,
C            IER .lt. 0 on an unrecoverable error, meaning that the
C                   solver is to stop immediately.
C          PSOL may not alter Y, FTY, or HRL1.
C          PSOL must be declared external in the calling program.
C
C MF     = The method flag.  Used only for input.  The legal values of
C          MF are 10, 11, 19, 20, 21, 29 .
C          MF is a two-digit integer, MF = 10*METH + MITER .
C          METH indicates the basic linear multistep method..
C            METH = 1 means the implicit Adams method.
C            METH = 2 means the method based on backward
C                     differentiation formulas (BDF-s).
C          MITER indicates the corrector iteration method.  Currently,
C            the values of MITER have the following meanings..
C
C          0 means functional iteration is used (no Jacobian matrix
C            is involved).
C
C          1 means SPIGMR, a scaled, preconditioned, incomplete version
C            of GMRES, a generalized minimum residual method, is used.
C            This is the best choice in general.
C
C          9 means that only a user-supplied matrix P (approximating A)
C            will be used, with no Krylov iteration done internally to
C            VODPK.  This option allows the user to provide the complete
C            linear system solution algorithm, if desired.
C
C The user can apply preconditioning to the linear system A*x = b,
C by means of arbitrary matrices (the preconditioners).
C
C RPAR     User-specified array used to communicate real parameters
C          to user-supplied subroutines.  If RPAR is a vector, then
C          it must be dimensioned in the user's main program.  If it
C          is unused or it is a scalar, then it need not be
C          dimensioned.
C
C IPAR     User-specified array used to communicate integer parameter
C          to user-supplied subroutines.  The comments on dimensioning
C          RPAR apply to IPAR.
C-----------------------------------------------------------------------
C Optional Input.
C
C The following is a list of the optional input provided for in the
C call sequence.  (See also Part ii.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of this input requires IOPT = 1, and in that
C case all of this input is examined.  A value of zero for any
C of these optional input variables will cause the default value to be
C used.  Thus to use a subset of the optional input, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C NAME    LOCATION      MEANING AND DEFAULT VALUE
C
C H0      RWORK(5)  The step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  The maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  The minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C MAXORD  IWORK(5)  The maximum order to be allowed.  The default
C                   value is 12 if METH = 1, and 5 if METH = 2.
C                   If MAXORD exceeds the default value, it will
C                   be reduced to the default value.
C                   If MAXORD is changed during the problem, it may
C                   cause the current order to be reduced.
C
C MXSTEP  IWORK(6)  Maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  Maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C
C MAXL    IWORK(8)  maximum number of iterations in the SPIGMR
C                   algorithm (.le. NEQ).  The default is
C                   MAXL = min(5, NEQ).
C
C KMP     IWORK(9)  number of vectors on which orthogonalization
C                   is done in the SPIGMR algorithm (.le. MAXL).
C                   The default is KMP = MAXL (complete GMRES method).
C                   See Ref. 2 for details on incomplete GMRES.
C                   Note:  When KMP .lt. MAXL and MITER = 1, the length
C                          of RWORK must be defined accordingly.  See
C                          the definition of RWORK above.
C-----------------------------------------------------------------------
C Optional Output.
C
C As optional additional output from VODPK, the variables listed
C below are quantities related to the performance of VODPK
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C Except where stated otherwise, all of this output is defined
C on any successful return from VODPK, and on any return with
C ISTATE = -1, -2, -4, -5, ..., -8.  On an illegal input return
C (ISTATE = -3), they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, output relevant to the error will be defined,
C as noted below.
C
C NAME    LOCATION      MEANING
C
C HU      RWORK(11) The step size in t last used (successfully).
C
C HCUR    RWORK(12) The step size to be attempted on the next step.
C
C TCUR    RWORK(13) The current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  In the output,
C                   TCUR will always be at least as far from the
C                   initial value of t as the current argument T,
C                   but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) A tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C NST     IWORK(11) The number of steps taken for the problem so far.
C
C NFE     IWORK(12) The number of f evaluations for the problem so far.
C
C NPE     IWORK(13) The number of preconditioner evaluations (JAC calls)
C                   so far.
C
C NQU     IWORK(14) The method order last used (successfully).
C
C NQCUR   IWORK(15) The order to be attempted on the next step.
C
C IMXER   IWORK(16) The index of the component of largest magnitude in
C                   the weighted local error vector ( e(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) The length of RWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) The length of IWORK actually required.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C NNI     IWORK(20) The number of nonlinear iterations so far (each of
C                   which calls the Krylov iterative linear solver).
C
C NCFN    IWORK(21) The number of convergence failures of the nonlinear
C                   (Newton) iteration so far.
C                   Note.. A measure of success is the overall
C                   rate of nonlinear convergence failures, NCFN/NST.
C
C NETF    IWORK(22) The number of error test failures of the integrator
C                   so far.
C
C
C NLI     IWORK(23) The number of linear iterations so far.
C                   Note.. a measure of the success of algorithm is
C                   the average number of linear iterations per
C                   nonlinear iteration, given by NLI/NNI.
C                   If this is close to MAXL, MAXL may be too small.
C
C NPS     IWORK(24) The number of preconditioning solve operations
C                   (PSOL calls) so far.
C
C
C NCFL    IWORK(25) The number of convergence failures of the linear
C                   iteration so far.
C                   Note.. A measure of success is the overall
C                   rate of linear convergence failures, NCFL/NNI.
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional output.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C NAME    BASE ADDRESS      DESCRIPTION
C
C YH      21             The Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the
C                        solution, evaluated at t = TCUR.
C
C ACOR     LENRW-NEQ+1   Array of size NEQ used for the accumulated
C                        corrections on each step, scaled in the output
C                        to represent the estimated local error in Y
C                        on the last step.  This is the vector e in
C                        the description of the error control.  It is
C                        defined only on a successful return from VODPK.
C
C-----------------------------------------------------------------------
C Interrupting and Restarting
C
C If the integration of a given problem by VODPK is to be
C interrrupted and then later continued, such as when restarting
C an interrupted run or alternating between two or more ODE problems,
C the user should save, following the return from the last VODPK call
C prior to the interruption, the contents of the call sequence
C variables and internal COMMON blocks, and later restore these
C values before the next VODPK call for that problem.  To save
C and restore the COMMON blocks, use subroutine VKSRC, as
C described below in part ii.
C
C-----------------------------------------------------------------------
C Part ii.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with VODPK.
C
C     FORM OF CALL                  FUNCTION
C
C  CALL VKSRC(RSAV,ISAV,JOB) Saves and restores the contents of
C                             the internal COMMON blocks used by
C                             VODPK. (See Part iii below.)
C                             RSAV must be a real array of length 52
C                             or more, and ISAV must be an integer
C                             array of length 52 or more.
C                             JOB=1 means save COMMON into RSAV/ISAV.
C                             JOB=2 means restore COMMON from RSAV/ISAV.
C
C                                VKSRC is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with VODPK.
C
C  CALL VINDY(,,,,,)         Provide derivatives of y, of various
C        (See below.)         orders, at a specified point T, if
C                             desired.  It may be called only after
C                             a successful return from VODPK.
C
C The detailed instructions for using VINDY are as follows.
C The form of the call is..
C
C  CALL VINDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C The input parameters are..
C
C T         = Value of independent variable where answers are desired
C             (normally the same as the T last returned by VODPK).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional output for TCUR and HU.)
C K         = Integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional output).  The capability corresponding
C             to K = 0, i.e. computing y(T), is already provided
C             by VODPK directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with VINDY.
C RWORK(21) = The base address of the history array YH.
C NYH       = Column length of YH, equal to the initial value of NEQ.
C
C The output parameters are..
C
C DKY       = A real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = Integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part iii.  COMMON Blocks.
C If VODPK is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in..
C   (1) the call sequence to VODPK,
C   (2) the two internal COMMON blocks
C         /VOD001/  of length  81  (48 single precision words
C                         followed by 33 integer words),
C         /VOD002/  of length  9  (1 single precision word
C                         followed by 8 integer words),
C
C         /VPK001/  of length 14 (3 single precision words
C                         followed by 11 integer words)
C
C If VODPK is used on a system in which the contents of internal
C COMMON blocks are not preserved between calls, the user should
C declare the above three COMMON blocks in the calling (main) program
C to insure that their contents are preserved.
C
C-----------------------------------------------------------------------
C Part iv.  Optionally Replaceable Solver Routines.
C
C Below are descriptions of two routines in the VODPK package which
C relate to the measurement of errors.  Either routine can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note.. The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) EWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above..
C     SUBROUTINE EWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the VODPK call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by EWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparison with
C errors in Y(i).  The EWT array returned by EWSET is passed to the
C VNORML routine (see below), and also used by VODPK in the computation
C of the optional output IMXER, the diagonal Jacobian approximation,
C and the increments for difference quotient Jacobians.
C
C In the user-supplied version of EWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C Optional Output.  In EWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of h**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in EWSET the statements..
C     COMMON /VOD001/ RVOD(48), IVOD(33)
C     COMMON /VOD002/ HU, NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
C     NQ = IVOD(28)
C     H = RVOD(21)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C
C (b) VNORML.
C The following is a real function routine which computes the weighted
C root-mean-square norm of a vector v..
C     D = VNORML (N, V, W)
C where..
C   N = the length of the vector,
C   V = real array of length N containing the vector,
C   W = real array of length N containing weights,
C   D = sqrt( (1/N) * sum(V(i)*W(i))**2 ).
C VNORML is called with N = NEQ and with W(i) = 1.0/EWT(i), where
C EWT is as set by subroutine EWSET.
C
C If the user supplies this function, it should return a non-negative
C value of VNORML suitable for use in the error control in VODPK.
C None of the arguments should be altered by VNORML.
C For example, a user-supplied VNORML routine might..
C   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
C   -ignore some components of V in the norm, with the effect of
C    suppressing the error control on those components of Y.
C-----------------------------------------------------------------------
C Other Routines in the VODPK Package.
C
C In addition to subroutine VODPK, the VODPK package includes the
C following subroutines and function routines..
C  VHIN     computes an approximate step size for the initial step.
C  VINDY    computes an interpolated value of the y vector at t = TOUT.
C  VSTEP    is the core integrator, which does one step of the
C           integration and the associated error control.
C  VSET     sets all method coefficients and test constants.
C  VJUST    adjusts the history array on a change of order.
C  VNLSK    solves the underlying nonlinear system -- the corrector.
C  VSOLPK   manages solution of linear system in chord iteration.
C  VSPIG    performs the SPIGMR algorithm.
C  VATV     computes a scaled, preconditioned product (I-hrl1*J)*v.
C  SVRORTHOG   orthogonalizes a vector against previous basis vectors.
C  SHEQR    generates a QR factorization of a Hessenberg matrix.
C  SHELS    finds the least squares solution of a Hessenberg system.
C  VUSOL    interfaces to the user-s PSOL routine (MITER = 9).
C  EWSET    sets the error weight vector EWT before each step.
C  VNORML    computes the weighted r.m.s. norm of a vector.
C  VKSRC    is a user-callable routine to save and restore
C           the contents of the internal COMMON blocks.
C  SAXPY, SSCAL, SCOPY, SDOT, and SNRM2 are basic linear
C           algebra modules (BLAS) used by this package.
C  R1MACH9   computes the unit roundoff in a machine-independent manner.
C  XERRWV   handles the printing of all error messages and warnings.
C           XERRWV is machine-dependent.

C Note..  VNORML, SDOT, SNRM2, and R1MACH9 are
C function routines.  All the others are subroutines.
C
C The intrinsic and external routines used by the VODPK package are..
C ABS, float, max, min, SIGN, SQRT, ijmgetmr, xerrab, and write.
C
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block VOD001 --------------------
C
      real ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block VOD002 --------------------
C
      real HU
      INTEGER NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
C
C Type declarations for labeled COMMON block VPK001 --------------------
C
      real DELT, SQRTN, RSQRTN
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LVSAV, KMP, MAXL, MNEWT,
     1      NLI, NPS, NCFL
C Type declarations for local variables --------------------------------
C
      EXTERNAL VNLSK
      LOGICAL IHIT, LAVD, LCFN, LCFL, LWARN
      real ATOLI, AVDIM, BIG, EWTI, FOUR, H0, HMAX, HMX,
     1   HUN, ONE, PT05, PT2, PT9, RCFL, RCFN, RH, RTOLI, SIZE,
     2   TCRIT, TNEXT, TOLSF, TP, TWO, ZERO
      INTEGER I, IER, IFLAG, IMXER, KGO, LENIW, LENIWK, LENRW, LENWK,
     1   LENWM, LF0, LIWP, LWP, MORD, MXHNL0, MXSTP0, NCFL0, NCFN0,
     2   NITER, NLI0, NNI0, NNID, NSTD, NSLAST, NWARN, IFAIL
      CHARACTER*80 MSG
C
C Type declaration for function subroutines called ---------------------
C
      real R1MACH9, VNORML
C
      DIMENSION MORD(2)
C-----------------------------------------------------------------------
C The following Fortran-77 declarations are to cause the values of the
C listed (local) variables to be saved between calls to VODPK.
C-----------------------------------------------------------------------
      SAVE MORD, MXHNL0, MXSTP0
      SAVE ZERO, ONE, TWO, FOUR, HUN, PT05, PT2, PT9
C-----------------------------------------------------------------------
C The following internal COMMON blocks contain variables which are
C communicated between subroutines in the VODPK package, or which are
C to be saved between calls to VODPK.
C In each block, real variables precede integers.
C The block /VOD001/ appears in subroutines VODPK, VINDY, VSTEP, VSET,
C VJUST, VNLSK, VSOLPK, VATV, and VKSRC.
C The block /VOD002/ appears in subroutines VODPK, VINDY, VSTEP, VNLSK,
C VSOLPK, VATV, and VKSRC.
C The block /VPK001/ appears in subroutines VODPK, VNLSK, VSOLPK,
C and VKSRC.
C
C The variables stored in the internal COMMON blocks are as follows..
C
C ACNRM  = Weighted r.m.s. norm of accumulated correction vectors.
C CCMXJ  = Threshhold on DRC for updating the Jacobian. (See DRC.)
C CONP   = The saved value of TQ(5).
C CRATE  = Estimated corrector convergence rate constant.
C DRC    = Relative change in H*RL1 since last VJAC call.
C EL     = Real array of integration coefficients.  See VSET.
C ETA    = Saved tentative ratio of new to old H.
C ETAMAX = Saved maximum value of ETA to be allowed.
C H      = The step size.
C HMIN   = The minimum absolute value of the step size H to be used.
C HMXI   = Inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C HNEW   = The step size to be attempted on the next step.
C HSCAL  = Stepsize in scaling of YH array.
C PRL1   = The saved value of RL1.
C RC     = Ratio of current H*RL1 to value on last VJAC call.
C RL1    = The reciprocal of the coefficient EL(1).
C TAU    = Real vector of past NQ step sizes, length 13.
C TQ     = A real vector of length 5 in which VSET stores constants
C          used for the convergence test, the error test, and the
C          selection of H at a new order.
C TN     = The independent variable, updated on each step taken.
C UROUND = The machine unit roundoff.  The smallest positive real number
C          such that  1.0 + UROUND .ne. 1.0
C ICF    = Integer flag for convergence failure in VNLSK..
C            0 means no failures.
C            1 means convergence failure with out of date Jacobian
C                   (recoverable error).
C            2 means convergence failure with current Jacobian or
C                   singular matrix (unrecoverable error).
C INIT   = Saved integer flag indicating whether initialization of the
C          problem has been done (INIT = 1) or not.
C IPUP   = Saved flag to signal updating of Newton matrix.
C JCUR   = Output flag from VJAC showing Jacobian status..
C            JCUR = 0 means J is not current.
C            JCUR = 1 means J is current.
C JSTART = Integer flag used as input to VSTEP..
C            0  means perform the first step.
C            1  means take a new step continuing from the last.
C            -1 means take the next step with a new value of MAXORD,
C                  HMIN, HMXI, N, METH, MITER, and/or matrix parameters.
C          On return, VSTEP sets JSTART = 1.
C JSV    = Integer flag for Jacobian saving, = sign(MF).
C KFLAG  = A completion code from VSTEP with the following meanings..
C               0      the step was succesful.
C              -1      the requested error could not be achieved.
C              -2      corrector convergence could not be achieved.
C              -3, -4  fatal error in VNLS (can not occur here).
C              -5      corrector loop unable to avoid failure in F.
C KUTH   = Input flag to VSTEP showing whether H was reduced by the
C          driver.  KUTH = 1 if H was reduced, = 0 otherwise.
C L      = Integer variable, NQ + 1, current order plus one.
C LMAX   = MAXORD + 1 (used for dimensioning).
C LOCJS  = A pointer to the saved Jacobian, whose storage starts at
C          WM(LOCJS), if JSV = 1.
C LYH, LEWT, LACOR, LSAVF, LWM, LIWM = Saved integer pointers
C          to segments of RWORK and IWORK.
C MAXORD = The maximum order of integration method to be allowed.
C METH/MITER = The method flags.  See MF.
C MSBJ   = The maximum number of steps between J evaluations, = 50.
C MXHNIL = Saved value of optional input MXHNIL.
C MXSTEP = Saved value of optional input MXSTEP.
C N      = The number of first-order ODEs, = NEQ.
C NEWH   = Saved integer to flag change of H.
C NEWQ   = The method order to be used on the next step.
C NHNIL  = Saved counter for occurrences of T + H = T.
C NQ     = Integer variable, the current integration method order.
C NQNYH  = Saved value of NQ*NYH.
C NQWAIT = A counter controlling the frequency of order changes.
C          An order change is about to be considered if NQWAIT = 1.
C NSLJ   = The number of steps taken as of the last Jacobian update.
C NSLP   = Saved value of NST as of last Newton matrix update.
C NYH    = Saved value of the initial value of NEQ.
C
C HU     = The step size in t last used.
C NCFN   = Number of nonlinear convergence failures so far.
C NETF   = The number of error test failures of the integrator so far.
C NFE    = The number of f evaluations for the problem so far.
C NPE    = The number of preconditioner evaluations (JAC calls) so far.
C NLU    = The number of matrix LU decompositions so far.
C NNI    = Number of nonlinear iterations so far.
C NQU    = The method order last used.
C NST    = The number of steps taken for the problem so far.
C
C DELT   = Convergence test constant in Krylov iterations.
C SQRTN  = SQRT(NEQ), for use in weights in Krylov convergence tests.
C RSQRTN = 1.0/SQRTN, also for use in convergence weights.
C JPRE   = Preconditioner type flag.
C JACFLG = Indicator for presence of user-supplied JAC routine.
C LOCWP  = Location of start of user-s WP array in WM work array.
C LOCIWP = Location of start of user-s IWP array in IWM work array.
C LVSAV  = Saved pointer to VSAV array in RWORK.
C KMP    = Number of vectors on which orthogonalization is done in
C          Krylov iteration.
C MAXL   = Maximum dimension of Krylov subspace used.
C MNEWT  = Newton iteration index.
C NLI    = Number of linear (Krylov) iterations done.
C NPS    = Number of preconditioner solvers (PSOL calls) done.
C NCFL   = Number of convergence failures in Krylov iteration.
C-----------------------------------------------------------------------
      COMMON /VOD001/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /VOD002/ HU, NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
      COMMON /VPK001/ DELT, SQRTN, RSQRTN, JPRE, JACFLG, LOCIWP,
     1                LOCWP, LVSAV, KMP, MAXL, MNEWT, NLI, NPS, NCFL
C
      DATA  MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
      DATA ZERO/0.0E0/, ONE/1.0E0/, TWO/2.0E0/, FOUR/4.0E0/,
     1   PT05/0.05E0/, PT2/0.2E0/, PT9/0.9E0/, HUN/100.0E0/
      real scaleps
      common /comvatv/ scaleps
c
c     Special scaling factor for Newton iteration error TDR 11/19/93
      real EFACN
      COMMON /TDR1/ EFACN
c
      scaleps = scalpk
      EFACN = efacnin
C
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, jump to Block G and return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT, MF.
C-----------------------------------------------------------------------
 20   IF (NEQ .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ .GT. N) GO TO 605
 25   N = NEQ
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      JSV = 0
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0) GO TO 608
      IF (MITER .GT. 1 .AND. MITER .LT. 9) GO TO 608
      IF (MITER .GE. 1) JPRE = IWORK(3)
      JACFLG = 0
      IF (MITER .GE. 1) JACFLG = IWORK(4)
      IF (MITER .GE. 1 .AND. MITER .NE. 9 .AND. JPRE .EQ. 0) JACFLG = 0
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = ZERO
      HMXI = ZERO
      HMIN = ZERO
      MAXL = min(5,N)
      KMP = MAXL
      DELT = PT05
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = min(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. ZERO) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. ZERO) GO TO 615
      HMXI = ZERO
      IF (HMAX .GT. ZERO) HMXI = ONE/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. ZERO) GO TO 616
      MAXL = IWORK(8)
      IF (MAXL .EQ. 0) MAXL = 5
      MAXL = min(MAXL,N)
      KMP = IWORK(9)
      IF (KMP .EQ. 0 .OR. KMP .GT. MAXL) KMP = MAXL
      DELT = RWORK(8)
      IF (DELT .EQ. 0.0E0) DELT = PT05
C-----------------------------------------------------------------------
C Set work array pointers and check lengths lrw and liw.
C Pointers to segments of RWORK and iwork are named by prefixing l to
C the name of the segment.  e.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are  YH, WM, EWT, SAVF, VSAV, ACOR.
C Within WM, LOCWP is the location of the WP work array,
C and within IWM, LOCIWP is the location of the IWP work array.
C-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .EQ. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      IF (MITER .EQ. 0) LENWK = 0
      IF (MITER .EQ. 1)
     1   LENWK = N*(MAXL+2+min(1,MAXL-KMP)) + (MAXL+3)*MAXL + 1
      IF (MITER .EQ. 9) LENWK = 2*N
      LWP = 0
      IF (MITER .GE. 1) LWP = IWORK(1)
      LENWM = LENWK + LWP
      LOCWP = LENWK + 1
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LVSAV = LSAVF + N
      LACOR = LVSAV + N
      IF (MITER .EQ. 0) LACOR = LVSAV
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 31
      LENIWK = 0
      LIWP = 0
      IF (MITER .GE. 1) LIWP = IWORK(2)
      LENIW = 30 + LENIWK + LIWP
      LOCIWP = LENIWK + 1
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1, N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. ZERO) GO TO 619
        IF (ATOLI .LT. ZERO) GO TO 620
 70     CONTINUE
C Load SQRT(N) and its reciprocal in common. ---------------------------
      SQRTN = SQRT(float(N))
      RSQRTN = ONE/SQRTN
      IF (ISTATE .EQ. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to VSTEP. --------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 100
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      CALL SCOPY (N, RWORK(LWM), 1, RWORK(LSAVF), 1)
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = R1MACH9(4)
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. ZERO) GO TO 625
      IF (H0 .NE. ZERO .AND. (T + H0 - TCRIT)*H0 .GT. ZERO)
     1   H0 = TCRIT - T
 110  JSTART = 0
      CCMXJ = PT2
      MSBJ = 50
      NHNIL = 0
      NST = 0
      NSLAST = 0
      HU = ZERO
      NQU = 0
      NPE = 0
      NLI0 = 0
      NNI0 = 0
      NCFN0 = 0
      NCFL0 = 0
      NWARN = 0
      NNI = 0
      NLI = 0
      NPS = 0
      NETF = 0
      NCFN = 0
      NCFL = 0
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      IFAIL = 0
      CALL F (N, T, Y, RWORK(LF0), RPAR, IPAR, IFAIL)
      NFE = 1
      IF (IFAIL .NE. 0) THEN
        MSG = 'VODPK-- Initial call to F failed.       '
        CALL XERRWV(MSG, 40, 251, 0, 0, 0, 0, 0, ZERO, ZERO)
        ISTATE = -8
        GO TO 580
        ENDIF
C Load the initial value vector in YH. ---------------------------------
      CALL SCOPY (N, Y, 1, RWORK(LYH), 1)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = ONE
      CALL EWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1, N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 621
 120    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
      IF (H0 .NE. ZERO) GO TO 180
C Call VHIN to set initial step size H0 to be attempted. ---------------
      CALL VHIN (N, T, RWORK(LYH), RWORK(LF0), F, RPAR, IPAR, TOUT,
     1   UROUND, RWORK(LEWT), ITOL, ATOL, Y, RWORK(LACOR), H0,
     2   NITER, IER)
      NFE = NFE + NITER
      IF (IER .EQ. -1) GO TO 622
      IF (IER .EQ. -2) THEN
        MSG = 'VODPK--  calls to F failed repeatedly in initial H calc.'
        CALL XERRWV(MSG, 56, 251, 0, 0, 0, 0, 0, ZERO, ZERO)
        ISTATE = -8
        GO TO 580
        ENDIF
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. ONE) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      CALL SSCAL (N, H0, RWORK(LF0), 1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      KUTH = 0
      NLI0 = NLI
      NNI0 = NNI
      NCFN0 = NCFN
      NCFL0 = NCFL
      NWARN = 0
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL VINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(ONE + HUN*UROUND)
      IF ((TP - TOUT)*H .GT. ZERO) GO TO 623
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. ZERO) GO TO 625
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 245
      CALL VINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator VSTEP.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken,
C check for poor Newton/Krylov performance, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL EWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      NSTD = NST - NSLAST
      NNID = NNI - NNI0
      IF (NSTD .LT. 10 .OR. NNID .EQ. 0) GO TO 255
      AVDIM = float(NLI - NLI0)/float(NNID)
      RCFN = float(NCFN - NCFN0)/float(NSTD)
      RCFL = float(NCFL - NCFL0)/float(NNID)
      LAVD = AVDIM .GT. (float(MAXL) - PT05)
      LCFN = RCFN .GT. PT9
      LCFL = RCFL .GT. PT9
      LWARN = LAVD .OR. LCFN .OR. LCFL
      IF (.NOT.LWARN) GO TO 255
      NWARN = NWARN + 1
      IF (NWARN .GT. 10) GO TO 255
      IF (LAVD) THEN
        MSG = 'VODPK-- Warning. Poor iterative algorithm performance   '
        CALL XERRWV (MSG, 56, 111, 0, 0, 0, 0, 0, ZERO, ZERO)
        MSG = '      at T = R1. Average no. of linear iterations = R2  '
        CALL XERRWV (MSG, 56, 111, 0, 0, 0, 0, 2, TN, AVDIM)
        ENDIF
      IF (LCFN) THEN
        MSG = 'VODPK-- Warning. Poor iterative algorithm performance   '
        CALL XERRWV (MSG, 56, 112, 0, 0, 0, 0, 0, ZERO, ZERO)
        MSG = '      at T = R1. Nonlinear convergence failure rate = R2'
        CALL XERRWV (MSG, 56, 112, 0, 0, 0, 0, 2, TN, RCFN)
        ENDIF
      IF (LCFL) THEN
        MSG = 'VODPK-- Warning. Poor iterative algorithm performance   '
        CALL XERRWV (MSG, 56, 113, 0, 0, 0, 0, 0, ZERO, ZERO)
        MSG = '      at T = R1. Linear convergence failure rate = R2   '
        CALL XERRWV (MSG, 56, 113, 0, 0, 0, 0, 2, TN, RCFL)
        ENDIF
 255  CONTINUE
      DO 260 I = 1, N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 510
 260    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*VNORML (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. ONE) GO TO 280
      TOLSF = TOLSF*TWO
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'VODPK--  Warning..internal T (=R1) and H (=R2) are'
      CALL XERRWV (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      (H = step size). solver will continue anyway'
      CALL XERRWV (MSG, 50, 101, 1, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'VODPK--  Above warning has been issued I1 times.  '
      CALL XERRWV (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      it will not be issued again for this problem'
      CALL XERRWV (MSG, 50, 102, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
 290  CONTINUE
C-----------------------------------------------------------------------
C  CALL VSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR, WM, IWM,
C                                     F, JAC, PSOL, VNLSK, RPAR, IPAR)
C-----------------------------------------------------------------------
      CALL VSTEP (Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LVSAV), RWORK(LACOR), RWORK(LWM),
     2   IWORK(LIWM), F, JAC, PSOL, VNLSK, RPAR, IPAR)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540, 550, 552, 555), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      KUTH = 0
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  if TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL VINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. ZERO) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. ZERO) GO TO 345
      CALL VINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from VODPK.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  CONTINUE
      CALL SCOPY (N, RWORK(LYH), 1, Y, 1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NPE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      IWORK(23) = NLI
      IWORK(24) = NPS
      IWORK(25) = NCFL
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C if there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH, and T is set to TN.  The optional outputs
C are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'VODPK--  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWV (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 1, 1, MXSTEP, 0, 1, TN, ZERO)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'VODPK--  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWV (MSG, 50, 202, 1, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'VODPK--  At T (=R1), too much accuracy requested  '
      CALL XERRWV (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 1, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'VODPK--  At T(=R1) and step size H(=R2), the error'
      CALL XERRWV (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      test failed repeatedly or with abs(H) = HMIN'
      CALL XERRWV (MSG, 50, 204, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with abs(H) = HMIN. ----
 540  MSG = 'VODPK--  At T (=R1) and step size H (=R2), the    '
      CALL XERRWV (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      or with abs(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 560
C KFLAG = -3.  Unrecoverable error from JAC. ---------------------------
 550  MSG = 'VODPK--  at T (=R1) an unrecoverable error return '
      CALL XERRWV(MSG, 50, 206, 0, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      was made from subroutine JAC      '
      CALL XERRWV(MSG, 40, 206, 0, 0, 0, 0, 1, TN, ZERO)
      ISTATE = -7
      GO TO 580
C KFLAG = -4.  Unrecoverable error from PSOL. --------------------------
 552  MSG = 'VODPK--  at T (=R1) an unrecoverable error return '
      CALL XERRWV(MSG, 50, 207, 0, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      was made from subroutine PSOL     '
      CALL XERRWV(MSG, 40, 207, 0, 0, 0, 0, 1, TN, ZERO)
      ISTATE = -7
      GO TO 580
C KFLAG = -5.  F returned IFAIL .ne. 0 repeatedly in corrector loop. ---
 555  MSG = 'VODPK--  at T (=R1) and step size H (=R2), F      '
      CALL XERRWV(MSG, 50, 208, 0, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      failed repeatedly in the corrector loop.    '
      CALL XERRWV(MSG, 50, 208, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -8
      GO TO 580
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = ZERO
      IMXER = 1
      DO 570 I = 1, N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  CONTINUE
      CALL SCOPY (N, RWORK(LYH), 1, Y, 1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NPE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      IWORK(23) = NLI
      IWORK(24) = NPS
      IWORK(25) = NCFL
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C Call the error message routine and then return.
C-----------------------------------------------------------------------
 601  MSG = 'VODPK--  ISTATE (=I1) illegal '
      CALL XERRWV (MSG, 30, 1, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'VODPK--  ITASK (=I1) illegal  '
      CALL XERRWV (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
      GO TO 700
 603  MSG='VODPK--   ISTATE (=I1) .gt. 1 but VODPK not initialized     '
      CALL XERRWV (MSG, 60, 3, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      GO TO 700
 604  MSG = 'VODPK--  NEQ (=I1) .lt. 1     '
      CALL XERRWV (MSG, 30, 4, 1, 1, NEQ, 0, 0, ZERO, ZERO)
      GO TO 700
 605  MSG = 'VODPK--  ISTATE = 3 and NEQ increased (I1 to I2)  '
      CALL XERRWV (MSG, 50, 5, 1, 2, N, NEQ, 0, ZERO, ZERO)
      GO TO 700
 606  MSG = 'VODPK--  ITOL (=I1) illegal   '
      CALL XERRWV (MSG, 30, 6, 1, 1, ITOL, 0, 0, ZERO, ZERO)
      GO TO 700
 607  MSG = 'VODPK--  IOPT (=I1) illegal   '
      CALL XERRWV (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
      GO TO 700
 608  MSG = 'VODPK--  MF (=I1) illegal     '
      CALL XERRWV (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
      GO TO 700
 611  MSG = 'VODPK--  MAXORD (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 11, 1, 1, MAXORD, 0, 0, ZERO, ZERO)
      GO TO 700
 612  MSG = 'VODPK--  MXSTEP (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 12, 1, 1, MXSTEP, 0, 0, ZERO, ZERO)
      GO TO 700
 613  MSG = 'VODPK--  MXHNIL (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 13, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
      GO TO 700
 614  MSG = 'VODPK--  TOUT (=R1) behind T (=R2)      '
      CALL XERRWV (MSG, 40, 14, 1, 0, 0, 0, 2, TOUT, T)
      MSG = '      integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
      GO TO 700
 615  MSG = 'VODPK--  HMAX (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
      GO TO 700
 616  MSG = 'VODPK--  HMIN (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 16, 1, 0, 0, 0, 1, HMIN, ZERO)
      GO TO 700
 617  CONTINUE
      MSG='VODPK--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 17, 1, 2, LENRW, LRW, 0, ZERO, ZERO)
      GO TO 700
 618  CONTINUE
      MSG='VODPK--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWV (MSG, 60, 18, 1, 2, LENIW, LIW, 0, ZERO, ZERO)
      GO TO 700
 619  MSG = 'VODPK--  RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWV (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
      GO TO 700
 620  MSG = 'VODPK--  ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWV (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'VODPK--  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWV (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
      GO TO 700
 622  CONTINUE
      MSG='VODPK--  TOUT (=R1) too close to T(=R2) to start integration'
      CALL XERRWV (MSG, 60, 22, 1, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  CONTINUE
      MSG='VODPK--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWV (MSG, 60, 23, 1, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  CONTINUE
      MSG='VODPK--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWV (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  CONTINUE
      MSG='VODPK--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWV (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'VODPK--  At start of problem, too much accuracy   '
      CALL XERRWV (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'VODPK--  Trouble from VINDY. ITASK = I1, TOUT = R1'
      CALL XERRWV (MSG, 50, 27, 1, 1, ITASK, 0, 1, TOUT, ZERO)
C
 700  CONTINUE
      ISTATE = -3
      RETURN
C
 800  MSG = 'VODPK--  Run aborted.. apparent infinite loop     '
ccc      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)
      call xerrab('*** PROBABLY YOU NEED TO SET ts = 0. AND istate = 1')
      RETURN
C----------------------- End of Subroutine VODPK -----------------------
      END
      SUBROUTINE VNLSK (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM,
     1                 F, JAC, PSOL, NFLAG, RPAR, IPAR)
C
      EXTERNAL F, JAC, PSOL
      real Y, YH, VSAV, SAVF, EWT, ACOR, WM, RPAR
      INTEGER IWM, LDYH, NFLAG, IPAR
      DIMENSION Y(*), YH(LDYH, *),  SAVF(*), VSAV(*), EWT(*), ACOR(*),
     1          IWM(*), WM(*), RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C Call sequence input -- YH, LDYH, F, JAC, EWT, ACOR, IWM, WM,
C                        NFLAG, RPAR, IPAR
C Call sequence output -- Y, YH, VSAV, SAVF, ACOR, IWM, WM, NFLAG
C COMMON block variables accessed..
C        /VOD001/  ACNRM, CRATE, DRC, H, ICF, IPUP, JCUR, JSTART,
C                  METH, MITER, N, NSLP, RC, RL1, TN, TQ
C        /VOD002/  NFE, NNI, NPE, NST
C        /VPK001/  JACFLG, LOCIWP, LOCWP, MNEWT
C Subroutines called.. F, JAC, PSOL, SAXPY, SCOPY, SSCAL, VSOLPK
C Function subroutines called.. VNORML
C-----------------------------------------------------------------------
C Subroutine VNLSK is a nonlinear system solver, which uses either
C functional iteration (MITER = 0), or a combination of an inexact
C Newton method and preconditioned Krylov iteration (MITER .gt. 0)
C to solve the implicit system for the corrector y vector.
C It calls Subroutine JAC (user-supplied) for preprocessing the
C preconditioner, and Subroutine VSOLPK for the Krylov iteration.
C
C In addition to variables described elsewhere, communication with
C VNLSK uses the following variables..
C
C Y          = The dependent variable, a vector of length N, input.
C YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
C              and output.  On input, it contains predicted values.
C LDYH       = A constant .ge. N, the first dimension of YH, input.
C VSAV       = A work array of length N.
C SAVF       = A work array of length N.
C EWT        = An error weight vector of length N, input.
C ACOR       = A work array of length N, used for the accumulated
C              corrections to the predicted y vector.
C WM,IWM     = Real and integer work arrays associated with matrix
C              operations in Newton iteration (MITER .ne. 0).
C F          = Dummy name for user-supplied routine for f.
C JAC        = Dummy name for user-supplied routine for Jacobian data
C              and associated preconditioner matrix.
C PSOL       = Dummy name for user-supplied subroutine to solve
C              preconditioner linear system.
C NFLAG      = Input/output flag, with values and meanings as follows..
C              INPUT
C                  0 first call for this time step.
C                 -1 convergence failure in previous call to VNLSK.
C                 -2 error test failure in VSTEP.
C              OUTPUT
C                  0 successful completion of nonlinear solver.
C                 -1 convergence failure or failure in JAC.
C                 -2 unrecoverable error in matrix preprocessing
C                    (cannot occur here).
C                 -3 unrecoverable error in PSOL.
C                 -4 IFAIL .ne. 0 from F repeatedly.
C RPAR, IPAR = Dummy names for user's real and integer work arrays.
C
C IPUP       = Own variable flag with values and meanings as follows..
C              0,            do not update preconditioner.
C              MITER .ne. 0, update the preconditioner, because it is
C                            the initial step, user input changed,
C                            there was an error test failure, or an
C                            update is indicated by a change in the
C                            scalar RC or step counter NST.
C
C For more details, see comments in driver subroutine.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block VOD001 --------------------
C
      real ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block VOD002 --------------------
C
      real HU
      INTEGER NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
C
C Type declarations for labeled COMMON block VPK001 --------------------
C
      real DELT, SQRTN, RSQRTN
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LVSAV, KMP, MAXL, MNEWT,
     1      NLI, NPS, NCFL
C
      COMMON /VOD001/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /VOD002/ HU, NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
      COMMON /VPK001/ DELT, SQRTN, RSQRTN, JPRE, JACFLG, LOCIWP,
     1                LOCWP, LVSAV, KMP, MAXL, MNEWT, NLI, NPS, NCFL
c
c     Special scaling factor for Newton iteration error TDR 11/19/93
      real EFACN
      COMMON /TDR1/ EFACN
c
C
C Type declarations for local variables --------------------------------
C
      real CCMAX, CRDOWN, CSCALE, DEL, DCON, DELP, HRL1,RDIV
      real ONE, TWO, ZERO
      INTEGER I, IERPJ, IERSL, M, MAXCOR, MSBP, IFAIL
C
C Type declaration for function subroutines called ---------------------
C
      real VNORML
C-----------------------------------------------------------------------
C The following Fortran-77 declarations are to cause the values of the
C listed (local) variables to be saved between calls to VODPK.
      SAVE CCMAX, CRDOWN, MAXCOR, MSBP, RDIV
      SAVE ONE, TWO, ZERO
C-----------------------------------------------------------------------
      DATA CCMAX /0.3E0/, CRDOWN /0.3E0/, MAXCOR /3/, MSBP /20/,
     1     RDIV   /2.0E0/
      DATA ONE /1.0E0/, TWO /2.0E0/, ZERO /0.0E0/
C-----------------------------------------------------------------------
C Up to MAXCOR corrector iterations are taken.  A convergence test is
C made on the RMS norm of each correction, weighted by the error
C weight vector EWT.  The sum of the corrections is accumulated in the
C vector ACOR(*).  The YH array is not altered in the corrector loop.
C-----------------------------------------------------------------------
      IF (JSTART .EQ. 0) NSLP = 0
      IF (NFLAG .EQ. 0) ICF = 0
      IF (NFLAG .EQ. -2) IPUP = MITER
      IF ( (JSTART .EQ. 0) .OR. (JSTART .EQ. -1) ) IPUP = MITER
      IF (JACFLG .EQ. 0) THEN
        IPUP = 0
        CRATE = ONE
        GO TO 220
      ENDIF
      DRC = ABS(RC-ONE)
      IF (DRC .GT. CCMAX .OR. NST .GE. NSLP+MSBP) IPUP = MITER
 220  M = 0
      HRL1 = H*RL1
      DELP = ZERO
      MNEWT = 0
      CALL SCOPY (N, YH(1,1), 1, Y, 1 )
      IFAIL = 0
      CALL F (N, TN, Y, SAVF, RPAR, IPAR, IFAIL)
      NFE = NFE + 1
      IF (IFAIL .NE. 0) GO TO 420
      IF (IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
C If indicated, the preconditioner matrix is reevaluated and
C preprocessed before starting the corrector iteration.  IPUP is set
C to 0 as an indicator that this has been done.
C-----------------------------------------------------------------------
      JCUR = 1
      IERPJ = 0
      CALL JAC (F, N, TN, Y, YH, EWT, SAVF, ACOR, HRL1,
     1   WM(LOCWP), IWM(LOCIWP), IERPJ, RPAR, IPAR)
      NPE = NPE + 1
      IPUP = 0
      RC = ONE
      DRC = ZERO
      CRATE = ONE
      NSLP = NST
      IF (IERPJ .NE. 0) GO TO 420
 250  DO 260 I = 1, N
 260    ACOR(I) = ZERO
 270  IF (MITER .NE. 0) GO TO 350
C-----------------------------------------------------------------------
C In the case of functional iteration, update Y directly from
C the result of the last function evaluation.
C-----------------------------------------------------------------------
      DO 290 I = 1, N
        SAVF(I) = RL1*(H*SAVF(I) - YH(I,2))
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = VNORML (N, Y, EWT)
      DO 300 I = 1, N
 300    Y(I) = YH(I,1) + SAVF(I)
      CALL SCOPY (N, SAVF, 1, ACOR, 1)
      GO TO 400
C-----------------------------------------------------------------------
C In the case of the Newton method, compute the corrector error,
C and solve the linear system with that as right-hand side and
C A as coefficient matrix.  The correction is scaled by the factor
C 2/(1+RC) to account for changes in H*RL1 since the last JAC call.
C-----------------------------------------------------------------------
 350  DO 360 I = 1, N
 360    VSAV(I) = HRL1*SAVF(I) - (RL1*YH(I,2) + ACOR(I))
      CALL VSOLPK (Y, SAVF, VSAV, EWT, WM, IWM, F, PSOL, IERSL,
     1             RPAR, IPAR)
      NNI = NNI + 1
      IF (METH .EQ. 2 .AND. JACFLG .EQ. 1 
     1    .AND. MITER .EQ. 9 .AND. RC .NE. ONE) THEN
        CSCALE = TWO/(ONE + RC)
        CALL SSCAL (N, CSCALE, VSAV, 1)
      ENDIF
      IF (IERSL .LT. 0) GO TO 440
      IF (IERSL .GT. 0) GO TO 410
      DEL = VNORML (N, VSAV, EWT)
      CALL SAXPY (N, ONE, VSAV, 1, ACOR, 1)
      DO 380 I = 1, N
 380    Y(I) = YH(I,1) + ACOR(I)
C-----------------------------------------------------------------------
C Test for convergence.  If M.gt.0, an estimate of the convergence
C rate constant is stored in CRATE, and this is used in the test.
C-----------------------------------------------------------------------
 400  IF (M .NE. 0) CRATE = max(CRDOWN*CRATE,DEL/DELP)
      DCON = DEL*min(ONE,CRATE)/TQ(4)
      IF (DCON .LE. ONE*EFACN) GO TO 450
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. RDIV*DELP) GO TO 410
      MNEWT = M
      DELP = DEL
      IFAIL = 0
      CALL F (N, TN, Y, SAVF, RPAR, IPAR, IFAIL)
      NFE = NFE + 1
      IF (IFAIL .EQ. 0) GO TO 270
C
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1 .OR. JACFLG .EQ. 0) GO TO 420
      ICF = 1
      IPUP = MITER
      GO TO 220
C
 420  CONTINUE
      ICF = 2
      NFLAG = -1
      RETURN
 440  CONTINUE
      NFLAG = -3
      RETURN
C Return for successful step. ------------------------------------------
 450  NFLAG = 0
      JCUR = 0
      ICF = 0
      IF (M .EQ. 0) ACNRM = DEL
      IF (M .GT. 0) ACNRM = VNORML (N, ACOR, EWT)
      RETURN
C----------------------- End of Subroutine VNLSK -----------------------
      END
      SUBROUTINE VSOLPK (Y, SAVF, X, EWT, WM, IWM, F, PSOL, IERSL,
     1                   RPAR, IPAR)
      EXTERNAL F, PSOL
      real Y, SAVF, X, EWT, WM, RPAR
      INTEGER IWM, IERSL, IPAR
      DIMENSION Y(*), SAVF(*), X(*), EWT(*), WM(*), IWM(*),
     1   RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C Call sequence input -- Y, SAVF, X, EWT, F, PSOL, RPAR, IPAR
C Call sequence output -- Y, SAVF, X, WM, IWM, IERSL
C COMMON block variables accessed..
C        /VOD001/  H, RL1, TQ, TN, MITER, N
C        /VPK001/  DELT, SQRTN, RSQRTN, JPRE, LOCIWP, LOCWP,
C                  KMP, MAXL, MNEWT, NLI, NPS, NCFL
C Subroutines called.. F, PSOL, SCOPY, SSCAL, VSPIG, VUSOL
C-----------------------------------------------------------------------
C This routine interfaces with  VSPIG  or  VUSOL  for the solution of
C the linear system arising from a Newton iteration (MITER .ne. 0).
C
C In addition to variables described elsewhere, communication with
C VSOLPK uses the following variables..
C WM    = real work space containing data for the algorithm
C         (Krylov basis vectors, Hessenberg matrix, etc.)
C IWM   = integer work space containing data for the algorithm
C X     = the right-hand side vector on input, and the solution vector
C         on output, of length N.
C IERSL = output flag (in COMMON)..
C         IERSL =  0 means no trouble occurred.
C         IERSL =  1 means the iterative method failed to converge.
C                    If the preconditioner is out of date, the step
C                    is repeated with a new preconditioner.  Otherwise,
C                    the stepsize is reduced (forcing a new evalua-
C                    tion of the preconditioner) and the step is
C                    repeated.
C         IERSL = -1 means there was a nonrecoverable error in the
C                    iterative solver.  The stepsize is reduced in
C                    VSTEP and the step is repeated.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block VOD001 --------------------
C
      real ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block VPK001 --------------------
C
      real DELT, SQRTN, RSQRTN
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LVSAV, KMP, MAXL, MNEWT,
     1      NLI, NPS, NCFL
C
C Type declarations for local variables --------------------------------
C
      real DELTA, HRL1
      INTEGER IFLAG, LB, LDL, LGMR, LHES, LQ, LV, LWK, MAXLP1, NPSL
C-----------------------------------------------------------------------
      COMMON /VOD001/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /VPK001/ DELT, SQRTN, RSQRTN, JPRE, JACFLG, LOCIWP,
     1                LOCWP, LVSAV, KMP, MAXL, MNEWT, NLI, NPS, NCFL
C-----------------------------------------------------------------------
C
      IERSL = 0
      HRL1 = H*RL1
      DELTA = DELT*TQ(4)
      IF (MITER .EQ. 1) THEN
C-----------------------------------------------------------------------
C Use the VSPIG algorithm to solve the linear system A*x = -f.
C-----------------------------------------------------------------------
        MAXLP1 = MAXL + 1
        LV = 1
        LB = LV + N*MAXL
        LHES = LB + N + 1
        LQ = LHES + MAXL*MAXLP1
        LWK = LQ + 2*MAXL
        LDL = LWK + min(1,MAXL-KMP)*N
        CALL SCOPY (N, X, 1, WM(LB), 1)
        CALL SSCAL (N, RSQRTN, EWT, 1)
        CALL VSPIG (TN, Y, SAVF, WM(LB), EWT, N, MAXL, MAXLP1, KMP,
     1     DELTA, HRL1, JPRE, MNEWT, F, PSOL, NPSL, X, WM(LV), WM(LHES),
     2     WM(LQ), LGMR, WM(LOCWP), IWM(LOCIWP), WM(LWK), WM(LDL),
     3     RPAR, IPAR, IFLAG)
        NLI = NLI + LGMR
        NPS = NPS + NPSL
        CALL SSCAL (N, SQRTN, EWT, 1)
        IF (IFLAG .NE. 0) NCFL = NCFL + 1
        IF (IFLAG .GE. 2) IERSL = 1
        IF (IFLAG .LT. 0) IERSL = -1
        RETURN
      ELSE IF (MITER .EQ. 9) THEN
C-----------------------------------------------------------------------
C Use VUSOL, which interfaces to PSOL, to solve the linear system
C (No Krylov iteration).
C-----------------------------------------------------------------------
        LB = 1
        LWK = LB + N
        CALL SCOPY (N, X, 1, WM(LB), 1)
        CALL VUSOL (N, TN, Y, SAVF, WM(LB), EWT, DELTA, HRL1, JPRE,
     1     MNEWT, PSOL, NPSL, X, WM(LOCWP), IWM(LOCIWP), WM(LWK),
     2     RPAR, IPAR, IFLAG)
        NPS = NPS + NPSL
        IF (IFLAG .NE. 0) NCFL = NCFL + 1
        IF (IFLAG .EQ. 3) IERSL = 1
        IF (IFLAG .LT. 0) IERSL = -1
        RETURN
       ENDIF
C----------------------- End of Subroutine VSOLPK ----------------------
      END
      SUBROUTINE VSPIG (TN, Y, SAVF, B, WGHT, N, MAXL, MAXLP1,
     1  KMP, DELTA, HB0, JPRE, MNEWT, F, PSOL, NPSL, X, V, HES, Q,
     2  LGMR, WP, IWP, WK, DL, RPAR, IPAR, IFLAG)
      EXTERNAL F, PSOL
      real TN, Y, SAVF, B, WGHT, DELTA, HB0, X, V, HES,
     1   Q, WP, WK, DL, RPAR
      INTEGER N, MAXL, MAXLP1, KMP, JPRE, MNEWT, NPSL, LGMR, IWP,
     1   IFLAG, IPAR
      DIMENSION Y(*), SAVF(*), B(*), WGHT(*), X(*), V(N,*),
     1   HES(MAXLP1,*), Q(*), WP(*), IWP(*), WK(*), DL(*),
     2   RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C Call sequence input --  TN, Y, SAVF, B, WGHT, N, MAXL, MAXLP1, DELTA,
C                         HB0, JPRE, MNEWT, F, PSOL, RPAR, IPAR
C Call sequence output -- B, KMP, DELTA, NPSL, X, V, HES, Q, LGMR, WP,
C                         IWP, WK, DL, RPAR, IPAR, IFLAG
C COMMON block variables accessed.. None
C Subroutines called.. F, SVRORTHOG, PSOL, SAXPY, SCOPY, SHELS, SHEQR,
C                      SSCAL, VATV
C Function subroutines called.. SNRM2
C-----------------------------------------------------------------------
C This routine solves the linear system  A * x = b using SPIGMR,
C a scaled preconditioned incomplete version of the generalized
C minimum residual method GMRES.
C An initial guess of x = 0 is assumed.
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            B = the right hand side of the system A*x = b.
C                B is also used as work space when computing
C                the final approximation.
C                (B is the same as V(*,MAXL+1) in the call to VSPIG.)
C
C         WGHT = the vector of length N containing the nonzero
C                elements of the diagonal scaling matrix.
C
C            N = the order of the matrix A, and the lengths
C                of the vectors WGHT, B and X.
C
C         MAXL = the maximum allowable order of the matrix H.
C
C       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES .
C
C          KMP = the number of previous vectors the new vector VNEW
C                must be made orthogonal to.  KMP .le. MAXL.
C
C        DELTA = tolerance on residuals  b - A*x  in weighted RMS norm.
C
C          HB0 = current value of (step size h) * (coefficient beta0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array used by routine VATV and PSOL.
C
C           DL = real work array used for calculation of the residual
C                norm rho when the method is incomplete (KMP.lt.MAXL).
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         LGMR = the number of iterations performed and
C                the current order of the upper Hessenberg
C                matrix HES .
C
C         NPSL = the number of calls to PSOL.
C
C         V    = the N by (LGMR+1) array containing the LGMR
C                orthogonal vectors V(*,1) to V(*,LGMR).
C
C         HES  = the upper triangular factor of the QR decomposition
C                of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C                entries are the scaled inner-products of A*V(*,i)
C                and V(*,k).
C
C         Q    = real array of length 2*MAXL containing the components
C                of the Givens rotations used in the QR decomposition
C                of HES.  It is loaded in SHEQR and used in SHELS.
C
C        IFLAG = integer error flag..
C                0 means convergence in LGMR iterations, LGMR.le.MAXL.
C                1 means the convergence test did not pass in MAXL
C                  iterations, but the residual norm is .lt. 1,
C                  or .lt. norm(b) if MNEWT = 0, and so x is computed.
C                2 means the convergence test did not pass in MAXL
C                  iterations, residual .gt. 1, and x is undefined.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C                4 means ATV failed with repeated failures of F.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
C
C Type declarations for local variables --------------------------------
C
      real BNRM, BNRM0, C, DLNRM, PROD, RHO, S, SNORMW, TEM
      INTEGER I, IER, INFO, IP1, I2, J, K, LL, LLP1
C
C Type declarations for function subroutines called --------------------
C
      real SNRM2
C
      IFLAG = 0
      LGMR = 0
      NPSL = 0
C-----------------------------------------------------------------------
C The initial residual is the vector b.  Apply scaling to b, and test
C for an immediate return with x = 0 or x = b.
C-----------------------------------------------------------------------
      DO 10 I = 1, N
 10     V(I,1) = B(I)*WGHT(I)
      BNRM0 = SNRM2 (N, V, 1)
      BNRM = BNRM0
      IF (BNRM0 .GT. DELTA) GO TO 30
      IF (MNEWT .GT. 0) GO TO 20
      CALL SCOPY (N, B, 1, X, 1)
      RETURN
 20   DO 25 I = 1, N
 25     X(I) = 0.0E0
      RETURN
 30   CONTINUE
C Apply inverse of left preconditioner to vector b. --------------------
      IER = 0
      IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GO TO 40
      CALL PSOL (N, TN, Y, SAVF, WK, HB0, WP, IWP, B, 1,
     1           IER, RPAR, IPAR)
      NPSL = 1
      IF (IER .NE. 0) GO TO 300
 40   IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GO TO 55
C Calculate norm of scaled vector V(*, 1) and normalize it. ------------
      DO 50 I = 1, N
 50     V(I,1) = B(I)*WGHT(I)
      BNRM = SNRM2 (N, V, 1)
      DELTA = DELTA*(BNRM/BNRM0)
 55   TEM = 1.0E0/BNRM
      CALL SSCAL (N, TEM, V(1,1), 1)
C Zero out the HES array. ----------------------------------------------
      DO 65 J = 1, MAXL
        DO 60 I = 1, MAXLP1
 60        HES (I,J) = 0.0E0
 65     CONTINUE
C-----------------------------------------------------------------------
C Main loop to compute the vectors V(*,2) to V(*,MAXL).
C The running product PROD is needed for the convergence test.
C-----------------------------------------------------------------------
      PROD = 1.0E0
      DO 90 LL = 1, MAXL
        LGMR = LL
C-----------------------------------------------------------------------
C Call routine VATV to compute VNEW = Abar*v(ll), where Abar is
C the matrix A with scaling and inverse preconditioner factors applied.
C Call routine SVRORTHOG to orthogonalize the new vector VNEW = V(*,LL+1).
C Call routine SHEQR to update the factors of HES .
C-----------------------------------------------------------------------
        CALL VATV (Y, SAVF, V(1,LL), WGHT, X, F, PSOL, RPAR, IPAR,
     1        V(1,LL+1), WK, WP, IWP, HB0, JPRE, IER, NPSL)
        IF (IER .NE. 0) GO TO 310
        CALL SVRORTHOG (V(1, LL+1), V,  HES , N, LL, MAXLP1, KMP,
     1     SNORMW)
        HES (LL+1, LL) = SNORMW
        CALL SHEQR (HES , MAXLP1, LL, Q, INFO, LL)
        IF (INFO .EQ. LL) GO TO 120
C-----------------------------------------------------------------------
C Update RHO, the estimate of the norm of the residual b - A*xl.
C If KMP .lt. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
C necessarily orthogonal for LL .gt. KMP.  The vector DL must then
C be computed, and its norm used in the calculation of RHO
C-----------------------------------------------------------------------
        PROD = PROD*Q(2*LL)
        RHO = ABS(PROD*BNRM)
        IF ((LL.GT.KMP) .AND. (KMP.LT.MAXL)) THEN
          IF (LL .EQ. KMP+1) THEN
            CALL SCOPY (N, V(1,1), 1, DL, 1)
            DO 75 I = 1, KMP
              IP1 = I + 1
              I2 = I*2
              S = Q(I2)
              C = Q(I2-1)
              DO 70 K = 1, N
 70             DL(K) = S*DL(K) + C*V(K,IP1)
 75           CONTINUE
            ENDIF
          S = Q(2*LL)
          C = Q(2*LL-1)/SNORMW
          LLP1 = LL + 1
          DO 80 K = 1, N
 80         DL(K) = S*DL(K) + C*V(K,LLP1)
          DLNRM = SNRM2 (N, DL, 1)
          RHO = RHO*DLNRM
          ENDIF
C-----------------------------------------------------------------------
C Test for convergence.  If passed, compute approximation xl.
C If failed and LL .lt. MAXL, then continue iterating.
C-----------------------------------------------------------------------
        IF (RHO .LE. DELTA) GO TO 200
        IF (LL .EQ. MAXL) GO TO 100
C-----------------------------------------------------------------------
C Rescale so that the norm of V(1,LL+1) is one.
C-----------------------------------------------------------------------
        TEM = 1.0E0/SNORMW
        CALL SSCAL (N, TEM, V(1,LL+1), 1)
 90     CONTINUE
 100  CONTINUE
      IF (RHO .LE. 1.0E0) GO TO 150
      IF (RHO .LE. BNRM .AND. MNEWT .EQ. 0) GO TO 150
 120  CONTINUE
      IFLAG = 2
      RETURN
 150  IFLAG = 1
C-----------------------------------------------------------------------
C Compute the approximation xl to the solution.
C Since the vector X was used as work space, and the initial guess
C of the Newton correction is zero, x must be reset to zero.
C-----------------------------------------------------------------------
 200  CONTINUE
      LL = LGMR
      LLP1 = LL + 1
      DO 210 K = 1, LLP1
 210    B(K) = 0.0E0
      B(1) = BNRM
      CALL SHELS (HES , MAXLP1, LL, Q, B)
      DO 220 K = 1, N
 220    X(K) = 0.0E0
      DO 230 I = 1, LL
        CALL SAXPY (N, B(I), V(1,I), 1, X, 1)
 230    CONTINUE
      DO 240 I = 1, N
 240    X(I) = X(I)/WGHT(I)
      IF (JPRE .LE. 1) RETURN
      CALL PSOL (N, TN, Y, SAVF, WK, HB0, WP, IWP, X, 2,
     1          IER, RPAR, IPAR)
      NPSL = NPSL + 1
      IF (IER .NE. 0) GO TO 300
      RETURN
C Error returns forced by routine PSOL. --------------------------------
 300  CONTINUE
      IF (IER .LT. 0) IFLAG = -1
      IF (IER .GT. 0) IFLAG = 3
      RETURN
C Error returns forced by ATV (by failure of PSOL or F). ---------------
 310  CONTINUE
      IF (IER .EQ. -1) IFLAG = -1
      IF (IER .EQ. 1) IFLAG = 3
      IF (IER .EQ. 2) IFLAG = 4
      RETURN
C----------------------- End of Subroutine VSPIG -----------------------
      END
      SUBROUTINE VATV (Y, SAVF, V, WGHT, FTEM, F, PSOL, RPAR, IPAR,
     1                Z, VTEM, WP, IWP, HB0, JPRE, IER, NPSL)
      EXTERNAL F, PSOL
      real Y, SAVF, V, WGHT, FTEM, RPAR, Z, VTEM, WP, HB0
      INTEGER IPAR, IWP, JPRE, IER, NPSL
      DIMENSION Y(*), SAVF(*), V(*), WGHT(*), FTEM(*), Z(*),
     1   VTEM(*), WP(*), IWP(*), RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C Call sequence input -- Y, SAVF, V, WGHT, F, PSOL, RPAR, IPAR,
C                        WP, IWP, HB0, NPSL
C Call sequence output --Z, IER, NPSL
C COMMON block variables accessed..
C        /VOD001/  TN, N
C        /VOD002/  NFE
C Subroutines called.. F, PSOL, SCOPY
C Function subroutines called.. SNRM2
C-----------------------------------------------------------------------
C This routine computes the product
C
C   (D-inverse)*(P1-inverse)*(I - hb0*df/dy)*(P2-inverse)*(D*v),
C
C where D is a diagonal scaling matrix, and P1 and P2 are the
C left and right preconditioning matrices, respectively.
C v is assumed to have WRMS norm equal to 1.
C The product is stored in Z.  This is computed by a
C difference quotient, a call to F, and two calls to PSOL.
C If the F call fails, up to 4 more tries are made with a smaller 
C value for the increment.
C-----------------------------------------------------------------------
C
C      On entry
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            V = real array of length N (can be the same array as Z).
C
C         WGHT = array of length N containing scale factors.
C                1/WGHT(i) are the diagonal elements of the matrix D.
C
C         FTEM = work array of length N.
C
C         VTEM = work array of length N used to store the
C                unscaled version of v.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C          HB0 = current value of (step size h) * (coefficient beta0).
C
C         JPRE = preconditioner type flag.
C
C
C      On return
C
C            Z = array of length N containing desired scaled
C                matrix-vector product.
C
C          IER = error flag.
C                IER =  1 means PSOL failed recoverably.
C                IER = -1 means PSOL failed unrecoverably.
C                IER =  2 means calls to F failed repeatedly.
C
C         NPSL = the number of calls to PSOL.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block VOD001 --------------------
C
      real ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block VOD002 --------------------
C
      real HU
      INTEGER NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
C
C Type declarations for local variables --------------------------------
C
      real FAC, SIG, RSIG
      INTEGER I, IERP, NFAIL, IFAIL

C
C Type declaration for function subroutines called ---------------------
C
      real SNRM2
C
      COMMON /VOD001/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /VOD002/ HU, NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
      real scaleps
      common /comvatv/ scaleps
C-----------------------------------------------------------------------
      NFAIL = 0
C Set vtem = D * v. ----------------------------------------------------
      DO 10 I = 1, N
 10     VTEM(I) = V(I)/WGHT(I)
      IF (JPRE .GE. 2) GO TO 30
C
C JPRE = 0 or 1.  Save y in Z and increment Y by VTEM. -----------------
      CALL SCOPY (N, Y, 1, Z, 1)
      DO 20 I = 1, N
 20     Y(I) = Z(I) + VTEM(I)*scaleps
      SIG = scaleps
      RSIG = 1.0E0/scaleps
      FAC = HB0*RSIG
      GO TO 60
C
C JPRE = 2 or 3.  Apply inverse of right preconditioner to VTEM. -------
 30   CONTINUE
      CALL PSOL (N, TN, Y, SAVF, FTEM, HB0, WP, IWP, VTEM, 2,
     1          IERP, RPAR, IPAR)
      NPSL = NPSL + 1
      IF (IERP .NE. 0) GO TO 100
C Calculate l-2 norm of (D-inverse) * VTEM. ----------------------------
      DO 40 I = 1, N
 40     Z(I) = VTEM(I)*WGHT(I)
      RSIG = SNRM2 (N, Z, 1)/scaleps
      SIG = 1.0E0/RSIG
C Save y in Z and increment Y by VTEM/norm. ----------------------------
      CALL SCOPY (N, Y, 1, Z, 1)
 45   DO 50 I = 1, N
 50     Y(I) = Z(I) + VTEM(I)*SIG
      FAC = HB0*RSIG
C
C For all JPRE, call F with incremented Y argument, and restore Y. -----
 60   CONTINUE
      IFAIL = 0
      CALL F (N, TN, Y, FTEM, RPAR, IPAR, IFAIL)
      NFE = NFE + 1
      IF (IFAIL .NE. 0) THEN
        NFAIL = NFAIL + 1
        IF (NFAIL .GT. 5) GO TO 110
        SIG = 0.25E0*SIG
        RSIG = 4.0E0*RSIG
        GO TO 45
        ENDIF
      CALL SCOPY (N, Z, 1, Y, 1)
C Set Z = (I - HB0*Jacobian) * VTEM, using difference quotient. --------
      DO 70 I = 1, N
 70     Z(I) = FTEM(I) - SAVF(I)
      DO 80 I = 1, N
 80     Z(I) = VTEM(I) - FAC*Z(I)
C Apply inverse of left preconditioner to Z, if nontrivial. ------------
      IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GO TO 85
      CALL PSOL (N, TN, Y, SAVF, FTEM, HB0, WP, IWP, Z, 1,
     1           IERP, RPAR, IPAR)
      NPSL = NPSL + 1
      IF (IERP .NE. 0) GO TO 100
 85   CONTINUE
C Apply D-inverse to Z and return. -------------------------------------
      DO 90 I = 1, N
 90     Z(I) = Z(I)*WGHT(I)
      IER = 0
      RETURN
C
 100  IER = 1
      IF (IERP .LT. 0) IER = -1
      RETURN
 110  IER = 2
      RETURN
C----------------------- End of Subroutine VATV ------------------------
      END
      SUBROUTINE VUSOL (N, TN, Y, SAVF, B, WGHT, DELTA, HB0, JPRE,
     1   MNEWT, PSOL, NPSL, X, WP, IWP, WK, RPAR, IPAR, IFLAG)
      EXTERNAL PSOL
      real TN, Y, SAVF, B, WGHT, DELTA, HB0, X, WP, WK, RPAR
      INTEGER N, JPRE, MNEWT, NPSL, IWP, IPAR, IFLAG
      DIMENSION Y(*), SAVF(*), B(*), WGHT(*), X(*),
     1   WP(*), IWP(*), WK(*), RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C This routine solves the linear system A * x = b using only
C calls to the user-supplied routine PSOL (no Krylov iteration).
C If the norm of the right-hand side vector b is smaller than DELTA,
C the vector x returned is x = b (if MNEWT = 0) or x = 0 otherwise.
C PSOL is called with an LR argument of 1 (if JPRE = 1 or 3),
C then 2 (if JPRE = 2 or 3).
C-----------------------------------------------------------------------
C
C      On entry
C
C          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
C
C           TN = current value of t.
C
C            Y = array containing current dependent variable vector.
C
C         SAVF = array containing current value of f(t,y).
C
C            B = the right hand side of the system A*x = b.
C
C         WGHT = the vector of length N containing the nonzero
C                elements of the diagonal scaling matrix.
C
C            N = the order of the matrix A, and the lengths
C                of the vectors WGHT, b and x.
C
C        DELTA = tolerance on residuals  b - A*x  in weighted RMS norm.
C
C          HB0 = current value of (step size h) * (coefficient beta0).
C
C         JPRE = preconditioner type flag.
C
C        MNEWT = Newton iteration counter (.ge. 0).
C
C           WK = real work array used by PSOL.
C
C           WP = real work array used by preconditioner PSOL.
C
C          IWP = integer work array used by preconditioner PSOL.
C
C      On return
C
C         X    = the final computed approximation to the solution
C                of the system A*x = b.
C
C         NPSL = the number of calls to PSOL.
C
C        IFLAG = integer error flag..
C                0 means no trouble occurred.
C                3 means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C               -1 means there was a nonrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
C
C Type declarations for local variables --------------------------------
C
      real BNRM
      INTEGER I, IER
C
C Type declaration for function subroutines called ---------------------
C
      real VNORML
C
      IFLAG = 0
      NPSL = 0
C-----------------------------------------------------------------------
C Test for an immediate return with x = 0 or x = b.
C-----------------------------------------------------------------------
      BNRM = VNORML (N, B, WGHT)
      IF (BNRM .GT. DELTA) GO TO 30
      IF (MNEWT .GT. 0) GO TO 10
      CALL SCOPY (N, B, 1, X, 1)
      RETURN
 10   DO 20 I = 1, N
 20     X(I) = 0.0E0
      RETURN
C Apply inverse of left preconditioner to vector b. --------------------
 30   IER = 0
      IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GO TO 40
      CALL PSOL (N, TN, Y, SAVF, WK, HB0, WP, IWP, B, 1,
     1           IER, RPAR,IPAR)
      NPSL = 1
      IF (IER .NE. 0) GO TO 100
C Apply inverse of right preconditioner to result, and copy to X. ------
 40   IF (JPRE .LE. 1) GO TO 50
      CALL PSOL (N, TN, Y, SAVF, WK, HB0, WP, IWP, B, 2,
     1           IER, RPAR, IPAR)
      NPSL = NPSL + 1
      IF (IER .NE. 0) GO TO 100
 50   CALL SCOPY (N, B, 1, X, 1)
      RETURN
C-----------------------------------------------------------------------
C This block handles error returns forced by routine PSOL.
C-----------------------------------------------------------------------
 100  CONTINUE
      IF (IER .LT. 0) IFLAG = -1
      IF (IER .GT. 0) IFLAG = 3
      RETURN
C----------------------- End of Subroutine VUSOL -----------------------
      END
      SUBROUTINE VKSRC (RSAV, ISAV, JOB)
      real RSAV
      INTEGER ISAV, JOB
      DIMENSION RSAV(*), ISAV(*)
C-----------------------------------------------------------------------
C Call sequence input -- RSAV, ISAV, JOB
C Call sequence output -- RSAV, ISAV
C COMMON block variables accessed -- All of /VOD001/, /VOD002/, /VODPK1/
C
C Subroutines/functions called by VKSRC.. None
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of the
C COMMON blocks VOD001, VOD002, VPK001, used internally by VODPK.
C
C RSAV = real array of length 52 or more.
C ISAV = integer array of length 52 or more.
C JOB  = flag indicating to save or restore the COMMON blocks..
C        JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV).
C        JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV).
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block VOD001 --------------------
C
      real RVOD1
      INTEGER IVOD1
C
C Type declarations for labeled COMMON block VOD002 --------------------
C
      real RVOD2
      INTEGER IVOD2
C
C Type declarations for labeled COMMON block VPK001 --------------------
C
      real RVPK1
      INTEGER IVPK1
C
C Type declarations for local variables --------------------------------
C
C
      INTEGER I, IOFF, LENIV1, LENIV2, LENRV1, LENRV2, LRVK1, LIVK1
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE LENRV1, LENIV1, LENRV2, LENIV2, LRVK1, LIVK1
C-----------------------------------------------------------------------
      COMMON /VOD001/ RVOD1(48), IVOD1(33)
      COMMON /VOD002/ RVOD2(1), IVOD2(8)
      COMMON /VPK001/ RVPK1(3), IVPK1(11)
      DATA LENRV1 /48/, LENIV1 /33/, LENRV2 /1/, LENIV2 /8/,
     1   LRVK1 /3/, LIVK1 /11/
C
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1, LENRV1
 10     RSAV(I) = RVOD1(I)
      DO 12 I = 1, LENRV2
 12     RSAV(LENRV1+I) = RVOD2(I)
      IOFF = LENRV1 + LENRV2
      DO 14 I = 1, LRVK1
 14     RSAV(IOFF+I) = RVPK1(I)
C
      DO 20 I = 1, LENIV1
 20     ISAV(I) = IVOD1(I)
      DO 22 I = 1, LENIV2
 22     ISAV(LENIV1+I) = IVOD2(I)
      IOFF = LENIV1 + LENIV2
      DO 24 I = 1, LIVK1
 24     ISAV(IOFF+I) = IVPK1(I)
C
      RETURN
C
 100  CONTINUE
      DO 110 I = 1, LENRV1
 110     RVOD1(I) = RSAV(I)
      DO 112 I = 1, LENRV2
 112     RVOD2(I) = RSAV(LENRV1+I)
      IOFF = LENRV1 + LENRV2
      DO 114 I = 1, LRVK1
 114    RVPK1(I) = RSAV(IOFF+I)
C
      DO 120 I = 1, LENIV1
 120     IVOD1(I) = ISAV(I)
      DO 122 I = 1, LENIV2
 122     IVOD2(I) = ISAV(LENIV1+I)
      IOFF = LENIV1 + LENIV2
      DO 124 I = 1, LIVK1
 124    IVPK1(I) = ISAV(IOFF+I)
C
      RETURN
C----------------------- End of Subroutine VKSRC -----------------------
      END
      SUBROUTINE VHIN (N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
     1   EWT, ITOL, ATOL, Y, TEMP, H0, NITER, IER)
      EXTERNAL F
      real T0, Y0, YDOT, RPAR, TOUT, UROUND, EWT, ATOL, Y,
     1   TEMP, H0
      INTEGER N, IPAR, ITOL, NITER, IER
      DIMENSION Y0(*), YDOT(*), EWT(*), ATOL(*), Y(*),
     1   TEMP(*), RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
C                        EWT, ITOL, ATOL, Y, TEMP
C Call sequence output -- H0, NITER, IER
C COMMON block variables accessed -- None
C
C Subroutines called by VHIN.. F
C Function routines called by VHIN.. VNORML
C-----------------------------------------------------------------------
C This routine computes the step size, H0, to be attempted on the
C first step, when the user has not supplied a value for this.
C
C First we check that TOUT - T0 differs significantly from zero.  Then
C an iteration is done to approximate the initial second derivative
C and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
C A bias factor of 1/2 is applied to the resulting h.
C The sign of H0 is inferred from the initial values of TOUT and T0.
C
C Communication with VHIN is done with the following variables..
C
C N      = Size of ODE system, input.
C T0     = Initial value of independent variable, input.
C Y0     = Vector of initial conditions, input.
C YDOT   = Vector of initial first derivatives, input.
C F      = Name of subroutine for right-hand side f(t,y), input.
C RPAR, IPAR = Dummy names for user's real and integer work arrays.
C TOUT   = First output value of independent variable
C UROUND = Machine unit roundoff
C EWT, ITOL, ATOL = Error weights and tolerance parameters
C                   as described in the driver routine, input.
C Y, TEMP = Work arrays of length N.
C H0     = Step size to be attempted, output.
C NITER  = Number of iterations (and of f evaluations) to compute H0,
C          output.
C IER    = The error flag, returned with the value
C          IER = 0  if no trouble occurred, or
C          IER = -1 if TOUT and T0 are considered too close to proceed.
C          IER = -2 if F returned IFAIL .ne. 0, after 5 attempts to
C                   avoid this condition.
C-----------------------------------------------------------------------
C
C Type declarations for local variables --------------------------------
C
      real AFI, ATOLI, DELYI, HALF, HG, HLB, HNEW, HRAT,
     1     HUB, HUN, PT1, PT2, T1, TDIST, TROUND, TWO, YDDNRM
      INTEGER I, ITER, IFAIL
C
C Type declaration for function subroutines called ---------------------
C
      real VNORML
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE HALF, HUN, PT1, TWO
      DATA HALF /0.5E0/, HUN /100.0E0/, PT1 /0.1E0/, PT2 /0.2E0/,
     1     TWO /2.0E0/
C
      NITER = 0
      TDIST = ABS(TOUT - T0)
      TROUND = UROUND*max(ABS(T0),ABS(TOUT))
      IF (TDIST .LT. TWO*TROUND) GO TO 100
C
C Set a lower bound on h based on the roundoff level in T0 and TOUT. ---
      HLB = HUN*TROUND
C Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. -
      HUB = PT1*TDIST
      ATOLI = ATOL(1)
      DO 10 I = 1, N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        DELYI = PT1*ABS(Y0(I)) + ATOLI
        AFI = ABS(YDOT(I))
        IF (AFI*HUB .GT. DELYI) HUB = DELYI/AFI
 10     CONTINUE
C
C Set initial guess for h as geometric mean of upper and lower bounds. -
      ITER = 0
      HG = SQRT(HLB*HUB)
C If the bounds have crossed, exit with the mean value. ----------------
      IF (HUB .LT. HLB) THEN
        H0 = HG
        GO TO 90
      ENDIF
C
C Looping point for iteration. -----------------------------------------
 50   CONTINUE
C Estimate the second derivative as a difference quotient in f. --------
      T1 = T0 + HG
      DO 60 I = 1, N
 60     Y(I) = Y0(I) + HG*YDOT(I)
      IFAIL = 0
      CALL F (N, T1, Y, TEMP, RPAR, IPAR, IFAIL)
      IF (IFAIL .NE. 0) THEN
        HNEW = PT2*HG
        GO TO 75
        ENDIF
      DO 70 I = 1, N
 70     TEMP(I) = (TEMP(I) - YDOT(I))/HG
      YDDNRM = VNORML (N, TEMP, EWT)
C Get the corresponding new value of h. --------------------------------
      IF (YDDNRM*HUB*HUB .GT. TWO) THEN
        HNEW = SQRT(TWO/YDDNRM)
      ELSE
        HNEW = SQRT(HG*HUB)
      ENDIF
 75   ITER = ITER + 1
C-----------------------------------------------------------------------
C Test the stopping conditions.
C Stop if the new and previous h values differ by a factor of .lt. 2.
C Stop if four iterations have been done.  Also, stop with previous h
C if HNEW/HG .gt. 2 after first iteration, as this probably means that
C the second derivative value is bad because of cancellation error.
C-----------------------------------------------------------------------
      IF (ITER .GE. 4) GO TO 80
      HRAT = HNEW/HG
      IF ( (HRAT .GT. HALF) .AND. (HRAT .LT. TWO) ) GO TO 80
      IF ( (ITER .GE. 2) .AND. (HNEW .GT. TWO*HG) ) THEN
        HNEW = HG
        GO TO 80
      ENDIF
      HG = HNEW
      GO TO 50
C
C Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ----
 80   IF (IFAIL .NE. 0) GO TO 110
      H0 = HNEW*HALF
      IF (H0 .LT. HLB) H0 = HLB
      IF (H0 .GT. HUB) H0 = HUB
 90   H0 = SIGN(H0, TOUT - T0)
      NITER = ITER
      IER = 0
      RETURN
C Error return for TOUT - T0 too small. --------------------------------
 100  IER = -1
      RETURN
C Error return for IFAIL .ne. 0 from F. --------------------------------
 110  IER = -2
      RETURN
C----------------------- End of Subroutine VHIN ------------------------
      END
      SUBROUTINE VINDY (T, K, YH, LDYH, DKY, IFLAG)
      real T, YH, DKY
      INTEGER K, LDYH, IFLAG
      DIMENSION YH(LDYH,*), DKY(*)
C-----------------------------------------------------------------------
C Call sequence input -- T, K, YH, LDYH
C Call sequence output -- DKY, IFLAG
C COMMON block variables accessed..
C     /VOD001/ --  H, TN, UROUND, L, N, NQ
C     /VOD002/ --  HU
C
C Subroutines called by VINDY.. SSCAL, XERRWV
C Function routines called by VINDY.. None
C-----------------------------------------------------------------------
C VINDY computes interpolated values of the K-th derivative of the
C dependent variable vector y, and stores it in DKY.  This routine
C is called within the package with K = 0 and T = TOUT, but may
C also be called by the user for any K up to the current order.
C (See detailed instructions in the usage documentation.)
C-----------------------------------------------------------------------
C The computed values in DKY are gotten by interpolation using the
C Nordsieck history array YH.  This array corresponds uniquely to a
C vector-valued polynomial of degree NQCUR or less, and DKY is set
C to the K-th derivative of this polynomial at T.
C The formula for DKY is..
C              q
C  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
C             j=K
C where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
C The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
C communicated by COMMON.  The above sum is done in reverse order.
C IFLAG is returned negative if either K or T is out of bounds.
C
C Discussion above and comments in driver explain all variables.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block VOD001 --------------------
C
      real ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block VOD002 --------------------
C
      real HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
C Type declarations for local variables --------------------------------
C
      real C, HUN, R, S, TFUZZ, TN1, TP, ZERO
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
      CHARACTER*80 MSG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE HUN, ZERO
C
      COMMON /VOD001/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /VOD002/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
      DATA HUN /100.0E0/, ZERO /0.0E0/
C
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TFUZZ = HUN*UROUND*(TN + HU)
      TP = TN - HU - TFUZZ
      TN1 = TN + TFUZZ
      IF ((T-TP)*(T-TN1) .GT. ZERO) GO TO 90
C
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1, NQ
 10     IC = IC*JJ
 15   C = float(IC)
      DO 20 I = 1, N
 20     DKY(I) = C*YH(I,L)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1, JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1, J
 30       IC = IC*JJ
 35     C = float(IC)
        DO 40 I = 1, N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      CALL SSCAL (N, R, DKY, 1)
      RETURN
C
 80   MSG = 'VINDY--  K (=I1) illegal      '
      CALL XERRWV (MSG, 30, 51, 1, 1, K, 0, 0, ZERO, ZERO)
      IFLAG = -1
      RETURN
 90   MSG = 'VINDY--  T (=R1) illegal      '
      CALL XERRWV (MSG, 30, 52, 1, 0, 0, 0, 1, T, ZERO)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWV (MSG, 60, 52, 1, 0, 0, 0, 2, TP, TN)
      IFLAG = -2
      RETURN
C----------------------- End of Subroutine VINDY -----------------------
      END
      SUBROUTINE VSTEP (Y, YH, LDYH, YH1, EWT, SAVF, VSAV, ACOR,
     1                  WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR)
      EXTERNAL F, JAC, PSOL, VNLS
      real Y, YH, YH1, EWT, SAVF, VSAV, ACOR, WM, RPAR
      INTEGER LDYH, IWM, IPAR
      DIMENSION Y(*), YH(LDYH,*), YH1(*), EWT(*), SAVF(*), VSAV(*),
     1   ACOR(*), WM(*), IWM(*), RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C Call sequence input -- Y, YH, LDYH, YH1, EWT, SAVF, VSAV,
C                        ACOR, WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR
C Call sequence output -- YH, ACOR, WM, IWM
C COMMON block variables accessed..
C     /VOD001/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13),
C               TQ(5), TN, JCUR, JSTART, KFLAG, KUTH,
C               L, LMAX, MAXORD, MITER, N, NEWQ, NQ, NQWAIT
C     /VOD002/  HU, NCFN, NETF, NFE, NQU, NST
C
C Subroutines called by VSTEP.. F, SAXPY, SCOPY, SSCAL,
C                               VJUST, VNLS, VSET
C Function routines called by VSTEP.. VNORML
C-----------------------------------------------------------------------
C VSTEP performs one step of the integration of an initial value
C problem for a system of ordinary differential equations.
C VSTEP calls subroutine VNLS for the solution of the nonlinear system
C arising in the time step.  Thus it is independent of the problem
C Jacobian structure and the type of nonlinear system solution method.
C VSTEP returns a completion flag KFLAG (in COMMON).
C A return with KFLAG = -1, -2, or -5 means either ABS(H) = HMIN or 10
C consecutive failures occurred.  On a return with KFLAG negative,
C the values of TN and the YH array are as of the beginning of the last
C step, and H is the last step size attempted.
C
C Communication with VSTEP is done with the following variables..
C
C Y      = An array of length N used for the dependent variable vector.
C YH     = An LDYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C LDYH   = A constant integer .ge. N, the first dimension of YH.
C          N is the number of ODEs in the system.
C YH1    = A one-dimensional array occupying the same space as YH.
C EWT    = An array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = An array of working storage, of length N.
C          also used for input of YH(*,MAXORD+2) when JSTART = -1
C          and MAXORD .lt. the current order NQ.
C VSAV   = A work array of length N passed to subroutine VNLS.
C ACOR   = A work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = Real and integer work arrays associated with matrix
C          operations in VNLS.
C F      = Dummy name for the user supplied subroutine for f.
C JAC    = Dummy name for the user supplied Jacobian subroutine.
C PSOL   = Dummy name for the subroutine passed to VNLS, for
C          possible use there.
C VNLS   = Dummy name for the nonlinear system solving subroutine,
C          whose real name is dependent on the method used.
C RPAR, IPAR = Dummy names for user's real and integer work arrays.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block VOD001 --------------------
C
      real ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block VOD002 --------------------
C
      real HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
C Type declarations for local variables --------------------------------
C
      real ADDON, BIAS1,BIAS2,BIAS3, CNQUOT, DDN, DSM, DUP,
     1     ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF,
     2     ETAQ, ETAQM1, ETAQP1, FLOTL, ONE, ONEPSM,
     3     R, THRESH, TOLD, ZERO
      INTEGER I, I1, I2, IBACK, J, JB, KFC, KFH, MXNCF, NCF, NFLAG
      INTEGER IFAIL
C
C Type declaration for function subroutines called ---------------------
C
      real VNORML
C
C Type declaration for the additional diagnostics set up to inquire ----
C during the running of vodpk. -----------------------------------------
C
      character*80 msgjm
      integer nrcv, ierrjm, ijmgetmr, imxer2, imxer3
      real big2, big3, size2
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE ADDON, BIAS1, BIAS2, BIAS3,
     1     ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF,
     2     KFC, KFH, MXNCF, ONEPSM, THRESH, ONE, ZERO, ETAQ, ETAQM1
C-----------------------------------------------------------------------
      COMMON /VOD001/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /VOD002/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
C
      DATA KFC/-3/, KFH/-7/, MXNCF/10/
      DATA ADDON  /1.0E-6/,    BIAS1  /6.0E0/,     BIAS2  /6.0E0/,
     1     BIAS3  /10.0E0/,    ETACF  /0.25E0/,    ETAMIN /0.1E0/,
     2     ETAMXF /0.2E0/,     ETAMX1 /1.0E4/,     ETAMX2 /10.0E0/,
     3     ETAMX3 /10.0E0/,    ONEPSM /1.00001E0/, THRESH /1.5E0/
      DATA ONE/1.0E0/, ZERO/0.0E0/
C
      KFLAG = 0
      TOLD = TN
      NCF = 0
      JCUR = 0
      NFLAG = 0
      IF (JSTART .GT. 0) GO TO 20
      IF (JSTART .EQ. -1) GO TO 100
C-----------------------------------------------------------------------
C On the first call, the order is set to 1, and other variables are
C initialized.  ETAMAX is the maximum ratio by which H can be increased
C in a single step.  It is normally 1.5, but is larger during the
C first 10 steps to compensate for the small initial H.  If a failure
C occurs (in corrector convergence or error test), ETAMAX is set to 1
C for the next increase.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      NQNYH = NQ*LDYH
      TAU(1) = H
      PRL1 = ONE
      RC = ZERO
      ETAMAX = ETAMX1
      NQWAIT = 2
      HSCAL = H
      IFAIL = 0
      GO TO 200
C-----------------------------------------------------------------------
C Take preliminary actions on a normal continuation step (JSTART.GT.0).
C If the driver changed H, then ETA must be reset and NEWH set to 1.
C If a change of order was dictated on the previous step, then
C it is done here and appropriate adjustments in the history are made.
C On an order decrease, the history array is adjusted by VJUST.
C On an order increase, the history array is augmented by a column.
C On a change of step size H, the history array YH is rescaled.
C-----------------------------------------------------------------------
 20   CONTINUE
      IF (KUTH .EQ. 1) THEN
        ETA = min(ETA,H/HSCAL)
        NEWH = 1
        ENDIF
 50   IF (NEWH .EQ. 0) GO TO 200
      IF (NEWQ .EQ. NQ) GO TO 150
      IF (NEWQ .LT. NQ) THEN
        CALL VJUST (YH, LDYH, -1)
        NQ = NEWQ
        L = NQ + 1
        NQWAIT = L
        GO TO 150
        ENDIF
      IF (NEWQ .GT. NQ) THEN
        CALL VJUST (YH, LDYH, 1)
        NQ = NEWQ
        L = NQ + 1
        NQWAIT = L
        GO TO 150
      ENDIF
C-----------------------------------------------------------------------
C The following block handles preliminaries needed when JSTART = -1.
C If N was reduced, zero out part of YH to avoid undefined references.
C If MAXORD was reduced to a value less than the tentative order NEWQ,
C then NQ is set to MAXORD, and a new H ratio ETA is chosen.
C Otherwise, we take the same preliminary actions as for JSTART .gt. 0.
C In any case, NQWAIT is reset to L = NQ + 1 to prevent further
C changes in order for that many steps.
C The new H ratio ETA is limited by the input H if KUTH = 1,
C by HMIN if KUTH = 0, and by HMXI in any case.
C Finally, the history array YH is rescaled.
C-----------------------------------------------------------------------
 100  CONTINUE
      LMAX = MAXORD + 1
      IF (N .EQ. LDYH) GO TO 120
      I1 = 1 + (NEWQ + 1)*LDYH
      I2 = (MAXORD + 1)*LDYH
      IF (I1 .GT. I2) GO TO 120
      DO 110 I = I1, I2
 110    YH1(I) = ZERO
 120  IF (NEWQ .LE. MAXORD) GO TO 140
      FLOTL = float(LMAX)
      IF (MAXORD .LT. NQ-1) THEN
        DDN = VNORML (N, SAVF, EWT)/TQ(1)
        ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
        ENDIF
      IF (MAXORD .EQ. NQ .AND. NEWQ .EQ. NQ+1) ETA = ETAQ
      IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ+1) THEN
        ETA = ETAQM1
        CALL VJUST (YH, LDYH, -1)
        ENDIF
      IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ) THEN
        DDN = VNORML (N, SAVF, EWT)/TQ(1)
        ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
        CALL VJUST (YH, LDYH, -1)
        ENDIF
      ETA = min(ETA,ONE)
      NQ = MAXORD
      L = LMAX
 140  IF (KUTH .EQ. 1) ETA = min(ETA,ABS(H/HSCAL))
      IF (KUTH .EQ. 0) ETA = max(ETA,HMIN/ABS(HSCAL))
      ETA = ETA/max(ONE,ABS(HSCAL)*HMXI*ETA)
      NEWH = 1
      NQWAIT = L
      IF (NEWQ .LE. MAXORD) GO TO 50
C Rescale the history array for a change in H by a factor of ETA. ------
 150  R = ONE
      DO 180 J = 2, L
        R = R*ETA
        CALL SSCAL (N, R, YH(1,J), 1 )
 180    CONTINUE
      H = HSCAL*ETA
      HSCAL = H
      RC = RC*ETA
      NQNYH = NQ*LDYH
C-----------------------------------------------------------------------
C This section computes the predicted values by effectively
C multiplying the YH array by the Pascal triangle matrix.
C VSET is called to calculate all integration coefficients.
C RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
C-----------------------------------------------------------------------
 200  TN = TN + H
      I1 = NQNYH + 1
      DO 220 JB = 1, NQ
        I1 = I1 - LDYH
        DO 210 I = I1, NQNYH
 210      YH1(I) = YH1(I) + YH1(I+LDYH)
 220  CONTINUE
      CALL VSET
      RL1 = ONE/EL(2)
      RC = RC*(RL1/PRL1)
      PRL1 = RL1
C
C Call the nonlinear system solver. ------------------------------------
C
      CALL VNLS (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM,
     1           F, JAC, PSOL, NFLAG, RPAR, IPAR)
C
      IF (NFLAG .EQ. 0) GO TO 450
C-----------------------------------------------------------------------
C The VNLS routine failed to achieve convergence (NFLAG .NE. 0).
C The YH array is retracted to its values before prediction.
C The step size H is reduced and the step is retried, if possible.
C Otherwise, an error exit is taken.
C-----------------------------------------------------------------------
        NCF = NCF + 1
        NCFN = NCFN + 1
        ETAMAX = ONE
        TN = TOLD
        I1 = NQNYH + 1
        DO 430 JB = 1, NQ
          I1 = I1 - LDYH
          DO 420 I = I1, NQNYH
 420        YH1(I) = YH1(I) - YH1(I+LDYH)
 430      CONTINUE
        IF (NFLAG .LT. -1) GO TO 680
        IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 670
        IF (NCF .EQ. MXNCF) GO TO 670
        ETA = ETACF
        ETA = max(ETA,HMIN/ABS(H))
        NFLAG = -1
        GO TO 150
C-----------------------------------------------------------------------
C The corrector has converged (NFLAG = 0).  The local error test is
C made and control passes to statement 500 if it fails.
C-----------------------------------------------------------------------
 450  CONTINUE
      DSM = ACNRM/TQ(2)
C-----------------------------------------------------------------------
C This section is here to obtain information while the integrator ------
C is running,  The keywords are:
C              s       give information about the status
C              kaboom  returns to the user prompt, this option is useful
C                      when the solution failes to converge.
C This section was added by JLM 3/7/92
C-----------------------------------------------------------------------
ccTDR  Disable because PYTHON versions get stuck looking for signal 12/3/03?
ccTDR      ierrjm = ijmgetmr(msgjm,80,1,nrcv)
ccTDR      if (ierrjm .eq. 0) then
ccTDR         if (msgjm(1:nrcv).eq.'status' .or. msgjm(1:nrcv).eq.'s' ) then
ccTDR           big2 = 0.0e0
ccTDR           imxer2 = 1
ccTDR           do 452 i = 1,n
ccTDR              size2 = abs( acor(i)*ewt(i) )
ccTDR              if (big2 .ge. size2) go to 452
ccTDR              big2 = size2
ccTDR              imxer2 = i
ccTDR 452       continue
ccTDR           big3 = 0.0e0
ccTDR           imxer3 = 1
ccTDR           do 453 i = 1,n
ccTDR              size2 = abs( vsav(i)*ewt(i) )
ccTDR              if (big2 .ge. size2) go to 453
ccTDR              big3 = size2
ccTDR              imxer3 = i
ccTDR 453       continue
ccTDR            write(*,*) 'nfe = ', nfe, '   npe = ', nje, '   imxtstep = '
ccTDR     .           , imxer2,'   imxnewt = ', imxer3
ccTDR            write(*,*) 'time = ', tn-h, '   dt = ', h
ccTDR            write(*,*) 'bigts = ', big2, '   dsm = ', dsm
ccTDR         elseif(msgjm(1:nrcv).eq.'kaboom' .or. msgjm(1:nrcv).eq.'k')then
ccTDR            call xerrab("")
ccTDR         endif
ccTDR      endif
C-----------------------------------------------------------------------
      IF (DSM .GT. ONE) GO TO 500
C-----------------------------------------------------------------------
C After a successful step, update the YH and TAU arrays and decrement
C NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
C for use in a possible order increase on the next step.
C If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2.
C-----------------------------------------------------------------------
      KFLAG = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO 470 IBACK = 1, NQ
        I = L - IBACK
 470    TAU(I+1) = TAU(I)
      TAU(1) = H
      DO 480 J = 1, L
        CALL SAXPY (N, EL(J), ACOR, 1, YH(1,J), 1 )
 480    CONTINUE
      NQWAIT = NQWAIT - 1
      IF ((L .EQ. LMAX) .OR. (NQWAIT .NE. 1)) GO TO 490
      CALL SCOPY (N, ACOR, 1, YH(1,LMAX), 1 )
      CONP = TQ(5)
 490  IF (ETAMAX .NE. ONE) GO TO 560
      IF (NQWAIT .LT. 2) NQWAIT = 2
      NEWQ = NQ
      NEWH = 0
      ETA = ONE
      HNEW = H
      GO TO 690
C-----------------------------------------------------------------------
C The error test failed.  KFLAG keeps track of multiple failures.
C Restore TN and the YH array to their previous values, and prepare
C to try the step again.  Compute the optimum step size for the
C same order.  After repeated failures, H is forced to decrease
C more rapidly.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      NETF = NETF + 1
      NFLAG = -2
      TN = TOLD
      I1 = NQNYH + 1
      DO 520 JB = 1, NQ
        I1 = I1 - LDYH
        DO 510 I = I1, NQNYH
 510      YH1(I) = YH1(I) - YH1(I+LDYH)
 520  CONTINUE
      IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 660
      ETAMAX = ONE
      IF (KFLAG .LE. KFC) GO TO 530
C Compute ratio of new H to current H at the current order. ------------
      FLOTL = float(L)
      ETA = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
      ETA = max(ETA,HMIN/ABS(H),ETAMIN)
      IF ((KFLAG .LE. -2) .AND. (ETA .GT. ETAMXF)) ETA = ETAMXF
      GO TO 150
C-----------------------------------------------------------------------
C Control reaches this section if 3 or more consecutive failures
C have occurred.  It is assumed that the elements of the YH array
C have accumulated errors of the wrong order.  The order is reduced
C by one, if possible.  Then H is reduced by a factor of 0.1 and
C the step is retried.  After a total of 7 consecutive failures,
C an exit is taken with KFLAG = -1.
C-----------------------------------------------------------------------
 530  IF (KFLAG .EQ. KFH) GO TO 660
      IF (NQ .EQ. 1) GO TO 540
      ETA = max(ETAMIN,HMIN/ABS(H))
      CALL VJUST (YH, LDYH, -1)
      L = NQ
      NQ = NQ - 1
      NQWAIT = L
      GO TO 150
 540  ETA = max(ETAMIN,HMIN/ABS(H))
      H = H*ETA
      HSCAL = H
      TAU(1) = H
      IFAIL = 0
      CALL F (N, TN, Y, SAVF, RPAR, IPAR, IFAIL)
      NFE = NFE + 1
      IF (IFAIL .NE. 0) GO TO 670
      DO 550 I = 1, N
 550    YH(I,2) = H*SAVF(I)
      NQWAIT = 10
      GO TO 200
C-----------------------------------------------------------------------
C If NQWAIT = 0, an increase or decrease in order by one is considered.
C Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
C be multiplied at order q, q-1, or q+1, respectively.
C The largest of these is determined, and the new order and
C step size set accordingly.
C A change of H or NQ is made only if H increases by at least a
C factor of THRESH.  If an order change is considered and rejected,
C then NQWAIT is set to 2 (reconsider it after 2 steps).
C-----------------------------------------------------------------------
C Compute ratio of new H to current H at the current order. ------------
 560  FLOTL = float(L)
      ETAQ = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
      IF (NQWAIT .NE. 0) GO TO 600
      NQWAIT = 2
      ETAQM1 = ZERO
      IF (NQ .EQ. 1) GO TO 570
C Compute ratio of new H to current H at the current order less one. ---
      DDN = VNORML (N, YH(1,L), EWT)/TQ(1)
      ETAQM1 = ONE/((BIAS1*DDN)**(ONE/(FLOTL - ONE)) + ADDON)
 570  ETAQP1 = ZERO
      IF (L .EQ. LMAX) GO TO 580
C Compute ratio of new H to current H at current order plus one. -------
      CNQUOT = (TQ(5)/CONP)*(H/TAU(2))**L
      DO 575 I = 1, N
 575    SAVF(I) = ACOR(I) - CNQUOT*YH(I,LMAX)
      DUP = VNORML (N, SAVF, EWT)/TQ(3)
      ETAQP1 = ONE/((BIAS3*DUP)**(ONE/(FLOTL + ONE)) + ADDON)
 580  IF (ETAQ .GE. ETAQP1) GO TO 590
      IF (ETAQP1 .GT. ETAQM1) GO TO 620
      GO TO 610
 590  IF (ETAQ .LT. ETAQM1) GO TO 610
 600  ETA = ETAQ
      NEWQ = NQ
      GO TO 630
 610  ETA = ETAQM1
      NEWQ = NQ - 1
      GO TO 630
 620  ETA = ETAQP1
      NEWQ = NQ + 1
      CALL SCOPY (N, ACOR, 1, YH(1,LMAX), 1)
C Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ----
 630  IF (ETA .LT. THRESH .OR. ETAMAX .EQ. ONE) GO TO 640
      ETA = min(ETA,ETAMAX)
      ETA = ETA/max(ONE,ABS(H)*HMXI*ETA)
      NEWH = 1
      HNEW = H*ETA
      GO TO 690
 640  NEWQ = NQ
      NEWH = 0
      ETA = ONE
      HNEW = H
      GO TO 690
C-----------------------------------------------------------------------
C All returns are made through this section.
C On a successful return, ETAMAX is reset and ACOR is scaled.
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      IF (IFAIL .NE. 0) KFLAG = -5
      GO TO 720
 680  IF (NFLAG .EQ. -2) KFLAG = -3
      IF (NFLAG .EQ. -3) KFLAG = -4
      IF (NFLAG .EQ. -4) KFLAG = -5
      GO TO 720
 690  ETAMAX = ETAMX3
      IF (NST .LE. 10) ETAMAX = ETAMX2
 700  R = ONE/TQ(2)
      CALL SSCAL (N, R, ACOR, 1)
 720  JSTART = 1
      RETURN
C----------------------- End of Subroutine VSTEP -----------------------
      END
      SUBROUTINE VSET
C-----------------------------------------------------------------------
C Call sequence communication.. None
C COMMON block variables accessed..
C     /VOD001/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1),
C                 METH, NQ, NQWAIT
C
C Subroutines called by VSET.. None
C Function routines called by VSET.. None
C-----------------------------------------------------------------------
C VSET is called by VSTEP and sets coefficients for use there.
C
C For each order NQ, the coefficients in EL are calculated by use of
C  the generating polynomial lambda(x), with coefficients EL(i).
C      lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
C For the backward differentiation formulas,
C                                     NQ-1
C      lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .
C                                     i = 1
C For the Adams formulas,
C                              NQ-1
C      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,
C                              i = 1
C      lambda(-1) = 0,    lambda(0) = 1,
C where c is a normalization constant.
C In both cases, xi(i) is defined by
C      H*xi(i) = t sub n  -  t sub (n-i)
C              = H + TAU(1) + TAU(2) + ... TAU(i-1).
C
C
C In addition to variables described previously, communication
C with VSET uses the following..
C   TAU    = A vector of length 13 containing the past NQ values
C            of H.
C   EL     = A vector of length 13 in which vset stores the
C            coefficients for the corrector formula.
C   TQ     = A vector of length 5 in which vset stores constants
C            used for the convergence test, the error test, and the
C            selection of H at a new order.
C   METH   = The basic method indicator.
C   NQ     = The current order.
C   L      = NQ + 1, the length of the vector stored in EL, and
C            the number of columns of the YH array being used.
C   NQWAIT = A counter controlling the frequency of order changes.
C            An order change is about to be considered if NQWAIT = 1.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block VOD001 --------------------
C
      real ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for local variables --------------------------------
C
      real AHATN0, ALPH0, CNQM1, CORTES, CSUM, ELP, EM,
     1     EM0, FLOTI, FLOTL, FLOTNQ, HSUM, ONE, RXI, RXIS, S, SIX,
     2     T1, T2, T3, T4, T5, T6, TWO, XI, ZERO
      INTEGER I, IBACK, J, JP1, NQM1, NQM2
C
      DIMENSION EM(13)
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE CORTES, ONE, SIX, TWO, ZERO
C
      COMMON /VOD001/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
C
      DATA CORTES /0.1E0/
      DATA ONE  /1.0E0/, SIX /6.0E0/, TWO /2.0E0/, ZERO /0.0E0/
C
      FLOTL = float(L)
      NQM1 = NQ - 1
      NQM2 = NQ - 2
      GO TO (100, 200), METH
C
C Set coefficients for Adams methods. ----------------------------------
 100  IF (NQ .NE. 1) GO TO 110
      EL(1) = ONE
      EL(2) = ONE
      TQ(1) = ONE
      TQ(2) = TWO
      TQ(3) = SIX*TQ(2)
      TQ(5) = ONE
      GO TO 300
 110  HSUM = H
      EM(1) = ONE
      FLOTNQ = FLOTL - ONE
      DO 115 I = 2, L
 115    EM(I) = ZERO
      DO 150 J = 1, NQM1
        IF ((J .NE. NQM1) .OR. (NQWAIT .NE. 1)) GO TO 130
        S = ONE
        CSUM = ZERO
        DO 120 I = 1, NQM1
          CSUM = CSUM + S*EM(I)/float(I+1)
 120      S = -S
        TQ(1) = EM(NQM1)/(FLOTNQ*CSUM)
 130    RXI = H/HSUM
        DO 140 IBACK = 1, J
          I = (J + 2) - IBACK
 140      EM(I) = EM(I) + EM(I-1)*RXI
        HSUM = HSUM + TAU(J)
 150    CONTINUE
C Compute integral from -1 to 0 of polynomial and of x times it. -------
      S = ONE
      EM0 = ZERO
      CSUM = ZERO
      DO 160 I = 1, NQ
        FLOTI = float(I)
        EM0 = EM0 + S*EM(I)/FLOTI
        CSUM = CSUM + S*EM(I)/(FLOTI+ONE)
 160    S = -S
C In EL, form coefficients of normalized integrated polynomial. --------
      S = ONE/EM0
      EL(1) = ONE
      DO 170 I = 1, NQ
 170    EL(I+1) = S*EM(I)/float(I)
      XI = HSUM/H
      TQ(2) = XI*EM0/CSUM
      TQ(5) = XI/EL(L)
      IF (NQWAIT .NE. 1) GO TO 300
C For higher order control constant, multiply polynomial by 1+x/xi(q). -
      RXI = ONE/XI
      DO 180 IBACK = 1, NQ
        I = (L + 1) - IBACK
 180    EM(I) = EM(I) + EM(I-1)*RXI
C Compute integral of polynomial. --------------------------------------
      S = ONE
      CSUM = ZERO
      DO 190 I = 1, L
        CSUM = CSUM + S*EM(I)/float(I+1)
 190    S = -S
      TQ(3) = FLOTL*EM0/CSUM
      GO TO 300
C
C Set coefficients for BDF methods. ------------------------------------
 200  DO 210 I = 3, L
 210    EL(I) = ZERO
      EL(1) = ONE
      EL(2) = ONE
      ALPH0 = -ONE
      AHATN0 = -ONE
      HSUM = H
      RXI = ONE
      RXIS = ONE
      IF (NQ .EQ. 1) GO TO 240
      DO 230 J = 1, NQM2
C In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
        HSUM = HSUM + TAU(J)
        RXI = H/HSUM
        JP1 = J + 1
        ALPH0 = ALPH0 - ONE/float(JP1)
        DO 220 IBACK = 1, JP1
          I = (J + 3) - IBACK
 220      EL(I) = EL(I) + EL(I-1)*RXI
 230    CONTINUE
      ALPH0 = ALPH0 - ONE/float(NQ)
      RXIS = -EL(2) - ALPH0
      HSUM = HSUM + TAU(NQM1)
      RXI = H/HSUM
      AHATN0 = -EL(2) - RXI
      DO 235 IBACK = 1, NQ
        I = (NQ + 2) - IBACK
 235    EL(I) = EL(I) + EL(I-1)*RXIS
 240  T1 = ONE - AHATN0 + ALPH0
      T2 = ONE + float(NQ)*T1
      TQ(2) = ABS(ALPH0*T2/T1)
      TQ(5) = ABS(T2/(EL(L)*RXI/RXIS))
      IF (NQWAIT .NE. 1) GO TO 300
      CNQM1 = RXIS/EL(L)
      T3 = ALPH0 + ONE/float(NQ)
      T4 = AHATN0 + RXI
      ELP = T3/(ONE - T4 + T3)
      TQ(1) = ABS(ELP/CNQM1)
      HSUM = HSUM + TAU(NQ)
      RXI = H/HSUM
      T5 = ALPH0 - ONE/float(NQ+1)
      T6 = AHATN0 - RXI
      ELP = T2/(ONE - T6 + T5)
      TQ(3) = ABS(ELP*RXI*(FLOTL + ONE)*T5)
 300  TQ(4) = CORTES*TQ(2)
      RETURN
C----------------------- End of Subroutine VSET ------------------------
      END
      SUBROUTINE VJUST (YH, LDYH, IORD)
      real YH
      INTEGER LDYH, IORD
      DIMENSION YH(LDYH,*)
C-----------------------------------------------------------------------
C Call sequence input -- YH, LDYH, IORD
C Call sequence output -- YH
C COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N
C COMMON block variables accessed..
C     /VOD001/ -- HSCAL, TAU(13), LMAX, METH, N, NQ,
C
C Subroutines called by VJUST.. SAXPY
C Function routines called by VJUST.. None
C-----------------------------------------------------------------------
C This subroutine adjusts the YH array on reduction of order,
C and also when the order is increased for the stiff option (METH = 2).
C Communication with VJUST uses the following..
C IORD  = An integer flag used when METH = 2 to indicate an order
C         increase (IORD = +1) or an order decrease (IORD = -1).
C HSCAL = Step size H used in scaling of Nordsieck array YH.
C         (If IORD = +1, VJUST assumes that HSCAL = TAU(1).)
C See References 1 and 2 for details.
C-----------------------------------------------------------------------
C
C Type declarations for labeled COMMON block VOD001 --------------------
C
      real ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for local variables --------------------------------
C
      real ALPH0,ALPH1, HSUM, ONE, PROD, T1, XI,XIOLD, ZERO
      INTEGER I, IBACK, J, JP1, LP1, NQM1, NQM2, NQP1
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this integrator.
C-----------------------------------------------------------------------
      SAVE ONE, ZERO
C
      COMMON /VOD001/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
C
      DATA ONE /1.0E0/, ZERO /0.0E0/
C
      IF ((NQ .EQ. 2) .AND. (IORD .NE. 1)) RETURN
      NQM1 = NQ - 1
      NQM2 = NQ - 2
      GO TO (100, 200), METH
C-----------------------------------------------------------------------
C Nonstiff option...
C Check to see if the order is being increased or decreased.
C-----------------------------------------------------------------------
 100  CONTINUE
      IF (IORD .EQ. 1) GO TO 180
C Order decrease. ------------------------------------------------------
      DO 110 J = 1, LMAX
 110    EL(J) = ZERO
      EL(2) = ONE
      HSUM = ZERO
      DO 130 J = 1, NQM2
C Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
        HSUM = HSUM + TAU(J)
        XI = HSUM/HSCAL
        JP1 = J + 1
        DO 120 IBACK = 1, JP1
          I = (J + 3) - IBACK
 120      EL(I) = EL(I)*XI + EL(I-1)
 130    CONTINUE
C Construct coefficients of integrated polynomial. ---------------------
      DO 140 J = 2, NQM1
 140    EL(J+1) = float(NQ)*EL(J)/float(J)
C Subtract correction terms from YH array. -----------------------------
      DO 170 J = 3, NQ
        DO 160 I = 1, N
 160      YH(I,J) = YH(I,J) - YH(I,L)*EL(J)
 170    CONTINUE
      RETURN
C Order increase. ------------------------------------------------------
C Zero out next column in YH array. ------------------------------------
 180  CONTINUE
      LP1 = L + 1
      DO 190 I = 1, N
 190    YH(I,LP1) = ZERO
      RETURN
C-----------------------------------------------------------------------
C Stiff option...
C Check to see if the order is being increased or decreased.
C-----------------------------------------------------------------------
 200  CONTINUE
      IF (IORD .EQ. 1) GO TO 300
C Order decrease. ------------------------------------------------------
      DO 210 J = 1, LMAX
 210    EL(J) = ZERO
      EL(3) = ONE
      HSUM = ZERO
      DO 230 J = 1, NQM2
C Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
        HSUM = HSUM + TAU(J)
        XI = HSUM/HSCAL
        JP1 = J + 1
        DO 220 IBACK = 1, JP1
          I = (J + 4) - IBACK
 220      EL(I) = EL(I)*XI + EL(I-1)
 230    CONTINUE
C Subtract correction terms from YH array. -----------------------------
      DO 250 J = 3, NQ
        DO 240 I = 1, N
 240      YH(I,J) = YH(I,J) - YH(I,L)*EL(J)
 250    CONTINUE
      RETURN
C Order increase. ------------------------------------------------------
 300  DO 310 J = 1, LMAX
 310    EL(J) = ZERO
      EL(3) = ONE
      ALPH0 = -ONE
      ALPH1 = ONE
      PROD = ONE
      XIOLD = ONE
      HSUM = HSCAL
      IF (NQ .EQ. 1) GO TO 340
      DO 330 J = 1, NQM1
C Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
        JP1 = J + 1
        HSUM = HSUM + TAU(JP1)
        XI = HSUM/HSCAL
        PROD = PROD*XI
        ALPH0 = ALPH0 - ONE/float(JP1)
        ALPH1 = ALPH1 + ONE/XI
        DO 320 IBACK = 1, JP1
          I = (J + 4) - IBACK
 320      EL(I) = EL(I)*XIOLD + EL(I-1)
        XIOLD = XI
 330    CONTINUE
 340  CONTINUE
      T1 = (-ALPH0 - ALPH1)/PROD
C Load column L + 1 in YH array. ---------------------------------------
      LP1 = L + 1
      DO 350 I = 1, N
 350    YH(I,LP1) = T1*YH(I,LMAX)
C Add correction terms to YH array. ------------------------------------
      NQP1 = NQ + 1
      DO 370 J = 3, NQP1
        CALL SAXPY (N, EL(J), YH(1,LP1), 1, YH(1,J), 1 )
 370  CONTINUE
      RETURN
C----------------------- End of Subroutine VJUST -----------------------
      END
      real FUNCTION VNORML (N, V, W)
      real V, W
      INTEGER N
      DIMENSION V(N), W(N)
C-----------------------------------------------------------------------
C Call sequence input -- N, V, W
C Call sequence output -- None
C COMMON block variables accessed -- None
C
C Subroutines/functions called by VNORML.. None
C-----------------------------------------------------------------------
C This function routine computes the weighted root-mean-square norm
C of the vector of length N contained in the array V, with weights
C contained in the array W of length N..
C     VNORML = sqrt( (1/N) * sum( V(i)*W(i) )**2 )
CCC commented section can calculate the maximum of V(i)*W(i)  TDR 05/27/92
C-----------------------------------------------------------------------
      real SUM, big, size
      INTEGER I
      big = 0.0e0
      do 10 I = 1, N
         size = abs( V(I)*W(I) )
         if(big .ge. size) go to 10
         big = size
 10   continue
      VNORML = big
C
ccc      SUM = 0.0E0
ccc         DO 10 I = 1, N
ccc  10     SUM = SUM + (V(I)*W(I))**2
ccc      VNORML = SQRT(SUM/float(N))
      RETURN
C----------------------- End of Function VNORML -------------------------
      END
