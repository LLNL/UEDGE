      SUBROUTINE SYLD96(MATT,MATP,CION,CIZB,CRMB)
      IMPLICIT none
      INTEGER MATT,MATP,CION,CIZB
      real    CRMB

c############################################################################
c ** - modifications by TDR 2/20/98
c 
c ** - Coding received from David Elder, Feb. 6, 1998; reduced argument
c ** - list for subroutine (deleted cneutd,cbombf,cbombz,cebd); deleted 
c ** - params, Use statement for cyield instead of include statement
c ** - REAL -> real for mppl; ALOG10 -> log10
c 
c ** -Input variables:
c     cion  # integer atomic number of impurity, e.g., for carbon, cion=6
c     cizb  # integer max charge state of plasma ions; hydrogen cizb=1
c     crmb  # real mass number of plasma ions; =2 for deut., =3 for tritium
c
c ** -Output variables:
c     matt  # integer flag giving target material
c     matp  # integer flag giving plasma material
c
c###########################################################################   
C
C  *********************************************************************
C  *                                                                   *
C  *  SYLD96:  SETS UP SPUTTERING YIELD DATA IN COMMON AREA CYIELD.    *
C  *  THE DATA IS TAKEN FROM ECKSTEIN IPP 9/82 (FEB 1993)              *
C  *                                                                   *
C  *  ARGUMENTS :-                                                     *
C  *  MATT : CODE INTEGER GIVING MATERIAL USED IN TARGET               *
C  *  MATP : CODE INTEGER GIVING BACKGROUND MATERIAL                   *
C  *  CNEUTD : NEUT SPUTTER OPTION - deleted from this version 2/20/98 *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
ctdr      include    'params'
C     INCLUDE "CYIELD"
cdtr      include    'cyield'

      Use(Cyield)   # ceth,cetf,cq,ntars,cidata

      real ETH(7,12), ETF(7,12), Q(7,12), EBD(12)
      LOGICAL IDATA(7,12)
      INTEGER I,J,NSPEC
      CHARACTER*18 TARMAT(19)
      CHARACTER*6  PLAMAT(7)
C
      NSPEC=7
      NTARS=12
C
C  NSPEC = NUMBER OF IMPURITY SPECIES IN PLASMA.
C  NTARS = NUMBER OF TARGET MATERIALS.
C   CETH = THRESHOLD ENERGY FOR TARGET-ION INTERACTION CONSIDERED (EV)
C   CETF = THOMAS-FERMI INTERACTION POTENTIAL (EV)
C     CQ = YIELD FACTOR (ATOMS/ION)
C CIDATA = LOGICAL FLAG INDICATING WHETHER DATA IS AVAILABLE (T OR F)
C
C     DATA FROM ECKSTEINS LATEST REPORT IPP 9/82 HAS BEEN USED.
C     THE TWO COMPOUNDS (TITANIUM CARBIDE AND SILICON CARBIDE)
C     IN THE ORIGINAL REFERENCE ARE NOT AVAILABLE IN THE LATEST
C     REPORT AND HAVE BEEN REPLACED WITH LITHIUM AND CHROMIUM.
C     NOTE ALSO THAT ECKSTEIN HAS CHANGED THE DEFINITION OF HIS
C     NUCLEAR STOPPING CROSS SECTION.  HE HAS ALSO RECOMMENDED
C     THAT WE USE FITS TO EXPERIMENTAL DATA, WHERE AVAILABLE, IN
C     PREFERENCE TO THE CALCULATIONS BASED ON EMPIRICAL FORMULAE
C     (PP 335-338, WHICH ARE AVAILABLE FOR A LARGE RANGE OF
C     PROJECTILE-TARGET COMBINATIONS).  FOR THE TIME BEING I HAVE
C     REPLACED THE GENERAL TABLES WITH FITS TO EXPERIMENTAL DATA
C     ONLY FOR H, D, AND HE ON BE.  THIS RESULTS IN HIGHER YIELDS
C     (APPROX 2X).
C                       LORNE HORTON MAY 93
C
 
      DATA TARMAT/
     &  ' ALUMINIUM       ',' BERYLLIUM       ',' COPPER          ',
     &  ' GRAPHITE        ',' TITANIUM        ',' IRON            ',
     &  ' NICKEL          ',' MOLYBDENUM      ',' TUNGSTEN        ',
     &  ' BORON           ',' LITHIUM         ',' CHROMIUM        ',
     &  ' "DEUTERIUM"     ',' "HELIUM"        ',' "NEON"          ',
     &  ' "ARGON"         ',' "OXYGEN"        ',' "CHLORINE"      ',
     &  ' "NITROGEN"      ' /
 
      DATA PLAMAT/
     &  ' H    ',' D    ',' T    ',' HE4  ',' C    ',' SELF ',' O    '/
 
      DATA  ETH/
     &     23.87, 14.86, 12.91, 12.51, 16.32, 24.02, 18.55,
     &     12.2 , 10.0 , 14.69, 13.9 , 28.08, 24.17, 32.71,
     &     57.25, 28.90, 20.64, 17.07, 13.27, 25.17, 14.01,
     &     31.00, 27.00, 29.00, 32.15, 52.98, 52.98, 61.54,
     &     59.49, 31.51, 23.71, 20.56, 19.45, 34.96, 21.23,
     &     61.39, 31.63, 23.12, 19.54, 16.70, 31.03, 17.95,
     &     66.80, 34.12, 24.69, 20.67, 17.00, 31.89, 18.14,
     &    172.36, 83.30, 56.47, 44.28, 25.75, 48.83, 25.47,
     &    447.02,209.37,136.26,102.07, 41.20, 62.06, 35.92,
     &     23.14, 21.56, 23.46, 25.83, 43.25, 40.97, 50.30,
     &      6.22,  6.92,  8.03,  9.10, 15.94, 11.94, 18.61,
     &     54.47, 28.39, 21.01, 17.96, 16.07, 29.46, 17.40/
C
      DATA ETF/
     &    1059.,   1097.,   1135.,   2448.,  10297.,  34550.,  15720.,
     &     256.,    282.,    308.,    720.,   4153.,   2208.,   6971.,
     &    2926.,   2972.,   3017.,   6293.,  22701., 224652.,  32727.,
     &     415.,    447.,    479.,   1087.,   5688.,   5688.,   9298.,
     &    2054.,   2097.,   2139.,   4503.,  16949., 117915.,  24846.,
     &    2544.,   2590.,   2635.,   5517.,  20270., 174122.,  29437.,
     &    2799.,   2846.,   2893.,   6045.,  22014., 206991.,  31860.,
     &    4719.,   4768.,   4817.,   9945.,  34188., 533127.,  48329.,
     &    9871.,   9925.,   9978.,  20376.,  66517.,1998893.,  91993.,
     &     333.,    361.,    389.,    894.,   4856.,   3717.,   8021.,
     &     185.,    209.,    232.,    557.,   3506.,   1129.,   6014.,
     &    2296.,   2340.,   2383.,   5002.,  18577., 144458.,  27091./
C
      DATA Q/
     &     0.08,  0.14,  0.19,  0.37,  1.65,  4.21,  2.36,
     &     0.128, 0.220, 0.14,  0.707, 1.00,  0.67,  1.35,
     &     0.08,  0.14,  0.19,  0.38,  1.83, 14.23,  2.73,
     &     0.035, 0.10,  0.12,  0.20,  0.75,  0.75,  1.02,
     &     0.06,  0.11,  0.15,  0.30,  1.41,  7.44,  2.07,
     &     0.07,  0.12,  0.16,  0.33,  1.59, 10.44,  2.36,
     &     0.07,  0.12,  0.16,  0.33,  1.60, 11.51,  2.38,
     &     0.05,  0.09,  0.12,  0.24,  1.20, 16.27,  1.81,
     &     0.04,  0.07,  0.10,  0.20,  1.02, 33.47,  1.55,
     &     0.05,  0.08,  0.11,  0.21,  0.80,  0.67,  1.08,
     &     0.10,  0.16,  0.21,  0.40,  1.37,  0.69,  1.82,
     &     0.07,  0.12,  0.17,  0.34,  1.61,  9.54,  2.38/
C
      DATA IDATA/
     &            7*.true.,
     &            7*.true.,
     &            7*.true.,
     &            7*.true.,
     &            7*.true.,
     &            7*.true.,
     &            7*.true.,
     &            7*.true.,
     &            7*.true.,
     &            7*.true.,
     &            7*.true.,
     &            7*.true./
C
C  TABLE OF BINDING ENERGIES TO BE USED AS DEFAULT WHEN ZERO IS
C  SPECIFIED IN THE INPUT FILE.  FOR THE GASEOUS IMPURITIES I HAVE
C  SET EBD = 0
C
      DATA EBD/3.36,3.38,3.52,7.42,4.89,4.34,
     &         4.46,6.83,8.68,5.73,1.67,4.12/
C
C-----------------------------------------------------------------------
C INITIALISE COMMON BLOCK CYIELD.
C-----------------------------------------------------------------------
C
      DO 10 I=1,NTARS
        DO 20 J=1,NSPEC
         CETH(J,I)   = ETH(J,I)
         CETF(J,I)   = ETF(J,I)
         CQ(J,I)     = Q(J,I)
         CIDATA(J,I) = IDATA(J,I)
   20   CONTINUE
   10 CONTINUE
C
C-----------------------------------------------------------------------
C  ASSIGN TARGET AND BACKGROUND ION MATERIALS.
C  MATERIALS 13-18 ARE SPECIAL "GAS" TARGET CASES, WHERE THE YIELD IS
C  ALWAYS 1.0
C-----------------------------------------------------------------------
C
      MATT = 4
      IF (CION.EQ.13) MATT = 1
      IF (CION.EQ.4)  MATT = 2
      IF (CION.EQ.29) MATT = 3
      IF (CION.EQ.6)  MATT = 4
      IF (CION.EQ.22) MATT = 5
      IF (CION.EQ.26) MATT = 6
      IF (CION.EQ.28) MATT = 7
      IF (CION.EQ.42) MATT = 8
      IF (CION.EQ.74) MATT = 9
      IF (CION.EQ.5)  MATT = 10
      IF (CION.EQ.3)  MATT = 11
      IF (CION.EQ.24) MATT = 12
      IF (CION.EQ.1)  MATT = 13
      IF (CION.EQ.2)  MATT = 14
      IF (CION.EQ.10) MATT = 15
      IF (CION.EQ.18) MATT = 16
      IF (CION.EQ.8)  MATT = 17
      IF (CION.EQ.17) MATT = 18
      IF (CION.EQ.7)  MATT = 19
C
C---- TARGET BINDING ENERGY FROM INTERNAL DATA IF EXTERNAL INPUT ZERO
C
ctdr      IF (CEBD.EQ.0.0 .AND. MATT.LE.12) CEBD = EBD(MATT)
C
C---- PLASMA MATERIAL.  CAN BE SET EXPLICITLY FROM INPUT DATA
C
      MATP = 6
      IF (NINT(CRMB).LE.4) MATP = NINT (CRMB)
      IF (CIZB.EQ.6) MATP = 5
      IF (CIZB.EQ.8) MATP = 7
cdtr      IF (CNEUTD.EQ.1) MATP = CBOMBF
      WRITE (*,*) 'TARGET MATERIAL IS     ' , TARMAT(MATT)
      WRITE (*,*) 'BOMBARDING IONS ARE    ' , PLAMAT(MATP)
cdtr      IF (CNEUTD.EQ.1) CALL PRI ('         WITH ZIMP', CBOMBZ)
      RETURN
      END
C
C
C
c ------------------------------------------------------------------------c
       real FUNCTION YLD96(MATP,MATT,ENERGY)
       IMPLICIT none
       real ENERGY,X1,X12,X2
       INTEGER MATP,MATT
C
C  *********************************************************************
C  *                                                                   *
C  *  YLD96:   CALCULATES YIELD OF MATERIAL "MATP" HITTING TARGET      *
C  *  MADE OF MATERIAL "MATT" WITH AN ENERGY "ENERGY".                 *
C  *                                                                   *
C  *  ARGUMENTS :-                                                     *
C  *  MATP : CODE INTEGER GIVING MATERIAL OF IMPACTING PARTICLE        *
C  *  MATT : CODE INTEGER GIVING MATERIAL USED IN TARGET               *
C  *  ENERGY : ENERGY OF IMPACTING PARTICLE                            *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "CYIELD"
cdtr      include    'cyield'
      Use(Cyield)    # ceth,cetf,cq,ntars,cidata

C
      IF (MATT.EQ.13.OR.MATT.EQ.14.OR.MATT.EQ.15
     >      .OR.MATT.EQ.16.OR.MATT.EQ.17.OR.MATT.EQ.18
     >      .or.matt.eq.19) THEN
        YLD96 = 1.0
        RETURN
      ENDIF
C
      IF(ENERGY.LE.0.0) GOTO 100
      IF(CIDATA(MATP,MATT)) THEN
      IF((CETH(MATP,MATT)/ENERGY).GT.1.0E0) GOTO 100
         X1=ENERGY/CETF(MATP,MATT)
         X12=SQRT(X1)
         X2=CETH(MATP,MATT)/ENERGY
         YLD96=CQ(MATP,MATT)*(0.5*LOG(1.0+1.2288*X1))
     1         /(X1+0.1728*X12+0.008*X1**0.1504)
     2         *(1-X2)*(1-X2)*(1.0-X2**(2.0/3.0))
      ELSE
         YLD96=0.0
      ENDIF
      RETURN
C  ERROR TRAPPING, OUTSIDE RANGE SET YIELD=0,   E=0 OR E/E0 > 1
  100 YLD96=0.0
      END
c
c ------------------------------------------------------------------------c
      SUBROUTINE SPUTCHEM(IOPTCHEM,E0,TEMP,FLUX,YCHEM)

      IMPLICIT NONE
      INTEGER  IOPTCHEM
      real     E0,TEMP,FLUX,YCHEM,YGARCIA,YHAASZ,YROTH96,
     >         YHAASZ97,YHAASZ97M,FLUX_cgs
      INTRINSIC MAX
      Use(Cyield)     #redf_haas

c#########################################################################
c
c ** -Modified 2/20/98 to use flux in MKS [1/m**2 s] rather than previous
c ** -[1/cm**2 s]; TDR 
c 
c ** -Code obtained from David Elder, 2/6/98; originally written by
c ** -Houyang Guo at JET 
c
c#########################################################################
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FOR D --> C                                  *
C  *                                                                   *  
C  *  IOPTCHEM       -  Options for chemical sputtering:               *
C  *         1       -  Garcia-Rosales' formula (EPS94)                *
C  *         2       -  according to Pospieszczyk (EPS95)              *
C  *         3       -  Vietzke (in Physical processes of the inter-   *
C  *                    action of Fusion Plasmas with Solids)          *
C  *         4       -  Haasz (Submitted to J.Nucl.Mater.,Dec. 1995)   *
C  *         5       -  Roth & Garcia-Rosales (Submitted to Nucl.      *
C  *                    Fusion, March 1996)                            *
C  *         6       -  Haasz 1997 (Brian Mech's PhD Thesis)           *
C  *         7       -  Haasz 1997 + reduced 1/5 from 10->5 eV (Porter)*
C  *                                                                   *
C  *  E0   (eV)      -  Ion or neutral incident energy                 *
C  *  TEMP (K)       -  Temperature at target or wall                  *
C  *  FLUX (m-2s-1)  -  Ion or neutral flux  (changed to MKS 2/20/98)  *
C  *  YCHEM          -  Chemical Sputtering yield                      *
C  *                                                                   *
C  *********************************************************************
C
c
c      Change the flux from MKS to cgs units to use old cgs routines
c      conversion done on 2/20/98
      FLUX_cgs  = 1e4*FLUX
 
      IF      (IOPTCHEM.EQ.1) THEN
        YCHEM = YGARCIA(E0,TEMP,FLUX_cgs)
      ELSE IF (IOPTCHEM.EQ.2) THEN
        YCHEM = 0.04254*(MAX(5E18,FLUX_cgs)/5E18)**(-0.477)
      ELSE IF (IOPTCHEM.EQ.3) THEN
        YCHEM = 0.0215*(MAX(1E14,FLUX_cgs)/1E16)**(-0.1)
      ELSE IF (IOPTCHEM.EQ.4) THEN
        YCHEM = YHAASZ(E0,TEMP)
      ELSE IF (IOPTCHEM.EQ.5) THEN
        YCHEM = YROTH96(E0,TEMP,FLUX_cgs)
      ELSE IF (IOPTCHEM.EQ.6) THEN
        YCHEM = YHAASZ97(E0,TEMP)
      ELSE IF (IOPTCHEM.EQ.7) THEN
        YCHEM = YHAASZ97M(E0,TEMP,redf_haas)
      END IF

      RETURN
      END

C -------------
c ------------------------------------------------------------------------c
      FUNCTION YROTH96(E0,TEMP,FLUX)

      IMPLICIT NONE
      real    E0,TEMP,FLUX
      real    ETHC,ETFC,QC,SNC
      real    CSURF,CSP3
      real    YPHYS,YSURF,YTHERM,YROTH96
      INTRINSIC MIN
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING CALCULATED BY Garcia-Rosales FORMULA        *
C  *                                                                   *  
C  *  ETHC (eV)  -  Threshold energy for D -> C physical sputtering    *
C  *  ETFC (eV)  -  Thomas-Fermi energy                                *
C  *  SNC        -  Stopping power                                     *
C  *  QC         -  Fitting parameters                                 *
C  *                                                                   *
C  *  CSURF      -  Fitting parameters                                 *
C  *  CSP3       -  Carbon at surface                                  *
C  *                                                                   *
C  *  YPHYS      -  Physical sputtering yield                          *
C  *  YSURF      -  Sputtering due to SURFACE process                  *
C  *  YTHERM     -  Sputtering due to THERMAL process                  *
C  *  YROTH96    -  Total CHEMICAL sputtering yield                    *
C  *                                                                   *
C  *********************************************************************

C ---------------------------------------------------------
C Total Chemical Sputtering Yield:
C            Ychem = Ysurf+Ytherm*(1+D*Yphys)
C ---------------------------------------------------------

C 
C 1> PHYSICAL SPUTTERING YIELD
C
      ETHC = 27.0
      ETFC = 447.0      
      QC   = 0.1
C
C  - Stopping Power      
C
      SNC = 0.5*LOG(1.+1.2288*E0/ETFC)/(E0/ETFC
     >    + 0.1728*SQRT(E0/ETFC)
     >    + 0.008*(E0/ETFC)**0.1504)
C
C  - Physical Sputtering Yield
C
      IF (E0.GT.ETHC) THEN
         YPHYS = QC*SNC*(1-(ETHC/E0)**(2./3.))*(1-ETHC/E0)**2
      ELSE
         YPHYS = 0.0
      ENDIF
C
C 2> YSURF: Surface Process 
C
      CSURF  = 1/(1.+1E13*EXP(-2.45*11604/TEMP))
      CSP3   = CSURF*(2E-32*FLUX+EXP(-1.7*11604/TEMP))
     >         /(2E-32*FLUX+(1+2E29/FLUX*EXP(-1.8*11604/TEMP))
     >         *EXP(-1.7*11604/TEMP))

      IF (E0.GT.1.) THEN
         YSURF = CSP3*QC*SNC*(1-(1./E0)**(2./3.))*(1-1./E0)**2
     >           /(1.+EXP((MIN(90.0,E0)-90.)/50.))
      ELSE
         YSURF = 0.0
      ENDIF
C
C 3> YTHERM: Thermak Activated Process
C
      YTHERM = CSP3*0.033*EXP(-1.7*11604/TEMP)
     >         /(2E-32*FLUX+EXP(-1.7*11604/TEMP))
C
C 4> YCHEM: Total Chemical Sputtering Yield
C
      YROTH96 = YSURF + YTHERM * (1 + 125 * YPHYS)

CW    WRITE(6,*) 'YROTH96 = ',YROTH96

      RETURN
      END
      
      
C -------------
c ------------------------------------------------------------------------c
      FUNCTION YGARCIA(E0,TEMP,FLUX)

      IMPLICIT NONE
      real    E0,TEMP,FLUX
      real    ETHC,ETFC,QC,SNC
      real    YPHYS,YCHEM_TH,YCHEM_ATH,YGARCIA
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING CALCULATED BY Garcia-Rosales FORMULA        *
C  *                                                                   *  
C  *  ETHC (eV)  -  Threshold energy for D -> C physical sputtering    *
C  *  ETFC (eV)  -  Thomas-Fermi energy                                *
C  *  SNC        -  Stopping power                                     *
C  *  QC         -  Fitting parameters                                 *
C  *                                                                   *
C  *  YPHYS      -  Physical sputtering yield                          *
C  *  YCHEM_TH   -  Thermal activated mechanism                        *
C  *  YCHEM_ATH  -  Athermal mechanism                                 *
C  *                                                                   *
C  *********************************************************************
C            

      ETHC = 27.0
      ETFC = 447.0      
      QC   = 0.1
C
C Check for impact energies below threshold
C
      IF (E0.GT.ETHC) THEN
C
C Stopping Power      
C
      SNC = 0.5*LOG(1.+1.2288*E0/ETFC)/(E0/ETFC
     >    + 0.1728*SQRT(E0/ETFC)
     >    + 0.008*(E0/ETFC)**0.1504)
C 
C Physical Sputtering Yield
C
         YPHYS = QC*SNC*(1-(ETHC/E0)**(2./3.))*(1-ETHC/E0)**2
C
      ELSE
         YPHYS = 0.0
      ENDIF
C
C Chemical Sputtering Yield
C
      YCHEM_TH = 6E19*EXP(-1.*11604/TEMP)
     >         /(1E15+3E27*EXP(-2.*11604/TEMP))
     >         * (2.+200*YPHYS)*(MAX(1E16,FLUX)/1E16)**(-0.1)

      YCHEM_ATH = 0.05*EXP(E0*1E-3*(20.-1*11604./TEMP))
     >          / ((1.+EXP((E0-150.)/25.))*(1.+EXP((TEMP-740.)/25.)))

      YGARCIA = YCHEM_TH + YCHEM_ATH

CW    WRITE(6,*) 'YTH = ',YCHEM_TH,'YATH = ',YCHEM_ATH

      RETURN
      END
      
      
C -------------
c ------------------------------------------------------------------------c

      FUNCTION YHAASZ(E0,TEMP)
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz NEW DATA (Dec. 1995)            *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *                                                                   *
C  *********************************************************************
C            
      IMPLICIT NONE

      real    E0,TEMP
      real    FITC300(4),FITC350(4),FITC400(4),FITC450(4),FITC500(4),
     >        FITC550(4),FITC600(4),FITC650(4),FITC700(4),FITC750(4),
     >        FITC800(4),FITC850(4),FITC900(4),FITC950(4),FITC1000(4),
     >        FITC1050(4),FITC1100(4)
      real    POLY_C(4),YFIT,FITE0
      real    YHAASZ
      INTEGER I
C
C     Poly. fit c. /       a0,      a1,      a2,      a3
C
      DATA FITC300 / -0.01789, 0.02309, 0.00089,-0.00315/
      DATA FITC350 / -0.01691, 0.02020, 0.00451,-0.00407/
      DATA FITC400 / -0.01128, 0.01230, 0.00922,-0.00483/
      DATA FITC450 / -0.00401, 0.00453, 0.01226,-0.00493/
      DATA FITC500 / -0.01000, 0.02097,-0.00032,-0.00153/
      DATA FITC550 / -0.02022, 0.04019,-0.01430, 0.00253/
      DATA FITC600 /  0.00047,-0.00319, 0.00950,-0.00025/
      DATA FITC650 /  0.02921,-0.05657, 0.03467,-0.00226/
      DATA FITC700 /  0.00561,-0.00081,-0.01044, 0.00939/
      DATA FITC750 /  0.00225, 0.00205,-0.00949, 0.00800/
      DATA FITC800 /  0.00900,-0.02109, 0.01366, 0.00048/
      DATA FITC850 /  0.00483,-0.01691, 0.01513,-0.00152/
      DATA FITC900 /  0.00569,-0.02211, 0.02185,-0.00427/
      DATA FITC950 /  0.00317,-0.01827, 0.02081,-0.00482/
      DATA FITC1000/  0.00436,-0.02075, 0.02290,-0.00574/
      DATA FITC1050/  0.00463,-0.02082, 0.02285,-0.00601/
      DATA FITC1100/  0.00920,-0.02942, 0.02802,-0.00723/
C
C Find right polynomial fit coefficients for a given temperature
C
      IF      (TEMP.LE.300.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC300(I)
         ENDDO
      ELSE IF (TEMP.GT.300.0 .AND. TEMP.LE.350.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC350(I)
         ENDDO
      ELSE IF (TEMP.GT.350.0 .AND. TEMP.LE.400.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC400(I)
         ENDDO
      ELSE IF (TEMP.GT.400.0 .AND. TEMP.LE.450.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC450(I)
         ENDDO
      ELSE IF (TEMP.GT.450.0 .AND. TEMP.LE.500.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC500(I)
         ENDDO
      ELSE IF (TEMP.GT.500.0 .AND. TEMP.LE.550.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC550(I)
         ENDDO
      ELSE IF (TEMP.GT.550.0 .AND. TEMP.LE.600.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC600(I)
         ENDDO
      ELSE IF (TEMP.GT.600.0 .AND. TEMP.LE.650.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC650(I)
         ENDDO
      ELSE IF (TEMP.GT.650.0 .AND. TEMP.LE.700.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC700(I)
         ENDDO
      ELSE IF (TEMP.GT.700.0 .AND. TEMP.LE.750.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC750(I)
         ENDDO
      ELSE IF (TEMP.GT.750.0 .AND. TEMP.LE.800.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC800(I)
         ENDDO
      ELSE IF (TEMP.GT.800.0 .AND. TEMP.LE.850.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC850(I)
         ENDDO
      ELSE IF (TEMP.GT.850.0 .AND. TEMP.LE.900.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC900(I)
         ENDDO
      ELSE IF (TEMP.GT.900.0 .AND. TEMP.LE.950.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC950(I)
         ENDDO
      ELSE IF (TEMP.GT.950.0 .AND. TEMP.LE.1000.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC1000(I)
         ENDDO
      ELSE IF (TEMP.GT.1000.0 .AND. TEMP.LE.1050.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC1050(I)
         ENDDO
      ELSE
         DO I = 1,4
           POLY_C(I) = FITC1100(I)
         ENDDO
      ENDIF
C
C Calculate chemical yield according to the 3th poly. fit
C
      IF      (E0.LT.10.0)  THEN
         FITE0 = 10.
      ELSE IF (E0.GT.200.0) THEN
         FITE0 = 200.
      ELSE
         FITE0 = E0
      ENDIF 
C
      YFIT = 0.0
      DO I = 1,4
        YFIT = YFIT + POLY_C(I)*log10(FITE0)**(I-1)
      ENDDO

      YHAASZ = YFIT

CW    WRITE(6,*) 'YHAASZ = ',YHAASZ

      RETURN
      END
c
c
c
C -------------
c ------------------------------------------------------------------------c

      FUNCTION YHAASZ97(E0,TEMP)
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz NEW DATA (February 1997)        *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *                                                                   *
C  *********************************************************************
C            
      IMPLICIT NONE

      real    E0,TEMP
      real    FITC300(4),FITC350(4),FITC400(4),FITC450(4),FITC500(4),
     >        FITC550(4),FITC600(4),FITC650(4),FITC700(4),FITC750(4),
     >        FITC800(4),FITC850(4),FITC900(4),FITC950(4),FITC1000(4)
      real    POLY_C(4),YFIT,FITE0
      real    YHAASZ97
      INTEGER I
C
C     Poly. fit c. /       a0,      a1,      a2,      a3
C
      DATA FITC300 / -0.03882, 0.07432,-0.03470, 0.00486/
      DATA FITC350 / -0.05185, 0.10126,-0.05065, 0.00797/
      DATA FITC400 / -0.06089, 0.12186,-0.06240, 0.01017/
      DATA FITC450 / -0.08065, 0.16884,-0.09224, 0.01625/
      DATA FITC500 / -0.08872, 0.19424,-0.10858, 0.01988/
      DATA FITC550 / -0.08728, 0.20002,-0.11420, 0.02230/
      DATA FITC600 / -0.05106, 0.13146,-0.07514, 0.01706/
      DATA FITC650 /  0.07373,-0.13263, 0.09571,-0.01672/
      DATA FITC700 /  0.02722,-0.03599, 0.02064, 0.00282/
      DATA FITC750 /  0.09052,-0.18253, 0.12362,-0.02109/
      DATA FITC800 /  0.02604,-0.05480, 0.04025,-0.00484/
      DATA FITC850 /  0.03478,-0.08537, 0.06883,-0.01404/
      DATA FITC900 /  0.02173,-0.06399, 0.05862,-0.01380/
      DATA FITC950 / -0.00086,-0.01858, 0.02897,-0.00829/
      DATA FITC1000/ -0.01551, 0.01359, 0.00600,-0.00353/
C
C Find right polynomial fit coefficients for a given temperature
C
      IF      (TEMP.LE.300.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC300(I)
         ENDDO
      ELSE IF (TEMP.GT.300.0 .AND. TEMP.LE.350.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC350(I)
         ENDDO
      ELSE IF (TEMP.GT.350.0 .AND. TEMP.LE.400.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC400(I)
         ENDDO
      ELSE IF (TEMP.GT.400.0 .AND. TEMP.LE.450.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC450(I)
         ENDDO
      ELSE IF (TEMP.GT.450.0 .AND. TEMP.LE.500.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC500(I)
         ENDDO
      ELSE IF (TEMP.GT.500.0 .AND. TEMP.LE.550.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC550(I)
         ENDDO
      ELSE IF (TEMP.GT.550.0 .AND. TEMP.LE.600.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC600(I)
         ENDDO
      ELSE IF (TEMP.GT.600.0 .AND. TEMP.LE.650.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC650(I)
         ENDDO
      ELSE IF (TEMP.GT.650.0 .AND. TEMP.LE.700.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC700(I)
         ENDDO
      ELSE IF (TEMP.GT.700.0 .AND. TEMP.LE.750.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC750(I)
         ENDDO
      ELSE IF (TEMP.GT.750.0 .AND. TEMP.LE.800.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC800(I)
         ENDDO
      ELSE IF (TEMP.GT.800.0 .AND. TEMP.LE.850.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC850(I)
         ENDDO
      ELSE IF (TEMP.GT.850.0 .AND. TEMP.LE.900.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC900(I)
         ENDDO
      ELSE IF (TEMP.GT.900.0 .AND. TEMP.LE.950.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC950(I)
         ENDDO
      ELSE IF (TEMP.GT.950.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC1000(I)
         ENDDO
      ENDIF
C
C Calculate chemical yield according to the 3th poly. fit
C
      IF      (E0.LT.10.0)  THEN
         FITE0 = 10.
      ELSE IF (E0.GT.200.0) THEN
         FITE0 = 200.
      ELSE
         FITE0 = E0
      ENDIF 
C
      YFIT = 0.0
      DO I = 1,4
        YFIT = YFIT + POLY_C(I)*log10(FITE0)**(I-1)
      ENDDO

      YHAASZ97 = YFIT

CW    WRITE(6,*) 'YHAASZ97 = ',YHAASZ97

      RETURN
      END
c
c
c
C -------------
c ------------------------------------------------------------------------c

      FUNCTION YHAASZ97M(E0,TEMP,reducf)
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz NEW DATA (February 1997)        *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *  with the addition of a new fit below 10 eV as suggested by       *
C  *  J.Davis and parameterized by G. Porter; now interpolates between *
C  *  5 and 10 eV to lower value (YDAVIS98), and is fixed below 5 eV.  *
C  *  Default value for reducf=0.2; change redf_haas                   *
C  *                                                                   *
C  *********************************************************************
C            
      IMPLICIT NONE

      real E0, TEMP, reducf
      real YHAASZ97M, YDAVIS98, YHAASZ97
      real m1,m2,m3,FRAC

      DATA m1/602.39/, m2/202.24/, m3/43.561/

      IF (E0 .GE. 10) THEN
         YHAASZ97M = YHAASZ97(E0,TEMP)
      ELSEIF (E0 .LT. 10. .AND. E0 .GE. 5.) THEN
         FRAC = (E0-5.)/5.
         YDAVIS98 = reducf/(m2*((TEMP/m1)**2 - 1)**2 + m3)
         YHAASZ97M = FRAC*YHAASZ97(E0,TEMP)+ (1.-FRAC)*YDAVIS98
      ELSEIF (E0 .LT. 5.) THEN
         YDAVIS98 = reducf/(m2*((TEMP/m1)**2 - 1)**2 + m3)
         YHAASZ97M = YDAVIS98   
      ENDIF

      RETURN
      END
