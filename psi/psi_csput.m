c/*
c
c=======================================================
c       Plasma Surface Interaction (PSI) package
c  for simulation of plasma-wall sputtering processes.
c
c                 PSI Version 1.2   
c              Date: March 6, 2007
c
c
c Copyright 2007 by A.Yu. Pigarov and R.D. Smirnov
c
c
c                   Developers:
c          A.Yu. Pigarov and R.D. Smirnov
c
c
c         E-mails: apigarov@ucsd.edu
c                  rsmirnov@ucsd.edu
c
c    PSI package has been developed as part of
c   the Dust Transport in Tokamaks (DUSTT) code:
c         
c [1] A.Yu.Pigarov, S.I. Krasheninnikov, T.K. Soboleva,
c     and T.D. Rognlien "Dust-particle transport in
c     tokamak edge plasmas", Phys. Plasmas 12 (2005), 122508.
c [2] R.D. Smirnov, A.Yu. Pigarov, M. Rosenberg,
c     S.I. Krasheninnikov, and D.A. Mendis "Modelling of
c     dynamics and transport of carbon dust particles in 
c     tokamaks", Plasma Phys. Control. Fusion 49 (2007), 347-371.
c
c
c    The PSI package is distributed in the hope 
c that it will be useful to computational scientists in field
c of fusion plasma and other plasma related studies.
c
c    The PSI package is FREE SOFTWARE. 
c You can use, copy, and modify this software for any purpose
c and without fee provided that the above copyright
c notice appear in all copies.
c
c    The PSI package is distributed "AS IS", 
c i.e. without any warranty including all implied warranties
c of merchantability and fitness.
c   
c=======================================================
c
c*/
      subroutine psiyspc(izpt,ampt,izsf,amsf)
      implicit none
      
      Use (psitab)
      Use (psitab_s)
      
      integer izpt,izsf
      real*8 ampt,amsf
c
c  constants for chemical sputtering
c
      common/psicsptc/tdcspt,picsp,ethy,fcsp,cspf
      
      real*8 tdcspt,picsp,ethy,fcsp,cspf
      real*8 z1,z2,am1,am2,es,s
      real*8 psiwsub
      real*8 psi_hydchemfac,psi_cncspt
      external psi_hydchemfac,psi_cncspt

      z1=dble(izpt)
      z2=dble(izsf)
      am1=ampt
      am2=amsf
      
      es=psiwsub(izsf)

          s = -1.d0
      picsp = dacos(s)

        s  = 4.d0*am1+am2
        s  = s*s
        s  = s/(4.d0*am2*am1)
      ethy = es*s

      fcsp = dsqrt(am2*ethy/am1)
      fcsp = dsqrt(fcsp/(20.d0*z1))
      cspf = psi_hydchemfac(izsf)*psi_cncspt(izsf)

      end


      real*8 FUNCTION psiYHAASZ(E0,TEMP)
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz's NEW DATA (February 1997)        *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *  with the addition of a new fit below 10 eV as suggested by       *
C  *  J.Davis and parameterized by G. Porter; now interpolates between *
C  *  5 and 10 eV to lower value (YDAVIS98), and is fixed below 5 eV.  *
C  *  Default value for reducf=0.2; change redf_haas                   *
C  *  E0   (eV)      -  Ion or neutral incident energy                 *
C  *  TEMP (K)       -  Temperature at target or wall C                *                                                                   *
C  *********************************************************************
C            
      IMPLICIT NONE

      real*8 E0, TEMP
      real*8 reducf
      real*8 psiYHAASZ97, psiYDAVIS98
      real*8 m1,m2,m3,FRAC

      DATA reducf /1.0d0/
      DATA m1/602.39d0/, m2/202.24d0/, m3/43.561d0/

      IF (E0 .GE. 10.d0) THEN
       psiYHAASZ = psiYHAASZ97(E0,TEMP)
      ELSEIF (E0 .LT. 10.d0 .AND. E0 .GE. 5.d0) THEN
       FRAC = (E0-5.d0)/5.d0
       psiYDAVIS98 = reducf/(m2*((TEMP/m1)**2 - 1.d0)**2 + m3)
       psiYHAASZ = FRAC*psiYHAASZ97(E0,TEMP)+ (1.d0-FRAC)*psiYDAVIS98
      ELSEIF (E0 .LT. 5.d0) THEN
       psiYDAVIS98 = reducf/(m2*((TEMP/m1)**2 - 1.d0)**2 + m3)
       psiYHAASZ = psiYDAVIS98   
      ENDIF
      RETURN
      END


      real*8 FUNCTION psiYHAASZ97(E0,TEMP)
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz's NEW DATA (February 1997)        *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *                                                                   *
C  *********************************************************************
C            
      IMPLICIT NONE


      real*8  E0,TEMP
      real*8  FITC300(4),FITC350(4),FITC400(4),FITC450(4),FITC500(4),
     >        FITC550(4),FITC600(4),FITC650(4),FITC700(4),FITC750(4),
     >        FITC800(4),FITC850(4),FITC900(4),FITC950(4),FITC1000(4)
      real*8  POLY_C(4),YFIT,FITE0
      real*8  FITDRV(4),DRV,FNC
      real*8  M1,M2,M3
      DATA    M1 /602.39D0/
      INTEGER I
C
C     Poly. fit c. /       a0,      a1,      a2,      a3
C
      DATA FITC300  / -0.03882, 0.07432,-0.03470, 0.00486 /
      DATA FITC350  / -0.05185, 0.10126,-0.05065, 0.00797 /
      DATA FITC400  / -0.06089, 0.12186,-0.06240, 0.01017 /
      DATA FITC450  / -0.08065, 0.16884,-0.09224, 0.01625 /
      DATA FITC500  / -0.08872, 0.19424,-0.10858, 0.01988 /
      DATA FITC550  / -0.08728, 0.20002,-0.11420, 0.02230 /
      DATA FITC600  / -0.05106, 0.13146,-0.07514, 0.01706 /
      DATA FITC650  /  0.07373,-0.13263, 0.09571,-0.01672 /
      DATA FITC700  /  0.02722,-0.03599, 0.02064, 0.00282 /
      DATA FITC750  /  0.09052,-0.18253, 0.12362,-0.02109 /
      DATA FITC800  /  0.02604,-0.05480, 0.04025,-0.00484 /
      DATA FITC850  /  0.03478,-0.08537, 0.06883,-0.01404 /
      DATA FITC900  /  0.02173,-0.06399, 0.05862,-0.01380 /
      DATA FITC950  / -0.00086,-0.01858, 0.02897,-0.00829 /
      DATA FITC1000 / -0.01551, 0.01359, 0.00600,-0.00353 /
C
      DATA FITDRV   / -2.93d-4, 6.43d-4,-4.59d-4, 9.52d-5 /
C
C Find right polynomial fit coefficients for a given temperature
C
      IF      (TEMP.LE.300.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC300(I)
         ENDDO
      ELSE IF (TEMP.GT.300.d0 .AND. TEMP.LE.350.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC350(I)
         ENDDO
      ELSE IF (TEMP.GT.350.d0 .AND. TEMP.LE.400.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC400(I)
         ENDDO
      ELSE IF (TEMP.GT.400.d0 .AND. TEMP.LE.450.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC450(I)
         ENDDO
      ELSE IF (TEMP.GT.450.d0 .AND. TEMP.LE.500.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC500(I)
         ENDDO
      ELSE IF (TEMP.GT.500.d0 .AND. TEMP.LE.550.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC550(I)
         ENDDO
      ELSE IF (TEMP.GT.550.d0 .AND. TEMP.LE.600.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC600(I)
         ENDDO
      ELSE IF (TEMP.GT.600.d0 .AND. TEMP.LE.650.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC650(I)
         ENDDO
      ELSE IF (TEMP.GT.650.d0 .AND. TEMP.LE.700.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC700(I)
         ENDDO
      ELSE IF (TEMP.GT.700.d0 .AND. TEMP.LE.750.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC750(I)
         ENDDO
      ELSE IF (TEMP.GT.750.d0 .AND. TEMP.LE.800.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC800(I)
         ENDDO
      ELSE IF (TEMP.GT.800.d0 .AND. TEMP.LE.850.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC850(I)
         ENDDO
      ELSE IF (TEMP.GT.850.d0 .AND. TEMP.LE.900.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC900(I)
         ENDDO
      ELSE IF (TEMP.GT.900.d0 .AND. TEMP.LE.950.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC950(I)
         ENDDO
      ELSE IF (TEMP.GT.950.d0) THEN
         DO I = 1,4
           POLY_C(I) = FITC1000(I)
         ENDDO
      ENDIF
C
C Calculate chemical yield according to the 3th poly. fit
C
      IF      (E0.LT.10.d0)  THEN
         FITE0 = 10.d0
      ELSE IF (E0.GT.200.d0) THEN
         FITE0 = 200.d0
      ELSE
         FITE0 = E0
      ENDIF 
C
      IF (TEMP.LE.1000.D0) THEN
       YFIT = 0.d0
       DO I = 1,4
         YFIT = YFIT + POLY_C(I)*dlog10(FITE0)**(I-1)
       ENDDO
      ELSE
       DRV = 0.D0
       FNC = 0.D0
       DO I = 1,4
         DRV = DRV + FITDRV(I)*dlog10(FITE0)**(I-1)
         FNC = FNC + POLY_C(I)*dlog10(FITE0)**(I-1)
       ENDDO
       M2=-64.2147D0*DRV/(FNC*FNC)
       M3=1.D0/FNC-2.21144D0*M2
       YFIT = 1.D0/(M2*((TEMP/M1)**2 - 1.d0)**2 + M3)
      ENDIF

      psiYHAASZ97 = YFIT
      RETURN

      END


      real*8 FUNCTION psiYROTH96(E0,TEMP,FLUX)
      IMPLICIT NONE

      real*8    E0,TEMP,FLUX
      real*8    ETHC,ETFC,QC,SNC
      real*8    CSURF,CSP3
      real*8    YPHYS,YSURF,YTHERM
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING CALCULATED BY Garcia-Rosales' FORMULA        *
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
C  *  E0   (eV)      -  Ion or neutral incident energy                 *
C  *  TEMP (K)       -  Temperature at target or wall                  *
C  *  FLUX (10^20 cm-2s-1)  -  Ion or neutral flux                     *
C  *  YCHEM          -  Chemical Sputtering yield                      *
C  *********************************************************************
C
C ---------------------------------------------------------
C Total Chemical Sputtering Yield:
C            Ychem = Ysurf+Ytherm*(1+D*Yphys)
C ---------------------------------------------------------
C 
C 1> PHYSICAL SPUTTERING YIELD
C
      ETHC = 27.d0
      ETFC = 447.d0      
      QC   = 0.1
C
C  - Stopping Power      
C
      SNC = 0.5*DLOG(1.d0+1.2288*E0/ETFC)/(E0/ETFC
     >    + 0.1728*DSQRT(E0/ETFC)
     >    + 0.008*(E0/ETFC)**0.1504)
C
C  - Physical Sputtering Yield
C
      IF (E0.GT.ETHC) THEN
         YPHYS = QC*SNC*(1.d0-(ETHC/E0)**(2./3.))*(1.d0-ETHC/E0)**2
      ELSE
         YPHYS = 0.d0
      ENDIF
C
C 2> YSURF: Surface Process 
C
      CSURF  = 1.d0/(1.+1.d13*DEXP(-2.45*11604./TEMP))
      CSP3   = CSURF*(2.d-12*FLUX+DEXP(-1.7*11604./TEMP))
     >         /(2.d-12*FLUX+(1.+2.d9/FLUX*DEXP(-1.8*11604./TEMP))
     >         *DEXP(-1.7*11604./TEMP))


      IF (E0.GT.1.d0) THEN
         YSURF = CSP3*QC*SNC*(1.-(1./E0)**(2./3.))*(1.-1./E0)**2
     >           /(1.+DEXP((DMIN1(90.d0,E0)-90.)/50.))
      ELSE
         YSURF = 0.d0
      ENDIF
C
C 3> YTHERM: Thermak Activated Process
C
      YTHERM = CSP3*0.033*DEXP(-1.7*11604./TEMP)
     >         /(2.d-12*FLUX+DEXP(-1.7*11604./TEMP))
C
C 4> YCHEM: Total Chemical Sputtering Yield
C
      psiYROTH96 = YSURF + YTHERM * (1.d0 + 125.d0 * YPHYS)
      RETURN

      END


      real*8 function psicspyugl(e,amu)
      implicit none
      real*8 e,amu

      common/psicsptc/tdcspt,picsp,ethy,fcsp,csph
      real*8 tdcspt,picsp,ethy,fcsp,csph

      real*8 gmax
      data gmax/1.5d0/

      real*8 g,ss,sss,t,x

         x=e/ethy

         ss = amu
      if(ss .lt. 0.d0) ss = -ss
      if(ss .gt. 0.999) then
       g=1.d0
      else if(ss .lt. 1.d-6) then
         g=0.d0
      else if(x .lt. 4.00000001d0) then
         g=1.d0
      else
         sss= ss*ss
         t  = 1.d0-sss
         if(t .lt. 0.d0) then
         sss= -t 
         else 
         sss=t
       endif
         t  = dsqrt(sss)/ss
         sss= 2.d0/picsp
         t  = datan(t)*sss
         if(t .lt. 0.d0 ) t = -t
         if(t .gt. 1.d0 ) then
         t = 0.d0
         else
         t = 1.d0 - t
       endif
         t  = t*dsqrt(x - 4.d0)
         sss= fcsp * dsqrt(t)
         g  = dexp( (-sss) * dlog(ss) )
       endif

      if(g .gt. gmax) g=gmax

      psicspyugl=g
      return

      end


      subroutine psicsptd(td)
      implicit none
      real*8 td

      common/psicsptc/tdcspt,picsp,ethy,fcsp,csph
      real*8 tdcspt,picsp,ethy,fcsp,csph

      tdcspt=td

      end


      real*8 function psiycsptf(FLUX, t)
      implicit none
      real*8 FLUX,t

      real*8 cv1, cv2
      data cv1/1.7d0/, cv2/1.8d0/

      real*8 s,g

      s=11600.0d0/(t+0.2d0)
      g=2.d-12*FLUX+DEXP(-cv1*s)

      psiycsptf=g/(g+2.d9/FLUX*DEXP(-(cv2+cv1)*s))
      return

      end


      real*8 function psicspy(e,amu)
      implicit none
      
      real*8 e,amu

      common/psicsptc/tdcspt,picsp,ethy,fcsp,cspf
      real*8 tdcspt,picsp,ethy,fcsp,cspf
      real*8 psiYHAASZ,psicspyugl
      real*8 y,yu

      y =psiYHAASZ(e,tdcspt)
      yu=psicspyugl(e,amu)
      psicspy=y*yu*cspf
      return

      end


      double precision function psi_hydchemfac(iz)
      implicit none
      integer iz
c---------------------------------------
c  chemical sputtering reduction factor
c
c         The factor reduces 
c     the chemical sputtering yield
c     calculated by analytic formula
c   for those periodic system elements
c   which do not produce hydrides.
c --------------------------------------
      double precision cfac(92)
      data cfac /
     ,   1.d0 , 1.d-3, 1.d-2, 1.d-2, 4.d-1,
     ,   1.d0 , 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,
     ,   1.d-3, 1.d-3 /

      psi_hydchemfac=cfac(iz)
      return

      end
