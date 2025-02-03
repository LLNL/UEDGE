c!include "../mppl.h"
c!include "../sptodp.h"
c!include "api.h"
!	SEE FMOMBAL SUBROUTINE FOR INPUT PARAMETER SPECIFICATIONS
	subroutine coulfric(amu,denz2,dloglam,mntau,zi_api,capm,capn,
     >	ela,elab,tempa)
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
	real capm(KXA,miso,KXA,miso),denz2(miso),
     >		 capn(KXA,KXA,miso,miso),mntau(miso,miso),
     >		amu(miso),tempa(miso),
     >		ela(*),elab(*),zi_api(miso,nzch)
c
c	COMPUTES COULOMB COLLISION FRICTION TERMS
c	ELA AND ELAB, WHERE THE FRICTION FORCE IS
c
c	R sub (ALPHA,Z) sup (M) == FRICTION FOR SPECIES ALPHA,
c					   CHARGE STATE Z (M-th FORCE)
c
c	= ZI(ALPHA,Z)*[ SUM(K) { ELA(M,K,ALPHA)*U sub (ALPHA,Z) (K) }
c
c			  + LAMBDA(ALPHA,M) ]
c
c	LAMBDA(ALPHA,M) = SUM(BETA) SUM(K) [ ELAB(M,K,ALPHA,BETA) *
c
c				    U sub (ALPHA) sup(2) (K) ]
c
c	ELA(M,K,ALPHA) = SUM(BETA)[ <<MA NA/TAU-AB>> CAPM(M,K,ALPHA,BETA) ]
c
c	ELAB(M,K,ALPHA,BETA) = <<MA NA/TAU-AB>> CAPN(M,K,ALPHA,BETA)
c
c*****	Calculate M''s and N''s for each isotope.
	call neomn(amu,capm,capn,tempa)
c*****	Evaluate calculating n/tauab.
	const=(4.0/3.0/sqrt(pi0))*4.0*pi0*(coulom/(4.0*pi0*epsilo))**2
     >	*dloglam
	do 100 misa=1,miso
c*****	Thermal velocity-m/s (TEMP in Joules)
	  vtherm=sqrt(2.0*tempa(misa)/amu(misa)/promas)
	  dza = (denz2(misa)*coulom) * ( const
     >	     /(promas*amu(misa)*vtherm**3))
	  do 90 misb=1,miso
	    mntau(misa,misb) = dza * (denz2(misb)*coulom)
   90   continue	
  100 continue
c*****	Calculate la(=ela) and lab (=elab)
	call neolab(mntau,capm,capn,ela,elab)
	return
	end
c---- End of subroutine coulfric ---------------------------------------
c-----------------------------------------------------------------------
	subroutine fmombal(amu,den,dloglam,epar,friction,gradp,gradt,
     >	nuion,nurec,qcond,qneut,ucond,uneut,umass,
     >	parcurrent,tempa,natomic,misotope,nchstate,ldir,fricc)
  	implicit real (a-h,o-z), integer (i-n)
        Use(Reduced_ion_constants)
	Use(Reduced_ion_variables)
        Use(Timing)   # istimingon,ttimpc
        real(Size4) sec4, gettime, tsval
	real amu(misotope), den(misotope*(nchstate+1)),
     >	gradp(misotope*nchstate), gradt(misotope*nchstate),
     >	nurec(misotope*nchstate), nuion(misotope*(nchstate+1)),
     >	tempa(misotope), qcond(*),
     >	friction(*), ucond(*), uneut(*), qneut(*), fricc(*)
	real amat(KXA*(MXMISO1)*KXA*(MXMISO1)), denz2(MXMISO),
     >	caplams(KXA*MXMISO), uresp(NBA*KMXZ), usave(KXA*MXNZCH*MXMISO)
	integer natomic(misotope)
        integer ndima
	save uresp, usave, caplams
!
!	****	DESCRIPTION OF INPUT VARIABLES	****
!
!
!	ALL PHYSICAL QUANTITIES ARE IN MKS UNITS, INCLUDING
!	TEMPERATURE IN JOULES (1 eV = 1.6022E-19 Joule)
!	DENSITY IN 1/METERS**3
!	ELECTRIC FIELD IN VOLTS/METER
!	VELOCITY IN METERS/SEC
!	HEAT FLUX IN JOULE/(METER**2 - SEC)
!
!	INTEGERS:
!
!	  MISOTOPE
!		NUMBER OF DISTINCT ISOTOPES (DIFFERENT MASSES).
!		THE ORDERING OF ALL ARRAYS A(ISO) BELOW MUST BE AS FOLLOWS:
!		FIRST ISOTOPE (ISO=1) <==> ELECTRONS
!		SUBSEQUENT ISOTOPES (2 <= ISO <= MISOTOPE)
!
!	  NATOMIC(1:MISOTOPE)
!		FOR IONS, THE NUCLEAR CHARGE OF ISOTOPE (ISO). FOR ELECTRONS,
!		NATOMIC(1) == 1.
!
!	  NCHSTATE
!		THE MAXIMUM OF NATOMIC(ISO) OVER ALL MISOTOPE ISOTOPES
!
!	  LDIR
!		IF LDIR >= 2, THEN A COMPLETE (NEW) MATRIX INVERSION
!		IS DONE TO COMPUTE THE FLOWS. FOR LDIR < 2, AN INCREMENTAL
!		MATRIX INVERSION IS DONE, BASED ON A PREVIOUS CALL TO
!		FMOMBAL WITH LDIR >= 2. THIS SAVES CPU TIME AND IS USEFUL
!		WHEN SMALL VARIABLE PERTURBATIONS
!		ARE MADE TO COMPUTE THE NUMERICAL JACOBIAN.
!
!
!	REAL VARIABLES:
!
!	  AMU(1:MISOTOPE)
!		ATOMIC MASS OF ISOTOPE (ISO), RELATIVE TO THE PROTON MASS.
!		FOR EXAMPLE, AMU(1) = 5.45E-4.
!
!	  DEN(1:MISOTOPE,0:NCHSTATE)	[1/M**3]
!		DEN(ISO,NZ) IS THE DENSITY OF THE ISOTOPE
!		(ISO) AND CHARGE STATE (NZ). NOTE THAT NZ=0 CORRESPONDS
!		TO THE NEUTRAL DENSITY (NEEDED FOR ATOMIC PHYSICS), AND
!		NZ = NATOM(ISO) TO THE FULLY STRIPPED ION. DEN
!		SHOULD BE PRESCRIBED AT THE SAME SPATIAL GRID POINT
!		REQUIRED FOR THE VELOCITIES AND HEAT FLUXES.
!
!	  DLOGLAM
!		LOG(LAMBDA) REQUIRED FOR THE COULOMB COLLISION RATE
!
!	  EPAR		[VOLTS/METER]
!		PARALLEL ELECTRIC FIELD (INPUT)
!
!	  FRICTION(1:MISOTOPE,1:NCHSTATE)		[JOULES/METER**4]
!		FRICTION(ISO,NZ) IS THE PARALLEL COMPONENT OF
!		THE FRICTION FORCE (COULOMB+ATOMIC COLLISIONS)
!		FOR ISOTOPE ISO, CHARGE STATE NZ (NOT INCLUDING NEUTRALS)
!
!         FRICC(1:MISOPTE,1:NCHSTATE,1:4)                [JOULES/METER**4]
!               COMPONENTS OF FRICTION WITH 
!               FRICC(,,1) = ela(1,)*usol(1,,) for velocity
!               FRICC(,,2) = ela(2,)*usol(2,,) for qcond
!               FRICC(,,3) = ela(3,)*usol(3,,) for h
!               FRICC(,,4) = caplam ? same for all isotope charge states
!               FRICC(,,5) = ioniz and recomb totals
!
!	  GRADP,GRADT(1:MISOTOPE,1:NCHSTATE)		[JOULES/METER**4]
!		GRADP,GRADT(ISO,NZ) IS THE PARALLEL COMPONENT OF
!		THE PRESSURE(TEMPERATURE) GRADIENT,
!
!		GRADP,T  =  B dot Grad (P,T) / |B|
!
!		EVALUATED AT THE SAME SPATIAL GRID POINT WHERE THE
!		PARALLEL SPEEDS AND HEAT FLUXES ARE TO BE COMPUTED, OF
!		CHARGE STATE NZ OF ISOTOPE ISO. NEUTRAL GRADIENTS
!		ARE NOT INCLUDED HERE
!
!	  NUION(1:MISOTOPE,0:NCHSTATE)	[1/SEC]
!		NUION(ISO,NZ) IS THE IONIZATION RATE (DUE TO
!		ELECTRON IMPACT) FOR ISOTOPE ISO, STATE NZ, FOR NZ = 0
!		TO NZ = NATOM(ISO)-1
!
!	  NUREC(1:MISOTOPE,1:NCHSTATE)	[1/SEC]
!		NUREC(ISO,NZ) IS THE RECOMBINATION RATE (DUE TO
!		ELECTRON IMPACT) FOR ISOTOPE ISO, STATE NZ, FOR NZ = 1
!		TO NZ = NATOM(ISO)
!
!	  PARCURRENT			[A/METER**2]
!		PARALLEL CURRENT AT PRESENT GRID POINT (OUTPUT, COMPUTED FROM OHMS LAW)
!
!	  QCOND(1:MISOTOPE,1:NCHSTATE)	[JOULE/(METER**2-SEC)]
!		QCOND(ISO,NZ) IS THE PARALLEL RANDOM HEAT FLUX FOR CHARGE
!		STATE NZ OF ISOTOPE ISO IN RESPONSE TO THE ELECTRIC
!		AND PRESSURE FORCES. (NEUTRALS NOT INCLUDED HERE)
!		THIS, TOGETHER WITH UCOND AND EPAR, COMPRISES THE
!		OUTPUT FROM THIS CALCULATION
!
!	  UCOND(1:MISOTOPE,1:NCHSTATE)	[METER/SEC]
!		PARALLEL FLOW SPEED (ABSOLUTE, NOT RELATIVE TO UMASS)
!		OF EACH CHARGE STATE FOR EACH ISOTOPE (NOT
!		INCLUDING NEUTRALS)
!
!	  QNEUT(1:MISOTOPE), UNEUT(1:MISOTOPE)
!		HEAT FLUX AND PARALLEL FLOW SPEED OF NEUTRAL
!		SPECIES FOR EACH ISOTOPE. IF QNEUT UNKNOWN, USE 0 FOR NOW
!
!	  TEMPA(1:MISOTOPE)		[JOULES]
!		TEMPA(ISO) IS THE TEMPERATURE AT GRID POINT
!		FOR ISOTOPE ISO
!
!	  UMASS		[METERS/SEC]
!		THE (EXTERNALLY PROVIDED) MEAN MASS FLOW ALONG
!		MAGNETIC FIELD LINES, WHICH IS THE SUM OVER ALL
!		IONIZED SPECIES (AND ELECTRONS) OF THE INDIVIDUAL
!		MASS-DENSITY WEIGHTED FLOWS:
!
!		UMASS = SUM ( MASS*DENSITY*U|| )
! 		        ________________________
!		        SUM ( MASS*DENSITY )
!
!		THE EQUATION FOR EVOLVING UMASS IN TIME IS OBTAINED
!		BY SUMMING THE INDIVIDUAL MOMENTUM EQUATIONS FOR ALL
!		IONS. THIS SUMMATION ANNIHILATES THE DOMINANT COULOMB
!		FRICTION TERMS AND INCLUDES RESIDUAL COLLISIONAL
!		COUPLING (DUE TO ATOMIC PHYSICS) TO THE NEUTRALS
!		POPULATION (THROUGH UNEUT,QNEUT).
!

c ... Start timing impurity calculations.
      if (istimingon .eq. 1 .and. misotope .gt. 1) tsval = gettime(sec4)

c ... Initialize amat array.
        ndima=KXA*(MXMISO1)*KXA*(MXMISO1)
        do i=1,ndima
          amat(i)=zero
        enddo

c****	Setup zi_api, mass-density, Z*e*density arrays
	call setden(amu,den,denmass,denz,denz2,zi_api)

c****	Setup force arrays (gradp, gradt stuff)
	call setforce(den,denz,denmass,epar,gradp,gradt,
     >	tempa,uneut,qneut,umass,fmomenta,nuion)

c*****	Compute Coulomb friction coefficients la, lab
	call coulfric(amu,denz2,dloglam,mntau,zi_api,capm,capn,
     >	ela,elab,tempa)

c*****	Compute response matrices for inversion of momentum equations
c*****	in charge-state space
	call zrespond(den,denmass,zi_api,ela,nuion,nurec,uresp,
     >	usave,caplams,fmomenta,ldir)

c*****	Solve average-ion equations in mass-isotope space
	call mrespond(elab,zi_api,den,denmass,uresp,sbar,
     >	KXA*(miso+1),amat,nurec,umass,ldir)

c*****	Compute individual ion flows (for all charge states)
c
c	NOTE: USOL(m=1,Z,mass) = PARTICLE FLOW RELATIVE TO UMASS
c	
c	      USOL(m=2,Z,mass) = -(2/5)*CONDUCTIVE HEAT FLUX / PRESSURE
c
c	      USOL(m=3,Z,mass) = HIGHER-ORDER FLOW (H)
c
	call mzrespond(elab,uresp,sbar,caplam,caplams,usol,
     >	usave,den,denz,parcurrent,qcond,ucond,tempa,umass,ldir)
c
c	COMPUTE FRICTION FORCES
c	
	call getfrict(friction,fricc,caplam,denmass,ela,nuion,
     .  nurec,usol,zi_api)

c ... End timing impurity calculations.
        if (istimingon .eq. 1 .and. misotope .gt. 1)
     .     ttimpc = ttimpc + gettime(sec4) - tsval

	return
	end
c---- End of subroutine fmombal ----------------------------------------
c-----------------------------------------------------------------------
	subroutine getfrict(friction,fricc,caplam,denmass,ela,
     >	nuion,nurec,usol,zi_api)
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
	real friction(miso,nzch), usol(KXA,nzch,miso), 
     .       fricc(miso,nzch,5)
	real nuion(miso,0:nzch), nurec(miso,nzch)
	real denmass(miso,0:nzch), caplam(KXA,miso),
     >	ela(3,3,miso), zi_api(miso,nzch)
	m = 1		!ONLY FRICTION FORCE
	do misa = 1,miso
	nzmax = natom(misa)
	  do nz = 1,nzmax
	  friction(misa,nz) =
     >	  zi_api(misa,nz) * ( ela(m,1,misa)*usol(1,nz,misa)
     >	  + ela(m,2,misa)*usol(2,nz,misa)
     >	  + ela(m,3,misa)*usol(3,nz,misa) + caplam(m,misa) )
            fricc(misa,nz,1) = zi_api(misa,nz)*ela(m,1,misa)*usol(1,nz,misa)
            fricc(misa,nz,2) = zi_api(misa,nz)*ela(m,2,misa)*usol(2,nz,misa)
            fricc(misa,nz,3) = zi_api(misa,nz)*ela(m,3,misa)*usol(3,nz,misa)
            fricc(misa,nz,4) = zi_api(misa,nz)*caplam(m,misa)
	  friction(misa,nz) = friction(misa,nz)
     >	  - denmass(misa,nz)*usol(m,nz,misa)*
     >	  (nuion(misa,nz) + nurec(misa,nz))*al32(m)
            fricc(misa,nz,5) = - denmass(misa,nz)*usol(m,nz,misa)*
     >	    (nuion(misa,nz) + nurec(misa,nz))*al32(m)
	  if( nz.gt.1 ) then
      	    friction(misa,nz) = friction(misa,nz) + usol(m,nz-1,misa)
     >	    *denmass(misa,nz-1) * nuion(misa,nz-1) * al32(m)
            fricc(misa,nz,5) = fricc(misa,nz,5) + usol(m,nz-1,misa)
     >	    *denmass(misa,nz-1) * nuion(misa,nz-1) * al32(m)
          endif
	  if( nz.lt.nzmax ) then
      	    friction(misa,nz) = friction(misa,nz) + usol(m,nz+1,misa)
     >	    *denmass(misa,nz+1) * nurec(misa,nz+1) * al32(m)
            fricc(misa,nz,5) = fricc(misa,nz,5) + usol(m,nz+1,misa)
     >	    *denmass(misa,nz+1) * nurec(misa,nz+1) * al32(m)
          endif
	  enddo
	enddo
	return
	end
c---- End of subroutine getfrict ---------------------------------------
c-----------------------------------------------------------------------
	subroutine inicon
c*****	******************************************************************
c*****	INICON loads the common block containing physical constants and  *
c*****	conversion constants.                                           *
c*****	Last Revision:12/85 W.A.Houlberg and S.E.Attenberger ORNL.       *
c*****	******************************************************************
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
c*****	Physical constants.
c*****	Elementary charge - (coulomb).
	coulom=1.6022e-19
c*****	Permittivity of free space - (farad/meter).
	epsilo=8.8419e-12
c*****	Proton mass - (kilogram).
	promas=1.6726e-27
c*****	Unit conversions.
c*****	Joules per keV.
	xj7kv=1.6022e-16
c****	Constants
	one = 1.0e0
c*****	Calculate pi0 for greatest accuracy.
	pi0=acos(-one)
c*****	Set zero for initializi_aping arrays.
	zero=0.e0
c*****	Integer constants
	ilam1 = 1
	ilam2 = 2
	ilam3 = 3
	iacci = 4
	iforc = 5
	
	al32(1) = 1.
	al32(2) = 5./2.
	al32(3) = 35./8.
	return
	end
c---- End of subroutine inicon -----------------------------------------
c-----------------------------------------------------------------------
        subroutine initmombal(misotope,natomic,nchstate)
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
	integer natomic(misotope)

	mise = 1		!Electrons MUST be 1st element

c****	Load common block elements
	miso = misotope
	nzch = nchstate
	if( miso.gt.MXMISO ) call xerrab('MISO > MXMISO')
	if( nzch.gt.MXNZCH ) call xerrab('NZCH > MXNZCH')
	do misa = 1,miso
	   natom(misa) = natomic(misa)
        enddo

	return
	end
c---- End of subroutine initmombal -------------------------------------
c-----------------------------------------------------------------------
	subroutine mrespond(elab,zi_api,den,denmass,uresp,sbar,
     >	nmat,amat,nurec,umass,ldir)
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
        integer nmat
	real elab(KXA,miso,KXA*miso),uresp(KXA,nzch,miso,NBA),
     >	sbar(KXA,miso+1), sbar0(KXA*(MXMISO1)), ubar0(KXA*(MXMISO1)),
     >	zi_api(miso,*),amat(nmat,nmat),denmass(miso,0:nzch)
	real den(miso,0:nzch), nurec(miso,nzch)
	real usum(MXNZCH), bmat(KXA*(MXMISO1)*KXA*(MXMISO1))
	integer nrow(KXA*(MXMISO1)), z1
        integer i, j
	save nrow, bmat, sbar0, ubar0
	macc = miso+1	!Space for acci = acceleration force
	mflow = 1	!Flow
	z1 = 1
	fanorm = one / (totmass * anorm)	!Norm for AI force
c
c	COMPUTE ZI,DENMASS-WEIGHTED AVERAGES OF RESPONSE FUNCTIONS
c	AND FORM INTER-ISOTOPIC COUPLING MATRIX
c	NOTE: FOR MISA = 1 (ELECTRONS), M = 1 (FLOW), REPLACE
c	WITH NO NET FLOW CONDITION SUM M(ALPHA)*DEN(ALPHA,Z)U(ALPHA,Z) = 0
c	TO AVOID SINGULAR MATRIX (IN REF. FRAME OF MASS FLOW)
c
c	AUGMENT (FOR MISA=MACC, K=1) SYSTEM WITH EQUATIONS FOR MEAN
c	ACCELERATION TERM AI [=ubar(1,macc)]
c
c	RHO*AI + SUM( M N NUrec U(1,Z=1,alpha) ) = NEUTRAL SOURCE
c	AND U(1,1,alpha) MUST BE EXPRESSED IN TERMS OF AI, UBAR
c
c	JPAR = SUM( DENZ * U(1,Z,alpha) )
c
c
c	NOTE THAT U(m,Z,alpha) = -BIGLAM(m,Z,alpha;mp,beta)*ubar(mp,beta)
c	WHERE BIGLAM = SUM(over kp) elab(kp,alpha,mp,beta)*uresp(m,Z,alpha;kp)
c
        j = KXA*(macc-1)
	do mmisb = 1,nmat
	  amat(mflow,mmisb) = zero   # electron "isotope"
	  amat(1+j,mmisb) = zero
	  amat(2+j,mmisb) = zero
	  amat(3+j,mmisb) = zero
	enddo
	sbar(mflow,mise) = zero
	sbar(mflow,macc) = sumforce  * fanorm
	sbar(2,macc) = zero
	sbar(3,macc) = zero
	amat(1+j,1+j) = one   	!totmass*anorm * fanorm
	amat(2+j,2+j) = one	!NOT USED FOR ANYTHING PHYSICAL
	amat(3+j,3+j) = one	!NOT USED FOR ANYTHING PHYSICAL
	do 10 misa = 1,miso
	  nz = natom(misa)
	  dnur1 = (fanorm*denmass(misa,z1))*nurec(misa,z1)
c
c	COMPUTE ACCELERATION SOURCE, MATRIX ELEMENTS
c
	  sbar(1,macc) = sbar(1,macc) - dnur1*
     >	  uresp(mflow,z1,misa,iforc)
	  amat(1+j,1+j) = amat(1+j,1+j) + dnur1*
     >	  uresp(mflow,z1,misa,iacci)
	  do nzsum = 1,nz
c
c     	REPLACE ELECTRON MOMENTUM WITH NO NET FLOW CONDITION
c     	(IN REFERENCE FRAME MOVING WITH MASS FLOW)
c
	    rhomass = denmass(misa,nzsum)/totmass
	    sbar(mflow,mise) = sbar(mflow,mise) +
     >	    rhomass*uresp(mflow,nzsum,misa,iforc)
	    amat(mflow,1+j) = amat(mflow,1+j) -
     >	    rhomass*uresp(mflow,nzsum,misa,iacci)
	  enddo
          i = KXA*(misa-1)
	  do 10 m = 1,KXA
	    if( (misa.ne.mise) .or. (m.ne.mflow) )then
	      sbar(m,misa) =
     >	      sdot(nz,zi_api(misa,1),miso,uresp(m,z1,misa,iforc),KXA)
	      amat(m+i,1+j) =
     >	      -sdot(nz,zi_api(misa,1),miso,uresp(m,z1,misa,iacci),KXA)
	    endif
	  do 10 mpmisb = 1,KXA*miso
	    do nzb = 1,nz
	    usum(nzb) = uresp(m,nzb,misa,ilam1)*elab(ilam1,misa,mpmisb)
     >	              + uresp(m,nzb,misa,ilam2)*elab(ilam2,misa,mpmisb)
     >	              + uresp(m,nzb,misa,ilam3)*elab(ilam3,misa,mpmisb)
	    enddo
	    if( (misa.ne.mise) .or. (m.ne.mflow) )
     >	    amat(m+i,mpmisb) = sdot(nz,zi_api(misa,z1),miso,usum,1)
c
c	PERFORM SUMS OVER MISA (ISOTOPE MASS A)
c
	    if( m.eq.mflow )then
	      amat(1+j,mpmisb) = amat(1+j,mpmisb)
     >	      - dnur1*usum(z1)
	      do nzsum = 1,nz
	        amat(mflow,mpmisb) = amat(mflow,mpmisb)
     >	        + (denmass(misa,nzsum)/totmass)*usum(nzsum)
	      enddo
	    endif
 10	continue
c
c	ADD DIAGONAL ELEMENTS
c
	mmin = 1		!Ignore electron momentum balance here
	do 20 mmisa = mmin+1,KXA*miso
	  amat(mmisa,mmisa) = amat(mmisa,mmisa) + one
 20	continue
c
c	INVERT MATRIX TO OBTAIN REDUCED FLOWS
c
	if( ldir.gt.1 )then
	  call sgefa(amat,nmat,nmat,nrow,info)
	  nmat2 = nmat*nmat
	  call scopy(nmat2,amat,1,bmat,1)
	  call scopy(nmat,sbar,1,sbar0,1)
	  if( info.ne.0 )then
	    call xerrab('mrespond:  Condition No. = 0 in Solver')
	  endif
	else
	  do i = 1,nmat
	  sbar(i,1) = sbar(i,1) + sbar0(i) -
     >	      sdot(nmat,amat(i,1),nmat,ubar0,1)
	  enddo
	endif
        call sgesl(bmat,nmat,nmat,nrow,sbar,0)
	if( ldir.gt.1 )call scopy(nmat,sbar,1,ubar0,1)
	return
	end
c---- End of subroutine mrespond ---------------------------------------
c-----------------------------------------------------------------------
	subroutine mzrespond(elab,uresp,ubar,caplam,caplams,usol,
     >	usave,den,denz,parcurrent,qcond,ucond,tempa,umass,ldir)
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
	real elab(KXA*miso,KXA*miso),uresp(KXA*nzch,miso,NBA),
     >	ubar(KXA,miso+1), usol(KXA*nzch,miso), denz(miso,nzch),
     >	den(miso,0:nzch), qcond(miso,nzch), ucond(miso,nzch)
	real caplam(KXA,miso), caplams(KXA*miso), usave(*), tempa(miso)
	integer m,misa
	acci1 = ubar(1,miso+1)
	parcurrent = zero
	uscale = one	!/csnorm
	qscale = one	!/(csnorm * tempnorm * dennorm)
	mtotal = KXA*miso
	do mmisa = 1,mtotal
          misa = 1 + (mmisa-1)/KXA
          m = mmisa - (misa-1)*KXA
	  caplam(m,misa) =
     >	  sdot(mtotal,ubar,1,elab(mmisa,1),mtotal)
	enddo
	
	do misa = 1,miso
	  do knz = 1,KXA*natom(misa)
	    usol(knz,misa) = uresp(knz,misa,iforc)
     >  	- uresp(knz,misa,1)*caplam(1,misa)
     >  	- uresp(knz,misa,2)*caplam(2,misa)
     >  	- uresp(knz,misa,3)*caplam(3,misa)
     >  	+ uresp(knz,misa,iacci)*acci1
	  enddo
c
c	CONVERT BACK TO ABSOLUTE (NOT RELATIVE TO UMASS) FLOW
c	CONVERT TO RANDOM HEAT FLUX (MULT BY -2.5) FOR K=2
c	AND STORE NORMED FORM FOR EXTERNAL USE
c
	  k1 = 1
	  do nz = 1,natom(misa)
	    k2 = k1+1
	    ucond(misa,nz) = (usol(k1,misa) + umass) * uscale
	    qcond(misa,nz) = -2.5*usol(k2,misa)*den(misa,nz)
     >	                *tempa(misa) * qscale
	    parcurrent = parcurrent + denz(misa,nz)*ucond(misa,nz)
	    k1 = k1 + KXA
	  enddo
	enddo
c
c	STORE SOLUTION FOR INCREMENTAL MODE
c
	if( ldir.gt.1 )then
	mtot2 = mtotal*nzch
	call scopy(mtot2,usol,1,usave,1)
	call scopy(mtotal,caplam,1,caplams,1)
	acci0 = acci1
	endif
	return
	end
c---- End of subroutine mzrespond --------------------------------------
c-----------------------------------------------------------------------
      subroutine neolab(mntau,capm,capn,ela,elab)
c*****	******************************************************************
c*****	NEOLAB calculates the neoclassical friction coefficients.        *
c*****	References:                                                      *
c*****	Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079 (H&S).              *
c*****	Last Revision: 5/94 W.A.Houlberg and S.P.Hirshman ORNL.          *
c*****	Input:                                                           *
c*****	miso-number of isotopes.                                         *
c*****	mntau(misb)-ma<<na/tauab>>-kg/m**3/s.                           *
c*****	     -Eqn 5.25 of H&S.                                    *
c*****	capm(3,miso,3,miso)-matrix of test particle elements.            *
c*****	capn(3,miso,3,miso)-matrix of field particle elements.           *
c*****	Output:                                                          *
c*****	ela(3,3,misa)-test particle component-kg/m**3/s.                 *
c*****	elab(3,misa,3,misb)-field particle component-kg/m**3/s.           *
c*****	******************************************************************
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
      real      mntau(miso,miso)
      real      ela(KXA,KXA,miso), elab(KXA,miso,KXA,miso)
      real      capm(KXA,miso,KXA,miso), capn(KXA,miso,KXA,miso)

c ... Initialize ela array.
      do i=1,KXA
        do j=1,KXA
          do k=1,miso
            ela(i,j,k)=zero
          enddo
        enddo
      enddo

      do 40 misb=1,miso
        do 30 misa=1,miso
          do 20 k=1,KXA
	    do 10 j=1,KXA
c*****	Sum over isotopes misb for test particle component.
              ela(j,k,misa)=ela(j,k,misa)
     >                  +mntau(misa,misb)*capm(j,misa,k,misb)
c*****	Field particle component.
              elab(j,misa,k,misb)=mntau(misa,misb)*capn(j,misa,k,misb)
 10	      continue
 20	    continue
 30	  continue
 40	continue
      return
      end
c---- End of subroutine neolab -----------------------------------------
c-----------------------------------------------------------------------
      subroutine neomn(amu,capm,capn,tempa)
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
c*****	******************************************************************
c*****	NEOMN calculates the test particle (M) and field particle (N)    *
c*****	matrix elements of the collision operator using the Laguerre    *
c*****	polynomials of order 3/2 as basis functions for each isotopic   *
c*****	species combination.                                            *
c*****	References:                                                      *
c*****	Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079 (H&S).              *
c*****	Hirshman, Phys Fluids 20 (1977) 589.                             *
c*****	Last Revision: 3/94 W.A.Houlberg and S.E.Attenberger ORNL.       *
c*****	Input:                                                           *
c*****	miso-number of isotopes.                                         *
c*****	amu(miso)-atomic mass of isotopes-arbitrary.                     *
c*****	tempa(miso)-temperature of isotopes-arbitrary.                      *
c*****	Output:                                                          *
c*****	capm(3,miso,3,miso)-matrix of test particle elements.            *
c*****	capn(3,miso,3,miso)-matrix of field particle elements.           *
c*****	Comments:                                                        *
c*****	The indices on the matrices are one greater than the notation in *
c*****	the review article so as to avoid 0 as an index.                *
c*****	******************************************************************
      real      capm(KXA,miso,KXA,miso),  capn(KXA,miso,KXA,miso)
	real amu(miso), tempa(miso)
c*****	Loop over isotope a.
      do 20 misa=1,miso
c*****	Loop over isotope b.
        do 10 misb=1,miso
c*****	Ratio of masses.
          xmab=amu(misa)/amu(misb)
c*****	Ratio of temperatures.
          xtab=tempa(misa)/tempa(misb)
c*****	Ratio of thermal velocities, vtb/vta.
          xab=sqrt(xmab/xtab)
c*****	Elements of M.
          xab2=xab**2
          yab32=(one+xab2)*sqrt(one+xab2)
          yab52=(one+xab2)*yab32
          yab72=(one+xab2)*yab52
          yab92=(one+xab2)*yab72
c*****	Eqn 4.11 for M00 H&S.
          capm(1,misa,1,misb)=-(one+xmab)/yab32
c*****	Eqn 4.12 for M01 H&S.
          capm(1,misa,2,misb)=-(3.0/2.0)*(one+xmab)/yab52
c*****	Eqn 4.15 for M02 H&S.
          capm(1,misa,3,misb)=-(15.0/8.0)*(one+xmab)/yab72
c*****	Eqn 4.8 for M10 H&S.
          capm(2,misa,1,misb)=capm(1,misa,2,misb)
c*****	Eqn 4.13 for M11 H&S.
          capm(2,misa,2,misb)=-(13.0/4.0+4.0*xab2+(15.0/2.0)*xab2**2)
     >                         /yab52
c*****	Eqn 4.16 for M12 H&S.
          capm(2,misa,3,misb)=-(69.0/16.0+6.0*xab2+(63.0/4.0)*xab2**2)
     >                         /yab72
c*****	Eqn 4.8 for M20 H&S.
          capm(3,misa,1,misb)=capm(1,misa,3,misb)
c*****	Eqn 4.8 for M21 H&S.
          capm(3,misa,2,misb)=capm(2,misa,3,misb)
c*****	Eqn 5.21 for M22 H&S.
          capm(3,misa,3,misb)=-(433.0/64.0+17.0*xab2+(459.0/8.0)*xab2**2
     >                        +28.0*xab2**3+(175.0/8.0)*xab2**4)/yab92
c*****	Elements of N.
c*****	Momentum conservation, Eqn 4.11 for N00 H&S.
          capn(1,misa,1,misb)=-capm(1,misa,1,misb)
c*****	Eqn 4.9 and 4.12 for N01 H&S.
          capn(1,misa,2,misb)=-xab2*capm(1,misa,2,misb)
c*****	Eqn 4.15 for N02 H&S - corrected by multiplying rhs by (Ta/Tb)
          capn(1,misa,3,misb)=-xab2**2*capm(1,misa,3,misb)
c*****	Momentum conservation, Eqn 4.12 for N10 H&S.
          capn(2,misa,1,misb)=-capm(2,misa,1,misb)
c*****	Eqn 4.14 for N11 H&S.
          capn(2,misa,2,misb)=(27.0/4.0)*xtab*xab2/yab52
c*****	Eqn 4.17 for N12 H&S.
          capn(2,misa,3,misb)=(225.0/16.0)*xtab*xab2**2/yab72
c*****	Momentum conservation for N20.
          capn(3,misa,1,misb)=-capm(3,misa,1,misb)
c*****	Eqn 4.9 and 4.17 for N21 H&S.
          capn(3,misa,2,misb)=(225.0/16.0)*xab2/yab72
c*****	Eqn 5.22 for N22 H&S.
          capn(3,misa,3,misb)=(2625.0/64.0)*xtab*xab2**2/yab92
   10   continue
   20 continue
      return
      end
c---- End of subroutine neomn ------------------------------------------
c-----------------------------------------------------------------------
	subroutine printit(amu,den,epar,zi_api,ela,
     >  caplam,sbar,usol,fmom,denmass,nuion,nurec,
     >	gradp,gradt,parcurrent,xirl,tempa,uneut,umass)
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
	real sbar(KXA,miso+1), usol(KXA,nzch,miso)
	real den(miso,0:nzch), caplam(KXA,miso)
	real denmass(miso,0:nzch)
	real nuion(miso,0:nzch), nurec(miso,nzch),
     >	gradp(miso,*), gradt(miso,*), ela(3,3,miso),
     >	       zi_api(miso,nzch), fmom(KXA,nzch,miso), amu(miso),
     >	       tempa(miso), uneut(miso)
	
c*****	Print out average-ion quantities
	write(*,490)xirl,umass
 490	format(' XI(IGRID) = ',1pe10.3,' UMASS = ',1pe10.3,/,
     >	' NOTE: RESULTS ARE PHYSICALLY VALID IF |UBAR| < VTHERM AND',
     >	' |QBAR/P| < VTHERM')
	acci = acci0*anorm
	sumflow = 0.
	sumaflow= 0.
	do 55 misa = 1,miso
	  vtherm=sqrt(2.0*tempa(misa)/amu(misa)/promas)
	do 55 m = 1,2		!Only look at flow, heat flow
	if( m.eq.1 )then
	write(*,500)m,amu(misa),sbar(m,misa),
     >	1000.*tempa(misa)/xj7kv,vtherm,caplam(m,misa)
	  write(*,510)
	  if( misa.ne.1 )
     >	  write(*,520)0,uneut(misa),den(misa,0),nuion(misa,0)
	else if(m.eq.2)then
	  write(*,505)m,amu(misa),sbar(m,misa)
	  write(*,530)
	endif
	nzmax = natom(misa)
	do 55 nz = 1,nzmax
	force = fmom(m,nz,misa)
	friction = zi_api(misa,nz) * ( ela(m,1,misa)*usol(1,nz,misa)
     >	+ ela(m,2,misa)*usol(2,nz,misa)
     >	+ ela(m,3,misa)*usol(3,nz,misa) + caplam(m,misa) )
	friction = friction - denmass(misa,nz)*usol(m,nz,misa)*
     >	(nuion(misa,nz) + nurec(misa,nz))*al32(m)
	if( nz.gt.1 )friction = friction + usol(m,nz-1,misa)
     >	*denmass(misa,nz-1) * nuion(misa,nz-1) * al32(m)
	if( nz.lt.nzmax )friction = friction + usol(m,nz+1,misa)
     >	*denmass(misa,nz+1) * nurec(misa,nz+1) * al32(m)
	if( m.eq.1 )then
	  force = force + denmass(misa,nz)*acci
	  sumflow = sumflow + denmass(misa,nz)*usol(1,nz,misa)
	  sumaflow= sumaflow+ denmass(misa,nz)*abs(usol(1,nz,misa))
	  write(*,550)nz,usol(m,nz,misa),force,friction,
     >	  den(misa,nz),zi_api(misa,nz),nuion(misa,nz),nurec(misa,nz)
	else if( m.eq.2 )then
	  pres = den(misa,nz)*tempa(misa)
	  write(*,550)nz,-2.5*usol(m,nz,misa)*pres,
     >	  force,friction,gradp(misa,nz),den(misa,nz)*gradt(misa,nz)
	endif
 55	continue
	sumaflow = max(sumaflow,totmass*abs(umass))
 500	format(/' UBAR (m=',i1,', mass=',1pe9.2,') = ',1pe10.3,
     >	' (m/sec)  TEMP = ',1pe10.3,' (eV)  VTHERM = ',1pe10.3,
     >	' (m/sec)','  LAMBDA =',1pe10.3)
 505	format(/' QBAR/P(m=',i1,', mass=',1pe9.2,') = ',1pe10.3,
     >	' (m/sec)')
 510	format('  Z  U||-UMASS',
     >	'     FORCE    FRICTION',
     >	'       DEN         ZI     NU-ION     NU-REC')
 520	format(i3,1pe11.3,2(11x),1pe11.3,11x,1pe11.3)
 530	format('  Z      Q||',
     >	'       FORCE    FRICTION',
     >	'   GRAD||P   N*GRAD||T')
 550	format(i3,1p8e11.3)
	write(*,560)sumflow/sumaflow,parcurrent,
     >	totmass*acci,epar
 560	format(/' NORMALIZED NET MASS FLOW, SUM(RHO*U)/SUM(RHO*|U|)',
     >	' -> 0 (U = U||-UMASS) = ',1pe11.3,
     >	/' J||(calc) = ',1pe11.3,/,
     >	' MEAN INERTIAL TERM RHO*AI   = ',1pe11.3,/,
     >	' E||     = ',1pe11.3/)
	return	
	end
c---- End of subroutine printit ----------------------------------------
c-----------------------------------------------------------------------
	subroutine setden(amu,den,denmass,denz,denz2,zi_api)
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
	real amu(miso), denz2(miso)
	real den(miso,0:nzch), denmass(miso,0:nzch), denz(miso,*)
	real zi_api(miso,*)
c*****	Sum over isotopes for electron density-isotope 1.
	do misa = 2,miso
	  do nza = 1,natom(misa)
	    denz(misa,nza) = den(misa,nza)*real(nza)*coulom
	  enddo
	enddo
	denz(1,1) = -den(1,1)*coulom
c*****	Compute mass-density, ZI, SUM(ZI), and total mass density
	totmass = zero
	do misa = 1,miso
	  bmass = amu(misa)*promas
	  denz2(misa)=zero
	  do nza = 1,natom(misa)
	    denmass(misa,nza) = bmass*den(misa,nza)
	    zi_api(misa,nza) = den(misa,nza)*real(nza)**2
	    denz2(misa)=denz2(misa) + zi_api(misa,nza)
	    totmass = totmass + denmass(misa,nza)
	  enddo
	enddo
c*****	Compute neutral mass-density
	do misa = 2,miso
	  denmass(misa,0) = amu(misa)*promas*den(misa,0)
	enddo
c*****	Provide pedestal for zi_api to avoid matrix inversion singularity
	zi_apiped = 1.E-4
	do misa=1,miso
	  zisum = zero
	  do nza=1,natom(misa)
	    zi_api(misa,nza)=zi_api(misa,nza)/denz2(misa) + zi_apiped
	    zisum = zisum + zi_api(misa,nza)
	  enddo
	  do nza = 1,natom(misa)	!zi_api normed to 1
	    zi_api(misa,nza) = zi_api(misa,nza)/zisum
	  enddo
	enddo
	return
	end
c---- End of subroutine setden -----------------------------------------
c-----------------------------------------------------------------------
	subroutine setforce(den,denz,denmass,epar,gradp,gradt,
     >	tempa,uneut,qneut,umass,fmom,nuion)
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
	real gradp(miso,*), gradt(miso,*), nuion(miso,0:nzch),
     >	denmass(miso,0:nzch), den(miso,0:nzch), denz(miso,nzch),
     >	fmom(KXA,nzch,miso), tempa(miso), uneut(miso), qneut(miso)
c
c	COMPUTE FORCES PROPORTIONAL TO GRAD(P), GRAD(T)
c
c	FMOMENTA(1,ALPHA) = GRAD(P sub alpha) -e*N*E|| - M N NU(ion)( UN - U ) delta(Z,1)
c	FMOMENTA(2,ALPHA) =-2.5 N GRAD T
c	FMOMENTA(3,ALPHA) = 0
c
	sumforce = zero
	do misa = 1,miso
	  do nz = 1,natom(misa)
	    fmom(1,nz,misa) = gradp(misa,nz) - denz(misa,nz)*epar
	    fmom(2,nz,misa) = -2.5*den(misa,nz)*gradt(misa,nz)
	    fmom(3,nz,misa) = 0.
	    if( nz.eq.1 .and. misa.ne.mise )then
	    fmom(1,nz,misa)
     >	    = fmom(1,nz,misa) - al32(1)*denmass(misa,0)*	!(Neutral mass-density)
     >	    nuion(misa,0)*(uneut(misa) - umass)
	    presa = den(misa,nz)*tempa(misa)
	    fmom(2,nz,misa)
     >	    = fmom(2,nz,misa) - al32(2)*denmass(misa,0)*	!(Neutral mass-density)
     >	    nuion(misa,0)*(-2.*qneut(misa)/(5.*presa))
	   endif
	    sumforce = sumforce - fmom(1,nz,misa)
	  enddo
	enddo
	return
	end
c---- End of subroutine setforce ---------------------------------------
c-----------------------------------------------------------------------
	subroutine uinvm2(k,n,a,b,x,nrow,scas,lss,msbs,
     .	nbrt,istore,iflag)
c	ALSO, BE SURE TO LINK IN BLAS LIBRARY (-LBLAS)
	implicit real (a-h,o-z), integer(i-n)
c*****	******************************************************************
c*****	UINVM2 performs the matrix inversion for a system of k equations *
c*****	each having n points.                                           *
c*****	Previous Last Revision: 3/90 W.A.Houlberg ORNL.                  *
c*****	Last Revision: 5/94 S.P.Hirshman ORNL (multiple right sides)     *
c*****	        and BLAS routines                                 *
c*****	Input:                                                           *
c*****	k-number of equations per block.                                 *
c*****	n-number of nodes in each equation (number of blocks)            *
c*****	nbrt-number of right hand sides (independent source terms)       *
c*****	a(n*k*(3*k))-coefficient array.                                  *
c*****	b(n*k,nbrt)-right sides (sources) loaded as 1-D vector b(k,n,nbrt)   *
c*****	loaded as one dimensional x(n,k,nbrt).                           *
c*****	iflag: if non-zero, routine is being called with previously
c*****	istore: if non-zero, store for restart scas,lss,msb vectors
c*****	       computed values of nrow, a
c*****	       use iflag != 0 for incremental update
c*****	Output:                                                          *
c*****	x(n*k,nbrt)-solution array loaded as one dimensional x(k,n,nbrt).*
c*****	nrow(i)-record of row permutations for i-th row of matrix a      *
c*****	 during the gaussian elimination with partial pivoting.   *
c*****	 dimension in calling routine must be >= n*k              *
c*****	iflag-error flag for solution.                                   *
c*****	=0 normal solution.                                         *
c*****	=1 singular matrix found in inverting a.                    *
c*****	Comments:                                                        *
c*****	We wish to approximate the solution of the system of linear     *
c*****	equations fx=b where f is a block tridiagonal matrix of order n *
c*****	and each block is of order k. thus there n*k equations in the   *
c*****	system. Let f(i,j) represent the (i,j) block of the matrix f.   *
c*****	Let b(i) denote the corresponding k by 1 block of the vector b. *
c*****	Then the 3k X n*k array a contains on entry into UINVM the   *
c*****	nontrivial information to describe the coefficient matrix f.    *
c*****	The matrix "A"           is described more easily by displaying *
c*****	the transpose of "A" which follows (note: "A" is stacked            *
c*****	so that the diagonal elements are aligned in the same column)   *
c*****	0 block    f(1,1)      f(1,2)   b(1)   1st k eqns..             *
c*****	f(2,1)     f(2,2)     f(2,3)    b(2)   next k eqns..            *
c*****	...        ...        ...      ...    ...                       *
c*****	f(n-1,n-2) f(n-1,n-1) f(n-1,n)  b(n-1) next to last k eqns..    *
c*****	f(n,n-1)   f(n,n)     0 block   b(n)   last k eqns..            *
*	INDEXING OF MATRIX A
*	F(I,J,M,N)*X(N,J) = B(M,I)    ,   AMAT(L) = LMAT(I,J,M,N) strung out
*
*	HERE, 1 <= (M,N) <= J,  1 <= I <= N,  I-1 <= J <= I+1
*
*	THUS, A(L) = LMAT(I,J,M,N) FOR
*
*	L = N + (M-1)*KX3 + (J-I+1)*K + (I-1)*KSQ3
*
*	WHERE KX3 = 3*K (OR 3*K+1 for original version), KSQ3 = K * KX3
*
c*****	Partial pivoting is within a subblock of k rows.                *
c*****	For full partial pivoting we would have to search over 2k rows  *
c*****	for the pivot element.                                          *
c*****	Double precision may be required for s,sca,anor,raia,raig,x,a.  *
c*****	******************************************************************
	real a(3*n*k*k), b(n*k,nbrt), x(n*k,nbrt), scas(*)
	integer lss(*), msbs(*)
	integer nrow(n*k)
	nk=n*k
	km1=k-1
	k2m1=k+km1
	k2=k+k
	k3=k2+k
c*****	Previously computed a-matrix inversion
	if( iflag.ne.0 )goto 1000
c*****	Shift first row blocks 1 block to left (on input, 1st block = 0)
	na=0
	do 30 j=1,k
	  napk=na+k
	  call scopy(k2,a(napk+1),1,a(na+1),1)
	  napk2=na+k2
	  do i=1,k
	    a(napk2+i)=0.0
	  enddo
	  na=na+k3
 30	continue
	na=0
	do i=1,nk
	  nrow(i)=na
	  na=na + k3
	enddo
 1000	continue
c*****	Begin loop through n blocks.
c*****	nl=(iblock-1)*k+1,  nh=(iblock+1)*k  for iblock=1,...,n-1.
	nl=1
	icount = 0
	do 50 iblock = 1,n
	  nh=nl+k2m1
	  nhs=nl+km1
	  kblock = k
	  krow = k3
	  if( iblock.eq.n )then
	    nh = nk
	    nhs= nk
	    kblock = km1
	    krow = k
	  endif
	  do 90 j=1,kblock
c*****	Determine element of maximum modulus in next column.
c*****	Only allow permutations within blocks of k rows.
	    if( iflag.ne.0 )goto 1010
	    s=0.0
	    do l=nl,nhs
	      sa=a(nrow(l)+j)
	      if(sa.lt.0.0) sa=-sa
	      if(s.lt.sa) then
	        ni=l   	    	    !Pivot row
	        s=sa
	      endif
	    enddo
	    if(s.eq.0.0) then
c*****	Singular matrix found.
	      iflag=1
	      return
	    endif
c*****	Row exchange (pivot)
	    ms=nrow(ni)
	    nrow(ni)=nrow(nl)
	    nrow(nl)=ms
	    anor=a(ms+j)
	    jp1=j+1
 1010	    continue
	    nl=nl+1
c*****	Gaussian elimination within block
	    do l=nl,nh
	      icount = icount+1
	      if( iflag.eq.0 )then
	        sca=a(nrow(l)+j)
	        if(sca.ne.0.0) then
	          sca=sca/anor
c*****	  Subtract pivot row elements (ms) from row (l) elements
	          call saxpy(krow-jp1+1,-sca,a(ms+jp1),1,
     >	          a(nrow(l)+jp1),1)
	        endif
	        ls = nrow(l)/k3 + 1
	        msb= ms/k3 + 1
	        if( istore.ne.0 )then
	          lss(icount) = ls
	          msbs(icount) = msb
	          scas(icount) = sca
	        endif
	      else
	        ls = lss(icount)
	        msb= msbs(icount)
	        sca = scas(icount)
	      endif
c*****	  Subtract pivot row right side from other right-hand side(s)
	      if(sca.ne.0.0)call saxpy(nbrt,-sca,b(msb,1),nk,b(ls,1),nk)
	    enddo
 90	  continue
c*****	Move block 2 to 1, block 3 to 2, 0 out block 3 of rows used in
c*****	this stage that will also be used in next.
	  if( iblock.ne.n .and. iflag.eq.0 )then
	    do l=nl,nh
	      call scopy(k2,a(nrow(l)+k+1),1,a(nrow(l)+1),1)
	      do j=1,k
	        a(nrow(l)+j+k2)=0.0
	      enddo
	    enddo
	  endif
	  if(nh.gt.nk) call xerrab('uinvm2:  nh > nk')
 50	continue
	if( icount .gt. 2*k*k*n ) call xerrab('uinvm2:  icount > 2*k*k*n')
c*****	Process back substitution in blocks of k equations.  1st 2 blocks
c*****	are special, but have special coding only for the 1st block.
	if(a(nrow(nk)+k).eq.0.0) then
c*****	Singular matrix found.
	  iflag=1
	  return
	endif
	do 100 lb = 1,nbrt
	  nkb = nrow(nk)/k3 + 1
	  x(nk,lb)=b(nkb,lb)/a(nrow(nk)+k)
	  neq=nk-1
	  if(k.ne.1) then
	    ntrm=1
	    do lk=km1,1,-1
	      neqb = nrow(neq)/k3 + 1
	      s = b(neqb,lb) -
     >	      sdot(ntrm,x(neq+1,lb),1,a(nrow(neq)+lk+1),1)
	      x(neq,lb)=s/a(nrow(neq)+lk)
	      ntrm=ntrm+1
	      neq=neq-1
	    enddo
	  endif
	  if( neq.le.0 )goto 100
 180	  ntrm=k
	  do lk=k,1,-1
	    neqb = nrow(neq)/k3 + 1
	    s = b(neqb,lb) -
     >	    sdot(ntrm,x(neq+1,lb),1,a(nrow(neq)+lk+1),1)
	    x(neq,lb)=s/a(nrow(neq)+lk)
	    ntrm=ntrm+1
	    neq=neq-1
	  enddo
	  if(neq.gt.0) goto 180
 100	continue
	iflag=0
	return
	end
c---- End of subroutine uinvm2 -----------------------------------------
c-----------------------------------------------------------------------
	subroutine zrespond(den,denmass,zi_api,ela,nuion,nurec,
     >	uresp,usave,caplams,fmom,ldir)
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
	parameter( NSIZE = 2*KXA*KXA*MXNZCH )
	real zmat(KNX,MXMISO), source(KXA*MXNZCH,MXMISO)
	real amat(KNX), asource(KXA*MXNZCH*NBA)
	real den(miso,0:nzch), zi_api(miso,nzch), ela(KXA,KXA,miso)
	real nuion(miso,0:nzch), nurec(miso,nzch)
	real uresp(KXA,nzch,miso,NBA), usave(KXA,nzch,miso)
	real fmom(KXA,nzch,miso),denmass(miso,0:nzch)
	real caplams(KXA,miso), xsol(NBA*KXA*MXNZCH)
	real scaz(NSIZE,MXMISO)
	integer nrowz(KXA*MXNZCH,MXMISO), lsz(NSIZE,MXMISO),
     >	msz(NSIZE,MXMISO)
	save zmat, source, scaz, lsz, msz, nrowz
c
c	FOR LDIR > 1, DO COMPLETE NEW MATRIX INVERSION
c	FOR LDIR <=1, DO INCREMENTAL MATRIX INVERSION, USING
c
c	[ ZMAT(old) + del-ZMAT ][USOL(old) + del-U] = SOURCE + del-SOURCE
c
c	==>   ZMAT(old)*USOL(new) = S(new) + [S(old) - ZMAT(new)*USOL(old)]
c
c	HERE, ZMAT(new)*USOL(old) = ZMAT*USAVE - S(IACCI)*ACCI(old) -
c
c	 - S(1)*CAPLAM(1) - S(2)*CAPLAM(2) - S(3)*CAPLAM(3)
c
c
c	COMPUTE NORMS FOR AI
c	AI = ANORM * aI, WHERE
c	ANORM = SUM(ela(1,1,a))/SUM(DENMASS)
c
	elsum = zero
	do misa = 1,miso
	  elsum = elsum + ela(1,1,misa)
	enddo
	anorm = elsum / totmass
c
c	FOR EACH FIXED ISOTOPE, COMPUTE RESPONSE FUNCTIONS WHICH SOLVE
c	THE EQUATION (FOR ALPHA FIXED)
c
c	ZMAT(m,Z,mm,ZZ) u(mm,ZZ;k) = SOURCE(m,Z;k)
c
c	WHERE FOR
c
c	k = 1, 2, 3 [ uresp(m,Z;k) == ULAMk ]
c
c	SOURCE(m,Z;k)  = +ZI(alpha,Z) delta(m,k)
c
c	k = 4 [uresp(m,Z;k) == UACCI] SOURCE = DEN(alpha,Z)*MASS(alpha) delta(m,1)
c
c	k = 5 [uresp(m,Z;k) == UFORC] SOURCE = FORCE(alpha,m,Z)   (see paper...)
c
c	THEN, THE SOLUTION (IN TERMS OF THE INTER-SPECIES COUPLING TERMS
c	LAMBDA AND OTHER DRIVING TERMS) IS:
c
c	u(alpha,m,Z) =  SUM(k=1,2,3) [-ULAMk * LAMBDA(alpha,k) ]
c
c	          + UACCI * AI + UFORC
c
c	NOTE THAT THE ORDERING (STACKING) OF SOURCE AND URESP IS THE
c	SAME (SEE SUBROUTINE ZSOURCE FOR DETAILS)
c
	kx3 = 3*KXA
	kxsq3 = kx3 * KXA
	do 100 misa = 1,miso
	nz = natom(misa)
	nkz = nz * KXA
	nresp = NBA*nkz		!Number of response array elements per isotope
	nknz = kxsq3*nz
c
c	SET UP SOURCE TERMS: SOURCE(M,Z,k=1,NBA)
c	NOTE: Z = 1,NATOM IS DIFFERENT SIZE FOR EACH ISOTOPE
c

c ... Initalize asource array.
        do i=1,nresp
          asource(i)=zero
        enddo

	call zsource(asource,zi_api,den,fmom(1,1,misa),denmass,misa,nz)
	nforc = nkz*(iforc-1)
	nacci = nkz*(iacci-1)
	if( ldir.le.1 )then
	do i = 1,nkz
	asource(i+nforc) = asource(i+nforc) + ( source(i,misa) +
     >  asource(i+nacci)*acci0 -
     >	(asource(i)*caplams(1,misa) + asource(i+nkz)*caplams(2,misa)
     >	+asource(i+2*nkz)*caplams(3,misa)) )
	enddo
	endif
c
c	SET UP Z-MATRIX ARRAY ELEMENTS
c

c ... Initialize amat array.
        do i=1,nknz
          amat(i)=zero
        enddo

	do 10 m = 1,KXA
	do 10 mp= 1,KXA
	i1 = mp + kx3*(m-1)
	i2 = m + nforc
	do 10 jz= 1,nz
	index = i1 + (jz-1)*kxsq3	!consecutive block-by-block index
	is = i2 + (jz-1)*KXA
c
c	LOWER DIAGONAL ELEMENTS (IONIZATION FROM LOWER CHARGE STATE)
c
	ozi_api = al32(m)/zi_api(misa,jz)
	if (jz.gt.1)then
	  jp = jz-1
	  if( (m.eq.mp) )then
	  amat(index) = ozi_api*denmass(misa,jp)*nuion(misa,jp)
	  if( ldir.le.1 )
     >	  asource(is) = asource(is) - amat(index)*usave(mp,jp,misa)
	  endif
	endif
c
c	DIAGONAL ELEMENTS
c
	index = index + KXA
	amat(index) = ela(m,mp,misa)
	if( (m.eq.mp) )
     >	amat(index) = amat(index) - denmass(misa,jz)*
     >	( nuion(misa,jz) + nurec(misa,jz) )*ozi_api
	if( ldir.le.1 )
     >	asource(is) = asource(is) - amat(index)*usave(mp,jz,misa)
c
c	UPPER DIAGONAL ELEMENTS (RECOMBINATION FROM HIGHER CHARGE STATE)
c
	if (jz.lt.nz)then
	  index = index + KXA
	  jp = jz+1
	  if( (m.eq.mp) )then
	  amat(index) = ozi_api*denmass(misa,jp)*nurec(misa,jp)
	  if( ldir.le.1 )
     >	  asource(is) = asource(is) - amat(index)*usave(mp,jp,misa)
	  endif
	endif
 10	continue
c
c	SOLVE TRIDIAGONAL SYSTEM WITH NBA RIGHT SIDES STORED IN SOURCE
c	STORE FOR FAST (INCREMENTAL) RECALCULATION (LDIR = 0,1)
c
	noff = 1 + nforc
	istore = 1
	if( ldir.gt.1 )then
	  call scopy(nkz,asource(noff),1,source(1,misa),1)
	  lflag = 1
	  noff = 1
	  nr = NBA
	  iflag = 0
	else
	  lflag = iforc
	  nr = 1
	  iflag = 1
	  call scopy(nknz,zmat(1,misa),1,amat,1)
	endif
	call uinvm2(KXA,nz,amat,asource(noff),xsol(noff),
     >	nrowz(1,misa),scaz(1,misa),lsz(1,misa),msz(1,misa),
     >	nr,istore,iflag)
	if( ldir.gt.1 )call scopy(nknz,amat,1,zmat(1,misa),1)
	if( iflag.ne.0 )then
	call xerrab(" CALL TO UINVM2 FAILED! ")
	endif
c
c	STORE SOLUTION VECTOR IN RESPONSE ARRAY URESP FOR EACH ISOTOPE
c	URESP(KXA=1,2,3,NATOM,n=1,..,NBA).
c	Note response arrays are arranged consecutively
c
	do ntype = lflag,NBA
	noff = 1 + nkz*(ntype-1)
	call scopy(nkz,xsol(noff),1,uresp(1,1,misa,ntype),1)
	enddo
 100	continue
	return
	end
c---- End of subroutine zrespond ---------------------------------------
c-----------------------------------------------------------------------
	subroutine zsource(source,zi_api,den,fmom,denmass,misa,nz)
  	implicit real (a-h,o-z), integer (i-n)
	Use(Reduced_ion_constants)
	real zi_api(miso,*), den(miso,0:nz), source(KXA,nz,NBA), fmom(KXA,*)
	real denmass(miso,0:nz)
	do 20 jz= 1,nz
	do 20 m = 1,KXA
	if( m.eq.1 )then
	  source(m,jz,ilam1) = one
	  source(m,jz,iacci) = denmass(misa,jz) * anorm / zi_api(misa,jz)
	else if( (m.eq.ilam2) .or. (m.eq.ilam3) )then
	  source(m,jz,m) = one
	endif
	source(m,jz,iforc) = fmom(m,jz)/zi_api(misa,jz)
 20	continue
	return
	end
c---- End of subroutine zsource ----------------------------------------
c-----------------------------------------------------------------------
