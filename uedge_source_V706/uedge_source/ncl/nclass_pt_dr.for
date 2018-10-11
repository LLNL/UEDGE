      PROGRAM NCLASS_PT_DR
!***********************************************************************
!NCLASS_PT_DR is a stand-alone driver to call NCLASS for a single radial
!  point with simplified geometry
!NCLASS calculates the neoclassical quantities using the reduced charge
!  state formalism of Hirshman, the friction coefficients of Hirshman,
!  and the velocity dependent viscosities of Shaing
!References:
!  NCLASS equations:
!    Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  Multiple species formulation of neoclassical transport:
!    Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079
!  Friction coefficients, pitch angle diffusion:
!    Hirshman, Sigmar, Phys Fluids 19 (1976) 1532
!    Hirshman, Phys Fluids 20 (1977) 589
!  Reduced charge state formalism:
!    Boley, Gelbard, Hirshman, Phys Fluids 22 (1979) 1280
!  Orbit squeezing:
!    Shaing, Hazeltine, Phys Fluids B 4 (1992) 2547
!    Shaing, Hsu, Hazeltine, Phys Plasmas 1 (1994) 3365
!  Viscosity coefficients:
!    Shaing, Hsu, Yokoyama, Wakatani,Phys Plasmas 2 (1995) 349
!    Shaing, Yokoyama, Wakatani, Hsu, Phys Plasmas 3 (1996) 965
!  Trapped particle fraction:
!    Lin-Liu, Miller, Phys Plasmas 2 (1995) 1666
!  Other implementations of Hirshman formalism:
!    Tani, Azumi, Devoto, J. Comput Phys 98 (1992) 332
!    Kessel, Nucl Fusion 34 (1994) 1221
!  NCLASS applications:
!    Bush, et al, J Nucl Mater 241-243 (1997) 892
!    Efthimion, et al, Phys Plasmas 5 (1998) 1832
!    Houlberg, Baylor, Fusion Technol 34 (1998) 591
!    Efthimion, et al, 1998 Yokoyama IAEA to be published
!    Hillis, et al, Phys Plasmas 6 (1999) 1985
!    Wade, Houlberg, Baylor, Phys Rev Lett submitted
!  W.A.Houlberg 6/99
!***********************************************************************
      IMPLICIT NONE
      INCLUDE 'pamx_mi.inc'
      INCLUDE 'pamx_ms.inc'
      INCLUDE 'pamx_mz.inc'
!Declaration of input to NCLASS
      INTEGER        k_order,                 k_potato
      INTEGER        m_i,                     m_z
      REAL           c_den,                   c_potb,
     #               c_potl
      REAL           p_b2,                    p_bm2,
     #               p_eb,                    p_fhat,
     #               p_fm(3),                 p_ft,
     #               p_grbm2,                 p_grphi,
     #               p_gr2phi,                p_ngrth
      REAL           amu_i(mx_mi),            grt_i(mx_mi),
     #               temp_i(mx_mi)
      REAL           den_iz(mx_mi,mx_mz),     fex_iz(3,mx_mi,mx_mz),
     #               grp_iz(mx_mi,mx_mz)
!Declaration of output from NCLASS
      INTEGER        iflag,                   m_s
      INTEGER        jm_s(mx_ms),             jz_s(mx_ms)
      REAL           p_bsjb,                  p_etap,
     #               p_exjb
      REAL           calm_i(3,3,mx_mi)
      REAL           caln_ii(3,3,mx_mi,mx_mi),capm_ii(3,3,mx_mi,mx_mi),
     #               capn_ii(3,3,mx_mi,mx_mi)
      REAL           bsjbp_s(mx_ms),          bsjbt_s(mx_ms),
     #               dn_s(mx_ms),             gfl_s(5,mx_ms),
     #               qfl_s(5,mx_ms),          sqz_s(mx_ms),
     #               upar_s(3,3,mx_ms),       utheta_s(3,3,mx_ms),
     #               vn_s(mx_ms),             veb_s(mx_ms),
     #               qeb_s(mx_ms),            xi_s(mx_ms),
     #               ymu_s(3,3,mx_ms)
      REAL           chip_ss(mx_ms,mx_ms),    chit_ss(mx_ms,mx_ms),
     #               dp_ss(mx_ms,mx_ms),      dt_ss(mx_ms,mx_ms)
!Declaration of local variables
      CHARACTER      label*120
      INTEGER        i,                       im,
     #               iz,                      iza,
     #               j,                       jm,
     #               jza,                     k,
     #               k_out,                   l,
     #               nin,                     nout
      INTEGER        idum(8)
      REAL           a0,                      bt0,
     #               bpol,                    btor,
     #               btot,                    ds,
     #               e0,                      p_eps,
     #               p_q,                     ppr,
     #               q0,                      r0,
     #               uthai
      REAL           dq_s(mx_ms),             vq_s(mx_ms)
      REAL           z_coulomb,               z_electronmass,
     #               z_j7kv,                  z_mu0,
     #               z_pi,                    z_protonmass
      REAL           dum(8),                  edum(8),
     #               rdum(8)
!Declaration of functions
      REAL           RARRAY_SUM
      NAMELIST/indata/k_out,k_order,k_potato,
     #                c_den,p_eps,p_grphi,p_gr2phi,p_q,
     #                r0,a0,e0,bt0,q0,
     #                amu_i,temp_i,den_iz
!Physical and conversion constants
      z_coulomb=1.6022e-19
      z_electronmass=9.1095e-31
      z_j7kv=1.6022e-16
      z_mu0=1.2566e-06
      z_pi=ACOS(-1.0)
      z_protonmass=1.6726e-27
!Set default data
!  k_out-option for output to nout (-)
!       =1 errors only
!       =2 errors and results 
!       =else no output
      k_out=2
!  k_order-order of v moments to be solved (-)
!         =2 u and p_q
!         =3 u, p_q, and u2
!         =else error
      k_order=2
!  k_potato-option to include potato orbits (-)
!          =0 off
!          =else on
      k_potato=0
!  c_den-density cutoff below which species is ignored (/m**3)
      c_den=1.0e10
!  p_eps-inverse aspect ratio (-)
      p_eps=0.1
!  p_grphi-radial electric field Phi' (V/rho)
      p_grphi=-5.0e3
!  p_gr2phi-radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
      p_gr2phi=-2.0e4
!  p_q-safety factor (-)
      p_q=-3.0
!  r0-major radius (m)
      r0=3.0
!  a0-minor radius, scale length for rho (m)
      a0=1.0
!  e0-axial elongation (-)
      e0=1.0
!  bt0-axial toroidal field (T)
      bt0=-3.0
!  q0-axial safety factor (-)
      q0=-3.0
!Set the input unit and open file
      nin=1
      OPEN(unit=nin,
     #     file='in_nclass_pt.nml',
     #     status='old',
     #     access='sequential')
!Overwrite default data with namelist input and close file
      READ(nin,indata)
      CLOSE(unit=nin)
!Calculate input for NCLASS
!  m_i-number of isotopes (1<mi<mx_mi+1)
!  m_z-highest charge state of all species (0<mz<mx_mz+1)
      m_i=0
      m_z=0
	DO i=1,mx_mi
        ds=0.0
	  DO j=1,mx_mz
	    IF(den_iz(i,j).gt.c_den) THEN
            ds=den_iz(i,j)
	      IF(j.gt.m_z) m_z=j
	    ENDIF
	  ENDDO
	  IF(ds.gt.c_den) m_i=m_i+1
	ENDDO
!  grt_i(i)-temperature gradient of i (keV/rho)
!  grp_iz(i,z)-pressure gradient of i,z (keV/m**3/rho)
      DO i=1,m_i
        grt_i(i)=-temp_i(i)/a0/2.0
        DO j=1,m_z
          grp_iz(i,j)=-temp_i(i)*den_iz(i,j)/a0
        ENDDO   
      ENDDO
!  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
      c_potb=e0*bt0/2/q0**2
!  c_potl-q(0)*R(0) (m)
      c_potl=q0*r0
!  p_b2-<B**2> (T**2)
      p_b2=bt0**2*(1.0+0.5*p_eps**2)
!  p_bm2-<1/B**2> (/T**2)
      p_bm2=(1.0+1.5*p_eps**2)/bt0**2
!  p_eb-<E.B> (V*T/m)
      p_eb=0.1*bt0/(2.0*z_pi*r0)
!  p_fhat-mu_0*F/(dPsi/dr) (rho/m)
      p_fhat=p_q/p_eps
!  p_fm(3)-poloidal moments of geometric factor for PS viscosity (-)
      DO i=1,3
        p_fm(i)=0.0
	ENDDO
      IF(p_eps.gt.0.0) THEN
        DO i=1,3
          p_fm(i)=i*((1.0-SQRT(1.0-p_eps**2))/p_eps)**(2.0*i)
     #            *(1.0+i*SQRT(1.0-p_eps**2))/((1.0-p_eps**2)**1.5
     #            *(p_q*r0)**2)
        ENDDO   
      ENDIF
!  p_ft-trapped fraction (-)
      p_ft=1.46*SQRT(p_eps)
!  p_grbm2-<grad(rho)**2/B**2> (rho**2/m**2/T**2)
      p_grbm2=1.0/bt0**2
!  p_ngrth-<n.grad(Theta)> (1/m)
      p_ngrth=1.0/(p_q*r0)
!Set the output unit and open file
      nout=10
      OPEN(unit=nout,
     #     file='out_nclass_pt.dat',
     #     status='unknown')
!     Echo out updated namelist data
      WRITE(nout,indata)
!Evaluate friction and viscosity coefficients and flows
      CALL NCLASS(k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2,
     #            p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi,
     #            p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz,m_s,
     #            jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii,capm_ii,
     #            capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,sqz_s,upar_s,
     #            utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss,chit_ss,
     #            dp_ss,dt_ss,iflag)
!Check warning flags
      IF(iflag.eq.-1) THEN
        label='WARNING:NCLASS-no potato orbit viscosity'
        CALL WRITE_LINE(nout,label,0,0)
      ELSEIF(iflag.eq.-2) THEN
        label='WARNING:NCLASS-Pfirsch-Schluter viscosity'
        CALL WRITE_LINE(nout,label,0,0)
      ELSEIF(iflag.eq.-3) THEN
        label='WARNING:NCLASS-no banana viscosity'
        CALL WRITE_LINE(nout,label,0,0)
      ELSEIF(iflag.eq.-4) THEN
        label='WARNING:NCLASS-no viscosity'
        CALL WRITE_LINE(nout,label,0,0)
      ENDIF
!Check error flags
      IF(iflag.gt.0) THEN
        IF(k_out.gt.0) THEN
          IF(iflag.eq.1) THEN
            label='ERROR:NCLASS-k_order must be 2 or 3, k_order='
            idum(1)=k_order
            CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag.eq.2) THEN
            label='ERROR:NCLASS-require 1<m_i<mx_mi, m_i='
            idum(1)=m_i
            CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag.eq.3) THEN
            label='ERROR:NCLASS-require 0<m_z<mx_mz, m_z='
            idum(1)=m_z
            CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag.eq.4) THEN
            label='ERROR:NCLASS-require 0<m_s<mx_ms, m_s='
            idum(1)=m_s
            CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag.eq.5) THEN
            label='ERROR:NCLASS-inversion of flow matrix failed'
            CALL WRITE_LINE(nout,label,0,0)
          ENDIF
        ENDIF
        GOTO 1000
      ENDIF
!Check for optional output
      IF(k_out.gt.1) THEN
!  Species identification
        label='     *** Species Identification ***'
        CALL WRITE_LINE(nout,label,2,1)
        label='     Isotope     Species      Charge        Mass'//
     #        '     Density Temperature  Chg Factor'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -           -     Coulomb         AMU'//
     #        '       /m**3         keV           -'
        CALL WRITE_LINE(nout,label,0,0)
        DO i=1,m_s
          im=jm_s(i)
          iz=jz_s(i)
          iza=IABS(iz)
          idum(1)=im
          idum(2)=i
          idum(3)=iz
          rdum(1)=amu_i(im)
          rdum(2)=den_iz(im,iza)
          rdum(3)=temp_i(im)
          rdum(4)=xi_s(i)
          CALL WRITE_IR(nout,3,idum,4,rdum,2)
        ENDDO
!  Friction coefficients
        label='     *** Friction Coefficients ***'
        CALL WRITE_LINE(nout,label,2,0)
        DO im=1,m_i
          DO jm=1,m_i
            CALL WRITE_LINE(nout,' ',0,0)            
            label='  Isotopes ='
            idum(1)=im
            idum(2)=jm
            CALL WRITE_LINE_IR(nout,label,2,idum,0,rdum,0)           
!           Mkl
            label='  Mkl(-) ='
            CALL WRITE_LINE(nout,label,0,0)           
            DO k=1,k_order
              DO l=1,k_order
                rdum(l)=capm_ii(k,l,im,jm)
              ENDDO
              CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
            ENDDO            
!           Nkl
            label='  Nkl(-) ='
            CALL WRITE_LINE(nout,label,0,0)           
            DO k=1,k_order
              DO l=1,k_order
                rdum(l)=capn_ii(k,l,im,jm)
              ENDDO
              CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
            ENDDO   
          ENDDO    
        ENDDO
!  Reduced friction coefficients
        label='     *** Reduced Friction Coefficients ***'
        CALL WRITE_LINE(nout,label,2,0)
!       cal(M)
        DO im=1,m_i
          CALL WRITE_LINE(nout,' ',0,0)
          label='  calMij (kg/m**3/s) for Isotope ='
          idum(1)=im
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)           
          DO j=1,k_order
            DO k=1,k_order
              rdum(k)=calm_i(j,k,im)
            ENDDO
            CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO            
        ENDDO   
!       cal(N)
        DO im=1,m_i
          CALL RARRAY_ZERO(3,dum)
          DO jm=1,m_i
            CALL WRITE_LINE(nout,' ',0,0)
            label='  calNij (kg/m**3/s) for Isotopes ='
            idum(1)=im
            idum(2)=jm
            CALL WRITE_LINE_IR(nout,label,2,idum,0,rdum,0)           
            DO k=1,k_order
              dum(k)=dum(k)-caln_ii(k,1,im,jm)
              DO l=1,k_order
                rdum(l)=caln_ii(k,l,im,jm)
              ENDDO
              CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
            ENDDO            
          ENDDO   
!         Momentum check
          CALL WRITE_LINE(nout,' ',0,0)
          label='Momentum check for Isotope ='
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          label='  sum_b[-calNk1] (kg/m**3/s) ='
          CALL WRITE_LINE_IR(nout,label,0,idum,k_order,dum,2)           
          label='  calMk1 (kg/m**3/s)         ='
          CALL WRITE_LINE_IR(nout,label,0,idum,k_order,calm_i(1,1,im),2)           
          label='  Difference ='
          DO k=1,k_order
            dum(k+3)=dum(k)-calm_i(k,1,im)
          ENDDO
          CALL WRITE_LINE_IR(nout,label,0,idum,k_order,dum(4),2)           
        ENDDO
!  Normalized viscosities
        label='     *** Normalized Viscosities ***'
        CALL WRITE_LINE(nout,label,2,0)
        DO i=1,m_s
          CALL WRITE_LINE(nout,' ',0,0)
          label='  muij (kg/m**3/s) for Species ='
          idum(1)=i
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          DO k=1,k_order
            DO l=1,k_order
              rdum(l)=ymu_s(k,l,i)
            ENDDO
            CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO            
        ENDDO
!  Normalized parallel flows
        label='     *** Normalized Parallel Flows ***'
        CALL WRITE_LINE(nout,label,2,0)
        DO i=1,m_s
          CALL WRITE_LINE(nout,' ',0,0)
          label='  upar_ij (T*m/s) for Species ='
          idum(1)=i
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          DO j=1,k_order
            DO k=1,k_order
              rdum(k)=upar_s(j,k,i)
            ENDDO
            CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO
        ENDDO
!  Normalized poloidal flows
        label='     *** Normalized Poloidal Flows ***'
        CALL WRITE_LINE(nout,label,2,0)
        DO i=1,m_s
          CALL WRITE_LINE(nout,' ',0,0)
          label='  utheta_ij (m/s/T) for Species ='
          idum(1)=i
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          DO j=1,k_order
            DO k=1,k_order
              rdum(k)=utheta_s(j,k,i)
            ENDDO
            CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO
        ENDDO
!  Radial particle fluxes
        label='     *** Radial Particle Fluxes ***'
        CALL WRITE_LINE(nout,label,2,1)
        label='     Species          BP          PS          CL'//
     #        '       <E.B>         src       total'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -     /m**2/s     /m**2/s     /m**2/s'//
     #        '     /m**2/s     /m**2/s     /m**2/s'
        CALL WRITE_LINE(nout,label,0,0)
        CALL RARRAY_ZERO(6,dum)
        CALL RARRAY_ZERO(6,edum)
!       Load into rdum the five flux components and total
!       Load into dum the z-weighted + charge (ion) components
!       Load into edum the z-weighted - charge (electron) components 
        DO i=1,m_s
          iz=jz_s(i)
          idum(1)=i
          CALL RARRAY_COPY(5,gfl_s(1,i),1,rdum,1)
          rdum(6)=RARRAY_SUM(5,rdum,1)
          CALL WRITE_IR(nout,1,idum,6,rdum,2)           
          DO k=1,5
            IF(iz.gt.0) THEN
              dum(k)=dum(k)+iz*gfl_s(k,i)
            ELSE
              edum(k)=edum(k)+iz*gfl_s(k,i)
            ENDIF
          ENDDO
          IF(iz.gt.0) THEN
            dum(6)=dum(6)+iz*rdum(6)
          ELSE
            edum(6)=edum(6)+iz*rdum(6)
          ENDIF
        ENDDO
!  Ambipolarity check
        label='     *** Ambipolarity Check ***'
        CALL WRITE_LINE(nout,label,2,1)
        label='  Sum_i Z_i*Gamma_i ='
        CALL WRITE_LINE(nout,label,0,0)
        idum(1)=0
        CALL WRITE_IR(nout,1,idum,6,dum,2)           
        label='  Sum_e Z_e*Gamma_e ='
        CALL WRITE_LINE(nout,label,0,0)
        CALL WRITE_IR(nout,1,idum,6,edum,2)           
        label='  Difference ='
        CALL WRITE_LINE(nout,label,0,0)
        DO k=1,6
          rdum(k)=edum(k)+dum(k)
        ENDDO
        CALL WRITE_IR(nout,1,idum,6,rdum,2)
!  Particle transport is returned in three forms
!  Particle flux consistency check
        label='     *** Particle Flux Consistency Check ***'           
        CALL WRITE_LINE(nout,label,2,1)
        label='     Species       gfl_s   dn_s,vn_s   dp_s,dt_s'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -     /m**2/s     /m**2/s     /m**2/s'
        CALL WRITE_LINE(nout,label,0,0)
        DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          idum(1)=i
!    Gamma total is sum of components
          rdum(1)=RARRAY_SUM(4,gfl_s(1,i),1)
!    Gamma total is diffusive plus pinch
          rdum(2)=-dn_s(i)*(grp_iz(im,iza)-den_iz(im,iza)*grt_i(im))
     #            /temp_i(im)+den_iz(im,iza)*(vn_s(i)+veb_s(i))
!    Gamma total is sum over T' and p' of all species
          rdum(3)=0.0
          DO j=1,m_s
            jm=jm_s(j)
            jza=IABS(jz_s(j))
            rdum(3)=rdum(3)-dt_ss(j,i)*grt_i(jm)/temp_i(jm)
     #              -dp_ss(j,i)*grp_iz(jm,jza)/den_iz(jm,jza)/temp_i(jm)
          ENDDO
          rdum(3)=(rdum(3)+veb_s(i))*den_iz(im,iza)
          CALL WRITE_IR(nout,1,idum,3,rdum,2)
        ENDDO
!  Particle diffusion, velocity
        label='     *** Particle Diffusion, Velocity ***'           
        CALL WRITE_LINE(nout,label,2,1)
        label='     Species         D_s         V_s'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -      m**2/s         m/s'
        CALL WRITE_LINE(nout,label,0,0)
        DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          idum(1)=i
          rdum(1)=dn_s(i)
          rdum(2)=vn_s(i)+veb_s(i)+gfl_s(5,i)/den_iz(im,iza)
          CALL WRITE_IR(nout,1,idum,2,rdum,2)
        ENDDO
!  Particle diffusivity matrices
!    On p'/p           
        label='     *** Particle Diffusivity Matrices ***'           
        CALL WRITE_LINE(nout,label,2,1)
        label='     Species    on dp/dr ... by Species'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -      m**2/s ...'
        CALL WRITE_LINE(nout,label,0,0)
        DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,dp_ss(1,i),2)
        ENDDO
!    On T'/T
         label='     Species    on dT/dr ... by Species'
        CALL WRITE_LINE(nout,label,1,0)
        label='           -      m**2/s ...'
        CALL WRITE_LINE(nout,label,0,0)
        DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,dt_ss(1,i),2)
        ENDDO
!  Radial conduction fluxes
        label='     *** Radial Conduction Fluxes ***'
        CALL WRITE_LINE(nout,label,2,1)
        label='     Species          BP          PS          CL'//
     #        '       <E.B>         src       total'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -      w/m**2      w/m**2      w/m**2'//
     #        '      w/m**2      w/m**2      w/m**2'
        CALL WRITE_LINE(nout,label,0,0)
        CALL RARRAY_ZERO(6,dum)
        CALL RARRAY_ZERO(6,edum)
        DO i=1,m_s
          iz=jz_s(i)
          idum(1)=i
          CALL RARRAY_COPY(5,qfl_s(1,i),1,rdum,1)
          rdum(6)=RARRAY_SUM(5,qfl_s(1,i),1)
          CALL WRITE_IR(nout,1,idum,6,rdum,2)           
          DO k=1,5
            IF(iz.gt.0) THEN
              dum(k)=dum(k)+qfl_s(k,i)
            ELSE
              edum(k)=edum(k)+qfl_s(k,i)
            ENDIF
          ENDDO
          IF(iz.gt.0) THEN
            dum(6)=dum(6)+rdum(6)
          ELSE
            edum(6)=edum(6)+rdum(6)
          ENDIF
        ENDDO
!     Total electron radial conduction flux
        label='  Sum of electron conduction fluxes ='
        CALL WRITE_LINE(nout,label,1,0)
        idum(1)=0
        CALL WRITE_IR(nout,1,idum,6,edum,2)
!     Total ion radial conduction flux
        label='  Sum of ion conduction fluxes ='
        CALL WRITE_LINE(nout,label,1,0)
        idum(1)=0
        CALL WRITE_IR(nout,1,idum,6,dum,2)
!  Heat conduction is returned in two forms, add diagonal form
!  Heat flux consistency check
        label='     *** Heat Flux Consistency Check ***'           
        CALL WRITE_LINE(nout,label,2,1)
        label='     Species       qfl_s   dq_s,vq_s,  cp_s,ct_s'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -      w/m**2      w/m**2      w/m**2'
        CALL WRITE_LINE(nout,label,0,0)
        DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          idum(1)=i
!    Conduction total is sum of components
          rdum(1)=RARRAY_SUM(5,qfl_s(1,i),1)
!    Diagonal conductivity plus convective velocity
          dq_s(i)=chit_ss(i,i)+chip_ss(i,i)
          vq_s(i)=rdum(1)/(den_iz(im,iza)*z_j7kv*temp_i(im))
     #            +dq_s(i)*grt_i(im)/temp_i(im)
          rdum(2)=den_iz(im,iza)*z_j7kv*(-dq_s(i)*grt_i(im)
     #                                   +vq_s(i)*temp_i(im))
!    Conduction total is sum over T' and p' of all species
          rdum(3)=0.0
          DO j=1,m_s
            jm=jm_s(j)
            jza=IABS(jz_s(j))
            rdum(3)=rdum(3)-chit_ss(j,i)*grt_i(jm)/temp_i(jm)
     #                     -chip_ss(j,i)*grp_iz(jm,jza)
     #                      /den_iz(jm,jza)/temp_i(jm)
          ENDDO
          rdum(3)=(rdum(3)+qeb_s(i))*den_iz(im,iza)*temp_i(im)*z_j7kv
          CALL WRITE_IR(nout,1,idum,3,rdum,2)
        ENDDO
!  Heat conduction, velocity
        label='     *** Heat Conduction, Velocity ***'           
        CALL WRITE_LINE(nout,label,2,1)
        label='     Species         X_s         V_s'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -      m**2/s         m/s'
        CALL WRITE_LINE(nout,label,0,0)
        DO i=1,m_s
          idum(1)=i
          rdum(1)=dq_s(i)
          rdum(2)=vq_s(i)
          CALL WRITE_IR(nout,1,idum,2,rdum,2)
        ENDDO
!     Effective ion heat conduction, velocity
        label='  Effective ion conduction, velocity ='
        CALL WRITE_LINE(nout,label,1,0)
        idum(1)=0
        CALL RARRAY_ZERO(4,rdum)
        DO i=1,m_s
          im=jm_s(i)
          iz=jz_s(i)
          iza=IABS(iz)
          IF(iz.gt.0) THEN
            rdum(1)=rdum(1)+(chit_ss(i,i)+chip_ss(i,i))*den_iz(im,iza)
            rdum(2)=rdum(2)+RARRAY_SUM(5,qfl_s(1,i),1)/temp_i(im)/z_j7kv
            rdum(3)=rdum(1)+den_iz(im,iza)
            rdum(4)=rdum(4)+(chit_ss(i,i)+chip_ss(i,i))*den_iz(im,iza)
     #                      *grt_i(im)/temp_i(im)
          ENDIF
        ENDDO
        rdum(1)=rdum(1)/rdum(3)
        rdum(2)=(rdum(2)+rdum(4))/rdum(3)
        CALL WRITE_IR(nout,1,idum,2,rdum,2)
!  Thermal diffusivity matrices
!    On p'/p           
        label='     *** Thermal Diffusivity Matrices ***'           
        CALL WRITE_LINE(nout,label,2,1)
        label='     Species    on dp/dr ... by Species'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -      m**2/s ...'
        CALL WRITE_LINE(nout,label,0,0)
        DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,chip_ss(1,i),2)
        ENDDO
!    On T'/T
         label='     Species    on dT/dr ... by Species'
        CALL WRITE_LINE(nout,label,1,0)
        label='           -      m**2/s ...'
        CALL WRITE_LINE(nout,label,0,0)
        DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,chit_ss(1,i),2)
        ENDDO
!  Radial energy (conduction+convection) fluxes
        label='     *** Radial Energy (Cond+Conv) Fluxes ***'
        CALL WRITE_LINE(nout,label,2,1)
        label='     Species    ban-plat          PS   classical'//
     #        '       <E.B>  extern src       total'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -      w/m**2      w/m**2      w/m**2'//
     #        '      w/m**2      w/m**2      w/m**2'
        CALL WRITE_LINE(nout,label,0,0)
        CALL RARRAY_ZERO(6,dum)
        DO i=1,m_s
          idum(1)=i
          im=jm_s(i)
          iz=jz_s(i)
          DO k=1,5
            rdum(k)=qfl_s(k,i)+2.5*gfl_s(k,i)*temp_i(im)*z_j7kv
          ENDDO
          rdum(6)=RARRAY_SUM(5,rdum,1)
          CALL WRITE_IR(nout,1,idum,6,rdum,2)           
          IF(iz.gt.0) THEN
            DO k=1,5
              dum(k)=dum(k)+qfl_s(k,i)+2.5*gfl_s(k,i)*temp_i(im)*z_j7kv
            ENDDO
            dum(6)=dum(6)+rdum(6)
          ENDIF
        ENDDO
!     Total ion radial energy flux
        label='  Sum of ion energy fluxes ='
        CALL WRITE_LINE(nout,label,1,0)
        idum(1)=0
        CALL WRITE_IR(nout,1,idum,6,dum,2)           
!  Bootstrap current is returned in two forms
!  Bootstrap current consistency check
        label='     *** Bootstrap Current Consistency Check ***'           
        CALL WRITE_LINE(nout,label,2,1)
        label='      p_bsjb   dp_s,dt_s'
        CALL WRITE_LINE(nout,label,0,0)
        label='    A*T/m**2    A*T/m**2'
        CALL WRITE_LINE(nout,label,0,0)
!    Total bootstrap current
        rdum(1)=p_bsjb
!    Bootstrap current summed over T' and p' components
        rdum(2)=0.0
        DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          rdum(2)=rdum(2)-bsjbt_s(i)*grt_i(im)/temp_i(im)
     #                   -bsjbp_s(i)*grp_iz(im,iza)
     #                    /den_iz(im,iza)/temp_i(im)
        ENDDO
        CALL WRITE_IR(nout,0,idum,2,rdum,2)
!  Bootstrap current arrays
!    On p'/p           
        label='     *** Bootstrap Current Arrays ***'           
        CALL WRITE_LINE(nout,label,2,1)
        label='     Species    on dp/dr    on dT/dr'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -    A*T/m**2    A*T/m**2'
        CALL WRITE_LINE(nout,label,0,0)
        DO i=1,m_s
          idum(1)=i
          rdum(1)=bsjbp_s(i)
          rdum(2)=bsjbt_s(i)
          CALL WRITE_IR(nout,1,idum,2,rdum,2)
        ENDDO
!  Current response to external source
        label='     *** Bootstrap and External Source Current ***'           
        CALL WRITE_LINE(nout,label,2,1)
        label='   Bootstrap    External       Total'
        CALL WRITE_LINE(nout,label,0,0)
        label='    A*T/m**2    A*T/m**2    A*T/m**2'
        CALL WRITE_LINE(nout,label,0,0)
        rdum(1)=p_bsjb
        rdum(2)=p_exjb
        rdum(3)=rdum(1)+rdum(2)
        CALL WRITE_IR(nout,0,idum,3,rdum,2)
!  Flow velocities on outside midplane
        label='     *** Flow Velocities on Outside Midplane ***'
        CALL WRITE_LINE(nout,label,2,1)
        label='     Species       v-tor       v-pol       v-par'//
     #        '      v-perp'
        CALL WRITE_LINE(nout,label,0,0)
        label='           -         m/s         m/s         m/s'//
     #        '         m/s'
        CALL WRITE_LINE(nout,label,0,0)                        
        btor=bt0/(1.0+p_eps)
        bpol=btor/p_fhat
        btot=SQRT(btor**2+bpol**2)*btor/ABS(btor)
        DO i=1,m_s
          idum(1)=i
          im=jm_s(i)
          iz=jz_s(i)
          iza=IABS(iz)
          ppr=p_fhat*grp_iz(im,iza)*z_j7kv
     #        /(z_coulomb*iz*den_iz(im,iza))+p_fhat*p_grphi
          uthai=utheta_s(1,1,i)+utheta_s(1,2,i)+utheta_s(1,3,i)
!         Toroidal
          rdum(1)=uthai*btor-ppr/btor
!         Poloidal
          rdum(2)=uthai*bpol
!         Parallel
          rdum(3)=uthai*btot-ppr/btot
!         Perpendicular
          rdum(4)=ppr*bpol/btot/btor
          CALL WRITE_IR(nout,1,idum,4,rdum,2)           
        ENDDO
!  Miscellaneous parameters
        label='     *** Miscellaneous Parameters ***'
        CALL WRITE_LINE(nout,label,2,1)
        label='<J_bs.B>/Bt0 (A/m**2) ='
        rdum(1)=p_bsjb/bt0
        CALL WRITE_LINE_IR(nout,label,0,idum,1,rdum,2)
        label='<J_ex.B>/Bt0 (A/m**2) ='
        rdum(1)=p_exjb/bt0
        CALL WRITE_LINE_IR(nout,label,0,idum,1,rdum,2)
        label='eta parallel (Ohm*m) ='
        rdum(1)=p_etap
        CALL WRITE_LINE_IR(nout,label,0,idum,1,rdum,2)
      ENDIF
 1000 STOP
      END
      
