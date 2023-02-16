      subroutine nclassb(dens0, ddens0_drho, temp0, dtemp0_drho)

c!   Mimics PROGRAM NCLASS_PT_DR in NCLASS_PT tar file from NTCC website
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
      Use(Vars_nclass)

c!...Declaration of local variables
      integer i,j
      real ds

c!...Input arguments     
      real dens0         #-local density, 1/m**3
      real ddens0_drho   #-local density gradient, (1/m**3)/rho
      real temp0         #-local temperature, keV
      real dtemp0_drho   #-local temperature gradient, keV/rho

c       print *, dens0, ddens0_drho, temp0, dtemp0_drho    
c       print *, 'mx_mi:', mx_mi
c       print *, 'mx_ms:', mx_ms
c       print *, 'mx_mz:', mx_mz


c!...Calculate input for NCLASS

c!...only electrons and protons    
      amu_i(1)=5.4463e-4 
      amu_i(2)=1.0

c!  grt_i(i)-temperature gradient of i (keV/rho)
c!  grp_iz(i,z)-pressure gradient of i,z (keV/m**3/rho)
       DO i=1,2
         grt_i(i)=dtemp0_drho
         temp_i(i)=temp0
         DO j=1,1
           grp_iz(i,j)=dtemp0_drho*dens0 + ddens0_drho*temp0
           den_iz(i,j)=dens0
         ENDDO   
       ENDDO

c!  m_i-number of isotopes (1<mi<mx_mi+1)
c!  m_z-highest charge state of all species (0<mz<mx_mz+1)
      m_i=0
      m_z=0
      DO i=1,2
         ds=0.0
       DO j=1,1
         IF(den_iz(i,j).gt.c_den) THEN
             ds=den_iz(i,j)
           IF(j.gt.m_z) m_z=j
         ENDIF
       ENDDO
       IF(ds.gt.c_den) m_i=m_i+1
      ENDDO

c!  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
       c_potb=e0*bt0/2/q0**2
c!  c_potl-q(0)*R(0) (m)
       c_potl=q0*r0
c!  p_b2-<B**2> (T**2)
       p_b2=bt0**2*(1.0+0.5*p_eps**2)
c!  p_bm2-<1/B**2> (/T**2)
       p_bm2=(1.0+1.5*p_eps**2)/bt0**2
c!  p_eb-<E.B> (V*T/m)
       p_eb=0.1*bt0/(2.0*z_pi*r0)
c!  p_fhat-mu_0*F/(dPsi/dr) (rho/m)
       p_fhat=p_q/p_eps
c!  p_fm(3)-poloidal moments of geometric factor for PS viscosity (-)
      DO i=1,3
         p_fm(i)=0.0
      ENDDO
       IF(p_eps.gt.0.0) THEN
         DO i=1,3
           p_fm(i)=i*((1.0-SQRT(1.0-p_eps**2))/p_eps)**(2.0*i)
     .            *(1.0+i*SQRT(1.0-p_eps**2))/((1.0-p_eps**2)**1.5
     .            *(p_q*r0)**2)
         ENDDO   
      ENDIF
c!  p_ft-trapped fraction (-)
      DO i=1,3 
        p_ft(i)=1.46*SQRT(p_eps)
      ENDDO
c!  p_grbm2-<grad(rho)**2/B**2> (rho**2/m**2/T**2)
      p_grbm2=1.0/bt0**2
c!  p_ngrth-<n.grad(Theta)> (1/m)
      p_ngrth=1.0/(p_q*r0)


!Evaluate friction and viscosity coefficients and flows
      CALL NCLASS(k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2,
     .            p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi,
     .            p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz,m_s,
     .            jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii,capm_ii,
     .            capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,sqz_s,upar_s,
     .            utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss,chit_ss,
     .            dp_ss,dt_ss,iflag)

      print *, "...finished call to NCLASS"

      return 
      end
