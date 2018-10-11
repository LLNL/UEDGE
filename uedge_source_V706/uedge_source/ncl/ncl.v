ncl
{
mx_mi=9  #!  mx_mi-max isotopes (-)
mx_ms=40 #!  mx_ms-max species (-)
mx_mz=18 #!  mx_mz-max charge of any species (-)
}

***** Vars_nclass:
k_out integer /2/
k_order integer /2/
k_potato integer /0/
c_den real /1.0e10/
p_eps real /0.1/
p_grphi real /-5.0e3/
p_gr2phi real /-2.0e4/
p_q real /-3.0/
r0 real /3.0/
a0 real /1.0/
e0 real /1.0/
bt0 real /-3.0/
q0 real /-5.0/

m_i integer
m_z integer
c_potb real
c_potl real
p_b2 real
p_bm2 real
p_eb real
p_fhat real
p_fm(3) real
p_ft(3) real
p_grbm2 real
p_ngrth real

amu_i(mx_mi) real
temp_i(mx_mi) real
den_iz(mx_mi,mx_mz) real

grt_i(mx_mi) real
fex_iz(3,mx_mi,mx_mz) real
grp_iz(mx_mi,mx_mz) real


z_coulomb real /1.6022e-19/
z_electronmass real /9.1095e-31/
z_j7kv real /1.6022e-16/
z_mu0 real /1.2566e-06/
z_pi real /-1.57/
z_protonmass real /1.6726e-27/

iflag integer
m_s integer
jm_s(mx_ms) integer
jz_s(mx_ms) integer
p_bsjb real
p_etap real
p_exjb real
calm_i(3,3,mx_mi) real
caln_ii(3,3,mx_mi,mx_mi) real
capm_ii(3,3,mx_mi,mx_mi) real
capn_ii(3,3,mx_mi,mx_mi) real
bsjbp_s(mx_ms) real
bsjbt_s(mx_ms) real
dn_s(mx_ms) real          
gfl_s(5,mx_ms) real       
qfl_s(5,mx_ms) real         
sqz_s(mx_ms) real
upar_s(3,3,mx_ms) real     
utheta_s(3,3,mx_ms) real
vn_s(mx_ms) real            
veb_s(mx_ms) real
qeb_s(mx_ms) real           
xi_s(mx_ms) real
ymu_s(3,3,mx_ms) real
chip_ss(mx_ms,mx_ms) real   
chit_ss(mx_ms,mx_ms) real
dp_ss(mx_ms,mx_ms) real     
dt_ss(mx_ms,mx_ms) real

nclassb(dens0:real, ddens0_drho:real, temp0:real, dtemp0_drho:real) subroutine
	#Input variables for NCLASS driver routine nclassb()
        #dens0         :local density, 1/m**3
        #ddens0_drho   :local density gradient, (1/m**3)/rho
        #temp0         :local temperature, keV
        #dtemp0_drho   :local temperature gradient, keV/rho
