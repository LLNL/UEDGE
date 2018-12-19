psi

{
Nmat_ptm = integer 7			# max number of materials
Nprj_ptm = integer 8			# max number of projectiles
Npsi_ptm = integer Nmat_ptm*Nprj_ptm

#   mesh dimensions
Nt_ptm   = integer 100			# max number of nodes in Ti mesh
Ntd_ptm  = integer 100			# max number of nodes in Twall mesh
Ntt_ptm  = integer Nt_ptm*Ntd_ptm
Nphi_ptm = integer 100			# max number of nodes in Phi/T mesh
Ntde_ptm = integer 30			# max number of nodes in Td mesh
Nrd_ptm  = integer 30			# max number of nodes in rd mesh

#   arrays of deduced dimensions
Ntpsi_ptm  = integer Nt_ptm*Npsi_ptm
Nttpsi_ptm = integer Ntpsi_ptm*Ntd_ptm
Ntpm_ptm   = integer Nt_ptm*Nphi_ptm*Nmat_ptm
Ntdmat_ptm = integer Nt_ptm*Nmat_ptm
Nrtem_ptm  = integer Nrd_ptm*Ntde_ptm*Nmat_ptm
}


***** const :
c		real		/2.99792458d+08/
hplanck		real		/6.62606876d-34/
kboltz		real		/1.3806503d-23/
sigma		real		/5.6704d-08/
eps0		real		/8.854187817d-12/
pi		real		/3.141592653589793d+00/
one		complex		/(1.d0,0.d0)/
ione		complex		/(0.d0,1.d0)/


***** psiusr :
N_teta		integer		/99/
E_eps		real		/1.d-3/
Teps		real		/1.d-6/
Ti_eps		real		/1.d-2/


***** psitab :

IZmatter(Nmat_ptm)	integer		# array of atomic numbers (in periodic system of chemical elements) of wall materials. Note, order is important
Namfpsi (Nmat_ptm)	character*32	# array of file names into which the PSI data are stored for each material
IZprjctl(Nprj_ptm)	integer		# array of atomic numbers (in periodic system of chemical elements) of projectiles
AZprjctl(Nprj_ptm)	real	[amu]	# array for the mass of projectiles

# some physical parameters of wall materials
mat_mas(Nmat_ptm)	real	[amu]	# matter atom mass
mat_wsub(Nmat_ptm)	real	[eV]	# sublimation heat
mat_a0(Nmat_ptm)	real	[A]	# lattice constant
mat_Tmelt(Nmat_ptm)	real	[K]	# melting temperature
mat_Hmelt(Nmat_ptm)	real	[eV]	# melting enthalpy
mat_Tsub(Nmat_ptm)	real	[K]	# sublimation temperature
mat_Wf(Nmat_ptm)	real	[eV]	# work function

#   Table of chemical parameters of solid matter
psi_ro(Ntd_ptm,Nmat_ptm)	real	[g/cm^3]		# matter number density
psi_cp(Ntd_ptm,Nmat_ptm)	real	[J/(g grad K)]		# heat capacity
psi_kp(Ntd_ptm,Nmat_ptm)	real	[W/(m grad K)]		# heat conductivity
psi_ep(Ntd_ptm,Nmat_ptm)	real				# Black-body emissivity coefficient
psi_hs(Ntd_ptm,Nmat_ptm)	real	[eV]			# thermal vaporization heat
psi_gs(Ntd_ptm,Nmat_ptm)	real	[10^20 atoms/(cm^2 s)]	# (LOG) thermal vaporization flux
psi_pbtd(Ntd_ptm,Nmat_ptm)	real				# coefficient of thermal vaporization

#   Table of PSI parameters of hydrogen+solid matter (self-sputtering including)
rn_mx(Nt_ptm,Npsi_ptm)			real		# <Rn(Teff)>
re_mx(Nt_ptm,Npsi_ptm)			real		# <Re(Teff)>
ypsp_mx(Nt_ptm,Npsi_ptm)		real		# <Yphsput(Teff)>
epsp_mx(Nt_ptm,Npsi_ptm)		real		# <Ephsput(Teff)>
ycsp_mx(Nt_ptm,Ntd_ptm,Npsi_ptm)	real		# <Ychsput(Teff,Twall)>
yres_mx(Nt_ptm,Ntd_ptm,Npsi_ptm)	real		# <Yres(Teff,Twall)>

dse_mx(Nt_ptm,Nphi_ptm,Nmat_ptm)	real
dtheps_mx(Nrd_ptm,Ntde_ptm,Nmat_ptm)	real

rdk(Nrd_ptm)			real
tkkk(Ntde_ptm)			real

#   limits of meshes over T_plasma and T_wall

# log mesh on [Ti_min,Ti_max]
Tmin_pt		real	[eV]	/1.d-1/
Tmax_pt		real	[eV]	/1.d+4/

# log mesh on [Td_min,Td_max] 
Tdmin_pt	real	[K]	/3.d+2/
Tdmax_pt	real	[K]	/7.d+3/

# mesh on [phi_min,phi_max]
Phimin_pt	real	/-1.d+1/	# Phi/T
Phimax_pt	real	/+1.d+1/	# Phi/T

# log mesh on [rd_min,rd_max]
rdmin_pt	real	[micron]	/1.d-2/
rdmax_pt	real	[micron]	/1.d+1/

# log mesh on [Td_mine,Td_maxe]
Tdmine_pt	real	[K]	/3.d+2/
Tdmaxe_pt	real	[K]	/5.d+3/


***** psitab_s local hidden :

psi_ro_s(Ntdmat_ptm)		real
psi_cp_s(Ntdmat_ptm)		real
psi_kp_s(Ntdmat_ptm)		real
psi_ep_s(Ntdmat_ptm)		real
psi_gs_s(Ntdmat_ptm)		real
psi_hs_s(Ntdmat_ptm)		real
psi_pbtd_s(Ntdmat_ptm)		real

rn_mx_s(Ntpsi_ptm)		real
re_mx_s(Ntpsi_ptm)		real
ypsp_mx_s(Ntpsi_ptm)		real
epsp_mx_s(Ntpsi_ptm)		real
ycsp_mx_s(Nttpsi_ptm)		real
yres_mx_s(Nttpsi_ptm)		real
dse_mx_s(Ntpm_ptm)		real
dtheps_mx_s(Nrtem_ptm)		real



***** psitabint :

Nmat_pt		integer		# number of wall materials
Nprj_pt		integer		# number of projectiles
Npsi_pt		integer		# Nmat_pt*Nprj_pt
Nphi_pt		integer		# number of nodes in Phi/T mesh
Nt_pt		integer		# number of nodes in Tplasma mesh
Ntd_pt		integer		# number of nodes in Twall mesh
Nttm		integer		# Nt_ptm*Nt_ptm
Nttdm		integer		# Ntd_ptm*Nt_ptm
Ntphm		integer		# Nphi_ptm*Nt_ptm
Nttdmc		integer		# Ntd_ptm*(Nt_ptm+1)
Ntm		integer		# Nt_ptm
Ntde_pt		integer		# number of nodes in Td mesh
Nrd_pt		integer		# number of nodes in radius
Nrtdem		integer		# Nrd_ptm*Ntde_ptm 

I_pwi(Nprj_ptm,Nmat_ptm)	integer		# packed indices
I_prj(Npsi_ptm)			integer		# indices of projectiles
I_wal(Npsi_ptm)			integer		# indiced of wall matterials

atmi		real
atma		real
datm		real
patm		real		# [Tmin,Tmax]
atdmi		real
atdma		real
datdm		real
patdm		real		# [Tdmin,Tdmax]
aphmi		real
aphma		real
daphm		real
paphm		real		# [Phimin,Phimax]
atdmie		real
atdmae		real
datdme		real
patdme		real		# [Tdmine,Tdmaxe]
armi		real
arma		real
darm		real
parm		real		# [rdmin,rdmax]

iatm		integer
iatdm		integer
iaphm		integer
iatdme		integer
iarm		integer

sml		real	/+1.d-97/
asml		real	/-2.24d2/
big		real	/+1.d+97/
abig		real	/+2.24d2/


***** psitabmxx :

grn_mx(Nt_ptm)			real		# <Rn(Teff)>
gre_mx(Nt_ptm)			real		# <Re(Teff)>
gypsp_mx(Nt_ptm)		real		# <Yphsput(Teff)>
gepsp_mx(Nt_ptm)		real		# <Ephsput(Teff)>
gycsp_mx(Nt_ptm,Ntd_ptm)	real		# <Ychsput(Teff,Twall)>
gyres_mx(Nt_ptm,Ntd_ptm)	real		# <Yres(Teff,Twall)>
gdse_mx(Nt_ptm,Nphi_ptm)	real		# <Ysee(Teff,Phi/T)>
gytheps_mx(Nrd_ptm,Ntde_ptm)	real		# <Ytheps(rd,Td)>


***** psiparloc :

imat		integer

rotd		real
cptd		real
kptd		real
eptd		real
hstd		real
gstd		real
pbtd		real
fysee		real
theps		real
frn(Nprj_ptm)	real
fre(Nprj_ptm)	real
fypsp(Nprj_ptm)	real
fepsp(Nprj_ptm)	real
fycsp(Nprj_ptm)	real
fyres(Nprj_ptm)	real


***** psi_sub :

psidotab		subroutine
