
# Boundary conditions
bbb.ncore[0] = 2.e19    #hydrogen ion density on core
bbb.iflcore = 0         #flag; =0, fixed Te,i; =1, fixed power on core
bbb.tcoree = 250.       #core Te
bbb.tcorei = 250.       #core Ti
bbb.tedge = 2.          #fixed wall,pf Te,i if istewcon=1, etc
bbb.recycp[0] = 0.95    #hydrogen recycling coeff at plates
bbb.recycm = 0.1

# Transport coefficients (m**2/s)
bbb.difni[0] = 1.       #D for radial hydrogen diffusion
bbb.kye = 1.            #chi_e for radial elec energy diffusion
bbb.kyi = 1.            #chi_i for radial ion energy diffusion
bbb.travis[0] = 1.      #eta_a for radial ion momentum diffusion

# Flux limits
bbb.flalfe = 0.21       #electron parallel thermal conduct. coeff
bbb.flalfi = 0.21       #ion parallel thermal conduct. coeff
bbb.flalfv = 0.5        #ion parallel viscosity coeff
bbb.flalfgx = 1.        #neut. gas part. flux in poloidal direction
bbb.flalfgy = 1.        #neut. gas part. flux in radial direction
bbb.flalfgxy = 1.       #neut. gas part. flux in  mixed derivatives
bbb.flalftgx = 1.       #neut. gas thermal flux, poloidal direction
bbb.flalftgy = 1.       #neut. gas thermal flux, radial direction
bbb.lgmax = 0.1		#max scale length for flalfgx,y
bbb.lgtmax = 0.1	#max scale length for flalftgx,y
bbb.lgvmax = 0.1	#max scale length for flalfvgx,y


# Solver package
bbb.svrpkg = "nksol"    #Newton solver using Krylov method
bbb.premeth = "ilut"    #Solution method for precond. Jacobian matrix

# Parallel neutral momentum equation
bbb.isupgon[0] = 1
bbb.ineudif = 2             #=2 for evolving pg=ng*tg variable
bbb.isngon[0] = 0
com.ngsp = 1
com.nhsp = 2

bbb.ziin[1] = 0

# Atomic data switches
com.istabon = 10            #=10 specifics hydrogen data file ehr2.dat
