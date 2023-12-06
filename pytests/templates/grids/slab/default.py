bbb.isfixlb[0]=0    # =1 fixes values on left boundary (nib, upb, teb, tib, yylb0)
                    # =2 for symmetry point
grd.radx = 4.5e-2       #outer "radial" wall
grd.rad0 = 0.
grd.radm = -5.0e-4      #minimum "radial" position
grd.alfyt=-2.5          #radial nonuniformity factor <0 => expanding
grd.za0 = 0.            #poloidal symmetry plane location
grd.zaxpt = 1.0         #poloidal location of flx.x-point
grd.zax = 1.5           #poloidal location of divertor plate
grd.alfxt=5.5

grd.btfix = 2.          #constant total B-field
grd.bpolfix = .2        #constant poloidal B-field
bbb.ngrid = 1
com.nxleg[0,0]=0
com.nxcore[0,0]=0
com.nxcore[0,1]=8
com.nxleg[0,1]=24
com.nycore[0]=1
com.nysol=16
com.isgrdsym = 1        #symmeterize the grid about z=0
bbb.isnfmiy = 1         #diff bbb.fmiy to preserve symmetry about X-pt

# Finite-difference algorithms [upwind, central diff, etc.]
bbb.methn = 33		#ion continuty eqn
bbb.methu = 33		#ion parallel momentum eqn
bbb.methe = 33		#electron energy eqn
bbb.methi = 33		#ion energy eqn
bbb.methg = 33		#neutral gas continuity eqn