bbb.mhdgeo = 1 		#use MHD equilibrium
bbb.gengrid=0		#read mesh from gridue file
com.geometry = "snull" # Defines a LSN geometry for the problem
bbb.GridFileName = "gridue_8x4" # Defines the gridue file name
com.isnonog = 1		#nonorthogal differencing used

bbb.methn = 33		#ion continuty eqn
bbb.methu = 33		#ion parallel momentum eqn
bbb.methe = 33		#electron energy eqn
bbb.methi = 33		#ion energy eqn
bbb.methg = 66		#neutral gas continuity eqn