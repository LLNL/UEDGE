bbb.matwso[0] = 1		#recycle on main-chamber wall
bbb.isnwcono = 1		#if=1, set bbb.ni[,,com.ny+1]=bbb.nwallo
bbb.isnwconi = 1		#if=1, set PF wall bbb.ni=bbb.nwalli
bbb.allocate()		#bbb.allocate() space of bbb.nwallo,i
bbb.nwallo = 1.e18
bbb.nwalli = 1.e18
bbb.igaso=100