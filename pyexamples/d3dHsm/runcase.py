# import uedge
from uedge import bbb, rundt, uedgeplots

# read UEDGE settings
exec(open("rd_d3dHsm_in.py").read())

# do a quick preliminary run to set all internals
bbb.restart=0; bbb.ftol=1e10; bbb.dtreal = 1e-6; bbb.exmain()

# Plot the mesh
uedgeplots.plotmesh()

# Take a small time step
bbb.restart=1; bbb.ftol=1e-8; 
bbb.isbcwdt=1
bbb.dtreal = 1e-14; bbb.itermx=30; bbb.exmain()

# run to steady state
case = rundt.UeRun()
case.converge(savedir=".")

# Check that we've really reached steady state
bbb.dtreal=1e20; bbb.isbcwdt=0; bbb.exmain()

# export the solution in hdf5 file
case.save_intermediate(".", "d3dHsm")

# Plot the ion density
uedgeplots.plotmeshval(bbb.ni, title='Ion density', units=r'm$^-3$')

###-refine the grid, interpolate to new grid, and restart:
com.nxleg[0,0]=20; bbb.newgeo=1; bbb.icntnunk=0
bbb.dtreal = 1e-14; bbb.isbcwdt=1; bbb.itermx=30; bbb.exmain()

# Plot the new mesh
uedgeplots.plotmesh()

# Converge to steady state
case.converge(savedir=".")

bbb.dtreal=1e20; bbb.isbcwdt=0; bbb.exmain()

# Plot electron temperature
uedgeplots.plotmeshval(bbb.te/bbb.ev, title='Electron temperature', units='eV')

input("Wait for figures to appear, then press return to exit.")
