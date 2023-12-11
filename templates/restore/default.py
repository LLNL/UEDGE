bbb.restart = 1			#Begin from savefile
bbb.allocate()			#allocate plasma/neutral arrays
hdf5_restore("solution.h5")    #read in the solution from pfb file