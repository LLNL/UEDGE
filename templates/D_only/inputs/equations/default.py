
if setni:
    bbb.isnion = 1
if setup:
    bbb.isupon = 1
if sette:
    bbb.isteon = 1
if setti:
    bbb.istion = 1
if setphi:
    bbb.isphion = 1
if setphiofft:
    bbb.isphiofft = 0
if settg:
    bbb.istgon = 0
if setng:
    bbb.isngon = 0
    bbb.isngon[0] = 0
if setupg:
    bbb.isupgon = 0
    bbb.isupgon[0] = 1


# Catch-all for turning off potential equaiton in slab geometry
if bbb.mhdgeo == -1:
	bbb.isphion = 0
	bbb.isphiofft = 1
