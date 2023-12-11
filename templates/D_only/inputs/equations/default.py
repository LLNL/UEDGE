bbb.isnion = 1
bbb.isupon = 1
bbb.isteon = 1
bbb.istion = 1
bbb.isphion = 1
bbb.isphiofft = 0
bbb.istgon = 0
bbb.isngon = 0
bbb.isngon[0] = 0
bbb.isupgon = 0
bbb.isupgon[0] = 1


# Catch-all for turning off potential equaiton in slab geometry
if bbb.mhdgeo == -1:
	bbb.isphion = 0
	bbb.isphiofft = 1