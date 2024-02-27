bbb.isupwo[1] = 0
bbb.ineudif = 2
com.ngsp=1
com.nhsp=2
bbb.ziin[1]=0
bbb.travis[1] = 0.	#shouldn't be used for neutrals - but to be sure

bbb.istgon[0] = 1   # Turn on D0 temperature equation
bbb.cftiexclg = 0.  # Remove the Tg part in the Ti equation
bbb.cfdiss = 1.0
bbb.istgcore[0] = 2 #..albedo
bbb.istgpfc[0] = 4#2  #..zml flow change to 4 later
bbb.istgwc[0] = 3#2   #..flow
bbb.istglb[0] = 4
bbb.istgrb[0] = 4

bbb.cfnidhgy = 1.0
bbb.cfnidhg2 = 1.0
bbb.cfnidhdis = 1.0 #..consider drift heating for dissociation
bbb.cfnidhmol = 1.0 #..ignore molecular v in dissociation drift heating

bbb.recyce = 0.3
bbb.recycm = -0.9