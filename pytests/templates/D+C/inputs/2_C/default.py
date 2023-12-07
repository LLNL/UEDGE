# This part sets up carbon.
# Should be split into smaller parts, each part tested
## Impurity gas basics
com.ngsp = 2                #total number of gas species
bbb.isngon[1] = 1           #turns on impurity gas
bbb.ngbackg[1] = 1.e9	    #neutral impurity background for added source
bbb.ingb = 2		    #exponent for strength of ngbackg turn-on
bbb.istgcon[1] = 1	    #=1 for constant tg(2) at tgas(2)
bbb.tgas[1] = 1.	    #value for tg when istgcon=1
bbb.rcxighg = 0.            # best value; ratio of imp cx to hyd cx
bbb.kelighi[1] = 5.e-16     #elastic sig_v for imp_gas/h_ion
bbb.kelighg[1] = 5.e-16     #elastic sig_v for imp_gas/h_gas
bbb.n0g[1] = 1.e16          #imp. gas density normalization

# Impurity gas boundary conditions
bbb.recycp[1] = 0.01        #plate recycling of impurities
bbb.recycw[1] = 1e-4        #wall recycling; matwsi,o set above for hyd
bbb.isch_sput[1]=7          # Haasz/Davis chemical sputtering model
bbb.isph_sput[1]=3 	        # physical sputtering model

# Reduce C isputtering owing to boronization
bbb.fphysylb[1,0] = 0.5
bbb.fphysyrb[1,0] = 0.5
bbb.fchemylb[1,0] = 0.1
bbb.fchemyrb[1,0] = 0.1

bbb.t_wall = 300.
bbb.t_plat = 500.

## Impurity ions
bbb.isimpon = 6             #Use force-balance only
com.nzsp[0] = 6             #number chrg states impurity isotope #1

bbb.csfaclb[2:8,0] = 2.191
bbb.csfacrb[2:8,0] = 2.191
bbb.minu[2:8] = 12.
bbb.ziin[:2] = [1, 0]
bbb.ziin[2:8] = list(range(1,7))
bbb.znuclin[:2] = 1
bbb.znuclin[2:8] = 6
bbb.n0[2:8] = 1.e17

bbb.nzbackg = 1.e9          #background density for impurities
bbb.inzb = 2		    #exponent for switching on nzbackg
bbb.ismctab = 2             # use Braams" rate tables
com.mcfilename[0] = "C_rates.strahl"     # Imp rate file name

bbb.isnicore[7] = 3
bbb.curcore[7] = 0.

bbb.isnwcono[2:8] = 3
bbb.isnwconi[2:8] = 3
bbb.nwomin[2:8] = 1.e7
bbb.nwimin[2:8] = 1.e7

# Scale factor converting (upi-upg)**2 energy to thermal energy
bbb.cfnidh = 0.2