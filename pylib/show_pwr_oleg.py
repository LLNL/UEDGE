#-show power fluxes in the outer divertor leg domain


pwrx=bbb.feex+bbb.feix
pwry=bbb.feey+bbb.feiy

#-calculate photon fluxes on material surfaces
bbb.pradpltwl()



print "Convection and conduction power entering outer leg = ", \
    sum(pwrx[com.ixpt2[0],com.iysptrx+1:com.ny])/1e3, " [kW]"

print "Recombination power entering outer leg = ", \
    sum((bbb.fnix[com.ixpt2[0],com.iysptrx+1:com.ny,0])*bbb.ebind*bbb.ev)/1e3, " [kW]"

print ""



print "Convection and conduction power to outer target plate = ", \
    sum(pwrx[com.nx,com.iysptrx+1:com.ny])/1e3, " [kW]"

print "Recombination power to outer target plate = ", \
    sum((bbb.fnix[com.nx,com.iysptrx+1:com.ny,0])*bbb.ebind*bbb.ev)/1e3, " [kW]"

print "Hydrogen photon power to outer target plate = ", \
    sum(bbb.pwr_plth[:,1]*com.sx[com.nx,:])/1e3, " [kW]"

print "Impurity photon power to outer target plate = ", \
    sum(bbb.pwr_pltz[:,1]*com.sx[com.nx,:])/1e3, " [kW]"

print ""


print "Hydrogen radiated power lost in outer leg = ", \
    sum((bbb.erliz+bbb.erlrc)[com.ixpt2[0]:com.nx,0:com.ny])/1e3, " [kW]"

if (bbb.isimpon > 0):
    imprad=sum((bbb.prad*com.vol)[com.ixpt2[0]:com.nx,0:com.ny])/1e3, " [kW]"
else:
    imprad=0.0
print "Impurity radiated power lost in outer leg = ", imprad

print ""



print "Convection and conduction power to outer leg outer side = ", \
    sum(pwry[com.ixpt2[0]:com.nx,com.ny])/1e3, " [kW]"

print "Recombination power to outer leg outer side = ", \
    sum((bbb.fniy[com.ixpt2[0]:com.nx,com.ny,0])*bbb.ebind*bbb.ev)/1e3, " [kW]"

print "Hydrogen photon power to outer leg outer side = ", \
    sum(bbb.pwr_wallh[com.ixpt2[0]:com.nx]*com.sy[com.ixpt2[0]:com.nx,com.ny])/1e3, " [kW]"

print "Impurity photon power to outer leg outer side = ", \
    sum(bbb.pwr_wallz[com.ixpt2[0]:com.nx]*com.sy[com.ixpt2[0]:com.nx,com.ny])/1e3, " [kW]"

print ""


print "Convection and conduction power to outer leg inner side = ", \
    sum(pwry[com.ixpt2[0]:com.nx,0])/1e3, " [kW]"

print "Recombination power to outer leg inner side = ", \
    sum((bbb.fniy[com.ixpt2[0]:com.nx,0,0])*bbb.ebind*bbb.ev)/1e3, " [kW]"
