bbb.matwso[0] = 1		#recycle on main-chamber wall
bbb.isnwcono = 1		#if=1, set bbb.ni[,,com.ny+1]=bbb.nwallo
bbb.isnwconi = 1		#if=1, set PF wall bbb.ni=bbb.nwalli
bbb.allocate()		#bbb.allocate() space of bbb.nwallo,i
bbb.nwallo = 1.e18
bbb.nwalli = 1.e18


# gas pumping
bbb.nwsor = 2  #.. number of sources (sinks).  defalut: =1
  #.. PRF
bbb.igasi[0] = 0.   #.. pumping. default =0.
bbb.xgasi[0] = 0.   #.. position from inner. default =0.
bbb.wgasi[0] = 100. #.. width. default =100.(m)
bbb.igspsori[0] = 1 #.. species. default =1
bbb.albdsi[0] = 0.25 #.. pumping flux = (1.-albdsi)*1/4*n*sqrt(8T/pi/me)
  #.. wall
bbb.igaso[0] = 0.   #.. default =0.
bbb.xgaso[0] = 0.   #.. position from inner. default =0.
bbb.wgaso[0] = 100. #.. width. default =100.(m)
bbb.igspsoro[0] = 1 #.. species. default =1
bbb.albdso[0] = 0.25 #.. pumping flux = (1.-albdsi)*1/4*n*sqrt(8T/pi/me)

#..new pumping for atoms
  #.. PRF
bbb.igasi[1] = 0.   #.. pumping. default =0.
bbb.xgasi[1] = 0.   #.. position from inner. default =0.
bbb.wgasi[1] = 100. #.. width. default =100.(m)
bbb.igspsori[1] = 2 #.. species. default =1
bbb.albdsi[1] = 0.9999 #.. pumping flux = (1.-albdsi)*1/4*n*sqrt(8T/pi/me)
  #.. wall
bbb.igaso[1] = 0.   #.. default =0.
bbb.xgaso[1] = 0.   #.. position from inner. default =0.
bbb.wgaso[1] = 100. #.. width. default =100.(m)
bbb.igspsoro[1] = 2 #.. species. default =1
bbb.albdso[1] = 0.9999 #.. pumping flux = (1.-albdsi)*1/4*n*sqrt(8T/pi/me)
  #.. recycling
bbb.matwsi[1] = 1   #.. default =0
#bbb.recycw[0] = 0.9    #.. wall recycling coef. for igsp=1 (for PRF and wall)
bbb.matwso[1] = 1   #.. coefficient is recycw