# Currents and potential parameters
bbb.b0 = 1.		        # =1 for normal direction of B-field
bbb.rsigpl=1.e-9		#anomalous cross-field conductivity
bbb.cfjhf=1.	   	    #turn-on heat flow from current [bbb.fqp]
bbb.cfjve = 1.		    #makes bbb.vex = vix - bbb.cfjve*bbb.fqx
bbb.cfjpy = 0.		    #diamag. cur. in flx.y-direction
bbb.cfjp2 = 0.		    #diamag. cur. in 2-direction

bbb.newbcl = 1		    #Sheath BC [bee,i] from current equation
bbb.newbcr = 1		    #Sheath BC [bee,i] from current equation
bbb.isfdiax =0.		    #Factor to turn on diamag. contrib. to sheath

bbb.cfyef = 0.6	    	#ExB drift in flx.y-dir.
bbb.cf2ef = 0.6		    #ExB drift in 2-dir.
bbb.cfybf = 1.0e-20     #turns on bbb.vycb - radial grad_B drift
bbb.cf2bf = 1.0e-20     #turns on bbb.v2cb - perp grad_B drift [nearly pol]
bbb.cfqybf = 1.0e-20    #turns on bbb.vycb contrib to radial current
bbb.cfq2bf = 1.0e-20    #turns on bbb.v2cb contrib to perp["2"] current
bbb.cftef = 1.0e-20		#turns on bbb.v2ce for toroidal velocity
bbb.cftdd = 1.0e-20		#turns on bbb.v2dd [diagm vel] for toroidal vel
bbb.cfqym = 1.0e-20		#turns on inertial correction to bbb.fqy current
bbb.cfydd = 0.		    #Diamag. drift in flx.y-dir. [always=0]
bbb.cf2dd = 0.		    #Diamag. drift in 2-dir. [always=0]
bbb.cfrd = 0.		    #Resistive drift in flx.y and 2 dirs.
bbb.cfbgt = 0.	    	#Diamag. energy drift [always=0]
bbb.isnewpot = 0
bbb.rnewpot = 0.