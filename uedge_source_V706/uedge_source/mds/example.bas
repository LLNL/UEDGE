

echo = 1
integer sock = mdsconnect("atlas.gat.com")
chameleon ip = mdsvalue(sock,"ptdata('ip',110493)")
chameleon ip_time = mdsvalue(sock,"DIM_OF(ptdata('ip',110493))")
chameleon alt_ip1 = mdsvalue(sock,"ptdata('ip',$)",110493)
chameleon alt_ip2 = mdsvalue(sock,"ptdata($,110493)","ip")
chameleon alt_ip3 = mdsvalue(sock,"ptdata($,$)","ip",110493)

#
# EFIT01 is the inter-shot "jt" efit results
# EFIT02 is the inter-shot "mse" efit results
#
mdsopen(sock,"EFIT02",110493)
#
#
#
chameleon gnids  = mdsvalue(sock,"GETNCI($,'NID_NUMBER')","\TOP.RESULTS.GEQDSK:*:READA_NAME")
chameleon gparents  = mdsvalue(sock,"GETNCI(GETNCI($,'PARENT'),'FULLPATH')",gnids)
chameleon gtags  = mdsvalue(sock,"GETNCI($,'RECORD')",gnids)
#
# get psi(r,z) of equilibria (get index by looking at gtags)
#
chameleon psirz = mdsvalue(sock,gparents(16))
#
#
#
#
chameleon anids  = mdsvalue(sock,"GETNCI($,'NID_NUMBER')","\TOP.RESULTS.AEQDSK:*:READA_NAME")
chameleon aparents  = mdsvalue(sock,"GETNCI(GETNCI($,'PARENT'),'FULLPATH')",anids)
chameleon atags  = mdsvalue(sock,"GETNCI($,'RECORD')",anids)
#
# get time of equilibria (get index by looking at atags)
#
chameleon atime = mdsvalue(sock,aparents(4))
chameleon qmin = mdsvalue(sock,aparents(75))
mdsopen(sock,"EFIT02",110493)

#
# Thomson Data
#
mdsopen(sock,"ELECTRONS",110493)
chameleon date_loaded = mdsvalue(sock,"\TOP.TS.BLESSED.HEADER:DATE_LOADED")
mdssetdefault("\TOP.TS.BLESSED.CORE")
chameleon ts_r  = mdsvalue("R")


mdsclose()

mdsdisconnect()
