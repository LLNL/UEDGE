integer saveecho
saveecho = echo
echo = 0
#
# establish connection with dtm server/ncsa tool
#   port = init_dtm("alvin.llnl.gov:7777")
#
function init_dtm(server)
integer port

f_init_dtm(server,&port)
return port

endf
#
# terminate connection with dtm server/ncsa tool
#   
#

function close_dtm(port)
integer status

f_close_dtm(port,&status)
return status

endf


#
# re-establish connection with rpc server
# dtm ports are lost after this
#
function init_dce(server)
integer status
f_init_dce(server,&status)
return status
endf

#
# define rpc timeout in seconds
#
#
function dce_timeout(timeout)
integer status

f_dce_timeout(timeout,&status)
return status

endf

#
# run xgraph, see man page on xgraph for params
#
#
function xgraph($display,$title,x,$xunits,y,$yunits,$params)

integer dims(5)
integer xlen,ylen
integer status
if((length(myshape(x)) > 1) | (length(myshape(y)) > 1))then
	remark "xgraph only supports 1-D signals "
	return
endif
if((type(x) <> "real") | (type(y) <> "real"))then
	remark "xgraph only supports floating point arrays"
	return
endif

dims = myshape(x)
xlen = dims(1)
dims = myshape(y)
ylen = dims(1)
if(xlen .ne. ylen)then
	remark "Signals must be of equal length"
	return
endif

f_xgraph("-display "//$display//" -x """//$xunits//""""\
	  //" -y """//$yunits//""""//" -t """//$title//""" "//$params ,\
	  x,y,xlen,&status)

return status

endf

#
# send 2d data to an ncsa tool via dtm
#
#
function senddtm(port,$title,y)

integer dims(5)
integer xlen,ylen
integer status
if((length(myshape(y)) > 2) | (length(myshape(y)) < 2))then
	remark "senddtm only supports 2-D signals "
	return
endif
if(type(y) <> "real")then
	remark "senddtm only supports floating point arrays"
	return
endif

dims = myshape(y)
xlen = dims(1)
ylen = dims(2)

f_senddtm(port,$title,y,xlen,ylen,&status)

return status

endf

#
# transform ledge format data into r-z space and then
# send it to an ncsa tool
#

function senddtm_rzxform(port,$title,y,r,z,outxlen,outylen,rmin,rmax,zmin,zmax)

integer dims(5)
integer xlen,ylen
integer status
if((length(myshape(y)) > 2) | (length(myshape(y)) < 2))then
	remark "senddtm_rzxform only supports 2-D signals "
	return
endif
if(type(y) <> "real")then
	remark "senddtm_rzxform only supports floating point arrays"
	return
endif

dims = myshape(y)
xlen = dims(1)
ylen = dims(2)


f_senddtm_rzxform(port,$title,y,r,z,xlen,ylen,outxlen,outylen,rmin,rmax,zmin,zmax,&status)


return status
endf

#
# transform ledge format data into r-z space and then
# send it to an ncsa tool
#

function rzxform(y,r,z,outxlen,outylen,rmin,rmax,zmin,zmax)
integer dims(5)
integer xlen,ylen
integer status
real out(outxlen,outylen)
if((length(myshape(y)) > 2) | (length(myshape(y)) < 2))then
	remark "rzxform only supports 2-D signals on input "
	return
endif
if((length(myshape(out)) > 2) | (length(myshape(out)) < 2))then
	remark "rzxform only supports 2-D signals on output"
	return
endif
if(type(y) <> "real")then
	remark "senddtm_rzxform only supports floating point arrays"
	return
endif
if(type(out) <> "real")then
	remark "senddtm_rzxform only supports floating point arrays"
	return
endif

dims = myshape(y)
xlen = dims(1)
ylen = dims(2)

f_rzxform(&out,y,r,z,xlen,ylen,outxlen,outylen,rmin,rmax,zmin,\
zmax,&status)
return out

endf

function myshape(y)
integer dims(5),i
integer numdims(5),j
if(length(shape(y)) .eq. 1) return shape(y)

dims = shape(y)
j = 0;
do i = 1,5
	if(dims(i) .gt. 1)then
		j = j + 1
		numdims(j) = dims(i)
	endif
enddo
return numdims(1:j)
endf


#init_dce("sas-hp.nersc.gov")
init_dce("alvin.llnl.gov")
echo = saveecho

