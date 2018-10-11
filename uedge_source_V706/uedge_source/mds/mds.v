mds
#skeleton variable descripion file written by basiskit
#comments look like this
#
#
#
COMMENTSIZE = 72
DIMSTRSIZE = 500
UNITSIZE = 16
*******MDSVars:
# put your variables here in place of these
mdsplus_socket integer 
mdsplus_error integer 
mdsplus_diiid_server character*13 /'atlas.gat.com'/

*******MDSRoutines:
# describe your routines here like this
mdsvalue()  builtin [1-100]
mdswfevent()  builtin [0-100]
mdssetupevent()  builtin [1-100]
mdsput()  builtin [1-100]
mdsconnect(server:string) integer function
mdsdisconnect() builtin [0-100]
mdsopen() builtin [1-100]
mdssetdefault() builtin [1-100]
mdsclose() builtin [0-100]
mdslogin(s:integer,user:string,passwd:string)  builtin [2-100]
mdssave(attr:string)	subroutine
mdssavevar(varname:string;pkg:string) integer function
mdscopyvar(svarname:string,spkg:string,dvarname:string,dpkg:string;ilow:integer,ihi:integer,icol:integer) integer function


*******MDSRestored:
