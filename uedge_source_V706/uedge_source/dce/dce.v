dce
# This is the post-processor for the UEDGE code.
# Impurity radiation diagnostic capability is included.
{
}

***** DCE_vars:
dce_server character*14 /'alvin.llnl.gov'/   # primary dce server hostname

***** Runroutines: 
dceinit subroutine # Initializes dce package using dce_server
f_init_dce(server:string,status:integer) subroutine	
f_dce_timeout(timeout:integer,status:integer) subroutine  
f_init_dtm(server:string,port:integer) subroutine  
f_close_dtm(port:integer,status:integer) subroutine  
f_senddtm(port:integer,title:string,x:real,xlen:integer,ylen:integer,status:integer) subroutine
f_xgraph(sparams:string,x:real,f:real,xlen:integer,status:integer) subroutine  
f_senddtm_rzxform(port:integer,title:string,f:real,r:real,z:real,\
      xlen:integer,ylen:integer,outxlen:integer,outylen:integer,\
      rmin:real,rmax:real,zmin:real,zmax:real,status:integer) subroutine
f_rzxform(out:real,f:real,r:real,z:real,xlen:integer,ylen:integer,
	  outxlen:integer,outylen:integer,rmin:real,rmax:real,
	   zmin:real,zmax:real,status:integer) subroutine
f_interprzxform(out:real,f:real,r:real,z:real,vint:real,rint:real,zint:real,\
		numint:integer,xlen:integer,ylen:integer,\
	  outxlen:integer,outylen:integer,rmin:real,rmax:real,\
	   zmin:real,zmax:real,status:integer) subroutine
f_interpdtm_rzxform(port:integer,title:string,f:real,r:real,z:real,\
			vint:real,rint:real,zint:real,numint:integer,\
      xlen:integer,ylen:integer,outxlen:integer,outylen:integer,\
      rmin:real,rmax:real,zmin:real,zmax:real,status:integer) subroutine
