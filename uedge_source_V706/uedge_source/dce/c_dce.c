
#include <rpc/rpc.h>
#include "basis_dce.h"
#include <stdio.h>
#include <time.h>

/*
 undef WELL_KNOWN_PORT to use portmapper only
 usually defined/(or not) in the PackageMakefile
*/

#ifdef WELL_KNOWN_PORT 
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netdb.h>
static int dcesock;
#endif


#define PORT 12100

static char current_server[256] = "alvin.llnl.gov";
static char revchk[] =  "@(#)c_dce.c	1.9\t4/10/95 Basis package DCE linked";

static CLIENT *cl = 0;
#ifdef linux
c_dceinit_(server,s)	/* 	package initialization     */
#endif
#ifdef sun
c_dceinit_(server,s)	/* 	package initialization     */
#endif
#ifdef __alpha
c_dceinit_(server,s)	/* 	package initialization     */
#endif
#ifdef sgi
c_dceinit_(server,s)	/* 	package initialization     */
#endif
#ifdef hpux
c_dceinit(server,s)	/* 	package initialization     */
#endif
#ifdef cray
C_DCEINIT(server,s)	/* 	package initialization     */
#endif
char *server;
int *s;	/* 	status return  */
{
	int to;
	dtmserver id;
	char *getlogin();
	int status;
#ifdef WELL_KNOWN_PORT
	struct sockaddr_in saddr;
	long addrlen;
	struct hostent *host;
#endif
	to = 60;	/* 	default timeout     */
	*s = 0;
	if(cl != NULL){
		my_clnt_sperror();
		clnt_destroy(cl);
		cl = NULL;
	}
        if(server)strncpy(current_server,server,256);



	printf("\nIf using xgraph make sure that host\n");
	printf("%s has access to specified Xservers\n\n",current_server);

#ifdef WELL_KNOWN_PORT
        host = gethostbyname(current_server);
	if(host != NULL){
           dcesock = socket(AF_INET,SOCK_STREAM,0);
	   if(dcesock > 0){
              saddr.sin_family = AF_INET;
              saddr.sin_addr.s_addr = INADDR_ANY;
              memcpy(&saddr.sin_addr,host->h_addr_list[0],host->h_length);
              saddr.sin_port = htons(PORT);
              addrlen = sizeof(struct sockaddr);
              if(connect(dcesock, (struct sockaddr *) &saddr, sizeof(saddr)) >= 0)
              {
                  cl = clnttcp_create(&saddr,DCEPROG,DCEVERS,&dcesock,0,0);
                  if(cl == NULL){ 
                    perror("DCE clnttcp to well known port:");
	            cl = clnt_create(current_server,DCEPROG,DCEVERS,"tcp");
                  }
	      } else {
	          perror("DCE connect to well known port:");
                  cl = clnt_create(current_server,DCEPROG,DCEVERS,"tcp");
              }
	   } else {
	      perror("DCE socket to well known port:");
              cl = clnt_create(current_server,DCEPROG,DCEVERS,"tcp");
           }
	} else {
	   perror("DCE gethostbyname to well known port:");
           cl = clnt_create(current_server,DCEPROG,DCEVERS,"tcp");
        }
#else
        cl = clnt_create(current_server,DCEPROG,DCEVERS,"tcp");

#endif




	if(cl == NULL){
		clnt_pcreateerror(current_server);
		return;
	} else {
#ifdef cray
		C_DCE_TIMEOUT(&to,&status);
#endif
#ifdef hpux
		c_dce_timeout(&to,&status);
#endif
#ifdef sun
		c_dce_timeout_(&to,&status);
#endif
#ifdef __alpha
		c_dce_timeout_(&to,&status);
#endif
#ifdef sgi
		c_dce_timeout_(&to,&status);
#endif
		*s = 1;
	}

}


/* cause rpc server to open dtm connection with specified dtm server     */
#ifdef cray
C_INIT_DTM(server,p,s)
#endif
#ifdef hpux
c_init_dtm(server,p,s)
#endif
#ifdef linux
c_init_dtm_(server,p,s)
#endif
#ifdef sun
c_init_dtm_(server,p,s)
#endif
#ifdef __alpha
c_init_dtm_(server,p,s)
#endif
#ifdef sgi
c_init_dtm_(server,p,s)
#endif
char *server;	/* 	host:port     */
int *p;	/* 	port return     */
{

	int *tmp;
	dtmserver rpcserver;
	*p = 0;	/* 	error return     */
#ifdef INIT_DTM
	if(cl != NULL){
#endif
		rpcserver.name = server;
		tmp = init_dtm_1(&rpcserver,cl);
		if(tmp == NULL){
			my_clnt_sperror();
			clnt_destroy(cl);
			cl = NULL;
		} else {
			*p = *tmp;
			if(*p == -1){
				fprintf(stderr,"DTM error\n");
				*p = 0;	/* 	dtm errors are -1     */
			}
		}
#ifdef INIT_DTM
	} else {
		fprintf(stderr,"Not connected to server %s\n",current_server);
	}
#endif
	return;
}

/* cause rpc server to destroy dtm connection with dtm server     */
#ifdef cray
C_CLOSE_DTM(p,s)
#endif
#ifdef hpux
c_close_dtm(p,s)
#endif
#ifdef sun
c_close_dtm_(p,s)
#endif
#ifdef linux
c_close_dtm_(p,s)
#endif
#ifdef __alpha
c_close_dtm_(p,s)
#endif
#ifdef sgi
c_close_dtm_(p,s)
#endif
int *p;	/* 	port from C_INIT_DTM     */
int *s;		/*    status return     */
{

	int *tmp;
	*s = 0;
#ifdef CLOSE_DTM
	if(cl != NULL){
#endif
		tmp = close_dtm_1(p,cl);
		if(tmp == NULL){
			my_clnt_sperror();
			clnt_destroy(cl);
			cl = NULL;
		} else if(*tmp == -1){
			fprintf(stderr,"DTM error\n");
			*s = 0;	/* 	dtm errors are -1     */
		} else{
			*s = 1;
		}
#ifdef CLOSE_DTM
	} else {
		fprintf(stderr,"Not connected to server %s\n",current_server);
	}
#endif
	return;
}
/* set rpc timeout      */
#ifdef cray
C_DCE_TIMEOUT(i,s)
#endif
#ifdef hpux
c_dce_timeout(i,s)
#endif
#ifdef linux
c_dce_timeout_(i,s)
#endif
#ifdef sun
c_dce_timeout_(i,s)
#endif
#ifdef __alpha
c_dce_timeout_(i,s)
#endif
#ifdef sgi
c_dce_timeout_(i,s)
#endif
int *i;	/* 	seconds     */
int *s;		/*    status return     */
{
	struct timeval timeout;
	*s = 0;
	if(cl){
		timeout.tv_sec = *i;
		timeout.tv_usec = 0;
		clnt_control(cl,CLSET_TIMEOUT,&timeout);
		*s = 1;
	} else {
		fprintf(stderr,"Not connected to server %s\n",current_server);
	}
}

/* send data to program that supports dtm on specified port    */
#ifdef cray
C_SENDDTM(p,t,f,ix,iy,s)
#endif
#ifdef hpux
c_senddtm(p,t,f,ix,iy,s)
#endif
#ifdef linux
c_senddtm_(p,t,f,ix,iy,s)
#endif
#ifdef sun
c_senddtm_(p,t,f,ix,iy,s)
#endif
#ifdef __alpha
c_senddtm_(p,t,f,ix,iy,s)
#endif
#ifdef sgi
c_senddtm_(p,t,f,ix,iy,s)
#endif
int *p;	/* 	port     */
char *t;	/* 	dtm title     */
float *f;	/* 	2d array      */
int *ix,*iy;	/* 	x and y dimensions of f     */
int *s;		/*    status return     */
{

	ledgeimage i;
	int *tmp;
	*s = 0;
#ifdef DTMOUT
	if(cl != NULL){
#endif
		i.dtmport = *p;
		i.title = t;
		i.imagedata.imagedata_val = f;
		i.imagedata.imagedata_len = *ix * *iy;
		i.xlen = *ix;
		i.ylen = *iy;
		tmp = dtmout_1(&i,cl);
		if(tmp == NULL){
			my_clnt_sperror();
			clnt_destroy(cl);
			cl = NULL;
		} else if(*tmp == -1){
			fprintf(stderr,"DTM error\n");
			*s = 0;	/* 	dtm errors are -1     */
		} else{
			*s = 1;
		}
#ifdef DTMOUT
	} else {
		fprintf(stderr,"Not connected to server %s\n",current_server);
	}
#endif

	return;
}
/* send ledge data to rpc server
   ledge data is an array of values and a set of polygon
   vertices that describe region in space that the values
   describe. The server does a polygon fill and then
   send the 2d image to an ncsa tool via dtm.

   */

#ifdef cray
C_SENDDTM_RZXFORM(p,t,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef hpux
c_senddtm_rzxform(p,t,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef linux
c_senddtm_rzxform_(p,t,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef sun
c_senddtm_rzxform_(p,t,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef __alpha
c_senddtm_rzxform_(p,t,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef sgi
c_senddtm_rzxform_(p,t,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
int *p;
char *t;
float *f,*r,*z,*rmin,*rmax,*zmin,*zmax;
int *ix,*iy,*ox,*oy;
int *s;		/*    status return     */
{

	ledge_pr_image i;
	int *tmp;
	*s = 0;
#ifdef DTMRZOUT
	if(cl != NULL){
#endif
		i.dtmport = *p;
		i.title = t;
		i.prdata.prdata_val = f;
		i.prdata.prdata_len = *ix * *iy;
		i.r_verts.r_verts_val = r;
		i.r_verts.r_verts_len = *ix * *iy * 5;
		i.z_verts.z_verts_val = z;
		i.z_verts.z_verts_len = *ix * *iy * 5;
		i.xlen = *ix;
		i.ylen = *iy;
		i.outxlen = *ox;
		i.outylen = *oy;
		i.rmin = *rmin;
		i.rmax = *rmax;
		i.zmin = *zmin;
		i.zmax = *zmax;
		tmp = dtmrzout_1(&i,cl);
		if(tmp == NULL){
			my_clnt_sperror();
			clnt_destroy(cl);
			cl = NULL;
		} else if(*tmp == -1){
			fprintf(stderr,"DTM error\n");
			*s = 0;	/* 	dtm errors are -1     */
		}  else {
			*s = 1;
		}
#ifdef DTMRZOUT
	} else {
		fprintf(stderr,"Not connected to server %s\n",current_server);
	}
#endif

	return;
}

#ifdef cray
C_INTERPDTM_RZXFORM(p,t,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef hpux
c_interpdtm_rzxform(p,t,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef linux
c_interpdtm_rzxform_(p,t,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef sun
c_interpdtm_rzxform_(p,t,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef __alpha
c_interpdtm_rzxform_(p,t,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef sgi
c_interpdtm_rzxform_(p,t,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
int *p;
char *t;
float *f,*r,*z,*rmin,*rmax,*zmin,*zmax;
float *rint,*zint,*vint;
int *nint,*ix,*iy,*ox,*oy;
int *s;		/*    status return     */
{

	ledge_pr_intimage i;
	int *tmp;
	*s = 0;
#ifdef DTMRZOUT_INT
	if(cl != NULL){
#endif
		i.dtmport = *p;
		i.title = t;
		i.prdata.prdata_val = f;
		i.prdata.prdata_len = *ix * *iy;
		i.r_verts.r_verts_val = r;
		i.r_verts.r_verts_len = *ix * *iy * 5;
		i.z_verts.z_verts_val = z;
		i.z_verts.z_verts_len = *ix * *iy * 5;
		i.r_int.r_int_val  = rint;
		i.r_int.r_int_len  = *nint;
		i.z_int.z_int_val  = zint;
		i.z_int.z_int_len  = *nint;
		i.v_int.v_int_val  = vint;
		i.v_int.v_int_len  = *nint;
		i.xlen = *ix;
		i.ylen = *iy;
		i.outxlen = *ox;
		i.outylen = *oy;
		i.rmin = *rmin;
		i.rmax = *rmax;
		i.zmin = *zmin;
		i.zmax = *zmax;
		tmp = dtmrzout_int_1(&i,cl);
		if(tmp == NULL){
			my_clnt_sperror();
			clnt_destroy(cl);
			cl = NULL;
		} else if(*tmp == -1){
			fprintf(stderr,"DTM error\n");
			*s = 0;	/* 	dtm errors are -1     */
		}  else {
			*s = 1;
		}
#ifdef DTMRZOUT_INT
	} else {
		fprintf(stderr,"Not connected to server %s\n",current_server);
	}
#endif

	return;
}
/* send 1d data to rpc server which then runs xgraph     */
#ifdef cray
C_XGRAPH(p,x,f,ix,s)
#endif
#ifdef hpux
c_xgraph(p,x,f,ix,s)
#endif
#ifdef linux
c_xgraph_(p,x,f,ix,s)
#endif
#ifdef sun
c_xgraph_(p,x,f,ix,s)
#endif
#ifdef __alpha
c_xgraph_(p,x,f,ix,s)
#endif
#ifdef sgi
c_xgraph_(p,x,f,ix,s)
#endif
char *p;	/* 	xgraph parameters     */
float *f,*x;	/* 	x and y axis arrays     */
int *ix;	/* 	number of data points     */
int *s;		/*    status return     */
{

	ledgesignal i;
	int *tmp;
	int j;
	j = 1;
	*s = 0;
#ifdef XGRAPHOUT
	if(cl != NULL){
#endif
		i.params = p;
		i.signaldata.signaldata_val = f;
		i.signaldata.signaldata_len = *ix;
		i.xdata.xdata_val = x;
		i.xdata.xdata_len = *ix;
		i.xlen = *ix;
		tmp = xgraphout_1(&i,cl);
		if(tmp == NULL){
			my_clnt_sperror();
			clnt_destroy(cl);
			cl = NULL;
		}  else {
			*s = 1;
		}
#ifdef XGRAPHOUT
	} else {
		fprintf(stderr,"Not connected to server %s\n",current_server);
	}
#endif

	return;
}
/* print out rpc error message     */
my_clnt_sperror()
{


	struct rpc_err err;

	clnt_geterr(cl,&err);
	clnt_perrno(err.re_status);
	fprintf(stderr,"\n");
}
/* send ledge data to rpc server
   ledge data is an array of values and a set of polygon
   vertices that describe region in space that the values
   describe. The server does a polygon fill and then
   send the 2d image to an ncsa tool via dtm.

   */

#ifdef cray
C_RZXFORM(o,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef hpux
c_rzxform(o,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef linux
c_rzxform_(o,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef sun
c_rzxform_(o,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef __alpha
c_rzxform_(o,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef sgi
c_rzxform_(o,f,r,z,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
float *o,*f,*r,*z,*rmin,*rmax,*zmin,*zmax;
int *ix,*iy,*ox,*oy;
int *s;		/*    status return     */
{

	ledge_pr_image i;
	rawimage *ret;
	int *tmp;
	*s = 0;
#ifdef RZXFORM
	if(cl != NULL){
#endif
		i.dtmport = 0;
		i.title = "none";
		i.prdata.prdata_val = f;
		i.prdata.prdata_len = *ix * *iy;
		i.r_verts.r_verts_val = r;
		i.r_verts.r_verts_len = *ix * *iy * 5;
		i.z_verts.z_verts_val = z;
		i.z_verts.z_verts_len = *ix * *iy * 5;
		i.xlen = *ix;
		i.ylen = *iy;
		i.outxlen = *ox;
		i.outylen = *oy;
		i.rmin = *rmin;
		i.rmax = *rmax;
		i.zmin = *zmin;
		i.zmax = *zmax;
		ret = rzxform_1(&i,cl);
		if(ret == NULL){
			my_clnt_sperror();
			clnt_destroy(cl);
			cl = NULL;
		}  else {
		    memcpy(o,ret->rawimage_val,ret->rawimage_len * sizeof(float));
		    *s = 1;
		}
#ifdef RZXFORM
	} else {
		fprintf(stderr,"Not connected to server %s\n",current_server);
	}
#endif

	return;
}

#ifdef cray
C_INTERPRZXFORM(o,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef hpux
c_interprzxform(o,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef linux
c_interprzxform_(o,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef sun
c_interprzxform_(o,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef __alpha
c_interprzxform_(o,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif
#ifdef sgi
c_interprzxform_(o,f,r,z,vint,rint,zint,nint,ix,iy,ox,oy,rmin,rmax,zmin,zmax,s)
#endif

float *o,*f,*r,*z,*rmin,*rmax,*zmin,*zmax;
float *rint,*zint,*vint;
int *nint,*ix,*iy,*ox,*oy;
int *s;		/*    status return     */
{

	ledge_pr_intimage i;
	rawimage *ret;
	int *tmp;
	*s = 0;
#ifdef RZXFORM_INT
	if(cl != NULL){
#endif
		i.dtmport = 0;
		i.title = "none";
		i.prdata.prdata_val = f;
		i.prdata.prdata_len = *ix * *iy;
		i.r_verts.r_verts_val = r;
		i.r_verts.r_verts_len = *ix * *iy * 5;
		i.z_verts.z_verts_val = z;
		i.z_verts.z_verts_len = *ix * *iy * 5;
		i.r_int.r_int_val  = rint;
		i.r_int.r_int_len  = *nint;
		i.z_int.z_int_val  = rint;
		i.z_int.z_int_len  = *nint;
		i.v_int.v_int_val  = vint;
		i.v_int.v_int_len  = *nint;
		i.xlen = *ix;
		i.ylen = *iy;
		i.outxlen = *ox;
		i.outylen = *oy;
		i.rmin = *rmin;
		i.rmax = *rmax;
		i.zmin = *zmin;
		i.zmax = *zmax;
		ret = rzxform_int_1(&i,cl);
		if(ret == NULL){
			my_clnt_sperror();
			clnt_destroy(cl);
			cl = NULL;
		}  else {
		    memcpy(o,ret->rawimage_val,ret->rawimage_len * sizeof(float));
		    *s = 1;
		}
#ifdef RZXFORM_INT
	} else {
		fprintf(stderr,"Not connected to server %s\n",current_server);
	}
#endif

	return;
}
