#define MDSPLUS_FILE_VERSION 1.1

#define DEFAULT -999999
#include <stdio.h>
#include <score.h>
#include <Basis.h>

/* 
   This package makes use of some basis entry points which
   have been changing as basis converts from fortran to C.
   The following ifdefs allows this source to compile with
   different versions of basis. 
*/
#if !defined(BASIS_MAJOR_VERSION) || BASIS_MAJOR_VERSION < 12
#error "MDS package only compilable starting with basis 12"
#elif BASIS_MAJOR_VERSION == 12 && BASIS_MINOR_VERSION == 1 && \
      BASIS_MICRO_VERSION <= 999999
#define BASIS_ENTRIES_1
#elif BASIS_MAJOR_VERSION == 12 && BASIS_MINOR_VERSION == 2 && \
      BASIS_MICRO_VERSION > 1
#define BASIS_ENTRIES_2
#else
#error "I do not know what to do with this version of basis"
#endif


/*#include <ipdesc.h>*/
#include <mdsip.h>
#include <time.h>
#include <par_.h>
#ifndef __hpux
#define setmdserror setmdserror_
#define setmdssocket setmdssocket_
#define getmdssocket getmdssocket_
#define getmdsfileversion getmdsfileversion_
#endif
/* 	basis types  */
#include <par_.h>
/*
#define NULL_TYPE  0
#define INTEGER  1
#define REAL  2
#define DOUBLE  3
#define COMPLEX  4
#define DOUBLE_C  5
#define LOGICAL  6
#define RANGE  7
#define STRUCTURE  8
#define GROUP  9
#define NAME  10
#define FUNCTION_TYPE  11
#define LHS12
#define FWA  13
#define STRING_ADDRESS  14
#define LASTTYPE = STRING_ADDRESS
*/
/* TYPE IDENTIFIERS FOR MDSPLUS DESCRIPTORS */
#define IDTYPE_UCHAR 2
#define IDTYPE_USHORT 3
#define IDTYPE_ULONG 4
#define IDTYPE_ULONGLONG 5
#define IDTYPE_CHAR 6
#define IDTYPE_SHORT 7
#define IDTYPE_LONG 8
#define IDTYPE_LONGLONG 9
#define IDTYPE_FLOAT 10
#define IDTYPE_DOUBLE 11
#define IDTYPE_COMPLEX 12
#define IDTYPE_COMPLEX_DOUBLE 13
#define IDTYPE_CSTRING 14
#define OK 1
#define TRUE 1
#define MAXSTRING 500
#define MAXDIMS 7
extern struct descrip *MakeDescrip(struct descrip *in_descrip, char dtype, char ndims, int *dims, void *ptr);
extern struct descrip *MakeDescripWithLength(struct descrip *in_descrip, char dtype, int length,char ndims, int *dims, void *ptr);




typedef struct s_f_entry f_entry;

typedef Integer (*PFHand)(Integer nargs, char *fname,BA_stack_elem *sxi);
struct s_f_entry
{char *fname;                /* Actual name of the function  */
  Integer arg1typ;            /* Allowable types for arg 1    */
  PFHand handler;
  Integer ni;
  void **info;};


Integer Mdsbfcn(Integer nargs,char *fname, BA_stack_elem *sx)
{

  static char msg1[] = "Mds: Function not implemented yet.";
  int status;
  Integer Mdsvalue();
  Integer Mdsput();
  Integer Mdsopen();
  Integer Mdsclose();
  Integer Mdswfevent();
  Integer Mdssetupevent();
  Integer Mdssetdefault();
  Integer Mdslogin();
  Integer Mdsdisconnect();



  sx->ndim = 0;
  sx->type_code = 0;
  status = OK;
  /*printf("mdsbfcn called:%s \n",fname);*/
  if(strncmp(fname,"mdsvalue",7) == 0){
    Mdsvalue(nargs,fname,sx);
  }else if(strncmp(fname,"mdsput",6) == 0){
    Mdsput(nargs,fname,sx);
  }else if(strncmp(fname,"mdsopen",7) == 0){
    Mdsopen(nargs,fname,sx);
  }else if(strncmp(fname,"mdsclose",8) == 0){
    Mdsclose(nargs,fname,sx);
  }else if(strncmp(fname,"mdslogin",8) == 0){
    Mdslogin(nargs,fname,sx);
  }else if(strncmp(fname,"mdswfevent",10) == 0){
    Mdswfevent(nargs,fname,sx);
  }else if(strncmp(fname,"mdssetupevent",13) == 0){
    Mdssetupevent(nargs,fname,sx);
  }else if(strncmp(fname,"mdsdisconnect",13) == 0){
    Mdsdisconnect(nargs,fname,sx);
  }else if(strncmp(fname,"mdssetdefault",13) == 0){
    Mdssetdefault(nargs,fname,sx);
  } else {
    BasError(msg1); 
    BasError(fname); 
    status = ERR;
  }
  
  return status;


}

Integer setstatusreturnp(BA_stack_elem *sxi,int **ix)
{
  sxi->ndim = 0;
  sxi->extent[0] = 1;
  sxi->lower[0] = 1;
  sxi->stride[0] = 1;
  sxi->type_code = INTEGER;
  Parget(sxi);
  *ix = (int *) Paraddr(sxi);
}
Integer setstatusreturn(BA_stack_elem *sxi,int status)
{

  int *ix;
  setstatusreturnp(sxi,&ix);
  *ix = status;
}

/* 
   Connecto to tcp server.
   For DIII-D as of this writing the server is atlas.gat.com.
*/
Integer Mdsconnect(char *server)
{
  int status;

  status = ConnectToMds(server);
  setmdssocket(&status);
  setmdserror(&status);
  return(status);
}


/* 
MDSPlus is hierarchical, to access data the proper "tree" must
be opened.
*/
Integer Mdsopen(Integer nargs, char *fname, BA_stack_elem *sxi)
{
  BA_stack_elem syi_,*syi=&syi_;
  char errmsg[1024];
  Integer sock;
  Integer shot;
  Integer argind;
  char tree[1025];
  int i,j,n;
  int *status,mdsstatus;
  int getmdssocket();

     
  
  setstatusreturnp(sxi,&status);

  *status = OK;
  argind = 1;
     
  /* get first argument which may be the socket returned by mdsconnect  */
  Fcnargb(argind++,syi);
  if (syi->type_code != INTEGER){
    getmdssocket(&sock); 
  } else {
    sock = *(Integer *)Paraddr(syi);
    Parrel(syi);
    Fcnargb(argind++,syi);
  }

  if (sock == 0){
    sprintf(errmsg,"%s: sock = 0, use mdsconnect to connet to server.",fname);
    BasError(errmsg);
    *status =  ERR;
    setmdserror(status);
    return(*status);
  }

  if (Parlen(syi) != 1){
    sprintf(errmsg,"%s: tree expected.",fname);
    BasError(errmsg);
    *status =  ERR;
    setmdserror(status);
    return(*status);
  }
  memset(tree,0L,sizeof(tree));
  strncpy(tree,(char *)Paraddr(syi),((sizeof(tree)-1) < -syi->type_code)?(sizeof(tree)-1):-syi->type_code);
  /*printf("sock %d:    tree %s\n",sock,tree);*/
  Parrel(syi);

  Fcnargb(argind++,syi);
  if (syi->type_code != INTEGER){
    sprintf(errmsg,"%s: argument %d has illegal type. Must be integer shot number.",fname,argind-1);
    BasError(errmsg);
    *status = ERR;
    setmdserror(status);
    return(*status);
  }
  shot = *(Integer *)Paraddr(syi);
  Parrel(syi);
  /*printf("open tree %s using socket %d on shot %d\n",tree,sock,shot);*/
     

  *status = MdsOpen(sock,tree,shot);
  setmdserror(status);
  return(*status);
}

/* 
This routine closes the opened tree
*/
Integer Mdsclose(Integer nargs, char *fname, BA_stack_elem *sxi)
{
  BA_stack_elem syi_,*syi=&syi_;
  char errmsg[1024];
  Integer sock;
  Integer shot;
  Integer argind;
  int i,j,n;
  int *status,mdsstatus;
  int getmdssocket();

     

  setstatusreturnp(sxi,&status);


  *status = OK;
  argind = 1;
     
  /* get first argument which may be the socket returned by mdsconnect  */
  if(argind <= nargs){
     Fcnargb(argind++,syi);
     if (syi->type_code != INTEGER){
       getmdssocket(&sock); 
     } else {
       sock = *(Integer *)Paraddr(syi);
       Parrel(syi);
       Fcnargb(argind++,syi);
     }

     if (sock == 0){
       sprintf(errmsg,"%s: sock = 0, use mdsconnect to connet to server.",fname);
       BasError(errmsg);
       *status =  ERR;
       setmdserror(status);
       return(*status);
     }
  } else {
     getmdssocket(&sock); 
  }

  if (sock == 0){
     sprintf(errmsg,"%s: sock = 0, use mdsconnect to connet to server.",fname);
     BasError(errmsg);
     *status =  ERR;
     setmdserror(status);
     return(*status);
  }


  *status = MdsClose(sock);
  setmdserror(status);
  return(*status);
}

/* 
This routine closes the socket connection 
*/
Integer Mdsdisconnect(Integer nargs, char *fname, BA_stack_elem *sxi)
{
  BA_stack_elem syi_,*syi=&syi_;
  char errmsg[1024];
  Integer sock,dsock;
  Integer shot;
  Integer argind;
  int i,j,n;
  int *status,mdsstatus;
  int getmdssocket();

     

  setstatusreturnp(sxi,&status);


  *status = OK;
  argind = 1;
     
  if(argind <= nargs){
  /* get first argument which may be the socket returned by mdsconnect  */
     Fcnargb(argind++,syi);
     if (syi->type_code != INTEGER){
       getmdssocket(&sock); 
     } else {
       sock = *(Integer *)Paraddr(syi);
       Parrel(syi);
       Fcnargb(argind++,syi);
     }
  } else {
     getmdssocket(&sock); 
  }

  if (sock == 0){
    sprintf(errmsg,"%s: sock = 0, use mdsconnect to connet to server.",fname);
    BasError(errmsg);
    *status =  ERR;
    setmdserror(status);
    return(*status);
  }

  *status = DisconnectFromMds(sock);
  getmdssocket(&dsock); 
  if ( dsock == sock){
     sock = 0;
     setmdssocket(&sock);
  } else {
     sock = 0;
  }
  setmdserror(status);
  return(*status);
}

/* 
MDSPlus is hierarchical, this sets the default tree, something like unix cd
*/
Integer Mdssetdefault(Integer nargs, char *fname, BA_stack_elem *sxi)
{
  BA_stack_elem syi_,*syi=&syi_;
  char errmsg[1024];
  Integer sock;
  Integer shot;
  Integer argind;
  char tree[1025];
  int i,j,n;
  int *status,mdsstatus;
  int getmdssocket();

  setstatusreturnp(sxi,&status);


  *status = OK;
  argind = 1;
     
  /* get first argument which may be the socket returned by mdsconnect  */
  Fcnargb(argind++,syi);
  if (syi->type_code != INTEGER){
    getmdssocket(&sock); 
  } else {
    sock = *(Integer *)Paraddr(syi);
    Parrel(syi);
    Fcnargb(argind++,syi);
  }

  if (sock == 0){
    sprintf(errmsg,"%s: sock = 0, use mdsconnect to connet to server.",fname);
    BasError(errmsg);
    *status =  ERR;
    setmdserror(status);
    return(*status);
  }

  if (Parlen(syi) != 1){
    sprintf(errmsg,"%s: tree expected.",fname);
    BasError(errmsg);
    *status =  ERR;
    setmdserror(status);
    return(*status);
  }
  memset(tree,0L,sizeof(tree));
  strncpy(tree,(char *)Paraddr(syi),((sizeof(tree)-1) < -syi->type_code)?(sizeof(tree)-1):-syi->type_code);
  /*printf("sock %d:    tree %s\n",sock,tree);*/
  Parrel(syi);


  *status = MdsSetDefault(sock,tree);
  setmdserror(status);
  return(*status);
}


/* 
This is the main data access routine for MDSPlus.
*/
Integer Mdsvalue(Integer nargs, char *fname, BA_stack_elem *sxi)
{
  BA_stack_elem syi_,*syi=&syi_;
  BA_stack_elem mdsarg[256];
  struct descrip *mdsdarg[256];
  struct descrip ans;
  float *rx,*cx,*ry,*cy;
  double *mdsdrx,*drx;
  float *mdsrx,*mdscx,*mdsry,*mdscy;
  real8 *dx, *dy, *dcx, *dcy;
  Integer npoints,*ix,*iy,dim;
  Integer *mdsix,*mdsiy;
  unsigned short *usmdsix;
  short *smdsix;
  char errmsg[1024];
  Integer sock;
  char com[1025];
  char *chx;
  unsigned char *uchx;
  int i,j,n;
  int status,mdsstatus;
  Integer argind;
  int getmdssocket();

     
  sxi->ndim = 0;
  sxi->extent[0] = 1;
  sxi->lower[0] = 1;
  sxi->stride[0] = 1;
  sxi->type_code = INTEGER;


  status = OK;
  /*printf("Mdsvalue called with %d args\n",nargs);*/
     
  /* get first argument which may be the socket returned by mdsconnect  */
  Fcnargb(1,syi);
  if (syi->type_code != INTEGER){
    getmdssocket(&sock); 
    /*printf("first argument not socket, using socket %d\n",sock);*/
    argind = 2;
  } else {
    sock = *(Integer *)Paraddr(syi);
    Parrel(syi);
    Fcnargb(2,syi);
    argind = 3;
  }

  if (sock == 0){
    sprintf(errmsg,"%s: sock = 0, use mdsconnect to connet to server.",fname);
    BasError(errmsg);
    status =  ERR;
    setstatusreturn(sxi,status);
    setmdserror(&status);
    return(status);
  }

  if (Parlen(syi) != 1){
    sprintf(errmsg,"%s: tdi expression expected.",fname);
    BasError(errmsg);
    status =  ERR;
    setstatusreturn(sxi,status);
    setmdserror(&status);
    return(status);
  }
  memset(com,0L,sizeof(com));
  strncpy(com,(char *)Paraddr(syi),((sizeof(com)-1) < -syi->type_code)?(sizeof(com)-1):-syi->type_code);
  /*printf("sock %d:    comm %s\n",sock,com);*/
  Parrel(syi);

  for(i = 0; i < 256; ++i) mdsdarg[i] = NULL;
  for(i = argind,n=0; i <= nargs; ++i,++n){
    Fcnargb(i,&mdsarg[n]);
    for(npoints = 1, j = 0;j<mdsarg[n].ndim;++j){
      npoints *= mdsarg[n].extent[j];
    }
    /*printf("arg %d, input arg type = %d, npoints=%d\n",i,mdsarg[n].type_code,npoints);*/
    /* convert input arguments into mds descriptors  */
    if(mdsarg[n].type_code >= 0){
      switch(mdsarg[n].type_code){
      case DOUBLE:
	if(npoints > 1)
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_DOUBLE,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	else
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_DOUBLE,NULL,NULL,Paraddr(&mdsarg[n]));
	mdsdarg[n]->length= ArgLen(mdsdarg[n]);
	break;
      case REAL:
	if(npoints > 1)
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_FLOAT,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	else
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_FLOAT,NULL,NULL,Paraddr(&mdsarg[n]));
	mdsdarg[n]->length= ArgLen(mdsdarg[n]);
	break;
      case COMPLEX:
	if(npoints > 1)
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_COMPLEX,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	else
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_COMPLEX,NULL,NULL,Paraddr(&mdsarg[n]));
	mdsdarg[n]->length= ArgLen(mdsdarg[n]);
	break;
      case DOUBLE_C:
	if(npoints > 1)
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_COMPLEX_DOUBLE,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	else
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_COMPLEX_DOUBLE,NULL,NULL,Paraddr(&mdsarg[n]));
	mdsdarg[n]->length= ArgLen(mdsdarg[n]);
	break;
      case INTEGER:
      case LOGICAL:
	if(npoints > 1)
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_LONG,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	else
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_LONG,NULL,NULL,Paraddr(&mdsarg[n]));
	/*printf("%d\n",*(int *)mdsdarg[n]->ptr);*/
	mdsdarg[n]->length= ArgLen(mdsdarg[n]);
	break;
      }
    } else if(mdsarg[n].type_code < 0){
	if(npoints > 1){
      mdsdarg[n] = MakeDescripWithLength(calloc(1,sizeof(struct descrip)),
	 IDTYPE_CSTRING,-mdsarg[n].type_code,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	} else
      mdsdarg[n] = MakeDescripWithLength(calloc(1,sizeof(struct descrip)),
	 IDTYPE_CSTRING,-mdsarg[n].type_code,NULL,NULL,Paraddr(&mdsarg[n]));
    }
  }
  memset(&ans,0L,sizeof(struct descrip));
  mdsdarg[nargs-argind+1] = &ans;
  mdsdarg[nargs-argind+2] = NULL;

     
  /*fprintf(stderr,"command = %s ans in %d\n",com,nargs-argind+1);*/
  if((mdsstatus = MdsValue(sock,com,mdsdarg[0],mdsdarg[1],mdsdarg[2],
			   mdsdarg[3],mdsdarg[4],mdsdarg[5],mdsdarg[6],mdsdarg[7],mdsdarg[8],mdsdarg[9],mdsdarg[10],
			   mdsdarg[11],mdsdarg[12],mdsdarg[13],mdsdarg[14],mdsdarg[15],mdsdarg[16],mdsdarg[17],mdsdarg[18],mdsdarg[19],mdsdarg[20],
			   mdsdarg[21],mdsdarg[22],mdsdarg[23],mdsdarg[24],mdsdarg[25],mdsdarg[26],mdsdarg[27],mdsdarg[28],mdsdarg[29],mdsdarg[30],
			   mdsdarg[31],mdsdarg[32],mdsdarg[33],mdsdarg[34],mdsdarg[35],mdsdarg[36],mdsdarg[37],mdsdarg[38],mdsdarg[39],mdsdarg[40],
			   NULL)) & 1){
    sxi->ndim = ans.ndims;
    for(npoints = 1,i = 0; i < sxi->ndim; ++i){
      sxi->lower[i] = 1;
      sxi->extent[i] = ans.dims[i];
      sxi->stride[i] = 1;
      npoints *= ans.dims[i];
    }

    /*fprintf(stderr,"now return the results\n");*/
    /* translate mds descriptor returned into a basis data type  */
    switch (ans.dtype){
    case IDTYPE_DOUBLE:
      sxi->type_code = DOUBLE;
      Parget(sxi);
      memcpy(Paraddr(sxi),ans.ptr,npoints*sizeof(double));
      break;
    case IDTYPE_FLOAT:
      sxi->type_code = REAL;
      Parget(sxi);
      memcpy(Paraddr(sxi),ans.ptr,npoints*sizeof(float));
      break;
    case IDTYPE_COMPLEX:
      sxi->type_code = COMPLEX;
      Parget(sxi);
      memcpy(Paraddr(sxi),ans.ptr,npoints*sizeof(float)*2);
      break;
    case IDTYPE_COMPLEX_DOUBLE:
      sxi->type_code = DOUBLE_C;
      Parget(sxi);
      memcpy(Paraddr(sxi),ans.ptr,npoints*sizeof(double)*2);
      break;
    case IDTYPE_LONG:
      sxi->type_code = INTEGER;
      Parget(sxi);
      memcpy(Paraddr(sxi),ans.ptr,npoints*sizeof(int));
      break;
    case IDTYPE_SHORT:
      sxi->type_code = INTEGER;
      Parget(sxi);
      ix = (int *) Paraddr(sxi);
      smdsix = (short *) ans.ptr;
      for(i = 0; i < npoints; ++i){
	ix[i] = smdsix[i];
      }
      break;
    case IDTYPE_USHORT:
      sxi->type_code = INTEGER;
      Parget(sxi);
      ix = (int *) Paraddr(sxi);
      usmdsix = (unsigned short *) ans.ptr;
      for(i = 0; i < npoints; ++i){
	ix[i] = usmdsix[i];
      }
      break;
    case IDTYPE_CHAR:
      sxi->type_code = INTEGER;
      Parget(sxi);
      ix = (int *) Paraddr(sxi);
      chx = (char *) ans.ptr;
      for(i = 0; i < npoints; ++i){
	ix[i] = chx[i];
      }
      break;
    case IDTYPE_UCHAR:
      sxi->type_code = INTEGER;
      Parget(sxi);
      ix = (int *) Paraddr(sxi);
      uchx = (char *) ans.ptr;
      for(i = 0; i < npoints; ++i){
	ix[i] = uchx[i];
      }
      break;
    case IDTYPE_CSTRING:
      sxi->type_code = -ans.length;
      Parget(sxi);
      memcpy(Paraddr(sxi),ans.ptr,ans.length*npoints);
      break;
    default:
      Parget(sxi);
      ix = (int *) Paraddr(sxi);
      sprintf(errmsg,"%s: Unsupported data type returned = %d.",fname,ans.dtype);
      BasError(errmsg);
      status = ERR;
      *ix = ERR;
      break;
    }

  } else {
    Parget(sxi);
    ix = (int *) Paraddr(sxi);
    sprintf(errmsg,"%s: MdsValue error = %s.",fname,ans.ptr);
    BasError(errmsg);
    status = ERR;
    *ix = ERR;
  }
  /*fprintf(stderr,"now free all the args\n");*/
  for(i = argind,n = 0; i <= nargs; ++i,++n){
    free(mdsdarg[n]);
    Parrel(&mdsarg[n]);
  }

  setmdserror(&status);
  return status;
}
/* 
This is the main data put routine for MDSPlus.
*/
Integer Mdsput(Integer nargs, char *fname, BA_stack_elem *sxi)
{
  BA_stack_elem syi_,*syi=&syi_;
  BA_stack_elem mdsarg[256];
  struct descrip *mdsdarg[256];
  struct descrip ans;
  float *rx,*cx,*ry,*cy;
  float *mdsrx,*mdscx,*mdsry,*mdscy;
  real8 *dx, *dy, *dcx, *dcy;
  Integer npoints,*ix,*iy,dim;
  Integer *mdsix,*mdsiy;
  char errmsg[1024];
  Integer sock;
  Integer argind;
  char com[1025],node[2049];
  int i,j,n;
  int status,mdsstatus;

     

  sxi->ndim = 0;
  sxi->extent[0] = 1;
  sxi->lower[0] = 1;
  sxi->stride[0] = 1;
  sxi->type_code = INTEGER;

  status = OK;
  /*printf("MdsPut called with %d args\n",nargs);*/
     
  /* get first argument which may be the socket returned by mdsconnect  */
  Fcnargb(1,syi);
  if (syi->type_code != INTEGER){
    getmdssocket(&sock); 
    /*printf("first argument not socket, using socket %d\n",sock);*/
    argind = 2;
  } else {
    sock = *(Integer *)Paraddr(syi);
    Parrel(syi);
    Fcnargb(2,syi);
    argind = 3;
  }

  if (sock == 0){
    sprintf(errmsg,"%s: sock = 0, use mdsconnect to connet to server.",fname);
    BasError(errmsg);
    status =  ERR;
    setstatusreturn(sxi,status);
    setmdserror(&status);
    return(status);
  }

  /* get second argument which is the tree node */
  if (Parlen(syi) != 1){
    sprintf(errmsg,"%s: tree node expected.",fname);
    BasError(errmsg);
    status =  ERR;
    setstatusreturn(sxi,status);
    setmdserror(&status);
    return(status);
  }
  memset(node,0L,sizeof(node));
  strncpy(node,(char *)Paraddr(syi),((sizeof(node)-1) < -syi->type_code)?(sizeof(node)-1):-syi->type_code);
  Parrel(syi);

  /* 	this is a valid mdsplus data command   */
  Fcnargb(argind,syi);
  if (Parlen(syi) != 1){
    sprintf(errmsg,"%s: tdi expression expected.",fname);
    BasError(errmsg);
    status =  ERR;
    setstatusreturn(sxi,status);
    setmdserror(&status);
    return(status);
  }
  memset(com,0L,sizeof(com));
  strncpy(com,(char *)Paraddr(syi),((sizeof(com)-1) < -syi->type_code)?(sizeof(com)-1):-syi->type_code);
  Parrel(syi);
  ++argind;

  for(i = 0; i < 256; ++i) mdsdarg[i] = NULL;
  for(i = argind,n=0; i <= nargs; ++i,++n){
    Fcnargb(i,&mdsarg[n]);
    for(npoints = 1, j = 0;j<mdsarg[n].ndim;++j){
      npoints *= mdsarg[n].extent[j];
    }
    /*printf("arg %d, input arg type = %d, npoints=%d\n",i,mdsarg[n].type_code,npoints);*/
    /* convert input arguments into mds descriptors  */
    if(mdsarg[n].type_code >= 0){
      /*printf("input type=%d\n",mdsarg[n].type_code);*/
      switch(mdsarg[n].type_code){
      case DOUBLE:
	if(npoints > 1)
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_DOUBLE,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	else
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_DOUBLE,NULL,NULL,Paraddr(&mdsarg[n]));
	mdsdarg[n]->length= ArgLen(mdsdarg[n]);
	break;
      case REAL:
	if(npoints > 1)
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_FLOAT,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	else
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_FLOAT,NULL,NULL,Paraddr(&mdsarg[n]));
	mdsdarg[n]->length= ArgLen(mdsdarg[n]);
	break;
      case COMPLEX:
	if(npoints > 1)
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_COMPLEX,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	else
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_COMPLEX,NULL,NULL,Paraddr(&mdsarg[n]));
	mdsdarg[n]->length= ArgLen(mdsdarg[n]);
	break;
      case DOUBLE_C:
	if(npoints > 1)
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_COMPLEX_DOUBLE,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	else
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_COMPLEX_DOUBLE,NULL,NULL,Paraddr(&mdsarg[n]));
	mdsdarg[n]->length= ArgLen(mdsdarg[n]);
	break;
      case INTEGER:
      case LOGICAL:
	if(npoints > 1)
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_LONG,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	else
	  mdsdarg[n] = MakeDescrip(calloc(1,sizeof(struct descrip)),
				   IDTYPE_LONG,NULL,NULL,Paraddr(&mdsarg[n]));
	/*printf("%d\n",*(int *)mdsdarg[n]->ptr);*/
	mdsdarg[n]->length= ArgLen(mdsdarg[n]);
	break;
      }
    } else if(mdsarg[n].type_code < 0){
	if(npoints > 1){
      mdsdarg[n] = MakeDescripWithLength(calloc(1,sizeof(struct descrip)),
	 IDTYPE_CSTRING,-mdsarg[n].type_code,mdsarg[n].ndim,mdsarg[n].extent,Paraddr(&mdsarg[n]));
	} else
      mdsdarg[n] = MakeDescripWithLength(calloc(1,sizeof(struct descrip)),
	 IDTYPE_CSTRING,-mdsarg[n].type_code,NULL,NULL,Paraddr(&mdsarg[n]));
    }
  }
  memset(&ans,0L,sizeof(struct descrip));
  mdsdarg[nargs-argind+1] = &ans;
  mdsdarg[nargs-argind+2] = NULL;


     
  mdsstatus = MdsPut(sock,node,com,mdsdarg[0],mdsdarg[1],
                     mdsdarg[2],mdsdarg[3],mdsdarg[4],mdsdarg[5],
                     mdsdarg[6],mdsdarg[7],mdsdarg[8],mdsdarg[9],mdsdarg[10],
		     mdsdarg[11],mdsdarg[12],mdsdarg[13],mdsdarg[14],mdsdarg[15],mdsdarg[16],mdsdarg[17],mdsdarg[18],mdsdarg[19],mdsdarg[20],
		     mdsdarg[21],mdsdarg[22],mdsdarg[23],mdsdarg[24],mdsdarg[25],mdsdarg[26],mdsdarg[27],mdsdarg[28],mdsdarg[29],mdsdarg[30],
		     mdsdarg[31],mdsdarg[32],mdsdarg[33],mdsdarg[34],mdsdarg[35],mdsdarg[36],mdsdarg[37],mdsdarg[38],mdsdarg[39],mdsdarg[40],
		     NULL);
  Parget(sxi);
  ix = (int *) Paraddr(sxi);
  *ix = mdsstatus;

  if((mdsstatus & 1) == 0){
    sprintf(errmsg,"%s: MdsPut error = %s.",com,ans.ptr);
    BasError(errmsg);
  }
  status = ERR;
  for(i = argind,n=0; i <= nargs; ++i,++n){
    free(mdsdarg[n]);
    Parrel(&mdsarg[n]);
  }

  setmdserror(&status);
  return(status);
}
/* 
Log into the mdsplus server as a particular user
*/
Integer Mdslogin(Integer nargs, char *fname, BA_stack_elem *sxi)
{
  BA_stack_elem syi_,*syi=&syi_;
  BA_stack_elem mdsarg[256];
  struct descrip *mdsdarg[256];
  struct descrip ans;
  char errmsg[1024];
  Integer sock;
  Integer argind;
  char user[1025],pass[2049];
  int i,j,n;
  int status,mdsstatus,*ix;

     
  setstatusreturnp(sxi,&ix);

  status = OK;

  /* get first argument which may be the socket returned by mdsconnect  */
  Fcnargb(1,syi);
  if (syi->type_code != INTEGER){
    getmdssocket(&sock); 
    printf("first argument not socket, using socket %d\n",sock);
    argind = 2;
  } else {
    sock = *(Integer *)Paraddr(syi);
    Parrel(syi);
    Fcnargb(2,syi);
    argind = 3;
  }

  if (sock == 0){
    sprintf(errmsg,"%s: sock = 0, use mdsconnect to connet to server.",fname);
    BasError(errmsg);
    status =  ERR;
    *ix = status;
    setmdserror(&status);
    return(status);
  }

  /* get second argument which is the user name */
  if (Parlen(syi) != 1){
    sprintf(errmsg,"%s: username expected.",fname);
    BasError(errmsg);
    status =  ERR;
    *ix = status;
    setmdserror(&status);
    return(status);
  }
  memset(user,0L,sizeof(user));
  strncpy(user,(char *)Paraddr(syi),((sizeof(user)-1) < -syi->type_code)?(sizeof(user)-1):-syi->type_code);
  Parrel(syi);

  /* users passwd for the mdsplus server   */
  Fcnargb(argind++,syi);
  if (Parlen(syi) != 1){
    sprintf(errmsg,"%s: password expected.",fname);
    BasError(errmsg);
    status =  ERR;
    setmdserror(&status);
    *ix = status;
    return(status);
  }
  memset(pass,0L,sizeof(pass));
  strncpy(pass,(char *)Paraddr(syi),((sizeof(pass)-1) < -syi->type_code)?(sizeof(pass)-1):-syi->type_code);
  Parrel(syi);

  mdsdarg[nargs-argind+1] = NULL;


     
  mdsstatus = MdsLogin(sock,user,pass);
  *ix = mdsstatus;
  status = mdsstatus;


  setmdserror(&status);
  return(status);
}


/* 
block and wait for mdsplus event
*/

static int id = 0;
static Integer caltrans_mdsplus_eventast(dumm)
     int *dumm;
{
  printf("mdsevent handler called with arg %d\n",*dumm);
  id = *dumm;
}
     
Integer Mdswfevent(Integer nargs, char *fname, BA_stack_elem *sxi)
{
  BA_stack_elem syi_,*syi=&syi_;
  char errmsg[1024];
  Integer sock;
  Integer shot;
  Integer argind;
  char event[1025];
  int i,j,n;
  int *status,mdsstatus;
  Integer caltrans_mdsplus_eventast();
  int *received,*ix;
  void MdsDispatchEvent();


  setstatusreturnp(sxi,&status);   


  *status = OK;
  argind = 1;
     
  /* get first argument which may be the socket returned by mdsconnect  */
  if(argind <= nargs){
    Fcnargb(argind,syi);
    if (syi->type_code != INTEGER){
      getmdssocket(&sock); 
    } else {
      sock = *(Integer *)Paraddr(syi);
      Parrel(syi);
      ++argind;
    }

    if (sock == 0){
      sprintf(errmsg,"%s: sock = 0, use mdsconnect to connet to server.",fname);
      BasError(errmsg);
      *status =  ERR;
      setmdserror(status);
      return(*status);
    }

    if (argind <= nargs) {
      Fcnargb(argind++,syi);
      if( Parlen(syi) == 1){
        memset(event,0L,sizeof(event));
        strncpy(event,(char *)Paraddr(syi),((sizeof(event)-1) < -syi->type_code)?(sizeof(event)-1):-syi->type_code);
	/*printf("sock %d:    event %s\n",sock,event);*/
	Parrel(syi);
   
	if(argind == nargs){
	  Fcnargb(argind,syi);
	  if (syi->type_code == INTEGER){
	    received = (int *)malloc(sizeof(int));
	    *received = *(int *)Paraddr(syi);
	    id = *received;
	    Parrel(syi);
	    MdsEventCan(sock,id);
	    MdsEventAst(sock,event,caltrans_mdsplus_eventast,received,&id);
	    printf("ast addy = 0x%x\n",caltrans_mdsplus_eventast);
	    printf("astprm addy = 0x%x\n",received);
	  } 
	} else {
          sprintf(errmsg,"%s: event and id expected.",fname);
          BasError(errmsg);
          *status =  ERR;
          setmdserror(status);
          return(*status);
	}
      }
    }
  } else {
    getmdssocket(&sock); 
  }
  if (sock == 0){
    sprintf(errmsg,"%s: sock = 0, use mdsconnect to connet to server.",fname);
    BasError(errmsg);
    *status =  ERR;
    setmdserror(status);
    return(*status);
  }

  /*printf("id %d received %d\n",id,*received);*/
  MdsDispatchEvent(sock);
  *status = id;

  setmdserror(status);
  return(*status);
}
Integer Mdssetupevent(Integer nargs, char *fname, BA_stack_elem *sxi)
{
  BA_stack_elem syi_,*syi=&syi_;
  char errmsg[1024];
  Integer sock;
  Integer shot;
  Integer argind;
  char event[1025];
  int i,j,n;
  int *status,mdsstatus;
  Integer caltrans_mdsplus_eventast();
  int *received,*ix;


  setstatusreturnp(sxi,&status);   

  *status = OK;
  argind = 1;
     
  /* get first argument which may be the socket returned by mdsconnect  */
  Fcnargb(argind,syi);
  if (syi->type_code != INTEGER){
    getmdssocket(&sock); 
  } else {
    sock = *(Integer *)Paraddr(syi);
    Parrel(syi);
    ++argind;
  }

  if (sock == 0){
    sprintf(errmsg,"%s: sock = 0, use mdsconnect to connet to server.",fname);
    BasError(errmsg);
    *status =  ERR;
    setmdserror(status);
    return(*status);
  }

  Fcnargb(argind++,syi);
  if( Parlen(syi) == 1){
    memset(event,0L,sizeof(event));
    strncpy(event,(char *)Paraddr(syi),((sizeof(event)-1) < -syi->type_code)?(sizeof(event)-1):-syi->type_code);
    /*printf("sock %d:    event %s\n",sock,event);*/
    Parrel(syi);
   
    if(argind == nargs){
      Fcnargb(argind,syi);
      if (syi->type_code == INTEGER){
        received = (int *)malloc(sizeof(int));
	*received = *(int *)Paraddr(syi);
	id = *received;
	Parrel(syi);
	MdsEventCan(sock,id);
	MdsEventAst(sock,event,caltrans_mdsplus_eventast,received,&id);
	printf("ast addy = 0x%x\n",caltrans_mdsplus_eventast);
      }
    } else {
      sprintf(errmsg,"%s: event and id expected.",fname);
      BasError(errmsg);
      *status =  ERR;
      setmdserror(status);
      return(*status);
    }
  } else {
    sprintf(errmsg,"%s: event and id expected.",fname);
    BasError(errmsg);
    *status =  ERR;
    setmdserror(status);
    return(*status);
  }

  setmdserror(status);
  return(*status);
}

void Mdsplusserv(Integer *pndb, Integer *pjvar, char *name,
		 Integer *ptc, Address *fwa, Integer *ndim,
		 Integer *ilow, Integer *ihi, Integer *icol,
		 char *attr, Integer param, Integer *moreargs,
		 Integer *pistage)

{
  Integer jvar, njvar, ndb, tc, istage;
  char lname[MAXSTRING];
  Integer _mdstcl();
  char comm[MAXSTRING];
  struct descrip ans;
  short lans;
  struct descrip *last;
  int sock;
  void WriteCodeInfo();

  istage = *pistage;
  ndb = *pndb;
  tc = *ptc;
  last = NULL;
  memset(&ans,0L,sizeof(struct descrip));

  getmdssocket(&sock);


  if(istage == 0){
    WriteCodeInfo();

  } else if (istage == 1){
    if(ndb == 0){
      sprintf(comm,"MDSSave: macros not implemented\n");
      BasError(comm);
    }  else {
         WriteByJvar(pndb,pjvar,name,ptc,fwa,ndim,ilow,ihi,icol);
    }
  } else if (istage == 2){
    _mdstcl("write");
  }

}
Integer WriteByJvar(Integer *pndb, Integer *pjvar, char *name,
		 Integer *ptc, Address *fwa, Integer *pndim,
		 Integer *ilow, Integer *ihi, Integer *icol)

{
  char comm[MAXSTRING];
  char dimst[MAXSTRING];
  char errmsg[1024];
  char lname[MAXSTRING];
  Integer jvar, ndb, tc, istage,ndim;
  char comment[1025];
  Integer _mdstcl();
  int len,icom;
  Integer igrp,gsize;
  Integer Rtgnvar1();
  char units[MAXSTRING];
  Integer limstr;
  Integer lilow[MAXDIMS],lihi[MAXDIMS],licol[MAXDIMS];
  Integer ailow[MAXDIMS],aihi[MAXDIMS],aicol[MAXDIMS];
  integer i;
  int limitscol[MAXDIMS];
  Integer size;
  char *ftext;
  Integer status;
  char nodename[13];
  Address av;
#if defined(BASIS_ENTRIES_2)
  BA_str *ba_com;
#endif


  jvar = *pjvar;
  ndb = *pndb;
  tc = *ptc;
  ndim = *pndim;

  if(tc != FUNCTION_TYPE && Isalloc(*fwa) == 0){
      printf("unallocated %s\n",name);
      return ERR;
  }

  limstr = Rtglimst(jvar,ndb);
  if( limstr != 0) {
    Rtxdb(jvar,ndb,fwa,&ndim,ailow,aihi,aicol,-1);
    memcpy(lilow,ilow,sizeof(lilow));
    memcpy(licol,icol,sizeof(licol));
    memcpy(lihi,ihi,sizeof(lihi));
    lihi[ndim-1] = aihi[ndim-1];
    licol[ndim-1] = lihi[ndim-1] - lilow[ndim-1] + 1;
  } else {
    memcpy(lihi,ihi,sizeof(lihi));
    memcpy(lilow,ilow,sizeof(lilow));
    memcpy(licol,icol,sizeof(licol));
    memcpy(aicol,icol,sizeof(licol));
  }
  for(i=0;i<ndim;++i){
    if(aicol[i] <= 0 || licol[i] <= 0)return;
  }
  for(i=ndim;i<MAXDIMS;++i){
    lilow[i]  = 0;
    lihi[i]  = 0;
    licol[i]  = 0;
    ailow[i]  = 0;
    aihi[i]  = 0;
    aicol[i]  = 0;
  }
   
  memset(comm,0L,sizeof(comm));
  sprintf(units,"                   ");
  Rtgunits(jvar,ndb,units);
  nullterm(units,sizeof(units)-1);
  if(strcmp(units," ") == 0) strcpy(units,"none");

  nullterm(name,MAXSTRING);

  memset(nodename,0L,sizeof(nodename));
  if(strlen(name) > 12) {
     sprintf(nodename,"bas%09d",jvar);
     sprintf(errmsg,"Variable name %s  > 12 characters, writing into node name %s", name,nodename);
     Remark(errmsg, strlen(errmsg));
  } else {
     strncpy(nodename,name,12);
  }
  sprintf(lname,":%s",nodename);
  if(tc >= 0){
    switch(tc){
    case REAL:
      sprintf(comm,"add node %s /usage=signal",lname);
      _mdstcl(comm);
      if((status = _buildsignal(lname,*fwa,ndim,licol,units,IDTYPE_FLOAT) & 1) == 0)
          return status;
      break;
    case DOUBLE:
      sprintf(comm,"add node %s /usage=signal",lname);
      _mdstcl(comm);
      if((status = _buildsignal(lname,*fwa,ndim,licol,units,IDTYPE_DOUBLE) & 1) == 0)
          return status;
      break;
    case COMPLEX:
      sprintf(comm,"add node %s /usage=signal",lname);
      _mdstcl(comm);
      if((status = _buildsignal(lname,*fwa,ndim,licol,units,IDTYPE_COMPLEX) & 1) == 0)
          return status;
      break;
    case DOUBLE_C:
      sprintf(comm,"add node %s /usage=signal",lname);
      _mdstcl(comm);
      if((status = _buildsignal(lname,*fwa,ndim,licol,units,IDTYPE_COMPLEX_DOUBLE) & 1) == 0)
          return status;
      break;
    case INTEGER:
      sprintf(comm,"add node %s /usage=signal",lname);
      _mdstcl(comm);
      if((status = _buildsignal(lname,*fwa,ndim,licol,units,IDTYPE_LONG) & 1) == 0)
          return status;
      break;
    case LOGICAL:
      sprintf(comm,"add node %s /usage=signal",lname);
      _mdstcl(comm);
      if((status = _buildsignal(lname,*fwa,ndim,licol,units,IDTYPE_LONG) & 1) == 0)
          return status;
      break;
    case FUNCTION_TYPE:
      sprintf(comm,"add node %s /usage=text",lname);
      _mdstcl(comm);
      ftext = (char *)Rtgftext(jvar,ndb,&size);
      _buildtext(lname,ftext,size);
      break;
    default:
      sprintf(comm,"MDSSave: type %d for %s not yet implemented.",tc,name);
      BasError(comm);
      return;
      break;
    }
  } else {
     sprintf(comm,"add node %s /usage=signal",lname);
     _mdstcl(comm);
    _buildtextarray(lname,*fwa,-tc,ndim,licol,units);
  }



  sprintf(lname,":%s:comment",nodename);
  sprintf(comm,"add node %s /usage=text",lname);
  _mdstcl(comm);
  len = 0;
  icom = 0;
  memset(comm,32L,sizeof(comm));
  comm[MAXSTRING-1] = '\0';
  memset(comment,0L,sizeof(comment));
#if defined(BASIS_ENTRIES_1)
  while(Rtxcomd(jvar,ndb,++icom,comm) >= 0){
#elif defined(BASIS_ENTRIES_2)
  while((ba_com = _BA_rtxcom(jvar,ndb,++icom)) != NULL){
      strncpy(comm,ba_com,sizeof(comm)-1);
#else
#error "Can't get comments for variable"
#endif
    if(icom == 2) sprintf(&comment[strlen(comment)],"\n");
    if((strlen(comment) + strlen(comm)) > 1024){
      break;
    } else {
      strcpy(&comment[strlen(comment)],comm);
    }
    memset(comm,32L,sizeof(comm));
    comm[MAXSTRING-1] = '\0';
  }
  _buildtext(lname,comment,strlen(comment));

  memset(comm,0L,sizeof(comm));
  sprintf(lname,":%s:group",nodename);
  sprintf(comm,"add node %s /usage=text",lname);
  _mdstcl(comm);
  for(igrp = 1;igrp <= 1000; ++igrp) if (Rtgnvar1(igrp, ndb) > jvar) break;
#if defined(BASIS_ENTRIES_1)
  memset(comm,32L,sizeof(comm));
  comm[MAXSTRING-1] = '\0';
  Rtggname(igrp-1, ndb, comm);
#elif defined(BASIS_ENTRIES_2)
  strncpy(comm,BA_rtggname(igrp-1,ndb),sizeof(comm)-1);
#else
#error "Can't get group name for variable"
#endif
  _buildtext(lname,comm,strlen(comm));

  memset(comm,0L,sizeof(comm));
  sprintf(lname,":%s:package",nodename);
  sprintf(comm,"add node %s /usage=text",lname);
  _mdstcl(comm);
#if defined(BASIS_ENTRIES_1)
  Glbpknam(ndb, comm,sizeof(comm)-1);
#elif defined(BASIS_ENTRIES_2)
  strncpy(comm,BA_glbpknam(ndb),sizeof(comm)-1);
#else
#error "Can't get package name for variable"
#endif
  nullterm(comm,sizeof(comm)-1);
  _buildtext(lname,comm,strlen(comm));

  sprintf(lname,":%s:units",nodename);
  sprintf(comm,"add node %s /usage=text",lname);
  _mdstcl(comm);
  _buildtext(lname,units,strlen(units));

  sprintf(lname,":%s:basisname",nodename);
  sprintf(comm,"add node %s /usage=text",lname);
  _mdstcl(comm);
  _buildtext(lname,name,strlen(name));

  if(limstr != 0){
    memset(comm,0L,sizeof(comm));
    sprintf(lname,":%s:limited",nodename);
    sprintf(comm,"add node %s /usage=text",lname);
    _mdstcl(comm);
    memset(comm,0L,sizeof(comm));
    strncpy(comm,(char *)StrgAddr(limstr),sizeof(comm)-1);
    nullterm(comm,sizeof(comm)-1);
    _buildtext(lname,comm,strlen(comm));
  }

  av = BA_db_var(ndb,jvar);
  BA_get_dim_expr(av,dimst,sizeof(dimst));
  nullterm(dimst,MAXSTRING);
  if(strncmp("scalar",dimst,6) != 0 && ndim > 0){
     memset(comm,0L,sizeof(comm));
     memset(dimst,0L,sizeof(dimst));
     sprintf(dimst,"(");
     for(i=0;i<ndim;++i){
        sprintf(&dimst[strlen(dimst)],"%d:%d",lilow[i],lihi[i]);
        if(i < (ndim-1)) sprintf(&dimst[strlen(dimst)],",");
     }
     sprintf(&dimst[strlen(dimst)],")");
     sprintf(lname,":%s:dim_string",nodename);
     sprintf(comm,"add node %s /usage=text",lname);
     _mdstcl(comm);
     _buildtext(lname,dimst,strlen(dimst));
  }

  memset(limitscol,0L,sizeof(limitscol));
  limitscol[0] = 1;
  sprintf(lname,":%s:basistype",nodename);
  sprintf(comm,"add node %s /usage=signal",lname);
  _mdstcl(comm);
  if((status = _buildsignal(lname,ptc,1,limitscol,NULL,IDTYPE_LONG) & 1) == 0)
          return status;

  memset(limitscol,0L,sizeof(limitscol));
  limitscol[0] = ndim;

  sprintf(lname,":%s:ilow",nodename);
  sprintf(comm,"add node %s /usage=signal",lname);
  _mdstcl(comm);
  if((status = _buildsignal(lname,ilow,1,limitscol,NULL,IDTYPE_LONG) & 1) == 0)
          return status;
  sprintf(lname,":%s:ihi",nodename);
  sprintf(comm,"add node %s /usage=signal",lname);
  _mdstcl(comm);
  if((status = _buildsignal(lname,ihi,1,limitscol,NULL,IDTYPE_LONG) & 1) == 0)
          return status;
  sprintf(lname,":%s:icol",nodename);
  sprintf(comm,"add node %s /usage=signal",lname);
  _mdstcl(comm);
  if((status = _buildsignal(lname,icol,1,limitscol,NULL,IDTYPE_LONG) & 1) == 0)
          return status;

  if(limstr != 0){
    sprintf(lname,":%s:limited_ilow",nodename);
    sprintf(comm,"add node %s /usage=signal",lname);
    _mdstcl(comm);
    if((status = _buildsignal(lname,ailow,1,limitscol,NULL,IDTYPE_LONG) & 1) == 0)
          return status;
    sprintf(lname,":%s:limited_ihi",nodename);
    sprintf(comm,"add node %s /usage=signal",lname);
    _mdstcl(comm);
    if((status = _buildsignal(lname,aihi,1,limitscol,NULL,IDTYPE_LONG) & 1) == 0)
          return status;
    sprintf(lname,":%s:limited_icol",nodename);
    sprintf(comm,"add node %s /usage=signal",lname);
    _mdstcl(comm);
    if((status = _buildsignal(lname,aicol,1,limitscol,NULL,IDTYPE_LONG) & 1) == 0)
          return status;
  }

}
void WriteCodeInfo()
{
  time_t t;
  char comm[MAXSTRING];
  float x;
  _mdstcl("add node :code_info /usage=text");

  _mdstcl("add node :code_info:mdsversion /usage=numeric");
  x = MDSPLUS_FILE_VERSION;
  _buildr4numeric(":code_info:mdsversion",&x);

  _mdstcl("add node :code_info:date_saved /usage=text");
  t = time(NULL);
  memset(comm,0L,sizeof(comm));
  strncpy(comm,ctime(&t),sizeof(comm));
  nlterm(comm,sizeof(comm));
  sprintf(&comm[strlen(comm)]," %s",tzname[0]);
  _buildtext(":code_info:date_saved",comm,strlen(comm));

  _mdstcl("add node :code_info:codename /usage=text");
  memset(comm,0L,sizeof(comm));
  strncpy(comm,codename,sizeof(codename));
  _buildtext(":code_info:codename",comm,strlen(comm));

  _mdstcl("add node :code_info:codedate /usage=text");
  memset(comm,0L,sizeof(comm));
  strncpy(comm,codedate,sizeof(codedate));
  _buildtext(":code_info:codedate",comm,strlen(comm));

  _mdstcl("add node :code_info:codetime /usage=text");
  memset(comm,0L,sizeof(comm));
  strncpy(comm,codetime,sizeof(codetime));
  _buildtext(":code_info:codetime",comm,strlen(comm));

  _mdstcl("add node :code_info:runtime /usage=text");
  memset(comm,0L,sizeof(comm));
  strncpy(comm,cruntime,sizeof(cruntime));
  _buildtext(":code_info:runtime",comm,strlen(comm));

  _mdstcl("add node :code_info:rundate /usage=text");
  memset(comm,0L,sizeof(comm));
  strncpy(comm,crundate,sizeof(crundate));
  _buildtext(":code_info:rundate",comm,strlen(comm));

  _mdstcl("add node :code_info:runmach /usage=text");
  memset(comm,0L,sizeof(comm));
  strncpy(comm,crunmach,sizeof(crunmach));
  _buildtext(":code_info:runmach",comm,strlen(comm));

  _mdstcl("add node :code_info:probname /usage=text");
  memset(comm,0L,sizeof(comm));
  strncpy(comm,probname,sizeof(probname));
  _buildtext(":code_info:probname",comm,strlen(comm));
}


static Integer _mdstcl(char *exp)
{
  struct descrip dexp;
  char errmsg[1024];
  char com[1024];
  struct descrip ans;
  struct descrip *last;
  integer status;
  int getmdssocket();
  int sock;
  getmdssocket(&sock);
  last = NULL;
  memset(&ans,0L,sizeof(struct descrip));
  sprintf(com,"tcl(\"%s\")",exp);
  /*printf("tcl command: %s\n",com);*/
  status = MdsValue(sock,com,&ans,last,NULL);
  /*printf("status = %d\n",status);*/
  if((status & 1) == 0){
    sprintf(errmsg,"%s: MdsTcl error = %s.",com,ans.ptr);
    BasError(errmsg);
  }
}

static Integer _buildsignal(char *node, void *data, int ndim, int *icol,char *units,char tc)
{
  struct descrip *pexp;
  char errmsg[1024];
  char com[1024];
  struct descrip ans;
  struct descrip *last;
  integer status;
  integer npoints;
  int getmdssocket();
  char lndim;
  int sock;
  int i;
  getmdssocket(&sock);
  lndim = ndim;
  last = NULL;
  memset(&ans,0L,sizeof(struct descrip));
  if(units == NULL || strncmp(units,"none",4) == 0 || strlen(units) == 0){
    sprintf(com,"build_signal($,,)");
  } else {
    sprintf(com,"build_signal(build_with_units($,\"%s\"),,)",units);
  }
  for(npoints = 1,i = 0; i < ndim; ++i){
    npoints *= icol[i];
  }
  if(npoints > 1)
    pexp = MakeDescrip(calloc(1,sizeof(struct descrip)),
		       tc,lndim,icol,data);
  else
    pexp = MakeDescrip(calloc(1,sizeof(struct descrip)),
		       tc,NULL,NULL,data);
  pexp->length= ArgLen(pexp);
  /*printf("%s %s %d %d\n",node,com,npoints,*data);*/
  status = MdsPut(sock,node,com,pexp,&ans,last,NULL);
  /*printf("status = %d\n",status);*/
  if((status & 1) == 0){
    sprintf(errmsg,"%s: BuildIData error = %s on node %s.",com,ans.ptr,node);
    BasError(errmsg);
    return status;
  }
  return OK;
  free(pexp);
}
static Integer _buildtextarray(char *node, void *data, int len, int ndim, int *icol,char *units)
{
  struct descrip *pexp;
  char errmsg[1024];
  char com[1024];
  struct descrip ans;
  struct descrip *last;
  integer status;
  integer npoints;
  int getmdssocket();
  char lndim;
  int sock;
  int i;
  getmdssocket(&sock);
  lndim = ndim;
  last = NULL;
  memset(&ans,0L,sizeof(struct descrip));
  if(units == NULL || strncmp(units,"none",4) == 0 || strlen(units) == 0){
    sprintf(com,"build_signal($,,)");
  } else {
    sprintf(com,"build_signal(build_with_units($,\"%s\"),,)",units);
  }
  for(npoints = 1,i = 0; i < ndim; ++i){
    npoints *= icol[i];
  }
  if(npoints > 1)
    pexp = MakeDescripWithLength(calloc(1,sizeof(struct descrip)),
		       IDTYPE_CSTRING,len,lndim,icol,data);
  else
    pexp = MakeDescripWithLength(calloc(1,sizeof(struct descrip)),
		       IDTYPE_CSTRING,len,NULL,NULL,data);
  pexp->length= ArgLen(pexp);
  /*printf("%s %s %d %d\n",node,com,npoints,*data);*/
  status = MdsPut(sock,node,com,pexp,&ans,last,NULL);
  /*printf("status = %d\n",status);*/
  if((status & 1) == 0){
    sprintf(errmsg,"%s: BuildIData error = %s.",com,ans.ptr);
    BasError(errmsg);
    return status;
  }
  return OK;
  free(pexp);
}
static Integer _buildr4numeric(char *node, float *rdata)
{
  struct descrip *pexp,dexp;
  char errmsg[1024];
  char com[1024];
  struct descrip ans;
  struct descrip *last;
  integer status;
  int getmdssocket();
  char ndim;
  int dim[MAX_DIMS];
  int sock;
  int i;
  float x;
  getmdssocket(&sock);
  last = NULL;
  memset(&ans,0L,sizeof(struct descrip));
  sprintf(com,"data($)");
  ndim = 1;
  dim[0] = 1;
  pexp = MakeDescrip(calloc(1,sizeof(struct descrip)),IDTYPE_FLOAT,NULL,NULL,rdata);
  pexp->length= ArgLen(pexp);
  /*printf("just before mdsput\n");*/
  /*printf("%s %s %g 0x%x\n",node,com,*rdata,pexp);*/
  status = MdsPut(sock,node,com,pexp,&ans,last,NULL);
  /*printf("status = %d\n",status);*/
  if((status % 2) == 0){
    sprintf(errmsg,"%s: Buildr4numeric error = %s.",com,ans.ptr);
    BasError(errmsg);
  }
  free(pexp);
}

static Integer _buildtext(char *node, char *text,int len)
{
  struct descrip dexp,*pexp;
  char errmsg[1024];
  char com[1024];
  char *ltext;
  struct descrip ans;
  struct descrip *last;
  integer status;
  int getmdssocket();
  int sock;
  getmdssocket(&sock);
  last = NULL;
  memset(&ans,0L,sizeof(struct descrip));
  sprintf(com,"data($)");
  ltext = calloc(1,len+1);
  memcpy(ltext,text,len);
  pexp = MakeDescripWithLength(calloc(1,sizeof(struct descrip)),

					 IDTYPE_CSTRING,len,NULL,NULL,ltext);
  status = MdsPut(sock,node,com,pexp,&ans,last,NULL);
  /*printf("status = %d\n",status);*/
  if((status % 2) == 0){
    sprintf(errmsg,"%s: BuildText error = %s status = %d. text:%s ",com,ans.ptr,status,ltext);
    BasError(errmsg);
    return status;
  }
  return OK;
  free(pexp);
  free(ltext);
}
void Mdssave(char *attr)
{
  Rtserv(attr,Mdsplusserv,NULL,"v serve:info"," ");
  Rtserv(attr,Mdsplusserv,NULL,"f db:global serve:info"," ");
}
Integer Mdssavevar(char *varname,char *pkg)
{
  Integer jvar, ndb, tc, ndim,istage,npack;
  Address fwa;
  Integer ilow[MAXDIMS],ihi[MAXDIMS],icol[MAXDIMS];
  char errmsg[1024];
  Integer status;
  char lname[MAXSTRING];


  status = OK;
  if(pkg && strlen(pkg) > 0){
     npack = Glbpknum(pkg);
#if defined(BASIS_ENTRIES_1)
     if(Parfind(npack,varname,&jvar,&ndb,&tc,-1) == ERR){
  	return ERR;
     }
#elif defined(BASIS_ENTRIES_2)
     if(BA_parfind(npack,varname,&jvar,&ndb,&tc,-1) == ERR){
  	return ERR;
     }
#else
#error "Can't call parfind for variable"
#endif
  }else{
     if(Rtfinder(varname,strlen(varname),&jvar,&ndb,&tc,"Mdssavevar") == ERR){
  	return ERR;
	}
  }

  Rtxdb(jvar,ndb,&fwa,&ndim,ilow,ihi,icol,0);
  if(ndb == 0){
    sprintf(errmsg,"MDSSave: macros not implemented\n");
    BasError(errmsg);
    status = ERR;
  } else {
       status = WriteByJvar(&ndb,&jvar,varname,&tc,&fwa,&ndim,ilow,ihi,icol);
  }

  return status;
}
static Integer nullterm(char *s,Integer len)
{
  char *c;
  c = (char *)index(s,(int)32);
  if(c)*c = '\0';
}
static Integer nlterm(char *s,Integer len)
{
  char *c;
  c = (char *)index(s,(int)'\n');
  if(c)*c = '\0';
}


/* This was written to copy a variable read into the mds package to the
package it originally came from. Borrowed alot from the pfb restore 
function. (see pfbrst.c in the pfb package contained in the basis source) */

Integer Mdscopyvar(char *svarname,char *spkg,char *dvarname, char *dpkg, Integer *ailow, Integer *aihi, Integer *aicol)
{
  Integer sjvar, sndb, stc, sndim,sistage,snpack;
  Integer djvar, dndb, dtc, dndim,distage,dnpack;
  Address sfwa;
  Address dfwa;
  Integer silow[MAXDIMS],sihi[MAXDIMS],sicol[MAXDIMS];
  Integer dilow[MAXDIMS],dihi[MAXDIMS],dicol[MAXDIMS];
  char errmsg[1024];
  Integer status;
  char lname[MAXSTRING];
  Integer kslot,isnew,doassign;
  Integer npoints;
  integer i;
  Address av;
  BA_stack_elem sx;
  Integer ns,ni[MAXDIMS], mi[MAXDIMS], ii[MAXDIMS];


  isnew = 0;
  doassign = 0;

  fprintf(stderr,"copy %s.%s to %s.%s\n",spkg,svarname,dpkg,dvarname);
  status = OK;
  if(spkg && strlen(spkg) > 0){
     snpack = Glbpknum(spkg);
#if defined(BASIS_ENTRIES_1)
     if(Parfind(snpack,svarname,&sjvar,&sndb,&stc,-1) == ERR){
        sprintf(errmsg,"Mdscopyvar:source variable %s.%s not found\n",spkg,svarname);
        BasError(errmsg);
  	return ERR;
     }
#elif defined(BASIS_ENTRIES_2)
     if(BA_parfind(snpack,svarname,&sjvar,&sndb,&stc,-1) == ERR){
        sprintf(errmsg,"Mdscopyvar:source variable %s.%s not found\n",spkg,svarname);
        BasError(errmsg);
  	return ERR;
     }
#else
#error "Can't call parfind for variable"
#endif
  }else{
     if(Rtfinder(svarname,strlen(svarname),&sjvar,&sndb,&stc," ") == ERR){
        sprintf(errmsg,"Mdscopyvar:source variable %s not found\n",svarname);
        BasError(errmsg);
  	return ERR;
	}
  }

  Rtxdb(sjvar,sndb,&sfwa,&sndim,silow,sihi,sicol,2);
  if(sndb == 0){
    sprintf(errmsg,"MDSSave: macros not implemented\n");
    BasError(errmsg);
    status = ERR;
  } else {
       if(dpkg && strlen(dpkg) > 0){
          dnpack = Glbpknum(dpkg);
#if defined(BASIS_ENTRIES_1)
          Parfind(dnpack,dvarname,&djvar,&dndb,&dtc,-1);
#elif defined(BASIS_ENTRIES_2)
          BA_parfind(dnpack,dvarname,&djvar,&dndb,&dtc,-1);
#else
#error "Can't call parfind for variable"
#endif
       }else{
          if(Rtfinder(dvarname,strlen(dvarname),&djvar,&dndb,&dtc," ") == ERR){
	      dndb = Glbpknum("global");
	  }
       }
       /*printf("djvar = %d dndb = %d\n",djvar,dndb);*/
       if(djvar == 0){
         isnew = 1;
         sprintf(errmsg,"new variable %s\n",dvarname);
         Remark(errmsg,strlen(errmsg));
	 kslot = StrgFind(dvarname,strlen(dvarname));
	 sx.type_code = Utstrcod("name");
	 sx.ptr = kslot;
	 sx.extent[0] = 0;
	 sx.naml = Utgetcl(dvarname);
	 sx.ns = 0;
	 dtc = stc;
	 djvar = Rtinname(dndb,dtc,&sx);
	 if(djvar <= 0){
	   sprintf(errmsg,"Mdscopyvar: variable create failed on %s", dvarname);
           Remark(errmsg, strlen(errmsg));
	   Kaboom(0);
	 }
	 
       } else if (djvar > 0){
         Rtxdb(djvar,dndb,&dfwa,&dndim,dilow,dihi,dicol,2);
/*         sprintf(errmsg,"variable exists %s\n",dvarname);
         Remark(errmsg,strlen(errmsg));*/
       } else {
         sprintf(errmsg,"unknown variable %s\n",dvarname);
	 Remark(errmsg,strlen(errmsg));
	 return ERR;
       }


       av = BA_db_var(dndb,djvar);

       if(BA_is_dynamic(av)){
	 if(Isalloc(dfwa) == 0){
            sprintf(errmsg,"variable is dynamic and not allocated %s\n",dvarname);
            Remark(errmsg,strlen(errmsg));
	 }
           if(Rtischam(djvar,dndb)){
	      doassign = 1;
	      isnew = 1;
	      Rtstc(djvar,dndb,stc);
	   }
           Rtxdb(sjvar,sndb,&sfwa,&dndim,dilow,dihi,dicol,2);
           Rtstgget(djvar,dndb,dndim,dicol,dilow,0);
           Rtxdb(djvar,dndb,&dfwa,&dndim,dilow,dihi,dicol,2);
           if(*ailow != (Integer *)DEFAULT) memcpy(dilow,ailow,sndim*sizeof(Integer));
           if(*aihi != (Integer *)DEFAULT) memcpy(dihi,aihi,sndim*sizeof(Integer));
           if(*aicol != (Integer *)DEFAULT) memcpy(dicol,aicol,sndim*sizeof(Integer));
       }




       sx.type_code = stc;
       Rtxdb(sjvar,sndb,&sfwa,&sx.ndim,sx.lower,ni,sx.extent,2);
       Pargetc(&sx,sfwa);
       if(dtc != FUNCTION_TYPE){
          ns = 0;
          status = Parasgn2(dtc,dfwa,dndim,dilow,dihi,dicol,ns,ni,mi,ii,&sx);
          if(status == ERR){
            sprintf(errmsg,"mdscopyvar:Could not copy  %s.%s to %s.%s\n",spkg,svarname,dpkg,dvarname);
            Remark(errmsg,strlen(errmsg));
            status = -1;
          }
       } else {
          ns = -stc;
          Parsestr_((Address *) &sfwa,&ns);
       }



    
  }

  return status;

}

