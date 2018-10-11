#include <stdio.h>
#include <score.h>
#include <Basis.h>
#include <par_.h>
#include <dlfcn.h>
#include "links.h"



#define LIBLEN 1024
#define SYMLEN 256
#define MAXARG 1024


/* Struct to create linked list of opened libraries and symbols  */
struct dlentry_ {
	struct linkage link;
	char libpath[LIBLEN];
        void *handle;
	char symname[SYMLEN];
        int (*sym)();
};

typedef struct dlentry_ dlentry;




/* Set basis return for single integer returns
 * Don't do this more than once per return to basis as
 * it will create a memory leak. */
static Integer setstatusreturnp(BA_stack_elem *sxi,int **ix)
{
	sxi->ndim = 0;
	sxi->extent[0] = 1;
	sxi->lower[0] = 1;
	sxi->stride[0] = 1;
	sxi->type_code = INTEGER;
	Parget(sxi);
	*ix = (int *) Paraddr(sxi);
}
static Integer setstatusreturn(BA_stack_elem *sxi,int status)
{
	                                                                                
	  int *ix;
          setstatusreturnp(sxi,&ix);
	  *ix = status;
}


Integer SHLbfcn(Integer nargs,char *fname, BA_stack_elem *sx)
{

  static char msg1[] = "SHL: Function not implemented yet.";
  Integer call_external();
  int status;
	char *c;

  sx->ndim = 0;
  sx->type_code = 0;
  status = OK;
  if(strncmp(fname,"call_external",13) == 0){
	status  = call_external(nargs,fname);
	setstatusreturn(sx,status);
  } else if(strncmp(fname,"exec_external",13) == 0){
	status  = exec_external(nargs,fname);
	setstatusreturn(sx,status);
  } else {
    BasError(msg1); 
    BasError(fname); 
    status = ERR;
  }
  
  return status;


}





/* Get the library path and symbol name from the basis arguments.  */
static int get_lib_sym(Integer nargs, char *fname, char *libpath, char *symname)
{
    BA_stack_elem syi_,*syi=&syi_;
    int status;
    char errmsg[1024];
    Integer argind;

    status = OK;
    argind = 1;

    Fcnargb(argind++,syi);
    if (syi->type_code >= 0){
      sprintf(errmsg,"%s: Path to shared library expected.",fname);
      BasError(errmsg);
      status =  ERR;
      return(status);
    }
    memset(libpath,0L,LIBLEN);
    strncpy(libpath,(char *)Paraddr(syi),((LIBLEN-1) < -syi->type_code)?(LIBLEN-1):-syi->type_code);
    Parrel(syi);


    Fcnargb(argind++,syi);
    if (syi->type_code >= 0){
      sprintf(errmsg,"%s: Symbol name in shared library expected.",fname);
      BasError(errmsg);
      status =  ERR;
      return(status);
    }
    memset(symname,0L,SYMLEN);
    strncpy(symname,(char *)Paraddr(syi),((SYMLEN-1) < -syi->type_code)?(SYMLEN-1):-syi->type_code);
    Parrel(syi);

    return(status);


}


/* Meant to be basically the same as the IDL call_external in 
 * the way arguments are passed (argc,argv). It's completely
 * up the user to ensure that the arguments match up. */
Integer call_external(Integer nargs, char *fname)
{
    BA_stack_elem syi_,*syi=&syi_;
    BA_stack_elem barg[MAXARG];
    char libpath[LIBLEN];
    char symname[SYMLEN];
    int status;
    char errmsg[1024],*error;
    void *ce_args[MAXARG];
    int i,argind;
    dlentry *e, *find_sym(),*add_sym();

    if( get_lib_sym(nargs,fname,libpath,symname) != OK) return(ERR);

    if( (e = find_sym(fname,libpath,symname)) != (dlentry *) NULL){
	 status = OK;
    } else if ( (e = add_sym(fname,libpath,symname)) != (dlentry *) NULL){
	 status = OK;
    } else {
         status =  ERR;
	 return(status);
    }

    argind = 3;	/* the first two args are found with get_lib_sym  */
    i = 0;
    while( i < MAXARG &&  argind <= nargs){
        Fcnargb(argind,&barg[i]);
	ce_args[i] = (void *) Paraddr(&barg[i]);
	++i;
	++argind;
    }
    status = (*(e->sym))(i,ce_args);
    argind = 3;
    i = 0;
    while( i < MAXARG &&  argind <= nargs){
	Parrel(&barg[i]);
	++i;
	++argind;
    }
    return(status);
}

/* Similar to call_external except arguments passed directly
 * instead of in an (argc,argv) style. It's completely
 * up the user to ensure that the arguments match up. */
Integer exec_external(Integer nargs, char *fname)
{
    BA_stack_elem syi_,*syi=&syi_;
    BA_stack_elem barg[MAXARG];
    char libpath[LIBLEN];
    char symname[SYMLEN];
    int status;
    char errmsg[1024],*error;
    void *ca[MAXARG];
    int i,argind;
    dlentry *e, *find_sym(),*add_sym();

    if(nargs > (31)){
         sprintf(errmsg,"%s: Too many arguments, exec_external only supports 30 agruments at the moment.",fname);
         BasError(errmsg);
	 return(ERR);
    }
    if( get_lib_sym(nargs,fname,libpath,symname) != OK) return(ERR);

    if( (e = find_sym(fname,libpath,symname)) != (dlentry *) NULL){
	 status = OK;
    } else if ( (e = add_sym(fname,libpath,symname)) != (dlentry *) NULL){
	 status = OK;
    } else {
         status =  ERR;
	 return(status);
    }

    argind = 3;	/* the first two args are found with get_lib_sym  */
    i = 0;
    while( i < MAXARG &&  argind <= nargs){
        Fcnargb(argind,&barg[i]);
	ca[i] = (void *) Paraddr(&barg[i]);
	++i;
	++argind;
    }
    status = (*(e->sym))(ca[0],ca[1],ca[2],ca[3],ca[4],ca[5],ca[6],ca[7],
		    ca[8],ca[9],ca[10],ca[11],ca[12],ca[13],ca[14],ca[15],
		    ca[16],ca[17],ca[18],ca[19],ca[20],ca[21],ca[22],ca[23],
		    ca[24],ca[25],ca[26],ca[27],ca[28],ca[29]);
    argind = 3;
    i = 0;
    while( i < MAXARG &&  argind <= nargs){
	Parrel(&barg[i]);
	++i;
	++argind;
    }
    return(status);
}

static struct linkage *dllhead = NULL;
static struct linkage *dlltail = NULL;


/* Run the linked list and return a pointer to the first entry structure
 * that matches both library path and symbol name.  */
dlentry *find_sym(fname,libpath,symname)
char *fname,*libpath,*symname;
{
	dlentry *e;
	for(e = (dlentry *)dllhead; e != NULL; e = (dlentry *)e->link.next){
		if(strncmp(libpath,e->libpath,LIBLEN) == 0 &&
		   strncmp(symname,e->symname,SYMLEN) == 0)
		   	return e;
	}

	return (dlentry *)NULL;

}



/* If not already in linked list open the library, find the symbol,
 * and put an entry into the linked list. */
dlentry *add_sym(fname,libpath,symname)
char *fname,*libpath,*symname;
{
    dlentry *e,*find_sym();
    int (*userfunc)();
    char errmsg[1024],*error;

    if((e = find_sym(fname,libpath,symname)) != (dlentry *)NULL)return e;

    e = (dlentry *)calloc(1,sizeof(dlentry));

    if(e){
       e->handle = dlopen(libpath,RTLD_LAZY|RTLD_GLOBAL);/* dlopen is a system call  */
       if(!e->handle){
         sprintf(errmsg,"%s: Couldn't open shared library %s.",fname,libpath);
         BasError(errmsg);
         sprintf(errmsg,"%s",dlerror());
         BasError(errmsg);
	 free(e);
	 return (dlentry *)NULL;
       }
       strncpy(e->libpath,libpath,LIBLEN);
       strncpy(e->symname,symname,SYMLEN);
       e->sym = dlsym(e->handle,symname);
       if ((error = dlerror()) != NULL)  {
         sprintf(errmsg,"%s: Couldn't find sym %s in shared library %s.",fname,symname,libpath);
         BasError(errmsg);
         sprintf(errmsg,"%s",error);
         BasError(errmsg);
	 free(e);
	 return (dlentry *)NULL;
       }
       add_link(e,&dllhead,&dlltail);
    }
    /*free(error);*/	/* It's unclear whether the string returned by dlerror
			 needs to be freed or not.*/

    return e;

}



/* Run the linked list and dlclose all handles of matching libraries. After
 * the last one is closed the library should be closed. A subsequent 
 * call_external should get any new library if changes were made. */
int Close_library(libpath)
char *libpath;
{
	dlentry *e;
	dlentry t;
	int found;
	found = ERR;
	for(e = (dlentry *)dllhead; e != NULL; e = (dlentry *)e->link.next){
		if(strncmp(libpath,e->libpath,LIBLEN) == 0){
			found = OK;
                        rem_link(e,&dllhead,&dlltail);
			dlclose(e->handle);
			memcpy(&t,e,sizeof(dlentry));
			free(e);
			e = &t;
		}
	}

	return found;

}
