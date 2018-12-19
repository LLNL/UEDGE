/*

=======================================================
       Import of Data from File (IDF) package
     for time-saving programming of data input
              to scientific codes.

                 IDF Version 1.1   
              Date: February 1, 2000


 Copyright 1999 by A. Yu. Pigarov and I.V. Saltanova


                   Developers:
         A. Yu. Pigarov and I.V. Saltanova


         E-mails: apigarov@psfc.mit.edu
                  apigarov@pppl.gov
                  iras@rex.pfc.mit.edu

         IDF Documentation is available at:
  http://www2.psfc.mit.edu/library/preprints.html

          A.Yu. Pigarov, I.V. Saltanova, 
      "Import Data from File (IDF) utilities
  for programming data input to scientific codes",
         Plasma Science and Fusion Center, 
       Massachusetts Institute of Technology,
    Cambridge, MA, PSFC/JA-00-01, January 2000.



    The IDF package is distributed in the hope 
that it will be useful to many computational scientists.

    The IDF package is FREE SOFTWARE. 
You can use, copy, and modify this software for any purpose
and without fee provided that the above copyright
notice appear in all copies.

    The IDF package is distributed "AS IS", 
i.e without any warranty including all implied warranties
of merchantability and fitness.
   
=======================================================

*/


#include "idflib.h"


#include "idf.h"
#include "idfusr.h"



#if IDF_CPP_ == 1
int idf_c(char* name, char* val)
#else
int idf_c(name, val)
char *name;
char *val;
#endif
{
int ierr;
unsigned int uval;
void *cval;

cval = &uval;
ierr = idf_get_one(name, 'c' , 1, cval);
*val = (char) uval;

return ierr;
}


#if IDF_CPP_ == 1
int idf_uc(char* name, unsigned char* val)
#else
int idf_uc(name, val)
char *name;
unsigned char *val;
#endif
{
int ierr;
unsigned int uval;
void *cval;

cval = &uval;
ierr = idf_get_one(name, 'c' , 1, cval);
*val = (unsigned char) uval;

return ierr;
}


#if IDF_CPP_ == 1
int idf_u(char* name, unsigned int* val)
#else
int idf_u(name, val)
char *name;
unsigned int *val;
#endif
{
int ierr;
void *cval;

cval = val;
ierr = idf_get_one(name, 'c' , 1, cval);

return ierr;
}


#if IDF_CPP_ == 1
int idf_s(char* name, char* val, int n)
#else
int idf_s(name, val, n)
char *name;
char *val;
int   n;
#endif
{
int ierr;
void *cval;

cval = val;
ierr = idf_get_one(name, 's' , n, cval);

return ierr;
}


#if IDF_CPP_ == 1
int idf_t(char* name, char* val, int n)
#else
int idf_t(name, val, n)
char *name;
char *val;
int   n;
#endif
{
int ierr;
void *cval;

cval = val;
ierr = idf_get_one(name, 't' , n, cval);

return ierr;
}


#if IDF_CPP_ == 1
int idf_i(char* name, int* val)
#else
int idf_i(name, val)
char *name;
int  *val;
#endif
{
int ierr;
void *ival;

ival = val;
ierr = idf_get_one(name, 'i' , 1, ival);

return ierr;
}


#if IDF_CPP_ == 1
int idf_l(char* name, long* val)
#else
int idf_l(name, val)
char *name;
long *val;
#endif
{
int ierr;
void *lval;

lval = val;
ierr = idf_get_one(name, 'l', 1, lval);

return ierr;
}


#if IDF_CPP_ == 1
int idf_f(char* name, float* val)
#else
int idf_f(name, val)
char  *name;
float *val;
#endif
{
int ierr;
void *fval;

fval = val;
ierr = idf_get_one(name, 'f' , 1, fval);

return ierr;
}


#if IDF_CPP_ == 1
int idf_d(char* name, double* val)
#else
int idf_d(name, val)
char   *name;
double *val;
#endif
{
int ierr;
void *dval;

dval = val;
ierr = idf_get_one(name, 'd' , 1, dval);

return ierr;
}


#if IDF_CPP_ == 1
int idf_z(char* name, float* val)
#else
int idf_z(name, val)
char  *name;
float *val;
#endif
{
int ierr;
void *fval;

fval = val;
ierr = idf_get_one(name, 'z' , 1, fval);

return ierr;
}


#if IDF_CPP_ == 1
int idf_w(char* name, double* val)
#else
int idf_w(name, val)
char   *name;
double *val;
#endif
{
int ierr;
void *dval;

dval = val;
ierr = idf_get_one(name, 'w' , 1, dval);

return ierr;
}


#if IDF_CPP_ == 1
int idf_iarr(char* name, int* val, int n)
#else
int idf_iarr(name, val, n)
char *name;
int  *val;
int   n;
#endif
{
int ierr;
void *ival;

ival = val;
ierr = idf_get_1Darray(name, 'i' , n, ival);

return ierr;
}


#if IDF_CPP_ == 1
int idf_larr(char* name, long* val, int n)
#else
int idf_larr(name, val, n)
char *name;
long *val;
int   n;
#endif
{
int ierr;
void *lval;

lval = val;
ierr = idf_get_1Darray(name, 'l' , n, lval);

return ierr;
}



#if IDF_CPP_ == 1
int idf_farr(char* name, float* val, int n)
#else
int idf_farr(name, val, n)
char  *name;
float *val;
int    n;
#endif
{
int ierr;
void *fval;

fval = val;
ierr = idf_get_1Darray(name, 'f' , n, fval);

return ierr;
}



#if IDF_CPP_ == 1
int idf_darr(char* name, double* val, int n)
#else
int idf_darr(name, val, n)
char   *name;
double *val;
int     n;
#endif
{
int ierr;
void *dval;

dval = val;
ierr = idf_get_1Darray(name, 'd' , n, dval);

return ierr;
}



#if IDF_CPP_ == 1
int idf_zarr(char* name, float* val, int n)
#else
int idf_zarr(name, val, n)
char  *name;
float *val;
int    n;
#endif
{
int ierr;
void *fval;

fval = val;
ierr = idf_get_1Darray(name, 'z' , n, fval);

return ierr;
}



#if IDF_CPP_ == 1
int idf_warr(char* name, double* val, int n)
#else
int idf_warr(name, val, n)
char   *name;
double *val;
int     n;
#endif
{
int ierr;
void *dval;

dval = val;
ierr = idf_get_1Darray(name, 'w' , n, dval);

return ierr;
}
