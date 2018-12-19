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
int IDF_FORTRAN(idfinit_,idfinit,IDFINIT)
               (int* dummy)
#else
int IDF_FORTRAN(idfinit_,idfinit,IDFINIT)
               (dummy)
int* dummy;
#endif
{
 return ( idf_init() );
}



void IDF_FORTRAN(idfinish_,idfinish,IDFINISH)  ()
{
 idf_finish();
}



#if IDF_CPP_ == 1
int IDF_FORTRAN(idfini_,idfini,IDFINI)
               (char* IDFfile_name, int* n)
#else
int IDF_FORTRAN(idfini_,idfini,IDFINI)
               (IDFfile_name,n)
char *IDFfile_name;
int  *n;
#endif
{
int ierr=0;
int jerr;
int nn;
char *sbuf;

idf_err_clear();

sbuf = NULL;
 nn = *n;


/*----------------------
  check input file name
*-----------------------*/

   jerr = idf_file_name_(IDFfile_name,nn);
if(jerr) {ierr=1; goto err;} /*err:1-5*/

   sbuf = idf_file_Cname();
 
/*-------------------
  init internal data
 --------------------*/

   jerr = idf_file_ini();
if(jerr) {ierr=2; jerr=5+jerr; goto err;} /*err:6,7*/

   jerr = idf_frmt_ini();
if(jerr) {ierr=3; jerr=7+jerr; goto err;} /*err:8,9,10*/

   jerr = idf_name_ini();
if(jerr) {ierr=4; jerr=10+jerr; goto err;}/*err:11,12*/

          idf_list_nul();
   jerr = idf_list_ini();
if(jerr) {ierr=5; jerr=13; goto err;} /*err:13*/

   jerr = idf_data_ini();
if(jerr) {ierr=6; jerr=13+jerr; goto err;}/*err:14,15*/

   jerr = idf_dbuf_ini();
if(jerr) {ierr=7; jerr=15+jerr; goto err;}/*err:16,17*/

   jerr = idf_val_ini();
if(jerr) {ierr=8; jerr=17+jerr; goto err;}/*err:18,19*/

   jerr = idf_fbuf_ini();
if(jerr) {ierr=7; jerr=25+jerr; goto err;}/*err:26,27*/

   jerr = idf_form_ini();
if(jerr) {ierr=8; jerr=27+jerr; goto err;}/*err:28,29*/

#if IDF_FLOAT_DEFAULT == 0
   jerr = idf_form_ini();
if(jerr) {ierr=7; jerr=30; goto err;}/*err:30*/
#endif


/*-------------------
     open IDFfile
 --------------------*/

   jerr = idf_file_open(sbuf);
if(jerr) {ierr=9; jerr=20; goto err;} /*err;20*/

idf_file_Cname_free();

/*------------------------------
  create catalogs of named data
 -------------------------------*/

   jerr = idf_catalog();
if(jerr) {ierr=11; jerr=21+jerr;} /*22-25*/

idf_order_put( IDF_ORDER_OPEN );

idf_jus_fort();

err:
if(ierr)
 {
  idf_err_put(ierr,jerr);

#if IDF_MISTAKE == 1
  idf_err_prn();
#endif
 }

return ierr;
}



void IDF_FORTRAN(idfend_,idfend,IDFEND) ()
{
 idf_end();
}



#if IDF_CPP_ == 1
int IDF_FORTRAN(idfopen_,idfopen,IDFOPEN)
               (char* IDFfile_name, int* n)
#else
int IDF_FORTRAN(idfopen_,idfopen,IDFOPEN)
               (IDFfile_name,n)
char *IDFfile_name;
int  *n;
#endif
{
/* to be called from FORTRAN */

int ierr=0;
char *sbuf;
int   jerr;
int   nn;

sbuf = NULL;
  nn = *n;

   jerr = idf_order_cmp( IDF_ORDER_INIT );
if(jerr)
 {
  ierr=10;
  jerr=21;
  goto err;
 }

   jerr = idf_file_name_(IDFfile_name,nn);
if(jerr) {ierr=1;goto err;} /*err:1-5*/

   sbuf = idf_file_Cname();

   jerr = idf_file_open(sbuf);
if(jerr) {ierr=9; jerr=20; goto err;}/*err:20*/

idf_file_Cname_free();


   jerr = idf_catalog();
if(jerr) {ierr=11; jerr=22;} /*22*/

idf_order_put( IDF_ORDER_OPEN );

idf_jus_fort();

err:
if(ierr)
 {
  idf_err_put(ierr,jerr);

#if IDF_MISTAKE == 1
  idf_err_prn();
#endif
 }

return ierr;
}



void IDF_FORTRAN(idfclose_,idfclose,IDFCLOSE) ()
{
 idf_close();
}



void IDF_FORTRAN(idferprn_,idferprn,IDFERPRN) ()
{
 idf_err_prn();
}



#if IDF_CPP_ == 1
int IDF_FORTRAN(idfc_,idfc,IDFC)
               (char* name, int* nam, void* val)
#else
int IDF_FORTRAN(idfc_,idfc,IDFC)
               (name, nam, val)
char *name;
int  *nam;
void *val;
#endif
{
int ierr;
int jus,nm;
unsigned int uval;
void *cval;
unsigned char uc,*puc;
char           c,*pc;

  nm = *nam;
cval = &uval;

ierr = idf_get_one_(name, nm, 'c' , 1, cval);

   jus = idf_jus_get();
if(jus)
 {
    uc = (unsigned char) uval;
   puc = (unsigned char*) val;
  *puc = uc;
 }
else
 {
    c = (char) uval;
   pc = (char*) val;
  *pc = c;  
 }

return ierr;
}



#if IDF_CPP_ == 1
int IDF_FORTRAN(idfs_,idfs,IDFS)
               (char* name, int* nam, 
                void*  val, int* n)
#else
int IDF_FORTRAN(idfs_,idfs,IDFS)
                (name, nam, val, n)
char *name;
int  *nam;
void *val;
int  *n;
#endif
{
int ierr;
int nm,nn;

   nm = *nam;
   nn = *n;

 idf_jus_fort();

 ierr = idf_get_one_(name, nm, 's' , nn, val);

return ierr;
}


#if IDF_CPP_ == 1
int IDF_FORTRAN(idft_,idft,IDFT)
               (char* name, int* nam, 
                void*  val, int* n)
#else
int IDF_FORTRAN(idft_,idft,IDFT)
               (name, nam, val, n)
char *name;
int  *nam;
void *val;
int  *n;
#endif
{
int ierr;
int nm,nn;

  nm = *nam;
  nn = *n;

idf_jus_fort();

ierr = idf_get_one_(name, nm, 't' , nn, val);

return ierr;
}


#if IDF_CPP_ == 1
int IDF_FORTRAN(idfi_,idfi,IDFI)
                (char* name, int* nam, int* val)
#else
int IDF_FORTRAN(idfi_,idfi,IDFI)
               (name, nam, val)
char *name;
int  *nam;
int  *val;
#endif
{
int ierr,nm;
void *ival;

  nm = *nam;
ival = val;
ierr = idf_get_one_(name, nm, 'i' , 1, ival);

return ierr;
}


#if IDF_CPP_ == 1
int IDF_FORTRAN(idff_,idff,IDFF)
               (char* name, int* nam, float* val)
#else
int IDF_FORTRAN(idff_,idff,IDFF)
(name, nam, val)
char  *name;
int   *nam;
float *val;
#endif
{
int ierr,nm;
void *fval;

  nm = *nam;
fval = val;
ierr = idf_get_one_(name, nm,'f' , 1, fval);

return ierr;
}


#if IDF_CPP_ == 1
int IDF_FORTRAN(idfd_,idfd,IDFD)
               (char* name, int* nam, double* val)
#else
int IDF_FORTRAN(idfd_,idfd,IDFD)
               (name, nam, val)
char   *name;
int    *nam;
double *val;
#endif
{
int ierr,nm;
void *dval;

  nm = *nam;
dval = val;
ierr = idf_get_one_(name, nm, 'd' , 1, dval);

return ierr;
}


#if IDF_CPP_ == 1
int IDF_FORTRAN(idfz_,idfz,IDFZ)
               (char* name, int* nam, float* val)
#else
int IDF_FORTRAN(idfz_,idfz,IDFZ)
               (name, nam, val)
char  *name;
int   *nam;
float *val;
#endif
{
int ierr,nm;
void *fval;

  nm = *nam;
fval = val;
ierr = idf_get_one_(name, nm, 'z' , 1, fval);

return ierr;
}


#if IDF_CPP_ == 1
int IDF_FORTRAN(idfw_,idfw,IDFW)
               (char* name, int* nam, double* val)
#else
int IDF_FORTRAN(idfw_,idfw,IDFW)
               (name, nam, val)
char   *name;
int    *nam;
double *val;
#endif
{
int ierr,nm;
void *dval;

  nm = *nam;
dval = val;
ierr = idf_get_one_(name, nm, 'w' , 1, dval);

return ierr;
}


#if IDF_CPP_ == 1
int IDF_FORTRAN(idfiarr_,idfiarr,IDFIARR)
               (char* name, int* nam,
                int* val, int* n)
#else
int IDF_FORTRAN(idfiarr_,idfiarr,IDFIARR)
               (name, nam, val, n)
char *name;
int  *nam;
int  *val;
int  *n;
#endif
{
int ierr,nm,nn;
void *ival;

  nm = *nam;
  nn = *n;
ival = val;
ierr = idf_get_1Darray_(name, nm, 'i' , nn, ival);

return ierr;
}


#if IDF_CPP_ == 1
int IDF_FORTRAN(idffarr_,idffarr,IDFFARR)
               (char* name, int* nam,
                float* val, int* n)
#else
int IDF_FORTRAN(idffarr_,idffarr,IDFFARR)
               (name, nam, val, n)
char  *name;
int   *nam;
float *val;
int   *n;
#endif
{
int ierr,nm,nn;
void *fval;

  nm = *nam;
  nn = *n;
fval = val;
ierr = idf_get_1Darray_(name, nm, 'f' , nn, fval);

return ierr;
}


#if IDF_CPP_ == 1
int IDF_FORTRAN(idfdarr_,idfdarr,IDFDARR)
               (char* name, int* nam,
                double* val, int* n)
#else
int IDF_FORTRAN(idfdarr_,idfdarr,IDFDARR)
               (name, nam, val, n)
char   *name;
int    *nam;
double *val;
int    *n;
#endif
{
int ierr,nm,nn;
void *dval;

  nm = *nam;
  nn = *n;
dval = val;
ierr = idf_get_1Darray_(name, nm, 'd' , nn, dval);

return ierr;
}


#if IDF_CPP_ == 1
int IDF_FORTRAN(idfzarr_,idfzarr,IDFZARR)
               (char* name, int* nam,
                float* val, int* n)
#else
int IDF_FORTRAN(idfzarr_,idfzarr,IDFZARR)
               (name, nam, val, n)
char  *name;
int   *nam;
float *val;
int   *n;
#endif
{
int ierr,nm,nn;
void *fval;

  nm = *nam;
  nn = *n;
fval = val;
ierr = idf_get_1Darray_(name, nm, 'z' , nn, fval);

return ierr;
}


#if IDF_CPP_ == 1
int IDF_FORTRAN(idfwarr_,idfwarr,IDFWARR)
               (char* name, int* nam,
                double* val, int* n)
#else
int IDF_FORTRAN(idfwarr_,idfwarr,IDFWARR)
               (name, nam, val, n)
char   *name;
int    *nam;
double *val;
int    *n;
#endif
{
int ierr,nm,nn;
void *dval;

  nm = *nam;
  nn = *n;
dval = val;
ierr = idf_get_1Darray_(name, nm, 'w' , nn, dval);

return ierr;
}



#if IDF_CPP_ == 1
int IDF_FORTRAN(idfget_,idfget,IDFGET)
               (char* name, int* nam,
                char* frmt, int* nfrmt, int* align,
                void* val)
#else
int IDF_FORTRAN(idfget_,idfget,IDFGET)
               (name,nam, frmt,nfrmt, align, val)
char   *name;
int    *nam;
char   *frmt;
int    *nfrmt;
int    *align;
void   *val;
#endif
{
int ierr;
int nm,nn,algn;

  nm = *nam;
  nn = *nfrmt;
algn = *align;

idf_jus_fort();

ierr = idf_get_(name,nm, frmt,nn, algn, val);

return ierr;
}



#if IDF_CPP_ == 1
int IDF_FORTRAN(idfarray_,idfarray,IDFARRAY)
               (char* name, int* nam,
                char* frmt, int* nfrmt,
                void* val)
#else
int IDF_FORTRAN(idfarray_,idfarray,IDFARRAY)
               (name,nam, frmt,nfrmt, val)
char   *name;
int    *nam;
char   *frmt;
int    *nfrmt;
void   *val;
#endif
{
static int align = 0;
int ierr;
int nm,nn;

  nm = *nam;
  nn = *nfrmt;

idf_jus_fort();

ierr = idf_get_(name,nm, frmt,nn, align, val);

return ierr;
}
