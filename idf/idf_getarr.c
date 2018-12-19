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
int idf_get_1Darray(char* name, char frmt, int Mdim,
                    void* value)
#else
int idf_get_1Darray(name, frmt, Mdim, value)
char *name;
char  frmt;
int   Mdim;
void *value;
#endif
{
int ierr;
int nam;

nam = strlen(name);

ierr = idf_get_1Darray_(name,nam, frmt, Mdim, value);

return ierr;
}



#if IDF_CPP_ == 1
int idf_get_1Darray_(char* name, int nam,
                     char frmt, int Mdim,
                     void* value)
#else
int idf_get_1Darray_(name,nam, frmt, Mdim, value)
char *name;
int   nam;
char  frmt;
int   Mdim;
void *value;
#endif
{
static int   Ndim = 1;
static int   style = 1;
int ierr;
int dim;

dim = Mdim;

ierr = idf_getarray(name,nam, frmt, Ndim,&Mdim, style, value);

return ierr;
}



#if IDF_CPP_ == 1
int idf_getarray(char* name, int nam, char frmt,
                 int Ndim, int* Mdim, int style,
                 void* value)
#else
int idf_getarray(name,nam, frmt, Ndim,Mdim, style, value)
char *name;
int   nam;
char  frmt;
int   Ndim;
int  *Mdim;
int   style;
void *value;
#endif
{
/*
-----------------------------------------------
import of single format 1D array of data units

'name'  is the string for named data
'nam'   is the number of symbols in name
'frmt'  is the format symbol
'Ndim'  number of dimensions assigned to array
'Mdim'  dimensions of array
'style' type of storage for array:
        0 = C/C++ and 1 = Fortran convention
        That is, if 0 then the last index varies faster
        and if 1 then the first index varies faster.
'val'   is the pointer to the target
-----------------------------------------------

In success it returns the number of imported elements
and returns zero at an error;
Caution: because of offset option, array can contain gaps.
*/

int ierr=0;
int nval=1;
int Nlist;
register int i;
unsigned int cfrmt;
int nmax,ofs,jelem,kup,kkup;
int jsize,jfrmt,kfrmt;
int ndim;
int *mdim;
int krec,ibrk,jerr,flag,lvl;
IDF_NAME_ *Anam;
char *val;
void *tag;

idf_frmt_clean();
idf_frmtreat_nul();

lvl   = 0;
cfrmt = frmt;
nmax  = 0;
jsize = 0;
val   = (char*) value;
jelem = 0;
kup   = 0;

   jerr = idf_order_cmp( IDF_ORDER_OPEN );
if(jerr)
 {
  /* improper order of function call*/
  ierr=28;
  goto err;
 }

if(nam<1)
 {
  /*empty name*/
  ierr=21;
  goto err;
 }

   jerr = idf_extname_check(name,nam);
if(jerr)
 {
  /* improper input name*/
  if(jerr==4)
   {
    /* improper size of name*/
    ierr=22;
   }
  else
   {
    ierr = jerr; /*1-3*/
   }
  goto err;
 }

    jfrmt = idf_frmt_type(cfrmt);
if(!jfrmt)
 {
  /* improper format symbol*/
  ierr = 4;
  goto err;
 }

jsize = idf_frmt_size(jfrmt);


   nmax = idf_dim_nmax(Ndim,Mdim);
if(nmax<0)
 {
  /* improper dimensionality */
  ierr = 5-nmax; /*1-2*/
  goto err;
 }

   Nlist = idf_list_num();
if(Nlist < 1)
 {
  /* empty name list */
  ierr = 5;
  goto err;
 }

krec=0;
kup =0;
for(i=0;i<Nlist;i++)
 {
     Anam = idf_list_get(name, nam, i); 
  if(Anam != NULL)
   {
    krec++;
    kkup=0;

       jerr = idf_data_begin(Anam);
    if(jerr)
     {
      /* unable to move to data field*/
      lvl = 1;
      ierr = 11;
      break;
     }

      ndim = Anam->ndim;
      mdim = Anam->Mshft;
       ofs = idf_dim_offset(Ndim,Mdim, style, ndim,mdim);
    if(ofs<0)
     {
      /* offset*/
      lvl = 1;
      ierr = 7-ofs; /*1-3*/
      break;
     }
    else
     {
      jelem = ofs;
      ibrk = ofs*jsize;
      val = val + ibrk;
     }

        kfrmt = idf_keyw_frmt(jfrmt,Anam->jtype);
    if(!kfrmt)
     {
      /* format type inconsistent with that assigned to name*/
      lvl = 1;
      ierr = 19;
      break;
     }

          ibrk=0;
    while(ibrk==0)
     {

           flag = idf_data_value();

        if(flag<0)
         {
          lvl  = 1-flag;
          ierr = 11-flag; /*1-3*/
          ibrk=2;
         }
        else if(flag==1)
         {
          /* no more data */
          ibrk=1;
         }
        else
         {
          /*data unit up */
             tag = val;
             jerr = idf_data_up(jfrmt,tag,nval);
          if(jerr)
           {
            /* conversion error*/
            lvl = 5;
            ierr=20;
            ibrk=2;
           }
          else
           {
            kup++;
            kkup++;

               jelem++;
            if(jelem >= nmax)
             {
              /* end of array is reached */
              ibrk=1;
             }
            else
             {
              val = val + jsize;
             }
           }
         }

     }/*ibrk search over data units */


    if(ierr)
     {
      flag=1;
     }
    else if(!kkup)
     {
      /* no data unit */
      lvl=1;
      ierr=16;
      flag=1;
     }
    else
     {
      flag=0;
     }

    if(flag) break;
   }

 } /* search over name list records */


if(!ierr)
 {
  if(krec==0)
   {
    /*no such name*/
    lvl=0;
    ierr = 15;
   }
  else if(kup==0)
   {
    /* no data units */
    lvl=0;
    ierr = 18;
   }
 }

err:
if(ierr)
 {
  idf_get_err_prn(name,nam,ierr,lvl);
  flag=0;
 }
else
 {
  flag = kup;
 }

return flag;
}
