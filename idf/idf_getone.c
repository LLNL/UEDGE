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
int idf_get_one(char* name, char frmt,
                int nval, void* val)
#else
int idf_get_one(name, frmt, nval, val)
char *name;
char  frmt;
int   nval;
void *val;
#endif
{
int ierr,nam;

nam = strlen(name);

ierr = idf_get_one_(name,nam, frmt, nval, val);

return ierr;
}



#if IDF_CPP_ == 1
int idf_get_one_(char* name, int nam, char frmt,
                 int nval, void* val)
#else
int idf_get_one_(name,nam, frmt, nval, val)
char *name;
int   nam;
char  frmt;
int   nval;
void *val;
#endif
{
/*
--------------------------------------------
   import of single data unit

'name' is the string for named data
'frmt' is the format symbol
'nval' is the maximal size of bite field
       if target is string or text
'val'  is the pointer to the target
--------------------------------------------

returns zero in success.
*/

int ierr=0;
int Nlist;
register int i;
unsigned int cfrmt;
int jfrmt,kfrmt;
int krec,ibrk,jerr,flag,lvl,kup;
IDF_NAME_ *Anam;


idf_frmt_clean();
idf_frmtreat_nul();

lvl=0;
cfrmt = frmt;
kfrmt = 0;

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
  /* improper input format symbol*/
  ierr = 4;
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

       jerr = idf_data_begin(Anam);
    if(jerr)
     {
      /* unable to move to data field*/
      lvl = 1;
      ierr = 11;
      break;
     }

       jerr = idf_dim(Anam);
    if(jerr)
     {
      /* dimensions are prescribed to single data unit*/
      lvl = 1;
      ierr = 17;
      break;
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
             jerr = idf_data_up(jfrmt,val,nval);
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
           }
         }

     }/*search over data units */

    if(ierr) break;
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
    ierr = 16;
   }
 }

err:
if(ierr)
 {
  idf_get_err_prn(name,nam,ierr,lvl);
 }

return ierr;
}
