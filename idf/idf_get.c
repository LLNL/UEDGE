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




#define JTRC_S_IDF 0




#if IDF_CPP_ == 1
int idf_get(char* name, char* frmt, int align, void* value)
#else
int idf_get(name, frmt, align, value)
char *name;
char *frmt;
int   align;
void *value;
#endif
{
register int i;
int ierr,ibrk;
int nam,nfrmt;
char c;

  nam = strlen(name);

  i=0;
      ibrk=0;
while(ibrk==0)
 {
          c = frmt[i];
       if(c == IDF_EOS) ibrk=1;
  else if(c == ';'    ) ibrk=2;
  else                  i++;
 }

nfrmt = i;

#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get start: %s %d %s %d\n", name,nam,frmt,nfrmt);
#endif

    ierr = idf_get_(name,nam, frmt,nfrmt, align, value);

return ierr;
}



#if IDF_CPP_ == 1
int idf_get_(char* name, int nam,
             char* frmt, int nfrmt, int align,
             void* value)
#else
int idf_get_(name,nam, frmt,nfrmt, align, value)
char *name;
int   nam;
char *frmt;
int   nfrmt;
int   align;
void *value;
#endif
{
/*
-----------------------------------------------
import of packed structure of data units

'name'  is the string for named data
'frmt'  is the format symbol
'align' is the alignment rule code
'rwnd'  is the flag to rewind the format string
        if it get exhausted
'val'   is the pointer to the target
-----------------------------------------------

In success it returns the number of imported elements
and returns negative number at an error;
Caution: offset mode for array import and 
free memory allocation are unacceptable.
*/
static int rwnd  = 0;
static int style = 0;

int ierr=0;
int Nlist;
register int i;
int nval,jsize,jfrmt,kfrmt,kup,kkup,krec,jalloc;
int Ndim,ndim;
int *Mdim,*mdim;
int jgt,jelem,nmax,ofs;
int ibrk,jerr,flag,fflag,lvl;
IDF_NAME_ *Anam;
IDF_FRMT_TREAT *FrmtS;
void *val;
char *cval;

idf_frmt_clean();

FrmtS = idf_frmtreat_address();
        idf_frmtreat_nul();

lvl   = 0;
jsize = 0;
jfrmt = 0;
kup   = 0;
kkup  = 0;
krec  = 0;
nmax  = 0;
jelem = 0;

 val  = NULL;
cval  = (char*) value;

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
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=1\n");
#endif

   Nlist = idf_list_num();
if(Nlist < 1)
 {
  /* empty name list */
  ierr = 5;
  goto err;
 }
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=2\n");
#endif

   jerr = idf_frmt_prp(frmt,nfrmt);
if(jerr)
 {
  /* improper input format string*/
  ierr = 23;
  goto err;
 }
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=3\n");
#endif
#if IDF_FRMT_TRACE == 1
 idf_frmt_prn();
#endif



   jgt = idf_frmtreat_array(value);
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=4 jgt=%d\n",jgt);
#endif
if(jgt < 0)
{
 /* improper input format string*/
 ierr = 25;
 goto err;
}

else if(jgt>0)
{
/* 
---------------------------------
        array target 
  with array subscripts option 
---------------------------------
*/

jfrmt = FrmtS->jtype;
jsize = FrmtS->jsize;
nval  = jsize;
Ndim  = FrmtS->Ndimen;
Mdim  = FrmtS->Mdimen;

if(Ndim)
 {
     nmax = idf_dim_nmax(Ndim,Mdim);
  if(nmax<0)
   {
    /* improper dimensionality */
    ierr = 5-nmax; /*1-2*/
    goto err;
   }
 }
else
 {
  nmax=1;
 }
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=5\n");
#endif

krec=0;
kup =0;
for(i=0;i<Nlist;i++)
 {
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=6 krec=%d kup=%d\n",krec,kup);
#endif

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
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=7\n");
#endif
    cval = FrmtS->Target;

        ndim = Anam->ndim;
        mdim = Anam->Mshft;
       jelem = 0;

    if(Ndim)
     {
         ofs = idf_dim_offset(Ndim,Mdim, style, ndim,mdim);
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=70 ofs=%d\n",ofs);
#endif
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
        ibrk  = ofs*jsize;
        cval  = cval + ibrk;
       }
     }
    else
     {
      if(ndim)
       {
        ierr=9;
        break;
       }
     }
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=8\n");
#endif

        kfrmt = idf_keyw_frmt(jfrmt,Anam->jtype);
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=9 kfrmt=%d\n",kfrmt);
#endif
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
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=10 flag=%d\n",flag);
#endif
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
              val = cval;
             jerr = idf_data_up(jfrmt,val,nval);
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=100 jerr=%d\n",jerr);
#endif
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
              cval = cval + jsize;
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

}/*jgt>0*/ 



else
{
/* 
---------------------------------
      arbitrary target 
 without array subscripts option 
---------------------------------
*/

for(i=0;i<Nlist;i++)
 {
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=11\n");
#endif
     Anam = idf_list_get(name, nam, i);
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=12\n");
#endif 
  if(Anam != NULL)
   {
    krec++;
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=13 krec=%d\n",krec);
#endif
    if(krec>1)
     {
      lvl=1;
      ierr=27;
     }

    else if( idf_data_begin(Anam) )
     {
      /* unable to move to data field*/
      lvl = 1;
      ierr = 11;
     }

    else if( idf_frmt_begin(align,rwnd,cval) )
     {
      lvl = 1;
      ierr = 25;
     }

    else
     {
            ibrk=0;
      while(ibrk==0)
       {

           flag = idf_frmt_get();
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=14 flag=%d\n",flag);
#endif
        if(flag<0)
         {
          /* format getting error*/
          ierr=26;
          ibrk=2;
         }

        else if(flag==2)
         {
          /* no more format units*/
          ibrk=1;
         }

        else
         {
          jfrmt = FrmtS->jtype;
          nval  = FrmtS->jsize;
           val  = FrmtS->Target;
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=150 jfrmt=%d nval=%d\n",jfrmt,nval);
#endif
             fflag = idf_data_value();
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=15 fflag=%d\n",fflag);
#endif
          if(fflag<0)
           {
            lvl  = 1-fflag;
            ierr = 11-fflag; /*1-3*/
            ibrk=2;
           }
          else if(fflag==1)
           {
            /* no more data */
            ibrk=1;
           }
          else
           {
            /*data unit up */
               jerr = idf_data_up(jfrmt,val,nval);
#if JTRC_S_IDF == 1
fprintf(stdout,"idf_get trace, point=16 jerr=%d\n",jerr);
#endif
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
         }

       }/*ibrk search over data units */

     }

   }/*name matches*/


  if(ierr) break;

 } /* search over name list records */


if(!ierr)
 {
  if(krec==0)
   {
    /*no such name*/
    lvl=0;
    ierr = 15;
   }
 }

}/*jgt*/


err:
if(ierr)
 {
  idf_get_err_prn(name,nam,ierr,lvl);
  fflag = -ierr;
 }
else
 {
  fflag = kup;
 }

return fflag;
}
