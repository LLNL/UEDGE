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



#define IDF_PREFER_INT_IN_VOID 1




#if IDF_CPP_ == 1
int idf_str_to_v(char* str, int nn,
                 double* val, int* type)
#else
int idf_str_to_v(str,nn, val, type)
char   *str;
int     nn;
double *val;
int    *type;
#endif
{
/*
 -----------------------------------
    converts string to a number.
     the number type is guessed
  according to C language convention
 -----------------------------------
returns:
  0  OK
  1  no data symbols
  2  no digits
  3  improper symbol
  4  overflow
  5  improper number

type = 5 long
       7 double
       9 double complex
*/
int n,nk;
char c;
unsigned int ci;
int ic;
long l;
int ierr;
double dvl[2];
int typ;
register int i;
char *Str;

 dvl[0] = 0.0;
 dvl[1] = 0.0;

ierr = 0;
 typ = 0;

/* single digit */
if(nn == 1)
 {
  c = str[0];

  if( isdigit(c) )
   {
    ic = c - '0';
    dvl[0] = (double)ic;
    typ = IDF_FRMT_L;
    goto fin;
   }
  else
   {
    /* no digits */
    ierr=2;
    goto fin;
   }
 }


/* skip white insignificant spaces */
  i=nn;
while(i>0)
 {
  i--;

     c = str[i];
  if(c != ' ')
   {
    i++;
    break;
   }
 }
n = i;

if(i==0)
 {
  /* no digits */
  ierr=2;
  goto fin;
 }

if(i==1)
 {
  c = str[0];

  if( isdigit(c) )
   {
    ic = c - '0';
    dvl[0] = (double)ic;
    typ = IDF_FRMT_L;
    goto fin;
   }
  else
   {
    /* no digits */
    ierr=2;
    goto fin;
   }
 }

  i=0;
while(i<n)
 {
     c = str[i];
  if(c != ' ') break;
  i++;
 }

if(i==n)
 {
  /* no data symbols */
  ierr=1;
  goto fin;
 }


Str = str + i;
 nk = n-i;


/* single digit again */
if(nk == 1)
 {
  c = Str[0];

  if( isdigit(c) )
   {
    ic = c - '0';
    dvl[0] = (double)ic;
    typ = IDF_FRMT_L;
    goto fin;
   }
  else
   {
    /* no digits */
    ierr=2;
    goto fin;
   }
 }


/* check for possible complex */
ic = 0;
for(i=0;i<nk;i++)
 {
     c = Str[i];
  if(c == ',')
   {
    ic=1;
    break;
   }
 }
if(ic)
 {
   typ = IDF_FRMT_W;
  ierr = idf_str_to_z(Str,nk, dvl);
  goto fin; 
 }


/* check for possible floating point */
ic = 0;
for(i=0;i<nk;i++)
 {
     c = Str[i];
  if(c == '.')
   {
    ic=1;
    break;
   }
 }
if(ic)
 {
   typ = IDF_FRMT_D;
  ierr = idf_str_to_d(Str,nk, dvl);
  goto fin; 
 }


/* check for octal */
   c = Str[0];
if(c == '0')
 {
  ic = 0;
  for(i=1;i<nk;i++)
   {
       ci = (unsigned int) Str[i];
    if( !idf_octal_symbol(ci) )
     {
      ic=1;
      break;
     }
   }
  if(!ic)
   {
     typ = IDF_FRMT_L;
    ierr = idf_str_to_l(Str,nk, &l);
    dvl[0] = (double)l;
    goto fin;
   }

  /* check for hex */
     c = Str[1];
  if(c=='x' || c=='X')
   {
    if(nk==2)
     {
      /* hex zero*/
      typ = IDF_FRMT_L;
      goto fin;      
     }
    else
     {
      ic=0;
      for(i=2;i<nk;i++)
       {
        ci = (unsigned int) Str[i];
        if( !idf_hex_symbol(ci) )
         {
          ic=1;
          break;
         }
       }
      if(!ic)
       {
         typ = IDF_FRMT_L;
        ierr = idf_str_to_l(Str,nk, &l);
        dvl[0] = (double)l;
        goto fin;
       }
     }
   }
 }


/* check for possible exponent */
ic = 0;
for(i=0;i<nk;i++)
 {
  ci = (unsigned int) Str[i];

  if( idf_exp_symbol(ci) )
   {
    ic=1;
    break;
   }
 }
if(ic)
 {
   typ = IDF_FRMT_D;
  ierr = idf_str_to_d(Str,nk, dvl);
  goto fin;
 }

#if IDF_PREFER_INT_IN_VOID == 1
     typ = IDF_FRMT_L;
    ierr = idf_str_to_l(Str,nk, &l);
 if(ierr==4)
  {
   /*if overflow then get double */
    typ = IDF_FRMT_D;
   ierr = idf_str_to_d(Str,nk, dvl);
  }
 else
  {
   dvl[0] = (double)l;
  }
#else
  typ = IDF_FRMT_D;
 ierr = idf_str_to_d(Str,nk, dvl);
#endif

fin:
 *type = typ;
val[0] = dvl[0];
val[1] = dvl[1];

return ierr;
}
