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




#if IDF_CPP_ == 1
int idf_str_to_d(char* str, int nn, double* dval)
#else
int idf_str_to_d(str,nn, dval)
char   *str;
int     nn;
double *dval;
#endif
{
/*
 -----------------------------
   converts string to double
 -----------------------------
returns:
  0  OK
  1  no data symbols
  2  no digits
  3  improper symbol
  4  overflow
  5  improper number
*/

static long lmax = IDF_LONG_MAX;

long lmin,lmscl;
long iscl,jscl,kscl;
long iexp,jexp,ig;
double dvl,g;
char c;
int ierr,sgn,sge,ip,ic,k,n;
register int i;


 dvl = 0.0;
ierr = 0;
sgn  = 0;
  c  = '\0';

if(nn == 1)
 {
  c = str[0];

  if( isdigit(c) )
   {
    ic = c - '0';
    dvl = dvl + (double)ic;
    goto fin;
   }
  else
   {
    /* no digits */
    ierr=2;
    goto err;
   }
 }


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
  goto err;
 }

if(i==1)
 {
  c = str[0];

  if( isdigit(c) )
   {
    ic = c - '0';
    dvl = dvl + (double)ic;
    goto fin;
   }
  else
   {
    /* no digits */
    ierr=2;
    goto err;
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
  goto err;
 }


if(c == '+')
 {
     i++;
  if(i<n)
   {
    c = str[i];
   }
  else
   {
    ierr=2;
    goto err;
   }
 }
else if(c == '-')
 {
     i++;
  if(i<n)
   {
    sgn=1;
    c = str[i];
   }
  else
   {
    /* no digits */
    ierr=2;
    goto err;
   }
 }

  ip = 0;
if(c == '.')
 {
     i++;
  if(i<n)
   {
    c = str[i];
    ip++;
   }
  else
   {
    ierr=2;
    goto err;
   }
 }

if( !isdigit(c) )
 {
  /* should have at least one digit */
  ierr=2;
  goto err;
 }


 lmin = -lmax;
lmscl = (lmax - 10)/10;

 iscl = 0;
 jscl = 1;
 kscl = 0;
 iexp = 0;

    k = 0;
do {
    if(iscl < lmscl)
     {
      ic = c - '0';
      iscl = iscl*10 + ic;
      if(ip)
       {
        if(iexp == lmin)
         {
          /*improper number*/
          ierr=5;
          k=3;
         }
        else
         {
          iexp--;
         }
       }
     }
    else if(jscl < lmscl)
     {
      ic = c - '0';
      kscl = kscl*10 + ic;
      jscl = jscl*10;
      if(ip)
       {
        if(iexp == lmin)
         {
          /*improper number*/
          ierr=5;
          k=3;
         }
        else
         {
          iexp--;
         }
       }
     }
    else
     {
      if(ip==0)
       {
        if(iexp == lmax)
         {
          /* overflow*/
          ierr=4;
          k=3;
         }
        else
         {
          iexp++;
         }
       }
     }

       i++;
    if(i<n)
     {
      c = str[i];

      if(c == '.')
       {
        if(ip)
         {
          /* multiple dot */
          ierr=3;
          k=3;
         }
        else
         {
          ip++;

             i++;
          if(i<n)
           {
            c = str[i];
           }
          else
           {
            /* no exponent */
            k=2;
            break;
           }
         }
       }

      if( !isdigit(c) )
       {
             if(c == 'e') k=1;
        else if(c == 'E') k=1;
        else if(c == 'd') k=1;
        else if(c == 'D') k=1;
        else if(c == 'g') k=1;
        else if(c == 'G') k=1;
        else
         {
          /* improper symbol */
          ierr=3;
          k=3;
         }

        if(k==1)
         {
             i++;
          if(i<n)
           {
            c = str[i];
           }
          else
           {
            /* empty Exponent */
            k=2;
           }
         }
       }
     }

    else
     {
      /* number finished without Exp*/
      k=2;
     }

   } while(k==0);


if(k==3) goto err;


if(k==1)
 {
  /* get exponent transform */
  sge=0;

  if(c == '+')
   {
       i++;
    if(i<n)
     {
      c = str[i];
     }
    else
     {
      ierr=5;
      goto err;
     }
   }
  else if(c == '-')
   {
       i++;
    if(i<n)
     {
      sge=1;
      c = str[i];
     }
    else
     {
      ierr=5;
      goto err;
     }
   }

 jexp = 0;

      k=0;
  do {
      if( isdigit(c) )
       {
        ic = c - '0';
        ig = (lmax-ic)/10;

        if(jexp > ig)
         {
          /* improper number */
          ierr=5;
          k=2;
         }
        else
         {
          jexp = jexp*10 + ic;

             i++;
          if(i < n)
           {
            c = str[i];
           }
          else
           {
            k=1;
           }
         }
       }
      else
       {
        ierr=5;
        k=2;
       }
       
     } while(k==0);

  if(k==2) goto err;

  if(sge) jexp = -jexp;
  iexp = iexp + jexp;  

 }/* get Exp*/


            k=0;
if(iscl==0) k++;
if(jscl==0) k++;
if(kscl==0) k++;
if(k != 3)
 {
  /* non zero */

    dvl = (double)iscl;
  if(jscl != 1)
   {
    dvl = dvl*(double)jscl + (double)kscl;
   } 

    ip = idf_float_compose(dvl, iexp, &g);
 if(ip)
  {
   ierr = 4;
  }
  dvl = g;
 }


err: if(sgn) dvl = -dvl;

fin: *dval = dvl;

return ierr;
}
