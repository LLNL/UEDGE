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
int idf_str_to_l(char* str, int nn, long* lval)
#else
int idf_str_to_l(str,nn, lval)
char *str;
int   nn;
long *lval;
#endif
{
/*
 -----------------------------------
   converts string to long integer
 -----------------------------------
returns:
  0  OK
  1  no data symbols
  2  no digits
  3  improper symbol
  4  overflow
  5  improper number
*/
static long lmax = IDF_LONG_MAX;

long ig,ival;
char c;
int ierr,sgn,jtyp,base,ic,n;
register int i;

ival = 0L;
ierr = 0;
sgn  = 0;
  c  = '\0';

if(nn == 1)
 {
  c = str[0];

  if( isdigit(c) )
   {
    ic = c - '0';
    ival = ival + ic;
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

if(i == 0)
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
    ival = ival + ic;
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


if( isdigit(c) )
 {
  if(c == '0')
   {
       i++;
    if(i<n)
     {
      c = str[i];

      if((c=='x')||(c=='X'))
       {
        /* hex */
        jtyp=2;
        base=16;

           i++;
        if(i<n)
         {
          c = str[i];
         }
        else
         {
          /* hex zero */
          goto fin;
         }
       }
      else
       {
        /* octal */
        jtyp=1;
        base=8;
       }      
     }
    else
     {
      /* octal zero */
      goto fin;
     }
   }

  else
   {
    /* decimal */
    jtyp = 0;
    base = 10;

    ic = c - '0';
    ival = ival + ic;

       i++;
    if(i<n)
     {
      c = str[i];
     }
    else
     {
      /* decimal zero*/
      goto fin;
     }
   }/*decimal*/
 }

else
 {
  /* improper symbol */
  ierr=3;
  goto err;
 }



if(jtyp == 0)
 {
  /* decimal */

  do {
      if( isdigit(c) )
       {
        ic = c - '0';
        ig = (lmax-ic)/base;

        if(ival > ig)
         {
          /*overflow */
          ierr = 4;
          ival = lmax;
          break;
         }

        else
         {
          ival = ival*base + ic;
         }
       }

      else
       {
        ierr=3;
        break;
       }

              i++;
      c = str[i];

     } while(i<n);

 } /* end of decimal */


else if(jtyp==2)
 {
  /* hex */

  do {
      if( isdigit(c) )
       {
        ic = c - '0';
       }
      else if( islower(c) )
       {
        ic = c - 'a';
        ic = ic + 10;
       }
      else if( isupper(c) )
       {
        ic = c - 'A';
        ic = ic + 10;
       }
      else
       {
        ierr=3;
        break;
       }

      if(ic < base)
       {
        ig = (lmax-ic)/base;

        if(ival > ig)
         {
          /*overflow */
          ierr = 4;
          ival = lmax;
          break;
         }

        else
         {
          ival = ival*base + ic;
         }
       }
      else
       {
        /* out of range */
        ierr=5;
        break;
       }

              i++;
      c = str[i];

     } while(i<n);

 }/* end of hex */


else
 {
  /* octal */

  do {
      if( isdigit(c) )
       {
        ic = c - '0';

        if(ic < base)
         {
          ig = (lmax-ic)/base;

          if(ival > ig)
           {
            /*overflow */
            ierr = 4;
            ival = lmax;
            break;
           }

          else
           {
            ival = ival*base + ic;
           }
         }

        else
         {
          /* out of range */
          ierr=5;
          break;
         }
       }

      else
       {
        ierr=3;
        break;
       }

              i++;
      c = str[i];

     } while(i<n);

 } /* end of octal */



err: if(sgn) ival = -ival;

fin: *lval = ival;

return ierr;
}
