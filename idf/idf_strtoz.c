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
int idf_str_to_z(char* str, int nn, double* dval)
#else
int idf_str_to_z(str,nn, dval)
char   *str;
int     nn;
double *dval;
#endif
{
/*
 -----------------------------------
  converts string to complex double
 -----------------------------------
returns:
  0  OK
  1  no data symbols
  2  no digits
  3  improper symbol
  4  overflow
  5  improper number
*/
char c;
int k,kk,kc,flag;
register int i;
double d;
char *Str;

   flag = 0;
dval[0] = (double)0.0;
dval[1] = (double)0.0;

  k=0;
  i=0;
while(i<nn)
 {
     c = str[i];
  if(c == ',')
   {
    k=i;
    break;
   }
  i++;
 }

if(k==0||k==nn)
 {
  flag = 5;
 }
else
 {
  Str = str;

     kk = idf_str_to_d(Str,k, &d);
  if(kk)
   {
    flag = kk;
   }
  else
   {
    dval[0] = d;

    k++;
    Str = str + k;
    k = nn-k;

       kc = idf_str_to_d(Str,k, &d);
    if(kc)
     {
      flag = kc;
     }
    else
     {
      dval[1] = d;
     }
   }
 }

return flag;
}
