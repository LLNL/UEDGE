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
int idf_dim(IDF_NAME_* A)
#else
int idf_dim(A)
IDF_NAME_ *A;
#endif
{
int n;

       n = A->ndim;
return n;
}



#if IDF_CPP_ == 1
int idf_dim_nmax(int Ndim, int* Mdim)
#else
int idf_dim_nmax(Ndim,Mdim)
int Ndim;
int *Mdim;
#endif
{
/*
  >0  nmax
  -1  no dimensions for target array
  -2  improper dimension in target array
*/
int k,flag;
register int i;

flag = 0;

if(Ndim < 1)
 {
  flag = -1;
 }

else
 {
  k=1;
  for(i=0;i<Ndim;i++)
   {
    if(Mdim[i] < 1)
     {
      flag=-2;
      break;
     }
    else
     {
      k = k * Mdim[i];
     }
   }/*i*/

  if(!flag) flag = k;
 }

return flag;
}



#if IDF_CPP_ == 1
int idf_dim_offset(int Ndim, int* Mdim, int style,
                   int ndim, int* mdim)
#else
int idf_dim_offset(Ndim,Mdim, style, ndim,mdim)
int  Ndim;
int *Mdim;
int  style;
int  ndim;
int *mdim;
#endif
{
/*
 >=0  offset
  -1  inconsistency in number of dimensions
  -2  improper dimension
  -3  offset is equal or exceed the maximal number of elements
*/
int k,flag,of,nmax;
register int i;

flag = 0;

if(ndim)
 {
  if(ndim != Ndim)
   {
    flag = -1;
   }

  else
   {
    for(i=0;i<ndim;i++)
     {
         k = mdim[i];
      if(k < 0 || k>=Mdim[i])
       {
        flag = -2;
        break;
       }
     }/*i*/

    if(!flag)
     {
      of  =0;
      k   =1;
      nmax=1;

      if(style)
       {
        /* Fortran*/
        for(i=0;i<ndim;i++)
         {
          of = of + mdim[i]*k;
          k = Mdim[i];
          nmax = nmax*k;
         }
       }

      else
       {
        /* C/C++ */
        i=ndim;
        while(i>0)
         {
          i--;
          of = of + mdim[i]*k;
          k = Mdim[i];
          nmax = nmax*k;
         }
       }

         k = of + 1;
      if(k >= nmax)
       {
        /* offset is equal or exceed the maximal number of elements*/
        flag=-3;
       }
      else
       {
        /* offset */
        flag = of;
       }
     }/*!flag*/

   }
 }

return flag;
}
