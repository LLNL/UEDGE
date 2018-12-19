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



#if IDF_FLOAT_DEFAULT == 1

 static int     idf_float_iexp   = IDF_EXPOWTAB_MAX;
 static double  idf_float_Ntab[] = IDF_NTAB ;
 static double  idf_float_Ptab[] = IDF_PTAB ;


#else /* IDF_FLOAT_DEFAULT */

/*
 Warning: this option may lead
 to the loss of accuracy!
*/

 static int     idf_float_iexp = 0;

 static double  idf_float_Ntab = NULL;
 static double  idf_float_Ptab = NULL;



int idf_float_ini()
{
int ierr=0;
register int i;
double d;
int iexp,l,n;

if(idf_float_Ptab != NULL || idf_float_Ntab != NULL)
 {
  ierr=1;
  goto err;
 }

iexp = 1;
   l = 1;
   i = 0;
while(l < IDF_DBL_10_EXP_MAX)
 {
  iexp = l;
  l = iexp << 1;
  i++; 
 }

if(iexp < 2) {ierr=2; goto err;}

l = sizeof(double);

      n = i*l;
   idf_float_Ptab = (double*) malloc(n);
if(idf_float_Ptab == NULL) {ierr=1; goto err;}

      n = (i+1)*l;
   idf_float_Ntab = (double*) malloc(n);
if(idf_float_Ntab == NULL) {ierr=1; goto err;}


idf_float_iexp = iexp;

  n = i;
  d = 10.0;
while(i>0)
 {
  i--;
  idf_float_Ptab[i] = d;
  if(i) d = d*d;
 }

  i = n;
idf_float_Ntab[i] = 1.0;
  d = 0.1;
while(i>0)
 {
  i--;
  idf_float_Ntab[i] = d;
  if(i) d = d*d;
 }

err:
return ierr;
}


void idf_float_end()
{
 if(idf_float_Ptab != NULL)
  {
   free(idf_float_Ptab);
   idf_float_Ptab = NULL;
  }

 if(idf_float_Ntab != NULL)
  {
   free(idf_float_Ntab);
   idf_float_Ntab = NULL;
  }

 idf_float_iexp = 0;
}


#endif /* IDF_FLOAT_DEFAULT */



#if IDF_CPP_ == 1
void idf_float_decompose(double d, int* isig,
                         double* base, int* iexp)
#else
void idf_float_decompose(d, isig, base, iexp)
double  d;
int    *isig;
double *base;
int    *iexp;
#endif
{
/* 
    normalize double as
  (base=X.XXXX) (e+-) iexp
*/ 

double val;
register int i;
int iex,ip,sig;
int     PowM;
double *Ptab,*Ntab;

PowM = idf_float_iexp;
Ptab = idf_float_Ptab;
Ntab = idf_float_Ntab;

iex = 0;

  if(d == (double)0.0)
   {
    sig = 1;
    val = 0.0;
   }

  else if(d < 0.0)
   {
    sig = -1;
    val = -d;
   }

  else
   {
    sig = 1;
    val = d;
   }

 if(val < IDF_DOUBLE_MIN) goto fin;

 if(val < 1.0)
  {
   ip = PowM;
    i = 0;
   while(val < 0.1)
    {
     if(val <= Ntab[i])
      {
       val = val * Ptab[i];
       iex = iex - ip;
      }
     else
      {
       ip = ip >> 1;
       i++;
      }
    }

   if(val < 1.0)
    {
     val = val*10.0;
     iex = iex - 1;
    }
  }

 else if(val >= 10.0)
  {
   ip = PowM;
    i = 0;
   while(val >= 10.0)
    {
     if(val >= Ptab[i])
      {
       val = val * Ntab[i];
       iex = iex + ip;
      }
     else
      {
       ip = ip >> 1;
       i++;
      }
    }
/*
   if(val >= 10.0)
    {
     val = val*0.1;
     iex = iex + 1;
    }
*/
  }

fin:
*iexp = iex;
*base = val;
*isig = sig;
}



#if IDF_CPP_ == 1
int idf_float_compose(double base, int iex, double* d)
#else
int idf_float_compose(base, iex, d)
double  base;
int     iex;
double *d;
#endif
{
/* 
    normalize double as
  (base=X.XXXX) (e+-) iexp
*/ 

double val,g;
register int i;
int iexp,ip,iofl;
int     PowM;
double *Ptab,*Ntab;

PowM = idf_float_iexp;
Ptab = idf_float_Ptab;
Ntab = idf_float_Ntab;

 val = base;
iexp = iex;
iofl = 0;

if(iexp>0)
 {
  i=0;
  ip = PowM;
  while(iexp > 0)
   {
    if(iexp >= ip)
     {
               g = IDF_DOUBLE_MAX/Ptab[i];
      if(val < g)
       {
        val = val * Ptab[i];
        iexp = iexp - ip;
       }
      else
       {
        val = IDF_DOUBLE_MAX;
        iofl=1;
        break;
       }
     }
    else
     {
      ip = ip >> 1;
      i++;
     }
   }
 }

else if(iexp < 0)
 {
  i=0;
  ip = PowM;
  while(iexp < 0)
   {
    if(iexp <= -ip)
     {
      val = val * Ntab[i];
      iexp = iexp + ip;
     }
    else
     {
      ip = ip >> 1;
      i++;
     }
   }
 }

*d = val;

return iofl;
}
