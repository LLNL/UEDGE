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



#define IDF_EXACT_MATH  1



#if IDF_CPP_ == 1
int idf_mathd_add(double a, double b, double* c)
#else
int idf_mathd_add(a,b, c)
double  a;
double  b;
double *c;
#endif
{
int asig,bsig,csig;
int aexp,bexp,cexp;
double abase,bbase,cbase;
int iofl,jofl;
double val;

idf_float_decompose(a, &asig, &abase, &aexp);
idf_float_decompose(b, &bsig, &bbase, &bexp);

if(aexp > bexp)
 {
     csig  = asig*bsig;
     cexp  = bexp - aexp;
      iofl = idf_float_compose(bbase, cexp, &cbase);
  if(!iofl)
   {
       cbase = abase + cbase * (double)csig;
    if(cbase < 0.0)
     {
      cbase = -cbase;
      csig  = -asig;
     }
    else
     {
      csig = asig;
     }
    cexp  = aexp;
   }
 }

else if(aexp==bexp)
 {
  csig  = asig*bsig;
  iofl = 0;
       cbase = abase + bbase * (double)csig;
    if(cbase < 0.0)
     {
      cbase = -cbase;
      csig  = -asig;
     }
    else
     {
      csig = asig;
     }
  cexp  = aexp;
 }

else
 {
      csig  = asig*bsig;
      cexp  = aexp - bexp;
      iofl  = idf_float_compose(abase, cexp, &cbase);
  if(!iofl)
   {
       cbase = bbase + cbase * (double)csig;
    if(cbase < 0.0)
     {
      cbase = -cbase;
      csig  = -bsig;
     }
    else
     {
      csig = bsig;
     }
    cexp  = bexp;
   }
 }

if(!iofl)
 {
     jofl = idf_float_compose(cbase, cexp, &val);
  if(jofl)
   {
    if(csig < 0) val = -IDF_DOUBLE_MAX;
    else         val =  IDF_DOUBLE_MAX;    
   }
  else
   {
#if IDF_EXACT_MATH == 1
    val = a+b;
#else
    if(csig < 0) val = -val;
#endif
   }
 }
else
 {
  jofl = 1;
  if(csig < 0) val = -IDF_DOUBLE_MAX;
  else         val =  IDF_DOUBLE_MAX;
 }

 *c = val;

    if(jofl) jofl=IDF_MATH_OVERFLOW;
return jofl;
}



#if IDF_CPP_ == 1
int idf_mathd_sub(double a, double b, double* c)
#else
int idf_mathd_sub(a,b, c)
double  a;
double  b;
double *c;
#endif
{
int iofl;
double d;

d = -b;

       iofl = idf_mathd_add(a,b, c);
return iofl;
}



#if IDF_CPP_ == 1
int idf_mathd_mul(double a, double b, double* c)
#else
int idf_mathd_mul(a,b, c)
double  a;
double  b;
double *c;
#endif
{
int asig,bsig,csig;
int aexp,bexp,cexp;
double abase,bbase,cbase;
int iofl;
double val;

idf_float_decompose(a, &asig, &abase, &aexp);
idf_float_decompose(b, &bsig, &bbase, &bexp);

     csig  = asig*bsig;
     cexp  = aexp + bexp;
     cbase = abase * bbase;

     iofl = idf_float_compose(cbase, cexp, &val);
  if(iofl)
   {
    if(csig < 0) val = -IDF_DOUBLE_MAX;
    else         val =  IDF_DOUBLE_MAX;    
   }
  else
   {
#if IDF_EXACT_MATH == 1
    val = a*b;
#else
    if(csig < 0) val = -val;
#endif
   }

 *c = val;
    if(iofl) iofl=IDF_MATH_OVERFLOW;
return iofl;
}



#if IDF_CPP_ == 1
int idf_mathd_div(double a, double b, double* c)
#else
int idf_mathd_div(a,b, c)
double  a;
double  b;
double *c;
#endif
{
int asig,bsig,csig;
int aexp,bexp,cexp;
double abase,bbase,cbase;
int iofl;
double val;

idf_float_decompose(a, &asig, &abase, &aexp);
idf_float_decompose(b, &bsig, &bbase, &bexp);

if(bbase < IDF_DOUBLE_MIN)
 {
  val = IDF_DOUBLE_MAX;
  iofl=IDF_MATH_EXCEPTION;
 }

else
 {
     csig  = asig*bsig;
     cexp  = aexp - bexp;
     cbase = abase / bbase;

     iofl = idf_float_compose(cbase, cexp, &val);
  if(iofl)
   {
    if(csig < 0) val = -IDF_DOUBLE_MAX;
    else         val =  IDF_DOUBLE_MAX;
    iofl=IDF_MATH_OVERFLOW;    
   }
  else
   {
#if IDF_EXACT_MATH == 1
    val = a/b;
#else
    if(csig < 0) val = -val;
#endif
   }
 }

 *c = val;
return iofl;
}



#if IDF_CPP_ == 1
int idf_mathd_pow(double a, double b, double* c)
#else
int idf_mathd_pow(a,b, c)
double  a;
double  b;
double *c;
#endif
{
double xa,xb, x,y;
int iofl;
int sig,k;
double val;

iofl=0;

   xa = 1.0-a;
if(xa < 0.0) xa = -xa;

if(xa < IDF_DOUBLE_MIN)
 {
  val = 1.0;
  goto fin;
 }

   xa = 1.0-b;
if(xa < 0.0) xa = -xa;

if(xa < IDF_DOUBLE_MIN)
 {
  val = a;
  goto fin;
 }

if(a < 0.0) xa=-a; 
else        xa=a;

if(b < 0.0) xb=-b; 
else        xb=b;

if(xa < IDF_DOUBLE_MIN)
 {
  if(xb < IDF_DOUBLE_MIN)
   {
    val = 0.0;
   }
  else if(b < 0.0)
   {
    val = IDF_DOUBLE_MAX;
    iofl=IDF_MATH_DOMAIN;
   }
  else
   {
    val = 0.0;
   }
 }

else if(xb < IDF_DOUBLE_MIN)
 {
  val = 1.0;
 }

else if(a < 0.0)
 {
  if(xb > IDF_LONG_MAX)
   {
    val = IDF_DOUBLE_MAX;
    iofl=IDF_MATH_DOMAIN;
   }
  else
   {
       k = (int)b;
       y = (double)k;
       y = y-b;
    if(y < 0.0) y=-y;
    if(y < IDF_DOUBLE_MIN)
     {
      if(k==0)
       {
        val = 1.0;
       }
      else if(k==1)
       {
        val = a;
       }
      else if(k==-1)
       {
        val = 1.0/a;
       }
      else
       {
           xa = -a;
        if(xa < 1.0)
         {
          x = 1./xa;
          k = -k;
         }
        else
         {
          x = xa;
         }

         y = (double)k;
        xa = log(x);

        if(y > 0.0)
         {
             iofl = idf_mathd_mul(xa,y, &xb);
          if(iofl)
           {
            val = IDF_DOUBLE_MAX;
            iofl=IDF_MATH_RANGE;
           }
          else if(xb > IDF_LNMAXDOUBLE)
           {
            val = IDF_DOUBLE_MAX;
            iofl=IDF_MATH_RANGE;
           }
          else
           {
            val = pow(x,y);
            iofl=0;
           }    
         }
        else
         {
             iofl = idf_mathd_mul(xa,-y, &xb);
          if(iofl)
           {
            val = 0.0;
           }
          else if(xb > IDF_LNMAXDOUBLE)
           {
            val = 0.0;
           }
          else
           {
            val = pow(x,y);
           }
          iofl=0;
         }

           sig = k & 1;
        if(sig) val = -val;
       }
     }

    else
     {
      val = IDF_DOUBLE_MAX;
      iofl=IDF_MATH_RANGE;
     } 
   }
 }

else
 {
  /* a>0 */
  if(a < 1.0)
   {
    x = 1.0/a;
    y = -b;
   }
  else
   {
    x = a;
    y = b;
   }

  xa = log(x);

  if(y > 0.0)
   {
       iofl = idf_mathd_mul(xa,y, &xb);
    if(iofl)
     {
      val = IDF_DOUBLE_MAX;
      iofl=IDF_MATH_RANGE;
     }
    else if(xb > IDF_LNMAXDOUBLE)
     {
      val = IDF_DOUBLE_MAX;
      iofl=IDF_MATH_RANGE;
     }
    else
     {
      val = pow(x,y);
      iofl=0;
     }    
   }
  else
   {
       iofl = idf_mathd_mul(xa,-y, &xb);
    if(iofl)
     {
      val = 0.0;
     }
    else if(xb > IDF_LNMAXDOUBLE)
     {
      val = 0.0;
     }
    else
     {
      val = pow(x,y);
     }
    iofl=0;
   }
 }

fin:
 *c = val;
return iofl;
}
