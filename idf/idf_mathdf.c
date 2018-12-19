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
int idf_mathd_sin(double x, double* y)
#else
int idf_mathd_sin(x,y)
double x;
double *y;
#endif
{
int iofl,sig;
double xx,f;

 if(x < 0.0)
  {
   sig=1;
   xx = -x;
  }
 else
  {
   sig=0;
   xx=x;
  }

 if(xx > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   f = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }
 else if(xx < IDF_DOUBLE_MIN)
  {
   f = 0.0;
   iofl=0;
  }
 else
  {
   f = sin(xx);
   if(sig) f = -f;
   iofl=0;
  }

*y = f;
return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_cos(double x, double* y)
#else
int idf_mathd_cos(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx,f;

 if(x < 0.0) xx = -x; else xx=x;
 if(xx > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   f = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }
 else if(xx < IDF_DOUBLE_MIN)
  {
   f = 1.0;
   iofl=0;
  }
 else
  {
   f = cos(xx);
   iofl=0;
  }

*y = f;
return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_tan(double x, double* y)
#else
int idf_mathd_tan(x,y)
double x;
double *y;
#endif
{
int iofl,sig;
double xx,f;

 if(x < 0.0)
  {
   sig=1;
   xx = -x;
  }
 else
  {
   sig=0;
   xx=x;
  }

 if(xx > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   f = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }
 else if(xx < IDF_DOUBLE_MIN)
  {
   f = 0.0;
   iofl=0;
  }
 else
  {
      xx = cos(x);
   if(xx < 0.0) x = -xx;
   if(xx < IDF_DOUBLE_MIN)
    {
     /* improper argument */
     f = IDF_DOUBLE_MAX;
     iofl=IDF_MATH_DOMAIN;
    }
   else
    {
     f = tan(xx);
     if(sig) f = -f;
     iofl=0;
    }
  }

*y = f;
return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_sqrt(double x, double* y)
#else
int idf_mathd_sqrt(x,y)
double x;
double *y;
#endif
{
int iofl;

 if(x < 0.0)
  {
   /* improper argument */
   *y = 0.0;
   iofl=IDF_MATH_DOMAIN;
  }
 else if(x > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   *y = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }
 else
  {
   *y = sqrt(x);
   iofl=0;
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_log(double x, double* y)
#else
int idf_mathd_log(x,y)
double x;
double *y;
#endif
{
int iofl;

 if(x < IDF_DOUBLE_MIN)
  {
   /* improper argument */
   *y = -IDF_DOUBLE_MAX;
   iofl=IDF_MATH_DOMAIN;
  }
 else if(x > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   *y = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }
 else
  {
   *y = log(x);
   iofl=0;
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_log10(double x, double* y)
#else
int idf_mathd_log10(x,y)
double x;
double *y;
#endif
{
int iofl;

 if(x < IDF_DOUBLE_MIN)
  {
   /* improper argument */
   *y = -IDF_DOUBLE_MAX;
   iofl=IDF_MATH_DOMAIN;
  }
 else if(x > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   *y = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }
 else
  {
   *y = log10(x);
   iofl=0;
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_asin(double x, double* y)
#else
int idf_mathd_asin(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx;

 if(x < 0.0) xx = -x; else xx=x;
 if(xx > 1.0)
  {
   /* improper argument */
   if(x < 0.0) xx=-1.0; else xx=1.0;
   *y = xx;
   iofl=IDF_MATH_DOMAIN;
  }
 else
  {
   *y = asin(x);
   iofl=0;
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_acos(double x, double* y)
#else
int idf_mathd_acos(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx;

 if(x < 0.0) xx = -x; else xx=x;
 if(xx > 1.0)
  {
   /* improper argument */
   if(x < 0.0) xx=-1.0; else xx=1.0;
   *y = xx;
   iofl=IDF_MATH_DOMAIN;
  }
 else
  {
   *y = acos(x);
   iofl=0;
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_atan(double x, double* y)
#else
int idf_mathd_atan(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx;

 if(x < 0.0) xx = -x; else xx=x;
 if(xx > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   *y = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }
 else
  {
   *y = atan(x);
   iofl=0;
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_sinh(double x, double* y)
#else
int idf_mathd_sinh(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx,f,g,sig;

 if(x < 0.0)
  {
   sig = -1.0;
    xx = -x;
  }
 else
  {
   sig = 1.0;
    xx = x;
  }

 if(xx > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   f = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }

 else if(xx < IDF_DOUBLE_MIN)
  {
   f = 0.0;
   iofl=0;
  }

 else if(xx < 1.0)
  {
   f = exp(xx);
   g = f - 1.0/f;
   f = g*0.5*sig;

   iofl=0;
  }

 else
  {
      iofl = idf_mathd_add(xx,xx, &f);
   if(iofl)
    {
     f = 0.5;
     iofl=0;
    }

   else
    {
        g = xx + xx;
     if(g > IDF_LNMAXDOUBLE)
      {
       f = 0.5;
      }
     else
      {
       f = 1.0 - exp(-g);
       f = f*0.5;
      }

        g = xx + log(f);
     if(g > IDF_LNMAXDOUBLE)
      {
       f = IDF_DOUBLE_MAX;
       iofl=IDF_MATH_RANGE;
      }
     else
      {
       f = f*exp(xx);
       f = f*sig;
      }
    }
  }

*y = f;
return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_cosh(double x, double* y)
#else
int idf_mathd_cosh(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx,f,g;

 if(x < 0.0) xx = -x;
 else        xx =  x;

 if(xx > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   f = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }

 else if(xx < IDF_DOUBLE_MIN)
  {
   f = 0.0;
   iofl=0;
  }

 else if(xx < 1.0)
  {
   f = exp(xx);
   g = f + 1.0/f;
   f = g*0.5;
   iofl=0;
  }

 else
  {
      iofl = idf_mathd_add(xx,xx, &f);
   if(iofl)
    {
     f = 0.5;
     iofl=0;
    }

   else
    {
        g = xx + xx;
     if(g > IDF_LNMAXDOUBLE)
      {
       f = 0.5;
      }
     else
      {
       f = 1.0 + exp(-g);
       f = f*0.5;
      }

        g = xx + log(f);
     if(g > IDF_LNMAXDOUBLE)
      {
       f = IDF_DOUBLE_MAX;
       iofl=IDF_MATH_RANGE;
      }
     else
      {
       f = f*exp(xx);
      }
    }
  }

*y = f;
return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_tanh(double x, double* y)
#else
int idf_mathd_tanh(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx,f,g,sig;

 if(x < 0.0)
  {
   sig = -1.0;
    xx = -x;
  }
 else
  {
   sig = 1.0;
    xx = x;
  }

 if(xx > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   f = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }

 else if(xx < IDF_DOUBLE_MIN)
  {
   f = 0.0;
   iofl=0;
  }

 else if(xx < 1.0)
  {
   f = xx + xx;
   g = exp(-f);
   f = 1.0 - g;
   f = f/(1.0+g);
   f = f*sig;
   iofl=0;
  }

 else
  {
      iofl = idf_mathd_add(xx,xx, &f);
   if(iofl)
    {
     f = 0.5;
     iofl=0;
    }

   else
    {
        g = xx + xx;
     if(g > IDF_LNMAXDOUBLE)
      {
       f = 0.0;
      }
     else
      {
       f = exp(-g);
       g = f / (1.0 + f);
       f = 0.5 - g;
       g = f+f;
       f = g*sig;
      }
    }
  }

*y = f;
return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_atanh(double x, double* y)
#else
int idf_mathd_atanh(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx,sig;
double f,g,h;

 if(x < 0.0)
  {
    xx = -x;
   sig = -1.0;
  } 
 else
  { 
    xx =  x;
   sig = 1.0;
  }

 if(xx > 1.0)
  {
   /* improper argument */
   *y = x;
   iofl=IDF_MATH_DOMAIN;
  }

 else
  {
   /* atanh(x0 = 0.5 ln [(1+x)/(1-x)] */

   /* *y = atanh(x); */

   f = 1.0 - xx;
   g = 1.0 + xx;

      iofl = idf_mathd_div(g,f, &h);
   if(iofl)
    {
     *y = 0.0;
     iofl = IDF_MATH_RANGE;
    }
   else
    {
     g = log(h)*0.5;
     f = g*sig;

     *y = f;
     iofl=0;
    }
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_asinh(double x, double* y)
#else
int idf_mathd_asinh(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx,sig;
double f,g,h;

 if(x < 0.0)
  {
    xx = -x;
   sig = -1.0;
  } 
 else
  { 
    xx =  x;
   sig = 1.0;
  }

 if(xx > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   
   *y = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_DOMAIN;
  }

 else
  {
   /* asinh(x) = ln[x + sqrt(x^2+1)] */

   /* *y = asinh(x); */

   if(xx > 1.0)
    {
     f = 1.0 + 1.0/xx;
     g = f*f;
     f = 1.0 + sqrt(g);
     g = log(f);
     f = log(xx);

        iofl = idf_mathd_add(f,g, &h);
     if(iofl)
      {
       *y = 0.0;
       iofl = IDF_MATH_RANGE;
      }
     else
      {
        f = h*sig;
       *y = f;
      }
    }

   else
    {
     f = xx*xx + 1.0;
     g = x + sqrt(f);
     h = log(g);
     f = h*sig;

     *y = f;
     iofl = 0;
    }
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_acosh(double x, double* y)
#else
int idf_mathd_acosh(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx,f,g;


 if(x < 1.0)
  {
   *y = 0.0;
   iofl=IDF_MATH_DOMAIN;
  }

 else
  {
      x = -x;
   if(x > IDF_DOUBLE_MAX)
    {
     *y = IDF_DOUBLE_MAX;
     iofl=IDF_MATH_DOMAIN;
    }
   else
    {
     /* acosh(x) = ln[x+sqrt(x^2-1)] */

     /*  *y = acosh(x); */

     xx = 1.0/x;
      f = (1.0-xx)*(1.0+xx);
     xx = 1.0 + sqrt(f);
      f = log(xx);
     xx = log(x);

        iofl = idf_mathd_add(xx,f, &g);
     if(iofl)
      {
       *y = 0.0;
       iofl = IDF_MATH_RANGE;
      }
     else
      {
       *y = g;
      }
    }
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_abs(double x, double* y)
#else
int idf_mathd_abs(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx;

 if(x < 0.0) xx = -x; else xx=x;
 if(xx > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   *y = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }
 else
  {
   *y = xx;
   iofl=0;
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_int(double x, double* y)
#else
int idf_mathd_int(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx;

 if(x < 0.0) xx = -x; else xx=x;
 if(xx > IDF_INT_MAX)
  {
   /* improper argument */
   *y = IDF_INT_MAX;
   iofl=IDF_MATH_DOMAIN;
  }
 else
  {
   *y = (int)x;
   iofl=0;
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_hypot(double xx, double yy, double* d)
#else
int idf_mathd_hypot(xx,yy, d)
double xx,yy;
double *d;
#endif
{
int iofl;
double x,y,t;

iofl = 0;

if(xx < 0.0) x=-xx; else x=xx;
if(yy < 0.0) y=-yy; else y=yy;

if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *d = 0.0;
 }

else if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *d = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_RANGE;
 }

else if(x > y)
 {
  if(x < 1.0)
   {
     t = x*x + y*y;
    *d = sqrt(t);
   }
  else
   {
     t = y/x;
    *d = x * sqrt(1.0 + t*t);
   }
 }
else
 {
  if(y < 1.0)
   {
     t = x*x + y*y;
    *d = sqrt(t);
   }
  else
   {
     t = x/y;
    *d = y * sqrt(1.0 + t*t);
   }
 }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_min(double x, double y, double* z)
#else
int idf_mathd_min(x,y, z)
double x,y;
double *z;
#endif
{
int iofl;
double xx,yy;

 if(x < 0.0) xx = -x; else xx=x;
 if(y < 0.0) yy = -y; else yy=y;

 if(xx > IDF_DOUBLE_MAX || yy > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   *z = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }
 else
  {
   if(x < y) xx=x; else xx=y;
   *z = xx;
   iofl=0;
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_max(double x, double y, double* z)
#else
int idf_mathd_max(x,y, z)
double x,y;
double *z;
#endif
{
int iofl;
double xx,yy;

 if(x < 0.0) xx = -x; else xx=x;
 if(y < 0.0) yy = -y; else yy=y;

 if(xx > IDF_DOUBLE_MAX || yy > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   *z = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_RANGE;
  }
 else
  {
   if(x > y) xx=x; else xx=y;
   *z = xx;
   iofl=0;
  }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathd_exp(double x, double* y)
#else
int idf_mathd_exp(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx;

 if(x < 0.0)
  {
      xx = -x;
   if(xx > IDF_LNMAXDOUBLE)
    {
     *y = 0.0;
     iofl=0;
    }
   else
    {
     *y = exp(x);
     iofl=0;
    }
  }

 else
  {
   if(x > IDF_LNMAXDOUBLE)
    {
     *y = IDF_DOUBLE_MAX;
     iofl=IDF_MATH_DOMAIN;
    }
   else
    {
     *y = exp(x);
     iofl=0;
    }
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_sq1mz2(double x, double* y)
#else
int idf_mathd_sq1mz2(x,y)
double x;
double *y;
#endif
{
/* sqrt(1-z^2) */
int iofl;
double xx,z;

 if(x > 1.0)
  {
   /* improper argument */
   *y = 0.0;
   iofl=IDF_MATH_DOMAIN;
  }

 else
  {
   if(x < 0.0) xx=-x;
   else        xx=x;

   if(xx > IDF_DOUBLE_MAX)
    {
     /* improper argument */
     *y = IDF_DOUBLE_MAX;
     iofl=IDF_MATH_RANGE;
    }

   else if(xx < IDF_DOUBLE_MIN)
    {
     *y = 1.0;
     iofl=0;
    }

   else
    {
        z = 1.0 - xx*xx;
     if(z < 0.0) z = 0.0;

     *y = sqrt(z);
     iofl=0;
    }
  }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_sq1pz2(double x, double* y)
#else
int idf_mathd_sq1pz2(x,y)
double x;
double *y;
#endif
{
/* sqrt(1+z^2) */
int iofl;
double xx,z;

   if(x < 0.0) xx=-x;
   else        xx=x;

   if(xx > IDF_DOUBLE_MAX)
    {
     /* improper argument */
     *y = IDF_DOUBLE_MAX;
     iofl=IDF_MATH_RANGE;
    }

   else if(xx < IDF_DOUBLE_MIN)
    {
     *y = 1.0;
     iofl=0;
    }

   else if(xx < 1.0)
    {
      z = 1.0 + xx*xx;
     *y = sqrt(z);
     iofl=0;
    }

   else
    {
        iofl = idf_mathd_div(1.0,xx,&z);
     if(iofl)
      {
       /* improper argument */
       *y = IDF_DOUBLE_MAX;
       iofl=IDF_MATH_RANGE;
      }
     else
      {
       xx = 1.0 + z*z;
        z = z/sqrt(xx);

       *y = z;
       iofl=0;
      }
    }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathd_pow10(double x, double* y)
#else
int idf_mathd_pow10(x,y)
double x;
double *y;
#endif
{
int iofl;
double xx;

 if(x < 0.0)
  {
      xx = -x;
   if(xx > IDF_LNMAXDOUBLE)
    {
     *y = 0.0;
     iofl=0;
    }
   else
    {
     xx = x*IDF_LN10;
     *y = exp(xx);
     iofl=0;
    }
  }

 else
  {
   if(x > IDF_LNMAXDOUBLE)
    {
     *y = IDF_DOUBLE_MAX;
     iofl=IDF_MATH_DOMAIN;
    }
   else
    {
        xx = x*IDF_LN10;
     if(xx > IDF_LNMAXDOUBLE)
      {
       *y = IDF_DOUBLE_MAX;
       iofl=IDF_MATH_DOMAIN;
      }
     else
      {
       *y = exp(xx);
       iofl=0;
      }
    }
  }

return iofl;
}
