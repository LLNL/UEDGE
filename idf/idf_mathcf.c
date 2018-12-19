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
int idf_mathc_cmplx(double xr, double xi,
                    double* yr, double* yi)
#else
int idf_mathc_cmplx(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,t,d;

iofl = 0;

 if(xr < 0.0) x = -xr; else x=xr;
 if(xi < 0.0) y = -xi; else y=xi;

if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  /* improper argument */
  iofl=IDF_MATH_DOMAIN;
 }

else if(x > y)
 {
  if(x < 1.0)
   {
    t = x*x + y*y;
    d = sqrt(t);
   }
  else
   {
    d = y/x;
    t = sqrt(1.0 + d*d);
    iofl = idf_mathd_mul(x,t, &d);
   }
 }

else
 {
  if(y < 1.0)
   {
    t = x*x + y*y;
    d = sqrt(t);
   }
  else
   {
    d = x/y;
    t = sqrt(1.0 + d*d);
    iofl = idf_mathd_mul(x,t, &d);
   }
 }

if(iofl)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl=IDF_MATH_DOMAIN;
 }
else
 {
  *yr = xr;
  *yi = xi;
  iofl=0;
 } 

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_real(double xr, double xi,
                   double* yr,double* yi)
#else
int idf_mathc_real(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double xx;

 if(xr < 0.0) xx = -xr; else xx=xr;
 if(xx > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   *yr = IDF_DOUBLE_MAX;
   *yi = 0.0;
   iofl=IDF_MATH_DOMAIN;
  }
 else
  {
   *yr = xr;
   *yi = 0.0;
   iofl=0;
  }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_imag(double xr, double xi,
                   double* yr, double* yi)
#else
int idf_mathc_imag(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double xx;

 if(xi < 0.0) xx = -xi; else xx=xi;
 if(xx > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   *yr = IDF_DOUBLE_MAX;
   *yi = 0.0;
   iofl=IDF_MATH_DOMAIN;
  }
 else
  {
   *yr = xi;
   *yi = 0.0;
   iofl=0;
  }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_polar(double mod, double arg,
                    double* yr, double* yi)
#else
int idf_mathc_polar(mod,arg, yr,yi)
double mod,arg;
double *yr,*yi;
#endif
{
int iofl;
double xx,yy;

   *yr = 0.0;
   *yi = 0.0;
   iofl= 0;

 if(mod < 0.0) xx = -mod; else xx=mod;
 if(arg < 0.0) yy = -arg; else yy=arg;

 if(xx > IDF_DOUBLE_MAX || yy > IDF_DOUBLE_MAX)
  {
   /* improper argument */
   *yr = IDF_DOUBLE_MAX;
   *yi = IDF_DOUBLE_MAX;
   iofl=IDF_MATH_DOMAIN;
  }

 else
  {
      iofl = idf_mathd_cos(arg, &xx);
   if(iofl) goto err;

      iofl = idf_mathd_sin(arg, &yy);
   if(iofl) goto err;

   *yr = mod * xx;
   *yi = mod * yy;

   err:
   if(iofl)
    {
     *yr = IDF_DOUBLE_MAX;
     *yi = IDF_DOUBLE_MAX;
     iofl = IDF_MATH_RANGE;
    }
  }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_abs(double xr, double xi, double* d)
#else
int idf_mathc_abs(xr,xi, d)
double xr,xi;
double *d;
#endif
{
int iofl;
double x,y,t;

iofl = 0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *d = 0.0;
 }

else if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *d = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
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
    y = y/x;
    t = sqrt(1.0 + y*y);
    iofl = idf_mathd_mul(x,t, d);
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
    x = x/y;
    t = sqrt(1.0 + x*x);
    iofl = idf_mathd_mul(y,t, d);
   }
 }
 
return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_arg(double xr, double xi, double* a)
#else
int idf_mathc_arg(xr,xi, a)
double xr,xi;
double *a;
#endif
{
int iofl;
double x,y,t,d;

iofl = 0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(y < IDF_DOUBLE_MIN)
 {
  d = 0.0;
 }

else if(y > IDF_DOUBLE_MAX || x > IDF_DOUBLE_MAX)
 {
  d = IDF_DOUBLE_MAX;
  iofl=IDF_MATH_DOMAIN;
 }

else if(x > y)
 {
  if(x < 1.0)
   {
    t = sqrt(x*x + y*y);
   }
  else
   {
    t = y/x;
    d = sqrt(1.0 + t*t);
    iofl = idf_mathd_mul(x,d, &t);
   }
  
  if(iofl)
   {
    d = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_RANGE;
   }
  else
   {
    d = acos(xr/t);
    if(xi < 0.0) d = -d;
   }
 }

else
 {
  if(y < 1.0)
   {
    t = sqrt(x*x + y*y);
   }
  else
   {
    t = x/y;
    d = sqrt(1.0 + t*t);
    iofl = idf_mathd_mul(y,d, &t);
   }

  if(iofl)
   {
    d = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_RANGE;
   }
  else
   {
    d = acos(xr/t);
    if(xi < 0.0) d = -d;
   }
 }

*a = d;
return iofl;
}


#if IDF_CPP_ == 1
int idf_mathc_sin(double xr, double xi,
                  double* yr, double* yi)
#else
int idf_mathc_sin(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d,cx,cy;


   iofl = idf_mathd_sin(xr, &x);
if(iofl) goto err;
   iofl = idf_mathd_cosh(xi, &y);
if(iofl) goto err;
   iofl = idf_mathd_mul(x,y, &cx);
if(iofl) goto err;


   iofl = idf_mathd_cos(xr, &x);
if(iofl) goto err;
   iofl = idf_mathd_sinh(xi, &y);
if(iofl) goto err;
   iofl = idf_mathd_mul(x,y, &cy);
if(iofl) goto err;

iofl = idf_mathc_abs(cx,cy, &d);

err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr = cx;
    *yi = cy;
   }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_cos(double xr, double xi,
                  double* yr, double* yi)
#else
int idf_mathc_cos(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d,cx,cy;

   iofl = idf_mathd_cos(xr, &x);
if(iofl) goto err;
   iofl = idf_mathd_cosh(xi, &y);
if(iofl) goto err;
   iofl = idf_mathd_mul(x,y, &cx);
if(iofl) goto err;


   iofl = idf_mathd_sin(xr, &x);
if(iofl) goto err;
   iofl = idf_mathd_sinh(xi, &y);
if(iofl) goto err;
      x = -x;
   iofl = idf_mathd_mul(x,y, &cy);
if(iofl) goto err;

iofl = idf_mathc_abs(cx,cy, &d);

err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr = cx;
    *yi = cy;
   }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_tan(double xr, double xi,
                  double* yr, double* yi)
#else
int idf_mathc_tan(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d,cx,cy;


   iofl = idf_mathd_mul(xr,2.0, &x);
if(iofl) goto err;
   iofl = idf_mathd_mul(xi,2.0, &y);
if(iofl) goto err;


   iofl = idf_mathd_cos(x, &cx);
if(iofl) goto err;
   iofl = idf_mathd_cosh(y, &cy);
if(iofl) goto err;
   iofl = idf_mathd_add(cx,cy, &d);
if(iofl) goto err;


   iofl = idf_mathd_sin(x, &cx);
if(iofl) goto err;
   iofl = idf_mathd_sinh(y, &cy);
if(iofl) goto err;

   iofl = idf_mathd_div(cx,d, &x);
if(iofl) goto err;
   iofl = idf_mathd_div(cy,d, &y);
if(iofl) goto err;

iofl = idf_mathc_abs(x,y, &d);

err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr = x;
    *yi = y;
   }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_sinh(double xr, double xi,
                   double* yr, double* yi)
#else
int idf_mathc_sinh(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d,cx,cy;


   iofl = idf_mathd_sinh(xr, &x);
if(iofl) goto err;
   iofl = idf_mathd_cos(xi, &y);
if(iofl) goto err;
   iofl = idf_mathd_mul(x,y, &cx);
if(iofl) goto err;


   iofl = idf_mathd_cosh(xr, &x);
if(iofl) goto err;
   iofl = idf_mathd_sin(xi, &y);
if(iofl) goto err;
   iofl = idf_mathd_mul(x,y, &cy);
if(iofl) goto err;

iofl = idf_mathc_abs(cx,cy, &d);

err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr = cx;
    *yi = cy;
   }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_cosh(double xr, double xi,
                   double* yr, double* yi)
#else
int idf_mathc_cosh(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d,cx,cy;


   iofl = idf_mathd_cosh(xr, &x);
if(iofl) goto err;
   iofl = idf_mathd_cos(xi, &y);
if(iofl) goto err;
   iofl = idf_mathd_mul(x,y, &cx);
if(iofl) goto err;


   iofl = idf_mathd_sinh(xr, &x);
if(iofl) goto err;
   iofl = idf_mathd_sin(xi, &y);
if(iofl) goto err;
   iofl = idf_mathd_mul(x,y, &cy);
if(iofl) goto err;

iofl = idf_mathc_abs(cx,cy, &d);

err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr = cx;
    *yi = cy;
   }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_tanh(double xr, double xi,
                   double* yr, double* yi)
#else
int idf_mathc_tanh(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d,cx,cy;


   iofl = idf_mathd_mul(xr,2.0, &x);
if(iofl) goto err;
   iofl = idf_mathd_mul(xi,2.0, &y);
if(iofl) goto err;


   iofl = idf_mathd_cos(x, &cx);
if(iofl) goto err;
   iofl = idf_mathd_cosh(y, &cy);
if(iofl) goto err;
   iofl = idf_mathd_add(cx,cy, &d);
if(iofl) goto err;


   iofl = idf_mathd_sin(x, &cx);
if(iofl) goto err;
   iofl = idf_mathd_sinh(y, &cy);
if(iofl) goto err;

   iofl = idf_mathd_div(cy,d, &x);
if(iofl) goto err;
   iofl = idf_mathd_div(cx,d, &y);
if(iofl) goto err;

iofl = idf_mathc_abs(x,y, &d);

err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr = x;
    *yi = y;
   }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_log(double xr, double xi,
                  double* yr, double* yi)
#else
int idf_mathc_log(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,t,d;

iofl = 0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *yr = -IDF_DOUBLE_MAX;
  *yi = 0.0;
  iofl = IDF_MATH_DOMAIN;
 }

else if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }

else if(x > y)
 {
  if(x < 1.0)
   {
    t = x*x + y*y;
    d = sqrt(t);
   }
  else
   {
    d = y/x;
    t = sqrt(1.0 + d*d);
    iofl = idf_mathd_mul(x,t, &d);
   }


  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
       iofl = idf_mathd_log(d, &x);
    if(iofl)
     {
      *yr = IDF_DOUBLE_MAX;
      *yi = IDF_DOUBLE_MAX;
      iofl = IDF_MATH_DOMAIN;
     }
    else
     {
      y = acos(xr/d);
      if(xi < 0.0) y = -y;

      *yr = x;
      *yi = y;
     }
   }
 }

else
 {
  if(y < 1.0)
   {
    t = x*x + y*y;
    d = sqrt(t);
   }
  else
   {
    d = x/y;
    t = sqrt(1.0 + d*d);
    iofl = idf_mathd_mul(y,t, &d);
   }

  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
       iofl = idf_mathd_log(d, &x);
    if(iofl)
     {
      *yr = IDF_DOUBLE_MAX;
      *yi = IDF_DOUBLE_MAX;
      iofl = IDF_MATH_DOMAIN;
     }
    else
     {
      y = acos(xr/d);
      if(xi < 0.0) y = -y;

      *yr = x;
      *yi = y;
     }
   }
 }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_log10(double xr, double xi,
                    double* yr, double* yi)
#else
int idf_mathc_log10(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y;


   iofl = idf_mathc_log(xr,xi, &x,&y);
if(iofl)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }
else
 {
  x = x / IDF_LN10;
  y = y / IDF_LN10;

  *yr = x;
  *yi = y;
 }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_sqrt(double xr, double xi,
                   double* yr, double* yi)
#else
int idf_mathc_sqrt(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,t,d;

  iofl = 0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *yr = 0.0;
  *yi = 0.0;
 }

else if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }

else if(x > y)
 {
  if(x < 1.0)
   {
    t = x*x + y*y;
    d = sqrt(t);
   }
  else
   {
    d = y/x;
    t = sqrt(1.0 + d*d);
    iofl = idf_mathd_mul(x,t, &d);
   }

  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }

  else
   {
    t = sqrt(d);

    if(xr > 0.0)
     {
        d = (t + xr) * 0.5;
        x = sqrt(d);
      *yr = x;
        y = xi * 0.5 / x;
      *yi = y;
     }
    else
     {
         d = (t - xr) * 0.5;
         y = sqrt(d);
      if(xi < 0.0) y = -y;
      *yi = y;
        x = xi * 0.5 / y;
      *yr = x;
     }
   }
 }

else
 {
  if(x < 1.0)
   {
    t = x*x + y*y;
    d = sqrt(t);
   }
  else
   {
    d = x/y;
    t = sqrt(1.0 + d*d);
    iofl = idf_mathd_mul(y,t, &d);
   }

  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }

  else
   {
    t = sqrt(d);

    if(xr > 0.0)
     {
        d = (t + xr) * 0.5;
        x = sqrt(d);
      *yr = x;
        y = xi * 0.5 / x;
      *yi = y;
     }
    else
     {
         d = (t - xr) * 0.5;
         y = sqrt(d);
      if(xi < 0.0) y = -y;
      *yi = y;
        x = xi * 0.5 / y;
      *yr = x;
     }
   }
 }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_exp(double xr, double xi,
                  double* yr, double* yi)
#else
int idf_mathc_exp(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d;

   iofl = idf_mathd_exp(xr, &d);
if(iofl) goto err;


   iofl = idf_mathd_cos(xi, &x);
if(iofl) goto err;
   iofl = idf_mathd_sin(xi, &y);

err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr = x * d;
    *yi = y * d;
   }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_sq1mz2(double xr, double xi,
                     double* yr, double* yi)
#else
int idf_mathc_sq1mz2(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
/* sqrt(1-z^2) */
int iofl;
double x,y,t,d,cx,cy;

iofl = 0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }

else if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *yr = 1.0;
  *yi = 0.0;
  iofl = 0;
 }

else
 {
  /* 1-z^2 */
  if(x > y)
   {
    if(x < 1.0)
     {
       t = 1.0 + x*x + y*y;
      cy = (-2.0) * xr*xi / t;
      cx = (1.0 - xr*xr + xi*xi) / t;
       d = sqrt(t);
     }

    else
     {
       t = xi/xr;
       d = 1.0/x;
       d = d/x;
       d = d + t*t + 1.0;
      cx = 1.0 - 2.0/d;
      cy = (-2.0)*t/d;
       t = sqrt(d);

      iofl = idf_mathd_mul(x,t, &d);
     }
   }

  else
   {
    if(y < 1.0)
     {
       t = 1.0 + x*x + y*y;
      cy = (-2.0) * xr*xi / t;
      cx = (1.0 - xr*xr + xi*xi) / t;
       d = sqrt(t);
     }

    else
     {
       t = xr/xi;
       d = 1.0/y;
       d = d/y;
       d = d + t*t + 1.0;
      cy = (-2.0)*t / d;
      cx = 1.0 + cy*t;
       t = sqrt(d);

      iofl = idf_mathd_mul(y,t, &d);
     }
   }


  /* sqrt(1-z^2) */

  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }

  else
   {
       iofl = idf_mathc_sqrt(cx,cy, &x,&y);
    if(iofl) goto err;

       iofl = idf_mathd_mul(x,d, &cx);
    if(iofl) goto err;
       iofl = idf_mathd_mul(y,d, &cy);
    if(iofl) goto err;

    iofl = idf_mathc_abs(cx,cy, &d);

    err:
    if(iofl)
     {
      *yr = IDF_DOUBLE_MAX;
      *yi = IDF_DOUBLE_MAX;
      iofl = IDF_MATH_DOMAIN;
     }
    else
     {
      *yr = cx;
      *yi = cy;
     }

   }
 }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_sq1pz2(double xr, double xi,
                     double* yr, double* yi)
#else
int idf_mathc_sq1pz2(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
/* sqrt(1+z^2) */
int iofl;
double x,y,t,d,cx,cy;

iofl = 0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }

else if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *yr = 1.0;
  *yi = 0.0;
  iofl = 0;
 }

else
 {
  /* 1+z^2 */
  if(x > y)
   {
    if(x < 1.0)
     {
      t = 1.0 + x*x + y*y;
      cy = 2.0 * xr*xi / t;
      cx = (1.0 + xr*xr - xi*xi) / t;
      d = sqrt(t);
     }

    else
     {
      t = xi/xr;
      d = 1.0/x;
      d = d/x;
      d = d + t*t + 1.0;
      cy = 2.0*t/d;
      cx = 1.0 - cy*t;
       t = sqrt(d);

      iofl = idf_mathd_mul(x,t, &d);
     }
   }

  else
   {
    if(y < 1.0)
     {
      t = 1.0 + x*x + y*y;
      cy = 2.0 * xr*xi / t;
      cx = (1.0 + xr*xr - xi*xi) / t;
      d = sqrt(t);
     }

    else
     {
      t = xr/xi;
      d = 1.0/y;
      d = d/y;
      d = d + t*t + 1.0;
      cy = 2.0*t/d;
      cx = (d - 2.0)/d;
       t = sqrt(d);

      iofl = idf_mathd_mul(y,t, &d);
     }
   }

  /* sqrt(1+z^2) */
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }

  else
   {
       iofl = idf_mathc_sqrt(cx,cy, &x,&y);
    if(iofl) goto err;

       iofl = idf_mathd_mul(x,d, &cx);
    if(iofl) goto err;
       iofl = idf_mathd_mul(y,d, &cy);
    if(iofl) goto err;

    iofl = idf_mathc_abs(cx,cy, &d);

    err:
    if(iofl)
     {
      *yr = IDF_DOUBLE_MAX;
      *yi = IDF_DOUBLE_MAX;
      iofl = IDF_MATH_DOMAIN;
     }
    else
     {
      *yr = cx;
      *yi = cy;
     }

   }
 }


return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_acos(double xr, double xi,
                   double* yr, double* yi)
#else
int idf_mathc_acos(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d,cx,cy;

/* -i log(z + i sqrt(1-z^2) ] */

iofl=0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }

else if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *yr = IDF_PIOVER2;
  *yi = 0.0;
  iofl = 0;
 }

else
 {
  /* w = sqrt(1-z^2) */
     iofl = idf_mathc_sq1mz2(xr,xi, &x,&y);
  if(iofl) goto err;

  /* c = z + iw */
     iofl = idf_mathd_sub(xr,y, &cx);
  if(iofl) goto err;
     iofl = idf_mathd_add(xi,x, &cy);
  if(iofl) goto err;

  /* -i log(c) */
     iofl = idf_mathc_log(cx,cy, &x,&y);
  if(iofl) goto err;

     iofl = idf_mathc_abs(x,y, &d);

  err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr =  y;
    *yi = -x;
   }
 }


return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_asin(double xr, double xi,
                   double* yr, double* yi)
#else
int idf_mathc_asin(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d,cx,cy;

/* -i log(iz + sqrt(1-z^2) ] */

iofl=0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }

else if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *yr = 0.0;
  *yi = 0.0;
  iofl = 0;
 }

else
 {
  /* w = sqrt(1-z^2) */
     iofl = idf_mathc_sq1mz2(xr,xi, &x,&y);
  if(iofl) goto err;

  /* c = iz + w */
     iofl = idf_mathd_sub(x,xi, &cx);
  if(iofl) goto err;
     iofl = idf_mathd_add(xr,y, &cy);
  if(iofl) goto err;

  /* -i log(c) */
     iofl = idf_mathc_log(cx,cy, &x,&y);
  if(iofl) goto err;

     iofl = idf_mathc_abs(x,y, &d);

  err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr =  y;
    *yi = -x;
   }
 }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_atan(double xr, double xi,
                   double* yr, double* yi)
#else
int idf_mathc_atan(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d,t,cx,cy;

/* -i/2  log[(1+iz)/(1-iz)] */

iofl=0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }

else if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *yr = 0.0;
  *yi = 0.0;
  iofl = 0;
 }

else
 {
  /* w = (1+iz)/(1-iz) = {(1-x^2-y^2) + i 2x}/ (x^2+(1+y)^2) */

     d = xi + 1.0;
  if(d < 0.0) d=-d;
  if(x < IDF_DOUBLE_MIN && d < IDF_DOUBLE_MIN)
   {
    /* pole: x=0 y=-1;*/
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl=IDF_MATH_DOMAIN;
   }

  else if(x > y)
   {
    if(x < 1.0)
     {
       d = xi + 1.0;
       t = x*x + d*d;
      cx = 2.0*d/t - 1.0;
      cy = 2.0 * xr / t;
     }
    else
     {
       t = (xi + 1.0)/xr;
       d = 1.0 + t*t;
      cy = (2.0/xr)/d;
      cx = cy*t - 1.0;
     }
   }

  else
   {
    if(y < 1.0)
     {
       d = xi + 1.0;
       t = x*x + d*d;
      cx = 2.0*d/t - 1.0;
      cy = 2.0 * xr / t;
     }
    else
     {
       t = (xi + 1.0)/y;
       x = xr/y;
       d = x*x + t*t;
      cy = (2.0*x/y)/d;
      cx = (2.0*t/y)/d - 1.0;
     }
   }

  /* log(w) */
  if(!iofl)
   {
       iofl = idf_mathc_log(cx,cy, &x,&y);
    if(iofl)
     {
      *yr = IDF_DOUBLE_MAX;
      *yi = IDF_DOUBLE_MAX;
      iofl=IDF_MATH_DOMAIN;
     }
    else
     {
      /* i/2 log(w) */
      *yr =  y*0.5;
      *yi = -x*0.5;
      iofl = 0;
     }
   }
 }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_acosh(double xr, double xi,
                    double* yr, double* yi)
#else
int idf_mathc_acosh(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,cx,cy;

/* i acos(z) = log[z + i sqrt(1-z^2)] */

iofl=0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }

else if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *yr = 0.0;
  *yi = IDF_PIOVER2;
  iofl = 0;
 }

else
 {
  /* w = sqrt(1-z^2) */
     iofl = idf_mathc_sq1mz2(xr,xi, &x,&y);
  if(iofl) goto err;

  /* c = z + iw */
     iofl = idf_mathd_sub(xr,y, &cx);
  if(iofl) goto err;
     iofl = idf_mathd_add(xi,x, &cy);
  if(iofl) goto err;

  /* log(c) */
     iofl = idf_mathc_log(cx,cy, &x,&y);

  err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr = x;
    *yi = y;
   }

 }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_asinh(double xr, double xi,
                    double* yr, double* yi)
#else
int idf_mathc_asinh(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,cx,cy;

/* -i asin(iz) = log[z + sqrt(1+z^2)] */

iofl=0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }

else if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *yr = 0.0;
  *yi = 0.0;
  iofl = 0;
 }

else
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;

  /* w = sqrt(1+z^2) */
     iofl = idf_mathc_sq1pz2(xr,xi, &x,&y);
  if(iofl) goto err;

  /* c = z + w */
     iofl = idf_mathd_add(x,xr, &cx);
  if(iofl) goto err;
     iofl = idf_mathd_add(y,xi, &cy);
  if(iofl) goto err;

  /* log(c) */
     iofl = idf_mathc_log(cx,cy, &x,&y);

  err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr = x;
    *yi = y;
   }
 }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_atanh(double xr, double xi,
                    double* yr, double* yi)
#else
int idf_mathc_atanh(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d,t,cx,cy;

/* 1/2  log[(1+z)/(1-z)] */

iofl = 0;

if(xr < 0.0) x=-xr; else x=xr;
if(xi < 0.0) y=-xi; else y=xi;

if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *yr = IDF_DOUBLE_MAX;
  *yi = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }

else if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *yr = 0.0;
  *yi = 0.0;
  iofl = 0;
 }

else
 {
  /* w = (1+z)/(1-z) = {(1-x^2-y^2) + i 2y}/ ((x-1)^2+y^2) */

     d = xr - 1.0;
  if(d < 0.0) d=-d;
  if(y < IDF_DOUBLE_MIN && d < IDF_DOUBLE_MIN)
   {
    /* pole: x=1 y=0;*/
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl=IDF_MATH_DOMAIN;
   }

  else if(x > y)
   {
    if(x < 1.0)
     {
       d = 1.0 - xr;
       t = y*y + d*d;
      cx = 2.0*d/t - 1.0;
      cy = 2.0 * xi / t;
     }
    else
     {
       t = (1.0 - xr)/xr;
       y = xi/xr;
       d = y*y + t*t;
      cx = (2.0*t/xr)/d - 1.0;
      cy = (2.0*y/xr)/d;
     }
   }

  else
   {
    if(y < 1.0)
     {
       d = 1.0 - xr;
       t = y*y + d*d;
      cx = 2.0*d/t - 1.0;
      cy = 2.0 * xi / t;
     }
    else
     {
       t = (1.0 - xr)/xi;
       d = 1.0 + t*t;
      cy = (2.0/xi)/d;
      cx = cy*t - 1.0;
     }
   }

  /* log(w) */
  if(!iofl)
   {
       iofl = idf_mathc_log(cx,cy, &x,&y);
    if(iofl)
     {
      *yr = IDF_DOUBLE_MAX;
      *yi = IDF_DOUBLE_MAX;
      iofl=IDF_MATH_DOMAIN;
     }
    else
     {
      *yr = x*0.5;
      *yi = y*0.5;
     }
   }
 }

return iofl;
}


#if IDF_CPP_ == 1
int idf_mathc_min(double xr, double xi,
                  double yr, double yi,
                  double* zr, double* zi)
#else
int idf_mathc_min(xr,xi, yr,yi, zr,zi)
double xr,xi;
double yr,yi;
double *zr,*zi;
#endif
{
int iofl,ido;
double xd,yd,d;

 
   iofl = idf_mathc_abs(xr,xi, &xd);
if(iofl) goto err;
   iofl = idf_mathc_abs(yr,yi, &yd);
if(iofl) goto err;

   d = xd-yd;
if(d < 0.0) d = -d;
if(d < IDF_DOUBLE_MIN)
 {
     iofl = idf_mathc_arg(xr,xi, &xd);
  if(iofl) goto err;
     iofl = idf_mathc_arg(yr,yi, &yd);
  if(iofl) goto err;

  if(xd > yd)
   {
    ido=0;
   }
  else
   {
    ido=1;
   }
 }

else if(xd > yd)
 {
  ido=0;
 }
else
 {
  ido=1;
 }

if(ido)
 {
  *zr = xr;
  *zi = xi;
 }
else
 {
  *zr = yr;
  *zi = yi;
 }

err:
  if(iofl)
   {
    *zr = IDF_DOUBLE_MAX;
    *zi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_max(double xr, double xi,
                  double yr, double yi,
                  double* zr, double* zi)
#else
int idf_mathc_max(xr,xi, yr,yi, zr,zi)
double xr,xi;
double yr,yi;
double *zr,*zi;
#endif
{
int iofl,ido;
double xd,yd,d;

 
   iofl = idf_mathc_abs(xr,xi, &xd);
if(iofl) goto err;
   iofl = idf_mathc_abs(yr,yi, &yd);
if(iofl) goto err;

   d = xd-yd;
if(d < 0.0) d = -d;
if(d < IDF_DOUBLE_MIN)
 {
     iofl = idf_mathc_arg(xr,xi, &xd);
  if(iofl) goto err;
     iofl = idf_mathc_arg(yr,yi, &yd);
  if(iofl) goto err;

  if(xd < yd)
   {
    ido=0;
   }
  else
   {
    ido=1;
   }
 }

else if(xd < yd)
 {
  ido=0;
 }
else
 {
  ido=1;
 }

if(ido)
 {
  *zr = xr;
  *zi = xi;
 }
else
 {
  *zr = yr;
  *zi = yi;
 }

err:
  if(iofl)
   {
    *zr = IDF_DOUBLE_MAX;
    *zi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_RANGE;
   }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_pow10(double xr, double xi,
                    double* yr, double* yi)
#else
int idf_mathc_pow10(xr,xi, yr,yi)
double xr,xi;
double *yr,*yi;
#endif
{
int iofl;
double x,y,d,z;

   iofl = idf_mathd_pow10(xr, &d);
if(iofl) goto err;

if(xi < 0.0) z=-xi; else z=xi;
if(z < 1.0)
 {
  z = xi*IDF_LN10;
 }
else
 {
     iofl = idf_mathd_mul(xi,IDF_LN10, &z);
  if(iofl) goto err;
 }

     iofl = idf_mathd_cos(z, &x);
  if(iofl) goto err;
     iofl = idf_mathd_sin(z, &y);


err:
  if(iofl)
   {
    *yr = IDF_DOUBLE_MAX;
    *yi = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_DOMAIN;
   }
  else
   {
    *yr = x * d;
    *yi = y * d;
   }

return iofl;
}

