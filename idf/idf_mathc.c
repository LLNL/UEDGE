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
int idf_mathc_add(double ar, double ai,
                  double br, double bi,
                  double* cr, double* ci)
#else
int idf_mathc_add(ar,ai,br,bi, cr,ci)
double ar,ai,br,bi;
double *cr,*ci;
#endif
{
int iofl;

*cr=0.0;
*ci=0.0;

    iofl = idf_mathd_add(ar,br, cr);
if(!iofl)
 {
  iofl = idf_mathd_add(ai,bi, ci);
 }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_sub(double ar, double ai,
                  double br, double bi,
                  double* cr, double* ci)
#else
int idf_mathc_sub(ar,ai,br,bi, cr,ci)
double ar,ai,br,bi;
double *cr,*ci;
#endif
{
int iofl;
double d;

*cr=0.0;
*ci=0.0;

       d = -br;
    iofl = idf_mathd_add(ar,d, cr);
if(!iofl)
 {
     d = -bi;
  iofl = idf_mathd_add(ai,d, ci);
 }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_mul(double ar, double ai,
                  double br, double bi,
                  double* cr, double* ci)
#else
int idf_mathc_mul(ar,ai,br,bi, cr,ci)
double ar,ai,br,bi;
double *cr,*ci;
#endif
{
int iofl;
double x,y,d,xx,yy,dd,cx,cy;

if(ar < 0.0) x=-ar; else x=ar;
if(ai < 0.0) y=-ai; else y=ai;

if(br < 0.0) xx=-br; else xx=br;
if(bi < 0.0) yy=-bi; else yy=bi;

if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *cr = 0.0;
  *ci = 0.0;
 iofl = 0;
 }

else if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *cr = IDF_DOUBLE_MAX;
  *ci = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_OVERFLOW;
 }

if(xx < IDF_DOUBLE_MIN && yy < IDF_DOUBLE_MIN)
 {
  *cr = 0.0;
  *ci = 0.0;
 iofl = 0;
 }

else if(xx > IDF_DOUBLE_MAX || yy > IDF_DOUBLE_MAX)
 {
  *cr = IDF_DOUBLE_MAX;
  *ci = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_OVERFLOW;
 }

else
 {
  iofl = 0;

  if(x > y)
   {
     d = x;
    ar = ar/x;
    ai = ai/x;
   }
  else
   {
     d = y;
    ar = ar/y;
    ai = ai/y;
   }

  if(xx > yy)
   {
    dd = xx;
    br = br/xx;
    bi = bi/xx;
   }
  else
   {
    dd = yy;
    br = br/yy;
    bi = bi/yy;
   }

  cx = ar*br - ai*bi;
  cy = ar*bi + ai*br;

  if(d<1.0 && dd>1.0)
   {
    d = d*dd;
    x = cx * d;
    y = cy * d;
   }

  else if(d>1.0 && dd<1.0)
   {
    d = d*dd;
    x = cx * d;
    y = cy * d;
   }

  else if(d > dd)
   {
    cx = cx * d;
    cy = cy * d;

        iofl = idf_mathd_mul(cx,dd, &x);
    if(!iofl)
      {
       iofl = idf_mathd_mul(cy,dd, &y);
      }
   }
  else
   {
    cx = cx * dd;
    cy = cy * dd;

        iofl = idf_mathd_mul(cx,d, &x);
    if(!iofl)
      {
       iofl = idf_mathd_mul(cy,d, &y);
      }
   }

  if(iofl)
   {
    *cr = IDF_DOUBLE_MAX;
    *ci = IDF_DOUBLE_MAX;
    iofl = IDF_MATH_OVERFLOW;
   }
  else
   {
    *cr = x;
    *ci = y;
   }
 }

return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_div(double ar, double ai,
                  double br, double bi,
                  double* cr, double* ci)
#else
int idf_mathc_div(ar,ai,br,bi, cr,ci)
double ar,ai,br,bi;
double *cr,*ci;
#endif
{
int iofl;
double x,y,d,xx,yy,dd,cx,cy;

if(br < 0.0) x=-br; else x=br;
if(bi < 0.0) y=-bi; else y=bi;

if(ar < 0.0) xx=-ar; else xx=ar;
if(ai < 0.0) yy=-ai; else yy=ai;

if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
 {
  *cr = IDF_DOUBLE_MAX;
  *ci = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_EXCEPTION;
 }

else if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *cr = IDF_DOUBLE_MAX;
  *ci = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_OVERFLOW;
 }

else if(xx > IDF_DOUBLE_MAX || yy > IDF_DOUBLE_MAX)
 {
  *cr = IDF_DOUBLE_MAX;
  *ci = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_OVERFLOW;
 }

else
 {
  if(x > y)
   {
     d = bi/br;
    dd = 1.0 + d*d;

    x = (br/x)/dd;
    y = d/dd;
    d = x;
   }
  else
   {
     d = br/bi;
    dd = 1.0 + d*d;

    x = d/dd;
    y = (bi/y)/dd;
    d = y;
   }

  if(xx > yy)
   {
    dd = xx;
   }
  else
   {
    dd = yy;
   }

  xx = ar/dd;
  yy = ai/dd;

  cx = xx*x + yy*y;
  cy = xx*y - yy*x;

  if(dd>1.0 && d<1.0)
   {
    cx = cx * dd;
    cy = cy * dd;

        iofl = idf_mathd_div(cx,d, &x);
    if(!iofl)
     {
      iofl = idf_mathd_div(cy,d, &y);
     }

    if(iofl)
     {
      *cr = IDF_DOUBLE_MAX;
      *ci = IDF_DOUBLE_MAX;
      iofl = IDF_MATH_OVERFLOW;
     }
    else
     {
      *cr = x;
      *ci = y;
      iofl = 0;
     }
   }

  else
   {
    dd = dd/d;
    *cr = cx * dd;
    *ci = cy * dd;
    iofl = 0;
   }

 }


return iofl;
}



#if IDF_CPP_ == 1
int idf_mathc_pow(double ar, double ai,
                  double br, double bi,
                  double* cr, double* ci)
#else
int idf_mathc_pow(ar,ai,br,bi, cr,ci)
double ar,ai,br,bi;
double *cr,*ci;
#endif
{
int iofl;
double x,y,xx,yy,t,tt,d,dd;

iofl = 0;

if(ar < 0.0) x=-ar; else x=ar;
if(ai < 0.0) y=-ai; else y=ai;

if(br < 0.0) xx=-br; else xx=br;
if(bi < 0.0) yy=-bi; else yy=bi;


if(x > IDF_DOUBLE_MAX || y > IDF_DOUBLE_MAX)
 {
  *cr = IDF_DOUBLE_MAX;
  *ci = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }

else if(xx > IDF_DOUBLE_MAX || yy > IDF_DOUBLE_MAX)
 {
  *cr = IDF_DOUBLE_MAX;
  *ci = IDF_DOUBLE_MAX;
  iofl = IDF_MATH_DOMAIN;
 }


else
 {
  if(x < IDF_DOUBLE_MIN && y < IDF_DOUBLE_MIN)
   {
    d = 0.0;
    t = 0.0;
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

    if(!iofl)
     {
      if(d > IDF_DOUBLE_MIN)
       {
        t = acos(ar/d);
        if(ai < 0.0) t = -t;
       }
      else
       {
        t = 0.0;
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

    if(!iofl)
     {
      if(d > IDF_DOUBLE_MIN)
       {
        t = acos(ar/d);
        if(ai < 0.0) t = -t;
       }
      else
       {
        t = 0.0;
       }
     }
   }


  if(iofl)
   {
    dd = 0.0;
    tt = 0.0;
   }

  else if(xx < IDF_DOUBLE_MIN && yy < IDF_DOUBLE_MIN)
   {
    dd = 0.0;
    tt = 0.0;
   }

  else if(xx > yy)
   {
    if(xx < 1.0)
     {
      tt = xx*xx + yy*yy;
      dd = sqrt(tt);
     }
    else
     {
      dd = yy/xx;
      tt = sqrt(1.0 + dd*dd);
      iofl = idf_mathd_mul(xx,tt, &dd);
     }

    if(!iofl)
     {
      if(dd > IDF_DOUBLE_MIN)
       {
        tt = acos(br/dd);
        if(bi < 0.0) tt = -tt;
       }
      else
       {
        tt = 0.0;
       }
     }
   }

  else
   {
    if(yy < 1.0)
     {
      tt = xx*xx + yy*yy;
      dd = sqrt(tt);
     }
    else
     {
      dd = xx/yy;
      tt = sqrt(1.0 + dd*dd);
      iofl = idf_mathd_mul(yy,tt, &dd);
     }

    if(!iofl)
     {
      if(dd > IDF_DOUBLE_MIN)
       {
        tt = acos(br/dd);
        if(bi < 0.0) tt = -tt;
       }
      else
       {
        tt = 0.0;
       }
     }
   }

  if(iofl)
   {
    *cr = IDF_DOUBLE_MAX;
    *ci = IDF_DOUBLE_MAX;
    iofl=IDF_MATH_DOMAIN;
   }

  else if(d < IDF_DOUBLE_MIN)
   {
    if(dd < IDF_DOUBLE_MIN)
     {
      *cr = 0.0;
      *ci = 0.0;
      iofl = 0;
     }
    else if(br < 0.0)
     {
      *cr = IDF_DOUBLE_MAX;
      *ci = IDF_DOUBLE_MAX;
      iofl=IDF_MATH_DOMAIN;
     }
    else
     {
      *cr = 0.0;
      *ci = 0.0;
      iofl=0;
     }
   }

  else if(dd < IDF_DOUBLE_MIN)
   {
    *cr = 1.0;
    *ci = 0.0;
    iofl=0;    
   }

  else
   {
       iofl = idf_mathd_log(d, &tt);
    if(iofl) goto err;

    if(tt < 0.0) x=-tt; else x=tt;
    if(t  < 0.0) y=-t;  else y=t;
    if(x > y) d = x; else d = y;

     x = tt/d;
     y = t/d;
    xx = br/dd;
    yy = bi/dd;

     t = x*xx - y*yy;
    tt = y*xx + yy*x;

    if(d<1.0 && dd>1.0)
     {
      d = d*dd;
      x = t*d;
      y = tt*d;
     }
    else if(d>1.0 && d<1.0)
     {
      d = d*dd;
      x = t*d;
      y = tt*d;
     }
    else if(d > dd)
     {
       t = t*d;
      tt = tt*d;

         iofl = idf_mathd_mul(t,dd, &x);
     if(!iofl)
      {
       iofl = idf_mathd_mul(tt,dd, &y);
      }
     }
    else
     {
       t = t*dd;
      tt = tt*dd;

         iofl = idf_mathd_mul(t,d, &x);
     if(!iofl)
      {
       iofl = idf_mathd_mul(tt,d, &y);
      }
     }

    if(iofl)
     {
      *cr = IDF_DOUBLE_MAX;
      *ci = IDF_DOUBLE_MAX;
      iofl=IDF_MATH_DOMAIN;
     }
    else
     {
      /* exp(x + iy) */

         iofl = idf_mathd_exp(x, &d);
      if(iofl) goto err;

         iofl = idf_mathd_cos(y, &xx);
      if(iofl) goto err;
         iofl = idf_mathd_sin(y, &yy);

      err:
      if(iofl)
       {
        *cr = IDF_DOUBLE_MAX;
        *ci = IDF_DOUBLE_MAX;
        iofl=IDF_MATH_DOMAIN;
       }
      else
       {
        *cr = xx * d;
        *ci = yy * d;
       }
     }
   }
   
 }

return iofl;
}
