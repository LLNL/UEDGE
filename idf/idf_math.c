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



static int idf_math_sub=0;
static int idf_math_err=0;



void idf_matherr_nul()
{
 idf_math_sub = 0;
 idf_math_err = 0;
}


#if IDF_CPP_ == 1
void idf_matherr_put(int sub, int err)
#else
void idf_matherr_put(sub,err)
int sub;
int err;
#endif
{
 idf_math_sub = sub;
 idf_math_err = err;
}



#if IDF_CPP_ == 1
int idf_mathd(char typ, double a, double b, double* c)
#else
int idf_mathd(typ, a, b, c)
char    typ;
double  a;
double  b;
double *c;
#endif
{
int flag;

 if(typ == '+')
  {
   flag = idf_mathd_add(a,b, c);
  }
 else if(typ == '-')
  {
   flag = idf_mathd_sub(a,b, c);
  }
 else if(typ == '*')
  {
   flag = idf_mathd_mul(a,b, c);
  }
 else if(typ == '/')
  {
   flag = idf_mathd_div(a,b, c);
  }
 else if(typ == '^')
  {
   flag = idf_mathd_pow(a,b, c);
  }
 else
  {
   flag = IDF_MATH_OPERAND;
  }

return flag;
}



#if IDF_CPP_ == 1
int idf_mathc(char typ,
              double ar, double ai,
              double br, double bi,
              double* cr, double* ci)
#else
int idf_mathc(typ, ar,ai, br,bi, cr,ci)
char    typ;
double  ar,ai;
double  br,bi;
double *cr,*ci;
#endif
{
int flag;

 if(typ == '+')
  {
   flag = idf_mathc_add(ar,ai,br,bi, cr,ci);
  }
 else if(typ == '-')
  {
   flag = idf_mathc_sub(ar,ai,br,bi, cr,ci);
  }
 else if(typ == '*')
  {
   flag = idf_mathc_mul(ar,ai,br,bi, cr,ci);
  }
 else if(typ == '/')
  {
   flag = idf_mathc_div(ar,ai,br,bi, cr,ci);
  }
 else if(typ == '^')
  {
   flag = idf_mathc_pow(ar,ai,br,bi, cr,ci);
  }
 else
  {
   flag = IDF_MATH_OPERAND;
  }

return flag;
}


#if IDF_CPP_ == 1
int idf_math(char math,
             int typ1, double vr1, double vi1,
             int typ2, double vr2, double vi2,
             int* jtyp, double* vr, double* vi)
#else
int idf_math(math, typ1,vr1,vi1, typ2,vr2,vi2,
                   jtyp,vr,vi)
char    math;
int     typ1;
double  vr1,vi1;
int     typ2;
double  vr2,vi2;
int    *jtyp;
double *vr,*vi;
#endif
{
/* 
  0 - OK
 >0 - see matherr
*/
int ierr=0;
int typ;

if(typ1 < typ2) typ=typ2; 
else            typ=typ1;

if(typ == IDF_FRMT_W)
 {
  /* complex */
  ierr = idf_mathc(math, vr1,vi1, vr2,vi2, vr,vi);
 }

else if(typ == IDF_FRMT_D)
 {
  /* real */
  ierr = idf_mathd(math, vr1, vr2, vr);
  *vi = 0.0;
 }

else
 {
  ierr=IDF_MATH_OPERARG;
 }

if(ierr)
 {
  *vr=0.0;
  *vi=0.0;

  idf_math_sub = 0;
  idf_math_err = ierr;
 }

  *jtyp=typ;

 return ierr;
}



#if IDF_CPP_ == 1
int idf_mathdf(int jtype, double x1, double x2,
               double* y, double* yy)
#else
int idf_mathdf(jtype, x1,x2, y,yy)
int     jtype;
double  x1,x2;
double *y,*yy;
#endif
{
int iofl=0;

*yy = 0.0;

if(jtype==1)
 {
  iofl = idf_mathd_abs(x2,y);
 }
else if(jtype==2)
 {
  iofl = idf_mathd_int(x2,y);
 }
else if(jtype==3)
 {
  iofl = idf_mathd_sqrt(x2,y);
 }
else if(jtype==4)
 {
  iofl = idf_mathd_log(x2,y);
 }
else if(jtype==5)
 {
  iofl = idf_mathd_log10(x2,y);
 }
else if(jtype==6)
 {
  iofl = idf_mathd_exp(x2,y);
 }
else if(jtype==7)
 {
  iofl = idf_mathd_cos(x2,y);
 }
else if(jtype==8)
 {
  iofl = idf_mathd_sin(x2,y);
 }
else if(jtype==9)
 {
  iofl = idf_mathd_tan(x2,y);
 }
else if(jtype==10)
 {
  iofl = idf_mathd_acos(x2,y);
 }
else if(jtype==11)
 {
  iofl = idf_mathd_asin(x2,y);
 }
else if(jtype==12)
 {
  iofl = idf_mathd_atan(x2,y);
 }
else if(jtype==13)
 {
  iofl = idf_mathd_pow(x1,x2,y);
 }
else if(jtype==14)
 {
  iofl = idf_mathd_sinh(x2,y);
 }
else if(jtype==15)
 {
  iofl = idf_mathd_cosh(x2,y);
 }
else if(jtype==16)
 {
  iofl = idf_mathd_tanh(x2,y);
 }
else if(jtype==17)
 {
  iofl = idf_mathd_asinh(x2,y);
 }
else if(jtype==18)
 {
  iofl = idf_mathd_acosh(x2,y);
 }
else if(jtype==19)
 {
  iofl = idf_mathd_atanh(x2,y);
 }
else if(jtype==20)
 {
  iofl = idf_mathc_cmplx(x1,x2, y,yy);
 }
else if(jtype==21)
 {
  iofl = idf_mathc_cmplx(0.0,x2, y,yy);
 }
else if(jtype==22)
 {
  *y  = x2;
  *yy = 0.0;
 }
else if(jtype==23)
 {
  *y  = 0.0;
  *yy = 0.0;
 }
else if(jtype==24)
 {
  iofl = idf_mathc_polar(x1,x2, y,yy);
 }
else if(jtype==26)
 {
  iofl = idf_mathd_min(x1,x2,y);
 }
else if(jtype==27)
 {
  iofl = idf_mathd_max(x1,x2,y);
 }
else if(jtype==28)
 {
  iofl = idf_mathd_hypot(x1,x2,y);
 }
else if(jtype==29)
 {
  *y  = x1;
  *yy = 0.0;
 }
else if(jtype==31)
 {
  iofl = idf_mathd_sq1pz2(x2,y);
 }
else if(jtype==32)
 {
  iofl = idf_mathd_sq1mz2(x2,y);
 }
else if(jtype==33)
 {
  iofl = idf_mathd_pow10(x2,y);
 }
else
 {
  iofl=IDF_MATH_DFUN;
 }

return iofl;
} 



#if IDF_CPP_ == 1
int idf_mathcf(int jtype,
               double xr1, double xi1,
               double xr2, double xi2,
               double* yr, double* yi)
#else
int idf_mathcf(jtype, xr1,xi1, xr2,xi2, yr,yi)
int     jtype;
double  xr1,xi1;
double  xr2,xi2;
double *yr,*yi;
#endif
{
int iofl=0;

if(jtype==1)
 {
  iofl = idf_mathc_abs(xr2,xi2, yr);
  *yi = 0.0;
 }
else if(jtype==3)
 {
  iofl = idf_mathc_sqrt(xr2,xi2, yr,yi);
 }
else if(jtype==4)
 {
  iofl = idf_mathc_log(xr2,xi2, yr,yi);
 }
else if(jtype==5)
 {
  iofl = idf_mathc_log10(xr2,xi2, yr,yi);
 }
else if(jtype==6)
 {
  iofl = idf_mathc_exp(xr2,xi2, yr,yi);
 }
else if(jtype==7)
 {
  iofl = idf_mathc_cos(xr2,xi2, yr,yi);
 }
else if(jtype==8)
 {
  iofl = idf_mathc_sin(xr2,xi2, yr,yi);
 }
else if(jtype==9)
 {
  iofl = idf_mathc_tan(xr2,xi2, yr,yi);
 }
else if(jtype==10)
 {
  iofl = idf_mathc_acos(xr2,xi2, yr,yi);
 }
else if(jtype==11)
 {
  iofl = idf_mathc_asin(xr2,xi2, yr,yi);
 }
else if(jtype==12)
 {
  iofl = idf_mathc_atan(xr2,xi2, yr,yi);
 }
else if(jtype==13)
 {
  iofl = idf_mathc_pow(xr1,xi1,xr2,xi2, yr,yi);
 }
else if(jtype==14)
 {
  iofl = idf_mathc_sinh(xr2,xi2, yr,yi);
 }
else if(jtype==15)
 {
  iofl = idf_mathc_cosh(xr2,xi2, yr,yi);
 }
else if(jtype==16)
 {
  iofl = idf_mathc_tanh(xr2,xi2, yr,yi);
 }
else if(jtype==17)
 {
  iofl = idf_mathc_asinh(xr2,xi2, yr,yi);
 }
else if(jtype==18)
 {
  iofl = idf_mathc_acosh(xr2,xi2, yr,yi);
 }
else if(jtype==19)
 {
  iofl = idf_mathc_atanh(xr2,xi2, yr,yi);
 }
else if(jtype==20)
 {
  iofl = idf_mathc_cmplx(xr1,xr2, yr,yi);
 }
else if(jtype==21)
 {
  *yr = xr2;
  *yi = 0.0;
  iofl=0;
 }
else if(jtype==22)
 {
  iofl = idf_mathc_real(xr2,xi2, yr,yi);
 }
else if(jtype==23)
 {
  iofl = idf_mathc_imag(xr2,xi2, yr,yi);
 }
else if(jtype==24)
 {
  iofl = idf_mathc_polar(xr1,xr2, yr,yi);
 }
else if(jtype==25)
 {
  iofl = idf_mathc_arg(xr2,xi2, yr);
  *yi = 0.0;
 }
else if(jtype==26)
 {
  iofl = idf_mathc_min(xr1,xi1,xr2,xi2, yr,yi);
 }
else if(jtype==27)
 {
  iofl = idf_mathc_max(xr1,xi1,xr2,xi2, yr,yi);
 }
else if(jtype==30)
 {
  *yr =  xr2;
  *yi = -xi2;
  iofl=0;
 }
else if(jtype==31)
 {
  iofl = idf_mathc_sq1pz2(xr2,xi2, yr,yi);
 }
else if(jtype==32)
 {
  iofl = idf_mathc_sq1mz2(xr2,xi2, yr,yi);
 }
else if(jtype==33)
 {
  iofl = idf_mathc_pow10(xr2,xi2, yr,yi);
 }
else
 {
  iofl=IDF_MATH_CFUN;
 }

return iofl;
} 



#if IDF_CPP_ == 1
int idf_math_func(int ftype, int narg,
                  int typ1, double vr1, double vi1,
                  int typ2, double vr2, double vi2,
                  int* typ, double* vr, double* vi)
#else
int idf_math_func(ftype,narg,
                  typ1,vr1,vi1, typ2,vr2,vi2,
                  typ,vr,vi)
int     ftype;
int     narg;
int     typ1;
double  vr1,vi1;
int     typ2;
double  vr2,vi2;
int    *typ;
double *vr,*vi;
#endif
{
/* 
  0 - OK
 >0 - see matherr
*/
int ierr=0;
int iofl,jarg,targ,jret,atyp,jtyp,jrt;

jtyp=0;
jret=0;
 jrt=0;

 if(ftype < 1)
  {
   ierr=IDF_MATH_FUN;
   goto err;
  }

 /* number of arguments */
 jarg = idf_keywF_narg(ftype);

 /* argument type: 0-void D-double W-complex*/
 targ = idf_keywF_targ(ftype);

 /* return type: 0-void D-double W-complex*/
  jrt = idf_keywF_tfun(ftype);

 /* real type of passed arguments */
 if(typ1 < typ2) atyp=typ2; 
 else            atyp=typ1;

 if(targ)
  {
   if(targ != atyp)
    {
     /* incompatible arguments */
     ierr=IDF_MATH_FARG;
    }
   else
    {
     jtyp=targ;
     if(jarg != narg)
      {
       /* improper number of arguments*/
       ierr=IDF_MATH_NARG;
      }
    }
  }
 else
  {
   jtyp = atyp;
   if(jarg != narg)
    {
     /* improper number of arguments*/
     ierr=IDF_MATH_NARG;
    }
  }

 if(jtyp == IDF_FRMT_W)
  {
   /* complex */
      iofl = idf_mathcf(ftype, vr1,vi1, vr2,vi2, vr,vi);
   if(iofl)
    {
     ierr=iofl;
    }
   else
    {
     if(jrt)
      {
       jret = jrt;
      }
     else
      {
       jret = jtyp;
      }
    }
  }
 else if(jtyp == IDF_FRMT_D)
  {
   /* real */
      iofl = idf_mathdf(ftype, vr1, vr2, vr,vi);
   if(iofl)
    {
     ierr=iofl;
    }
   else
    {
     if(jrt)
      {
       jret = jrt;
      }
     else
      {
       jret = jtyp;
      }
    }
  }

err:
if(ierr)
 {
  *vr=0.0;
  *vi=0.0;

  idf_math_sub = 0;
  idf_math_err = ierr;
 }

  *typ=jret;

return ierr;
}


void idf_matherr_prn()
{
#define IDF_MATHERR_NMES 12

static char *ErrMath[IDF_MATHERR_NMES] = {

/* 1*/ "overflow in arithmetics",
/* 2*/ "division by zero",
/* 3*/ "out of domain",
/* 4*/ "out of range",
/* 5*/ "improper math operand",
/* 6*/ "improper math function",
/* 7*/ "improper complex math function",
/* 8*/ "inconsistent types of binary math token",
/* 9*/ "unknown math function",
/*10*/ "incompatible types of arguments in function",
/*11*/ "imprpoper number of arguments in function",
/*  */ "unknown"
                          };
char *ssub;
int isub,imes;

   imes = idf_math_err;
if(imes > 0)
 {
  if(imes > IDF_MATHERR_NMES)
     imes = IDF_MATHERR_NMES;
     imes--;

  fprintf(stdout,"math error: %s\n",ErrMath[imes]);

     isub = idf_math_sub;
  if(isub)
   {
       ssub = idf_keywF_name(isub);
    if(ssub != NULL)
     {
      fprintf(stdout,"math function=%s\n",ssub);
     }
   }
 }

}
