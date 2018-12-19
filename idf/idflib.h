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


#ifndef IDFLIB_H
#define IDFLIB_H



#include <math.h>

/* 
IDF uses the following macros normally defined
in the standard C header files:
limits.h : LONG_MAX, INT_MAX, CHAR_MAX.
values.h : LN_MAXDOUBLE
float.h  : DBL_MAX, DBL_MIN, FLT_MAX, FLT_MIN, DBL_MAX_10_EXP.

If limits.h, values.h, or float.h header file does not exist,
delete the corresponding #include command and 
either define the remainder macros for your machine
or try to use their default values declared further in this file.
 */
#include <limits.h>
#include <float.h>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

/*I am using 'malloc' and 'free' C functions,
  the prototypes of which should be historically in malloc.h.
  Recently malloc merged to <stdlib.h>. 
  However ansi C might expect malloc.h in /usr/include.
  GNU can automatically detect if system needs this header
   and defines HAVE_MALLOC_H macro 

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
*/



/*----------------------------------------
      define macro IDF_IEEE as 1
 if your software supports IEEE standard.
------------------------------------------*/
#define IDF_IEEE 1



/*------------------------------------
  Default settings instead of float.h
*-------------------------------------*/

/* DBL_MAX is the largest NORMALIZED double number */
#ifndef DBL_MAX
#if IDF_IEEE == 1
#define DBL_MAX 1.79769313486231570e+308
#else
#define DBL_MAX 1.70141183460469229e+38
#endif
#endif

/* DBL_MIN is the smallest NORMALIZED double number */
#ifndef DBL_MIN
#if IDF_IEEE == 1
#define DBL_MIN 2.2250738585072014e-308
#else
#define DBL_MIN 2.93873587705571877e-39
#endif
#endif

/* DBL_MAX_10_EXP is the largest base 10 exponent of double number */
#ifndef DBL_MAX_10_EXP
#if IDF_IEEE == 1
#define DBL_MAX_10_EXP 308
#else
#define DBL_MAX_10_EXP 38
#endif
#endif

/* FLT_MAX is the largest NORMALIZED float number */
#ifndef FLT_MAX
#if IDF_IEEE == 1
#define FLT_MAX 3.40282347e+38f
#else
#define FLT_MAX 1.701411835e+38f
#endif
#endif

/* FLT_MIN is the smallest NORMALIZED float number */
#ifndef FLT_MIN
#if IDF_IEEE == 1
#define FLT_MIN 1.17549435e-38f
#else
#define FLT_MIN 2.938735877e-39f
#endif
#endif

/* INT_MAX is the largest int number */
#ifndef INT_MAX
#define INT_MAX 2147483647
#endif

/* LONG_MAX is the largest long int number */
#ifndef LONG_MAX
#define LONG_MAX 2147483647L
#endif

/* CHAR_MAX is the largest char number */
#ifndef CHAR_MAX
#define CHAR_MAX 255
#endif

/* pi */
#ifndef M_PI
#define M_PI 3.1415926535897932385e0
#endif

/* pi/2 */
#ifndef M_PI_2
#define M_PI_2 1.5707963267948966192e0
#endif

/* log(10) */
#ifndef M_LN10
#define M_LN10 2.3025850929940456840e0
#endif

/* the natural log of the largest double */
#ifndef LN_MAXDOUBLE
#if IDF_IEEE == 1
#define LN_MAXDOUBLE 7.09782712893383973e2
#else
#define LN_MAXDOUBLE 8.80296919311130495e1
#endif
#endif



/*---------------------------------------------------- 
   macro IDF_CPP_STYLE controls the declaration 
of C function arguments as well as function prototype

If IDF_CPP_STYLE is defined as 1
then function prototypes will be declared in 
1989 Standard C style:

type fun(type Arg1, type Arg2)
{ function body }

If IDF_CPP_STYLE is defined as 0
then function prototypes will be declared in 
Pre-Standard C style:

type fun(Arg1,Arg2)
     type Arg1; type Arg2;
{ function body } 
--------------------------------------------------------*/
#define IDF_CPP_STYLE 1



/*------------------------------------------
  define IDF_CPP as 0 if you use a C compiler
      to create IDF object library 
--------------------------------------------*/
#define IDF_CPP 0



/*
------------------------------------------
 Specification of naming conventions for
   C functions to be called from FORTRAN
------------------------------------------
 Define the macros IDF_FORTRAN_TYPE:
    for   ANSI_F77  as   1
    for       CRAY  as   1
    for        AIX  as   2
    for       HPUX  as   2
    for        SVS  as   2
    for      other  as   3
*/
#define IDF_FORTRAN_TYPE 3


#if IDF_FORTRAN_TYPE == 1
#define IDF_FORTRAN(fid_, fid, FID) FID 
#endif
#if IDF_FORTRAN_TYPE == 2
#define IDF_FORTRAN(fid_, fid, FID) fid 
#endif
#if IDF_FORTRAN_TYPE == 3
#define IDF_FORTRAN(fid_, fid, FID) fid_ 
#endif

#ifndef IDF_FORTRAN
#define IDF_FORTRAN(fid_, fid, FID) fid_ 
#endif



#endif   /* IDFLIB_H */
