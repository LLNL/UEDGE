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


#include "idf.i"


/*maximal number of standard math constants */
#define IDF_MLIST_MAX 27

/*
 -----------------------
   numerical constants
 -----------------------

IDF_RAD      1 radian         5.7295779513082320877E+1
IDF_DEG      1 degree         1.7453292519943295769E-2
IDF_EU       Euler number     0.5772156649015328606
IDF_PI       pi number        3.1415926535897932385
IDF_2PI      2*pi             6.2831853071795864769
IDF_PIPI     pi*pi            9.8696044010893586188
IDF_PIOVER2  pi/2             1.5707963267948966192
IDF_PIOVER3  pi/3             1.0471975511965977462
IDF_PIOVER4  pi/4             0.78539816339744880962
IDF_1_PI     1/pi             0.31830988618379067154
IDF_2_PI     2/pi             0.63661977236758134308
IDF_EXP      exp(1)           2.7182818284590453254
IDF_EXP2     exp(2)           7.3890560989306502272
IDF_EXPI     exp(pi)          2.3140692632779269006
IDF_EXPU     exp(Eu)          1.7810724179901979852
IDF_EXPE     exp(exp(1))      1,5154262241479264190E+1
IDF_SEXP     sqrt(exp(1))     1.6487212707001281468
IDF_SQRT2    sqrt(2)          1.4142135623730950488
IDF_SQRT3    sqrt(3)          1.7320508075688772935
IDF_SQRT5    sqrt(5)          2.2360679774997896964
IDF_SPI      sqrt(pi)         1.7724538509055160273
IDF_S2PI     sqrt(2*pi)       2.5066282746310005024
IDF_LN2      ln(2)            0.69314718055994530942
IDF_LN3      ln(3)            1.0986122886681096914
IDF_LN10     ln(10)           2.30258509299404568402
IDF_LNPI     ln(pi)           1.1447298858494001741
IDF_LNU      ln(Eu)          -0.54953931298164482234
*/

     
#if IDF_CPP_ == 1
int idf_mlist(char* name, int nn, double* val)
#else
int idf_mlist(name,nn, val)
char   *name;
int     nn;
double *val;
#endif
{

static char  *Clist[IDF_MLIST_MAX] = {
               "RAD",      "DEG",        "EU",
                "PI",      "2PI",      "PIPI",
           "PIOVER2",  "PIOVER3",   "PIOVER4",
              "1_PI",     "2_PI",       "EXP",
              "EXP2",     "EXPI",      "EXPU",
              "EXPE",     "SEXP",     "SQRT2",
             "SQRT3",    "SQRT5",       "SPI",
              "S2PI",      "LN2",       "LN3",
              "LN10",     "LNPI",       "LNU"
                                     };

static int   NClist[IDF_MLIST_MAX] = {
                  3,    3,    2,
                  2,    3,    4,
                  6,    6,    6,
                  4,    4,    3,
                  4,    4,    4,
                  4,    4,    5,
                  5,    5,    3,
                  4,    3,    3,
                  4,    4,    3
                                     };


static double Vlist[IDF_MLIST_MAX] = {
 5.7295779513082320877E+1,  1.7453292519943295769E-2,
 0.5772156649015328606   ,  3.1415926535897932385   ,
 6.2831853071795864769   ,  9.8696044010893586188   ,
 1.5707963267948966192   ,  1.0471975511965977462   ,
 0.78539816339744880962  ,  0.31830988618379067154  ,
 0.63661977236758134308  ,  2.7182818284590453254   ,
 7.3890560989306502272   ,  2.3140692632779269006   ,
 1.7810724179901979852   ,  1.5154262241479264190E+1,
 1.6487212707001281468   ,  1.4142135623730950488   ,
 1.7320508075688772935   ,  2.2360679774997896964   ,
 1.7724538509055160273   ,  2.5066282746310005024   ,
 0.69314718055994530942  ,  1.0986122886681096914   ,
 2.30258509299404568402  ,  1.1447298858494001741   ,
-0.54953931298164482234
                                     };

int k,n;
register int i,j;
char *s;


  i=0;
while(i<IDF_MLIST_MAX)
 {
  s = Clist[i];
  n = NClist[i];

  if(n == nn)
   {
    k=1;
    for(j=0;j<n;j++)
     {
      if(s[j] != name[j])
       {
        k=0;
        break;
       }
     }
    if(k) break;
   }

  i++;
 }

if(i==IDF_MLIST_MAX)
 {
  *val = 0.0;
  k=0;
 }
else
 {
  k=1;
  *val = Vlist[i];
 }

return k;
}
