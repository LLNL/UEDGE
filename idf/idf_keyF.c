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



#define IDF_KEYW_MAX 33

#define IDF_KEYW_NUM 33
#define IDF_KEYW_NUM_1 34




static char *KeyW[IDF_KEYW_MAX] = {
               "abs",                "int", 
              "sqrt",                "log", 
             "log10",                "exp", 
               "cos",                "sin", 
               "tan",               "acos", 
              "asin",               "atan", 
               "pow",               "sinh", 
              "cosh",               "tanh",
             "asinh",              "acosh", 
             "atanh",            "complex",
             "cmplx",               "real",
              "imag",              "polar",
               "arg",                "min",
               "max",              "hypot",
            "double",               "conj",
          "sqrt1pz2",           "sqrt1mz2",
             "pow10"
                                  };

static int KeyL[IDF_KEYW_MAX] = {
                   3,                   3,
                   4,                   3,
                   5,                   3,
                   3,                   3,
                   3,                   4,
                   4,                   4,
                   3,                   4,
                   4,                   4,
                   5,                   5,
                   5,                   7,
                   5,                   4,
                   4,                   5,
                   3,                   3,
                   3,                   5,
                   6,                   4,
                   8,                   8,
                   5
                                };



#if IDF_CPP_ == 1
char* idf_keywF_name(int inum)
#else
char* idf_keywF_name(inum)
int inum;
#endif
{
char *c;

 if(inum<1 || inum>IDF_KEYW_MAX)
  {
   c = NULL;
  }
 else
  {
   c = KeyW[inum-1];
  }

return c;
}



#if IDF_CPP_ == 1
int idf_keywF_lname(int inum)
#else
int idf_keywF_lname(inum)
int inum;
#endif
{
int c;

 if(inum<1 || inum>IDF_KEYW_MAX)
  {
   c = 0;
  }
 else
  {
   c = KeyL[inum-1];
  }

return c;
}



#if IDF_CPP_ == 1
int idf_keywF(char* name, int nn)
#else
int idf_keywF(name,nn)
char *name;
int   nn;
#endif
{
static int nK = IDF_KEYW_MAX;

int k,kk,n;
register int i,j;
char *s;
char c;

 kk=0;
  i=0;
while(i<nK)
 {
  s = KeyW[i];
  n = KeyL[i];

  if(n == nn)
   {
    k=1;
    for(j=0;j<n;j++)
     {
         c = s[j];
      if(c != name[j])
       {
        k=0;
        break;
       }
     }

    if(k==0)
     {
      k=1;
      for(j=0;j<n;j++)
       {
           c = toupper(s[j]);
        if(c != name[j])
         {
          k=0;
          break;
         }
       }
     }

    if(k)
     {
      kk = i + 1;
      break;
     }
   }

  i++;
 }

return kk;
}



#if IDF_CPP_ == 1
int idf_keywF_narg(int jtyp)
#else
int idf_keywF_narg(jtyp)
int jtyp;
#endif
{
/* number of arguments */

static int nK = IDF_KEYW_NUM;

static int KeyV[IDF_KEYW_NUM_1] = {
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 2, 1, 1,
        1, 1, 1, 1, 2,
        1, 1, 1, 2, 1,
        2, 2, 2, 1, 1,
        1, 1, 1, 
        0
                      };

int k;
register int j;

        j = jtyp;
     if(j < 1 ) j=nK;
else if(j > nK) j=nK;
else            j--;

       k = KeyV[j];
return k;
}



#if IDF_CPP_ == 1
int idf_keywF_targ(int jtyp)
#else
int idf_keywF_targ(jtyp)
int jtyp;
#endif
{
/* argument type: 0-void D-double W-complex*/

static int nK = IDF_KEYW_NUM;

static int KeyT[IDF_KEYW_NUM_1] = {
         0 ,    IDF_FRMT_D ,
         0 ,             0 ,
         0 ,             0 ,
         0 ,             0 ,
         0 ,             0 ,
         0 ,             0 ,
         0 ,             0 ,
         0 ,             0 ,
         0 ,             0 ,
         0 ,    IDF_FRMT_D ,
IDF_FRMT_D ,             0 ,
         0 ,    IDF_FRMT_D ,
         0 ,             0 ,
         0 ,    IDF_FRMT_D ,
IDF_FRMT_D ,    IDF_FRMT_W ,
         0 ,             0 ,
         0
                      };

int k;
register int j;

        j = jtyp;
     if(j < 1 ) j=nK;
else if(j > nK) j=nK;
else            j--;

       k = KeyT[j];
return k;
}



#if IDF_CPP_ == 1
int idf_keywF_tfun(int jtyp)
#else
int idf_keywF_tfun(jtyp)
int jtyp;
#endif
{
/* return type: 0-void D-double W-complex*/

static int nK = IDF_KEYW_NUM;

static int KeyT[IDF_KEYW_NUM_1] = {
          0 ,        IDF_FRMT_D ,
          0 ,                 0 ,
          0 ,                 0 ,
          0 ,                 0 ,
          0 ,                 0 ,
          0 ,                 0 ,
          0 ,                 0 ,
          0 ,                 0 ,
          0 ,                 0 ,
          0 ,        IDF_FRMT_W ,
 IDF_FRMT_W ,        IDF_FRMT_D ,
 IDF_FRMT_D ,        IDF_FRMT_W ,
 IDF_FRMT_D ,                 0 ,
          0 ,        IDF_FRMT_D ,
 IDF_FRMT_D ,        IDF_FRMT_W ,
          0,                  0 ,
          0,
          0
                      };

int k;
register int j;

        j = jtyp;
     if(j < 1 ) j=nK;
else if(j > nK) j=nK;
else            j--;

       k = KeyT[j];
return k;
}
