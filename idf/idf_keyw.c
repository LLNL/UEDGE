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



#if IDF_CPP_ == 1
int idf_keyw(char* name)
#else
int idf_keyw(name)
char *name;
#endif
{
static int nK = 14;

static char *KeyW[14] = {
              "char",          "character",
               "int",            "integer",
              "long",             "float",
            "double",              "void",
            "static",              "text",
            "string",           "complex",
            "extern",             "local" 
                        };

static int KeyL[14] = {
                   4,                   8,
                   3,                   7,
                   4,                   5,
                   6,                   4,
                   6,                   4,
                   6,                   7,
                   6,                   5
                      };

/* char  string text int  long  float  double  complex  void  static extern*/
/*    1       2    3   4     5      6       7        8    10     11   12   */
static int KeyT[15] = {
                   1,                   1,
                   4,                   4,
                   5,                   6,
                   7,                  10,
                  11,                   3,
                   2,                   8,
                  12,                  11,
                   0 
                      };
int k,n,nn;
register int i,j;
unsigned int ci;
char *s;
char c;

nn = strlen(name);

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
      ci = name[j];
       c = tolower(ci);

      if(s[j] != c)
       {
        k=0;
        break;
       }
     }
    if(k) break;
   }

  i++;
 }

       k = KeyT[i];
return k;
}



#if IDF_CPP_ == 1
int idf_keyw_void(int jtyp)
#else
int idf_keyw_void(jtyp)
int jtyp;
#endif
{
static int nK = 13;

static int KeyVv[13] = {

/* char string text int long float double complex complex void  static extern*/
     0,     0,   0,  0,   0,    0,     0,      0,      0,   1,      1,     1,
     0
                      };

int k;
register int j;


     if(jtyp < 1 ) j=nK;
else if(jtyp > nK) j=nK;
else               j=jtyp;
                   j--;

       k = KeyVv[j];
return k;
}



#if IDF_CPP_ == 1
int idf_keyw_local(int jtyp)
#else
int idf_keyw_local(jtyp)
int jtyp;
#endif
{
static int nK = 13;

static int KeyVl[13] = {

/* char string text int long float double complex complex void  local  extern*/
     0,     0,   0,  0,   0,    0,     0,      0,      0,   0,      1,     0,
     0
                      };

int k;
register int j;

     if(jtyp < 1 ) j=nK;
else if(jtyp > nK) j=nK;
else               j=jtyp;
                   j--;

       k = KeyVl[j];
return k;
}



#if IDF_CPP_ == 1
int idf_keyw_digit(int jtyp)
#else
int idf_keyw_digit(jtyp)
int jtyp;
#endif
{
static int nK = 13;

static int KeyVV[13] = {

/* char string text int long float double complex complex void  local  extern*/
     0,     0,   0,  1,   1,    1,     1,      1,       1,   1,     1,     0,
     0
                      };

int k;
register int j;

     if(jtyp < 1 ) j=nK;
else if(jtyp > nK) j=nK;
else               j=jtyp;
                   j--;

       k = KeyVV[j];
return k;
}



#if IDF_CPP_ == 1
int idf_keyw_cnv(int k1, int k2)
#else
int idf_keyw_cnv(k1,k2)
int k1;
int k2;
#endif
{
/*
 convolution of two data types assigned to name
*/
static int nK = 12;

static char KeyVs[145] = {

/*           char  str txt  int  long float double cmplf cmpld void local ext*/
/*  char */  '1', '0', '0', '0', '0',  '0',  '0',   '0',  '0', '1',  '1',  '1',
/*string */  '0', '2', '0', '0', '0',  '0',  '0',   '0',  '0', '2',  '2',  '2',
/*  text */  '0', '0', '3', '0', '0',  '0',  '0',   '0',  '0', '3',  '3',  '3',
/*   int */  '0', '0', '0', '4', '0',  '0',  '0',   '0',  '0', '4',  '4',  '4',
/*  long */  '0', '0', '0', '0', '5',  '0',  '0',   '0',  '0', '5',  '5',  '5',
/* float */  '0', '0', '0', '0', '0',  '6',  '0',   '8',  '9', '6',  '6',  '6',
/*double */  '0', '0', '0', '0', '0',  '0',  '7',   '9',  '9', '7',  '7',  '7',
/*fcmplx */  '0', '0', '0', '0', '0',  '0',  '0',   '8',  '9', '8',  '8',  '8',
/*dcmplx */  '0', '0', '0', '0', '0',  '0',  '0',   '9',  '9', '9',  '9',  '9',
/*  void */  '1', '2', '3', '4', '5',  '6',  '7',   '9',  '9', 'v',  'v',  'v',
/* local */  '1', '2', '3', '4', '5',  '6',  '7',   '9',  '9', 'v',  'v',  '0',
/*extern */  '1', '2', '3', '4', '5',  '6',  '7',   '9',  '9', 'v',  '0',  'v',
             '0'
                      };

int k;
char c;
register int j;

     if((k1<1)||(k1>nK)) j=145;
else if((k2<1)||(k2>nK)) j=145;
else j = (k1-1)*nK + k2;

                  j--;
        c = KeyVs[j];
     if(c == 'v')
      {
       k = IDF_FRMT_VOID;
      }
     else
      {
       k = c - '0';
      }

return k;
}



#if IDF_CPP_ == 1
int idf_keyw_frmt(int kf,int kd)
#else
int idf_keyw_frmt(kf,kd)
int kf;
int kd;
#endif
{
/*
 convolution of format type and data type assigned to name
*/
static int nK = 12;
static int nF = 12;

static char KeyVf[145] = {

/*           char  str txt  int  long float double cmplf cmpld void local ext*/
/*  char */  '1', 'V', 'V', '0', '0',  '0',  '0',   '0',  '0', '1',  '0', '0',
/*string */  'V', '2', 'V', '0', '0',  '0',  '0',   '0',  '0', '2',  '0', '0',
/*  text */  'V', 'V', '3', '0', '0',  '0',  '0',   '0',  '0', '3',  '0', '0',
/*   int */  '0', '0', '0', '4', 'v',  'v',  'v',   'v',  'v', 'v',  '0', '0',
/*  long */  '0', '0', '0', '5', '5',  'v',  'v',   'v',  'v', 'v',  '0', '0',
/* float */  '0', '0', '0', 'v', 'v',  '6',  'v',   'v',  'v', 'v',  '0', '0',
/*double */  '0', '0', '0', 'v', 'v',  '7',  '7',   'v',  'v', 'v',  '0', '0',
/*fcmplx */  '0', '0', '0', 'v', 'v',  'v',  'v',   '8',  'v', 'v',  '0', '0',
/*dcmplx */  '0', '0', '0', 'v', 'v',  'v',  'v',   '9',  '9', '9',  '0', '0',
/*  void */  '1', '2', '3', '4', '5',  '6',  '7',   '8',  '9', 'v',  '0', '0',
/* local */  '0', '0', '0', '0', '0',  '0',  '0',   '0',  '0', '0',  '0', '0',
/*extern */  '0', '0', '0', '0', '0',  '0',  '0',   '0',  '0', '0',  '0', '0',
             '0'
                      };

int k;
char c;
register int j;

     if((kf<1)||(kf>nF)) j=145;
else if((kd<1)||(kd>nK)) j=145;
else j = (kf-1)*nK + kd;

                  j--;
        c = KeyVf[j];
     if(c=='v' || c=='V')
      {
       k = IDF_FRMT_VOID;
      }
     else
      {
       k = c - '0';
      }

return k;
}
