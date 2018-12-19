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


#define IDF_DTOS_ACC_MAX 15
#define IDF_DTOS_NBUF 27



#if IDF_CPP_ == 1
char* idf_dtos(double d, char* cntr)
#else
char* idf_dtos(d, cntr)
double d;
char *cntr;
#endif
{
static char DigStr[] = "0123456789*";
static char Buf[IDF_DTOS_NBUF];

int iacc;
char ce,c;
register int i,j;
int isig,iexp,ic;
double base;
char *str;

if(cntr==NULL)
 {
  /* no input format string */
  iacc=1;
  ce = 'e';
 }
else
 {
   i=0;
  ce=' ';
  iacc=0;
  ic=0;
  while(ic==0)
   {
    c = cntr[i];

    if(c==' ')
     {
      ;
     }

    else if( isdigit(c) )
     {
      iacc = iacc*10 + (c-'0');
     }

    else if(c=='e'||c=='d'||c=='g')
     {
      ce = c;
      if(iacc) ic=1;
     }

    else
     {
      ic=2;
     }

    i++;
   }

  if(ce==' ') ce='e';

  if(iacc > IDF_DTOS_ACC_MAX)
   {
    iacc = IDF_DTOS_ACC_MAX;
   }
  else if(iacc < 1)
   {
    iacc=1;
   }
 }

  i = 0;
str = Buf;

idf_float_decompose(d, &isig, &base, &iexp);

if(base < IDF_DOUBLE_MIN)
 {
  str[i] = '0';
      i++;
  str[i] = '.';
      i++;
  str[i] = '0';
      i++;
  str[i] = IDF_EOS;
  goto fin;
 }


if(isig < 0)
 {
  str[0] = '-';
  i++;
 }

if(base >= 1.0)
 {
     ic = (int)base;
  if(ic > 9) ic=10;
 }
else
 {
  ic = 0;
 }

  str[i] = DigStr[ic];
      i++;

str[i] = '.';
    i++;

base = base - (double)ic;

   isig=0;
      j=0;
while(j<iacc)
 {
  base = base*10.0;

     ic = (int)base;
  if(ic) isig++;

  str[i] = DigStr[ic];
      i++;

  base = base - (double)ic;
  j++;
 }

if(isig==0 && iacc!=1)
 {
  j=i;
  while(j>0)
   {
    j--;
    if(str[j] != '0') break;
   }
  i=j+2;
 }

if(iexp==0)
 {
  str[i] = IDF_EOS;
  goto fin;
 }

str[i] = ce;
    i++;

if(iexp<0)
 {
  str[i] = '-';
      i++;
  iexp = -iexp;
 }

if(iexp>=1000)
 {
  ic=iexp/1000;
  if(ic>9) j=10; else j=ic;

  str[i] = DigStr[j];
      i++;

  iexp = iexp - ic*1000;
 }

if(iexp>=100)
 {
  j=iexp/100;

  str[i] = DigStr[j];
      i++;

  iexp = iexp - j*10;
 }

if(iexp>=10)
 {
  j=iexp/10;

  str[i] = DigStr[j];
      i++;

  iexp = iexp - j*10;
 }

  str[i] = DigStr[iexp];
      i++;

str[i] = IDF_EOS;

fin:
return str;
}
