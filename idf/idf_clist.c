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
#include "idf.i"


#if IDF_CPP_ == 1
int idf_clist(char* name, int nn, double* val)
#else
int idf_clist(name,nn, val)
char *name;
int nn;
double *val;
#endif
{
static char Seq[] = "IDF_";

int k,l,n;
register int i;
char *s;

if(nn > 4)
 {
  /* check for name started with IDF_ */
  l=0;
  for(i=0;i<4;i++)
   {
    if(name[i] != Seq[i])
     {
      l=1;
      break;
     }
   }

  if(l==0)
   {
     n = nn - 4;
     s = name + 4;

       /* search in math constant list */
       k = idf_mlist(s,n, val);
    if(k==0)
     {
      /* search in physical constant list */
      k = idf_plist(s,n, val);
     }
   }
  else
   {
    *val = 0.0;
    k=0;
   }
 }

else
 {
  *val = 0.0;
  k=0;
 }

return k;
}
