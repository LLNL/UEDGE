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
int idf_extname_check(char* name, int n)
#else
int idf_extname_check(name,n)
char *name;
int   n;
#endif
{
int k;
register int i;
unsigned int c;

  k = -1;
  i = 0;
while(i<n)
 {
  c = (unsigned int) name[i];


  if(c == IDF_EOS)
   {
    k=4;
   }

  else if(i >= IDF_NAME_LENGTH)
   {
    k=3;
   }

  else if( !isascii(c) )
   {
#if IDF_FOREIGN_LANGUAGE == 1
    if(k == -1) k=0;
#else
    k=2;
#endif
   }

  else if(c=='_' || c=='%')
   {
    if(k != 0) k=2;
   }

  else if(c == ' ')
   {
    if(k != -1) k=1;
   }

  else if( isalpha(c) )
   {
    if(k == -1) k=0;
   }

  else if( isdigit(c) )
   {
    if(k == -1) k=2;
   }

  else
   {
    k=2;
   }

  if(k>0) break;

  i++;
 }/*i*/

if(k == -1) k=4;

return k;
}
