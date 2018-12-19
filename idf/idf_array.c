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
#include "idfusr.h"



#if IDF_CPP_ == 1
int idf_get_array(char* name, char* frmt, void* value)
#else
int idf_get_array(name, frmt, Mdim, value)
char *name;
char *frmt;
void *value;
#endif
{
static int align=0;
int ierr;

ierr = idf_get(name, frmt, align, value);

return ierr;
}



#if IDF_CPP_ == 1
int idf_get_array_(char* name, int nam,
                   char* frmt, int nfrmt,
                   void* value)
#else
int idf_get_array_(name,nam, frmt,nfrmt, value)
char *name;
int   nam;
char *frmt;
int   nfrmt;
void *value;
#endif
{
static int align = 0;
int ierr;

ierr = idf_get_(name,nam, frmt,nfrmt, align, value);

return ierr;
}


