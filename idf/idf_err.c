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


static int idf_level=0;
static int idf_error=0;
static int idf_prnt =0;



void idf_err_clear()
{
 idf_level = 0;
 idf_error = 0;
 idf_prnt  = 0;

 idf_matherr_nul  ();
 idf_txt_err_nul  ();
 idf_file_clearerr();
 idf_list_clearerr();
}



#if IDF_CPP_ == 1
void idf_err_put(int lvl, int err)
#else
void idf_err_put(lvl,err)
int lvl,err;
#endif
{
 idf_level = lvl;
 idf_error = err;
}


void idf_err_prn()
{
#define IDF_ERR_MAIN_N 31

static char *ErrMAIN[IDF_ERR_MAIN_N] = {
/* 1, 1*/ "empty file name",
/* 2, 1*/ "too long file name",
/* 3, 1*/ "improper symbol in file name",
/* 4, 1*/ "unable to update input file name",
/* 5, 1*/ "no memory to allocate file C-name",
/* 6, 2*/ "file stucture init already",
/* 7, 2*/ "no memory for file structure",
/* 8, 3*/ "format stucture init already",
/* 9, 3*/ "no memory for format structure",
/*10, 3*/ "improper format initialization",
/*11, 4*/ "name stucture init already",
/*12, 4*/ "no memory for name structure",
/*13, 5*/ "unable to init list structure",
/*14, 6*/ "data stucture init already",
/*15, 6*/ "no memory for data structure",
/*16, 7*/ "buf stucture init already",
/*17, 7*/ "no memory for buf structure",
/*18, 8*/ "val stucture init already",
/*19, 8*/ "no memory for val structure",
/*20, 9*/ "unable to open input file",
/*21,10*/ "improper usage of 'open', call 'init' first",
/*22,11*/ "unable to create catalog of names",
/*23,11*/ "improper names in catalog",
/*24,11*/ "unable to rewind file",
/*25,11*/ "unable to fulfil catalog with local values",
/*26, 7*/ "formular unit stucture init already",
/*27, 7*/ "no memory for formular unit structure",
/*28, 8*/ "formular stucture init already",
/*29, 8*/ "no memory for formular structure",
/*30, 7*/ "error in float structure",
/*     */ "unknown error"
                                 }; 
int jerr,lvl;
register int j;
char *s;

if(idf_prnt==0)
{
    lvl = idf_level;
   jerr = idf_error;
if(jerr>0)
 {
     j = jerr;
  if(j > IDF_ERR_MAIN_N) j = IDF_ERR_MAIN_N;
   {
     j--;
    s = ErrMAIN[j];

    fprintf(stdout,"IDF error=%d level=%d:\n%s\n",
            jerr,lvl,s);
   }

  if(jerr==10)
   {
    fprintf(stdout,"format table(%2d)=%s\n",
            IDF_FRMT_NTYP,IDF_FRMT_LIST);
   }

  if(lvl==11)
   {
    /* error message for catalog creating */
    idf_list_err_prn();
   }

  if(lvl > 10)
   {
    /* error message associated with file operation*/
    idf_file_errprn();
   }

 }

idf_prnt = 1;
}

}
