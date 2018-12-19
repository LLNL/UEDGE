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




int idf_list_ffl()
{
/*
--------------------------------------
 calculate values of local variables
--------------------------------------
*/
int ierr;
int n;
register int i;
int l,jtype;
int flag,k,ibrk,kup;
char *s;
IDF_NAME_L *A;
IDF_NAME_LIST *NameList;

       A = NULL;
NameList = NULL;
    ierr = 0;
     kup = 0;

   n = idf_list_numL();
if(n < 1)
 {
  /* nothing to do, empty local list */
  goto fin;
 }

   NameList = idf_list_address();
if(NameList == NULL)
 {
  fprintf(stdout,"list structure was lost\n");
  ierr = 1;
  goto fin;
 }

    i=0;
  while(i<n)
   {
    A = NameList->ListL[i];

    /* name attributes*/
    s = A->name;
    l = A->length;

    /* data type assigned in file*/
    jtype = A->jtype;


       flag = idf_data_beginL(A);
    if(flag)
     {
      /* unable to move to data field*/
      NameList->lvl  = 3;
      NameList->kerr = 16;
      NameList->irec = i;
      ierr=1;
      break;
     }

          ibrk=0;
           kup=0;
    while(ibrk==0)
     {

           flag = idf_data_value();
        if(flag<0)
         {
          NameList->lvl  = 2-flag;
          NameList->kerr = 17-flag; /*1-3*/
          NameList->irec = i;
          ierr=1;
          ibrk=2;
         }
        else if(flag==1)
         {
          /* no more data */
          ibrk=1;
         }
        else
         {
          /*data unit up */
             k = idf_val_lcns(A);
          if(k)
           {
            NameList->lvl  = 6;
            NameList->kerr = 20+k; /*1-2*/
            NameList->irec = i;
            ierr=1;
            ibrk=2;
           }
          else
           {
            kup++;
           }
         }

     }/*search over data units */

    if(!ierr)
     {
      if(!kup)
       {
        /* no data unit */
        NameList->lvl  = 3;
        NameList->kerr = 17;
        NameList->irec = i;
        ierr=1;
        flag=1;
       }
      else
       {
        flag=0;
       }
     }
    else
     {
      flag=1;
     }

    if(flag) break;

    i++;
   } /* nonloc list search */

fin:
return ierr;
}
