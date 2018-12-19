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
void idf_get_err_prn(char* name, int n,
                     int jerr, int lvl)
#else
void idf_get_err_prn(name,n, jerr,lvl)
char *name;
int   n;
int   jerr;
int   lvl;
#endif
{
#define IDF_ERR_M_N 29

static char *ErrM[IDF_ERR_M_N] = {

/* 1 ,0*/ "space symbol inside input name of data",
/* 2 ,0*/ "improper symbol inside input name of data",
/* 3 ,0*/ "too long input name for data or no EoS",
/* 4 ,0*/ "improper input format symbol",
/* 5 ,0*/ "empty list of named datas",
/* 6 ,0*/ "no dimensions for target array",
/* 7 ,0*/ "improper dimension for target array",
/* 8 ,0*/ "inconsistency in number of dimensions",
/* 9 ,0*/ "improper offset index in dimension list",
/* 10,0*/ "offset is equal or exceed the maximal number of elements",
/* 11,0*/ "unable to move to the beginning of data unit",
/* 12,2*/ "error in treating data unit",
/* 13,3*/ "error in conversion of data unit body",
/* 14,4*/ "error in conversion of formular results",
/* 15,0*/ "no such name in file",
/* 16,0*/ "no data units for this name",
/* 17,1*/ "dimensions are prescribed to single data unit in file",
/* 18,1*/ "no data units found for this name",
/* 19,1*/ "format type inconsistent with that assigned to name",
/* 20,5*/ "unable to convert",
/* 21,0*/ "empty input name of data",
/* 22,0*/ "improper size of name of data",
/* 23,0*/ "improper input format string",
/* 24,0*/ "memory allocation mode is not permitted",
/* 25,0*/ "improper input format parameters",
/* 26,0*/ "format unit error",
/* 27,0*/ "duplicated names are not allowed in memory allocation mode",
/* 28,0*/ "improper order for function call, you must first init and open file",
/*   ,0*/ "unknown error"
                      };
register int j;
char *s=NULL;

if(n)
 {
  fprintf(stdout,"input name(%3d)=%s\n",n,name);  
 }

if(jerr>0)
 {
     j = jerr;
  if(j > IDF_ERR_M_N) j = IDF_ERR_M_N;
   {
     j--;
    s = ErrM[j];

    fprintf(stdout,"%s\n",s);
   }
 }


if(lvl>0)
 {
  /* name was obtained*/
  idf_data_prn0();
 }

if(lvl>1)
 {
  /* data unit was obtained*/
  idf_data_prn1();

  if(lvl==3)
   {
    /* data buffer was obtained*/
    idf_dbuf_prn();
   }

  else if(lvl==4)
   {
    /* formular was obtained*/
    idf_form_prn();
    idf_matherr_prn();
   }

  else if(lvl==5)
   {
    /* value was obtained*/
    idf_val_prn();
   }

  idf_txt_err_prn();
 }

idf_frmt_err_prn();
idf_frmtreat_err_prn();
}



void idf_list_err_prn()
{
#define IDF_ERR_L_N 23

static char *ErrL[IDF_ERR_L_N] = {

/* 1,0*/ "list structure has been already init",
/* 2,0*/ "not enough memory to init name list",
/* 3,1*/ "improper name,unable to put it into list",
/* 4,1*/ "unable to put local name to nonlocal list",
/* 5,0*/ "no more memory to place name in nonlocal list",
/* 6,0*/ "non-local list is full",
/* 7,1*/ "unable to put nonlocal name to local list",
/* 8,0*/ "no more memory to place name in local list",
/* 9,0*/ "local list is full",
/*10,2*/ "error in name search",
/*11,0*/ "unexpected overflow of name list",
/*12,0*/ "empty non-local list",
/*13,7*/ "data name redefines the standard name",
/*14,7*/ "data name redefines the standard function name",
/*15,7*/ "data name is in non-local and local lists",
/*16,3*/ "unable to move to local data field",
/*17,3*/ "no data units for local name",
/*18,3*/ "error in treating local data unit",
/*19,4*/ "error in conversion of local data unit body",
/*20,5*/ "error in conversion of local formular results",
/*21,6*/ "local type inconsistent with that in file",
/*22,6*/ "unable to convert to local",
/*  ,0*/ "unknown list error"
                                  };

register int j;
int jerr,lvl,irec;
char *s=NULL;
IDF_NAME_LIST *NameList;

   NameList = idf_list_address();
if(NameList != NULL)
{
   jerr = NameList->kerr;
   lvl  = NameList->lvl;
   irec = NameList->irec;

if(jerr>0)
 {
     j = jerr;
  if(j > IDF_ERR_L_N) j = IDF_ERR_L_N;
   {
     j--;
    s = ErrL[j];

    fprintf(stdout,"catalog: %s\n",s);
   }
 }

if(lvl < 1)
 {
  /* no additional info*/;
 }

else if(lvl==1)
 {
  idf_name_errprn();
 }

else if(lvl==2)
 {
  /* error in obtaing the name*/
  idf_name_errprn();
  idf_txt_err_prn();
 }

else if(lvl<6)
 {
  /*name was obtained*/
  idf_list_prnL(irec);

  /* data unit */
  idf_data_prn1();

  if(lvl==4)
   {
    /* data buffer was obtained*/
    idf_dbuf_prn();
   }

  else if(lvl==5)
   {
    /* formular was obtained*/
    idf_form_prn();
    idf_matherr_prn();
   }

  idf_txt_err_prn();  
 }

else if(lvl==6)
 {
  /*name was obtained*/
  idf_list_prnL(irec);

  /* data unit */
  idf_data_prn1();

  /* value was obtained*/
  idf_val_prn(); 
 }

else if(lvl==7)
 {
  /*improper local name was obtained*/
  idf_list_prnL(irec);
 }

else if(lvl==8)
 {
  /*improper name was obtained*/
  idf_list_prn(irec);
 }

}

}
