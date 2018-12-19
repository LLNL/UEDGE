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

 

IDF_VAL *FtxVal = NULL;



IDF_VAL* idf_val_address()
{
 return ( FtxVal );
}



int idf_val_ini()
{
int ierr=0;
long k;

if(FtxVal != NULL)
 {
  ierr=1;
 }

else
 {
        k = sizeof(IDF_VAL);
     FtxVal = (IDF_VAL*) malloc(k);
  if(FtxVal == NULL)
   {
    ierr=2;
   }
 }

return ierr;
}


void idf_val_end()
{
if(FtxVal != NULL)
 {
  free(FtxVal);
  FtxVal = NULL;
 }
}


void idf_val_nul()
{
 FtxVal->jtype  = 0;
}


void idf_val_prn()
{
int typ;

   typ = FtxVal->jtype;
if(typ == IDF_FRMT_C)
 {
  fprintf(stdout,"Cval=%c\n",FtxVal->Cval);
 }
else if(typ == IDF_FRMT_S || typ == IDF_FRMT_T)
 {
  fprintf(stdout,"Sval(%d)=%s\n",FtxVal->len,FtxVal->Sval);
 }
else if(typ == IDF_FRMT_I)
 {
  fprintf(stdout,"Ival=%d\n",FtxVal->Ival);
 }
else if(typ == IDF_FRMT_L)
 {
  fprintf(stdout,"Lval=%ld\n",FtxVal->Lval);
 }
else if(typ == IDF_FRMT_F)
 {
  fprintf(stdout,"Fval=%e\n",FtxVal->Fval[0]);
 }
else if(typ == IDF_FRMT_D)
 {
  fprintf(stdout,"Dval=%e\n",FtxVal->Dval[0]);
 }
else if(typ == IDF_FRMT_Z)
 {
  fprintf(stdout,"Zval=(%e , %e)\n",FtxVal->Fval[0],FtxVal->Fval[1]);
 }
else if(typ == IDF_FRMT_W)
 {
  fprintf(stdout,"Wval=(%e , %e)\n",FtxVal->Dval[0],FtxVal->Dval[1]);
 }

}



#if IDF_CPP_ == 1
int idf_val_put(int kfrmt)
#else
int idf_val_put(kfrmt)
int kfrmt;
#endif
{
/*-----------------------------
   convert buffer into value
  according to external format
 -----------------------------*/

unsigned int cerr = ' ';
int jerr,flag;
int Jtype,typ;
IDF_DBUF *Db;
int k;
char   *pcv;
int    *piv;
long   *plv;
float  *pfv;
double *pdv;

  Db = idf_dbuf_address();
jerr = 0;
 typ = 0;
flag = 0;

   Jtype = Db->jtype;
if(Jtype < 1)
 {
  /* no symbols in data unit buffer*/
  jerr=37;
 }

else if(kfrmt==IDF_FRMT_C)
 {
  pcv = &FtxVal->Cval;
     k = idf_dbuf_to_c(pcv);
  if(k) jerr=k;
  typ = kfrmt;
 }

else if(kfrmt==IDF_FRMT_S)
 {
  pcv = FtxVal->Sval;
  piv = &FtxVal->len;
     k = idf_dbuf_to_s(pcv,piv);
  if(k) jerr=6+k;
  typ = kfrmt;
 }

else if(kfrmt==IDF_FRMT_T)
 {
  pcv = FtxVal->Sval;
  piv = &FtxVal->len;
     k = idf_dbuf_to_t(pcv,piv);
  if(k) jerr=12+k;
  typ = kfrmt;
 }

else if(kfrmt==IDF_FRMT_I)
 {
  piv = &FtxVal->Ival;
#if IDF_DBUF_FUN == 1
     k = idf_dbuf_to_i_(piv);
#else
     k = idf_dbuf_to_i (piv);
#endif
  if(k) jerr=14+k;
  typ = kfrmt;
 }

else if(kfrmt==IDF_FRMT_L)
 {
  plv = &FtxVal->Lval;
#if IDF_DBUF_FUN == 1
     k = idf_dbuf_to_l_(plv);
#else
     k = idf_dbuf_to_l (plv);
#endif
  if(k) jerr=14+k;
  typ = kfrmt;
 }

else if(kfrmt==IDF_FRMT_F)
 {
  pfv = FtxVal->Fval;
#if IDF_DBUF_FUN == 1
     k = idf_dbuf_to_f_(pfv);
#else
     k = idf_dbuf_to_f (pfv);
#endif
  if(k) jerr=21+k;
  typ = kfrmt;
 }

else if(kfrmt==IDF_FRMT_D)
 {
  pdv = FtxVal->Dval;
#if IDF_DBUF_FUN == 1
     k = idf_dbuf_to_d_(pdv);
#else
     k = idf_dbuf_to_d (pdv);
#endif
  if(k) jerr=21+k;
  typ = kfrmt;
 }

else if(kfrmt==IDF_FRMT_Z)
 {
  pfv = FtxVal->Fval;
     k = idf_dbuf_to_z(pfv);
  if(k) jerr=28+k;
  typ = kfrmt;
 }

else if(kfrmt==IDF_FRMT_W)
 {
  pdv = FtxVal->Dval;
     k = idf_dbuf_to_w(pdv);
  if(k) jerr=28+k;
  typ = kfrmt;
 }

else if(kfrmt==IDF_FRMT_VOID)
 {
  pdv = FtxVal->Dval;
     k = idf_dbuf_to_v(pdv,&typ);
  if(k)
   {
    if(typ==IDF_FRMT_L)
     {
      jerr=14+k;
      FtxVal->Lval = (long)FtxVal->Dval[0];
     }
    else if(typ==IDF_FRMT_W)
     {
      jerr=28+k;
     }
    else
     {
      jerr=21+k;
     }
   }
  else
   {
    if(typ==IDF_FRMT_L)
     {
      FtxVal->Lval = (long)FtxVal->Dval[0];
     }    
   }
 }

else
 {
  /* improper format for buffer conversion */
  jerr=36;
 }

  FtxVal->jtype = typ;

if(jerr)
 {
  k = jerr + 115;
  idf_txt_err_put(k, Db->kstr,Db->kurs,Db->pos,cerr);
  flag = 1;
 }

return flag;
}



int idf_val_buf()
{
/*-----------------------------
   convert buffer into value
  according to buffer format
 -----------------------------*/
int kfrmt;
int k;
IDF_DBUF *Db;

  Db = idf_dbuf_address();

   kfrmt = Db->jtype;
if(kfrmt == IDF_FRMT_DIGIT)
 {
  kfrmt = IDF_FRMT_VOID;
 }

       k = idf_val_put(kfrmt);
return k;
}



int idf_val_form()
{
/*--------------------
  copy formular value 
  --------------------*/

unsigned int cerr = ' ';
int frmt;
int k;
IDF_FORM *Df;

  Df = idf_form_address();
   k = 0;

   frmt = Df->jtype;
if(frmt == IDF_FRMT_DIGIT)
 {
  FtxVal->jtype   = IDF_FRMT_D;
  FtxVal->Dval[0] = Df->Val[0];
 }
else if(frmt == IDF_FRMT_CMPLX)
 {
  FtxVal->jtype   = IDF_FRMT_W;
  FtxVal->Dval[0] = Df->Val[0];
  FtxVal->Dval[1] = Df->Val[1];
 }
else
 {
  FtxVal->jtype   = 0;
  k = 153;
  idf_txt_err_put(k, Df->kstr,Df->kurs,Df->kpos, cerr);
 }

return k;
}



#if IDF_CPP_ == 1
int idf_val_get(int jfrmt, void* val, int nval)
#else
int idf_val_get(jfrmt, val,nval)
int   jfrmt;
void *val;
int   nval;
#endif
{
/*
 -----------------------------
   convert buffer into value
 -----------------------------

returns 0 - OK
       -1 - inconsistency between input format and
            data type in file
       -2 - overflow in data unit conversion 
            according to input format 
*/
int k,flag;
int kfrmt;
int l;
float f;
double d;

flag  = 0;
k     = 0;
kfrmt = FtxVal->jtype;
if(kfrmt==IDF_FRMT_C)
 {
  k = idf_cnv_c_to_val(jfrmt,FtxVal->Cval, val,nval);
 }
else if(kfrmt==IDF_FRMT_S)
 {
  l = FtxVal->len;
  k = idf_cnv_s_to_val(jfrmt,FtxVal->Sval,l, val,nval);
 }
else if(kfrmt==IDF_FRMT_T)
 {
  l = FtxVal->len;
  k = idf_cnv_s_to_val(jfrmt,FtxVal->Sval,l, val,nval);
 }
else if(kfrmt==IDF_FRMT_I)
 {
  k = idf_cnv_i_to_val(jfrmt,FtxVal->Ival, val);
 }
else if(kfrmt==IDF_FRMT_L)
 {
  k = idf_cnv_l_to_val(jfrmt,FtxVal->Lval, val);
 }
else if(kfrmt==IDF_FRMT_F)
 {
  f = FtxVal->Fval[0];
  k = idf_cnv_f_to_val(jfrmt,f, val);
 }
else if(kfrmt==IDF_FRMT_D)
 {
  d = FtxVal->Dval[0];
  k = idf_cnv_d_to_val(jfrmt,d, val);
 }
else if(kfrmt==IDF_FRMT_Z)
 {
  k = idf_cnv_z_to_val(jfrmt,FtxVal->Fval, val);
 }
else if(kfrmt==IDF_FRMT_W)
 {
  k = idf_cnv_w_to_val(jfrmt,FtxVal->Dval, val);
 }
else
 {
  k=1;
 }

return flag;
}



#if IDF_CPP_ == 1
int idf_val_lcns(IDF_NAME_L *A)
#else
int idf_val_lcns(A)
IDF_NAME_L *A;
#endif
{
int kfrmt,jfrmt,k;

   kfrmt = FtxVal->jtype;
if(kfrmt==IDF_FRMT_Z || kfrmt==IDF_FRMT_W)
 {
  jfrmt = IDF_FRMT_W;
 }
else
 {
  jfrmt = IDF_FRMT_D;
 }

  A->jest = jfrmt;

  k = idf_val_get(jfrmt, A->val,kfrmt);

return k;
}
