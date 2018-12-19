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
int  idf_cnv_c_to_val(int jfrmt, char cv,
                      void* sval, int nval)
#else
int  idf_cnv_c_to_val(jfrmt,cv, sval,nval)
int   jfrmt;
char  cv;
void *sval;
int   nval;
#endif
{
char *val;
int k=0;

val = (char*)sval;

if(jfrmt==IDF_FRMT_C)
 {
  *val = cv;
 }
else if(jfrmt==IDF_FRMT_S)
 {
  val[0] = cv;
  val[1] = IDF_EOS;
 }
else if(jfrmt==IDF_FRMT_T)
 {
  val[0] = cv;
  val[1] = IDF_EOS;
 }
else
 {
  k=1;
 }


return k;
}



#if IDF_CPP_ == 1
int  idf_cnv_s_to_val(int jfrmt, char* sv, int nv,
                      void* sval, int nval)
#else
int  idf_cnv_s_to_val(jfrmt,sv,nv, sval,nval)
int   jfrmt;
char *sv;
int   nv;
void *sval;
int   nval;
#endif
{
int jus;
int ierr;
 
   jus = idf_jus_get();
if(jus)
 {
  ierr = idf_cnv_s_to_val_uc(jfrmt,sv,nv, sval,nval);
 }
else
 {
  ierr = idf_cnv_s_to_val_c(jfrmt,sv,nv, sval,nval);
 }

return ierr;
}



#if IDF_CPP_ == 1
int  idf_cnv_s_to_val_c(int jfrmt, char* sv, int nv,
                      void* sval, int nval)
#else
int  idf_cnv_s_to_val_c(jfrmt,sv,nv, sval,nval)
int   jfrmt;
char *sv;
int   nv;
void *sval;
int   nval;
#endif
{
char *val;
int k=0;
int n;
register int i;

val = (char*)sval;

if(jfrmt==IDF_FRMT_C)
 {
  *val = sv[0];
 }

else if(jfrmt==IDF_FRMT_S)
 {
  if(nval)
   {
    if(nv)
     {
      if(nv > nval) k=nval;
      else          k=nv;

        i=0;
      while(i<k)
       {
        val[i] = sv[i];
            i++;
       }
      val[i] = IDF_EOS;
      k=0;
     }
    else
     {
      val[0] = IDF_EOS;
     }
   }
  else
   {
    val[0] = IDF_EOS;
   } 
 }

else if(jfrmt==IDF_FRMT_T)
 {
  if(nval)
   {
    if(nv)
     {
      if(nval <= nv)
       {
          i=0;
        while(i<nval)
         {
          val[i] = sv[i];
              i++;
         }
       }
      else
       {
          i=0;
        while(i<nv)
         {
          val[i] = sv[i];
              i++;
         }
#if IDF_TEXT_FULFIL == 1
        while(i<nval)
         {
          val[i] = IDF_TEXT_FULFIL_SYMBOL;
              i++;
         }
#else
        val[i] = IDF_EOS;
#endif
       }
     }
    else
     {
      i=0;
#if IDF_TEXT_FULFIL == 1
      while(i<nval)
       {
        val[i] = IDF_TEXT_FULFIL_SYMBOL;
            i++;
       }
#else
      val[i] = IDF_EOS;
#endif
     }
   }
 }
else
 {
  k=1;
 }


return k;
}



#if IDF_CPP_ == 1
int  idf_cnv_s_to_val_uc(int jfrmt, char* sv, int nv,
                      void* sval, int nval)
#else
int  idf_cnv_s_to_val_uc(jfrmt,sv,nv, sval,nval)
int   jfrmt;
char *sv;
int   nv;
void *sval;
int   nval;
#endif
{
unsigned char *val;
int k=0;
int n;
register int i;

val = (unsigned char*)sval;

if(jfrmt==IDF_FRMT_C)
 {
  *val = (unsigned char) sv[0];
 }

else if(jfrmt==IDF_FRMT_S)
 {
  if(nval)
   {
    if(nv)
     {
      if(nv > nval) k=nval;
      else          k=nv;

        i=0;
      while(i<k)
       {
        val[i] = (unsigned char) sv[i];
            i++;
       }
      val[i] = (unsigned char) IDF_EOS;
      k=0;
     }
    else
     {
      val[0] = (unsigned char) IDF_EOS;
     }
   }
  else
   {
    val[0] = (unsigned char) IDF_EOS;
   } 
 }


else if(jfrmt==IDF_FRMT_T)
 {
  if(nval)
   {
    if(nv)
     {
      if(nval <= nv)
       {
          i=0;
        while(i<nval)
         {
          val[i] = (unsigned char) sv[i];
              i++;
         }
       }
      else
       {
          i=0;
        while(i<nv)
         {
          val[i] = (unsigned char) sv[i];
              i++;
         }
#if IDF_TEXT_FULFIL == 1
        while(i<nval)
         {
          val[i] = (unsigned char) IDF_TEXT_FULFIL_SYMBOL;
              i++;
         }
#else
        val[i] = (unsigned char) IDF_EOS;
#endif
       }
     }
    else
     {
      i=0;
#if IDF_TEXT_FULFIL == 1
      while(i<nval)
       {
        val[i] = (unsigned char) IDF_TEXT_FULFIL_SYMBOL;
            i++;
       }
#else
      val[i] = (unsigned char) IDF_EOS;
#endif
     }
   }
 }
else
 {
  k=1;
 }


return k;
}



#if IDF_CPP_ == 1
int  idf_cnv_i_to_val(int jfrmt, int iv, void* val)
#else
int  idf_cnv_i_to_val(jfrmt,iv, val)
int   jfrmt;
int   iv;
void *val;
#endif
{
int *i;
long *l;
float *f;
double *d;
int k=0;

if(jfrmt==IDF_FRMT_I)
 {
   i = (int*)val;
  *i = iv;  
 }
else if(jfrmt==IDF_FRMT_L)
 {
   l = (long*)val;
  *l = (long)iv;  
 }
else if(jfrmt==IDF_FRMT_F)
 {
   f = (float*)val;
  *f = (float)iv;  
 }
else if(jfrmt==IDF_FRMT_D)
 {
   d = (double*)val;
  *d = (double)iv;  
 }
else if(jfrmt==IDF_FRMT_Z)
 {
  f    = (float*)val;
  f[0] = (float)iv;
  f[1] = (float)0.0;  
 }
else if(jfrmt==IDF_FRMT_W)
 {
  d    = (double*)val;
  d[0] = (double)iv;
  d[1] = (double)0.0;  
 }
else
 {
  k=1;
 }

return k;
}


#if IDF_CPP_ == 1
int  idf_cnv_l_to_val(int jfrmt, long iv, void* val)
#else
int  idf_cnv_l_to_val(jfrmt,iv, val)
int   jfrmt;
long  iv;
void *val;
#endif
{
int *i;
long *l;
float *f;
double *d;
int k=0;
double g;

if(jfrmt==IDF_FRMT_I)
 {
     g = (double)iv;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_INT_MAX)
   {
    k=2;
   }
  else
   {
     i = (int*)val;
    *i = (int)iv;
   }  
 }
else if(jfrmt==IDF_FRMT_L)
 {
   l = (long*)val;
  *l = iv;  
 }
else if(jfrmt==IDF_FRMT_F)
 {
     g = (double)iv;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_FLOAT_MAX)
   {
    k=2;
   }
  else
   {
     f = (float*)val;
    *f = (float)iv;
   }  
 }
else if(jfrmt==IDF_FRMT_D)
 {
   d = (double*)val;
  *d = (double)iv;  
 }
else if(jfrmt==IDF_FRMT_Z)
 {
     g = (double)iv;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_FLOAT_MAX)
   {
    k=2;
   }
  else
   {
    f    = (float*)val;
    f[0] = (float)iv;
    f[1] = (float)0.0;
   }  
 }
else if(jfrmt==IDF_FRMT_W)
 {
  d    = (double*)val;
  d[0] = (double)iv;
  d[1] = (double)0.0;  
 }
else
 {
  k=1;
 }

return k;
}


#if IDF_CPP_ == 1
int  idf_cnv_f_to_val(int jfrmt, float iv, void* val)
#else
int  idf_cnv_f_to_val(jfrmt,iv, val)
int   jfrmt;
float iv;
void *val;
#endif
{
int *i;
long *l;
float *f;
double *d;
int k=0;
double g;

if(jfrmt==IDF_FRMT_I)
 {
     g = (double)iv;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_INT_MAX)
   {
    k=2;
   }
  else
   {
     i = (int*)val;
    *i = (int)iv;
   }  
 }
else if(jfrmt==IDF_FRMT_L)
 {
     g = (double)iv;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_LONG_MAX)
   {
    k=2;
   }
  else
   {
     l = (long*)val;
    *l = (long)iv;
   }  
 }
else if(jfrmt==IDF_FRMT_F)
 {
   f = (float*)val;
  *f = iv;  
 }
else if(jfrmt==IDF_FRMT_D)
 {
   d = (double*)val;
  *d = (double)iv;  
 }
else if(jfrmt==IDF_FRMT_Z)
 {
  f    = (float*)val;
  f[0] = iv;
  f[1] = (float)0.0;  
 }
else if(jfrmt==IDF_FRMT_W)
 {
  d    = (double*)val;
  d[0] = (double)iv;
  d[1] = (double)0.0;  
 }
else
 {
  k=1;
 }

return k;
}


#if IDF_CPP_ == 1
int  idf_cnv_d_to_val(int jfrmt, double iv, void* val)
#else
int  idf_cnv_d_to_val(jfrmt,iv, val)
int    jfrmt;
double iv;
void  *val;
#endif
{
int *i;
long *l;
float *f;
double *d;
int k=0;
double g;

if(jfrmt==IDF_FRMT_I)
 {
     g = (double)iv;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_INT_MAX)
   {
    k=2;
   }
  else
   {
     i = (int*)val;
    *i = (int)iv;
   }  
 }
else if(jfrmt==IDF_FRMT_L)
 {
     g = (double)iv;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_LONG_MAX)
   {
    k=2;
   }
  else
   {
     l = (long*)val;
    *l = (long)iv;
   }  
 }
else if(jfrmt==IDF_FRMT_F)
 {
     g = (double)iv;
  if(g > (double)IDF_FLOAT_MAX)
   {
    k=2;
   }
  else
   {
     f = (float*)val;
    *f = (float)iv;
   }  
 }
else if(jfrmt==IDF_FRMT_D)
 {
   d = (double*)val;
  *d = (double)iv;
 }
else if(jfrmt==IDF_FRMT_Z)
 {
     g = (double)iv;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_FLOAT_MAX)
   {
    k=2;
   }
  else
   {
    f    = (float*)val;
    f[0] = (float)iv;
    f[1] = (float)0.0;
   }  
 }
else if(jfrmt==IDF_FRMT_W)
 {
  d    = (double*)val;
  d[0] = (double)iv;
  d[1] = (double)0.0;  
 }
else
 {
  k=1;
 }

return k;
}


#if IDF_CPP_ == 1
int  idf_cnv_z_to_val(int jfrmt, float* iv, void* val)
#else
int  idf_cnv_z_to_val(jfrmt,iv, val)
int    jfrmt;
float *iv;
void  *val;
#endif
{
float ff;
int *i;
long *l;
float *f;
double *d;
int k=0;
double g;

ff = iv[0];

if(jfrmt==IDF_FRMT_I)
 {
     g = (double)ff;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_INT_MAX)
   {
    k=2;
   }
  else
   {
     i = (int*)val;
    *i = (int)ff;
   }  
 }
else if(jfrmt==IDF_FRMT_L)
 {
     g = (double)ff;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_LONG_MAX)
   {
    k=2;
   }
  else
   {
     l = (long*)val;
    *l = (long)ff;
   }  
 }
else if(jfrmt==IDF_FRMT_F)
 {
   f = (float*)val;
  *f = ff;  
 }
else if(jfrmt==IDF_FRMT_D)
 {
   d = (double*)val;
  *d = (double)ff;  
 }
else if(jfrmt==IDF_FRMT_Z)
 {
  f    = (float*)val;
  f[0] = ff;
  f[1] = iv[1];  
 }
else if(jfrmt==IDF_FRMT_W)
 {
  d    = (double*)val;
  d[0] = (double)ff;
    ff = iv[1];
  d[1] = (double)ff;  
 }
else
 {
  k=1;
 }

return k;
}


#if IDF_CPP_ == 1
int  idf_cnv_w_to_val(int jfrmt, double* iv, void* val)
#else
int  idf_cnv_w_to_val(jfrmt,iv, val)
int     jfrmt;
double *iv;
void   *val;
#endif
{
int *i;
long *l;
float *f;
double *d;
int k=0;
double g,gg;

gg = iv[0];

if(jfrmt==IDF_FRMT_I)
 {
     g = (double)gg;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_INT_MAX)
   {
    k=2;
   }
  else
   {
     i = (int*)val;
    *i = (int)gg;
   }  
 }
else if(jfrmt==IDF_FRMT_L)
 {
     g = (double)gg;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_LONG_MAX)
   {
    k=2;
   }
  else
   {
     l = (long*)val;
    *l = (long)gg;
   }  
 }
else if(jfrmt==IDF_FRMT_F)
 {
     g = (double)gg;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_FLOAT_MAX)
   {
    k=2;
   }
  else
   {
     f = (float*)val;
    *f = (float)gg;
   }  
 }
else if(jfrmt==IDF_FRMT_D)
 {
   d = (double*)val;
  *d = (double)gg;  
 }
else if(jfrmt==IDF_FRMT_Z)
 {
     g = (double)gg;
  if(g < 0.0) g = -g;
  if(g > (double)IDF_FLOAT_MAX)
   {
    k=2;
   }
  else
   {
    f    = (float*)val;
    f[0] = (float)gg;

      gg = iv[1];
       g = (double)gg;
    if(g < 0.0) g = -g;
    if(g > (double)IDF_FLOAT_MAX)
     {
      k=2;
     }
    else
     {
      f[1] = (float)gg;
     }    
   }  
 }
else if(jfrmt==IDF_FRMT_W)
 {
  d    = (double*)val;
  d[0] = iv[0];
  d[1] = iv[1];  
 }
else
 {
  k=1;
 }

return k;
}
