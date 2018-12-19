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


/* sizeof pointer in bytes */
#define IDF_POINTER_SIZE sizeof(void*)


static int  idf_FrmtSize_List   [IDF_FRMT_NTYP_1];
static int  idf_FrmtWord_List   [IDF_FRMT_NTYP_1];
static int  idf_FrmtPointer_List[IDF_FRMT_NTYP_1];



#if IDF_CPP_ == 1
int idf_frmt_type(unsigned int cc)
#else
int idf_frmt_type(cc)
unsigned int cc;
#endif
{
static char FrmtType_List[IDF_FRMT_NTYP_1] = {
                                              IDF_FRMT_LIST
                                             };
register int i;
char c;
int k;

c = tolower(cc);

k=0;
for(i=0;i<IDF_FRMT_NTYP;i++)
 {
  if(c == FrmtType_List[i])
   {
    k = i+1;
    break;
   }
 }

 return k;
}



#if IDF_CPP_ == 1
int idf_frmt_ctype(char cc)
#else
int idf_frmt_ctype(cc)
char cc;
#endif
{
char c;
int k;

c = tolower(cc);

     if(c == 'c') k=IDF_FRMT_C;
else if(c == 's') k=IDF_FRMT_S;
else if(c == 't') k=IDF_FRMT_T;
else if(c == 'i') k=IDF_FRMT_I;
else if(c == 'l') k=IDF_FRMT_L;
else if(c == 'f') k=IDF_FRMT_F;
else if(c == 'd') k=IDF_FRMT_D;
else if(c == 'z') k=IDF_FRMT_Z;
else if(c == 'w') k=IDF_FRMT_W;
else              k=IDF_FRMT_UNKNOWN;

 return k;
}



int idf_frmt_table_prp()
{
static char nfu[] = "idf_frmt_table_prp";

static char FrmtType_List_Orig[] = IDF_FRMT_LIST ;

int ierr=0;
int Nf;
int kp,kk;
char c;
register int j,k;

   Nf = strlen(FrmtType_List_Orig);
if(Nf<1 || Nf>IDF_FRMT_NTYP)
 {
  fprintf(stdout,"%s improper size=%d of format type table\n",
          nfu,Nf);
  ierr=1;
  goto err;
 }


for(j=0;j<Nf;j++)
 {
  c = FrmtType_List_Orig[j];

     if(c == 'c') k=IDF_FRMT_C;
else if(c == 's') k=IDF_FRMT_S;
else if(c == 't') k=IDF_FRMT_T;
else if(c == 'i') k=IDF_FRMT_I;
else if(c == 'l') k=IDF_FRMT_L;
else if(c == 'f') k=IDF_FRMT_F;
else if(c == 'd') k=IDF_FRMT_D;
else if(c == 'z') k=IDF_FRMT_Z;
else if(c == 'w') k=IDF_FRMT_W;
else              k=IDF_FRMT_UNKNOWN;

     kk = j+1;
  if(kk != k)
   {
    ierr=2;
    break;
   }
 }

if(ierr)
 {
  fprintf(stdout,"%s inconsistency in format type=%c and table list=%d\n",
          nfu,c,k);
  goto err;
 }


for(j=0;j<IDF_FRMT_NTYP;j++)
 {
  c = FrmtType_List_Orig[j];

  if(c == 'c') 
   {
    kk=sizeof(char);
     k=kk;
    kp=sizeof(char*);
   }
  else if(c == 's')
   {
    kk=sizeof(char);
     k=kk;
    kp=sizeof(char*);
   }
  else if(c == 't')
   {
    kk=sizeof(char); 
     k=kk;
    kp=sizeof(char*);
   }
  else if(c == 'i')
   {
    kk=sizeof(int);
     k=kk;
    kp=sizeof(int*);
   }
  else if(c == 'l')
   {
    kk=sizeof(long);
     k=kk;
    kp=sizeof(long*);
   }
  else if(c == 'f')
   {
    kk=sizeof(float);
     k=kk;
    kp=sizeof(float*);
   }
  else if(c == 'd')
   {
    kk=sizeof(double);
     k=kk;
    kp=sizeof(double*);
   }
  else if(c == 'z')
   {
    kk=sizeof(float);
     k=2*kk;
    kp=sizeof(float*);
   }
  else if(c == 'w')
   {
    kk=sizeof(double);
     k=2*kk;
    kp=sizeof(double*);
   }
  else
   {
    kk=0;
     k=kk;
    kp=sizeof(void*);
   }

  idf_FrmtSize_List   [j] = k;
  idf_FrmtWord_List   [j] = kk;
  idf_FrmtPointer_List[j] = kp;
 }

k=0;
j=0;
while(j<IDF_FRMT_NTYP)
 {
  if(idf_FrmtSize_List[j] < 1)
   {
    k=1;
    break;
   }

  if(idf_FrmtWord_List[j] < 1)
   {
    k=1;
    break;
   }

  if(idf_FrmtPointer_List[j] < 1)
   {
    k=1;
    break;
   }

  j++;
 }

if(k)
 {
  fprintf(stdout,"%s improper size for format type=%d in table list\n",
          nfu,j);
  ierr=3;
  goto err;
 }
else
 {
  j = IDF_FRMT_NTYP;

  idf_FrmtSize_List   [j] = 0;
  idf_FrmtWord_List   [j] = 0;
  idf_FrmtPointer_List[j] = sizeof(void*);
 }

/*
k = IDF_POINTER_SIZE;
j = 0;
   kk = sizeof(char*);
if(kk != k) j=kk;
   kk = sizeof(int*);
if(kk != k) j=kk;
   kk = sizeof(long*);
if(kk != k) j=kk;
   kk = sizeof(float*);
if(kk != k) j=kk;
   kk = sizeof(double*);
if(kk != k) j=kk;

if(j)
 {
  fprintf(stdout,"%s improper size for pointer=%d %d\n",
          nfu,k,j);
  ierr=4;
 }
*/

err:
return ierr;
}



#if IDF_CPP_ == 1
int idf_frmt_size(int jf)
#else
int idf_frmt_size(jf)
int jf;
#endif
{
register int j;

     if(jf <               1) j = IDF_FRMT_NTYP_1;
else if(jf > IDF_FRMT_NTYP_1) j = IDF_FRMT_NTYP_1;
else                          j = jf;
                              j--;

return ( idf_FrmtSize_List[j] );
}


#if IDF_CPP_ == 1
int idf_frmt_word(int jf)
#else
int idf_frmt_word(jf)
int jf;
#endif
{
register int j;

     if(jf <               1) j = IDF_FRMT_NTYP_1;
else if(jf > IDF_FRMT_NTYP_1) j = IDF_FRMT_NTYP_1;
else                          j = jf;
                              j--;

return ( idf_FrmtWord_List[j] );
}


#if IDF_CPP_ == 1
int idf_frmt_psize(int jf)
#else
int idf_frmt_psize(jf)
int jf;
#endif
{
register int j;

     if(jf <               1) j = IDF_FRMT_NTYP_1;
else if(jf > IDF_FRMT_NTYP_1) j = IDF_FRMT_NTYP_1;
else                          j = jf;
                              j--;

return ( idf_FrmtPointer_List[j] );
}


#if IDF_CPP_ == 1
int idf_frmt_cns(int jf)
#else
int idf_frmt_cns(jf)
int jf;
#endif
{
/* format for local constants */
static int  FrmtCns[IDF_FRMT_NTYP_1] = {
        /* c */ IDF_FRMT_UNKNOWN ,
        /* s */ IDF_FRMT_UNKNOWN ,
        /* t */ IDF_FRMT_UNKNOWN ,
        /* i */ IDF_FRMT_D       ,
        /* l */ IDF_FRMT_D       ,
        /* f */ IDF_FRMT_D       ,
        /* d */ IDF_FRMT_D       ,
        /* z */ IDF_FRMT_W       ,
        /* w */ IDF_FRMT_W       ,
        /* v */ IDF_FRMT_VOID 
                                       };
register int j;

     if(jf <               1) j = IDF_FRMT_NTYP_1;
else if(jf > IDF_FRMT_NTYP_1) j = IDF_FRMT_NTYP_1;
else                          j = jf;
                              j--;

return ( FrmtCns[j] );
}


#if IDF_CPP_ == 1
int idf_frmt_amode(int jf)
#else
int idf_frmt_amode(jf)
int jf;
#endif
{
/* alignment mode flag: */
static int  FrmtA[IDF_FRMT_NTYP_1] = {
        /* c */ IDF_AMODE_TEXT   ,
        /* s */ IDF_AMODE_TEXT   ,
        /* t */ IDF_AMODE_TEXT   ,
        /* i */ IDF_AMODE_BINARI ,
        /* l */ IDF_AMODE_BINARI ,
        /* f */ IDF_AMODE_BINARI ,
        /* d */ IDF_AMODE_BINARI ,
        /* z */ IDF_AMODE_BINARI ,
        /* w */ IDF_AMODE_BINARI ,
        /* v */ IDF_AMODE_BINARI 
                                       };
register int j;

     if(jf <               1) j = IDF_FRMT_NTYP_1;
else if(jf > IDF_FRMT_NTYP_1) j = IDF_FRMT_NTYP_1;
else                          j = jf;
                              j--;

return ( FrmtA[j] );
}
