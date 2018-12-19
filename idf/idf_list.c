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




/* sizeof IDF_NAME_ */
static int NameList_SIZE = 0;

/* sizeof IDF_NAME_L */
static int NameList_SIZE_L = 0;

/* list of names */
static IDF_NAME_LIST NameList;




IDF_NAME_LIST* idf_list_address()
{
 return ( &NameList );
}



void idf_list_nul()
{
register int i;

for(i=0;i<IDF_NAME_LIST_N;i++)
 {
  NameList.List[i] = NULL;
 }

for(i=0;i<IDF_NAME_LISTL_N;i++)
 {
  NameList.ListL[i] = NULL;
 }

NameList.Num  = 0;
NameList.NumL = 0;

NameList.kerr = 0;
NameList.lvl  = 0;
NameList.irec = 0;
}



int idf_list_ini()
{
int ierr=0;
register int i;
IDF_NAME_  *A;
IDF_NAME_L *AL;


NameList_SIZE   = sizeof(IDF_NAME_ );
NameList_SIZE_L = sizeof(IDF_NAME_L);


for(i=0;i<IDF_NAME_LIST_N;i++)
 {
     A = NameList.List[i];
  if(A != NULL)
   {
    /*list structure No=%d has been already init*/
    NameList.kerr = 1;
    NameList.lvl  = 0;
    NameList.irec = i;
    ierr=1;
    break;
   }

     A = (IDF_NAME_*) malloc(NameList_SIZE);
  if(A == NULL)
   {
    /*not enough memory to init name list*/
    NameList.kerr = 2;
    NameList.lvl  = 0;
    NameList.irec = i;
    ierr=2;
    break;
   }

  A->name   = NULL;
  A->length = 0;
  A->kstr   = 0;
  A->kurs   = 0;
  A->dpos   = 0;
  A->nd     = 0;
  A->jtype  = 0;
  A->ndim   = 0;
  A->Mshft  = NULL;

  NameList.List[i] = A;
 }

NameList.Num = 0;

if(ierr) goto err;



for(i=0;i<IDF_NAME_LISTL_N;i++)
 {
     AL = NameList.ListL[i];
  if(AL != NULL)
   {
    /*list structure No=%d has been already init*/
    NameList.kerr = 1;
    NameList.lvl  = 0;
    NameList.irec = i;
    ierr=1;
    break;
   }

     AL = (IDF_NAME_L*) malloc(NameList_SIZE_L);
  if(AL == NULL)
   {
    /*not enough memory to init name list*/
    NameList.kerr = 2;
    NameList.lvl  = 0;
    NameList.irec = i;
    ierr=2;
    break;
   }

  AL->name   = NULL;
  AL->length = 0;
  AL->kstr   = 0;
  AL->kurs   = 0;
  AL->dpos   = 0;
  AL->nd     = 0;
  AL->jtype  = 0;
  AL->jest   = 0;
  AL->val[0] = (double)0.0;
  AL->val[1] = (double)0.0;

  NameList.ListL[i] = AL;
 }

NameList.NumL = 0;

err:
return ierr;
}


void idf_list_clean()
{
register int i;
IDF_NAME_ *A;
IDF_NAME_L *AL;

for(i=0;i<IDF_NAME_LIST_N;i++)
 {
     A = NameList.List[i];
  if(A != NULL)
   {
    if(A->name != NULL)
     {
      free(A->name);
     }

    if(A->Mshft != NULL)
     {
      free(A->Mshft);
     }

    A->name   = NULL;
    A->length = 0;
    A->kstr   = 0;
    A->kurs   = 0;
    A->dpos   = 0;
    A->nd     = 0;
    A->jtype  = 0;
    A->ndim   = 0;
    A->Mshft  = NULL;
   }

 }

 NameList.Num = 0;


for(i=0;i<IDF_NAME_LISTL_N;i++)
 {
     AL = NameList.ListL[i];
  if(AL != NULL)
   {
    if(AL->name != NULL)
     {
      free(AL->name);
     }

    AL->name   = NULL;
    AL->length = 0;
    AL->kstr   = 0;
    AL->kurs   = 0;
    AL->dpos   = 0;
    AL->nd     = 0;
    AL->jtype  = 0;
    AL->jest   = 0;
    AL->val[0] = (double)0.0;
    AL->val[1] = (double)0.0;
   }

 }

 NameList.NumL   = 0;

}


void idf_list_end()
{
register int i;
IDF_NAME_ *A;
IDF_NAME_L *AL;

for(i=0;i<IDF_NAME_LIST_N;i++)
 {
     A = NameList.List[i];
  if(A != NULL)
   {
    if(A->name != NULL)
     {
      free(A->name);
      A->name = NULL;
     }

    if(A->Mshft != NULL)
     {
      free(A->Mshft);
      A->Mshft = NULL;
     }

    free(A);

    NameList.List[i] = NULL;
   }

 }

 NameList_SIZE = 0;
 NameList.Num = 0;


for(i=0;i<IDF_NAME_LISTL_N;i++)
 {
     AL = NameList.ListL[i];
  if(AL != NULL)
   {
    if(AL->name != NULL)
     {
      free(AL->name);
      AL->name = NULL;
     }

    free(AL);

    NameList.ListL[i] = NULL;
   }

 }

 NameList_SIZE_L = 0;
 NameList.NumL   = 0;

}



int idf_list_num()
{
 return ( NameList.Num );
}

int idf_list_numL()
{
 return ( NameList.NumL );
}


#if IDF_CPP_ == 1
IDF_NAME_* idf_list_get(char* s, int n, int irec)
#else
IDF_NAME_* idf_list_get(s,n,irec)
char *s;
int   n;
int   irec;
#endif
{
int k;
int nn;
char *ss;
register int j;
IDF_NAME_ *A;

A = NULL;

     if(NameList.Num <  1) goto fin;
else if(irec         <  0) goto fin;
else if(irec         >= IDF_NAME_LIST_N) goto fin;

A = NameList.List[irec];


    ss = A->name;
    nn = A->length;

    if(n == nn)
     {
      k=1;
      for(j=0;j<n;j++)
       {
        if(s[j] != ss[j])
         {
          k=0;
          break;
         }
       }
     }

    else
     {
      k=0;
     }


  if(k==0)
   {
    A = NULL;
   }

fin:
return A;
}


#if IDF_CPP_ == 1
IDF_NAME_L* idf_list_getL(char* s, int n, int irec)
#else
IDF_NAME_L* idf_list_getL(s,n,irec)
char *s;
int   n;
int   irec;
#endif
{
int k;
int nn;
char *ss;
register int j;
IDF_NAME_L *A;

A = NULL;

     if(NameList.NumL <  1) goto fin;
else if(irec          <  0) goto fin;
else if(irec          >= IDF_NAME_LISTL_N) goto fin;

A = NameList.ListL[irec];


    ss = A->name;
    nn = A->length;

    if(n == nn)
     {
      k=1;
      for(j=0;j<n;j++)
       {
        if(s[j] != ss[j])
         {
          k=0;
          break;
         }
       }
     }

    else
     {
      k=0;
     }


  if(k==0)
   {
    A = NULL;
   }

fin:
return A;
}


#if IDF_CPP_ == 1
int idf_list_getLV(char* s, int n,
                   int* vtyp, double* val)
#else
int idf_list_getLV(s,n, vtyp,val)
char   *s;
int     n;
int    *vtyp;
double *val;
#endif
{
/*
  2 no value 
  1 no name
  0 got it
 <0 error
*/

int k,flag;
int nn,nl;
char *ss;
register int i,j;
IDF_NAME_L *A;

A = NULL;

   nl = idf_list_numL();
if(nl == 0)
 {
  flag=1;
 }

else if(n<1)
 {
  flag = -1;
 }

else
 {
  flag = 1;
  for(i=0;i<nl;i++)
   {
     A = NameList.ListL[i];
    ss = A->name;
    nn = A->length;

    if(n == nn)
     {
      ss = A->name;
      nn = A->length;

      k=1;
      if(n == nn)
       {
        for(j=0;j<n;j++)
         {
          if(s[j] != ss[j])
           {
            k=0;
            break;
           }
         }
       }
      if(k)
       {
        val[0] = A->val[0];
        val[1] = A->val[1];
         *vtyp = A->jtype;

        if(A->jest)
         {
          flag=0;
         }
        else
         {
          flag = 2;
         }
       }
     }
   }/*i*/
 }

return flag;
}



#if IDF_CPP_ == 1
int idf_list_put(IDF_NAME* Name)
#else
int idf_list_put(Name)
IDF_NAME *Name;
#endif
{
int ierr=0;
char *s,*ss;
int n;
int ndim;
int *mdim;
int irec,k;
register int i;
IDF_NAME_ *A;

   A = NULL;
  ss = NULL;
mdim = NULL;

   irec = NameList.Num;
if(irec < IDF_NAME_LIST_N)
 {
  s = Name->name;
  n = Name->Nlength;

  if(n < 1)
   {
    /*improper name,unable to make a record*/
    NameList.kerr = 3;
    NameList.lvl  = 1;
    NameList.irec = irec;
    ierr=2;
    goto err;
   }

   k = Name->jlocal;
if(k)
 {
  /*unable to put local name to nonlocal list */
  NameList.kerr = 4;
  NameList.lvl  = 1;
  NameList.irec = irec;
  ierr=2;
  goto err;
 }

  A = NameList.List[irec];

      k = (n+1)*sizeof(char);
     ss = (char*) malloc(k);
  if(ss == NULL)
   {
    /* no more memory*/
    NameList.kerr = 5;
    NameList.lvl  = 0;
    NameList.irec = irec;
    ierr=3;
    goto err;
   }

    i=0;
  while(i<n)
   {
    ss[i] = s[i];
    i++;
   }
  ss[i] = IDF_EOS;

  A->name   = ss;
  A->length = n;
  A->kstr   = Name->le;
  A->kurs   = Name->ke;
  A->dpos   = Name->dspos;
  A->nd     = Name->depos - Name->dspos;  
  A->jtype  = Name->jtype;

     ndim = Name->ndim;
  A->ndim = ndim;

  if(ndim > 0)
   {
          k = ndim*sizeof(int);
       mdim = (int*) malloc(k);
    if(mdim == NULL)
     {
      /* no more memory*/
      NameList.kerr = 5;
      NameList.lvl  = 0;
      NameList.irec = irec;
      ierr=3;
      goto err;
     }

    for(i=0;i<ndim;i++)
     {
      mdim[i] = Name->Mshft[i];
     }

    A->Mshft = mdim;
   }

  NameList.Num = irec + 1;
 }

else
 {
  /* non-local list is full */
  NameList.kerr = 6;
  NameList.lvl  = 0;
  NameList.irec = irec;
  ierr=1;
 }

err:
return ierr;
}



#if IDF_CPP_ == 1
int idf_list_putL(IDF_NAME* Name)
#else
int idf_list_putL(Name)
IDF_NAME *Name;
#endif
{
int ierr=0;
char *s,*ss;
int n,nn;
int irec,jrec;
int k,jest;
register int i,j;
IDF_NAME_L *A;

   A = NULL;
  ss = NULL;


   irec = NameList.NumL;
if(irec >= IDF_NAME_LISTL_N)
 {
  /* local list is full */
  NameList.kerr = 9;
  NameList.lvl  = 0;
  NameList.irec = irec;
  ierr=3;
  goto err;
 }

   s = Name->name;
   n = Name->Nlength;

   k = Name->jlocal;
if(k==0)
 {
  /*unable to put nonlocal name to local list */
  NameList.kerr = 7;
  NameList.lvl  = 1;
  NameList.irec = irec;
  ierr=1;
  goto err;
 }

if(n < 1)
 {
  /*improper name,unable to make a record*/
  NameList.kerr = 3;
  NameList.lvl  = 1;
  NameList.irec = irec;
  ierr=2;
  goto err;
 }


 jest = 0;

if(irec)
 {
  jrec = 0;
  for(i=0;i<irec;i++)
   {
    A = NameList.ListL[i];

    ss = A->name;
    nn = A->length;

    k = 1;

    if(n != nn)
     {
      k=0;
     }
    else
     {
      for(j=0;j<n;j++)
       {
        if(ss[j] != s[j])
         {
          k=0;
          break;
         }
       }
     }

    if(k)
     {
      /* names coinside */
      jest = 1;
      jrec = i;
      break;
     }

   }/* seach list */

 }


if(jest)
 {
  /*-------------- 
    update record
   ---------------*/

  A->kstr   = Name->le;
  A->kurs   = Name->ke;
  A->dpos   = Name->dspos;
  A->nd     = Name->depos - Name->dspos;  
  A->jtype  = Name->jtype;

  A->val[0]  = 0.0;
  A->val[1]  = 0.0;
 }

else
 {
  /*----------- 
    add record
   ------------*/

  A = NameList.ListL[irec];

      k = (n+1)*sizeof(char);
     ss = (char*) malloc(k);
  if(ss == NULL)
   {
    /* no more memory*/
    NameList.kerr = 8;
    NameList.lvl  = 0;
    NameList.irec = irec;
    ierr=4;
    goto err;
   }

    i=0;
  while(i<n)
   {
    ss[i] = s[i];
    i++;
   }
  ss[i] = IDF_EOS;

  A->name   = ss;
  A->length = n;
  A->kstr   = Name->le;
  A->kurs   = Name->ke;
  A->dpos   = Name->dspos;
  A->nd     = Name->depos - Name->dspos;  
  A->jtype  = Name->jtype;

  A->val[0] = 0.0;
  A->val[1] = 0.0;

  NameList.NumL = irec + 1;
 }


err:
return ierr;
}



#if IDF_CPP_ == 1
void idf_list_prn_(IDF_NAME_* A)
#else
void idf_list_prn_(A)
IDF_NAME_ *A;
#endif
{
 register int i;

 if(A->length > 0)
  {
   fprintf(stdout,"LIST Name(%3d)=%s\n",
           A->length,A->name);

   if(A->ndim > 0)
    {
     fprintf(stdout,"Mdim(%1d)=",A->ndim);
     for(i=0;i<A->ndim;i++)
      {
       fprintf(stdout," %5d",A->Mshft[i]);
      }
     fprintf(stdout,"\n");
    }
   else
    {
     fprintf(stdout,"ndim=%d\n",A->ndim);
    }

   fprintf(stdout,"jtype=%2d\n",A->jtype);
   fprintf(stdout,"Data position=%ld string=%d kurs=%d length=%d\n", 
           A->dpos,A->kstr,A->kurs,A->nd);
  }
}



#if IDF_CPP_ == 1
void idf_list_prnL_(IDF_NAME_L* A)
#else
void idf_list_prnL_(A)
IDF_NAME_L *A;
#endif
{
 register int i;

 if(A->length > 0)
  {
   fprintf(stdout,"LIST LOCAL Name(%3d)=%s\n",
           A->length,A->name);

   fprintf(stdout,"jtype=%2d jest=%1d\n",
           A->jtype,A->jest);
   fprintf(stdout,"Data position=%ld string=%d kurs=%d length=%d\n", 
           A->dpos,A->kstr,A->kurs,A->nd);
   if(A->jest)
    {
     fprintf(stdout,"val=%24.16e %24.16e\n",
             A->val[0],A->val[1]);
    }
  }
}



#if IDF_CPP_ == 1
void idf_list_prn(int irec)
#else
void idf_list_prn(irec)
int irec;
#endif
{
int k;
IDF_NAME_ *A;

A = NULL;
                           k=1;
     if(NameList.Num <  1) k=0;
else if(irec         <  0) k=0;
else if(irec         >= IDF_NAME_LIST_N) k=0;
A = NameList.List[irec];

if(k)
 {
  idf_list_prn_(A);
 }

}



#if IDF_CPP_ == 1
void idf_list_prnL(int irec)
#else
void idf_list_prnL(irec)
int irec;
#endif
{
int k;
IDF_NAME_L *A;

A = NULL;
                            k=1;
     if(NameList.NumL <  1) k=0;
else if(irec          <  0) k=0;
else if(irec          >= IDF_NAME_LISTL_N) k=0;
A = NameList.ListL[irec];

if(k)
 {
  idf_list_prnL_(A);
 }

}


int idf_list_do()
{
#if IDF_NAME_TRACE == 1
static char nfu[] = "idf_list_do:";
#endif

int ierr=0;
int n,nl;
int flag,k,jerr,jlocal,ibrk;
int kst,kurs;
long pos;
register int i,j;
IDF_NAME *Name;

     n = IDF_NAME_LIST_N  + 1;
    nl = IDF_NAME_LISTL_N + 1;
  flag = 0;
jlocal = 0;

 kst = 0;
kurs = 0;
 pos = 0L;

Name = idf_name_address();

      i=0;
      j=0;

   ibrk=0;
do
 {

#if IDF_NAME_TRACE == 1
    fprintf(stdout,"%s record No=%d\n",nfu,i);
#endif

   /*-------------------------
     get named data structure
    -------------------------- */

     k = idf_name_get(&kst,&kurs,&pos, &jlocal);

#if IDF_NAME_TRACE == 1
     idf_name_prn();
#endif

  if(k == -1)
   {
#if IDF_NAME_TRACE == 1
      fprintf(stdout,"%s named data with eof-flag\n",nfu);
#endif

    if(jlocal)
     {
      /* add record to local name list */
         jerr = idf_list_putL(Name);
      if(jerr)
       {
        flag = -1;
       }
      else
       {
#if IDF_NAME_TRACE == 1
         idf_list_prnL(j);
#endif
        flag = 1;
       }
      j++;
     }
    else
     {
      /* add record to non-local name list */
         jerr = idf_list_put(Name);
      if(jerr)
       {
        flag = -1;
       }
      else
       {
#if IDF_NAME_TRACE == 1
          idf_list_prn(i);
#endif
        flag = 1;
       }
      i++;
     }
   }

  if(k == -2)
   {
#if IDF_NAME_TRACE == 1
     fprintf(stdout,"%s empty name with eof-flag\n",nfu);
#endif
    flag = 1;
   }

  else if(k)
   {
    /* error in name search */
    NameList.kerr = 10;
    NameList.lvl  = 2;
    NameList.irec = -1;
/*
    idf_name_errprn();
    idf_txt_err_prn();
*/
    flag = -1;
   }

  else
   {
#if IDF_NAME_TRACE == 1
      fprintf(stdout,"%s standard named data\n",nfu);
#endif

    if(jlocal)
     {
      /* add record to local name list */
         jerr = idf_list_putL(Name);
      if(jerr)
       {
        flag = -1;
       }
      else
       {
#if IDF_NAME_TRACE == 1
          idf_list_prnL(j);
#endif
       }
      j++;
     }
    else
     {
      /* add record to non-local name list */
         jerr = idf_list_put(Name);
      if(jerr)
       {
        flag = -1;
       }
      else
       {
#if IDF_NAME_TRACE == 1
          idf_list_prn(i);
#endif
       }
      i++;
     }
   }


  if(flag)
   {
    ibrk=1;
   }
  else if((i > n)||(j>nl))
   {
    /*unexpected lists overflow */
    NameList.kerr = 11;
    NameList.lvl  = 0;

    if(jlocal)
     {
      NameList.irec = j;
     }
    else
     {
      NameList.irec = i;
     }

    flag = -1;
    ibrk=2;
   }

 } while(ibrk==0);


if(flag < 0)
 {
  ierr=1;
 }

return ierr;
}



int idf_list_check()
{
int ierr=0;
int n;
IDF_NAME_ *A;
int nl;
IDF_NAME_L *AL;
int l,ll;
char *s,*ss;
register int i,j,k;
int jest,m;
double val;

   n = idf_list_num();
if(n < 1)
 {
  /*empty non-local list*/
  NameList.kerr = 12;
  NameList.lvl  = 0;
  NameList.irec = -1;

  ierr=1;
  goto err;
 }

else
 {
  jest=0;
     i=0;
  while(i<n)
   {
    A = NameList.List[i];
    s = A->name;
    l = A->length;

       jest = idf_clist(s,l, &val);
    if(jest) break;

    i++;
   }

  if(jest)
   {
    /*data name redefines the standard name*/
    NameList.kerr = 13;
    NameList.lvl  = 8;
    NameList.irec = i;

    ierr=2;
    goto err;
   }


  jest=0;
     i=0;
  while(i<n)
   {
    A = NameList.List[i];
    s = A->name;
    l = A->length;

       jest = idf_keywF(s,l);
    if(jest) break;

    i++;
   }

  if(jest)
   {
    /*data name redefines the standard function name*/
    NameList.kerr = 14;
    NameList.lvl  = 8;
    NameList.irec = i;

    ierr=6;
    goto err;
   }
 }


   nl = idf_list_numL();
if(nl)
 {
  jest=0;

    i=0;
  while(i<n)
   {
    A = NameList.List[i];
    s = A->name;
    l = A->length;

      j=0;
    while(j<nl)
     {
      AL = NameList.ListL[j];
      ss = AL->name;
      ll = AL->length;

      m = 1;

      if(l != ll)
       {
        m=0;
       }
      else
       {
        for(k=0;k<n;k++)
         {
          if(ss[k] != s[k])
           {
            m=0;
            break;
           }
         }
       }

      if(m)
       {
        /* records name coinside */
        jest = j+1;
        break;
       }

      j++;
     } /* loc list */

    if(jest) break;

    i++;
   } /* nonloc list */

  if(jest)
   {
    /*data name is in non-local and local lists*/
    NameList.kerr = 15;
    NameList.lvl  = 7;
    NameList.irec = jest-1;

    ierr=3;
    goto err;
   }


  jest=0;
     i=0;
  while(i<nl)
   {
    AL = NameList.ListL[i];
    ss = AL->name;
    ll = AL->length;

       jest = idf_clist(ss,ll, &val);
    if(jest) break;

    i++;
   }

  if(jest)
   {
    /*data name redefines the standard name*/
    NameList.kerr = 13;
    NameList.lvl  = 7;
    NameList.irec = i;

    ierr=4;
    goto err;
   }


  jest=0;
     i=0;
  while(i<nl)
   {
    AL = NameList.ListL[i];
    ss = AL->name;
    ll = AL->length;

       jest = idf_keywF(ss,ll);
    if(jest) break;

    i++;
   }

  if(jest)
   {
    /*data name redefines the standard function name*/
    NameList.kerr = 14;
    NameList.lvl  = 7;
    NameList.irec = i;

    ierr=5;
    goto err;
   }

 }

err:
return ierr;
}


void idf_list_catalog()
{
register int i;
int n;

   n = NameList.Num;
if(n)
 {
  for(i=0;i<n;i++)
   {
    idf_list_prn(i);
   }
 }
}


void idf_list_catalogL()
{
register int i;
int n;

   n = NameList.NumL;
if(n)
 {
  for(i=0;i<n;i++)
   {
    idf_list_prnL(i);
   }
 }
}


void idf_list_clearerr()
{
 NameList.kerr = 0;
 NameList.lvl  = 0;
 NameList.irec = 0;
}
