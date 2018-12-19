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




IDF_FRMT_TREAT FrmtS;




IDF_FRMT_TREAT* idf_frmtreat_address()
{
 return ( &FrmtS );
}


void idf_frmtreat_nul()
{
 register int i;

 FrmtS.ielem   = -1;
 FrmtS.iunit   =  0;
 FrmtS.irep    =  0;
 FrmtS.jrep    =  0;
 FrmtS.offset  =  0;
 FrmtS.pointer =  0;
 FrmtS.jalloc  =  0;
 FrmtS.jmem    =  0;
 FrmtS.jtype   =  0;
 FrmtS.jsize   =  0;
 FrmtS.jwrd    =  0;
 FrmtS.swrd    =  0;
 FrmtS.amode   =  0;
 FrmtS.align   =  0;
 FrmtS.blk     =  0;
 FrmtS.jerr    =  0;

 FrmtS.Ndimen  =  0;
 for(i=0;i<IDF_FRMT_DIM_MAX;i++)
  {
   FrmtS.Mdimen[i] = 0;
  }

 FrmtS.Otarget =  NULL;
 FrmtS.Atarget =  NULL;
 FrmtS.Target  =  NULL;
}


void idf_frmtreat_prn()
{
 int n;
 register int i;

 fprintf(stdout,"Format: elem=%d unit=%d\n",
         FrmtS.ielem,FrmtS.iunit);
 fprintf(stdout,"irept=%d jrept=%d\n",
         FrmtS.irep,FrmtS.jrep);
 fprintf(stdout,"offset=%d gap=%d\n",
         FrmtS.offset,FrmtS.gap);
 fprintf(stdout,"pointer=%d jalloc=%d\n",
         FrmtS.pointer,FrmtS.jalloc);
 fprintf(stdout,"jmem=%d jtype=%d amode=%d\n",
         FrmtS.jmem,FrmtS.jtype,FrmtS.amode);
 fprintf(stdout,"jsize=%d wsize=%d swrd=%d\n",
         FrmtS.jsize,FrmtS.jwrd,FrmtS.swrd);
 fprintf(stdout,"align=%d blk=%d jerr=%d\n",
         FrmtS.align,FrmtS.blk,FrmtS.jerr);

    n = FrmtS.Ndimen;
 if(n>0)
  {
   fprintf(stdout,"Mdim[%1d] =",n);
   for(i=0;i<n;i++)
    {
     fprintf(stdout," %d", FrmtS.Mdimen[i]);
    }
   fprintf(stdout,"\n");
  }
}



#if IDF_CPP_ == 1
int idf_frmt_begin(int align, int blk, char* target)
#else
int idf_frmt_begin(align,blk,target)
int   align;
int   blk;
char *target;
#endif
{
int ierr=0;
int swrd,wrd;
IDF_FRMT *Frmt;

 idf_frmtreat_nul();

 FrmtS.align   =  align;
 FrmtS.blk     =  blk;
 FrmtS.Otarget =  target;

   Frmt = idf_frmt_address();
if(Frmt == NULL)
 {
  /* format structure was lost */
  swrd=0;
  ierr=1;
 }

else
 {
     wrd  = Frmt->swrd;
  if(wrd  < 1)
   {
    swrd=0;
    ierr=2;
   }
  else
   {
        swrd = idf_align_bnd(align,wrd);
    if(!swrd)
     {
      /* unspecified alignment rule */
      ierr=3;
     }
   }
 }

 FrmtS.swrd  =  swrd;
 FrmtS.jerr  =  ierr;

return ierr;
}



#if IDF_CPP_ == 1
int idf_frmt_set_pntr(void** p, int jsize)
#else
int idf_frmt_set_pntr(p,jsize)
void **p;
int  jsize;
#endif
{
int ierr=0;
void *kv;
char *k;

k = NULL;

if(jsize>0)
 {
     k = (char*) malloc(jsize);
  if(k==NULL) ierr=1;
 }

kv = k;
*p = kv;

return ierr;
}


#if IDF_CPP_ == 1
void idf_frmt_put_pntr(void** p, void* ptr)
#else
void idf_frmt_put_pntr(p,ptr)
void **p;
void *ptr;
#endif
{
*p = ptr;
}


#if IDF_CPP_ == 1
char* idf_frmt_get_pntr(void** p)
#else
char* idf_frmt_get_pntr(p)
void **p;
#endif
{
void *k;

k = *p;

return ( (char*)k );
}



#if IDF_CPP_ == 1
int idf_frmtreat_array(void* val)
#else
int idf_frmtreat_array(val)
void *val;
#endif
{
/* 
check if the target is an array or a single object 
*/
int ierr;
int yes;
register int i;
int junit,ipntr,jaloc,jtype,jsize,jmem,ndim;
IDF_FRMT *Frmt;
char *Atrgt;
char *Tag;
void *tag;

ierr  = 0;
Atrgt = NULL;
Tag   = NULL;

   Frmt  = idf_frmt_address();

if(Frmt == NULL)
 {
  yes=0;
 }

else if(Frmt->NfrmtE != 1)
 {
  yes=0;
 }

else if(Frmt->NfrmtU != 1)
 {
  yes=0;
 }

else if(Frmt->Rfrmt[0] != '%')
 {
  yes=0; 
 }

else
 {
  junit = Frmt->Mrpt[0];

  if(Frmt->Ufrmt[junit].Nrep != 1)
   {
    yes=0;
   }

  else
   {
    ipntr = Frmt->Ufrmt[junit].pointer;
    jaloc = Frmt->Ufrmt[junit].jalloc;
    jtype = Frmt->Ufrmt[junit].type;
    jsize = Frmt->Ufrmt[junit].jsize;
    jmem  = Frmt->Ufrmt[junit].jmem;

    if(ipntr != 0)
     {
      if(jaloc != 0)
       {
        /* allocate memory here */

           Atrgt = (char*) malloc(jmem);
        if(Atrgt == NULL)
         {
          /* unable to allocate memory for pointer target */
          ierr=5;
         }

        else
         {
          tag = Atrgt;
          idf_frmt_put_pntr(val,tag);
          Tag = Atrgt;
         }

        yes=3;
       }

      else
       {
        /* memory has been allocated to pointer */
           Atrgt = idf_frmt_get_pntr(val);
        if(Atrgt==NULL)
         {
          /* improper target pointer */
          ierr=6;
         }
        else
         {
          Tag = Atrgt;
         }

        yes=2;
       }
     }

    else
     {
      /* non-pointer*/
      Tag = (char*) val;
      yes=1; 
     }
   }
 }


if(yes)
 {
  FrmtS.ielem  = 0;
  FrmtS.iunit  = junit;
  FrmtS.Otarget= (char*)val;
  FrmtS.Atarget= Atrgt;
  FrmtS.Target = Tag;
  FrmtS.irep   = 0;
  FrmtS.jrep   = jmem/jsize;
  FrmtS.pointer= ipntr;
  FrmtS.jalloc = jaloc;
  FrmtS.jtype  = jtype;
  FrmtS.jsize  = jsize;
  FrmtS.jmem   = jmem;

  ndim = Frmt->Ufrmt[junit].Ndim;
  FrmtS.Ndimen = ndim;

  if(ndim>0)
   {
    for(i=0;i<ndim;i++)
     {
      FrmtS.Mdimen[i] = Frmt->Ufrmt[junit].Mdim[i];
     }
   }

  FrmtS.jerr   = ierr;

  if(ierr) yes = -1;
 }


return yes;
}



int idf_frmt_get()
{
/*
  <0 - error
   2 - format string exhausted
   1 - next format unit up
*/
int ierr=0;
int flag,krept,jcurs,irep,jrep;
char ctype;
int ofs,gap,align,swrd,blk;
int jtype,jsize,ipntr,jalloc,jmem;
int icurs,junit;
int l,ksize,wsize,amode;
register int i;
IDF_FRMT *Frmt;
char *trgt;
char *Atrgt;
void *tag,*tagg;

 Frmt  = idf_frmt_address();

 irep  = FrmtS.irep;
 jrep  = FrmtS.jrep;

flag=0;
if(irep>0 || jrep>0)
 {
  jrep--;

  if(jrep>0)
   {
    Atrgt = FrmtS.Target;
    jsize = FrmtS.jsize;

    trgt = Atrgt + jsize;

    FrmtS.Target = trgt;
    FrmtS.irep   = irep;
    FrmtS.jrep   = jrep;

    flag=1;
   }

  else
   {
    irep--;

    if(irep>0)
     {
      align = FrmtS.align;
      ofs   = FrmtS.offset;
      gap   = FrmtS.gap;
      swrd  = FrmtS.swrd;
      blk   = FrmtS.blk;
      ipntr = FrmtS.pointer;
      jalloc= FrmtS.jalloc;
      wsize = FrmtS.jwrd;
      amode = FrmtS.amode;
      jmem  = FrmtS.jmem;
      jsize = FrmtS.jsize;
      jrep  = jmem/jsize;

      if(ipntr)
       {
         /* target is a pointer */
         ksize = ipntr;
         idf_align(align,swrd,amode, ksize,wsize, &gap,&ofs);

         if(jalloc)
          {
           /* allocate memory here */

              Atrgt = (char*) malloc(jmem);
           if(Atrgt == NULL)
            {
             /* unable to allocate memory for pointer target */
             ierr=5;
            }

           else
            {
             trgt = FrmtS.Otarget + ofs;

             tag  = trgt;
             tagg = Atrgt;
             idf_frmt_put_pntr(tag,tagg);

             FrmtS.Atarget = Atrgt;
             trgt = Atrgt;
            }
          }

         else
          {
           /* memory has been allocated to pointer */
           trgt = FrmtS.Otarget + ofs;

                tag = trgt;
              Atrgt = idf_frmt_get_pntr(tag);
           if(Atrgt==NULL)
            {
             /* improper target pointer */
             ierr=6;
            }
           else
            {
             FrmtS.Atarget = Atrgt;
             trgt = Atrgt;
            }
          }
        }

       else
        {
         /* non-pointer target */
         ksize = jmem;
         idf_align(align,swrd,amode, ksize,wsize, &gap,&ofs);

         trgt = FrmtS.Otarget + ofs;
        }

      if(ierr)
       {
        flag=-1;
       }
      else
       {
        FrmtS.Target = trgt;
        FrmtS.gap    = gap;
        FrmtS.offset = ofs + ksize;
        FrmtS.irep   = irep;
        FrmtS.jrep   = jrep;

        flag=1;
       }
     }

   }
 }


if(!flag)
 {
  icurs = FrmtS.ielem;
  align = FrmtS.align;
  ofs   = FrmtS.offset;
  gap   = FrmtS.gap;
  swrd  = FrmtS.swrd;
  blk   = FrmtS.blk;

  icurs++;

  while(flag==0)
   {
    /*------------------------------
       rewind if the end is reached
      ------------------------------*/

    if(icurs >= Frmt->NfrmtE)
     {
      if(blk)
       {
        icurs=0;
       }
      else
       {
        /*----------------------
           no more format units
          ---------------------*/
        flag = 2;
       }
     }

    else    
     {
        ctype = Frmt->Rfrmt[icurs];

     if(ctype == '{')
      {
       /*----------------------
          new period openning
        -----------------------*/

       krept = Frmt->Mrpt[icurs];

       /* insert repeatition counter */
       Frmt->Krpt[icurs] = krept;

       icurs = icurs + 1;
      }


     else if(ctype == '}')
      {
       /*----------------------
          continue the period
          return to nearest "{" 
         ----------------------*/

       jcurs = Frmt->Mrpt[icurs];

          krept = Frmt->Krpt[jcurs];
          krept--;
       if(krept)
        {
         /* reduce the counter and go inside "{" */
         Frmt->Krpt[jcurs] = krept;
         icurs = jcurs + 1;
        }
       else
        {
         /* nearest period is exhausted */
         Frmt->Krpt[jcurs] = 0;
         icurs = icurs + 1;
        }
      }


     else if(ctype == '#')
      {
       /*----------------------------
          new format unit is found, 
            but must be skipped
         ----------------------------*/

       junit = Frmt->Mrpt[icurs];

          ipntr = Frmt->Ufrmt[junit].pointer;
       if(ipntr)
        {
         /* pointer target */

         wsize = ipntr;
         jsize = ipntr;
         amode = IDF_AMODE_BINARI;
         irep  = Frmt->Ufrmt[junit].Nrep;

            jalloc = Frmt->Ufrmt[junit].jalloc;
         if(jalloc)
          {
           /* allocate memory here */

           jmem  = Frmt->Ufrmt[junit].jmem;
           ksize = jsize;

           while(irep>0)
            {
             irep--;

             idf_align(align,swrd,amode, ksize,wsize, &gap,&ofs);

             trgt = FrmtS.Otarget + ofs;

              tag = trgt;
                l = idf_frmt_set_pntr(tag,jmem);
             if(l)
              {
               /* unable to allocate memory for pointer target */
               ierr=6;
               break;
              }
             else
              {
               ofs = ofs + ksize;
              }
            }/*irep*/
          }

         else
          {
           /* memory has been allocated to pointer */

           ksize = jsize;

           while(irep>0)
            {
             irep--;
             idf_align(align,swrd,amode, ksize,wsize, &gap,&ofs);
             ofs = ofs + ksize;
            }
          }
        }

       else
        {
         /* non-pointer target */

         irep  = Frmt->Ufrmt[junit].Nrep;
         jtype = Frmt->Ufrmt[junit].type;
         jsize = Frmt->Ufrmt[junit].jsize;
         jmem  = Frmt->Ufrmt[junit].jmem;

         ksize = jmem;
         wsize = idf_frmt_word (jtype);
         amode = idf_frmt_amode(jtype);

         while(irep>0)
          {
           irep--;
           idf_align(align,swrd,amode, ksize,wsize, &gap,&ofs);
           ofs = ofs + ksize;
          }

        }

       icurs = icurs + 1;
      }


     else if(ctype == '%')
      {
       /*---------------------------------
          new 'real' format unit is found
         ---------------------------------*/

       junit = Frmt->Mrpt[icurs];

       irep  = Frmt->Ufrmt[junit].Nrep;
       ipntr = Frmt->Ufrmt[junit].pointer;
       jalloc= Frmt->Ufrmt[junit].jalloc;
       jtype = Frmt->Ufrmt[junit].type;
       jsize = Frmt->Ufrmt[junit].jsize;
       jmem  = Frmt->Ufrmt[junit].jmem;
       jrep  = jmem/jsize;

       if(ipntr)
        {
         /* target is a pointer */

         wsize = ipntr;
         amode = IDF_AMODE_BINARI;
         ksize = ipntr;

         idf_align(align,swrd,amode, ksize,wsize, &gap,&ofs);

         if(jalloc)
          {
           /* allocate memory here */

              Atrgt = (char*) malloc(jmem);
           if(Atrgt == NULL)
            {
             /* unable to allocate memory for pointer target */
             ierr=5;
            }

           else
            {
             trgt = FrmtS.Otarget + ofs;

              tag = trgt;
             tagg = Atrgt;
             idf_frmt_put_pntr(tag,tagg);

             FrmtS.Atarget = Atrgt;
             trgt = Atrgt;
            }
          }

         else
          {
           /* memory has been allocated to pointer */

           trgt = FrmtS.Otarget + ofs;

                tag = trgt;
              Atrgt = idf_frmt_get_pntr(tag);
           if(Atrgt==NULL)
            {
             /* improper target pointer */
             ierr=6;
            }
           else
            {
             FrmtS.Atarget = Atrgt;
             trgt = Atrgt;
            }
          }
        }

       else
        {
         /* non-pointer target */
         ksize = jmem;
         wsize = idf_frmt_word (jtype);
         amode = idf_frmt_amode(jtype);

         idf_align(align,swrd,amode, ksize,wsize, &gap,&ofs);

         FrmtS.Atarget = NULL;
         trgt = FrmtS.Otarget + ofs;
        }

       if(ierr)
        {
         flag = -1;
        }
       else
        {
         FrmtS.iunit  = junit;
         FrmtS.Target = trgt;
         FrmtS.gap    = gap;
         FrmtS.offset = ofs + ksize;
         FrmtS.irep   = irep;
         FrmtS.jrep   = jrep;
         FrmtS.pointer= ipntr;
         FrmtS.jalloc = jalloc;
         FrmtS.jtype  = jtype;
         FrmtS.jsize  = jsize;
         FrmtS.jwrd   = wsize;
         FrmtS.jmem   = jmem;
         FrmtS.amode  = amode;

         flag = 1;
        }
      }


     else
      {
       /*------------------------------------
         improper symbol in format structure
        -------------------------------------*/
       ierr=4;
       flag=-1;
      }

     }/*!rewind*/

    }/*while*/

  FrmtS.ielem = icurs;
  FrmtS.jerr  = ierr;

 }/*flag*/


return flag;
}



void idf_frmtreat_err_prn()
{
#define IDF_ERR_FRMT_TREAT_N 7

static char *ErrFRMT_T[IDF_ERR_FRMT_TREAT_N] = {

/* 1*/ "format structure was lost",
/* 2*/ "improper format max word",
/* 3*/ "improper alignment rule",
/* 4*/ "improper symbol in format structure",
/* 5*/ "unable to allocate memory for pointer target",
/* 6*/ "improper target pointer",
/*  */ "unknown error"
                                               };

int jerr;
register int j;
char *s;

   jerr = FrmtS.jerr;
if(jerr>0)
 {
     j = jerr;
  if(j > IDF_ERR_FRMT_TREAT_N) j = IDF_ERR_FRMT_TREAT_N;
   {
     j--;
    s = ErrFRMT_T[j];

    fprintf(stdout,"format search error=%d element=%d unit=%d:\n%s\n",
            jerr,FrmtS.ielem,FrmtS.iunit,s);

    /*
    if(jerr>4)
     {
      idf_frmtreat_prn();
     }
    */
   }
 }

}
