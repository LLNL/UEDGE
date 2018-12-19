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



IDF_FBUF *Fbuf = NULL;




IDF_FBUF* idf_fbuf_address()
{
 return ( Fbuf );
}



int idf_fbuf_ini()
{
int ierr=0;
int k;

if(Fbuf != NULL)
 {
  ierr=1;
 }

else
 {
        k = sizeof(IDF_FBUF);
     Fbuf = (IDF_FBUF*) malloc(k);
  if(Fbuf == NULL)
   {
    ierr=2;
   }
 }

return ierr;
}


void idf_fbuf_end()
{
if(Fbuf != NULL)
 {
  free(Fbuf);
  Fbuf = NULL;
 }
}


void idf_fbuf_prn()
{
 if(Fbuf != NULL)
 {
  fprintf(stdout,"formular unit buffer has size=%d and type=%1d\n",
          Fbuf->length,Fbuf->jtype);

  fprintf(stdout,"it starts at position=%ld string=%d kursor=%d\n",
          Fbuf->pos,Fbuf->kstr,Fbuf->kurs);

  if(Fbuf->length)
   {
    fprintf(stdout,"content of character buffer:\n%s\n",
            Fbuf->buf);
   }

  if(Fbuf->jconv==1)
   {
    fprintf(stdout,"buffer was converted to constant=(%e , %e) of type=%d\n",
            Fbuf->val[0],Fbuf->val[1],Fbuf->jf);    
   }
  else if(Fbuf->jconv==2)
   {
    fprintf(stdout,"buffer was converted to function=%d\n",
            Fbuf->jf);    
   }
 }
}


#if IDF_CPP_ == 1
void idf_fbuf_begin(int* Kstr, int* Kurs, long* Kpos, int Jtype)
#else
void idf_fbuf_begin(Kstr,Kurs,Kpos, Jtype)
int  *Kstr,*Kurs;
long *Kpos;
int   Jtype;
#endif
{
 Fbuf->buf[0] = IDF_EOS;
 Fbuf->length = 0;

 Fbuf->kstr   = *Kstr;
 Fbuf->kurs   = *Kurs;
 Fbuf->pos    = *Kpos;

 Fbuf->jtype  = Jtype;

 Fbuf->jconv  = 0;
 Fbuf->jf     = 0;
 Fbuf->val[0] = 0.0;
 Fbuf->val[1] = 0.0;
}


void idf_fbuf_nul()
{
 Fbuf->buf[0] = IDF_EOS;
 Fbuf->length = 0;

 Fbuf->kstr   = 0;
 Fbuf->kurs   = 0;
 Fbuf->pos    = 0;

 Fbuf->jtype  = 0;

 Fbuf->jconv  = 0;
 Fbuf->jf     = 0;
 Fbuf->val[0] = 0.0;
 Fbuf->val[1] = 0.0;
}



#if IDF_CPP_ == 1
int idf_fbuf_number(int* Kstr, int* Kurs, long* Pos,
                    unsigned int* Symb)
#else
int idf_fbuf_number(Kstr,Kurs,Pos, Symb)
int  *Kstr,*Kurs;
long *Pos;
unsigned int *Symb;
#endif
{
int ierr=0;
int ibrk,istp,k,kd;
unsigned int c,ci,cntr;
register int i;
int kstr,kurs;
long pos;
char *Buf;

 Buf = Fbuf->buf;

   i = 0;
  kd = 0;
   c = *Symb;
kstr = *Kstr;
kurs = *Kurs;
 pos = *Pos;
cntr = ' ';

if(c=='.')
 {
  istp=1;
  Buf[i] = (char)c;
      i++;
 }
else if( isdigit(c) )
 {
  istp=0;
  kd++;
  Buf[i] = (char)c;
      i++;
 }
else
 {
  /* number in formular starts with improper symbol*/
  Buf[i] = IDF_EOS;
  ierr=1;
  goto err;
 }


      ibrk=0;
while(ibrk==0)
 {
  if(i > IDF_FBUF_LENGTH)
   {
    /*too many symbols for a formular number*/
    ierr=5;
    ibrk=2;
   }

  else
   {
    k = idf_file_getc(&ci);
    kurs++;
    pos++;

    if(k==2)
     {
      /* I/O error */
      ierr=2;
      ibrk=2;
     }

    else if(k==1)
     {
      /* unexpected EOF */
      kstr++;
      kurs=0;

      ierr=3;
      ibrk=2;
     }

    else if( idf_string_over(ci) )
     {
      /* it is not a comment */
      kstr++;
      kurs=0;
     }

    else if(ci=='.')
     {
      if(istp==0)
       {
        istp=1;
        Buf[i] = (char)ci;
            i++;
       }
      else
       {
        /*improper floating point in formular*/
        ierr=4;
        ibrk=2;
       }
     }

    else if( idf_exp_symbol(ci) )
     {
      if(istp>2)
       {
        /*improper float number style in formular*/
        ierr=6;
        ibrk=2;
       }
      else
       {
        istp=3;
        Buf[i] = 'e';
            i++;
       }
     }

    else if(ci=='+' || ci=='-')
     {
      if(istp<3)
       {
        /* number has no exponent*/
        cntr = ci;
        ibrk=1;
       }
      else if(istp==3)
       {
        istp=4;
        Buf[i] = (char)ci;
            i++;
       }
      else if(istp==4)
       {
        /*improper float number style in formular*/
        ierr=6;
        ibrk=2;
       }
      else
       {
        /* number with exponent*/
        cntr = ci;
        ibrk=1;
       }
     }

    else if( isdigit(ci) )
     {
      kd++;
      if(istp==1)
       {
        istp=2;
       }
      else if(istp==3 || istp==4)
       {
        istp=5;
       }

      Buf[i] = (char)ci;
          i++;
     }

    else if(ci == ' ')
     {
      ;
     }

    else if(ci == '\t')
     {
#if IDF_TXT_TABULATORS == 1
      kurs+=7;
#else
      ;
#endif
     }

    else if(ci == '\v')
     {
#if IDF_TXT_VERTICES == 1
      kstr+=8;
#else
      ;
#endif
     }

    else if(ci==')' || ci==',')
     {
      if(istp==3 || istp==4)
       {
        /*improper symbol in formular number*/
        ierr=7;
        ibrk=2;
       }
      else
       {
        /* end of number*/
        cntr=ci;
        ibrk=1;
       }
     }

    else if(ci=='*' || ci=='^' || ci=='/')
     {
      if(istp==3 || istp==4)
       {
        /*improper symbol in formular number*/
        ierr=7;
        ibrk=2;
       }
      else
       {
        /* end of number*/
        cntr=ci;
        ibrk=1;
       }
     }

    else if(ci=='$')
     {
      /* end of number and formular*/
      cntr=ci;
      ibrk=1;
     }

    else if(ci==';')
     {
      /* unexpected end of data*/
      ierr=8;
      ibrk=2;
     }

    else
     {
      /*improper symbol in formular number*/
      ierr=7;
      ibrk=2;
     }

   }/*i<Nbuf*/

 } /* while*/


if(!ierr)
 {
  if(i==0 || kd==0)
   {
    /* empty number*/
    ierr=9;
   }
 }

 Buf[i] = IDF_EOS;


err:
Fbuf->length = i;

*Kstr = kstr;
*Kurs = kurs;
*Pos  = pos;
*Symb = cntr;

return ierr;
}



#if IDF_CPP_ == 1
int idf_fbuf_name(int* Kstr, int* Kurs, long* Pos, 
                  unsigned int* Symb)
#else
int idf_fbuf_name(Kstr,Kurs,Pos, Symb)
int  *Kstr,*Kurs;
long *Pos;
unsigned int *Symb;
#endif
{
int ierr=0;
int ibrk,k,kd;
unsigned int c,ci,cntr;
register int i;
int kstr,kurs;
long pos;
char *Buf;

 Buf = Fbuf->buf;

   i = 0;
  kd = 0;
   c = *Symb;
kstr = *Kstr;
kurs = *Kurs;
 pos = *Pos;
cntr = ' ';


if(
   ( isalpha(c) )
#if IDF_FOREIGN_LANGUAGE == 1
   || ( !isascii(c) )
#endif
  )
 {
  kd++;
  Buf[i] = (char)c;
      i++;
 }

else
 {
  /* name in formular starts with improper symbol*/
  Buf[i] = IDF_EOS;
  ierr=10;
  goto err;
 }


      ibrk=0;
while(ibrk==0)
 {
  if(i > IDF_FBUF_LENGTH)
   {
    /*too many symbols for a formular name*/
    ierr=11;
    ibrk=2;
   }

  else
   {
    k = idf_file_getc(&ci);
    kurs++;
    pos++;

    if(k==2)
     {
      /* I/O error */
      ierr=2;
      ibrk=2;
     }

    else if(k==1)
     {
      /* unexpected EOF */
      kstr++;
      kurs=0;

      ierr=3;
      ibrk=2;
     }

    else if( idf_string_over(ci) )
     {
      /* it is not a comment */
      kstr++;
      kurs=0;

#if IDF_EOS_BREAK_FNAME == 1
      /* EOS is a separator for name*/
      cntr = ' ';
      ibrk=1;
#endif
     }

    else if(
             ( isalpha(ci) )
#if IDF_FOREIGN_LANGUAGE == 1
            || ( !isascii(ci) )
#endif
           )
     {
      kd++;
      Buf[i] = (char)ci;
          i++;
     }

    else if( isdigit(ci) )
     {
      if(kd==0)
       {
        ierr=10;
        ibrk=2;
       }
      else
       {
        Buf[i] = (char)ci;
            i++;
       }
     }

    else if(ci == '_')
     {
      if(kd==0)
       {
        ierr=10;
        ibrk=2;
       }
      else
       {
        Buf[i] = (char)ci;
            i++;
       }
     }

    else if( idf_form_symbol(ci) )
     {
      /* end of name */
      cntr = ci;
      ibrk=1;
     }

    else if(ci == ' ')
     {
      /* white space is a separator for name*/
      cntr = ' ';
      ibrk=1;
     }

    else if(ci == '\t')
     {
#if IDF_TXT_TABULATORS == 1
      kurs+=7;
#endif
      /* white space is a separator for name*/
      cntr = ' ';
      ibrk=1;
     }

    else if(ci == '\v')
     {
#if IDF_TXT_VERTICES == 1
      kstr+=8;
#endif
      /* white space is a separator for name*/
      cntr = ' ';
      ibrk=1;
     }

    else if(ci=='$')
     {
      /* end of name and formular*/
      cntr=ci;
      ibrk=1;
     }

    else if(ci==';')
     {
      /* unexpected end of data*/
      ierr=8;
      ibrk=2;
     }

    else
     {
      /*improper symbol in formular name*/
      ierr=12;
      ibrk=2;
     }

   }/*i<Nbuf*/

 } /* while*/


if(!ierr)
 {
  if(i==0 || kd==0)
   {
    /* empty number*/
    ierr=13;
   }
 }


Buf[i] = IDF_EOS;


err:
Fbuf->length = i;

*Kstr = kstr;
*Kurs = kurs;
*Pos  = pos;
*Symb = cntr;

return ierr;
}



int idf_fbuf_conv()
{
int k,jcns,jfun,jlc;
int flag;
int typ;
double d[2];
char *buf;
int   n;

flag = 0;
   k = Fbuf->jtype;
 buf = Fbuf->buf;
   n = Fbuf->length;

if(k == IDF_FBUF_TYPE_NUMBER)
 {
     flag = idf_str_to_d(buf,n, d);
  if(flag)
   {
    flag = flag + 14;
   }
  else
   {
    /* digital constant*/
    Fbuf->jconv  = 1;
    Fbuf->jf     = IDF_FRMT_D;
    Fbuf->val[0] = d[0];
    Fbuf->val[1] = 0.0;
   }
 }

else if(k == IDF_FBUF_TYPE_NAME)
 {
     jcns = idf_clist(buf,n, d);
  if(jcns)
   {
    /* constant from standard list*/
    Fbuf->jconv  = 1;
    Fbuf->jf     = IDF_FRMT_D;
    Fbuf->val[0] = d[0];
    Fbuf->val[1] = 0.0;
   }
  else
   {
       jfun = idf_keywF(buf,n);
    if(jfun)
     {
      /* function from standard list*/
      Fbuf->jconv = 2;
      Fbuf->jf    = jfun;
     }
    else
     {
         jlc = idf_list_getLV(buf,n, &typ, d);
      if(jlc==0)
       {
        if(typ==IDF_FRMT_D || typ==IDF_FRMT_W)
         {
          /* constant from local name list*/
          Fbuf->jconv  = 1;
          Fbuf->jf     = typ;
          Fbuf->val[0] = d[0];
          Fbuf->val[1] = d[1];
         }
        else
         {
          /*inconsistency in local name list type*/
          flag=20;
         }
       }
      else if(jlc==1)
       {
        /* unknown name*/
        flag=21;
       }
      else if(jlc==2)
       {
        /* no value*/
        flag=22;
       }
      else
       {
        /* error in list name*/
        flag=23;
       }
     }
   }
 }

else
 {
  /* improper buffer type */
  flag=14;
 }

return flag;
}
