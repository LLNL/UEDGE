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

 

IDF_DATA *Dtp = NULL;




IDF_DATA* idf_data_address()
{
 return ( Dtp );
}


int idf_data_ini()
{
int ierr=0;
long k;

if(Dtp != NULL)
 {
  ierr=1;
 }

else
 {
        k = sizeof(IDF_DATA);
     Dtp = (IDF_DATA*) malloc(k);
  if(Dtp == NULL)
   {
    ierr=2;
   }
 }

return ierr;
}


void idf_data_end()
{
if(Dtp != NULL)
 {
  free(Dtp);
  Dtp = NULL;
 }
}


void idf_data_prn0()
{
 fprintf(stdout,"data field starts at position=%ld string=%d kursor=%d\n",
         Dtp->kpos0,Dtp->kstr0,Dtp->kurs0); 

 fprintf(stdout,"data type assigned to name in file jtype=%d\n",
         Dtp->jtype); 
}


void idf_data_prn1()
{
 fprintf(stdout,"data unit=%d starts at position=%ld string=%d kursor=%d\n",
        Dtp->iunt+1, Dtp->kpos,Dtp->kstr,Dtp->kurs);
}


#if IDF_CPP_ == 1
int idf_data_begin(IDF_NAME_* A)
#else
int idf_data_begin(A)
IDF_NAME_ *A;
#endif
{
int k;
register int i;
long kpos;

 Dtp->kstr0  = A->kstr;
 Dtp->kurs0  = A->kurs;
      kpos   = A->dpos;
 Dtp->kpos0  = kpos;

 Dtp->kstr   = A->kstr;
 Dtp->kurs   = A->kurs;
 Dtp->kpos   = kpos;

 Dtp->cntr   = '[';
 Dtp->ctor   = ' ';

 Dtp->flag   = 0;
 Dtp->nrpt   = 0;
 Dtp->Jbuf   = 0;
 Dtp->DTPn   = 0;

 Dtp->jtype  = A->jtype;

 Dtp->iunt   = 0;


 /* zeroing the formular result */
 idf_form_nul();

 /* zeroing the data buffer result */
 idf_dbuf_nul();
 /* idf_val_nul ();*/


/* put cursor at the beginning of data field */

   k = idf_file_fseek(kpos);

/*
for(i=0;i<IDF_DTP_MAX;i++)
 {
  Dtp->DTPpos[i] = 0L;
 }

for(i=0;i<IDF_DTP_MAX;i++)
 {
  Dtp->DTPstr[i] = 0;
  Dtp->DTPkur[i] = 0;
  Dtp->DTPrpt[i] = 0;
 }
*/

return k;
}


#if IDF_CPP_ == 1
int idf_data_beginL(IDF_NAME_L* A)
#else
int idf_data_beginL(A)
IDF_NAME_L *A;
#endif
{
int k;
register int i;
long kpos;

 Dtp->kstr0  = A->kstr;
 Dtp->kurs0  = A->kurs;
      kpos   = A->dpos;
 Dtp->kpos0  = kpos;

 Dtp->kstr   = A->kstr;
 Dtp->kurs   = A->kurs;
 Dtp->kpos   = kpos;

 Dtp->cntr   = '[';
 Dtp->ctor   = ' ';

 Dtp->flag   = 0;
 Dtp->nrpt   = 0;
 Dtp->Jbuf   = 0;
 Dtp->DTPn   = 0;

 Dtp->jtype  = A->jtype;

 Dtp->iunt   = 0;


 /* zeroing the formular result */
 idf_form_nul();

 /* zeroing the data buffer result */
 idf_dbuf_nul();
 /*idf_val_nul ();*/


/* put cursor at the beginning of data field */

   k = idf_file_fseek(kpos);

/*
for(i=0;i<IDF_DTP_MAX;i++)
 {
  Dtp->DTPpos[i] = 0L;
 }

for(i=0;i<IDF_DTP_MAX;i++)
 {
  Dtp->DTPstr[i] = 0;
  Dtp->DTPkur[i] = 0;
  Dtp->DTPrpt[i] = 0;
 }
*/

return k;
}


int idf_data_unit()
{
/*
----------------------- 
 get data unit symbols
 or value if formular
-----------------------

 1    no more data
 0    data unit up
-1    error

'[]' are assigned to fictitious start and end of a number
*/

char cntr;
unsigned int c,ci,cerr;
int ibrk,k,kk,nrpt,flag,jerr,irpt,Jtype,Jbuf;
int kstr,kurs;
long kpos,lpos;
int *pKstr,*pKurs;
long *pKpos;
int lp;

kstr  = Dtp->kstr;
kurs  = Dtp->kurs;
kpos  = Dtp->kpos;

pKstr = &kstr;
pKurs = &kurs;
pKpos = &kpos;

cntr = Dtp->cntr;
flag = Dtp->flag;
ci   = Dtp->ctor;
irpt = Dtp->DTPn;

Jtype = 0;
 Jbuf = 0;
 nrpt = 1;
 jerr = 0;
 cerr = ' ';


if(cntr == ';')
 {
  /* end of data field has been reached */
  ibrk=3;
  goto fin;
 }


      ibrk=0;
while(ibrk==0)
 {

 /*--------------------
   get the next symbol 
  ---------------------*/

  if(flag)
   {
    /* return stored symbol */
    c = ci;

    flag = 0;
   }

  else
   {
    /* get next symbol */

      k = idf_file_getc(&c);
      kurs++;
      kpos++;

      if(k==2)
       {
        /* I/O error */
        jerr=1;
	ibrk=2;
       }

      else if(k==1)
       {
	/* unexpected IDF_EOF */
        jerr=2;
        ibrk=2;
       }
   }

      if(jerr)
       {
        ;
       }

    /*---------------------
      comma delimiter case
     ----------------------*/

      else if(c == ',')
       {
             if(cntr == ',') kk=1;
        else if(cntr == ';') kk=1;
        else if(cntr == '*') kk=1;
        else if(cntr == '[') kk=1;
        else                 kk=0;

        if(kk)
         {
          /* extra , */
          cerr=cntr;
          jerr=3;
          ibrk=2;
         }
        else
         {
          if(cntr == ']')
           {
            /* data up*/
            cntr = ',';
            ibrk=1;
           }
          else
           {
            /* skip it */
            cntr = ',';
           }
         }

       }/*end of comma */


    /*------------------------
      simple data period case
     -------------------------*/

      else if(c == '*')
       {
             if(cntr == ',') kk=1;
        else if(cntr == '*') kk=1;
        else if(cntr == '{') kk=1;
        else if(cntr == '}') kk=1;
        else if(cntr == '[') kk=1;
        else                 kk=0;

        if(kk)
         {
          /* improper star */
          cerr=cntr;
          jerr=4;
          ibrk=2;
         }
        else
         {
          if(cntr == ']')
           {
            /* number of repeatitions */
#if IDF_DBUF_FUN == 1
               k = idf_dbuf_to_i_(&lp);
#else
               k = idf_dbuf_to_i(&lp);
#endif
            if(k)
             {
              /* repeatition number conversion error*/
              jerr=4+k;
              ibrk=2;
             }
            else if(lp < 0)
             {
              /* improper repeatition number*/
              jerr=12;
              ibrk=2;
             }
            else
             {
              nrpt = lp;
              cntr = '*';
             }
           }
          else
           {
            /* missing number of repeations */
            jerr=13;
            ibrk=2;
           }
         }

       }/*end of star */


    /*-----------------
      left brace case
     ------------------*/

      else if(c == '{')
       {
        if(cntr == ']')
         {
          /*data up*/
          cntr = ' ';
          ci = c;
          flag=1;
          ibrk=1;
         }
        else if(cntr == '*')
         {
          /* put data period */
          if(nrpt)
           {
            if(irpt < IDF_DTP_MAX_NP)
             {
                 kk = idf_file_ftell(&lpos);
              if(kk==1)
               {
                /*unexpected EoF after data period*/
                jerr=14;
                ibrk=2;
               }
              else if(kk<0)
               {
                /* data period get position error*/
                jerr=15;
                ibrk=2;
               }
              else
               {
                Dtp->DTPstr[irpt] = kstr;
                Dtp->DTPkur[irpt] = kurs;
                Dtp->DTPpos[irpt] = lpos;
                Dtp->DTPrpt[irpt] = nrpt;
                            irpt++;
                cntr = '{';
               }
             }
            else
             {
              /* too many periods */
              jerr=16;
              ibrk=2;
             }
           }

          else
           {
            /* skip zero period */
            ci = c;
               k = idf_skip_data_period(&kstr,&kurs,&kpos, &ci);
            if(k)
             {
              /* period skipping error*/
              ibrk=2;
             }
            else
             {
              cntr = '}';
             }
           }
         }
        else
         {
          /* unity repeatition period */
            if(irpt < IDF_DTP_MAX_NP)
             {
                 kk = idf_file_ftell(&lpos);
              if(kk==1)
               {
                /*unexpected EoF after data period*/
                jerr=14;
                ibrk=2;
               }
              else if(kk<0)
               {
                /* data period get position error*/
                jerr=15;
                ibrk=2;
               }
              else
               {
                Dtp->DTPstr[irpt] = kstr;
                Dtp->DTPkur[irpt] = kurs;
                Dtp->DTPpos[irpt] = lpos;
                Dtp->DTPrpt[irpt] = 1;
                            irpt++;
                cntr = '{';
               }
             }
            else
             {
              /* too many periods */
              jerr=16;
              ibrk=2;
             }
         }
       }/* end of { */


    /*-----------------
      right brace case
     ------------------*/

      else if(c == '}')
       {
             if(cntr==',') k=1;
        else if(cntr=='*') k=1;
        else if(cntr=='[') k=1;
        else               k=0;

        if(k)
         {
          /* unexpected } */
          cerr=cntr;
          jerr=17;
          ibrk=2;
         }
        else if(cntr == ']')
         {
          /* data up*/
          cntr = ' ';
          ci = c;
          flag=1;
          ibrk=1;          
         }
        else if(irpt < 1)
         {
          /* extra } */
          jerr=18;
          ibrk=2;
         }
        else
         {
          kk = irpt - 1;

             lp = Dtp->DTPrpt[kk];
             lp--;
          if(lp < 1)
           {
            /* nearest period is over*/
            irpt = kk;
            cntr = '}';
           }
          else
           {
            /* continue the nearest period*/
            kstr = Dtp->DTPstr[kk];
            kurs = Dtp->DTPkur[kk];
            lpos = Dtp->DTPpos[kk];

            Dtp->DTPrpt[kk] = lp;

               k = idf_file_fseek(lpos);
            if(k)
             {
              /* return i/o error */
              jerr=19;
              ibrk=2;
             }
            else
             {
              /* return to period beginning */
              cntr = '{';
             }
           }
         }

       }/* end of } */


    /*-------------
      comment case
     --------------*/

      else if(c == '/')
       {
        if(cntr == ']')
         {
          cntr = ' ';
          ci = c;
          flag=1;
          ibrk=1;
         }
        else
         {
          /* skip comment if appropriate */

          ci = c;

          k = idf_stxt_skip_cmnt(pKstr,pKurs,pKpos, &ci);

          if(k ==  2)
           {
            /* improper symbol '/' */
            jerr = 20;
            cerr = ci;
            ibrk=2;
           }
          else if(k ==  1)
           {
            /* unexpected IDF_EOF */
            jerr = 2;
            ibrk=2;
           }
          else if(k < 0)
           {
            /* comment error */
            ibrk=2;
           }
          else
           {
            /* store symbol */
            flag=1;
           }
         }

       } /* '/' */


    /*-------------
      string  case 
     --------------*/
 
      else if(c == '\"')
       {
        if(cntr == ']')
         {
          if(nrpt==0)
           {
            /*skip data unit with zero repeatition*/
            nrpt=1;
            cntr = ' ';
            ci = c;
            flag=1;
           }
          else
           {
            /*data up*/
            cntr = ' ';
            ci = c;
            flag=1;
            ibrk=1;
           }
         }
        else
         {
             ci = c;
          Jtype = IDF_FRMT_S;
          Jbuf  = 1;

          idf_dbuf_begin(pKstr,pKurs,pKpos, Jtype);

          k = idf_dbuf_str(pKstr,pKurs,pKpos, &ci);

          if(k < 0)
           {
            /* string search error */
            ibrk=2;
           }
          else
           {
            /* string placed inside buffer */
            cntr = ']';
            flag = 1;
           }
         }

       } /* "" */


    /*-------------
       text  case 
     --------------*/
 
      else if(c == '#')
       {
        if(cntr == ']')
         {
          if(nrpt==0)
           {
            /*skip data unit with zero repeatition*/
            nrpt=1;
            cntr = ' ';
            ci = c;
            flag=1;
           }
          else
           {
            /*data up*/
            cntr = ' ';
            ci = c;
            flag=1;
            ibrk=1;
           }
         }
        else
         {
             ci = c;
          Jtype = IDF_FRMT_T;
          Jbuf  = 1;

          idf_dbuf_begin(pKstr,pKurs,pKpos, Jtype);

          k = idf_dbuf_txt(pKstr,pKurs,pKpos, &ci);

          if(k < 0)
           {
            /* text search error */
            ibrk=2;
           }
          else
           {
            /* text placed inside buffer */
            cntr = ']';
            flag=1;
           }
         }

       } /* # # */


    /*---------------
      character case 
     ----------------*/

      else if(c == '\'')
       {
        if(cntr == ']')
         {
          if(nrpt==0)
           {
            /*skip data unit with zero repeatition*/
            nrpt=1;
            cntr = ' ';
            ci = c;
            flag=1;
           }
          else
           {
            /*data up*/
            cntr = ' ';
            ci = c;
            flag=1;
            ibrk=1;
           }
         }
        else
         {
             ci = c;
          Jtype = IDF_FRMT_C;
          Jbuf  = 1;

          idf_dbuf_begin(pKstr,pKurs,pKpos, Jtype);

          k = idf_dbuf_char(pKstr,pKurs,pKpos, &ci);

          if(k < 0)
           {
            /* character search error */
            ibrk=2;
           }
          else
           {
            /* character placed inside buffer */
            cntr = ']';
            flag=1;
           }
         }

       } /* '' */


    /*---------------
      numerical case 
     ----------------*/

      else if( idf_data_symbol(c) )
       {
        if(cntr == ']')
         {
          if(nrpt==0)
           {
            /*skip data unit with zero repeatition*/
            nrpt=1;
            cntr = ' ';
            ci = c;
            flag=1;
           }
          else
           {
            /*data up*/
            cntr = ' ';
            ci = c;
            flag=1;
            ibrk=1;
           }
         }
        else
         {
             ci = c;
          Jtype = IDF_FRMT_DIGIT;
          Jbuf  = 1;

          idf_dbuf_begin(pKstr,pKurs,pKpos, Jtype);

          k = idf_dbuf_number(pKstr,pKurs,pKpos, &ci);

          if(k == -1)
           {
            /* I/O error */
            jerr = 1;
            ibrk=2;
           }
          else if(k == -2)
           {
            /* unexpected IDF_EOF */
            jerr = 2;
            ibrk=2;
           }
          else if(k == -3)
           {
            /* too many symbols for number */
            jerr = 21;
            ibrk=2;
           }
          else
           {
            /* number symbols placed inside buffer */
            cntr = ']';
            flag=1;
           }

         }

       } /* numerics */


    /*---------------
      complex case 
     ----------------*/

      else if(c == '(' )
       {
        if(cntr == ']')
         {
          if(nrpt==0)
           {
            /*skip data unit with zero repeatition*/
            nrpt=1;
            cntr = ' ';
            ci = c;
            flag=1;
           }
          else
           {
            /*data up*/
            cntr = ' ';
            ci = c;
            flag=1;
            ibrk=1;
           }
         }
        else
         {
             ci = c;
          Jtype = IDF_FRMT_CMPLX;
          Jbuf  = 1;

          idf_dbuf_begin(pKstr,pKurs,pKpos, Jtype);

          k = idf_dbuf_cmplx(pKstr,pKurs,pKpos, &ci);

          if(k < 0)
           {
            /* complex search error */
            ibrk=2;
           }
          else
           {
            /* complex is placed inside buffer */
            cntr = ']';
            flag=1;
           }
         }

       } /* complex */


    /*---------------
      formular case 
     ----------------*/

      else if(c == '$')
       {
        if(cntr == ']')
         {
          if(nrpt==0)
           {
            /*skip data unit with zero repeatition*/
            nrpt=1;
            cntr = ' ';
            ci = c;
            flag=1;
           }
          else
           {
            /*data up*/
            cntr = ' ';
            ci = c;
            flag=1;
            ibrk=1;
           }
         }
        else
         {
          ci   = c;
          Jtype = IDF_FRMT_VOID;
          Jbuf = 2;

          k = idf_form_do(pKstr,pKurs,pKpos, &ci);

          if(k)
           {
            /* formular search error */
            ibrk=2;
           }
          else
           {
            /* formular results is placed inside buffer */
            cntr = ']';
            flag=1;
           }
         }

       } /* $$ */


    /*-------------
       EOS case
     --------------*/

      else if( idf_string_over(c) )
       {
        kstr++;
        kurs=0;
       }


    /*----------------------------
      white space delimiters case
     -----------------------------*/

      else if(c == ' ' )
       {
        ;
       }

      else if(c == '\t')
       {
#if IDF_TXT_TABULATORS == 1
        kurs+=7;
#else
        ;
#endif
        }

      else if(c == '\v')
       {
#if IDF_TXT_VERTICES == 1
        kstr+=8;
#else
        ;
#endif
        }



    /*-------------------
       end of datas case
     --------------------*/

      else if(c == ';')
       {
             if(cntr == ',') kk=1;
        else if(cntr == '*') kk=1;
        else if(cntr == '{') kk=1;
        else if(cntr == '[') kk=1;
        else                 kk=0;
        if(kk)
         {
          /* unexpected end of data field */
          cerr=cntr;
          jerr=22;
          ibrk=2;
         }
        else
         {
          if(cntr == ']')
           {
            if(nrpt==0)
             {
              /*skip data unit with zero repeatition*/
              nrpt=1;
              cntr = ';';
              ci = c;
              flag=1;
             }
            else
             {
              /*data up*/
              cntr = ';';
              ci = c;
              flag=1;
              ibrk=1;
             }
           }
          else if(irpt)
           {
            /* unclosed data period { before data end*/
            cerr=';';
            jerr=23;
            ibrk=2; 
           }
          else
           {
            /*no more data units*/
            ibrk=3;
           }
         }
       } /* ; */


    /*-----------------------------------------
      all other symbols should not appear here
     ------------------------------------------ */

      else
       {
        if(cntr == ']')
         {
            if(nrpt==0)
             {
              /*skip data unit with zero repeatition*/
              nrpt=1;
              cntr = ' ';
              ci = c;
              flag=1;
             }
            else
             {
              /*data up*/
              cntr = ' ';
              ci = c;
              flag=1;
              ibrk=1;
             }          
         }
        else
         {
          /* improper symbol in data unit field */
          cerr = c;
          jerr = 24;
          ibrk=2;
         }
       }

 }/*ibrk*/


 Dtp->kstr = kstr;
 Dtp->kurs = kurs;
 Dtp->kpos = kpos;

 Dtp->cntr = cntr;
 Dtp->nrpt = nrpt;
 Dtp->flag = flag;
 Dtp->ctor = ci;
 Dtp->DTPn = irpt;
 Dtp->Jbuf = Jbuf;


if(jerr)
 {
  k = jerr + 90;
  idf_txt_err_put(k, kstr,kurs,kpos, cerr);
 }

fin:
     if(ibrk==1) flag=0;  /* data unit up */
else if(ibrk==2) flag=-1; /* error        */
else             flag=1;  /* no more data */

return flag;
}




int idf_data_value()
{
/*
----------------------- 
  get data unit value
-----------------------

 1    no more data
 0    data unit value up
-1    error in getting data unit
-2    error in buffer coversion
-3    error in formular
*/
unsigned int cerr = ' ';
int flag,jerr,k,kk;
int Jbuf,Nrpt,Jtype;
IDF_FORM *Df;
IDF_DBUF *Db;

Db = idf_dbuf_address();
Df = idf_form_address();

Jbuf = Dtp->Jbuf;
Nrpt = Dtp->nrpt;

flag=0;
jerr=0;


/*-----------------------
    get next data unit
  if repeatition is over
 -------------------------*/

if(Nrpt < 1)
 {
  Dtp->Jbuf = 0;
  Df->jtype = 0;
  Db->jtype = 0;

  /* get data unit */

     k = idf_data_unit();

  if(k==1)
   {
    /* no more data */
    flag=1;
    goto fin;
   }

  else if(k==-1)
   {
    /* data unit getting error*/

    if(Dtp->Jbuf==2)
     {
      flag = -3;
     }
    else
     {
      flag = -1;
     }

    goto err;
   }

  else
   {
    Jbuf = Dtp->Jbuf;

#if IDF_DATA_TRACE == 1
    if(Jbuf==1)
     {
      idf_dbuf_prn();
     }
    else if(Jbuf==2)
     {
      idf_form_prn();
     }
#endif

       Nrpt = Dtp->nrpt;
    if(Nrpt < 1)
     {
      /* improper repeatition number */
      jerr = 1;
      flag = -1;
      goto err;
     }
   }
 }


if(Jbuf == 1)
{
 /*-----------------------------------
  convert buffer according to name type
     or data type in file itself
  ------------------------------------*/

     k = idf_val_buf();
  if(k)
   {
    /* buffer conversion error */
    flag = -2;
   }

#if IDF_DATA_TRACE == 1
   idf_val_prn();
#endif
}


else if(Jbuf==2)
{
/*------------------------
  convert formular value
    into output value 
 -------------------------*/

    k = idf_val_form();
 if(k)
   {
    /* formular conversion error */
    flag = -3;
   }

#if IDF_DATA_TRACE == 1
   idf_val_prn();
#endif
}

else
{
 /* improper raw data unit storage type*/
 jerr=2;
 flag = -1;
}


err:
if(jerr)
 {
  k = jerr + 153;
  idf_txt_err_put(k, Dtp->kstr,Dtp->kurs,Dtp->kpos, cerr);
 }

if(Nrpt > 0)
 {
  Nrpt = Nrpt - 1;
  Dtp->nrpt = Nrpt;
 }

if(flag==0)
 {
  Dtp->iunt = Dtp->iunt + 1;
 }

fin:
return flag;
}



#if IDF_CPP_ == 1
int idf_data_up(int jfrmt, void* val, int nval)
#else
int idf_data_up(jfrmt, val,nval)
int   jfrmt;
void *val;
int   nval;
#endif
{
/*
-------------------------------------- 
 convert and send the data unit value 
    according to external format
--------------------------------------
*/
unsigned int cerr=' ';
int k;
int jerr;

jerr = 0;

    /* convert according to input format */
       k = idf_val_get(jfrmt, val,nval);
    if(k == -1)
     {
      /* incompatible types for conversion */
      jerr=1;
     }
    else if(k == -2)
     {
      /* overflow in conversion */
      jerr=2;
     }

if(jerr)
 {
  k = jerr + 155;
  idf_txt_err_put(k, Dtp->kstr,Dtp->kurs,Dtp->kpos, cerr);
 }

return jerr;
}
