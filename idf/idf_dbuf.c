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




IDF_DBUF *Dbuf = NULL;




IDF_DBUF* idf_dbuf_address()
{
 return ( Dbuf );
}



int idf_dbuf_ini()
{
int ierr=0;
int k;

if(Dbuf != NULL)
 {
  ierr=1;
 }

else
 {
        k = sizeof(IDF_DBUF);
     Dbuf = (IDF_DBUF*) malloc(k);
  if(Dbuf == NULL)
   {
    ierr=2;
   }
 }

return ierr;
}


void idf_dbuf_end()
{
if(Dbuf != NULL)
 {
  free(Dbuf);
  Dbuf = NULL;
 }
}


void idf_dbuf_nul()
{
 Dbuf->jtype  = 0;
}


void idf_dbuf_prn()
{
 if(Dbuf != NULL)
 {
  fprintf(stdout,"data buffer has size=%d and type=%1d\n",
          Dbuf->length,Dbuf->jtype);

  fprintf(stdout,"it starts at position=%ld string=%d kursor=%d\n",
          Dbuf->pos,Dbuf->kstr,Dbuf->kurs);

  if(Dbuf->length)
   {
    fprintf(stdout,"content of data buffer:\n%s\n",
            Dbuf->Data_buf);
   }
 }
}


#if IDF_CPP_ == 1
void idf_dbuf_begin(int* Kstr, int* Kurs, long* Kpos, int Jtype)
#else
void idf_dbuf_begin(Kstr,Kurs,Kpos, Jtype)
int  *Kstr,*Kurs;
long *Kpos;
int   Jtype;
#endif
{
 Dbuf->Data_buf[0] = IDF_EOS;
 Dbuf->length = 0;

 Dbuf->kstr   = *Kstr;
 Dbuf->kurs   = *Kurs;
 Dbuf->pos    = *Kpos;

 Dbuf->jtype  = Jtype;
}



#if IDF_CPP_ == 1
int idf_dbuf_number(int* Kstr,int* Kurs,long* Kpos,
                    unsigned int* Symb)
#else
int idf_dbuf_number(Kstr,Kurs,Kpos, Symb)
int  *Kstr,*Kurs;
long *Kpos;
unsigned int *Symb;
#endif
{
/*
------------------------
 put 'number' symbols 
     into buffer
------------------------

  1 - empty data
  0 - OK, number buffer done
 -1 - I/O error
 -2 - unexpected IDF_EOF
 -3 - too long number
*/

unsigned int ci;
int k,kk,ibrk,flag;
register int i;
int  kstr,kurs;
long kpos;
char *Buf;

Buf = Dbuf->Data_buf;

    i  = 0;
   ci = *Symb;
if(ci != ' ')
 {
  Buf[i] = (char)ci;
      i++;
 }
Buf[i] = IDF_EOS;

Dbuf->length = i;


kpos = *Kpos;
kstr = *Kstr;
kurs = *Kurs;

Dbuf->kstr   = kstr;
Dbuf->kurs   = kurs;
Dbuf->pos    = kpos;

flag = 0;

      ibrk=0;
while(ibrk==0)
  {
      /* get next symbol */
       k = idf_file_getc(&ci);
       kurs++;
       kpos++;

       if(k==2)
        {
	 /* I/O error */
	 flag = -1;
         ibrk=1;
	}

       else if(k==1)
	{
	 /* unexpected IDF_EOF */
	 flag = -2;
         ibrk=1;
        }

       else if( idf_string_over(ci) )
        {
         kstr++;
         kurs=0;
#if IDF_EOS_BREAK_NUMBER == 1
         if(i)
          {
           ci = ' ';
           ibrk=2;
          }
         else
          {
           /* empty data buffer */
           flag=1;
           ibrk=2;
          }
#endif
	}

       else if(i >= IDF_DBUF_LENGTH)
        {
	 /* too long number */
	 flag = -3;
         ibrk=1;
        }

       else if( idf_data_symbol(ci) )
	 {
          Buf[i] = (char)ci;
              i++;
         }

       else
        {
#if IDF_TXT_TABULATORS == 1
         if(ci == '\t')
          {
           kurs+=7;
           ci = ' ';
          }
#endif

#if IDF_TXT_VERTICES == 1
         if(ci == '\v')
          {
           kstr+=8;
           ci = ' ';
          }
#endif

         if(i)
          {
           ibrk=2;
          }
         else
          {
           /* empty data buffer */
           flag=1;
           ibrk=2;
          }
	}

  }/*ibrk*/


Buf[i] = IDF_EOS;
Dbuf->length = i;

*Kstr = kstr;
*Kurs = kurs;
*Kpos = kpos;

*Symb = ci;

return flag;
}



#if IDF_CPP_ == 1
int idf_dbuf_char(int* Kstr, int* Kurs, long* Kpos,
                  unsigned int* Symb)
#else
int idf_dbuf_char(Kstr,Kurs,Kpos, Symb)
int *Kstr,*Kurs;
long *Kpos;
unsigned int *Symb;
#endif
{
/*
-------------------------------------
 put 'character' symbols (between '')
        into buffer
-------------------------------------

  0 - OK
 -1 - I/O error
 -2 - unexpected IDF_EOF
 -3 - improper symbol in character list
 -4 - improper char
 -5 - empty char
 -6 - improper char case
*/

unsigned int ci;
int k,kc,ibrk,flag,is,io;
register int i;
int kstr,kurs;
long kpos;
char *Buf;


Buf = Dbuf->Data_buf;
Buf[0] = IDF_EOS;

           i = 0;
Dbuf->length = i;

kpos = *Kpos;
kstr = *Kstr;
kurs = *Kurs;

Dbuf->kstr   = kstr;
Dbuf->kurs   = kurs;
Dbuf->pos    = kpos;

   ci = *Symb;
if(ci != '\'')
 {
  flag = -6;
  goto err;
 }

flag = 0;
  kc = 0;
  is = 0;
  io = 0;

      ibrk=0;
while(ibrk==0)
  {
      /* get next symbol */
       k = idf_file_getc(&ci);
       kurs++;
       kpos++;

       if(k==2)
        {
	 /* I/O error */
	 flag = -1;
         ibrk=2;
	}

       else if(k==1)
	{
	 /* unexpected IDF_EOF */
          kstr++;
          kurs=0;

	  flag = -2;
          ibrk=2;
	 }

        else if( idf_string_over(ci) )
	 {
	  kstr++;
          kurs=0;
	 }

        else if(kc)
         {
          if(ci == '\'')
           {
            /* concatenation case */
            kc = 0;
           }
          else
           {
            /* char finished */
            ibrk=1;
           }
         }

        else if(ci == '\'')
         {
          if(io>0)
           {
            /* end of char */
            kc++;
           }
          else if(is==1)
           {
            /* escape sequence */
            io++;
            is--;

            Buf[i] = ci;
                i++;
           }
          else
           {
            /* empty */
            flag = -5;
            ibrk=2;
           }
         }

        else if(ci == '\\')
         {
          if(io)
           {
            flag = -4;
            ibrk=2;
           }
          else if(is)
           {
            /*escape sequence */
            io++;
            is--;
            Buf[i] = ci;
                i++;
           }
          else
           {
            is++;

            Buf[i] = ci;
                i++;
           }
         }

        else
         {
          if(ci == '\t')
	   {
#if IDF_TXT_TABULATORS == 1
            kurs+=7;
#endif
            flag = -3;
            ibrk=2;
	   }

          else if(ci == '\v')
	   {
#if IDF_TXT_VERTICES == 1
	    kstr+=8;
#endif
            flag = -3;
            ibrk=2;
	   }

          else if(is)
           {
            if(io>2)
             {
              flag = -4;
              ibrk=2;
             }
            else
             {
              io++;

              Buf[i] = ci;
                  i++;
             }
           }

          else
           {
            if(io)
             {
              flag = -3;
              ibrk=2;
             }
            else
             {
              io++;

              Buf[i] = ci;
                  i++;
             }
           }
         }

  } /* ibrk*/


err:
if(flag < 0)
 {
  ibrk = 8 - flag;
  idf_txt_err_put(ibrk, kstr,kurs,kpos, ci);
 }

Buf[i] = IDF_EOS;
Dbuf->length = i;

*Kstr = kstr;
*Kurs = kurs;
*Kpos = kpos;

*Symb = ci;

return flag;
}



#if IDF_CPP_ == 1
int idf_dbuf_str(int* Kstr, int* Kurs, long* Kpos,
                 unsigned int* Symb)
#else
int idf_dbuf_str(Kstr,Kurs,Kpos, Symb)
int *Kstr,*Kurs;
long *Kpos;
unsigned int *Symb;
#endif
{
/*
-----------------------------------
 put 'string' symbols (between "")
        into buffer
-----------------------------------

  0 - OK
 -1 - I/O error
 -2 - unexpected IDF_EOF
 -3 - unexpected IDF_EOS
 -4 - empty string
 -5 - too long string
 -6 - improper string case
*/

unsigned int ci,c,crpt;
int k,kk,blk,ibrk,jbrk,flag,kb,ks,ic;
int kstr,kurs;
long kpos;
register int i,j;
char *Buf;

Buf = Dbuf->Data_buf;
Buf[0] = IDF_EOS;

           i = 0;
Dbuf->length = i;

kpos = *Kpos;
kstr = *Kstr;
kurs = *Kurs;

Dbuf->kstr   = kstr;
Dbuf->kurs   = kurs;
Dbuf->pos    = kpos;

   ci = *Symb;
if(ci != '\"')
 {
  flag = -6;
  goto err;
 }

flag = 0;
  kk = 0;
  kb = 0;

crpt = ' ';
 blk = 0;

      ibrk=0;
while(ibrk==0)
 {

  if(i < IDF_DBUF_LENGTH)
   {
    if(blk)
     {
      /* return stored symbol */
      ci = crpt;
      blk=0;
     }

    else
     {
       /* get next symbol */
       k = idf_file_getc(&ci);
       kurs++;
       kpos++;

       if(k==2)
        {
	 /* I/O error */
	 flag = -1;
	 ibrk=1;
	}

       else if(k==1)
	{
	 /* unexpected IDF_EOF */
          kstr++;
          kurs=0;

	  flag = -2;
	  ibrk=1;
	 }

     }/*blk*/


        if(flag)
         {
          ;
         }
        else if(kb)
         {
          if(ci == '\"')
           {
            kb=0;
           }
          else
           {
            /*end of string */
            ibrk=2;
           }
         }

        else if(ci == '\"')
         {
          /* string is over */
          if(i)
           {
            kb++;
           }
          else
           {
            /* empty string*/
            flag = -4;
            ibrk=1;
           }
         }

        else if(ci == '\\')
         {
          /* escape sequence or concatenation mark */
          kk=0;
          ks=1;

                  jbrk=0;
            while(jbrk==0)
             {
              k = idf_file_getc(&c);
              kurs++;
              kpos++;

              if(k==2)
               {
	        /* I/O error */
	        flag = -1;
	        ibrk=1;
                jbrk=-1;
	       }

              else if(k==1)
	       {
	        /* unexpected IDF_EOF */
                kstr++;
                kurs=0;

	        flag = -2;
	        ibrk=1;
                jbrk=-1;
	       }

              else if(c == '\\')
	       {
	        if(kk)
                 {
                  /* // are separated by white space */
                  crpt = c;
                  blk=1;
                  jbrk=1;
                 }
                else
                 {
                  ks++;
                 }
	       }

              else if(c == ' ')
	       {
	        kk++;
	       }

              else if( idf_string_over(c) )
	       {
                /* IDF_EOS is found, skip symbols */
                kk=0;
                ks--;

	        kstr++;
                kurs=0;

                jbrk=1;
	       }

              else if(c == '\t')
	       {
                kk = kk + 8;
#if IDF_TXT_TABULATORS == 1
                kurs+=7;
#endif
	       }

              else if(c == '\v')
	       {
                if(kurs) kk = kk + (kurs-1)*8;
#if IDF_TXT_VERTICES == 1
	        kstr+=8;
#endif
	       }

              else if(c == '\"')
               {
                ic = ks/2;
                ic = ks - ic*2;
                if((kk>0)||(ic==0))
                 {
                  crpt = c;
                  blk=1;
                  jbrk=1;
                 }
                else
                 {
                  /* skip escape sequence \" */
                  jbrk=2;
                 }
               }

              else
               {
                /* escape sequence */
                jbrk=2;
               }

             }/*jbrk*/

            if(ks)
             {
              for(j=0;j<ks;j++)
               {
                if(i < IDF_DBUF_LENGTH)
                 {
                  Buf[i] = '\\';
                      i++;
                 }
                else
                 {
                  flag = -5;
                  ibrk=1;
                 }
               }
             }

            if(kk)
             {
              for(j=0;j<kk;j++)
               {
                if(i < IDF_DBUF_LENGTH)
                 {
                  Buf[i] = ' ';
                      i++;
                 }
                else
                 {
                  flag = -5;
                  ibrk=1;
                 }
               }
             }

            if(jbrk==2)
             {            
                if(i < IDF_DBUF_LENGTH)
                 {
                  Buf[i] = c;
                      i++;
                 }
                else
                 {
                  flag = -5;
                  ibrk=1;
                 }
             }

         } /* '/' */

        else if( idf_string_over(ci) )
	 {
	  kstr++;
          kurs=0;

          flag = -3;
          ibrk=1;
	 }

        else if(ci == '\t')
	 {
#if IDF_TXT_TABULATORS == 1
          kurs+=7;
#endif
              for(j=0;j<8;j++)
               {
                if(i < IDF_DBUF_LENGTH)
                 {
                  Buf[i] = ' ';
                      i++;
                 }
                else
                 {
                  flag = -5;
                  ibrk=1;
                 }
               }
	 }

        else if(ci == '\v')
	 {
#if IDF_TXT_VERTICES == 1
	  kstr+=8;
#endif
             jbrk = (kurs-1)*8;
          if(jbrk>0)
           {
              for(j=0;j<jbrk;j++)
               {
                if(i < IDF_DBUF_LENGTH)
                 {
                  Buf[i] = ' ';
                      i++;
                 }
                else
                 {
                  flag = -5;
                  ibrk=1;
                 }
               }
            }
	 }

        else
	 {
          /* any other symbol*/
          Buf[i] = ci;
              i++;
         }

   } /*i*/


  else
   {
    /* too large input */
    flag = -5;
    ibrk=1;
   }

 } /*ibrk*/


err:
if(flag < 0)
 {
  ibrk = 14 - flag;
  idf_txt_err_put(ibrk, kstr,kurs,kpos, ci);
 }

Buf[i] = IDF_EOS;
Dbuf->length = i;

*Kstr = kstr;
*Kurs = kurs;
*Kpos = kpos;

*Symb = ci;

return flag;
}



#if IDF_CPP_ == 1
int idf_dbuf_txt(int* Kstr, int* Kurs, long* Kpos,
                 unsigned int* Symb)
#else
int idf_dbuf_txt(Kstr,Kurs,Kpos, Symb)
int *Kstr,*Kurs;
long *Kpos;
unsigned int *Symb;
#endif
{
/*
-----------------------------------
 put 'string' symbols (between # #)
        into buffer
-----------------------------------

  0 - text unit done
 -1 - I/O error
 -2 - unexpected IDF_EOF
 -3 - empty text
 -4 - incorrect text
 -5 - improper trigraph case
 -6 - too long text input
 -7 - improper text case
*/

unsigned int c,ci;
int k,ibrk,flag,blk,kd;
int kstr,kurs;
long kpos;
register int i;
char *Buf;

Buf = Dbuf->Data_buf;
Buf[0] = IDF_EOS;

           i = 0;
Dbuf->length = i;

kpos = *Kpos;
kstr = *Kstr;
kurs = *Kurs;

Dbuf->kstr   = kstr;
Dbuf->kurs   = kurs;
Dbuf->pos    = kpos;


   ci = *Symb;
if(ci != '#')
 {
  flag = -7;
  goto err;
 }

flag = 0;
  kd = 0;

 blk = 0;
   c = ' ';

      ibrk=0;
while(ibrk==0)
 {
  if(i < IDF_DBUF_LENGTH)
   {
   if(blk)
    {
     /* return stored symbol */
     ci = c;
     blk=0;
    }
   else
    {
      /* get next symbol */
       k = idf_file_getc(&ci);
       kurs++;
       kpos++;

       if(k==2)
        {
	 /* I/O error */
	 flag = -1;
	 ibrk=1;
	}

       else if(k==1)
	{
	 /* unexpected IDF_EOF */
          kstr++;
          kurs=0;

	  flag = -2;
	  ibrk=1;
	 }
    }/*blk*/

        if(flag)
         {
          ;
         }

        else if(kd)
         {
          if(ci == '#')
           {
            /* concatenation case */
            kd=0;
           }
          else
           {
            /* text unit is over */
            ibrk=1;
           }
         }

        else if(ci == '?')
         {
          c = ci;

          k = idf_stxt_skip_trig(&kstr,&kurs,&kpos, &c);
          
          if(k == -1)
           {
	    /* I/O error */
	    flag = -1;
	    ibrk=1;
           }
          else if(k == -2)
	   {
	    /* unexpected IDF_EOF */
            flag = -2;
	    ibrk=1;
	   }
          else if(k == -3)
	   {
	    /* improper trigraph */
            flag = -5;
	    ibrk=1;
	   }
          else if(k==1)
	   {
            Buf[i] = '?';
                i++;
            blk=1;
           }
          else if(k==2)
	   {
            Buf[i] = '?';
                i++;
            Buf[i] = '?';
                i++;
            blk=1;
           }

          else
           {
            Buf[i] = c;
                i++;
           }
         }

        else if(ci == '#')
         {
          if(i==0)
           {
            flag = -3;
           }

          kd++;
         }

        else if( idf_string_over(ci) )
	 {
	  kstr++;
          kurs=0;
	 }

        else
         {
          Buf[i] = ci;
              i++;
         }

   } /*i*/


  else
   {
    /* too large input */
    flag = -6;
    ibrk=1;
   }

 } /*ibrk*/

err:
if(flag < 0)
 {
  ibrk = 20 - flag;
  idf_txt_err_put(ibrk, kstr,kurs,kpos, ci);
 }

Buf[i] = IDF_EOS;
Dbuf->length = i;

*Kstr = kstr;
*Kurs = kurs;
*Kpos = kpos;

*Symb = ci;

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_cmplx(int* Kstr, int* Kurs, long* Kpos,
                   unsigned int* Symb)
#else
int idf_dbuf_cmplx(Kstr,Kurs,Kpos, Symb)
int *Kstr,*Kurs;
long *Kpos;
unsigned int *Symb;
#endif
{
/*
-----------------------------------------
  skip 'complex unit' symbols between ()
-----------------------------------------

  0 - complex skipped
 -1 - I/O error
 -2 - unexpected IDF_EOF
 -3 - empty complex
 -4 - incorrect complex
 -5 - improper symbol in complex
 -6 - too long complex
 -7 - improper complex case
*/
unsigned int ci;
int k,ibrk,flag;
int kstr,kurs;
long kpos;
register int i;
int n,nn,id;
char *Buf;

Buf = Dbuf->Data_buf;
Buf[0] = IDF_EOS;

           i = 0;
Dbuf->length = i;

kpos = *Kpos;
kstr = *Kstr;
kurs = *Kurs;

Dbuf->kstr   = kstr;
Dbuf->kurs   = kurs;
Dbuf->pos    = kpos;

   ci = *Symb;
if(ci != '(')
 {
  flag = -7;
  goto err;
 }

flag = 0;
   n = 0;
  nn = 0;
  id = 0;

      ibrk=0;
while(ibrk==0)
 {
  if(i < IDF_DBUF_LENGTH)
   {
      /* get next symbol */
       k = idf_file_getc(&ci);
       kurs++;
       kpos++;

       if(k==2)
        {
	 /* I/O error */
	 flag = -1;
	 ibrk=1;
	}

       else if(k==1)
	{
	 /* unexpected IDF_EOF */
          kstr++;
          kurs=0;

	  flag = -2;
	  ibrk=1;
	 }

        else if(idf_string_over(ci))
	 {
	  kstr++;
          kurs=0;
	 }

        else if(ci == ')')
         {
          if(id != 1)
           {
            /* incorrect complex*/
            flag = -4;
            ibrk=1;
           }
          else if(n==0||nn==0)
           {
            /* empty complex*/
            flag = -3;
            ibrk=1;
           }
          else
           {
            /* complex is over */
            ci = ' ';
            ibrk=1;
           }
         }

        else if(ci == ',')
         {
          if(id)
           {
            /* incoorect complex number*/
            flag = -4;
            ibrk=1;
           }

          Buf[i] = ci;
              i++;
          id++;
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

        else if( idf_data_symbol(ci) )
         {
          if(id)
           {
            nn++;
           }
          else
           {
            n++;
           }

          Buf[i] = ci;
              i++;
         }

        else
         {
          /* improper symbol*/
          flag = -5;
          ibrk=1;
         }

   } /*i*/

  else
   {
    /* too large input */
    flag = -6;
    ibrk=1;
   }

 } /*ibrk*/

err:
if(flag < 0)
 {
  ibrk = 83 - flag;
  idf_txt_err_put(ibrk, kstr,kurs,kpos, ci);
 }

Buf[i] = IDF_EOS;
Dbuf->length = i;

*Kstr = kstr;
*Kurs = kurs;
*Kpos = kpos;

*Symb = ci;

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_add(unsigned int c)
#else
int idf_dbuf_add(c)
unsigned int c;
#endif
{
/*
--------------------------
 add character to buffer
--------------------------

  0 - OK
 -1 - buffer overflow
*/

int flag=0;
register int i;

   i = Dbuf->length;
if(i < IDF_DBUF_LENGTH)
 {
  Dbuf->Data_buf[i] = (char)c;
                 i++;
  Dbuf->Data_buf[i] = IDF_EOS;
  Dbuf->length = i;
 }

else
 {
  flag = -1;
 }

return flag;
}



#if IDF_CPP_ == 1
int idf_dbuf_to_d_(double* dp)
#else
int idf_dbuf_to_d_(dp)
double *dp;
#endif
{
int k;
double d;
char *buf;
char *endptr;
int flag=0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_DIGIT)
{
 /* improper number type */
 d = (double)0.0;
 flag=6;
}

else
{
 k = 0;
 buf = Dbuf->Data_buf;

 d = strtod(buf, &endptr);


 if(buf == endptr)
  {
   /* error in conversion */
   d = 0.0;
   flag=7;
  }
}

 *dp = d;

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_to_f_(float* fp)
#else
int idf_dbuf_to_f_(fp)
float *fp;
#endif
{
int k;
double dp;
char *buf;
char *endptr;
int flag=0;

 *fp = (float)0.0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_DIGIT)
{
 /* improper number type */
 flag=6;
}

else
{
 k = 0;
 buf = Dbuf->Data_buf;

 dp = strtod(buf, &endptr);

 if(buf == endptr)
  {
   /* error in conversion */
   flag=7;
  }
 else
  {
   if(dp > IDF_FLOAT_MAX)
    {
     flag=4;
    }
   else
    {
     *fp = (float)dp;
    }
  }
}

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_to_l_(long* lp)
#else
int idf_dbuf_to_l_(lp)
long *lp;
#endif
{
static int base=0; 

int k,flag;
long l;
char *buf;
char *endptr;

flag= 0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_DIGIT)
{
 /* improper number type */
 l = 0L;
 flag=6;
}

else
{
 buf = Dbuf->Data_buf;

 l = strtol(buf, &endptr ,base);

 if(buf == endptr)
  {
   /* error in conversion */
   l = 0L;
   flag=7;
  }
}

 *lp = l;

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_to_i_(int* lp)
#else
int idf_dbuf_to_i_(lp)
int *lp;
#endif
{
static long li = IDF_INT_MAX;
static int base=0; 

long l,lc;
int k,kk;
char *buf;
char *endptr;
int flag=0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_DIGIT)
{
 /* improper number type */
 l = 0;
 flag=6;
}

else
{
 buf = Dbuf->Data_buf;

 l = strtol(buf, &endptr, base);

 if(buf == endptr)
  {
   /* error in conversion */
   l = 0L;
   flag=7;
  }

 else
  {
   if(l < 0)
    {
     lc=-l;
     kk=1;
    }
   else
    {
     lc=l;
     kk=0;
    }
 
   if(lc > li)
    {
     if(kk) l = -li; else l=li;
     flag=4;
    }
  }
}
 
*lp = (int)l;

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_to_d(double* dp)
#else
int idf_dbuf_to_d(dp)
double *dp;
#endif
{
int k;
int flag=0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_DIGIT)
{
 /* improper number type */
 *dp = (double)0.0;
 flag=6;
}

else
{
 flag = idf_str_to_d(Dbuf->Data_buf, Dbuf->length, dp);
}

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_to_f(float* fp)
#else
int idf_dbuf_to_f(fp)
float *fp;
#endif
{
int k;
double dp;
int flag=0;

 *fp = (float)0.0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_DIGIT)
{
 /* improper number type */
 flag=6;
}

else
{
  k = idf_str_to_d(Dbuf->Data_buf, Dbuf->length, &dp);

if(k)
 {
  flag = k;
 }
else
 {
  if(dp > IDF_FLOAT_MAX)
   {
    flag=4;
   }
  else
   {
    *fp = (float)dp;
   }
 }
}

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_to_l(long* lp)
#else
int idf_dbuf_to_l(lp)
long *lp;
#endif
{
int k;
int flag=0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_DIGIT)
{
 /* improper number type */
 *lp = 0L;
 flag=6;
}

else
{
 flag = idf_str_to_l(Dbuf->Data_buf, Dbuf->length, lp);
}

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_to_i(int* lp)
#else
int idf_dbuf_to_i(lp)
int *lp;
#endif
{
static long li = IDF_INT_MAX;

long l,lc;
int k,kk;
int flag=0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_DIGIT)
{
 /* improper number type */
 l = 0;
 flag=6;
}

else
{
  k = idf_str_to_l(Dbuf->Data_buf, Dbuf->length, &l);

if(k)
 {
  flag=k;
 }
else
 {
  if(l < 0)
   {
    lc=-l;
    kk=1;
   }
  else
   {
    lc=l;
    kk=0;
   }
 
  if(lc > li)
   {
    if(kk) l = -li; else l=li;
    flag=4;
   }
 }
}
 
*lp = (int)l;

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_to_t(char* s, int* ns)
#else
int idf_dbuf_to_t(s,ns)
char *s;
int  *ns;
#endif
{
int k;
register int i;
char *buf;
int nb;
int flag=0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_T)
{
 /* improper number type */
 s[0] = IDF_EOS;
 *ns = 0;
 flag=2;
}

else
{
  buf = Dbuf->Data_buf;
   nb = Dbuf->length;    
if(nb)
 {
     i=0;
    while(i<nb)
     {
      s[i] = buf[i];
        i++;
     }
 s[i] = IDF_EOS;
 *ns = nb;
 }

else
 {
  /* no data symbols*/
  s[0] = IDF_EOS;
  *ns = 0;
  flag=1;
 }
}

return flag;
}



#if IDF_CPP_ == 1
int idf_dbuf_to_c(char* c)
#else
int idf_dbuf_to_c(c)
char *c;
#endif
{
int k;
int flag=0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_C)
 {
  /* improper number type */
  *c = IDF_EOS;
  flag=6;
 }

else
 {
  flag = idf_str_to_c(Dbuf->Data_buf, Dbuf->length, c);
 }

return flag;
}



#if IDF_CPP_ == 1
int idf_dbuf_to_s(char* str, int* ns)
#else
int idf_dbuf_to_s(str,ns)
char *str;
int  *ns;
#endif
{
int k;
int flag=0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_S)
 {
  /* improper number type */
  str[0] = IDF_EOS;
  *ns = 0;
  flag=6;
 }

else
 {
  flag = idf_str_to_s(Dbuf->Data_buf, Dbuf->length, str,IDF_DBUF_LENGTH);
  *ns = strlen(str);
 }

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_to_z(float* cf)
#else
int idf_dbuf_to_z(cf)
float *cf;
#endif
{
int k;
int flag;
double cd[2];
double d;

flag = 0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_CMPLX)
{
 /* improper number type */
 cf[0] = (float)0.0;
 cf[1] = (float)0.0;
 flag=6;
}

else
{
    k = idf_str_to_z(Dbuf->Data_buf, Dbuf->length, cd);
 if(k)
  {
   flag = k;
  }

 else
  {
      d = cd[0];
   if(d > IDF_FLOAT_MAX)
    {
     flag=4;
    }
   else
    {
     cf[0] = (float)d;
    }

      d = cd[1];
   if(d > IDF_FLOAT_MAX)
    {
     flag=4;
    }
   else
    {
     cf[1] = (float)d;
    }
  }
}

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_to_w(double* cd)
#else
int idf_dbuf_to_w(cd)
double *cd;
#endif
{
int k;
int flag=0;

   k = Dbuf->jtype;
if(k != IDF_FRMT_CMPLX)
 {
  /* improper number type */
  cd[0] = (double)0.0;
  cd[1] = (double)0.0;
  flag=6;
 }

else
 {
  flag = idf_str_to_z(Dbuf->Data_buf, Dbuf->length, cd);
 }

return flag;
}


#if IDF_CPP_ == 1
int idf_dbuf_to_v(double* cd, int* typ)
#else
int idf_dbuf_to_v(cd,typ)
double *cd;
int    *typ;
#endif
{
int k;
int flag=0;

   k = Dbuf->jtype;
if(k==IDF_FRMT_C||k==IDF_FRMT_S||k==IDF_FRMT_T)
 {
  /* improper number type */
  cd[0] = (double)0.0;
  cd[1] = (double)0.0;
  *typ=0;
  flag=6;
 }

else
 {
  flag = idf_str_to_v(Dbuf->Data_buf, Dbuf->length, cd, typ);

  /*
  printf("buf=%s val=%e %e typ=%d\n", Dbuf->Data_buf,cd[0],cd[1],*typ);
  */
 }

return flag;
}

