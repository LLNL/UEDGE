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
int idf_stxt_skip_cmnt(int* Kstr, int* Kurs, long* Pos,
                       unsigned int* Symb)
#else
int idf_stxt_skip_cmnt(Kstr,Kurs,Pos, Symb)
int *Kstr,*Kurs;
long *Pos;
unsigned int *Symb;
#endif
{
/*
----------------
  skip comments
----------------

  2 - it is not a comment (additional symbol is in Symb)
  1 - comment finished with eof
  0 - comment skipped
 -1 - I/O error
 -2 - unexpected EOF
 -3 - improper comment case
*/
unsigned int c,ci;
int k,kk,ibrk,flag;
int kstr,kurs;
long pos;

kstr = *Kstr;
kurs = *Kurs;
 pos = *Pos;

   c = *Symb;
if(c != '/')
 {
  /* improper comment case */
  flag = -3;
  goto err;
 }

flag = 0;

      /*----------------
        get next symbol 
       -----------------*/

           k = idf_file_getc(&c);
           kurs++;
           pos++;

           if(k==2)
            {
	     /* I/O error */
	     flag = -1;
	    }

           else if(k==1)
	    {
	     /* unexpected EOF */
              kstr++;
              kurs=0;

	      flag = -2;
	     }


            /*-----------
               EoS case 
              -----------*/

            else if( idf_string_over(c) )
	     {
              /* it is not a comment */
	      kstr++;
              kurs=0;

              flag=2;
	     }


	     /*-------------
                 case: //   
              --------------*/

	     else if(c == '/')
		  {
		   /* skip symbols upto the end-of-string or EOF */

                         ibrk=0;
		   while(ibrk==0)
		      {
		       k = idf_file_getc(&ci);
                       kurs++;
                       pos++;

		       if(k==2)
			 {
			  /* I/O error */
			  flag = -1;
                          ibrk=2;
			 }

		       else if(k==1)
			 {
			  kstr++;
                          kurs=0;
                          flag = 1;
                          ibrk=1;
			 }

		       else if(idf_string_over(ci))
			 {
			  kstr++;
                          kurs=0;
                          flag = 0;
			  ibrk=1;
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

		      }/*ibrk*/
		  }


	     /* ----------------
                   case: */ /* 
                ----------------*/

	     else if (c == '*')
		  {
		   /* skip symbols upto the */

                   ibrk=0;
                   kk=0;
		   while(ibrk==0)
		      {
		       k = idf_file_getc(&ci);
                       kurs++;
                       pos++;

		       if(k==2)
			 {
			  /* I/O error */
			  flag = -1;
                          ibrk=1;
			 }

		       else if(k==1)
			 {
			  kstr++;
                          kurs=0;

			  flag = -2;
			  ibrk=1;
			 }

		       else if( idf_string_over(ci) )
			 {
			  kstr++;
                          kurs=0;
			 }

		       else if(kk==0)
			 {
			  if(ci == '*') kk=1;
			 }

		       else if(kk==1)
			 {
			   if(ci == '/')
			     {
			      ibrk=1;
			     }
			   else
			     {
			      kk=0;
			     }
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

		      }/*ibrk*/
		  }


             /* ------------------
                  not a comment 
                ------------------*/

             else if(c == '\t')
	          {
#if IDF_TXT_TABULATORS == 1
                   kurs+=7;
#endif
                   flag=2;
	          }

             else if(c == '\v')
	          {
#if IDF_TXT_VERTICES == 1
	           kstr+=8;
#endif
                   flag=2;
	          }

	     else
		  {
		   /* not a comment */
                   flag = 2;
		  }

err:
*Kstr = kstr;
*Kurs = kurs;
*Pos  = pos;

/* replace comment by single white space */
if(flag < 0)
 {
  ibrk = -flag;
  idf_txt_err_put(ibrk, kstr,kurs,pos, c);

  *Symb = c;
 }
else if(flag==2)
 {
  *Symb = c;
 }
else
 {
  *Symb = ' ';
 }

return flag;
}



#if IDF_CPP_ == 1
int idf_stxt_skip_cmplx(int* Kstr, int* Kurs, long* Pos,
                        unsigned int* Symb)
#else
int idf_stxt_skip_cmplx(Kstr,Kurs,Pos, Symb)
int  *Kstr,*Kurs;
long *Pos;
unsigned int *Symb;
#endif
{
/*
-----------------------------------------
  skip 'complex unit' symbols between ()
-----------------------------------------

  0 - complex skipped
 -1 - I/O error
 -2 - unexpected EOF
 -3 - empty complex
 -4 - incorrect complex
 -5 - improper symbol in complex
 -6 - too long complex
 -7 - improper complex case
*/
unsigned int ci;
int k,ibrk,flag;
int kstr,kurs;
long pos;
int n,nn,nk,id;

kstr = *Kstr;
kurs = *Kurs;
 pos = *Pos;

   ci = *Symb;
if(ci != '(')
 {
  flag = -7;
  goto err;
 }

flag = 0;
   n = 0;
  nn = 0;
  nk = 0;
  id = 0;

      ibrk=0;
while(ibrk==0)
 {
  if(nk < IDF_DBUF_LENGTH)
   {
      /* get next symbol */
       k = idf_file_getc(&ci);
       kurs++;
       pos++;

       if(k==2)
        {
	 /* I/O error */
	 flag = -1;
	 ibrk=1;
	}

       else if(k==1)
	{
	 /* unexpected EOF */
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
          id++;
          nk++;
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
          nk++;
         }

        else
         {
          /* improper symbol*/
          flag = -5;
          ibrk=1;
         }
   }

  else
   {
    /* too long input */
    flag = -6;
    ibrk=1;
   }/*nk*/

 } /*ibrk*/

err:
if(flag < 0)
 {
  ibrk = 83 - flag;
  idf_txt_err_put(ibrk, kstr,kurs,pos, ci);
 }

*Kstr = kstr;
*Kurs = kurs;
*Pos  = pos;

*Symb = ci;

return flag;
}



#if IDF_CPP_ == 1
int idf_stxt_skip_form(int* Kstr, int* Kurs, long* Pos,
                       unsigned int* Symb)
#else
int idf_stxt_skip_form(Kstr,Kurs,Pos, Symb)
int  *Kstr,*Kurs;
long *Pos;
unsigned int *Symb;
#endif
{
/*
-------------------------------------
  skip 'formular' symbols between $$
   and concatenation if appropriate
-------------------------------------

  0 - formular skipped
 -1 - I/O error
 -2 - unexpected EOF
 -3 - empty formular
 -4 - incorrect formular
 -5 - improper formular case
*/
unsigned int ci;
int k,ibrk,flag,is;
int kstr,kurs;
long pos;
int n,kd;

kstr = *Kstr;
kurs = *Kurs;
 pos = *Pos;

   ci = *Symb;
if(ci != '$')
 {
  flag = -5;
  goto err;
 }

flag = 0;
   n = 0;
  kd = 0;
  is = 0;

      ibrk=0;
while(ibrk==0)
 {
      /* get next symbol */
       k = idf_file_getc(&ci);
       kurs++;
       pos++;

       if(k==2)
        {
	 /* I/O error */
	 flag = -1;
	 ibrk=1;
	}

       else if(k==1)
	{
	 /* unexpected EOF */
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

        else if(kd)
         {
          if(ci == '$')
           {
            /* concatenation case */
            kd=0;
           }

          else
           {
            /* formular is over */
            if(n)
             {
              if(is)
               {
                flag = -4;
               }
              else
               {
                flag = 0;
               }
             }
            else
             {
              flag = -3;
             }

            ibrk=1;
           }
         }

        else if(ci == '$')
         {
          kd++;
         }

        else if(ci == '(')
         {
          is++;
         }

        else if(ci == ')')
         {
          is--;
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

        else
         {
          n++;
         }

 } /*ibrk*/

err:
if(flag < 0)
 {
  ibrk = 3 - flag;
  idf_txt_err_put(ibrk, kstr,kurs,pos, ci);
 }

*Kstr = kstr;
*Kurs = kurs;
*Pos  = pos;

*Symb = ci;

return flag;
}



#if IDF_CPP_ == 1
int idf_stxt_skip_char(int* Kstr, int* Kurs, long* Pos,
                       unsigned int* Symb)
#else
int idf_stxt_skip_char(Kstr,Kurs,Pos, Symb)
int  *Kstr,*Kurs;
long *Pos;
unsigned int *Symb;
#endif
{
/*
-------------------------------------
  skip 'character' symbols between ''
   and concatenation if appropriate
-------------------------------------

  0 - characters skipped
 -1 - I/O error
 -2 - unexpected EOF
 -3 - improper symbol in character list
 -4 - improper char
 -5 - empty char
 -6 - improper char case
*/

unsigned int ci;
int k,ibrk,flag;
int kc,is,io;
int kstr,kurs;
long pos;

kstr = *Kstr;
kurs = *Kurs;
 pos = *Pos;

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
       pos++;

       if(k==2)
        {
	 /* I/O error */
	 flag = -1;
         ibrk=2;
	}

       else if(k==1)
	{
	 /* unexpected EOF */
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
          if(io)
           {
            /* end of char */
            kc++;
           }
          else if(is==1)
           {
            /* escape sequence */
            io++;
            is--;
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
           }
          else
           {
            is++;
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
            if(io>4)
             {
              flag = -4;
              ibrk=2;
             }
            else
             {
              io++;
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
             }
           }
         }

  } /* ibrk*/


err:
if(flag < 0)
 {
  ibrk = 8 - flag;
  idf_txt_err_put(ibrk, kstr,kurs,pos, ci);
 }

*Kstr = kstr;
*Kurs = kurs;
*Pos  = pos;

*Symb = ci;

return flag;
}



#if IDF_CPP_ == 1
int idf_stxt_skip_str(int* Kstr, int* Kurs, long* Pos,
                      unsigned int* Symb)
#else
int idf_stxt_skip_str(Kstr,Kurs,Pos, Symb)
int  *Kstr,*Kurs;
long *Pos;
unsigned int *Symb;
#endif
{
/*
-----------------------------------
  skip 'string' symbols between ""
  and concatenation if appropriate
-----------------------------------

  0 - string skipped
 -1 - I/O error
 -2 - unexpected EOF
 -3 - unexpected EOS
 -4 - empty string
 -5 - too long string
 -6 - improper string case
*/

unsigned int ci,c,crpt;
int k,kk,blk,ibrk,jbrk,flag,kb,ks,n,ic;
int kstr,kurs;
long pos;

kstr = *Kstr;
kurs = *Kurs;
 pos = *Pos;

   ci = *Symb;
if(ci != '\"')
 {
  flag = -6;
  goto err;
 }

flag = 0;
   n = 0;
  kk = 0;
  kb = 0;

crpt = ' ';
 blk = 0;

      ibrk=0;
while(ibrk==0)
 {

  if(n < IDF_DBUF_LENGTH)
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
       pos++;

       if(k==2)
        {
	 /* I/O error */
	 flag = -1;
	 ibrk=1;
	}

       else if(k==1)
	{
	 /* unexpected EOF */
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
          if(n)
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
              pos++;

              if(k==2)
               {
	        /* I/O error */
	        flag = -1;
	        ibrk=1;
                jbrk=1;
	       }

              else if(k==1)
	       {
	        /* unexpected EOF */
                kstr++;
                kurs=0;

	        flag = -2;
	        ibrk=1;
                jbrk=1;
	       }

              else if(c == '\\')
	       {
	        if(kk)
                 {
                  /* \\ are separated by white space */
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
                /* EOS is found, skip symbols */
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
                  kk++;
                  jbrk=1;
                 }
                else
                 {
                  /* skip escape sequence \" */
                  kk++;
                  jbrk=2;
                 }
               }

              else
               {
                /* return all symbols */
                kk++;
                jbrk=2;
               }

             }/*jbrk*/
            
            n = n + kk + ks;

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
          n = n + 8;
#if IDF_TXT_TABULATORS == 1
          kurs+=7;
#endif
	 }

        else if(ci == '\v')
	 {
          if(kurs) n = n + (kurs-1)*8;
#if IDF_TXT_VERTICES == 1
	  kstr+=8;
#endif
	 }

        else
	 {
          n++;
         }

   } /*n*/


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
  idf_txt_err_put(ibrk, kstr,kurs,pos, ci);
 }

*Kstr = kstr;
*Kurs = kurs;
*Pos  = pos;

*Symb = ci;

return flag;
}



#if IDF_CPP_ == 1
int idf_stxt_skip_crpt(int* Kstr, int* Kurs, long* Pos,
                       unsigned int* Symb)
#else
int idf_stxt_skip_crpt(Kstr,Kurs,Pos, Symb)
int  *Kstr,*Kurs;
long *Pos;
unsigned int *Symb;
#endif
{
/*
 1 - if Symb is doubled in a stream
 0 - else (next symbol is stored in Symb)
<0 - error
*/

unsigned int c,ci;
int k,flag;

c = *Symb;

       k = idf_file_getc(&ci);

       *Kurs = *Kurs + 1;
       *Pos  = *Pos + 1;

       if(k==2)
        {
	 /* I/O error */
	 flag = -1;
	}

       else if(k==1)
	{
	 /* unexpected EOF */
          *Kstr = *Kstr + 1;
          *Kurs = 0;

	  flag = -2;
	 }

        else if(ci == c)
         {
          flag = 1;
         }

        else
         {
          flag = 0;
         }

*Symb = ci;

return flag;
}



#if IDF_CPP_ == 1
int idf_stxt_skip_trig(int* Kstr, int* Kurs, long* Pos,
                       unsigned int* Symb)
#else
int idf_stxt_skip_trig(Kstr,Kurs,Pos, Symb)
int *Kstr,*Kurs;
long *Pos;
unsigned int *Symb;
#endif
{
/*
------------------------------------------
   replace trigraphs symbols ??xxx
------------------------------------------

  2 - not a trigraph: ??x
  1 - not a trigraph: ?x
  0 - trigraph replaced (result in Symb)
 -1 - I/O error
 -2 - unexpected EOF
 -3 - improper trigraph case
*/

unsigned int c,ci;
int k,flag;
int kstr,kurs;
long pos;

kstr = *Kstr;
kurs = *Kurs;
 pos = *Pos;

   ci = *Symb;
if(ci != '<')
 {
  c = ci;
  flag = -3;
  goto err;
 }

flag = 0; 

 /* get next symbol */
  k = idf_file_getc(&ci);
  kurs++;
  pos++;

  if(k==2)
   {
    /* I/O error */
    c = ci;
    flag = -1;
   }

  else if(k==1)
   {
    /* unexpected EOF */
    kstr++;
    kurs=0;

    c = ci;
    flag = -2;
   }

  else if(ci == '?')
   {
    /* double ? */

      /* get next symbol */
       k = idf_file_getc(&ci);
       kurs++;
       pos++;

       if(k==2)
        {
	 /* I/O error */
         c = ci;
	 flag = -1;
	}

       else if(k==1)
	{
	 /* unexpected EOF */
          kstr++;
          kurs=0;

          c = ci;
	  flag = -2;
	}

       /* trigraphs */

       else if(ci == '=')
	{
         c = '#';
        }
       else if(ci == '(')
	{
         c = '[';
        }
       else if(ci == '/')
	{
         c = '\\';
        }
       else if(ci == ')')
	{
         c = ']';
        }
       else if(ci == '\'')
	{
         c = '^';
        }
       else if(ci == '<')
	{
         c = '{';
        }
       else if(ci == '!')
	{
         c = '|';
        }
       else if(ci == '>')
	{
         c = '}';
        }
       else if(ci == '-')
	{
         c = '~';
        }

       else
	{
         /* double ? */
         c = ci;
         flag=2;
        }
   }

  else
   {
    /* single ? */
    c = ci;
    flag=1;
   }


err:
*Kstr = kstr;
*Kurs = kurs;
*Pos  = pos;

*Symb = c;

return flag;
}



#if IDF_CPP_ == 1
int idf_stxt_skip_txt(int* Kstr, int* Kurs, long* Pos,
                      unsigned int* Symb)
#else
int idf_stxt_skip_txt(Kstr,Kurs,Pos, Symb)
int *Kstr,*Kurs;
long *Pos;
unsigned int *Symb;
#endif
{
/*
-------------------------------------
  skip 'text' symbols between # #
   and concatenation if appropriate
-------------------------------------

  0 - text unit skipped
 -1 - I/O error
 -2 - unexpected EOF
 -3 - empty text
 -4 - incorrect text
 -5 - improper trigraph case
 -6 - too long text input
 -7 - improper text case
*/

unsigned int c,ci;
int k,ibrk,flag,n,blk,kd;
int kstr,kurs;
long pos;

kstr = *Kstr;
kurs = *Kurs;
 pos = *Pos;

   ci = *Symb;
if(ci != '#')
 {
  flag = -7;
  goto err;
 }

flag = 0;
   n = 0;
  kd = 0;

 blk = 0;
   c = ' ';

      ibrk=0;
while(ibrk==0)
 {
  if(n < IDF_DBUF_LENGTH)
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
       pos++;

       if(k==2)
        {
	 /* I/O error */
	 flag = -1;
	 ibrk=1;
	}

       else if(k==1)
	{
	 /* unexpected EOF */
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

          k = idf_stxt_skip_trig(&kstr,&kurs,&pos, &c);
          
          if(k == -1)
           {
	    /* I/O error */
	    flag = -1;
	    ibrk=1;
           }
          else if(k == -2)
	   {
	    /* unexpected EOF */
            flag = -2;
	    ibrk=1;
	   }
          else if(k == -3)
	   {
	    /* improper trigraph */
            flag = -5;
	    ibrk=1;
	   }
          else if(k>0)
	   {
            n = n + k;
            blk=1;
           }
          else
           {
            n++;
           }
         }

        else if(ci == '#')
         {
          if(n==0)
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
          n++;
         }

   }/* n*/

  else
   {
    /* too long text */
    flag = -6;
    ibrk=1;
   }

 } /*ibrk*/

err:
if(flag < 0)
 {
  ibrk = 20 - flag;
  idf_txt_err_put(ibrk, kstr,kurs,pos, ci);
 }

*Kstr = kstr;
*Kurs = kurs;
*Pos  = pos;

*Symb = ci;

return flag;
}

