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
int idf_skip_data(int* Nstr, int* Nkurs, long* Npos,
                  unsigned int* Symb)
#else
int idf_skip_data(Nstr,Nkurs,Npos, Symb)
int  *Nstr,*Nkurs;
long *Npos;
unsigned int *Symb;
#endif
{
/*
====================================
  skip symbols belong to data field

  flag =  1 - empty data with EOF flag
          0 - OK
         <0 - error
====================================
*/
    unsigned int c,ci,cerr,crpt;
    int flag,ibrk,blk,is,n,k,jerr;
    int kstr,kurs;
    long pos;
    int Kstr,Kurs;
    long Kpos;

    kstr  = *Nstr;
    kurs  = *Nkurs;
     pos  = *Npos;

    flag  = 0;

    jerr  = 0;
    cerr  = ' ';

ci = *Symb;
if( !idf_name_seps(ci) )
 {
  jerr=7;
  cerr=ci;
  goto err;
 }

     k    = 0;
    is    = 0;
     n    = 0;

     crpt = ' ';
     blk  = 0;

       ibrk=0;
 while(ibrk==0)
   {

    /*--------------------
      get the next symbol 
     ---------------------*/

   if(blk)
    {
     /* stored symbol up */
     ci = crpt;
     blk=0;
    }

   else
    {
      k = idf_file_getc(&ci);
      kurs++;
      pos++;

      if(k==2)
       {
        /* I/O error */
        jerr = 1;
        flag = -1;
	ibrk=2;
       }

      else if(k==1)
       {
        kstr++;
        kurs=0;

        if(n)
         {
	  /* unexpected EOF */
          jerr = 2;
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* no data & EoF */
          flag = 1;
          ibrk=1;
         }
       }

    }/*blk*/


    /*-------------
       EOS case
     --------------*/

      if(flag)
       {
        ;
       }

      else if( idf_string_over(ci) )
       {
        /* file string is over*/
	kstr++;
        kurs=0;
       }


    /*-------------
      comment case
     --------------*/

      else if(ci == '/')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_cmnt(&Kstr,&Kurs,&Kpos, &c);

        if(k ==  2)
         {
          /* improper symbol '/' */
          jerr = 3;
          cerr = ci;
          flag = -1;
          ibrk=1;
         }
        else if(k ==  1)
         {
          /* unexpected EOF */
          jerr = 2;
          flag = -1;
          ibrk=1;
         }
        else if(k < 0)
         {
          /* comment error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          blk=1;
          crpt=c;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* '/' */


    /*-------------
      string  case 
     --------------*/

      else if(ci == '\"')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_str(&Kstr,&Kurs,&Kpos, &c);

        if(k < 0)
         {
          /* string error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* string skipped */
          blk=1;
          crpt = c;
          n++;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* "" */


    /*---------------
      character case 
     ----------------*/

      else if(ci == '\'')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_char(&Kstr,&Kurs,&Kpos, &c);

        if(k < 0)
         {
          /* char error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* char skipped */
          blk=1;
          crpt=c;
          n++;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* '' */


    /*---------------
      formular case 
     ----------------*/

      else if(ci == '$')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_form(&Kstr,&Kurs,&Kpos, &c);

        if(k == -1)
         {
          /* formular error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* formular skipped */
          blk=1;
          crpt=c;
          n++;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* $$ */


    /*---------------
        text case 
     ----------------*/

      else if(ci == '#')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_txt(&Kstr,&Kurs,&Kpos, &c);

        if(k < 0)
         {
          /* text error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* text unit skipped */
          blk=1;
          crpt=c;
          n++;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* # # */


    /*---------------
      complex case 
     ----------------*/

      else if(ci == '(')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_cmplx(&Kstr,&Kurs,&Kpos, &c);

        if(k < 0)
         {
          /* complex error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* complex skipped */
          blk=1;
          crpt=c;
          n++;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* () */


       /*---------------------------
         search for data separator 
        ---------------------------*/

      else if(ci == ';')
       {
        if(is)
         {
          /* unclosed { */
          jerr = 4;
          cerr = ci;
          flag = -1;
	  ibrk=1;
         }
        else if(n)
	 {
          /* data done */
          flag    = 0;
	  ibrk=1;
	 }
        else
	 {
          /* empty data */
          jerr = 5;
          cerr = ci;
          flag = -1;
	  ibrk=1;
	 }

       } /* separator */


      /*------------------
         improper symbols
        ------------------*/

      else if(ci == '='||ci == ':')
       {
        cerr = ci;
        jerr = 3;
        flag = -1;
	ibrk=1;
       }


      /*----------
         case {}
        ----------*/

      else if(ci == '{' )
       {
        is++;
       }

      else if(ci == '}' )
       {
        if(is<1)
         {
          /* extra } */
          jerr = 6;
          flag = -1;
	  ibrk=1;
         }
        is--;
       }


      /*----------------
         data attributes
        ----------------*/

      else if(ci == ','||ci=='*'||ci==' ')
       {
        ;
       }


      /*-----------------
           tabulators
       -----------------*/

      else if(ci == '\t')
       {
#if IDF_TXT_TABULATORS == 1
        kurs+=7;
#else
        ;
#endif
       }


       /*-----------
           vertices  
         -----------*/

      else if(ci == '\v')
       {
#if IDF_TXT_VERTICES == 1
        kstr+=8;
#else
        ;
#endif
       }


       /*-------------------------------
         add ordinary character to data
        -------------------------------- */

      else
       {
        n++;
       }


   } /* end of search */

err:
*Nstr  = kstr;
*Nkurs = kurs;
*Npos  = pos;

*Symb = ci;

  if(jerr)
   {
    ibrk = jerr + 27;
    idf_txt_err_put(ibrk, kstr,kurs,pos, cerr);
   }

 return flag;
}



#if IDF_CPP_ == 1
int idf_skip_data_period(int* Nstr, int* Nkurs, long* Npos,
                         unsigned int* Symb)
#else
int idf_skip_data_period(Nstr,Nkurs,Npos, Symb)
int  *Nstr,*Nkurs;
long *Npos;
unsigned int *Symb;
#endif
{
/*
====================================
  skip data units inside the period

          0 - OK
         <0 - error
====================================
*/
    unsigned int c,ci,cerr,crpt;
    int flag,ibrk,blk,is,n,k,jerr;
    int kstr,kurs;
    long pos;
    int Kstr,Kurs;
    long Kpos;

    kstr  = *Nstr;
    kurs  = *Nkurs;
     pos  = *Npos;

    flag  = 0;

    jerr  = 0;
    cerr  = ' ';

   ci = *Symb;
if(ci != '{')
 {
  jerr=8;
  cerr=ci;
  goto err;
 }

     k    = 0;
    is    = 1;
     n    = 0;

     crpt = ' ';
     blk  = 0;

       ibrk=0;
 while(ibrk==0)
   {

    /*--------------------
      get the next symbol 
     ---------------------*/

   if(blk)
    {
     /* stored symbol up */
     ci = crpt;
     blk=0;
    }

   else
    {
      k = idf_file_getc(&ci);
      kurs++;
      pos++;

      if(k==2)
       {
        /* I/O error */
        jerr = 1;
        flag = -1;
	ibrk=2;
       }

      else if(k==1)
       {
        kstr++;
        kurs=0;

        if(n)
         {
	  /* unexpected EOF */
          jerr = 2;
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* no data & EoF */
          flag = 1;
          ibrk=1;
         }
       }

    }/*blk*/


    /*-------------
       EOS case
     --------------*/

      if(flag)
       {
        ;
       }

      else if( idf_string_over(ci) )
       {
        /* file string is over*/
	kstr++;
        kurs=0;
       }


    /*-------------
      comment case
     --------------*/

      else if(ci == '/')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_cmnt(&Kstr,&Kurs,&Kpos, &c);

        if(k ==  2)
         {
          /* improper symbol '/' */
          jerr = 3;
          cerr = ci;
          flag = -1;
          ibrk=1;
         }
        else if(k ==  1)
         {
          /* unexpected EOF */
          jerr = 2;
          flag = -1;
          ibrk=1;
         }
        else if(k < 0)
         {
          /* comment error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          blk=1;
          crpt=c;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* '/' */


    /*-------------
      string  case 
     --------------*/

      else if(ci == '\"')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_str(&Kstr,&Kurs,&Kpos, &c);

        if(k < 0)
         {
          /* string error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* string skipped */
          blk=1;
          crpt = c;
          n++;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* "" */


    /*---------------
      character case 
     ----------------*/

      else if(ci == '\'')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_char(&Kstr,&Kurs,&Kpos, &c);

        if(k < 0)
         {
          /* char error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* char skipped */
          blk=1;
          crpt=c;
          n++;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* '' */


    /*---------------
      formular case 
     ----------------*/

      else if(ci == '$')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_form(&Kstr,&Kurs,&Kpos, &c);

        if(k == -1)
         {
          /* formular error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* formular skipped */
          blk=1;
          crpt=c;
          n++;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* $$ */


    /*---------------
        text case 
     ----------------*/

      else if(ci == '#')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_txt(&Kstr,&Kurs,&Kpos, &c);

        if(k < 0)
         {
          /* text error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* text unit skipped */
          blk=1;
          crpt=c;
          n++;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* # # */


    /*---------------
      complex case 
     ----------------*/

      else if(ci == '(')
       {
        Kstr = kstr;
        Kurs = kurs;
        Kpos = pos;
           c = ci;

        k = idf_stxt_skip_cmplx(&Kstr,&Kurs,&Kpos, &c);

        if(k < 0)
         {
          /* complex error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* complex skipped */
          blk=1;
          crpt=c;
          n++;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* () */


       /*---------------------------
         search for data separator 
        ---------------------------*/

      else if(ci == ';')
       {
        if(is)
         {
          /* unclosed { */
          jerr = 4;
          cerr = ci;
          flag = -1;
	  ibrk=1;
         }
        else
	 {
          /* unexpected end of period */
          jerr = 7;
          cerr = ci;
          flag = -1;
	  ibrk=1;
	 }

       } /* separator */


      /*------------------
         improper symbols
        ------------------*/

      else if(ci == '='||ci == ':')
       {
        cerr = ci;
        jerr = 3;
        flag = -1;
	ibrk=1;
       }


      /*----------
         case {
        ----------*/

      else if(ci == '{' )
       {
        is++;
       }


      /*----------
         case }
        ----------*/

      else if(ci == '}' )
       {
           is--;
        if(is==0)
         {
          if(n)
           {
            /* data period finished */
            flag=0;
            ibrk=2;
           }
          else
           {
            /* empty data period */
            jerr = 5;
            cerr = ci;
            flag = -1;
	    ibrk=1;
           }
         }
        else if(is<0)
         {
          /* extra } */
          jerr = 6;
          flag = -1;
	  ibrk=1;
         }
       }


      /*----------------
         data attributes
        ----------------*/

      else if(ci == ','||ci=='*'||ci==' ')
       {
        ;
       }


      /*-----------------
           tabulators
       -----------------*/

      else if(ci == '\t')
       {
#if IDF_TXT_TABULATORS == 1
        kurs+=7;
#else
        ;
#endif
       }


       /*-----------
           vertices  
         -----------*/

      else if(ci == '\v')
       {
#if IDF_TXT_VERTICES == 1
        kstr+=8;
#else
        ;
#endif
       }


       /*-------------------------------
         add ordinary character to data
        -------------------------------- */

      else
       {
        n++;
       }


   } /* end of search */


err:
*Nstr  = kstr;
*Nkurs = kurs;
*Npos  = pos;

*Symb = ci;

  if(jerr)
   {
    ibrk = jerr + 34;
    idf_txt_err_put(ibrk, kstr,kurs,pos, cerr);
   }

 return flag;
}
