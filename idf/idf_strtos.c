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
int idf_str_to_s(char* str, int n,
                 char* sval, int nval)
#else
int idf_str_to_s(str,n, sval,nval)
char *str;
int   n;
char *sval;
int   nval;
#endif
{
/*
 -----------------------------------
  converts string to another string
     according to C standards
 -----------------------------------
returns:
  0  OK
  1  no data symbols
  2  unknown escape sequence
  3  improper symbol
  4  empty escape sequence
  5  improper string
*/
static int icmax = IDF_CHAR_MAX;
static int hexbase = 16;
static int octbase = 8;

int flag,ic,iv,k,ibrk;
char c,cc;
register int i,j;

  flag = 0;
     j = 0;
sval[j] = IDF_EOS;

if(nval < 1)
 {
  /* nothing to do */
  goto fin;
 }

if(n < 1)
 {
  /* no data symbols */
  flag=1;
  goto fin;
 }


i = 0;
c = str[i];
        i++;

if(n == 1)
 {
  /* single char */
  if(c == '\\')
   {
    /*empty escape sequence */
    flag=4;
   }
  else
   {
    sval[j] = c;
         j++;
    sval[j] = IDF_EOS;
   }

  goto fin; 
 }


while(i<=n)
{
 if(c == '\\')
  {
   /* escape sequence */
    c = str[i];
            i++;

    if(c == 'a')
     {
      cc = '\a';
     }
    else if(c == 'b')
     {
      cc = '\b';
     }
    else if(c == 'f')
     {
      cc = '\f';
     }
    else if(c == 'n')
     {
      cc = '\n';
     }
    else if(c == 'r')
     {
      cc = '\r';
     }
    else if(c == 't')
     {
      cc = '\t';
     }
    else if(c == 'v')
     {
      cc = '\v';
     }
    else if(c == '\\')
     {
      cc = '\\';
     }
    else if(c == '\"')
     {
      cc = '\"';
     }
    else if(c == '\'')
     {
      cc = '\'';
     }
    else if(c == '\?')
     {
      cc = '\?';
     }
    else if(c == '$')
     {
      cc = '$';
     }


    else if(c == 'x')
     {
      /* two hex symbols expected */
      if(i==n)
       {
        /* empty */
        flag=4;
       }

      else
       {
        iv = 0;

        c = str[i];
                i++;

        if( isdigit(c) )
         {
          ic = c - '0';
         }
        else if( islower(c) )
         {
          ic = c - 'a';
          ic = ic + 10;
         }
        else if( isupper(c) )
         {
          ic = c - 'A';
          ic = ic + 10;
         }
        else
         {
          flag=2;
         }

        if(flag)
         {
          ;
         }
        else if(ic < hexbase)
         {
             iv = iv*hexbase + ic;
          if(iv > icmax)
           {
            flag=2;
           }
          else
           {
            if(i<n)
             {
              c = str[i];
                      i++;

              if( isdigit(c) )
               {
                ic = c - '0';
               }
              else if( islower(c) )
               {
                ic = c - 'a';
                ic = ic + 10;
               }
              else if( isupper(c) )
               {
                ic = c - 'A';
                ic = ic + 10;
               }
              else
               {
                flag=2;
               }

              if(flag)
               {
                ;
               }
              else if(ic < hexbase)
               {
                   iv = iv*hexbase + ic;
                if(iv > icmax)
                 {
                  flag=2;
                 }
                else
                 {
                  cc = (char)iv;
                 }
               }
              else
               {
                /* out of range */
                flag=2;
               }
             }
            else
             {
              /* end of input */
              flag=2;
             }
           }
         }
        else
         {
          /* out of range */
          flag=2;
         }

       }/*i<n*/

     }/*hex*/


    else if( isdigit(c) )
     {
      /* three octal symbols expected */
          iv   = 0;
           k   = 0;
          ibrk = 0;
        do
         {
             ic = c - '0';
          if(ic < octbase)
           {
               iv = iv*octbase + ic;
            if(iv > icmax)
             {
              flag=2;
              ibrk=1;
             }
           }
          else
           {
            /* out of range */
            flag=2;
            ibrk=1;
           }

          if(i<n)
           {
            c = str[i];
                    i++;

               k++;
            if(k > 2)
             {
              ibrk=1;
             }
            else if( !isdigit(c) )
             {
              flag=3;
              ibrk=1;
             }
           }
          else
           {
            /* end of input */
               k++;
            if(k>2)
             {
              ibrk=1;
             }
            else
             {
              flag=2;
              ibrk=1;
             }
           }

         } while(ibrk==0);

      if(flag==0) cc = (char)iv;

     }/*oct*/


    else
     {
      /* improper escape symbol */
      flag=3;
     }


    if(flag==0)
     {
      if(j<nval)
       {
        sval[j] = cc;
             j++;
       }
      else
       {
        /* end of output */
        break;
       }
     }

   }/* esc */


  else
   {
    /* ordinary symbol */

    if(j<nval)
     {
      /* add symbol */
      sval[j] = c;
           j++;

      c = str[i];
              i++;
     }
    else
     {
      /* end of output */
      break;
     }
   }

} /* i*/


if(j==0)
 {
  /* empty */
  flag=1;
 }
else
 {
  sval[j] = IDF_EOS;
 }

fin:
return flag;
}
