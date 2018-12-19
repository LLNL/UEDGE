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



#define IDF_STRTOC_SKIP_W 0




/*
     This is a map of the ASCII character  set.
     It contains octal values for each character.

     |000 NUL|001 SOH|002 STX|003 ETX|004 EOT|005 ENQ|006 ACK|007 BEL|
     |010 BS |011 HT |012 NL |013 VT |014 NP |015 CR |016 SO |017 SI |
     |020 DLE|021 DC1|022 DC2|023 DC3|024 DC4|025 NAK|026 SYN|027 ETB|
     |030 CAN|031 EM |032 SUB|033 ESC|034 FS |035 GS |036 RS |037 US |
     |040 SP |041  ! |042  " |043  # |044  $ |045  % |046  & |047  ' |
     |050  ( |051  ) |052  * |053  + |054  , |055  - |056  . |057  / |
     |060  0 |061  1 |062  2 |063  3 |064  4 |065  5 |066  6 |067  7 |
     |070  8 |071  9 |072  : |073  ; |074  < |075  = |076  > |077  ? |
     |100  @ |101  A |102  B |103  C |104  D |105  E |106  F |107  G |
     |110  H |111  I |112  J |113  K |114  L |115  M |116  N |117  O |
     |120  P |121  Q |122  R |123  S |124  T |125  U |126  V |127  W |
     |130  X |131  Y |132  Z |133  [ |134  \ |135  ] |136  ^ |137  _ |
     |140  ` |141  a |142  b |143  c |144  d |145  e |146  f |147  g |
     |150  h |151  i |152  j |153  k |154  l |155  m |156  n |157  o |
     |160  p |161  q |162  r |163  s |164  t |165  u |166  v |167  w |
     |170  x |171  y |172  z |173  { |174  | |175  } |176  ~ |177 DEL|


     This is a map of the ASCII character  set.
     It contains hexadecimal values for each character.

     | 00 NUL| 01 SOH| 02 STX| 03 ETX| 04 EOT| 05 ENQ| 06 ACK| 07 BEL|
     | 08 BS | 09 HT | 0A NL | 0B VT | 0C NP | 0D CR | 0E SO | 0F SI |
     | 10 DLE| 11 DC1| 12 DC2| 13 DC3| 14 DC4| 15 NAK| 16 SYN| 17 ETB|
     | 18 CAN| 19 EM | 1A SUB| 1B ESC| 1C FS | 1D GS | 1E RS | 1F US |
     | 20 SP | 21  ! | 22  " | 23  # | 24  $ | 25  % | 26  & | 27  ' |
     | 28  ( | 29  ) | 2A  * | 2B  + | 2C  , | 2D  - | 2E  . | 2F  / |
     | 30  0 | 31  1 | 32  2 | 33  3 | 34  4 | 35  5 | 36  6 | 37  7 |
     | 38  8 | 39  9 | 3A  : | 3B  ; | 3C  < | 3D  = | 3E  > | 3F  ? |
     | 40  @ | 41  A | 42  B | 43  C | 44  D | 45  E | 46  F | 47  G |
     | 48  H | 49  I | 4A  J | 4B  K | 4C  L | 4D  M | 4E  N | 4F  O |
     | 50  P | 51  Q | 52  R | 53  S | 54  T | 55  U | 56  V | 57  W |
     | 58  X | 59  Y | 5A  Z | 5B  [ | 5C  \ | 5D  ] | 5E  ^ | 5F  _ |
     | 60  ` | 61  a | 62  b | 63  c | 64  d | 65  e | 66  f | 67  g |
     | 68  h | 69  i | 6A  j | 6B  k | 6C  l | 6D  m | 6E  n | 6F  o |
     | 70  p | 71  q | 72  r | 73  s | 74  t | 75  u | 76  v | 77  w |
     | 78  x | 79  y | 7A  z | 7B  { | 7C  | | 7D  } | 7E  ~ | 7F DEL|


     This is a map of the ASCII character  set.
     It contains decimal values for each character.

     |  0 NUL|  1 SOH|  2 STX|  3 ETX|  4 EOT|  5 ENQ|  6 ACK|  7 BEL|
     |  8 BS |  9 HT | 10 NL | 11 VT | 12 NP | 13 CR | 14 SO | 15 SI |
     | 16 DLE| 17 DC1| 18 DC2| 19 DC3| 20 DC4| 21 NAK| 22 SYN| 23 ETB|
     | 24 CAN| 25 EM | 26 SUB| 27 ESC| 28 FS | 29 GS | 30 RS | 31 US |
     | 32 SP | 33  ! | 34  " | 35  # | 36  $ | 37  % | 38  & | 39  ' |
     | 40  ( | 41  ) | 42  * | 43  + | 44  , | 45  - | 46  . | 47  / |
     | 48  0 | 49  1 | 50  2 | 51  3 | 52  4 | 53  5 | 54  6 | 55  7 |
     | 56  8 | 57  9 | 58  : | 59  ; | 60  < | 61  = | 62  > | 63  ? |
     | 64  @ | 65  A | 66  B | 67  C | 68  D | 69  E | 70  F | 71  G |
     | 72  H | 73  I | 74  J | 75  K | 76  L | 77  M | 78  N | 79  O |
     | 80  P | 81  Q | 82  R | 83  S | 84  T | 85  U | 86  V | 87  W |
     | 88  X | 89  Y | 90  Z | 91  [ | 92  \ | 93  ] | 94  ^ | 95  _ |
     | 96  ` | 97  a | 98  b | 99  c |100  d |101  e |102  f |103  g |
     |104  h |105  i |106  j |107  k |108  l |109  m |110  n |111  o |
     |112  p |113  q |114  r |115  s |116  t |117  u |118  v |119  w |
     |120  x |121  y |122  z |123  { |124  | |125  } |126  ~ |127 DEL|

*/




#if IDF_CPP_ == 1
int idf_str_to_c(char* str, int nn, char* cval)
#else
int idf_str_to_c(str,nn, cval)
char *str;
int   nn;
char *cval;
#endif
{
/*
 -----------------------------------
   converts string to character
 -----------------------------------
returns:
  0  OK
  1  no data symbols
  2  unknown escape sequence
  3  improper symbol
  4  empty escape sequence
  5  improper character
*/
static int icmax = IDF_CHAR_MAX;
static int hexbase = 16;
static int octbase = 8;

int flag,n,ic,iv;
char c;
register int i;

 flag = 0;
*cval = ' ';

if(nn < 1)
 {
  /* no data symbols*/
  flag=1;
  goto fin;
 }


#if IDF_STRTOC_SKIP_W == 1

/* skip white spaces */

  i=nn;
while(i>0)
 {
  i--;

     c = str[i];
  if(c != ' ') break;
 }
n = i+1;

if(n == 1)
 {
  /* white space */
  goto fin;
 }

  i=0;
while(i<n)
 {
     c = str[i];
  if(c != ' ') break;
  i++;
 }

if(i==n)
 {
  /* white space */
  goto fin;
 }

#else

n = nn;

c = str[0];
i = 1;

#endif


if(c == '\\')
 {
  if(i==n)
   {
    /*empty escape sequence */
    flag=4;
   }
  else
   {
    /* escape sequence */
       c = str[i];
               i++;

    if(c == 'a')
     {
      *cval = '\a';
     }
    else if(c == 'b')
     {
      *cval = '\b';
     }
    else if(c == 'f')
     {
      *cval = '\f';
     }
    else if(c == 'n')
     {
      *cval = '\n';
     }
    else if(c == 'r')
     {
      *cval = '\r';
     }
    else if(c == 't')
     {
      *cval = '\t';
     }
    else if(c == 'v')
     {
      *cval = '\v';
     }
    else if(c == '\\')
     {
      *cval = '\\';
     }
    else if(c == '\"')
     {
      *cval = '\"';
     }
    else if(c == '\'')
     {
      *cval = '\'';
     }
    else if(c == '\?')
     {
      *cval = '\?';
     }
    else if(c == '$')
     {
      *cval = '$';
     }


    else if(c == 'x')
     {
      /* hex symbols expected */
      if(i<n)
      {
       iv = 0;
      while(i<n)
       {
        c = str[i];

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
          break;
         }

        if(ic < hexbase)
         {
             iv = iv*hexbase + ic;
          if(iv > icmax)
           {
            flag=2;
            break;
           }
         }
        else
         {
          /* out of range */
          flag=2;
          break;
         }
 
        i++;
       }

      if(flag==0)
       {
        *cval = (char)iv;
       }
      }/*i<n*/

      else
      {
       /* empty es */
       flag=4;
      }

     }/* hex */


    else if( isdigit(c) )
     {
      /* octal symbols expected */
      iv = 0;

         ic = c - '0';
      if(ic >= octbase)
       {
        /* out of range */
        flag=2;
       }
      else
       {
        iv = ic;
        while(i<n)
         {
          c = str[i];

          if( !isdigit(c) )
           {
            /* improper symbol */
            flag=3;
            break;
           }

          else
           {
               ic = c - '0';
            if(ic < octbase)
             {
                 iv = iv*octbase + ic;
              if(iv > icmax)
               {
                flag=2;
                break;
               }
             }
            else
             {
              /* out of range */
              flag=2;
              break;
             }
           }

          i++;
         }

        if(flag==0)
         {
          *cval = (char)iv;
         }
       }

     }/* oct */


    else
     {
      if(c == ' ') ic=0; 
      else         ic=1;

      while(i<n)
       {
        if(str[i] != ' ')
         {
          if(ic==0)
           {
            c = str[i];
           }
          ic++;
         }

        i++;
       }

      if(ic == 1)
       {
        *cval = c;
       }
      else
       {
        /* improper symbol */
        flag=3;
       }
     }
   }
 }


else
 {
  if(i==n)
   {
    /* single symbol*/
    *cval = c;
   }
  else
   {
    /* improper char */
    flag = 5;
   }
 }

fin:
return flag;
}
