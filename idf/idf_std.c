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




#if IDF_CPP_ == 1
int idf_string_over(unsigned int c)
#else
int idf_string_over(c)
unsigned int c;
#endif
{
int k;

     if(c == '\n'    ) k=1;
else if(c == '\r'    ) k=1;
else if(c == IDF_EOS ) k=2;
else if(c == IDF_EOF ) k=2;
else                   k=0;

return k;
}



#if IDF_CPP_ == 1
int idf_name_seps(unsigned int c)
#else
int idf_name_seps(c)
unsigned int c;
#endif
{
int k;

       if(c == ':') k = 1;
  else if(c == '=') k = 2;
  else              k = 0;

return k;
}


#if IDF_CPP_ == 1
int idf_number_over(unsigned int c)
#else
int idf_number_over(c)
unsigned int c;
#endif
{
int k;

       if( isspace(c) ) k = 1;
  else if(c ==    ',' ) k = 2;
  else if(c ==    '*' ) k = 3;
  else if(c ==    ';' ) k = 4;
  else if(c ==    '}' ) k = 5;
  else if(c ==    '{' ) k = 6;
  else                  k = 0;

return k;
}


#if IDF_CPP_ == 1
int idf_data_seps(unsigned int c)
#else
int idf_data_seps(c)
unsigned int c;
#endif
{
int k;

       if(c == ';') k = 2;
  else if(c == ',') k = 1;
  else              k = 0;

return k;
}


#if IDF_CPP_ == 1
int idf_name_symbol(unsigned int c)
#else
int idf_name_symbol(c)
unsigned int c;
#endif
{
int k;

 if(c == '_')
  {
   k=1;
  }
 else if(c == '%')
  {
   k=2;
  }
 else if(c == '(')
  {
   k=3;
  }
 else if(c == ')')
  {
   k=3;
  }
 else if(c == '[')
  {
   k=4;
  }
 else if(c == ']')
  {
   k=4;
  }
 else if( isalpha(c) )
  {
   k=5;
  }
 else if( isdigit(c) )
  {
   k=6;
  }
 else
  {
   k=0;
  }

return k;
}


#if IDF_CPP_ == 1
int idf_data_symbol(unsigned int c)
#else
int idf_data_symbol(c)
unsigned int c;
#endif
{
int k;

 if(c == '+')
  {
   k=1;
  }
 else if(c == '-')
  {
   k=2;
  }
 else if(c == '.')
  {
   k=3;
  }
 else if( isalpha(c) )
  {
   k=4;
  }
 else if( isdigit(c) )
  {
   k=5;
  }
 else
  {
   k=0;
  }

return k;
}


#if IDF_CPP_ == 1
int idf_hex_symbol(unsigned int c)
#else
int idf_hex_symbol(c)
unsigned int c;
#endif
{
/* 
hex symbols: xX 0123456789 abcdef ABCDEF

returns:
   >0 hex symbol index
    0 if x or X
   <0 if non hex
*/
static int base = 16;
int k;

if( isdigit(c) )
 {
  k = c - '0';
  if(k < base) k=k+1; else k=0;
 }
else if((c == 'x')||(c == 'X'))
 {
  k=0;
 }
else if( islower(c) )
 {
  k = c - 'a';
  if(k < base) k=k+11; else k=0;
 }
else if( isupper(c) )
 {
  k = c - 'A';
  if(k < base) k=k+17; else k=0;
 }
else
 {
  k=-1;
 }

return k;
}


#if IDF_CPP_ == 1
int idf_octal_symbol(unsigned int c)
#else
int idf_octal_symbol(c)
unsigned int c;
#endif
{
/* 
octal symbols: 01234567

returns:
   >=0 octal symbol index
    <0 if non octal
*/
int k;

if( isdigit(c) )
 {
  k = c - '0';
  if(k > 8) k=0;
 }
else
 {
  k=-1;
 }

return k;
}


#if IDF_CPP_ == 1
int idf_exp_symbol(unsigned int c)
#else
int idf_exp_symbol(c)
unsigned int c;
#endif
{
/*
 floating point exponent symbols
*/
int k;

     if(c == 'e') k=1;
else if(c == 'E') k=1;
else if(c == 'd') k=1;
else if(c == 'D') k=1;
else if(c == 'g') k=1;
else if(c == 'G') k=1;
else              k=0;

return k;
}


#if IDF_CPP_ == 1
int idf_form_symbol(unsigned int c)
#else
int idf_form_symbol(c)
unsigned int c;
#endif
{
/*
 additional symbols allowed in formular
*/
int k;

     if(c == '*')  k=1;
else if(c == '/')  k=1;
else if(c == '(')  k=1;
else if(c == ')')  k=1;
else if(c == '^')  k=1;
else if(c == ',')  k=1;
else if(c == '+')  k=1;
else if(c == '-')  k=1;
else               k=0;

return k;
}


#if IDF_CPP_ == 1
int idf_file_symbol(unsigned int c)
#else
int idf_file_symbol(c)
unsigned int c;
#endif
{
/*
 additional symbols allowed for file name
*/
int k;

     if(c == '.')  k=1;
else if(c == ']')  k=1;
else if(c == '[')  k=1;
else if(c == ':')  k=1;
else if(c == '/')  k=1;
else if(c == '\\') k=1;
else if(c == '_')  k=1;
else if(c == '#')  k=1;
else if(c == '$')  k=1;
else if(c == '%')  k=1;
else if(c == ',')  k=1;
else               k=0;

return k;
}
