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


static IDF_TXT_ERR TxT_Err;



void idf_txt_err_nul()
{
TxT_Err.jerr = 0;

TxT_Err.kstr = 0;
TxT_Err.kurs = 0;
TxT_Err.kpos = 0L;

TxT_Err.cerr = ' ';
}



#if IDF_CPP_ == 1
void idf_txt_err_put(int jerr, int kstr, int kurs, long kpos,
                     unsigned int c)
#else
void idf_txt_err_put(jerr, kstr,kurs,kpos, c)
int  jerr;
int  kstr,kurs;
long kpos;
unsigned int c;
#endif
{
int n;
char *s;
register int i;

TxT_Err.jerr = jerr;

TxT_Err.kstr = kstr;
TxT_Err.kurs = kurs;
TxT_Err.kpos = kpos;

TxT_Err.cerr = c;
}



void idf_txt_err_prn()
{
#define IDF_ERR_MESSAGE_N 205

static char *ErrMes[IDF_ERR_MESSAGE_N] = {

/*  1*/ "file reading error in comment",
/*  2*/ "unexpected EoF in comment",
/*  3*/ "improper comment treating algorithm",
/*  4*/ "file reading error in formular",
/*  5*/ "unexpected EoF in formular",
/*  6*/ "empty formular",
/*  7*/ "improper formular",
/*  8*/ "improper formular treating algorithm",
/*  9*/ "file reading error in char",
/* 10*/ "unexpected EoF in char",
/* 11*/ "improper symbol in char field",
/* 12*/ "improper char sequence",
/* 13*/ "empty char",
/* 14*/ "improper char treating algorithm",
/* 15*/ "file reading error in string",
/* 16*/ "unexpected EoF in string",
/* 17*/ "unexpected EoS in string",
/* 18*/ "empty string",
/* 19*/ "suspiciously long string input",
/* 20*/ "improper string treating algorithm",
/* 21*/ "file reading error in text",
/* 22*/ "unexpected EoF in text",
/* 23*/ "empty text",
/* 24*/ "incorrect text unit",
/* 25*/ "improper trigraph in text",
/* 26*/ "suspiciously long text input",
/* 27*/ "improper text treating algorithm",
/* 28*/ "file reading error in data field",
/* 29*/ "unexpected EoF in data field",
/* 30*/ "improper symbol in data field",
/* 31*/ "unclosed { in data field",
/* 32*/ "empty data field",
/* 33*/ "extra } in data field",
/* 34*/ "error in skipping data field",
/* 35*/ "file reading error in data period field",
/* 36*/ "unexpected EoF in data period field",
/* 37*/ "improper symbol in data period field",
/* 38*/ "unclosed { in data period field",
/* 39*/ "empty data period field",
/* 40*/ "algorithm error in skipping data period field",
/* 41*/ "unexpected end of data period field",
/* 42*/ "error in skipping data period field",
/* 43*/ "file reading error in name field",
/* 44*/ "unexpected EoF in name field",
/* 45*/ "improper symbol in name field",
/* 46*/ "suspiciously long name input",
/* 47*/ "empty name field",
/* 48*/ "too long name field",
/* 49*/ "empty raw name buffer",
/* 50*/ "improper raw name buffer",
/* 51*/ "letter in name dimension field",
/* 52*/ "name starts with a digit",
/* 53*/ "too large number in name dimension field",
/* 54*/ "name starts with underscore or percent",
/* 55*/ "too many dimensions specified in file",
/* 56*/ "name starts with ()[]",
/* 57*/ "in name a digit before [ or (",
/* 58*/ "improper appearance of ] or ) in name",
/* 59*/ "extra [ or ( in name",
/* 60*/ "improper usage of [ or ) in name",
/* 61*/ "unclosed [ or ( detected in name",
/* 62*/ "empty name",
/* 63*/ "mixture [] and () in name",
/* 64*/ "improper symbol in raw name buffer",
/* 65*/ "empty name",
/* 66*/ "improper name",
/* 67*/ "no significant symbols in name",
/* 68*/ "white symbol inside a name",
/* 69*/ "name treating algorith error",
/* 70*/ "name is a keyword",
/* 71*/ "no name itself",
/* 72*/ "name symbols are separated by white space",
/* 73*/ "extra white space inside name",
/* 74*/ "empty sub-name",
/* 75*/ "algorithm error=9 for sub-names",
/* 76*/ "algorithm error=10 for sub-names",
/* 77*/ "improper sequence of keywords in name",
/* 78*/ "improperly placed keyword in name",
/* 79*/ "local names must not be of 'c,s,t' types",
/* 80*/ "improper position of keyword in name",
/* 81*/ "locals must not have dimensions",
/* 82*/ "name of data is followed by EoF",
/* 83*/ "empty named data near EoF",
/* 84*/ "file reading error in complex",
/* 85*/ "unexpected EoF in complex",
/* 86*/ "empty complex",
/* 87*/ "incorrect complex",
/* 88*/ "improper symbol in complex",
/* 89*/ "suspiciously long complex input",
/* 90*/ "improper complex treating algorithm",
/* 91*/ "error in reading symbol in data unit",
/* 92*/ "unexpected EoF in getting data unit",
/* 93*/ "extra comma separator in data unit",
/* 94*/ "improper usage of '*' in data unit",
/* 95*/ "no symbols in conversion of data unit repeatition number",
/* 96*/ "no digits in conversion of data unit repeatition number",
/* 97*/ "improper symbol in conversion of data unit repeatition number",
/* 98*/ "overflow in conversion of data unit repeatition number",
/* 99*/ "improper integer in conversion of data unit repeatition number",
/*100*/ "improper data type in conversion of data unit repeatition number",
/*101*/ "unable to convert into data unit repeatition number",
/*102*/ "negative data unit repeatition number",
/*103*/ "missing data unit repeatition number",
/*104*/ "unexpected EoF after data period {",
/*105*/ "unable to get position of data period {",
/*106*/ "too many data periods {",
/*107*/ "unexpected right brace in data unit",
/*108*/ "extra } in data unit",
/*109*/ "unable to continue the nearest period",
/*110*/ "improper symbol / in data unit",
/*111*/ "too many symbols in digital data unit",
/*112*/ "unexpected end of data field",
/*113*/ "unclosed data period { before data end",
/*114*/ "improper symbol in data unit field",
/*115*/ "unable to get a data unit",
/*116*/ "no symbols in character unit",
/*117*/ "unkwnown escape sequence in character unit",
/*118*/ "improper symbol in character unit",
/*119*/ "empty escape sequence in character unit",
/*120*/ "improper character unit",
/*121*/ "improper data type for conversion into character",
/*122*/ "no symbols in string unit",
/*123*/ "unkwnown escape sequence in string unit",
/*124*/ "improper symbol in string unit",
/*125*/ "empty escape sequence in string unit",
/*126*/ "improper string unit",
/*127*/ "improper data type for conversion into string",
/*128*/ "no symbols in text unit",
/*129*/ "improper data type for conversion into text",
/*130*/ "no symbols in integer/long unit",
/*131*/ "no digits in integer/long unit",
/*132*/ "improper symbol in integer/long unit",
/*133*/ "overflow in conversion integer/long unit",
/*134*/ "improper integer/long unit",
/*135*/ "improper data type for conversion into integer/long",
/*136*/ "unable to convert buffer into integer/long",
/*137*/ "no symbols in float/double unit",
/*138*/ "no digits in float/double unit",
/*139*/ "improper symbol in float/double unit",
/*140*/ "overflow in conversion float/double unit",
/*141*/ "improper style of float/double unit",
/*142*/ "improper data type for conversion into float/double", 
/*143*/ "unable to convert buffer into float/double",
/*144*/ "no symbols in complex unit",
/*145*/ "no digits in complex unit",
/*146*/ "improper symbol in complex unit",
/*147*/ "overflow in conversion complex unit",
/*148*/ "improper style of complex unit",
/*149*/ "improper data type for conversion into complex", 
/*150*/ "unable to convert buffer into complex",
/*151*/ "improper format for buffer conversion",
/*152*/ "no data buffer available for conversion",
/*153*/ "no data on formular result available for conversion",
/*154*/ "data unit got improper repeatition number",
/*155*/ "improper raw data unit storage type",
/*156*/ "inconsistency between input format and data type",
/*157*/ "overflow in data conversion according to input format",
/*158*/ "improper formular call",
/*159*/ "file reading error in formular",
/*160*/ "unexpected EoF formular case",
/*161*/ "too many elements accumulated in formular",
/*162*/ "too many constants accumulated in formular",
/*163*/ "parse error in formular",
/*164*/ "error in comment in formular",
/*165*/ "number in formular starts with improper symbol",
/*166*/ "i/o error in formular",
/*167*/ "unexpected EoF in formular",
/*168*/ "improper floating point in formular",
/*169*/ "too many symbols for a formular number",
/*170*/ "improper float number style in formular",
/*171*/ "improper symbol in formular number",
/*172*/ "unexpected ; in formular",
/*173*/ "expected number in formular, got no digits",
/*174*/ "name in formular starts with improper symbol",
/*175*/ "too many symbols for a formular name",
/*176*/ "improper symbol in formular name",
/*177*/ "expected name in formular, got no letters",
/*178*/ "improper formular buffer type for conversion",
/*179*/ "no symbols to convert formular buffer into double",
/*180*/ "no digits to convert formular buffer into double",
/*181*/ "improper symbol in conversion of formular buffer into double",
/*182*/ "overflow in conversion of formular buffer into double",
/*183*/ "unable to convert formular buffer into double",
/*184*/ "inconsistency in formular constant and local name list type",
/*185*/ "unknown name for constant or function in formular",
/*186*/ "local named data has not been calculated for formular constant",
/*187*/ "unable to search in local name list for formular constant",
/*188*/ "improper usage of function name in formular",
/*189*/ "too many accumulated 'brackets' in formular",
/*190*/ "unable to distinguish number in formular",
/*191*/ "more than 2 arguments for function in formular",
/*192*/ "improper binary token",
/*193*/ "math error in binary token execution",
/*194*/ "no number before operand in binary token",
/*195*/ "extra bracket in formular",
/*196*/ "parse error in accumulated expression",
/*197*/ "math error in function execution",
/*198*/ "argument is missing for function in formular",
/*199*/ "too many arguments for function in formular",
/*200*/ "improper logic in accumulated expression",
/*201*/ "improper accumulated expression",
/*202*/ "unclosed bracket in formular",
/*203*/ "improper symbol in formular body",
/*204*/ "improper formular",
/*   */ "unknown error"
 };

register int j;
int jerr;
char *s=NULL;

   jerr = TxT_Err.jerr;
if(jerr>0)
 {
     j = jerr;
  if(j > IDF_ERR_MESSAGE_N) j = IDF_ERR_MESSAGE_N;
   {
     j--;
    s = ErrMes[j];

    fprintf(stdout,"error=%4d at position=%ld string=%d kursor=%d\n",
            jerr,TxT_Err.kpos,TxT_Err.kstr,TxT_Err.kurs);
    fprintf(stdout,"message:%s  ,cerr=%c\n",
            s,TxT_Err.cerr);
   }
 }

}
