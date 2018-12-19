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


#ifndef IDF_I
#define IDF_I

/* if =1 then non-ascii symbols are treated as letters,
  otherwise appearence of non-ascii symbol elsewhere
  than in text type data unit will result in a fatal mistake */
#define IDF_FOREIGN_LANGUAGE 0

/* In C language, the text type data is declared as CHAR 
   (not as UNSIGNED CHAR).
   In FORTRAN language, character may correspond to 
   either CHAR or UNSIGNED CHAR (depending on a compiler).
   Define the following macro as 1 if FORTRAN character
   corresponds to UNSIGNED CHAR, and define macro as zero otherwise*/
#define IDF_TEXT_CHAR_FORTRAN 1

/* if =1 then print warning messages */
#define IDF_WARNING 0

/* if =1 then print error messages
  otherwise use fuctions to get messages */
#define IDF_MISTAKE 1

/* if =1 include \v in calculating strings number in file */
#define IDF_TXT_VERTICES 0

/* if =1 include \t in calculating chars number in file string */
#define IDF_TXT_TABULATORS 0

/* if =1 then EOS symbol is a blank space separator
 the arithmetic type data unit,
 otherwise EOS is simply ignored */
#define IDF_EOS_BREAK_NUMBER 1

/* if =1 then EOS symbol is a blank space separator 
  for a name in formular, otherwise EOS is simply ignored */
#define IDF_EOS_BREAK_FNAME 1

/*
If '0' then all white spaces neglected in the data 
If '1' then white space is an error if appeares
between letters or digits.
In all other cases, name is allowed to be composed
from the keyword specified data type (ie. int,float,double etc)
and the name itself which should not contain white spaces.
*/
#define IDF_NAME_SPACE 2

/*
If 0 then all white spaces neglected in the name 
If != 0 then white space is an error if appeares
between letters or digits. 
*/
#define IDF_DATA_SPACE 1


/* if != 0 the trace the name, data, and format units */
#define IDF_NAME_TRACE 0
#define IDF_DATA_TRACE 0
#define IDF_FRMT_TRACE 0


/* the following macros are used to control the order
   in which idf functions are called */ 
#define IDF_ORDER_START 0
#define IDF_ORDER_INIT 1
#define IDF_ORDER_OPEN 2


/* maximal length of input-file-name string */
#define IDF_FSTR_MAX 125
#define IDF_FSTR_MAX_1 126

/* if !=0 store file name */
#define IDF_FSTR 1

typedef struct {                
                /* stream */
                FILE *fd;

                /* if !=0 then file is open */
                int   iopen;

                /* if !=0 then error */
                int   ferr;

                /* if !=0 then at EoF */
                int   fend;

#if IDF_FSTR == 1
                /* name of a file */
                char  name[IDF_FSTR_MAX_1];
#endif

               } IDF_FILE;


/* format types for data import
   1. 'c' char
   2. 's' char string
   3. 't' array of characters
   4. 'i' int
   5. 'l' long int
   6. 'f' float
   7. 'd' double
   8. 'z' complex float
   9. 'w' complex double

  10. 'v' void   
*/
#define IDF_FRMT_NTYP 9
#define IDF_FRMT_NTYP_1 10
#define IDF_FRMT_LIST "cstilfdzw"

#define IDF_FRMT_C 1
#define IDF_FRMT_S 2
#define IDF_FRMT_T 3
#define IDF_FRMT_I 4
#define IDF_FRMT_L 5
#define IDF_FRMT_F 6
#define IDF_FRMT_D 7
#define IDF_FRMT_Z 8
#define IDF_FRMT_W 9
#define IDF_FRMT_VOID 10
#define IDF_FRMT_LOCAL 11
#define IDF_FRMT_UNKNOWN 0

#define IDF_FRMT_CHAR 3
#define IDF_FRMT_DIGIT 7
#define IDF_FRMT_CMPLX 9


/* maximal number of units in a format string*/
#define IDF_FRMT_UNIT_NMAX 50

/* maximal number of elements in a format string*/
#define IDF_FRMT_ELEM_NMAX 50
#define IDF_FRMT_ELEM_NMAX_1 51

/* maximal number of dimensions
which can be assigned to format unit */
#define IDF_FRMT_DIM_MAX 3

/* maximal number for format unit repeatition*/
#define IDF_FRMT_MAX_IRPT 65535

/* maximal number in format unit dimensions*/
#define IDF_FRMT_MAX_IDIM 1000000

typedef struct {
                /* repeatition number */
                int Nrep;

                /* type= 1-9: cstilfdzw */
                int type;

                /* number of bytes in format element:
                it is sizeof (c,i,l,f,d,z,w) or
                maximal mumber of bytes in text or string*/
                int jsize;

                /* if pointer=0 then target is an object.
                If >0, the target is a pointer and 
                the value of pointer is the sizeof(type) */
                int pointer;

                /* if >0 then allocate memory to pointer */
                int jalloc;

                /* memory (in bytes) associated with format unit */
                int jmem;

                /* number of dimensions*/
                int Ndim;

                /* number of elements in a dimension*/
                int Mdim[IDF_FRMT_DIM_MAX];

               } IDF_FRMT_UNIT;


typedef struct {
                /* current number of format units*/
                int NfrmtU;

                /* current number of real format elements*/
                int NfrmtE;

                /* repeatition numbers for each element */
                int   Mrpt[IDF_FRMT_ELEM_NMAX];
                int   Krpt[IDF_FRMT_ELEM_NMAX];

                /* type of format element:%{}# */ 
                char  Rfrmt[IDF_FRMT_ELEM_NMAX_1];

                /* data on format unit */
                IDF_FRMT_UNIT Ufrmt[IDF_FRMT_UNIT_NMAX];

                /* maximal size of word in bytes */
                int swrd;

                /* error flag and position in format string*/
                int jerr;
                int kpos;

               } IDF_FRMT;

/*
 alignment rules:
     0 => byte boundary;
   1-6 => the boundary which is smaller of the selfsize 
          or the specified boundary:
               word     (1),
               longword (2),
               quadword (3),
               octaword (4),
               hexword  (5),
           and page     (6);
     7 => natural boundary appropriate to the type;
     8 => word packing size, pragma pack(1);
     9 => longword packing size, pragma pack(2);
    10 => quadword packing size, pragma pack(3);
    11 => octaword packing size, pragma pack(4);
    12 => hexword packing size,  pragma pack(5);
    13 => page packing size,     pragma pack(9);
*/
#define IDF_ALIGN_NATURAL 7
#define IDF_ALIGN_NMAX 14

#define IDF_ALIGN_BYTE 1
#define IDF_ALIGN_WORD 2
#define IDF_ALIGN_LONGWORD 4
#define IDF_ALIGN_QUADWORD 8
#define IDF_ALIGN_OCTAWORD 16
#define IDF_ALIGN_HEXWORD 32
#define IDF_ALIGN_PAGE 512

#define IDF_AMODE_TEXT 1
#define IDF_AMODE_BINARI 2

typedef struct {
                /* format element index */
                int ielem;

                /* format unit index */
                int iunit;

                /* repeatition counters */
                int irep;
                int jrep;

                /* offset and unfilled gap*/
                int offset;
                int gap;

                /* pointer flag */
                int pointer;

                /* memory allocation flag*/
                int jalloc;

                /* memory size for a format unit*/
                int jmem;

                /* pointer onto Original target */
                char *Otarget;

                /* pointer onto Allocated memory */
                char *Atarget;

                /* pointer onto current target */
                char *Target;

                /* format type */
                int jtype;

                /* format and its base word sizes */
                int jsize;
                int jwrd;

                /* max word size in format string*/
                int swrd;

                /* current alignment mode */
                int amode;

                /* alignment rule */
                int align;

                /* if >0 rewind if format string is over*/
                int blk;

                /* error flag */
                int jerr;

                /* number of dimensions*/
                int Ndimen;

                /* number of elements in a dimension*/
                int Mdimen[IDF_FRMT_DIM_MAX];
               
               } IDF_FRMT_TREAT;


/* maximal length of data-name buffer string */
#define IDF_NAME_LENGTH  255
#define IDF_NAME_LENGTH_1 256

/* maximal number of dimensions */
#define IDF_NAME_MAX_DIM 3

/* maximal number of digits denoting size in one of dimensions*/
#define IDF_NAME_MAX_ID 6

/* maximal number of symbols in 'real' name */
#define IDF_NAME_SIZE 125
#define IDF_NAME_SIZE_1 126

/*
keywords assigned to name:
 1. char character
 2. string 
 3. text
 4. int integer
 5. long
 6. float
 7. double
 8. complex
 9. complex
10. void
11. local
*/

/* max number of keywords assigned to a name */
#define IDF_KEYW_MAX_IN_NAME 5

typedef struct {
                /* buffer for name symbols */
		char name_buf[IDF_NAME_LENGTH_1];

                /* number of symbols put into buffer*/
		int length;

                /* name separator symbol */
		unsigned int eflag;

                /* starting position of name */
		int lb;
		int kb;
                long posb;

                /* ending position of name */
		int le;
		int ke;
                long pose;

                /* additional buffer for composed name analysis */
		char name_raw[IDF_NAME_LENGTH_1];
		int  nlength;

                /* resulting name (EoS included)*/
		char name[IDF_NAME_LENGTH_1];
                /* number of characters in name (without EoS)*/
		int  Nlength;

                /* number of dimensions assigned to name [][][]...*/
                int  ndim;
                /* numbers corresponding to [][][]...*/
                int  Mshft[IDF_NAME_MAX_DIM];

                /* data type keyword assigned to name */
                int  jtype;

                /* flag assigned to name. If flag=1 then 
                   data unit is local (ie for internal usage)
                   and it is external if flag=0 */
                int jlocal;

                /* ending position of data field
                   which is associated with that name */
		int ld;
		int kd;
                long posd;

                /* 'real' starting position of data 
                    from beginning of the file */
                long dspos;

                /* 'real' ending position of data 
                    from beginning of the file */
                long depos;

	       } IDF_NAME;




typedef struct {
                /* name of non-local data */
		char *name;
		int  length;

                /* starting position and length of data */
                int kstr;
                int kurs;
                long dpos;
                int  nd;

                /* data type keyword given in file */
                int  jtype;

                /* dimensions assigned to data */
                int  ndim;
                int  *Mshft;

	       } IDF_NAME_;


typedef struct {
                /* name of local data */
		char *name;
		int  length;

                /* starting position and size of data field */
                int kstr;
                int kurs;
                long dpos;
                int  nd;

                /* data type keyword given in file */
                int  jtype;

                /* If 0 then no value exists yet. If != then
                 value has been calculated and jest is its type
                 This is because locals may appear not consequently*/
                 int jest;

                /* data value is nominally strored as complex double
                  and will be converted if nessassary*/
                double  val[2];

	       } IDF_NAME_L;



/* maximal number of non-local named data in file */
#define IDF_NAME_LIST_N 200

/* maximal number of local named data in file */
#define IDF_NAME_LISTL_N 50

typedef struct {
                /* number of non-local (external) named data stored */
                int       Num;

                /* list of non-local names */
                IDF_NAME_ *List[IDF_NAME_LIST_N];

                /* number of local (internal) named data stored */
                int       NumL;

                /* list of local names */
                IDF_NAME_L *ListL[IDF_NAME_LISTL_N];

                /* error flag, level, and record index */
                int kerr;
                int lvl;
                int irec;

                } IDF_NAME_LIST;



/* maximal number of data periods */
#define IDF_DTP_MAX 50

/* maximal number repeatitions for a period */
#define IDF_DTP_MAX_NP 10000

typedef struct {
                /* starting position of named data*/
                int  kstr0;
                int  kurs0;
                long kpos0;

                /* type of data assigned to name in file
                  c=1 s=2 t=3 i=4 l=5 f=6 d=7 z=8 w=9 v=10   */
                int jtype;

                /* data unit number */
                int iunt;

                /* control symbol */
                char cntr;

                /* stored symbol */
                char ctor;

                /* unget flag */
                int flag;

                /* position of current data unit*/
                int  kstr;
                int  kurs;
                long kpos;

                /*  repeatition number for a unit */
                int nrpt;

                /* if 0 no data unit yet
                   if 1 data unit is stored in character buffer
                   if 2 data value is stored as formular 
                   execution result
                */
                int Jbuf;

                /* number of opened data periods */
                int  DTPn;

                /* position of { in file */
                int   DTPstr[IDF_DTP_MAX];
                int   DTPkur[IDF_DTP_MAX];
                long  DTPpos[IDF_DTP_MAX];

                /* length of period */
                int  DTPrpt[IDF_DTP_MAX];

                } IDF_DATA;



/* maximal length for data unit buffer */
#define IDF_DBUF_LENGTH  327
#define IDF_DBUF_LENGTH_1  328

/* if !=0 then use standard C functions
   to convert character string to a number */
#define IDF_DBUF_FUN 0

typedef struct {
                /* length of current record */
		int length;

                /* starting position in file */
                int  kstr;
                int  kurs;
                long pos;

                /* type of data: 0-nothing, C,S,T,DIGIT,CMPLX */
                int jtype;

		char Data_buf[IDF_DBUF_LENGTH_1];

	       } IDF_DBUF;


typedef struct {
                /* type of data: 0-nothing, C,S,T,I,L,F,D,Z,W */
                int jtype;

                char Cval;

		char Sval[IDF_DBUF_LENGTH_1];
                int  len;

                int  Ival;

                long Lval;

                float Fval[2];

                double Dval[2];                

	       } IDF_VAL;


/* maximal length for formular unit buffer */
#define IDF_FBUF_LENGTH  125
#define IDF_FBUF_LENGTH_1  126

#define IDF_FBUF_TYPE_NUMBER 1
#define IDF_FBUF_TYPE_NAME 2

typedef struct {
                /* length of current record */
		int length;

                /* starting position in file */
                int  kstr;
                int  kurs;
                long pos;

                /* type of data: 0-nothing, 1-number 2-name */
                int jtype;

                /* result of conversion: 0 - nothing, 
                  1 - numerical constant value, 2 - function */
                int jconv;

                /* it is function type if jconv=2,
                 and it is constant type if jconv=1 */
                int jf;

                /* numerical constant value*/
                double val[2];

                /* character buffer*/
		char buf[IDF_FBUF_LENGTH_1];

	       } IDF_FBUF;



/*maximal number of accumulated constants in formular*/
#define IDF_FORM_CNS_MAX 50

/*maximal number of accumulated elements in formular*/
#define IDF_FORM_ELM_MAX 99
#define IDF_FORM_ELM_MAX_1 100

/*maximal number of accumulated brackets in formular*/
#define IDF_FORM_BRACE_MAX 32

typedef struct {
                /* starting position in file */
                int  kstr;
                int  kurs;
                long kpos;

                /* brace sequence                            */
                int  iBRC;                       /* counter  */
                char  Brace[IDF_FORM_BRACE_MAX]; /* type: ([{*/

                /* list of constants accumulated in formular  */
                int   iCNS;                      /* counter   */
                int    CNStyp[IDF_FORM_CNS_MAX]; /* type      */
                double CNSvr [IDF_FORM_CNS_MAX]; /* real part */
                double CNSvi [IDF_FORM_CNS_MAX]; /* image part*/

                /* list of elements accumulated in formular        */
                int   iELM;                        /* counter      */
                char   ELMtyp[IDF_FORM_ELM_MAX_1]; /* type:#(){[],%*/
                int    ELMdat[IDF_FORM_ELM_MAX_1]; /* parameter    */

                /* formular data: 0 - no result, DIGIT, CMPLX */
                int jtype;

                /* formular execution result is stored 
                   as complex double */                
		double Val[2];

	       } IDF_FORM;


/* macros for math expression treating */
#define IDF_HIER_NOTHING -1
#define IDF_HIER_NUMBER   0
#define IDF_HIER_PLUS     1
#define IDF_HIER_MULT     2
#define IDF_HIER_POW      3


/* maximal base-10-exponent value */
#define IDF_DBL_10_EXP_MAX DBL_MAX_10_EXP

/* max and min double value */
#define IDF_DOUBLE_MAX DBL_MAX
#define IDF_DOUBLE_MIN DBL_MIN

/* max and min float value */
#define IDF_FLOAT_MAX FLT_MAX
#define IDF_FLOAT_MIN FLT_MIN

/* max long value */
#define IDF_LONG_MAX INT_MAX

/* max int value */
#define IDF_INT_MAX INT_MAX

/* max char value */
#define IDF_CHAR_MAX CHAR_MAX


#if IDF_IEEE == 1

#define IDF_EXPOWTAB_MAX 256
#define IDF_PTAB {1.e256,1.e128,1.e64,1.e32,1.e16,1.e8,1.e4,1.e2,1.e1}
#define IDF_NTAB {1.e-256,1.e-128,1.e-64,1.e-32,1.e-16,1.e-8,1.e-4,0.01,0.1,1.0}

#else

#define IDF_EXPOWTAB_MAX 32
#define IDF_PTAB {1.e32,1.e16,1.e8,1.e4,1.e2,1.e1}
#define IDF_NTAB {1.e-32,1.e-16,1.e-8,1.e-4,0.01,0.1,1.0}

#endif

/* initialization option for float structure*/
#define IDF_FLOAT_DEFAULT 1

/* natural logarithm of maxiaml double*/
#define IDF_LNMAXDOUBLE LN_MAXDOUBLE

/* log(10.0) */
#define IDF_LN10     M_LN10

/* pi*0.5 */
#define IDF_PIOVER2  M_PI_2
#define IDF_PI       M_PI

/* math errors */
#define IDF_MATH_OVERFLOW  1
#define IDF_MATH_EXCEPTION 2
#define IDF_MATH_DOMAIN    3
#define IDF_MATH_RANGE     4
#define IDF_MATH_OPERAND   5
#define IDF_MATH_DFUN      6
#define IDF_MATH_CFUN      7
#define IDF_MATH_OPERARG   8
#define IDF_MATH_FUN       9
#define IDF_MATH_FARG     10
#define IDF_MATH_NARG     11

/* text symbols */
#define IDF_EOS '\0'
#define IDF_EOF EOF


/* In the case when less symbols are imported
   from data file, fulfil the remaider positions with
   the specified symbol (IDF_TEXT_FULFIL_SYMBOL)
   if macro IDF_TEXT_FULFIL=1 
*/
#define IDF_TEXT_FULFIL 1
#define IDF_TEXT_FULFIL_SYMBOL IDF_EOS


/* if != 0 print additional error messages */
#define IDF_ERRMESSAGE_DETAILED 0

typedef struct {
                /* error number */
                int jerr;

                /* ending position */
                int kstr,kurs;
                long kpos;

                /* error symbol */
                unsigned int cerr;

               } IDF_TXT_ERR;

#endif /* IDF_I */
