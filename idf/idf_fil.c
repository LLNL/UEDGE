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


static char *Cname=NULL;

IDF_FILE  *FtxF=NULL;



IDF_FILE* idf_file_address()
{
 return ( FtxF );
}


char* idf_file_Cname()
{
 return Cname;
}


int idf_file_ini()
{
int ierr=0;
int k;

if(FtxF != NULL)
 {
  ierr=1;
  goto err;
 }

      k = sizeof(IDF_FILE);
   FtxF = (IDF_FILE*) malloc(k);
if(FtxF == NULL)
 {
  ierr=2;
 }

else
 {
  FtxF->iopen = 0;
  FtxF->ferr  = 0;
  FtxF->fend  = 0;
  FtxF->fd    = NULL;

#if IDF_FSTR == 1
  FtxF->name[0] = IDF_EOS;
#endif
 }

err:
return ierr;
}


void idf_file_end()
{
if(FtxF != NULL)
 {
  free(FtxF);
  FtxF = NULL;
 }

if(Cname != NULL)
 {
  free(Cname);
  Cname = NULL;
 }
}


void idf_file_prn()
{
if(FtxF != NULL)
 {
#if IDF_FSTR == 1
  fprintf(stdout,"File=%s\n",FtxF->name);
#endif

  fprintf(stdout,"File fopen=%d ferr=%d fend=%d\n",
          FtxF->iopen,FtxF->ferr,FtxF->fend);
 }
}


#if IDF_CPP_ == 1
int idf_file_name(char* Name)
#else
int idf_file_name(Name)
char *Name;
#endif
{
/*
-----------------------------------
  Check the C/C++ file name string
-----------------------------------
*/
int ierr=0;
unsigned int c;
register int i;
int l;

   l = strlen(Name);
if(l<1)
 {
  ierr=1;
 }

else if(l>IDF_FSTR_MAX)
 {
  ierr=2;
 }

else
 {
        i=0;
  while(i<l)
   {
    c = Name[i];

    if( !isascii(c) )
     {
#if IDF_FOREIGN_LANGUAGE == 1
      ;
#else
      ierr=3;
      break;
#endif
     }

    else if( !isalnum(c) )
     {
      if( !idf_file_symbol(c) )
       {
        /* improper symbol in file name*/
        ierr=3;
        break;
       }
     }

    i++;
   }
 }

if(ierr)  FtxF->ferr=ierr;

return ierr;
}



void idf_file_Cname_free()
{
if(Cname != NULL)
 {
  free(Cname);
  Cname = NULL;
 }
}


#if IDF_CPP_ == 1
int idf_file_name_(char* Name, int l)
#else
int idf_file_name_(Name,l)
char *Name;
int   l;
#endif
{
/*
----------------------------------------
 convert the non-C/C++ file name string
----------------------------------------
*/
int ierr=0;
unsigned int c;
register int i;
int n;

if(l<1)
 {
  ierr=1;
 }

else if(l>IDF_FSTR_MAX)
 {
  ierr=2;
 }

else
 {
        i=0;
  while(i<l)
   {
    c = Name[i];

    if( !isascii(c) )
     {
#if IDF_FOREIGN_LANGUAGE == 1
      ;
#else
      ierr=3;
      break;
#endif
     }

    else if( !isalnum(c) )
     {
      if( !idf_file_symbol(c) )
       {
        /* improper symbol in file name*/
        ierr=3;
        break;
       }
     }
    i++;
   }
 }

if(!ierr)
 {
  n = l+1;

  if(Cname != NULL)
   {
    /* unable to update input file name */
    ierr=4;
   }
  else
   {
       Cname = (char*) malloc(n);
    if(Cname == NULL)
     {
      /* no memory*/
      ierr=5;
     }
    else
     {
            i=0;
      while(i<l)
       {
        Cname[i] = Name[i];
              i++;
       }
      Cname[i] = IDF_EOS;
     }
   }
 }

return ierr;
}



#if IDF_CPP_ == 1
int idf_file_open(char* name)
#else
int idf_file_open(name)
char *name;
#endif
{
/*
-----------------------------
   open target idf_file
-----------------------------
*/
#if IDF_FSTR == 1
register int i;
char *s;
char c;
#endif

int ierr=0;
FILE *fil;

 fil = NULL;

  if(FtxF->iopen)
   {
    FtxF->ferr = 1;
    ierr = 1;
   }

  else
   {
     /*-----------------------*/
       fil = fopen(name,"r");
     /*-----------------------*/

    if(fil == NULL)
     {
      FtxF->ferr = 2;
      FtxF->fd   = NULL;

#if IDF_FSTR == 1
      FtxF->name[0] = IDF_EOS;
#endif
      ierr = 2;
     }

    else
     {
      /* rewind(fil); */

      FtxF->fd    = fil;
      FtxF->iopen = 1;
      FtxF->ferr  = 0;
      FtxF->fend  = 0;

#if IDF_FSTR == 1
        s = FtxF->name;
        i=0;
        while(i<IDF_FSTR_MAX)
         {
             c = name[i];
          if(c == IDF_EOS) break;
          s[i] = c;
            i++;
         }
        s[i] = IDF_EOS;
#endif
     }
   }

 return ierr;
}



void idf_file_close()
{
/*
-----------------------------
   close target idf_file
-----------------------------
*/
FILE *fil;
int ferr;

    fil = FtxF->fd;

   ferr = FtxF->ferr;
if(ferr)
 {
  fclose(fil);
 }

else
 {
     if(FtxF->iopen)
      {
       fflush(fil);
       fclose(fil);
      }

     else
      {
       ferr = 3;
      }
 }

	FtxF->fd    = NULL;
        FtxF->ferr  = ferr;
	FtxF->iopen = 0;
        FtxF->fend  = 0;

#if IDF_FSTR == 1
        FtxF->name[0] = IDF_EOS;
#endif

if(Cname != NULL)
 {
  free(Cname);
  Cname = NULL;
 }

}



int idf_file_opened()
{
/*
-----------------------------
   status of target idf_file
-----------------------------
*/

 return ( FtxF->iopen );
}



int idf_file_ferr()
{
/*
-----------------------------
   error in target idf_file
-----------------------------
*/

 return ( FtxF->ferr );
}



#if IDF_CPP_ == 1
int idf_file_getc(unsigned int *symb)
#else
int idf_file_getc(symb)
unsigned int *symb;
#endif
{
/*
-----------------------------
   get nex symbol from file

	returns:
	   0 - o'kay
	   1 - EOF
	   2 - ERROR
-----------------------------
*/

  unsigned int ci;
  int k=0;
  int  l;
  FILE *fil;

fil = FtxF->fd;
 ci = IDF_EOS;

  if( !FtxF->iopen )
   {
    /* file is not open */
    FtxF->ferr = 5;
    *symb = IDF_EOS;
    k=2;
   }

  else if( feof(fil) )
   {
    /* EoF is reached */
    FtxF->fend = 1;

     k=1;
     *symb = IDF_EOS;
   }

  else
   {
    ci = fgetc(fil);

    /*--------search for EOF & I/O error */

    if( ferror(fil) )
     {
      /* ERROR */
      k=2;
      FtxF->ferr = 4;
      *symb = IDF_EOS;
     }

    else if(ci==IDF_EOF)
     {
      if( feof(fil) )
       {
        FtxF->fend = 1;
        k=1;
       }
      else
       {
        k=2;
       }

      *symb = IDF_EOS;
     }

    else
     {
      *symb = ci;
     }
   }

 return k;
}



#if IDF_CPP_ == 1
int idf_file_gets(char* str, int nstr, int* flag)
#else
int idf_file_gets(str,nstr,flag)
char *str;
int   nstr;
int  *flag;
#endif
{
/*
------------------------------
   get nstr chars from file

	returns:
           2 - finished with EOF
	   1 - finished 
	   0 - not finished
	  -1 - file is not opened
	  -2 - ERROR
------------------------------
*/
  int is;
  register int i,k;
  unsigned int c;

 if(FtxF->iopen)
  {
      i=0;
      is = 0;
    while(i<nstr)
     {
         k = idf_file_getc(&c);

      if(k == 1)
       {
        is = 2;
        str[i] = IDF_EOS;
        break;
       }

      else if(k == 2)
       {
        is = -2;
        str[i] = IDF_EOS;
        break;
       }

      else if( idf_string_over(c) )
       {
        is = 1;
        str[i] = c;
        break;
       }

      else
       {
        str[i] = c;
       }

      i++;

     }/*i*/
  }

 else
  {
   FtxF->ferr = 5;
   i = 0;
   is = -1;
   str[0] = IDF_EOS;
  }

   *flag = is;

return i;
}


#if IDF_CPP_ == 1
int idf_file_ftell(long* lp)
#else
int idf_file_ftell(lp)
long *lp;
#endif
{
/*
   1 - EOF
   0 - OK
  <0 - ERROR
*/
int k;
long int l;
FILE *fil;

fil = FtxF->fd;
  k = 0;

if(FtxF->iopen)
 {
  if( feof(fil) )
   {
    FtxF->fend = 1;
    k = 1;
    l = 0L;
   }
  else
   {
       l = ftell(fil);
    if(l < 0L)
     {
      FtxF->ferr = 7;
      k = -2;
      l = 0L;
     }
   }
 }

else
 {
  /* file is not open */
  FtxF->ferr = 6;
  k = -1;
  l = 0L;
 }

*lp = l;
 return k;
}



#if IDF_CPP_ == 1
int idf_file_fseek(long k)
#else
int idf_file_fseek(k)
long k;
#endif
{
int l=0;

if(FtxF->iopen)
 {
     l = (int) fseek(FtxF->fd, k, 0);
  if(l)
   {
    l=2;
    FtxF->ferr = 9;
   }
 }

else
 {
  FtxF->ferr = 8;
  l = 1;
 }

 return l;
}



int idf_file_rewind()
{
int l=0;
FILE *fil;

fil = FtxF->fd;

if(FtxF->iopen)
 {
  rewind(fil);

  if( ferror(fil) )
   {
    clearerr(fil);
    rewind(fil);
    FtxF->fend = 0;
   }
  else
   {
    FtxF->fend = 0;
   }
 }

else
 {
  FtxF->ferr = 10;
  l = 1;
 }

 return l;
}



void idf_file_errprn()
{
#define IDF_FILE_ERR_M 11

static char *ErrF[IDF_FILE_ERR_M] = {
/* 1*/ "unable to open new file, previous one is not closed",
/* 2*/ "unable to open the file",
/* 3*/ "unable to close the unopen file",
/* 4*/ "i/o error",
/* 5*/ "unable to get symbol, file is not opened",
/* 6*/ "unable to determine position, file is not opened",
/* 7*/ "unable to determine position in file",
/* 8*/ "unable to move, file is not opened",
/* 9*/ "unable to move",
/*10*/ "unable to rewind",
/*  */ "unknown error"
                                    };

register int j;
int jerr;
char *s=NULL;

if(FtxF != NULL)
{
   jerr = FtxF->ferr;
if(jerr>0)
 {
     j = jerr;
  if(j > IDF_FILE_ERR_M) j = IDF_FILE_ERR_M;
   {
     j--;
    s = ErrF[j];

    fprintf(stdout,"I/O error=%4d %s\n",
            jerr,s);

#if IDF_FSTR == 1
    fprintf(stdout,"in file=%s\n",FtxF->name);
#endif
   }
 }
}

}


void idf_file_clearerr()
{
 if(FtxF != NULL)
  {
   FtxF->ferr = 0;
  }
}
