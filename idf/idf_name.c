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



#define IDF_SKIP_DATA_IN_NAME 0
#define IDF_VOID_UNTYPE 0


static IDF_NAME *idf_Name = NULL;



IDF_NAME* idf_name_address()
{
 return ( idf_Name );
}


int idf_name_ini()
{
int ierr=0;
int k;

if(idf_Name != NULL)
 {
  ierr=1;
 }

else
 {
            k = sizeof(IDF_NAME);
     idf_Name = (IDF_NAME*) malloc(k);
   if(idf_Name == NULL)
   {
    ierr=2;
   }
 }

return ierr;
}


void idf_name_end()
{
 if(idf_Name != NULL)
  {
   free(idf_Name);
   idf_Name = NULL;
  }
}


#if IDF_CPP_ == 1
void idf_name_begin(int* kst, int* kurs, long* pos)
#else
void idf_name_begin(kst,kurs,pos)
int  *kst,*kurs;
long *pos;
#endif
{
/*
prepare buffer for next name of data in file
*/
register int i;

    idf_Name->name_buf[0] = IDF_EOS;
    idf_Name->length      = 0;

    idf_Name->eflag       = ' ';

    idf_Name->lb          = *kst;
    idf_Name->kb          = *kurs;
    idf_Name->posb        = *pos;
    idf_Name->le          = *kst;
    idf_Name->ke          = *kurs;
    idf_Name->pose        = *pos;

    idf_Name->name_raw[0] = IDF_EOS;
    idf_Name->nlength     = 0;

    idf_Name->name[0]     = IDF_EOS;
    idf_Name->Nlength     = 0;

    idf_Name->jtype       = 0;
    idf_Name->jlocal      = 0;

    idf_Name->ndim        = 0;

    for(i=0;i<IDF_NAME_MAX_DIM;i++)
     {
    idf_Name->Mshft[i] = 0;
     }

    idf_Name->ld          = *kurs;
    idf_Name->kd          = *kst;
    idf_Name->posd        = *pos;

    idf_Name->dspos       = 0;
    idf_Name->depos       = 0;
}


void idf_name_prn()
{
 register int i;

 if(idf_Name->nlength > 0)
  {
   fprintf(stdout,"Name_RAW(%3d)=%s\n",
           idf_Name->nlength,idf_Name->name_raw);

   if(idf_Name->ndim > 0)
    {
     fprintf(stdout,"Mdim(%1d)=",idf_Name->ndim);
     for(i=0;i<idf_Name->ndim;i++)
      {
       fprintf(stdout," %5d",idf_Name->Mshft[i]);
      }
     fprintf(stdout,"\n");
    }
   else
    {
     fprintf(stdout,"ndim=%d\n",idf_Name->ndim);
    }

   fprintf(stdout,"ef=%c sb=%d kb=%d pb=%ld  se=%d ke=%d pe=%ld\n",
           idf_Name->eflag,
           idf_Name->lb,idf_Name->kb,idf_Name->posb,
           idf_Name->le,idf_Name->ke,idf_Name->pose);

   if(idf_Name->Nlength > 0)
    {
     fprintf(stdout,"Name(%3d)=%s\n",
             idf_Name->Nlength,idf_Name->name);

     fprintf(stdout,"Jtype=%d Jlocal=%d\n",
             idf_Name->jtype,idf_Name->jlocal);
    }

   fprintf(stdout,"Data: dspos=%ld depos=%ld  se=%d ke=%d pe=%ld\n",
           idf_Name->dspos,idf_Name->depos,
           idf_Name->ld,idf_Name->kd,idf_Name->posd);
  }
}



void idf_name_errprn()
{
 register int i;

 if(idf_Name->length > 0)
  {
   fprintf(stdout,"Name_BUFFER(%3d)=%s\n",
           idf_Name->length,idf_Name->name_buf);
  }

 fprintf(stdout,"Starting position=%ld string=%d kursor=%d\n",
         idf_Name->posb,idf_Name->lb,idf_Name->kb);
 fprintf(stdout,"Ending   position=%ld string=%d kursor=%d eflag=%c\n",
         idf_Name->pose,idf_Name->le,idf_Name->ke,idf_Name->eflag);


 if(idf_Name->nlength > 0)
  {
   fprintf(stdout,"Name_RAW(%3d)=%s\n",
           idf_Name->nlength,idf_Name->name_raw);

   if(idf_Name->ndim > 0)
    {
     fprintf(stdout,"Mdim(%1d)=",idf_Name->ndim);
     for(i=0;i<idf_Name->ndim;i++)
      {
       fprintf(stdout," %5d",idf_Name->Mshft[i]);
      }
     fprintf(stdout,"\n");
    }
   else
    {
     fprintf(stdout,"ndim=%d\n",idf_Name->ndim);
    }

   if(idf_Name->Nlength > 0)
    {
     fprintf(stdout,"Name(%3d)=%s\n",
             idf_Name->Nlength,idf_Name->name);

     fprintf(stdout,"Jtype=%d Jlocal=%d\n",
             idf_Name->jtype,idf_Name->jlocal);

     fprintf(stdout,"Starting Data positions=%ld %ld\n",
           idf_Name->dspos,idf_Name->depos);

     fprintf(stdout,"Ending   Data position=%ld %d %d\n",
           idf_Name->posd,idf_Name->ld,idf_Name->kd);
    }
  }
}



#if IDF_CPP_ == 1
int idf_name_buf(int* Kstr_, int* Kurs_, long* Pos_)
#else
int idf_name_buf(Kstr_,Kurs_,Pos_)
int *Kstr_,*Kurs_;
long *Pos_;
#endif
{
/*
====================================
  fulfill the buffer for data name

  flag =  1 - empty name with EOF flag
          0 - OK
         -1 - error
====================================
*/
static int jspc = IDF_NAME_SPACE;

unsigned int c,ci;
unsigned int eflag,cerr,crpt;
register int i;
int flag,ibrk,blk,is,n,k,jerr;
int kstr,kurs;
long pos;
int Kstr,Kurs;
long Kpos;
char *Name;

    Name  = idf_Name->name_buf;

    kstr  = *Kstr_;
    kurs  = *Kurs_;
     pos  = *Pos_;

    flag  = 0;

    jerr  = 0;
    cerr  = ' ';
    eflag = ' ';

    crpt  = ' ';
    blk   = 0;

    i   = 0;
    k   = 0;
    is  = 0;
    n   = 0;

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
          /* no more names */
          flag = 1;
          ibrk=1;
         }
       }
     }/*blk*/


    /*-------------
       IDF_EOS case
     --------------*/

      if( idf_string_over(ci) )
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

        if(k == 2)
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


#if IDF_SKIP_DATA_IN_NAME == 1

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
          /* improper symbol \" */
          jerr = 3;
          cerr = ci;
          flag = -1;
          ibrk=1;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* "" */


    /*-------------
       text  case 
     --------------*/

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
          /* improper symbol \# */
          jerr = 3;
          cerr = ci;
          flag = -1;
          ibrk=1;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* ## */


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
          /* improper symbol \' */
          jerr = 3;
          cerr = '\'';
          flag = -1;
          ibrk=1;
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

        if(k < 0)
         {
          /* formular error */
          flag = -1;
          ibrk=1;
         }
        else
         {
          /* improper symbol \$ */
          jerr = 3;
          cerr = ci;
          flag = -1;
          ibrk=1;
         }

        kstr = Kstr;
        kurs = Kurs;
         pos = Kpos;

       } /* $$ */

#else /* IDF_SKIP_DATA_IN_NAME */

      /*--------------------------
         improper symbols in name
        --------------------------*/

      else if(ci=='\"'||ci=='\''||ci=='#'||ci=='$')
       {
        cerr = ci;
        jerr = 3;
        flag = -1;
	ibrk=1;
       }

#endif /* IDF_SKIP_DATA_IN_NAME */


       /*---------------------------
         search for name separator 
        ---------------------------*/

      else if( idf_name_seps(ci) )
       {
        if(i > IDF_NAME_LENGTH)
	 {
          /* too long name */
	  jerr = 4;
          cerr = ci;
          flag = -1;
	  ibrk=2;
	 }
        else if(n)
	 {
          /* name done */
	  Name[i] = IDF_EOS;
	  eflag   = ci;
          flag    = 0;
	  ibrk=1;
	 }
        else
	 {
          /* empty name */
	  Name[i] = IDF_EOS;
	  eflag   = ci;

          jerr = 5;
          cerr = ci;
          flag = -1;
	  ibrk=1;
	 }

       } /* separator */


      /*-----------------
           tabulators
       -----------------*/

      else if(ci == '\t')
       {
#if IDF_TXT_TABULATORS == 1
        kurs+=7;
#endif
        if(jspc)
         {
          if(!is)
           {
	    if(i>IDF_NAME_LENGTH)
	     {
	      jerr = 6;
              flag = -1;
	      ibrk=2;
             }
	    else
             {
	      Name[i] = ' ';
	           i++;
	     }
           }
          is++;
         }
       }


       /*-----------
           vertices  
         -----------*/

      else if(ci == '\v')
       {
#if IDF_TXT_TABULATORS == 1
        kstr+=8;
#endif
        if(jspc)
         {
          if(!is)
           {
	    if(i>IDF_NAME_LENGTH)
	     {
	      jerr = 6;
              flag = -1;
	      ibrk=2;
             }
	    else
             {
	      Name[i] = ' ';
	           i++;
	     }
           }
          is++;
         }
       }


       /*--------------
           white space  
         --------------*/

      else if(ci == ' ')
       {
        if(jspc)
         {
          if(!is)
           {
	    if(i>IDF_NAME_LENGTH)
	     {
	      jerr = 6;
              flag = -1;
	      ibrk=2;
             }
	    else
             {
	      Name[i] = ' ';
	           i++;
	     }
           }
          is++;
         }
       }


       /*-------------------------------
         add ordinary character to name
        -------------------------------- */

      else if( idf_name_symbol(ci) )
       {
        n++;
        is=0;

	if(i>IDF_NAME_LENGTH)
         {
	  jerr = 6;
          flag = -1;
          ibrk=2;
	 }
	else
	 {
	  Name[i] = (char) ci;
	       i++;
	 }
       }


       /*-------------------------------
           if foreign language allowed
         add non-ASCII character to name
        -------------------------------- */

#if IDF_FOREIGN_LANGUAGE == 1
      else if( !isascii(ci) )
       {
        n++;
        is=0;

	if(i>IDF_NAME_LENGTH)
         {
	  jerr = 6;
          flag = -1;
          ibrk=2;
	 }
	else
	 {
	  Name[i] = (char) ci;
	       i++;
	 }
       }
#endif


      /*--------------------------
         improper symbols in name
        --------------------------*/

      else
       {
        cerr = ci;
        jerr = 3;
        flag = -1;
	ibrk=1;
       }

   } /* end of search */


     Name[i]          = IDF_EOS;
     idf_Name->le     = kstr;
     idf_Name->ke     = kurs;
     idf_Name->pose   = pos;
     idf_Name->length = i;
     idf_Name->eflag  = eflag;

  *Kstr_ = kstr;
  *Kurs_ = kurs;
  *Pos_  = pos;

  if(jerr)
   {
    ibrk = jerr + 42;
    idf_txt_err_put(ibrk, kstr,kurs,pos, cerr);
   }

 return flag;
}


#if IDF_CPP_ == 1
int idf_name_raw(int* kstr, int* kurs, long* kpos)
#else
int idf_name_raw(kstr,kurs,kpos)
int  *kstr,*kurs;
long *kpos;
#endif
{
static int jspc = IDF_NAME_SPACE;

int ierr;
int n,nb,iu,iuu;
char *Buf,*Name;
int *Mshft;
register int i,iname;
char c,cerr;
int kdim,idim;
int idd,id,kp,ip,is,jerr;
int Mdim[IDF_NAME_MAX_DIM];

   nb = idf_Name->length;
  Buf = idf_Name->name_buf;
 Name = idf_Name->name_raw;
Mshft = idf_Name->Mshft;

ierr = 0;
jerr = 0;
cerr = ' ';

if(nb < 1)
 {
  /*empty raw name buffer*/
  ierr=1;
  goto err;
 }

   i = strlen(Buf);
if(i != nb)
 {
  /* improper raw name buffer size */
  ierr=2;
  goto err;
 }

kdim =  0;
idim =  0;
iname=  0;
id   =  0;
kp   =  0;
ip   = -1;
is   =  0;
iu   =  0;
iuu  =  0;

for(i=0;i<nb;i++) 
 {
  c = Buf[i];

  if( isspace(c) )
   {
    if(jspc)
     {
      if(iname)
       {
        if(kp==0)
         {
          if(!is)
           {
            Name[iname] = ' ';
            iname++;
            is++;
           }
         }
       }
     }/*jspc*/
   }


  else if( isalpha(c) )
   {
    is=0;
    if(kp)
     {
      /* letter after '[' */
      cerr=c;
      jerr=3;
     }
    else if(iname)
     {
      Name[iname] = c;
      iname++;
     }
    else
     {
      if(id)
       {
        /* name starts with a digit */
        cerr=c;
        jerr=4;
       }
      else
       {
        Name[iname] = c;
        iname++;
       }
     }
   }


#if IDF_FOREIGN_LANGUAGE == 1
  else if( !isascii(c) )
   {
    is=0;
    if(kp)
     {
      /* letter after '[' */
      cerr=c;
      jerr=3;
     }
    else if(iname)
     {
      Name[iname] = c;
      iname++;
     }
    else
     {
      if(id)
       {
        /* name starts with a digit */
        cerr=c;
        jerr=4;
       }
      else
       {
        Name[iname] = c;
        iname++;
       }
     }
   }
#endif


  else if(c=='_'||c=='%')
   {
    is=0;
    if(kp)
     {
      /* '_' or '%' after '[' */
      jerr=16;
     }
    else if(iname)
     {
      Name[iname] = c;
      iname++;
     }
    else
     {
      /* name starts with '_' */
      jerr=6;
     }
   }


  else if( isdigit(c) )
   {
    is=0;
    if(kp)
     {
      /* calculate number inside [] */
      if(ip == -1)
       {
        /* improper digit position*/
        jerr=9;
       }
      else if(id)
       {
        idd = c - '0';

           kdim = kdim*10 + idd;
        if(kdim > IDF_FRMT_MAX_IDIM)
         {
          /* too large number in name dimension*/
          jerr=5;
         }
       }
      else
       {
        idd = c - '0';
        kdim = idd;
       }
      id++;
     }
    else if(iname)
     {
      /* digit is a part of name */
      Name[iname] = c;
      iname++;
     }
    else
     {
      /* name starts with a digit */
      cerr=c;
      jerr=4;
     } 
   }


  else if(c == '[')
   {
    is=0;
    iu++;
    kp++;
    if(iuu)
     {
      /* mixture [] and () */
      jerr=15;
     }
    else if(kp > IDF_NAME_MAX_DIM)
     {
      /* too many dimensions*/
      jerr=7;
     }
    else if(iname==0)
     {
      /* name starts with '[' */
      jerr=8;
     }
    else if(id)
     {
      /* digit before '[' */
      jerr=9;
     }
    else if(ip != -1)
     {
      /* improper use of '[' */
      jerr=12;
     }
    else
     {
      ip++;
     }
   }/* '[' */


  else if(c == '(')
   {
    is=0;
    iuu++;
    kp++;
    if(iu)
     {
      /* mixture [] and () */
      jerr=15;
     }
    else if(kp > IDF_NAME_MAX_DIM)
     {
      /* too many dimensions*/
      jerr=7;
     }
    else if(iname==0)
     {
      /* name starts with '[' */
      jerr=8;
     }
    else if(id)
     {
      /* digit before '[' */
      jerr=9;
     }
    else if(ip != -1)
     {
      /* improper use of '[' */
      jerr=12;
     }
    else
     {
      ip++;
     }
   }/* '(' */


  else if(c == ']')
   {
    is=0;
    iu++;
    if(iuu)
     {
      /* mixture [] and () */
      jerr=15;
     }
    else if(iname==0)
     {
      /* name starts with  ']' */
      jerr=8;
     }
    else if(kp)
     {
      if(ip==0)
       {
        if(idim < IDF_NAME_MAX_DIM)
         {
          if(id==0) kdim=0;
          Mdim[idim] = kdim;
               idim++;

          ip--;
          id=0;
         }
        else
         {
          /* too many dims*/
          jerr=7;
         }
       }
      else
       {
        /* extra ']' detected */
        jerr=11;
       }
     }
    else
     {
      /* improper appearence of  ']' */
      jerr=10;
     }
   } /* ] */


  else if(c == ')')
   {
    is=0;
    iuu++;
    if(iu)
     {
      /* mixture [] and () */
      jerr=15;
     }
    else if(iname==0)
     {
      /* name starts with  ']' */
      jerr=8;
     }
    else if(kp)
     {
      if(ip==0)
       {
        if(idim < IDF_NAME_MAX_DIM)
         {
          if(id==0) kdim=0;
          Mdim[idim] = kdim;
               idim++;

          ip--;
          id=0;
         }
        else
         {
          /* too many dims*/
          jerr=7;
         }
       }
      else
       {
        /* extra ']' detected */
        jerr=11;
       }
     }
    else
     {
      /* improper appearence of  ']' */
      jerr=10;
     }
   }/* ) */


  else if(c == ',')
   {
    is=0;
    if(iu>0 || iname==0)
     {
      /* improper comma position */
      cerr=c;
      jerr=16;
     }
    else
     {
      if(kp==0 || ip==-1)
       {
        /* improper appearence of  ',' */
        cerr=c;
        jerr=16;
       }
      else
       {
        if(idim < IDF_NAME_MAX_DIM)
         {
          if(id==0) kdim=0;
          Mdim[idim] = kdim;
               idim++;
          id=0;
         }
        else
         {
          /* too many dims*/
          jerr=7;
         }
       }
     }
   }/* , */


  else
   {
    /* improper symbol in name */
    cerr=c;
    jerr=16;
   }

  if(jerr) break;

 } /*i */


if(jerr==0)
 {
  if(ip != -1)
   {
    /* unclosed '[' */
    jerr=13;
   }
  else if(jspc)
   {
    if(iname)
    {
     i=iname;
     while(i>0)
      {
       i--;
       if(Name[i] != ' ') break;
      }
          iname = i+1;
     Name[iname] = IDF_EOS;
    }

    if(iname<1)
     {
      /*empty name */
      jerr=14;
     }
   }

  else
   {
    if(iname<1)
     {
      jerr=14;
     }
    else if(ip != -1)
     {
      jerr=13;
     }
   }
 }


 Name[iname] = IDF_EOS;


if(jerr==0)
 {
  if(idim)
   {

    if(iuu)
     {
      /* convert dims FORTRAN style to C */
      i=idim;
      iname=0;
      while(i>0)
       {
        i--;

        Mshft[iname] = Mdim[i];
              iname++;
       }
     }

    else
     {
      /* C style dims */
      for(i=0;i<idim;i++)
       {
        Mshft[i] = Mdim[i];
       }
     }
   }
 }

  idf_Name->nlength = iname;
  idf_Name->ndim    = idim;

if(jerr)
 {
  ierr=jerr+2;
 }

err:
if(ierr)
 {
  kp = ierr + 48;
  idf_txt_err_put(kp, *kstr,*kurs,*kpos, cerr);
 }

return ierr;
}



#if IDF_CPP_ == 1
int idf_name_unit(char* Buf, int ksbuf, int nend,
                  char* Name, int* nam, int* kebuf)
#else
int idf_name_unit(Buf,ksbuf,nend, Name,nam,kebuf)
char *Buf;
int   ksbuf;
int   nend;
char *Name;
int  *nam;
int  *kebuf;
#endif
{
/*
----------------------------------------
          get name unit 
if name divided by space into sub-names
----------------------------------------
returns: 2 - no name
         1 - last name
         0 - sub-name
        -1 - alg error
*/
int nn,nl;
int k,kk,ke;
register int i,j;
char c;
int flag;

flag = 0;
   k = 0;
   i = ksbuf;
  ke = ksbuf;
  nn = nend;
   j = 0;

if(ke >= nn)
 {
  /* no name */
  Name[0] = IDF_EOS;
  nl=0;
  ke=i;
  flag = 2;
 }

else
 {
  /* skip spaces before the name */
  while(i<nn)
   {
       c = Buf[i];
    if(c != ' ') break;
    i++;
   }

  if(i >= nn)
   {
    /* empty name */
    Name[0] = IDF_EOS;
    nl=0;
    ke=i;
    flag = 2;
   }

  else
   {
    while(i<nn)
     {
         c = Buf[i];
      if(c == ' ')
       {
        k=1;
       }
      else
       {
        Name[j] = c;
             j++;
       }
      
      i++;

      if(k) break;
     }/*i*/

    if(j==0)
     {
      /* algorithm error*/
      Name[0] = IDF_EOS;
      nl=0;
      ke=i;
      flag = -1;
     }

    else if(k)
     {
      if(i == nn)
       {
        /* last name */
        Name[j] = IDF_EOS;
        nl = j;
        ke=i;
        flag=1;
       }

      else
       {
        /* composition of sub-names is posible*/

        Name[j] = IDF_EOS;
        nl = j;
        ke=i;

        kk=0;
        for(j=i;j<nn;j++)
         {
             c = Buf[j];
          if(c != ' ')
           {
            kk=1;
            break;
           }
         }

        if(kk)
         {
          flag=0;
         }
        else
         {
          ke = nn;
          flag=1;
         }
       }
     }

    else
     {
      /* last name */
      Name[j] = IDF_EOS;
      nl = j;
      ke=i;
      flag=1;
     }

   }/*i<nn*/
 }/*ke<nn*/

*nam = nl;
*kebuf = ke;

return flag;
}



#if IDF_CPP_ == 1
int idf_name_compose(char* Buf, int n, char* Name, int* nam,
                     int* Jtype, int* Jlocal)
#else
int idf_name_compose(Buf,n, Name,nam, Jtype,Jlocal)
char *Buf;
int   n;
char *Name;
int  *nam;
int  *Jtype;
int  *Jlocal;
#endif
{
/*
--------------------------------
   get name and its attributes
--------------------------------
returns:
        0 - OK
        1 - alg error #1
        2 - empty name
        3 - name contains only keywords
        4 - too many keywords assigned to name
        5 - white space splits name
        6 - incompatible keywords in name
*/
int ibrk,k,flag,kk,kt;
int ib,ie;
int jlocal,jtype,ktype;

  flag = 0;
jlocal = 0;
jtype  = IDF_FRMT_VOID;

ib=0;
ie=0;
kk=0;
      ibrk=0;
while(ibrk==0)
 {
  /* get unit */
     k = idf_name_unit(Buf,ib,n, Name,nam,&ie);

  if(k < 0)
   {
    /* error in getting name unit*/
    flag = 1;
    ibrk=2;
   }

  else if(k == 2)
   {
    /* empty name */
    flag = 2;
    ibrk=2;
   }

  else if(k == 1)
   {
    /* last name */
       ktype = idf_keyw(Name);
    if(ktype != 0)
     {
      /* name contains only keywords*/
      flag = 3;
      ibrk=2;
     }
    else
     {
      /* name received */
      ibrk=1;
     }
   }

  else if(kk > IDF_KEYW_MAX_IN_NAME)
   {
    flag = 4;
    ibrk=2;
   }

  else
   {
    /* sub-name */
       ktype = idf_keyw(Name);
    if(ktype==IDF_FRMT_LOCAL)
     {
      jlocal=1;
     }
    else if(ktype==0)
     {
      /* white space splits name */
      flag = 5;
      ibrk=2;
     }
    else
     {
      /* type*/
         kt = idf_keyw_cnv(jtype,ktype);
      if(kt)
       {
        jtype = kt;
       }
      else
       {
        /* improper convolution of type keywords */
        flag = 6;
        ibrk=2;
       }

      kk++;
     }

    ib = ie;
   }

 }/* ibrk*/

*Jtype = jtype;
*Jlocal = jlocal;

return flag;
}



#if IDF_CPP_ == 1
int idf_name(int* kstr, int* kurs, long* kpos)
#else
int idf_name(kstr,kurs,kpos)
int  *kstr;
int  *kurs;
long *kpos;
#endif
{
static int jspc = IDF_NAME_SPACE;

int ierr;
int n,nn,nl;
char *Buf,*Name;
register int i,j;
char c;
unsigned int cerr;
int k,kk;
int jtype,jlocal;


    n = idf_Name->nlength;
  Buf = idf_Name->name_raw;
 Name = idf_Name->name;

  ierr = 0;
    nl = 0;
 jtype = IDF_FRMT_VOID;
jlocal = 0;
  cerr = ' ';

if(n < 1)
 {
  /* empty name*/
  ierr=1;
  goto err;
 }

   i = strlen(Buf);
if(i != n)
 {
  /* improper name */
  ierr=2;
  goto err;
 }


if(jspc == 0)
 {
  /* ignore white spaces */
    j=0;
  for(i=0;i<n;i++)
   {
       c = Buf[i];
    if(c != ' ')
     {
      Name[j] = c;
           j++;
     }
   }/*i*/

  if(j)
   {
    /* single name*/
    Name[j] = IDF_EOS;
    nl = j;

       jtype = idf_keyw(Name);
    if(jtype)
     {
      /* name is a keyword */
      ierr=6;
     }
   }
  else
   {
    /* empty name */
    Name[j] = IDF_EOS;
    nl = j;
    ierr=3;
   }

 }/*jspc=0*/



else if(jspc == 1)
 {
  k=0;
  i=0;
  j=0;

  while(i<n)
   {
       c = Buf[i];
    if(c != ' ') break;
    i++;
   }

  if(i < n)
   {
    while(i<n)
     {
         c = Buf[i];
      if(c == ' ')
       {
        k=1;
        break;
       }
      else
       {
        Name[j] = c;
             j++;
       }
      
      i++;
     }
   }/*i<n*/

  if(k)
   {
    /* white space inside the name */
    Name[0] = IDF_EOS;
    nl=0;
    ierr=4;
   }
  else if(j)
   {
    /* single name */
    Name[j] = IDF_EOS;
    nl = j;

       jtype = idf_keyw(Name);
    if(jtype)
     {
      /* name is a keyword */
      ierr=6;
     }
   }
  else
   {
    /* no name */
    Name[0] = IDF_EOS;
    nl=0;
    ierr=7;
   }

 }/*jspc=1*/



else
 {
  /* skip multiple spaces */
  k=0;
  j=0;
  for(i=0;i<n;i++)
   {
       c = Buf[i];
    if(c == ' ')
     {
      if(k==0)
       {
        Buf[j] = c;
            j++;
        k++;
       }
     }
    else
     {
      Buf[j] = c;
          j++;
      k=0;
     }
   }
  nn = j;
  Buf[j] = IDF_EOS;

          kk = idf_name_compose(Buf,nn, Name,&nl,
                               &jtype,&jlocal);
       if(kk==1) ierr=11;
  else if(kk==2) ierr=12;
  else if(kk==3) ierr=7;
  else if(kk==4) ierr=14;
  else if(kk==5) ierr=8;
  else if(kk==6) ierr=13;

 }/*jspc*/



  idf_Name->Nlength = nl;
  idf_Name->jtype   = jtype;
  idf_Name->jlocal  = jlocal;


if(ierr==0)
 {
  if(jlocal)
   {
       kk = idf_Name->ndim;
    if(kk)
     {
      /* locals must not have dimensions */
      ierr=17;
     }

    else
     {
         k = idf_keyw_digit(jtype);
      if(k==0)
       {
        /* locals must be digital, not char and string */
        ierr=15;
       }
     }
   }
 }

err:
if(ierr)
 {
  k = ierr + 64;
  idf_txt_err_put(k, *kstr,*kurs,*kpos, cerr);
 }

return ierr;
}



#if IDF_CPP_ == 1
int idf_name_dpos(int kse)
#else
int idf_name_dpos(kse)
int kse;
#endif
{
long l;
int k;

    k = idf_file_ftell(&l);
if(!k)
 {
  if(kse)
   {
    idf_Name->depos = l;
   }
  else
   {
    idf_Name->dspos = l;
   }
 }

return k;
}



#if IDF_CPP_ == 1
int idf_name_get(int* Kstr, int* Kurs, long* Pos, int* Local)
#else
int idf_name_get(Kstr,Kurs,Pos,Local)
int  *Kstr,*Kurs;
long *Pos;
int  *Local;
#endif
{
/*
 Create the IDF_NAME structure
     >0 - error
      0 - OK
     -1 - named data followed by EoF
     -2 - empty name with EoF flag
*/

int ierr;
int jerr,flag,jlocal;
unsigned int cf,cerr;

flag = 0;
ierr = 0;
cerr = ' ';

 idf_name_begin(Kstr,Kurs,Pos);
 jlocal=0;

/*---------------
   get the name 
  --------------*/

   jerr = idf_name_buf(Kstr,Kurs,Pos);
if(jerr==1)
 {
  /* empty name with EOF flag */
  flag = -2;
  goto err;
 }
else if(jerr<0)
 {
  flag=1;
  goto err;
 }
 
   jerr = idf_name_raw(Kstr,Kurs,Pos);
if(jerr)
 {
  flag=2;
  goto err;
 }

   jerr = idf_name(Kstr,Kurs,Pos);
if(jerr)
 {
  flag=3;
  goto err;
 }


/*-------------------------------
  real starting position of data 
       from file beginning 
  -------------------------------*/

   jerr = idf_name_dpos(0);
if(jerr==1)
 {
  /* empty named data */
  ierr=1;
  flag=4;
 }
else if(jerr<0)
 {
  flag=5;
  goto err;
 }


/*-------------------------------
  search and check the data field 
  -------------------------------*/

  cf = idf_Name->eflag;

   jerr = idf_skip_data(Kstr,Kurs,Pos, &cf);

  idf_Name->ld    = *Kstr;
  idf_Name->kd    = *Kurs;
  idf_Name->posd  = *Pos;

if(jerr<0)
 {
  flag=6;
  goto err;
 }
else if(jerr==1)
 {
  /* empty named data near EoF*/
  ierr=2;
  flag=7;
  goto err;
 }


/*---------------------
   local/nonlocal flag 
 ----------------------*/ 

  jlocal = idf_Name->jlocal;


/*-------------------------------
  real ending position of data 
       from file beginning 
  -------------------------------*/

   jerr = idf_name_dpos(1);
if(jerr==1)
 {
  /* empty named data followed by eOF*/
  ierr=2;
  flag = -1;
 }
else if(jerr<0)
 {
  flag=8;
 }


err:
 *Local = jlocal;

if(ierr)
 {
  jerr = ierr + 81;
  idf_txt_err_put(jerr, *Kstr,*Kurs,*Pos, cerr);
 }

return flag;
}
