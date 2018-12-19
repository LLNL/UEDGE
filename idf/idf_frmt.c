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




IDF_FRMT *Frmt=NULL;



IDF_FRMT* idf_frmt_address()
{
 return ( Frmt );
}


int idf_frmt_ini()
{
int ierr=0;
long k;

if(Frmt != NULL)
 {
  ierr=1;
  goto err;
 }

      k = sizeof(IDF_FRMT);
   Frmt = (IDF_FRMT*) malloc(k);
if(Frmt == NULL)
 {
  ierr=2;
  goto err;
 }


   k = idf_frmt_table_prp();
if(k)
 {
  ierr=3;
 }

err:
return ierr;
}


void idf_frmt_end()
{
if(Frmt != NULL)
 {
  idf_frmt_clean();

  free(Frmt);
  Frmt = NULL;
 }
}


void idf_frmt_clean()
{
register int i;
int *Iarr;

if(Frmt != NULL)
 {
  Iarr = Frmt->Mrpt;
  for(i=0;i<IDF_FRMT_ELEM_NMAX;i++)
   {
    Iarr[i] = 0;
   }

  Iarr = Frmt->Krpt;
  for(i=0;i<IDF_FRMT_ELEM_NMAX;i++)
   {
    Iarr[i] = 0;
   }

  for(i=0;i<IDF_FRMT_ELEM_NMAX;i++)
   {
    Frmt->Rfrmt[i] = ' ';
   }
  Frmt->Rfrmt[0] = IDF_EOS;

  Frmt->NfrmtE = 0;

  for(i=0;i<IDF_FRMT_UNIT_NMAX;i++)
   {
    Frmt->Ufrmt[i].type   = 0;
    Frmt->Ufrmt[i].jsize  = 0;
    Frmt->Ufrmt[i].jalloc = 0;
    Frmt->Ufrmt[i].pointer= 0;
    Frmt->Ufrmt[i].jmem   = 0;

    Frmt->Ufrmt[i].Ndim = 0;    
   }

  Frmt->NfrmtU = 0;

  Frmt->jerr = 0;
  Frmt->kpos = 0;
 }
}


#if IDF_CPP_ == 1
int idf_frmt_prp(char* str, int n)
#else
int idf_frmt_prp(str,n)
char *str;
int   n;
#endif
{
int ierr=0;
int jbrk;
int pos;
int jtype,jsize,jalloc,jmem,ipntr;
int idigt,irpt,jrpt,idim,dimtype;
int ic,k,flag;
char c,cc,cntr,cdo;
register int i,ibrk;
int j,ju;
int iopen,jopen;
int Mdim[IDF_FRMT_DIM_MAX];
char *Rfrmt;
int  *Mrpt,*Krpt;

Rfrmt = Frmt->Rfrmt;
 Mrpt = Frmt->Mrpt;
 Krpt = Frmt->Krpt;

ipntr = 0;
jalloc= 0;
jtype = 0;
idigt = 0;
ic    = 0;
cntr  = '{';

    i = 0;
    j = 0;
   ju = 0;
iopen = 0;
  pos = 0;
 flag = 0;
 irpt = 0;
 jrpt = 0;
 idim = 0;

   k=0;
cntr= ' ';
cdo = ' ';

if(n < 1)
 {
  /* empty format string*/
  ierr=46;
  goto err;
 }

      jbrk=0;
while(jbrk==0)
 {
  /*-----------------------
     seach for format unit
    ----------------------*/

  if(i<n)
   {
    c = str[i];
            i++;
            pos++;

    if( !isascii(c) )
     {
      /* improper format symbol*/
      ierr=10;
     }

    else if( isdigit(c) )
     {
      if(idigt)
       {
           ic = ic*10 + (c - '0');
        if(ic > IDF_FRMT_MAX_IRPT)
         {
          /* large number */
          ierr=8;
         }
       }
      else
       {
        ic = c - '0';
        idigt=1;
        cntr = '7';
       }
     }/* digit */


    else if( isalpha(c) )
     {
      /* can be format symbol without % or # */

           jtype = idf_frmt_ctype(c);
        if(jtype)
         {
          /* format unit up !!! */
          if(idigt)
           {
            /* repeatition number*/
            if(ic > 0)
             {
              irpt=ic;
              jrpt=1;
              flag=1;
             }
            else
             {
              /*improper repeatition number for format symbol*/
              ierr=11;
             }
           }
          else
           {
            /* single target*/
            irpt=1;
            jrpt=1;
            flag=1;
           }

          idigt=0;
          cntr='%';
          cdo = '%';
         }

        else
         {
          /* improper format symbol*/
          ierr=10;
         }
     } /* letter*/


    else if(c == '{')
     {
      k=1;
      if(idigt)
       {
        if(ic < 1)
         {
          /* zero repeatition number */
          ierr=16;
          k=0;
         }
        else
         {
          idigt=0;
          irpt=ic;
         }
       }
      else
       {
        irpt=1;
       }

      if(k)
       {
        if(j < IDF_FRMT_ELEM_NMAX)
         {
          iopen++;

          Rfrmt[j] = '{';
          Mrpt [j] = irpt;
          Krpt [j] = iopen;
                j++;

          cntr = '{';
         }
        else
         {
          /* too many format elements*/
          ierr=17;
         }
       }
     }/* { */


    else if(c == '}')
     {
           if(cntr=='7') k=1;
      else if(cntr==',') k=1;
      else if(cntr=='{') k=1;
      else if(iopen==0)  k=2;
      else               k=0;

      if(k==1)
       {
        /* improper } */
        ierr=18;
       }
      else if(k==2)
       {
        /* extra } */
        ierr=19;
       }
      else
       {
        ibrk=j;
        k=0;
        while(ibrk>0)
         {
          ibrk--;
             cc = Rfrmt[ibrk];
          if(cc == '{')
           {
            if(Krpt[ibrk] == iopen)
             {
              k=1;
              break;
             }
           }
         }

        if(k)
         {
          if(j < IDF_FRMT_ELEM_NMAX)
           {
            Rfrmt[j] = '}';
            Mrpt [j] = ibrk;
            Krpt [j] = -iopen;
                  j++;

            iopen--;
            cntr = '}';
           }

          else
           {
            ierr=17;
           }
         }
        else
         {
          /* return point was not found */
          ierr=20;
         }
       }              

     }/* } */


    else if(c=='%' || c=='#' || c=='&' || c=='@')
     {
      /* format unit body*/
      if(cntr=='7')
       {
        /* repeatition number for a unit*/
        irpt=ic;
       }
      else
       {
        irpt=1;
       }

      if(c=='&')
       {
        cntr = '%';
        cdo  = '%';
        ipntr=1;
       }
      else if(c=='@')
       {
        cntr = '%';
        cdo  = '%';
        ipntr=1;
        jalloc=1;
       }
      else if(c=='#')
       {
        cntr = '#';
        cdo  = '#';
        ipntr=0;
       }
      else
       {
        cntr = '%';
        cdo  = '%';
        ipntr=0;
       }

      idigt=0;
      jrpt =1;

           ibrk=0;
     while(ibrk==0)
      {
       if(i<n)
        {
          c = str[i];
                  i++;
                  pos++;

          if(c == '&')
           {
            /* pointer flag received*/
            if(ipntr)
             {
              /* multiple pointer symbol in format body */
              ierr=13;
             }
            else if(idigt)
             {
              /* improper position of pointer symbol */
              ierr=7;
             }
            else
             {
              ipntr=1;
             }
           }

          else if(c == '@')
           {
            /* pointer flag received*/
            if(ipntr)
             {
              /* multiple pointer symbol in format body */
              ierr=13;
             }
            else if(idigt)
             {
              /* improper position of pointer symbol */
              ierr=7;
             }
            else
             {
              ipntr=1;
              jalloc=1;
             }
           }

          else if(c == '#')
           {
            if(cntr == '%')
             {
              cdo = '#';
             }
            else if(idigt)
             {
              /* improper position of # symbol */
              ierr = 7;
             }
           }

          else if( isalpha(c) )
           {
               jtype = idf_frmt_ctype(c);
            if(jtype)
             {
              /* format unit up */
              if(idigt)
               {
                if(ic > 0)
                 {
                  /* ic is repeatition or sizeof number*/
                  ibrk=1;
                  jrpt=ic;
                 }
                else
                 {
                  /*zero number in format unit body*/
                  ierr=12; 
                 }
               }
              else
               {
                /* single target*/
                ic=1;
                ibrk=1;
                jrpt=1;
               }
             }
            else
             {
              /* improper letter*/
              ierr=1;
             }
           }/* letter*/


          else if( isdigit(c) )
           {
            if(idigt)
             {
              if(idigt<IDF_FRMT_MAX_IRPT)
               {
                ic = ic*10 + (c - '0');
                idigt++;
               }
              else
               {
                /* too many digits in format unit body*/
                ierr=6;
               }
             }
            else
             {
              ic = c -'0';
              idigt=1;
             }
           }/*digit*/


          else if(c == ',')
           {
            /* comma in format unit body */
            ierr=2;
           }

          else if(c=='{' || c=='}')
           {
            /* { or } in format unit body */
            ierr=4;
           }

          else if(c=='[' || c==']' || c=='(' || c==')')
           {
            /* [] or () in format unit body */
            ierr=21;
           }

          else if(c == ' ')
           {
            ;
           }

          else if(c == '\t')
           {
            /* skip tabulator */
            pos+=7;
           }

          else if(c == '\v')
           {
            /* improper symbol*/
            ierr=5;
           }

          else
           {
            /* improper symbol */
            ierr=5;
           }

         }/*i<n*/

        else
         {
          /* encomplete format unit at the end of string*/
          ierr=3;
         }/*i>=n*/

        if(ierr) ibrk=2;

       }/*ibrk*/

      if(ibrk==1)
       {
        if(irpt<1 || jrpt<1)
         {
          /*improper repeatition number for unit*/
          ierr=9;
         }
        else
         {
          flag=1;
         }
       }

     }/* % # & @ */


    else if(c==',')
     {
           if(cntr=='{') k=1;
      else if(cntr==',') k=1;
      else if(cntr=='7') k=1;
      else               k=0;
      if(k)
       {
        /* improper appearence of comma*/
        ierr=15;
       }
      else
       {
        cntr = ',';
       }
     }

    else if(c=='[' || c==']' || c=='(' || c==')')
     {
      /* [] or () in format string */
      ierr=22;
     }

    else if(c == ' ')
     {
      ;
     }

    else if(c == '\t')
     {
      pos+=7;
     }

    else if(c == '\v')
     {
      /* improper symbol in format string*/
      ierr=14;
     }

    else if(c == ';')
     {
      /* end of format string*/

      if(cntr=='7')
       {
        /* improper digit position*/
        ierr=23;
       }
      else if(cntr==',')
       {
        ierr=15;
       }
      else if(cntr=='{')
       {
        if(i==1)
         {
          /* empty format string */
          ierr=46;
         }
        else
         {
          /* unclosed { */
          ierr=24;
         }
       }
      else if(iopen)
       {
        /* unclosed { period*/
        ierr=41;
       }
      else
       {
        /* no more elements*/
        jbrk=1;
       }
     }

    else
     {
      /*improper symbol in format string*/
      ierr=14;
     }

   }/*i<n*/

  else
   {
    if(cntr=='7')
     {
      /* improper digit position*/
      ierr=23;
     }
    else if(cntr==',')
     {
      ierr=15;
     }
    else if(cntr=='{')
     {
      /* unclosed { */
      ierr=24;
     }
    else if(iopen)
     {
      /* unclosed { period*/
      ierr=41;
     }
    else
     {
      /* no more elements*/
      jbrk=1;
     }

   }/*i>=n*/


    /*-----------------------
       search for dimensions 
      -----------------------*/

    if(jbrk==0 && ierr==0)
     {
      if(flag)
       {
        if(i>=n)
         {
          /* no dims*/
          jmem  = jrpt;

          idim  = 1;
          Mdim[0] = jrpt;

          if(jtype==IDF_FRMT_S || jtype==IDF_FRMT_T)
           {
            jsize = jrpt;
           }
          else
           {
            jsize = 1;
           }

          if(irpt<1 || jsize<1 || jmem<1)
           {
            /* improper dimensions*/
            ierr=42;
           }

          else if(j<IDF_FRMT_ELEM_NMAX)
           {
            Rfrmt[j] = cdo;
            Mrpt [j] = ju;
            Krpt [j] = -1;
                  j++;

            if(ju<IDF_FRMT_UNIT_NMAX)
             {
              Frmt->Ufrmt[ju].Nrep    = irpt;
              Frmt->Ufrmt[ju].type    = jtype;
              Frmt->Ufrmt[ju].pointer = ipntr;
              Frmt->Ufrmt[ju].jsize   = jsize;
              Frmt->Ufrmt[ju].jmem    = jmem;
              Frmt->Ufrmt[ju].jalloc  = jalloc;
              Frmt->Ufrmt[ju].Ndim    = idim;
              Frmt->Ufrmt[ju].Mdim[0] = Mdim[0];
                          ju++;
             }
            else
             {
              /* too many format units */
              ierr=25;
             }
           }
          else
           {
            ierr=17;
           }
         }/*i>=n*/

        else
         {
          idim   =0;
          dimtype=0;
          jopen  =0;
          idigt  =0;

                ibrk=0;
          while(ibrk==0)
           {
            if(i<n)
             {
              c = str[i];
                      i++;
                      pos++;

              if(c=='[')
               {
                if(idigt)
                 {
                  /*improper position of [*/
                  ierr=29;
                 }
                else if(dimtype==1)
                 {
                  if(jopen)
                   {
                    /*improper sequence of [] of () */
                    ierr=27;
                   }
                  else
                   {
                    jopen++;
                   }
                 }
                else if(dimtype==2)
                 {
                  /* mixed [] with () */
                  ierr=26;
                 }
                else
                 {
                  dimtype=1;
                  jopen++;
                 }
               } /* [ */


              else if(c=='(')
               {
                if(idigt)
                 {
                  /*improper position of (*/
                  ierr=30;
                 }
                else if(dimtype==2)
                 {
                  if(jopen)
                   {
                    /*improper sequence of [] of () */
                    ierr=27;
                    ibrk=2;
                   }
                  else
                   {
                    jopen++;
                   }
                 }
                else if(dimtype==1)
                 {
                  /* mixed [] with () */
                  ierr=26;
                 }
                else
                 {
                  dimtype=2;
                  jopen++;
                 }
               }/* ( */


              else if(c==']')
               {
                if(dimtype==2)
                 {
                  ierr=26;
                 }
                else if(jopen)
                 {
                  if(idigt)
                   {
                    if(idim < IDF_FRMT_DIM_MAX)
                     {
                      Mdim[idim] = ic;
                           idim++;

                      idigt=0;
                      jopen--;
                     }
                    else
                     {
                      /*too many dimensions*/
                      ierr=32;
                     }
                   }
                  else
                   {
                    /*empty dimension*/
                    ierr=31;
                   }
                 }
                else
                 {
                  /* extra ] */
                  ierr=28;
                 }
               }/* ] */


              else if(c==')')
               {
                if(dimtype==1)
                 {
                  ierr=26;
                 }
                else if(jopen)
                 {
                  if(idigt)
                   {
                    if(idim < IDF_FRMT_DIM_MAX)
                     {
                      Mdim[idim] = ic;
                           idim++;

                      idigt=0;
                      jopen--;
                     }
                    else
                     {
                      /*too many dimensions*/
                      ierr=32;
                     }
                   }
                  else
                   {
                    /*empty dimension*/
                    ierr=31;
                   }
                 }
                else
                 {
                  /* extra ) */
                  ierr=33;
                 }
               }/* ) */


              else if(c==',')
               {
                if(jopen)
                 {
                  if(dimtype==1)
                   {
                    /*improper comma*/
                    ierr=34;
                   }
                  else if(idigt)
                   {
                    if(idim < IDF_FRMT_DIM_MAX)
                     {
                      Mdim[idim] = ic;
                           idim++;

                      idigt=0;
                     }
                    else
                     {
                      /*too many dimensions*/
                      ierr=32;
                     }
                   }
                  else
                   {
                    /*improper comma*/
                    ierr=34;
                   }
                 }
                else if(idigt)
                 {
                  /* improper comma */
                  ierr=35;
                 }
                else
                 {
                  /* format unit ends */
                  cntr = ',';
                  ibrk=1;
                 }
               }/*,*/


              else if( isdigit(c) )
               {
                if(jopen)
                 {
                  if(idigt)
                   {
                       ic = ic*10 + (c - '0');
                    if(ic > IDF_FRMT_MAX_IDIM)
                     {
                      /*too large number in format unit dimension*/
                      ierr=36;
                     }
                   }
                  else
                   {
                    ic = c - '0';
                    idigt=1;
                   }
                 }
                else
                 {
                  /* new unit starts with a digit */
                  i--; /* unget*/
                  pos--;
                  ibrk=1;
                 }
               }/* digit*/


              else if(c == ' ')
               {
                ;
               }

              else if(c=='\t')
               {
                pos+=7;
               }

              else if(c=='\v')
               {
                ierr=37;
               }


              else
               {
                if(jopen)
                 {
                  /*improper symbol in dimensions*/
                  ierr=37;
                 }
                else
                 {
                  /* new unit */
                  i--; /*unget*/
                  pos--;
                  ibrk=1;
                 }
               }

             }/*i<n*/

            else
             {
              if(jopen)
               {
                /*unfinished dimensions*/
                ierr=38;
               }
              else if(idigt)
               {
                /* encomplete format unit */
                ierr=39;
               }
              else if(iopen)
               {
                /* unclosed period*/
                ierr=41;
               }
              else
               {
                /* end of format*/
                ibrk=1;
               }
             }/*i>=n*/

            if(ierr) ibrk=2;

           }/*ibrk*/


          if(ierr==0)
           {
            jsize = 0;
            jmem  = 0;

            if(idim==1)
             {
              if(jtype==IDF_FRMT_T || jtype==IDF_FRMT_S)
               {
                jsize = Mdim[0];
                jmem  = jrpt*jsize;
               }
              else
               {
                jsize = 1;
                jmem  = jrpt*Mdim[0];
               }
             }

            else if(idim>1)
             {
              if(jtype==IDF_FRMT_T || jtype==IDF_FRMT_S)
               {
                if(dimtype==1)
                 {
                  Mdim[0] = Mdim[0]*jrpt;
                             ibrk=idim-1;
                  jsize=Mdim[ibrk];
                  jmem = jsize;
                 }
                else if(dimtype==2)
                 {
                  jsize = Mdim[0];
                              ibrk = idim-1;
                  jmem = Mdim[ibrk]*jrpt;
                  Mdim[ibrk] = jmem;
                 }
               }
              else
               {
                jsize = 1;
                if(dimtype==1)
                 {
                  Mdim[0] = Mdim[0]*jrpt;
                            ibrk=idim-1;
                  jmem=Mdim[ibrk];
                 }
                else if(dimtype==2)
                 {
                              ibrk = idim-1;
                  jmem = Mdim[ibrk]*jrpt;
                  Mdim[ibrk] = jmem;
                 }
               }

              while(ibrk>0)
               {
                ibrk--;
                jmem = jmem*Mdim[ibrk];
               }
             }

            else
             {
              /* no dims*/
              dimtype = 1;
                jmem  = jrpt;

                 idim = 1;
              Mdim[0] = jrpt;

              if(jtype==IDF_FRMT_T || jtype==IDF_FRMT_S)
               {
                jsize = jrpt;
               }
              else
               {
                jsize = 1;
               }
             }


            if(irpt<1 || jsize<1 || jmem<1)
             {
              /* improper dimensions*/
              ierr=42;
             }

            else if(j<IDF_FRMT_ELEM_NMAX)
             {
              Rfrmt[j] = cdo;
              Mrpt [j] = ju;
              Krpt [j] = -1;
                    j++;

              if(ju<IDF_FRMT_UNIT_NMAX)
               {
                Frmt->Ufrmt[ju].Nrep   = irpt;
                Frmt->Ufrmt[ju].type   = jtype;
                Frmt->Ufrmt[ju].pointer= ipntr;
                Frmt->Ufrmt[ju].jsize  = jsize;
                Frmt->Ufrmt[ju].jmem   = jmem;
                Frmt->Ufrmt[ju].jalloc = jalloc;
                Frmt->Ufrmt[ju].Ndim   = idim;

                if(dimtype==1)
                 {
                  for(ibrk=0;ibrk<idim;ibrk++)
                   {
                    Frmt->Ufrmt[ju].Mdim[ibrk] = Mdim[ibrk];
                   }
                 }
                else if(dimtype==2)
                 {
                  ibrk=idim;
                  k=0;
                  while(ibrk>0)
                   {
                    ibrk--;
                    Frmt->Ufrmt[ju].Mdim[k] = Mdim[ibrk];
                    k++;
                   }
                 }
                else
                 {
                  /* improper dimtype*/
                  ierr=43;
                 }
                ju++;
               }

              else
               {
                /* too many format units */
                ierr=25;
               }
             }
            else
             {
              ierr=17;
             }
           }

         }/*i<n*/

        if(ierr==0)
         {
          cntr == ' ';
          idigt = 0;
          ipntr = 0;
          jalloc= 0;
          idim  = 0;
          jtype = 0;
         }

       }/*flag*/

      irpt = 0;
      jrpt = 0;
      flag = 0;
     }

  if(ierr) jbrk=2;

 }/*jbrk*/


if(!ierr)
 {
  k = 0;

  if(ju<1 || j<1)
   {
    /*no format units in string */
    ierr=40;
   }
 
  else
   {
    for(ibrk=0;ibrk<ju;ibrk++)
     {
         jtype = Frmt->Ufrmt[ibrk].type;
      if(jtype==0)
       {
        ierr=44;
        break;
       }

      else
       {
           ic = idf_frmt_size(jtype);
        if(ic < 1)
         {
          /* improper format */
          ierr=44;
          break;
         }

        else
         {
          /* number of bytes to allocate */
          jmem = Frmt->Ufrmt[ibrk].jmem;
          jmem = jmem*ic;
          Frmt->Ufrmt[ibrk].jmem = jmem;

          /* number of bytes in format element */
          jsize = Frmt->Ufrmt[ibrk].jsize;
          jsize = jsize*ic;
          Frmt->Ufrmt[ibrk].jsize = jsize;

             ipntr = Frmt->Ufrmt[ibrk].pointer;
          if(ipntr)
           {
            jrpt = idf_frmt_psize(jtype);
            Frmt->Ufrmt[ibrk].pointer = jrpt;
           }
          else
           {
            /* maximal word size */
            jrpt = idf_frmt_word(jtype);
           }

          if(jrpt > k) k=jrpt;
         }
       }

     }/*ibrk*/

    if(ierr==0)
     {
      if(k==0)
       {
        /* improper format word size */
        ierr=45;
       }
     }
   }
 
 }

err:
    Frmt->Rfrmt[j] = IDF_EOS;

    Frmt->NfrmtU = ju;
    Frmt->NfrmtE = j;
    Frmt->swrd   = k;

    Frmt->jerr = ierr;
    Frmt->kpos = pos;

return ierr;
}



void idf_frmt_prn()
{
int n,m,idim;
register int i,j;

if(Frmt != NULL)
 {
     n = Frmt->NfrmtE;
  if(n)
   {
    m = Frmt->NfrmtU;

    fprintf(stdout,"INPUT FORMAT:\n");

    fprintf(stdout,"Number of elements=%d\n",n);
    fprintf(stdout,"Number of units=%d\n",m);
    fprintf(stdout,"Maximal word=%d\n",Frmt->swrd);

    fprintf(stdout,"Elements:%s\n",Frmt->Rfrmt);

    fprintf(stdout,"Repeatitions:\n");
    j=0;
    for(i=0;i<n;i++)
     {
      fprintf(stdout," %5d", Frmt->Mrpt[i]);
         j++;
      if(j==5)
       {
        fprintf(stdout,"\n");
        j=0;
       }
     }
    if(j)
     {
      fprintf(stdout,"\n");
     }

    if(m)
     {
      for(i=0;i<m;i++)
       {
        fprintf(stdout,"\nFormat unit=%d\n",i+1);

        fprintf(stdout,"Repeatition number=%d\n",
                       Frmt->Ufrmt[i].Nrep);

        if(Frmt->Ufrmt[i].pointer != 0)
         {
          fprintf(stdout,"Pointer Type=%d\n",
                       Frmt->Ufrmt[i].type);
         }
        else
         {
          fprintf(stdout,"Type=%d\n",
                       Frmt->Ufrmt[i].type);
         }

        fprintf(stdout,"Number of bytes=%d\n",
                       Frmt->Ufrmt[i].jsize);

        fprintf(stdout,"Allocation=%d Bytes of memory=%d\n",
                       Frmt->Ufrmt[i].jalloc,
                       Frmt->Ufrmt[i].jmem);

        idim = Frmt->Ufrmt[i].Ndim;
        fprintf(stdout,"Number of dimensions=%d\n",idim);
        if(idim)
         {
          for(j=0;j<idim;j++)
           {
            fprintf(stdout," %8d", Frmt->Ufrmt[i].Mdim[j]);
           }
          fprintf(stdout,"\n");
         }
       }
     }

   }/*n*/
 }

}


int idf_frmt_noalloc()
{
 int flag,k;
 int n;
 register int i;

flag = 0;

if(Frmt != NULL)
 {
     n = Frmt->NfrmtU;
  if(n)
   {
    for(i=0;i<n;i++)
     {
         k = Frmt->Ufrmt[i].jalloc;
      if(k)
       {
        flag=1;
        break;
       }
     }
   }
 }

return flag;
}


void idf_frmt_err_prn()
{
#define IDF_ERR_FRMT_N 47

static char *ErrFRMT[IDF_ERR_FRMT_N] = {

/* 1*/ "improper letter in format unit body",
/* 2*/ "comma in format unit body",
/* 3*/ "encomplete format unit at the end of string",
/* 4*/ "{ or } in format unit body",
/* 5*/ "improper symbol in format unit body",
/* 6*/ "too many digits in format unit body",
/* 7*/ "improper position of &, @, or # in format unit body",
/* 8*/ "large repeatition number in format string",
/* 9*/ "improper repeatition number for format unit",
/*10*/ "improper format symbol in a string",
/*11*/ "improper repeatition number for format symbol",
/*12*/ "zero number in format unit body",
/*13*/ "multiple pointer symbol in format body",
/*14*/ "improper symbol in format string",
/*15*/ "improper comma separator in format",
/*16*/ "zero repeatition number in format string",
/*17*/ "string contains too many format elements",
/*18*/ "improper } in format string",
/*19*/ "extra } in format string",
/*20*/ "return point was not found in format string",
/*21*/ "improperly placed [] or () in format unit body",
/*22*/ "improperly placed [] or () in format string",
/*23*/ "improper digits at the end of format string",
/*24*/ "unclosed { in format string",
/*25*/ "too many format units in a string",
/*26*/ "mixing [] with () in format string",
/*27*/ "improper sequence of [] or () in format unit dimensions",
/*28*/ "unclosed [ in format unit dimensions",
/*29*/ "improper position of [ in format unit dimensions",
/*30*/ "improper position of ( in format unit dimensions",
/*31*/ "empty [] in format unit",
/*32*/ "too many dimensions for format unit",
/*33*/ "extra ) in format unit dimensions",
/*34*/ "improper comma separator in format unit dimensions",
/*35*/ "improper comma in format string",
/*36*/ "too large number in format unit dimension",
/*37*/ "improper symbol in format unit dimension",
/*38*/ "encomplete dimension at the end of format string",
/*39*/ "encomplete format unit at the end of string",
/*40*/ "no format units in string",
/*41*/ "unclosed { period in format string",
/*42*/ "improper dimensions for format unit",
/*43*/ "improper dimension type for format unit",
/*44*/ "impoper format in structure",
/*45*/ "improper format word size",
/*46*/ "empty format string",
/*  */ "unknown error"
                                       };
int jerr;
register int j;
char *s;

if(Frmt != NULL)
{
   jerr = Frmt->jerr;
if(jerr>0)
 {
     j = jerr;
  if(j > IDF_ERR_FRMT_N) j = IDF_ERR_FRMT_N;
   {
     j--;
    s = ErrFRMT[j];

    fprintf(stdout,"format string error=%d position=%d:\n%s\n",
            jerr,Frmt->kpos,s);
   }
 }
}

}


