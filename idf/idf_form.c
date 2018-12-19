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



#define IDF_FORM_TRACE_IT 0

/* to view the formular, it is max number
of characters in the output strings */
#define IDF_FORM_VIEW_MAXL 76



static IDF_FORM *Dform=NULL;




IDF_FORM* idf_form_address()
{
 return ( Dform );
}



int idf_form_ini()
{
int ierr=0;
int k;

if(Dform != NULL)
 {
  ierr=1;
 }

else
 {
        k = sizeof(IDF_FORM);
     Dform = (IDF_FORM*) malloc(k);
  if(Dform == NULL)
   {
    ierr=2;
   }
 }

return ierr;
}


void idf_form_end()
{
if(Dform != NULL)
 {
  free(Dform);
  Dform = NULL;
 }
}


void idf_form_nul()
{
 Dform->kstr = 0;
 Dform->kurs = 0;
 Dform->kpos = 0;

 Dform->iBRC = -1;
 Dform->iCNS = -1;
 Dform->iELM = -1;

 Dform->jtype  = 0;
 Dform->Val[0] = (double)0.0;
 Dform->Val[1] = (double)0.0;
}



void idf_form_prn()
{
register int i,j;
int k;

 if(Dform != NULL)
 {
  fprintf(stdout,"formular data has type=%1d\n",
          Dform->jtype);

  fprintf(stdout,"it starts at position=%ld string=%d kursor=%d\n",
          Dform->kpos,Dform->kstr,Dform->kurs);

  if(Dform->jtype==0)
   {
    idf_form_view();
   }
  else
   {
    fprintf(stdout,"value: (%e , %e)\n",
            Dform->Val[0],Dform->Val[1]);
   }
 }
}


void idf_form_view()
{
static char frmt[] = "3d";

register int i,j;
int k,il,ic,n;
double vr,vi;
char c;
char *s;

 if(Dform != NULL)
  {
       k = Dform->iELM;
       k++;
    if(k>1)
     {
      il=0;
      for(i=0;i<k;i++)
       {
        c = Dform->ELMtyp[i];

        if(c=='{'||c=='}'||c==',')
         {
          fprintf(stdout,"%c",c);
          il++;
         }

        else if(c==')'||c==']'||c=='('||c=='[')
         {
          fprintf(stdout,"%c",c);
          il++;
         }

        else if(c=='%')
         {
          ic = Dform->ELMdat[i];

           s = idf_keywF_name(ic);
           n = idf_keywF_lname(ic);

          if(s!=NULL && n>0)
           {
               j = il+(n/2)+1;
            if(j > IDF_FORM_VIEW_MAXL)
             {
              fprintf(stdout,"\n");
              il=0;
             }

            for(j=0;j<n;j++)
             {
              fprintf(stdout,"%c",s[j]);
             }
            il += n;
           }

          else
           {
               j = il+9;
            if(j > IDF_FORM_VIEW_MAXL)
             {
              fprintf(stdout,"\n");
              il=0;
             }

            fprintf(stdout,"unknown[]");
            il += 9;
           }
         }

        else if(c=='#')
         {
             j = Dform->ELMdat[i];
          if(j>=0 && j<=Dform->iCNS)
           {
               ic = Dform->CNStyp[j];
               vr = Dform->CNSvr [j];
               vi = Dform->CNSvi [j];

            if(ic==IDF_FRMT_D)
             {
                 j = il+4;
              if(j > IDF_FORM_VIEW_MAXL)
               {
                fprintf(stdout,"\n");
                il=0;
               }

              fprintf(stdout,"dbl(");
              il += 4;

              s = idf_dtos(vr,frmt);
              n = strlen(s);

                 j = il+(2*n/3)+1;
              if(j > IDF_FORM_VIEW_MAXL)
               {
                fprintf(stdout,"\n");
                il=0;
               }

              fprintf(stdout,"%s)",s);
              il += n+1;
             }

            else if(ic==IDF_FRMT_D)
             {
                 j = il+6;
              if(j > IDF_FORM_VIEW_MAXL)
               {
                fprintf(stdout,"\n");
                il=0;
               }

              fprintf(stdout,"complx(");
              il += 7;

              s = idf_dtos(vr,frmt);
              n = strlen(s);

                 j = il+(2*n/3)+1;
              if(j > IDF_FORM_VIEW_MAXL)
               {
                fprintf(stdout,"\n");
                il=0;
               }

              fprintf(stdout,"%s,",s);
              il += n+1;

              s = idf_dtos(vi,frmt);
              n = strlen(s);

                 j = il+(2*n/3)+1;
              if(j > IDF_FORM_VIEW_MAXL)
               {
                fprintf(stdout,"\n");
                il=0;
               }

              fprintf(stdout,"%s)",s);
              il += n+1;
             }

            else
             {
                 j = il+9;
              if(j > IDF_FORM_VIEW_MAXL)
               {
                fprintf(stdout,"\n");
                il=0;
               }

              fprintf(stdout,"unknown()");
              il += 9;
             }
           }

          else
           {
            fprintf(stdout,"%c",c);
            il++;
           }
         }

        else
         {
          fprintf(stdout,"%c",c);
          il++;
         }

        if(il>IDF_FORM_VIEW_MAXL)
         {
          fprintf(stdout,"\n");
          il=0;
         }
       }

      if(il)
       {
        fprintf(stdout,"\n");
       }
     }
  }
}



#if IDF_CPP_ == 1
int idf_form_do(int* Kstr, int* Kurs, long* Kpos,
                unsigned int* Symb)
#else
int idf_form_do(Kstr,Kurs,Kpos, Symb)
int  *Kstr,*Kurs;
long *Kpos;
unsigned int *Symb;
#endif
{
int ierr=0;
char cntr,cc,cf;
char math,funa;
unsigned int c,ci,cerr;
int ibrk,jbrk,k,kk;
int flag,jerr,Jtype,sig;
int ielem,icns,hier,hhier,iopen,kend;
int iexpr,itok,iofl,jmath,jtype;
int kstr,kurs;
long kpos;
int typ,typ1,typ2;
double vr,vi,vr1,vi1,vr2,vi2;
int  *pKstr,*pKurs;
long *pKpos;
char   *Brace;
int    *CNStyp;
double *CNSvr,*CNSvi;
char   *ELMtyp;
int    *ELMdat;
IDF_FBUF *Fbuf;

Fbuf = idf_fbuf_address();
       idf_fbuf_nul();

kstr  = *Kstr;
kurs  = *Kurs;
kpos  = *Kpos;

 Dform->kstr = kstr;
 Dform->kurs = kurs;
 Dform->kpos = kpos;

 Dform->iBRC = -1;
 Dform->iCNS = -1;
 Dform->iELM = -1;

 Dform->jtype  = 0;
 Dform->Val[0] = (double)0.0;
 Dform->Val[1] = (double)0.0;

 vr=0.0; vi=0.0;
vr1=0.0;vi1=0.0;
vr2=0.0;vi2=0.0;
typ=0;typ1=0;typ2=0;

cerr=' ';

   ci = *Symb;
if(ci != '$')
 {
  /* improper call*/
  ierr=1;
  goto err;
 }

 pKstr = &kstr;
 pKurs = &kurs;
 pKpos = &kpos;

 Brace = Dform->Brace;

  icns = -1;
CNStyp = Dform->CNStyp;
CNSvr  = Dform->CNSvr;
CNSvi  = Dform->CNSvi;

ielem  = -1;
ELMtyp = Dform->ELMtyp;
ELMdat = Dform->ELMdat;

 cntr  = '{';

       iopen  = 0;
 Brace[iopen] = cntr;

       ielem++;
ELMtyp[ielem] = cntr;
ELMdat[ielem] = 0;

 flag  = 0;
hhier  = IDF_HIER_NOTHING;
 kend  = 0;
 cerr  = ' ';

      jbrk=0;
while(jbrk==0)
{

/*=================================================
  accumulate the sequence of formula elemements
      until it can be partially executed
  =================================================*/

 sig   = 0;
 hier  = hhier;
 jmath = 0;

      ibrk=0;
while(ibrk==0)
 {

 /*--------------------
   get the next symbol 
  ---------------------*/

  if(flag)
   {
    /* return stored symbol */
    c = ci;

    flag = 0;
    jerr = 0;
   }

  else
   {
    /* get next symbol */

      k = idf_file_getc(&c);
      kurs++;
      kpos++;

      if(k==2)
       {
        /* I/O error */
        jerr=2;
       }
      else if(k==1)
       {
	/* unexpected EOF */
        jerr=3;
       }
      else
       {
        jerr=0;
       }
   }

      if(jerr)
       {
        ierr=jerr;
       }


    /*-------------
       EOS case
     --------------*/

      else if( idf_string_over(c) )
       {
        kstr++;
        kurs=0;
       }


    /*-----------------
      end of formular
     ------------------*/

      else if(kend)
       {
        if(c == '$')
         {
          /* concatenation */
          kend=0;
         }

        else
         {
          /* truly formular end */
          ci = c;

               if(cntr=='+') k=1;
          else if(cntr=='-') k=1;
          else if(cntr=='*') k=1;
          else if(cntr=='/') k=1;
          else if(cntr=='^') k=1;
          else if(cntr=='(') k=1;
          else if(cntr=='{') k=1;
          else if(cntr=='%') k=3;
          else if(cntr=='&') k=2;
          else if(hier==IDF_HIER_NOTHING) k=1;
          else if(iopen!=0 ) k=1;
          else               k=0;

          if(k)
           {
            if(k==3)
             {
              ierr=31;
             }
            else if(k==2)
             {
              ierr=33;
             }
            else
             {
              /* parse error*/
              ierr=6;
             }
            cerr=c;
           }

          else
           {
            if(ielem < IDF_FORM_ELM_MAX_1)
             {
              cntr = '}';
                     ielem++;
              ELMtyp[ielem] = cntr;
              ELMdat[ielem] = 0;

              /* must be executed*/
              hhier = IDF_HIER_NUMBER;
              jmath=0;
              ibrk=1;
             }
            else
             {
              /* too many elements accumulated in formular*/
              ierr=4;
             }
           }
         }
       }


      else if(c == '$' )
       {
        /* concatenation is possible*/
        kend++;
       }


     /*--------------
        math symbols 
       --------------*/

      else if(c=='+' || c=='-')
       {
             if(cntr=='+') k=1;
        else if(cntr=='-') k=1;
        else if(cntr=='*') k=1;
        else if(cntr=='/') k=1;
        else if(cntr=='^') k=1;
        else if(cntr=='%') k=4;
        else if(cntr=='&') k=3;
        else if(cntr=='(') k=2;
        else if(cntr=='{') k=2;
        else if(cntr=='[') k=2;
        else if(cntr==',') k=2;
        else               k=0;

        if(k==2)
         {
          /* should be sign of a number,
            consider it as 0+ or 0-*/

          if(hier != IDF_HIER_NOTHING)
           {
            ierr=33;
           }
          else
           {
            cntr = '&';
            if(c=='-') sig=1; else sig=0;
           }
         }

        else if(k)
         {
          if(k==4)
           {
            ierr=31;
           }
          else if(k==3)
           {
            ierr=33;
           }
          else
           {
            /* parse error*/
            ierr=6;
           }
          cerr=c;
         }

        else
         {
          /* math symbol */
          if(ielem < IDF_FORM_ELM_MAX_1)
           {
            cntr = (char)c;
                   ielem++;
            ELMtyp[ielem] = cntr;
            ELMdat[ielem] = IDF_HIER_PLUS;

            if(hier >= IDF_HIER_PLUS)
             {
              /* can be executed*/
              hhier = IDF_HIER_PLUS;
              jmath=1;
              ibrk=1;
             }
            else
             {
              hier = IDF_HIER_PLUS;
             }
           }

          else
           {
            /* too many elements accumulated in formular*/
            ierr=4;
           }
         }
       }


      else if(c=='*')
       {
             if(cntr=='+') k=1;
        else if(cntr=='-') k=1;
        else if(cntr=='*') k=1;
        else if(cntr=='/') k=1;
        else if(cntr=='^') k=1;
        else if(cntr=='%') k=3;
        else if(cntr=='&') k=2;
        else if(cntr=='(') k=1;
        else if(cntr=='[') k=1;
        else if(cntr=='{') k=1;
        else if(cntr==',') k=1;
        else               k=0;

        if(k)
         {
          if(k==2)
           {
            ierr=33;
           }
          else if(k==3)
           {
            ierr=31;
           }
          else
           {
            /* parse error*/
            ierr=6;
           }
          cerr=c;
         }

        else
         {
          /* math symbol */
          kk = idf_file_getc(&ci);
          kurs++;
          kpos++;

          if(kk==2)
           {
            /* I/O error */
            ierr=2;
           }
          else if(kk==1)
           {
	    /* unexpected EOF */
            ierr=3;
           }
          else
           {
            if(ci=='*')
             {
              /* sequence ** is equivalent to ^ */
              cntr = '^';
              kk = IDF_HIER_POW;
             }
            else
             {
              cntr = '*';
              kk = IDF_HIER_MULT;
              flag=1;
             }

            if(ielem < IDF_FORM_ELM_MAX_1)
             {
                     ielem++;
              ELMtyp[ielem] = cntr;
              ELMdat[ielem] = kk;

              if(hier >= kk)
               {
                /* can be executed*/
                hhier = kk;
                jmath=1;
                ibrk=1;
               }
              else
               {
                hier = kk;
               }
             }

            else
             {
              /* too many elements accumulated in formular*/
              ierr=4;
             }
           }
         }
       }


      else if(c=='^')
       {
             if(cntr=='+') k=1;
        else if(cntr=='-') k=1;
        else if(cntr=='*') k=1;
        else if(cntr=='/') k=1;
        else if(cntr=='^') k=1;
        else if(cntr=='%') k=3;
        else if(cntr=='&') k=2;
        else if(cntr=='(') k=1;
        else if(cntr=='[') k=1;
        else if(cntr=='{') k=1;
        else if(cntr==',') k=1;
        else               k=0;

        if(k)
         {
          if(k==2)
           {
            ierr=33;
           }
          else if(k==3)
           {
            ierr=31;
           }
          else
           {
            /* parse error*/
            ierr=6;
           }
          cerr=c;
         }

        else
         {
          /* math symbol */
          if(ielem < IDF_FORM_ELM_MAX_1)
           {
            cntr = '^';
                   ielem++;
            ELMtyp[ielem] = cntr;
            ELMdat[ielem] = IDF_HIER_POW;

              if(hier >= IDF_HIER_POW)
               {
                /* can be executed*/
                hhier = IDF_HIER_POW;
                jmath=1;
                ibrk=1;
               }
              else
               {
                hier = IDF_HIER_POW;
               }
           }
          else
           {
            /* too many elements accumulated in formular*/
            ierr=4;
           }
         }
       }


    /*--------------------
      comment or division
     ---------------------*/

      else if(c == '/')
       {
        /* skip comment if appropriate */

          ci = c;

          k = idf_stxt_skip_cmnt(pKstr,pKurs,pKpos, &ci);

          if(k == 2)
           {
            /* division symbol */
                 if(cntr=='+') k=1;
            else if(cntr=='-') k=1;
            else if(cntr=='*') k=1;
            else if(cntr=='/') k=1;
            else if(cntr=='^') k=1;
            else if(cntr=='%') k=3;
            else if(cntr=='&') k=2;
            else if(cntr=='(') k=1;
            else if(cntr=='[') k=1;
            else if(cntr=='{') k=1;
            else if(cntr==',') k=1;
            else               k=0;

            if(k)
             {
              if(k==2)
               {
                ierr=33;
               }
              else if(k==3)
               {
                ierr=31;
               }
              else
               {
                /* parse error*/
                ierr=6;
               }
              cerr=c;
             }

            else
             {
              /* math symbol */
              if(ielem < IDF_FORM_ELM_MAX_1)
               {
                cntr = '/';
                       ielem++;
                ELMtyp[ielem] = cntr;
                ELMdat[ielem] = IDF_HIER_MULT;

                if(hier >= IDF_HIER_MULT)
                 {
                  /* can be executed*/
                  hhier = IDF_HIER_MULT;
                  jmath=1;
                  ibrk=1;
                 }
                else
                 {
                  hier = IDF_HIER_MULT;
                 }

                flag = 1;
               }

              else
               {
                /* too many elements accumulated in formular*/
                ierr=4;
               }
             }
           }
          else if(k == 1)
           {
            /* unexpected EOF */
            ierr = 3;
           }
          else if(k < 0)
           {
            /* comment error */
                 if(k==-1) ierr=2;
            else if(k==-2) ierr=3;
            else           ierr=7;
            cerr=c;
           }
          else
           {
            /* store symbol */
            flag=1;
           }

       } /* '/' */


    /*---------------
        name case 
     ----------------*/

      else if(
              ( isalpha(c) )
#if IDF_FOREIGN_LANGUAGE == 1
              || ( !isascii(c) )
#endif
             )
       {
        if(cntr=='%')
         {
          ierr=31;
          cerr=c;
         }
        else if(cntr==')' || cntr==']')
         {
          /* parse error*/
          ierr=6;
          cerr=c;
         }
        else
         {
             ci = c;
          Jtype = IDF_FBUF_TYPE_NAME;

          idf_fbuf_begin(pKstr,pKurs,pKpos, Jtype);

             k = idf_fbuf_name(pKstr,pKurs,pKpos, &ci);
          if(k)
           {
            /* name search error: 17-20 */
            ierr=7+k;
           }
          else
           {
            /* name is placed inside buffer */
               kk = idf_fbuf_conv();
            if(kk)
             {
              /*conversion error: 21-30*/
              ierr=7+kk;
             }
            else
             {
                 k = Fbuf->jconv;
              if(k==1)
               {
                /* name was converted into number*/
                if(icns < IDF_FORM_CNS_MAX)
                 {
                  vr = Fbuf->val[0];
                  vi = Fbuf->val[1];

                  if(cntr=='&')
                   {
                    /* hidden sign*/
                    if(sig)
                     {
                      vr = -vr;
                      vi = -vi;
                      sig = 0;
                     }
                   }

                         icns++;
                  CNStyp[icns] = Fbuf->jf;
                  CNSvr [icns] = vr;
                  CNSvi [icns] = vi;

                  if(ielem < IDF_FORM_ELM_MAX_1)
                   {
                    cntr = '#';
                           ielem++;
                    ELMtyp[ielem] = cntr;
                    ELMdat[ielem] = icns;

                    if(hier < IDF_HIER_NUMBER)
                     {
                      hier = IDF_HIER_NUMBER;
                     }

                    flag = 1;
                   }
                  else
                   {
                    /* too many elements accumulated in formular*/
                    ierr=4;
                   }
                 }
                else
                 {
                  /* too many constants accumulated in formular*/
                  ierr=5;
                 }
               }

              else
               {
                /* name was converted into function*/

                if(cntr=='&')
                 {
                  /* hidden sign: 0+fun, 0-fun*/
                  if(icns < IDF_FORM_CNS_MAX)
                   {
                           icns++;
                    CNStyp[icns] = IDF_FRMT_D;
                    CNSvr [icns] = 0.0;
                    CNSvi [icns] = 0.0;

                    if(ielem < IDF_FORM_ELM_MAX)
                     {
                             ielem++;
                      ELMtyp[ielem] = '#';
                      ELMdat[ielem] = icns;

                      if(sig) cntr = '-';
                      else    cntr = '+';

                      hier = IDF_HIER_PLUS;
                      sig  = 0;

                      /* math symbol */
                             ielem++;
                      ELMtyp[ielem] = cntr;
                      ELMdat[ielem] = hier;
                     }
                    else
                     {
                      /* too many elements accumulated in formular*/
                      ierr=4;
                     }
                   }
                  else
                   {
                    /* too many constants accumulated in formular*/
                    ierr=5;
                   }
                 }

                if(ierr)
                 {
                  ;
                 }
                else if(ielem < IDF_FORM_ELM_MAX_1)
                 {
                  cntr = '%';
                         ielem++;
                  ELMtyp[ielem] = cntr;
                  ELMdat[ielem] = Fbuf->jf;

                  flag=1;
                 }
                else
                 {
                  /* too many elements accumulated in formular*/
                  ierr=4;
                 }
               }
             }
           }
         }
       } /* name */


    /*---------------
       digital case 
     ----------------*/

      else if( isdigit(c) || c=='.')
       {
             if(cntr==']') k=1;
        else if(cntr==')') k=1;
        else if(cntr=='%') k=2;
        else if(cntr=='#') k=1;
        else               k=0;
        if(k)
         {
          if(k==2)
           {
            ierr=31;
           }
          else
           {
            /* parse error*/
            ierr=6;
           }
          cerr=c;
         }
        else
         {
             ci = c;
          Jtype = IDF_FBUF_TYPE_NUMBER;

          idf_fbuf_begin(pKstr,pKurs,pKpos, Jtype);

             k = idf_fbuf_number(pKstr,pKurs,pKpos, &ci);
          if(k)
           {
            /* number search error: 8-16 */
            ierr=7+k;
           }
          else
           {
            /* number is placed inside buffer */
               kk = idf_fbuf_conv();
            if(kk)
             {
              /*conversion error: 21-30*/
              ierr=7+kk;
             }
            else
             {
              if(icns < IDF_FORM_CNS_MAX)
               {
                vr = Fbuf->val[0];
                vi = Fbuf->val[1];

                if(cntr=='&')
                 {
                  /* hidden sign*/
                  if(sig)
                   {
                    vr = -vr;
                    vi = -vi;
                    sig = 0;
                   }
                 }

                       icns++;
                CNStyp[icns] = Fbuf->jf;
                CNSvr [icns] = vr;
                CNSvi [icns] = vi;

                if(ielem < IDF_FORM_ELM_MAX_1)
                 {
                  cntr = '#';
                         ielem++;
                  ELMtyp[ielem] = cntr;
                  ELMdat[ielem] = icns;

                  if(hier < IDF_HIER_NUMBER)
                   {
                    hier = IDF_HIER_NUMBER;
                   }

                  flag = 1;
                 }
                else
                 {
                  /* too many elements accumulated in formular*/
                  ierr=4;
                 }
               }
              else
               {
                /* too many constants accumulated in formular*/
                ierr=5;
               }
             }
           }
         }

       } /* digit */


    /*-------------
       brackets
     --------------*/

      else if(c == '(' )
       {
             if(cntr==')') k=1;
        else if(cntr==']') k=1;
        else if(cntr=='#') k=1;
        else               k=0;

        if(k)
         {
          /* parse error*/
          ierr=6;
          cerr=c;
         }

        else
         {
          if(cntr=='&')
           {
            /* hidden sign: 0+(, 0-(*/
            if(icns < IDF_FORM_CNS_MAX)
             {
                     icns++;
              CNStyp[icns] = IDF_FRMT_D;
              CNSvr [icns] = 0.0;
              CNSvi [icns] = 0.0;

              if(ielem < IDF_FORM_ELM_MAX)
               {
                       ielem++;
                ELMtyp[ielem] = '#';
                ELMdat[ielem] = icns;

                if(sig) cntr = '-';
                else    cntr = '+';

                hier = IDF_HIER_PLUS;
                sig  = 0;

                /* math symbol */
                       ielem++;
                ELMtyp[ielem] = cntr;
                ELMdat[ielem] = hier;

                cc='(';
               }
              else
               {
                /* too many elements accumulated in formular*/
                ierr=4;
               }
             }
            else
             {
              /* too many constants accumulated in formular*/
              ierr=5;
             }
           }

          else if(cntr=='%')
           {
            cc='[';
           }

          else
           {
            cc='(';
           }

          if(ierr)
           {
            ;
           }
          else if(iopen < IDF_FORM_BRACE_MAX)
           {
                  iopen++;
            Brace[iopen] = cc;

            if(ielem < IDF_FORM_ELM_MAX_1)
             {
                     ielem++;
              ELMtyp[ielem] = cc;
              ELMdat[ielem] = iopen;

              hier = IDF_HIER_NOTHING;
              cntr = cc;
             }
            else
             {
              /* too many elements accumulated in formular*/
              ierr=4;
             }
           }

          else
           {
            /* too many accumulated braces */
            ierr=32;
            cerr=cc;
           }
         }
       }


      else if(c == ')' )
       {
             if(cntr=='+') k=1;
        else if(cntr=='-') k=1;
        else if(cntr=='*') k=1;
        else if(cntr=='/') k=1;
        else if(cntr=='^') k=1;
        else if(cntr=='%') k=3;
        else if(cntr=='&') k=2;
        else if(cntr=='[') k=1;
        else if(cntr=='(') k=1;
        else if(cntr=='{') k=1;
        else if(cntr==',') k=1;
        else if(hier==IDF_HIER_NOTHING) k=1;
        else if(iopen==0 ) k=1;
        else               k=0;

        if(k)
         {
          if(k==2)
           {
            ierr=33;
           }
          else if(k==3)
           {
            ierr=31;
           }
          else
           {
            /* parse error*/
            ierr=6;
           }
          cerr=c;
         }

        else
         {
          if(ielem < IDF_FORM_ELM_MAX_1)
           {
                    cc = Brace[iopen];
                 if(cc=='[') cntr=']';
            else if(cc=='(') cntr=')';
            else             cntr=' ';

                   ielem++;
            ELMtyp[ielem] = cntr;
            ELMdat[ielem] = iopen;

            /* must be executed*/
            hhier = IDF_HIER_NUMBER;
            jmath=0;
            ibrk=1;
           }
          else
           {
            /* too many elements accumulated in formular*/
            ierr=4;
           }
         }
       }


    /*-----------------
       comma separator
     ------------------*/

      else if(c == ',' )
       {
             if(cntr=='+') k=1;
        else if(cntr=='-') k=1;
        else if(cntr=='*') k=1;
        else if(cntr=='/') k=1;
        else if(cntr=='^') k=1;
        else if(cntr=='%') k=3;
        else if(cntr=='&') k=2;
        else if(cntr=='(') k=1;
        else if(cntr=='[') k=1;
        else if(cntr=='{') k=1;
        else if(cntr==',') k=1;
        else if(hier==IDF_HIER_NOTHING) k=1;
        else               k=0;

        if(k)
         {
          if(k==2)
           {
            ierr=33;
           }
          else if(k==3)
           {
            ierr=31;
           }
          else
           {
            /* parse error*/
            ierr=6;
           }
          cerr=c;
         }

        else
         {
          if(ielem < IDF_FORM_ELM_MAX_1)
           {
            cntr = ',';
                   ielem++;
            ELMtyp[ielem] = cntr;
            ELMdat[ielem] = 0;

            /* must be executed*/
            if(hier > IDF_HIER_NUMBER)
             {
              hhier = IDF_HIER_NUMBER;
              jmath=0;
              ibrk=1;
             }
           }
          else
           {
            /* too many elements accumulated in formular*/
            ierr=4;
           }
         }
       }


    /*-----------------
      white space case
     ------------------*/

      else if(c == ' ')
       {
        ;
       }

      else if(c == '\t')
       {
#if IDF_TXT_TABULATORS == 1
        kurs+=7;
#else
        ;
#endif
       }

      else if(c == '\v')
       {
#if IDF_TXT_VERTICES == 1
        kstr+=8;
#else
        ;
#endif
       }


    /*-----------------------------------------
      all other symbols should not appear here
     ------------------------------------------ */

      else
       {
        ierr=46;
        cerr=c;
       }


  if(ierr)
   {
    ci = '$';
    ibrk=2;
   }


#if IDF_FORM_TRACE_IT == 1
fprintf(stdout,"c=%c hier=%d ie=% ic=%d io=%d ib=%d cntr=%c ierr=%d\n",
        c,hier,ielem,icns,iopen,ibrk,cntr,ierr);
#endif

 }/*ibrk*/


    Dform->iBRC = iopen;
    Dform->iELM = ielem;
    Dform->iCNS = icns;


 if(ierr) break;


/*=================================================
  if no errors, execute the accumulated expression
  =================================================*/


iexpr = 0;
itok  = 0;
funa  = ' ';
math  = ' ';
hier  = IDF_HIER_NOTHING;


      ibrk=0;
while(ibrk==0)
 {
  ielem--;

  cf = ELMtyp[ielem];

    /*------------------
       numerical target
      ------------------*/

    if(cf == '#')
     {
      if(hier==IDF_HIER_NOTHING)
       {
        /* first number */
        icns = ELMdat[ielem];

        typ2 = CNStyp[icns];
         vr2 = CNSvr [icns];
         vi2 = CNSvi [icns];

        itok++;
        hier = IDF_HIER_NUMBER;
       }

      else if(hier==IDF_HIER_NUMBER)
       {
        /*second number*/
        if(funa != ',')
         {
          /* case:##, improper binary token */
          ierr=35;
          cerr=cf;
         }
        else if(itok>1)
         {
          /* more then 2 arguments */
          ierr=34;
          cerr=cf;
         }
        else
         {
          /* next argument */
          icns = ELMdat[ielem];

          typ1 = CNStyp[icns];
           vr1 = CNSvr [icns];
           vi1 = CNSvi [icns];

          itok++;
          hier = IDF_HIER_NUMBER;
         }
       }

      else
       {
        icns = ELMdat[ielem];

        if(jmath)
         {
          if(hier >= hhier)
           {
            k=1;
           }
          else if(iexpr)
           {
            /* can not be executed since could be 'a+b*' */
                   ielem++;
            ELMtyp[ielem] = math;
            ELMdat[ielem] = hier;

            icns++;
                   ielem++;
            ELMtyp[ielem] = '#';
            ELMdat[ielem] = icns;

            CNStyp[icns] = typ2;
            CNSvr [icns] = vr2;
            CNSvi [icns] = vi2;

                   ielem++;
            ELMtyp[ielem] = cntr;
            ELMdat[ielem] = hhier;

            k=0;
            ibrk=1;
           }
          else
           {
            /* should be some portion of math*/
            ierr=35;
            cerr=cntr;
            k=0;
           }
         }
        else
         {
          k=1;
         }

        if(k)
         {
          typ1 = CNStyp[icns];
           vr1 = CNSvr [icns];
           vi1 = CNSvi [icns];

          /* execution of binary token*/
          iofl = idf_math(math,
                          typ1,vr1,vi1, typ2,vr2,vi2,
                          &typ,&vr,&vi);

#if IDF_FORM_TRACE_IT == 1
fprintf(stdout,"math: v1=%d %e %e v2=%d %e %e rez=%d %e %e\n",
typ1,vr1,vi1, typ2,vr2,vi2, typ,vr,vi);
#endif
          if(iofl)
           {
            /*  math error*/
            ierr=36;
            cerr=math;
           }
          else
           {
            typ2 = typ;
             vr2 = vr;
             vi2 = vi;

            itok = 1;
            hier = IDF_HIER_NUMBER;
            math = ' ';
            iexpr++;
           }
         }
       }
     }


    /*-------------
       arithmetic
      -------------*/

    else if(cf=='+' || cf=='-')
     {
      if(hier != IDF_HIER_NUMBER)
       {
        /* no number*/
        ierr=37;
        cerr=cf;
       }
      else
       {
        math = cf;
        hier = IDF_HIER_PLUS;
       }
     }

    else if(cf=='*' || cf=='/')
     {
      if(hier != IDF_HIER_NUMBER)
       {
        /* no number*/
        ierr=37;
        cerr=cf;
       }
      else
       {
        math = cf;
        hier = IDF_HIER_MULT;
       }
     }

    else if(cf=='^')
     {
      if(hier != IDF_HIER_NUMBER)
       {
        /* no number*/
        ierr=37;
        cerr=cf;
       }
      else
       {
        math = cf;
        hier = IDF_HIER_POW;
       }
     }


    /*-------------
        bracket
      -------------*/

    else if(cf == '(')
     {
      if(hier != IDF_HIER_NUMBER) k=1;
      else if(iopen < 1  )        k=1;
      else if(cntr != ')')        k=1;
      else                        k=0;
      if(k)
       {
        ierr=3;
        cerr=cf;
       }
      else
       {
           kk = ELMdat[ielem];
        if(kk != iopen)
         {
          /* extra bracket */
          ierr=38;
          cerr=cf;
         }
        else
         {
                         ielem--;
             cc = ELMtyp[ielem];
          if(cc=='[' || cc=='{')
           {
            /* {(#) or [(#)  */
            hier = IDF_HIER_NUMBER;
            kk=1;
           }
          else if(cc=='+'||cc=='-'||cc=='*'||cc=='/'||cc=='^')
           {
            /*  +*-/^(#) */
            hier = ELMdat[ielem];
            kk=1;
           }
          else
           {
            /* 2nd run parse error*/
            ierr=39;
            cerr=cf;
            kk=0;
           }

          if(kk)
           {
            cntr = '#';
                   ielem++;
            ELMtyp[ielem] = cntr;
            ELMdat[ielem] = icns;

            CNStyp[icns] = typ2;
            CNSvr [icns] = vr2;
            CNSvi [icns] = vi2;

            iopen--;
            hhier = hier;
            math=' ';
            ibrk=1;
           }
         }
       }

     }/* ( */


    /*-------------
        function
      -------------*/

    else if(cf == '[')
     {
      if(hier != IDF_HIER_NUMBER) k=1;
      else if(iopen < 1  )        k=1;
      else if(cntr != ']')        k=1;
      else if(itok == 0  )        k=1;
      else                        k=0;
      if(k)
       {
        ierr=3;
        cerr=cf;
       }
      else
       {
           kk = ELMdat[ielem];
        if(kk != iopen)
         {
          ierr=38;
          cerr=cf;
         }
        else
         {
                         ielem--;
             cc = ELMtyp[ielem];
          if(cc!='%')
           {
            /* parse error */
            ierr=39;
            cerr=cf;
           }
          else
           {
            jtype = ELMdat[ielem];

               iofl = idf_math_func(jtype,itok,
                                    typ1,vr1,vi1, typ2,vr2,vi2,
                                    &typ,&vr,&vi);
#if IDF_FORM_TRACE_IT == 1
fprintf(stdout,"func: v1=%d %e %e v2=%d %e %e rez=%d %e %e\n",
typ1,vr1,vi1, typ2,vr2,vi2, typ,vr,vi);
#endif
            if(iofl)
             {
              /* func math error*/
              ierr=40;
              cerr=cf;
             }
            else
             {
              typ2 = typ;
               vr2 = vr;
               vi2 = vi;

                             ielem--;
                 cc = ELMtyp[ielem];
              if(cc=='['||cc=='('||cc=='{'||cc==',')
               {
                /*  (%[#] [%[#]*/
                kk=1;
                hier = IDF_HIER_NUMBER;
               }
              else if(cc=='+'||cc=='-'||cc=='*'||cc=='/'||cc=='^')
               {
                /*  +*-/^ %[#]  */
                hier = ELMdat[ielem];
                kk=1;
               }
              else
               {
                kk=0;
                ierr=39;
                cerr=cc;
               }

              if(kk)
               {
                cntr = '#';
                       ielem++;
                ELMtyp[ielem] = cntr;
                ELMdat[ielem] = icns;

                CNStyp[icns] = typ2;
                CNSvr [icns] = vr2;
                CNSvi [icns] = vi2;

                iopen--;
                itok=1;
                hhier=hier;
                math=' ';
                funa=' ';
                ibrk=1;                
               }
             }
           }
         }
       }
     } /* [ */


    else if(cf == ',')
     {
      if(cntr==']')
       {
        if(hier != IDF_HIER_NUMBER)
         {
          /* argument is missing */
          ierr=41;
          cerr=cntr;
         }
        else
         {
          if(funa==',')
           {
            ierr=42;
            cerr=cntr;
           }
          else
           {
            funa = ',';
           }
         }
       }

      else if(cntr=='}'||cntr==','||cntr==')')
       {
        /* improper logic*/
        ierr=43;
        cerr=cntr;
       }

      else
       {
        if(hier != IDF_HIER_NUMBER)
         {
          /* improper expression */
          ierr=44;
          cerr=cf;
         }
        else
         {
          /* expression finished */
                 ielem++;
          ELMtyp[ielem] = '#';
          ELMdat[ielem] = icns;

          CNStyp[icns] = typ2;
          CNSvr [icns] = vr2;
          CNSvi [icns] = vi2;

                 ielem++;
          ELMtyp[ielem] = cntr;
          ELMdat[ielem] = hhier;

          ibrk=1;
         }
       }

     }/* , */


    /*----------------
       starting point
      ----------------*/

    else if(cf == '{')
     {
      if(hier == IDF_HIER_NUMBER)
       {
        if(cntr==']'||cntr==')'||cntr==',')
         {
          ierr=43;
          cerr=cf;
         }

        else if(iopen)
         {
          /* unclosed brace */
          ierr=45;
          cerr=Brace[iopen];
         }

        else if(cntr=='}')
         {
          /*formular finished*/
          icns=0;
                 ielem++;
          ELMtyp[ielem] = '#';
          ELMdat[ielem] = icns;

          CNStyp[icns] = typ2;
          CNSvr [icns] = vr2;
          CNSvi [icns] = vi2;

                 ielem++;
          ELMtyp[ielem] = cntr;
          ELMdat[ielem] = 0;

          ibrk=2;
          Dform->jtype  = typ2;
          Dform->Val[0] = vr2;
          Dform->Val[1] = vi2;

          /* end of formular flag*/
          jbrk=2;
         }

        else
         {
          if(iexpr==0)
           {
            ierr=44;
            cerr=cf;
           }
          else
           {
            /*expression done*/
            icns=0;
                   ielem++;
            ELMtyp[ielem] = '#';
            ELMdat[ielem] = icns;

            CNStyp[icns] = typ2;
            CNSvr [icns] = vr2;
            CNSvi [icns] = vi2;

                   ielem++;
            ELMtyp[ielem] = cntr;
            ELMdat[ielem] = hhier;

            ibrk=1;
           }
         }
       }

      else
       {
        /*alg error */
        ierr=43;
        cerr=cf;
       }
     }


    else if(cf==')' || cf==']')
     {
      /* improper formular decoding alg*/
      ierr=43;
      cerr=cf;
     }

    else
     {
      /* improper symbol*/
      ierr=44;
      cerr=cf;
     }

  if(ierr) ibrk=1;

#if IDF_FORM_TRACE_IT == 1
fprintf(stdout,"cf=%c hier=%d ie=% ic=%d io=%d ib=%d cntr=%c ierr=%d\n",
        cf,hier,ielem,icns,iopen,ibrk,cntr,ierr);
#endif

 }/*ibrk*/


 if(ierr)
  {
   jbrk=1;
   ci = '$';
  }

 else
  {
   Dform->iBRC = iopen;
   Dform->iELM = ielem;
   Dform->iCNS = icns;
  }

}/*jbrk*/


err:
*Kstr = kstr;
*Kurs = kurs;
*Kpos = kpos;
*Symb = ci;

if(ierr)
 {
  k = ierr + 157;
  idf_txt_err_put(k, kstr,kurs,kpos, cerr);
 }

return ierr;
}



#if IDF_CPP_ == 1
int idf_form_up(int jfrmt, void* val)
#else
int idf_form_up(jfrmt,val)
int   jfrmt;
void *val;
#endif
{
/*
  0 - OK
  1 - improper format type for formular data
  2 - overflow in conversion formular into int
  3 - overflow in conversion formular into long
  4 - overflow in conversion formular into float
  5 - overflow in conversion formular into complex float
  6 - no formular data stored
*/
int flag,jtype;
int *i;
long *l;
float *f;
double *d;
double g,gg;

 flag = 0;
jtype = Dform->jtype;
 
if(jtype < 1)
 {
  /* no formular data stored */
  flag=6;
 }

else if(jfrmt==4)
 {
  i = (int*)val;
  g = Dform->Val[0];
         gg = (double)IDF_INT_MAX;
  if(g > gg)
   {
    flag = 2;
    *i = (int)IDF_INT_MAX;
   }
  else
   {
    *i = (int)g;
   }
 }

else if(jfrmt==5)
 {
  l = (long*)val;
  g = Dform->Val[0];
         gg = (double)IDF_LONG_MAX;
  if(g > gg)
   {
    flag = 3;
    *l = (long)IDF_LONG_MAX;
   }
  else
   {
    *l = (long)g;
   }
 }

else if(jfrmt==6)
 {
  f = (float*)val;
  g = Dform->Val[0];
         gg = (double)IDF_FLOAT_MAX;
  if(g > gg)
   {
    flag = 4;
    *f = (float)IDF_FLOAT_MAX;
   }
  else
   {
    *f = (float)g;
   }
 }

else if(jfrmt==7)
 {
  d = (double*)val;
  g = Dform->Val[0];
  *d = g;
 }

else if(jfrmt==8)
 {
  f = (float*)val;
  g = Dform->Val[0];
         gg = (double)IDF_FLOAT_MAX;
  if(g > gg)
   {
    flag = 5;
    f[0] = (float)IDF_FLOAT_MAX;
   }
  else
   {
    f[0] = (float)g;

       g = Dform->Val[1];
    if(g > gg)
     {
      flag = 5;
      f[1] = (float)IDF_FLOAT_MAX;
     }
    else
     {
      f[1] = (float)g;
     }
   }
 }

else if(jfrmt==9)
 {
  d = (double*)val;
  g = Dform->Val[0];
  d[0] = g;
  g = Dform->Val[1];
  d[1] = g;
 }

else
 {
  /* improper format type for formular data*/
  flag=1;
 }

return flag;
}
