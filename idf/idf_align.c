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
int idf_align_bnd(int align,int swrd)
#else
int idf_align_bnd(align,swrd)
int align;
int swrd;
#endif
{
static int Abnd[] = {
     IDF_ALIGN_BYTE     , IDF_ALIGN_WORD     , IDF_ALIGN_LONGWORD ,
     IDF_ALIGN_QUADWORD , IDF_ALIGN_OCTAWORD , IDF_ALIGN_HEXWORD  ,
     IDF_ALIGN_PAGE     , IDF_ALIGN_BYTE     , IDF_ALIGN_WORD     ,
     IDF_ALIGN_LONGWORD , IDF_ALIGN_QUADWORD , IDF_ALIGN_OCTAWORD ,
     IDF_ALIGN_HEXWORD  , IDF_ALIGN_PAGE     , 0
                    };

register int j;
int wrd;

     if(align < 0             ) j = IDF_ALIGN_NMAX;
else if(align > IDF_ALIGN_NMAX) j = IDF_ALIGN_NMAX;
else                            j = align;

if(align == IDF_ALIGN_NATURAL) wrd = swrd;
else                           wrd = Abnd[j];

return wrd;
}


#if IDF_CPP_ == 1
void idf_align(int align, int swrd, int amode,
               int ksize, int wsize,
               int* Gap, int* Ofs)
#else
void idf_align(align,swrd,amode, ksize,wsize, Gap,Ofs)
int align,swrd,amode;
int ksize,wsize;
int *Gap,*Ofs;
#endif
{
int gap,ofs;
int l,lw;

gap = *Gap;
ofs = *Ofs;

  if(align == 0)
   {
    /* no gaps */
    ofs = ofs + gap;
    gap = 0;
   }

  else if(align < 8)
   {
     if(amode==IDF_AMODE_BINARI)
      {
       /* t->b  b->b */
       l = ksize;
       ofs = ofs + gap;

       if(wsize <= gap)
        {
         l = l - wsize;
         ofs = ofs - wsize;
        }

       lw = l/swrd;
       lw = lw*swrd;
       if(l > lw) lw = lw+swrd;
       gap = lw - l;
      }

     else
      {
       /* t->t  b->t */
       if(ksize <= gap)
        {
         gap = gap - ksize;
        }

       else
        {
          l = ksize - gap;
         lw = l/swrd;
         lw = lw*swrd;
         if(l > lw) lw = lw+swrd;
         gap = lw-l;
        }
      }
   }

  else
   {
    /* prescribed packing size */
    l = ksize;
    ofs = ofs + gap;

    lw = l/swrd;
    lw = lw*swrd;
    if(l > lw) lw = lw+swrd;
    gap = lw - l;
   }

*Ofs = ofs;
*Gap = gap;
}
