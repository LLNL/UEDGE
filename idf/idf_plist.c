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


/*maximal number of standard physical constants */
#define IDF_PLIST_MAX 30

/*
 --------------------------- 
   physical constants (CGS)
 ---------------------------
     
IDF_C      c         Speed of light              2.99792458e10    cm/sec
IDF_H      h         Planck constant             6.6260196e-27    erg sec
IDF_HB     hbar      Planck constant/2pi         1.054592e-27     erg sec
IDF_K      k         Boltzmann constant          1.380658e-16     erg/ K
IDF_K_EV   k         Boltzmann constant in eV    8.617385e-5      k/|e| eV/K
IDF_E      e         Elementary charge           4.806532e-10     statcoul
IDF_ME     m_e       Electron mass               9.1093897e-28    g
IDF_MP     m_p       Proton mass                 1.6726231e-24    g
IDF_MPME   m_p/m_e                               1836.152755
IDF_MEMP   m_e/m_p                               5.446169971e-4
IDF_G      g         Gravitational constant      6.6732e-8        dyne cm^2/g^2
IDF_RY     Ry        Rydberg constant            109737.31534     cm^(-1)
IDF_RY_EV  Ry        Rydberg constant            13.6056981       eV
IDF_RB     rB        Bohr radius                 0.529177249e-8   cm
IDF_CS     pi rB^2   Atomic cross section        0.87973567e-16   cm^2
IDF_RE     rE        Classical electron radius   2.8179e-13       cm
IDF_ALP    alpha     Fine-structure constant     7.297351e-3
IDF_IALP  1/alpha                                137.036
IDF_MUB    muB       Bohr magneton               5.78838263e-5    eV/Tesla
IDF_CW   hbar/m_e/c  Compton electron wavelength 2.4263e-10       cm
IDF_FA     F         Faraday constant            2.892599e14      statcoul/mol
IDF_CR1   8pi hc     First radiation costant     4.992579e-15     erg/mol
IDF_CR2   hc/k       Second radiation constant   1.43883          cm K
IDF_S     sigma      Stefan-Boltzmann constant   5.66961e-5       erg/cm2/s/K^4
IDF_NA     NA        Avogadro number             6.022169e23      1/mol
IDF_R     R=kNA      Gas constant                8.31434e7        erg/deg/mol
IDF_TO     To        Standard temperature        273.15           K
IDF_PO   Po=LokTo    Atmospheric pressure        1.0133e6         dyne/cm2
IDF_LO     Lo        Loschmidt's number          2.6868e19        cm-3
IDF_VO   Vo=RTo/Po   Normal volume perfect gas   2.24136e4        cm3/mol

*/


#if IDF_CPP_ == 1
int idf_plist(char* name, int nn, double* val)
#else
int idf_plist(name,nn, val)
char   *name;
int     nn;
double *val;
#endif
{

static char  *Clist[IDF_PLIST_MAX] = {
            "C",      "H",     "HB",
            "K",   "K_EV",      "E",
           "ME",     "MP",   "MPME",
         "MEMP",      "G",     "RY",
        "RY_EV",     "RB",     "CS",
           "RE",    "ALP",   "IALP",
          "MUB",     "CW",     "FA",
          "CR1",    "CR2",      "S",
           "NA",      "R",     "TO",
           "PO",     "LO",     "VO"
                                     };

static int   NClist[IDF_PLIST_MAX] = {
                  1,    1,    2,
                  1,    4,    1,
                  2,    2,    4,
                  4,    1,    2,
                  5,    2,    2,
                  2,    3,    4,
                  3,    2,    2,
                  3,    3,    1,
                  2,    1,    2,
                  2,    2,    2
                                     };


static double Vlist[IDF_PLIST_MAX] = {
   2.99792458e10 ,  6.6260196e-27 ,  1.054592e-27  ,
   1.380658e-16  ,  8.617385e-5   ,  4.806532e-10  ,
   9.1093897e-28 ,  1.6726231e-24 ,  1836.152755   ,
   5.446169971e-4,  6.6732e-8     ,  109737.31534  ,
   13.6056981    ,  0.529177249e-8,  0.87973567e-16,
   2.8179e-13    ,  7.297351e-3   ,  137.036       ,
   5.78838263e-5 ,  2.4263e-10    ,  2.892599e14   ,
   4.992579e-15  ,  1.43883       ,  5.66961e-5    ,
   6.022169e23   ,  8.31434e7     ,  273.15        ,
   1.0133e6      ,  2.6868e19     ,  2.24136e4
                                     };

int k,n;
register int i,j;
char *s;


  i=0;
while(i<IDF_PLIST_MAX)
 {
  s = Clist[i];
  n = NClist[i];

  if(n == nn)
   {
    k=1;
    for(j=0;j<n;j++)
     {
      if(s[j] != name[j])
       {
        k=0;
        break;
       }
     }
    if(k) break;
   }

  i++;
 }

if(i==IDF_PLIST_MAX)
 {
  *val = 0.0;
  k=0;
 }
else
 {
  k=1;
  *val = Vlist[i];
 }

return k;
}
