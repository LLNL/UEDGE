#include <stdio.h>

#define IDF_CPP 0
#define IDF_CPP_STYLE 1
#define IDF_FORTRAN_TYPE 3
#include "idfuser.h"



int main()
/* The first sample program in C (see IDF manual) */
{
/*initialization of input IDF data file String*/
static char Fname[12] = "couette.dat";

/* declaration of input parameters */
double Lz,Pr[3];
float Knu,Uw,Tw;
int Nz,Nv[2];
int ier;

/* start the IDF package */
   ier=idf_init();
if(ier) goto er;

/* open the data file */
   ier=idf_open(Fname);
if(ier) goto err;

/* reading the data from file */
   ier=idf_d   ("Lz" ,    &Lz );
if(ier) goto err;
   ier=idf_f   ("Knu",    &Knu);
if(ier) goto err;
   ier=idf_f   ("Uw" ,    &Uw );
if(ier) goto err;
   ier=idf_f   ("Tw" ,    &Tw );
if(ier) goto err;
   ier=idf_darr("Pr" ,  Pr , 3  );
if(ier!=3) goto err;
   ier=idf_i   ("Nz" ,    &Nz );
if(ier) goto err;
   ier=idf_iarr("Nv" ,  Nv , 2 );
if(ier!=2) goto err;

/* print the imported data sets*/
printf("Lz=%12.3e Knu=%12.3e\n",   Lz,Knu);
printf("Uw=%12.3e  Tw=%12.3e\n",   Uw,Tw);
printf("Pr=%12.3e %12.3e %12.3e\n",Pr[0],Pr[1],Pr[2]);
printf("Nz=%5d Nv=%5d, %5d\n",     Nz,Nv[0],Nv[1]);

/* close the data file */
err: idf_close();

/* finish the IDF */
er: idf_finish();

return ier;
}
