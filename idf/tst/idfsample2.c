#include <stdio.h>

#define IDF_CPP 0
#define IDF_CPP_STYLE 1
#define IDF_FORTRAN_TYPE 3
#include "idfuser.h"


typedef struct {
                double Lz;
                float  Knu,Uw,Tw;
                double Pr[3];
                int    Nz;
                int    Nv[2];
               } COUETTE;


int Couette_input(Cinp)
COUETTE *Cinp;
{
/*initialization of input IDF data file String*/
char Fname[13] = "couetteS.dat";

/* initialization of structure data name */
char StructName[8] = "COUETTE";

/*initialization of alignment mode */
int amode = 7;  /*assume the `natural boundary' rule*/

/*initialization of Target String */
char Target_String[14] = "dfffd[3]ii[2]";

int ier;

/* start the IDF package */
   ier=idf_init();
if(ier) goto er;

/* open the data file */
   ier=idf_open(Fname);
if(ier) goto err;

/* reading the data from file to the Cinp structure*/
   ier = idf_get(StructName, Target_String, amode, Cinp);
if(ier>0) ier=0;

/* close the data file */
err: idf_close();

/* finish the IDF */
er:  idf_finish();

 return ier;
}


int main()
/* The second sample program in C (see IDF manual) */
{
/* declaration of input structure */
COUETTE Cinput;
int ier=0;

Cinput.Lz    = 0.0;
Cinput.Knu   = 0.0;
Cinput.Uw    = 0.0;
Cinput.Tw    = 0.0;
Cinput.Pr[0] = 0.0;
Cinput.Pr[1] = 0.0;
Cinput.Pr[2] = 0.0;
Cinput.Nz    = 0;
Cinput.Nv[0] = 0;
Cinput.Nv[1] = 0;

/* import data from file */
    ier=Couette_input(&Cinput);
 if(ier) goto err;

/* print the imported data sets*/
printf("Lz=%12.3e Knu=%12.3e\n", Cinput.Lz,Cinput.Knu);
printf("Uw=%12.3e  Tw=%12.3e\n", Cinput.Uw,Cinput.Tw);
printf("Pr=%12.3e, %12.3e, %12.3e\n",
       Cinput.Pr[0],Cinput.Pr[1],Cinput.Pr[2]);
printf("Nz=%5d Nv=%5d, %5d\n",
       Cinput.Nz,Cinput.Nv[0],Cinput.Nv[1]);

err: return ier;
}
