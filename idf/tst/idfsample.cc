#include <stream.h>
#include <strstream.h>

#define IDF_CPP 1
#define IDF_CPP_STYLE 1
#define IDF_FORTRAN_TYPE 3
#include "idfuser.h"


class COUET {
            int ierror;
     public:
            double Lz;
            float  Knu,Uw,Tw;
            double Pr[3];
            int    Nz;
            int    Nv[2];

            COUET(char*,char*);
       void COUETprn();
        int COUETerr();      
            };


COUET :: COUET(char* Fname, char* Sname)
{
//assume the `byte boundary' rule
int amode = 7;

//initialization of Target String
char Target_String[14] = "dfffd[3]ii[2]";

// start the IDF package
   ierror = idf_init();
if(ierror) goto er;

// open the data file
   ierror = idf_open(Fname);
if(ierror) goto err;

// reading the data from file to the Cinp structure
   ierror = idf_get(Sname, Target_String, amode, &Lz);
if(ierror>0) ierror=0;

err: idf_close(); //close the data file

er: idf_finish(); //finish the IDF
}

int COUET :: COUETerr()
{
 return ierror; // transfer the error flag
}

void COUET :: COUETprn()
{
// function prints the imported data sets
cout << "Lz="<<Lz<<" Knu="<<Knu<<"\n";
cout << "Uw="<<Uw<<"  Tw="<<  Tw<<"\n";
cout << "Pr="<<Pr[0]<<" "<<Pr[1]<<" , "<<Pr[2]<<"\n";
cout << "Nz="<<Nz<<" Nv=" <<Nv[0]<<" , " <<Nv[1] <<"\n";
}


int main()
/* The first sample program in C++ (see the IDF manual) */
{
// initialization of input IDF data file String
static char File_Name[13] = "couetteS.dat";

// initialization of structure data name
char Class_Name[8] = "COUETTE";

int ier;

// create the class object with the IDF based constructor
COUET Cinp(File_Name , Class_Name);

// if no errors print input data
    ier = Cinp.COUETerr();
if(!ier) Cinp.COUETprn();

return ier;
}

