#include "idflib.h"


#include "idf.h"
#include "idfusr.h"


int main()
{
int ierr=0;
static char FN[] = "idftst.dat";
int k,l;
register int i;

char   named[] = "a___b___c";
double abc=0.0;

char   namei[] = "iabc";
int    iabc=0;

char   names[] = "Sabc";
char   stroka[30];

char   namet[] = "Tabc";
char   text[30];

char   namedarr[] = "D_arr";
int    ndarr=5;
double Darr[5]={0.0,0.0,0.0,0.0,0.0};

char   nameStru[] = "Stru";
struct {char c; int i; float f; char cc[3]; double g; float h[3];} St;

char   nameStra[] = "Stra";
struct {char c; int i; float *f; char cc[3]; double *g; float h[3];} Sta;



   ierr = idf_ini(FN);
if(ierr) goto err;

/* idf_list_catalogL(); */
/* idf_list_catalog (); */


ierr = idf_d(named, &abc);
     printf("%s=%24.16e flag=%d\n",named,abc,ierr);

ierr = idf_i(namei, &iabc);
     printf("%s=%d flag=%d\n",namei,iabc,ierr);

ierr = idf_s(names, stroka, 30);
     printf("%s=%s flag=%d\n",names,stroka,ierr);

ierr = idf_t(namet, text, 30);
     printf("%s=%s flag=%d\n",namet,text,ierr);

ierr = idf_darr(namedarr, Darr, ndarr);
     printf("%s flag=%d\n",namedarr,ierr);
     for(i=0;i<ndarr;i++) printf("d(%d)=%e\n",i,Darr[i]);

ierr = idf_get_array(namedarr, "d[5]", Darr);
     printf("%s flag=%d\n",namedarr,ierr);
     for(i=0;i<ndarr;i++) printf("d(%d)=%e\n",i,Darr[i]);

     St.c=' '; St.i=0; St.f=-1.0; St.cc[0]='\0'; St.g=0.0;
     St.h[0]=-1; St.h[1]=-2; St.h[2]=-3;
ierr = idf_get(nameStru, "ci#f%s[3]d#3f", 7, &St);
     printf("%s flag=%d\n",nameStru,ierr);
     printf("c=%c i=%d f=%e cc=%s\n d=%e h=%e %e %e\n",
     St.c,St.i,St.f,St.cc,St.g,St.h[0],St.h[1],St.h[2]);
   
     Sta.c=' '; Sta.i=0; Sta.f=NULL; Sta.cc[0]='\0';
     Sta.g=(double*)malloc(2*sizeof(double));
     Sta.h[0]=-1; Sta.h[1]=-2; Sta.h[2]=-3;
ierr = idf_get(nameStra, "ci@f[2]%s[3],&2d#3f", 7, &Sta);
     printf("%s flag=%d\n",nameStra,ierr);
     printf("c=%c i=%d f=%e %e cc=%s\n d=%e %e h=%e %e %e\n",
     Sta.c,Sta.i,Sta.f[0],Sta.f[1],Sta.cc,Sta.g[0],Sta.g[1],
     Sta.h[0],Sta.h[1],Sta.h[2]);

if(Sta.g != NULL) free(Sta.g);
if(Sta.f != NULL) free(Sta.f);

err:
idf_end();

return ierr;
}
