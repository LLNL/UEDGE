
#include <stdio.h>
static int ic,iex,nd,nxi,nyi;
static float tp;
static float *xd,*yd,*zd;
static float xi,yi,zi;
static int ndim,nout;
static int *iwk;
static float *wk;
static float srelpr;
int first_interp;
static float max,min;
struct point {
    float r;
    float z;
};

float interp(cents,inimage,innum,p)
struct point *cents;
float *inimage,*p;
int innum;
{

    int i,j;
		int num;

    xi = p[0];
    yi = p[1];
    zi = 0;
    if(first_interp == 0){
        iex = 0;
        ic = 1;
        nd = innum;
        min = inimage[0];
        max = inimage[0];
        for(i = 0; i < innum ; ++i){
                if(max < inimage[i]) max = inimage[i];
                if(min > inimage[i]) min = inimage[i];
        }
        xd = malloc((nd + 1) * sizeof(float));
        yd = malloc((nd + 1)* sizeof(float));
        zd = malloc((nd + 1) * sizeof(float));
				num = 0;
        for(i = 0; i < innum ; ++i){
                xd[i] = cents[i].r;
                yd[i] = cents[i].z;
								/*printf("%g %g\n",xd[num],yd[num]);*/
								zd[num] = inimage[i];
								++num;
        }
				/*fflush(stdout);*/
				nd = num;
        tp = 10.0;
        nxi = 1;
        nyi = 1;
        ndim = 1;
        nout = 6;
        first_interp = 1;
				iwk = malloc(23 * nd * sizeof(int));
				wk = malloc(19 * nd * sizeof(float));
				maceps_(&srelpr);
        masub_(&srelpr,&ic,&iex,&nd,xd,yd,zd,&tp,&nxi,&nyi,&xi,&yi,&zi,&ndim,&nout,iwk,wk);

    } else {
        ic = 4;
        masub_(&srelpr,&ic,&iex,&nd,xd,yd,zd,&tp,&nxi,&nyi,&xi,&yi,&zi,&ndim,&nout,iwk,wk);
    }
    if(zi < min)zi = min;
    if(zi > max) zi = min;




    return(zi);
}
