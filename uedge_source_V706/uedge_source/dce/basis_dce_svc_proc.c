#include <rpc/rpc.h>
#include "basis_dce.h"
#include <stdio.h>
#include <string.h>
#include <malloc.h>

extern int first_interp;
extern int statistics;
int num_of_calls[100];
static int port;
static int status;
struct point {
    float r;
    float z;
};

rawimage return_image;

rawimage *rzxform_1(f,c)
ledge_pr_image *f;
CLIENT *c;
{
    struct point *cells;
    int i,j,k;
    int xlen,ylen;
    int cellnum,index;
    float *outimage;
    float *inimage;
    float tmp;
    int rank, dims[3];    /* the rank and dimensions of the data */
    float *data;                /* the data */


    status = 1;
    xlen = f->xlen;
    ylen = f->ylen;
    cells = malloc(xlen * ylen * 8 * sizeof(struct point));
    if(cells == NULL)
        printf("cells is NULL\n");
    outimage = malloc(f->outxlen * f->outylen * sizeof(float) * 2);
    if(outimage == NULL)
        printf("outimage is NULL\n");
    inimage = malloc(xlen * ylen * sizeof(float) * 2);
    if(inimage == NULL)
        printf("inimage is NULL\n");
    cellnum  = 1;
    for(j = 0 ; j < ylen; ++j){
        for(i = 0 ; i < xlen; ++i){
            inimage[j * xlen + i] = f->prdata.prdata_val[j * xlen + i];
        }
    }
    for(j = 0 ; j < ylen; ++j){
        for(i = 0 ; i < xlen; ++i){
            for(k = 1; k <= 4; ++k){
                index = i + j * xlen + k * xlen * ylen;
                cells[cellnum].r = f->r_verts.r_verts_val[index];
                cells[cellnum].z = f->z_verts.z_verts_val[index];
                ++cellnum;
            }
            ++cellnum;
            ++cellnum;

        }
    }
		first_interp = 0;
    ledgecvt(xlen,ylen,f->outxlen,f->outylen,f->rmin,f->rmax,f->zmin,
        f->zmax,cells,NULL,NULL,0,inimage,outimage,0);


    rank = 2;
		if(return_image.rawimage_len)free(return_image.rawimage_val);
		return_image.rawimage_len = f->outxlen * f->outylen;
		return_image.rawimage_val = outimage;

    free(cells);
    free(inimage);
		/*++num_of_calls[RZXFORM];*/
    return(&return_image);
}
rawimage *rzxform_int_1(f,c)
ledge_pr_intimage *f;
CLIENT *c;
{
    struct point *cells,*int_pts;
		float *int_vals;
    int i,j,k;
    int cellnum;
    float *outimage;
    float *inimage;
    float tmp;
    int rank, dims[3];    /* the rank and dimensions of the data */
    float *data;                /* the data */

    status = 1;
    cells = malloc(f->xlen * f->ylen * 8 * sizeof(struct point));
    int_pts = malloc(f->r_int.r_int_len * sizeof(struct point));
    int_vals = malloc(f->v_int.v_int_len * sizeof(struct point));
    if(cells == NULL)
        printf("cells is NULL\n");
    if(int_pts == NULL)
        printf("int_pts is NULL\n");
    if(int_vals == NULL)
        printf("int_vals is NULL\n");
    outimage = malloc(f->outxlen * f->outylen * sizeof(float) * 2);
    if(outimage == NULL)
        printf("outimage is NULL\n");
    inimage = malloc(f->xlen * f->ylen * sizeof(float) * 2);
    if(inimage == NULL)
        printf("inimage is NULL\n");
    cellnum  = 1;
    for(i = 0 ; i < f->ylen; ++i){
        for(j = 0 ; j < f->xlen; ++j){
            inimage[i * f->xlen + j] = f->prdata.prdata_val[j * f->ylen + i];
        }
    }
		for(i = 0; i < f->r_int.r_int_len; ++i){
			 int_pts[i].r = f->r_int.r_int_val[i];
			 int_pts[i].z = f->z_int.z_int_val[i];
			 int_vals[i] = f->v_int.v_int_val[i];
		}
    for(i = 0 ; i < f->ylen; ++i){
        for(j = 0 ; j < f->xlen; ++j){
            for(k = 1; k <= 4; ++k){
                cells[cellnum].r = f->r_verts.r_verts_val[i + j * f->ylen +
                                                                                                                 k * f->xlen * f->ylen];
                cells[cellnum].z = f->z_verts.z_verts_val[i + j * f->ylen +
                                                                                                                 k * f->xlen * f->ylen];
                ++cellnum;
            }
            ++cellnum;
            ++cellnum;

        }
    }
		first_interp = 0;
    ledgecvt(f->xlen,f->ylen,f->outxlen,f->outylen,f->rmin,f->rmax,f->zmin,
        f->zmax,cells,int_pts,int_vals,f->v_int.v_int_len,inimage,outimage,1);


    rank = 2;
		if(return_image.rawimage_len)free(return_image.rawimage_val);
		return_image.rawimage_len = f->outxlen * f->outylen;
		return_image.rawimage_val = outimage;

    free(cells);
		free(int_pts);
		free(int_vals);
    free(inimage);
		/*++num_of_calls[RZXFORM];*/
    return(&return_image);
}

freemal(a,len)
char **a;
int len;
{
    if(*a){
			 free((char *)*a);
    }
		*a = (char *)calloc(1,len);
}
