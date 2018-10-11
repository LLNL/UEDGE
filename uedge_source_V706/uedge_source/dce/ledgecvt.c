/* 
         This code takes output from ledge, which is a 2d array
         in psi,r space, and outputs a 2d array in r,z space. 
*/
#include <stdio.h>

struct point {
    float r;
    float z;
};
struct line {
    struct point p1,
    p2;
};
static float rstep,zstep;
static float rmin,rmax,zmin,zmax;

ledgecvt(inx,iny,outx,outy,prmin,prmax,pzmin,pzmax,cells,int_pts,int_vals,intnum,inimage,outimage,interp_fl)
int intnum,inx,iny,outx,outy;
float prmin,prmax,pzmin,pzmax;
float *outimage,*inimage;
struct point *cells,*int_pts,*int_vals;
int interp_fl;
{

    struct point maxrpoint(),minrpoint();
    struct point maxzpoint(),minzpoint();
    float interp();
    int i,j;    /*    indices into input array  */
    int k;
    int m,n;    /*  indices into output array  */
    int mstart,nstart;
    int mstop,nstop;
    struct point p;
    struct point Mr,mr;
    struct point Mz,mz;


    rmax = prmax;
    zmax = pzmax;
    rmin = prmin;
    zmin = pzmin;
    rstep = (rmax - rmin) / (outx - 1);
    zstep = (zmax - zmin) / (outy - 1);


    for(j = 0; j < iny; ++j){
        for(i = 0; i < inx; ++i){
            forder(&cells[j * inx * 6 + i * 6]);
        }
    }
    memset(outimage,0,outx * outy * sizeof(float));
    for(j = 0; j < iny; ++j){
        for(i = 0; i < inx ; ++i){
            Mr = maxrpoint(&cells[j * inx * 6 + i * 6]);
            mr = minrpoint(&cells[j * inx * 6 + i* 6 ]);
            Mz = maxzpoint(&cells[j * inx * 6 + i* 6 ]);
            mz = minzpoint(&cells[j * inx * 6 + i* 6 ]);
						/*printf("%g %g %g %g\n",mr.r,Mr.r,mz.z,Mz.z);*/
						if(Mr.r <= rmin || mr.r >= rmax ||
						   Mz.z <= zmin || mz.z >= zmax)continue;
            nstart = (mz.z - zmin) / zstep ;
            nstop = (Mz.z - zmin) / zstep;
            mstart = (mr.r - rmin) / rstep ;
            mstop = (Mr.r - rmin) / rstep;
            if(nstart < 0)nstart = 0;
            if(nstop >= outy) nstop = outy - 1;
            if(mstart < 0)mstart = 0;
            if(mstop >= outx) mstop = outx - 1;
            for(n = nstart; n <= nstop; ++n){
                for(m = mstart; m <= mstop; ++m){
                    p.r = (float)m * rstep + rmin;
                    p.z = (float)n * zstep + zmin;
                    if(inside(p,&cells[j * inx * 6+ i* 6],4,4 * rmax) ||
                       ontheline(p,&cells[j * inx * 6+ i* 6],4)){
                        if(interp_fl)
                            outimage[n * outx + m] = 
                                interp(int_pts,int_vals,intnum,&p);
                        else
                            outimage[n * outx + m] = inimage[j * inx + i];
                    }
                }
            }
        }
    }
}

struct point maxrpoint(p)
struct point *p;
{
    int i;
    struct point mp;
    mp.r = rmin;
    for(i = 1; i <= 4; ++i)
        if(mp.r < p[i].r){
            mp.z = p[i].z;
            mp.r = p[i].r;
        }

    if(mp.r > rmax)mp.r = rmax;
    if(mp.z > zmax)mp.z = zmax;
    return mp;
}
struct point minrpoint(p)
struct point *p;
{
    int i;
    struct point mp;
    mp.r = rmax;
    for(i = 1; i <= 4; ++i)
        if(mp.r > p[i].r){
            mp.z = p[i].z;
            mp.r = p[i].r;
        }

    if(mp.r < rmin)mp.r = rmin;
    if(mp.z < zmin)mp.z = zmin;
    return mp;
}
struct point maxzpoint(p)
struct point *p;
{
    int i;
    struct point mp;
    mp.z = zmin;
    for(i = 1; i <= 4; ++i)
        if(mp.z < p[i].z){
            mp.z = p[i].z;
            mp.r = p[i].r;
        }

    if(mp.r > rmax)mp.r = rmax;
    if(mp.z > zmax)mp.z = zmax;
    return mp;
}
struct point minzpoint(p)
struct point *p;
{
    int i;
    struct point mp;
    mp.z = zmax;
    for(i = 1; i <= 4; ++i)
        if(mp.z > p[i].z){
            mp.z = p[i].z;
            mp.r = p[i].r;
        }

    if(mp.r < rmin)mp.r = rmin;
    if(mp.z < zmin)mp.z = zmin;
    return mp;
}


outrind(r)
float r;
{

    return (int)((r - rmin)/ rstep);

}
outzind(z)
float z;
{

    return (int)((z - zmin)/ zstep);

}
forder(p)
struct point *p;
{
    struct point t;
    int i;


    if(!order(p)){
        t = p[3];
        p[3]  = p[2];
        p[2] = t;
        order(p);
    }
}


order(p)
struct point *p;
{
    struct point t;
    if(ccw(p[1],p[2],p[3]) != ccw(p[2],p[3],p[4])){
        /*fprintf(stderr,"1:vertices not in order \n");*/
    } else {
        return 1;
    }

    t = p[3];
    p[3] = p[4];
    p[4] = t;

    if(ccw(p[1],p[2],p[3]) != ccw(p[2],p[3],p[4])){
        /*fprintf(stderr,"2:vertices not in order \n");*/
    } else {
        return 1;
    }

    /*fprintf(stderr,"3:vertices not in order \n");*/
    return 0;


}
