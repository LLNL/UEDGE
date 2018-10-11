/*

The following code is based on algorithms ccw (p350) ,intersect (p351) and
inside (p354) from:
            Algorithms in C 
            by Robert Sedgewick
            Copyright 1990 by Addison-Wesley Publishing Co., Inc.
            ISBN 0-201-51425-7

Keyed in and adapted by:
            Bill Meyer
            Lawrence Livermore National Lab
            meyer8@llnl.gov
            (510)423-6039

Changes from the text:
            Made max an input parameter to inside. 

*/
#include <stdio.h>

#define VMAX 8  /*  maximum number of polygon vertices  */

struct point {
    float x;
    float y;
};
struct line {
    struct point p1,
    p2;
};

int ccw(p0,p1,p2)
struct point p0,p1,p2;
/* when traveling from p0 to p1 to p2 does one
     turn clockwise (return -1) , turn counterclockwise (return 1), 
     or are the points colinear (return 0) */
{

    float dx1,dx2,dy1,dy2;
    dx1 = p1.x - p0.x;
    dy1 = p1.y - p0.y;
    dx2 = p2.x - p0.x;
    dy2 = p2.y - p0.y;
    if(dx1*dy2 > dy1*dx2) return 1;        /*    counterclockwise  */
    if(dx1*dy2 < dy1*dx2) return -1;        /*    clockwise  */
    if((dx1*dx2 < 0) || (dy1*dy2 < 0)) return -1;   /*p0 between p2 and p1*/
    if((dx1*dx1+dy1*dy1) < (dx2*dx2+dy2*dy2)) return 1;/*p1 between p0 and p2*/
    return 0;   /*  else p2 between p0 and p1  */
}


/*
Given two lines, if the endpoints of one line are on different
sides (different values of ccw) of the other line they must intersect.
*/
int intersect(l1,l2)
struct line l1,l2;
{

    int status;
    status = ((ccw(l1.p1,l1.p2,l2.p1)
        *ccw(l1.p1,l1.p2,l2.p2)) <= 0)
        && ((ccw(l2.p1,l2.p2,l1.p1)
        *ccw(l2.p1,l2.p2,l1.p2)) <= 0);
    return(status);

}
int inside(t,p,N,max)
struct point t,*p;
int N;
float max;
{

    int i,count=0,j=0;
    struct line lt,lp;
    p[0] = p[N];
    p[N+1] = p[1];
    lt.p1 = t;
    lt.p2 = t;
    lt.p2.x = max;
    for( i = 1; i <= N; i++){
        lp.p1 = p[i];
        lp.p2 = p[i];
        if(!intersect(lp,lt)){
            lp.p2 = p[j];
            j = i;
            if(intersect(lp,lt))count++;
        }
    }
    return count & 1;
}
int ontheline(t,p,N)
struct point t,*p;
int N;
{

    int i,count=0,j=0;
    p[0] = p[N];
    p[N+1] = p[1];
    for( i = 1; i <= N; i++){
        if(p[i].x == t.x && p[i].y == t.y){
						fprintf(stderr,"direct hit\n");
						return(1);
				}
				if(ccw(p[i],p[i+1],t) == 0){
					/*fprintf(stderr,"on the line\n");*/
					return(1);
				}
    }
		return(0);
}
