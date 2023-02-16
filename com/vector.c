
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include "llnltyps.h"
#include <stdlib.h>

#define ZERO   RCONST(0.0)  /* real 0.0   */

/* typedef double real;
typedef  int integer; */

/* Part I: Real N_Vector Kernel Prototypes (Machine Environment-Independent) */

/***************************************************************
 *                                                             *
 * Type: N_Vector                                              *
 *-------------------------------------------------------------*
 * The type N_Vector is an abstract vector type. The fields of *
 * its concrete representation should not be accessed          *
 * directly, but rather through the macros given below.        *
 *                                                             *
 * A user may assume that the N components of an N_Vector      *
 * are stored contiguously. A pointer to the first component   *
 * can be obtained via the macro N_VRDATA.                      *
 *                                                             *
 ***************************************************************/

/*typedef struct {
  integer length;
  real   *rdata;
  integer   *idata;
} *N_Vector; */


/*************************************************************************
 *                                                                       *
 * Memory Allocation and Deallocation: NVNEWR, NVNEWI, NVFREE         *
 *                                                                       *
 *************************************************************************/


/***************************************************************
 *                                                             *
 * Function : NVNEWR                                          *
 * Usage    : x = NVNEWR(N);                                  *
 * Function : NVNEWI                                          *
 * Usage    : x = NVNEWI(N);                                  *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns a new N_Vector of length N. The parameter machEnv   *
 * is a pointer to machine environment-specific information.   *
 * It is ignored in the sequential machine environment and the *
 * user in this environment should simply pass NULL for this   *
 * argument. If there is not enough memory for a new N_Vector, *
 * then N_VNew returns NULL.                                   *
 *                                                             *
 ***************************************************************/

   void * nvnewr_(integer *n);
   void * nvnewi_(integer *n);


/***************************************************************
 *                                                             *
 * Function : NVFREER                                          *
 * Usage    : NVFREER(x);                                      *
 * Function : NVFREEI                                          *
 * Usage    : NVFREEI(x);                                      *
 *-------------------------------------------------------------*
 *                                                             *
 * Frees the N_Vector x. It is illegal to use x after the call *
 * N_VFREER(x) and  N_VFREEI(x).                               *
 *                                                             *
 ***************************************************************/

void nvfreer_(real x);
void nvfreei_(void *x);

/***************************************************************/
void * nvnewr_(integer *NLENGTH)
{
  integer N, i;
  real *vr;
 
  N = *NLENGTH;

  if (N <= 0) return(NULL);

  vr = (real *) malloc(N * sizeof(real));
  if (vr == NULL) {
    free(vr);
    return(NULL);
  }

    for (i=0; i < N; i++)
      vr[i] =  ZERO;

  return(vr);
}
/***************************************************************/

void nvfreer_(real x)
{
  /*  free(x); */
}
/***************************************************************/
void * nvnewi_(integer *NLENGTH)
{
  integer N, i;
  integer *vr;
 
  N = *NLENGTH;

  if (N <= 0) return(NULL);

  vr = (integer *) malloc(N * sizeof(integer));
  if (vr == NULL) {
    free(vr);
    return(NULL);
  }

    for (i=0; i < N; i++)
      vr[i] =  0;

  return(vr);
}
/***************************************************************/

void nvfreei_(void *x)
{
  free(x);
}

