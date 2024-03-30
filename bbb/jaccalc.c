#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <string.h>
#if defined(UEDGE_WITH_OMP)
#include <omp.h>
#endif
typedef long Int;
typedef double real;

#define DO_TIMING       1
#define PRINT_LEVEL     1
#define uedge_min(a,b)  (((a)<(b)) ? (a) : (b))

/* Fortran protos */
extern void jac_calc_iv_(Int *iv, Int *neq, real *t, real* yl, real *yldot00, Int *ml, Int *mu,
                         real *wk, Int *nnzmx, real *jac, Int *ja, Int *ia, real *yldot_pert, Int *nnz);

int
jac_calc_seq_c_(Int  *neq,
                real *t,
                real *yl,
                real *yldot00,
                Int  *ml,
                Int  *mu,
                real *wk,
                Int  *nnzmx,
                real *jac,
                Int  *ja,
                Int  *ia,
                real *yldot_pert,
                Int  *nnz)
{
   if (PRINT_LEVEL > 0)
   {
      printf(" =============================================\n"
             " Jac_calc C SEQ version                       \n"
             "  ** n = %ld, nnzmx = %ld **                  \n"
             " =============================================\n",
             *neq, *nnzmx);
   }

   Int iv;

   *nnz = 1;
   for (iv = 1; iv <= *neq; iv++)
   {
      jac_calc_iv_(&iv, neq, t, yl, yldot00, ml, mu, wk, nnzmx, jac, ja, ia, yldot_pert, nnz);
   }
   ia[*neq] = *nnz;

   return 0;
}

#if defined(UEDGE_WITH_OMP)

#define USE_OMP_VERSION 1
#define SEQ_CHECK       0

typedef struct
{
   Int   num_threads;
   Int   neq;
   Int   nnzmx;
   Int   nnzmx_t;
   real  nnzmx_f;
   real *yl;
   real *yldot00;
   real *wk;
   real *jac;
   Int  *ja;
   real *yldot_pert;
} jac_calc_omp_data;

jac_calc_omp_data *jcod = NULL;

extern void omp_copy_module_(void);

/* partition range [0, 1, ..., n) */
void
jac_calc_1dpartition(Int  num_threads,
                     Int  thread_id,
                     Int *begin,
                     Int *end,
                     Int  n)
{
   Int n_per_thread = (n + num_threads - 1) / num_threads;
   *begin = uedge_min(n_per_thread * thread_id, n);
   *end = uedge_min(*begin + n_per_thread, n);
}

jac_calc_omp_data *
jac_calc_omp_data_create(void)
{
   jac_calc_omp_data *data = (jac_calc_omp_data *) calloc(1, sizeof(jac_calc_omp_data));

   return data;
}

/* allocate private memory for threads */
int
jac_calc_omp_data_init(Int                num_threads,
                       Int                neq,
                       Int                nnzmx,
                       real               nnzmx_f,
                       jac_calc_omp_data *data)
{
   if (num_threads == data -> num_threads &&
       neq         == data -> neq         &&
       nnzmx       == data -> nnzmx       &&
       nnzmx_f     == data -> nnzmx_f)
   {
      return 0;
   }

   const Int nnzmx_0 = (nnzmx + num_threads - 1) / num_threads;
   const Int nnzmx_t = (Int) (nnzmx * nnzmx_f + nnzmx_0 * (1.0 - nnzmx_f));

   data -> num_threads = num_threads;
   data -> neq         = neq;
   data -> nnzmx       = nnzmx;
   data -> nnzmx_t     = nnzmx_t;
   data -> nnzmx_f     = nnzmx_f;

   Int nt = num_threads - 1;

   data -> yl         = (real *) realloc(data -> yl,         nt * (neq + 2) * sizeof(real));
   data -> yldot00    = (real *) realloc(data -> yldot00,    nt * (neq + 2) * sizeof(real));
   data -> wk         = (real *) realloc(data -> wk,         nt * neq       * sizeof(real));
   data -> jac        = (real *) realloc(data -> jac,        nt * nnzmx_t   * sizeof(real));
   data -> ja         = (Int  *) realloc(data -> ja,         nt * nnzmx_t   * sizeof(Int));
   data -> yldot_pert = (real *) realloc(data -> yldot_pert, nt * neq       * sizeof(real));

   return 0;
}

int
jac_calc_omp_data_destroy(jac_calc_omp_data *data)
{
   if (!data)
   {
      return 0;
   }

   free(data -> yl);
   free(data -> yldot00);
   free(data -> wk);
   free(data -> jac);
   free(data -> ja);
   free(data -> yldot_pert);

   return 0;
}

int
jac_calc_omp_get_thread_data(Int                thread_id,
                             real              *yl,
                             real              *yldot00,
                             real              *wk,
                             real              *jac,
                             Int               *ja,
                             real              *yldot_pert,
                             jac_calc_omp_data *data,
                             jac_calc_omp_data *thread_data)
{
   thread_data -> num_threads = data -> num_threads;
   thread_data -> neq         = data -> neq;
   thread_data -> nnzmx       = data -> nnzmx;
   thread_data -> nnzmx_t     = data -> nnzmx_t;

   if (thread_id)
   {
      Int j, tid = thread_id - 1;
      thread_data -> yl         = data -> yl         + tid * (data -> neq + 2);
      thread_data -> yldot00    = data -> yldot00    + tid * (data -> neq + 2);
      thread_data -> wk         = data -> wk         + tid * (data -> neq);
      thread_data -> jac        = data -> jac        + tid * (data -> nnzmx_t);
      thread_data -> ja         = data -> ja         + tid * (data -> nnzmx_t);
      thread_data -> yldot_pert = data -> yldot_pert + tid * (data -> neq);

      for (j = 0; j < data -> neq + 2; j++)
      {
         thread_data -> yl[j]      = yl[j];
         thread_data -> yldot00[j] = yldot00[j];
      }
   }
   else
   {
      thread_data -> yl         = yl;
      thread_data -> yldot00    = yldot00;
      thread_data -> wk         = wk;
      thread_data -> jac        = jac;
      thread_data -> ja         = ja;
      thread_data -> yldot_pert = yldot_pert;
   }

   return 0;
}

int
jac_calc_omp_init(void)
{
   omp_set_dynamic(0);
   omp_copy_module_();

   return 0;
}

int
jac_calc_omp_c_(Int  *neq,
                real *t,
                real *yl,
                real *yldot00,
                Int  *ml,
                Int  *mu,
                real *wk,
                Int  *nnzmx,
                real  nnzmx_f,
                real *jac,
                Int  *ja,
                Int  *ia,
                real *yldot_pert,
                Int  *nnz)
{
   Int num_threads = omp_get_max_threads();

   if (PRINT_LEVEL > 0)
   {
      printf(" =============================================\n"
             " Jac_calc OpenMP C version, Num. Threads = %ld\n"
             "  ** n = %ld, nnzmx = %ld, nnzmx_f = %.2f **  \n"
#if SEQ_CHECK
             "  ** Check with serial version is ON **       \n"
#endif
             " =============================================\n",
             num_threads, *neq, *nnzmx, nnzmx_f);
   }

#if SEQ_CHECK
   Int nnz_s = 0;
   Int *ia_s = (Int *) calloc(*neq + 1, sizeof(Int));
   Int *ja_s = (Int *) calloc(*nnzmx, sizeof(Int));
   real *jac_s = (real *) calloc(*nnzmx, sizeof(real));
   jac_calc_seq_c_(neq, t, yl, yldot00, ml, mu, wk, nnzmx, jac_s, ja_s, ia_s, yldot_pert, &nnz_s);
#endif

   /* copy global variables in Fortran modules into thread privates*/
   jac_calc_omp_init();

   if (!jcod)
   {
      jcod = jac_calc_omp_data_create();
   }
   jac_calc_omp_data_init(num_threads, *neq, *nnzmx, nnzmx_f, jcod);

   Int *neq_all = (Int *) malloc((num_threads + 1) * sizeof(Int));
   Int *nnz_all = (Int *) malloc((num_threads + 1) * sizeof(Int));

   #pragma omp parallel
   {
      Int iv, iv_start, iv_end, thread_id = omp_get_thread_num();
      Int ml_t = *ml, mu_t = *mu, nnz_t = 1;
      real t_t = *t;
      jac_calc_omp_data jcod_t;

      jac_calc_omp_get_thread_data(thread_id, yl, yldot00, wk, jac, ja, yldot_pert, jcod, &jcod_t);
      jac_calc_1dpartition(num_threads, thread_id, &iv_start, &iv_end, *neq);

      for (iv = iv_start + 1; iv <= iv_end; iv++)
      {
         jac_calc_iv_(&iv, &jcod_t.neq, &t_t, jcod_t.yl, jcod_t.yldot00, &ml_t, &mu_t, jcod_t.wk,
                      &jcod_t.nnzmx_t, jcod_t.jac, jcod_t.ja, ia, jcod_t.yldot_pert, &nnz_t);
      }

      if (PRINT_LEVEL > 1)
      {
         printf("thread %ld: [%ld, %ld], nnz %ld\n", thread_id, iv_start+1, iv_end, nnz_t);
      }

      neq_all[thread_id + 1] = iv_end - iv_start;
      nnz_all[thread_id + 1] = nnz_t - 1;

      #pragma omp barrier
      #pragma omp master
      {
         neq_all[0] = nnz_all[0] = 0;
         for (iv = 1; iv < num_threads; iv ++)
         {
            neq_all[iv + 1] += neq_all[iv];
            nnz_all[iv + 1] += nnz_all[iv];
         }
         assert(neq_all[num_threads] == *neq);
      }
      #pragma omp barrier

      if (thread_id)
      {
         for (iv = neq_all[thread_id]; iv < neq_all[thread_id + 1]; iv ++)
         {
            ia[iv] = ia[iv] + nnz_all[thread_id];
         }

         for (iv = nnz_all[thread_id]; iv < nnz_all[thread_id + 1]; iv ++)
         {
            ja[iv] = jcod_t.ja[iv - nnz_all[thread_id]];
            jac[iv] = jcod_t.jac[iv - nnz_all[thread_id]];
         }
      }
      else
      {
         ia[*neq] = *nnz = nnz_all[num_threads] + 1;
      }
   } /* #pragma omp parallel */

#if SEQ_CHECK
   Int i;
   if (nnz_s != *nnz) { printf("SEQ_CHECK: nnz error %ld %ld\n", nnz_s, *nnz); exit(0); }
   for (i = 0; i <= *neq; i++) { if (ia_s[i] != ia[i]) { printf("SEQ_CHECK: ia error, i %ld, %ld %ld\n", i, ia_s[i], ia[i]); exit(0); } }
   for (i = 0; i <= *nnz; i++) { if (ja_s[i] != ja[i]) { printf("SEQ_CHECK: ja error, \n"); exit(0); } }
   for (i = 0; i <= *nnz; i++) { if (jac_s[i] != jac[i]) { printf("SEQ_CHECK: jac error, i = %ld: %e %e\n", i, jac_s[i], jac[i]); exit(0);} }
   printf(" SEQ_CHECK: OK...\n");
   free(ia_s);
   free(ja_s);
   free(jac_s);
#endif

   free(neq_all);
   free(nnz_all);

   jac_calc_omp_data_destroy(jcod); free(jcod); jcod = NULL;

   return 0;
}
#endif

int
jac_calc_c_(Int  *neq,
            real *t,
            real *yl,
            real *yldot00,
            Int  *ml,
            Int  *mu,
            real *wk,
            Int  *nnzmx,
            real *nnzmx_f,
            real *jac,
            Int  *ja,
            Int  *ia,
            real *yldot_pert,
            Int  *nnz)
{
#if DO_TIMING
   struct timeval t1, t2;
   gettimeofday(&t1, NULL);
#endif

#if USE_OMP_VERSION
   int ret = jac_calc_omp_c_(neq, t, yl, yldot00, ml, mu, wk, nnzmx, *nnzmx_f, jac, ja, ia, yldot_pert, nnz);
#else
   int ret = jac_calc_seq_c_(neq, t, yl, yldot00, ml, mu, wk, nnzmx, jac, ja, ia, yldot_pert, nnz);
#endif

#if DO_TIMING
   gettimeofday(&t2, NULL);
   double elapsedTime = (t2.tv_sec - t1.tv_sec);
   elapsedTime += (t2.tv_usec - t1.tv_usec) / 1e6;
   printf(" @@Time jac_calc_c@@    %e s\n", elapsedTime);
#endif

   return ret;
}

int jacmm_c_(Int *neq, real *jac, Int *ja, Int *ia)
{
   Int i, j;
   FILE *fp = fopen("Jacobian.IJ", "w");

   fprintf(fp, "%% %ld %ld %ld\n", *neq, *neq, ia[*neq] - 1);

   for (i = 0; i < *neq; i++)
   {
      for (j = ia[i]; j < ia[i + 1]; j++)
      {
         fprintf(fp, "%ld %ld %.15e\n", i + 1, ja[j - 1], jac[j - 1]);
      }
   }
   fclose(fp);

   return 0;
}
