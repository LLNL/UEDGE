#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <assert.h>
#include <sys/time.h>

#define uedge_min(a,b)  (((a)<(b)) ? (a) : (b))

#define SEQ_CHECK

typedef long Int;
typedef double real;

extern void pandf1_(Int *xc, Int *yc, Int *ieq, Int *neq, real *time, real *yl, real *yldot);

extern void pandf_(Int *xc, Int *yc, Int *neq, real *time, real *yl, real *yldot);

extern void jac_calc_iv_(Int *iv, Int *neq, real *t, real* yl, real *yldot00, Int *ml, Int *mu,
                         real *wk, Int *nnzmx, real *jac, Int *ja, Int *ia, real *yldot_pert, Int *nnz,
                         real *tpandf, Int *npandf);

void uedge_GetThreadPartition(Int num_threads, Int thread_id, Int *begin, Int *end, Int n)
{
   Int n_per_thread = (n + num_threads - 1) / num_threads;
   *begin = uedge_min(n_per_thread * thread_id, n);
   *end = uedge_min(*begin + n_per_thread, n);
}

int pandf_time_(Int *neq, real *t, real *yl, real *yldot00, Int *ml, Int *mu, real *wk, Int *nnzmx, real *yldot_pert)
{
   Int *ia0 = (Int *) calloc(*neq + 1, sizeof(Int));
   Int *ja0 = (Int *) calloc(*nnzmx, sizeof(Int));
   real *jac0 = (real *) calloc(*nnzmx, sizeof(real));
   Int nnz0 = 1;
   real tpandf = 0.0;
   Int npandf = 0;

   Int iv;
   for (iv = 1; iv <= *neq; iv++)
   {
      jac_calc_iv_(&iv, neq, t, yl, yldot00,
                   ml, mu, wk,
                   nnzmx,
                   jac0, ja0, ia0,
                   yldot_pert, &nnz0, &tpandf, &npandf);
   }
   ia0[*neq] = nnz0;

   Int mone = -1;
   real time = 0.0;

   struct timeval t1, t2;
   gettimeofday(&t1, NULL);

   pandf_(&mone, &mone, neq, &time, yl, wk);

   gettimeofday(&t2, NULL);

   double elapsedTime = (t2.tv_sec - t1.tv_sec);
   elapsedTime += (t2.tv_usec - t1.tv_usec) / 1e6;

   printf("%e seconds elapsed\n", elapsedTime);

   free(ia0);
   free(ja0);
   free(jac0);

   return 0;
}

int jac_calc_c_(Int *neq, real *t, real *yl, real *yldot00, Int *ml, Int *mu,
                real *wk, Int *nnzmx, real *jac, Int *ja, Int *ia, real *yldot_pert, Int *nnz)
{
   printf("=================\n Jac calc C\n ===============\n");
   Int i;

   omp_set_num_threads(2);
   Int nt = omp_get_max_threads();
   printf("nt %ld\n", nt);

   real **jac_thread = (real **) calloc(nt, sizeof(real *));
   Int **ja_thread = (Int **) calloc(nt, sizeof(Int *));
   Int **ia_thread = (Int **) calloc(nt, sizeof(Int *));
   real **wk_thread = (real **) calloc(nt, sizeof(real *));
   real **yl_thread = (real **) calloc(nt, sizeof(real *));
   real **yldot00_thread = (real **) calloc(nt, sizeof(real *));
   Int *nnz_thread = (Int *) calloc(nt, sizeof(Int));
   real *tpandf_thread = (real *) calloc(nt, sizeof(real));
   Int *npandf_thread = (Int *) calloc(nt, sizeof(Int));
   Int *n_thread = (Int *) calloc(nt, sizeof(Int));
   real *t_thread = (real *) calloc(nt, sizeof(real));

   for (i = 0; i < nt; i++)
   {
      jac_thread[i] = (real *) calloc(*nnzmx, sizeof(real));
      ja_thread[i] = (Int *) calloc(*nnzmx, sizeof(Int));
      ia_thread[i] = (Int *) calloc(*neq + 1, sizeof(Int));
      wk_thread[i] = (real *) calloc(*neq, sizeof(real));
      yl_thread[i] = (real *) calloc(*neq + 2, sizeof(real));
      yldot00_thread[i] = (real *) calloc(*neq + 2, sizeof(real));
      Int j;
      for (j = 0; j < *neq; j++)
      {
         wk_thread[i][j] = wk[j];
      }
      for (j = 0; j < *neq + 2; j++)
      {
         yl_thread[i][j] = yl[j];
         yldot00_thread[i][j] = yldot00[j];
      }
      nnz_thread[i] = 1;
      t_thread[i] = *t;
   }

#ifdef SEQ_CHECK
   Int *ia0 = (Int *) calloc(*neq + 1, sizeof(Int));
   Int *ja0 = (Int *) calloc(*nnzmx, sizeof(Int));
   real *jac0 = (real *) calloc(*nnzmx, sizeof(real));
   Int nnz0 = 1;
   Int npandf_seq = 0;
   real tpandf_seq = 0.;

   Int iv;
   for (iv = 1; iv <= *neq; iv++)
   {
      jac_calc_iv_(&iv, neq, t, yl, yldot00,
                   ml, mu, wk,
                   nnzmx,
                   jac0, ja0, ia0,
                   yldot_pert, &nnz0, &tpandf_seq, &npandf_seq);
   }
   ia0[*neq] = nnz0;
#endif

   #pragma omp parallel
   {
      Int iv;
      Int iv_start, iv_end;
      Int thread_id = omp_get_thread_num();

      uedge_GetThreadPartition(nt, thread_id, &iv_start, &iv_end, *neq);

      n_thread[thread_id] = iv_end - iv_start;

      {
         for (iv = iv_start + 1; iv <= iv_end; iv++)
         {
            #pragma omp critical
            {
               jac_calc_iv_(&iv, neq, &t_thread[thread_id], yl_thread[thread_id], yldot00_thread[thread_id],
                            ml, mu, wk_thread[thread_id],
                            nnzmx,
                            jac_thread[thread_id], ja_thread[thread_id], ia_thread[thread_id],
                            yldot_pert, &nnz_thread[thread_id], &tpandf_thread[thread_id], &npandf_thread[thread_id]);
            }
         }
         printf("thread %ld: [%ld, %ld], nnz %ld\n", thread_id, iv_start+1, iv_end, nnz_thread[thread_id]);
      }
   }

   Int nnz_tot = 0;
   for (i = 0; i < nt; i++)
   {
      nnz_tot += nnz_thread[i] - 1;
   }
   *nnz = nnz_tot + 1;

   Int k = 0, s = 0, s2 = 0;
   Int *ia_i = ia;
   for (i = 0; i < nt; i++)
   {
      Int ni = n_thread[i];
      Int nnzi = nnz_thread[i] - 1;
      Int j;
      for (j = 0; j < nnzi; j++)
      {
         jac[k] = jac_thread[i][j];
         ja[k] = ja_thread[i][j];
         k++;
      }
      for (j = 0; j < ni; j++)
      {
         ia_i[j] = ia_thread[i][j+s2] + s;
      }
      ia_i += ni;
      s += nnzi;
      s2 += ni;
   }
   ia[*neq] = nnz_tot + 1;

   if (k != nnz_tot)
   {
      printf("k != nnz error\n");
      exit(0);
   }

#ifdef SEQ_CHECK
   if (nnz0 != *nnz) { printf("nnz error %ld %ld\n", nnz0, *nnz); exit(0); }

   for (i = 0; i <= *neq; i++)
   {
      if (ia0[i] != ia[i]) { printf("error ia. i %ld, %ld %ld\n", i, ia0[i], ia[i]); exit(0);}
   }
   for (i = 0; i <= *nnz; i++)
   {
      if (ja0[i] != ja[i]) { printf("error ja\n"); exit(0);}
   }
   for (i = 0; i <= *nnz; i++)
   {
      if (jac0[i] != jac[i]) { printf("error jac\n"); exit(0);}
   }
   free(ia0);
   free(ja0);
   free(jac0);
#endif

   for (i = 0; i < nt; i++)
   {
      free(jac_thread[i]);
      free(ja_thread[i]);
      free(ia_thread[i]);
      free(wk_thread[i]);
      free(yl_thread[i]);
      free(yldot00_thread[i]);
   }
   free(jac_thread);
   free(ja_thread);
   free(ia_thread);
   free(wk_thread);
   free(nnz_thread);
   free(tpandf_thread);
   free(npandf_thread);
   free(n_thread);
   free(t_thread);


   return 0;
}

void pandf1_c_(Int *xc, Int *yc, Int *iv, Int *neq, real *t, real *yl, real *wk)
{
   //Int nt_max = omp_get_max_threads();
   //Int nt = omp_get_num_threads();
   //printf("%ld %ld\n", nt_max, nt);
//#pragma omp critical
   {
      pandf1_(xc, yc, iv, neq, t, yl, wk);
   }
}

#if 0

int rltest1_(Int *xcptr, Int *ycptr, Int *ivptr, Int *neqptr, real *tptr, real *yl, real *wk)
{
   Int i, neq = *neqptr;

   omp_set_num_threads(6);
   Int nt = omp_get_max_threads();

   real **wk_thread = (real **) calloc(nt, sizeof(real *));
   real **yl_thread = (real **) calloc(nt, sizeof(real *));

   for (i = 0; i < nt; i++)
   {
      wk_thread[i] = (real *) calloc(neq, sizeof(real));
      yl_thread[i] = (real *) calloc(neq + 2, sizeof(real));
      Int j;
      for (j = 0; j < neq; j++)
      {
         wk_thread[i][j] = wk[j];
      }
      for (j = 0; j < neq + 2; j++)
      {
         yl_thread[i][j] = yl[j];
      }
   }

   pandf1_(xcptr, ycptr, ivptr, neqptr, tptr, yl, wk);

#pragma omp parallel
   {
      Int thread_id = omp_get_thread_num();
      Int xc = *xcptr, yc = *ycptr, iv = *ivptr, n = *neqptr;
      real t = *tptr;
      #pragma omp critical
      {
         pandf1_(&xc, &yc, &iv, &n, &t, yl_thread[thread_id], wk_thread[thread_id]);
      }
   }

   Int count = 0;
   real diff = 0.0;
   for (i = 0; i < neq; i++)
   {
      Int j;
      for (j = 0; j < nt; j++)
      {
         if (wk_thread[j][i] != wk[i])
         {
            count++;
            diff += fabs(wk_thread[j][i] - wk[i]);
            //printf("i %d, %e %e\n", i, wk_thread[0][i], wk_thread[1][i]);
         }
      }
   }
   printf("count %ld diff %e\n", count, diff);

   exit(0);

   return 0;
}

#endif

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
