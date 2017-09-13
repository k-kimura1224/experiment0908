/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_forward.c
 * @brief  heuristic for this problem
 * @author Keiji Kimura
 *
 * This file implements heuristic algorithm for this problem. It uses the stepwise method with forward
 * selection, and find a feasible solution of a given subproblem. The forward selection involves
 * starting with no variables in the model. This algorithm calculates AIC/BIC of the model that
 * increased each variable, selects adding the variable that improves the model the most, and repeat
 * this process until none improves the model.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>
#include <string.h>
#include <math.h>

#include "heur_forward.h"
#include "probdata_linereg.h"
#include "convenient_tool.h"
#include "set_myparameter.h"
#include "call_cblas.h"

#define HEUR_NAME             "forward"
#define HEUR_DESC             "primal heuristic using forward selection"
#define HEUR_DISPCHAR         'f'
#define HEUR_PRIORITY         1000
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         10
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define MYPARA_LOG            0

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   int   a;
};


/*
 * Local methods
 */

/** calculate residual sum of squares */
static
SCIP_Real calcRSSvalue(
   int                   dim,                /**< dimension */
   SCIP_Real*            a,                  /**< solution */
   SCIP_Real*            xy,                 /**< (submatrix X)' y */
   SCIP_Real             r                   /**< constant term */
   )
{
   assert(dim > 0);
   assert(a != NULL);
   assert(xy != NULL);

   return - SCIPcblasDdot(a, xy, dim) + r;
}


/** calculate objective value  */
static
SCIP_Real calcObjval(
   int                   n,                  /**< the number of data points */
   int                   dim,                /**< dimension */
   SCIP_Real             para_regterm,       /**< parameter to control the regularization term */
   SCIP_Real             rss                 /**< residual sum of squares */
   )
{
   assert(n > 0);
   assert(dim > 0);
   assert(EPSEQ(para_regterm, 2.0, 1e-08) || EPSEQ(para_regterm, log((SCIP_Real) n), 1e-08));
   assert(rss > 0.0);

   return (SCIP_Real) n * log(rss) + para_regterm * (SCIP_Real) dim;
}


/** compute a solution and set new inverse matrix */
static
SCIP_RETCODE  computeSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   n,                  /**< the number of data points */
   int                   dim,                /**< dimension */
   SCIP_Real*            sub_x,              /**< array [n*(dim-1)] */
   SCIP_Real*            x,                  /**< array [n] */
   SCIP_Real*            inv_xx,             /**< array [(dim-1)*(dim-1)] */
   SCIP_Real*            sub_xy_new,         /**< array [dim] */
   SCIP_Real*            a,                  /**< array [dim] */
   SCIP_Real*            inv_xx_new          /**< array [dim*dim] */
   )
{
   SCIP_Real* vec_b;
   SCIP_Real* vec_c;
   SCIP_Real* vec_d;
   SCIP_Real u;
   SCIP_Real* vec_v;
   SCIP_Real* mat_v;
   SCIP_Real* mat_z;

   int j;
   int t;
   int dim_1 = dim - 1;
   int dimdim_1 = dim_1 * dim_1;
   int dimdim = dim * dim;

   assert(scip != NULL);
   assert(n > 0);
   assert(dim > 0);
   assert(sub_x != NULL);
   assert(x != NULL);
   assert(inv_xx != NULL);
   assert(sub_xy_new != NULL);
   assert(a != NULL);
   assert(inv_xx_new != NULL);

   /* alloc */
   SCIP_CALL( SCIPallocBufferArray(scip, &vec_b, dim_1));
   SCIP_CALL( SCIPallocBufferArray(scip, &vec_c, dim_1));
   SCIP_CALL( SCIPallocBufferArray(scip, &vec_d, n));
   SCIP_CALL( SCIPallocBufferArray(scip, &vec_v, dim_1));
   SCIP_CALL( SCIPallocBufferArray(scip, &mat_v, dimdim_1));
   SCIP_CALL( SCIPallocBufferArray(scip, &mat_z, dimdim));

   /* 1. b <- sub_x^t x_i */
   SCIP_CALL( SCIPcblasDgemv3(sub_x, n, dim_1, x, vec_b));

   /* 2. c <- inv_xx b */
   SCIP_CALL( SCIPcblasDgemv2(inv_xx, dim_1, dim_1, vec_b, vec_c));

   /* 3. d <- - sub_x c + x_i */
   SCIP_CALL( SCIPcblasDgemv1(sub_x, n, dim_1, vec_c, x, -1.0, 1.0, vec_d));

   /* 4. u <- 1/<x_i, d> */
   u = 1.0 / SCIPcblasDdot(x, vec_d, n);

   /* 5. v <- - u c */
   SCIP_CALL( SCIPcblasDscal(vec_c, dim_1, -u, vec_v));

   /* 6. V <- inv_xx + u c c^t */
   SCIP_CALL( SCIPcblasDger(inv_xx, vec_c, vec_c, dim_1, dim_1, u, mat_v));

   /* 7.
    * Z = ( V  v )
    *     ( v' u )
    */
   /* V */
   for( j = 0; j <dim_1; j++ )
   {
      for( t = 0; t < dim_1; t++ )
         *(mat_z + j + (t*dim) ) = SCIPmatColMajor(mat_v, dim_1, j, t);
   }
   /* v */
   for( j = 0; j < dim_1; j++ )
   {
      *(mat_z + dim_1 + (j*dim) ) = vec_v[j];
      *(mat_z + j + (dim_1*dim)) = vec_v[j];
   }

   /* u */
   *(mat_z + dim_1 + (dim_1*dim)) = u;

   /* 8. a_old <- Z (sub_xy_new) */
   SCIP_CALL( SCIPcblasDgemv2(mat_z, dim, dim, sub_xy_new, a));

   for( j = 0; j < dimdim; j++ )
      inv_xx_new[j] = mat_z[j];

   SCIPfreeBufferArray(scip, &vec_b);
   SCIPfreeBufferArray(scip, &vec_c);
   SCIPfreeBufferArray(scip, &vec_d);
   SCIPfreeBufferArray(scip, &vec_v);
   SCIPfreeBufferArray(scip, &mat_v);
   SCIPfreeBufferArray(scip, &mat_z);

   return SCIP_OKAY;
}

/** check feasiblity by using linear dependence */
static
SCIP_Bool checkFeasiblity(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   p,                  /**< the number of explanatory variables */
   int                   ndep,               /**< the number of linerly dependent sets */
   int                   ubnumex,            /**< parameter to control the number of selected explanatory variables */
   SCIP_PROBDATA*        probdata,           /**< user problem data */
   int*                  list,               /**< variable information */
   int                   index               /**< index of variable */
   )
{
   int* indexdepsets;
   int* sizedepsets;

   int i;
   int j;
   int m;
   int size;
   int buf;
   int* ip;
   SCIP_Bool memo;
   SCIP_Bool in;

   assert(scip != NULL);
   assert(p > 0);
   assert(ndep > 0);
   assert(probdata != NULL);
   assert(list != NULL);
   assert(index >= 0 && index < p && list[index] == 0);
   assert(ubnumex == -1 || (ubnumex > 0 && ubnumex < p) );

   indexdepsets = SCIPprobdataGetindexdepsets(probdata);
   sizedepsets = SCIPprobdataGetsizedepsets(probdata);

   assert(indexdepsets != NULL);
   assert(sizedepsets != NULL);

   buf = 0;
   ip = &indexdepsets[0];

   for( i = 0; i < ndep; i++ )
   {
      ip += buf;
      size = sizedepsets[i];
      if( ubnumex == -1 || size <= ubnumex )
      {
         memo = FALSE;
         in = FALSE;
         for( j = 0; j < size; j++ )
         {
            m = *(ip+j);

            if( index == m )
               in = TRUE;
            else if( index < m && in == FALSE )
            {
               memo = TRUE;
               break;
            }

            if( index != m && list[m] <= 0 )
            {
               memo = TRUE;
               break;
            }
         }

         if( memo == FALSE )
         {
            return FALSE;
         }
      }

      buf = size;
   }

   return TRUE;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyForward)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurForward(scip) );

   return SCIP_OKAY;
}


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeForward)
{   /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecForward)
{  /*lint --e{715}*/

   /* probdata */
   SCIP_PROBDATA* probdata;
   int n;
   int p;
   int ndep;
   int ubnumex;
   SCIP_Real para_regterm;
   SCIP_Real* orig_x;
   SCIP_Real* orig_xy;
   SCIP_Real* orig_xx;
   SCIP_Real r;
   SCIP_VAR** var_z;
   SCIP_Bool savememory;

   int dim;
   int* list;
   SCIP_Bool checkld;

   SCIP_Real* a_old;
   SCIP_Real* sub_x;
   SCIP_Real* inv_xx;
   SCIP_Real* sub_xy;

   SCIP_Real* a;
   SCIP_Real* a_new;
   SCIP_Real* sub_xy_new;
   SCIP_Real* inv_xx_store;
   SCIP_Real* inv_xx_new;

   SCIP_Real rss;
   SCIP_Real rss_new;
   SCIP_Real objval = 1e+06;
   SCIP_Real objval_new;

   /* set solution */
   int nsols;
   int store;
   SCIP_SOL** sols;
   SCIP_VAR* var;

   int ndim;
   int ni;
   int dimdim;
   int dim_1;
   int memo;

   int i;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);

#if MYPARA_LOG
   printf("forward selection");
   SCIPprintLongLine();
#endif

   /* get heuristic data */
   /*
   SCIP_HEURDATA* heurdata;
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(lastsolindices != NULL);
   */

   /* get values from probdata */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   n = SCIPprobdataGetNdatas(probdata);
   p = SCIPprobdataGetNexvars(probdata);
   para_regterm = SCIPprobdataGetPara_RegTerm(probdata);
   ndep = SCIPprobdataGetNdep(probdata);
   ubnumex = SCIPprobdataGetPara_NumEx(probdata);

   orig_x = SCIPprobdataGetx(probdata);
   orig_xx = SCIPprobdataGetxx(probdata);
   orig_xy = SCIPprobdataGetxy(probdata);
   r = SCIPprobdataGetyy(probdata);
   savememory = SCIPprobdataGetSaveMemory(probdata);

   var_z = SCIPprobdataGetVars_z(probdata);

   assert(n > 0);
   assert(p > 0);
   assert(EPSEQ(para_regterm, 2.0, 1e-08) || EPSEQ(para_regterm, log((SCIP_Real) n), 1e-08));
   assert(ndep >= 0);
   assert(ubnumex == -1 || (ubnumex > 0 && ubnumex < p) );

   assert(orig_x != NULL);
   assert((savememory == FALSE && orig_xx != NULL)
         || (savememory == TRUE && orig_xx == NULL));
   assert(orig_xy != NULL);

#if MYPARA_LOG
   {
      /* branching info */
      int ublb;
      int* branchinfo;

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &branchinfo, 3*p));
      /* get branching info */
      SCIP_CALL( SCIPinitIntArrayZero(3*p, branchinfo));

      for(i = 0; i < p; ++i )
      {
         ublb = SCIPround(scip, SCIPcomputeVarUbLocal(scip, var_z[i])
               + SCIPcomputeVarLbLocal(scip, var_z[i]));
         *(branchinfo+(ublb*p)+i) = 1;
      }

      for( i = 0; i < 3; i++ )
      {
         for( j = 0; j < p; j++ )
         {
            printf("%d, ", *(branchinfo + (i*p) + j));
         }
         SCIPstartNewLine();
      }

      /* free */
      SCIPfreeBufferArray(scip, &branchinfo);
   }
#endif


#if MYPARA_LOG
   printf("step1:\n");
#endif

   /* alloc */
   SCIP_CALL( SCIPallocBufferArray(scip, &list, p));
   SCIP_CALL( SCIPallocBufferArray(scip, &sub_x, n));
   SCIP_CALL( SCIPallocBufferArray(scip, &inv_xx, 1));
   SCIP_CALL( SCIPallocBufferArray(scip, &sub_xy, 1));
   SCIP_CALL( SCIPallocBufferArray(scip, &a_old, 1));

   SCIP_CALL( SCIPinitIntArrayZero( p, list));
   dim = 0;
   for( i = 0; i < p; i++ )
   {
      var = var_z[i];

      assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );

      if( (int)SCIPround(scip, SCIPcomputeVarUbLocal(scip, var)) == 0 )   /* then z_i is fixed to 0 */
         list[i] = -1;
      else if( (int)SCIPround(scip, SCIPcomputeVarLbLocal(scip, var)) == 0 ) /* then z_i is not fixed */
         list[i] = 0;
      else /* then z_i is fixed to 1 */
      {
         dim++;
         list[i] = dim;

         if( dim==1 )
         {
            /* set Y */
            if( savememory )
            {
               ni = n * i;
               inv_xx[0] = 1.0 / SCIPcblasDdot(&orig_x[ni], &orig_x[ni], n);
            }
            else
               inv_xx[0] = 1.0 / SCIPmatColMajor(orig_xx, p, i, i);

            a_old[0] = orig_xy[i] * inv_xx[0];

            rss = calcRSSvalue(1, a_old, &orig_xy[i], r);
            objval = calcObjval(n, dim, para_regterm, rss);

            /* update sub_x and sub_xy */
            SCIP_CALL( SCIPcblasCopy(&orig_x[n * i], &sub_x[0], n));
            sub_xy[0] = orig_xy[i];
         }
         else
         {
            ndim = n * dim;
            ni = n * i;
            /* realloc */
            SCIP_CALL( SCIPreallocBufferArray(scip, &a_old, dim));
            SCIP_CALL( SCIPreallocBufferArray(scip, &inv_xx, dim*dim));
            SCIP_CALL( SCIPreallocBufferArray(scip, &sub_xy, dim));

            /* update sub_xy */
            sub_xy[dim-1] = orig_xy[i];

            /* calculate a_old and update inv_xx */
            SCIP_CALL( computeSolution(scip, n, dim, sub_x, &orig_x[ni], inv_xx, sub_xy, a_old, inv_xx));

            rss = calcRSSvalue(dim, a_old, sub_xy, r);

            if( rss < 0.0 )
            {
               /* free */
               SCIPfreeBufferArray(scip, &list);
               SCIPfreeBufferArray(scip, &sub_x);
               SCIPfreeBufferArray(scip, &inv_xx);
               SCIPfreeBufferArray(scip, &sub_xy);
               SCIPfreeBufferArray(scip, &a_old);

               return SCIP_OKAY;
            }

            objval = calcObjval(n, dim, para_regterm, rss);

            /* update sub_x */
            SCIP_CALL( SCIPreallocBufferArray(scip, &sub_x, ndim));
            SCIP_CALL( SCIPcblasCopy(&orig_x[n*i], &sub_x[ndim-n], n));
         }
#if MYPARA_LOG
         printf("---> %dth variable, objval:%f\n", i, objval);
#endif
      }
   }

#if MYPARA_LOG
   printf("list:");
   for( i = 0; i < p; i++ )
   {
      if( i % 10 == 0 )
         SCIPstartNewLine();

      printf("%d, ", list[i]);
   }
   SCIPstartNewLine();
   printf("step2:\n");
#endif

   if( dim == 0 )
   {
      memo = -1;
      rss = 1e+06;
      dim++;
      for( i = 0; i < p; i++ )
      {
         if( list[i] == 0 )
         {
            if( savememory )
            {
               ni = n * i;
               a_old[0] = orig_xy[i] / SCIPcblasDdot(&orig_x[ni], &orig_x[ni], n);
            }
            else
               a_old[0] = orig_xy[i] / SCIPmatColMajor(orig_xx, p, i, i);

            rss_new = calcRSSvalue(1, a_old, &orig_xy[i], r);
            if( rss_new < rss )
            {
               rss = rss_new;
               memo = i;
            }
         }
#if MYPARA_LOG
         printf("%d: rss = %f\n", i, rss_new);
#endif
      }

      assert( memo >= 0 && memo < p );

      objval = calcObjval(n, dim, para_regterm, rss);
      list[memo] = dim;

      /* update sub_x and sub_xy */
      SCIP_CALL( SCIPcblasCopy(&orig_x[n * memo], &sub_x[0], n));
      sub_xy[0] = orig_xy[memo];

      /* set inverse */
      if( savememory )
      {
         ni = n * memo;
         inv_xx[0] = 1.0 / SCIPcblasDdot(&orig_x[ni], &orig_x[ni], n);
      }
      else
         inv_xx[0] = 1.0 / SCIPmatColMajor(orig_xx, p, memo, memo);

#if MYPARA_LOG
   printf("---> %dth variable, objval:%f\n", memo, objval);
#endif
   }

#if MYPARA_LOG
   printf("list:");
   for( i =0; i < p; i++ )
   {
      if( i % 10 == 0 )
         SCIPstartNewLine();
      printf("%d, ", list[i]);
   }
   SCIPstartNewLine();
#endif

   if( objval >= 1e+06 )
   {
      printf("error:heur_forward.c\n");
      SCIPexit();
   }

   while( 1 )
   {
      if( dim == ubnumex )
      {
         assert( dim > 0 && ubnumex > 0 && ubnumex < p && ubnumex < n );
         break;
      }

      dim++;
      memo = -1;
      rss = 1e+06;

      dimdim = dim * dim;
      dim_1 = dim - 1;
      ndim = n * dim;

#if MYPARA_LOG
      printf("(dim=%d) ", dim);
      SCIPprintLongLine();
#endif

      /* alloc */
      SCIP_CALL( SCIPallocBufferArray(scip, &a_new, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &a, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &inv_xx_store, dimdim));
      SCIP_CALL( SCIPallocBufferArray(scip, &inv_xx_new, dimdim));
      SCIP_CALL( SCIPallocBufferArray(scip, &sub_xy_new, dim));
      SCIP_CALL( SCIPcblasCopy(sub_xy, sub_xy_new, dim_1));

      for( i = 0; i < p; i++ )
      {
         if( list[i] == 0 )
         {
            if( ndep == 0 )
               checkld = TRUE;
            else
               checkld = checkFeasiblity(scip, p, ndep, ubnumex, probdata, list, i);

            if( checkld == FALSE )
               continue;

            /* define sub_xy_new */
            sub_xy_new[dim_1] = orig_xy[i];

            /* calculate a and store inv_xx_store for updating inv_xx */
            SCIP_CALL( computeSolution(scip, n, dim, sub_x, &orig_x[n*i], inv_xx, sub_xy_new, a_new, inv_xx_new));

            /* test */
            rss_new = calcRSSvalue(dim, a_new, sub_xy_new, r);
            if( rss_new < rss && rss_new > 0.0 )
            {
               rss = rss_new;
               memo = i;
               SCIP_CALL( SCIPcblasCopy(inv_xx_new, inv_xx_store, dimdim));
               SCIP_CALL( SCIPcblasCopy(a_new, a, dim));
            }
#if MYPARA_LOG
            printf("%d: rss = %f\n", i, rss_new);
#endif
         }
      }

      objval_new = calcObjval(n, dim, para_regterm, rss);
      if( objval_new < objval )
      {
         assert( memo >= 0 && memo < p );

         objval = objval_new;
         list[memo] = dim;

#if MYPARA_LOG
         printf("---> %dth variable, objval:%f\n", memo, objval);
#endif
         /* realloc */
         SCIP_CALL( SCIPreallocBufferArray(scip, &inv_xx, dimdim));
         SCIP_CALL( SCIPreallocBufferArray(scip, &a_old, dim));
         SCIP_CALL( SCIPreallocBufferArray(scip, &sub_x, ndim));
         SCIP_CALL( SCIPreallocBufferArray(scip, &sub_xy, dim));

         /* update inv_xx */
         SCIP_CALL( SCIPcblasCopy(inv_xx_store, inv_xx, dimdim));

         /* update a_old */
         SCIP_CALL( SCIPcblasCopy(a, a_old, dim));

         /* update sub_x */
         SCIP_CALL( SCIPcblasCopy(&orig_x[n*memo], &sub_x[ndim-n], n));

         /* update sub_xy */
         sub_xy[dim_1] = orig_xy[memo];

      }
      else
      {
         memo = -1;
#if MYPARA_LOG
         printf("--> no selection, (objval:%f)\n", objval_new);
#endif
      }

      /* free */
      SCIPfreeBufferArray(scip, &a_new);
      SCIPfreeBufferArray(scip, &a);
      SCIPfreeBufferArray(scip, &inv_xx_store);
      SCIPfreeBufferArray(scip, &inv_xx_new);
      SCIPfreeBufferArray( scip, &sub_xy_new);

      if( memo == -1 )
      {
         dim--;
         break;
      }

   }

   /* check object value of solution */
   nsols = SCIPgetNSols(scip);

   if( nsols < MP_NUM_SOL )
   {
      store = 1;
   }
   else
   {
      sols = SCIPgetSols(scip);
      nsols = MP_NUM_SOL;

      if( objval < SCIPgetSolOrigObj(scip, sols[nsols-1]) )
         store = 1;
      else
         store = 0;
   }

   if( store )
   {

      SCIP_Real* y;
      int nvars;

      SCIP_VAR** var_a;
      SCIP_VAR** var_ep;
      SCIP_VAR* var_rss;
      SCIP_VAR* var_log;

      SCIP_SOL* sol;
      SCIP_Real* solvals;
      SCIP_Bool success;
      SCIP_VAR** vars;
      SCIP_Real* ep;

      int ct;

      nvars = SCIPprobdataGetNvars(probdata);
      y = SCIPprobdataGety(probdata);

      assert(y != NULL);
      assert(nvars == 2 * p + n + 2);

      /* variables */
      var_a = SCIPprobdataGetVars_a(probdata);
      var_ep = SCIPprobdataGetVars_ep(probdata);
      var_rss = SCIPprobdataGetVar_rss(probdata);
      var_log = SCIPprobdataGetVar_log(probdata);

      assert(var_a != NULL);
      assert(var_z != NULL);
      assert(var_ep != NULL);
      assert(var_rss != NULL);
      assert(var_log != NULL);

      /* alloc */
      SCIP_CALL( SCIPallocBufferArray(scip, &ep, n));
      SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars));
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars));

      SCIP_CALL( SCIPcblasDgemv1(sub_x, n, dim, a_old, y, -1.0, 1.0, ep));

      ct = 0;

      /* set solution */
      /* a */
      for( i = 0; i < p; ++i )
      {
         vars[ct] = var_a[i];

         if( list[i] > 0 )
            solvals[ct] = a_old[list[i]-1];
         else
            solvals[ct] = 0.0;

         ct++;
      }

      /* z */
      for( i = 0; i < p; i++ )
      {
         vars[ct] = var_z[i];

         if( list[i] > 0 )
            solvals[ct] = 1.0;
         else
            solvals[ct] = 0.0;

         ct++;
      }

      /* ep */
      for( i = 0; i < n; ++i )
      {
         vars[ct] = var_ep[i];
         solvals[ct] = ep[i];
         ct++;
      }

      vars[ct] = var_rss;
      solvals[ct] = SCIPcblasDdot(ep, ep, n);
      ct++;

      vars[ct] = var_log;
      solvals[ct] = log(solvals[ct-1]);
      ct++;

      assert(nvars == 2 * p + n + 2);
      assert(ct == nvars);

      SCIP_CALL( SCIPcreateSol(scip, &sol, heur));
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, solvals));
      SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, TRUE, FALSE, TRUE, TRUE, &success));

      /* free */
      SCIPfreeBufferArray(scip, &ep);
      SCIPfreeBufferArray(scip, &solvals);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* free */
   SCIPfreeBufferArray(scip, &list);
   SCIPfreeBufferArray(scip, &sub_x);
   SCIPfreeBufferArray(scip, &inv_xx);
   SCIPfreeBufferArray(scip, &sub_xy);
   SCIPfreeBufferArray(scip, &a_old);

   *result = SCIP_FOUNDSOL;
   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the local primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurForward(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Local primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecForward, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyForward) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeForward) );

   return SCIP_OKAY;
}
