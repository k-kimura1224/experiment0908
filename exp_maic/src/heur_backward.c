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

/**@file   heur_backward.c
 * @brief  heuristic for this problem
 * @author Keiji Kimura
 *
 * This file implements heuristic algorithm for this problem. It uses the stepwise method with
 * backward elimination, and find a feasible solution of a given subproblem. The backward elimination
 * involves starting with all variables in the model. This algorithm calculates AIC/BIC of the model
 * that decreased each variable, selects removing the variable that improves the model the most, and
 * repeat this process until none improves the model.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>
#include <string.h>
#include <math.h>

#include "heur_backward.h"
#include "probdata_linereg.h"
#include "convenient_tool.h"
#include "set_myparameter.h"
#include "call_cblas.h"

#define HEUR_NAME             "backward"
#define HEUR_DESC             "primal heuristic using backward elimination"
#define HEUR_DISPCHAR         'b'
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


/** calculate AIC */
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


/** fix variable and check feasiblity by using linear dependence */
static
SCIP_Bool fixVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   p,                  /**< the number of explanatory variables */
   int                   ndep,               /**< the number of linerly dependent sets */
   SCIP_PROBDATA*        probdata,           /**< user problem data */
   int*                  branchinfo          /**< branching information */
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
   int memo;

   assert(scip != NULL);
   assert(p > 0);
   assert(ndep > 0);
   assert(probdata != NULL);
   assert(branchinfo != NULL);

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
      memo = -1;
      for( j = 0; j < size; j++ )
      {
         m = *(ip+j);

         if( branchinfo[m] == 1 )
         {
            memo = -2;
            break;
         }
         else if( branchinfo[p+m] == 1 )
            memo = m;
      }

      assert( memo >= -2 && memo < p );
      if( memo > -1 )
      {
         *(branchinfo + p + memo ) = 0;
         *(branchinfo + memo ) = 1;
      }
      else if( memo == -1 )
      {
         return FALSE;
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
SCIP_DECL_HEURCOPY(heurCopyBackward)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurBackward(scip) );

   return SCIP_OKAY;
}


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeBackward)
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
SCIP_DECL_HEUREXEC(heurExecBackward)
{  /*lint --e{715}*/

   /* probdata */
   SCIP_PROBDATA* probdata;
   int n;
   int p;
   int ndep;
   SCIP_Real* y;
   SCIP_Real* orig_x;
   SCIP_Real* orig_xx;
   SCIP_Real* orig_xy;
   SCIP_Real r;
   int ubnumex;
   int nvars;
   SCIP_Real para_regterm;
   SCIP_Bool savememory;

   int dim;
   int* list;

   SCIP_Real rss;
   SCIP_Real rss_new;
   SCIP_Real objval = 1e+06;
   SCIP_Real objval_new;

   SCIP_Real* a_old;
   SCIP_Real* a;
   SCIP_Real* a_new;
   SCIP_Real* sub_xx;   /* sub matrix of orig_xx */
   SCIP_Real* sub_xy;   /* sub vector of orig_xy */

   /* variables */
   SCIP_VAR** var_a;
   SCIP_VAR** var_z;
   SCIP_VAR** var_ep;
   SCIP_VAR* var_rss;
   SCIP_VAR* var_log;

   /* branching info */
   int ublb;
   int* branchinfo;
   SCIP_Real sum_branchinfo[3];

   /* set solution */
   SCIP_Real* sub_x;
   SCIP_Real* ep;
   int nsols;
   int store;
   SCIP_SOL** sols;
   SCIP_SOL* sol;
   SCIP_Real* solvals;
   SCIP_Bool success;
   SCIP_VAR** vars;

   int i;
   int j;
   int t;
   int ct;
   int ct_a;
   int memo;
   int dpv;

   assert(heur != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);

#if MYPARA_LOG
   printf("backward selection!\n");
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
   ndep = SCIPprobdataGetNdep(probdata);
   nvars = SCIPprobdataGetNvars(probdata);
   ubnumex = SCIPprobdataGetPara_NumEx(probdata);
   para_regterm = SCIPprobdataGetPara_RegTerm(probdata);

   assert(n > 0);
   assert(p > 0);
   assert(ndep >= 0);
   assert(nvars == 2 * p + n + 2);
   assert(EPSEQ(para_regterm, 2.0, 1e-08) || EPSEQ(para_regterm, log((SCIP_Real) n), 1e-08));
   assert(ubnumex == -1);
   assert(n > p);

   if( ubnumex != -1 )
   {
      SCIPerrorMessage("heuristics/backward/freq should be set to -1\n");
      return SCIP_ERROR;
   }

   y = SCIPprobdataGety(probdata);
   orig_x = SCIPprobdataGetx(probdata);
   orig_xx = SCIPprobdataGetxx(probdata);
   orig_xy = SCIPprobdataGetxy(probdata);
   r = SCIPprobdataGetyy(probdata);
   savememory = SCIPprobdataGetSaveMemory(probdata);

   assert(y != NULL);
   assert(orig_x != NULL);
   assert((savememory == FALSE && orig_xx != NULL)
         || (savememory == TRUE && orig_xx == NULL));
   assert(orig_xy != NULL);

   /* variables */
   var_a = SCIPprobdataGetVars_a(probdata);
   var_z = SCIPprobdataGetVars_z(probdata);
   var_ep = SCIPprobdataGetVars_ep(probdata);
   var_rss = SCIPprobdataGetVar_rss(probdata);
   var_log = SCIPprobdataGetVar_log(probdata);

   assert(var_a != NULL);
   assert(var_z != NULL);
   assert(var_ep != NULL);
   assert(var_rss != NULL);
   assert(var_log != NULL);

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

#if MYPARA_LOG
   for( i = 0; i < 3; i++ )
   {
      for( j = 0; j < p; j++ )
      {
         printf("%d, ", *(branchinfo + (i*p) + j));
      }
      SCIPstartNewLine();
   }
#endif

   if( ndep )
   {
      if( fixVariable(scip, p, ndep, probdata, branchinfo) == FALSE )
      {
         SCIPfreeBufferArray(scip, &branchinfo);
         return SCIP_OKAY;
      }
   }

#if MYPARA_LOG
   printf("linear dependence\n");
   if( ndep )
   {
      for( i = 0; i < 3; i++ )
      {
         for( j = 0; j < p; j++ )
         {
            printf("%d, ", *(branchinfo + (i*p) + j));
         }
         SCIPstartNewLine();
      }
   }
#endif

   for( i = 0; i < 3; i++ )
      sum_branchinfo[i] = SCIPcalcIntSum(&branchinfo[i*p], p);

   /* alloc */
   SCIP_CALL( SCIPallocBufferArray(scip, &list, p));

   /* list */
   for( i = 0; i < p; i++ )
      list[i] = 1 - branchinfo[i];

   dim = p - sum_branchinfo[0];
   objval = 1e+06;
   SCIP_CALL( SCIPallocBufferArray(scip, &a_old, dim));

   while( 1 )
   {
      dim--;
      memo = -1;
      rss = 1e+06;

#if MYPARA_LOG
      printf("(dim=%d) ", dim);
      SCIPprintLongLine();
#endif

      /* alloc */
      SCIP_CALL( SCIPallocBufferArray(scip, &a_new, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &sub_xx, dim*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &sub_xy, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &a, dim));

      for( i = 0; i < p; i++ )
      {
         /*
          * solve
          *    Q a_new = sub_xy
          *
          */

         if( *(branchinfo + p + i) == 1 && list[i] == 1 )
         {
            /* generate Q and sub_xy */
            /* Q */
            ct = 0;
            if( savememory )
            {
               for( j = 0; j < p; j++ )
               {
                  if( list[j]==1 && j != i )
                  {
                     for( t = 0; t < p; t++ )
                     {
                        if( list[t]==1 && t != i )
                           sub_xx[ct++] = SCIPcblasDdot(&orig_x[(j*n)], &orig_x[(t*n)], n);
                     }
                  }
               }
            }
            else
            {
               for( j = 0; j < p; j++ )
               {
                  if( list[j]==1 && j != i )
                  {
                     for( t = 0; t < p; t++ )
                     {
                        if( list[t]==1 && t != i )
                           sub_xx[ct++] = SCIPmatColMajor(orig_xx, p, j, t);
                     }
                  }
               }
            }

            assert(ct == dim * dim);

            /* sub_xy */
            ct = 0;
            for( j = 0; j < p; j++ )
            {
               if( list[j]==1 && j != i )
                  sub_xy[ct++] = orig_xy[j];
            }

            assert(ct == dim);

            dpv = SCIPclapackDposv(scip, sub_xx, sub_xy, dim, a_new );

            if( dpv == 0 )
            {
               /* test */
               rss_new = calcRSSvalue(dim, a_new, sub_xy, r);
               if( rss_new < rss ){
                  rss = rss_new;
                  memo = i;
                  SCIP_CALL( SCIPcblasCopy(a_new, a, dim));
               }
#if MYPARA_LOG
               printf("%d: rss = %f\n", i, rss_new);
#endif
            }
         }
      }

      assert(memo >= -1 && memo < p);
      if( memo == -1 )
      {
         /* free */
         SCIPfreeBufferArray(scip, &a_new);
         SCIPfreeBufferArray(scip, &sub_xx);
         SCIPfreeBufferArray(scip, &sub_xy);
         SCIPfreeBufferArray(scip, &a);
         dim++;
         break;
      }

      objval_new = calcObjval(n, dim, para_regterm, rss);
      if( objval_new < objval )
      {
         objval = objval_new;
         list[memo] = 0;
#if MYPARA_LOG
         printf("---> %dth variable, objval:%f\n", memo, objval);
#endif
         /* copy */
         SCIPfreeBufferArray(scip, &a_old);
         SCIP_CALL( SCIPallocBufferArray(scip, &a_old, dim));
         SCIP_CALL( SCIPcblasCopy(a, a_old, dim));
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
      SCIPfreeBufferArray(scip, &sub_xx);
      SCIPfreeBufferArray(scip, &sub_xy);
      SCIPfreeBufferArray(scip, &a);

      if( memo == -1 )
      {
         dim++;
         break;
      }
      else if( SCIPcalcIntSum(list, p) == sum_branchinfo[2] )
      {
         break;
      }
   }

   if( objval >= 1e+06 )
   {
      /* free */
      SCIPfreeBufferArray(scip, &branchinfo);
      SCIPfreeBufferArray(scip, &list);
      SCIPfreeBufferArray(scip, &a_old);
      return SCIP_OKAY;
   }

   nsols = SCIPgetNSols(scip);

   if( nsols < MP_NUM_SOL )
   {
      store = 1;
   }
   else
   {
      sols = SCIPgetSols(scip);
      nsols = MP_NUM_SOL;

      if( objval < SCIPgetSolOrigObj(scip,sols[nsols-1]) )
         store = 1;
      else
         store = 0;
   }

   if( store )
   {
      /* alloc */
      SCIP_CALL( SCIPallocBufferArray(scip, &ep, n));
      SCIP_CALL( SCIPallocBufferArray(scip, &sub_x, n*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars));
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars));

      SCIP_CALL( SCIPcreateSubmatColMajor(orig_x, n, p, list, sub_x));
      SCIP_CALL( SCIPcblasDgemv1( sub_x, n, dim, a_old, y, -1.0, 1.0, ep));

      ct=0;
      ct_a=0;

      /* set solution */
      /* a */
      for( i = 0; i < p; ++i )
      {
         vars[ct] = var_a[i];

         if( list[i] == 1 )
            solvals[ct] = a_old[ct_a++];
         else
            solvals[ct] = 0.0;

         ct++;
      }

      assert(ct_a == dim);

      /* z */
      for( i = 0; i < p; i++ )
      {
         vars[ct] = var_z[i];
         solvals[ct] = (SCIP_Real)list[i];
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
      solvals[ct] = SCIPcblasDdot( ep, ep, n);
      ct++;

      vars[ct] = var_log;
      solvals[ct] = log(solvals[ct-1]);
      ct++;

      assert(nvars == 2 * p + n + 2);
      assert(ct == nvars);

      SCIP_CALL( SCIPcreateSol(scip, &sol, heur));
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, solvals));
      SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, FALSE, FALSE, TRUE, TRUE, &success));

      /* free */
      SCIPfreeBufferArray(scip, &ep);
      SCIPfreeBufferArray(scip, &sub_x);
      SCIPfreeBufferArray(scip, &solvals);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* free */
   SCIPfreeBufferArray(scip, &branchinfo);
   SCIPfreeBufferArray(scip, &list);
   SCIPfreeBufferArray(scip, &a_old);

   *result = SCIP_FOUNDSOL;
   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the local primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurBackward(
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
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecBackward, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyBackward) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeBackward) );

   return SCIP_OKAY;
}
