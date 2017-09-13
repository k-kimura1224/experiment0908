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

/**@file   relax_dposv.c
 * @brief  computing a lower bound
 * @author Keiji Kimura
 *
 * This file implements a method to compute a lower bound
 *
 * Although the relaxation problem of this problem is nonconvex, we can compute the lower bound
 * by solving a linear system. This file calls DPOSV implemented in LAPACK.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <math.h>
#include <assert.h>
#include <string.h>

#include "relax_dposv.h"
#include "probdata_linereg.h"
#include "convenient_tool.h"
#include "call_cblas.h"
#include "get_lineardependence.h"
#include "get_branchvar.h"
#include "set_myparameter.h"

#define RELAX_NAME             "myrelaxator_dposv"
#define RELAX_DESC             "relaxator dposv"
#define RELAX_PRIORITY         1000
#define RELAX_FREQ             1
#define RELAX_INCLUDESLP       TRUE    /* What does this parameter mean? */

#define RELAX_TIMELOG          0

#if RELAX_TIMELOG
#include <time.h>
#endif

/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   int test;
};

/*
 * Local methods
 */

/** solve subproblem */
static
SCIP_RETCODE solveSubproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< user problem data */
   int                   n,                  /**< the number of data pooints */
   int                   p,                  /**< the number of explantory variables */
   int                   ubnumex,            /**< upper bound of the number of explantory variables */
   SCIP_Real             para_regterm,       /**< parameter to control the regularization term */
   SCIP_VAR**            var_z               /**< variable z */
   )
{
   /* probdata */
   SCIP_Real* data_x;
   SCIP_Real* data_y;
   SCIP_Real* orig_xx;
   SCIP_Real* orig_xy;
   SCIP_Real yy;
   SCIP_Bool savememory;

   int* index;
   int dim;
   int dimdim;
   SCIP_Real* sub_xx;
   SCIP_Real* sub_xy;
   SCIP_Real* sub_sol;
   SCIP_Real rss;
   SCIP_Real log_rss;
   SCIP_Real dualbound;
   int info;

   /* for solution */
   int nsols;
   int maxsols = MP_NUM_SOL;
   int store;

   int ct;
   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(probdata != NULL);
   assert(n > 0);
   assert(p > 0);
   assert(ubnumex >= 1 && ubnumex <= p);
   assert(EPSEQ(para_regterm, 2.0, 1.0e-08) || EPSEQ(para_regterm, log((SCIP_Real) n), 1.0e-08));
   assert(var_z != NULL);
   assert(maxsols > 0);

   /* get values from the probdata */
   data_x = SCIPprobdataGetx(probdata);
   data_y = SCIPprobdataGety(probdata);
   orig_xy = SCIPprobdataGetxy(probdata);
   yy = SCIPprobdataGetyy(probdata);
   orig_xx = SCIPprobdataGetxx(probdata);
   savememory = SCIPprobdataGetSaveMemory(probdata);

   assert(data_x != NULL);
   assert(data_y != NULL);
   assert((savememory == FALSE && orig_xx != NULL)
         || (savememory == TRUE && orig_xx == NULL));
   assert(orig_xy != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &index, ubnumex));

   /* set index and dimension */
   dim = 0;
   for( i = 0; i < p; ++i )
   {
      assert(SCIPvarGetType(var_z[i]) == SCIP_VARTYPE_BINARY);
      if( SCIPround(scip, SCIPcomputeVarLbLocal(scip, var_z[i])) == 1 )
      {
         index[dim] = i;
         dim++;

         if( dim == ubnumex )
            break;
      }
   }

   assert(dim == ubnumex);
   assert(dim < n);

   dimdim = dim * dim;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &sub_xx, dimdim));
   SCIP_CALL( SCIPallocBufferArray(scip, &sub_xy, dim));
   SCIP_CALL( SCIPallocBufferArray(scip, &sub_sol, dim));

   if( savememory )
   {
      int* nindex;
      SCIP_Real* dp1;
      SCIP_Real* dp2;

      SCIP_CALL( SCIPallocBufferArray(scip, &nindex, dim));
      for( i = 0; i < dim; i++ )
         nindex[i] = n * index[i];


      for( i = 0; i < dim; i++ )
      {
         k = nindex[i];
         dp1 = &sub_xx[i*dim];
         for( j = i; j < dim; j++ )
            *(dp1 + j) = SCIPcblasDdot(&data_x[k], &data_x[nindex[j]], n);
      }

      for( i = 0; i < dim; i++ )
      {
         dp1 = &sub_xx[(i*dim)];
         dp2 = &sub_xx[i];
         for( j = 0; j < i; j++ )
            *(dp1 + j) = *(dp2 + j * dim);
      }

      SCIPfreeBufferArray(scip, &nindex);
   }
   else
   {
      ct = 0;
      for( i = 0; i < dim; i++ )
      {
         k = index[i] * p;
         for( j = 0; j < dim; j++ )
         {
            sub_xx[ct++] = orig_xx[k + index[j]];
         }
      }
   }

   for( i = 0; i < dim; i++ )
      sub_xy[i] = orig_xy[index[i]];

   /* solve linear system */
   info = SCIPclapackDposv(scip, sub_xx, sub_xy, dim, sub_sol);
   assert( info == 0 );
   if( info != 0 )
   {
      SCIPerrorMessage("info:%d\n", info);
      SCIPexit();
   }

   /* calculate dualbound */
   rss = yy - SCIPcblasDdot(sub_xy, sub_sol, dim);
   assert( rss > 1e-08 );
   log_rss = log(rss);
   dualbound = (SCIP_Real) n * log_rss;
   dualbound += para_regterm * (SCIP_Real) dim;

   /* update lower bound */
   SCIP_CALL( SCIPupdateLocalLowerbound(scip, dualbound));

   /* get the number of stored solutions */
   nsols = SCIPgetNSols(scip);
   if( nsols < maxsols )
   {
     store = 1;
   }
   else
   {
      SCIP_Real primalbound;
      SCIP_SOL** sols;

      assert( nsols == maxsols );

      primalbound = dualbound;
      sols = SCIPgetSols(scip);

      if( primalbound < SCIPgetSolOrigObj(scip, sols[nsols-1]) )
         store = 1;
      else
         store = 0;
   }

   if( store )
   {
      SCIP_SOL* sol;
      SCIP_Real* vals_ep;
      SCIP_Real* sub_data;
      SCIP_Real* solvals;
      SCIP_HEUR* heur;
      SCIP_Bool success;
      SCIP_VAR** vars;

      /* probdata */
      int nvars;
      SCIP_VAR** var_a;
      SCIP_VAR** var_ep;
      SCIP_VAR* var_rss;
      SCIP_VAR* var_log;

      nvars = SCIPprobdataGetNvars(probdata);
      var_a = SCIPprobdataGetVars_a(probdata);
      var_ep = SCIPprobdataGetVars_ep(probdata);
      var_rss = SCIPprobdataGetVar_rss(probdata);
      var_log = SCIPprobdataGetVar_log(probdata);

      assert(nvars == 2 * p + n + 2);
      assert(var_a != NULL);
      assert(var_ep != NULL);
      assert(var_rss != NULL);
      assert(var_log != NULL);

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &vals_ep, n));
      SCIP_CALL( SCIPallocBufferArray(scip, &sub_data, n*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars));
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars));

      /* set values of solution */
      ct = 0;
      for( i = 0; i < p; i++ )
      {
         vars[ct] = var_a[i];
         solvals[ct] = 0.0;
         ct++;
      }

      for( i = 0; i < dim; i++ )
         solvals[index[i]] = sub_sol[i];

      /* set sub-matrix */
      for( i = 0; i < dim; i++ )
         SCIP_CALL( SCIPcblasCopy(&data_x[n*index[i]], &sub_data[n*i], n));

      /* ep = y - Ax */
      SCIP_CALL( SCIPcblasDgemv1(sub_data, n, dim, sub_sol, data_y, -1.0, 1.0, vals_ep));

      for( i = 0; i < p; i++ )
      {
         vars[ct] = var_z[i];
         solvals[ct] = 0.0;
         ct++;
      }

      for( i = 0; i < dim; i++ )
         solvals[p+index[i]] = 1.0;

      for( i = 0; i < n; i++ )
      {
         vars[ct] = var_ep[i];
         solvals[ct] = vals_ep[i];
         ct++;
      }

      vars[ct] = var_rss;
      solvals[ct] = rss;
      ct++;

      vars[ct] = var_log;
      solvals[ct] = log_rss;
      ct++;

      assert(ct == nvars);

      heur = SCIPfindHeur(scip, "trysol");
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur));
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, solvals));
      SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, FALSE, FALSE, TRUE, TRUE, &success));

      /*
      if( success == FALSE )
      {
         SCIPprintIntVal(dim);
         assert(0);
      }
      */

      /* free */
      SCIPfreeBufferArray(scip, &vals_ep);
      SCIPfreeBufferArray(scip, &sub_data);
      SCIPfreeBufferArray(scip, &solvals);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* free */
   SCIPfreeBufferArray(scip, &index);
   SCIPfreeBufferArray(scip, &sub_xx);
   SCIPfreeBufferArray(scip, &sub_xy);
   SCIPfreeBufferArray(scip, &sub_sol);

   return SCIP_OKAY;
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


/** calculate sum_branchinfo */
static
SCIP_RETCODE calcSumBranchinfo(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   p,                  /**< the number of explantory variables */
   int                   pp,                 /**< p * p */
   int*                  branchinfo,         /**< branching information */
   int*                  sum_branchinfo
   )
{
#if 1
   int* ip0 = &branchinfo[0];
   int* ip1 = &branchinfo[pp];
   int i;

   assert(scip != NULL);
   assert(p > 0);
   assert(pp == 2 * p);
   assert(branchinfo != NULL);
   assert(sum_branchinfo != NULL);

   sum_branchinfo[0] = 0;
   sum_branchinfo[2] = 0;

   for( i = 0; i < p; i++ )
   {
      sum_branchinfo[0] += *(ip0 + i);
      sum_branchinfo[2] += *(ip1 + i);
   }
   sum_branchinfo[1] = p - sum_branchinfo[0] - sum_branchinfo[2];

#else /* These are the same result */
   int i;

   assert(scip != NULL);
   assert(p > 0);
   assert(p == p * p);
   assert(branchinfo != NULL);
   assert(sum_branchinfo != NULL);

   for( i = 0; i < 3; i++ )
      sum_branchinfo[i] = SCIPcalcIntSum(&branchinfo[i*p], p);
#endif

   return SCIP_OKAY;
}

/*
 * Callback methods of relaxator
 */

/** copy method for relaxator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_RELAXCOPY(relaxCopyDposv)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(relax != NULL);
   assert(strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0);

   /* call inclusion method of relaxator */
   SCIP_CALL( SCIPincludeRelaxDposv(scip));

   return SCIP_OKAY;
}


/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreeDposv)
{   /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;

   assert(relax != NULL);
   assert(strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0);
   assert(scip != NULL);

   /* free relaxator data */
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   SCIPfreeMemory(scip, &relaxdata);
   SCIPrelaxSetData(relax, NULL);

   return SCIP_OKAY;
}


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecDposv)
{
   /* for probdata */
   SCIP_PROBDATA* probdata;
   int n;
   int p;
   int ubnumex;
   SCIP_Real para_regterm;
   int ndep;
   SCIP_Real* y;
   SCIP_Real* data_x;
   SCIP_Real* orig_xx;
   SCIP_Real* orig_xy;
   SCIP_Real yy;
   SCIP_VAR** var_z;
   SCIP_Bool savememory;

   SCIP_NODE* node;
   SCIP_VAR* lastbranchvar;
   int *branchinfo;                 /* arraty for branching information */
   int sum_branchinfo[3];
   int dim;
   int dimdim;
   int* index;

   /* for while */
   SCIP_Real* sub_xx;
   SCIP_Real* sub_xy;
   SCIP_Real* sub_sol;
   SCIP_Real* sub_x;
   int* ldindex;
   int* sub_maxdep;
   int* sub_depsets;
   int sub_ndep;
   int check;
   int info;

   /* for updating local dual bound */
   SCIP_Real rss;                   /* residual sum of square */
   SCIP_Real log_rss;               /* log ( rss ) */
   SCIP_Real dualbound;             /* local dual bound */

   /* for setting the primal solution */
   int store;
   int nsols;

   int i;
   int j;
   int k;
   int pp;
   int ct;

   /* pointers */
   int* ip;
   int* ip0;
   int* ip1;
   int* ip2;

#if RELAX_TIMELOG
   clock_t start;
   clock_t end;
   SCIP_Real runtime = 0.0;
   printf("start [%d] ---------------------------------\n", __LINE__);
   start = clock();
#endif

   assert(relax != NULL);
   assert(scip != NULL);
   assert(strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0);
   assert(result != NULL);

   /* get relaxator data */
   /*
   SCIP_RELAXDATA* relaxdata;
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   */

   /* get the probdata */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* get values from the probdata */
   n = SCIPprobdataGetNdatas(probdata);
   p = SCIPprobdataGetNexvars(probdata);
   ubnumex = SCIPprobdataGetPara_NumEx(probdata);
   para_regterm = SCIPprobdataGetPara_RegTerm(probdata);
   var_z = SCIPprobdataGetVars_z(probdata);
   pp = 2 * p;

   assert(n > 0);
   assert(p > 0);
   assert(ubnumex == -1 || (ubnumex > 0 && ubnumex < p) );
   assert(EPSEQ(para_regterm, 2.0, 1e-08) || EPSEQ(para_regterm, log((SCIP_Real) n), 1e-08));
   assert(var_z != NULL);

   /* get current node */
   node = SCIPgetCurrentNode(scip);

   /* get the last branching variable */
   SCIP_CALL( SCIPvarGetLastBranching(scip, node, &lastbranchvar));
   assert( lastbranchvar == NULL || SCIPvarGetType(lastbranchvar) == SCIP_VARTYPE_BINARY );

   /* this means that lastbranchvar is fixed to 1 */
   if( lastbranchvar != NULL &&
         (int)SCIPround(scip, SCIPcomputeVarLbLocal(scip, lastbranchvar)) == 1 )
   {
      SCIP_NODE* parent;
      SCIP_Real parentdual;

      assert( (int)SCIPround(scip, SCIPcomputeVarUbLocal(scip, lastbranchvar)) == 1);

      if( ubnumex != -1 )
      {
         int numfix1 = 0;

         assert(ubnumex > 0 && ubnumex < p);

         /* count the number of fixed binary variables z = 1 */
         for( i = 0; i < p && numfix1 < ubnumex; ++i )
            numfix1 += (int)SCIPround(scip, SCIPcomputeVarLbLocal(scip, var_z[i]));

         assert(numfix1 >= 1 && numfix1 <= ubnumex);
         if( numfix1 == ubnumex )
         {
            if( p >= n )
            {
               /* update the lower bound and find a feasible solution */
               SCIP_CALL( solveSubproblem(scip, probdata, n, p, ubnumex, para_regterm, var_z));
            }
            else
            {
               /* get parent node */
               parent = SCIPnodeGetParent(node);
               /* get local lower bound */
               parentdual = SCIPgetNodeLowerbound(scip, parent);
               /* update */
               SCIP_CALL( SCIPupdateLocalLowerbound(scip, parentdual + para_regterm));
            }

#if RELAX_TIMELOG
            end = clock();
            runtime = SCIPcalcTime(start, end) - runtime;
            printf("|[%d] ->%f\n", __LINE__, runtime);
            printf("|total->%f\n", SCIPcalcTime(start, end));
            printf("end [%d] -----------------------------------\n\n", __LINE__);
#endif
            /* cut-off */
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
      }

      /* get parent node */
      parent = SCIPnodeGetParent(node);
      /* get local lower bound */
      parentdual = SCIPgetNodeLowerbound(scip, parent);
      /* update */
      SCIP_CALL( SCIPupdateLocalLowerbound(scip, parentdual + para_regterm));

#if RELAX_TIMELOG
      end = clock();
      runtime = SCIPcalcTime(start, end) - runtime;
      printf("|[%d] ->%f\n", __LINE__, runtime);
      printf("|total->%f\n", SCIPcalcTime(start, end));
      printf("end [%d] -----------------------------------\n\n", __LINE__);
#endif

      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* allocate memory for branchinfo */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchinfo, 3*p));

   /* get branching information */
   /* branchinfo[i] = 1 if z_i is fixed to 0 */
   ip0 = &branchinfo[0];
   for( i = 0; i < p; i++ )
   {
      *(ip0 + i) = 1 - SCIPround(scip, SCIPcomputeVarUbLocal(scip, var_z[i]));
      assert( *(ip0 + i) == 1 || *(ip0 + i) == 0 );
   }

   /* branchinfo[2p+i] = 1 if z_i is fixed to 1 */
   ip1 = &branchinfo[pp];
   for( i = 0; i < p; i++ )
   {
      *(ip1 + i) = SCIPround(scip, SCIPcomputeVarLbLocal(scip, var_z[i]));
      assert( *(ip1 + i) == 1 || *(ip1 + i) == 0 );
   }

   /* branchinfo[p+i] = 1 if z_i is not fixed */
   ip = &branchinfo[p];
   for( i = 0; i < p; i++ )
   {
      *(ip + i) = 1 - *(ip0 + i) - *(ip1 + i);
      assert( *(ip + i) == 1 || *(ip + i) == 0 );
   }

   /* get values from the probdata */
   ndep = SCIPprobdataGetNdep(probdata);
   assert(ndep >= 0);

   if( ndep )
   {
      /* fix non-fixed variable z and check feasiblity by using linear dependence */
      if( fixVariable(scip, p, ndep, probdata, branchinfo) == FALSE )
      {
         /* cut-off */
         SCIPfreeBufferArray(scip, &branchinfo);
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   /* calculate sub_branchinfo */
   SCIP_CALL( calcSumBranchinfo(scip, p, pp, branchinfo, sum_branchinfo));

   dim = sum_branchinfo[1] + sum_branchinfo[2];
   assert(dim > 0);

   if( dim >= n )
   {
#if 0
      ip = &branchinfo[p];
      for( i = 0; i < p; i++ )
      {
         if( *(ip + i) == 1 )
         {
            branchinfo[i] = 1;
            *(ip + i) = 0;

            sum_branchinfo[0]++;
            sum_branchinfo[1]--;
            dim--;

            if( dim == n - 1 )
               break;
         }
      }
      assert(dim == n - 1);
#else
      SCIPfreeBufferArray(scip, &branchinfo);
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
#endif
   }


   /* get values from the probdata */
   y = SCIPprobdataGety(probdata);
   data_x = SCIPprobdataGetx(probdata);
   orig_xx = SCIPprobdataGetxx(probdata);
   orig_xy = SCIPprobdataGetxy(probdata);
   yy = SCIPprobdataGetyy(probdata);
   savememory = SCIPprobdataGetSaveMemory(probdata);

   assert(y != NULL);
   assert(data_x != NULL);
   assert((savememory == FALSE && orig_xx != NULL)
         || (savememory == TRUE && orig_xx == NULL));
   assert(orig_xy != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &index, dim));

   while( 1 )
   {
      assert(dim > 0);
      dimdim = dim * dim;

      ct = 0;
      for( i = 0; i < p; i++ )
      {
         if( branchinfo[i] == 0 )
         {
            index[ct++] = i;

            if( dim == ct )
               break;
         }
      }

      assert( dim == ct );

      /* alloc */
      SCIP_CALL( SCIPallocBufferArray(scip, &sub_xx, dimdim));
      SCIP_CALL( SCIPallocBufferArray(scip, &sub_xy, dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &sub_sol, dim));

      if( savememory )
      {
         int* nindex;
         SCIP_Real* dp1;
         SCIP_Real* dp2;

         SCIP_CALL( SCIPallocBufferArray(scip, &nindex, dim));
         for( i = 0; i < dim; i++ )
            nindex[i] = n * index[i];


         for( i = 0; i < dim; i++ )
         {
            k = nindex[i];
            dp1 = &sub_xx[i*dim];
            for( j = i; j < dim; j++ )
               *(dp1 + j) = SCIPcblasDdot(&data_x[k], &data_x[nindex[j]], n);
         }

         for( i = 0; i < dim; i++ )
         {
            dp1 = &sub_xx[(i*dim)];
            dp2 = &sub_xx[i];

            for( j = 0; j < i; j++ )
               *(dp1 + j) = *(dp2 + j * dim);
         }

         SCIPfreeBufferArray(scip, &nindex);
      }
      else
      {
         ct = 0;
         for( i = 0; i < dim; i++ )
         {
            k = index[i] * p;
            for( j = 0; j < dim; j++ )
               sub_xx[ct++] = orig_xx[k + index[j]];
         }
      }

      for( i = 0; i < dim; i++ )
         sub_xy[i] = orig_xy[index[i]];

      /* solve with dposv of clapack */
      info = SCIPclapackDposv(scip, sub_xx, sub_xy, dim, sub_sol);

      if( info==0 )
         break;   /* solving is successful */

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &sub_x, n * dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &ldindex, dim));

      /* set sub-matrix */
      for( i = 0; i < dim; i++ )
         SCIP_CALL( SCIPcblasCopy(&data_x[n*index[i]], &sub_x[n*i], n));

      SCIP_CALL( SCIPgetNLineDependSet(scip, sub_x, n, dim, ldindex));
      sub_ndep = SCIPcalcIntSum(ldindex, dim);

      if( sub_ndep == 0 )
      {
         SCIPfreeBufferArray(scip, &sub_x);
         SCIPfreeBufferArray(scip, &ldindex);

         SCIPfreeBufferArray(scip, &sub_xx);
         SCIPfreeBufferArray(scip, &sub_xy);
         SCIPfreeBufferArray(scip, &sub_sol);

         SCIPfreeBufferArray(scip, &branchinfo);
         SCIPfreeBufferArray(scip, &index);

         *result = SCIP_SUCCESS;
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &sub_maxdep, sub_ndep));
      SCIP_CALL( SCIPallocBufferArray(scip, &sub_depsets, sub_ndep*dim));
      SCIP_CALL( SCIPinitIntArrayZero(sub_ndep*dim, sub_depsets));

      ct = 0;
      for( i = 0; i < dim; i++ )
      {
         if( ldindex[i] == 1 )
            sub_maxdep[ct++] = i;
      }

      assert(sub_ndep == ct);

      for( i = 0; i < sub_ndep; i++ )
         *(sub_depsets + (i * dim) + sub_maxdep[i]) = 1;

      SCIP_CALL( SCIPgetLineDependSet(scip, sub_xx, dim, sub_ndep, sub_maxdep, ldindex, sub_depsets));

      check = 0;

      for( i = 0; i < sub_ndep; i++ )
      {
         ip = &branchinfo[p];
         ip2 = &sub_depsets[i*dim];
         for( j = 0; j < dim; j++ )
         {
            if( *(ip2 + j) == 1 )
            {
               if( *(branchinfo + index[j]) == 1 )
                  break;

               if( *(ip + index[j]) == 1 )
               {
                  *(ip + index[j]) = 0;
                  *(branchinfo + index[j]) = 1;
                  check = 1;
                  break;
               }
            }
         }
      }

      /* free */
      SCIPfreeBufferArray(scip, &sub_x);
      SCIPfreeBufferArray(scip, &ldindex);
      SCIPfreeBufferArray(scip, &sub_maxdep);
      SCIPfreeBufferArray(scip, &sub_depsets);

      SCIPfreeBufferArray(scip, &sub_xx);
      SCIPfreeBufferArray(scip, &sub_xy);
      SCIPfreeBufferArray(scip, &sub_sol);

      if( !check )
      {
         /* free */
         SCIPfreeBufferArray(scip, &branchinfo);
         SCIPfreeBufferArray(scip, &index);

#if RELAX_TIMELOG
         end = clock();
         runtime = SCIPcalcTime(start, end) - runtime;
         printf("|[%d] ->%f\n", __LINE__, runtime);
         printf("|total->%f\n", SCIPcalcTime(start, end));
         printf("end [%d] -----------------------------------\n\n", __LINE__);
#endif

         *result=SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* calculate sub_branchinfo */
      SCIP_CALL( calcSumBranchinfo(scip, p, pp, branchinfo, sum_branchinfo));

      dim = sum_branchinfo[1] + sum_branchinfo[2];

   }

   /* calculate dualbound */
   rss = yy - SCIPcblasDdot(sub_xy, sub_sol, dim);


   if( rss < 0.0 || EPSEQ(rss, 0.0, 1.0e-05) )
   {
      SCIP_Real* alphaep;
      SCIP_Real alpha = 1e+8;
      SCIP_Real* ep;
      SCIP_Real* x;

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &ep, n));
      SCIP_CALL( SCIPallocBufferArray(scip, &x, n*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &alphaep, n));

      /* set sub-matrix */
      for( i = 0; i < dim; i++ )
         SCIP_CALL( SCIPcblasCopy(&data_x[n*index[i]], &x[n*i], n));

      /* ep = y - Ax */
      SCIP_CALL( SCIPcblasDgemv1( x, n, dim, sub_sol, y, -1.0, 1.0, ep));

      /* alphaep := alpah ep */
      SCIP_CALL( SCIPcblasDscal(ep, n, alpha, alphaep));

      rss = SCIPcblasDdot(alphaep, alphaep, n);

      assert( rss > 1.0e-12 );

      log_rss = log(rss);
      log_rss -= log(alpha);
      log_rss -= log(alpha);

      dualbound = (SCIP_Real) n * log_rss;
      dualbound += para_regterm * (SCIP_Real)sum_branchinfo[2];

      /* update the lower bound */
      SCIP_CALL( SCIPupdateLocalLowerbound(scip, dualbound));

      /* free */
      SCIPfreeBufferArray(scip, &branchinfo);
      SCIPfreeBufferArray(scip, &index);
      SCIPfreeBufferArray(scip, &sub_xx);
      SCIPfreeBufferArray(scip, &sub_xy);
      SCIPfreeBufferArray(scip, &sub_sol);

      SCIPfreeBufferArray(scip, &ep);
      SCIPfreeBufferArray(scip, &x);
      SCIPfreeBufferArray(scip, &alphaep);

#if RELAX_TIMELOG
      end = clock();
      runtime = SCIPcalcTime(start, end) - runtime;
      printf("|[%d] ->%f\n", __LINE__, runtime);
      printf("|total->%f\n", SCIPcalcTime(start, end));
      printf("end [%d] -----------------------------------\n\n", __LINE__);
#endif
      *result=SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   log_rss = log(rss);
   dualbound = (SCIP_Real) n * log_rss;
   dualbound += para_regterm * (SCIP_Real)sum_branchinfo[2];

   /* update the lower bound */
   SCIP_CALL( SCIPupdateLocalLowerbound(scip, dualbound));


   /* get the number of stored solutions */
   nsols = SCIPgetNSols(scip);

   /* get the number of stored solutions */
   nsols = SCIPgetNSols(scip);

   assert( ubnumex >= -1 && ubnumex < p );
   if( ubnumex == -1 )
   {
      if( nsols < MP_NUM_SOL )
      {
         store = 1;
      }
      else
      {
         SCIP_SOL** sols;
         SCIP_Real primal;

         sols = SCIPgetSols(scip);
         primal = dualbound;
         primal += para_regterm * (SCIP_Real) sum_branchinfo[1];

         if( primal < SCIPgetSolOrigObj(scip, sols[nsols-1]) )
            store = 1;
         else
            store = 0;
      }
   }
   else
   {
      if( dim > ubnumex )
         store = 0;
      else
      {
         if( nsols < MP_NUM_SOL )
         {
            store = 1;
         }
         else
         {
            SCIP_SOL** sols;
            SCIP_Real primal;

            sols = SCIPgetSols(scip);
            primal = dualbound;
            primal += para_regterm * (SCIP_Real) sum_branchinfo[1];

            if( primal < SCIPgetSolOrigObj(scip, sols[nsols-1]) )
               store = 1;
            else
               store = 0;
         }
      }
   }

   if( store )
   {
      SCIP_SOL* sol;
      SCIP_Real* ep;
      SCIP_Real* x;
      SCIP_Real* solvals;
      SCIP_HEUR* heur;
      SCIP_Bool success;

      /* probdata */
      int nvars;
      SCIP_VAR** vars;
      SCIP_VAR** var_a;
      SCIP_VAR** var_ep;
      SCIP_VAR* var_rss;
      SCIP_VAR* var_log;

      nvars = SCIPprobdataGetNvars(probdata);
      var_a = SCIPprobdataGetVars_a(probdata);
      var_ep = SCIPprobdataGetVars_ep(probdata);
      var_rss = SCIPprobdataGetVar_rss(probdata);
      var_log = SCIPprobdataGetVar_log(probdata);

      assert(nvars == (2 * p) + n + 2);
      assert(var_a != NULL);
      assert(var_ep != NULL);
      assert(var_rss != NULL);
      assert(var_log != NULL);

      /* alloc */
      SCIP_CALL( SCIPallocBufferArray(scip, &ep, n));
      SCIP_CALL( SCIPallocBufferArray(scip, &x, n*dim));
      SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars));
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars));

      /* set values of solution */
      ct = 0;
      for( i = 0; i < p; i++ )
      {
         vars[ct] = var_a[i];
         solvals[ct] = 0.0;
         ct++;
      }

      for( i = 0; i < dim; i++ )
         solvals[index[i]] = sub_sol[i];

      /* set sub-matrix */
      for( i = 0; i < dim; i++ )
         SCIP_CALL( SCIPcblasCopy(&data_x[n*index[i]], &x[n*i], n));

      /* ep = y - Ax */
      SCIP_CALL( SCIPcblasDgemv1(x, n, dim, sub_sol, y, -1.0, 1.0, ep));

      for( i = 0; i < p; i++ )
      {
         vars[ct] = var_z[i];
         solvals[ct] = 0.0;
         ct++;
      }

      for( i = 0; i < dim; i++ )
         solvals[p+index[i]] = 1.0;

      for( i = 0; i < n; i++ )
      {
         vars[ct] = var_ep[i];
         solvals[ct] = ep[i];
         ct++;
      }

      vars[ct] = var_rss;
      solvals[ct] = rss;
      ct++;

      vars[ct] = var_log;
      solvals[ct] = log_rss;
      ct++;

      assert(ct == nvars);

      heur = SCIPfindHeur(scip, "trysol");
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur));
      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, solvals));
      SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, FALSE, FALSE, TRUE, TRUE, &success));

      /* free */
      SCIPfreeBufferArray(scip, &ep);
      SCIPfreeBufferArray(scip, &x);
      SCIPfreeBufferArray(scip, &solvals);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* free */
   SCIPfreeBufferArray(scip, &branchinfo);
   SCIPfreeBufferArray(scip, &index);
   SCIPfreeBufferArray(scip, &sub_xx);
   SCIPfreeBufferArray(scip, &sub_xy);
   SCIPfreeBufferArray(scip, &sub_sol);

#if RELAX_TIMELOG
   end = clock();
   runtime = SCIPcalcTime(start, end) - runtime;
   printf("|[%d] ->%f\n", __LINE__, runtime);
   printf("|total ->%f\n", SCIPcalcTime(start, end));
   printf("end [%d] -----------------------------------\n\n", __LINE__);
#endif

   *result=SCIP_SUCCESS;
   return SCIP_OKAY;
}


/** create this relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxDposv(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;

   /* create myrelaxator data */
   SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );


   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, RELAX_INCLUDESLP,
         relaxExecDposv, relaxdata) );

   assert(relax != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetRelaxCopy(scip, relax, relaxCopyDposv) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeDposv) );


   return SCIP_OKAY;
}
