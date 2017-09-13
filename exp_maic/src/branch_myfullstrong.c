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

/**@file   branch_myfullstrong.c
 * @brief  full strong branching
 * @author Keiji Kimura
 *
 * This file implements full strong branching for minimization of AIC/BIC
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>

#include "branch_myfullstrong.h"
#include "probdata_linereg.h"
#include "set_myparameter.h"
#include "convenient_tool.h"
#include "call_cblas.h"

#define BRANCHRULE_NAME            "fullstrong_maic"
#define BRANCHRULE_DESC            "my full strong branching"
#if 0
#define BRANCHRULE_PRIORITY        200000
#else
#define BRANCHRULE_PRIORITY        100000
#endif
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define MYPARA_LOG   0

/*
 * Data structures
 */
/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
   int a;
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


/** check whether variable z_k can be fixed to 1 */
static
SCIP_Bool checkBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   p,                  /**< the number of explanatory variables */
   int                   ndep,               /**< the number of linerly dependent sets */
   SCIP_PROBDATA*        probdata,           /**< user problem data */
   int*                  branchinfo,         /**< branching information */
   int                   k                   /**< index of variable */
   )
{
   SCIP_Bool result = TRUE;

   int* indexdepsets;
   int* sizedepsets;

   int i;
   int j;
   int m;
   int ct;
   int size;
   int buf;
   int* ip;
   int check;

   assert(scip != NULL);
   assert(p > 0);
   assert(ndep >= 0);
   assert(probdata != NULL);
   assert(branchinfo != NULL);
   assert(k >= 0 && k < p);

   if( ndep == 0 )
      return TRUE;

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
      ct = 0;
      check = 0;
      for( j = 0; j < size; j++ )
      {
         m = *(ip+j);

         if( m == k )
            check = 1;
         else if( m > k && check == 0 )
            break;

         if( branchinfo[m] == 1 )
         {
            ct = 0;
            break;
         }
         else if( branchinfo[p+m] == 1 )
         {
            ct++;
            if( ct == 2 )
               break;
         }

      }
      assert( ct == 0 || ct == 1 || ct == 2 );
      if( ct == 1 && check == 1 )
      {
         result = FALSE;
         break;
      }
      buf = size;
   }

   return result;
}


/** fix variable by using linear dependence */
static
SCIP_RETCODE fixVariable(
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
            break;
         else if( branchinfo[p+m] == 1 )
            memo = m;
      }

      assert( memo >= -1 && memo < p );
      if( memo != -1 )
      {
         *(branchinfo + p + memo ) = 0;
         *(branchinfo + memo ) = 1;
      }

      buf = size;
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of branching rule
 */

/* TODO: Implement all necessary branching rule methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyMyfullstrong)
{  /*lint --e{715}*/

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleMyfullstrong(scip) ) ;

   return SCIP_OKAY;
}


/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeMyfullstrong)
{  /*lint --e{715}*/

   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitMyfullstrong)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
#if 1
static
SCIP_DECL_BRANCHEXIT(branchExitMyfullstrong)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}
#else
#define branchExitMyfullstrong NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BRANCHINITSOL(branchInitsolMyfullstrong)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myfullstrong branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolMyfullstrong NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolMyfullstrong)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myfullstrong branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolMyfullstrong NULL
#endif


/** branching execution method for fractional LP solutions */
#if 0
static
SCIP_DECL_BRANCHEXECLP(branchExeclpMyfullstrong)
{  /*lint --e{715}*/

   printf("mybranching rule!\n");
   return SCIP_OKAY;
}
#else
#define branchExeclpMyfullstrong NULL
#endif


/** branching execution method for external candidates */
#if 0
static
SCIP_DECL_BRANCHEXECEXT(branchExecextMyfullstrong)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myfullstrong branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextMyfullstrong NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#if 1
static
SCIP_DECL_BRANCHEXECPS(branchExecpsMyfullstrong)
{  /*lint --e{715}*/

   SCIP_VAR** cands;
   int ncands;
   SCIP_NODE* childnode_0;
   SCIP_NODE* childnode_1;

   /* probdata */
   SCIP_PROBDATA* probdata;
   int n;
   int p;
   int ndep;
   SCIP_VAR** var_z;
   SCIP_Real* orig_xx;
   SCIP_Real* orig_xy;
   SCIP_Real r;

   int dim;
   SCIP_Real rss;
   SCIP_Real rss_new;
   SCIP_Real* a;
   SCIP_Real* sub_xx;
   SCIP_Real* sub_xy;

   int ublb;
   int* branchinfo;
   int* branchinfo_new;
   int* list;

   int i;
   int j;
   int t;
   int ct;
   int ind;
   int dpv;
   int p3;

#if MYPARA_LOG
   printf("[myfullstrong brnaching]");
   SCIPprintLongLine();
#endif

   /* get branching rule data */
   /*
   SCIP_BRANCHRULEDATA* branchruledata;
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   */

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   n = SCIPprobdataGetNdatas(probdata);
   p = SCIPprobdataGetNexvars(probdata);
   ndep = SCIPprobdataGetNdep(probdata);
   p3 = 3 * p;

   orig_xx = SCIPprobdataGetxx(probdata);
   orig_xy = SCIPprobdataGetxy(probdata);
   r = SCIPprobdataGetyy(probdata);
   var_z = SCIPprobdataGetVars_z(probdata);

   assert(n > 0);
   assert(p > 0);
   assert(ndep >= 0);
   assert(n > p);
   assert(orig_xx != NULL);
   assert(orig_xy != NULL);
   assert(var_z != NULL);

   if( p >= n  || SCIPprobdataGetSaveMemory(probdata) == TRUE )
   {
      SCIPerrorMessage("branching/frequent/priority should be set to 200000\n");
      SCIPerrorMessage("branching/fullstrong_maic/priority should be set to 100000\n");
      return SCIP_ERROR;
   }

   /* alloc */
   SCIP_CALL( SCIPallocBufferArray(scip, &list, p));
   SCIP_CALL( SCIPallocBufferArray(scip, &branchinfo, p3));
   SCIP_CALL( SCIPallocBufferArray(scip, &branchinfo_new, p3));

   SCIP_CALL( SCIPinitIntArrayZero(p, list));
   SCIP_CALL( SCIPinitIntArrayZero(p3, branchinfo));

   /* get pseudo candidates (non-fixed integer variables) */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, NULL, &ncands));

   for( i = 0; i < ncands; i++ )
   {
      for( j = 0; j < p; j++ )
      {
         if( cands[i] == var_z[j] )
         {
            list[j] = 1;
            break;
         }
      }
   }

#if MYPARA_LOG
   printf("list:");
   SCIPprintIntVec( p, list);
#endif

   /* get branching info */
   for( i = 0; i < p; ++i )
   {
      ublb = SCIPround(scip, SCIPcomputeVarUbLocal(scip, var_z[i])
            + SCIPcomputeVarLbLocal(scip, var_z[i]));
      *(branchinfo + (ublb * p) + i) = 1;
   }

#if MYPARA_LOG
   for( i = 0; i < 3; i++ )
   {
      for( j = 0; j < p; j++ )
      {
         printf("%d, ", *(branchinfo + (i * p) + j));
      }
      SCIPstartNewLine();
   }
#endif

   rss = -1.0;
   ind = -1;

   for( i = 0; i < p; i++ )
   {
      /*
       * solve
       *   sub_xx a = sub_xy
       */

      if( list[i] == 1 )
      {
         /* check whether z_i can be fixed to 1*/
         if( checkBranching(scip, p, ndep, probdata, branchinfo, i) == TRUE )
         {
            /* copy */
            for( j = 0; j < p3; j++ )
               branchinfo_new[j] = branchinfo[j];

            branchinfo_new[i] = 1;
            branchinfo_new[p+i] = 0;

            if( ndep )
            {
               /* fix variable by using linear dependence */
               SCIP_CALL( fixVariable(scip, p, ndep, probdata, branchinfo_new));
            }

            dim = p - SCIPcalcIntSum(&branchinfo_new[0], p);

            /* alloc */
            SCIP_CALL( SCIPallocBufferArray(scip, &a, dim));
            SCIP_CALL( SCIPallocBufferArray(scip, &sub_xx, dim*dim));
            SCIP_CALL( SCIPallocBufferArray(scip, &sub_xy, dim));

            /* generate Q and sub_xy */
            /* Q */
            ct = 0;
            for( j = 0; j < p; j++ )
            {
               if( branchinfo_new[j] == 0 )
               {
                  assert(i != j);

                  for( t = 0; t < p; t++ )
                  {
                     if( branchinfo_new[t] == 0 && t != i  )
                        sub_xx[ct++] = SCIPmatColMajor(orig_xx, p, j, t);
                  }
               }
            }

            assert( ct == dim * dim );

            /* sub_xy */
            ct = 0;
            for( j = 0; j < p; j++ )
            {
               if( branchinfo_new[j] ==0)
                  sub_xy[ct++] = orig_xy[j];
            }

            assert( ct == dim );

            dpv = SCIPclapackDposv(scip, sub_xx, sub_xy, dim, a);

            if( dpv == 0 )
            {
               /* test */
               rss_new = calcRSSvalue(dim, a, sub_xy, r);
               if( rss_new > rss )
               {
                  rss = rss_new;
                  ind = i;
               }
#if MYPARA_LOG
               printf("%d: rss = %f\n", i, rss_new);
#endif
            }

            /* free */
            SCIPfreeBufferArray(scip, &sub_xx);
            SCIPfreeBufferArray(scip, &sub_xy);
            SCIPfreeBufferArray(scip, &a);
         }
      }
   }

#if MYPARA_LOG
   printf("max->%dth var. \n", ind);
#endif

   if( ind == -1 )
   {
      /* free */
      SCIPfreeBufferArray(scip, &list);
      SCIPfreeBufferArray(scip, &branchinfo);
      SCIPfreeBufferArray(scip, &branchinfo_new);

      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPbranchVar( scip, var_z[ind], &childnode_0, NULL, &childnode_1));

   /* free */
   SCIPfreeBufferArray(scip, &list);
   SCIPfreeBufferArray(scip, &branchinfo);
   SCIPfreeBufferArray(scip, &branchinfo_new);

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}
#else
#define branchExecpsMyfullstrong NULL
#endif


/*
 * branching rule specific interface methods
 */

/** creates the myfullstrong branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleMyfullstrong(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create myfullstrong branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );

   /* TODO: (optional) create branching rule specific data here */

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsMyfullstrong) );
#if 0
   SCIP_CALL( SCIPsetBranchruleInitsol(scip, branchrule, branchInitsolMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextMyfullstrong) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpMyfullstrong) );
#endif

   return SCIP_OKAY;
}
