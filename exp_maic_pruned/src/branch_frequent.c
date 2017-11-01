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

/**@file   branch_frequent.c
 * @brief  branching rule for this problem
 * @author Keiji Kimura
 *
 * This file implements the branching rule by using stored solutions
 *
 * The score of the binary variable z is calculated as follows:
 * score(z) = \sum_{solution \in pool} solution(z).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "branch_frequent.h"

#include "set_myparameter.h"
#include "convenient_tool.h"
#include "probdata_linereg.h"

#define BRANCHRULE_NAME            "frequent"
#define BRANCHRULE_DESC            "frequent variable on good models"
#if 1
#define BRANCHRULE_PRIORITY        200000
#else
#define BRANCHRULE_PRIORITY        100000
#endif
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0


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

/* put your local methods here, and declare them static */


/*
 * Callback methods of branching rule
 */

/* TODO: Implement all necessary branching rule methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyFrequent)
{  /*lint --e{715}*/

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleFrequent(scip) ) ;

   return SCIP_OKAY;
}


/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeFrequent)
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
SCIP_DECL_BRANCHINIT(branchInitFrequent)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
#if 1
static
SCIP_DECL_BRANCHEXIT(branchExitFrequent)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}
#else
#define branchExitFrequent NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BRANCHINITSOL(branchInitsolFrequent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myrule branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolFrequent NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolFrequent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myrule branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolFrequent NULL
#endif


/** branching execution method for fractional LP solutions */
#if 0
static
SCIP_DECL_BRANCHEXECLP(branchExeclpFrequent)
{  /*lint --e{715}*/

   printf("mybranching rule!\n");
   return SCIP_OKAY;
}
#else
#define branchExeclpFrequent NULL
#endif


/** branching execution method for external candidates */
#if 0
static
SCIP_DECL_BRANCHEXECEXT(branchExecextFrequent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of myrule branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextFrequent NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#if 1
static
SCIP_DECL_BRANCHEXECPS(branchExecpsFrequent)
{  /*lint --e{715}*/

   SCIP_NODE* childnode_0;
   SCIP_NODE* childnode_1;

   int nsols;
   SCIP_SOL** sols;
   int score;
   int maxscore;
   int indexbranching;

   int i;
   int j;

   /* get branching rule data */
   /*
   SCIP_BRANCHRULEDATA* branchruledata;
   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   */

#if 0
   SCIP_VAR**  cands;
   int ncands;

   /* get pseudo candidates (non-fixed integer variables) */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, NULL, &ncands) );
   assert(ncands > 0);

   nsols = SCIPgetNSols(scip);
   assert(nsols >= 0);

   if( nsols == 0 )
   {
      assert(cands[0] != NULL && SCIPvarGetType(cands[0]) == SCIP_VARTYPE_BINARY );
      SCIP_CALL( SCIPbranchVar(scip, cands[0], &childnode_0, NULL, &childnode_1));


      *result = SCIP_BRANCHED;
      return SCIP_OKAY;
   }

   sols = SCIPgetSols(scip);

   if( nsols > MP_NUM_SOL )
      nsols = MP_NUM_SOL;

   maxscore = -1;
   indexbranching = -1;
   for( i = 0; i < ncands; i++ )
   {
      score = 0;
      assert(SCIPvarGetType(cands[i]) == SCIP_VARTYPE_BINARY);

      for( j = 0; j < nsols; j++ )
         score += SCIPround(scip, SCIPgetSolVal(scip, sols[j], cands[i]));

      if( maxscore < score )
      {
         maxscore = score;
         indexbranching = i;

         if( maxscore == nsols )
            break;
      }
   }

   assert(indexbranching >= 0 && indexbranching < ncands);
   assert(SCIPvarGetType(cands[indexbranching]) == SCIP_VARTYPE_BINARY);

   SCIP_CALL( SCIPbranchVar(scip, cands[indexbranching], &childnode_0, NULL, &childnode_1));
   *result = SCIP_BRANCHED;
   return SCIP_OKAY;
#else
   SCIP_PROBDATA* probdata;
   SCIP_VAR** var_z;
   SCIP_VAR** var;
   int p;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   p = SCIPprobdataGetNexvars(probdata);
   var_z = SCIPprobdataGetVars_z(probdata);
   assert(p > 0);
   assert(var_z != NULL);

   nsols = SCIPgetNSols(scip);
   assert(nsols >= 0);

   if( nsols == 0 )
   {
      for( i = 0; i < p; i++ )
      {
         if( (int)SCIPround(scip, SCIPcomputeVarLbLocal(scip, var_z[i])
                  + SCIPcomputeVarUbLocal(scip, var_z[i])) == 1 )
            break;
      }
      SCIP_CALL( SCIPbranchVar(scip, var_z[i], &childnode_0, NULL, &childnode_1));

      *result = SCIP_BRANCHED;
      return SCIP_OKAY;
   }

   sols = SCIPgetSols(scip);

   if( nsols > MP_NUM_SOL )
      nsols = MP_NUM_SOL;

   maxscore = -1;
   indexbranching = -1;
   for( i = 0; i < p; i++ )
   {
      var = &var_z[i];
      if( (int)SCIPround(scip, SCIPcomputeVarLbLocal(scip, *var)
               + SCIPcomputeVarUbLocal(scip, *var)) == 1 )
      {
         score = 0;
         assert(SCIPvarGetType(*var) == SCIP_VARTYPE_BINARY);

         for( j = 0; j < nsols; j++ )
            score += SCIPround(scip, SCIPgetSolVal(scip, sols[j], *var));

         if( maxscore < score )
         {
            maxscore = score;
            indexbranching = i;

            if( maxscore == nsols )
               break;
         }
      }
   }

   assert(indexbranching >= 0 && indexbranching < p);
   assert(SCIPvarGetType(var_z[indexbranching]) == SCIP_VARTYPE_BINARY);

   SCIP_CALL( SCIPbranchVar(scip, var_z[indexbranching], &childnode_0, NULL, &childnode_1));
   *result = SCIP_BRANCHED;
   return SCIP_OKAY;
#endif
}
#else
#define branchExecpsFrequent NULL
#endif

/*
 * branching rule specific interface methods
 */

/** creates the myrule branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleFrequent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create myrule branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );

   /* TODO: (optional) create branching rule specific data here */


   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyFrequent) );
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeFrequent) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitFrequent) );
   SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitFrequent) );
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsFrequent) );
#if 0
   SCIP_CALL( SCIPsetBranchruleInitsol(scip, branchrule, branchInitsolFrequent) );
   SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolFrequent) );
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextFrequent) );
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpFrequent) );
#endif

   return SCIP_OKAY;
}
