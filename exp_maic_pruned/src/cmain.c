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

/**@file   cmain.c
 * @brief  Main file for C compilation
 * @author Keiji Kimura
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <stdio.h>
#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"
#include "convenient_tool.h"
#include "reader_linereg.h"
#include "probdata_linereg.h"
#include "relax_dposv.h"
#include "set_myparameter.h"
#include "heur_forward.h"
#include "heur_backward.h"
#include "branch_frequent.h"
#include "branch_myfullstrong.h"

/** print some information about the best solution
 *  a: continuous variables
 *  z: binary varibales
 *  rss: residual sum of squares
 *  AIC: AIC value
 *  k: the number of selected explanatory variables, i.e., \sum z
 */
static
SCIP_RETCODE outputSolution(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_SOL* sol;
   SCIP_VAR* var;
   SCIP_Real solval;
   SCIP_Real abic;
   SCIP_Real para_reg;
   char varname[SCIP_MAXSTRLEN];
   int n;
   int p;
   int solvalint;
   int ct;
   int i;

   assert(scip != NULL);

   printf("\nSolution:\n");
   if( SCIPgetNSols(scip) == 0 )
   {
      printf("Not found\n");
      return SCIP_OKAY;
   }

   /* get the best solution */
   sol = SCIPgetBestSol(scip);

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   n = SCIPprobdataGetNdatas(probdata);
   p = SCIPprobdataGetNexvars(probdata);
   para_reg = SCIPprobdataGetPara_RegTerm(probdata);
   assert(n > 0);
   assert(p > 0);
   assert((SCIPprobdataGetMode(probdata) == 'a' && EPSEQ(para_reg, 2.0, 1e-08))
         || (SCIPprobdataGetMode(probdata) == 'b' && EPSEQ(para_reg, log((SCIP_Real) n), 1e-08)) );

   abic = (double) n * log( 2.0 * M_PI ) + (double) n + para_reg;

   /* print a part of the best solution */
   ct = 0;
   printf("\tz\t\t\ta\n\n");
   for( i = 0; i < p; i++ )
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "a%d", i + 1);
      var = SCIPfindVar(scip, varname);
      solval = SCIPgetSolVal(scip, sol, var);
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z%d", i + 1);
      var = SCIPfindVar(scip, varname);
      solvalint = SCIPgetSolVal(scip, sol, var);
      printf("%2d:\t%d\t%20.15g\n", i + 1, solvalint, solval);
      if( solvalint == 1 )
         ct++;

      abic += para_reg * (double) solvalint;
   }

   (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "rss", NULL);
   var = SCIPfindVar(scip, varname);
   solval = SCIPgetSolVal(scip, sol, var);
   abic += (double) n * log( solval / (double) n );

   /* print rss, AIC/BIC and k */
   printf("\nrss : %20.15g\n", solval);
   if( SCIPprobdataGetMode(probdata) == 'a' )
      printf("AIC : %20.15g\n", abic);
   else
      printf("BIC : %20.15g\n", abic);
   printf("k   :\t  %d\n\n", ct);

   return SCIP_OKAY;
}


/** print the frequency of selection in the top k solutions */
static
SCIP_RETCODE outputFrequency(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_SOL** sols;
   SCIP_VAR** var_z;
   int p;
   int z_val;
   int* n_z; /* frequency */
   int nsols;
   int k;
   int i;
   int j;

   assert(scip != NULL);

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   p = SCIPprobdataGetNexvars(probdata);
   var_z = SCIPprobdataGetVars_z(probdata);
   nsols = SCIPgetNSols(scip);
   assert(p > 0);

   printf("\nSolution:\n");

   k = 100;
   if( k > nsols )
      k = nsols;

   /* if there are no feasible solutions, return */
   if( k == 0 )
      return SCIP_OKAY;

   sols = SCIPgetSols(scip);

   /* alloc */
   SCIP_CALL( SCIPallocBufferArray(scip, &n_z, p));
   SCIP_CALL( SCIPinitIntArrayZero( p, n_z));

   /* count */
   for( i = 0; i < k; i++ )
   {
      for( j = 0; j < p; j++ )
      {
         z_val = SCIPgetSolVal( scip, sols[i], var_z[j]);
         if( z_val == 1 )
            n_z[j]++;
      }
   }

   /* print */
   for( i = 0; i < p; i++ )
   {
      printf("%d -- %d \n", i + 1, n_z[i]);
   }

   /* free */
   SCIPfreeBufferArray(scip, &n_z);

   return SCIP_OKAY;
}


static
SCIP_RETCODE runShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array containing shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include reader */
   SCIP_CALL( SCIPincludeReaderLinereg(scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxDposv(scip) );

   /* include primal heuristics based on forward selection */
   SCIP_CALL( SCIPincludeHeurForward(scip) );

   /* include primal heuristics based on backward elimination */
   SCIP_CALL( SCIPincludeHeurBackward(scip) );

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleFrequent(scip) );

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleMyfullstrong(scip) );

   /* set LINEREG-specific default parameters */
   SCIP_CALL( SCIPsetMyParameter(scip) );

   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   /* output the best solution */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_CALL( outputSolution(scip) );
      if( 0 )
         SCIP_CALL( outputFrequency(scip) );
   }

   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                   argc,               /**< number of shell parameters */
   char**                argv                /**< array containing shell parameters */
   )
{
   SCIP_RETCODE retcode;

   retcode = runShell(argc, argv, "scip.set");
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
