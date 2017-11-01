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

/**@file   set_myparameter.c
 * @brief  setting parameters for this problem
 * @author Keiji Kimura
 *
 * This file implements a function to set parameters for this problem
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "set_myparameter.h"

/** set parameters */
SCIP_RETCODE SCIPsetMyParameter(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* propagating */
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/pseudoobj/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/vbounds/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/dualfix/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/linear/propfreq", -1));

   /* branching rule */
   /*
   SCIP_CALL( SCIPsetBoolParam(scip, "branching/preferbinary", TRUE));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/inference/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/mostinf/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/allfullstrong/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/cloud/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/distribution/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/fullstrong/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/leastinf/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/multaggr/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/nodereopt/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/random/maxdepth", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/relpscost/maxdepth", 0));
   */

   /* no presolving */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0));

   /* if -1, disable the LP relaxation and only use my custom relaxation */
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1));

   /**
    * maximal number of solutions candidates to store
    * in the solution storage of the original problem [default:10]
   **/
   SCIP_CALL( SCIPsetIntParam(scip, "limits/maxorigsol", MP_NUM_SOL));
   SCIP_CALL( SCIPsetIntParam(scip, "limits/maxsol", MP_NUM_SOL));

   /* if -1, disable the all primal heuristic */
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/actconsdiving/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/clique/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/coefdiving/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/crossover/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/dins/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/dualval/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/feaspump/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/fixandinfer/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/fracdiving/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/guideddiving/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/intdiving/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/intshifting/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/linesearchdiving/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/localbranching/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/mutation/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/nlpdiving/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/objpscostdiving/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/octane/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/oneopt/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/proximity/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/pscostdiving/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/randrounding/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rens/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rins/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rootsoldiving/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rounding/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/shiftandpropagate/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/shifting/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/simplerounding/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/subnlp/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/trivial/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/twoopt/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/undercover/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/vbounds/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/veclendiving/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/zeroobj/freq", -1));

   return SCIP_OKAY;
}
