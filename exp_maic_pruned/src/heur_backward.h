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

/**@file   heur_backward.h
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

#ifndef __SCIP_HEUR_BACKWARD_H__
#define __SCIP_HEUR_BACKWARD_H__

#include "scip/scip.h"
#ifdef __cplusplus
extern "C" {
#endif

/** creates the local primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurBackward(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
