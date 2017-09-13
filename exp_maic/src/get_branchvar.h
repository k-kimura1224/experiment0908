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

/**@file   get_branchvar.h
 * @brief  function to return the last branching variable
 * @author Keiji Kimura
 *
 * This file implements the function to return the last branching variable
 *
 * If the last branching (binary) variable z_j = 1, the relaxator can compute the lower bound easily.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_BRANCHVAR_H__
#define __SCIP_BRANCHVAR_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** retrun the last branching variable */
extern
SCIP_RETCODE SCIPvarGetLastBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< the current node */
   SCIP_VAR**            result              /**< pointer to store the last branching variable */
   );

#ifdef __cplusplus
}
#endif

#endif
