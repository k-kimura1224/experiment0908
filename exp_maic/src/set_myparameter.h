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

/**@file   set_myparameter.h
 * @brief  setting parameters for this problem
 * @author Keiji Kimura
 *
 * This file implements a function to set parameters for this problem
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_SET_MYPARAMETER_H__
#define __SCIP_SET_MYPARAMETER_H__

#include "scip/scip.h"

#define MP_NUM_SOL      100      /* number of stored solution (df:10) */
#define MP_PRINTORIGP   0        /* if nonzero, print original problem.  */
#define MP_PRINTREFOP   0        /* if nonzero, print reformulated problem */

#ifdef __cplusplus
extern "C" {
#endif

/** set parameters */
extern
SCIP_RETCODE SCIPsetMyParameter(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
