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

/**@file   relax_dposv.h
 * @brief  computing a lower bound
 * @author Keiji Kimura
 *
 * This file implements a method to compute a lower bound
 *
 * Although the relaxation problem of this problem is nonconvex, we can compute the lower bound
 * by solving a linear system. This file calls DPOSV implemented in LAPACK.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_RELAX_DPOSV_H__
#define __SCIP_RELAX_DPOSV_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the myrelaxator relaxator and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeRelaxDposv(
	SCIP*						scip
   );

#ifdef __cplusplus
}
#endif

#endif
