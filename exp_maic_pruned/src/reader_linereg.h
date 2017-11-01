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

/**@file   reader_linereg.h
 * @brief  file reader for linereg regression data
 * @author Keiji Kimura
 *
 * This file implements the reader for linear regression data.
 *
 * This reader reads from the second line of a file.
 * The format of the file is as follows:
 * - [2nd line] the number of data
 * - [3rd line] the number of the explanatory variables
 * - [4th line] the index of the response variable(Begining is 1)
 * - [from 5th line] data points
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_READER_LINEREG_H__
#define __SCIP_READER_LINEREG_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the stp file reader in SCIP */
extern
SCIP_RETCODE SCIPincludeReaderLinereg(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
