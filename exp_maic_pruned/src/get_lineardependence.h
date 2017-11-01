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

/**@file   get_lineardependence.h
 * @brief  handler linear dependence
 * @author Keiji Kimura
 *
 * This file implements some functions for linear dependence
 *
 * If the column vectors of data are linearly dependent, we say that the data has linearly dependent
 * variables. we can efficiently compute the optimal value of this problem by using linear
 * dependence. SCIPgetNLineDependSet() and SCIPgetLineDependSet() find sets of linearly dependent
 * variables from the data.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GET_LINEARDEPENDENCE_H__
#define __GET_LINEARDEPENDENCE_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** retrun the number of linearly dependent sets from among colmun vectors */
extern
SCIP_RETCODE SCIPgetNLineDependSet(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real*      origmatrix,         /**< array for matrix */
   const int             n,                  /**< the number of rows */
   const int             m,                  /**< the number of columns */
   int*                  result              /**< array to return reslut */
   );


/** find linearly dependent sets from among colmun vectors */
extern
SCIP_RETCODE SCIPgetLineDependSet(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real*      matrix,             /**< symmetric matrix */
   const int             m,                  /**< dimension of matrix */
   const int             ndep,               /**< the number of linearly dependent sets */
   const int*            maxdep,             /**< maximal index at each set */
   const int*            ldindex,            /**< index of linearly dependent sets */
   int*                  depsets             /**< linearly dependent sets */
   );


/** print linearly dependent sets */
extern
SCIP_RETCODE SCIPprintLineDependSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   ndep,               /**< the number of linearly dependent set */
   int                   p,                  /**< the number of explanatory variables */
   int*                  depsets             /**< linearly dependent sets */
   );


#ifdef __cplusplus
}
#endif

#endif
