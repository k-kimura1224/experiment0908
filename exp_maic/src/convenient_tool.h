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

/**@file   convenient_tool.h
 * @brief  convenient tools for programming
 * @author Keiji Kimura
 *
 * This file implements some convenient functions.
 *
 * A list of the implemented functions is as follows:
 * - SCIPprintMat: print real matrix
 * - SCIPprintIntMat: print integer matrix
 * - SCIPstartNewLine: start a new line
 * - SCIPprintLongLine: print a long line
 * - SCIPprintVal: print a real value
 * - SCIPprintIntVal: print a integer value
 * - SCIPprintVec: print a real vector
 * - SCIPprintIntVec: print a integer vector
 * - SCIPcalcSum: calculate the sum of elements in real vector
 * - SCIPcalcIntSum: calculate the sum of elements in integer vector
 * - SCIPinitArrayZero: initialize all array elements to zero
 * - SCIPinitIntArrayZero: initialize all array elements to zero
 * - SCIPexit: exit a program
 * - SCIPsortBubble: execute bubble sort
 * - SCIPtransMat: transpose a matrix
 * - SCIPcalcTime: calculate time
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_CONVENIENT_TOOL_H__
#define __SCIP_CONVENIENT_TOOL_H__

#include <time.h>
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** print matrix */
extern
void SCIPprintMat(
   int                   n,                  /**< the number of rows */
   int                   m,                  /**< the number of columns */
   SCIP_Real*            Mat                 /**< matrix to print */
   );


/** print integer matrix */
extern
void SCIPprintIntMat(
   int                   n,                  /**< the number of rows */
   int                   m,                  /**< the number of columns */
   int*                  intmat              /**< integer matrix to print */
   );


/** start a new line */
extern
void SCIPstartNewLine(
   void
   );


/** print a long line */
extern
void SCIPprintLongLine(
   void
   );


/** print a real value */
extern
void SCIPprintVal(
   SCIP_Real             a
   );


/** print a integer value */
extern
void SCIPprintIntVal(
   int                   a
   );


/** print a real vector */
extern
void SCIPprintVec(
   int                   n,                  /**< dimension */
   SCIP_Real*            vec
   );


/** print a integer vector */
extern
void SCIPprintIntVec(
   int                   n,                  /**< dimension */
   int*                  vec
   );


/** calculate the sum of elements in real vector */
extern
SCIP_Real SCIPcalcSum(
   const SCIP_Real*      x,                  /**< real vector to calculate the sum of elements */
   int                   n                   /**< dimension */
   );


/** calculate the sum of elements in integer vector */
extern
int SCIPcalcIntSum(
   const int*            x,                  /**< integer vector to calculate the sum of elements */
   int                   n                   /**< dimension */
   );


/** initialize all array elements to zero */
extern
SCIP_RETCODE SCIPinitArrayZero(
   int                   n,                  /**< dimension */
   SCIP_Real*            x                   /**< array to initialize all elements to zero */
   );


/** initialize all array elements to zero */
extern
SCIP_RETCODE SCIPinitIntArrayZero(
   int                   n,                  /**< dimension */
   int*                  x                   /**< array to initialize all elements to zero */
   );


/** exit a program */
extern
void SCIPexit(
   void
   );


/** execute bubble sort */
extern
SCIP_RETCODE SCIPsortBubble(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   n,                  /**< dimension */
   SCIP_Real*            x,                  /**< array to sort elements */
   SCIP_Real*            y,                  /**< array to return sorted array */
   int*                  z                   /**< array to return indices of sorted array */
   );


/** transpose a matrix */
extern
SCIP_RETCODE SCIPtransMat(
   int                   n,                  /**< the number of rows */
   int                   m,                  /**< the number of columns */
   SCIP_Real*            inputmat,           /**< matrix to transpose elements */
   SCIP_Real*            outputmat           /**< matrix to return */
   );


/** calculate time */
extern
SCIP_Real SCIPcalcTime(
   clock_t               start,
   clock_t               end
   );


#ifdef __cplusplus
}
#endif

#endif
