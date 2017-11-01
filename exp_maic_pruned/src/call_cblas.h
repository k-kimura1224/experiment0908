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

/**@file   call_cblas.h
 * @brief  functions for calling cblas
 * @author Keiji Kimura
 *
 * This file implements some functions for calling cblas.
 * All matrices in this file are ColMajor.
 *
 * A list of the implemented functions is as follows:
 * - SCIPmatColMajor: return the element of matrix
 * - SCIPprintMatColMajor: print matrix with ColMajor
 * - SCIPcblasCopy: copy array
 * - SCIPcreateSubmatColMajor: create submatrix
 * - SCIPclapackDposv: compute a linear system with Cholesky decomposition
 * - SCIPcblasDgemv1: call dgemv function
 * - SCIPcblasDgemv2: call dgemv function
 * - SCIPcblasDgemv3: call dgemv function
 * - SCIPcblasDgemm1: call dgemm function
 * - SCIPcblasDdot: call ddot function
 * - SCIPcblasDscal: call dscal function
 * - SCIPcblasDaxpy: call daxpy function
 * - SCIPcblasDnrm: call dnrm2 function
 * - SCIPcblasDger: call dger function
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_CALL_CBLAS_H__
#define __SCIP_CALL_CBLAS_H__

#include "scip/scip.h"


/** return the element (i,j) of matrix \in R^{n times m} with ColMajor */
extern
SCIP_Real SCIPmatColMajor(
   const SCIP_Real*      matrix,             /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             i,                  /**< index of row */
   const int             j                   /**< index of column */
   );


/** print matrix with ColMajor */
extern
void  SCIPprintMatColMajor(
   const SCIP_Real*      matrix,             /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             m                   /**< the number of columns */
   );


/** copy array */
extern
SCIP_RETCODE SCIPcblasCopy(
   const SCIP_Real*      x,                  /**< original array */
   SCIP_Real*            y,                  /**< array for copy */
   const int             n                   /**< size of array */
   );


/** create submatrix */
extern
SCIP_RETCODE SCIPcreateSubmatColMajor(
   const SCIP_Real*      origmat,            /**< original matrix with ColMajor*/
   const int             n,                  /**< the number of rows original matrix */
   const int             m,                  /**< the number of columns of original matrix */
   const int*            x,                  /**< array to pick up column vectos */
   SCIP_Real*            submat              /**< array to store submatrix */
   );


/** compute a linear system with Cholesky decomposition;
 *
 *  This function assumes that the matrix is positive definite and symmetric.
 *
 *  @return Return 0 if solving process is successful
 */
extern
int SCIPclapackDposv(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real*      matrix,             /**< matrix with ColMajor */
   const SCIP_Real*      vector,             /**< vector */
   const int             dim,                /**< dimension */
   SCIP_Real*            solution            /**< array to store solution */
   );


/** call dgemv function of cblas
 *
 *  This function computes
 *  z[n] = (alpha) A[n*m] x[m] + (beta) y[n]
 */
extern
SCIP_RETCODE SCIPcblasDgemv1(
   const SCIP_Real*      mtrix_A,            /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             m,                  /**< the number of columns */
   const SCIP_Real*      vector_x,
   const SCIP_Real*      vector_y,
   const SCIP_Real       alpha,
   const SCIP_Real       beta,
   SCIP_Real*            vector_z            /**< array to store vector z */
   );


/** call dgemv function of cblas
 *
 *  This function computes
 *  z[n] = A[n*m] x[m]
 */
extern
SCIP_RETCODE SCIPcblasDgemv2(
   const SCIP_Real*      matrix_A,           /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             m,                  /**< the number of columns */
   const SCIP_Real*      vector_x,
   SCIP_Real*            vector_z            /**< array to store vector z*/
   );


/** call dgemv function of cblas
 *
 *  This function computes
 *  z[m] = (A[n*m])^T x[n]
 */
extern
SCIP_RETCODE SCIPcblasDgemv3(
   const SCIP_Real*      matrix_A,           /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             m,                  /**< the number of columns */
   const SCIP_Real*      vector_x,
   SCIP_Real*            vector_z            /**< array to store vector z */
   );


/** call dgemm function of cblas
 *
 *  This function computes
 *  B[m*m] = (A[n*m])^T A[n*m]
 */
extern
SCIP_RETCODE SCIPcblasDgemm1(
   const SCIP_Real*      matrix_A,           /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             m,                  /**< the number of columns */
   SCIP_Real*            matrix_B            /**< array to stor matrix B */
   );


/** call ddot function of cblas
 *
 *  This function computes the dot product
 */
extern
SCIP_Real SCIPcblasDdot(
   const SCIP_Real*      x,                  /**< array to compute vector */
   const SCIP_Real*      y,                  /**< array to compute vector */
   const int             n                   /**< dimention */
   );


/** call dscal function of cblas
 *
 *  This function scales a vector by a constant
 */
extern
SCIP_RETCODE SCIPcblasDscal(
   const SCIP_Real*      x,                  /**< array */
   const int             n,                  /**< size */
   const SCIP_Real       alpha,              /**< constant */
   SCIP_Real*            y                   /**< array to store scaled vector */
   );


/** call daxpy function of cblas
 *
 *  This function computes
 *  z[n] = (alpha) x[n] + (beta) y[n]
 */
extern
SCIP_RETCODE SCIPcblasDaxpy(
   const SCIP_Real*      x,                  /**< array for vector x */
   const SCIP_Real*      y,                  /**< array for vector y */
   const int             n,                  /**< size of array */
   const SCIP_Real       alpha,              /**< scalar for vector x */
   const SCIP_Real       beta,               /**< scalar for vector y */
   SCIP_Real*            z                   /**< array to store vector z */
   );


/** call dnrm2 of cblas
 *
 * This function returns the euclidean norm of a vector
 */
extern
SCIP_Real SCIPcblasDnrm(
   const SCIP_Real*      x,                  /**< array for a vector */
   const int             n                   /**< size of a vector */
   );


/** call dger function of cblas
 *
 *  This function computes
 *  B[n*m] = (alpha) x[n] (y[m])^t + A[n*m]
 */
extern
SCIP_RETCODE SCIPcblasDger(
   const SCIP_Real*      matrix_A,           /**< array for matrix A */
   const SCIP_Real*      vector_x,           /**< array for vector x */
   const SCIP_Real*      vector_y,           /**< array for vector y */
   const int             n,                  /**< size of vextor x */
   const int             m,                  /**< size of vector y */
   const SCIP_Real       alpha,              /**< scalar */
   SCIP_Real*            matrix_B            /**< array to store matrix B */
   );


#endif
