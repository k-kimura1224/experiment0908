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

/**@file   call_cblas.c
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


#if defined(__APPLE__)
#include "cblas.h"
#include "clapack.h"
#else
#include <cblas.h>
#ifdef __cpluscplus
extern "C" {
#endif
int dposv_( char* uplo, int* n, int* nrhs, double* A,
      int* lda, double* x, int* ldb, int* info);
#ifdef __cpluscplus
}
#endif
#endif

#include "call_cblas.h"
#include "convenient_tool.h"


/** return the element (i,j) of matrix \in R^{n times m} with ColMajor */
SCIP_Real SCIPmatColMajor(
   const SCIP_Real*      matrix,             /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             i,                  /**< index of row */
   const int             j                   /**< index of column */
   )
{
   assert(matrix != NULL);
   assert(n > 0);
   assert(i >= 0);
   assert(j >= 0);

   return *(matrix + i + (j * n));
}


/** print matrix with ColMajor */
void  SCIPprintMatColMajor(
   const SCIP_Real*      matrix,             /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             m                   /**< the number of columns */
   )
{
   int i;
   int j;

   assert(matrix != NULL);
   assert(n > 0);
   assert(m > 0);

   printf("\n");

   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < m; j++ )
      {
         printf(" %f,", SCIPmatColMajor(matrix, n, i, j));
      }
      printf("\n");
   }
}


/** copy array */
SCIP_RETCODE SCIPcblasCopy(
   const SCIP_Real*      x,                  /**< original array */
   SCIP_Real*            y,                  /**< array for copy */
   const int             n                   /**< size of array */
   )
{
   int inc = 1;

   assert(x != NULL);
   assert(y != NULL);
   assert(n > 0);

   cblas_dcopy(n, x, inc, y, inc);

   return SCIP_OKAY;
}


/** create submatrix */
SCIP_RETCODE SCIPcreateSubmatColMajor(
   const SCIP_Real*      origmat,            /**< original matrix with ColMajor*/
   const int             n,                  /**< the number of rows original matrix */
   const int             m,                  /**< the number of columns of original matrix */
   const int*            x,                  /**< array to pick up column vectos */
   SCIP_Real*            submat              /**< array to store submatrix */
   )
{
   int i;
   int ct = 0;

   assert(origmat != NULL);
   assert(n > 0);
   assert(m > 0);
   assert(x != NULL);
   assert(SCIPcalcIntSum(x, m) > 0);
   assert(submat != NULL);

   for( i = 0; i < m; i++ )
   {
      assert(x[i] == 1 || x[i] == 0);
      if( *(x + i) == 1 )
      {
         SCIP_CALL( SCIPcblasCopy(origmat + (i * n), submat + (ct * n), n));
         ct++;
      }
   }

   return SCIP_OKAY;
}


/** compute a linear system with Cholesky decomposition;
 *
 *  This function assumes that the matrix is positive definite and symmetric.
 *
 *  @return Return 0 if solving process is successful
 */
int SCIPclapackDposv(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real*      matrix,             /**< matrix with ColMajor */
   const SCIP_Real*      vector,             /**< vector */
   const int             dim,                /**< dimension */
   SCIP_Real*            solution            /**< array to store solution */
   )
{
   char uplo[2] = "L";
   int one = 1;
   int info;
   int n = dim;
   SCIP_Real* copymatrix;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(vector != NULL);
   assert(dim > 0);
   assert(solution != NULL);

   /* allocate memory for copymatrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &copymatrix, n * n) );

   SCIP_CALL( SCIPcblasCopy(matrix, copymatrix, n*n));
   SCIP_CALL( SCIPcblasCopy(vector, solution, n));

   dposv_(uplo, &n, &one, copymatrix, &n, solution, &n, &info);

   /* free */
   SCIPfreeBufferArray(scip, &copymatrix);

   return info;
}


/** call dgemv function of cblas
 *
 *  This function computes
 *  z[n] = (alpha) A[n*m] x[m] + (beta) y[n]
 */
SCIP_RETCODE SCIPcblasDgemv1(
   const SCIP_Real*      matrix_A,            /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             m,                  /**< the number of columns */
   const SCIP_Real*      vector_x,
   const SCIP_Real*      vector_y,
   const SCIP_Real       alpha,
   const SCIP_Real       beta,
   SCIP_Real*            vector_z            /**< array to store vector z */
   )
{
   int one = 1;

   assert(matrix_A != NULL);
   assert(n > 0);
   assert(m > 0);
   assert(vector_x != NULL);
   assert(vector_y != NULL);
   assert(vector_x != vector_z);
   assert(alpha > 0 || beta > 0);
   assert(vector_z != NULL);

   /* copy */
   SCIP_CALL( SCIPcblasCopy(vector_y, vector_z, n));

   cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, alpha, matrix_A, n, vector_x, one, beta, vector_z, one);

   return SCIP_OKAY;
}


/** call dgemv function of cblas
 *
 *  This function computes
 *  z[n] = A[n*m] x[m]
 */
SCIP_RETCODE SCIPcblasDgemv2(
   const SCIP_Real*      matrix_A,           /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             m,                  /**< the number of columns */
   const SCIP_Real*      vector_x,
   SCIP_Real*            vector_z            /**< array to store vector z*/
   )
{
   int one = 1;

   assert(matrix_A != NULL);
   assert(n > 0);
   assert(m > 0);
   assert(vector_x != NULL);
   assert(vector_z != NULL);
   assert(vector_x != vector_z);

   /* initialize */
   SCIP_CALL( SCIPinitArrayZero(n, vector_z));

   cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, 1.0, matrix_A, n, vector_x, one, 0.0, vector_z, one);

   return SCIP_OKAY;
}


/** call dgemv function of cblas
 *
 *  This function computes
 *  z[m] = (A[n*m])^T x[n]
 */
SCIP_RETCODE SCIPcblasDgemv3(
   const SCIP_Real*      matrix_A,           /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             m,                  /**< the number of columns */
   const SCIP_Real*      vector_x,
   SCIP_Real*            vector_z            /**< array to store vector z */
   )
{
   int one = 1;

   assert(matrix_A != NULL);
   assert(n > 0);
   assert(m > 0);
   assert(vector_x != NULL);
   assert(vector_z != NULL);
   assert(vector_x != vector_z);

   /* initialize */
   SCIP_CALL( SCIPinitArrayZero(m, vector_z));

   cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, matrix_A, n, vector_x, one, 0.0, vector_z, one);

   return SCIP_OKAY;
}


/** call dgemm function of cblas
 *
 *  This function computes
 *  B[m*m] = (A[n*m])^T A[n*m]
 */
SCIP_RETCODE SCIPcblasDgemm1(
   const SCIP_Real*      matrix_A,           /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             m,                  /**< the number of columns */
   SCIP_Real*            matrix_B            /**< array to stor matrix B */
   )
{
   SCIP_Real   zero=0.0;
   SCIP_Real   one=1.0;

   assert(matrix_A != NULL);
   assert(matrix_B != NULL);
   assert(matrix_A != matrix_B);
   assert(n > 0);
   assert(m > 0);

   /* initialize */
   SCIP_CALL( SCIPinitArrayZero(m * m, matrix_B));

   cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, m, n, one, matrix_A, n, matrix_A, n, zero, matrix_B, m);

   return SCIP_OKAY;
}


/** call ddot function of cblas
 *
 *  This function computes the dot product
 */
SCIP_Real SCIPcblasDdot(
   const SCIP_Real*      x,                  /**< array to compute vector */
   const SCIP_Real*      y,                  /**< array to compute vector */
   const int             n                   /**< dimention */
   )
{
   int one = 1;
   SCIP_Real xy;

   assert(x != NULL);
   assert(y != NULL);
   assert(n > 0);

   xy = cblas_ddot(n, x, one, y, one);

   return xy;
}


/** call dscal function of cblas
 *
 *  This function scales a vector by a constant
 */
SCIP_RETCODE SCIPcblasDscal(
   const SCIP_Real*      x,                  /**< array */
   const int             n,                  /**< size */
   const SCIP_Real       alpha,              /**< constant */
   SCIP_Real*            y                   /**< array to store scaled vector */
   )
{
   int one = 1;

   assert(x != NULL);
   assert(n > 0);
   assert(y != NULL);

   if( x != y )
      SCIP_CALL( SCIPcblasCopy(x, y, n));

   cblas_dscal(n, alpha, y, one);

   return SCIP_OKAY;
}


/** call daxpy function of cblas
 *
 *  This function computes
 *  z[n] = (alpha) x[n] + (beta) y[n]
 */
SCIP_RETCODE SCIPcblasDaxpy(
   const SCIP_Real*      x,                  /**< array for vector x */
   const SCIP_Real*      y,                  /**< array for vector y */
   const int             n,                  /**< size of array */
   const SCIP_Real       alpha,              /**< scalar for vector x */
   const SCIP_Real       beta,               /**< scalar for vector y */
   SCIP_Real*            z                   /**< array to store vector z */
   )
{
   int one = 1;

   assert(x != NULL);
   assert(y != NULL);
   assert(z != NULL);
   assert(n > 0);

   if( beta != 1.0 )
      SCIP_CALL( SCIPcblasDscal(y, n, beta, z));

   cblas_daxpy(n, alpha, x, one, z, one);

   return SCIP_OKAY;
}


/** call dnrm2 of cblas
 *
 * This function returns the euclidean norm of a vector
 */
SCIP_Real SCIPcblasDnrm(
   const SCIP_Real*      x,                  /**< array for a vector */
   const int             n                   /**< size of a vector */
   )
{
   int one = 1;

   assert(x != NULL);
   assert(n > 0);

   return cblas_dnrm2(n, x, one);
}


/** call dger function of cblas
 *
 *  This function computes
 *  B[n*m] = (alpha) x[n] (y[m])^t + A[n*m]
 */
SCIP_RETCODE SCIPcblasDger(
   const SCIP_Real*      matrix_A,           /**< array for matrix A */
   const SCIP_Real*      vector_x,           /**< array for vector x */
   const SCIP_Real*      vector_y,           /**< array for vector y */
   const int             n,                  /**< size of vextor x */
   const int             m,                  /**< size of vector y */
   const SCIP_Real       alpha,              /**< scalar */
   SCIP_Real*            matrix_B            /**< array to store matrix B */
   )
{
   int one = 1;

   assert(matrix_A != NULL);
   assert(vector_x != NULL);
   assert(vector_y != NULL);
   assert(n > 0);
   assert(m > 0);
   assert(matrix_B != NULL);

   /* copy */
   if( matrix_A != matrix_B )
      SCIP_CALL( SCIPcblasCopy(matrix_A, matrix_B, n * m));

   cblas_dger(CblasColMajor, n, m, alpha, vector_x, one, vector_y, one, matrix_B, n);

   return SCIP_OKAY;
}
