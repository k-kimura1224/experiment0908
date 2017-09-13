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

/**@file   get_lineardependence.c
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


#include "get_lineardependence.h"
#include "convenient_tool.h"
#include "call_cblas.h"


/** calclulate the number of linearly dependent sets from among colmun vectors */
SCIP_RETCODE SCIPgetNLineDependSet(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real*      origmatrix,         /**< array with ColMajor */
   const int             n,                  /**< the number of rows */
   const int             m,                  /**< the number of columns */
   int*                  result              /**< array to return reslut */
   )
{
   int i;
   int j;
   int memo = 0;
   int add  = 0;
   int rank = 1;
   SCIP_Real norm;
   SCIP_Real ep = 10e-06;
   SCIP_Real norm_x0;

   SCIP_Real *submatrix = NULL;
   SCIP_Real *invmatrix = NULL;
   SCIP_Real *vector_a = NULL;
   SCIP_Real *vector_b = NULL;
   SCIP_Real *vector_c = NULL;

   SCIP_Real *newsubmatrix = NULL;
   SCIP_Real *newinvmatrix = NULL;
   SCIP_Real *vector_newa = NULL;
   SCIP_Real *vector_newb = NULL;

   SCIP_Real *matrix_P = NULL;
   SCIP_Real *vector_q = NULL;
   SCIP_Real scalar_r;

   assert(scip != NULL);
   assert(origmatrix != NULL);
   assert(n > 0);
   assert(m > 1);
   assert(result != NULL);

   norm_x0 = SCIPcblasDdot(origmatrix, origmatrix, n);
   assert(norm_x0 > ep);

   *result = 0;

   /* allocate memory */
   assert(rank == 1);
   SCIP_CALL( SCIPallocMemoryArray(scip, &vector_a, rank));
   SCIP_CALL( SCIPallocMemoryArray(scip, &vector_b, n));
   SCIP_CALL( SCIPallocMemoryArray(scip, &submatrix, n * rank));
   SCIP_CALL( SCIPallocMemoryArray(scip, &invmatrix, rank * rank));

   /* initialize submatrix, inverse matrix, vector a and b
    * submatrix := the first column vector (\neq 0 ) of origmatrix
    * invmatrix := (submatrix'submatrix)^-1
    * a := invmatrix submatrix' x_k,
    * b := x_k - submatrix a
    * where x_k is the k-th column vector of origmatrix and
    * k = min{ i : 1<=i<=(n-1), x_0 and x_i are linearly independent }
    */
   for( i = 1; i < m; i++ )
   {
      memo = i;
      /* calclulate vector a */
      *vector_a = SCIPcblasDdot(origmatrix, origmatrix + (i * n), n) / norm_x0;
      /* calclulate vector b */
      SCIP_CALL( SCIPcblasDaxpy(origmatrix + (i * n), origmatrix, n, 1.0, - (*vector_a), vector_b));

      norm = SCIPcblasDnrm(vector_b, n);

      if( norm < ep )
      {
         *(result + i) = 1;
      }
      else
      {
         *(result + i) = 0;
         break;
      }
   }

   /* set submatrix */
   SCIP_CALL( SCIPcblasCopy(origmatrix, submatrix, n));
   /* set invmatrix */
   *invmatrix = 1.0 / norm_x0;

   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &vector_newb, n));


   add = 1;

   /* find linearly dependent setes */
   for( i = memo + 1; i < m; i++ )
   {
      if( add )
      {
         /* allocate memory */
         SCIP_CALL( SCIPallocMemoryArray(scip, &matrix_P, rank * rank));
         SCIP_CALL( SCIPallocMemoryArray(scip, &vector_q, rank));
         rank++;
         SCIP_CALL( SCIPallocMemoryArray(scip, &newsubmatrix, n * rank));
         SCIP_CALL( SCIPallocMemoryArray(scip, &newinvmatrix, rank * rank));
         SCIP_CALL( SCIPallocMemoryArray(scip, &vector_newa, rank));
         SCIP_CALL( SCIPallocMemoryArray(scip, &vector_c, rank));

         /* set newsubmatrix */
         SCIP_CALL( SCIPcblasCopy(submatrix, newsubmatrix, n * (rank - 1)));
         SCIP_CALL( SCIPcblasCopy(origmatrix + (memo * n), newsubmatrix + (n * (rank - 1)), n));
         add = 0;

         /* compute r, q and P
          * r = 1 / x_{memo}'
          * q = - r a
          * P = invmatrix - a q'
          */
         /* compute r */
         scalar_r = 1.0 / SCIPcblasDdot(origmatrix + (memo * n), vector_b, n);
         /* compute q */
         SCIP_CALL( SCIPcblasDscal(vector_a, rank - 1, - scalar_r, vector_q));
         /* compute P */
         SCIP_CALL( SCIPcblasDger(invmatrix, vector_a, vector_q, rank - 1, rank - 1, - 1.0, matrix_P));

         /* set newinvmatrix from r, q and P
          * newinvmatrix = ( P  q )
          *                ( q' r )
          */
         for( j = 0; j < rank - 1; j++ )
         {
            SCIP_CALL( SCIPcblasCopy(matrix_P + (j * (rank - 1)), newinvmatrix + (j * rank), rank - 1));
            newinvmatrix[rank - 1 + (j * rank)] = vector_q[j];
         }
         SCIP_CALL( SCIPcblasCopy(vector_q, newinvmatrix + (rank * (rank - 1)), rank - 1));
         newinvmatrix[(rank * rank) - 1] = scalar_r;
      }

      /* compute c, new a and new b
       * c = newsubmatrix' x_i
       * newa = newinvmatrix c
       * newb = x_i - newsubmatrix newa
       */
      /* compute c */
      SCIP_CALL( SCIPcblasDgemv3(newsubmatrix, n, rank, origmatrix + (i * n), vector_c));
      /* compute new a */
      SCIP_CALL( SCIPcblasDgemv1(newinvmatrix, rank, rank, vector_c, vector_c, 1.0, 0.0, vector_newa));
      /* compute new b */
      SCIP_CALL( SCIPcblasDgemv1(newsubmatrix, n, rank, vector_newa, origmatrix + (i * n), - 1.0, 1.0, vector_newb));

      /* compute norm of vector_newb */
      norm = SCIPcblasDnrm(vector_newb, n);

      /* if the norm of new b is zero, linear dependence is found */
      if( norm < ep )
      {
         *(result + i) = 1;
      }
      else
      {
         *(result + i) = 0;
         memo = i;
         add = 1;

         if( i < m - 1 )
         {
            /* free */
            SCIPfreeMemoryArrayNull(scip, &invmatrix);
            SCIPfreeMemoryArrayNull(scip, &vector_a);
            SCIPfreeMemoryArrayNull(scip, &submatrix);

            /* allocate */
            SCIP_CALL( SCIPallocMemoryArray(scip, &invmatrix, rank * rank));
            SCIP_CALL( SCIPallocMemoryArray(scip, &vector_a, rank));
            SCIP_CALL( SCIPallocMemoryArray(scip, &submatrix, n * rank));

            /* copy */
            SCIP_CALL( SCIPcblasCopy(newinvmatrix, invmatrix, rank * rank));  /* invmatrix <- newinvmatrix */
            SCIP_CALL( SCIPcblasCopy(vector_newa, vector_a, rank));           /* a <- new a */
            SCIP_CALL( SCIPcblasCopy(newsubmatrix, submatrix, n * rank));     /* submatrix <- newsubmatrix */
            SCIP_CALL( SCIPcblasCopy(vector_newb, vector_b, n));              /* b <- new b */

            /* free */
            SCIPfreeMemoryArrayNull(scip, &newinvmatrix);
            SCIPfreeMemoryArrayNull(scip, &vector_newa);
            SCIPfreeMemoryArrayNull(scip, &newsubmatrix);
            SCIPfreeMemoryArrayNull(scip, &matrix_P);
            SCIPfreeMemoryArrayNull(scip, &vector_q);
            SCIPfreeMemoryArrayNull(scip, &vector_c);
         }
      }
   }

   SCIPfreeMemoryArrayNull(scip, &submatrix);
   SCIPfreeMemoryArrayNull(scip, &invmatrix);
   SCIPfreeMemoryArrayNull(scip, &vector_a);
   SCIPfreeMemoryArrayNull(scip, &vector_b);
   SCIPfreeMemoryArrayNull(scip, &newsubmatrix);
   SCIPfreeMemoryArrayNull(scip, &newinvmatrix);
   SCIPfreeMemoryArrayNull(scip, &vector_newa);
   SCIPfreeMemoryArrayNull(scip, &vector_newb);
   SCIPfreeMemoryArrayNull(scip, &matrix_P);
   SCIPfreeMemoryArrayNull(scip, &vector_q);
   SCIPfreeMemoryArrayNull(scip, &vector_c);

   return SCIP_OKAY;
}


/** find linearly dependent sets from among colmun vectors */
SCIP_RETCODE SCIPgetLineDependSet(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real*      matrix,             /**< symmetric matrix */
   const int             m,                  /**< dimension of matrix */
   const int             ndep,               /**< the number of linearly dependent sets */
   const int*            maxdep,             /**< maximal index at each set */
   const int*            ldindex,            /**< index of linearly dependent sets */
   int*                  depsets             /**< linearly dependent sets */
   )
{
   int i;
   int j;
   int k;
   int t;
   int ct1 = 0;
   int ct2 = 0;
   int rank = 1;
   int dpv;
   SCIP_Real ep = 10E-06;
   SCIP_Real *submatrix = NULL;
   SCIP_Real *vector_d = NULL;
   SCIP_Real *vector_x = NULL;
   SCIP_Real *sort = NULL;
   int *num = NULL;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(m > 1);
   assert(ndep > 0);
   assert(maxdep != NULL);
   assert(ldindex != NULL);
   assert(depsets != NULL);

   for( i = 1; i < m; i++ )
   {
      if( i == maxdep[ct1] && rank > 1 )
      {
         /*
          * solve this linear system: submatrix x = d
          */

         /* allocate memory */
         SCIP_CALL( SCIPallocMemoryArray(scip, &submatrix, rank * rank));
         SCIP_CALL( SCIPallocMemoryArray(scip, &vector_d, rank));
         SCIP_CALL( SCIPallocMemoryArray(scip, &vector_x, rank));
         SCIP_CALL( SCIPallocMemoryArray(scip, &sort, rank));
         SCIP_CALL( SCIPallocMemoryArray(scip, &num, rank));

         /* set submatrix */
         ct2 = 0;
         for( j = 0; j < m; j++ )
         {
            if( ldindex[j] == 0 )
            {
               for( k = 0; k < m; k++ )
               {
                  if( ldindex[k] == 0 )
                     submatrix[ct2++] = SCIPmatColMajor(matrix, m, j, k);
                  else if( k == maxdep[ct1] )
                     break;
               }
            }
            else if( j==maxdep[ct1] )
               break;
         }

         assert(ct2 == rank * rank);

         /* set d */
         ct2 = 0;
         for( j = 0; j < m; j++ )
         {
            if( ldindex[j] == 0 )
               vector_d[ct2++] = SCIPmatColMajor(matrix, m, j, i);
            else if( j == maxdep[ct1] )
               break;
         }

         assert(ct2 == rank);

         /* solve submatrix x = d by using dposv function */
         dpv = SCIPclapackDposv(scip, submatrix, vector_d, rank, vector_x);

         assert(dpv == 0);

         if( dpv != 0 )
         {
            SCIPerrorMessage("dpv = %d\n", dpv);
            return SCIP_ERROR;
         }

         /* square elements of vector x */
         for( j = 0; j < rank; j++ )
            vector_x[j] *= vector_x[j];

         /* sort */
         SCIP_CALL( SCIPsortBubble(scip, rank, vector_x, sort, num));

         /* check whether elements are 0 */
         assert( !EPSEQ(sort[0], 0.0, ep) );
         for( j = 0; j < rank - 1; j++ )
         {
            if( EPSEQ(sort[j+1], 0.0, ep) )
               break;
         }

         /* find linearly dependent sets */
         for( k = 0; k < j + 1; k++ )
         {
            ct2 = 0;
            for( t = 0; t < m; t++ )
            {
               if( ldindex[t] == 0 )
               {
                  if( num[k] == ct2 )
                  {
                     *(depsets + (ct1 * m) + t) = 1;
                     break;
                  }
                  ct2++;
               }
            }
         }

         /* free */
         SCIPfreeMemoryArrayNull(scip, &submatrix);
         SCIPfreeMemoryArrayNull(scip, &vector_d);
         SCIPfreeMemoryArrayNull(scip, &vector_x);
         SCIPfreeMemoryArrayNull(scip, &sort);
         SCIPfreeMemoryArrayNull(scip, &num);
         ct1++;
         if( ct1 == ndep )
            break;
      }
      else if( i==maxdep[ct1] && rank==1 )
      {
         *(depsets+(ct1*m)) = 1;
         ct1++;
      }
      else
      {
         rank++;
      }
   }

   assert(ct1 == ndep);

   return SCIP_OKAY;
}


/** print linearly dependent sets */
SCIP_RETCODE SCIPprintLineDependSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   ndep,               /**< the number of linearly dependent set */
   int                   p,                  /**< the number of explanatory variables */
   int*                  depsets             /**< linearly dependent sets */
   )
{
   int i;
   int j;

   assert(scip != NULL);
   assert(p > 0);
   assert((ndep > 0 && depsets != NULL) || (ndep == 0 && depsets == NULL));

   if( ndep == 0 )
      SCIPinfoMessage(scip, NULL, "The given data has linear independence\n");
   else
   {
      SCIPinfoMessage(scip, NULL, "The given data has linear dependence\n");
      SCIPinfoMessage(scip, NULL, "%d sets of column vector are linearly dependent:\n", ndep);

      if( ndep <= 10 )
      {
         for( i = 0; i < ndep; i++ )
         {
            SCIPinfoMessage(scip, NULL, "{");
            for( j = 0; j < p; j++ )
            {
               if( *(depsets + (i * p) + j) == 1 )
                  SCIPinfoMessage(scip, NULL, " %d", j + 1);
            }
            SCIPinfoMessage(scip, NULL, "}, ");
         }
         SCIPinfoMessage(scip, NULL, "\n");
      }
   }

   return SCIP_OKAY;
}
