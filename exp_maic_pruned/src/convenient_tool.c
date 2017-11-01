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

/**@file   convenient_tool.c
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "convenient_tool.h"

/** print matrix */
void SCIPprintMat(
   int                   n,                  /**< the number of rows */
   int                   m,                  /**< the number of columns */
   SCIP_Real*            mat                 /**< matrix to print */
   )
{
   int i;
   int j;

   assert(mat != NULL);

   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < m; j++ )
      {
         printf(" %f,", *( mat + ( i * m ) + j ));
      }
      printf("\n");
   }
}


/** print integer matrix */
void SCIPprintIntMat(
   int                   n,                  /**< the number of rows */
   int                   m,                  /**< the number of columns */
   int*                  intmat              /**< integer matrix to print */
   )
{
   int i;
   int j;

   assert(intmat != NULL);

   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < m; j++ )
      {
         printf(" %d,", *( intmat + ( i * m ) + j ));
      }
      printf("\n");
   }
}


/** start a new line */
void SCIPstartNewLine(
   void
   )
{
   printf("\n");
}


/** print a long line */
void SCIPprintLongLine(
   void
   )
{
   printf("*------------------------------------------------------------*\n");
}


/** print a real value */
void SCIPprintVal(
   SCIP_Real             a
   )
{
   printf("%f", a);
   SCIPstartNewLine();
}


/** print a integer value */
void SCIPprintIntVal(
   int                   a
   )
{
   printf("%d", a);
   SCIPstartNewLine();
}


/** print a real vector */
void SCIPprintVec(
   int                   n,                  /**< dimension */
   SCIP_Real*            vec
   )
{
   int i;

   assert(vec != NULL);

   for( i = 0; i < n; i++ )
   {
      printf("%f,", vec[i]);
   }

   SCIPstartNewLine();
}


/** print a integer vector */
void SCIPprintIntVec(
   int                   n,                  /**< dimension */
   int*                  vec
   )
{
   int i;

   assert(vec != NULL);

   for( i = 0; i < n; i++)
   {
      printf("%d , ", vec[i]);
   }

   SCIPstartNewLine();
}


/** calculate the sum of elements in real vector */
SCIP_Real SCIPcalcSum(
   const SCIP_Real*      x,                  /**< real vector to calculate the sum of elements */
   int                   n                   /**< dimension */
   )
{
   int i;
   SCIP_Real y;

   assert(x != NULL);
   assert(n > 0);

   y = 0.0;
   for( i = 0; i < n; i++ )
   {
      y = y + x[i];
   }

   return y;
}


/** calculate the sum of elements in integer vector */
int SCIPcalcIntSum(
   const int*            x,                  /**< integer vector to calculate the sum of elements */
   int                   n                   /**< dimension */
   )
{
   int i;
   int y;

   assert(x != NULL);
   assert(n > 0);

   y = 0;
   for( i = 0; i < n; i++ )
   {
      y = y + x[i];
   }

   return y;
}


/** initialize all array elements to zero */
SCIP_RETCODE SCIPinitArrayZero(
   int                   n,                  /**< dimension */
   SCIP_Real*            x                   /**< array to initialize all elements to zero */
   )
{
   int i;

   assert(x != NULL);
   assert(n > 0);

   for( i = 0; i < n; i++ )
   {
      x[i] = 0.0;
   }

   return SCIP_OKAY;
}


/** initialize all array elements to zero */
SCIP_RETCODE SCIPinitIntArrayZero(
   int                   n,                  /**< dimension */
   int*                  x                   /**< array to initialize all elements to zero */
   )
{
   int i;

   assert(x != NULL);
   assert(n > 0);

   for( i = 0; i < n; i++ )
   {
      x[i] = 0;
   }

   return SCIP_OKAY;
}


/** exit a program */
void SCIPexit(
   void
   )
{
   SCIPstartNewLine();
   printf("SCIPexit!");
   SCIPstartNewLine();
   exit(1);
}


/** execute bubble sort */
SCIP_RETCODE SCIPsortBubble(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   n,                  /**< dimension */
   SCIP_Real*            x,                  /**< array to sort elements */
   SCIP_Real*            y,                  /**< array to return sorted array */
   int*                  z                   /**< array to return indices of sorted array */
   )
{
   int i;
   int j;
   SCIP_Real* a;
   SCIP_Real* b;
   SCIP_Real buf1;
   int buf2;

   assert(scip != NULL);
   assert(n > 0);
   assert(x != NULL);
   assert(y != NULL);
   assert(z != NULL);

   /* alloc */
   SCIP_CALL( SCIPallocBufferArray(scip, &a, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &b, n) );

   assert(a != NULL);
   assert(b != NULL);

   /* copy */
   for( i = 0; i < n; i++ )
   {
      a[i] = x[i];
      b[i] = i;
   }

   /* bubble sort */
   for( i = 0; i < ( n - 1 ); i++ )
   {
      for( j = ( n - 1 ); j > i; j-- )
      {
         if( a[j] > a[j-1] )
         {
            buf1 = a[j];
            a[j] = a[j-1];
            a[j-1] = buf1;

            buf2 = b[j];
            b[j] = b[j-1];
            b[j-1] = buf2;
         }
      }
   }

   /* copy */
   for( i = 0; i < n; i++ )
   {
      *(y + i) = a[i];
      *(z + i) = b[i];
   }

   /* free */
   SCIPfreeBufferArray(scip, &a);
   SCIPfreeBufferArray(scip, &b);

   return SCIP_OKAY;
}


/** transpose a matrix */
SCIP_RETCODE SCIPtransMat(
   int                   n,                  /**< the number of rows */
   int                   m,                  /**< the number of columns */
   SCIP_Real*            inputmat,           /**< matrix to transpose elements */
   SCIP_Real*            outputmat           /**< matrix to return */
   )
{
   int i;
   int j;

   assert(n > 0);
   assert(m > 0);
   assert(inputmat != NULL);
   assert(inputmat != NULL);

   for( i = 0; i < m; i++ )
   {
      for( j = 0; j < n; j++ )
      {
         *( outputmat + ( n * i ) + j ) = *( inputmat + ( j * m ) + i );
      }
   }

   return SCIP_OKAY;
}


/** calculate time */
SCIP_Real SCIPcalcTime(
   clock_t               start,
   clock_t               end
   )
{
   return (SCIP_Real)(end - start)/CLOCKS_PER_SEC;
}
