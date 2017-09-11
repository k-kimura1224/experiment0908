/* main.c */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "scip/scip.h"

/** exit a program if a input file has error */
static
SCIP_RETCODE readFail(
   int                   status
   )
{
   if( !status )
   {
      printf("Reading failed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/** read dimension from the input file */
static
SCIP_RETCODE readDataDim(
   const char            *filename,          /**< name of the input file */
   int*                  n,                  /**< pointer to store the number of data points*/
   int*                  p,                  /**< pointer to store the  number of explanatory variables */
   int*                  i_ex,               /**< pointer to store the index of explained varaible */
   char*                 firstline
   )
{
   FILE *file;
   int status;

   assert(filename != NULL);
   assert(n != NULL);
   assert(p != NULL);
   assert(i_ex != NULL);

   /* open file */
   file = fopen(filename, "r");
   if( file == NULL )
   {
      printf("Could not open file <%s>.\n", filename);
      return SCIP_ERROR;
   }

   /* skip one line */
   if( fgets(firstline, 1000, file) == NULL )
   {
      printf("Error reading file <%s>.\n", filename);
      return SCIP_ERROR;
   }

   /* read n */
   status = fscanf(file, "%d", n);
   SCIP_CALL( readFail(status));

   /* read p */
   status = fscanf(file, "%d", p);
   SCIP_CALL( readFail(status));

   /* read i_ex */
   status = fscanf(file, "%d", i_ex);
   SCIP_CALL( readFail(status));

   /* close file */
   fclose(file);

   return SCIP_OKAY;
}


/** read data points */
static
SCIP_RETCODE readData(
   const char*           filename,           /**< name of the input file */
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory varaibles */
   SCIP_Real*            data                /**< array to store data points */
   )
{
   int i;
   FILE *file;
   int buf;
   int status;
   char s[1000];

   /* open file */
   file = fopen(filename, "r");
   if( file==NULL )
   {
      printf("Could not open file <%s>.\n", filename);
      return SCIP_ERROR;
   }

   /* skip one line */
   if( fgets(s, 1000, file) == NULL )
   {
      printf("Error reading file <%s>.\n", filename);
      return SCIP_ERROR;
   }

   /* skip n */
   status = fscanf(file, "%d", &buf);
   SCIP_CALL( readFail(status));

   /* skip p */
   status = fscanf(file, "%d", &buf);
   SCIP_CALL( readFail(status));

   /* skip i_ex */
   status = fscanf(file, "%d", &buf);
   SCIP_CALL( readFail(status));

   /* read data points */
   for( i = 0; i < n * ( p + 1 ); i++ )
   {
      status = fscanf(file, "%lf", (data + i));
      SCIP_CALL( readFail(status));
   }

   /* close file */
   fclose(file);

   return SCIP_OKAY;
}


int
main(int argc,char *argv[])
{


   int n;
   int p;
   int i_ex;
   char firstline[1000];
   SCIP_Real* data;
   int* list;
   int* counter;
   int dim;

   SCIP* scip = NULL;
   int i;
   int j;
   int k;

   int s;
   int max;

   int buf;

   if( argc != 5 )
   {
      printf("Error: commandline arguments \n");
      return -1;
   }

   s = atoi(argv[2]);

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* read dimension of data points */
   SCIP_CALL( readDataDim(argv[1], &n, &p, &i_ex, firstline));

   if( s < 0 || s > n )
   {
      printf("Error: commandline arguments \n");
      return -1;
   }

   /* allocate memory for data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &data, n * (p+1)));
   /* read data points */
   SCIP_CALL( readData(argv[1], n, p, data));

   SCIP_CALL( SCIPallocMemoryArray(scip, &list, n));
   SCIP_CALL( SCIPallocMemoryArray(scip, &counter, s));

   for( i = 0; i < s; i++ )
      counter[i] = 0;

   max = floor((double)n / (double)s + 1.0e-08);

   assert( ( max * s ) + ( n % s ) == n );

   for( j = 0; j < n; j++ )
      list[j] = -1;

   for( j = 0; j < s; j++ )
      counter[j] = 0;

   srand(time(NULL));

   for( j = 0; j < max * s; j++ )
   {
      while( 1 )
      {
         buf = rand() % s;

         if( counter[buf] < max )
         {
            list[j] = buf;
            counter[buf]++;
            break;
         }
      }
   }

   for( j = max * s; j < n; j++ )
   {
      while( 1 )
      {
         buf = rand() % s;

         assert( counter[buf] == max || counter[buf] == max + 1 );
         if( counter[buf] == max )
         {
            list[j] = buf;
            counter[buf]++;
            break;
         }
      }
   }

   for( j = 0; j < n; j++ )
      printf("%d ", list[j]);

   printf("\n");

   for( i = 0; i < s; i++ )
   {
      FILE* file_sample;
      FILE* file_predict;
      char filename[SCIP_MAXSTRLEN];

      dim = 0;
      for( j = 0; j < s; j++ )
      {
         if( i != j )
            dim += counter[j];
      }

      (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%d_%s", i, argv[3]);
      printf("sample file:     %s (%d)\n", filename, dim);
      file_sample = fopen(filename, "w");

      fprintf(file_sample, "%s", firstline);
      fprintf(file_sample, "%d\n", dim);
      fprintf(file_sample, "%d\n", p);
      fprintf(file_sample, "%d\n", i_ex);

      for( j = 0; j < n; j++ )
      {
         if( list[j] != i )
         {
            for( k = 0; k < p; k++ )
            {
               fprintf(file_sample, "%d ", (int) data[j*(p+1) + k]);
            }
            fprintf(file_sample, "%d\n", (int) data[j*(p+1) + p]);
         }
      }

      fclose(file_sample);

      (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%d_%s", i, argv[4]);
      printf("prediction file: %s (%d)\n", filename, counter[i]);
      file_predict = fopen(filename, "w");

      assert( dim + counter[i] == n );
      fprintf(file_predict, "%s", firstline);
      fprintf(file_predict, "%d\n", counter[i]);
      fprintf(file_predict, "%d\n", p);
      fprintf(file_predict, "%d\n", i_ex);

      for( j = 0; j < n; j++ )
      {
         if( list[j] == i )
         {
            for( k = 0; k < p; k++ )
            {
               fprintf(file_predict, "%d ", (int) data[j*(p+1) + k]);
            }
            fprintf(file_predict, "%d\n", (int) data[j*(p+1) + p]);
         }
      }

      fclose(file_predict);

      printf("\n");
   }

   sleep(1);

   SCIPfreeMemoryArrayNull(scip, &data);
   SCIPfreeMemoryArrayNull(scip, &list);
   SCIPfreeMemoryArrayNull(scip, &counter);
   SCIP_CALL( SCIPfree(&scip) );
   BMScheckEmptyMemory();

   return 0;
}

