/* main.c */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

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
   int*                  i_ex                /**< pointer to store the index of explained varaible */
   )
{
   FILE *file;
   int status;
   char s[1000];

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
   if( fgets(s, 1000, file) == NULL )
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


/** calculate the mean values for normalization */
static
SCIP_RETCODE calcMean(
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory varaibles */
   SCIP_Real*            data,               /**< data points */
   SCIP_Real*            mean                /**< array to store the mean values */
   )
{
   int i;
   int j;

   assert(n > 0);
   assert(p > 0);
   assert(data != NULL);
   assert(mean != NULL);

   for( i = 0; i < p + 1; i++ )
   {
      *(mean + i) = 0.0;
      for( j = 0; j < n; j++ )
      {
         *(mean + i) += *(data + (j * (p + 1)) + i);
      }
      *(mean + i) *= 1.0 / (double) n;
   }

   return SCIP_OKAY;
}


/** calculate the variance values for normalization */
static
SCIP_RETCODE calcVariance(
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory varaibles */
   SCIP_Real*            data,               /**< data points */
   SCIP_Real*            mean,               /**< the mean values */
   SCIP_Real*            variance            /**< array to store the variance values */
   )
{
   int i;
   int j;

   assert(n > 1);
   assert(p > 0);
   assert(data != NULL);
   assert(mean != NULL);
   assert(variance != NULL);

   for( i = 0; i < p + 1; i++ )
   {
      *(variance + i) = 0.0;
      for( j = 0; j < n; j++ )
      {
         *(variance + i) += pow(*(data + (j * (p + 1)) + i) - mean[i], 2.0);
      }
      *(variance + i) *= 1.0 / ((double) n - 1.0);

   }

   return SCIP_OKAY;
}


/** normalize data */
static
SCIP_RETCODE normalization(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory variables */
   SCIP_Real*            data                /**< data points */
   )
{
   int i;
   int j;
   SCIP_Real* mean;
   SCIP_Real* variance;

   assert(scip != NULL);
   assert(n > 0);
   assert(p > 0);
   assert(data != NULL);

   /* allocate memory for mean and variance */
   SCIP_CALL( SCIPallocMemoryArray(scip, &mean, p + 1));
   SCIP_CALL( SCIPallocMemoryArray(scip, &variance, p + 1));

   /* calculate the mean values and the variance values */
   SCIP_CALL( calcMean(n, p, data, mean));
   SCIP_CALL( calcVariance(n, p, data, mean, variance));

   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < p + 1; j++ )
      {
         assert(!EPSEQ(variance[j], 0.0, 1e-06));
         *(data + (i * (p + 1)) + j) = (*(data + (i * ( p + 1)) + j) - mean[j]) / sqrt(variance[j]);
      }
   }

   /* free */
   SCIPfreeMemoryArrayNull(scip, &mean);
   SCIPfreeMemoryArrayNull(scip, &variance);

   return SCIP_OKAY;
}


/** divide data into explained variable and explanatory variables */
static
SCIP_RETCODE divideData(
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory variables */
   int                   i_ex,               /**< index of the explained variable */
   SCIP_Real*            data,               /**< data points */
   SCIP_Real*            explained,          /**< array to store the explained variable */
   SCIP_Real*            explanatory         /**< array to store the explanatory variables */
   )
{
   int i;
   int j;
   int ct = 0;

   assert(n > 0);
   assert(p > 0);
   assert(i_ex > 0 && i_ex <= p);
   assert(data != NULL);
   assert(explained != NULL);
   assert(explanatory != NULL);

   for( i = 0; i < n; i++ )
   {
      *(explained + i) = *(data + (i * (p + 1)) + (i_ex - 1));
   }

   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < p + 1; j++ )
      {
         if( j != i_ex - 1 )
         {
            *(explanatory + ct) = *(data + (i * (p + 1)) + j);
            ct++;
         }
      }
   }

   return SCIP_OKAY;
}


static
void objectfunction(
   int n
   )
{
   printf("minimize\n");
   printf("[ 2 ep1^2 ");

   int   i;
   for(i=1; i<n; ++i)
      printf("+ 2 ep%d^2 ", i+1);

   printf("] /2\n");

}


static
void c(
   int   ct
   )
{
   printf("c%d:  ", ct);
}


// ep_i + b + sum{ a_j*x_ij ) = y_i
static
void epsilon(
   int      i,
   int      p,
   double   y,
   double*  x
   )
{
    int   j;

    printf("ep%d + a%d", i+1, p+1);

    for(j=0; j<p; ++j){
       if( x[j]*x[j] > 0.00001 ){
          if( x[j] > 0.0 )
             printf(" + %.8f a%d", x[j], j+1);
          else
             printf(" - %.8f a%d", -x[j], j+1);
       }
    }

    if( y*y > 0.00001 ){
       if( y > 0.0 )
          printf(" = %.8f\n", y);
       else
          printf(" = - %.8f\n", -y);
    }else{
       printf(" = 0\n");
    }

}


// z_j = 0 -> a_j = 0
static
void indicator(
   int j
   )
{
   printf("z%d = 0 -> a%d = 0\n", j+1, j+1);
}


// sum{ z_j } = k
static
void zconst(
   int   p,
   int   k
   )
{
   int j;

   for(j=0; j<(p-1); ++j)
      printf("z%d + ", j+1);

   printf("z%d ", j+1);

   //printf("= %d\n", k);
   printf("<= %d\n", k);

}


static
void constraints(
   int      n,
   int      p,
   double*  y,
   double*  X,
   int      k
   )
{
   int   i,j;
   int   ct = 1;

   printf("subject to\n");

   // epsilon const
   for(i=0; i<n; ++i){
      c(ct); ct++;
      epsilon( i, p, y[i], &X[i*p]);
   }

   // indicator constraints
   for(j=0; j<p+1; j++){
      c(ct); ct++;
      indicator(j);
   }

   // zconst
   c(ct); ct++;
   zconst(p, k);

}


static
void boundary(
   int   n,
   int   p,
   int   k
   )
{
   int i,j;

   printf("bounds\n");

   // free
   for(j=0; j<(p+1); ++j)
      printf("a%d free\n", j+1);
   for(i=0; i<n; ++i)
      printf("ep%d free\n", i+1);

   // binary
   printf("binary\n");
   for(j=0; j<(p+1); ++j)
      printf("z%d\n", j+1);

}


int
main(int argc,char *argv[])
{


   int n;
   int p;
   int i_ex;
   SCIP_Real* data;
   SCIP_Real* explained;
   SCIP_Real* explanatory;

   SCIP* scip = NULL;
   int k;

   if( argc != 3 )
   {
      printf("Error: commandline arguments \n");
      return -1;
   }

   k = atoi(argv[2]);

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* read dimension of data points */
   SCIP_CALL( readDataDim(argv[1], &n, &p, &i_ex));

   if( k < 0 || k > p )
   {
      printf("Error: commandline arguments \n");
      return -1;
   }

   /* allocate memory for data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &data, n * (p+1)));
   /* read data points */
   SCIP_CALL( readData(argv[1], n, p, data));

   /* normalize data */
   SCIP_CALL( normalization(scip, n, p, data));

   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &explained, n));
   SCIP_CALL( SCIPallocMemoryArray(scip, &explanatory, n*p));

   /* divide data into explained variable and explanatory variables */
   SCIP_CALL( divideData(n, p, i_ex, data, explained, explanatory));

   // output lp format {{
   objectfunction(n);
   constraints(n, p, explained, explanatory, k);
   boundary(n, p, k);
   // }} output lp format

   printf("end\n");

   SCIPfreeMemoryArrayNull(scip, &data);
   SCIPfreeMemoryArrayNull(scip, &explained);
   SCIPfreeMemoryArrayNull(scip, &explanatory);
   SCIP_CALL( SCIPfree(&scip) );
   BMScheckEmptyMemory();

   return 0;
}

