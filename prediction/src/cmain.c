
#include <stdio.h>
#include <stdlib.h>
#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"
#include "convenient_tool.h"


/*
 * Local methods
 */

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


static
SCIP_RETCODE calcMean_and_Variance(
   const char*           sampledatafile,
   int                   predict_p,
   SCIP_Real*            mean,
   SCIP_Real*            variance
   )
{

   /* for sampledatafile */
   int n;
   int p;
   int i_ex;
   SCIP_Real* data;                  /**< array to store data */

   /* read sampledatafile */
   SCIP_CALL( readDataDim(sampledatafile, &n, &p, &i_ex));

   assert( p == predict_p );

   /* allocate memory for data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &data, n * (p + 1)));

   /* read data points */
   SCIP_CALL( readData(sampledatafile, n, p, data));

   /* calculate the mean values and the variance values */
   SCIP_CALL( calcMean(n, p, data, mean));
   SCIP_CALL( calcVariance(n, p, data, mean, variance));

   SCIPfreeMemoryArrayNull(scip, &data);

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


/** create variables */
static
SCIP_RETCODE createVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   n,
   int                   p
   )
{
   int i;
   char varname[SCIP_MAXSTRLEN];

   SCIP_VAR** a;
   SCIP_VAR** z;
   SCIP_VAR** ep;
   SCIP_VAR* rss;
   SCIP_VAR* log_rss;

   SCIP_Real para_regterm = 2.0;

   assert(scip != NULL);
   assert(n > 0);
   assert(p > 0);

   /* create variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &a, p));
   SCIP_CALL( SCIPallocMemoryArray(scip, &z, p));
   SCIP_CALL( SCIPallocMemoryArray(scip, &ep, n));

   SCIP_CALL( SCIPcreateVarBasic(scip, &rss, "rss", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS));
   SCIP_CALL( SCIPcreateVarBasic(scip, &log_rss, "log_rss", - SCIPinfinity(scip), SCIPinfinity(scip), (double) n, SCIP_VARTYPE_CONTINUOUS));

   for( i = 0; i < p; i++ )
   {
      /* create continuous variables a_i */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "a%d", i+1);
      SCIP_CALL( SCIPcreateVarBasic( scip, &a[i], varname, - SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS));

      /* create binary variables z_i */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z%d", i+1);
      SCIP_CALL( SCIPcreateVarBasic( scip, &z[i], varname, 0.0, 1.0, para_regterm, SCIP_VARTYPE_BINARY));
   }

   for( i = 0; i < n; i++ )
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "ep%d", i+1);
      SCIP_CALL( SCIPcreateVarBasic( scip, &ep[i], varname, - SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS));
   }

   /* add variables to problem */
   SCIP_CALL( SCIPaddVar(scip, rss));
   SCIP_CALL( SCIPreleaseVar(scip, &rss));
   SCIP_CALL( SCIPaddVar(scip, log_rss));
   SCIP_CALL( SCIPreleaseVar(scip, &log_rss));
   for( i = 0; i < p; i++ )
   {
      SCIP_CALL( SCIPaddVar(scip, a[i]));
      SCIP_CALL( SCIPaddVar(scip, z[i]));
      SCIP_CALL( SCIPreleaseVar(scip, &a[i]));
      SCIP_CALL( SCIPreleaseVar(scip, &z[i]));
   }
   for( i = 0; i < n; i++ )
   {
      SCIP_CALL( SCIPaddVar(scip, ep[i]));
      SCIP_CALL( SCIPreleaseVar(scip, &ep[i]));
   }

   return SCIP_OKAY;
}


static
SCIP_RETCODE writePred(
   const char* outputfile,
   int         TP,
   int         FP,
   int         FN,
   int         TN,
   SCIP_Real   Sensitivity,
   SCIP_Real   Specificity
   )
{
   FILE* file;

   printf("\noutput: %s\n", outputfile);
   file = fopen(outputfile, "w");
   fprintf(file, "%d %d %d %d %f %f", TP, FP, FN, TN, Sensitivity, Specificity);
   fclose(file);

   return SCIP_OKAY;
}


static
SCIP_RETCODE run(
   int                   argc,               /**< number of shell parameters */
   const char*           sampledatafile,
   const char*           predictdatafile,
   const char*           solutionfile,
   const char*           outputfile
   )
{
   SCIP* scip = NULL;

   /* for predictdatafile */
   int n;
   int p;
   int i_ex;
   SCIP_Real* data;                  /**< array to store data */
   SCIP_Real* x;                     /**< array to store data of explanatory variables */
   SCIP_Real* y;

   /* for normalization */
   SCIP_Real* mean;
   SCIP_Real* variance;

   int i,j;

   printf("\n");
   printf("sampledata:  %s\n", sampledatafile);
   printf("predictdata: %s\n", predictdatafile);
   printf("solution:    %s\n", solutionfile);

   /*********
    * Setup *
    *********/
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* read predictdatafile */
   SCIP_CALL( readDataDim(predictdatafile, &n, &p, &i_ex));

   printf("\n");
   printf("predictdata: n = %d, p = %d \n", n, p);

   /* allocate memory for data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &data, n * (p + 1)));

   /* read data points */
   SCIP_CALL( readData(predictdatafile, n, p, data));

   /* calculate mean and variance from sampledata */
   SCIP_CALL( SCIPallocBufferArray(scip, &mean, p+1));
   SCIP_CALL( SCIPallocBufferArray(scip, &variance, p+1));

   SCIP_CALL( calcMean_and_Variance(sampledatafile, p, mean, variance));

   /* normalization */
   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < p + 1; j++ )
      {
         if( !EPSEQ(variance[j], 0.0, 1e-06) )
         {
            *(data + (i * (p + 1)) + j) = (*(data + (i * ( p + 1)) + j) - mean[j]) / sqrt(variance[j]);
         }
      }
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &y, n));
   SCIP_CALL( SCIPallocMemoryArray(scip, &x, n*p));

   /* divide data into explained variable and explanatory variables */
   SCIP_CALL( divideData(n, p, i_ex, data, y, x));

   /* create a problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem"));

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

   /* create and add variables */
   SCIP_CALL( createVariables(scip, n, p));

   {
      SCIP_SOL* readsol;
      SCIP_Bool error;
      SCIP_HEUR* heur;
      char varname[SCIP_MAXSTRLEN];
      SCIP_VAR* var;

      SCIP_Real* solval;
      int* solvalint;

      SCIP_Real sum;
      int TP, FP;
      int FN, TN;
      SCIP_Real Sensitivity;
      SCIP_Real Specificity;

      heur = SCIPfindHeur(scip, "trysol");
      SCIP_CALL( SCIPcreateSol(scip, &readsol, heur));

      /* read solution */
      SCIP_CALL( SCIPreadSolFile(scip, solutionfile, readsol, FALSE, NULL, &error));
      assert( error == FALSE );

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &solval, p));
      SCIP_CALL( SCIPallocBufferArray(scip, &solvalint, p));

      for( i = 0; i < p; i++ )
      {
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "a%d", i + 1);
         var = SCIPfindVar(scip, varname);
         solval[i] = SCIPgetSolVal(scip, readsol, var);

         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z%d", i + 1);
         var = SCIPfindVar(scip, varname);
         solvalint[i] = SCIPgetSolVal(scip, readsol, var);
      }


      TP = 0;
      FP = 0;
      FN = 0;
      TN = 0;

      for( i = 0; i < n; i++ )
      {
         sum = 0.0;
         for( j = 0; j < p; j++ )
         {
            if( solvalint[j] == 1 )
            {
               sum += solval[j] * x[i*p + j];
            }
         }

         if( y[i] < 0.0 )
         {
            /* POSITVE */

            if( sum < 0.0 )
               TP++;
            else
               FN++;
         }
         else
         {
            /* NEGATIVE */

            if( sum > 0.0 )
               TN++;
            else
               FP++;
         }
      }

      printf("\n");
      printf("TP:\t%d\tFP:\t%d\n", TP, FP);
      printf("FN:\t%d\tTN:\t%d\n", FN, TN);

      Sensitivity = (double)TP/(TP+FN);
      Specificity = (double)TN/(FP+TN);

      printf("\n");
      printf("Sensitivity:\t%f\n", Sensitivity);
      printf("Specificity:\t%f\n", Specificity);

      SCIP_CALL( writePred(outputfile, TP, FP, FN, TN, Sensitivity, Specificity));

      SCIPfreeBufferArray(scip, &solval);
      SCIPfreeBufferArray(scip, &solvalint);

   }
   /********************
    * Deinitialization *
    ********************/

   SCIPfreeMemoryArrayNull(scip, &data);
   SCIPfreeMemoryArrayNull(scip, &y);
   SCIPfreeMemoryArrayNull(scip, &x);
   SCIPfreeBufferArray(scip, &mean);
   SCIPfreeBufferArray(scip, &variance);

   //SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                   argc,               /**< number of shell parameters */
   char**                argv                /**< array containing shell parameters */
   )
{
   SCIP_RETCODE retcode;

   if( argc != 5 )
   {
      printf("Usage:prediction sampledata predictdata solution output\n");
      return -1;
   }

   retcode = run(argc, argv[1], argv[2], argv[3], argv[4]);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
