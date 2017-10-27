#include <stdio.h>
#include <stdlib.h>
#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"
#include "convenient_tool.h"


/*
 * Local methods
 */

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
SCIP_RETCODE readSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   p,
   int*                  counter,
   const char*           solutionfile
   )
{
   SCIP_SOL* readsol;
   SCIP_Bool error;
   SCIP_HEUR* heur;
   char varname[SCIP_MAXSTRLEN];
   SCIP_VAR* var;
   int ct;
   int i;

   heur = SCIPfindHeur(scip, "trysol");
   SCIP_CALL( SCIPcreateSol(scip, &readsol, heur));

   /* read solution */
   SCIP_CALL( SCIPreadSolFile(scip, solutionfile, readsol, FALSE, NULL, &error));
   assert( error == FALSE );

   ct = 0;
   for( i = 0; i < p; i++ )
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z%d", i + 1);
      var = SCIPfindVar(scip, varname);
      if( SCIPisEQ(scip, SCIPgetSolVal(scip, readsol, var), 1.0) )
      {
         counter[i]++;
         ct++;
      }

      if( ct == 10 )
         break;
   }

   assert( ct == 10 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE writePred(
   const char*           outputfile,
   int                   p,
   int*                  counter
   )
{
   FILE* file;
   int i;

   printf("\noutput: %s\n", outputfile);
   file = fopen(outputfile, "w");

   for( i = 0; i < p; i++ )
   {
      fprintf(file, "%d %d\n", i, counter[i]);
   }

   fclose(file);

   return SCIP_OKAY;
}


static
SCIP_RETCODE run(
   int                   argc,               /**< number of shell parameters */
   const char*           sampledatafile,
   const char*           solutiondir,
   const char*           outputfile
   )
{
   SCIP* scip = NULL;

   /* for sampledatafile */
   int n;
   int p;
   int i_ex;

   int i,j;
   int* counter;

   printf("\n");
   printf("sampledata:  %s\n", sampledatafile);

   /*********
    * Setup *
    *********/
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* read predictdatafile */
   SCIP_CALL( readDataDim(sampledatafile, &n, &p, &i_ex));

   printf("\n");
   printf("data: n = %d, p = %d \n", n, p);

   SCIP_CALL( SCIPallocBufferArray(scip, &counter, p));
   for( i = 0; i < p; i++ )
      counter[i] = 0;

   /* create a problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem"));

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

   /* create and add variables */
   SCIP_CALL( createVariables(scip, n, p));

   /* read solutions */
   {
      char solutionfile[SCIP_MAXSTRLEN];
      for( i = 0; i < 10; i++ )
      {
         for( j = 0; j < 5; j++ )
         {
            FILE *File;
            (void) SCIPsnprintf(solutionfile, SCIP_MAXSTRLEN, "%s/%d/%d_sample_%s.sol", solutiondir, i, j, solutiondir);

            if( (File = fopen(solutionfile, "r")) != NULL)
               SCIP_CALL( readSolution( scip, p, counter, solutionfile));

            fclose(File);
         }
      }
   }

   SCIP_CALL( writePred(outputfile, p, counter));
   /********************
    * Deinitialization *
    ********************/

   SCIPfreeBufferArray(scip, &counter);
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

   if( argc != 4 )
   {
      printf("Usage:datafile dir_solution outputfile\n");
      return -1;
   }

   retcode = run(argc, argv[1], argv[2], argv[3]);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
