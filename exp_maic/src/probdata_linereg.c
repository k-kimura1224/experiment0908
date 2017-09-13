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

/**@file   probdata_linereg.c
 * @brief  problem data for minimization of AIC/BIC
 * @author Keiji Kimura
 *
 * This file implements the problem data for minimization of AIC/BIC
 *
 * The problem data contains the number of data, the number of explanatory variables, and data points.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "probdata_linereg.h"
#include "convenient_tool.h"
#include "call_cblas.h"
#include "get_lineardependence.h"

struct SCIP_ProbData
{
   /* for FiberSCIP */
   SCIP_Bool             ug;                 /**< indicates if this ug dual bound is set or not */
   int                   nSolvers;           /**< the number of solvers */
   SCIP_Real             ugDual;             /**< dual bound set by ug */

   /* data obtained from the input file */
   int                   n;                  /**< the number of data points */
   int                   p;                  /**< the number of explanatory variables */
   SCIP_Real*            y;                  /**< array for explained variable */
   SCIP_Real*            x;                  /**< array for explanatory variables */

   /* for computation */
   SCIP_Real*            xx;                 /**< this is defined as x'x */
   SCIP_Real*            xy;                 /**< this is defined as x'y */
   SCIP_Real             yy;                 /**< this is defined as y'y */
   int                   ndep;               /**< the number of linearly dependent sets */
   int*                  indexdepsets;       /**< array for linearly dependent sets */
   int                   sizealldepsets;     /**< size of indexdepsets */
   int*                  sizedepsets;        /**< sizes of linearly dependent sets */

   /* variables */
   SCIP_VAR**            a;                  /**< continuous variables */
   SCIP_VAR**            z;                  /**< binary variables */
   SCIP_VAR**            ep;                 /**< continuous variables */
   SCIP_VAR*             rss;                /**< continuous variable, residual sum of squares */
   SCIP_VAR*             log_rss;            /**< continuous variable, log(rss) */
   int                   nvars;              /**< the number of variables */

   /* parameters to control the regularization term */
   char                  mode;               /**< range: {ab} default: a */
   SCIP_Real             para_regterm;       /**< If mode is aic, this is 2. If bic, log(n). */

   /* parameter to control the number of selected explanatory variables */
   int                   ubnumex;            /**< range: [-1,100], default: -1 */

   /* parameter to control branching rules */
   char                  bmode;              /**< range: {asm} default: a */

   /* parameter to save cpu memory */
   SCIP_Bool             savememory;
};


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


/** check data */
static
SCIP_RETCODE checkData(
   int                   n,                  /**< the number of data points */
   int                   p,                  /**< the number of explanatory varaibles */
   int                   i_ex,               /**< index of explained varaible */
   SCIP_Real*            data                /**< data points */
   )
{
   int i;
   int buf;
   int p1;

   assert(n > 0);
   assert(p > 0);
   assert(data != NULL);
   assert(i_ex > 0 && i_ex <= p + 1);

   buf = i_ex - 1;
   p1 = p + 1;

   for( i = 1; i < n; i++ )
   {
      if( !EPSEQ(data[buf + i * p1], data[buf + (i-1) * p1], 1e-06) )
            break;
   }

   if( i == n )
   {
      SCIPerrorMessage("error in checkData\n");
      return SCIP_ERROR;
   }

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
         if( !EPSEQ(variance[j], 0.0, 1e-06 ) )
         {
            *(data + (i * (p + 1)) + j) = (*(data + (i * ( p + 1)) + j) - mean[j]) / sqrt(variance[j]);
         }
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


/** get a problem name */
static
SCIP_RETCODE getProblemName(
   const char*           filename,           /**< filename */
   char*                 probname,           /**< array to store probname */
   int                   maxsize             /* maximum size of probname */
   )
{
   int i;
   int j;
   int l;

   assert(filename != NULL);
   assert(probname != NULL);

   i = 0;
   j = 0;

   /* first find end of string */
   while( filename[i]!=0 )
      ++i;
   l = i;

   /* go back until '.' or '/' or '\' appears */
   while( (i > 0) && (filename[i] != '.') && (filename[i] != '/') && (filename[i] != '\\'))
      --i;

   /* if we found '.', search for '/' or '\\' */
   if( filename[i] == '.' )
   {
      l = i;
      while( (i > 0) && (filename[i] != '/') && (filename[i] != '\\') )
         --i;
   }

   /* crrect counter */
   if( (filename[i] == '/') || (filename[i] == '\\') )
      ++i;

   /* copy name */
   while( (i < l) && (filename[i] != 0) )
   {
      probname[j++] = filename[i++];
      if( j>maxsize-1)
         return SCIP_ERROR;
   }
   probname[j] = 0;

   return SCIP_OKAY;
}


/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   assert(scip != NULL);
   assert(probdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, probdata) );

   /* set parameters to control the regularization term */
   SCIP_CALL( SCIPgetCharParam(scip, "linereg/mode", &((*probdata)->mode)) );
   SCIP_CALL( SCIPgetCharParam(scip, "linereg/bmode", &((*probdata)->bmode)) );
   SCIP_CALL( SCIPgetIntParam(scip, "linereg/ubnumex", &((*probdata)->ubnumex)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "linereg/savememory", &((*probdata)->savememory)) );
   assert((*probdata)->mode == 'a' || (*probdata)->mode == 'b');
   assert((*probdata)->bmode == 'a' || (*probdata)->bmode == 's' || (*probdata)->bmode == 'm' );
   assert((*probdata)->ubnumex >= -1 );

   /* set parameters for ug */
   (*probdata)->ug = FALSE;
   (*probdata)->nSolvers =0;
   (*probdata)->ugDual = 0.0;

   return SCIP_OKAY;
}


/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   int i;
   int n;
   int p;
   n = (*probdata)->n;
   p = (*probdata)->p;

   SCIPdebugMessage("probdataFree \n");
   assert(scip != NULL);
   assert(probdata != NULL);
   assert(n > 0);
   assert(p > 0);

   /* free memory */
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->y);
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->x);
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->xy);
   if( (*probdata)->savememory == FALSE )
      SCIPfreeMemoryArrayNull(scip, &(*probdata)->xx);

   if( (*probdata)->ndep )
   {
      SCIPfreeMemoryArrayNull(scip, &(*probdata)->indexdepsets);
      SCIPfreeMemoryArrayNull(scip, &(*probdata)->sizedepsets);
   }

   /* release variables */
   for( i = 0; i < p; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->a[i]));
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->z[i]));
   }
   for( i = 0; i < n; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->ep[i]));
   }
   SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->rss));
   SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->log_rss));

   SCIPfreeMemoryArrayNull(scip, &(*probdata)->a);
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->z);
   SCIPfreeMemoryArrayNull(scip, &(*probdata)->ep);

   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


/** create constraints */
static
SCIP_RETCODE createConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   int n;
   int p;
   int ndep;
   int ubnumex;
   int* sizedepsets;
   int* indexdepsets;

   SCIP_VAR** var_a;
   SCIP_VAR** var_z;
   SCIP_VAR** var_ep;
   SCIP_VAR* var_rss;
   SCIP_VAR* var_log;

   SCIP_CONS* defrss;
   SCIP_CONS* deflog;
   SCIP_CONS** defep;
   SCIP_CONS** amcons1;
   SCIP_CONS** amcons2;
   char consname[SCIP_MAXSTRLEN];
   SCIP_Real one;
   SCIP_Real minusone;
   SCIP_Real bigm;

   int i;
   int j;

   n = probdata->n;
   p = probdata->p;
   ndep = probdata->ndep;
   ubnumex = probdata->ubnumex;

   one = 1.0;
   minusone = - 1.0;
   bigm = 100.0;

   if( ndep > 0 )
   {
      sizedepsets = probdata->sizedepsets;
      indexdepsets = probdata->indexdepsets;
   }
   else
   {
      sizedepsets = NULL;
      indexdepsets = NULL;
   }

   var_a = probdata->a;
   var_z = probdata->z;
   var_ep = probdata->ep;
   var_rss = probdata->rss;
   var_log = probdata->log_rss;

   SCIPdebugMessage("createConstraints \n");
   assert(scip != NULL);
   assert(probdata != NULL);
   assert(n > 0);
   assert(p > 0);
   assert(ndep >= 0);
   assert( n > p || (p >= n && ubnumex > 0 && ubnumex < p) );
   assert( (ndep == 0 && sizedepsets == NULL) || (ndep > 0 && sizedepsets != NULL) );
   assert(var_a != NULL);
   assert(var_z != NULL);
   assert(var_ep != NULL);
   assert(var_rss != NULL);
   assert(var_log != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &defep, n));
   SCIP_CALL( SCIPallocMemoryArray(scip, &amcons1, p));
   SCIP_CALL( SCIPallocMemoryArray(scip, &amcons2, p));

   /**
    * create the quadratic constraint defrss:
    * ep_1^2 + .. + ep_n^2 - rss = 0.0
    */
   {
      SCIP_Real* onen;

      /* allocate memory */
      SCIP_CALL( SCIPallocMemoryArray(scip, &onen, n));

      for( i = 0; i < n; i++ )
         onen[i] = 1.0;

      /* create the quadratic constraint */
      SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &defrss, "def_rss", 1, &var_rss, &minusone, n, var_ep, var_ep, onen, 0.0, 0.0));

      /* free */
      SCIPfreeMemoryArrayNull(scip, &onen);
   }

   /**
    * create teh nonlinear constraint deflog:
    * log( rss ) - log_rss = 0.0
    */
   {
      SCIP_EXPR* rssexpr;
      SCIP_EXPR* expr;
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Real coef;
      SCIP_EXPRTREE* exprtree;

      var1 = var_rss;
      var2 = var_log;
      coef = -1;

      /* setup expression */
      SCIP_CALL( SCIPexprCreate( SCIPblkmem(scip), &rssexpr, SCIP_EXPR_VARIDX, 0));

      /* expression for expr : log( rss ) */
      SCIP_CALL( SCIPexprCreate( SCIPblkmem(scip), &expr, SCIP_EXPR_LOG, rssexpr));

      /* expression tree from expr */
      SCIP_CALL( SCIPexprtreeCreate( SCIPblkmem(scip), &exprtree, expr, 1, 0, NULL));
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, 1, &var1));

      /* create the nonlinear constraint for exprtree - log_rss = 0.0 */
      SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &deflog, "def_log", 1, &var2, &coef, 1, &exprtree, &one, 0.0, 0.0));

      SCIP_CALL( SCIPexprtreeFree(&exprtree));

   }

   /**
    * create the linear constraints defep:
    * y_i - \sum_{j=1}^p a_j x_ij - ep_i = 0
    */
   {
      SCIP_Real* coef;

      /* allocate memory */
      SCIP_CALL( SCIPallocMemoryArray(scip, &coef, p));

      for( i = 0; i < n; i++ )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "def_ep%d", i+1);

         for( j = 0; j < p; j++ )
            coef[j] = probdata->x[ i + ( n * j ) ];

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &defep[i], consname, p, var_a, coef, probdata->y[i], probdata->y[i]));
         SCIP_CALL( SCIPaddCoefLinear(scip, defep[i], var_ep[i], one));
      }

      /* free */
      SCIPfreeMemoryArrayNull(scip, &coef);
   }

   /**
    * create the linear constraints amcons1:
    * 0 <= a_j + ( bigm ) z_j
    */
   {
      for( i = 0; i < p; i++ )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "amcons1_%d", i+1);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &amcons1[i], consname, 0, NULL, NULL, 0.0, SCIPinfinity(scip)));
         SCIP_CALL( SCIPaddCoefLinear(scip, amcons1[i], var_a[i], one));
         SCIP_CALL( SCIPaddCoefLinear(scip, amcons1[i], var_z[i], bigm));
      }
   }

   /**
    * create the linear constraints amcons2:
    * a_j - ( bigm ) z_j <= 0
    */
   {
      for( i = 0; i < p; i++ )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "amcons2_%d", i+1);
         SCIP_CALL( SCIPcreateConsBasicLinear( scip, &amcons2[i], consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0));
         SCIP_CALL( SCIPaddCoefLinear( scip, amcons2[i], var_a[i], one));
         SCIP_CALL( SCIPaddCoefLinear( scip, amcons2[i], var_z[i], - bigm));
      }
   }

   /* add constraints to problem and relsease them */
   SCIP_CALL( SCIPaddCons(scip, defrss));
   SCIP_CALL( SCIPaddCons(scip, deflog));
   SCIP_CALL( SCIPreleaseCons(scip, &defrss));
   SCIP_CALL( SCIPreleaseCons(scip, &deflog));

   for( i = 0; i < n; i++ )
   {
      SCIP_CALL( SCIPaddCons(scip, defep[i]));
      SCIP_CALL( SCIPreleaseCons(scip, &defep[i]));
   }

   for( i = 0; i < p; i++ )
   {
      SCIP_CALL( SCIPaddCons(scip, amcons1[i]));
      SCIP_CALL( SCIPreleaseCons(scip, &amcons1[i]));

      SCIP_CALL( SCIPaddCons(scip, amcons2[i]));
      SCIP_CALL( SCIPreleaseCons(scip, &amcons2[i]));
   }

   /**
    * create linearly dependent constraints
    * ex.
    *  linearly dependent set = { 1, 2, 3, 4}
    *  add this constraint
    *    0 <= z1 + z2 + z3 + z4 <= 3
    */
   if( ndep )
   {
      int* ip = &indexdepsets[0];
      int size;
      int buf = 0;
      for( i = 0; i < ndep; i++ )
      {
         size = sizedepsets[i];
         ip += buf;
         if( (ubnumex == -1) || (ubnumex > 0 && size - 1 < ubnumex) )
         {
            SCIP_CONS* ld;
            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "ld_%d", i+1);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &ld, consname, 0, NULL, NULL, 0.0, size-1));

            for( j = 0; j < size; j++ )
               SCIP_CALL( SCIPaddCoefLinear(scip, ld, var_z[*(ip+j)], one));

            SCIP_CALL( SCIPaddCons(scip, ld));
            SCIP_CALL( SCIPreleaseCons(scip, &ld));
         }
         buf = size;
      }
   }

   /* add a linear constraint to control the number of explanatory variables */
   if( ubnumex != -1 )
   {
      SCIP_CONS* ubnumexcons;
      assert(ubnumex >= 1 || ubnumex <= p);

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &ubnumexcons, "ubnumex", 0, NULL, NULL, 0.0, (SCIP_Real) ubnumex));

      for( i = 0; i < p; i++ )
         SCIP_CALL( SCIPaddCoefLinear(scip, ubnumexcons, var_z[i], one));

      SCIP_CALL( SCIPaddCons(scip, ubnumexcons));
      SCIP_CALL( SCIPreleaseCons(scip, &ubnumexcons));

   }

   /* free */
   SCIPfreeMemoryArrayNull(scip, &defep);
   SCIPfreeMemoryArrayNull(scip, &amcons1);
   SCIPfreeMemoryArrayNull(scip, &amcons2);

   return SCIP_OKAY;
}


/** create variables */
static
SCIP_RETCODE createVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   int i;
   int j;
   int in;
   int n;
   int p;
   SCIP_Real para_regterm;
   char varname[SCIP_MAXSTRLEN];
   SCIP_Real* data_x;
   SCIP_Bool check;

   n = probdata->n;
   p = probdata->p;
   para_regterm = probdata->para_regterm;
   data_x = probdata->x;

   assert(scip != NULL);
   assert(probdata != NULL);
   assert(n > 0);
   assert(p > 0);
   assert(EPSEQ(para_regterm, 2.0, 1e-08) || EPSEQ(para_regterm, log((SCIP_Real) n), 1e-08));

   /* create variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->a, p));
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->z, p));
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->ep, n));

   SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->rss, "rss", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS));
   SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->log_rss, "log_rss", - SCIPinfinity(scip), SCIPinfinity(scip), (double) n, SCIP_VARTYPE_CONTINUOUS));

   for( i = 0; i < p; i++ )
   {
      /* create continuous variables a_i */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "a%d", i+1);
      SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->a[i], varname, - SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS));

      /* check data */
      check = TRUE;
      in = i * n;
      for( j = 1; j < n; j++ )
      {
         if( !EPSEQ(data_x[in + j], data_x[in + (j-1)], 1e-06) )
            break;
      }

      if( j == n )
         check = FALSE;

      /* create binary variables z_i */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "z%d", i+1);
      if( check == TRUE )
      {
         SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->z[i], varname, 0.0, 1.0, para_regterm, SCIP_VARTYPE_BINARY));
      }
      else
      {
         SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->z[i], varname, 0.0, 0.0, para_regterm, SCIP_VARTYPE_BINARY));
      }
   }

   for( i = 0; i < n; i++ )
   {
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "ep%d", i+1);
      SCIP_CALL( SCIPcreateVarBasic( scip, &probdata->ep[i], varname, - SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS));
   }

   /* add variables to problem */
   SCIP_CALL( SCIPaddVar(scip, probdata->rss));
   SCIP_CALL( SCIPaddVar(scip, probdata->log_rss));
   for( i = 0; i < p; i++ )
   {
      SCIP_CALL( SCIPaddVar(scip, probdata->a[i]));
      SCIP_CALL( SCIPaddVar(scip, probdata->z[i]));
   }
   for( i = 0; i < n; i++ )
      SCIP_CALL( SCIPaddVar(scip, probdata->ep[i]));

   probdata->nvars = 2 + p + p + n;

   /* change priorities of branching variables */
   for( i = 0; i < p; i++ )
   {
      SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->z[i], 1000));
      SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->a[i], -1000));
   }
   for( i = 0; i < n; i++ )
      SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->ep[i], -1000));
   SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->rss, -1000));
   SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->log_rss, -1000));

   return SCIP_OKAY;
}


/*
 * Callback methods of probdata
 */

/** copies user data of source SCIP for the target SCIP */
static
SCIP_DECL_PROBCOPY(probcopyLinereg)
{
   int n;
   int p;
   int ndep;
   int size_all;
   int np;
   int pp;
   SCIP_Bool success;
   int i;

   n = sourcedata->n;
   p = sourcedata->p;
   ndep = sourcedata->ndep;
   size_all = sourcedata->sizealldepsets;

   np = n * p;
   pp = p * p;

   assert(scip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);
   assert(n > 0);
   assert(p > 0);
   assert( (ndep == 0 && size_all == 0) || (ndep > 0 && size_all > 0) );

   SCIPdebugMessage("########################## probcopy ###########################\n");

   /* allocate memory */
   SCIP_CALL( probdataCreate(scip, targetdata) );

   /* copy */
   (*targetdata)->n = sourcedata->n;
   (*targetdata)->p = sourcedata->p;
   (*targetdata)->ndep = sourcedata->ndep;
   (*targetdata)->nvars = sourcedata->nvars;
   (*targetdata)->sizealldepsets = sourcedata->sizealldepsets;

   /* copy parameters */
   (*targetdata)->mode = sourcedata->mode;
   (*targetdata)->para_regterm = sourcedata->para_regterm;
   (*targetdata)->ubnumex = sourcedata->ubnumex;
   (*targetdata)->bmode = sourcedata->bmode;
   (*targetdata)->savememory = sourcedata->savememory;
   assert( (sourcedata->savememory == TRUE && sourcedata->savememory)
         || sourcedata->savememory == FALSE );

   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->y, n));
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->x, np));
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->xy, p));

   if( sourcedata->savememory == FALSE )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->xx, pp));
   else
      (*targetdata)->xx = NULL;


   if( ndep )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->indexdepsets, size_all));
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->sizedepsets, ndep));
   }

   /* copy */
   (*targetdata)->yy = sourcedata->yy;

   for( i = 0; i < n; i++ )
      (*targetdata)->y[i] = sourcedata->y[i];

   for( i = 0; i < np; i++ )
      (*targetdata)->x[i] = sourcedata->x[i];

   if( sourcedata->savememory == FALSE )
   {
      for( i = 0; i < pp; i++ )
         (*targetdata)->xx[i] = sourcedata->xx[i];
   }

   for( i = 0; i < p; i++ )
      (*targetdata)->xy[i] = sourcedata->xy[i];

   if( ndep )
   {
      for( i = 0; i < size_all; i++ )
         (*targetdata)->indexdepsets[i] = sourcedata->indexdepsets[i];

      for( i = 0; i < ndep; i++ )
         (*targetdata)->sizedepsets[i] = sourcedata->sizedepsets[i];
   }
   else
   {
      (*targetdata)->indexdepsets = NULL;
      (*targetdata)->sizedepsets = NULL;
   }

   /* allocate memory for variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->a, p));
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->z, p));
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->ep, n));

   /* copy variables */
   for( i = 0; i < p; i++ )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->a[i], &((*targetdata)->a[i]), varmap, consmap, global, &success));
      assert(success);
      SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->a[i]));

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->z[i], &((*targetdata)->z[i]), varmap, consmap, global, &success));
      assert(success);
      SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->z[i]));
   }

   for( i = 0; i < n; i++ )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->ep[i], &((*targetdata)->ep[i]), varmap, consmap, global, &success));
      assert(success);
      SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->ep[i]));
   }

   SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->rss, &((*targetdata)->rss), varmap, consmap, global, &success));
   assert(success);
   SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->rss));

   SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->log_rss, &((*targetdata)->log_rss), varmap, consmap, global, &success));
   assert(success);
   SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->log_rss));

   /* branching variables */
   for( i = 0; i < p; i++ )
      SCIP_CALL( SCIPchgVarBranchPriority(scip, (*targetdata)->z[i], 1000));

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigLinereg)
{
   SCIPdebugMessage("probdelorigLinereg \n");
   SCIPdebugMessage("free original problem data\n");

   /* free the (original) probdata */
   SCIP_CALL( probdataFree(scip, probdata));

   return SCIP_OKAY;
}


/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransLinereg)
{
   SCIP_Real timelimit;
   SCIP_Bool update;

   int n;
   int p;
   int ndep;
   int size_all;
   int np;
   int pp;

   int i;

   n = sourcedata->n;
   p = sourcedata->p;
   ndep = sourcedata->ndep;
   size_all = sourcedata->sizealldepsets;

   np = n * p;
   pp = p * p;

   SCIPdebugMessage("probtransLinereg \n");
   assert(scip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);
   assert(n > 0);
   assert(p > 0);
   assert(ndep >= 0);
   assert( (ndep == 0 && size_all == 0) || (ndep > 0 && size_all > 0) );

   SCIP_CALL( SCIPgetBoolParam(scip, "linereg/countpresoltime", &update));

   /* adjust time limit to take into account reading time */
   if( update )
   {
      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit));
      timelimit -= SCIPgetReadingTime(scip);
      timelimit = MAX(0.0,timelimit);
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timelimit));
   }

   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata));

   /* copy */
   (*targetdata)->n = sourcedata->n;
   (*targetdata)->p = sourcedata->p;
   (*targetdata)->ndep = sourcedata->ndep;
   (*targetdata)->nvars = sourcedata->nvars;
   (*targetdata)->sizealldepsets = sourcedata->sizealldepsets;

   /* copy parameters */
   (*targetdata)->mode = sourcedata->mode;
   (*targetdata)->para_regterm = sourcedata->para_regterm;
   (*targetdata)->ubnumex = sourcedata->ubnumex;
   (*targetdata)->bmode = sourcedata->bmode;
   (*targetdata)->savememory = sourcedata->savememory;
   assert( (sourcedata->savememory == TRUE && sourcedata->savememory)
         || sourcedata->savememory == FALSE );

   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->y, n));
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->x, np));
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->xy, p));
   if( sourcedata->savememory == FALSE )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->xx, pp));
   else
      (*targetdata)->xx = NULL;

   if( ndep )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->indexdepsets, size_all));
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->sizedepsets, ndep));
   }

   /* copy */
   (*targetdata)->yy = sourcedata->yy;

   for( i = 0; i < n; i++ )
      (*targetdata)->y[i] = sourcedata->y[i];

   for( i = 0; i < np; i++ )
      (*targetdata)->x[i] = sourcedata->x[i];

   if( sourcedata->savememory == FALSE )
   {
      for( i = 0; i < pp; i++ )
         (*targetdata)->xx[i] = sourcedata->xx[i];
   }

   for( i = 0; i < p; i++ )
     (*targetdata)->xy[i] = sourcedata->xy[i];

   if( ndep )
   {
      for( i = 0; i < size_all; i++ )
         (*targetdata)->indexdepsets[i] = sourcedata->indexdepsets[i];

      for( i = 0; i < ndep; i++ )
         (*targetdata)->sizedepsets[i] = sourcedata->sizedepsets[i];
   }
   else
   {
      (*targetdata)->indexdepsets = NULL;
      (*targetdata)->sizedepsets = NULL;
   }

   /* alloc memory for variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->a, p));
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->z, p));
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->ep, n));

   SCIP_CALL( SCIPtransformVars(scip, p, sourcedata->a, (*targetdata)->a));
   SCIP_CALL( SCIPtransformVars(scip, p, sourcedata->z, (*targetdata)->z));
   SCIP_CALL( SCIPtransformVars(scip, n, sourcedata->ep, (*targetdata)->ep));
   SCIP_CALL( SCIPtransformVars(scip, 1, &sourcedata->rss, &(*targetdata)->rss));
   SCIP_CALL( SCIPtransformVars(scip, 1, &sourcedata->log_rss, &(*targetdata)->log_rss));

   /* branching variables */
   for( i = 0; i < p; i++ )
      SCIP_CALL( SCIPchgVarBranchPriority(scip, (*targetdata)->z[i], 1000));

   return SCIP_OKAY;
}


static
SCIP_DECL_PROBEXITSOL(probexitsolLinreg)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(probdata != NULL);

   return SCIP_OKAY;
}


/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransLinereg)
{
   SCIPdebugMessage("free transformed problem data\n");
   SCIP_CALL( probdataFree(scip, probdata) );


   return SCIP_OKAY;
}


/*
 * Interface methods
 */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   )
{
   SCIP_PROBDATA* probdata;
   int n;
   int p;
   int i_ex;
   int ndep;
   SCIP_Real* data;                  /**< array to store data */
   SCIP_Real* data_explanatoryvars;  /**< array to store data of explanatory variables */
   int* ldindex;
   int* maxdep;
   int* depsets;
   int i;
   int j;
   int ct;
   char probname[SCIP_MAXSTRLEN];
   int np;
   int pp;

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata));

   /* read dimension of data points */
   SCIP_CALL( readDataDim(filename, &n, &p, &i_ex));

   probdata->n = n;
   probdata->p = p;

   np = n * p;
   pp = p * p;

   if( probdata->mode == 'a' )
      probdata->para_regterm = 2.0;
   else if( probdata->mode == 'b' )
      probdata->para_regterm = log((SCIP_Real) n);
   else
   {
      SCIPerrorMessage("mode is %c\n", probdata->mode);
      assert(0);
      return SCIP_ERROR;
   }

   /* allocate memory for data */
   SCIP_CALL( SCIPallocMemoryArray(scip, &data, np + n));

   /* read data points */
   SCIP_CALL( readData(filename, n, p, data));

   /* check data */
   SCIP_CALL( checkData(n, p, i_ex, data));

   /* normalize data */
   SCIP_CALL( normalization(scip, n, p, data));

   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->y, n));
   SCIP_CALL( SCIPallocMemoryArray(scip, &data_explanatoryvars, np));

   /* divide data into explained variable and explanatory variables */
   SCIP_CALL( divideData(n, p, i_ex, data, probdata->y, data_explanatoryvars));

   /* free memory for data */
   SCIPfreeMemoryArrayNull(scip, &data);

   /* get a probdata's name */
   SCIP_CALL( getProblemName(filename, probname, SCIP_MAXSTRLEN));

   /* output information of the problem */
   SCIPinfoMessage(scip, NULL, "File name\t:\t%s\n", filename);
   SCIPinfoMessage(scip, NULL, "Problem name \t:\t%s\n", probname);
   SCIPinfoMessage(scip, NULL, "Number of data\t:\t%d\n", n);
   SCIPinfoMessage(scip, NULL, "Number of var\t:\t%d\n", p);

   /* define arrays with ColMajor for cblas and clpack
    * x[n*p] := data_explanatoryvars[n][p]
    * xx[p*p] := data_explanatoryvars' data_explanatoryvars
    * xy[p] := data_explanatoryvars' y
    * yy := y'y
    */

   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->x, np));
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->xy, p));

   /* calculate values */
   SCIP_CALL( SCIPtransMat(n, p, data_explanatoryvars, probdata->x));
   SCIP_CALL( SCIPcblasDgemv3(probdata->x, n, p, probdata->y, probdata->xy));
   probdata->yy = SCIPcblasDdot(probdata->y, probdata->y, n);

   if( probdata->savememory == FALSE )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->xx, pp));
      SCIP_CALL( SCIPcblasDgemm1(probdata->x, n, p, probdata->xx));
   }
   else
      probdata->xx = NULL;

   /* free for data_explanatoryvars */
   SCIPfreeMemoryArrayNull(scip, &data_explanatoryvars);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &ldindex, p));

   /* calculate the number of linearly dependent sets */
   SCIP_CALL( SCIPgetNLineDependSet(scip, probdata->x, n, p, ldindex));
   ndep = SCIPcalcIntSum(ldindex, p);
   probdata->ndep = ndep;

   /******************/
   if( ndep )
   {
      int ndep_p = ndep * p;
      int buf1;
      int buf2;
      int size;
      int sizealldepsets;

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &maxdep, ndep));
      SCIP_CALL( SCIPallocBufferArray(scip, &depsets, ndep_p));

      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->sizedepsets, ndep));

      SCIP_CALL( SCIPinitIntArrayZero(ndep_p, depsets));

      ct = 0;

      for(i = 0; i < p; i++)
      {
         if( ldindex[i] == 1 )
            maxdep[ct++] = i;
      }


      for( i = 0; i < ndep; i++ )
         *(depsets + (i * p) + maxdep[i]) = 1;

      /* find linearly dependent sets */
      if( probdata->savememory == FALSE )
      {
         assert(probdata->xx != NULL);
         SCIP_CALL( SCIPgetLineDependSet(scip, probdata->xx, p, ndep, maxdep, ldindex, depsets));
      }
      else
      {
         SCIP_Real* xx;
         assert(probdata->xx == NULL);
         SCIP_CALL( SCIPallocBufferArray(scip, &xx, pp));
         SCIP_CALL( SCIPcblasDgemm1(probdata->x, n, p, xx));
         SCIP_CALL( SCIPgetLineDependSet(scip, xx, p, ndep, maxdep, ldindex, depsets));
         SCIPfreeBufferArray(scip, &xx);
      }

      /* calculate sizedepsets */
      sizealldepsets = 0;
      for( i = 0; i < ndep; i++ )
      {
         probdata->sizedepsets[i] = SCIPcalcIntSum(&depsets[i * p], p);
         sizealldepsets += probdata->sizedepsets[i];
      }

      assert(sizealldepsets >= 2 * ndep);

      /* allocate memory */
      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->indexdepsets, sizealldepsets));

      /* set indexdepsets */
      ct = 0;
      for( i = 0; i < ndep; i++ )
      {
         buf1 = i * p;
         buf2 = 0;
         size = probdata->sizedepsets[i];
         for( j = 0; j < p; j++ )
         {
            if( depsets[buf1 + j] == 1 )
            {
               probdata->indexdepsets[ct++] = j;
               buf2++;

               if( buf2 == size )
                  break;
            }
         }
         assert(buf2 == size);
      }
      assert(ct == sizealldepsets);

      /* set sizealldepsets */
      probdata->sizealldepsets = sizealldepsets;
   }
   else
   {
      maxdep = NULL;
      depsets = NULL;

      probdata->indexdepsets = NULL;
      probdata->sizealldepsets = 0;
      probdata->sizedepsets = NULL;
   }

   /* print linearly dependent sets */
   SCIP_CALL( SCIPprintLineDependSet(scip, ndep, p, depsets));

   /*
    * if( ndep )
    *     free for depsets and maxdep
    */
   if( ndep )
   {
      /* free */
      SCIPfreeBufferArray(scip, &maxdep);
      SCIPfreeBufferArray(scip, &depsets);
   }

   /* free */
   SCIPfreeMemoryArrayNull(scip, &ldindex);

   /* create a problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname));
   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigLinereg));
   SCIP_CALL( SCIPsetProbTrans(scip, probtransLinereg));
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransLinereg));
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolLinreg));
   SCIP_CALL( SCIPsetProbCopy(scip, probcopyLinereg));

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata));

   /* disable sub-SCIP heuristics */
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rens/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rins/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/dins/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/crossover/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/mutation/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/clique/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/vbounds/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/bound/freq", -1));
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/zeroobj/freq", -1));

   /* create and add variables */
   SCIP_CALL( createVariables(scip, probdata));

   /* create and add constraints */
   SCIP_CALL( createConstraints(scip, probdata));

   return SCIP_OKAY;
}


/** return the number of data points */
int SCIPprobdataGetNdatas(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->n;
}


/** return the number of explanatory variables */
int SCIPprobdataGetNexvars(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->p;
}


/** return the number of variables in thip optimization */
int SCIPprobdataGetNvars(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->nvars;
}


/** return the number of linearly dependent sets */
int SCIPprobdataGetNdep(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->ndep;
}


/** return pointer to y */
SCIP_Real* SCIPprobdataGety(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->y;
}


/** return pointer to x */
SCIP_Real* SCIPprobdataGetx(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->x;
}


/** return pointer to xx */
SCIP_Real* SCIPprobdataGetxx(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->xx;
}


/** return pointer to xy */
SCIP_Real* SCIPprobdataGetxy(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->xy;
}


/** return yy */
SCIP_Real SCIPprobdataGetyy(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->yy;
}


/** return pointer to indexdepsets */
int* SCIPprobdataGetindexdepsets(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->indexdepsets;
}


/** return sizealldepsets */
int SCIPprobdataGetsizealldepsets(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->sizealldepsets;
}


/** return pointer to sizedepsets */
int* SCIPprobdataGetsizedepsets(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->sizedepsets;
}


/** return pointer to variable a */
SCIP_VAR** SCIPprobdataGetVars_a(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->a;
}


/** return pointer to variable z */
SCIP_VAR** SCIPprobdataGetVars_z(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->z;
}


/** return pointer to vairable ep */
SCIP_VAR** SCIPprobdataGetVars_ep(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->ep;
}


/** return variable rss */
SCIP_VAR* SCIPprobdataGetVar_rss(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->rss;
}


/** return variable log */
SCIP_VAR* SCIPprobdataGetVar_log(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->log_rss;
}


/** return solving mode */
char SCIPprobdataGetMode(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->mode;
}


/** return parameter to the regularization term */
SCIP_Real SCIPprobdataGetPara_RegTerm(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->para_regterm;
}


/** return parameter to the number of explanatory variables */
int SCIPprobdataGetPara_NumEx(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->ubnumex;
}


/** return branching mode */
char SCIPprobdataGetBranchMode(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->bmode;
}


/** return savememory */
SCIP_Bool SCIPprobdataGetSaveMemory(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   )
{
   assert(probdata != NULL);
   return probdata->savememory;
}


/** set dual bound by ug */
void SCIPprobdataSetDualBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             dual                /**< dual bound */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);

   probdata = SCIPgetProbData(scip);
   probdata->ug = TRUE;
   probdata->ugDual = dual;
}


/** set the number of solvers */
void SCIPprobdataSetNSolvers(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nSolvers            /**< the number of solvers */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   probdata->nSolvers = nSolvers;
}
