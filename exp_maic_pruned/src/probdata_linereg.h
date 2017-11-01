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

/**@file   probdata_linereg.h
 * @brief  problem data for minimization of AIC/BIC
 * @author Keiji Kimura
 *
 * This file implements the problem data for minimization of AIC/BIC
 *
 * The problem data contains the number of data, the number of explanatory variables, and data points.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_PROBDATA_LINEREG__
#define __SCIP_PROBDATA_LINEREG__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** sets up the problem data */
extern
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   );


/** return the number of data points */
extern
int SCIPprobdataGetNdatas(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return the number of explanatory variables */
extern
int SCIPprobdataGetNexvars(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return the number of variables in thip optimization */
extern
int SCIPprobdataGetNvars(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return the number of linearly dependent sets */
extern
int SCIPprobdataGetNdep(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return pointer to y */
extern
SCIP_Real* SCIPprobdataGety(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return pointer to x */
extern
SCIP_Real* SCIPprobdataGetx(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return pointer to xx */
extern
SCIP_Real* SCIPprobdataGetxx(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return pointer to xy */
extern
SCIP_Real* SCIPprobdataGetxy(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return yy */
extern
SCIP_Real SCIPprobdataGetyy(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );



/** return pointer to indexdepsets */
extern
int* SCIPprobdataGetindexdepsets(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return sizealldepsets */
extern
int SCIPprobdataGetsizealldepsets(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return pointer to sizedepsets */
extern
int* SCIPprobdataGetsizedepsets(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return pointer to variable a */
extern
SCIP_VAR** SCIPprobdataGetVars_a(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return pointer to variable z */
extern
SCIP_VAR** SCIPprobdataGetVars_z(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return pointer to vairable ep */
extern
SCIP_VAR** SCIPprobdataGetVars_ep(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return variable rss */
extern
SCIP_VAR* SCIPprobdataGetVar_rss(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return variable log */
extern
SCIP_VAR* SCIPprobdataGetVar_log(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return solving mode */
extern
char SCIPprobdataGetMode(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return parameter to the regularization term */
extern
SCIP_Real SCIPprobdataGetPara_RegTerm(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return parameter to the number of explanatory variables */
extern
int SCIPprobdataGetPara_NumEx(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return branching mode */
extern
char SCIPprobdataGetBranchMode(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** return savememory */
extern
SCIP_Bool SCIPprobdataGetSaveMemory(
   SCIP_PROBDATA*        probdata            /**< user problem data */
   );


/** set dual bound by ug */
void SCIPprobdataSetDualBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             dual                /**< dual bound */
   );


/** set the number of solvers */
void SCIPprobdataSetNSolvers(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nSolvers            /**< the number of solvers */
   );


#ifdef __cplusplus
}
#endif

#endif
