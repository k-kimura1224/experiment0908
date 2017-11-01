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

/**@file   reader_linereg.c
 * @brief  file reader for linereg regression data
 * @author Keiji Kimura
 *
 * This file implements the reader for linear regression data.
 *
 * This reader reads from the second line of a file.
 * The format of the file is as follows:
 * - [2nd line] the number of data
 * - [3rd line] the number of the explanatory variables
 * - [4th line] the index of the response variable(Begining is 1)
 * - [from 5th line] data points
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>
#include <string.h>

#include "probdata_linereg.h"
#include "reader_linereg.h"
#include "set_myparameter.h"


#define READER_NAME           "lineregreader"
#define READER_DESC           "file reader for linear regression data format"
#define READER_EXTENSION      "linereg"

#define DEFAULT_COUNTPRESOLTIME  TRUE      /**< count presolving time as part of overall solution time? */


/**@name Callback methods
 *
 * @{
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyLinereg)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderLinereg(scip) );

   return SCIP_OKAY;
}

/** problem reading method of the reader */
static
SCIP_DECL_READERREAD(readerReadLinereg)
{  /*lint --e{715}*/
   SCIP_RETCODE          retcode;
   SCIP_PROBDATA*        probdata;

   *result = SCIP_DIDNOTRUN;

   retcode = SCIPprobdataCreate(scip, filename);

   if( retcode == SCIP_READERROR )
      return SCIP_READERROR;

   SCIP_CALL( retcode );

   probdata = SCIPgetProbData(scip);
   if( SCIPgetStage(scip) == SCIP_STAGE_INIT ||  probdata == NULL )
      return SCIP_READERROR;

#if MP_PRINTORIGP
      SCIPinfoMessage(scip, NULL, "Original problem:\n");
      SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

#if MP_PRINTREFOP
      SCIP_CALL( SCIPpresolve(scip) );
      SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
#endif

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** problem writing method of the reader */
static
SCIP_DECL_READERWRITE(readerWriteLinereg)
{  /*lint --e{715}*/

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** includes the stp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderLinereg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   readerdata = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyLinereg) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadLinereg) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteLinereg) );

   /* include user parameters */

   SCIP_CALL( SCIPaddBoolParam(scip,
         "linereg/countpresoltime",
         "count presolving time to solving time?",
         NULL, FALSE, DEFAULT_COUNTPRESOLTIME, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip,
         "linereg/mode",
         "Solving mode: 'a'ic, 'b'ic",
         NULL, FALSE, 'a', "ab", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "linereg/ubnumex",
         "Upper bound of the number of explanatory variable. -1 means p.",
         NULL, FALSE, -1, -1, 100, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip,
         "linereg/bmode",
         "Branching mode: 'a'utomatic, 's'trong 'm'ost frequent",
         NULL, FALSE, 'a', "asm", NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
            "linereg/savememory",
            "save memory?",
            NULL, FALSE, FALSE, NULL, NULL));

   return SCIP_OKAY;
}

/**@} */
