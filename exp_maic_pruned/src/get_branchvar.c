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

/**@file   get_branchvar.c
 * @brief  function to return the last branching variable
 * @author Keiji Kimura
 *
 * This file implements the function to return the last branching variable
 *
 * If the last branching (binary) variable z_j = 1, the relaxator can compute the lower bound easily.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "convenient_tool.h"
#include "get_branchvar.h"


/** retrun the last branching variable */
SCIP_RETCODE SCIPvarGetLastBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< the current node */
   SCIP_VAR**            result              /**< pointer to store the last branching variable */
   )
{
   SCIP_VAR** branchvars;                    /* array of variables on which the branchings has been performed in all ancestors */
   SCIP_Real* branchbounds;                  /* array of bounds which the branchings in all ancestors set */
   SCIP_BOUNDTYPE* boundtypes;               /* array of boundtypes which the branchings in all ancestors set */
   int* nodeswitches;                        /* marks, where in the arrays the branching decisions of the next node on the path start
                                              * branchings performed at the parent of node always start at position 0. For single variable branching,
                                              * nodeswitches[i] = i holds */
   int nbranchvars;                          /* number of variables on which branchings have been performed in all ancestors
                                              *   if this is larger than the array size, arrays should be reallocated and method should be called again */
   int branchvarssize;                       /* available slots in arrays */
   int nnodes;                               /* number of nodes in the nodeswitch array */
   int nodeswitchsize;                       /* available slots in node switch array */

   SCIP_VAR* lastbranchvar;

   branchvarssize = SCIPnodeGetDepth(node);
   nodeswitchsize = branchvarssize;

   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchvars, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &branchbounds, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodeswitches, nodeswitchsize) );

   SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize, nodeswitches, &nnodes, nodeswitchsize );

   /* if the arrays were to small, we have to reallocate them and recall SCIPnodeGetAncestorBranchingPath */
   if( nbranchvars > branchvarssize || nnodes > nodeswitchsize )
   {
      branchvarssize = nbranchvars;
      nodeswitchsize = nnodes;

      /* memory reallocation */
      SCIP_CALL( SCIPreallocBufferArray(scip, &branchvars, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &branchbounds, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &boundtypes, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &nodeswitches, nodeswitchsize) );

      SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize, nodeswitches, &nnodes, nodeswitchsize);
      assert(nbranchvars == branchvarssize);
   }

   /* get the last branching */
   if( nbranchvars >= 1 )
      lastbranchvar = branchvars[nodeswitches[0]];
   else
      lastbranchvar = NULL;

#if 0
   if( lastbranchvar != NULL )
   {
      if( SCIPvarGetType(lastbranchvar) != SCIP_VARTYPE_BINARY )
      {
         SCIP_CALL( SCIPprintNodeRootPath(scip, node, NULL));
         printf("%d\n", SCIPvarGetType(lastbranchvar));
         SCIPexit();
      }
   }
#endif

   assert( lastbranchvar == NULL || SCIPvarGetType(lastbranchvar) == SCIP_VARTYPE_BINARY );


   /* free all local memory */
   SCIPfreeBufferArray(scip, &nodeswitches);
   SCIPfreeBufferArray(scip, &boundtypes);
   SCIPfreeBufferArray(scip, &branchbounds);
   SCIPfreeBufferArray(scip, &branchvars);

   *result = lastbranchvar;
   return SCIP_OKAY;
}
