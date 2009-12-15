/* 
   zoltanParams.h

   prototypes for Zoltan parameters from a file, call Zoltan to set
   the parameters

   zoltanParams library

   Jim Teresco

   Department of Computer Science
   Williams College

   and

   Computer Science Research Institute
   Sandia National Laboratories

   $Id: zoltanParams.h,v 1.4 2004/04/13 02:39:49 terescoj Exp $

   Modification History:
   2004/02/02 JDT Created
   2004/03/10 JDT Removed obsolete hier_num_partitions function
   2004/03/22 JDT Added zoltanParams_using_graph
   2004/04/12 JDT Added C++ extern "C" guard
*/

#ifndef __ZOLTANPARAMS_H
#define __ZOLTANPARAMS_H

#include <mpi.h>
#include "zoltan.h"

#ifdef __cplusplus
extern "C" {
#endif

void zoltanParams_hier_free();
void zoltanParams_hier_set_num_levels(int levels);
void zoltanParams_hier_set_partition(int level, int partition);
void zoltanParams_hier_set_param(int level, char *param, char *value);
int zoltanParams_hier_get_num_levels();
int zoltanParams_hier_get_partition(int level);
void zoltanParams_hier_use_params(int level, struct Zoltan_Struct *zz, 
				  int *ierr);
void zoltanParams_set_comm(MPI_Comm thecomm);
void zoltanParams_hier_setup(struct Zoltan_Struct *zz);
void zoltanParams_read_file(struct Zoltan_Struct *lb, char *file, 
			    MPI_Comm thecomm);
int zoltanParams_using_graph();

#ifdef __cplusplus
}
#endif

#endif
