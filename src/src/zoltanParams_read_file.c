/* 
   zoltanParams_read_file.c

   Read Zoltan parameters from a file, call Zoltan to set the parameters

   zoltanParams library

   Jim Teresco

   Department of Computer Science
   Williams College

   and

   Computer Science Research Institute
   Sandia National Laboratories

   $Id: zoltanParams_read_file.c,v 1.8 2006/06/19 20:00:13 terescoj Exp $

   Modification History:
   2004/01/30  JDT  Adapted from PMDB's Zoltan interface implementation
   2004/02/23  JDT  Preliminary fortran interface
   2004/03/10  JDT  Removed obsolete num_partitions callback things
   2004/03/22  JDT  Added ability to check if graph info is needed
   2006/06/16  JDT  New Fortran-callable version
   2006/06/19  JDT  Fixes and cleanup for new Fortran interface
*/

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include "zoltan.h"

#define SAFE_MALLOC(v,type,size) \
 {  v = (type) malloc(size) ; \
    if ( v == NULL) { \
       fflush(stdout); \
       fprintf(stderr,"in file %s, line %d, failed to allocate %ld bytes",\
               __FILE__,__LINE__,size); \
       MPI_Abort(MPI_COMM_WORLD,1); \
    } \
 }

struct zoltanParams_list_entry {
  char *param;
  char *value;
  struct zoltanParams_list_entry *next;
};

struct zoltanParams_hier_struct {
  int partition;
  struct zoltanParams_list_entry *first;
};

static struct zoltanParams_hier_struct **zph = NULL;
static int num_levels = 0;
static MPI_Comm comm;
static int using_graph = 0;

static void check_level(int level) {

  if (!zph) {
    fprintf(stderr,"check_level: must set number of levels first\n");
    return;
  }

  if (level >= num_levels) {
    fprintf(stderr,"check_level: invalid level\n");
  }
}

void zoltanParams_hier_free() {
  int i;

  if (!zph) {
    fprintf(stderr, "zoltanParams_hier_free warning: not allocated\n");
    return;
  }

  for (i=0; i<num_levels; i++) {
    free(zph[i]);
  }

  free(zph);
}

void zoltanParams_hier_set_num_levels(int levels) {
  int i;

#ifdef DEBUG
  printf("(zoltanParams_hier_set_num_levels) setting to %d\n", levels);  
  fflush(stdout);
#endif

  if (zph) {
    fprintf(stderr,"zoltanParams_hier_set_num_levels warning: already initialized, reinitializing\n");
    zoltanParams_hier_free();
  }

  if (levels <= 0) {
    fprintf(stderr, "(zoltanParams_hier_set_num_levels) num levels must be positive\n");
    return;
  }

  num_levels = levels;

  SAFE_MALLOC(zph, struct zoltanParams_hier_struct **, 
	      sizeof(struct zoltanParams_hier_struct *) * levels);

  for (i=0; i<levels; i++) {
    SAFE_MALLOC(zph[i],  struct zoltanParams_hier_struct *, 
		sizeof (struct zoltanParams_hier_struct));
    zph[i]->partition = 0;
    zph[i]->first = NULL;
  }
}

void zoltanParams_hier_set_partition(int level, int partition) {

#ifdef DEBUG
  int mypid;
  MPI_Comm_rank(comm, &mypid);

  printf("[%d] will compute partition %d at level %d\n", 
	 mypid, partition, level); fflush(stdout);
#endif

  check_level(level);

  zph[level]->partition = partition;
}

void zoltanParams_hier_set_param(int level, char *param, char *value) {
  struct zoltanParams_list_entry *newparam, *nextparam;

#ifdef DEBUG
  int mypid;
  MPI_Comm_rank(comm, &mypid);
  printf("[%d] will set param <%s> to <%s> at level %d\n", 
	 mypid, param, value, level); fflush(stdout);
#endif

  check_level(level);

  SAFE_MALLOC(newparam, struct zoltanParams_list_entry *,
	      sizeof(struct zoltanParams_list_entry));

  newparam->param = strdup(param);
  newparam->value = strdup(value);
  newparam->next = NULL;

  if ((strcasecmp(param, "LB_METHOD") == 0) &&
      ((strcasecmp(value, "PARMETIS") == 0)
       || (strcasecmp(value, "JOSTLE") == 0)))
    using_graph = 1;      
  
  if (!zph[level]->first) {
    zph[level]->first = newparam;
    return;
  }

  nextparam = zph[level]->first;
  while (nextparam->next) nextparam=nextparam->next;
  nextparam->next = newparam;    
}

int zoltanParams_hier_get_num_levels() {

  return num_levels;
}

int zoltanParams_hier_get_partition(int level) {

  check_level(level);

  return zph[level]->partition;
}

void zoltanParams_hier_use_params(int level, struct Zoltan_Struct *zz, int *ierr) {
  struct zoltanParams_list_entry *nextparam;

  *ierr = ZOLTAN_OK;
  check_level(level);
  
  nextparam = zph[level]->first;

  while (nextparam) {
    *ierr = Zoltan_Set_Param(zz, nextparam->param, nextparam->value);
    if (*ierr != ZOLTAN_OK) return;
    nextparam = nextparam->next;
  }
  
}

static int get_num_levels(void *data, int *ierr) {

  *ierr = ZOLTAN_OK;
  return zoltanParams_hier_get_num_levels();
}

static int get_partition(void *data, int level, int *ierr) {

  *ierr = ZOLTAN_OK;

  return zoltanParams_hier_get_partition(level);
}

static void get_method(void *data, int level, struct Zoltan_Struct *zz,
		       int *ierr) {

  zoltanParams_hier_use_params(level, zz, ierr);
}

void zoltanParams_set_comm(MPI_Comm thecomm) {

  /* remember the comm passed in */
  MPI_Comm_dup(thecomm, &comm);
}

void zoltanParams_hier_setup(struct Zoltan_Struct *zz) {

  /* make sure the hierarchical balancing callbacks are in place */
  if (Zoltan_Set_Fn(zz, ZOLTAN_HIER_NUM_LEVELS_FN_TYPE, 
		    (void (*)()) get_num_levels, NULL) == ZOLTAN_FATAL) {
    fprintf(stderr,"zoltanParams_hier_setup: set NUM_LEVELS callback failed\n");
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_HIER_PARTITION_FN_TYPE, 
		    (void (*)()) get_partition, NULL) == ZOLTAN_FATAL) {
    fprintf(stderr,"zoltanParams_hier_setup: set PARTITION callback failed\n");
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_HIER_METHOD_FN_TYPE, 
		    (void (*)()) get_method, NULL) == ZOLTAN_FATAL) {
    fprintf(stderr,"zoltanParams_hier_setup: set METHOD callback failed\n");
  }   
}

/*

  zoltanParams_read_file

  Set up the given Zoltan_Struct with parameters as specified
  in the given file.

  File format:

  Lines of the format:
  ZOLTAN_PARAM PARAM_VALUE

  If the parameter is LB_METHOD set to HIER, the next part of the file
  is interpreted as hierarchical balancing parameters:

  num_levels
  level 0 partitions for each proc
  level 0 parameters
  end with LEVEL END
  level 1 partitions for each proc
  level 1 parameters
  end with LEVEL END
  ...

  End file with EOF

*/
void zoltanParams_read_file(struct Zoltan_Struct *lb, char *file, 
			    MPI_Comm thecomm) {
  FILE *fp;
  char str1[500], str2[500], str3[1000];
  int numlevels, level, partition, proc;
  int ierr;
  int mypid, numprocs;

  /* remember the comm passed in */
  MPI_Comm_dup(thecomm, &comm);

  MPI_Comm_rank(comm, &mypid);
  MPI_Comm_size(comm, &numprocs);

  fp = fopen(file, "r");
  if (!fp) {
    fprintf(stderr,"Cannot open file %s for reading", file);
    return;
  } 

#ifdef DEBUG
  if (mypid == 0) {
    printf("Reading Zoltan parameters from file %s\n", file);
  }
#endif

  while (fscanf(fp, "%s %s\n", str1, str2) == 2) {
    ierr = Zoltan_Set_Param(lb, str1, str2);
    if (ierr != ZOLTAN_OK) {
      fprintf(stderr,"Zoltan_Set_Param failed to set param <%s> to <%s>",str1,str2);
    }
#ifdef DEBUG
    else {
      if (mypid == 0) {
	printf("Set Zoltan parameter <%s> to <%s>\n", str1, str2);
      }
    }
#endif
    if ((strcasecmp(str1, "LB_METHOD") == 0) &&
	((strcasecmp(str2, "PARMETIS") == 0) || 
	 (strcasecmp(str2, "JOSTLE") == 0)))
      using_graph = 1;
    
    if (strcasecmp(str1,"LB_METHOD") == 0 && strcasecmp(str2,"HIER") == 0) {

      zoltanParams_hier_setup(lb);

      /* the rest of the file contains hierarchical balancing parameters */
      fscanf(fp, "%d", &numlevels);
#ifdef DEBUG
      printf("[%d] read in numlevels=%d\n", mypid, numlevels);
#endif
      zoltanParams_hier_set_num_levels(numlevels);

      for (level=0; level<numlevels; level++) {
	/* first, a list of partitions for each proc should be in the file */
	for (proc=0; proc<numprocs; proc++) {
	  fscanf(fp, "%d", &partition);
	  if (proc == mypid) zoltanParams_hier_set_partition(level, partition);
	}
	/* then parameters until we get LEVEL END */
	while ((fscanf(fp, "%s %s\n", str1, str2) == 2) &&
	       (strcmp(str1, "LEVEL") != 0) &&
	       (strcmp(str2, "END") != 0)) {
	  
	  zoltanParams_hier_set_param(level, str1, str2);
	}
      }
    }
  }
  fclose(fp);
}

#ifdef FORTRAN_INTERFACE
/* Fortran callable interface -- tested only for Lahey compiler/Linux */
void zf90_zoltanparams_read_file_(int* addr_lb, int* nbyte, int* istring,
                                  int* nchar, int *f_communicator, int* ierr) {
  struct Zoltan_Struct *lb;
  extern MPI_Comm Zoltan_comm_f2c(int *);
  unsigned char *p;
  int i;
  MPI_Comm c_communicator;
  unsigned char* cstring;

  /* convert the pointer, byte-by-byte */
  p = (unsigned char *) &lb;
  for (i=0; i<(*nbyte); i++) {*p = (unsigned char)addr_lb[i]; p++;}

  /* convert the string, byte-by-byte */
  cstring = malloc((*nchar)+1);
  p = (unsigned char *) cstring;
  for (i=0; i<(*nchar); i++) {*p = (unsigned char)istring[i]; p++;}
  *p = '\0';

  /* convert to C communicator -- from Zoltan's comm library */
  c_communicator = Zoltan_comm_f2c(f_communicator);

  zoltanParams_read_file(lb, cstring, c_communicator);

  free(cstring);
}
/* other forms of name mangling */
void zf90_zoltanparams_read_file(int* addr_lb, int* nbyte, int* istring,
                                  int* nchar, int *f_communicator, int* ierr) {
  zf90_zoltanparams_read_file_(addr_lb,nbyte,istring,nchar,f_communicator,ierr);
}
void zf90_zoltanparams_read_file__(int* addr_lb, int* nbyte, int* istring,
                                  int* nchar, int *f_communicator, int* ierr) {
  zf90_zoltanparams_read_file_(addr_lb,nbyte,istring,nchar,f_communicator,ierr);
}
void ZF90_ZOLTANPARAMS_READ_FILE(int* addr_lb, int* nbyte, int* istring,
                                  int* nchar, int *f_communicator, int* ierr) {
  zf90_zoltanparams_read_file_(addr_lb,nbyte,istring,nchar,f_communicator,ierr);
}
#endif /* FORTRAN_INTERFACE */

/* return whether graph callbacks/other info is needed by parameters that
   have been set */
int zoltanParams_using_graph() {

  return using_graph;
}

#ifdef FORTRAN_INTERFACE
/* fortran-callable */
void zoltanparams_using_graph_(int *retval) {

  *retval = using_graph;
}
/* other forms of name mangling */
void zoltanparams_using_graph(int *retval) {
  *retval = using_graph;
}
void zoltanparams_using_graph__(int *retval) {
  *retval = using_graph;
}
void ZOLTANPARAMS_USING_GRAPH(int *retval) {
  *retval = using_graph;
}
#endif /* FORTRAN_INTERFACE */
