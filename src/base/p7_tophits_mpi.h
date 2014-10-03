#ifndef p7TOPHITS_MPI_INCLUDED
#define p7TOPHITS_MPI_INCLUDED
#include "p7_config.h"

#ifdef HAVE_MPI
#include <mpi.h>

#include "base/p7_tophits.h"

#define p7TOPHITS_MPI_HITBLOCKSIZE 100 /* when sending P7_HIT arrays, how many to send in one MPI message */

extern int p7_tophits_mpi_Send    (const P7_TOPHITS *th, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int p7_tophits_mpi_PackSize(const P7_TOPHITS *th, MPI_Comm comm, int *ret_n);
extern int p7_tophits_mpi_Pack    (const P7_TOPHITS *th, char *buf, int n, int *pos, MPI_Comm comm);
extern int p7_tophits_mpi_Unpack  (char *buf, int n, int *pos, MPI_Comm comm, P7_TOPHITS **ret_th, int *ret_nhits);
extern int p7_tophits_mpi_Recv    (int source, int tag, MPI_Comm comm, char **buf, int *nalloc, P7_TOPHITS **ret_th);

extern int p7_hit_mpi_Send    (const P7_HIT *hit, int nhit, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int p7_hit_mpi_PackSize(const P7_HIT *hit, int nhit, MPI_Comm comm, int *ret_n);
extern int p7_hit_mpi_Pack    (const P7_HIT *hit, int nhit, char *buf, int n, int *pos, MPI_Comm comm);
extern int p7_hit_mpi_Unpack  (char *buf, int n, int *pos, MPI_Comm comm, P7_HIT *hit, int nhit);
extern int p7_hit_mpi_Recv    (int source, int tag, MPI_Comm comm, char **buf, int *nalloc, P7_HIT *hit, int nhit);

#endif /*HAVE_MPI*/
#endif /*p7TOPHITS_MPI_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
