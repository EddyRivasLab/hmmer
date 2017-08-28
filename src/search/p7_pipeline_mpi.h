#ifndef p7_PIPELINE_MPI_INCLUDED
#define p7_PIPELINE_MPI_INCLUDED
#ifdef  HAVE_MPI

#include "p7_config.h"

#include <mpi.h>

#include "search/p7_pipeline.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
extern int p7_pipeline_stats_mpi_Send    (const P7_PIPELINE_STATS *stats, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int p7_pipeline_stats_mpi_PackSize(const P7_PIPELINE_STATS *stats, MPI_Comm comm, int *ret_n);
extern int p7_pipeline_stats_mpi_Pack    (const P7_PIPELINE_STATS *stats, char *buf, int n, int *pos, MPI_Comm comm);
extern int p7_pipeline_stats_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, P7_PIPELINE_STATS *stats);
extern int p7_pipeline_stats_mpi_Recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, P7_PIPELINE_STATS *stats);
#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*HAVE_MPI*/
#endif /*p7PIPELINE_MPI_INCLUDED*/

