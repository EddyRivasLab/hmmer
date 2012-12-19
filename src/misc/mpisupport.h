#ifdef  HAVE_MPI

#ifndef P7_MPISUPPORT_INCLUDED
#define P7_MPISUPPORT_INCLUDED

#include "p7_config.h"

#include <mpi.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "base/p7_tophits.h"

extern int p7_tophits_MPISend(P7_TOPHITS *th, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int p7_tophits_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, P7_TOPHITS **ret_th);

#endif /*P7_MPISUPPORT_INCLUDED*/
#endif /*HAVE_MPI*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
