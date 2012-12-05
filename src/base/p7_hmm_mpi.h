#ifndef P7_HMM_MPI_INCLUDED
#define P7_HMM_MPI_INCLUDED
#ifdef  HAVE_MPI

#include "p7_config.h"

#include <mpi.h>

#include "esl_alphabet.h"

#include "base/p7_hmm.h"

extern int p7_hmm_mpi_Send(P7_HMM *hmm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int p7_hmm_mpi_PackSize(P7_HMM *hmm, MPI_Comm comm, int *ret_n);
extern int p7_hmm_mpi_Pack(P7_HMM *hmm, char *buf, int n, int *pos, MPI_Comm comm);
extern int p7_hmm_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **abc, P7_HMM **ret_hmm);
extern int p7_hmm_mpi_Recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, P7_HMM **ret_hmm);

#endif /*HAVE_MPI*/
#endif /*P7_HMM_MPI_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
