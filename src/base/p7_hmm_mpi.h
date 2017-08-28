#ifndef p7HMM_MPI_INCLUDED
#define p7HMM_MPI_INCLUDED

#include "p7_config.h"

#ifdef  HAVE_MPI
#include <mpi.h>

#include "esl_alphabet.h"

#include "base/p7_hmm.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
extern int p7_hmm_mpi_Send    (const P7_HMM *hmm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int p7_hmm_mpi_PackSize(const P7_HMM *hmm, MPI_Comm comm, int *ret_n);
extern int p7_hmm_mpi_Pack    (const P7_HMM *hmm, char *buf, int n, int *pos, MPI_Comm comm);
extern int p7_hmm_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm,                 ESL_ALPHABET **byp_abc, P7_HMM **ret_hmm);
extern int p7_hmm_mpi_Recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **byp_abc, P7_HMM **ret_hmm);
#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*HAVE_MPI*/
#endif /*p7HMM_MPI_INCLUDED*/

