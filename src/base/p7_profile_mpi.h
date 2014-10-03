#ifndef p7PROFILE_MPI_INCLUDED
#define p7PROFILE_MPI_INCLUDED
#ifdef  HAVE_MPI

#include "p7_config.h"

#include <mpi.h>

#include "esl_alphabet.h"

#include "base/p7_bg.h"
#include "base/p7_profile.h"

extern int p7_profile_mpi_Send    (const P7_PROFILE *gm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int p7_profile_mpi_PackSize(const P7_PROFILE *gm, MPI_Comm comm, int *ret_n);
extern int p7_profile_mpi_Pack    (const P7_PROFILE *gm, char *buf, int n, int *pos, MPI_Comm comm);
extern int p7_profile_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm,                          ESL_ALPHABET **byp_abc, P7_PROFILE **ret_gm);
extern int p7_profile_mpi_Recv(int source, int tag,          MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **byp_abc, P7_PROFILE **ret_gm);

#endif /*HAVE_MPI*/
#endif /*p7HMM_MPI_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
