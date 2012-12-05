#ifndef P7_OPROFILE_MPI_INCLUDED
#define P7_OPROFILE_MPI_INCLUDED

#ifdef HAVE_MPI

#include "p7_config.h"

#include <mpi.h>

#include "esl_alphabet.h"

#include "dp_vector/p7_oprofile.h"


extern int p7_oprofile_mpi_Send(P7_OPROFILE *om, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int p7_oprofile_mpi_PackSize(P7_OPROFILE *om, MPI_Comm comm, int *ret_n);
extern int p7_oprofile_mpi_Pack(P7_OPROFILE *om, char *buf, int n, int *pos, MPI_Comm comm);
extern int p7_oprofile_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **abc, P7_OPROFILE **ret_om);
extern int p7_oprofile_mpi_Recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, P7_OPROFILE **ret_om);

#endif /*HAVE_MPI*/
#endif /*P7_OPROFILE_MPI_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
