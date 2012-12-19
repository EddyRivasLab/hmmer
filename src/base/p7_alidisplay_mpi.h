#ifndef P7_ALIDISPLAY_MPI_INCLUDED
#define P7_ALIDISPLAY_MPI_INCLUDED
#include "p7_config.h"

#ifdef  HAVE_MPI
#include <mpi.h>

#include "base/p7_alidisplay.h"

extern int p7_alidisplay_mpi_PackSize(const P7_ALIDISPLAY *ad, MPI_Comm comm, int *ret_n);
extern int p7_alidisplay_mpi_Pack    (const P7_ALIDISPLAY *ad, char *buf, int n, int *pos, MPI_Comm comm);
extern int p7_alidisplay_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, P7_ALIDISPLAY **ret_ad);

#endif /*HAVE_MPI*/
#endif /*P7_ALIDISPLAY_MPI_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
