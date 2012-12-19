#ifndef P7_DOMAIN_MPI_INCLUDED
#define P7_DOMAIN_MPI_INCLUDED
#include "p7_config.h"

#ifdef HAVE_MPI
#include <mpi.h>

#include "p7_domain.h"

extern int p7_domain_mpi_PackSize(const P7_DOMAIN *dcl, int ndom, MPI_Comm comm, int *ret_n);
extern int p7_domain_mpi_Pack    (const P7_DOMAIN *dcl, int ndom, char *buf, int n, int *pos, MPI_Comm comm);
extern int p7_domain_mpi_Unpack  (char *buf, int n, int *pos, MPI_Comm comm, P7_DOMAIN **ret_dcl, int ndom);

#endif /*HAVE_MPI*/
#endif /*P7_DOMAIN_MPI_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

