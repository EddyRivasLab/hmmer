/* Optional support for P7_DOMAIN communication under MPI.
 * 
 * Currently, one or more P7_DOMAINs are packed/unpacked as
 * part of one P7_HIT; we do not need a send/recv, only
 * pack/unpack. So, see p7_tophits_mpi for the send/recv.
 * 
 * Depends on p7_alidisplay_mpi, because it packs the alidisplay,
 * in addition to the domain coord info.
 * 
 * Contents:
 *    1. Pack/unpack of P7_DOMAIN, one domain hit and its alignment
 *    2. Unit tests
 *    3. Test driver
 */
#include "p7_config.h"
#ifdef HAVE_MPI
#include <mpi.h>

#include "easel.h"

#include "base/p7_domain.h"
#include "base/p7_domain_mpi.h"
#include "base/p7_alidisplay_mpi.h"

/*****************************************************************
 * 1. Pack/unpack of P7_DOMAIN
 *****************************************************************/

/* Function:  p7_domain_mpi_PackSize()
 * Synopsis:  Calculates size needed to pack a P7_DOMAIN array
 *
 * Purpose:   Calculate an upper bound on the number of bytes that
 *            <p7_domain_mpi_Pack()> will need to pack an array of
 *            <ndom> domain coords <dcl> (and the associated
 *            alidisplay for each domain) in a packed MPI message for
 *            MPI communicator <comm>. Return that number of bytes in
 *            <*ret_n>.
 *            
 *            <ndom> may be 0 (and <dcl> <NULL>), in which case
 *            <*ret_n> is returned as 0.
 *            
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is undefined.
 */
int
p7_domain_mpi_PackSize(const P7_DOMAIN *dcl, int ndom, MPI_Comm comm, int *ret_n)
{
  int n = 0;
  int sz;
  int d;
  int status;
  
  for (d = 0; d < ndom; d++)
    {
      if ( MPI_Pack_size( 1, MPI_INT,    comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz*10; /* iae,ibe;kae,kbe;ia,ib;ka,kb;is_reported,is_included*/
      if ( MPI_Pack_size( 1, MPI_FLOAT,  comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz*5;  /* envsc,domcorrection,dombias,oasc,bitscore*/
      if ( MPI_Pack_size( 1, MPI_DOUBLE, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed"); n += sz;    /* lnP */

      if ( (status = p7_alidisplay_mpi_PackSize(dcl[d].ad, comm, &sz)) != eslOK) return status;                 n += sz;
      // fixme: TJW is using scores_per_pos; SRE is not dealing with it
    }
  *ret_n = n;
  return eslOK;
}

/* Function:  p7_domain_mpi_Pack()
 * Synopsis:  Pack a P7_DOMAIN array into a buffer.
 *
 * Purpose:   Packs an array of <ndom> domain coords <dcl> (including
 *            the associated alidisplay for each domain) into an MPI
 *            packed message buffer <buf> of length <n> bytes,
 *            starting at byte position <*pos>, for MPI communicator
 *            <comm>. After packing the data, update <*pos> to the
 *            next byte in the buffer.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed data at
 *            position <*pos>. This typically requires a call to
 *            <p7_domain_mpi_PackSize()> first, and reallocation if
 *            needed.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <dcl> array, and <*pos> is set to the byte
 *            immediately following the last byte of the data
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <dcl> into <buf>. In either case, the state of
 *            <buf> and <*pos> is undefined, and both should
 *            be considered to be corrupted.
 */
int
p7_domain_mpi_Pack(const P7_DOMAIN *dcl, int ndom, char *buf, int n, int *pos, MPI_Comm comm)
{
  int d;
  int status;

  for (d = 0; d < ndom; d++)
    {
      /* (void *) casts are to suppress compiler warnings about dropping the const qualifier */
      if (MPI_Pack((void *) &(dcl[d].iae),           1, MPI_INT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].ibe),           1, MPI_INT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].kae),           1, MPI_INT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].kbe),           1, MPI_INT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].ia),            1, MPI_INT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].ib),            1, MPI_INT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].ka),            1, MPI_INT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].kb),            1, MPI_INT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].envsc),         1, MPI_FLOAT,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].domcorrection), 1, MPI_FLOAT,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].dombias),       1, MPI_FLOAT,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].oasc),          1, MPI_FLOAT,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].bitscore),      1, MPI_FLOAT,  buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].lnP),           1, MPI_DOUBLE, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].is_reported),   1, MPI_INT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");
      if (MPI_Pack((void *) &(dcl[d].is_included),   1, MPI_INT,    buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed");

      if (( status = p7_alidisplay_mpi_Pack(dcl[d].ad,      buf, n, pos, comm))!= eslOK)       return status;
      // fixme: TJW is using scores_per_pos; SRE is not dealing with it
    }
  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}

/* Function:  p7_domain_mpi_Unpack()
 * Synopsis:  Unpack an array of P7_DOMAINs from a buffer.
 *
 * Purpose:   Unpack an array of <ndom> domain coords (and associated
 *            alidisplays) from <buf>, starting from position <*pos>,
 *            where the total length of the buffer in bytes is
 *            <n>. The new domain array is allocated here, and
 *            returned in <*ret_dcl>, suitable for inclusion in a
 *            <P7_HIT> structure.
 *            
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_dcl>
 *            contains a newly allocated array of <P7_DOMAIN>, which the caller is
 *            responsible for free'ing. 
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_dcl> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
p7_domain_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, P7_DOMAIN **ret_dcl, int ndom)
{
  P7_DOMAIN *dcl = NULL;
  int        d;
  int        status;
  if(ndom == 0){  // Nick's testing code because H4 doesn't have real hits left
    *ret_dcl = NULL;
    return eslOK;
  }

  ESL_ALLOC(dcl, sizeof(P7_DOMAIN) * ndom);
  for (d = 0; d < ndom; d++) dcl[d].ad = NULL;

  for (d = 0; d < ndom; d++)
    {
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].iae),           1, MPI_INT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].ibe),           1, MPI_INT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].kae),           1, MPI_INT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].kbe),           1, MPI_INT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].ia),            1, MPI_INT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].ib),            1, MPI_INT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].ka),            1, MPI_INT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].kb),            1, MPI_INT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].envsc),         1, MPI_FLOAT,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].domcorrection), 1, MPI_FLOAT,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].dombias),       1, MPI_FLOAT,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].oasc),          1, MPI_FLOAT,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].bitscore),      1, MPI_FLOAT,  comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].lnP),           1, MPI_DOUBLE, comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].is_reported),   1, MPI_INT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");
      if ( MPI_Unpack( buf, n, pos, &(dcl[d].is_included),   1, MPI_INT,    comm) != MPI_SUCCESS) ESL_XEXCEPTION(eslESYS, "mpi unpack failed");

      dcl[d].scores_per_pos = NULL;       // fixme: TJW is using scores_per_pos; SRE is not dealing with it

      if (( status = p7_alidisplay_mpi_Unpack( buf, n, pos, comm, &(dcl[d].ad)))  != eslOK) goto ERROR;
    }
  *ret_dcl = dcl;
  return eslOK;
  
 ERROR:
  if (dcl) 
    {
      for (d = 0; d < ndom; d++) if (dcl[d].ad) p7_alidisplay_Destroy(dcl[d].ad);
      free(dcl);
    }
  *ret_dcl = NULL;
  return status;
}

/*****************************************************************
 * 2. Unit tests.
 *****************************************************************/
#ifdef p7DOMAIN_MPI_TESTDRIVE

#include "esl_random.h"

/* The pack/unpack test does no interprocess communication, so it can
 * run with any number of mpi processes, even just 1
 */
static void
utest_PackUnpack(ESL_RANDOMNESS *rng)
{
  char msg[]       = "utest_PackUnpack() failed";
  P7_DOMAIN *dcl1  = NULL;
  P7_DOMAIN *dcl2  = NULL;
  int        ndom  = 10;
  int        alen  = 100;
  int        d;
  int        n1;
  char      *buf   = NULL;
  int        pos   = 0;
  char       errbuf[eslERRBUFSIZE];
  int        status;

  if ( (dcl1 = p7_domain_Create(ndom)) == NULL) esl_fatal(msg);
  for (d = 0; d < ndom; d++)
    if (p7_domain_TestSample(rng, alen, &(dcl1[d])) != eslOK) esl_fatal(msg);
  
  if (p7_domain_mpi_PackSize(dcl1, ndom, MPI_COMM_WORLD, &n1) != eslOK) esl_fatal(msg);
  ESL_ALLOC(buf, sizeof(char) * n1);
  
  pos = 0;
  if (p7_domain_mpi_Pack(dcl1, ndom, buf, n1, &pos, MPI_COMM_WORLD) != eslOK) esl_fatal(msg);
  if (n1 != pos) esl_fatal(msg);

  pos = 0;
  if (p7_domain_mpi_Unpack(buf, n1, &pos, MPI_COMM_WORLD, &dcl2, ndom) != eslOK) esl_fatal(msg);
  if (n1 != pos) esl_fatal(msg);

  for (d = 0; d < ndom; d++)
    {
      if (p7_domain_Validate(&(dcl1[d]), errbuf) != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
      if (p7_domain_Validate(&(dcl2[d]), errbuf) != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
      if (p7_domain_Compare( &(dcl1[d]), &(dcl2[d]), 1e-6) != eslOK) esl_fatal(msg);
    }
  p7_domain_Destroy(dcl1, ndom);
  p7_domain_Destroy(dcl2, ndom);
  free(buf);
  return;

 ERROR:
  esl_fatal(msg);
}
#endif /*p7DOMAIN_MPI_TESTDRIVE*/
/*--------------- end, unit tests -------------------------------*/

/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef p7DOMAIN_MPI_TESTDRIVE
#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                        docgroup*/
  { (char *) "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, (char *) "show brief help on version and usage",            0 },
  { (char *) "-s",        eslARG_INT,     (char *)  "0", NULL, NULL,  NULL,  NULL, NULL, (char *) "set random number seed to <n>",                   0 },
  { (char *) "--stall",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, (char *) "arrest after start: for debugging MPI under gdb", 0 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for p7_domain_mpi.c domain coord MPI communication routines";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng      = NULL;
  int             stalling = esl_opt_GetBoolean(go,(char *)  "--stall");
  int             my_rank;
  int             nproc;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (my_rank == 0) {
    rng = esl_randomness_CreateFast(esl_opt_GetInteger(go,(char *) "-s"));
    fprintf(stderr, "## %s\n", argv[0]);
    fprintf(stderr, "#  rng seed  = %" PRIu32 "\n", esl_randomness_GetSeed(rng));
    fprintf(stderr, "#  MPI nproc = %d\n", nproc);
  }
#ifdef HAVE_GETPID
  fprintf(stderr, "#    %6d = %d\n", my_rank, getpid()); 
#endif
  while (stalling);	

  if (my_rank == 0) 
    utest_PackUnpack(rng);
  
  if (my_rank == 0) {
    fprintf(stderr, "#  status = ok\n");
    esl_randomness_Destroy(rng);
  }

  MPI_Finalize();
  esl_getopts_Destroy(go);
  exit(0); /* success */
}
#endif /*p7DOMAIN_MPI_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/



#else /*! HAVE_MPI*/
/* If we don't have MPI compiled in, provide some nothingness to:
 *   a. prevent Mac OS/X ranlib from bitching about .o file that "has no symbols" 
 *   b. prevent compiler from bitching about "empty compilation unit"
 *   c. automatically pass the automated tests.
 */
void p7_domain_mpi_DoAbsolutelyNothing(void) { return; }
#if defined p7DOMAIN_MPI_TESTDRIVE || p7DOMAIN_MPI_EXAMPLE || p7DOMAIN_MPI_BENCHMARK
int main(void) { return 0; }
#endif
#endif /*HAVE_MPI*/


