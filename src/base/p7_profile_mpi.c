/* Optional support for P7_PROFILE communication under MPI.
 * 
 * Contents:
 *    1. Communicating P7_PROFILE, a score profile.
 *    2. Unit tests.
 *    3. Test driver.
 *    4. Copyright and license information.
 */

#include "p7_config.h"		

#ifdef HAVE_MPI

#include <mpi.h>

#include "easel.h"
#include "esl_mpi.h"

#include "base/p7_profile.h"
#include "base/p7_profile_mpi.h"


/*****************************************************************
 * 1. Communicating P7_PROFILE, a score profile.
 *****************************************************************/

/* Function:  p7_profile_mpi_Send()
 * Synopsis:  Send a profile as an MPI message.
 *
 * Purpose:   Sends profile <gm> to MPI process <dest> (where
 *            <dest> ranges from 0..<nproc-1>), with MPI tag <tag> 
 *            for MPI communicator <comm>.
 *            
 *            In order to minimize alloc/free cycles in this routine,
 *            caller passes a pointer to a working buffer <*buf> of
 *            size <*nalloc> characters. If necessary (i.e. if <gm> is
 *            too big to fit), <*buf> will be reallocated and <*nalloc>
 *            increased to the new size. As a special case, if <*buf>
 *            is <NULL> and <*nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *            
 *            If <gm> is NULL, an end-of-data signal is sent, which
 *            <p7_profile_mpi_Recv()> knows how to interpret.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Note:      This was tested against a version that simply issues a series
 *            of MPI_Send()'s, rather than Pack()'ing into a buffer
 *            and issuing one MPI_Send(). The packed version seems to
 *            be significantly faster, although benchmarking MPI
 *            programs is difficult, and variance on the results is high.
 *            
 *            To optimize communication still further, one might try
 *            to avoid many or all of the MPI_Pack()'s. It might be
 *            feasible to change the allocation of a profile such that
 *            it is allocated in one contiguous chunk of memory. And
 *            once one embarks on that, the memory layout of the
 *            profile should also be optimized with respect to cache
 *            performance during DP alignment.
 *            
 *            rf, cs annotation is optional, but for simplicity, we always 
 *            transmit the two allocated strings, even if they were empty.
 */
int
p7_profile_mpi_Send(P7_PROFILE *gm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   sz, pos;
  int   code;
  int   n  = 0;
  int   Kp = 0;	/* alphabet size including degeneracies */
  int   M  = 0; /* model size in nodes */
  int   status;

  /* Figure out size. */
  if ( MPI_Pack_size(1, MPI_INT, comm, &sz)              != 0)     ESL_EXCEPTION(eslESYS, "mpi pack size failed");  n += sz; /* # of profiles in the transmission; 0=EOD */
  if (( status = p7_profile_mpi_PackSize(gm, comm, &sz)) != eslOK) return status;                                   n += sz;

  /* Allocate buffer */
  if (*buf == NULL || n > *nalloc) 
    {
      ESL_REALLOC(*buf, sizeof(char) * n);
      *nalloc = n; 
    }

  /* Pack the profile into the buffer */
  pos  = 0;
  code = (gm ? 1 : 0);
  if (MPI_Pack(&code, 1, MPI_INT, *buf, n, &pos,  comm)            != 0)     ESL_XEXCEPTION(eslESYS, "pack failed");
  if (gm && (status = p7_profile_mpi_Pack(gm, buf, n, &pos, comm)) != eslOK) ESL_XEXCEPTION(status,  "pack failed");

  /* Send the packed profile to destination  */
  MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm);
  return eslOK;
  
 ERROR:
  return status;
}


/* Function:  p7_profile_mpi_PackSize()
 * Synopsis:  Calculates size (in bytes) needed to pack a P7_PROFILE
 *
 * Purpose:   Calculate an upper bound on the number of bytes
 *            that <p7_profile_mpi_Pack()> will need to pack a
 *            <P7_PROFILE> <om> in a packed MPI message for
 *            MPI communicator <comm>. Return that number of 
 *            bytes in <*ret_n>.
 *            
 *            <gm> may be <NULL>, in which case <*ret_n> is 
 *            returned as 0.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 * 
 * Throws:    <eslESYS> if an MPI call fails; now <*ret_n> is 0.           
 */
int
p7_profile_mpi_PackSize(P7_PROFILE *gm, MPI_Comm comm, int *ret_n)
{
  int   n  = 0;	                        /* result: total number of bytes needed    */
  int   Kp = (gm ? gm->abc->Kp : 0);	/* alphabet size including degeneracies    */
  int   M  = (gm ? gm->M       : 0);    /* model size in nodes                     */
  int   sz;	                        /* size of some subpiece of data structure */
  int   status;

  /* Recall that MPI_Pack_size(x, ...) + MPI_Pack_size(x, ... ) != MPI_Pack_size(2x, ...).
   * Pack_size() units must exactly match to our Pack() calls.
   */
  if (gm)
    {
      if (MPI_Pack_size(                           1, MPI_INT,      comm, &sz) != 0)     ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz*4;               /* M,L,max_length,abc->type */
      if (MPI_Pack_size(         (M+1) *  p7P_NTRANS, MPI_FLOAT,    comm, &sz) != 0)     ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz;                 /* tsc             */
      if (MPI_Pack_size(         (M+1) * Kp * p7P_NR, MPI_FLOAT,    comm, &sz) != 0)     ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz;                 /* rsc[0]          */
      if (MPI_Pack_size(                 p7P_NXTRANS, MPI_FLOAT,    comm, &sz) != 0)     ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz*p7P_NXSTATES;    /* xsc[0..3]       */
      if (MPI_Pack_size(                           1, MPI_FLOAT,    comm, &sz) != 0)     ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz*2;               /* nj,pglocal      */
      if ((status = esl_mpi_PackOptSize(gm->name, -1, MPI_CHAR,     comm, &sz))!= eslOK) return status;                               n += sz;                 /* name (string)   */
      if ((status = esl_mpi_PackOptSize(gm->acc,  -1, MPI_CHAR,     comm, &sz))!= eslOK) return status;                               n += sz;                 /* acc (string)    */
      if ((status = esl_mpi_PackOptSize(gm->desc, -1, MPI_CHAR,     comm, &sz))!= eslOK) return status;                               n += sz;                 /* desc (string)   */
      if (MPI_Pack_size(                       (M+2), MPI_CHAR,     comm, &sz) != 0)     ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz*4;               /* rf,cs,mm,consensus */
      if (MPI_Pack_size(                 p7_NEVPARAM, MPI_FLOAT,    comm, &sz) != 0)     ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz;                 /* evparam         */
      if (MPI_Pack_size(                 p7_NCUTOFFS, MPI_FLOAT,    comm, &sz) != 0)     ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz;                 /* Pfam cutoffs    */
      if (MPI_Pack_size(                 p7_NOFFSETS, MPI_LONG_LONG,comm, &sz) != 0)     ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz;                 /* offs[]          */
      if (MPI_Pack_size(                           1, MPI_LONG_LONG,comm, &sz) != 0)     ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz*2;               /* roff, eoff      */
      /* <allocM> is not sent. */
      /* alphabet is only sent as <type>. */
    }
  *ret_n = n;
  return eslOK;
}



/* Function:  p7_profile_mpi_Pack()
 * Synopsis:  Pack a profile into an MPI buffer.
 *
 * Purpose:   Pack profile <gm> into an MPI packed message buffer <buf>
 *            of length <n> bytes, starting at byte position <*pos>,
 *            for MPI communicator <comm>.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed profile at
 *            position <*pos>. This typically requires a call to
 *            <p7_hmm_mpi_PackSize()> first, and reallocation if
 *            needed.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <gm>, and <*pos> is set to the byte
 *            immediately following the last byte of the HMM
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <gm> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 */
int
p7_profile_mpi_Pack(P7_PROFILE *gm, char *buf, int n, int *pos, MPI_Comm comm)
{
  int M  = (gm ? gm->M       : 0);
  int Kp = (gm ? gm->abc->Kp : 0);
  int x;

  if (gm) {
    if (MPI_Pack( &(gm->abc->type),              1, MPI_INT,       buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
    if (MPI_Pack( &(gm->M),                      1, MPI_INT,       buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack(   gm->tsc,    (M+1) * p7P_NTRANS, MPI_FLOAT,     buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack(   gm->rsc[0], (M+1)* Kp * p7P_NR, MPI_FLOAT,     buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    for (s = 0; s < p7P_NXSTATES; s++)
      if (MPI_Pack( gm->xsc[s],        p7P_NXTRANS, MPI_FLOAT,     buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack( &(gm->L),                      1, MPI_INT,       buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack( &(gm->nj),                     1, MPI_FLOAT,     buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack( &(gm->pglocal),                1, MPI_FLOAT,     buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if ((status = esl_mpi_PackOpt(gm->name,     -1, MPI_CHAR,      buf, n, pos, comm)) != eslOK) goto ERROR;
    if ((status = esl_mpi_PackOpt(gm->acc,      -1, MPI_CHAR,      buf, n, pos, comm)) != eslOK) goto ERROR; 
    if ((status = esl_mpi_PackOpt(gm->desc,     -1, MPI_CHAR,      buf, n, pos, comm)) != eslOK) goto ERROR;
    if (MPI_Pack(   gm->rf,                    M+2, MPI_CHAR,      buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack(   gm->mm,                    M+2, MPI_CHAR,      buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed");
    if (MPI_Pack(   gm->cs,                    M+2, MPI_CHAR,      buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack(   gm->consensus,             M+2, MPI_CHAR,      buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack(   gm->evparam,       p7_NEVPARAM, MPI_FLOAT,     buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack(   gm->cutoff,        p7_NCUTOFFS, MPI_FLOAT,     buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack(   gm->compo,          p7_MAXABET, MPI_FLOAT,     buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack(   gm->offs,          p7_NOFFSETS, MPI_LONG_LONG, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack( &(gm->roff),                   1, MPI_LONG_LONG, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack( &(gm->eoff),                   1, MPI_LONG_LONG, buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
    if (MPI_Pack( &(gm->max_length),             1, MPI_INT,       buf, n, pos, comm)  != 0)     ESL_EXCEPTION(eslESYS, "pack failed"); 
  }
  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}


int
p7_profile_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **abc, P7_PROFILE **ret_gm)
{

  


}


/* Function:  p7_profile_mpi_Recv()
 * Synopsis:  Receive a profile as an MPI message.
 *
 * Purpose:   Receive a profile from <source> (where <source> is usually
 *            process 0, the master) with tag <tag> from communicator <comm>,
 *            and return it in <*ret_gm>. 
 *            
 *            Caller must also provide the alphabet <abc> and the
 *            background model <bg> for this profile. (Of course, that means
 *            the caller already knows them, by an appropriate
 *            initialization.)
 *            
 *            To minimize alloc/free cycles in this routine, caller
 *            passes a pointer to a buffer <*buf> of size <*nalloc>
 *            characters. These are passed by reference because if
 *            necessary, <buf> will be reallocated and <nalloc>
 *            increased to the new size. As a special case, if <buf>
 *            is <NULL> and <nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *
 *            If the packed profile is an end-of-data signal, return
 *            <eslEOD>, and <*ret_gm> is <NULL>.
 *            
 * Returns:   <eslOK> on success. <*ret_gm> contains the new profile; it
 *            is allocated here, and the caller is responsible for
 *            free'ing it.  <*buf> may have been reallocated to a
 *            larger size, and <*nalloc> may have been increased.
 *            
 */
int
p7_profile_mpi_Recv(int source, int tag, MPI_Comm comm, const ESL_ALPHABET *abc, const P7_BG *bg, char **buf, int *nalloc,  P7_PROFILE **ret_gm)
{
  int         status;
  P7_PROFILE *gm    = NULL;
  int         n;
  int         position;
  MPI_Status  mpistatus;
  int         M;

  /* Probe first, because we need to know if our buffer is big enough.
   */
  MPI_Probe(source, tag, comm, &mpistatus);
  MPI_Get_count(&mpistatus, MPI_PACKED, &n);

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    ESL_REALLOC(*buf, sizeof(char) * n); 
    *nalloc = n; 
  }

  /* Receive the packed profile */
  MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus);

  /* Unpack it - watching out for the EOD signal of M = -1. */
  position = 0;
  if (MPI_Unpack(*buf, n, &position, &M,                            1, MPI_INT,   comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (M == -1) { *ret_gm = NULL; return eslEOD; }

  if ((gm = p7_profile_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR; }
  if (MPI_Unpack(*buf, n, &position, &(gm->mode),                   1, MPI_INT,   comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &position, &(gm->L),                      1, MPI_INT,   comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &position, gm->tsc,            p7P_NTRANS*M, MPI_FLOAT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &position, gm->rsc[0], p7P_NR*(M+1)*abc->Kp, MPI_FLOAT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &position, gm->xsc[0],          p7P_NXTRANS, MPI_FLOAT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");  
  if (MPI_Unpack(*buf, n, &position, gm->xsc[1],          p7P_NXTRANS, MPI_FLOAT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");  
  if (MPI_Unpack(*buf, n, &position, gm->xsc[2],          p7P_NXTRANS, MPI_FLOAT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");  
  if (MPI_Unpack(*buf, n, &position, gm->xsc[3],          p7P_NXTRANS, MPI_FLOAT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");  
  if (MPI_Unpack(*buf, n, &position, &(gm->nj),                     1, MPI_FLOAT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");  

  if ((status = esl_mpi_UnpackOpt(  *buf, n, &position,  (void**)&(gm->name),  NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if ((status = esl_mpi_UnpackOpt(  *buf, n, &position,  (void**)&(gm->acc),   NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;
  if ((status = esl_mpi_UnpackOpt(  *buf, n, &position,  (void**)&(gm->desc),  NULL, MPI_CHAR,  comm)) != eslOK) goto ERROR;

  if (MPI_Unpack(*buf, n, &position, gm->rf,                      M+2, MPI_CHAR,  comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &position, gm->mm,                      M+2, MPI_CHAR,  comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &position, gm->cs,                      M+2, MPI_CHAR,  comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &position, gm->consensus,               M+2, MPI_CHAR,  comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &position, gm->evparam,         p7_NEVPARAM, MPI_FLOAT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &position, gm->cutoff,          p7_NCUTOFFS, MPI_FLOAT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  if (MPI_Unpack(*buf, n, &position, gm->compo,           p7_NCUTOFFS, MPI_FLOAT, comm) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
  
  gm->abc = abc;
  gm->M   = M;
  *ret_gm = gm;
  return eslOK;

 ERROR:
  if (gm  != NULL) p7_profile_Destroy(gm);
  *ret_gm = NULL;
  return status;
}
/*--------------- end, P7_PROFILE communication -----------------*/


/*****************************************************************
 * 2. Unit tests.
 *****************************************************************/
#ifdef p7PROFILE_MPI_TESTDRIVE

static void
utest_SendRecv(ESL_RANDOMNESS *rng, int my_rank, int nproc)
{
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm  = NULL;
  P7_BG          *bg   = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_PROFILE     *gm2  = NULL;
  int             M    = 200;
  int             L    = 400;
  char           *wbuf = NULL;
  int             wn   = 0;
  int             i;
  char            errbuf[eslERRBUFSIZE];

  p7_hmm_Sample(rng, M, abc, &hmm); /* master and worker's sampled profiles are identical */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, L);
  p7_bg_SetLength  (bg, L);

  if (my_rank == 0)
    {
      for (i = 1; i < nproc; i++)
	{
	  ESL_DPRINTF1(("Master: receiving test profile\n"));
	  p7_profile_MPIRecv(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, abc, bg, &wbuf, &wn, &gm2);
	  ESL_DPRINTF1(("Master: test profile received\n"));

	  if (p7_profile_Validate(gm2, errbuf, 0.001) != eslOK) p7_Die("profile validation failed: %s", errbuf);
	  if (p7_profile_Compare(gm, gm2, 0.001) != eslOK) p7_Die("Received profile not identical to what was sent");

	  p7_profile_Destroy(gm2);
	}
    }
  else 
    {
      ESL_DPRINTF1(("Worker %d: sending test profile\n", my_rank));
      p7_profile_MPISend(gm, 0, 0, MPI_COMM_WORLD, &wbuf, &wn);
      ESL_DPRINTF1(("Worker %d: test profile sent\n", my_rank));
    }

  free(wbuf);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  return;
}
#endif /*p7PROFILE_MPI_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/


/*****************************************************************
 * 3. Test driver.
 *****************************************************************/
#ifdef p7PROFILE_MPI_TESTDRIVE
#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                        docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",            0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                   0 },
  { "--stall",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "arrest after start: for debugging MPI under gdb", 0 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for p7_profile_mpi.c core model MPI communication routines";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             stalling = esl_opt_GetBoolean(go, "--stall");
  int             my_rank;
  int             nproc;

  while (stalling);		/* wait for gdb to be attached, and operator sets stalling=FALSE */

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  utest_SendRecv(rng, my_rank, nproc);

  MPI_Finalize();
  fprintf(stderr, "#  status = ok\n");

  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  exit(0); /* success */
}
#endif /*p7PROFILE_MPI_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/

#else /*! HAVE_MPI*/
/* If we don't have MPI compiled in, provide some nothingness to:
 *   a. prevent Mac OS/X ranlib from bitching about .o file that "has no symbols" 
 *   b. prevent compiler from bitching about "empty compilation unit"
 *   c. automatically pass the automated tests.
 */
void p7_profile_mpi_DoAbsolutelyNothing(void) { return; }
#if defined p7PROFILE_MPI_TESTDRIVE || p7PROFILE_MPI_EXAMPLE || p7PROFILE_MPI_BENCHMARK
int main(void) { return 0; }
#endif
#endif /*HAVE_MPI*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
