/* Optional support for P7_PIPELINE_STATS communication under MPI.
 * 
 * (Right now the only thing we need to communicate in P7_PIPELINE is
 * the stats chunk of a pipeline in an MPI worker, which gets sent
 * back to be merged into the master at the end of the worker's work.)
 * 
 * This module follows a framework laid down for more complicated
 * structure communication, as in p7_hmm_mpi.c and
 * p7_profile_mpi.c. That design is overkill for a little fixed-size
 * structure like P7_PIPELINE_STATS, but we use it anyway, for
 * consistency and in case of future expansions.
 * 
 * Contents:
 *    1. Communicating P7_PIPELINE_STATS, accumulated accounting statistics.
 *    2. Unit tests.
 *    3. Test driver.
 */
#include "p7_config.h"		

#ifdef HAVE_MPI

#include <mpi.h>

#include "search/p7_pipeline.h"
#include "search/p7_pipeline_mpi.h"

/* Function:  p7_pipeline_stats_mpi_Send()
 * Synopsis:  Send a P7_PIPELINE_STATS back to the master.
 *
 * Purpose:   Sends a P7_PIPELINE_STATS as a packed unit to MPI process
 *            <dest> (where <dest> ranges from 0..<nproc-1>), tagged
 *            with MPI tag <tag>, for MPI communicator <comm>.
 *            
 *            In order to minimize alloc/free cycles in this routine,
 *            caller passes a pointer to a working buffer <*buf> of
 *            size <*nalloc> characters. If necessary (i.e. if <hmm> is
 *            too big to fit), <*buf> will be reallocated and <*nalloc>
 *            increased to the new size. As a special case, if <*buf>
 *            is <NULL> and <*nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *            
 * Returns:   <eslOK> on success; <*buf> may have been reallocated and
 *            <*nalloc> may have been increased.
 * 
 * Throws:    <eslESYS> if an MPI call fails; <eslEMEM> if a malloc/realloc
 *            fails. In either case, <*buf> and <*nalloc> remain valid and useful
 *            memory (though the contents of <*buf> are undefined). 
 */
int
p7_pipeline_stats_mpi_Send(const P7_PIPELINE_STATS *stats, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int n = 0;
  int code;
  int sz, pos;
  int status;

  /* Figure out size, including a status code (0=EOD=nothing sent) */
  if ( MPI_Pack_size(1, MPI_INT, comm, &sz)                       != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack size failed");  n += sz;
  if ((status = p7_pipeline_stats_mpi_PackSize(stats, comm, &sz)) != eslOK)       return status;                                   n += sz;

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) 
    {
      ESL_REALLOC(*buf, sizeof(char) * n);
      *nalloc = n; 
    }

  /* Pack the status code and structure into the buffer */
  /* The status code is the # of structures being sent as one MPI message; here 1 or 0 */
  pos  = 0;
  code = (stats ? 1 : 0);
  if (MPI_Pack(&code, 1, MPI_INT,                 *buf, n, &pos, comm)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi pack failed");
  if ((status = p7_pipeline_stats_mpi_Pack(stats, *buf, n, &pos, comm)) != eslOK)       return status;
  
  /* Send the packed buffer to the destination. */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != MPI_SUCCESS)  ESL_EXCEPTION(eslESYS, "mpi send failed");
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_pipeline_stats_mpi_PackSize()
 * Synopsis:  Calculates size (in bytes) needed to pack a P7_PIPELINE_STATS.
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the size in bytes.
 *
 * Throws:    <eslESYS> if MPI call fails, and <*ret_n> is undefined.
 */
int
p7_pipeline_stats_mpi_PackSize(const P7_PIPELINE_STATS *stats, MPI_Comm comm, int *ret_n)
{
  ESL_UNUSED(stats);		/* maybe we'll need <stats> someday, but for now, size is a predictable constant */
  int n = 0;
  int sz;
  
  if (stats) {
    if (MPI_Pack_size(1, MPI_UINT64_T, comm, &sz) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack size failed");  n += sz*14;
  }
  
  *ret_n = n;
  return eslOK;
}

/* Function:  p7_pipeline_stats_mpi_Pack()
 * Synopsis:  Pack a P7_PIPELINE_STATS into an MPI buffer.
 *
 * Purpose:   Packs <stats> into an MPI packed message buffer <buf>
 *            of length <n> bytes, starting at byte position <*pos>,
 *            for MPI communicator <comm>.
 *            
 *            The caller must know that <buf>'s allocation of <n>
 *            bytes is large enough to append the packed HMM at
 *            position <*pos>. This typically requires a call to
 *            <p7_pipeline_stats_mpi_PackSize()> first, and reallocation if
 *            needed.
 *            
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <stats>, and <*pos> is set to the byte
 *            immediately following the last byte of the stats
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <stats> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 */
int
p7_pipeline_stats_mpi_Pack(const P7_PIPELINE_STATS *stats, char *buf, int n, int *pos, MPI_Comm comm)
{
  if (stats) 
    {
      if (MPI_Pack((void *) &stats->nmodels,       1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &stats->nseqs,         1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &stats->nres,          1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &stats->nnodes,        1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 

      if (MPI_Pack((void *) &stats->n_past_ssv,    1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &stats->n_past_bias,   1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &stats->n_past_vit,    1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &stats->n_past_fwd,    1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 

      if (MPI_Pack((void *) &stats->n_output,      1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &stats->pos_past_msv,  1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &stats->pos_past_bias, 1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &stats->pos_past_vit,  1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &stats->pos_past_fwd,  1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
      if (MPI_Pack((void *) &stats->pos_output,    1, MPI_UINT64_T, buf, n, pos, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "pack failed"); 
    }
  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
}


/* Function:  p7_pipeline_stats_mpi_Unpack()
 * Synopsis:  Unpack one P7_PIPELINE_STATS struct from MPI buffer.
 *
 * Purpose:   Unpack one <P7_PIPELINE_STATS> structure from MPI packed buffer
 *            <buf>, starting from position <*pos>, where the total length
 *            of the buffer in bytes is <n>.
 *            
 *            Caller provides the space for <*stats>, which is
 *            probably just a reference to a structure on the stack
 *            (P7_PIPELINE_STATS is a fixed-size structure at
 *            present). If the structure becomes more complex in the
 *            future, we would allocate it here, as in our other MPI
 *            communication <_Unpack()> routines.
 *            
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*stats>
 *            is filled in with the communicated <P7_PIPELINE_STATS>.

 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, the contents of <*stats> is undefined; <*pos> 
 *            is undefined too, meaning that <*buf> is rendered useless.
 */
int
p7_pipeline_stats_mpi_Unpack(char *buf, int n, int *pos, MPI_Comm comm, P7_PIPELINE_STATS *stats)
{
  if (MPI_Unpack(buf, n, pos, &(stats->nmodels),       1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &(stats->nseqs),         1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &(stats->nres),          1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &(stats->nnodes),        1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 

  if (MPI_Unpack(buf, n, pos, &(stats->n_past_ssv),    1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &(stats->n_past_bias),   1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &(stats->n_past_vit),    1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &(stats->n_past_fwd),    1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 

  if (MPI_Unpack(buf, n, pos, &(stats->n_output),      1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &(stats->pos_past_msv),  1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &(stats->pos_past_bias), 1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &(stats->pos_past_vit),  1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &(stats->pos_past_fwd),  1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &(stats->pos_output),    1, MPI_UINT64_T, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "unpack failed"); 

  return eslOK;
} 

/* Function:  p7_pipeline_stats_mpi_Recv()
 * Synopsis:  Receive a P7_PIPELINE_STATS structure from an MPI sender.
 *
 * Purpose:   Receive a unit that consists of a single <P7_PIPELINE_STATS>
 *            structure sent by MPI <source> (<0..nproc-1>, or
 *            <MPI_ANY_SOURCE>) tagged as <tag> for MPI communicator <comm>.
 *            
 *            Caller provides space for <*stats>, usually just as a
 *            reference to a <P7_PIPELINE_STATS> structure on the
 *            stack, since this is currently a fixed-size structure.
 *            (If the <P7_PIPELINE_STATS> structure gets more complex
 *            in the future, this may change, to allocate the structure
 *            here as in MPI communication of <P7_HMM> or <P7_PROFILE>.)
 *            
 *            Units are prefixed by a status code that gives the
 *            number of structures to follow; here, 0 or 1.  If we
 *            receive a 1 code and we successfully unpack a stats
 *            structure, this routine will return <eslOK> with the
 *            stats in <*stats>.  If we receive a 0 code (a shutdown
 *            signal), this routine returns <eslEOD> and <*stats> is
 *            initialized to all zeros.
 *   
 *            Caller provides a working buffer <*buf> of size
 *            <*nalloc> characters. These are passed by reference, so
 *            that <*buf> can be reallocated and <*nalloc> increased
 *            if necessary. As a special case, if <*buf> is <NULL> and
 *            <*nalloc> is 0, the buffer will be allocated
 *            appropriately, but the caller is still responsible for
 *            free'ing it.
 *
 * Returns:   <eslOK> if stats are received in <*stats>.
 *  
 *            <eslEOD> if the communication unit was an end-of-data
 *            signal, with no communicated statistics; <*stats> is set
 *            to all zeros.
 *
 * Throws:    <eslEMEM> on allocation error; <eslESYS> on MPI communication
 *            error; in either case, contents of <*stats> are undefined.
 */
int
p7_pipeline_stats_mpi_Recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, P7_PIPELINE_STATS *stats)
{
  int        pos = 0;
  int        code;
  int        n;
  MPI_Status mpistatus;
  int        status;
  
  /* Probe first, because we need to know if our buffer is big enough. */
  if ( MPI_Probe(source, tag, comm, &mpistatus)  != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi probe failed");
  if ( MPI_Get_count(&mpistatus, MPI_PACKED, &n) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi get count failed");

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) 
    {
      ESL_REALLOC(*buf, sizeof(char) * n);
      *nalloc = n; 
    }

  /* Receive the entire packed work unit */
  if (MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi recv failed");

  /* Unpack the status code prefix */
  if (MPI_Unpack(*buf, n, &pos, &code, 1, MPI_INT, comm) != MPI_SUCCESS) ESL_EXCEPTION(eslESYS, "mpi unpack failed");

  if      (code == 0) { status = eslEOD; p7_pipeline_stats_Init(stats); }
  else if (code == 1)   status = p7_pipeline_stats_mpi_Unpack(*buf, *nalloc, &pos, comm, stats);
  else                  ESL_EXCEPTION(eslESYS, "bad mpi buffer transmission code");
  return status;

 ERROR: /* from ESL_REALLOC only */
  return status;
}
/*--------- end, P7_PIPELINE_STATS communication ----------------*/


/*****************************************************************
 * 2. Unit tests.
 *****************************************************************/
#ifdef p7PIPELINE_MPI_TESTDRIVE
#include "esl_random.h"

/* Solely for unit testing */
static int
pipeline_stats_Sample(ESL_RANDOMNESS *rng, P7_PIPELINE_STATS *stats)
{
  stats->nmodels       = esl_rnd_Roll(rng, 100);
  stats->nseqs         = esl_rnd_Roll(rng, 100);
  stats->nres          = esl_rnd_Roll(rng, 100);
  stats->nnodes        = esl_rnd_Roll(rng, 100);
  stats->n_past_ssv    = esl_rnd_Roll(rng, 100);
  stats->n_past_bias   = esl_rnd_Roll(rng, 100);
  stats->n_past_vit    = esl_rnd_Roll(rng, 100);
  stats->n_past_fwd    = esl_rnd_Roll(rng, 100);
  stats->n_output      = esl_rnd_Roll(rng, 100);
  stats->pos_past_msv  = esl_rnd_Roll(rng, 100);
  stats->pos_past_bias = esl_rnd_Roll(rng, 100);
  stats->pos_past_vit  = esl_rnd_Roll(rng, 100);
  stats->pos_past_fwd  = esl_rnd_Roll(rng, 100);
  stats->pos_output    = esl_rnd_Roll(rng, 100);
  return eslOK;
}

static int
pipeline_stats_Compare(P7_PIPELINE_STATS *s1, P7_PIPELINE_STATS *s2)
{
  if (s1->nmodels       != s2->nmodels)       return eslFAIL;
  if (s1->nseqs         != s2->nseqs)         return eslFAIL;
  if (s1->nres          != s2->nres)          return eslFAIL;
  if (s1->nnodes        != s2->nnodes)        return eslFAIL;
  if (s1->n_past_ssv    != s2->n_past_ssv)    return eslFAIL;
  if (s1->n_past_bias   != s2->n_past_bias)   return eslFAIL;
  if (s1->n_past_vit    != s2->n_past_vit)    return eslFAIL;
  if (s1->n_past_fwd    != s2->n_past_fwd)    return eslFAIL;
  if (s1->n_output      != s2->n_output)      return eslFAIL;
  if (s1->pos_past_msv  != s2->pos_past_msv)  return eslFAIL;
  if (s1->pos_past_bias != s2->pos_past_bias) return eslFAIL;
  if (s1->pos_past_vit  != s2->pos_past_vit)  return eslFAIL;
  if (s1->pos_past_fwd  != s2->pos_past_fwd)  return eslFAIL;
  if (s1->pos_output    != s2->pos_output)    return eslFAIL;
  return eslOK;
}

static void
utest_SendRecv(ESL_RANDOMNESS *rng, int my_rank, int nproc)
{
  char              msg[]   = "utest_SendRecv() failed";
  P7_PIPELINE_STATS mstats;	/* master's copy */
  P7_PIPELINE_STATS wstats;	/* worker's copy */
  char             *wbuf    = NULL;
  int               wn      = 0;
  uint32_t          rngseed;
  int               i;
  MPI_Status        mpistatus;
  
  if (my_rank == 0)
    {
      /* First we send our RNG seed to all workers */
      rngseed = esl_randomness_GetSeed(rng);
      for (i = 1; i < nproc; i++)
	if (MPI_Send( &rngseed, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD) != MPI_SUCCESS) esl_fatal(msg);

      /* Sample a made-up stats structure that will be identical to the workers' */
      if (pipeline_stats_Sample(rng, &mstats) != eslOK) esl_fatal(msg);

      /* Receive workers' structures, and compare. */
      for (i = 1; i < nproc; i++)
	{
	  if (p7_pipeline_stats_mpi_Recv(MPI_ANY_SOURCE, /*tag=*/0, MPI_COMM_WORLD, &wbuf, &wn, &wstats) != eslOK) esl_fatal(msg);
	  if (pipeline_stats_Compare(&mstats, &wstats)                                                   != eslOK) esl_fatal(msg);
	}
    }
  else
    {
      /* Worker(s) must first receive the exact same RNG seed that the master is using. */
      if (MPI_Recv(&rngseed, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &mpistatus) != MPI_SUCCESS) esl_fatal(msg);

      /* and then the worker(s) can create the exact same RNG (and random number sequence) that the master has */
      rng = esl_randomness_CreateFast(rngseed);                                     

      /* so when the worker samples its structure, the master has independently sampled an exact duplicate of it... */
      if (pipeline_stats_Sample(rng, &wstats) != eslOK) esl_fatal(msg);

      /* each worker sends the structure to the master (it's the same HMM for each worker. The test is intended for one master, one worker.) */
      if (p7_pipeline_stats_mpi_Send(&wstats, /*dest=*/0, /*tag=*/0, MPI_COMM_WORLD, &wbuf, &wn) != eslOK) esl_fatal(msg);

      /* worker's RNG is a private copy; destroy it. Master keeps its RNG, which the caller is responsible for. */
      esl_randomness_Destroy(rng);                                                  
     }
  free(wbuf);
  return;
}
#endif /*p7PIPELINE_MPI_TESTDRIVE*/


/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef p7PIPELINE_MPI_TESTDRIVE
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
static char banner[] = "unit test driver for p7_pipeline_mpi.c, pipeline stats MPI communication routines";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng      = NULL;
  int             stalling = esl_opt_GetBoolean(go, "--stall");
  int             my_rank;
  int             nproc;

  /* For more notes, see p7_hmm_mpi.c, our most-documented-model for an MPI unit test */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (my_rank == 0) {
    rng = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
    fprintf(stderr, "## %s\n", argv[0]);
    fprintf(stderr, "#  rng seed  = %" PRIu32 "\n", esl_randomness_GetSeed(rng));
    fprintf(stderr, "#  MPI nproc = %d\n", nproc);
  }
#ifdef HAVE_GETPID
  fprintf(stderr, "#    %6d = %d\n", my_rank, getpid()); 
#endif
  while (stalling);	

  utest_SendRecv(rng, my_rank, nproc);
  
  if (my_rank == 0) {
    fprintf(stderr, "#  status = ok\n");
    esl_randomness_Destroy(rng);
  }

  MPI_Finalize();
  esl_getopts_Destroy(go);
  exit(0); /* success */
}
#endif /*p7PIPELINE_MPI_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/

#else /*! HAVE_MPI*/
/* If we don't have MPI compiled in, provide some nothingness to:
 *   a. prevent Mac OS/X ranlib from bitching about .o file that "has no symbols" 
 *   b. prevent compiler from bitching about "empty compilation unit"
 *   c. automatically pass the automated tests.
 */
void p7_pipeline_mpi_DoAbsolutelyNothing(void) { return; }
#if defined p7PIPELINE_MPI_TESTDRIVE || p7PIPELINE_MPI_EXAMPLE || p7PIPELINE_MPI_BENCHMARK
int main(void) { return 0; }
#endif
#endif /*HAVE_MPI*/

