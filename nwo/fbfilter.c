/* The Forward/Backward filter (FB).
 *
 * FB is the third vectorized DP filter in the HMMER acceleration
 * pipeline, after SSV and VF.
 *
 * The Forward pass of the FB filter computes a local (only) ensemble
 * Forward score. If this score is deemed sufficient, a
 * Backward/Decoding pass does posterior decoding to identify i,k
 * lattice cells (for sequence positions i and profile positions k)
 * with more than a threshold amount of probability mass. These cells
 * are marked in a data structure (H4_SPARSEMASK) for subsequent
 * sparse dynamic programming with the full glocal/local dual-mode
 * model.
 *
 * The FB filter is memory-efficient, using checkpointed dynamic
 * programming. It requires $O(M \sqrt L)$ memory for a profile of
 * length $M$ and a sequence of length $L$.
 *
 * The code here is only the runtime dispatcher. The actual
 * implementations are in fbfilter_{sse,avx,avx512...}.c.
 *
 * See fwdfilter.md for notes.
 *
 * Contents:
 *    1. h4_fwdfilter(), h4_bckfilter()
 *    2. CPU dispatching to vector implementations.
 *    3. Benchmark driver
 *    4. Unit tests
 *    5. Test driver
 *    6. Example
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "h4_profile.h"
#include "h4_mode.h"
#include "h4_checkptmx.h"
#include "h4_sparsemask.h"

#include "fbfilter.h"

static int fwdfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc);
#if 0
static int bckfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh);
#endif

/*****************************************************************
 * 1. h4_fwdfilter(), h4_bckfilter() 
 *****************************************************************/

int
(*h4_fwdfilter)(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc) =
  fwdfilter_dispatcher;

#if 0
int
(*h4_bckfilter)(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh) =
  bckfilter_dispatcher;
#endif


/*****************************************************************
 * 2. CPU dispatching to vector implementations.
 *****************************************************************/

static int 
fwdfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc)
{
  /* SRE: for now */
  h4_fwdfilter = h4_fwdfilter_sse;
  return h4_fwdfilter_sse(dsq, L, hmm, mo, cpx, opt_sc);

  /* When a platform supports more than one vector implementation, 
   * put fastest one first to prefer enabling it over others.
   */
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512())
    {
      h4_fwdfilter = h4_fwdfilter_avx512;
      return h4_fwdfilter_avx512(dsq, L, hmm, mo, cpx, opt_sc);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      h4_fwdfilter = h4_fwdfilter_avx;
      return h4_fwdfilter_avx(dsq, L, hmm, mo, cpx, opt_sc);
    }
#endif

#ifdef eslENABLE_SSE4
  if (esl_cpu_has_sse4())
    {
      h4_fwdfilter = h4_fwdfilter_sse;
      return h4_fwdfilter_sse(dsq, L, hmm, mo, cpx, opt_sc);
    }
#endif
  
#ifdef eslENABLE_NEON
  h4_fwdfilter = h4_fwdfilter_neon;
  return h4_fwdfilter_neon(dsq, L, hmm, mo, cpx, opt_sc);
#endif

#ifdef eslENABLE_VMX
  h4_fwdfilter = h4_fwdfilter_vmx;
  return h4_fwdfilter_vmx(dsq, L, hmm, mo, cpx, opt_sc);
#endif

  esl_fatal("fwdfilter_dispatcher found no vector implementation - that shouldn't happen.");
}



#if 0
static int
bckfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh)
{
  /* SRE: for now */
  h4_bckfilter = h4_bckfilter_sse;
  return h4_bckfilter_sse(dsq, L, hmm, mo, cpx, opt_sc);

#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512())
    {
      h4_bckfilter = h4_bckfilter_avx512;
      return h4_bckfilter_avx512(dsq, L, hmm, mo, cpx, sm, sm_thresh);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      h4_bckfilter = h4_bckfilter_avx;
      return h4_bckfilter_avx(dsq, L, hmm, mo, cpx, sm, sm_thresh);
    }
#endif

#ifdef eslENABLE_SSE4
  if (esl_cpu_has_sse4())
    {
      h4_bckfilter = h4_bckfilter_sse;
      return h4_bckfilter_sse(dsq, L, hmm, mo, cpx, sm, sm_thresh);
    }
#endif
  
#ifdef eslENABLE_NEON
  h4_bckfilter = h4_bckfilter_neon;
  return h4_bckfilter_neon(dsq, L, hmm, mo, cpx, sm, sm_thresh);
#endif

#ifdef eslENABLE_VMX
  h4_bckfilter = h4_bckfilter_vmx;
  return h4_bckfilter_vmx(dsq, L, hmm, mo, cpx, sm, sm_thresh);
#endif

  esl_fatal("bckfilter_dispatcher found no vector implementation - that shouldn't happen.");
}

#endif


/*****************************************************************
 * 3. Benchmark driver
 *****************************************************************/


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/


/*****************************************************************
 * 5. Test driver
 *****************************************************************/



/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef h4FBFILTER_EXAMPLE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_hmmfile.h"
#include "h4_profile.h"
#include "h4_mode.h"
#include "h4_refmx.h"
#include "h4_checkptmx.h"

#include "general.h"
#include "fbfilter.h"
#include "reference_dp.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version info",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 2, argc, argv, "example of using the Fwd/Bck filter", "[-options] <hmmfile> <seqfile>");
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = NULL;
  H4_REFMX       *fwd     = NULL;
  H4_CHECKPTMX   *cpx     = NULL;
  float           fsc, fsc_ref;
  int             status;

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
  h4_hmmfile_Close(hfp);

  if ( esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, /*env=*/NULL, &sqfp) != eslOK)
    esl_fatal("couldn't open sequence file %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  cpx = h4_checkptmx_Create(hmm->M, 400, ESL_MBYTES(32));
  fwd = h4_refmx_Create(hmm->M, 400);

  h4_mode_SetLocal(mo);  // FB filter is implicitly local. We set local mode so reference Fwd scores match FB.
  
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      h4_mode_SetLength(mo, sq->n);

      h4_checkptmx_SetDumpMode(cpx, stdout);

      h4_fwdfilter(sq->dsq, sq->n, hmm, mo, cpx, &fsc);
      h4_reference_Forward(sq->dsq, sq->n, hmm, mo, fwd, &fsc_ref);

      h4_refmx_Dump(stdout, fwd);

      printf("fwdfilter raw score (bits): %.2f\n", fsc);
      printf("reference Fwd score (bits): %.2f\n", fsc_ref);

      h4_refmx_Reuse(fwd);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);


  h4_checkptmx_Destroy(cpx);
  h4_refmx_Destroy(fwd);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif // h4FBFILTER_EXAMPLE
