/* Reference implementation of envelope definition.
 * 
 * Envelope definition step comes after we define the anchor set, and
 * before we do the alignment.
 * 
 * Contents:
 *    1. Envelope definition API call.
 *    2. Internal functions: the steps of envelope definition.
 *    3. Copyright and license information.
 */
// TODO:
//       write an example driver
//       get reference anchor set definition coded as designed

#include "p7_config.h"

#include "easel.h"

#include "base/p7_envelopes.h" 
#include "base/p7_coords2.h"
#include "base/p7_profile.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_asc_fwdback.h"


static int envcoords(P7_ENVELOPE *env, int D, const P7_REFMX *apd, float threshold);
static int glocality(P7_ENVELOPE *env, int D, const P7_REFMX *apd);
static int outcoords(P7_ENVELOPE *env, int D, const P7_REFMX *apd, float epsilon);
static int approxsc (P7_ENVELOPE *env, int D, const P7_REFMX *afd, const P7_PROFILE *gm);


int
p7_reference_Envelopes(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_COORD2 *anch, int D,
		       const P7_REFMX *apu, const P7_REFMX *apd,
		       P7_REFMX *afu, P7_REFMX *afd, P7_ENVELOPE *env)
{
  int d;
  int status;

  for (d = 0; d < D; d++)
    {
      env[d].i0    = anch[d].n1;
      env[d].k0    = anch[d].n2;
      env[d].flags = 0;
    }
  
  if ((status = envcoords(env, D, apd, 0.5))   != eslOK) goto ERROR;
  if ((status = glocality(env, D, apd))        != eslOK) goto ERROR;
  if ((status = outcoords(env, D, apd, 0.005)) != eslOK) goto ERROR;

  if (D > 0)
    {
      if ((status = approxsc(env, D, afd, gm)) != eslOK) goto ERROR;

      p7_refmx_Reuse(afu);
      p7_refmx_Reuse(afd);
      
      for (d = 0; d < D; d++)
	if (! (env[d].flags & p7E_ENVSC_APPROX))
	  {
	    status = p7_ReferenceASCForward( dsq+env[d].oea-1, env[d].oeb-env[d].oea+1, gm,
					     &(anch[d]), 1, afu, afd, &(env[d].env_sc));
	    if (status != eslOK) goto ERROR;

	    p7_refmx_Reuse(afu);
	    p7_refmx_Reuse(afd);
	  }
    }
  else 
    env[0].env_sc = P7R_XMX(afd, L, p7R_C) + gm->xsc[p7P_C][p7P_MOVE];

  /* ka,kb are not set yet. alignment does that part...
   * even though we already know for glocal alis that ka=1,kb=M.
   */
  return eslOK;
  
 ERROR:
  return status;

}

/* all specials are in <apd> */
/* threshold = 0.5 */
static int
envcoords(P7_ENVELOPE *env, int D, const P7_REFMX *apd, float threshold)
{
  float  phomology;
  int    d, i, iz;
  int    L = apd->L;
  
  for (d = 0; d < D; d++)
    {
      /* Determine ia, envelope start position;
       *  min_{i <= i0} P(x_i in domain d) >= threshold
       */
      phomology = 1.0;
      iz        = (d == 0 ? 0 : env[d-1].i0);
      for (i = env[d].i0 - 1; i >= iz; i--)
	{
	  phomology -= P7R_XMX(apd, i, p7R_B);   // now phomology = P(res i in domain d)
	  if (phomology < threshold) break;      // if i is not in domain...
	}
      env[d].ia = i+1;                           // but i+1 was, so ia = i+1.

      /* Determine ib, envelope end position;
       *   max_{i >= i0} P(x_i in domain d) >= threshold
       */
      phomology = 1.0;
      iz        = (d == D-1 ? L+1 : env[d+1].i0);
      for (i = env[d].i0 + 1; i < iz; i++)
	{
	  phomology -= P7R_XMX(apd, i, p7R_E);    // now phomology = P(x_i in dom_d)
	  if (phomology < threshold) break;       // if i is not in domain...
	}
      env[d].ib = i-1;		                  // but i-1 was, so ib = i-1.
    }

  return eslOK;
}



/* 
 *  For each domain d,
 *    pL = \sum_i=i0(d-1)..i0(d)-1 P(i,L)
 *    pG = \sum_i=i0(d-1)..i0(d)-1 P(i,G)
 *  and if pG >= pL, it's a glocal alignment.
 *  
 *  We do need the i0(d)-1 row, though it may not look like it at
 *  first. From a B on i0(d)-1 row, we need to use B->G->DDDDk0-1 path
 *  to be able to reach the anchor Mk on row i0(d). In decoding, we
 *  wing unfold the G->DD->Mk paths, so this path exists in a decoding
 *  matrix. Only for G, though; P(L, i0(d)) is 0.
 * 
 */
static int
glocality(P7_ENVELOPE *env, int D, const P7_REFMX *apd)
{
  float pL = 0.0;
  float pG = 0.0;
  int   d, i, iz;

  for (d = 0; d < D; d++)
    {
      iz = (d == 0 ? 0 : env[d-1].i0);
      for (i = iz; i < env[d].i0; i++)
	{
	  pL += P7R_XMX(apd, i, p7R_L);
	  pG += P7R_XMX(apd, i, p7R_G);
	}
      if (pG >= pL) env[d].flags |= p7E_IS_GLOCAL;
    }
  return eslOK;
}


/* "outer envelope", oea..oeb, for env score purposes */
/* for epsilon = 0.005 (allowing 0.005 of the probability
 * mass to be in paths that extend beyond the envelope on 
 * each side), we lose a max of 0.01 of the mass, which
 * is a log-odds score of ~0.01 in nats; negligible.
 */
static int
outcoords(P7_ENVELOPE *env, int D, const P7_REFMX *apd, float epsilon)
{
  int d, i, iz;
  float phomology;
  int   s;
  int   L = apd->L;
  int   leftchoked; 
  int   rightchoked;

  for (d = 0; d < D; d++)
    {
      leftchoked = rightchoked = FALSE;

      phomology = 1.0;
      iz = (d == 0 ? 0 : env[d-1].i0);
      s  = (d == 0 ? p7R_N : p7R_J);
      for (i = env[d].i0 - 1; i >= iz; i--)
	{
	  if (P7R_XMX(apd, i, s) >= 1.0 - epsilon) { leftchoked = TRUE; break; }
	  phomology -= P7R_XMX(apd, i, p7R_B);
	  if (phomology < epsilon) break;
	}
      env[d].oea = i+1;

      phomology = 1.0;
      iz = (d == D-1 ? L : env[d+1].i0-1);
      s  = (d == D-1 ? p7R_C : p7R_J);
      for (i = env[d].i0; i < iz; i++)
	{
	  if (P7R_XMX(apd, i, s) >= 1.0 - epsilon) { rightchoked = TRUE; break; }
	  phomology -= P7R_XMX(apd, i+1, p7R_E); // now phomology = P(x_i+1 in dom_d)
	  if (phomology < epsilon) break;        // if i+1 is not in domain...
	}
      env[d].oeb = i;		                 // but i=1 was, so oeb = i.

      if (leftchoked && rightchoked) 
	env[d].flags |= p7E_ENVSC_APPROX;
    }
  return eslOK;
}



/* If t_NN = t_JJ = t_CC, then we can often calculate the envelope
 * score of a single domain by an accurate approximation from the ASC
 * Forward matrix for the whole sequence.
 * 
 * Assume that we know that Xi and Yj are "choke points" flanking our
 * domain, such that >= 1.0-epsilon of the posterior path probability
 * flows through each of them.
 * 
 * For d=0, X=N, Y=J; for d=D-1, X=J, Y=C; for internal d, X=Y=J.
 * 
 * Then F(j,X) - F(i,X) = the sum of all path suffixes that start in
 * Xi and end at Xj.
 * 
 * Then, because t_NN = t_JJ = tCC, we just need to add in the log
 * transitions to get us from S to Xi via N's, and get us from Xj to T
 * via C's to get the envelope score. That's t_SN + i * t_NN, and
 * (L-j) * t_CC + t_CT. (And t_SN = 0).
 */
static int
approxsc(P7_ENVELOPE *env, int D, const P7_REFMX *afd, const P7_PROFILE *gm) 
{
  int d;
  int L = afd->L;
  int s1, s2;


  for (d = 0; d < D; d++)
    if (env[d].flags & p7E_ENVSC_APPROX)
      {
	s1 = (d == 0   ? p7R_N : p7R_J);
	s2 = (d == D-1 ? p7R_C : p7R_J);

	env[d].env_sc = 
	  P7R_XMX(afd, env[d].oeb,   s2) -
	  P7R_XMX(afd, env[d].oea-1, s1);

	env[d].env_sc += 
	  gm->xsc[p7P_N][p7P_LOOP] * (env[d].oea-1) +
	  gm->xsc[p7P_C][p7P_LOOP] * (L-env[d].oeb) +
	  gm->xsc[p7P_C][p7P_MOVE];
      }
  return eslOK;
}


/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7REFERENCE_ENVELOPES_EXAMPLE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of envelope determination step";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_COORDS2     *anch    = p7_coords2_Create(0,0);
  P7_REFMX       *rxf     = NULL;
  P7_REFMX       *rxd     = NULL;
  P7_REFMX       *afu     = NULL;
  P7_REFMX       *afd     = NULL;
  P7_REFMX       *apu     = NULL;
  P7_REFMX       *apd     = NULL;
  P7_TRACE       *tr      = NULL;
  float          *wrk     = NULL;
  P7_ENVELOPES   *env     = p7_envelopes_Create(0,0);
  P7_COORDS2_HASH *hashtbl = p7_coords2_hash_Create(0,0,0);
  P7_MPAS_PARAMS  prm;
  P7_MPAS_STATS   stats;
  float           fsc, vsc, asc, asc_b;
  int             status;

  p7_Init();

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Read one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);
  esl_sqfile_Close(sqfp);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);

  /* Set the profile and null model's target length models */
  p7_bg_SetLength     (bg, sq->n);
  p7_profile_SetLength(gm, sq->n);

  /* Allocate DP matrices and tracebacks */
  rxf = p7_refmx_Create(gm->M, sq->n);
  rxd = p7_refmx_Create(gm->M, sq->n);
  tr  = p7_trace_Create();
  afu = p7_refmx_Create(gm->M, sq->n);
  afd = p7_refmx_Create(gm->M, sq->n);

  /* First pass analysis */
  p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf,  tr, &vsc);
  p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc);
  p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd, NULL);   
  p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd);   

  /* Customize MPAS parameters if you want; these are the defaults. */
  prm.max_iterations = 1000;
  prm.loss_threshold = 0.001;
  prm.be_verbose     = TRUE;

  /* MPAS algorithm gets us an anchor set */
  p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, tr, &wrk, hashtbl,
		       afu, afd, anch, &asc, &prm, &stats);

  
  printf("# ASC Forward UP:\n");    p7_refmx_Dump(stdout, afu);
  printf("# ASC Forward DOWN:\n"); p7_refmx_Dump(stdout, afd);

  /* We no longer need rxf and rxd. 
   * Use their space for apu/apd pair, which will briefly
   * hold ASC Backward matrices, then get used for ASC Decoding.
   */
  apu = rxf; p7_refmx_Reuse(apu);
  apd = rxd; p7_refmx_Reuse(apd);

  p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch->arr, anch->n, apu, apd, &asc_b);
  
  printf("# Backward score (raw, nats): %.2f\n", asc_b);
  printf("# ASC Backward UP:\n");   p7_refmx_Dump(stdout, apu);
  printf("# ASC Backward DOWN:\n"); p7_refmx_Dump(stdout, apd);

  /* ASC Decoding takes afu/afd and abu/abd as input;
   * overwrites abu/abd with decoding matrices
   */
  p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch->arr, anch->n, afu, afd, apu, apd, apu, apd);

  printf("# ASC Decoding UP matrix:\n");  p7_refmx_Dump(stdout, apu);
  printf("# ASC Decoding DOWN:\n");       p7_refmx_Dump(stdout, apu);


  /* Envelope calculation needs to get four matrices:
   * ASC Decoding pair, apu/apd, and it will leave these constant;
   * ASC Forward pair,  afu/afd, and it will overwrite these.
   */
  p7_envelopes_GrowTo(env, anch->n);
  p7_reference_Envelopes(sq->dsq, sq->n, gm, anch->arr, anch->n, apu, apd, afu, afd, env->arr);


  p7_envelopes_Destroy(env);
  p7_coords2_hash_Destroy(hashtbl);
  if (wrk) free(wrk);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(rxd);
  p7_refmx_Destroy(rxf);
  p7_coords2_Destroy(anch);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_ENVELOPES_EXAMPLE*/





/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
