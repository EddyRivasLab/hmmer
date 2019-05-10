/* Emitting (sampling) sequence from a profile HMM
 * 
 * Contents:
 *    1. Sequence sampling, h4_emit()
 *    2. Internal support functions
 *    x. Example
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"

static int sample_endpoints(ESL_RANDOMNESS *rng, const H4_PROFILE *hmm, int *ret_kstart, int *ret_kend);


/*****************************************************************
 * 1. Sequence sampling, h4_emit()
 *****************************************************************/

/* Function:  h4_emit()
 * Synopsis:  Sample a sequence from a profile
 * Incept:    SRE, Mon 06 May 2019 [Gurf Morlix, Series of Closing Doors]
 *
 * Purpose:   Sample a sequence from profile <hmm> in
 *            algorithm-dependent mode <mo>, using random number
 *            generator <rng>. Optionally return either the
 *            sampled sequence in <sq> or the sampled path in 
 *            <pi>, or both.
 *            
 *            <sq> and/or <pi>, if provided, are allocated by the
 *            caller. <sq> must be digital mode, and its alphabet must
 *            match the one in <hmm>. The function will <_Reuse()> and
 *            reinitialize them first, so caller doesn't have to call
 *            <_Reuse()> itself if it's reusing previously used
 *            structures.
 *            
 *            The minimum sampled sequence length is 1, because mute
 *            cycles are rejected.
 *            
 *
 * Args:      rng - random number generator
 *            hmm - profile to sample from
 *            mo  - algorithm-dependent profile configuration "mode"
 *            sq  - optRETURN: sampled digital sequence; or NULL
 *            pi  - optRETURN: sampled path
 *
 * Returns:   <eslOK> on success, and <sq> and <pi>, if provided, 
 *            contain the sampled sequence and path.
 *
 * Throws:    <eslEMEM> on allocation error. Now the states of <sq>
 *            <pi> are undefined, but they can be safely reused.
 */
int
h4_emit(ESL_RANDOMNESS *rng, const H4_PROFILE *hmm, const H4_MODE *mo, ESL_SQ *sq, H4_PATH *pi)
{
  int    k, kend;
  int8_t prv, st;
  float  xt[h4_NX][h4_NXT];
  int    i,x;
  int    nM;                  // number of M states used in B..E subpath. must be >0, else we reject the path (reject mute cycles)
  int    status;
  
  ESL_DASSERT1(( !sq || sq->abc == hmm->abc ));

  /* <mo> has the special transitions as bit scores. Backcalculate to
   * probabilities. 
   */
  for (i = 0; i < h4_NX; i++)
    for (x = 0; x < h4_NXT; x++)
      xt[i][x] = exp2f(mo->xsc[i][x]);
  
  do {
    if (pi) h4_path_Reuse(pi);
    if (sq) esl_sq_Reuse(sq);
    st = h4P_S;
    i  = 0;
    k  = 0;
    while (st != h4P_T)
      {
	/* Choose state transition. (including S|B|E|T that aren't explicit in path) */
	prv = st;
	switch (st) {
	case h4P_L:
	  if (( status = sample_endpoints(rng, hmm, &k, &kend)) != eslOK) goto ERROR;  // both k and kend have to be sampled, from implicit prob model's fragment distribution
	  st = h4P_ML;
	  break;
	case h4P_S:  st = h4P_N; break;
	case h4P_N:  st = (esl_rnd_FChoose(rng, xt[h4_N], h4_NXT)) == h4_MOVE ? h4P_B : h4P_N;                      break;
	case h4P_B:  st = (esl_rnd_FChoose(rng, xt[h4_B], h4_NXT)) == h4_MOVE ? h4P_G : h4P_L;       nM   = 0;      break;
	case h4P_G:  st = h4P_MG + esl_rnd_FChoose(rng, hmm->t[k], 3);                               kend = hmm->M; break;  // t[0][h4_TII] = 0 by convention, so this works
	case h4P_MG: st = ( k==kend ? h4P_E : h4P_MG + esl_rnd_FChoose(rng, hmm->t[k],          3));                break;  
	case h4P_IG: st = ( k==kend ? h4P_E : h4P_MG + esl_rnd_FChoose(rng, hmm->t[k] + h4_TIM, 3));                break;
	case h4P_DG: st = ( k==kend ? h4P_E : h4P_MG + esl_rnd_FChoose(rng, hmm->t[k] + h4_TDM, 3));                break;
	case h4P_ML: st = ( k==kend ? h4P_E : h4P_ML + esl_rnd_FChoose(rng, hmm->t[k],          3));                break;  // we're assuming h4P_ML-IL-DL order in path statetypes, and same for G types
	case h4P_IL: st = ( k==kend ? h4P_E : h4P_ML + esl_rnd_FChoose(rng, hmm->t[k] + h4_TIM, 3));                break;  // we're assuming h4_TMM(3) TIM(3) TDM(3) order in profile transition probs
	case h4P_DL: st = ( k==kend ? h4P_E : h4P_ML + esl_rnd_FChoose(rng, hmm->t[k] + h4_TDM, 3));                break;
	case h4P_E:  st = (esl_rnd_FChoose(rng, xt[h4_E], h4_NXT)) == h4_MOVE ? h4P_C : h4P_J;                      break;
	case h4P_J:  st = (esl_rnd_FChoose(rng, xt[h4_J], h4_NXT)) == h4_MOVE ? h4P_B : h4P_J;                      break;
	case h4P_C:  st = (esl_rnd_FChoose(rng, xt[h4_C], h4_NXT)) == h4_MOVE ? h4P_T : h4P_C;                      break;
	default: esl_fatal("inconceivable");
	}

	/* Based on the transition we just sampled, update k */
	if      (st == h4P_E)                  { k = 0; if (nM == 0) break; }  // nM=0 means a mute cycle; reject the path
	else if (st == h4P_ML && prv != h4P_L)   k++;                          // be careful about L->Mk, where we already set k 
	else if (st == h4P_MG)                   k++;
	else if (st == h4P_DG || st == h4P_DL)   k++;

	/* And based on the transition we just sampled, generate a residue. */
	if      (st == h4P_ML || st == h4P_MG)                           { x = esl_rnd_FChoose(rng, hmm->e[k], hmm->abc->K); nM++; i++; }
	else if (st == h4P_IL || st == h4P_IG)                           { x = esl_rnd_FChoose(rng, hmm->f,    hmm->abc->K);       i++; }
	else if ((st == h4P_N || st == h4P_C || st == h4P_J) && prv==st) { x = esl_rnd_FChoose(rng, hmm->f,    hmm->abc->K);       i++; }
	else    x = eslDSQ_SENTINEL;

	/* Add residue and state to growing sequence and path */
	if (sq && x != eslDSQ_SENTINEL)  { if ((status = esl_sq_XAddResidue(sq, x)) != eslOK) goto ERROR; }
	if (pi)                          { if ((status = h4_path_Append(pi, st))    != eslOK) goto ERROR; }
      }
  } while (nM == 0);  // outer loop is rejection sampling: reject any trace containing a mute cycle
    
  /* Terminate (optional) sq */
  if (sq && (status = esl_sq_XAddResidue(sq, eslDSQ_SENTINEL)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  return status;
}	       


/*****************************************************************
 * 2. Internal support functions
 *****************************************************************/

/* sample_endpoints()
 *
 * Purpose:   Given a profile <hmm> and random number source <rng>, sample
 *            a local begin transition from the _implicit_ probabilistic profile
 *            model, yielding a sampled start and end node; return these
 *            via <ret_kstart> and <ret_kend>.
 *            
 *            By construction, the entry at node <kstart> is into a
 *            match state, but the exit from node <kend> may turn
 *            out to be from either a match or delete state.
 *            
 *            The implicit probability model assumes that exits j are
 *            uniformly distributed conditional on a particular entry
 *            point i: $a_{ij} =$ constant $\forall j \geq i$.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            
 * Xref:      STL11/136-138           
 */
static int
sample_endpoints(ESL_RANDOMNESS *rng, const H4_PROFILE *hmm, int *ret_kstart, int *ret_kend)
{
  float *pstart = NULL;
  float  Z      = 0;    // normalization of L->Mk entry
  int    k;
  int    kstart, kend;

  int    status;

  /* We don't explicitly store the L->Mk probability distribution in a
   * model.  h4_profile.c::set_local_entry() sets LMk scores when we
   * convert probability model to profile scores. In emit.c, we want
   * to only depend on the probability model, for simplicity.  Of
   * various options, decided to just reproduce the set_local_entry()
   * calculation here in probability form. Can revisit if needed; this
   * is a little time consuming but we don't need to do it often. 
   * [SRE H6/142]
   */
  ESL_ALLOC(pstart, sizeof(float) * (hmm->M+1));
  h4_profile_CalculateOccupancy(hmm, pstart, NULL);                                       // this sets pstart[0] = 0; there's no M0 state.
  for (k = 1; k <= hmm->M; k++)  Z        +=  pstart[k]      * (float) (hmm->M + 1 - k);  // reproduce calculation in h4_profile::set_local_entry()
  for (k = 1; k <= hmm->M; k++)  pstart[k] = (pstart[k] / Z) * (float) (hmm->M + 1 - k);  //   (M+1-k) is weighting for how many ij fragments start at i                
  kstart = esl_rnd_FChoose(rng, pstart, hmm->M+1);                                        // sample the starting position from that distribution 
  kend   = kstart + esl_rnd_Roll(rng, hmm->M-kstart+1);                                   // and the exit uniformly from possible exits for it

  free(pstart);
  *ret_kstart = kstart;
  *ret_kend   = kend;
  return eslOK;
  
 ERROR:
  if (pstart) free(pstart);
  *ret_kstart = 0;
  *ret_kend   = 0;
  return status;
}


    
/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef h4EMIT_EXAMPLE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"

int
main(int argc, char **argv)
{
  char           *hmmfile = argv[1];
  ESL_RANDOMNESS *rng     = esl_randomness_Create(0);
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create(); // creating mode also puts it in default configuration except for setting length model.
  H4_PATH        *pi      = h4_path_Create();
  ESL_ALPHABET   *abc     = NULL;
  ESL_SQ         *sq      = NULL;
  int             L       = 400;
  int             nhmm    = 0;
  int             status;

  if (( status = h4_hmmfile_Open(hmmfile, NULL, &hfp)) != eslOK)
    esl_fatal("profile HMM input failed: %s", hfp->errmsg);

  h4_mode_SetLength(mo, L);

  while (( status = h4_hmmfile_Read(hfp, &abc, &hmm) ) == eslOK)
    {
      if (! sq) sq = esl_sq_CreateDigital(abc);

      h4_emit(rng, hmm, mo, sq, pi);

      esl_sq_SetName(sq, "sample");
      esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, FALSE);

      h4_path_Dump(stdout, pi);

      esl_sq_Reuse(sq);
      h4_path_Reuse(pi);
      h4_profile_Destroy(hmm);
      nhmm++;
    }
  if       (status == eslEFORMAT)   esl_fatal("Parse failed, bad profile HMM file format in %s:\n   %s", argv[1], hfp->errmsg);
  else if  (nhmm == 0)              esl_fatal("Empty input? no profiles read from %s", argv[1]);
  else if  (status != eslEOF )      esl_fatal("Unexpected error in profile HMM reader");

  esl_alphabet_Destroy(abc);
  esl_sq_Destroy(sq);
  h4_path_Destroy(pi);
  h4_mode_Destroy(mo);
  h4_hmmfile_Close(hfp);
  esl_randomness_Destroy(rng);
  return 0;
}
#endif // h4EMIT_EXAMPLE
