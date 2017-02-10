/* Null2 biased-composition score adjustment calculations.
 * Sparse implementation; production code.
 *
 * Contents:
 *   1. Null2 calculation routine.
 *   2. Example
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "base/p7_envelopes.h"
#include "base/p7_profile.h"

#include "misc/logsum.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/sparse_null2.h"


/*****************************************************************
 * 1. Null2 calculation routine
 *****************************************************************/

/* Function:  p7_sparse_Null2()
 * Synopsis:  Determine domain null2 scores, for composition bias correction.
 *
 * Purpose:   Calculate the ad hoc domain null2 scores for composition
 *            bias correction.  We have an anchor-set-constrained
 *            (ASC) posterior decoding matrix <apd>, for a comparison
 *            of profile <gm> to digital sequence <dsq> of length <L>,
 *            using anchors that are now in the envelope data <env>.
 *            Also in <env>, we have determined and stored the
 *            <ia..ib> and <oa..ob> envelope coordinates, and the
 *            envelope scores <env_sc>.
 *            
 *            Here we parameterize the ad hoc null2 model for each
 *            domain <d>; use it to calculate each domain null2 score;
 *            and store that in the envelope data <env>.
 *            
 *            We need two pieces of working memory. <wrk> is space for
 *            at least <M+1> floats.  <null2> is space for at least
 *            <Kp> floats for the alphabet size <gm->abc->Kp> including
 *            noncanonical symbols and degeneracies.
 *            
 *            <wrk[1..M]> gets used to collect expected usage counts
 *            $u_d(y)$. <null2[0..K-1;K..Kp-1]> gets used to hold the
 *            null2 model $W$'s log-odds score parameters.
 *            
 *            Upon successful return, <env->arr[d].null2_sc> fields
 *            have been set for each domain <d>. The only data in
 *            <env> that still needs to be set are the <ka..kb> model
 *            coords, which get set when we infer an alignment.
 *            
 *            <wrk[1..M]> contains the log usage frequency of each
 *            match state (marginalized over L/G), and <null2[]>
 *            contains the log-odds score parameters -- for the final
 *            domain <D>. It's hard to imagine this being of any use
 *            to the caller in production code, since it's only for
 *            the final domain. However, it's free information that
 *            might be useful to know in debugging/testing. This 
 *            workspace is handled by Easel's "bypass" idiom: usually
 *            you will initialize <wrk=NULL> and pass <&wrk> for
 *            the first and all successive calls, and 
 *            it will be reused while (re)allocating as needed.
 *            
 *            This routine uses <p7_FLogsum()> calculations, so caller
 *            must have initialized it with <p7_FLogsumInit()>.
 *            
 * Args:      dsq      : digital sequence 1..L
 *            L        : length of <dsq>
 *            gm       : search profile
 *            apd      : ASC posterior decoding matrix
 *            env      : envelope data for D domains, including i0/k0, ia..ib, and env_sc.
 *            wrk_byp  : ptr to working space for at least gm->M+1 floats; or to NULL and it'll be allocated (easel's "bypass" idiom)
 *            null2    : working space for at least gm->abc->Kp floats
 *
 * Returns:   <eslOK> on success and <env->arr[d].null2_sc> is set for each domain <d>.
 *
 * Throws:    <eslEMEM> on (re)allocation failure.
 */
int
p7_sparse_Null2(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *apd, P7_ENVELOPES *env, float **wrk_byp, float *null2)
{
  float               *wrk = NULL;       // M+1 floats of working space
  const P7_SPARSEMASK *sm  = apd->sm;    // sparse mask for <apd>
  const float         *ppp = apd->dp;    // <ppp> will step through sparse cells in <apd>
  int                  M   = gm->M;    
  float const         *rsc; 
  float                usum;             // Sum of expected counts, \sum_y u(y), for all but I/N/C/J.
  float                uincj;            // Deduced usage of I/N/C/J emissions, by Ld-usum
  float                Ld;               // Length of domain, as a float
  int                  need_degens;      // Flag gets raised if we need noncanonical residue scores in null2[] 
  float                base;             // We cap domain null2 score r_d to a max of s_d - base; above this, the domain call basically shouldn't be made at all
  int                  a,d,g,i,k,z;      // Indices: residues in alphabet, domains, segments, position in dsq, match states in model, sparse cells on row
  float                rd;		 // domain null2 score: this is what we're calculating for each domain <d>
  int                  status;
  
  /* Contract checks */
  ESL_DASSERT1(( apd->type == p7S_ASC_DECODE ));
  ESL_DASSERT1(( sm->L == L && env->L == L   ));
  ESL_DASSERT1(( sm->M == M && env->M == M   ));

  /* bypass idiom for the workspace */
  if      (esl_byp_IsInternal(wrk_byp)) { ESL_ALLOC  (wrk,      (gm->M + 1) * sizeof(float));                 }
  else if (esl_byp_IsReturned(wrk_byp)) { ESL_ALLOC  (*wrk_byp, (gm->M + 1) * sizeof(float)); wrk = *wrk_byp; }
  else if (esl_byp_IsProvided(wrk_byp)) { ESL_REALLOC(*wrk_byp, (gm->M + 1) * sizeof(float)); wrk = *wrk_byp; }

  /* One time initialization of noncanonical residue scores that we shouldn't be using for <dsq> anyway: */
  null2[gm->abc->K]    = 0.0;    /* gap character    */
  null2[gm->abc->Kp-2] = 0.0;	 /* nonresidue "*"   */
  null2[gm->abc->Kp-1] = 0.0;	 /* missing data "~" */  

  /* We have a cap = s_d - base calculation later, capping the domain null2 correction r_d to a maximum. */
  base = ( L ? (float) L * gm->xsc[p7P_C][p7P_LOOP] : 0.0f) + gm->xsc[p7P_C][p7P_MOVE];
  //        ^^^  must avoid 0 * -inf = NaN  ^^^

  /* here we go... a pass thru all of <apd>'s sparse cells, stopping
   * to count expected usage counts of M*k states in ia..ib cells for
   * each domain <d>.
   */
  g = 0;
  for (d = 1; d <= env->D; d++)
    {
      /* Because of anchor+envelope constrained sparsity, the mechanics of
       * simply finding relevant matrix cells ppp[i,k] are heavy. Don't
       * despair; all the following chunk of code is doing is counting
       * expected usage:
       *    u(d,k) = \sum_i=ia(d)..ib(d) P(\lambda_i = M*k)
       * for each emitting state M*k, for the current domain d.
       * These u_k counts are stored in wrk[1..M].
       * 
       * We don't need to count I/N/C/J expected emissions because all of
       * these are hardcoded to emit with 0 (background) score. All we
       * need is their total number. And since the total expected # of
       * counts is just the domain's length Ld, we can deduce that sum
       * just by Ld - \sum_k u(k).
       * 
       * At one time we did have inserts emitting with nonzero scores.
       * To ease any future return to such things, I've left stubs in the
       * code below that make it obvious where you'd do the insert counting.
       *              
       * We do this all in one block of code rather than breaking things
       * up into functions (count_expectations(), whatever) because of the
       * traversal mechanism over the sparse matrix. If we used functions,
       * we'd have to pass the *ppp variable around to keep state on where
       * we are in the traversal.
       */
      esl_vec_FSet(wrk, M+1, 0.0f);

      /* <ppp> maintenance (1): skip leading rows <i>, get <ppp> onto start of row <ia(d)> */
      while (env->arr[d].i0 > sm->seg[g].ib)  g++;                          // find which segment <g> that anchor <d> is in
      if (env->arr[d-1].i0 < sm->seg[g].ia)                                 // If <d> is 1st domain in seg <g>
	{                                                                   // then for ASC rows ia(g)..ia(d)-1, ASC is only in UP sector;
	  for (i = sm->seg[g].ia; i < env->arr[d].ia; i++)                  // advance <ppp> over UP(i) rows.
	    for (z = 0; z < sm->n[i] && sm->k[i][z] < env->arr[d].k0; z++)
	      ppp += p7S_NSCELLS;
	}
      else                                                                  // Else <d> isn't first,
	{                                                                   // and for ASC rows ib(d-1)..ia(d)-1, 
	  for (i = env->arr[d-1].ib+1; i < env->arr[d].ia; i++)             // advance <ppp> over both DOWN(i) and UP(i) supercells.
	    {
	      for (z = sm->n[i]-1; z >= 0       && sm->k[i][z] >= env->arr[d-1].k0; z--)  // Order of access doesn't matter here; only # of sparse supercells in DOWN(i)
		ppp += p7S_NSCELLS;  
	      for (z = 0;          z < sm->n[i] && sm->k[i][z] <  env->arr[d].k0;   z++)  // UP row includes any supercells 1..k0(d)-1
		ppp += p7S_NSCELLS; 
	    }
	}
      
      /* UP sector ia(d)..i0(d)-1; 1..k0(d)-1 */
      for (i = env->arr[d].ia; i < env->arr[d].i0; i++)
	{
	  /* <ppp> maintenance (2): skip leading DOWN supercells when i in ASC DOWN+UP row */
	  if (env->arr[d-1].i0 >= sm->seg[g].ia)
	    for (z = sm->n[i]-1; z >= 0 && sm->k[i][z] >= env->arr[d-1].k0; z--)
	      ppp += p7S_NSCELLS;
	  
	  for (z = 0; z < sm->n[i] && sm->k[i][z] < env->arr[d].k0; z++, ppp += p7S_NSCELLS)
	    {
	      // wrk[0]        += ppp[p7S_IL] + ppp[p7S_IG];  // If you had to count inserts, here's place 1 of 2 where you'd do it.
	      wrk[sm->k[i][z]] += ppp[p7S_ML] + ppp[p7S_MG];
	    }
	}

      /* DOWN sector i0(d)..ib(d); k0(d)..M */
      for (i = env->arr[d].i0; i <= env->arr[d].ib; i++)
	{
	  z = 0; while (z < sm->n[i] && sm->k[i][z] < env->arr[d].k0) z++;
	  for (; z < sm->n[i]; z++, ppp += p7S_NSCELLS)
	    {
	      // wrk[0]        += ppp[p7S_IL] + ppp[p7S_IG];  // If you had to count inserts, here's place 2 of 2 where you'd do it.
	      wrk[sm->k[i][z]] += ppp[p7S_ML] + ppp[p7S_MG];
	    }

	  /* <ppp> maintenance (3); skip trailing UP supercells when i is ASC DOWN/UP row */
	  if (env->arr[d+1].i0 <= sm->seg[g].ib)
	    for (z = 0; z < sm->n[i] && sm->k[i][z] < env->arr[d+1].k0; z++)   
	      ppp += p7S_NSCELLS;
	}
	
      /* <ppp> maintenance (4): if <d> is last anchor in seg, 
       * skip trailing rows ib(d)+1..ib(g), which must be
       * DOWN only rows.
       */
      if (env->arr[d+1].i0 > sm->seg[g].ib)
	for (i = env->arr[d].ib+1; i <= sm->seg[g].ib; i++)
	  for (z = sm->n[i]-1; z >= 0 && sm->k[i][z] >= env->arr[d].k0; z--)  // DOWN is k0(d)..M; order doesn't matter, so we can skip it backwards.
	    ppp += p7S_NSCELLS;

      /* Ta-da! Now the heavyweight sparse AEC matrix traversal is complete.
       * And now wrk[] contains expected usage counts per emitting state: u_d(y) in our notation.
       * ... but not including I/N/C/J emissions, which we'll deduce by difference from domain length Ld.
       */


      /* So let's do that. Determine I/N/C/J expected counts by difference.
       * Then convert all counts to log frequencies.
       */
      Ld    = (float) (env->arr[d].ib - env->arr[d].ia + 1);   // Length of domain d in residues
      usum  = esl_vec_FSum(wrk+1, M);                          // Sum of expected counts over all emitting Mk=1..M except I/N/C/J
      uincj = Ld - usum;                                       // Since the expected counts must add up to Ld, we can deduce N/C/J counts
      esl_vec_FScale(wrk+1, M, 1./Ld);                         // Now wrk[1..M] are normalized occurrence frequencies. Ld>0, so there's no div by zero.
      esl_vec_FLog  (wrk+1, M);                                // And now wrk[1..M] are log frequencies. -inf are possible, for wrk[k] = 0.
      uincj = log (uincj / Ld);                                // Now uincj is a log frequency too. It's possible to get -inf here.

      /* Now we're going to get to log-odds null2 residue scores log w(a)/f(a),
       * from our \sigma(k,a) = log [ e(k,a) / f(a) ] match emission scores.
       * Given *normalized frequencies* u(k), uincj, what we're aiming to do
       * looks like this:
       * 
       *     w(a) / f(a)  =    \sum_y u(y) [ e(y,a) / f(a) ]
       * log(w(a) / f(a)) = \logsum_y log(u(y)) + log [ e(s,a) / f(a) ]
       *                  = \logsum{ log(uincj) + 0), \logsum_k [ log(u(k)) + \sigma(k, a) ] 
       *
       * and we've already calculated log(uincj), log(u(k)) above, of course.
       */
      for (a = 0; a < gm->abc->K; a++)
	{
	  rsc      = gm->rsc[a] + p7P_NR + p7P_M; 
	  null2[a] = uincj;
	  for (k = 1; k <= M; k++)
	    {
	      null2[a] = p7_FLogsum(null2[a], wrk[k] + *rsc);   // -inf + -inf = -inf is possible here, and that's fine.
	      rsc += p7P_NR;
	    }
	}
      /* OK: now null2[a] = \log[ w(a) / f(a)] */


      /* If we need to, make valid scores for all degeneracies, by averaging the log odds ratios. 
       * Note H3 averaged odds ratios; H4 averages lod scores. 
       * Small optimization: see if x(ia..ib) has any noncanonical residues first.
       */
      need_degens = FALSE;
      for (i = env->arr[d].ia; i <= env->arr[d].ib; i++)
	if (! esl_abc_XIsCanonical(gm->abc, dsq[i])) { need_degens = TRUE; break; }
      if (need_degens)
	esl_abc_FAvgScVec(gm->abc, null2); // does not set gap, nonres, missing data; but we already initialized those to 0


      /* and at long last: calculate the domain null2 score for domain <d> */
      rd = 0.;
      for (i = env->arr[d].ia; i <= env->arr[d].ib; i++)
	rd += null2[dsq[i]];
      env->arr[d].null2_sc = ESL_MIN( rd, env->arr[d].env_sc - base );
    } // end loop over domains <d>

  if (esl_byp_IsInternal(wrk_byp)) { free(wrk); }
  return eslOK;

 ERROR:
  if (esl_byp_IsInternal(wrk_byp) && wrk) { free(wrk); }
  return status;
}


/* SRE: THE CODE BELOW IS LEGACY: WHEN P7_PIPELINE IS REWRITTEN, BURN IT */

/* Function:  p7_sparse_Null2ByExpectation()
 * Synopsis:  Calculate null2 lod scores, from sparse decoding matrix.
 *
 * Purpose:   Calculate $\log \frac{f'(x)}{f(x)}$, null2 lod scores for
 *            a given envelope, given a sparse posterior decoding matrix <sxd>, the 
 *            profile <gm>, and the envelope coords <iae,ibe,kae,kbe>.
 *            
 *            Requires a caller-allocated workspace <wrk> with room for at least
 *            <gm->M+1> floats.
 *            
 *            Stores the resulting lod scores in caller-allocated space <null2>,
 *            with room for at least <gm->abc->Kp> scores. 
 *            
 *            $f'(x)$ (the null2 emission probability for x) is defined (for <x=0..K-1>)
 *            the weighted average of the emission probabilities of profile states
 *            used to generate the envelope, weighted by the expected occurrence 
 *            frequency of those states (i.e. by using their posterior probabilities). 
 *            
 *            As for degenerate residues <x=K..Kp-1>: for residues (K+1..Kp-3), 
 *            f'(x) = the mean lod score over canonical residues included in the
 *            degeneracy. For not-residues (gap K, nonresidue *, missing data ~),
 *            set $\log \frac{f'(x)}{f(x)} = 0$, implying $f'(x)=f(x)$.
 *            
 *            This averaging behavior differs from H3.0, which averaged the odds
 *            ratios, not the log odds ratios; the difference is because of a switch
 *            to keeping null2[] in log space and saving some exponentiation calls.
 *            I don't believe that noncanonical residues have enough effect on the
 *            ad hoc null2[] calculation to make it worth worrying about the effect 
 *            of this change.
 *            
 * Args:      gm       - profile
 *            sxd      - sparse posterior probability matrix, already calculated
 *            iae,ibe  - envelope coords on sequence: 1..L
 *            kae,kbe  - envelope coords on profile: 1..M
 *            wrk      - caller provided tmp space for M+1 floats
 *            null2    - RETURN: caller provided space for null2[0..Kp-1] lod scores
 *            
 * Returns:   <eslOK> on success, and <null2> contains lod scores for a null2
 *            composition model.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_sparse_Null2ByExpectation(const P7_PROFILE *gm, const P7_SPARSEMX *sxd, 
			     int iae, int ibe, int kae, int kbe,
			     float *wrk, float *null2)
{
  const P7_SPARSEMASK *sm = sxd->sm;
  float *xc  = sxd->xmx;
  float *dpc = sxd->dp;
  float const *rsc;
  int    i,k,x,z;

  /* Skip ahead in <sxd>, to get <dpc> and <xpc> on start of sparse row <iae>
   * (actually, the next row that's stored >= iae)
   */
  for (i = 1; i < iae; i++)
    if (sm->n[i]) {
      if (sm->n[i-1] == 0) xc += p7R_NSCELLS;
      dpc += sm->n[i] * p7R_NSCELLS;
      xc  += p7R_NXCELLS;
    }

  /* wrk[k] is the expected # of uses of Mk (ML|MG). 
   * wrk[kae-1] is the expected # of uses of Ik (IL|IG) and NCJ emission.
   * Collect it, by summation over sparse posterior prob matrix <sxd>
   */
  esl_vec_FSet(wrk+kae-1, kbe-kae+2, 0.0f);
  for (i = iae; i <= ibe; i++)
    {
      if (sm->n[i])
	{
	  if (sm->n[i-1] == 0) xc += p7R_NSCELLS;

	  z = 0;
	  while (z < sm->n[i] && sm->k[i][z] <  kae) z++; 
	  for (; z < sm->n[i] && sm->k[i][z] <= kbe; z++) 
	    {
	      wrk[sm->k[i][z]] += dpc[p7S_ML] + dpc[p7S_MG];
	      wrk[kae-1]       += dpc[p7S_IL] + dpc[p7S_IG]; 
	      dpc += p7S_NSCELLS;
	    }
	  wrk[kae-1] += xc[p7S_N] + xc[p7S_JJ] + xc[p7S_CC];
	  xc         += p7S_NXCELLS;
	}
      else wrk[kae-1] += 1.;
    }

  /* Normalize those expected usage #'s to frequencies.
   * We use those as weights, to calculate null2[x] as 
   * a weighted average of the states used to generate the
   * iae..ibe envelope.
   * If kae..kbe=1..M, then \sum_{k=0}^{M} = Ld = (ibe-iae+1);
   * but because there might be posterior probability outside
   * the kae..kbe envelope, \sum_{k=0}^{M} <= Ld
   */
  esl_vec_FNorm(wrk+kae-1, kbe-kae+2);
  esl_vec_FLog (wrk+kae-1, kbe-kae+2);
  
  /* Collect null2's emission odds: 
   *  null2[x] = \sum_k wrk[k] [e(k,x) / f(x)]
   * but we do it in log space, because that's what the profile's emission scores are:
   *  log null2[x] = \logsum_k log wrk[k] + log [e(k,x)/f(x)]
   */
  for (x = 0; x < gm->abc->K; x++)
    {
      null2[x] = wrk[kae-1];	               /* wrk[0] * 1.0 : emission odds ratio for N,C,J,I emissions is 1.0 */
      rsc = gm->rsc[x] + p7P_NR * kae + p7P_M; /* initialize <rsc> ptr on MSC(kae) */
      for (k = kae; k <= kbe; k++)
	{
	  null2[x] = p7_FLogsum(null2[x], wrk[k] + *rsc); /* this is wrk[k] + MSC(k) for residue x: stepping thru rsc for efficiency */
	  rsc     += p7P_NR;
	}
    }
  /* Now null2[x] = \log \frac{f'(x)}{f(x)} for all canonical x=0..K-1 */
  
  /* make valid scores for all degeneracies, by averaging the odds ratios. */
  esl_abc_FAvgScVec(gm->abc, null2); /* does not set gap, nonres, missing  */ // note 3.0 averaged odds ratios; 3.1 averages lod scores
  null2[gm->abc->K]    = 0.0;        /* gap character    */
  null2[gm->abc->Kp-2] = 0.0;	     /* nonresidue "*"   */
  null2[gm->abc->Kp-1] = 0.0;	     /* missing data "~" */
  return eslOK;
}


/*****************************************************************
 * 2. Example
 *****************************************************************/
#ifdef p7SPARSE_NULL2_EXAMPLE

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
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "include all cells in sparse mx",                   0 },
  { "-s",         eslARG_INT,    "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
 
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of sparse envelope definition";

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
  P7_OPROFILE    *om      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_CHECKPTMX   *cx      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  P7_SPARSEMX    *sxf     = NULL;
  P7_SPARSEMX    *sxd     = NULL;
  P7_TRACE       *tr      = NULL;
  P7_SPARSEMX    *asf     = NULL;
  P7_SPARSEMX    *asb     = NULL;
  P7_SPARSEMX    *asd     = NULL;
  P7_ANCHORS     *anch    = p7_anchors_Create();
  P7_ANCHORS     *vanch   = p7_anchors_Create();
  P7_ANCHORHASH  *ah      = p7_anchorhash_Create();
  P7_ENVELOPES   *env     = p7_envelopes_Create();
  float          *wrk     = NULL;
  float          *wrk2    = NULL;
  float          *null2   = NULL;
  float           fsc, vsc, asc_f, asc_b;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
 
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Get one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create (hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config   (gm, hmm, bg);
  p7_oprofile_Convert (gm, om);

  p7_bg_SetLength     (bg, sq->n);
  p7_profile_SetLength(gm, sq->n);
  p7_oprofile_ReconfigLength(om, sq->n);

  /* We need a sparse mask <sm>.
   * To get it, run checkpointed Fwd/Bck/Decoding
   */
  cx = p7_checkptmx_Create(hmm->M, sq->n, ESL_MBYTES(32));
  sm = p7_sparsemask_Create(gm->M, sq->n, p7_VDEFAULT);
  if (esl_opt_GetBoolean(go, "-a")) 
    p7_sparsemask_AddAll(sm);
  else {
    p7_ForwardFilter (sq->dsq, sq->n, om, cx, /*fsc=*/NULL);
    p7_BackwardFilter(sq->dsq, sq->n, om, cx, sm, p7_SPARSIFY_THRESH);
  }

  /* Allocate DP matrices, traceback */
  sxf  = p7_sparsemx_Create(sm);
  sxd  = p7_sparsemx_Create(sm);
  asf  = p7_sparsemx_Create(sm);
  asb  = p7_sparsemx_Create(sm);
  asd  = p7_sparsemx_Create(sm);
  tr   = p7_trace_Create();

  /* First pass analysis */
  p7_SparseViterbi (sq->dsq, sq->n, gm, sm,  sxf, tr, &vsc);
  p7_SparseForward (sq->dsq, sq->n, gm, sm,  sxf,     &fsc);
  p7_SparseBackward(sq->dsq, sq->n, gm, sm,  sxd,     NULL);
  p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxd,     sxd);


  /* MPAS */
  p7_sparse_anchors_SetFromTrace(sxd, tr, vanch);
  p7_trace_Reuse(tr);
  p7_sparse_Anchors(rng, sq->dsq, sq->n, gm,
		    vsc, fsc, sxf, sxd, vanch,
		    tr, &wrk, ah,
		    asf, anch, &asc_f, 
		    NULL);

  /* Remaining ASC calculations; MPAS did <asf> for us and we need <asb>, <asd> */
  p7_sparse_asc_Backward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asb, &asc_b);
  p7_sparse_asc_Decoding(sq->dsq, sq->n, gm, anch->a, anch->D, asc_f, asf, asb, asd);
  // p7_spascmx_DumpSpecials(stdout, asd, anch->a, anch->D);

  /* Envelope determination: ia..ib, oa..ob, env_sc are set */
  p7_sparse_Envelopes(sq->dsq, sq->n, gm, anch->a, anch->D, asf, asd, env);

  /* Null2 calculation 
   * Ultimately we need to handle wrk2, null2 tmp memory space better.
   */
  ESL_ALLOC(wrk2,  (gm->M+1) * sizeof(float));
  ESL_ALLOC(null2, (gm->abc->Kp) * sizeof(float));
  p7_sparse_Null2(sq->dsq, sq->n, gm, asd, env, &wrk2, null2);

  p7_envelopes_Dump(stdout, env);

  /* flowthrough is deliberate */
 ERROR:
  if (wrk)   free(wrk);
  if (wrk2)  free(wrk2);
  if (null2) free(null2);
  p7_envelopes_Destroy(env);
  p7_anchorhash_Destroy(ah);
  p7_anchors_Destroy(anch);
  p7_anchors_Destroy(vanch);
  p7_sparsemx_Destroy(asd);
  p7_sparsemx_Destroy(asb);
  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemx_Destroy(sxf);
  p7_trace_Destroy(tr);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(cx);
  p7_profile_Destroy(gm); 
  p7_oprofile_Destroy(om);
  esl_sqfile_Close(sqfp); 
  esl_sq_Destroy(sq);
  p7_hmmfile_Close(hfp);  
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPARSE_NULL2_EXAMPLE*/
/*---------------------- end, example ---------------------------*/


