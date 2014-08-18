/* Reference implementation of envelope definition.
 * 
 * Envelope definition step comes after we define the anchor set, and
 * before we do the alignment.
 * 
 * Contents:
 *    1. Envelope definition API call.
 *    2. Internal functions: steps of envelope definition.
 *    3. Unit tests.
 *    4. Test driver.
 *    5. Example.
 *    6. Copyright and license information.
 */
#include "p7_config.h"

#include "easel.h"

#include "base/p7_envelopes.h" 
#include "base/p7_anchors.h"
#include "base/p7_profile.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_asc_fwdback.h"


static int envcoords(P7_ENVELOPES *env, int D, const P7_REFMX *apd, float threshold);
static int glocality(P7_ENVELOPES *env, int D, const P7_REFMX *apd);
static int outcoords(P7_ENVELOPES *env, int D, const P7_REFMX *apd, float epsilon);
static int approxsc (P7_ENVELOPES *env, int D, const P7_REFMX *afd, const P7_PROFILE *gm);

/*****************************************************************
 * 1. Envelope definition API call.
 *****************************************************************/

/* Function:  p7_reference_Envelopes()
 * Synopsis:  Define domain envelopes, given anchors.
 *
 * Purpose:   We're doing a comparison of profile <gm> and digital
 *            sequence <dsq> of length <L>. We have an anchor set
 *            <anch> defining <D> domains in the sequence (perhaps by
 *            using the most probable anchor set (MPAS) algorithm; see
 *            <p7_reference_Anchors()>). We've calculated the
 *            anchor-set-constrained (ASC) posterior decoding UP and
 *            DOWN matrices <apu> and <apd>, using ASC Forward
 *            matrices <afu> and <afd>. Now, calculate the envelope
 *            for each domain. 
 * 
 *            Return the envelope data in <env>. All the fields in the
 *            envelope structures are set except for the <ka> and <kb>
 *            fields, the alignment start/stop coords on the model.
 *            These two fields will be set when alignments are
 *            determined from the envelopes.
 *            
 *            <env> will be reallocated here as needed.
 *            
 *            The <afu> and <afd> matrices are probably re-used and
 *            overwritten here, for doing envelope score calculations
 *            -- Forward scores on isolated single domains. Caller
 *            must assume that their data has been invalidated.
 *            
 *            The envelope coords <ia>,<ib> are defined as the
 *            outermost coords that have p >= <threshold>
 *            (default=0.5) of generating x_i as part of the
 *            homologous domain.
 *            
 *            The envelope is defined as a local or glocal alignment
 *            by calculating the total marginal probability of using L
 *            or G in the ASC ensemble for this domain, and asking
 *            which is greater. If pG >= pL, the <p7E_IS_GLOCAL> flag
 *            is set in the envelope's <flags>.
 *            
 *            The 'outer envelope coords' <oea>, <oeb> are defined the
 *            same as the envelope coords, but for a lower p >=
 *            <epsilon> (default 0.005).
 *            
 *            The domain is considered to be 'well defined' (well
 *            separated from paths of adjacent domains) if
 *            D([N|J](oea-1)) >= 1-2*epsilon and D([J|C](oeb)) >=
 *            1-2*epsilon. If these conditions are met, we can
 *            calculate the envelope score (see below) by a fast
 *            approximation, just by arithmetic on the ASC Forward
 *            matrix, and we are guaranteed to obtain a score within
 *            ~4\epsilon nats (default=0.02) of the true envelope
 *            score. If these conditions are met, the
 *            <p7E_ENVSC_APPROX> flag is set in the envelope's
 *            <flags>, and the envelope score is calculated by the
 *            approximate method.
 *            
 *            The envelope score is the raw nat score of the ensemble
 *            of all paths through this domain's anchor, while erasing
 *            the influence of all other homology domains, as if this
 *            were the only domain in the entire sequence. For domains
 *            with the <p7E_ENVSC_APPROX> flag, it is calculated as
 *            above; for all other domains, it is calculated by an ASC
 *            Forward calculation on the isolated domain sequence
 *            <oea..oeb>, with the length model of <gm> left at <L>,
 *            and adding back the missing tNN/CC/JJ contributions for
 *            <L-Ld> residues outside the domain.
 *            
 * Args:      dsq  : digital sequence, 1..L
 *            L    : length of <dsq>
 *            gm   : profile, with length model already set to <L>
 *            anch : anchor set for <dsq>/<gm> comparison
 *            D    : number of domains defined by the anchor set
 *            apu  : ASC posterior decoding UP matrix
 *            apd  : ASC posterior decoding DOWN matrix
 *            afu  : ASC Forward UP matrix
 *            afd  : ASC Forward DOWN matrix
 *            env  : RETURN : envelope data for all <D> domains
 *
 * Returns:   <eslOK> on success. <env> may have been reallocated,
 *            and now contains the envelope data. <afu> and <afd> 
 *            must be assumed to be overwritten; caller can <_Reuse()>
 *            or <_Destroy()> them.
 *
 * Throws:    <eslEMEM> on allocation failure. State of data in <afu>, <afd>, 
 *            and <env> is undefined; they can be Reuse'd or Destroy'd.
 */
int
p7_reference_Envelopes(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D,
		       const P7_REFMX *apu, const P7_REFMX *apd, P7_REFMX *afu, P7_REFMX *afd, P7_ENVELOPES *env)
{
  P7_ANCHOR reanch[3];  // If we recalculate env score on isolated domain, must offset its anchor. [3] because of anchors.
  int       offset;     // Offset (in residues) to the domain start
  int       d;          // Index of domain
  int       status;   

  /* Reallocate <env> if needed, and initialize it */
  if (( status = p7_envelopes_Reinit(env, D)) != eslOK) goto ERROR;
  for (d = 1; d <= D; d++)
    {
      env->arr[d].i0    = anch[d].i0;
      env->arr[d].k0    = anch[d].k0;
      env->arr[d].alia  = 0;
      env->arr[d].alib  = 0;
      env->arr[d].ka    = 0;
      env->arr[d].kb    = 0;

      env->arr[d].flags = 0;
    }
  env->D = D;
  env->L = L;
  env->M = gm->M;
  p7_envelope_SetSentinels(env->arr, D, L, gm->M);
  
  if ((status = envcoords(env, D, apd, 0.5))   != eslOK) goto ERROR;
  if ((status = glocality(env, D, apd))        != eslOK) goto ERROR;
  if ((status = outcoords(env, D, apd, 0.005)) != eslOK) goto ERROR;

  if (D > 0)
    {
      /* All approximate scores must be calculated first, while we
       * still have the original <afd> matrix
       */
      if ((status = approxsc(env, D, afd, gm)) != eslOK) goto ERROR;

      /* Now, we can reuse the <afu>/<afd> matrix pair for any exact
       * env score calculations 
       */
      p7_refmx_Reuse(afu); p7_refmx_Reuse(afd);
      
      for (d = 1; d <= D; d++)
	if (! (env->arr[d].flags & p7E_ENVSC_APPROX))
	  {                                      
	    offset       = env->arr[d].oea - 1;  // oea..oeb is our subseq; its len = oeb-oea+1, or equivalently, oeb-offset
	    reanch[1].i0 = anch[d].i0 - offset;  // You can't use <anch> itself because of coord offset
	    reanch[1].k0 = anch[d].k0;

	    p7_anchor_SetSentinels(reanch, 1, env->arr[d].oeb-offset, gm->M);  // sentinel [D+1] has i0 = L+1, and our subseq L is new, so we must reinit sentinels

	    status = p7_ReferenceASCForward( dsq+offset, env->arr[d].oeb-offset, gm,
					     reanch, 1, afu, afd, &(env->arr[d].env_sc));
	    if (status != eslOK) goto ERROR;

	    /* Add back tNN/JJ/CC terms for <L-Ld> flanking residues, 
	     * as if this domain were alone in the sequence.
	     */
	    env->arr[d].env_sc += ( L - env->arr[d].oeb + offset ) * gm->xsc[p7P_N][p7P_LOOP];

	    p7_refmx_Reuse(afu); p7_refmx_Reuse(afd);
	  }
    }
  /* What if there's 0 domains in the sequence?  This will still all
   * work fine, and the returned <env> will be empty, with env->n = 0.
   */
  return eslOK;
  
 ERROR:
  return status;
}




/*****************************************************************
 * 2. Internal functions, implementing most of the steps
 *****************************************************************/

/* envcoords()
 * 
 * Sets <ia>,<ib> fields for all D domains.
 * 
 *   ia = min_{i <= i0) P(x_i in domain d) >= threshold
 *   ib = max_{i >= i0) P(x_i in domain d) >= threshold
 *   
 * where P(x_i in domain d) is the probability that x_i is emitted by
 * an M/I state on a path that includes the anchor in the same domain;
 * it is calculated cumulatively, by subtracting mass that leaves 
 * the homology domain:
 * 
 * P(x_i in d) = 
 *   i <= i0 :  1.0 - \sum_{j=i..i0-1} B(j)  
 *   i >= i0 :  1.0 - \sum_{j=i0..i-1} E(j)
 * 
 * where B(j) and E(j) are ASC decoding probabilities on row j.
 * 
 * All specials including B,E values are in the DOWN decoding matrix
 * <apd>, which is why we only need to see the DOWN matrix here.
 *
 * The default threshold is 0.5, but saying that is the caller's
 * job.
 */
static int
envcoords(P7_ENVELOPES *env, int D, const P7_REFMX *apd, float threshold)
{
  float  phomology;
  int    d, i;
  
  for (d = 1; d <= D; d++)
    {
      /* Determine ia, envelope start position;
       *  min_{i <= i0} P(x_i in domain d) >= threshold
       */
      phomology = 1.0;
      for (i = env->arr[d].i0 - 1; i >= env->arr[d-1].i0; i--)  // at d=1, i0(0)=0 sentinel makes this i0(1)-1 down to 0
	{
	  phomology -= P7R_XMX(apd, i, p7R_B);   // now phomology = P(res i in domain d)
	  if (phomology < threshold) break;      // if i is not in domain...
	}
      env->arr[d].ia = i+1;                      // but i+1 was, so ia = i+1.

      /* Determine ib, envelope end position;
       *   max_{i >= i0} P(x_i in domain d) >= threshold
       */
      phomology = 1.0;
      for (i = env->arr[d].i0; i < env->arr[d+1].i0; i++) // at d=D, i0(D+1)=L+1 sentinel makes this i0(d) to L
	{
	  phomology -= P7R_XMX(apd, i, p7R_E);    // now phomology = P(x_i+1 in dom_d)
	  if (phomology < threshold) break;       // if i+1 is not in domain...
	}
      env->arr[d].ib = i;		          // but i was, so ib = i.
    }
  return eslOK;
}



/* glocality()
 *
 * Determine whether each domain is glocal or local,
 * by marginalization over ensemble.
 *
 * Specifically, for each domain d,
 *    pL = \sum_i=i0(d-1)..i0(d)-1 P(i,L)
 *    pG = \sum_i=i0(d-1)..i0(d)-1 P(i,G)
 * 
 * If pG >= pL, it's a glocal alignment; set the p7E_IS_GLOCAL
 * flag for this domain envelope.
 *  
 * We do need the i0(d)-1 row, though it may not look like it at
 * first. From a B on i0(d)-1 row, we need to use B->G->DDDDk0-1 path
 * to be able to reach the anchor Mk on row i0(d). In decoding, we
 * wing unfold the G->DD->Mk paths, so this path exists in a decoding
 * matrix. Only for G, though; P(L, i0(d)) is 0.
 * 
 * (Because we're working in an ASC decoding matrix, not a standard
 * one, the mute partial cycle flaw is not relevant here.)
 */
static int
glocality(P7_ENVELOPES *env, int D, const P7_REFMX *apd)
{
  float pL, pG;
  int   d, i;

  for (d = 1; d <= D; d++)
    {
      pL = pG = 0.0;  

      for (i = env->arr[d-1].i0; i < env->arr[d].i0; i++)  // at d=1, i0(0)=0 sentinel makes this 0..i0(1)-1
	{
	  pL += P7R_XMX(apd, i, p7R_L);
	  pG += P7R_XMX(apd, i, p7R_G);
	}
      if (pG >= pL) env->arr[d].flags |= p7E_IS_GLOCAL;
    }
  return eslOK;
}


/* outcoords()
 * Determine the "outer envelope", oea..oeb, for env score purposes.
 *
 * Outer envelope coords are defined in the same way as the envelope
 * coords, but with a much lower threshold, \epsilon. The objective is
 * to identify coords that encompass essentially all of the probability
 * mass for this domain. 
 * 
 *   oea = min_{i <= i0} P(x_i in domain d) >= \epsilon (default 0.005)
 *   oeb = max_{i >= i0} P(x_i in domain d) >= \epsilon (default 0.005)
 *   
 * where P(x_i in domain d) =   
 *   i <= i0 :  1.0 - \sum_{j=i..i0-1} B(j)  
 *   i >= i0 :  1.0 - \sum_{j=i0..i-1} E(j)
 * 
 * i.e. a cumulative subtraction of the mass that's leaving the
 * homology domain as we move left or right from the anchor.
 * 
 * We are especially interested in whether we can identify flanking
 * "choke points" in the N/C/J states, through which passes all but a
 * negligible amount of path mass. If we can prove that oae and oeb
 * also define choke points, then we can use a fast approximation to
 * calculate the envelope score for this domain's ensemble, using
 * algebra in the Forward matrix of the entire ensemble. That lets us
 * avoid recalculating Forward on the isolated domain. 
 * 
 * Choke points don't have to exist, because path mass can flow into
 * (oea+1..i0) from a previous domain, and out of (i0..oeb-1) to a
 * next domain. If the domain is "well-defined", these flows are
 * negligible. 
 * 
 * Let's make that notion more formal now.
 * 
 * From the domain #, we know whether the domain is flanked on the
 * left by N or J, and on the right by J or C; call these X, Y
 * respectively. Let q be the total amount of probability mass on
 * paths into the previous (resp. next) domain.
 * 
 * For oea as our left bound:
 *   P(oea-1 in D) < epsilon              (1) // by construction; that's how we chose oae.
 *   X(oea-1) = 1 - P(oea-1 in D) - q     (2) // mass flows back, accumulates in state X; eventually flows back to prev domain
 *   1-X(oea-1) < epsilon + q             (3) // substitute (1) in (2)
 *   X(oea-1) >= 1 - epsilon - q.         (4) // rearrange (3)
 * If q <= epsilon:                       (5)
 *   X(oea-1) >= 1 - 2*epsilon.           (6) // substitute (5) in (4).
 * 
 * So, our definition of "well-definedness" of the left edge is that
 * if X(oea-1) >= 1-2*epsilon, we know that the amount of probability
 * mass in paths that enter to our right from a previous domain
 * (without bottlenecking through X(oea-1) is <= epsilon, and we know
 * that the amount of probability mass that's still in the homology
 * model is < epsilon. We've lost up to 2*epsilon of the probability
 * mass in paths at the left edge.
 * 
 * Similar analysis holds for right edge. 
 * 
 * Therefore, if the condition on X(oea-1) and Y(oeb) hold, we know
 * that if we only look at paths that pass through X(oea-1) and
 * Y(oeb), we are neglecting at most a mass of 4\epsilon in other
 * paths that use this domain's anchor.
 * 
 * For \epsilon = 0.005, we lose up to 0.02 in total mass of the ensemble
 * we count, corresponding to a maximum log-odds score loss of ~0.02 nats;
 * a negligible score difference.
 */
static int
outcoords(P7_ENVELOPES *env, int D, const P7_REFMX *apd, float epsilon)
{
  int   d, i, s;
  float phomology;
  int   leftchoked, rightchoked; 

  for (d = 1; d <= D; d++)
    {
      phomology = 1.0;
      s  = (d == 1 ? p7R_N : p7R_J);
      for (i = env->arr[d].i0 - 1; i >= env->arr[d-1].i0; i--)   // at d=1, i0(0)=0 sentinel makes this i0(1)-1 down to 0
	{
	  phomology -= P7R_XMX(apd, i, p7R_B);  // now phomology = P(x_i in domain d)
	  printf("%4d %.4f\n", i, phomology);
	  if (phomology < epsilon) break;       // if i is not in the domain...
	}
      env->arr[d].oea = i+1;                    // but i+1 was, so oea = i+1.
      leftchoked = ( P7R_XMX(apd, i, s) >= 1.0 - (2*epsilon) ? TRUE : FALSE );

      phomology = 1.0;
      s  = (d == D ? p7R_C : p7R_J);
      for (i = env->arr[d].i0; i < env->arr[d+1].i0; i++)   // at D=D, i0(D+1)=L+1 sentinel makes this i0(D)..L
	{
	  phomology -= P7R_XMX(apd, i, p7R_E);   // now phomology = P(x_i+1 in dom_d)
	  //	  printf("%4d %.4f\n", i, phomology);
	  if (phomology < epsilon) break;        // if i+1 is not in domain...
	}
      env->arr[d].oeb  = i;	                 // but i=1 was, so oeb = i.
      rightchoked = (P7R_XMX(apd, i, s) >= 1.0 - (2*epsilon) ? TRUE : FALSE );

      if (leftchoked && rightchoked)
	env->arr[d].flags |= p7E_ENVSC_APPROX;
    }
  return eslOK;
}



/* approxsc()
 * Where we can, calculate the env score by a fast approximation
 * in the main ASC Forward matrix.
 * 
 * If t_NN = t_JJ = t_CC, then we can often calculate the envelope
 * score of a single domain by an accurate approximation from the ASC
 * Forward matrix for the whole sequence.
 * 
 * Assume that we know that Xi and Yj are "choke points" flanking our
 * domain, such that >= 1.0-2*epsilon of the posterior path probability
 * flows through each of them.
 * 
 * For d=1, X=N, Y=J; for d=D, X=J, Y=C; for internal d, X=Y=J.
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
approxsc(P7_ENVELOPES *env, int D, const P7_REFMX *afd, const P7_PROFILE *gm) 
{
  int d;
  int L = afd->L;
  int s1, s2;

  for (d = 1; d <= D; d++)
    if (env->arr[d].flags & p7E_ENVSC_APPROX)
      {
	s1 = (d == 1 ? p7R_N : p7R_J);
	s2 = (d == D ? p7R_C : p7R_J);

	env->arr[d].env_sc = 
	  P7R_XMX(afd, env->arr[d].oeb,   s2) -
	  P7R_XMX(afd, env->arr[d].oea-1, s1);

	env->arr[d].env_sc += 
	  gm->xsc[p7P_N][p7P_LOOP] * (env->arr[d].oea-1) +
	  gm->xsc[p7P_C][p7P_LOOP] * (L-env->arr[d].oeb) +
	  gm->xsc[p7P_C][p7P_MOVE];
      }
  return eslOK;
}

/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7REFERENCE_ENVELOPES_TESTDRIVE
#include "hmmer.h"

/* "generation" test
 * Compare a randomly sampled profile to sequences sampled
 * from that profile.
 * 
 * This test is not very stringent, because we don't know the "true"
 * envelopes. Rather, this is more of a test that nothing obviously
 * bad happens, like a crash, or obviously incorrect data.
 * 
 * We test:
 *    1. Seq coordinates of each envelope are coherent:
 *       1 <= oea <= ia <= i0 <= ib <= oeb <= L
 *       
 *    2. Envelopes do not overlap (assuming default threshold of
 *       0.5 when defining them):
 *         ia(d) > ib(d-1)  for d = 2..D
 *       (Outer envelopes, in contrast, can overlap.)
 *       
 *    3. envsc(d) <= asc_sc <= fwdsc.
 *    
 *    4. If D=1 (single domain) in both the generated trace
 *       and the inferred envelopes, and the domain coords in 
 *       the trace are encompassed by the outer envelope,
 *       then envsc(d) >= generated trace score.
 */
static void
utest_generation(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int N)
{
  char             msg[] = "reference_envelopes:: generation unit test failed";
  ESL_SQ          *sq    = esl_sq_CreateDigital(abc);
  P7_BG           *bg    = p7_bg_Create(abc);
  P7_HMM          *hmm   = NULL;
  P7_PROFILE      *gm    = p7_profile_Create(M, abc);
  P7_TRACE        *gtr   = p7_trace_Create();            // generated trace
  P7_TRACE        *vtr   = p7_trace_Create();            // Viterbi trace
  P7_REFMX        *rxf   = p7_refmx_Create(M, 20);       // Fwd, Vit ~~> ASC Decode UP
  P7_REFMX        *rxd   = p7_refmx_Create(M, 20);       // Bck, Decode ~~> ASC Decode DOWN
  P7_REFMX        *afu   = p7_refmx_Create(M, 20);       // ASC Fwd UP
  P7_REFMX        *afd   = p7_refmx_Create(M, 20);       // ASC Fwd DOWN
  P7_REFMX        *apu   = rxf;                          // for 'clarity' we use two names for this mx
  P7_REFMX        *apd   = rxd;                          //   ... and this one too.
  float           *wrk   = NULL;
  P7_ANCHORS      *anch  = p7_anchors_Create();
  P7_ANCHORHASH   *ah    = p7_anchorhash_Create();
  P7_ENVELOPES    *env   = p7_envelopes_Create();
  float            tol   = 0.001;
  float  gsc, fsc, asc;
  int    idx;
  int    d;
  
  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);

  for (idx = 0; idx < N; idx++)
    {
      /* Emit sequence from model, using an arbitrary length model of <M>;
       * restrict the emitted sequence length to 6M, arbitrarily, to 
       * keep it down to something reasonable.
       */
      if ( p7_profile_SetLength(gm, M) != eslOK) esl_fatal(msg);
      do {
	esl_sq_Reuse(sq);
	if (p7_ProfileEmit(rng, hmm, gm, bg, sq, gtr) != eslOK) esl_fatal(msg);
      } while (sq->n > M * 6); 
      if (p7_trace_Index   (gtr)                      != eslOK) esl_fatal(msg);
      if (p7_trace_Score   (gtr, sq->dsq, gm, &gsc)   != eslOK) esl_fatal(msg);

      /* Reset the length model to the actual length sq->n, then
       * put it through the domain postprocessing analysis pipeline
       */
      if ( p7_profile_SetLength(gm, sq->n)                          != eslOK) esl_fatal(msg);
     
      /* First pass analysis */
      if ( p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf, vtr, NULL) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd,      NULL) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd)  != eslOK) esl_fatal(msg);

      /* Anchor determination (MPAS algorithm) */
      if ( p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, vtr, &wrk, ah,
				afu, afd, anch, &asc, NULL, NULL)  != eslOK) esl_fatal(msg);

      /* Reuse rxf,rxd as apu, apd; finish ASC analysis with Backward, Decoding */
      p7_refmx_Reuse(apu);  p7_refmx_Reuse(apd);
      if ( p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch->a, anch->D, apu, apd, NULL)               != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch->a, anch->D, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(msg);

      /* Envelope calculation */
      if ( p7_reference_Envelopes(sq->dsq, sq->n, gm, anch->a, anch->D, apu, apd, afu, afd, env) != eslOK) esl_fatal(msg);


      /* Test 1. Coords of each domain are coherent */
      if (anch->D != env->D) esl_fatal(msg);
      for (d = 1; d <= anch->D; d++)
	if (! (1 <= env->arr[d].oea &&
	       env->arr[d].oea <= env->arr[d].ia  &&
	       env->arr[d].ia  <= env->arr[d].i0  &&
	       env->arr[d].i0  <= env->arr[d].ib  &&
	       env->arr[d].ib  <= env->arr[d].oeb &&
	       env->arr[d].oeb <= sq->n)) esl_fatal(msg);

      /* Test 2. Envelopes do not overlap. */
      for (d = 1; d <= anch->D; d++)
	if (! (env->arr[d].ia > env->arr[d-1].ib)) esl_fatal(msg);

      /* Test 3. envsc(d) <= asc_sc <= fwdsc */
      for (d = 1; d <= anch->D; d++)
	if (! (env->arr[d].env_sc <= asc+tol && asc <= fsc+tol)) esl_fatal(msg);

      /* Test 4, only on D=1 case with generated trace's domain 
       * encompassed by the outer envelope 
       */
      if (gtr->ndom == 1 &&  anch->D   == 1 && 
	  gtr->sqfrom[0] >= env->arr[1].oea &&    // in <gtr>, domains are 0..D-1; in <env>, 1..D
	  gtr->sqto[0]   <= env->arr[1].oeb)
	if (! ( env->arr[1].env_sc >= gsc)) esl_fatal(msg);

      p7_envelopes_Reuse(env);
      p7_anchors_Reuse(anch);
      p7_anchorhash_Reuse(ah);
      p7_refmx_Reuse(rxf); p7_refmx_Reuse(rxd);
      p7_refmx_Reuse(afu); p7_refmx_Reuse(afd);
      p7_trace_Reuse(gtr); p7_trace_Reuse(vtr);
      esl_sq_Reuse(sq);
    }
      
  if (wrk) free(wrk);
  p7_envelopes_Destroy(env);
  p7_anchors_Destroy(anch);
  p7_anchorhash_Destroy(ah);
  p7_refmx_Destroy(afu); p7_refmx_Destroy(afd);
  p7_refmx_Destroy(rxf); p7_refmx_Destroy(rxd);
  p7_trace_Destroy(vtr); p7_trace_Destroy(gtr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_sq_Destroy(sq);
}

/* "singlemulti" test.
 * 
 * Use p7_modelsample_SinglePathedASC() to create a
 * profile/sequence/anchorset comparison that has only a single
 * possible path when anchor set constrained. Now the expected true
 * envelope(s) are known, from that single path, and we can compare.
 * 
 * In order to guarantee only one possible path, while allowing
 * multiple domains, the profile is limited to glocal-only.
 * 
 * We test:
 *     1. Trace and envelopes agree on number of domains.
 *     2. For each domain, oea==ia, oeb==ib, and these coords
 *        agree with the trace.
 *     3. In the case of a single domain (D=1), the envelope
 *        score == the trace score.
 */
static void
utest_singlemulti(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int N)
{
  char          msg[] = "reference_envelopes singlemulti unit test failed";
  P7_BG        *bg    = p7_bg_Create(abc);
  P7_HMM       *hmm   = NULL;
  P7_PROFILE   *gm    = NULL;
  ESL_DSQ      *dsq   = NULL;
  int           L;
  P7_TRACE     *gtr   = NULL;
  P7_ANCHOR    *anch  = NULL;
  int           D;
  P7_REFMX     *afu   = p7_refmx_Create(M, 20);
  P7_REFMX     *afd   = p7_refmx_Create(M, 20);
  P7_REFMX     *apu   = p7_refmx_Create(M, 20);
  P7_REFMX     *apd   = p7_refmx_Create(M, 20);
  P7_ENVELOPES *env   = p7_envelopes_Create();
  float         gsc;
  float         tol   = 0.001;
  int           idx;
  int           d;
  
  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_SinglePathedASC(rng, M, bg, &hmm, &gm, &dsq, &L, &gtr, &anch, &D, &gsc) != eslOK) esl_fatal(msg);

      if ( p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd, NULL)               != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCBackward(dsq, L, gm, anch, D, apu, apd, NULL)               != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCDecoding(dsq, L, gm, anch, D, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(msg);

      if ( p7_reference_Envelopes (dsq, L, gm, anch, D, apu, apd, afu, afd, env)      != eslOK) esl_fatal(msg);

      /* Test 1. Domain #'s agree */
      if (! (gtr->ndom == D && env->D == D)) esl_fatal(msg);

      /* Test 2. Envelope coords (and outer env coords) match trace.
       *         (Beware, trace domains are numbered 0..D-1, env domains are 1..D )
       */
      for (d = 1; d <= D; d++)
	{
	  if (! (env->arr[d].ia == gtr->sqfrom[d-1] &&
		 env->arr[d].ia == env->arr[d].oea)) esl_fatal(msg);
	  if (! (env->arr[d].ib == gtr->sqto[d-1] &&
		 env->arr[d].ib == env->arr[d].oeb)) esl_fatal(msg);
	}

      /* Test 3. If D == 1, envelope score == trace score. */
      if (D == 1 &&  esl_FCompare(env->arr[1].env_sc, gsc, tol) != eslOK) esl_fatal(msg);

      p7_envelopes_Reuse(env);
      p7_refmx_Reuse(afu); p7_refmx_Reuse(afd);
      p7_refmx_Reuse(apu); p7_refmx_Reuse(apd);

      free(dsq);
      free(anch);
      p7_trace_Destroy(gtr);
      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
    }
     
  p7_envelopes_Destroy(env);
  p7_refmx_Destroy(afu); p7_refmx_Destroy(afd);
  p7_refmx_Destroy(apu); p7_refmx_Destroy(apd);
  p7_bg_Destroy(bg);
}

#endif /*p7REFERENCE_ENVELOPES_TESTDRIVE*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7REFERENCE_ENVELOPES_TESTDRIVE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for reference envelopes inference";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  int             M    = 50;

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_generation (rng, M, abc, 10);  // test a bunch of seqs to try to make sure we exercise exact domain score recalculation
  utest_singlemulti(rng, M, abc, 10);

  fprintf(stderr, "#  status = ok\n");

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_ENVELOPES_TESTDRIVE*/


/*****************************************************************
 * 3. Example
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
  P7_ANCHORS    *anch     = p7_anchors_Create();
  P7_ANCHORHASH  *ah      = p7_anchorhash_Create();
  P7_ENVELOPES   *env     = p7_envelopes_Create();
  P7_REFMX       *rxf     = NULL;
  P7_REFMX       *rxd     = NULL;
  P7_REFMX       *afu     = NULL;
  P7_REFMX       *afd     = NULL;
  P7_REFMX       *apu     = NULL;
  P7_REFMX       *apd     = NULL;
  P7_TRACE       *tr      = NULL;
  float          *wrk     = NULL;
  P7_MPAS_PARAMS  prm;
  P7_MPAS_STATS   stats;
  float           fsc, vsc, asc, asc_b;
  int             status;

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
  p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, tr, &wrk, ah,
		       afu, afd, anch, &asc, &prm, &stats);

  
  //printf("# ASC Forward UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("# ASC Forward DOWN:\n"); p7_refmx_Dump(stdout, afd);

  /* We no longer need rxf and rxd. 
   * Use their space for apu/apd pair, which will briefly
   * hold ASC Backward matrices, then get used for ASC Decoding.
   */
  apu = rxf; p7_refmx_Reuse(apu);
  apd = rxd; p7_refmx_Reuse(apd);

  p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch->a, anch->D, apu, apd, &asc_b);
  
  //printf("# Backward score (raw, nats): %.2f\n", asc_b);
  //printf("# ASC Backward UP:\n");   p7_refmx_Dump(stdout, apu);
  //printf("# ASC Backward DOWN:\n"); p7_refmx_Dump(stdout, apd);

  /* ASC Decoding takes afu/afd and abu/abd as input;
   * overwrites abu/abd with decoding matrices
   */
  p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch->a, anch->D, afu, afd, apu, apd, apu, apd);

  //printf("# ASC Decoding UP matrix:\n");  p7_refmx_Dump(stdout, apu);
  //printf("# ASC Decoding DOWN:\n");       p7_refmx_Dump(stdout, apu);


  /* Envelope calculation needs to get four matrices:
   * ASC Decoding pair, apu/apd, and it will leave these constant;
   * ASC Forward pair,  afu/afd, and it will overwrite these.
   */
  p7_reference_Envelopes(sq->dsq, sq->n, gm, anch->a, anch->D, apu, apd, afu, afd, env);

  p7_envelopes_Dump(stdout, env);

  p7_envelopes_Destroy(env);
  p7_anchorhash_Destroy(ah);
  p7_anchors_Destroy(anch);
  if (wrk) free(wrk);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(rxd);
  p7_refmx_Destroy(rxf);
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
