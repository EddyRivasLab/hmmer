/* Envelope determination.
 * Reference implementation: debugging, development, pedagogy.
 * 
 * We've defined an anchor set, and we've done anchor-set-constrained
 * (ASC) decoding. Now we use marginalized ASC decoding data to
 * determine our "envelope" information for each domain:
 *     - whether domain d is glocal or local;
 *     - ha..hb, the bounds of the "homologous region";
 *     - ia..ib, the bounds of the "envelope": the alignable homologous region
 *     - oa..ob, the bound of the "outer envelope"
 *     - env_sc, the "envelope score".
 *     
 * After we have the envelope, the next step is anchor/envelope
 * constrained (AEC) alignment.
 *
 * Contents:
 *    1. Envelope definition API call.
 *    2. Internal functions: steps of envelope definition.
 *    3. Footnotes.
 *    4. Unit tests.
 *    5. Test driver.
 *    6. Example.
 */
#include "p7_config.h"

#include "easel.h"

#include "base/p7_envelopes.h" 
#include "base/p7_anchors.h"
#include "base/p7_profile.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_asc_fwdback.h"


static int envcoords(P7_ENVELOPES *env, int D, const P7_REFMX *apd);
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
 *            <anch> defining <D> domains in the sequence. We've
 *            calculated the anchor-set-constrained (ASC) posterior
 *            decoding UP and DOWN matrices <apu> and <apd>, using ASC
 *            Forward matrices <afu> and <afd>. Now, calculate the
 *            envelope for each domain.
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
 *
 * Notes:     See footnote [1] for definitions of the elements of the
 *            envelope.
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
      env->arr[d].ka    = 0;
      env->arr[d].kb    = 0;
      env->arr[d].flags = 0;
    }
  env->D = D;
  env->L = L;
  env->M = gm->M;
  p7_envelope_SetSentinels(env->arr, D, L, gm->M);
  
  if ((status = envcoords(env, D, apd))        != eslOK) goto ERROR;  // ia,ib (also ha,hb which aren't stored)
  if ((status = glocality(env, D, apd))        != eslOK) goto ERROR;  // is_glocal
  if ((status = outcoords(env, D, apd, 0.005)) != eslOK) goto ERROR;  // oa,ob

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
	    offset       = env->arr[d].oa - 1;   // oa..ob is our subseq; its len = ob-oa+1, or equivalently, ob-offset
	    reanch[1].i0 = anch[d].i0 - offset;  // You can't use <anch> itself because of coord offset
	    reanch[1].k0 = anch[d].k0;

	    p7_anchor_SetSentinels(reanch, 1, env->arr[d].ob-offset, gm->M);  // sentinel [D+1] has i0 = L+1, and our subseq L is new, so we must reinit sentinels

	    status = p7_ReferenceASCForward( dsq+offset, env->arr[d].ob-offset, gm,
					     reanch, 1, afu, afd, &(env->arr[d].env_sc));
	    if (status != eslOK) goto ERROR;

	    /* Add back tNN/JJ/CC terms for <L-Ld> flanking residues, 
	     * as if this domain were alone in the sequence.
	     */
	    env->arr[d].env_sc += ( L - env->arr[d].ob + offset ) * gm->xsc[p7P_N][p7P_LOOP];

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
 * First:
 *   ha = min_{i = i0(d-1)+1..i0(d)} P(x_i in domain d) > 0.5
 *   hb = max_{i = i0(d)..i0(d+1)-1} P(x_i in domain d) > 0.5
 *  equiv 
 *   hb = min_{i = i0(d)..i0(d+1)-1} P(x_{i+1} \in domain d) <= 0.5
 *            
 * Then: 
 *   ia = 1 + \argmax_{i=ha-1..i0(d)-1} \rho(i,B)
 *   ib =     \argmax_{i=i0(d)..hb} \rho(i,E)
 *  
 * All specials including B,E values are in the DOWN decoding matrix
 * <apd>, which is why we only need to see the DOWN matrix here.
 */
static int
envcoords(P7_ENVELOPES *env, int D, const P7_REFMX *apd)
{
  float  pB, pE;
  float  maxB, maxE;
  int    d, i;
  
  for (d = 1; d <= D; d++)
    {
      /* <ia> start position */
      maxB = 0.;
      pB   = 0.;
      for (i = env->arr[d-1].i0; i < env->arr[d].i0; i++)
	{
	  pB += P7R_XMX(apd, i, p7R_B);            // now pB = P(x_{i+1} \in d)
	  if ( P7R_XMX(apd, i, p7R_B) > maxB && pB > 0.5 )  
	  { 
	    env->arr[d].ia = i+1;
	    maxB = P7R_XMX(apd, i, p7R_B);
	  }
	}
      
      /* <ib> end positions */
      maxE = 0.;
      pE   = 0.;
      for (i = env->arr[d].i0; i < env->arr[d+1].i0; i++)
	{
	  if ( P7R_XMX(apd, i, p7R_E) >= maxE && pE <= 0.5)  
	    {
	      env->arr[d].ib = i;
	      maxE = P7R_XMX(apd, i, p7R_E);
	    }	  
	  pE += P7R_XMX(apd, i, p7R_E);          // now pE = 1.0 - P(x_{i+1} \in d)
	}
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
 * matrix. Only necessary for G, though; \rho_L(i0(d)) = 0
 * by construction.
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
 * Determine the "outer envelope", oa..ob, for env score purposes.
 *
 * Outer envelope coords are defined in the same way as the homologous
 * region in envcoords(), but with a much lower threshold,
 * \epsilon. The objective is to identify coords that encompass
 * essentially all of the probability mass for this domain.
 * 
 *   oa = min_{i <= i0} P(x_i in domain d) >= \epsilon
 *   ob = max_{i >= i0} P(x_i in domain d) >= \epsilon
 *   
 * where the default epsilon is 0.005,
 *
 * and where P(x_i in domain d) =   
 *  P(x_i \in d) = 
 *           \sum_{j=i0(d-1)..i-1} \rho_B(j)      for i0(d-1)+1 <= i <= i0(d)       
 *     1.0 - \sum_{j=i0(d)..i-1}   \rho_E(j)      for     i0(d) <= i <= i0(d+1)-1
 *
 * i.e. a cumulative addition of the mass that's entering the homology
 * domain as we move left from the anchor, or a cumulative subtraction
 * of mass leaving the domain as we move right.
 * 
 * See footnote [3] for why we think \epsilon = 0.005 is a reasonable
 * default, and for the proof of when the fast envelope score
 * approximation is valid.
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
	  if (phomology < epsilon) break;       // if i is not in the domain...
	}
      env->arr[d].oa = i+1;                     // but i+1 was, so oa = i+1.
      leftchoked = ( P7R_XMX(apd, i, s) >= 1.0 - (2*epsilon) ? TRUE : FALSE );

      phomology = 1.0;
      s  = (d == D ? p7R_C : p7R_J);
      for (i = env->arr[d].i0; i < env->arr[d+1].i0; i++)   // at D=D, i0(D+1)=L+1 sentinel makes this i0(D)..L
	{
	  phomology -= P7R_XMX(apd, i, p7R_E);   // now phomology = P(x_i+1 in dom_d)
	  //	  printf("%4d %.4f\n", i, phomology);
	  if (phomology < epsilon) break;        // if i+1 is not in domain...
	}
      env->arr[d].ob  = i;	                 // but i=1 was, so ob = i.
      rightchoked = (P7R_XMX(apd, i, s) >= 1.0 - (2*epsilon) ? TRUE : FALSE );

      if (leftchoked && rightchoked)
	env->arr[d].flags |= p7E_ENVSC_APPROX;
    }

  /* Do the oa = min(oa, ia); ob = max(ob, ib) part that
   * deals w/ rare pathological case
   */
  for (d = 1; d <= D; d++)
    {
      if (env->arr[d].oa > env->arr[d].ia) env->arr[d].oa = env->arr[d].ia;
      if (env->arr[d].ob > env->arr[d].ib) env->arr[d].ob = env->arr[d].ib;
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
	  P7R_XMX(afd, env->arr[d].ob,   s2) -
	  P7R_XMX(afd, env->arr[d].oa-1, s1);

	env->arr[d].env_sc += 
	  gm->xsc[p7P_N][p7P_LOOP] * (env->arr[d].oa-1) +
	  gm->xsc[p7P_C][p7P_LOOP] * (L-env->arr[d].ob) +
	  gm->xsc[p7P_C][p7P_MOVE];
      }
  return eslOK;
}

/*****************************************************************
 * 3. Footnotes
 *****************************************************************
 *
 * [1] DEFINITIONS
 * 
 * Whether envelope <d> is glocal or local is defined by the
 * marginalized probabilities of using the G vs L state in the path
 * ensemble for domain <d>:
 *      pL = \sum_{i=i0(d-1)..i0(d)-1} \rho_L(i)
 *      pG = \sum_{i=i0(d-1)..i0(d)-1} \rho_G(i)
 *      Domain d is glocal if pG >= pL 
 *      
 * The "homologous region" ha..hb for domain d is defined as the
 * region that is more likely homologous to the model than not, as
 * determined by posterior decoding:
 *     ha = min_{i = i0(d-1)+1..i0(d)} P(x_i in domain d) > 0.5
 *     hb = max_{i = i0(d)..i0(d+1)-1} P(x_i in domain d) > 0.5
 * The threshold of 0.5 is chosen to prevent adjacent inferred domains
 * from overlapping. The <ha..hb> coords are only used in defining the
 * alignable homologous region <ia..ib>; they are not recorded in the
 * envelope data structure.
 * 
 * The "envelope" ia..ib defines the bounds of the inferred alignment.
 * These bounds are defined as the most probable start/end
 * positions within the homologous region ha..hb:
 *      ia = 1 + \argmax_{i=ha-1..i0(d)-1} \rho(i,B)
 *      ib =     \argmax_{i=i0(d)..hb}     \rho(i,E)
 * This helps make sure that glocal alignments end at plausible
 * endpoints, which probability mass alone does not assure.     
 * 
 * The "outer envelope" oa..ob includes the ia..ib envelope, plus any
 * additional residues with non-negligible probability mass. The outer
 * envelope is used to calculate envelope scores (see below).  Outer
 * envelope coords <oa>,<ob> are defined as:
 *      oa = min_{i <= i0} P(x_i in domain d) >= epsilon 
 *      ob = max_{i >= i0} P(x_i in domain d) >= epsilon 
 * Since epsilon < 0.5 (default 0.005), oa <= ha and ob >= hb.
 *        
 * The "envelope score" <envsc> is the raw nat score of the
 * ensemble of all paths through this domain's anchor and
 * no other, as if this were the only domain in the entire
 * sequence.
 *      
 * A domain is considered to be "distinct" if its outer envelope
 * does not overlap the outer envelope of an adjacent domain. 
 * Then we know that:
 *    \rho_{N|J}(oa-1) \geq 1 - 2 \epsilon; and
 *    \rho_{J|C}(ob)   \geq 1 - 2 \epsilon 
 * on its left and right sides. When these conditions are met, we can
 * calculate the envelope score by a fast approximation using
 * arithmetic on existing ASC Forward matrix values, and we are
 * guaranteed to obtain a score within $4\epsilon$ nats (default =
 * 0.02) of the true envelope score. The <p7E_ENVSC_APPROX> flag is
 * set when this has been done. Otherwise, if the domain is not 
 * distinct, we recalculate the envelope score by ASC Forward on its
 * oa..ob outer envelope.
 * 
 * [2] TRICKSY OFF-BY-ONE ISSUES 
 * 
 * When we calculate <oa>,<ob>, we rearrange the definition of <ob>
 * so that both calculations are minima over i:
 *    oa(d) = min_{i <= ia(d)} P(x_i \in domain d) >= epsilon    
 *    ob(d) = min_{i >= ib(d)} P(x_{i+1} \in domain d) < epsilon 
 * Then the calculation can work in a single left to right pass
 * of increasing i, by identifying the first i that satisfies the
 * criterion. We do the same for <ha>,<hb>.  
 *
 * We could obtain P(x_i \in domain d) by marginalization over the
 * posterior probabilities in homologous M,I emitting states at i.
 * But equivalently -- and more conveniently -- we can obtain it by
 * summing up a cumulative distribution over B(i) and E(i) start and
 * endpoint posterior probabilities.
 * 
 * This can get confusing because of off-by-one issues. First think
 * about what coords are in play. Domain d must include its anchor
 * i0(d), and cannot include the anchors of adjacent domains, so:
 * 
 *    possible startpoints range [ i0(d-1)+1 .. i0(d) ];
 *    possible endpoints range                [ i0(d) .. i0(d+1)-1 ].
 *    
 * (Anchor sentinels i0(0)=0 and i0(D+1)=L+1 make this work even for
 * first, last domains.)
 * 
 * Next, let's look at exactly how we have to do our sums over B(i)
 * and E(i), because there's an off-by-one issue to keep in mind.  The
 * probability that you end exactly at residue <i> is E(i), but the
 * probability that you start exactly at <i> is B(i-1). Thus:
 * 
 *  P(x_i \in d) = 
 *           \sum_{j=i0(d-1)..i-1} \rho_B(j)      for i0(d-1)+1 <= i <= i0(d)       
 *     1.0 - \sum_{j=i0(d)..i-1}   \rho_E(j)      for     i0(d) <= i <= i0(d+1)-1
 *           
 * That is, as you move left to right in <i> on the left (starting)
 * side of the domain, the probability that x_i is in the left side of
 * the domain is the sum of B's up to i-1, including B(i-1)->Mk(x_i)
 * starts at x_i. On the right (ending) side, the probability that
 * you're still in the domain starts at 1.0 at the anchor i0(d), and
 * afterwards is 1.0 - the sum of E up to i-1, including
 * Mk(x_{i-1}->E) exits on the previous row i-1.
 * 
 * Next, think about writing that down in pseudocode as two independent
 * calculations:
 * 
 * For <oa> start point:
 *    pB = 0.0
 *    for i = i0(d-1) to i0(d)-1  // note off by one!
 *      pB += B(i)                // pB is now P(x_i+1 \in d)
 *      if (pB >= epsilon)        // if x_i+1 has met threshold
 *        oa(d) = i+1             //   then i+1 is the first x_i satisfying it
 *        break
 *
 * For <ob> end point:
 *    pE = 0.0                      // pE starts with P(x_i0(d) \in d)=1.0 at the anchor
 *    for i = i0(d) to i0(d+1) - 1: // then for each possible endpoint i:
 *      pE += E(i)                  //   pE is now 1.0 - P(x_i+1 \in d)
 *      if (1.0 - pE < epsilon)     //   if x_i+1 has dropped below threshold
 *        ob(d) = i                 //     then i was the last x_i satisfying it.
 *        break
 * 
 * After candidate oa..ob endpoints have been determined, we need to
 * go back through and do oa = min(oa, ia) and ob = max(ib, ob).  This
 * is to handle the pathological case that may arise if a maximum
 * likely endpoint has probability < epsilon. This final step makes
 * sure the outer envelope includes the alignable region, and assures
 * oa <= ia <= i0 <= ib <= oe coordinate order.
 * 
 * Then the last thing to do is make this work as a single
 * calculation, over all domains and i in a single pass, which is what
 * the implementation does.
 * 
 * 
 * 
 * [3] PROOF OF APPROXIMATION FOR ENVELOPE SCORES
 * 
 * To prove that a domain is "distinct", so we can use a fast
 * approximation to obtain the envelope score, we want to identify
 * that the outer envelope bounds oa/ob also serve as "choke points"
 * in the N/C/J states, through which passes all but a negligible
 * amount of path mass.
 * 
 * Choke points don't have to exist, because path mass can flow into
 * (oa+1..i0) from a previous domain, and out of (i0..ob-1) to a next
 * domain; but if the domain is "well-defined", these flows are
 * negligible.
 * 
 * Let's make the notion of "negligible" more formal now.
 * 
 * From the domain #, we know whether the domain is flanked on the
 * left by N or J, and on the right by J or C; we call these X, Y
 * respectively. Let q be the total amount of probability mass on
 * paths into the previous (resp. next) domain.
 * 
 * For oa as our left bound:
 *   P(oa-1 in D) < epsilon              (1) // by construction; that's how we chose oa.
 *   X(oa-1) = 1 - P(oa-1 in D) - q      (2) // mass flows back, accumulates in state X; eventually flows back to prev domain
 *   1-X(oa-1) < epsilon + q             (3) // substitute (1) in (2)
 *   X(oa-1) >= 1 - epsilon - q.         (4) // rearrange (3)
 * If q <= epsilon:                      (5)
 *   X(oa-1) >= 1 - 2*epsilon.           (6) // substitute (5) in (4).
 * 
 * So, our definition of "well-definedness" of the left edge is that
 * if X(oa-1) >= 1-2*epsilon, we know that the amount of probability
 * mass in paths that enter to our right from a previous domain
 * (without bottlenecking through X(oa-1)) is <= epsilon, and we know
 * that the amount of probability mass that's still in the homology
 * model (outside of the outer envelope) is < epsilon. So we've lost
 * up to 2*epsilon of the probability mass in paths at the left edge.
 * 
 * Analogous analysis holds for right edge. 
 * 
 * Therefore, if the condition on X(oa-1) and Y(ob) hold, we know that
 * if we only look at paths that pass through X(oa-1) and Y(ob), we
 * are neglecting at most a relative mass of 4\epsilon in other paths
 * that use this domain's anchor.
 * 
 * For \epsilon = 0.005, we lose up to 0.02 in total mass of the ensemble
 * we count, corresponding to a maximum log-odds score loss of ~0.02 nats;
 * we consider this to be a negligible score difference.
 */

/*****************************************************************
 * 4. Unit tests
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
 *       1 <= oa <= ia <= i0 <= ib <= ob <= L
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
	if (! (1 <= env->arr[d].oa &&
	       env->arr[d].oa <= env->arr[d].ia  &&
	       env->arr[d].ia <= env->arr[d].i0  &&
	       env->arr[d].i0 <= env->arr[d].ib  &&
	       env->arr[d].ib <= env->arr[d].ob &&
	       env->arr[d].ob <= sq->n)) esl_fatal(msg);

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
	  gtr->sqfrom[0] >= env->arr[1].oa &&    // in <gtr>, domains are 0..D-1; in <env>, 1..D
	  gtr->sqto[0]   <= env->arr[1].ob)
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
 *     2. For each domain, oa==ia, ob==ib, and these coords
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
		 env->arr[d].ia == env->arr[d].oa)) esl_fatal(msg);
	  if (! (env->arr[d].ib == gtr->sqto[d-1] &&
		 env->arr[d].ib == env->arr[d].ob)) esl_fatal(msg);
	}

      /* Test 3. If D == 1, envelope score == trace score. */
      if (D == 1 &&  esl_FCompare_old(env->arr[1].env_sc, gsc, tol) != eslOK) esl_fatal(msg);

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
 * 5. Test driver
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
 * 6. Example
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
  prm.nmax_sampling  = FALSE;
  prm.be_verbose     = FALSE;

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

