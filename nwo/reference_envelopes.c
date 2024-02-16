/* Envelope determination: reference implementation.
 *
 * After determining the most probable anchor set for a profile/sequence comparison,
 * which defines how many domains are in the sequence and where, the next step is to
 * define the bounds of each domain: the "envelope", and its envelope score.
 *
 * Once we know the envelope for each domain, the next step is to obtain an optimal
 * alignment, constrained by its envelope and anchor (AEC alignment, for
 * anchor/envelope-constrained).
 *
 * Contents:
 *   1. Envelope definition: h4_reference_Envelopes()
 *   2. Internal functions implementing envelope definition steps
 *   3. Unit tests
 *   4. Test driver
 *   5. Example
 *
 * See reference_envelopes.md for notes.  
 */
#include <h4_config.h>

#include "easel.h"

#include "h4_anchorset.h"
#include "h4_envset.h"
#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_refmx.h"

#include "reference_asc.h"

static int envcoords(H4_ENVSET *env, const H4_REFMX *apd);
static int glocality(H4_ENVSET *env, const H4_REFMX *apd);
static int outcoords(H4_ENVSET *env, const H4_REFMX *apd, float epsilon);
static int approxsc (H4_ENVSET *env, const H4_REFMX *afd, const H4_MODE *mo);
static int exactsc  (H4_ENVSET *env, const ESL_DSQ *dsq, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *afu, H4_REFMX *afd);


/*****************************************************************
 * 1. Envelope definition
 *****************************************************************/

/* Function:  h4_reference_Envelopes()
 * Synopsis:  Define domain envelopes, given anchors.
 *
 * Purpose:   We're comparing digital sequence <dsq> of length <L> to profile <hmm> in
 *            comparison mode <mo>. We have an anchor set <anch> that defines a set
 *            of domains in the sequence. We've calculated the anchor-set-constrained
 *            (ASC) posterior decoding UP and DOWN matrices <apu> and <apd>, using
 *            ASC Forward matrices <afu> and <afd>. Now, calculate the envelope for
 *            each domain.
 * 
 *            Return the envelope data in <env>. Caller provides an allocated <env>
 *            structure; it will be reinitialized, and reallocated if needed. All the
 *            fields in the envelope structures are set here except for the <ka> and
 *            <kb> fields, the alignment start/stop coords on the model.  These two
 *            fields are set later when alignments are determined from the envelopes.
 *            (They're set to 0 here.)
 *            
 *            The <afu> and <afd> matrices are (probably) re-used and overwritten
 *            here, for doing envelope score calculations -- Forward scores on
 *            isolated single domains. (Only probably, because it's possible that all
 *            domains are "well defined" and suitable for having their scores
 *            calculated by the fast approximation.) Caller must assume that the data
 *            in <afu/afd> has been invalidated.
 *            
 * Args:      dsq  : digital sequence, 1..L
 *            L    : length of <dsq>
 *            hmm  : profile
 *            mo   : comparison mode, with length model set to <L>
 *            anch : anchor set for <dsq>/<hmm> comparison
 *            apu  : ASC posterior decoding UP matrix
 *            apd  : ASC posterior decoding DOWN matrix
 *            afu  : ASC Forward UP matrix
 *            afd  : ASC Forward DOWN matrix
 *            env  : RETURN : envelope data for all <D> domains [caller-provided space; reused/resized as needed]
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
h4_reference_Envelopes(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_ANCHORSET *anch,
		       const H4_REFMX *apu, const H4_REFMX *apd, H4_REFMX *afu, H4_REFMX *afd, H4_ENVSET *env)
{
  int status;
  
  /* Reallocate <env> if needed, and initialize it using the anchorset */
  if (( status = h4_envset_Resize(env, anch->D))         != eslOK) return status;
  if (( status = h4_envset_CopyFromAnchorset(anch, env)) != eslOK) return status;

  if ((status = envcoords(env, apd))                     != eslOK) return status;  // ia,ib (also ha,hb which aren't stored)
  if ((status = glocality(env, apd))                     != eslOK) return status;  // is_glocal flag
  if ((status = outcoords(env, apd, 0.005))              != eslOK) return status;  // oa,ob; and envsc_approx flag
  if ((status = approxsc (env, afd, mo))                 != eslOK) return status;  // env_sc for domains w/ envsc_approx flag set
  if ((status = exactsc  (env, dsq, hmm, mo, afu, afd))  != eslOK) return status;  // env_sc for the rest

  return eslOK;
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
envcoords(H4_ENVSET *env, const H4_REFMX *apd)
{
  float  pB, pE;
  float  maxB, maxE;
  int    d, i;
  
  for (d = 1; d <= env->D; d++)
    {
      /* <ia> start position */
      maxB = 0.;
      pB   = 0.;
      for (i = env->e[d-1].i0; i < env->e[d].i0; i++)
	{
	  pB += H4R_XMX(apd, i, h4R_B);            // now pB = P(x_{i+1} \in d)
	  if ( H4R_XMX(apd, i, h4R_B) > maxB && pB > 0.5 )  
	  { 
	    env->e[d].ia = i+1;
	    maxB = H4R_XMX(apd, i, h4R_B);
	  }
	}
      
      /* <ib> end positions */
      maxE = 0.;
      pE   = 0.;
      for (i = env->e[d].i0; i < env->e[d+1].i0; i++)
	{
	  if ( H4R_XMX(apd, i, h4R_E) >= maxE && pE <= 0.5)  
	    {
	      env->e[d].ib = i;
	      maxE = H4R_XMX(apd, i, h4R_E);
	    }	  
	  pE += H4R_XMX(apd, i, h4R_E);          // now pE = 1.0 - P(x_{i+1} \in d)
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
 * If pG >= pL, it's a glocal alignment; set the h4E_IS_GLOCAL
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
glocality(H4_ENVSET *env, const H4_REFMX *apd)
{
  float pL, pG;
  int   d, i;

  for (d = 1; d <= env->D; d++)
    {
      pL = pG = 0.0;  
      for (i = env->e[d-1].i0; i < env->e[d].i0; i++)  // at d=1, i0(0)=0 sentinel makes this 0..i0(1)-1
	{
	  pL += H4R_XMX(apd, i, h4R_L);
	  pG += H4R_XMX(apd, i, h4R_G);
	}
      if (pG >= pL) env->e[d].flags |= h4E_IS_GLOCAL;
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
outcoords(H4_ENVSET *env, const H4_REFMX *apd, float epsilon)
{
  int   d, i, s;
  float phomology;
  int   leftchoked, rightchoked; 

  for (d = 1; d <= env->D; d++)
    {
      phomology = 1.0;
      s  = (d == 1 ? h4R_N : h4R_J);
      for (i = env->e[d].i0 - 1; i >= env->e[d-1].i0; i--)   // at d=1, i0(0)=0 sentinel makes this i0(1)-1 down to 0
	{
	  phomology -= H4R_XMX(apd, i, h4R_B);  // now phomology = P(x_i in domain d)
	  if (phomology < epsilon) break;       // if i is not in the domain...
	}
      env->e[d].oa = i+1;                     // but i+1 was, so oa = i+1.
      leftchoked = ( H4R_XMX(apd, i, s) >= 1.0 - (2*epsilon) ? TRUE : FALSE );

      phomology = 1.0;
      s  = (d == env->D ? h4R_C : h4R_J);
      for (i = env->e[d].i0; i < env->e[d+1].i0; i++)   // at D=D, i0(D+1)=L+1 sentinel makes this i0(D)..L
	{
	  phomology -= H4R_XMX(apd, i, h4R_E);   // now phomology = P(x_i+1 in dom_d)
	  if (phomology < epsilon) break;        // if i+1 is not in domain...
	}
      env->e[d].ob  = i;	                 // but i=1 was, so ob = i.
      rightchoked = (H4R_XMX(apd, i, s) >= 1.0 - (2*epsilon) ? TRUE : FALSE );

      if (leftchoked && rightchoked)
	env->e[d].flags |= h4E_ENVSC_APPROX;
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
approxsc(H4_ENVSET *env, const H4_REFMX *afd, const H4_MODE *mo)
{
  int d;
  int L = afd->L;
  int s1, s2;

  for (d = 1; d <= env->D; d++)
    if (env->e[d].flags & h4E_ENVSC_APPROX)
      {
	s1 = (d == 1      ? h4R_N : h4R_J);
	s2 = (d == env->D ? h4R_C : h4R_J);

	env->e[d].env_sc =
          H4R_XMX(afd, env->e[d].ob,   s2) -
	  H4R_XMX(afd, env->e[d].oa-1, s1);

	env->e[d].env_sc += 
	  mo->xsc[h4_N][h4_LOOP] * (env->e[d].oa-1) +
	  mo->xsc[h4_C][h4_LOOP] * (L-env->e[d].ob) +
	  mo->xsc[h4_C][h4_MOVE];

        env->e[d].env_sc -= mo->nullsc;
      }
  return eslOK;
}


static int
exactsc(H4_ENVSET *env, const ESL_DSQ *dsq, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *afu, H4_REFMX *afd)
{
  H4_ANCHOR    a[3];
  H4_ANCHORSET reanch;
  int          d;
  int          offset;
  int          status;
  
  /* hack a little temporary 1-anchor anchorset, bypassing malloc */
  reanch.a        = a;
  reanch.D        = 1;
  reanch.nalloc   = 3;
  reanch.nredline = 3; 

  for (d = 1; d <= env->D; d++)
    if (! (env->e[d].flags & h4E_ENVSC_APPROX))
      {
        h4_refmx_Reuse(afu);
        h4_refmx_Reuse(afd);

        offset         = env->e[d].oa - 1;        // oa..ob is our subseq; its len = ob-oa+1, or equivalently, ob-offset
        reanch.a[1].i0 = env->e[d].i0 - offset;   // can't use <anch> itself because of coord offset
        reanch.a[1].k0 = env->e[d].k0;
        
        if ((status = h4_anchorset_SetSentinels(&reanch, env->e[d].ob-offset, env->M)) != eslOK) return status;
        if ((status = h4_reference_asc_Forward(dsq+offset, env->e[d].ob-offset, hmm, mo, &reanch, afu, afd, &(env->e[d].env_sc))) != eslOK) return status;

        /* Add back tNN/JJ/CC terms for <L-Ld> flanking residues, 
         * as if this domain were alone in the sequence.
         */
        env->e[d].env_sc += ( env->L - env->e[d].ob + offset ) * mo->xsc[h4_N][h4_LOOP];
         
        /* ASC Forward already included mo->nullsc, don't subtract it again.
         */
      }
  return eslOK;
}


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef h4REFERENCE_ENVELOPES_TESTDRIVE

#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "h4_path.h"
#include "h4_anchorset.h"
#include "h4_anchorhash.h"

#include "emit.h"
#include "modelsample.h"
#include "reference_asc.h"
#include "reference_dp.h"
#include "reference_mpas.h"

static void
utest_generation(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int N)
{
  char             msg[] = "reference_envelopes:: generation unit test failed";
  H4_PROFILE      *hmm   = NULL;
  H4_MODE         *mo    = h4_mode_Create();
  ESL_SQ          *sq    = esl_sq_CreateDigital(abc);
  H4_PATH         *gpi   = h4_path_Create();              // generated path
  H4_PATH         *pi    = h4_path_Create();              
  H4_REFMX        *rxf   = h4_refmx_Create(M, 20);        // Vit, Fwd ... then ASC Decode UP
  H4_REFMX        *rxd   = h4_refmx_Create(M, 20);        // Bck, Decode ... then ASC Decode DOWN
  H4_REFMX        *afu   = h4_refmx_Create(M, 20);        // ASC Forward UP
  H4_REFMX        *afd   = h4_refmx_Create(M, 20);        // ASC Forward DOWN
  H4_REFMX        *apu   = rxf;                           // for "clarity", we use two names for one shared matrix
  H4_REFMX        *apd   = rxd;                           //  ... this one too
  H4_ANCHORSET    *anch  = h4_anchorset_Create(0,0,0);
  H4_ANCHORHASH   *ah    = h4_anchorhash_Create();
  H4_ENVSET       *env   = h4_envset_Create(0,0,0); 
  float           *wrk   = NULL;                          // workspace needed (and managed) by MPAS
  float            tol   = 0.001;
  int              idx;
  float            gsc, vsc, fsc, asc;
  int              ia,ib;
  int              d;

  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);

  for (idx = 1; idx <= N; idx++)
    {
      /* Emit sequence from model, using an arbitrary length model of <M>,
       * and restrict the emitted sequence length to 6M, arbitrarily, to 
       * keep it down to something reasonable. Then reset length model
       * to the actual length.
       */
      if ( h4_mode_SetLength(mo, M)                   != eslOK) esl_fatal(msg);
      do {
        if ( h4_emit(rng, hmm, mo, sq, gpi)           != eslOK) esl_fatal(msg);
      } while (sq->n > M * 6); 
      if ( h4_path_Score(gpi, sq->dsq, hmm, mo, &gsc) != eslOK) esl_fatal(msg);
      if ( h4_mode_SetLength(mo, sq->n)               != eslOK) esl_fatal(msg);

      /* First pass analysis */
      if ( h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxf, pi,  &vsc) != eslOK) esl_fatal(msg);   // use <rxf> matrix temporarily
      if ( h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc) != eslOK) esl_fatal(msg);   // overwrites Viterbi matrix with Forward in <rxf>
      if ( h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxd,      NULL) != eslOK) esl_fatal(msg);   // use <rxd> matrix temporarily
      if ( h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxd, rxd)  != eslOK) esl_fatal(msg);   // overwrites Backward with Decoding in <rxd>

      /* Determine most probable anchor set, followed by ASC decoding */
      if ( h4_reference_MPAS(rng, sq->dsq, sq->n, hmm, mo, rxf, rxd, pi, &wrk, ah, afu, afd, anch, &asc, NULL, NULL) != eslOK) esl_fatal(msg);
      if ( h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch, apu, apd, NULL)               != eslOK) esl_fatal(msg);
      if ( h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(msg);

      /* Envelope inference */
      if ( h4_reference_Envelopes(sq->dsq, sq->n, hmm, mo, anch, apu, apd, afu, afd, env) != eslOK) esl_fatal(msg);

      /* Test 1. Coords for each domain are coherent. */
      if (anch->D != env->D) esl_fatal(msg);
      for (d = 1; d <= anch->D; d++)
        {
          if (! (            1 <= env->e[d].oa &&
                  env->e[d].oa <= env->e[d].ia &&
                  env->e[d].ia <= env->e[d].i0 &&
                  env->e[d].i0 <= env->e[d].ib &&
                  env->e[d].ib <= env->e[d].ob &&
                  env->e[d].ob <= sq->n))
            esl_fatal(msg);
        }

      /* Test 2. Envelopes do not overlap (in ia/ib coords. Outer oa/ob coords can.) */
      for (d = 1; d <= env->D; d++)
        if (! (env->e[d].ia > env->e[d-1].ib)) esl_fatal(msg);  // boundary condition e[0].ib = 0

      /* Test 3. envsc(d) <= asc_sc <= fwdsc */
      for (d = 1; d <= env->D; d++)
        if (! (env->e[d].env_sc <= asc+tol && asc <= fsc+tol)) esl_fatal(msg);

      /* Test 4. For D=1 case (only), if outer envelope encompasses
       *         generated trace's domain, env_sc > trace score.
       */
      if (h4_path_GetDomainCount(gpi) == 1 && anch->D == 1)
        {
          h4_path_FetchDomainBounds(gpi, 1, &ia, &ib, /*ka=*/NULL, /*kb=*/NULL);
          if (env->e[1].oa <= ia && env->e[1].ob >= ib)
            {
              if (! (env->e[1].env_sc >= gsc)) esl_fatal(msg);
            }
        }

      // sq, gpi are reuse'd in h4_emit
      // pi, rxf, rxd are reuse'd in h4_reference_* DP routines
      // ah, afu, afd, anch reuse'd in MPAS
      // apu, apd reuse'd in asc_Backward
      // env reuse'd in Envelopes
    } 

  free(wrk);
  h4_envset_Destroy(env);
  h4_anchorhash_Destroy(ah);
  h4_anchorset_Destroy(anch);
  h4_refmx_Destroy(afu);  h4_refmx_Destroy(afd);    
  h4_refmx_Destroy(rxf);  h4_refmx_Destroy(rxd);    // apu,apd point at these and reuse them.
  h4_path_Destroy(pi);    h4_path_Destroy(gpi);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
}
#endif // h4REFERENCE_ENVELOPES_TESTDRIVE

/*****************************************************************
 * 4. Test driver
 *****************************************************************/

#ifdef h4REFERENCE_ENVELOPES_TESTDRIVE

#include "easel.h"
#include "esl_getopts.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                   docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show brief help summary",                   0 },
  { "-s",         eslARG_INT,     "0", NULL, NULL,  NULL,   NULL,  NULL, "set random number generator seed",          0 },
  { "-M",         eslARG_INT,    "10", NULL, NULL,  NULL,   NULL,  NULL, "set test profile length",                   0 },
  { "-N",         eslARG_INT,   "100", NULL, NULL,  NULL,   NULL,  NULL, "number of profile/seq comparisons to test", 0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show HMMER version number",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = h4_CreateDefaultApp(options, 0, argc, argv, "test driver for reference_envelopes", "[-options]");
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc = esl_alphabet_Create(eslAMINO);
  int             M   = esl_opt_GetInteger(go, "-M");
  int             N   = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng)); 

  utest_generation(rng, M, abc, N);

  fprintf(stderr, "#  status   = ok\n");

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif // h4REFERENCE_ENVELOPES_TESTDRIVE


/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef h4REFERENCE_ENVELOPES_EXAMPLE
#include <h4_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_refmx.h"

#include "general.h"
#include "reference_asc.h"
#include "reference_dp.h"
#include "reference_envelopes.h"
#include "reference_mpas.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show HMMER version info",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of reference envelope determination";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  H4_REFMX       *rxf     = NULL;
  H4_REFMX       *rxd     = NULL;
  H4_REFMX       *afu     = NULL;
  H4_REFMX       *afd     = NULL;
  H4_REFMX       *apu     = NULL;
  H4_REFMX       *apd     = NULL;
  H4_PATH        *pi      = h4_path_Create();
  H4_ANCHORSET   *anch    = h4_anchorset_Create(0,0,0);
  H4_ANCHORHASH  *ah      = h4_anchorhash_Create();
  H4_ENVSET      *env     = h4_envset_Create(0,0,0);
  float          *wrk     = NULL;
  float           fsc, vsc, asc;
  int             status;

  /* Read one profile */
  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);

  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);
 
  /* Read one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      esl_fatal("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  /* Set target length model */
  h4_mode_SetLength(mo, sq->n);

  /* Allocate DP matrices */
  rxf = h4_refmx_Create(hmm->M, sq->n);
  rxd = h4_refmx_Create(hmm->M, sq->n);
  afu = h4_refmx_Create(hmm->M, sq->n);
  afd = h4_refmx_Create(hmm->M, sq->n);
  apu = h4_refmx_Create(hmm->M, sq->n);
  apd = h4_refmx_Create(hmm->M, sq->n);

  /* First pass analysis */
  h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxf, pi,  &vsc);
  h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc);
  h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxd,      NULL);   
  h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxd, rxd);   

  /* Determine most probable anchor set */
  h4_reference_MPAS(rng, sq->dsq, sq->n, hmm, mo, rxf, rxd, pi, &wrk, ah,
                    afu, afd, anch, &asc, NULL, NULL);

  /* ASC decoding */
  h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch, apu, apd, NULL);
  h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, apu, apd, apu, apd);

  /* Envelope determination */
  h4_reference_Envelopes(sq->dsq, sq->n, hmm, mo, anch, apu, apd, afu, afd, env);

  h4_envset_Dump(stdout, env);

  h4_refmx_Destroy(rxf);   h4_refmx_Destroy(rxd);
  h4_refmx_Destroy(afu);   h4_refmx_Destroy(afd);
  h4_refmx_Destroy(apu);   h4_refmx_Destroy(apd);
  h4_envset_Destroy(env);
  h4_path_Destroy(pi);
  h4_anchorset_Destroy(anch);
  h4_anchorhash_Destroy(ah);
  free(wrk);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  h4_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif // h4REFERENCE_ENVELOPES_EXAMPLE
