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
#include "h4_config.h"

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
 *            afu  : ASC Forward UP matrix
 *            afd  : ASC Forward DOWN matrix
 *            apu  : ASC posterior decoding UP matrix
 *            apd  : ASC posterior decoding DOWN matrix
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
      }
  return eslOK;
}

/*****************************************************************
 * 3. Unit tests
 *****************************************************************/

// TK


/*****************************************************************
 * 4. Test driver
 *****************************************************************/

// TK


/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef h4REFERENCE_ENVELOPES_EXAMPLE
#include "h4_config.h"

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
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, -1, argc, argv, banner, usage);  // -1 means allow any number of cmdline args
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
