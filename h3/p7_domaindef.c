/* Definition of multidomain structure of a target sequence, and
 * rescoring as a sum of individual domains, with null2 correction.
 * 
 * Contents:
 *    1. The P7_DOMAINDEF object: allocation, reuse, destruction
 *    2. Routines inferring domain structure of a target sequence
 *    3. Internal routines 
 *    
 *
 * Exegesis:
 * 
 * The key function here is <p7_domaindef_ByPosteriorHeuristics()>.
 * Everything else is support structure. 
 * 
 * When you call <p7_domaindef_ByPosteriorHeuristics()>, you have a
 * per-sequence hit that's judged significant, and you've calculated
 * Forward/Backward score matrices for the complete sequence.  Thus,
 * the input data are the model <gm>, the sequence <sq>, and filled-in
 * forward and backward matrices <fwd>, <bck>.
 * 
 * The function then chews over this data, using posterior
 * probabilities and heuristics to define, score, and obtain
 * display-ready alignments for individual domains. When it's done,
 * your <fwd> and <bck> matrices have been effectively destroyed (they
 * get reused for individual domain alignment calculations), and
 * <ddef> contains all the per-domain results you need. It returns to
 * you the number of domains it's defined (in <ret_ndom>), and the
 * total per-sequence score derived by a sum of individual domain
 * scores (in <ret_sc>).
 * 
 * You can then use <p7_domaindef_Fetch()> to retrieve the start
 * position, end position, score, and display-ready alignment for each
 * domain, one at a time. 
 * 
 * The <P7_DOMAINDEF> structure is a reusable container that manages
 * all the necessary working memory and heuristic thresholds.
 *   
 * SRE, Thu Jan 24 09:28:01 2008 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <math.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_sqio.h"

#include "hmmer.h"

static int calculate_domain_posteriors(P7_PROFILE *gm, P7_GMX *fwd, P7_GMX *bck, P7_DOMAINDEF *ddef);
static int is_multidomain_region(P7_DOMAINDEF *ddef, int i, int j);
static int rescore_isolated_domain(P7_DOMAINDEF *ddef, P7_PROFILE *gm, const ESL_SQ *sq, P7_GMX *gx1, P7_GMX *gx2, int i, int j, 
				   P7_ALIDISPLAY **ret_ad, float *ret_sc);



/*****************************************************************
 * 1. The P7_DOMAINDEF object: allocation, reuse, destruction
 *****************************************************************/

/* Function:  p7_domaindef_Create()
 * Synopsis:  Creates a new <P7_DOMAINDEF> object.
 * Incept:    SRE, Fri Jan 25 13:21:31 2008 [Janelia]
 *
 * Purpose:   Creates a new <P7_DOMAINDEF> object, with <r> registered
 *            as its random number generator, using default settings
 *            for all thresholds.
 *
 * Returns:   a pointer to the new <P7_DOMAINDEF> object.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_DOMAINDEF *
p7_domaindef_Create(ESL_RANDOMNESS *r)
{
  P7_DOMAINDEF *ddef   = NULL;
  int           Lalloc = 512;	/* this initial alloc doesn't matter much; space is realloced as needed */
  int           nalloc = 32;
  int           status;

  /* level 1 alloc */
  ESL_ALLOC(ddef, sizeof(P7_DOMAINDEF));
  ddef->mocc = ddef->btot = ddef->etot = NULL;
  ddef->sp  = NULL;
  ddef->tr  = NULL;
  ddef->dcl = NULL;

  /* level 2 alloc: posterior prob arrays */
  ESL_ALLOC(ddef->mocc, sizeof(float) * (Lalloc+1));
  ESL_ALLOC(ddef->btot, sizeof(float) * (Lalloc+1));
  ESL_ALLOC(ddef->etot, sizeof(float) * (Lalloc+1));
  ddef->mocc[0] = ddef->etot[0] = ddef->btot[0] = 0.;
  ddef->Lalloc  = Lalloc;
  ddef->L       = 0;

  /* level 2 alloc: results storage */
  ESL_ALLOC(ddef->dcl, sizeof(struct p7_dom_s) * nalloc);
  ddef->nalloc = nalloc;
  ddef->ndom   = 0;

  /* default thresholds */
  ddef->rt1           = 0.25;
  ddef->rt2           = 0.10;
  ddef->rt3           = 0.20;
  ddef->nsamples      = 200;
  ddef->min_overlap   = 0.8;
  ddef->of_smaller    = TRUE;
  ddef->max_diagdiff  = 4;
  ddef->min_posterior = 0.25;
  ddef->min_endpointp = 0.02;

  /* allocate reusable, growable objects that domain def reuses for each seq */
  ddef->sp  = p7_spensemble_Create(1024, 64, 32); /* init allocs = # sampled pairs; max endpoint range; # of domains */
  ddef->tr  = p7_trace_Create();
  ddef->gtr = p7_trace_Create();

  /* keep a copy of ptr to the RNG */
  ddef->r = r;  
  return ddef;
  
 ERROR:
  p7_domaindef_Destroy(ddef);
  return NULL;
}


/* p7_domaindef_GrowTo()
 * Synopsis:  Reallocates a <P7_DOMAINDEF> for new seq length <L>
 * Incept:    SRE, Fri Jan 25 13:27:24 2008 [Janelia]
 *
 * Purpose:   Reallocates a <P7_DOMAINDEF> object <ddef> so that
 *            it can hold a sequence of up to <L> residues. 
 *
 *            (This might be a no-op, if <ddef> is already large
 *            enough.)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. In this case, the
 *            data in <ddef> are unaffected.
 */
static int
p7_domaindef_GrowTo(P7_DOMAINDEF *ddef, int L)
{
  void *p;
  int   status;

  if (L <= ddef->Lalloc) return eslOK;

  ESL_RALLOC(ddef->mocc, p, sizeof(float) * (L+1));
  ESL_RALLOC(ddef->btot, p, sizeof(float) * (L+1));
  ESL_RALLOC(ddef->etot, p, sizeof(float) * (L+1));
  ddef->Lalloc = L;
  return eslOK;

 ERROR:
  return status;
}

/* p7_domaindef_Add()
 * Synopsis:  Define a new domain.
 * Incept:    SRE, Fri Jan 25 13:31:12 2008 [Janelia]
 *
 * Purpose:   Define a new domain with start position <i>, end position <j>,
 *            score <sc>, and display alignment <ad>; store this information
 *            in <ddef> for future retrieval by the caller.
 *            
 *            <i>,<j> are in the <1..L> coordinate frame of the
 *            original <dsq> being analyzed. Likewise for the
 *            coordinates in the <ad> (and indeed, <ad->sqfrom = i>,
 *            <ad->sqto = j>).
 *            
 *            Once a non-<NULL> <ad> is passed into the <P7_DOMAINDEF>
 *            object, the object becomes responsible for <ad>'s
 *            memory. 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error. In this case, the domain
 *            is not added; the original data in <ddef> are unchanged.
 */
static int
p7_domaindef_Add(P7_DOMAINDEF *ddef, int i, int j, float sc, P7_ALIDISPLAY *ad)
{
  int   status;

  if (ddef->ndom == ddef->nalloc) {
    void *p;
    ESL_RALLOC(ddef->dcl, p, sizeof(struct p7_dom_s) * (ddef->nalloc*2));
    ddef->nalloc *= 2;    
  }

  ddef->dcl[ddef->n].i  = i;
  ddef->dcl[ddef->n].j  = j;
  ddef->dcl[ddef->n].sc = sc;
  ddef->dcl[ddef->n].ad = ad;
  ddef->ndom++;
  return eslOK;
  
 ERROR:
  return status;
}


/* Function:  p7_domaindef_Fetch()
 * Synopsis:  Retrieve one domain from a <P7_DOMAINDEF> object.
 * Incept:    SRE, Fri Jan 25 13:37:18 2008 [Janelia]
 *
 * Purpose:   Retrieve domain number <which> from domain definitions in <ddef>.
 *            <which> is numbered <0..ndomains-1>. The number of domains
 *            was returned by <p7_domaindef_ByPosteriorHeuristics()> (or 
 *            some other domain-defining function); alternatively, it may
 *            also be accessed as <ddef->ndom>.
 *            
 *            The information is retrieved in <opt_i> (start position
 *            of domain, <1..L>), <opt_j> (end position), <opt_sc>
 *            (score), and <opt_ad> (a pointer to the displayable
 *            alignment). All of these are optional; pass <NULL> for
 *            any that you aren't interested in retrieving.
 *            
 *            When you retrieve a <opt_ad> pointer, you become
 *            responsible for its memory; the <P7_DOMAINDEF> releases
 *            responsibility to you, and deletes its hold on the
 *            <ad>. You free alidisplays with
 *            <p7_alidisplay_Destroy()>. If you retrieve a domain
 *            coords and/or score without its <ad> -- i.e.\ passing
 *            <NULL> for <opt_ad> -- the <ad> is internally destroyed.
 *            
 *            A consequence of this memory management strategy for the
 *            <ad> is that although you may retrieve domain coords and
 *            scores more than once from the same <ddef>, you may only
 *            retrieve the <ad> for domain <which> once. If you need
 *            the <ad> for domain <which>, make sure you retrieve it
 *            the first time you ask for info on domain <which>.
 *            
 * Returns:   <eslOK> on success.
 */
int
p7_domaindef_Fetch(P7_DOMAINDEF *ddef, int which, int *opt_i, int *opt_j, float *opt_sc, P7_ALIDISPLAY *opt_ad)
{
  if (opt_i  != NULL) *opt_i  = ddef->dcl[which].i;
  if (opt_j  != NULL) *opt_j  = ddef->dcl[which].j;
  if (opt_sc != NULL) *opt_j  = ddef->dcl[which].sc;
  if (opt_ad != NULL) *opt_ad = ddef->dcl[which].ad; else { p7_alidisplay_Destroy(ddef->dcl[which].ad); }
  ddef->dcl[which].ad = NULL;
  return eslOK;
}


/* Function:  p7_domaindef_Reuse()
 * Synopsis:  Prepare to reuse a <P7_DOMAINDEF> on a new sequence.
 * Incept:    SRE, Fri Jan 25 13:48:36 2008 [Janelia]
 *
 * Purpose:   Prepare a <P7_DOMAINDEF> object <ddef> to be reused on
 *            a new sequence, reusing as much memory as possible.
 *            
 * Note:      Because of the way we handle alidisplays, handing them off to
 *            the caller, we don't reuse their memory; any unused
 *            alidisplays are destroyed. It's not really possible to
 *            reuse alidisplay memory. We need alidisplays to persist
 *            until all sequences have been processed and we're
 *            writing our final output to the user.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_domaindef_Reuse(P7_DOMAINDEF *ddef)
{
  int d;

  ddef->ndom = 0;
  ddef->L    = 0;
  for (d = 0; d < ddef->ndom; d++) {
    p7_alidisplay_Destroy(ddef->dcl[d].ad);
    ddef->dcl[d].ad = NULL;
  }

  p7_spensemble_Reuse(ddef->sp);
  p7_trace_Reuse(ddef->tr);	/* probable overkill; should already have been called */
  p7_trace_Reuse(ddef->gtr);	/* likewise */
  return eslOK;
}


/* Function:  p7_domaindef_Destroy()
 * Synopsis:  Destroys a <P7_DOMAINDEF>.
 * Incept:    SRE, Fri Jan 25 13:52:46 2008 [Janelia]
 *
 * Purpose:   Destroys a <P7_DOMAINDEF>.
 */
void
p7_domaindef_Destroy(P7_DOMAINDEF *ddef)
{
  int d;

  if (ddef == NULL) return;

  if (ddef->mocc != NULL) free(ddef->mocc);
  if (ddef->btot != NULL) free(ddef->btot);
  if (ddef->etot != NULL) free(ddef->etot);

  if (ddef->dcl  != NULL) {
    for (d = 0; d < ddef->ndom; d++)
      p7_alidisplay_Destroy(ddef->dcl[d]->ad);
    free(ddef->dcl);
  }

  p7_spensemble_Destroy(ddef->sp);
  p7_trace_Destroy(ddef->tr);
  p7_trace_Destroy(ddef->gtr);
  free(ddef);
  return;
}

/*****************************************************************
 * 2. Routines inferring domain structure of a target sequence
 *****************************************************************/

/* Function:  p7_domaindef_ByViterbi()
 * Synopsis:  Define domains in a sequence by maximum likelihood.
 * Incept:    SRE, Fri Jan 25 15:10:21 2008 [Janelia]
 *
 * Purpose:   Use a Viterbi (maximum likelihood) parse to determine
 *            the domain structure of sequence <sq> aligned to 
 *            model <gm>. Caller provides a filled Viterbi matrix
 *            in <gx1>, and a second matrix of at least the same
 *            size for scratch space in <gx2>.
 *            
 *            Upon return, <ddef> contains definitions of all the
 *            domains, bounds defined by Viterbi parse, individually
 *            scored by null2-corrected Forward, and aligned by
 *            optimal posterior accuracy. Caller can cycle over
 *            <ddef->ndom> domains and retrieve info on each one
 *            by <p7_domaindef_Fetch()>.
 *            
 * Returns:   <eslOK> on success.
 */
int
p7_domaindef_ByViterbi(P7_PROFILE *gm, ESL_SQ *sq, P7_GMX *gx1, P7_GMX *gx2, P7_DOMAINDEF *ddef) 
{
  P7_ALIDISPLAY *ad;
  float          dsc;
  float          vsc;
  int            d;
  int            status;

  p7_domaindef_GrowTo(ddef, sq->n);
  p7_GTrace  (sq->dsq, sq->n, gm, gx1, ddef->gtr);
  p7_trace_Index(ddef->gtr);
  
  for (d = 0; d < ddef->gtr->ndom; d++)
    {
      rescore_isolated_domain(ddef, gm, sq, gx1, gx2, ddef->gtr->sqfrom[d], ddef->gtr->sqto[d], &ad, &dsc);
      p7_domaindef_Add(ddef, ddef->gtr->sqfrom[d], ddef->gtr->sqto[d], dsc, ad);
    }
  return eslOK;
}


int
p7_domaindef_ByPosteriorHeuristics(P7_PROFILE *gm, ESL_SQ *sq, P7_GMX *fwd, P7_GMX *bck, P7_DOMAINDEF *ddef)
{
  P7_ALIDISPLAY *ad;
  float          dsc;
  int i;
  int lasti, j;
  int triggered;
  int d;
  int i2,j2;
  float prob;

  p7_domaindef_GrowTo(ddef, sq->n);                /* now ddef's btot,etot,mocc arrays are ready for seq of length n */
  calculate_domain_posteriors(gm, fwd, bck, ddef); /* now ddef->{btot,etot,mocc} are made.                           */

  /* That's all we need the original fwd, bck matrices for. 
   * Now we're free to reuse them for subsequent domain calculations.
   */

  lasti     = -1;
  triggered = FALSE;
  for (i = 1; i <= sq->n; i++)
    {
      if (! triggered) 
	{
	  if       (ddef->mocc[i] - (ddef->btot[i] - ddef->btot[i-1]) <  ddef->rt2) lasti = -1;
	  else if  (lasti == -1)                        lasti = i;
	  if       (ddef->mocc[i]                                     >= ddef->rt1) triggered = TRUE;  
	} 
      else if (ddef->mocc[i] - (ddef->etot[i] - ddef->etot[i-1])  <  ddef->rt2) 
	{

	  if (is_multidomain_region(ddef, i, j))
	    {
	      region_trace_ensemble(ddef, gm, sq, fwd, &nc);
	      for (d = 0; d < nc; d++) {
		p7_spensemble_GetClusterCoords(ddef->sp, d, &i2, &j2, NULL, NULL, NULL);
		rescore_isolated_domain(ddef, gm, sq, fwd, bck, i2, j2, &ad, &dsc);
		p7_domaindef_Add(ddef, i2, j2, dsc, ad);
	      }
	    }
	  else 
	    {
	      rescore_isolated_domain(ddef, gm, sq, fwd, bck, lasti, i, &ad, &dsc);
	      p7_domaindef_Add(ddef, lasti, i, dsc, ad);
	    }

	  lasti     = -1;
	  triggered = FALSE;
	}
    }
  return eslOK;
}



/*****************************************************************
 * 3. Internal routines 
 ***************************************************************** 

/* calculate_domain_posteriors()
 * SRE, Fri Feb  8 10:44:30 2008 [Janelia]
 *
 * The caller has calculated Forward and Backward matrices <fwd> and
 * <bck> for model <gm> aligned to a target sequence. (The target
 * sequence doesn't need to be provided, because all we need to know
 * is its length <L>, and that's available in either of the two DP
 * matrices.) 
 * 
 * We use this information to calculate the posterior probabilities
 * that we're in a begin state B, end state E, or any core model state
 * {M,D,I} at each target sequence position <i = 1..L>.
 * 
 * This information is stored in three arrays in <ddef>. This routine
 * expects that this storage has already been (re)allocated
 * appropriately for a target seq of length <L>.
 * 
 * <ddef->btot[i]> stores the cumulative expectation \sum_1^i of the
 * number of i's that were emitted (by an Mk state) immediately after
 * a B : i.e., the expected number of times domains have started at or
 * before position i.
 * 
 * <ddef->etot[i]> stores the cumulative expectation \sum_1^i of the
 * number of i's that were emitted (by an Mk or Dk state) and
 * immediately followed by an end transition to E : i.e., the expected
 * number of times domains have ended at or before position i.
 * 
 * <ddef->mocc[i]> stores the probability that residue i is emitted by
 * the core model, as opposed to the flanking N,C,J states : i.e., the
 * probability that i is in a domain.
 * 
 * Upon return, each of these arrays has been made, and <ddef->L> has
 * been set.
 * 
 * Ideas for future optimization:
 * 
 * - The calculations only need to access the xmx[CNJBE][i] special
 *   states in the fwd, bck matrices, so we could use streamlined
 *   (checkpointed?) matrices that only maintain this info over all
 *   i. This would be a step in letting us do domain parses in linear
 *   memory.
 *   
 * - indeed, the <btot>, <etot>, and <mocc> arrays could be made
 *   sparse; on long target sequences, we expect long stretches of
 *   negligible posterior probability that we're in the model or using
 *   a begin or end transition.
 *   
 * - indeed indeed, we don't really need to store the <btot>, <etot>,
 *   and <mocc> arrays at all. We can define regions in a single pass
 *   in O(1) extra memory, straight from the <fwd>, <bck> matrices, if
 *   we have to (xref J2/101). <p7_domaindef_ByPosteriorHeuristics()>
 *   is already implemented in a way to make this easy. We're not
 *   doing that for now, partly for clarity in the code, and partly
 *   because I think we'll want to output the <btot>, <etot>, and
 *   <mocc> arrays -- this view of the posterior decoding of the
 *   domain structure of a target sequence will be useful. Also, it's
 *   a lot easier to implement the <is_multidomain_region()> trigger
 *   if these arrays are available.
 */
static int
calculate_domain_posteriors(P7_PROFILE *gm, P7_GMX *fwd, P7_GMX *bck, P7_DOMAINDEF *ddef)
{
  int   L            = fwd->L;
  float overall_logp = fwd->xmx[p7G_NXCELLS*L + p7G_C] + gm->xsc[p7P_C][p7P_MOVE];
  float njcp;
  int   i;

  for (i = 1; i <= L; i++)
    {
      ddef->btot[i] = d->btot[i-1] + exp(fwd->xmx[(i-1)*p7G_NXCELLS+p7G_B] + bck->xmx[(i-1)*p7G_NXCELLS+p7G_B] - overall_logp);
      ddef->etot[i] = d->etot[i-1] + exp(fwd->xmx[i    *p7G_NXCELLS+p7G_E] + bck->xmx[i    *p7G_NXCELLS+p7G_E] - overall_logp);

      njc           = expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_N] + bck->xmx[p7G_NXCELLS*i + p7G_N] + gm->xsc[p7P_N][p7P_LOOP] - overall_logp);
      njc          += expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_J] + bck->xmx[p7G_NXCELLS*i + p7G_J] + gm->xsc[p7P_J][p7P_LOOP] - overall_logp);
      njc          += expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_C] + bck->xmx[p7G_NXCELLS*i + p7G_C] + gm->xsc[p7P_C][p7P_LOOP] - overall_logp);
      ddef->mocc[i] = 1. - njc;
    }
  ddef->L = gm->L;
  return eslOK;
}


/* is_multidomain_region()
 * SRE, Fri Feb  8 11:35:04 2008 [Janelia]
 *
 * This defines the trigger for when we need to hand a "region" off to
 * a deeper analysis (using stochastic tracebacks and clustering)
 * because there's reason to suspect it may encompass two or more
 * domains. 
 * 
 * The criterion is to find the split point z at which the expected
 * number of E occurrences preceding B occurrences is maximized, and
 * if that number is greater than the heuristic threshold <ddef->rt3>,
 * then return TRUE. In other words, we're checking to see if there's
 * any point in the region at which it looks like an E was followed by
 * a B, as expected for a multidomain interpretation of the region.
 * 
 * More precisely: return TRUE if  \max_z [ \min (B(z), E(z)) ]  >= rt3
 * where
 *   E(z) = expected number of E states occurring in region after z is emitted
 *        = \sum_{y=z}^{j} eocc[i]  =  etot[j] - etot[z-1]
 *
 *   B(z) = expected number of B states occurring in region before z is emitted
 *        = \sum_{y=i}^{z} bocc[i]  =  btot[z] - btot[i-1]               
 *        
 *        
 * Because this relies on the <ddef->etot> and <ddef->btot> arrays,
 * <calculate_domain_posteriors()> needs to have been called first.
 *
 * Xref:    J2/101.  
 */
static int
is_multidomain_region(P7_DOMAINDEF *ddef, int i, int j)
{
  int   z;
  float max;
  float expected_n;

  max = -1.0;
  for (z = i; z <= j; z++)
    {
      expected_n = ESL_MIN( (ddef->etot[j] - ddef->etot[z-1]), (ddef->btot[z] - ddef->btot[i-1]));
      max        = ESL_MAX(max, expected_n);
    }

  return ( (max >= ddef->rt3) ? TRUE : FALSE);
}


/* region_trace_ensemble()
 * SRE, Fri Feb  8 11:49:44 2008 [Janelia]
 *
 * Here, we've decided that region <i>..<j> in sequence <sq> might be
 * composed of more than one domain, and we're going to use clustering
 * of a posterior ensemble of stochastic tracebacks to sort it out.
 * 
 * To collect stochastic tracebacks, we need a Forward matrix
 * calculated over the region. Wastefully (see section on future
 * optimizations below), we calculate the Forward matrix here.  This
 * means we need the caller to provide the model <gm> and allocated
 * space for the Forward DP matrix in <gx>. 
 * 
 * The model will be reconfigured into multihit mode, with its target
 * length distribution set to <sq->n>: i.e., presumably identical to
 * the model configuration used to score the complete sequence (if it
 * weren't multihit, we wouldn't be worried about multiple domains).
 * (This also is potentially wasteful - see future optimizations notes
 * - and it's also a fracked up API, because we're changing the model
 * config and potentially confusing the caller.)
 * 
 * Caller also provides <ddef>, which defines heuristic parameters
 * that control the clustering, and provides working space for the
 * calculation and the answers.
 * 
 * Upon return, <*ret_nc> contains the number of clusters that were
 * defined.
 * 
 * The caller can retrieve info on each cluster by calling
 * <p7_spensemble_GetClusterCoords(ddef->sp...)> on the
 * <P7_SPENSEMBLE> object in <ddef>.
 * 
 * Other information on what's happened in working memory:
 * 
 * <ddef->sp> gets filled in, and upon return, it's holding the answers 
 *    (the cluster definitions). After the caller retrieves those answers,
 *    it needs to 
 * 
 * 
 * 
 * 
 * 
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */

static int
region_trace_ensemble(P7_DOMAINDEF *ddef, P7_PROFILE *gm, ESL_SQ *sq, int i, int j, P7_GMX *gx, int *ret_nc)
{
  float sc;
  int   t, d;

  p7_ReconfigLength(gm, L);
  p7_ReconfigMultihit(gm);
  p7_GForward(dsq, L, gm, gx, &sc);

  for (t = 0; t < ddef->nsamples; t++)
    {
      p7_StochasticTrace(ddef->r, dsq, L, gm, gx, ddef->tr);
      p7_trace_Index(ddef->tr);

      for (d = 0; d < ddef->tr->ndom; d++)
	p7_spensemble_Add(ddef->sp, t, ddef->tr->sqfrom[d], ddef->tr->sqto[d], ddef->tr->hmmfrom[d], ddef->tr->hmmto[d]);

      p7_trace_Reuse(ddef->tr);        
    }
  p7_spensemble_Cluster(ddef->sp, ddef->min_overlap, ddef->of_smaller, ddef->max_diagdiff, ddef->min_posterior, ddef->min_endpointp, ret_nc);
  return eslOK;
}


/* rescore_isolated_domain()
 * SRE, Fri Feb  8 09:18:33 2008 [Janelia]
 *
 * We have isolated a single domain from <i>..<j> in sequence <sq>,
 * and we want to score it in isolation and obtain an alignment
 * display for it.
 * 
 * (Later, we can add up all the individual domain scores from this
 * seq into a new per-seq score, to compare to the original per-seq
 * score).
 *  
 * The score is a Forward score for model <gm> configured in unihit
 * lobal mode (local wrt model, global wrt sequence). Lobal mode
 * forces the complete domain sequence from <i> to <j> into one
 * consistent alignment against the model.
 * 
 * The alignment is an optimal accuracy alignment (sensu IH Holmes),
 * also obtained in lobal mode.
 * 
 * The caller provides DP matrices <gx1> and <gx2> with sufficient
 * space to hold Forward and Backward calculations for this domain
 * against the model. (The caller will typically already have matrices
 * sufficient for the complete sequence lying around, and can just use
 * those.) The caller also provides a <P7_DOMAINDEF> object which is
 * (efficiently, we trust) managing any necessary temporary working
 * space and heuristic thresholds.
 * 
 * It isn't implemented yet, but this is eventually where null2
 * correction will go, since we obtain a posterior probability matrix
 * here; then the returned score will become a null2-corrected Forward
 * score.
 * 
 * Upon return, the <P7_ALIDISPLAY> <*ret_ad> has been newly created,
 * and the caller becomes responsible for freeing it. 
 * 
 * And here's what's happened to our working memory:
 * 
 * <ddef>: <ddef->tr> has been used, and possibly reallocated, for
 *         the OA trace of the domain. Before exit, we called
 *         <Reuse()> on it.
 * 
 * <gm>  : we've reconfigured it to lobal mode, and it's still in
 *         that mode upon return! watch out, caller.
 *            
 * <gx1> : happens to be holding a posterior probability matrix
 *         for the domain upon return, but we're not making that
 *         part of the spec, so caller shouldn't rely on this;
 *         spec just makes its contents "undefined".
 *         
 * <gx2> : happens to be holding OA score matrix for the domain
 *         upon return, but that's not part of the spec; officially
 *         its contents are "undefined".
 */
static int
rescore_isolated_domain(P7_DOMAINDEF *ddef, P7_PROFILE *gm, const ESL_SQ *sq, P7_GMX *gx1, P7_GMX *gx2, int i, int j, 
			P7_ALIDISPLAY **ret_ad, float *ret_sc)
{
  P7_ALIDISPLAY *ad = NULL;
  int             L = j-i+1;
  float           fsc, sc;

  p7_ReconfigLength(gm, 0);
  p7_ReconfigUnihit(gm);

  p7_GForward (sq->dsq + i, L, gm, gx1, &fsc);
  p7_GBackward(sq->dsq + i, L, gm, gx2, &sc);

  p7_PostProb(L, gm, gx1, gx2, fwd);          /* <gx1> is now overwritten with posterior probabilities */

  /* This is the info we need for a post-prob-oriented null2 correction,
   * and this is where the null2 correction call will go in the near future
   */

  p7_OptimalAccuracyDP(L, gm, gx1, gx2, &sc); /* <gx2> is now overwritten with OA scores */
  p7_OATrace(L, gm, gx1, gx2, ddef->tr);      /* <tr>'s seq coords are all offset by i-1 now, relative to original dsq  */
  
  ad = p7_alidisplay_Create(ddef->tr, 0, gm, sq);
  ad->sqfrom += i-1;
  ad->sqto   += i-1;

  *ret_ad = ad;
  *ret_sc = fsc;
  return eslOK;
}
  
    
