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
static int is_multidomain_region  (P7_DOMAINDEF *ddef, int i, int j);
static int region_trace_ensemble  (P7_DOMAINDEF *ddef, P7_PROFILE *gm, const ESL_SQ *sq, int i, int j, P7_GMX *gx, int *ret_nc);
static int rescore_isolated_domain(P7_DOMAINDEF *ddef, P7_PROFILE *gm, const ESL_SQ *sq, P7_GMX *gx1, P7_GMX *gx2, 
				   int i, int j, int noverlap,
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
  ESL_ALLOC(ddef->dcl, sizeof(P7_DOMAIN) * nalloc);
  ddef->nalloc = nalloc;
  ddef->ndom   = 0;

  ddef->nexpected  = 0.0;
  ddef->nregions   = 0;
  ddef->nclustered = 0;
  ddef->noverlaps  = 0;
  ddef->nenvelopes = 0;

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



/* Function:  p7_domaindef_Fetch()
 * Synopsis:  Retrieve one domain from a <P7_DOMAINDEF> object.
 * Incept:    SRE, Fri Jan 25 13:37:18 2008 [Janelia]
 *
 * Purpose:   Retrieve domain number <which> from domain definitions in <ddef>.
 *            <which> is numbered <0..ndomains-1>. The number of domains
 *            is in <ddef->ndom>.
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
p7_domaindef_Fetch(P7_DOMAINDEF *ddef, int which, int *opt_i, int *opt_j, float *opt_sc, P7_ALIDISPLAY **opt_ad)
{
  if (opt_i  != NULL) *opt_i  = ddef->dcl[which].i;
  if (opt_j  != NULL) *opt_j  = ddef->dcl[which].j;
  if (opt_sc != NULL) *opt_sc = ddef->dcl[which].sc;
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

  for (d = 0; d < ddef->ndom; d++) {
    p7_alidisplay_Destroy(ddef->dcl[d].ad);
    ddef->dcl[d].ad = NULL;
  }
  ddef->ndom = 0;
  ddef->L    = 0;

  ddef->nexpected  = 0.0;
  ddef->nregions   = 0;
  ddef->nclustered = 0;
  ddef->noverlaps  = 0;
  ddef->nenvelopes = 0;

  p7_spensemble_Reuse(ddef->sp);
  p7_trace_Reuse(ddef->tr);	/* probable overkill; should already have been called */
  p7_trace_Reuse(ddef->gtr);	/* likewise */
  return eslOK;
}


/* Function:  p7_domaindef_DumpPosteriors()
 * Synopsis:  Output posteriors that define domain structure to a stream.
 * Incept:    SRE, Fri Feb 29 08:32:14 2008 [Janelia]
 *
 * Purpose:   Output the posterior distributions from <ddef> that are used to
 *            define domain structure to a stream <ofp>, in xmgrace format.
 *            
 *            There are three distributions. The first set is
 *            <mocc[1..i..L]>, the probability that residue <i> is
 *            emitted by the core model (is in a domain). The second
 *            set is <btot[1..i..L]>, the cumulative expected number
 *            of times that a domain uses a B state (starts) at or
 *            before position <i>. The third set is <etot[1..i..L]>,
 *            the cumulative expected number of times that a domain
 *            uses an E state (ends) at or before position <i>. 
 *            
 *            These three fields will only be available after a call
 *            to domain definition by
 *            <p7_domaindef_ByPosteriorHeuristics()>.
 *
 * Returns:   <eslOK> on success
 *            
 * Xref:      J2/126
 */
int
p7_domaindef_DumpPosteriors(FILE *ofp, P7_DOMAINDEF *ddef)
{
  int i;

  for (i = 1; i <= ddef->L; i++)
    fprintf(ofp, "%d %f\n", i, ddef->mocc[i]);
  fprintf(ofp, "&\n");

  for (i = 1; i <= ddef->L; i++)
    fprintf(ofp, "%d %f\n", i, ddef->btot[i]);
  fprintf(ofp, "&\n");

  for (i = 1; i <= ddef->L; i++)
    fprintf(ofp, "%d %f\n", i, ddef->etot[i]);
  fprintf(ofp, "&\n");

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
      p7_alidisplay_Destroy(ddef->dcl[d].ad);
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
p7_domaindef_ByViterbi(P7_PROFILE *gm, const ESL_SQ *sq, P7_GMX *gx1, P7_GMX *gx2, P7_DOMAINDEF *ddef) 
{
  P7_ALIDISPLAY *ad;
  float          dsc;
  int            d;
  int            saveL = gm->L;

  p7_domaindef_GrowTo(ddef, sq->n);
  p7_GTrace  (sq->dsq, sq->n, gm, gx1, ddef->gtr);
  p7_trace_Index(ddef->gtr);
  
  p7_ReconfigUnihit(gm, 0);	  /* process each domain in unihit L=0 mode */

  for (d = 0; d < ddef->gtr->ndom; d++)
    rescore_isolated_domain(ddef, gm, sq, gx1, gx2, ddef->gtr->sqfrom[d], ddef->gtr->sqto[d], 0, &ad, &dsc);
  p7_ReconfigMultihit(gm, saveL); /* restore original profile configuration */
  return eslOK;
}


/* Function:  p7_domaindef_ByPosteriorHeuristics()
 * Synopsis:  Define domains in a sequence using posterior probs.
 * Incept:    SRE, Sat Feb 23 08:17:44 2008 [Janelia]
 *
 * Purpose:   Given a sequence <sq> and model <gm> for which we have 
 *            already calculated a Forward and Backward matrix
 *            <fwd> and <bck>; use posterior probability heuristics
 *            to determine an annotated domain structure. Caller
 *            provides a new or reused <ddef> object to hold these
 *            results.
 *            
 *            Upon return, <ddef> contains the definitions of all the
 *            domains: their bounds, their null-corrected Forward
 *            scores, and their optimal posterior accuracy alignments.
 *            Caller can cycle over <ddef->ndom> domains and retrieve
 *            info on each one by <p7_domain_Fetch()>.
 *            
 * Returns:   <eslOK> on success.           
 */
int
p7_domaindef_ByPosteriorHeuristics(P7_PROFILE *gm, const ESL_SQ *sq, P7_GMX *fwd, P7_GMX *bck, P7_DOMAINDEF *ddef)
{
  P7_ALIDISPLAY *ad;
  float          dsc;
  int i, j;
  int triggered;
  int d;
  int i2,j2;
  int last_j2;
  int saveL = gm->L;
  int nc;
  int noverlap;

  p7_domaindef_GrowTo(ddef, sq->n);                /* now ddef's btot,etot,mocc arrays are ready for seq of length n */
  calculate_domain_posteriors(gm, fwd, bck, ddef); /* now ddef->{btot,etot,mocc} are made.                           */
  ddef->nexpected = ddef->btot[sq->n];             /* posterior expectation for # of domains (same as etot[sq->n])   */

  /* That's all we need the original fwd, bck matrices for. 
   * Now we're free to reuse them for subsequent domain calculations.
   */

  p7_ReconfigUnihit(gm, saveL);	                   /* process each domain in unilocal mode                           */
  i     = -1;
  triggered = FALSE;
  for (j = 1; j <= sq->n; j++)
    {
      if (! triggered) 
	{			/* xref J2/101 for what the logic below is: */
	  if       (ddef->mocc[j] - (ddef->btot[j] - ddef->btot[j-1]) <  ddef->rt2) i = j;
	  else if  (i == -1)                                                        i = j;
	  if       (ddef->mocc[j]                                     >= ddef->rt1) triggered = TRUE;  
	} 
      else if (ddef->mocc[j] - (ddef->etot[j] - ddef->etot[j-1])  <  ddef->rt2) 
	{
	  /* We have a region i..j to evaluate. */
	  ddef->nregions++;

	  if (is_multidomain_region(ddef, i, j))
	    {
	      /* This region appears to contain more than one domain, so we have to 
               * resolve it by cluster analysis of posterior trace samples, to define
               * one or more domain envelopes.
	       */
	      ddef->nclustered++;

	      p7_ReconfigMultihit(gm, saveL);
	      region_trace_ensemble(ddef, gm, sq, i, j, fwd, &nc);
	      p7_ReconfigUnihit(gm, saveL);
	      last_j2 = 0;
	      for (d = 0; d < nc; d++) {
		p7_spensemble_GetClusterCoords(ddef->sp, d, &i2, &j2, NULL, NULL, NULL);
		if (i2 <= last_j2) ddef->noverlaps++;
		last_j2 = j2;

		noverlap = (last_j2 >= i ? last_j2-i+1 : 0);

		ddef->nenvelopes++;
		rescore_isolated_domain(ddef, gm, sq, fwd, bck, i2, j2, noverlap, &ad, &dsc);
	      }
	      p7_spensemble_Reuse(ddef->sp);
	      p7_trace_Reuse(ddef->tr);
	    }
	  else 
	    {
	      /* The region looks simple, single domain; convert the region to an envelope. */
	      ddef->nenvelopes++;
	      rescore_isolated_domain(ddef, gm, sq, fwd, bck, i, j, 0, &ad, &dsc);
	    }
	  i     = -1;
	  triggered = FALSE;
	}
    }

  p7_ReconfigMultihit(gm, saveL); /* restore original profile configuration */
  return eslOK;
}



/*****************************************************************
 * 3. Internal routines 
 *****************************************************************/

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
      ddef->btot[i] = ddef->btot[i-1] + exp(fwd->xmx[(i-1)*p7G_NXCELLS+p7G_B] + bck->xmx[(i-1)*p7G_NXCELLS+p7G_B] - overall_logp);
      ddef->etot[i] = ddef->etot[i-1] + exp(fwd->xmx[i    *p7G_NXCELLS+p7G_E] + bck->xmx[i    *p7G_NXCELLS+p7G_E] - overall_logp);

      njcp  = expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_N] + bck->xmx[p7G_NXCELLS*i + p7G_N] + gm->xsc[p7P_N][p7P_LOOP] - overall_logp);
      njcp += expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_J] + bck->xmx[p7G_NXCELLS*i + p7G_J] + gm->xsc[p7P_J][p7P_LOOP] - overall_logp);
      njcp += expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_C] + bck->xmx[p7G_NXCELLS*i + p7G_C] + gm->xsc[p7P_C][p7P_LOOP] - overall_logp);
      ddef->mocc[i] = 1. - njcp;
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
 *   E(z) = expected number of E states occurring in region before z is emitted
 *        = \sum_{y=i}^{z} eocc[i]  =  etot[z] - etot[i-1]
 *   B(z) = expected number of B states occurring in region after z is emitted
 *        = \sum_{y=z}^{j} bocc[i]  =  btot[j] - btot[z-1]               
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
      expected_n = ESL_MIN( (ddef->etot[z] - ddef->etot[i-1]), (ddef->btot[j] - ddef->btot[z-1]));
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
 * The caller provides model <gm> configured in multihit mode, with
 * its target length distribution set to <sq->n>: i.e., the model
 * configuration used to score the complete sequence (if it weren't
 * multihit, we wouldn't be worried about multiple domains).
 * 
 * Caller provides <ddef>, which defines heuristic parameters that
 * control the clustering, and provides working space for the
 * calculation and the answers. The <ddef->sp> object must have been
 * reused (i.e., it needs to be fresh; we're going to use it here);
 * the caller needs to Reuse() it specifically, because it can't just
 * Reuse() the whole <ddef>, when it's in the process of analyzing
 * regions.
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
 *    (the cluster definitions). When the caller is done retrieving those
 *    answers, it needs to <esl_spensemble_Reuse()> it before calling
 *    <region_trace_ensemble()> again.
 *    
 * <ddef->tr> is used as working memory for sampled traces.
 *    
 * <gm> gets reconfigured into multihit <L> mode!
 * 
 * <gx> is used for a Forward matrix on the domain <i..j>.
 * 
 */
static int
region_trace_ensemble(P7_DOMAINDEF *ddef, P7_PROFILE *gm, const ESL_SQ *sq, int i, int j, P7_GMX *gx, int *ret_nc)
{
  float sc;
  int   t, d;

  p7_GForward(sq->dsq+i-1, j-i+1, gm, gx, &sc);

  for (t = 0; t < ddef->nsamples; t++)
    {
      p7_StochasticTrace(ddef->r, sq->dsq+i-1, j-i+1, gm, gx, ddef->tr);
      p7_trace_Index(ddef->tr);

      for (d = 0; d < ddef->tr->ndom; d++)
	p7_spensemble_Add(ddef->sp, t, ddef->tr->sqfrom[d]+i-1, ddef->tr->sqto[d]+i-1, ddef->tr->hmmfrom[d], ddef->tr->hmmto[d]);

      p7_trace_Reuse(ddef->tr);        
    }
  p7_spensemble_Cluster(ddef->sp, ddef->min_overlap, ddef->of_smaller, ddef->max_diagdiff, ddef->min_posterior, ddef->min_endpointp, ret_nc);
  return eslOK;
}


/* rescore_isolated_domain()
 * SRE, Fri Feb  8 09:18:33 2008 [Janelia]
 *
 * We have isolated a single domain's envelope from <i>..<j> in
 * sequence <sq>, and now we want to score it in isolation and obtain
 * an alignment display for it.
 * 
 * (Later, we can add up all the individual domain scores from this
 * seq into a new per-seq score, to compare to the original per-seq
 * score).
 *  
 * The caller provides model <gm> configured in unilocal mode; by
 * using unilocal (as opposed to multilocal), we're going to force the
 * identification of a single domain in this envelope now.
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
 * The caller also passes us <noverlap>, the number of residues in
 * this domain envelope that overlap with a previously defined
 * envelope. (Overlapping envelopes may arise when stochastic
 * traceback clustering is used to define them.) We use <noverlap>
 * to make sure that per-seq null2 correction doesn't doublecount
 * any residues.
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
 * <gx1> : happens to be holding OA score matrix for the domain
 *         upon return, but that's not part of the spec; officially
 *         its contents are "undefined".
 *
 * <gx2> : happens to be holding a posterior probability matrix
 *         for the domain upon return, but we're not making that
 *         part of the spec, so caller shouldn't rely on this;
 *         spec just makes its contents "undefined".
 */
static int
rescore_isolated_domain(P7_DOMAINDEF *ddef, P7_PROFILE *gm, const ESL_SQ *sq, 
			P7_GMX *gx1, P7_GMX *gx2, int i, int j, int noverlap,
			P7_ALIDISPLAY **ret_ad, float *ret_sc)
{
  P7_DOMAIN     *dom = NULL;
  int            Ld = j-i+1;
  float          envsc, oasc;
  int            z;
  float          domcorrection, seqcorrection;


  p7_GForward (sq->dsq + i-1, Ld, gm, gx1, &envsc);
  p7_GBackward(sq->dsq + i-1, Ld, gm, gx2, NULL);

  p7_Null2Corrections(gm, sq->dsq + i-1, Ld, noverlap, gx1, gx2, &domcorrection, &seqcorrection);

  /* Find an optimal accuracy alignment */
  p7_PosteriorDecoding(Ld, gm, gx1, gx2, gx2);   /* <gx2> is now overwritten with post probabilities     */
  p7_OptimalAccuracyDP(Ld, gm, gx2, gx1, &oasc); /* <gx1> is now overwritten with OA scores              */
  p7_OATrace(Ld, gm, gx2, gx1, ddef->tr);        /* <tr>'s seq coords are offset by i-1, rel to orig dsq */
  
  /* hack the trace's sq coords to be correct w.r.t. original dsq */
  for (z = 0; z < ddef->tr->N; z++)
    if (ddef->tr->i[z] > 0) ddef->tr->i[z] += i-1;

  /* get ptr to next empty domain structure in domaindef's results */
  if (ddef->ndom == ddef->nalloc) {
    void *p;
    ESL_RALLOC(ddef->dcl, p, sizeof(P7_DOMAIN) * (ddef->nalloc*2));
    ddef->nalloc *= 2;    
  }
  dom = &(ddef->dcl[ddef->ndom]);

  /* store the results in it */
  dom->ienv          = i;
  dom->jenv          = j;
  dom->envsc         = fsc;
  dom->seqcorrection = seqcorrection;
  dom->domcorrection = domcorrection;
  dom->ad            = p7_alidisplay_Create(ddef->tr, 0, gm, sq);
  ddef->ndom++;

  p7_trace_Reuse(ddef->tr);
  *ret_ad = ad;
  *ret_sc = fsc;
  return eslOK;
}
  
    
/*****************************************************************
 * Example driver.
 *****************************************************************/

#ifdef p7DOMAINDEF_EXAMPLE

/* gcc -o domaindef_example -g -Wall -I../easel -L../easel -I. -L. -Dp7DOMAINDEF_EXAMPLE p7_domaindef.c -lhmmer -leasel -lm
 * ./domaindef_example <hmmfile> <seqfile>
 */


#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_alphabet.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "--occp",    eslARG_OUTFILE, NULL, NULL, NULL,  NULL,  NULL, NULL, "output posterior occupancies for xmgrace to <f>",  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of domain definition by posterior sampling";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *r       = esl_randomness_Create(42);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  char           *seqfile = esl_opt_GetArg(go, 2);
  int             format  = eslSQFILE_UNKNOWN;
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_ALPHABET   *abc     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *fwd     = NULL;
  P7_GMX         *bck     = NULL;
  P7_DOMAINDEF   *ddef    = NULL;
  P7_ALIDISPLAY  *ad      = NULL;
  char           *ofile   = NULL;
  FILE           *ofp     = NULL;
  float           overall_sc, sc;
  int             d;
  int             i,j;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
  if  (esl_sqio_Read(sqfp, sq) != eslOK) p7_Fail("Failed to read sequence");
  esl_sqfile_Close(sqfp);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);

  /* allocate the domain definition object */
  ddef = p7_domaindef_Create(r);

  /* allocate DP matrices for forward and backward */
  fwd = p7_gmx_Create(gm->M, sq->n);
  bck = p7_gmx_Create(gm->M, sq->n);

  /* run Forward and Backward */
  p7_FLogsumInit();
  p7_GForward (sq->dsq, sq->n, gm, fwd, &overall_sc); 
  p7_GBackward(sq->dsq, sq->n, gm, bck, &sc); 

  printf("Overall raw likelihood score: %.2f nats\n", overall_sc);

  /* define domain structure */
  p7_domaindef_ByPosteriorHeuristics(gm, sq, fwd, bck, ddef);

  /* retrieve and display results */
  for (d = 0; d < ddef->ndom; d++)
    {
      p7_domaindef_Fetch(ddef, d, &i, &j, &sc, &ad);
      printf("domain %-4d : %4d %4d  %6.2f\n", d+1, i, j, sc);

      p7_alidisplay_Print(stdout, ad, 50, 120);
      p7_alidisplay_Destroy(ad);
    }

  if ((ofile = esl_opt_GetString(go, "--occp")) != NULL) 
    {
      if ((ofp = fopen(ofile, "w")) == NULL) p7_Fail("Failed to open output file %s\n", ofile);
      p7_domaindef_DumpPosteriors(ofp, ddef);
      fclose(ofp);
    }

  p7_domaindef_Destroy(ddef);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7DOMAINDEF_EXAMPLE*/


#ifdef p7DOMAINDEF_EXAMPLE2
/* gcc -o domaindef_example2 -g -Wall -I../easel -L../easel -I. -L. -Dp7DOMAINDEF_EXAMPLE2 p7_domaindef.c -lhmmer -leasel -lm
 * ./domaindef_example2 <hmmfile> 
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_alphabet.h"
#include "esl_stopwatch.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-N",        eslARG_INT,   "1000", NULL, NULL,  NULL,  NULL, NULL, "number of sampled sequences",                      0 },
  { "-L",        eslARG_INT,    "400", NULL, NULL,  NULL,  NULL, NULL, "length config for the profile",                    0 },
  { "-V",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "find domains by Viterbi parsing",                  0 },
  { "-b",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "baseline timing - no domain processing",           0 },
  { "-v",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "be more verbose (show per sequence results)",      0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of domain definition by posterior sampling";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_RANDOMNESS *r       = esl_randomness_Create(42);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  char           *hmmfile = esl_opt_GetArg(go, 1);
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_ALPHABET   *abc     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_TRACE       *tr      = NULL;
  P7_GMX         *fwd     = NULL;
  P7_GMX         *bck     = NULL;
  P7_DOMAINDEF   *ddef    = NULL;
  int   N           = esl_opt_GetInteger(go, "-N");
  int   L0          = esl_opt_GetInteger(go, "-L");
  int   do_vit      = esl_opt_GetBoolean(go, "-V");
  int   do_baseline = esl_opt_GetBoolean(go, "-b");
  int   be_verbose  = esl_opt_GetBoolean(go, "-v");
  float           overall_sc, sc;
  int             idx;
  int             tot_true, tot_found;

  /* Read in one HMM */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L0);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L0, p7_LOCAL);

  /* Other initial allocations */
  sq   = esl_sq_CreateDigital(abc);
  ddef = p7_domaindef_Create(r);
  fwd  = p7_gmx_Create(gm->M, L0);
  bck  = p7_gmx_Create(gm->M, L0);
  tr   = p7_trace_Create();
  p7_FLogsumInit();

  tot_true = tot_found = 0;
  esl_stopwatch_Start(w);
  for (idx = 0; idx < N; idx++)
    {
      p7_ReconfigLength(gm, L0);
      p7_bg_SetLength(bg, L0);
      p7_ProfileEmit(r, hmm, gm, bg, sq, tr); /* sample a sequence from the profile */
      p7_trace_Index(tr);                      /* tr->ndom is the "true" domain number emitted */
      tot_true += tr->ndom;

      p7_ReconfigLength(gm, sq->n);
      p7_bg_SetLength(bg, sq->n);
      p7_gmx_GrowTo(fwd, gm->M, sq->n);
      p7_gmx_GrowTo(bck, gm->M, sq->n);

      if (do_vit) 
	{
	  p7_GViterbi (sq->dsq, sq->n, gm, fwd, &overall_sc); 
	  if (! do_baseline)
	    p7_domaindef_ByViterbi(gm, sq, fwd, bck, ddef);
	}
      else
	{
	  p7_GForward (sq->dsq, sq->n, gm, fwd, &overall_sc); 
	  if (! do_baseline) {
	    p7_GBackward(sq->dsq, sq->n, gm, bck, &sc);       
	    p7_domaindef_ByPosteriorHeuristics(gm, sq, fwd, bck, ddef);
	  }
	}

      tot_found += ddef->ndom;
      if (be_verbose) 
	printf("true: %d   found: %d\n", tr->ndom, ddef->ndom);

      p7_trace_Reuse(tr);
      p7_domaindef_Reuse(ddef);
    }
  esl_stopwatch_Stop(w);

  printf("Total domains: %d\n", tot_true);
  printf("Total found:   %d\n", tot_found);
  printf("Accuracy:      %.2f%%\n", 100. * (double) tot_found / (double) tot_true);
  esl_stopwatch_Display(stdout, w, "CPU time:   ");

  p7_domaindef_Destroy(ddef);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_profile_Destroy(gm);
  p7_trace_Destroy(tr);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_sq_Destroy(sq);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7DOMAINDEF_EXAMPLE2*/

