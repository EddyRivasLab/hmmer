/* P7_TOPHITS: implementation of ranked list of top-scoring hits
 * 
 * Contents:
 *    1. The P7_TOPHITS object.
 *    2. Standard (human-readable) output of pipeline results.
 *    3. Tabular (parsable) output of pipeline results.
 *    4. Benchmark driver.
 *    5. Test driver.
 *    6. Copyright and license information.
 * 
 * SRE, Fri Dec 28 07:14:54 2007 [Janelia] [Enigma, MCMXC a.D.]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "easel.h"

#include "hmmer.h"


/*****************************************************************
 * 1. The P7_TOPHITS object
 *****************************************************************/

/* Function:  p7_tophits_Create()
 * Synopsis:  Allocate a hit list.
 * Incept:    SRE, Fri Dec 28 07:17:51 2007 [Janelia]
 *
 * Purpose:   Allocates a new <P7_TOPHITS> hit list and return a pointer
 *            to it.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_TOPHITS *
p7_tophits_Create(void)
{
  P7_TOPHITS *h = NULL;
  int         default_nalloc = 256;
  int         status;

  ESL_ALLOC(h, sizeof(P7_TOPHITS));
  h->hit    = NULL;
  h->unsrt  = NULL;

  ESL_ALLOC(h->hit,   sizeof(P7_HIT *) * default_nalloc);
  ESL_ALLOC(h->unsrt, sizeof(P7_HIT)   * default_nalloc);
  h->Nalloc    = default_nalloc;
  h->N         = 0;
  h->nreported = 0;
  h->nincluded = 0;
  h->is_sorted = TRUE;       	/* but only because there's 0 hits */
  h->hit[0]    = h->unsrt;	/* if you're going to call it "sorted" when it contains just one hit, you need this */
  return h;

 ERROR:
  p7_tophits_Destroy(h);
  return NULL;
}


/* Function:  p7_tophits_Grow()
 * Synopsis:  Reallocates a larger hit list, if needed.
 * Incept:    SRE, Fri Dec 28 07:37:27 2007 [Janelia]
 *
 * Purpose:   If list <h> cannot hold another hit, doubles
 *            the internal allocation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; in this case,
 *            the data in <h> are unchanged.
 */
int
p7_tophits_Grow(P7_TOPHITS *h)
{
  void   *p;
  P7_HIT *ori    = h->unsrt;
  int     Nalloc = h->Nalloc * 2;	/* grow by doubling */
  int     i;
  int     status;

  if (h->N < h->Nalloc) return eslOK; /* we have enough room for another hit */

  ESL_RALLOC(h->hit,   p, sizeof(P7_HIT *) * Nalloc);
  ESL_RALLOC(h->unsrt, p, sizeof(P7_HIT)   * Nalloc);

  /* If we grow a sorted list, we have to translate the pointers
   * in h->hit, because h->unsrt might have just moved in memory. 
   */
  if (h->is_sorted) 
    {
      for (i = 0; i < h->N; i++)
	h->hit[i] = h->unsrt + (h->hit[i] - ori);
    }

  h->Nalloc = Nalloc;
  return eslOK;

 ERROR:
  return eslEMEM;
}


/* Function:  p7_tophits_CreateNextHit()
 * Synopsis:  Get pointer to new structure for recording a hit.
 * Incept:    SRE, Tue Mar 11 08:44:53 2008 [Janelia]
 *
 * Purpose:   Ask the top hits object <h> to do any necessary
 *            internal allocation and bookkeeping to add a new,
 *            empty hit to its list; return a pointer to 
 *            this new <P7_HIT> structure for data to be filled
 *            in by the caller.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_tophits_CreateNextHit(P7_TOPHITS *h, P7_HIT **ret_hit)
{
  P7_HIT *hit = NULL;
  int     status;

  if ((status = p7_tophits_Grow(h)) != eslOK) goto ERROR;
  
  hit = &(h->unsrt[h->N]);
  h->N++;
  if (h->N >= 2) h->is_sorted = FALSE;

  hit->name         = NULL;
  hit->acc          = NULL;
  hit->desc         = NULL;
  hit->sortkey      = 0.0;

  hit->score        = 0.0;
  hit->pre_score    = 0.0;
  hit->sum_score    = 0.0;

  hit->pvalue       = 0.0;
  hit->pre_pvalue   = 0.0;
  hit->sum_pvalue   = 0.0;

  hit->ndom         = 0;
  hit->nexpected    = 0.0;
  hit->nregions     = 0;
  hit->nclustered   = 0;
  hit->noverlaps    = 0;
  hit->nenvelopes   = 0;

  hit->flags        = p7_HITFLAGS_DEFAULT;
  hit->nreported    = 0;
  hit->nincluded    = 0;
  hit->best_domain  = -1;
  hit->dcl          = NULL;

  *ret_hit = hit;
  return eslOK;

 ERROR:
  *ret_hit = NULL;
  return status;
}



/* Function:  p7_tophits_Add()
 * Synopsis:  Add a hit to the top hits list.
 * Incept:    SRE, Fri Dec 28 08:26:11 2007 [Janelia]
 *
 * Purpose:   Adds a hit to the top hits list <h>. 
 * 
 *            <name>, <acc>, and <desc> are copied, so caller may free
 *            them if it likes.
 *            
 *            Only the pointer <ali> is kept. Caller turns over memory
 *            management of <ali> to the top hits object; <ali> will
 *            be free'd when the top hits structure is free'd.
 *
 * Args:      h        - active top hit list
 *            name     - name of target  
 *            acc      - accession of target (may be NULL)
 *            desc     - description of target (may be NULL) 
 *            sortkey  - value to sort by: bigger is better
 *            score    - score of this hit
 *            pvalue   - P-value of this hit 
 *            mothersc - score of parent whole sequence 
 *            motherp  - P-value of parent whole sequence
 *            sqfrom   - 1..L pos in target seq  of start
 *            sqto     - 1..L pos; sqfrom > sqto if rev comp
 *            sqlen    - length of sequence, L
 *            hmmfrom  - 0..M+1 pos in HMM of start
 *            hmmto    - 0..M+1 pos in HMM of end
 *            hmmlen   - length of HMM, M
 *            domidx   - number of this domain 
 *            ndom     - total # of domains in sequence
 *            ali      - optional printable alignment info
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> if reallocation failed.
 * 
 * Note:      Is this actually used anywhere? (SRE, 10 Dec 08) 
 *            I think it's not up to date.
 */
int
p7_tophits_Add(P7_TOPHITS *h,
	       char *name, char *acc, char *desc, 
	       double sortkey, 
	       float score,    double pvalue, 
	       float mothersc, double motherp,
	       int sqfrom, int sqto, int sqlen,
	       int hmmfrom, int hmmto, int hmmlen, 
	       int domidx, int ndom,
	       P7_ALIDISPLAY *ali)
{
  int status;

  if ((status = p7_tophits_Grow(h))                           != eslOK) return status;
  if ((status = esl_strdup(name, -1, &(h->unsrt[h->N].name))) != eslOK) return status;
  if ((status = esl_strdup(acc,  -1, &(h->unsrt[h->N].acc)))  != eslOK) return status;
  if ((status = esl_strdup(desc, -1, &(h->unsrt[h->N].desc))) != eslOK) return status;
  h->unsrt[h->N].sortkey  = sortkey;
  h->unsrt[h->N].score    = score;
  h->unsrt[h->N].pvalue   = pvalue;
  h->N++;

  if (h->N >= 2) h->is_sorted = FALSE;
  return eslOK;
}

/* hit_sorter(): qsort's pawn, below */
static int
hit_sorter(const void *vh1, const void *vh2)
{
  P7_HIT *h1 = *((P7_HIT **) vh1);  /* don't ask. don't change. Don't Panic. */
  P7_HIT *h2 = *((P7_HIT **) vh2);

  if      (h1->sortkey < h2->sortkey)  return  1;
  else if (h1->sortkey > h2->sortkey)  return -1;
  else                                 return  0;
}

/* Function:  p7_tophits_Sort()
 * Synopsis:  Sorts a hit list.
 * Incept:    SRE, Fri Dec 28 07:51:56 2007 [Janelia]
 *
 * Purpose:   Sorts a top hit list. After this call,
 *            <h->hit[i]> points to the i'th ranked 
 *            <P7_HIT> for all <h->N> hits.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_Sort(P7_TOPHITS *h)
{
  int i;

  if (h->is_sorted)  return eslOK;
  for (i = 0; i < h->N; i++) h->hit[i] = h->unsrt + i;
  if (h->N > 1)  qsort(h->hit, h->N, sizeof(P7_HIT *), hit_sorter);
  h->is_sorted = TRUE;
  return eslOK;
}

/* Function:  p7_tophits_Merge()
 * Synopsis:  Merge two top hits lists.
 * Incept:    SRE, Fri Dec 28 09:32:12 2007 [Janelia]
 *
 * Purpose:   Merge <h2> into <h1>. Upon return, <h1>
 *            contains the sorted, merged list. <h2>
 *            is effectively destroyed; caller should
 *            not access it further, and may as well free
 *            it immediately.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure, and
 *            both <h1> and <h2> remain valid.
 */
int
p7_tophits_Merge(P7_TOPHITS *h1, P7_TOPHITS *h2)
{
  void    *p;
  P7_HIT **new_hit = NULL;
  P7_HIT  *ori1    = h1->unsrt;	/* original base of h1's data */
  P7_HIT  *new2;
  int      i,j,k;
  int      Nalloc = h1->Nalloc + h2->Nalloc;
  int      status;

  /* Make sure the two lists are sorted */
  if ((status = p7_tophits_Sort(h1)) != eslOK) goto ERROR;
  if ((status = p7_tophits_Sort(h2)) != eslOK) goto ERROR;

  /* Attempt our allocations, so we fail early if we fail. 
   * Reallocating h1->unsrt screws up h1->hit, so fix it.
   */
  ESL_RALLOC(h1->unsrt, p, sizeof(P7_HIT) * Nalloc);
  ESL_ALLOC (new_hit, sizeof(P7_HIT *)    * Nalloc);
  for (i = 0; i < h1->N; i++)
    h1->hit[i] = h1->unsrt + (h1->hit[i] - ori1);

  /* Append h2's unsorted data array to h1. h2's data begin at <new2> */
  new2 = h1->unsrt + h1->N;
  memcpy(new2, h2->unsrt, sizeof(P7_HIT) * h2->N);

  /* Merge the sorted hit lists */
  for (i=0,j=0,k=0; i < h1->N && j < h2->N ; k++)
    new_hit[k] = (h2->hit[j]->sortkey > h1->hit[i]->sortkey) ? new2 + (h2->hit[j++] - h2->unsrt) : h1->hit[i++];
  while (i < h1->N) new_hit[k++] = h1->hit[i++];
  while (j < h2->N) new_hit[k++] = new2 + (h2->hit[j++] - h2->unsrt);

  /* h2 now turns over management of name, acc, desc memory to h1;
   * nullify its pointers, to prevent double free.  */
  for (i = 0; i < h2->N; i++)
    {
      h2->unsrt[i].name = NULL;
      h2->unsrt[i].acc  = NULL;
      h2->unsrt[i].desc = NULL;
    }

  /* Construct the new grown h1 */
  free(h1->hit);
  h1->hit    = new_hit;
  h1->Nalloc = Nalloc;
  h1->N     += h2->N;
  /* and is_sorted is TRUE, as a side effect of p7_tophits_Sort() above. */
  return eslOK;	
  
 ERROR:
  if (new_hit != NULL) free(new_hit);
  return status;
}


/* Function:  p7_tophits_GetMaxNameLength()
 * Synopsis:  Returns maximum name length.
 * Incept:    SRE, Fri Dec 28 09:00:13 2007 [Janelia]
 *
 * Purpose:   Returns the maximum name length of all the registered
 *            hits, in chars. This is useful when deciding how to
 *            format output.
 *            
 *            The maximum is taken over all registered hits. This
 *            opens a possible side effect: caller might print only
 *            the top hits, and the max name length in these top hits
 *            may be different than the max length over all the hits.
 *            
 *            If there are no hits in <h>, or none of the
 *            hits have names, returns 0.
 */
int
p7_tophits_GetMaxNameLength(P7_TOPHITS *h)
{
  int i, max, n;
  for (max = 0, i = 0; i < h->N; i++)
    if (h->unsrt[i].name != NULL) {
      n   = strlen(h->unsrt[i].name);
      max = ESL_MAX(n, max);
    }
  return max;
}

/* Function:  p7_tophits_Reuse()
 * Synopsis:  Reuse a hit list, freeing internals.
 * Incept:    SRE, Fri Jun  6 15:39:05 2008 [Janelia]
 *
 * Purpose:   Reuse the tophits list <h>; save as 
 *            many malloc/free cycles as possible,
 *            as opposed to <Destroy()>'ing it and
 *            <Create>'ing a new one.
 */
int
p7_tophits_Reuse(P7_TOPHITS *h)
{
  int i, j;

  if (h == NULL) return eslOK;
  if (h->unsrt != NULL) 
    {
      for (i = 0; i < h->N; i++)
	{
	  if (h->unsrt[i].name != NULL) free(h->unsrt[i].name);
	  if (h->unsrt[i].acc  != NULL) free(h->unsrt[i].acc);
	  if (h->unsrt[i].desc != NULL) free(h->unsrt[i].desc);
	  if (h->unsrt[i].dcl  != NULL) {
	    for (j = 0; j < h->unsrt[i].ndom; j++)
	      if (h->unsrt[i].dcl[j].ad != NULL) p7_alidisplay_Destroy(h->unsrt[i].dcl[j].ad);
	    free(h->unsrt[i].dcl);
	  }
	}
    }
  h->N         = 0;
  h->is_sorted = TRUE;
  h->hit[0]    = h->unsrt;
  return eslOK;
}

/* Function:  p7_tophits_Destroy()
 * Synopsis:  Frees a hit list.
 * Incept:    SRE, Fri Dec 28 07:33:21 2007 [Janelia]
 */
void
p7_tophits_Destroy(P7_TOPHITS *h)
{
  int i,j;
  if (h == NULL) return;
  if (h->hit   != NULL) free(h->hit);
  if (h->unsrt != NULL) 
    {
      for (i = 0; i < h->N; i++)
	{
	  if (h->unsrt[i].name != NULL) free(h->unsrt[i].name);
	  if (h->unsrt[i].acc  != NULL) free(h->unsrt[i].acc);
	  if (h->unsrt[i].desc != NULL) free(h->unsrt[i].desc);
	  if (h->unsrt[i].dcl  != NULL) {
	    for (j = 0; j < h->unsrt[i].ndom; j++)
	      if (h->unsrt[i].dcl[j].ad != NULL) p7_alidisplay_Destroy(h->unsrt[i].dcl[j].ad);
	    free(h->unsrt[i].dcl);
	  }
	}
      free(h->unsrt);
    }
  free(h);
  return;
}
/*---------------- end, P7_TOPHITS object -----------------------*/






/*****************************************************************
 * 2. Standard (human-readable) output of pipeline results
 *****************************************************************/

/* Function:  p7_tophits_Threshold()
 * Synopsis:  Apply score and E-value thresholds to a hitlist before output.
 * Incept:    SRE, Tue Dec  9 09:04:55 2008 [Janelia]
 *
 * Purpose:   After a pipeline has completed, go through it and mark all
 *            the targets and domains that are "significant" (satisfying
 *            the reporting thresholds set for the pipeline). 
 *            
 *            Also sets the final total number of reported targets
 *            <th->nreported> and the size of the search space for
 *            per-domain conditional E-value calculations,
 *            <pli->domZ>. By default, <pli->domZ> is the number of
 *            significant targets reported.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_Threshold(P7_TOPHITS *th, P7_PIPELINE *pli)
{
  int h, d;	/* counters over sequence hits, domains in sequences */
  
  /* First pass over sequences: flag all reportable, includable ones */
  th->nreported = 0;
  for (h = 0; h < th->N; h++)
    {
      if (p7_pli_TargetReportable(pli, th->hit[h]->score, th->hit[h]->pvalue))
	{
	  th->hit[h]->flags |= p7_IS_REPORTED;
	  th->nreported++;
	}
      if (p7_pli_TargetIncludable(pli, th->hit[h]->score, th->hit[h]->pvalue))
	{
	  th->hit[h]->flags |= p7_IS_INCLUDED;
	  th->nincluded++;
	}
    }
  
  /* Now we can determined domZ, the effective search space in which additional domains are found */
  if (pli->domZ_setby == p7_ZSETBY_NTARGETS) pli->domZ = (double) th->nreported;

  /* Second pass is over domains, flagging reportable/includable ones:
   * we always report or include at least the best single domain for each sequence, 
   * regardless of domain thresholding.
   */
  for (h = 0; h < th->N; h++)  
    {
      if (th->hit[h]->flags & p7_IS_REPORTED)
	for (d = 0; d < th->hit[h]->ndom; d++)
	  {
	    if (th->hit[h]->best_domain == d ||
		p7_pli_DomainReportable(pli, th->hit[h]->dcl[d].bitscore, th->hit[h]->dcl[d].pvalue))
	      {
		th->hit[h]->nreported++;
		th->hit[h]->dcl[d].is_reported = TRUE;
	      }
	  }
      if (th->hit[h]->flags & p7_IS_INCLUDED)
	for (d = 0; d < th->hit[h]->ndom; d++)
	  {
	    if (th->hit[h]->best_domain == d ||
		p7_pli_DomainIncludable(pli, th->hit[h]->dcl[d].bitscore, th->hit[h]->dcl[d].pvalue))
	      {
		th->hit[h]->nincluded++;
		th->hit[h]->dcl[d].is_included = TRUE;
	      }
	  }
    }
  return eslOK;
}


/* Function:  p7_tophits_CompareRanking()
 * Synopsis:  Compare current top hits to previous top hits ranking.
 * Incept:    SRE, Sun Feb  8 20:33:02 2009 [Janelia]
 *
 * Purpose:   Using a keyhash <kh> of the previous top hits and the
 *            their ranks, look at the current top hits list <th>
 *            and flag new hits that are included for the first time
 *            (by setting <p7_IS_NEW> flag) and hits that were 
 *            included previously, but are now below the inclusion
 *            threshold in the list (<by setting <p7_IS_DROPPED>
 *            flag). 
 *
 *            The <th> must already have been processed by
 *            <p7_tophits_Threshold()>. We assume the <is_included>,
 *            <is_reported> flags are set on the appropriate hits.
 * 
 *            Upon return, the keyhash <kh> is updated to hash the
 *            current top hits list and their ranks. 
 *            
 *            Optionally, <*opt_nnew> is set to the number of 
 *            newly included hits. jackhmmer uses this as part of
 *            its convergence criteria, for example.
 *            
 *            These flags affect output of top target hits from
 *            <p7_tophits_Targets()>. 
 *            
 *            It only makes sense to call this function in context of
 *            an iterative search.
 *            
 *            The <p7_IS_NEW> flag is comprehensive: all new hits
 *            are flagged (and counted in <*opt_nnew>). The <p7_WAS_DROPPED> 
 *            flag is not comprehensive: only those hits that still 
 *            appear in the current top hits list are flagged. If a 
 *            hit dropped entirely off the list, it isn't counted
 *            as "dropped". (This could be done, but we would want
 *            to have two keyhashes, one old and one new, to do the
 *            necessary comparisons efficiently.)
 *            
 *            If the target names in <th> are not unique, results may
 *            be strange.
 *
 * Args:      th         - current top hits list
 *            kh         - hash of top hits' ranks (in: previous tophits; out: <th>'s tophits)
 *            opt_nnew   - optRETURN: number of new hits above inclusion threshold
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> if <kh> needed to be reallocated but this failed.
 */
int
p7_tophits_CompareRanking(P7_TOPHITS *th, ESL_KEYHASH *kh, int *opt_nnew)
{
  int nnew = 0;
  int oldrank;
  int h;
  int status;

  /* Flag the hits in the list with whether they're new in the included top hits,
   * and whether they've dropped off the included list.
   */
  for (h = 0; h < th->N; h++)
    {
      esl_key_Lookup(kh, th->hit[h]->name, &oldrank);
      
      if (th->hit[h]->flags & p7_IS_INCLUDED) 
	{
	  if (oldrank == -1) { th->hit[h]->flags |= p7_IS_NEW; nnew++; }
	}
      else 
	{
	  if (oldrank >=  0) th->hit[h]->flags |= p7_IS_DROPPED;
	}
    }	

  /* Replace the old rank list with the new one */
  esl_keyhash_Reuse(kh);
  for (h = 0; h < th->N; h++)
    {
      if (th->hit[h]->flags & p7_IS_INCLUDED) 
	{
	  /* What happens when the same sequence name appears twice? It gets stored with higher rank */
	  status = esl_key_Store(kh, th->hit[h]->name, NULL);
	  if (status != eslOK && status != eslEDUP) goto ERROR;
	}
    }
  
  if (opt_nnew != NULL) *opt_nnew = nnew;
  return eslOK;

 ERROR:
  if (opt_nnew != NULL) *opt_nnew = 0;
  return status;
}


/* Function:  p7_tophits_Targets()
 * Synopsis:  Standard output format for a top target hits list.
 * Incept:    SRE, Tue Dec  9 09:10:43 2008 [Janelia]
 *
 * Purpose:   Output a list of the reportable top target hits in <th> 
 *            in human-readable ASCII text format to stream <ofp>, using
 *            final pipeline accounting stored in <pli>. 
 * 
 *            The tophits list <th> should already be sorted (see
 *            <p7_tophits_Sort()> and thresholded (see
 *            <p7_tophits_Threshold>).
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_Targets(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw)
{
  char   newness;
  int    h;
  int    d;
  int    namew = ESL_MAX(8,  p7_tophits_GetMaxNameLength(th));
  int    descw;
  int    have_printed_incthresh = FALSE;

  if (textw >  0) descw = ESL_MAX(32, textw - namew - 61); /* 61 chars excluding desc is from the format: 2 + 22+2 +22+2 +8+2 +<name>+1 */
  else            descw = INT_MAX;

  fprintf(ofp, "Scores for complete sequence%s (score includes all domains):\n", 
	  pli->mode == p7_SEARCH_SEQS ? "s" : "");

  /* The minimum width of the target table is 111 char: 47 from fields, 8 from min name, 32 from min desc, 13 spaces */
  fprintf(ofp, "  %22s  %22s  %8s\n",                              " --- full sequence ---",        " --- best 1 domain ---",   "-#dom-");
  fprintf(ofp, "  %9s %6s %5s  %9s %6s %5s  %5s %2s  %-*s %s\n", "E-value", " score", " bias", "E-value", " score", " bias", "  exp",  "N", namew, (pli->mode == p7_SEARCH_SEQS ? "Sequence":"Model"), "Description");
  fprintf(ofp, "  %9s %6s %5s  %9s %6s %5s  %5s %2s  %-*s %s\n", "-------", "------", "-----", "-------", "------", "-----", " ----", "--", namew, "--------", "-----------");

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)
      {
	d    = th->hit[h]->best_domain;

	if (! (th->hit[h]->flags & p7_IS_INCLUDED) && h < th->N-1 && ! have_printed_incthresh) {
	  fprintf(ofp, "  ------ inclusion threshold ------\n");
	  have_printed_incthresh = TRUE;
	}

	if      (th->hit[h]->flags & p7_IS_NEW)     newness = '+';
	else if (th->hit[h]->flags & p7_IS_DROPPED) newness = '-';
	else                                        newness = ' ';


	fprintf(ofp, "%c %9.2g %6.1f %5.1f  %9.2g %6.1f %5.1f  %5.1f %2d  %-*s %-.*s\n",
		newness,
		th->hit[h]->pvalue * pli->Z,
		th->hit[h]->score,
		th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
		th->hit[h]->dcl[d].pvalue * pli->Z,
		th->hit[h]->dcl[d].bitscore,
		th->hit[h]->dcl[d].dombias,
		th->hit[h]->nexpected,
		th->hit[h]->nreported,
		namew, th->hit[h]->name,
		descw, (th->hit[h]->desc == NULL ? "" : th->hit[h]->desc));
      }
  if (th->nreported == 0) fprintf(ofp, "\n   [No hits detected that satisfy reporting thresholds]\n");
  return eslOK;
}


/* Function:  p7_tophits_Domains()
 * Synopsis:  Standard output format for top domain hits and alignments.
 * Incept:    SRE, Tue Dec  9 09:32:32 2008 [Janelia]
 *
 * Purpose:   For each reportable target sequence, output a tabular summary
 *            of reportable domains found in it, followed by alignments of
 *            each domain.
 * 
 *            Similar to <p7_tophits_Targets()>; see additional notes there.
 */
int
p7_tophits_Domains(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw)
{
  int h, d;
  int nd;
  int namew, descw;

  fprintf(ofp, "Domain and alignment annotation for each %s:\n", pli->mode == p7_SEARCH_SEQS ? "sequence" : "model");

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)
      {
	namew = strlen(th->hit[h]->name);
	descw = (textw > 0 ?  ESL_MAX(32, textw - namew - 5) : INT_MAX);

	fprintf(ofp, ">> %s  %-.*s\n", th->hit[h]->name, descw, (th->hit[h]->desc == NULL ? "" : th->hit[h]->desc));

	/* The domain table is 101 char wide:
          #     score  bias  c-Evalue  i-Evalue hmmfrom   hmmto    alifrom  ali to    envfrom  env to     acc
         ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
           1 ?  123.4  23.1   9.7e-11    6.8e-9       3    1230 ..       1     492 []       2     490 .] 0.90
         123 ! 1234.5 123.4 123456789 123456789 1234567 1234567 .. 1234567 1234567 [] 1234567 1234568 .] 0.12
	*/
	fprintf(ofp, " %3s   %6s %5s %9s %9s %7s %7s %2s %7s %7s %2s %7s %7s %2s %4s\n",    "#",  "score",  "bias",  "c-Evalue",  "i-Evalue", "hmmfrom",  "hmm to", "  ", "alifrom",  "ali to", "  ", "envfrom",  "env to", "  ",  "acc");
	fprintf(ofp, " %3s   %6s %5s %9s %9s %7s %7s %2s %7s %7s %2s %7s %7s %2s %4s\n",  "---", "------", "-----", "---------", "---------", "-------", "-------", "  ", "-------", "-------", "  ", "-------", "-------", "  ", "----");
	nd = 0;
	for (d = 0; d < th->hit[h]->ndom; d++)
	  if (th->hit[h]->dcl[d].is_reported) 
	    {
	      nd++;
	      fprintf(ofp, " %3d %c %6.1f %5.1f %9.2g %9.2g %7d %7d %c%c %7ld %7ld %c%c %7d %7d %c%c %4.2f\n",
		      nd,
		      th->hit[h]->dcl[d].is_included ? '!' : '?',
		      th->hit[h]->dcl[d].bitscore,
		      th->hit[h]->dcl[d].dombias,
		      th->hit[h]->dcl[d].pvalue * pli->domZ,
		      th->hit[h]->dcl[d].pvalue * pli->Z,
		      th->hit[h]->dcl[d].ad->hmmfrom,
		      th->hit[h]->dcl[d].ad->hmmto,
		      (th->hit[h]->dcl[d].ad->hmmfrom == 1) ? '[' : '.',
		      (th->hit[h]->dcl[d].ad->hmmto   == th->hit[h]->dcl[d].ad->M) ? ']' : '.',
		      th->hit[h]->dcl[d].ad->sqfrom,
		      th->hit[h]->dcl[d].ad->sqto,
		      (th->hit[h]->dcl[d].ad->sqfrom == 1) ? '[' : '.',
		      (th->hit[h]->dcl[d].ad->sqto   == th->hit[h]->dcl[d].ad->L) ? ']' : '.',
		      th->hit[h]->dcl[d].ienv,
		      th->hit[h]->dcl[d].jenv,
		      (th->hit[h]->dcl[d].ienv == 1) ? '[' : '.',
		      (th->hit[h]->dcl[d].jenv == th->hit[h]->dcl[d].ad->L) ? ']' : '.',
		      (th->hit[h]->dcl[d].oasc / (1.0 + fabs((float) (th->hit[h]->dcl[d].jenv - th->hit[h]->dcl[d].ienv)))));
	    }

	
	fprintf(ofp, "\n  Alignments for each domain:\n");
	nd = 0;
	for (d = 0; d < th->hit[h]->ndom; d++)
	  if (th->hit[h]->dcl[d].is_reported) 
	    {
	      nd++;
	      fprintf(ofp, "  == domain %d    score: %.1f bits;  conditional E-value: %.2g\n",
		      nd, 
		      th->hit[h]->dcl[d].bitscore,
		      th->hit[h]->dcl[d].pvalue * pli->domZ);
	      p7_alidisplay_Print(ofp, th->hit[h]->dcl[d].ad, 40, textw);
	      fprintf(ofp, "\n");
	    }
      }
  if (th->nreported == 0) { fprintf(ofp, "\n   [No hits detected that satisfy reporting thresholds]\n"); return eslOK; }
  return eslOK;
}


/* Function:  p7_tophits_Alignment()
 * Synopsis:  Create a multiple alignment of all the included domains.
 * Incept:    SRE, Wed Dec 10 11:04:40 2008 [Janelia]
 *
 * Purpose:   Create a multiple alignment of all domains marked
 *            "includable" in the top hits list <th>, and return it in
 *            <*ret_msa>.
 *            
 *            Use of <optflags> is identical to <optflags> in <p7_MultipleAlignment()>.
 *            Possible flags include <p7_DIGITIZE>, <p7_ALL_CONSENSUS_COLS>,
 *            and <p7_TRIM>; they may be OR'ed together. Otherwise, pass
 *            <p7_DEFAULT> to set no flags.
 *
 *            Caller may optionally provide <inc_sqarr>, <inc_trarr>, and
 *            <inc_n> to include additional sequences in the alignment
 *            (the jackhmmer query, for example). Otherwise, pass <NULL, NULL, 0>.
 *
 * Returns:   <eslOK> on success, and <*ret_msa> points to a new MSA that
 *            the caller is responsible for freeing.
 *
 *            Returns <eslFAIL> if there are no reported domains that
 *            satisfy reporting thresholds, in which case <*ret_msa>
 *            is <NULL>.
 *
 * Throws:    <eslEMEM> on allocation failure; <eslECORRUPT> on 
 *            unexpected internal data corruption.
 *
 * Xref:      J4/29: incept.
 *            J4/76: added inc_sqarr, inc_trarr, inc_n, optflags 
 */
int
p7_tophits_Alignment(const P7_TOPHITS *th, const ESL_ALPHABET *abc, 
		     ESL_SQ **inc_sqarr, P7_TRACE **inc_trarr, int inc_n,
		     int optflags, ESL_MSA **ret_msa)
{
  ESL_SQ   **sqarr = NULL;
  P7_TRACE **trarr = NULL;
  ESL_MSA   *msa   = NULL;
  int        ndom  = 0;
  int        h, d, y;
  int        M;
  int        status;

  /* How many domains will be included in the new alignment? 
   * We also set model size M here; every alignment has a copy.
   */
  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_INCLUDED)
      {
	for (d = 0; d < th->hit[h]->ndom; d++)
	  if (th->hit[h]->dcl[d].is_included) 
	    ndom++;
      }
  if (inc_n+ndom == 0) { status = eslFAIL; goto ERROR; }

  if (inc_n)     M = inc_trarr[0]->M;          
  else           M = th->hit[0]->dcl[0].ad->M;
  
  /* Allocation */
  ESL_ALLOC(sqarr, sizeof(ESL_SQ *)   * (ndom + inc_n));
  ESL_ALLOC(trarr, sizeof(P7_TRACE *) * (ndom + inc_n));
  /* Inclusion of preexisting seqs, traces: make copy of pointers */
  for (y = 0; y < inc_n;        y++) { sqarr[y] = inc_sqarr[y];  trarr[y] = inc_trarr[y]; }
  for (;      y < (ndom+inc_n); y++) { sqarr[y] = NULL;          trarr[y] = NULL; }

  /* Make faux sequences, traces from hit list */
  y = inc_n;
  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_INCLUDED)
      {
	for (d = 0; d < th->hit[h]->ndom; d++)
	  if (th->hit[h]->dcl[d].is_included) 
	    {
	      if ((status = p7_alidisplay_Backconvert(th->hit[h]->dcl[d].ad, abc, &(sqarr[y]), &(trarr[y]))) != eslOK) goto ERROR;
	      y++;
	    }
      }
  
  /* Make the multiple alignment */
  if ((status = p7_tracealign_Seqs(sqarr, trarr, inc_n+ndom, M, optflags, &msa)) != eslOK) goto ERROR;

  /* Clean up */
  for (y = inc_n; y < ndom+inc_n; y++) esl_sq_Destroy(sqarr[y]);
  for (y = inc_n; y < ndom+inc_n; y++) p7_trace_Destroy(trarr[y]);
  free(sqarr);
  free(trarr);
  *ret_msa = msa;
  return eslOK;
  
 ERROR:
  if (sqarr != NULL) { for (y = inc_n; y < ndom+inc_n; y++) if (sqarr[y] != NULL) esl_sq_Destroy(sqarr[y]);   free(sqarr); }
  if (trarr != NULL) { for (y = inc_n; y < ndom+inc_n; y++) if (trarr[y] != NULL) p7_trace_Destroy(trarr[y]); free(trarr); }
  if (msa   != NULL) esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;
}
/*---------------- end, standard output format ------------------*/





/*****************************************************************
 * 3. Tabular (parsable) output of pipeline results.
 *****************************************************************/

/* Function:  p7_tophits_TabularTargets()
 * Synopsis:  Output parsable table of per-sequence hits.
 * Incept:    SRE, Wed Mar 18 15:26:17 2009 [Janelia]
 *
 * Purpose:   Output a parseable table of reportable per-sequence hits
 *            in sorted tophits list <th> in an easily parsed ASCII
 *            tabular form to stream <ofp>, using final pipeline
 *            accounting stored in <pli>.
 *            
 *            Designed to be concatenated for multiple queries and
 *            multiple top hits list.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_TabularTargets(FILE *ofp, char *queryname, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header)
{
  int qnamew = ESL_MAX(20, strlen(queryname));
  int tnamew = ESL_MAX(20, p7_tophits_GetMaxNameLength(th));
  int h,d;
  
  if (show_header)
    {  
      fprintf(ofp, "#%*s %22s %22s %33s\n", tnamew+qnamew, "", "--- full sequence ----", "--- best 1 domain ----", "--- domain number estimation ----");
      fprintf(ofp, "#%-*s %-*s %9s %6s %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n", 
	      tnamew-1, " target",             qnamew, "query",                "  E-value", " score", " bias", "  E-value", " score", " bias", "exp", "reg", "clu", " ov", "env", "dom", "rep", "inc", "description of target");
      fprintf(ofp, "#%*s %*s %9s %6s %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n", 
	      tnamew-1, "-------------------", qnamew, "--------------------", "---------", "------", "-----", "---------", "------", "-----", "---", "---", "---", "---", "---", "---", "---", "---", "---------------------");
     }

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)
      {
	d    = th->hit[h]->best_domain;

	fprintf(ofp, "%-*s %-*s %9.2g %6.1f %5.1f %9.2g %6.1f %5.1f %5.1f %3d %3d %3d %3d %3d %3d %3d %s\n", 
		tnamew, th->hit[h]->name,
		qnamew, queryname, 
		th->hit[h]->pvalue * pli->Z,
		th->hit[h]->score,
		th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
		th->hit[h]->dcl[d].pvalue * pli->Z,
		th->hit[h]->dcl[d].bitscore,
		th->hit[h]->dcl[d].dombias,
		th->hit[h]->nexpected,
		th->hit[h]->nregions,
		th->hit[h]->nclustered,
		th->hit[h]->noverlaps,
		th->hit[h]->nenvelopes,
		th->hit[h]->ndom,
		th->hit[h]->nreported,
		th->hit[h]->nincluded,
		(th->hit[h]->desc == NULL ? "-" : th->hit[h]->desc));
      }
  return eslOK;
}


/* Function:  p7_tophits_TabularDomains()
 * Synopsis:  Output parseable table of per-domain hits
 * Incept:    SRE, Wed Mar 18 16:57:58 2009 [Janelia]
 *
 * Purpose:   Output a parseable table of reportable per-domain hits
 *            in sorted tophits list <th> in an easily parsed ASCII
 *            tabular form to stream <ofp>, using final pipeline
 *            accounting stored in <pli>.
 *            
 *            Designed to be concatenated for multiple queries and
 *            multiple top hits list.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_TabularDomains(FILE *ofp, char *queryname, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header)
{
  int qnamew = ESL_MAX(20, strlen(queryname));
  int tnamew = ESL_MAX(20, p7_tophits_GetMaxNameLength(th));
  int h,d,nd;

  if (show_header)
    {
      fprintf(ofp, "#%*s %22s %40s %11s %11s %11s\n", tnamew+qnamew-1+13, "",                    "--- full sequence ---",        "-------------- this domain -------------",                "hmm coord",      "ali coord",     "env coord");
      fprintf(ofp, "#%-*s %5s %-*s %5s %9s %6s %5s %3s %3s %9s %9s %6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n",
	      tnamew-1, " target",            "tlen",   qnamew, "query",                "qlen",  "E-value",   "score",  "bias",  "#",   "of",  "c-Evalue",  "i-Evalue",  "score",  "bias",  "from",  "to",    "from",  "to",   "from",   "to",    "acc",  "description of target");
      fprintf(ofp, "#%*s %5s %*s %5s %9s %6s %5s %3s %3s %9s %9s %6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n", 
	      tnamew-1, "-------------------", "-----", qnamew, "--------------------", "-----", "---------", "------", "-----", "---", "---", "---------", "---------", "------", "-----", "-----", "-----", "-----", "-----", "-----", "-----", "----", "---------------------");
    }

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)
      {
	nd = 0;
	for (d = 0; d < th->hit[h]->ndom; d++)
	  if (th->hit[h]->dcl[d].is_reported) 
	    {
	      nd++;
	      fprintf(ofp, "%-*s %5ld %-*s %5d %9.2g %6.1f %5.1f %3d %3d %9.2g %9.2g %6.1f %5.1f %5d %5d %5ld %5ld %5d %5d %4.2f %s\n", 
		      tnamew, th->hit[h]->name, th->hit[h]->dcl[d].ad->L,
		      qnamew, queryname,        th->hit[h]->dcl[d].ad->M,
		      th->hit[h]->pvalue * pli->Z,
		      th->hit[h]->score,
		      th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
		      nd, 
		      th->hit[h]->nreported,
		      th->hit[h]->dcl[d].pvalue * pli->domZ,
		      th->hit[h]->dcl[d].pvalue * pli->Z,
		      th->hit[h]->dcl[d].bitscore,
		      th->hit[h]->dcl[d].dombias,
		      th->hit[h]->dcl[d].ad->hmmfrom,
		      th->hit[h]->dcl[d].ad->hmmto,
		      th->hit[h]->dcl[d].ad->sqfrom,
		      th->hit[h]->dcl[d].ad->sqto,
		      th->hit[h]->dcl[d].ienv,
		      th->hit[h]->dcl[d].jenv,
		      (th->hit[h]->dcl[d].oasc / (1.0 + fabs((float) (th->hit[h]->dcl[d].jenv - th->hit[h]->dcl[d].ienv)))),
		      (th->hit[h]->desc == NULL ? "-" : th->hit[h]->desc));		      
	    }
      }
  return eslOK;
}
/*------------------- end, tabular output -----------------------*/




/*****************************************************************
 * 4. Benchmark driver
 *****************************************************************/
#ifdef p7TOPHITS_BENCHMARK
/* 
  gcc -o benchmark-tophits -std=gnu99 -g -O2 -I. -L. -I../easel -L../easel -Dp7TOPHITS_BENCHMARK p7_tophits.c -lhmmer -leasel -lm 
  ./benchmark-tophits

  As of 28 Dec 07, shows 0.20u for 10 lists of 10,000 hits each (at least ~100x normal expectation),
  so we expect top hits list time to be negligible for typical hmmsearch/hmmscan runs.
  
  If needed, we do have opportunity for optimization, however - especially in memory handling.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_stopwatch.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-M",        eslARG_INT,     "10", NULL, NULL,  NULL,  NULL, NULL, "number of top hits lists to simulate and merge",   0 },
  { "-N",        eslARG_INT,  "10000", NULL, NULL,  NULL,  NULL, NULL, "number of top hits to simulate",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "benchmark driver for P7_TOPHITS";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_STOPWATCH  *w        = esl_stopwatch_Create();
  ESL_RANDOMNESS *r        = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             N        = esl_opt_GetInteger(go, "-N");
  int             M        = esl_opt_GetInteger(go, "-M");
  P7_TOPHITS    **h        = NULL;
  double         *sortkeys = NULL;
  char            name[]   = "not_unique_name";
  char            acc[]    = "not_unique_acc";
  char            desc[]   = "Test description for the purposes of making the benchmark allocate space";
  int             i,j;
  int             status;

  /* prep work: generate our sort keys before starting to time anything    */
  ESL_ALLOC(h,        sizeof(P7_TOPHITS *) * M); /* allocate pointers for M lists */
  ESL_ALLOC(sortkeys, sizeof(double) * N * M);   
  for (i = 0; i < N*M; i++) sortkeys[i] = esl_random(r);

  esl_stopwatch_Start(w);

  /* generate M "random" lists and sort them */
  for (j = 0; j < M; j++)
    {
      h[j] = p7_tophits_Create();
      for (i = 0; i < N; i++)
	p7_tophits_Add(h[j], name, acc, desc, sortkeys[j*N + i], 
		       (float) sortkeys[j*N+i], sortkeys[j*N+i],
		       (float) sortkeys[j*N+i], sortkeys[j*N+i],
		       i, i, N,
		       i, i, N,
		       i, N, NULL);
      p7_tophits_Sort(h[j]);
    }
  /* then merge them into one big list in h[0] */
  for (j = 1; j < M; j++)
    {
      p7_tophits_Merge(h[0], h[j]);
      p7_tophits_Destroy(h[j]);
    }      

  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  p7_tophits_Destroy(h[0]);
  status = eslOK;
 ERROR:
  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  if (sortkeys != NULL) free(sortkeys);
  if (h != NULL) free(h);
  return status;
}
#endif /*p7TOPHITS_BENCHMARK*/


/*****************************************************************
 * 5. Test driver
 *****************************************************************/

#ifdef p7TOPHITS_TESTDRIVE
/*
  gcc -o tophits_utest -std=gnu99 -g -O2 -I. -L. -I../easel -L../easel -Dp7TOPHITS_TESTDRIVE p7_tophits.c -lhmmer -leasel -lm 
  ./tophits_test
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_stopwatch.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of top hits to simulate",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options]";
static char banner[] = "test driver for P7_TOPHITS";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r        = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             N        = esl_opt_GetInteger(go, "-N");
  P7_TOPHITS     *h1       = NULL;
  P7_TOPHITS     *h2       = NULL;
  P7_TOPHITS     *h3       = NULL;
  char            name[]   = "not_unique_name";
  char            acc[]    = "not_unique_acc";
  char            desc[]   = "Test description for the purposes of making the test driver allocate space";
  double          key;
  int             i;

  h1 = p7_tophits_Create();
  h2 = p7_tophits_Create();
  h3 = p7_tophits_Create();
  
  for (i = 0; i < N; i++) 
    {
      key = esl_random(r);
      p7_tophits_Add(h1, name, acc, desc, key, (float) key, key, (float) key, key, i, i, N, i, i, N, 1, 1, NULL);
      key = 10.0 * esl_random(r);
      p7_tophits_Add(h2, name, acc, desc, key, (float) key, key, (float) key, key, i, i, N, i, i, N, 2, 2, NULL);
      key = 0.1 * esl_random(r);
      p7_tophits_Add(h3, name, acc, desc, key, (float) key, key, (float) key, key, i, i, N, i, i, N, 3, 3, NULL);
    }
  p7_tophits_Add(h1, "last",  NULL, NULL, -1.0, (float) key, key, (float) key, key, i, i, N, i, i, N, 1, 1, NULL);
  p7_tophits_Add(h1, "first", NULL, NULL, 20.0, (float) key, key, (float) key, key, i, i, N, i, i, N, 1, 1, NULL);

  p7_tophits_Sort(h1);
  if (strcmp(h1->hit[0]->name,   "first") != 0) esl_fatal("sort failed (top is %s = %f)", h1->hit[0]->name,   h1->hit[0]->sortkey);
  if (strcmp(h1->hit[N+1]->name, "last")  != 0) esl_fatal("sort failed (last is %s = %f)", h1->hit[N+1]->name, h1->hit[N+1]->sortkey);

  p7_tophits_Merge(h1, h2);
  if (strcmp(h1->hit[0]->name,     "first") != 0) esl_fatal("after merge 1, sort failed (top is %s = %f)", h1->hit[0]->name,     h1->hit[0]->sortkey);
  if (strcmp(h1->hit[2*N+1]->name, "last")  != 0) esl_fatal("after merge 1, sort failed (last is %s = %f)", h1->hit[2*N+1]->name, h1->hit[2*N+1]->sortkey);

  p7_tophits_Merge(h3, h1);
  if (strcmp(h3->hit[0]->name,     "first") != 0) esl_fatal("after merge 2, sort failed (top is %s = %f)", h3->hit[0]->name,     h3->hit[0]->sortkey);
  if (strcmp(h3->hit[3*N+1]->name, "last")  != 0) esl_fatal("after merge 2, sort failed (last is %s = %f)", h3->hit[3*N+1]->name,     h3->hit[3*N+1]->sortkey);
  
  if (p7_tophits_GetMaxNameLength(h3) != strlen(name)) esl_fatal("GetMaxNameLength() failed");

  p7_tophits_Destroy(h1);
  p7_tophits_Destroy(h2);
  p7_tophits_Destroy(h3);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7TOPHITS_TESTDRIVE*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/





