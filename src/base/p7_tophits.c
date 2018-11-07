/* P7_TOPHITS: implementation of ranked list of top-scoring hits
 * 
 * Contents:
 *    1. The P7_TOPHITS object.
 *    2. The P7_HIT object array in P7_TOPHITS.
 *    3. Debugging and development tools.
 *    4. Benchmark driver.
 *    5. Test driver.
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "base/p7_tophits.h"

/*****************************************************************
 *= 1. The P7_TOPHITS object
 *****************************************************************/

/* Function:  p7_tophits_Create()
 * Synopsis:  Allocate a hit list.
 *
 * Purpose:   Allocates a new <P7_TOPHITS> hit list, for an initial
 *            allocation of <int_hit_alloc> hits (this will be grown
 *            later as needed). Return a pointer to it.
 *            
 * Args:      init_hit_alloc  - initial allocation size, # of hits.
 *                              Often p7_TOPHITS_DEFAULT_INIT_ALLOC.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_TOPHITS *
p7_tophits_Create(int init_hit_alloc)
{
  P7_TOPHITS *h = NULL;
  int         status;

  ESL_ALLOC(h, sizeof(P7_TOPHITS));
  h->hit    = NULL;
  h->unsrt  = NULL;

  if (( h->unsrt = p7_hit_Create(init_hit_alloc) ) == NULL) goto ERROR;
  ESL_ALLOC(h->hit, sizeof(P7_HIT *) * init_hit_alloc);

  h->Nalloc               = init_hit_alloc;
  h->N                    = 0;
  h->nreported            = 0;
  h->nincluded            = 0;
  h->is_sorted_by_sortkey = TRUE;     /* but only because there's 0 hits */
  h->is_sorted_by_seqidx  = FALSE;
  h->hit[0]               = h->unsrt; /* if you're going to call it "sorted" when it contains just one hit, you need this */
  return h;

 ERROR:
  p7_tophits_Destroy(h);
  return NULL;
}


/* Function:  p7_tophits_Grow()
 * Synopsis:  Reallocates a larger hit list, if needed.
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
  void     *p;
  P7_HIT   *ori    = h->unsrt;
  int       Nalloc = h->Nalloc * 2;    /* grow by doubling */
  int       i;
  int       status;

  if (h->N < h->Nalloc) return eslOK; /* we have enough room for another hit */

  if (( status = p7_hit_Grow( &(h->unsrt), h->Nalloc, Nalloc)) != eslOK) goto ERROR;
  ESL_RALLOC(h->hit,   p, sizeof(P7_HIT *) * Nalloc);

  /* If we grow a sorted list, we have to translate the pointers
   * in h->hit, because h->unsrt might have just moved in memory. 
   */
  if (h->is_sorted_by_seqidx || h->is_sorted_by_sortkey)
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
  if (h->N >= 2) 
  {
      h->is_sorted_by_seqidx = FALSE;
      h->is_sorted_by_sortkey = FALSE;
  }

  hit->name         = NULL;
  hit->acc          = NULL;
  hit->desc         = NULL;
  hit->sortkey      = 0.0;

  hit->score        = 0.0;
  hit->pre_score    = 0.0;
  hit->sum_score    = 0.0;

  hit->lnP          = 0.0;
  hit->pre_lnP      = 0.0;
  hit->sum_lnP      = 0.0;

  hit->ndom         = 0;
  hit->noverlaps    = 0;
  hit->nexpected    = 0.0;

  hit->flags        = p7_HITFLAGS_DEFAULT;
  hit->nreported    = 0;
  hit->nincluded    = 0;
  hit->best_domain  = 0;
  hit->dcl          = NULL;
  hit->offset       = 0;

  *ret_hit = hit;
  return eslOK;

 ERROR:
  *ret_hit = NULL;
  return status;
}



/* hit_sorter(): qsort's pawn, below */
static int
hit_sorter_by_sortkey(const void *vh1, const void *vh2)
{
  P7_HIT *h1 = *((P7_HIT **) vh1);  
  P7_HIT *h2 = *((P7_HIT **) vh2);
  int     c;

  if      (h1->sortkey < h2->sortkey) return  1;
  else if (h1->sortkey > h2->sortkey) return -1;
  else {
    if ( (c = strcmp(h1->name, h2->name)) != 0) return c;

    /* if on different strand, the positive strand goes first, else use position */
    int dir1 = (h1->dcl[0].ia < h1->dcl[0].ib ? 1 : -1);
    int dir2 = (h2->dcl[0].ia < h2->dcl[0].ib ? 1 : -1);
    if (dir1 != dir2) return dir2; // so if dir1 is pos (1), and dir2 is neg (-1), this will return -1, placing h1 before h2;  otherwise, vice versa
    else {
      if     (h1->dcl[0].ia > h2->dcl[0].ia) return  1; 
      else if(h1->dcl[0].ia < h2->dcl[0].ia) return -1; 
      else                                   return  0;
    }
  }
}

static int
hit_sorter_by_seqidx_aliposition(const void *vh1, const void *vh2)
{
  P7_HIT *h1 = *((P7_HIT **) vh1);  
  P7_HIT *h2 = *((P7_HIT **) vh2);

  if      (h1->seqidx > h2->seqidx) return  1; /* first key, seq_idx (unique id for sequences), low to high */
  else if (h1->seqidx < h2->seqidx) return -1;
  // if on different strand, the positive strand goes first, else use position
  int dir1 = (h1->dcl[0].ia < h1->dcl[0].ib ? 1 : -1);
  int dir2 = (h2->dcl[0].ia < h2->dcl[0].ib ? 1 : -1);

  if (dir1 != dir2) return dir2; // so if dir1 is pos (1), and dir2 is neg (-1), this will return -1, placing h1 before h2;  otherwise, vice versa

  if      (h1->dcl[0].ia > h2->dcl[0].ia) return  1; 
  else if (h1->dcl[0].ia < h2->dcl[0].ia) return -1;
  else if (h1->dcl[0].ib < h2->dcl[0].ib) return  1; 
  else if (h1->dcl[0].ib > h2->dcl[0].ib) return -1; 
  return  0;
}

static int
hit_sorter_by_modelname_aliposition(const void *vh1, const void *vh2)
{
  P7_HIT *h1 = *((P7_HIT **) vh1);  
  P7_HIT *h2 = *((P7_HIT **) vh2);

  int res = esl_strcmp( h1->name, h2->name);

  if  ( res != 0 ) return  res; /* first key, seq_idx (unique id for sequences), low to high */

  // if on different strand, the positive strand goes first, else use position
  int dir1 = (h1->dcl[0].ia < h1->dcl[0].ib ? 1 : -1);
  int dir2 = (h2->dcl[0].ia < h2->dcl[0].ib ? 1 : -1);

  if (dir1 != dir2) return dir2; // so if dir1 is pos (1), and dir2 is neg (-1), this will return -1, placing h1 before h2;  otherwise, vice versa

  if      (h1->dcl[0].ia > h2->dcl[0].ia) return  1; 
  else if (h1->dcl[0].ia < h2->dcl[0].ia) return -1;
  else if (h1->dcl[0].ib < h2->dcl[0].ib) return  1; 
  else if (h1->dcl[0].ib > h2->dcl[0].ib) return -1; 
  return  0;
}


/* Function:  p7_tophits_SortBySortkey()
 * Synopsis:  Sorts a hit list.
 *
 * Purpose:   Sorts a top hit list. After this call,
 *            <h->hit[i]> points to the i'th ranked 
 *            <P7_HIT> for all <h->N> hits.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_SortBySortkey(P7_TOPHITS *th)
{
  int h;

  if (th->is_sorted_by_sortkey)  return eslOK;
  for (h = 0; h < th->N; h++) th->hit[h] = th->unsrt + h;
  if (th->N > 1)  qsort(th->hit, th->N, sizeof(P7_HIT *), hit_sorter_by_sortkey);
  th->is_sorted_by_seqidx  = FALSE;
  th->is_sorted_by_sortkey = TRUE;
  return eslOK;
}


/* Function:  p7_tophits_SortBySeqidxAndAlipos()
 * Synopsis:  Sorts a hit list by sequence index and position in that
 *            sequence at which the hit's first domain begins (used in nhmmer)
 *
 * Purpose:   Sorts a top hit list. After this call,
 *            <h->hit[i]> points to the i'th ranked
 *            <P7_HIT> for all <h->N> hits.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_SortBySeqidxAndAlipos(P7_TOPHITS *th)
{
  int h;

  if (th->is_sorted_by_seqidx)  return eslOK;
  for (h = 0; h < th->N; h++) th->hit[h] = th->unsrt + h;
  if (th->N > 1)  qsort(th->hit, th->N, sizeof(P7_HIT *), hit_sorter_by_seqidx_aliposition);
  th->is_sorted_by_sortkey = FALSE;
  th->is_sorted_by_seqidx  = TRUE;
  return eslOK;
}

/* Function:  p7_tophits_SortByModelnameAndAlipos()
 * Synopsis:  Sorts a hit list by model name and position in the query sequence
 *            sequence at which the hit's first domain begins (used in nhmmscan)
 *
 * Purpose:   Sorts a top hit list. After this call,
 *            <h->hit[i]> points to the i'th ranked
 *            <P7_HIT> for all <h->N> hits.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_SortByModelnameAndAlipos(P7_TOPHITS *th)
{
  int h;

  if (th->is_sorted_by_seqidx)  return eslOK;
  for (h = 0; h < th->N; h++) th->hit[h] = th->unsrt + h;
  if (th->N > 1)  qsort(th->hit, th->N, sizeof(P7_HIT *), hit_sorter_by_modelname_aliposition);
  th->is_sorted_by_sortkey = FALSE;
  th->is_sorted_by_seqidx  = TRUE;
  return eslOK;
}


/* Function:  p7_tophits_Merge()
 * Synopsis:  Merge two top hits lists.
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
  P7_HIT **new_hit = NULL;
  P7_HIT  *ori1    = h1->unsrt;    /* original base of h1's data */
  P7_HIT  *new2;
  int      i,j,k;
  int      Nalloc = h1->Nalloc + h2->Nalloc;
  int      status;

  /* Make sure the two lists are sorted */
  if ((status = p7_tophits_SortBySortkey(h1)) != eslOK) goto ERROR;
  if ((status = p7_tophits_SortBySortkey(h2)) != eslOK) goto ERROR;

  /* Attempt our allocations, so we fail early if we fail. 
   * Reallocating h1->unsrt screws up h1->hit, so fix it.
   */
  if (( status = p7_hit_Grow( &(h1->unsrt), h1->Nalloc, Nalloc)) != eslOK) goto ERROR;
  ESL_ALLOC (new_hit, sizeof(P7_HIT *)    * Nalloc);
  for (i = 0; i < h1->N; i++)
    h1->hit[i] = h1->unsrt + (h1->hit[i] - ori1);

  /* Append h2's unsorted data array to h1. h2's data begin at <new2> */
  new2 = h1->unsrt + h1->N;
  memcpy(new2, h2->unsrt, sizeof(P7_HIT) * h2->N);

  /* Merge the sorted hit lists */
  for (i=0,j=0,k=0; i < h1->N && j < h2->N ; k++)
    new_hit[k] = (hit_sorter_by_sortkey(&h1->hit[i], &h2->hit[j]) > 0) ? new2 + (h2->hit[j++] - h2->unsrt) : h1->hit[i++];
  while (i < h1->N) new_hit[k++] = h1->hit[i++];
  while (j < h2->N) new_hit[k++] = new2 + (h2->hit[j++] - h2->unsrt);

  /* h2 now turns over management of name, acc, desc memory to h1;
   * nullify its pointers, to prevent double free.  */
  for (i = 0; i < h2->N; i++)
  {
      h2->unsrt[i].name = NULL;
      h2->unsrt[i].acc  = NULL;
      h2->unsrt[i].desc = NULL;
      h2->unsrt[i].dcl  = NULL;
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

/* Function:  p7_tophits_GetMaxPositionLength()
 * Synopsis:  Returns maximum position length in hit list (targets).
 *
 * Purpose:   Returns the length of the longest hit location (start/end)
 *               of all the registered hits, in chars. This is useful when
 *               deciding how to format output.
 *
 *            The maximum is taken over all registered hits. This
 *            opens a possible side effect: caller might print only
 *            the top hits, and the max name length in these top hits
 *            may be different than the max length over all the hits.
 *
 *            Used specifically for nhmmer output, so expects only one
 *            domain per hit
 *
 *            If there are no hits in <h>, or none of the
 *            hits have names, returns 0.
 */
int
p7_tophits_GetMaxPositionLength(P7_TOPHITS *h)
{
  int i, max, n;
  char buffer [11];

  for (max = 0, i = 0; i < h->N; i++) {
    if (h->unsrt[i].dcl[0].ia > 0) {
      n = sprintf (buffer, "%d", h->unsrt[i].dcl[0].ia);
      max = ESL_MAX(n, max);
      n = sprintf (buffer, "%d", h->unsrt[i].dcl[0].ib);
      max = ESL_MAX(n, max);
    }
  }
  return max;
}

/* Function:  p7_tophits_GetMaxNameLength()
 * Synopsis:  Returns maximum name length in hit list (targets).
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

/* Function:  p7_tophits_GetMaxAccessionLength()
 * Synopsis:  Returns maximum accession length in hit list (targets).
 *
 * Purpose:   Same as <p7_tophits_GetMaxNameLength()>, but for
 *            accessions. If there are no hits in <h>, or none
 *            of the hits have accessions, returns 0.
 */
int
p7_tophits_GetMaxAccessionLength(P7_TOPHITS *h)
{
  int i, max, n;
  for (max = 0, i = 0; i < h->N; i++)
    if (h->unsrt[i].acc != NULL) {
      n   = strlen(h->unsrt[i].acc);
      max = ESL_MAX(n, max);
    }
  return max;
}

/* Function:  p7_tophits_GetMaxShownLength()
 * Synopsis:  Returns max shown name/accession length in hit list.
 *
 * Purpose:   Same as <p7_tophits_GetMaxNameLength()>, but 
 *            for the case when --acc is on, where
 *            we show accession if one is available, and 
 *            fall back to showing the name if it is not.
 *            Returns the max length of whatever is being
 *            shown as the reported "name".
 */
int
p7_tophits_GetMaxShownLength(P7_TOPHITS *h)
{
  int i, max, n;
  for (max = 0, i = 0; i < h->N; i++)
  {
    if (h->unsrt[i].acc != NULL && h->unsrt[i].acc[0] != '\0')
    {
      n   = strlen(h->unsrt[i].acc);
      max = ESL_MAX(n, max);
    }
    else if (h->unsrt[i].name != NULL)
    {
      n   = strlen(h->unsrt[i].name);
      max = ESL_MAX(n, max);
    }
  }
  return max;
}


/* Function:  p7_tophits_Reuse()
 * Synopsis:  Reuse a hit list, freeing internals.
 *
 * Purpose:   Reuse the tophits list <h>; save as 
 *            many malloc/free cycles as possible,
 *            as opposed to <Destroy()>'ing it and
 *            <Create>'ing a new one.
 */
int
p7_tophits_Reuse(P7_TOPHITS *th)
{
  int i;

  if (th == NULL) return eslOK;
  if (th->unsrt)
    {
      for (i = 0; i < th->N; i++)
	{
	  if (th->unsrt[i].name) { free(th->unsrt[i].name); th->unsrt[i].name = NULL; }
	  if (th->unsrt[i].acc)  { free(th->unsrt[i].acc);  th->unsrt[i].acc  = NULL; }
	  if (th->unsrt[i].desc) { free(th->unsrt[i].desc); th->unsrt[i].desc = NULL; }
	  if (th->unsrt[i].dcl)  { p7_domain_Destroy(th->unsrt[i].dcl, th->unsrt[i].ndom); th->unsrt[i].dcl = NULL; th->unsrt[i].ndom = 0; }
	}
    }
  th->N         = 0;
  th->is_sorted_by_seqidx = FALSE;
  th->is_sorted_by_sortkey = TRUE;  /* because there are 0 hits */
  th->hit[0]    = th->unsrt;
  return eslOK;
}

/* Function:  p7_tophits_Destroy()
 * Synopsis:  Frees a hit list.
 */
void
p7_tophits_Destroy(P7_TOPHITS *th)
{
  if (th) {
    if (th->hit)   free(th->hit);
    if (th->unsrt) p7_hit_Destroy(th->unsrt, th->N); /* destroy any allocated contents of the array of h->N hit objects we initialized/used */
    free(th);
  }
}
/*---------------- end, P7_TOPHITS object -----------------------*/


/*****************************************************************
 * 2. The P7_HIT object array in P7_TOPHITS
 *****************************************************************/

/* Function:  p7_hit_Create()
 * Synopsis:  Create a new array of P7_HITs
 *
 * Purpose:   Create a new array of <nhit_alloc> <P7_HIT> objects.
 *            Return a ptr to this new array.
 *
 *            Initialize all the memory ptrs in these structures, so
 *            it is safe to call <p7_hit_Destroy(ptr, nhit_alloc)>
 *            on the returned ptr. 
 */
P7_HIT *
p7_hit_Create(int nhit_alloc)
{
  P7_HIT *hit = NULL;
  int     h;
  int     status;

  ESL_ALLOC(hit, sizeof(P7_HIT) * nhit_alloc);
  for (h = 0; h < nhit_alloc; h++)
    {
      hit[h].name = NULL;
      hit[h].acc  = NULL;
      hit[h].desc = NULL;
      hit[h].dcl  = NULL;
      hit[h].ndom = 0;
    }
  return hit;

 ERROR:
  return NULL;
}
  

/* Function:  p7_hit_Grow()
 * Synopsis:  Change the allocation of a P7_HIT array.
 *
 * Purpose:   Given a ptr <*hitp> to a <P7_HIT> array,
 *            the old allocation size <oldalloc>, and
 *            a new allocation size <newalloc>;
 *            reallocate the array <*hitp>.
 *
 * Returns:   <eslOK> on success; <*hitp> may have moved.
 *
 * Throws:    <eslEMEM> on reallocation failure. <*hitp> is
 *            unchanged, as is the array's contents.
 */
int
p7_hit_Grow(P7_HIT **hitp, int oldalloc, int newalloc)
{
  int h;
  int status;

  ESL_REALLOC( (*hitp), sizeof(P7_HIT) * newalloc);
  for (h = oldalloc; h < newalloc; h++)
    {
      (*hitp)[h].name = NULL;
      (*hitp)[h].acc  = NULL;
      (*hitp)[h].desc = NULL;
      (*hitp)[h].dcl  = NULL;
      (*hitp)[h].ndom = 0;
    }
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_hit_Destroy()
 * Synopsis:  Free a P7_HIT array.
 *
 * Purpose:   Free array <hits> of <nhits> <P7_HIT> strucutures.
 */
void
p7_hit_Destroy(P7_HIT *hits, int nhits)
{
  int i;
  if (hits) {
    for (i = 0; i < nhits; i++)
      {
	if (hits[i].name) free(hits[i].name);
	if (hits[i].acc)  free(hits[i].acc);
	if (hits[i].desc) free(hits[i].desc);      
	if (hits[i].dcl)  p7_domain_Destroy(hits[i].dcl, hits[i].ndom);
      }
    free(hits);
  }
}
/*------------------ end, P7_HIT --------------------------------*/



/*****************************************************************
 * 3. Debugging and development tools.
 *****************************************************************/

/* Function:  p7_tophits_TestSample()
 * Synopsis:  Sample a random, bogus, mostly-syntactically-valid P7_TOPHITS
 *
 * Purpose:   Sample a random but syntactically valid <P7_TOPHITS>
 *            structure, using random number generator <rng>.  Return
 *            it thru <*ret_th>. It was allocated here; caller becomes
 *            responsible for freeing it with <p7_tophits_Destroy()>.
 *            
 *            <th->hit[]> 'sorted' array of ptrs is put in a
 *            randomized order. The <th->sortkey> value and
 *            <th->is_sorted*> flags are set randomly, and have
 *            nothing to do with the order of <th->hit[]>.  (Main use
 *            here is for testing faithful communication of the
 *            object, including its sorted ptrs.) 
 *
 * Returns:   <eslOK> on success, and <*ret_th> points to the sampled
 *            <P7_TOPHITS> object.
 *
 * Throws:    (no abnormal error conditions)
 * 
 * Notes:     Easel code spec requires that TestSample() generates an
 *            object that passes whatever Validate() looks for.
 */
int
p7_tophits_TestSample(ESL_RANDOMNESS *rng, P7_TOPHITS **ret_th)
{
  P7_TOPHITS *th    = NULL;
  int         nhits = 1049;	/* prime. don't make it divisible by any chunk size. */
  int         h,n;
  int         status;

  if (( th = p7_tophits_Create(nhits)) == NULL) { status = eslEMEM; goto ERROR; }
  th->nreported = 1+esl_rnd_Roll(rng, nhits);
  th->nincluded = 1+esl_rnd_Roll(rng, nhits);
  th->is_sorted_by_sortkey = esl_rnd_Roll(rng, 2);
  if (! th->is_sorted_by_sortkey) th->is_sorted_by_seqidx = esl_rnd_Roll(rng, 2);
  
  for (h = 0; h < nhits; h++)
    {
      if ( (status = p7_hit_TestSample(rng, &(th->unsrt[h]))) != eslOK) goto ERROR;
      th->N++;  /* keep th->N precisely up to date; p7_tophits_Destroy() must work w/ unfinished <th> on error  */
    }

  /* Put the hit[] "sorted" array in a random order w/ a Fisher-Yates shuffle */
  for (h = 0; h < th->N; h++)
    th->hit[h] = &(th->unsrt[h]);
  for (n = th->N; n > 1; n--) {	
    h = esl_rnd_Roll(rng, n);
    ESL_SWAP( th->hit[h], th->hit[n-1], P7_HIT *);
  }

  *ret_th = th;
  return eslOK;

 ERROR:
  if (th) p7_tophits_Destroy(th);
  *ret_th = NULL;
  return status;
}


int
p7_tophits_Validate(const P7_TOPHITS *th, char *errbuf)
{
  int i;
  int idx;
  int status;

  if (th->is_sorted_by_sortkey || th->is_sorted_by_seqidx)
    {
      for (i = 0; i < th->N; i++)
	{
	  idx = th->hit[i] - th->unsrt; /* i.e., by ptr arithmetic: #i in sorted list is #idx in unsorted list */
	  if (idx < 0 || idx >= th->N) ESL_FAIL(eslFAIL, errbuf, "sorted hit number %d points to bad address", i);
	  /* TestSample() currently doesn't sort its sampled hit array, so we don't test for proper sortedness */
	}
    }
  if (th->nreported < 0 || th->nreported > th->N) ESL_FAIL(eslFAIL, errbuf, "bad nreported field");
  if (th->nincluded < 0 || th->nincluded > th->N) ESL_FAIL(eslFAIL, errbuf, "bad nreported field");
  if (th->is_sorted_by_sortkey && th->is_sorted_by_seqidx) ESL_FAIL(eslFAIL, errbuf, "both sort type flags are up");
  for (i = 0; i < th->N; i++)
    if (( status = p7_hit_Validate( &(th->unsrt[i]), errbuf)) != eslOK) return status;

  return eslOK;
}

int
p7_tophits_Compare(const P7_TOPHITS *th1, const P7_TOPHITS *th2, float tol)
{
  int i;
  int idx1, idx2;
  int status;

  if (th1->N                    != th2->N)                    return eslFAIL;
  if (th1->nreported            != th2->nreported)            return eslFAIL;
  if (th1->nincluded            != th2->nincluded)            return eslFAIL;
  if (th1->is_sorted_by_sortkey != th2->is_sorted_by_sortkey) return eslFAIL;
  if (th1->is_sorted_by_seqidx  != th2->is_sorted_by_seqidx)  return eslFAIL;
  
  for (i = 0; i < th1->N; i++)
    if (( status = p7_hit_Compare( &(th1->unsrt[i]), &(th2->unsrt[i]), tol)) != eslOK) return status;

  if (th1->is_sorted_by_sortkey || th1->is_sorted_by_seqidx)
    {
      for (i = 0; i < th1->N; i++)
	{
	  idx1 = th1->hit[i] - th1->unsrt; /* i.e., by ptr arithmetic: #i in sorted list is #idx in unsorted list */
	  idx2 = th2->hit[i] - th2->unsrt; /* i.e., by ptr arithmetic: #i in sorted list is #idx in unsorted list */
	  if (idx1 != idx2) return eslFAIL;
	}
    }
  return eslOK;
}


/* Function:  p7_hit_TestSample()
 * Synopsis:  Sample a random, bogus, mostly syntactic P7_HIT.
 *
 * Purpose:   Sample a random but syntactically valid <P7_HIT>
 *            array, using random number generator <rng>, and 
 *            store it in <hit>, space provided by the caller
 *            (usually, one <P7_HIT> in an array that the caller
 *            has).
 */
int
p7_hit_TestSample(ESL_RANDOMNESS *rng, P7_HIT *hit)
{
  int d;
  int status;

  if ((status = esl_rsq_Sample(rng, eslRSQ_SAMPLE_GRAPH, 1+esl_rnd_Roll(rng, 30),  &(hit->name))) != eslOK) goto ERROR;
  if (esl_rnd_Roll(rng, 2)) { if ((status = esl_rsq_Sample(rng, eslRSQ_SAMPLE_ALNUM, 1+esl_rnd_Roll(rng, 10),  &(hit->acc)))  != eslOK) goto ERROR; }
  if (esl_rnd_Roll(rng, 2)) { if ((status = esl_rsq_Sample(rng, eslRSQ_SAMPLE_PRINT, 1+esl_rnd_Roll(rng, 120), &(hit->desc))) != eslOK) goto ERROR; }

  hit->window_length = 1 + esl_rnd_Roll(rng, 100000);
  hit->sortkey       = -1000. + 2000. * esl_random(rng);
  hit->score         = -1000. + 2000. * esl_random(rng);
  hit->pre_score     = -1000. + 2000. * esl_random(rng);
  hit->sum_score     = -1000. + 2000. * esl_random(rng);
  hit->lnP           = -1000. + 2000. * esl_random(rng);
  hit->pre_lnP       = -1000. + 2000. * esl_random(rng);
  hit->sum_lnP       = -1000. + 2000. * esl_random(rng);
  hit->ndom          = 1 + esl_rnd_Roll(rng, 10);
  hit->noverlaps     = esl_rnd_Roll(rng, hit->ndom);
  hit->nexpected     = esl_random(rng)*10;
  hit->flags         = p7_HITFLAGS_DEFAULT;
  if (esl_rnd_Roll(rng, 2)) hit->flags |= p7_IS_INCLUDED;
  if (esl_rnd_Roll(rng, 2)) hit->flags |= p7_IS_REPORTED;
  if (esl_rnd_Roll(rng, 2)) hit->flags |= p7_IS_NEW;
  if (esl_rnd_Roll(rng, 2)) hit->flags |= p7_IS_DROPPED;
  if (esl_rnd_Roll(rng, 2)) hit->flags |= p7_IS_DUPLICATE;
  hit->nreported     = 1 + esl_rnd_Roll(rng, hit->ndom);  
  hit->nincluded     = 1 + esl_rnd_Roll(rng, hit->ndom);
  hit->best_domain   =     esl_rnd_Roll(rng, hit->ndom);
  hit->seqidx        = 1 + esl_rnd_Roll(rng, 1000000);
  hit->subseq_start  = 1 + esl_rnd_Roll(rng, 1000000);
  hit->offset        = 1 + esl_rnd_Roll(rng, 1000000);

  if (( hit->dcl = p7_domain_Create(hit->ndom) ) == NULL) { status = eslEMEM; goto ERROR; }
  for (d = 0; d < hit->ndom; d++)
    if (( status = p7_domain_TestSample(rng, 1 + esl_rnd_Roll(rng, 100), &(hit->dcl[d]))) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  /* should free inside hit; caller has the shell of it though */
  return status;
}

int
p7_hit_Validate(const P7_HIT *hit, char *errbuf)
{
  int d;
  int status;

  if (hit->name == NULL) ESL_FAIL(eslFAIL, errbuf, "name cannot be NULL");
  if (isnan(hit->sortkey) ||
      isnan(hit->score) ||
      isnan(hit->pre_score) ||
      isnan(hit->sum_score) ||
      isnan(hit->lnP) ||
      isnan(hit->pre_lnP) ||
      isnan(hit->sum_lnP) ||
      isnan(hit->nexpected)) ESL_FAIL(eslFAIL, errbuf, "NaN found");
  
  if ( (hit->flags & (! ( p7_IS_REPORTED | p7_IS_INCLUDED | p7_IS_NEW | p7_IS_DROPPED | p7_IS_DUPLICATE))) != 0)
    ESL_FAIL(eslFAIL, errbuf, "unrecognized flag is up");
  
  if (hit->ndom < 0)  ESL_FAIL(eslFAIL, errbuf, "negative ndom");
  if (hit->noverlaps   < 0 || hit->noverlaps   >  hit->ndom) ESL_FAIL(eslFAIL, errbuf, "bad noverlaps");
  if (hit->nreported   < 0 || hit->nreported   >  hit->ndom) ESL_FAIL(eslFAIL, errbuf, "bad nreported");
  if (hit->nincluded   < 0 || hit->nincluded   >  hit->ndom) ESL_FAIL(eslFAIL, errbuf, "bad nincluded");
  if (hit->best_domain < 0 || hit->best_domain >= hit->ndom) ESL_FAIL(eslFAIL, errbuf, "bad best_domain");

  for (d = 0; d < hit->ndom; d++)
    if (( status = p7_domain_Validate(&(hit->dcl[d]), errbuf)) != eslOK) return status;
  
  return eslOK;
}

int
p7_hit_Compare(const P7_HIT *h1, const P7_HIT *h2, float tol)
{
  int d;
  int status;

  if (    strcmp(h1->name, h2->name) != 0) return eslFAIL;
  if (esl_strcmp(h1->acc,  h2->acc)  != 0) return eslFAIL;
  if (esl_strcmp(h1->desc, h2->desc) != 0) return eslFAIL;

  if ( h1->window_length != h2->window_length) return eslFAIL;
  if ( h1->ndom          != h2->ndom)          return eslFAIL;
  if ( h1->noverlaps     != h2->noverlaps)     return eslFAIL;
  if ( h1->flags         != h2->flags)         return eslFAIL;
  if ( h1->nreported     != h2->nreported)     return eslFAIL;
  if ( h1->nincluded     != h2->nincluded)     return eslFAIL;
  if ( h1->best_domain   != h2->best_domain)   return eslFAIL;
  if ( h1->seqidx        != h2->seqidx)        return eslFAIL;
  if ( h1->subseq_start  != h2->subseq_start)  return eslFAIL;
  if ( h1->offset        != h2->offset)        return eslFAIL;
  
  if ( esl_DCompare( h1->sortkey,   h2->sortkey,   tol ) != eslOK) return eslFAIL;
  if ( esl_FCompare( h1->score,     h2->score,     tol ) != eslOK) return eslFAIL;
  if ( esl_FCompare( h1->pre_score, h2->pre_score, tol ) != eslOK) return eslFAIL;
  if ( esl_FCompare( h1->sum_score, h2->sum_score, tol ) != eslOK) return eslFAIL;
  if ( esl_DCompare( h1->lnP,       h2->lnP,       tol ) != eslOK) return eslFAIL;
  if ( esl_DCompare( h1->pre_lnP,   h2->pre_lnP,   tol ) != eslOK) return eslFAIL;
  if ( esl_DCompare( h1->sum_lnP,   h2->sum_lnP,   tol ) != eslOK) return eslFAIL;
  if ( esl_DCompare( h1->nexpected, h2->nexpected, tol ) != eslOK) return eslFAIL;

  for (d = 0; d < h1->ndom; d++)
    if (( status = p7_domain_Compare(&(h1->dcl[d]), &(h2->dcl[d]), tol)) != eslOK) return status;
  return eslOK;
}


/*-------------- end, debug/devel tools ------------------------*/




/*****************************************************************
 * 4. Benchmark driver
 *****************************************************************/
#ifdef p7TOPHITS_BENCHMARK
/* 
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
  ESL_GETOPTS    *go       = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_STOPWATCH  *w        = esl_stopwatch_Create();
  ESL_RANDOMNESS *r        = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             N        = esl_opt_GetInteger(go, "-N");
  int             M        = esl_opt_GetInteger(go, "-M");
  P7_TOPHITS    **h        = NULL;
  P7_HIT         *hit      = NULL;
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
      h[j] = p7_tophits_Create(p7_TOPHITS_DEFAULT_INIT_ALLOC);
      for (i = 0; i < N; i++)
	{
	  p7_tophits_CreateNextHit(h[j], &hit);
	  esl_strdup(name, -1, &(hit->name));
	  esl_strdup(acc,  -1, &(hit->acc));
	  esl_strdup(desc, -1, &(hit->desc));
	  hit->sortkey    = sortkeys[j*N + i];
	  hit->score      = (float) sortkeys[j*N+i]; 
	  hit->pre_score  = 0.0;
	  hit->sum_score  = 0.0;
	  hit->lnP        = sortkeys[j*N+i];
	  hit->pre_lnP    = 0.0;
	  hit->sum_lnP    = 0.0;
	  hit->ndom       = N;
	  hit->noverlaps  = 0;
	  hit->nexpected  = 0;
	  hit->flags      = 0;
	  hit->nreported  = 0;
	  hit->nincluded  = 0;
	  hit->best_domain = 0;
	  hit->dcl         = NULL;
	}
      p7_tophits_SortBySortkey(h[j]);
  }
  /* then merge them into one big list in h[0] */
  for (j = 1; j < M; j++)
  {
      p7_tophits_Merge(h[0], h[j]);
      p7_tophits_Destroy(h[j]);
  }

  esl_stopwatch_Stop(w);

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

/* this is an obsolete derivative of old p7_tophits_Add(), which the
 * test driver happens to use; we should refactor.
 * The unit test only needs name and sort key;
 * we add the acc and desc too; and everything else is
 * at a default unused value.
 */
static int
tophits_Add(P7_TOPHITS *h, char *name, char *acc, char *desc, double sortkey)
{
  P7_HIT *hit = NULL;

  p7_tophits_CreateNextHit(h, &hit);
  esl_strdup(name, -1, &(hit->name));
  esl_strdup(acc,  -1, &(hit->acc));
  esl_strdup(desc, -1, &(hit->desc));
  hit->sortkey    = sortkey;
  return eslOK;
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
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

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  h1 = p7_tophits_Create(p7_TOPHITS_DEFAULT_INIT_ALLOC);
  h2 = p7_tophits_Create(p7_TOPHITS_DEFAULT_INIT_ALLOC);
  h3 = p7_tophits_Create(p7_TOPHITS_DEFAULT_INIT_ALLOC);
  
  for (i = 0; i < N; i++) 
  {
      key = esl_random(r);
      tophits_Add(h1, name, acc, desc, key);
      key = 10.0 * esl_random(r);
      tophits_Add(h2, name, acc, desc, key);
      key = 0.1 * esl_random(r);
      tophits_Add(h3, name, acc, desc, key);
  }
  tophits_Add(h1, "last",  NULL, NULL, -1.0);
  tophits_Add(h1, "first", NULL, NULL, 20.0);

  p7_tophits_SortBySortkey(h1);
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

  fprintf(stderr, "#  status = ok\n");
  return eslOK;
}
#endif /*p7TOPHITS_TESTDRIVE*/





