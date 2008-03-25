/* P7_TOPHITS: implementation of ranked list of top-scoring hits
 * 
 * Contents:
 *    1. The P7_TOPHITS object.
 *    2. Benchmark driver.
 *    3. Test driver.
 *    4. Copyright and license information.
 * 
 * SRE, Fri Dec 28 07:14:54 2007 [Janelia] [Enigma, MCMXC a.D.]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdlib.h>
#include <string.h>

#include "easel.h"

#include "hmmer.h"


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


/* Function:  
 * Synopsis:  
 * Incept:    SRE, Tue Mar 11 08:44:53 2008 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
p7_tophits_CreateNextHit(P7_TOPHITS *h, P7_HIT **ret_hit)
{
  P7_HIT *hit;
  int     status;

  if ((status = p7_tophits_Grow(h)) != eslOK) return status;
  
  hit = &(h->unsrt[h->N]);
  h->N++;
  if (h->N >= 2) h->is_sorted = FALSE;

  hit->name       = NULL;
  hit->acc        = NULL;
  hit->desc       = NULL;
  hit->sortkey    = 0.0;

  hit->score      = 0.0;
  hit->pre_score  = 0.0;
  hit->sum_score  = 0.0;

  hit->pvalue     = 0.0;
  hit->pre_pvalue = 0.0;
  hit->sum_pvalue = 0.0;

  hit->ndom       = 0;
  hit->nexpected  = 0.0;
  hit->nregions   = 0;
  hit->nclustered = 0;
  hit->noverlaps  = 0;
  hit->nenvelopes = 0;

  *ret_hit = hit;
  return eslOK;
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

/* Function:  p7_tophits_Destroy()
 * Synopsis:  Frees a hit list.
 * Incept:    SRE, Fri Dec 28 07:33:21 2007 [Janelia]
 */
void
p7_tophits_Destroy(P7_TOPHITS *h)
{
  int i;
  if (h == NULL) return;
  if (h->hit   != NULL) free(h->hit);
  if (h->unsrt != NULL) 
    {
      for (i = 0; i < h->N; i++)
	{
	  if (h->unsrt[i].name != NULL) free(h->unsrt[i].name);
	  if (h->unsrt[i].acc  != NULL) free(h->unsrt[i].acc);
	  if (h->unsrt[i].desc != NULL) free(h->unsrt[i].desc);
	}
      free(h->unsrt);
    }
  free(h);
  return;
}


/*****************************************************************
 * 2. Benchmark driver
 *****************************************************************/
#ifdef p7TOPHITS_BENCHMARK
/* 
  gcc -o benchmark-tophits -std=gnu99 -g -O2 -I. -L. -I../easel -L../easel -Dp7TOPHITS_BENCHMARK p7_tophits.c -lhmmer -leasel -lm 
  ./benchmark-tophits

  As of 28 Dec 07, shows 0.20u for 10 lists of 10,000 hits each (at least ~100x normal expectation),
  so we expect top hits list time to be negligible for typical hmmsearch/hmmpfam runs.
  
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
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                  0 },
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
  ESL_RANDOMNESS *r        = NULL;
  int             N        = esl_opt_GetInteger(go, "-N");
  int             M        = esl_opt_GetInteger(go, "-M");
  P7_TOPHITS    **h        = NULL;
  double         *sortkeys = NULL;
  char            name[]   = "not_unique_name";
  char            acc[]    = "not_unique_acc";
  char            desc[]   = "Test description for the purposes of making the benchmark allocate space";
  int             i,j;
  int             status;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

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
 * 3. Test driver
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
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                  0 },
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
  ESL_RANDOMNESS *r        = NULL;
  int             N        = esl_opt_GetInteger(go, "-N");
  P7_TOPHITS     *h1       = NULL;
  P7_TOPHITS     *h2       = NULL;
  P7_TOPHITS     *h3       = NULL;
  char            name[]   = "not_unique_name";
  char            acc[]    = "not_unique_acc";
  char            desc[]   = "Test description for the purposes of making the test driver allocate space";
  double          key;
  int             i;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

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





