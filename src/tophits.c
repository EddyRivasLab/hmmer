/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-1997 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the 
 *   GNU General Public License. See the files COPYING and 
 *   GNULICENSE for details.
 *    
 ************************************************************/

/* tophits.c
 * 
 * Routines for storing, sorting, displaying high scoring hits
 * and alignments.
 * 
 *****************************************************************************
 *
 * main API:
 * 
 * AllocTophits()       - allocation
 * FreeTophits()        - free'ing
 * RegisterHit()        - put information about a hit in the list
 * GetRankedHit()       - recovers information about a hit
 * FastSortTophits()    - puts top H hits in unsorted order in top H slots
 * FullSortTophits()    - sorts the top H hits.
 * 
 ***************************************************************************** 
 * Brief example of use:
 *
 *   struct tophit_s   *yourhits;   // list of hits
 *   struct fancyali_s *ali;        // (optional structure) alignment of a hit 
 *   int H = 1000;                  // max number of hits saved as scores
 *   int A = 100;                   // max number of alignments saved   
 *   
 *   yourhits = AllocTophits(H, A);
 *   (for every hit in a search) {
 *        if (do_alignments) 
 *           ali = Trace2FancyAli();  // You provide a function/structure here
 *        if (score > threshold)
 *           RegisterHit(yourhits, evalue, bitscore, name, desc, 
 *                       sqfrom, sqto, hmmfrom, hmmto, ali);
 *    }                   
 *      
 *   FastSortTophits(yourhits);      // Sort hits by evalue 
 *   for (i = 0; i < 100; i++)       // Recover hits out in ranked order
 *     {   
 *       GetRankedHit(yourhits, i, &evalue, &bitscore, &name, &desc,
 *                    &sqfrom, &sqto, &qfrom, &qto, &ali);
 *                                   // Presumably you'd print here...
 *     } 
 *   FreeTophits(yourhits);
 ***************************************************************************   
 * 
 * Output is controlled by the following parameters:
 *     T = threshold score (>=, bits)
 *     E = threshold score (<=, E-value)
 *     H = maximum number of hits kept as score and coords
 *     A = maximum number of hits kept as score, coords, and alignments
 *     
 * Estimated storage per hit: 
 *        coords:   16 bytes
 *        scores:    8 bytes
 *     name/desc:  192 bytes
 *     alignment: 1000 bytes   total = ~1200 bytes with alignment;
 *                                   = ~200 bytes without  
 *     Designed for: 10^5 hits (20 MB) or 10^4 alignments (10 MB)  
 */

#include <string.h>
#include <float.h>

#include "structs.h"
#include "funcs.h"


static void free_loser_alignments(struct tophit_s *h);


/* Function: AllocTophits()
 * 
 * Purpose:  Allocate a struct tophit_s, for maintaining
 *           a list of up to H top-scoring hits in a database search
 *           and up to A top-scoring hits with alignment info.
 *           
 * Args:     H - number of top hits to save          
 *           A - number of top hits to save with alignments
 *           
 * Return:   An allocated struct hit_s. Caller must free.
 */
struct tophit_s *
AllocTophits(int H, int A)
{
  int i;
  struct tophit_s *hitlist;
  
  hitlist = (struct tophit_s *) MallocOrDie (sizeof(struct tophit_s));
  hitlist->H     = H;
  hitlist->A     = A;
  hitlist->alloc = H * 2;         /* overallocation; memory/speed tradeoff */
  hitlist->sorts = 0;		  /* counter for # of sorts                */
  hitlist->pos   = 0;
  hitlist->best  = -DBL_MAX;	  /* current worst sort key */

  hitlist->hit = (struct hit_s **) MallocOrDie (sizeof(struct hit_s *)*
						(hitlist->alloc));
  for (i = 0; i < hitlist->alloc; i++)
    hitlist->hit[i] = NULL;

  return hitlist;
}
void
FreeTophits(struct tophit_s *hitlist)
{
  int pos;

  for (pos = 0; pos < hitlist->alloc; pos++)
    if (hitlist->hit[pos] != NULL)
      {
	if (hitlist->hit[pos]->ali != NULL)
	  FreeFancyAli(hitlist->hit[pos]->ali);
	if (hitlist->hit[pos]->name != NULL)
	  free(hitlist->hit[pos]->name);
	if (hitlist->hit[pos]->desc != NULL)
	  free(hitlist->hit[pos]->desc);
	free(hitlist->hit[pos]);
      }
  free(hitlist->hit);
  free(hitlist);
}

struct fancyali_s *
AllocFancyAli(void)
{
  struct fancyali_s *ali;

  ali = MallocOrDie (sizeof(struct fancyali_s));
  ali->rfline = ali->csline = ali->model = ali->mline = ali->aseq = NULL;
  ali->query  = ali->target = NULL;
  ali->sqfrom = ali->sqto   = 0;
  return ali;
}
void
FreeFancyAli(struct fancyali_s *ali)
{
  if (ali != NULL) {
    if (ali->rfline != NULL) free(ali->rfline);
    if (ali->csline != NULL) free(ali->csline);
    if (ali->model  != NULL) free(ali->model);
    if (ali->mline  != NULL) free(ali->mline);
    if (ali->aseq   != NULL) free(ali->aseq);
    if (ali->query  != NULL) free(ali->query);
    if (ali->target != NULL) free(ali->target);
    free(ali);
  }
}
    



/* Function: RegisterHit()
 * 
 * Purpose:  Examine a new hit and possibly add it to the
 *           list of top hits.
 *
 *           "ali", if provided, is a pointer to allocated memory
 *           for an alignment output structure.
 *           Management is turned over to the top hits structure.
 *           Caller should not free them; they will be free'd by
 *           the FreeTophits() call. 
 *
 *           In contrast, "name" and "desc" are copied, so the caller
 *           is still responsible for these.
 *           
 * Args:     hitlist  - active top hit list
 *           key      - value to sort by: bigger is better
 *           evalue   - evalue (or pvalue) of this hit 
 *           score    - score of this hit
 *           name     - name of target sequence 
 *           desc     - description of target sequence 
 *           sqfrom   - 1..L pos in target seq  of start
 *           sqto     - 1..L pos; sqfrom > sqto if rev comp
 *           hmmfrom  - 0..M+1 pos in HMM of start
 *           hmmto    - 0..M+1 pos in HMM of end
 *           ali      - optional printable alignment info
 *           
 * Return:   (void)
 *           hitlist may be modified, and may be resorted.          
 */
void
RegisterHit(struct tophit_s *hitlist, double key, double evalue, float score,
	    char *name, char *desc, int sqfrom, int sqto, int sqlen,
	    int hmmfrom, int hmmto, int hmmlen, struct fancyali_s *ali)
{
  /* Check if we have to add this hit to the current list.
   */
  if (key <= hitlist->best) { if (ali != NULL) FreeFancyAli(ali); return; }

  /* Check to see if list is full and we need to re-sort it.
   */
  if (hitlist->pos == hitlist->alloc)
    {
      FastSortTophits(hitlist);
			/* not a safe call yet; top A must be sorted */
      /*      free_loser_alignments(hitlist); */
      hitlist->best = hitlist->hit[hitlist->H-1]->sortkey;
      hitlist->sorts++;
    }

  /* Add the current hit to the list.
   * For efficiency, we reuse old hit_s structures in the
   * list. If we're reusing an existing structure, we
   * first must free its trace
   */
  if (hitlist->hit[hitlist->pos] == NULL)
    hitlist->hit[hitlist->pos] =
      (struct hit_s *) MallocOrDie (sizeof(struct hit_s));
  else 
    {
      if (hitlist->hit[hitlist->pos]->ali != NULL)
	FreeFancyAli(hitlist->hit[hitlist->pos]->ali);
      if (hitlist->hit[hitlist->pos]->name != NULL)
	free(hitlist->hit[hitlist->pos]->name);
      if (hitlist->hit[hitlist->pos]->desc != NULL)
	free(hitlist->hit[hitlist->pos]->desc);
    }

  hitlist->hit[hitlist->pos]->name    = Strdup(name);
  hitlist->hit[hitlist->pos]->desc    = Strdup(desc);
  hitlist->hit[hitlist->pos]->sortkey = key;
  hitlist->hit[hitlist->pos]->evalue  = evalue;
  hitlist->hit[hitlist->pos]->score   = score;
  hitlist->hit[hitlist->pos]->sqfrom  = sqfrom;
  hitlist->hit[hitlist->pos]->sqto    = sqto;
  hitlist->hit[hitlist->pos]->sqlen   = sqlen;
  hitlist->hit[hitlist->pos]->hmmfrom = hmmfrom;
  hitlist->hit[hitlist->pos]->hmmto   = hmmto;
  hitlist->hit[hitlist->pos]->hmmlen  = hmmlen;
  hitlist->hit[hitlist->pos]->ali     = ali;

  hitlist->pos++; 
  return;
}

/* Function: GetRankedHit()
 * Date:     SRE, Tue Oct 28 10:06:48 1997 [Newton Institute, Cambridge UK]
 * 
 * Purpose:  Recover the data from the i'th ranked hit.
 *           Any of the data ptrs may be passed as NULL for fields
 *           you don't want.
 *           
 *           name, desc, and ali are returned as pointers, not copies;
 *           don't free them.
 */ 
void
GetRankedHit(struct tophit_s *h, int rank, 
	     double *r_evalue, float *r_score, char **r_name, char **r_desc,
	     int *r_sqfrom, int *r_sqto, int *r_sqlen,
	     int *r_hmmfrom, int *r_hmmto, int *r_hmmlen,
	     struct fancyali_s **r_ali)
{
  if (r_evalue  != NULL) *r_evalue  = h->hit[rank]->evalue;
  if (r_score   != NULL) *r_score   = h->hit[rank]->score;
  if (r_name    != NULL) *r_name    = h->hit[rank]->name;
  if (r_desc    != NULL) *r_desc    = h->hit[rank]->desc;
  if (r_sqfrom  != NULL) *r_sqfrom  = h->hit[rank]->sqfrom;
  if (r_sqto    != NULL) *r_sqto    = h->hit[rank]->sqto;
  if (r_sqlen   != NULL) *r_sqlen   = h->hit[rank]->sqlen;
  if (r_hmmfrom != NULL) *r_hmmfrom = h->hit[rank]->hmmfrom;
  if (r_hmmto   != NULL) *r_hmmto   = h->hit[rank]->hmmto;
  if (r_hmmlen  != NULL) *r_hmmlen  = h->hit[rank]->hmmlen;
  if (r_ali     != NULL) *r_ali     = h->hit[rank]->ali;
}


/* Function: TophitsMaxName()
 * 
 * Purpose:  Returns the maximum name length in a top hits list.
 */
int
TophitsMaxName(struct tophit_s *h, int top_howmany)
{
  int i;
  int len, maxlen;
  
  maxlen = 0;
  for (i = 0; i < top_howmany && i < h->pos; i++)
    {
      len = strlen(h->hit[i]->name);
      if (len > maxlen) maxlen = len;
    }
  return maxlen;
}



#ifdef SRE_REMOVED
/* Function: PrintTopHits()
 * 
 * Purpose:  Format and print top hits list to a file pointer.
 * 
 * Args:     fp     - where to print to
 *           h      - top hits list
 *           n      - number of hits to print
 *           evd    - TRUE if we have a valid EVD fit
 *           mu     - EVD mu parameter
 *           lambda - EVD lambda parameter
 *           dbseqs - number of seqs in database (for calculating E-value)
 */
void
PrintTopHits(FILE *fp, struct tophit_s *h, int n, int evd, float mu, float lambda, int dbseqs)
{
  int i;
  int nmlen;			/* maximum name length */
  int desclen;			/* maximum description length */

  nmlen   = TophitsMaxName(h, n);
  if (nmlen < 12) nmlen = 12;
  if (nmlen > 61) nmlen = 61;
  desclen = 79 - (nmlen + 18); 
  printf("%-63.63s %6s %8s\n", "Sequence", "Score", "P-value");
  printf("%-63.63s %6s %8s\n", "--------", "-----", "-------");

  for (i = 0; i < n && i < h->pos; i++)
    {
      if (evd) 
	fprintf(fp, "%-*.*s  %-*.*s %6.1f %8.1g\n",
		nmlen, nmlen,     h->hit[i]->name,
		desclen, desclen, h->hit[i]->desc,
		h->hit[i]->score,
		ExtremeValueP2(h->hit[i]->score, mu, lambda, dbseqs));
      else
	fprintf(fp, "%-*.*s  %-*.*s %6.1f %8s\n",
		nmlen, nmlen,     h->hit[i]->name,
		desclen, desclen, h->hit[i]->desc,
		h->hit[i]->score,
		"?");
    }
}
#endif /*REMOVED*/
	      

/* Function: FastSortTophits()
 * 
 * Purpose:  Re-sort top hits list so that the top H are 
 *           in positions 0..H-1 and position H-1 has the H'th
 *           ranking hit. NOTE THAT THIS IS *NOT* A COMPLETE
 *           SORT! The top H are unsorted!
 *           
 *           Sorts by largest "sortkey" element in tophits list.
 * 
 *           This is an independent (i.e. my own, GPL'ed) implementation
 *           based on the select() function of Numerical
 *           Recipes in C (p. 342 of 2nd Edition, 1992, Cambridge
 *           University Press).
 *           
 * Args:     h  - hit list to sort.
 *                
 * Return:   sorted hit list.
 */                         
void
FastSortTophits(struct tophit_s *h)
{
  int left, right;
  int i, j;
  struct hit_s *swapfoo;	/* must declare this to use SWAP() */

  /* if we have <= H elements in the list, we're
   * already done.
   */
  if (h->pos <= h->H) return;

  left   = 0;
  right  = h->pos-1;
  if (h->H < h->pos) h->pos = h->H;
  /*CONSTCOND*/ while (1)
    {
      if (right - left <= 0) return;   /* case of single elem */
      else if (right - left == 1)      /* case of two elem */
	{
	  if (h->hit[left]->sortkey < h->hit[right]->sortkey)
	    SWAP(h->hit[left], h->hit[right]);
	  return;
	}
      else			/* general case */
	{
	  /* The purpose of the first three swaps is to set up
	   * guard positions. We use the leftmost guy as our
	   * partitioner. Right is guaranteed to be bigger (worse) than
	   * the partitioner. Left+1 is guaranteed to be smaller (better) than
	   * the partitioner. This way i can't pass right, and
	   * j can't pass left. This also means that left+1
           * and right are sorted, so we can start our next swaps
           * on left+2,right-1.
	   */
	  if (h->hit[right]->sortkey > h->hit[left]->sortkey)
	    SWAP(h->hit[left], h->hit[right]);
	  if (h->hit[right]->sortkey > h->hit[left+1]->sortkey)
	    SWAP(h->hit[left+1], h->hit[right]);
	  if (h->hit[left]->sortkey > h->hit[left+1]->sortkey)
	    SWAP(h->hit[left+1], h->hit[left]);

	  i = left+1;
	  j = right;
	  /*CONSTCOND*/ while (1)
	    {			/* i finds a worse score; j finds a better one */
	      i++; while (h->hit[left]->sortkey < h->hit[i]->sortkey) i++;
	      j--; while (h->hit[left]->sortkey > h->hit[j]->sortkey) j--;
	      if (j <= i) break;
	      SWAP(h->hit[i], h->hit[j]);
	    }
	  SWAP(h->hit[left], h->hit[j]);
	  if (h->H-1 >= j) left  = j;
	  if (h->H-1 <= j) right = j-1;
	}
    }
  /*NOTREACHED*/
  return;
}

/* Function: FullSortTophits()
 * 
 * Purpose:  Completely sort the top hits list. Calls
 *           qsort() to do the sorting, and uses 
 *           hit_comparison() to do the comparison.
 *           
 * Args:     h - top hits structure
 */          
int
hit_comparison(const void *vh1, const void *vh2)
{
				/* don't ask. */
  struct hit_s *h1 = *((struct hit_s **) vh1);
  struct hit_s *h2 = *((struct hit_s **) vh2);
  
  if      (h1->sortkey < h2->sortkey)  return  1;
  else if (h1->sortkey > h2->sortkey)  return -1;
  else if (h1->sortkey == h2->sortkey) return  0;

  /*NOTREACHED*/
  return 0;
}
void
FullSortTophits(struct tophit_s *h)
{
  /* Note: sort only top h->pos hits, not full allocation;
   * else we might try to sort nonexistent hits
   */
  qsort(h->hit, h->pos, sizeof(struct hit_s *), hit_comparison);
}


/* Function: free_loser_alignments()
 * 
 * Purpose:  Frees all the alignments below the A limit.
 *           Called only after a sort by FastSortTophits()
 *           or FullSortTophits().
 */           
static void
free_loser_alignments(struct tophit_s *h)
{
  int i;

  for (i = h->A+1; i < h->alloc && i < h->pos; i++)
    if (h->hit[i]->ali != NULL)
      {
	FreeFancyAli(h->hit[i]->ali);
	h->hit[i]->ali = NULL;
      } 
}
