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
 * FullSortTophits()    - sorts the top H hits.
 * 
 ***************************************************************************** 
 * Brief example of use:
 *
 *   struct tophit_s   *yourhits;   // list of hits
 *   struct fancyali_s *ali;        // (optional structure) alignment of a hit 
 *   
 *   yourhits = AllocTophits(200);
 *   (for every hit in a search) {
 *        if (do_alignments) 
 *           ali = Trace2FancyAli();  // You provide a function/structure here
 *        if (score > threshold)
 *           RegisterHit(yourhits, ...)
 *     }
 *      
 *   FullSortTophits(yourhits);      // Sort hits by evalue 
 *   for (i = 0; i < 100; i++)       // Recover hits out in ranked order
 *     {   
 *       GetRankedHit(yourhits, i, ...);
 *                                   // Presumably you'd print here...
 *     } 
 *   FreeTophits(yourhits);
 ***************************************************************************   
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
#include <limits.h>

#include "structs.h"
#include "funcs.h"

/* Function: AllocTophits()
 * 
 * Purpose:  Allocate a struct tophit_s, for maintaining
 *           a list of top-scoring hits in a database search.
 *           
 * Args:     lumpsize - allocation lumpsize
 *           
 * Return:   An allocated struct hit_s. Caller must free.
 */
struct tophit_s *
AllocTophits(int lumpsize)
{
  struct tophit_s *hitlist;
  
  hitlist        = MallocOrDie (sizeof(struct tophit_s));
  hitlist->hit   = NULL;
  hitlist->unsrt = MallocOrDie (lumpsize * sizeof(struct hit_s));
  hitlist->alloc = lumpsize;
  hitlist->num   = 0;
  hitlist->lump  = lumpsize; 
  return hitlist;
}
void
GrowTophits(struct tophit_s *h)
{
  h->unsrt = ReallocOrDie(h->unsrt,(h->alloc + h->lump) * sizeof(struct hit_s));
  h->alloc += h->lump;
}
void
FreeTophits(struct tophit_s *h)
{
  int pos;
  for (pos = 0; pos < h->num; pos++)
    {
      if (h->unsrt[pos].ali  != NULL) FreeFancyAli(h->unsrt[pos].ali);
      if (h->unsrt[pos].name != NULL) free(h->unsrt[pos].name);
      if (h->unsrt[pos].desc != NULL) free(h->unsrt[pos].desc);
    }
  free(h->unsrt);
  if (h->hit != NULL) free(h->hit);
  free(h);
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
 * Purpose:  Add a new hit to a list of top hits.
 *
 *           "ali", if provided, is a pointer to allocated memory
 *           for an alignment output structure.
 *           Management is turned over to the top hits structure.
 *           Caller should not free them; they will be free'd by
 *           the FreeTophits() call. 
 *
 *           In contrast, "name" and "desc" are copied, so caller
 *           is still responsible for these.
 *           
 *           Number of args is unwieldy.
 *           
 * Args:     h        - active top hit list
 *           key      - value to sort by: bigger is better
 *           pvalue   - P-value of this hit 
 *           score    - score of this hit
 *           motherp  - P-value of parent whole sequence
 *           mothersc - score of parent whole sequence 
 *           name     - name of target sequence 
 *           desc     - description of target sequence 
 *           sqfrom   - 1..L pos in target seq  of start
 *           sqto     - 1..L pos; sqfrom > sqto if rev comp
 *           sqlen    - length of sequence, L
 *           hmmfrom  - 0..M+1 pos in HMM of start
 *           hmmto    - 0..M+1 pos in HMM of end
 *           hmmlen   - length of HMM, M
 *           domidx   - number of this domain 
 *           ndom     - total # of domains in sequence
 *           ali      - optional printable alignment info
 *           
 * Return:   (void)
 *           hitlist is modified and possibly reallocated internally.
 */
void
RegisterHit(struct tophit_s *h, double key, 
	    double pvalue, float score, double motherp, float mothersc,
	    char *name, char *desc, 
	    int sqfrom, int sqto, int sqlen,
	    int hmmfrom, int hmmto, int hmmlen, 
	    int domidx, int ndom,
	    struct fancyali_s *ali)
{
  /* Check to see if list is full and we must realloc.
   */
  if (h->num == h->alloc) GrowTophits(h);

  h->unsrt[h->num].name    = Strdup(name);
  h->unsrt[h->num].desc    = Strdup(desc);
  h->unsrt[h->num].sortkey = key;
  h->unsrt[h->num].pvalue  = pvalue;
  h->unsrt[h->num].score   = score;
  h->unsrt[h->num].motherp = motherp;
  h->unsrt[h->num].mothersc= mothersc;
  h->unsrt[h->num].sqfrom  = sqfrom;
  h->unsrt[h->num].sqto    = sqto;
  h->unsrt[h->num].sqlen   = sqlen;
  h->unsrt[h->num].hmmfrom = hmmfrom;
  h->unsrt[h->num].hmmto   = hmmto;
  h->unsrt[h->num].hmmlen  = hmmlen;
  h->unsrt[h->num].domidx  = domidx;
  h->unsrt[h->num].ndom    = ndom;
  h->unsrt[h->num].ali     = ali;
  h->num++; 
  return;
}

/* Function: GetRankedHit()
 * Date:     SRE, Tue Oct 28 10:06:48 1997 [Newton Institute, Cambridge UK]
 * 
 * Purpose:  Recover the data from the i'th ranked hit.
 *           Any of the data ptrs may be passed as NULL for fields
 *           you don't want. hitlist must have been sorted first.
 *           
 *           name, desc, and ali are returned as pointers, not copies;
 *           don't free them!
 */ 
void
GetRankedHit(struct tophit_s *h, int rank, 
	     double *r_pvalue, float *r_score, 
	     double *r_motherp, float *r_mothersc,
	     char **r_name, char **r_desc,
	     int *r_sqfrom, int *r_sqto, int *r_sqlen,
	     int *r_hmmfrom, int *r_hmmto, int *r_hmmlen,
	     int *r_domidx, int *r_ndom,
	     struct fancyali_s **r_ali)
{
  if (r_pvalue  != NULL) *r_pvalue  = h->hit[rank]->pvalue;
  if (r_score   != NULL) *r_score   = h->hit[rank]->score;
  if (r_motherp != NULL) *r_motherp = h->hit[rank]->motherp;
  if (r_mothersc!= NULL) *r_mothersc= h->hit[rank]->mothersc;
  if (r_name    != NULL) *r_name    = h->hit[rank]->name;
  if (r_desc    != NULL) *r_desc    = h->hit[rank]->desc;
  if (r_sqfrom  != NULL) *r_sqfrom  = h->hit[rank]->sqfrom;
  if (r_sqto    != NULL) *r_sqto    = h->hit[rank]->sqto;
  if (r_sqlen   != NULL) *r_sqlen   = h->hit[rank]->sqlen;
  if (r_hmmfrom != NULL) *r_hmmfrom = h->hit[rank]->hmmfrom;
  if (r_hmmto   != NULL) *r_hmmto   = h->hit[rank]->hmmto;
  if (r_hmmlen  != NULL) *r_hmmlen  = h->hit[rank]->hmmlen;
  if (r_domidx  != NULL) *r_domidx  = h->hit[rank]->domidx;
  if (r_ndom    != NULL) *r_ndom    = h->hit[rank]->ndom;
  if (r_ali     != NULL) *r_ali     = h->hit[rank]->ali;
}

/* Function: TophitsMaxName()
 * 
 * Purpose:  Returns the maximum name length in a top hits list; 
 *           doesn't need to be sorted yet.
 */
int
TophitsMaxName(struct tophit_s *h)
{
  int i;
  int len, maxlen;

  maxlen = 0;
  for (i = 0; i < h->num; i++)
    {
      len = strlen(h->unsrt[i].name);
      if (len > maxlen) maxlen = len;
    }
  return maxlen;
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
           /* don't ask. don't change. and, Don't Panic. */
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
  int i;

  /* If we don't have /any/ hits, then don't
   * bother.
   */
  if (h->num == 0) return;

  /* Assign the ptrs in h->hit.
   */
  h->hit = MallocOrDie(h->num * sizeof(struct hit_s *));
  for (i = 0; i < h->num; i++)
    h->hit[i] = &(h->unsrt[i]);

  /* Sort the pointers. Don't bother if we've only got one.
   */
  if (h->num > 1)
    qsort(h->hit, h->num, sizeof(struct hit_s *), hit_comparison);
}



/* Function: TophitsReport()
 * Date:     Thu Dec 18 13:19:18 1997
 * 
 * Purpose:  Generate a printout summarizing how much
 *           memory is used by a tophits structure,
 *           how many hits are stored, and how much
 *           waste there is from not knowing nseqs.
 *           
 * Args:     h    - the sorted tophits list
 *           E    - the cutoff in Evalue      
 *           nseq - the final number of seqs used for Eval 
 *           
 * Return:   (void)
 *           Prints information on stdout
 */           
void
TophitsReport(struct tophit_s *h, double E, int nseq)
{
  int i;
  int memused;
  int x;
  int n;

  /* Count up how much memory is used
   * in the whole list.
   */
  memused = sizeof(struct hit_s) * h->alloc + sizeof(struct tophit_s);
  for (i = 0; i < h->num; i++)
    {
      if (h->unsrt[i].name != NULL)
	memused += strlen(h->unsrt[i].name) + 1;
      if (h->unsrt[i].desc != NULL)
	memused += strlen(h->unsrt[i].desc) + 1;
      if (h->unsrt[i].ali != NULL)
	{
	  memused += sizeof(struct fancyali_s);
	  x = 0;
	  if (h->unsrt[i].ali->rfline != NULL) x++;
	  if (h->unsrt[i].ali->csline != NULL) x++;
  	  if (h->unsrt[i].ali->model  != NULL) x++;
	  if (h->unsrt[i].ali->mline  != NULL) x++;
	  if (h->unsrt[i].ali->aseq   != NULL) x++;
	  memused += x * (h->unsrt[i].ali->len + 1);
	  
	  if (h->unsrt[i].ali->query  != NULL) 
	    memused += strlen(h->unsrt[i].ali->query) + 1;
	  if (h->unsrt[i].ali->target != NULL) 
	    memused += strlen(h->unsrt[i].ali->target) + 1;
	}
    }

  /* Count how many hits actually satisfy the E cutoff.
   */
  n = 0;
  for (i = 0; i < h->num; i++)
    {
      if (h->hit[i]->pvalue * (double) nseq >= E) break;
      n++;
    }

  /* Format and print a summary
   */
  printf("tophits_s report:\n");
  printf("     Total hits:           %d\n", h->num);
  printf("     Satisfying E cutoff:  %d\n", n);
  printf("     Total memory:         %dK\n", memused / 1000);
}
