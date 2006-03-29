/*****************************************************************
 * @LICENSE@
 *****************************************************************/

/* modelmakers.c
 * SRE, Fri Nov 15 10:00:04 1996
 * 
 * Construction of models from multiple alignments. Three versions:
 *    Handmodelmaker() -- use #=RF annotation to indicate match columns
 *    Fastmodelmaker() -- Krogh/Haussler heuristic
 *    Maxmodelmaker()  -- MAP model construction algorithm (Eddy, 
 *                          unpublished)
 *                          
 * The meat of the model construction code is in matassign2hmm().
 * The three model construction strategies simply label which columns
 * are supposed to be match states, and then hand this info to
 * matassign2hmm().
 * 
 * Two wrinkles to watch for:
 * 1) The alignment is assumed to contain sequence fragments. Look in
 *    fake_tracebacks() for how internal entry/exit points are handled.
 * 2) Plan7 disallows DI and ID transitions, but an alignment may
 *    imply these. Look in trace_doctor() for how DI and ID transitions
 *    are removed.
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#include "squid.h"
#include "msa.h"
#include "structs.h"
#include "funcs.h"



/* flags used for matassign[] arrays -- 
 *   assignment of aligned columns to match/insert states
 */
#define ASSIGN_MATCH      (1<<0) 
#define FIRST_MATCH       (1<<1)
#define LAST_MATCH        (1<<2)
#define ASSIGN_INSERT     (1<<3)
#define EXTERNAL_INSERT_N (1<<4)
#define EXTERNAL_INSERT_C (1<<5) 

static int build_cij(char **aseqs, int nseq, int *insopt, int i, int j,
		     float *wgt, float *cij);
static int estimate_model_length(MSA *msa);
static void matassign2hmm(MSA *msa, char **dsq,
			  int *matassign, struct plan7_s **ret_hmm,
			  struct p7trace_s ***ret_tr);
static void fake_tracebacks(char **aseq, int nseq, int alen, int *matassign,
			    struct p7trace_s ***ret_tr);
static void trace_doctor(struct p7trace_s *tr, int M, int *ret_ndi, 
			 int *ret_nid);
static void annotate_model(struct plan7_s *hmm, int *matassign, MSA *msa);
static void print_matassign(int *matassign, int alen);



/* Function: P7Handmodelmaker()
 * 
 * Purpose:  Manual model construction:
 *           Construct an HMM from an alignment, where the #=RF line
 *           of a HMMER alignment file is given to indicate
 *           the columns assigned to matches vs. inserts.
 *           
 *           NOTE: Handmodelmaker() will slightly revise the alignment
 *           if necessary, if the assignment of columns implies
 *           DI and ID transitions.
 *           
 *           Returns both the HMM in counts form (ready for applying
 *           Dirichlet priors as the next step), and fake tracebacks
 *           for each aligned sequence.
 *           
 * Args:     msa  - multiple sequence alignment          
 *           dsq  - digitized unaligned aseq's
 *           ret_hmm - RETURN: counts-form HMM
 *           ret_tr  - RETURN: array of tracebacks for aseq's
 *           
 * Return:   (void)
 *           ret_hmm and ret_tr alloc'ed here; FreeTrace(tr[i]), free(tr),
 *           FreeHMM(hmm).
 */            
void
P7Handmodelmaker(MSA *msa, char **dsq,
		 struct plan7_s **ret_hmm, struct p7trace_s ***ret_tr)
{
  int     *matassign;           /* MAT state assignments if 1; 1..alen */
  int      apos;                /* counter for aligned columns         */

  /* Make sure we have all the info about the alignment that we need */
  if (msa->rf == NULL)
    Die("Alignment must have RF annotation to hand-build an HMM");

				/* Allocation */
  matassign = (int *) MallocOrDie (sizeof(int) * (msa->alen+1));
  
  /* Determine match assignment from optional annotation
   */
  matassign[0] = 0;
  for (apos = 0; apos < msa->alen; apos++)
    {
      matassign[apos+1] = 0;
      if (!isgap(msa->rf[apos])) 
	matassign[apos+1] |= ASSIGN_MATCH;
      else 
	matassign[apos+1] |= ASSIGN_INSERT;
    }

  /* Hand matassign off for remainder of model construction
   */
  /*   print_matassign(matassign, msa->alen); */
  matassign2hmm(msa, dsq, matassign, ret_hmm, ret_tr);

  free(matassign);
  return;
}


/* Function: P7Fastmodelmaker()
 * 
 * Purpose:  Heuristic model construction:
 *           Construct an HMM from an alignment using the original
 *           Krogh/Haussler heuristic; any column with more 
 *           symbols in it than a given fraction is assigned to
 *           match.
 *           
 *           NOTE: Fastmodelmaker() will slightly revise the
 *           alignment if the assignment of columns implies
 *           DI and ID transitions.
 *           
 *           Returns the HMM in counts form (ready for applying Dirichlet
 *           priors as the next step). Also returns fake traceback
 *           for each training sequence.
 *           
 * Args:     msa     - multiple sequence alignment
 *           dsq     - digitized unaligned aseq's
 *           maxgap  - if more gaps than this, column becomes insert.
 *           ret_hmm - RETURN: counts-form HMM
 *           ret_tr  - RETURN: array of tracebacks for aseq's
 *           
 * Return:   (void)
 *           ret_hmm and ret_tr alloc'ed here; FreeTrace(tr[i]), free(tr),
 *           FreeHMM(hmm).       
 */
void
P7Fastmodelmaker(MSA *msa, char **dsq, float maxgap,
		 struct plan7_s **ret_hmm, struct p7trace_s ***ret_tr)
{
  int     *matassign;           /* MAT state assignments if 1; 1..alen */
  int      idx;                 /* counter over sequences              */
  int      apos;                /* counter for aligned columns         */
  int      ngap;		/* number of gaps in a column          */

  /* Allocations: matassign is 1..alen array of bit flags
   */
  matassign = (int *) MallocOrDie (sizeof(int) * (msa->alen+1));
  
  /* Determine match assignment by counting symbols in columns
   */
  matassign[0] = 0;
  for (apos = 0; apos < msa->alen; apos++) {
      matassign[apos+1] = 0;

      ngap = 0;
      for (idx = 0; idx < msa->nseq; idx++) 
	if (isgap(msa->aseq[idx][apos])) 
	  ngap++;
      
      if ((float) ngap / (float) msa->nseq > maxgap)
	matassign[apos+1] |= ASSIGN_INSERT;
      else
	matassign[apos+1] |= ASSIGN_MATCH;
  }

  /* Once we have matassign calculated, all modelmakers behave
   * the same; matassign2hmm() does this stuff (traceback construction,
   * trace counting) and sets up ret_hmm and ret_tr.
   */
  matassign2hmm(msa, dsq, matassign, ret_hmm, ret_tr);

  free(matassign);
  return;
}
  

/* Function: P7Maxmodelmaker()
 * 
 * Purpose:  The Unholy Beast of HMM model construction algorithms --
 *           maximum a posteriori construction. A tour de force and
 *           probably overkill. MAP construction for Krogh 
 *           HMM-profiles is fairly straightforward, but MAP construction of
 *           Plan 7 HMM-profiles is, er, intricate.
 *           
 *           Given a multiple alignment, construct an optimal (MAP) model
 *           architecture. Return a counts-based HMM.
 *           
 * Args:     msa     - multiple sequence alignment 
 *           dsq     - digitized, unaligned seqs
 *           maxgap  - above this, trailing columns are assigned to C           
 *           prior   - priors on parameters to use for model construction
 *           null    - random sequence model emissions
 *           null_p1 - random sequence model p1 transition 
 *           mpri    - prior on architecture: probability of new match node
 *           ret_hmm - RETURN: new hmm (counts form)
 *           ret_tr  - RETURN: array of tracebacks for aseq's
 *
 * Return:   (void)
 *           ret_hmm and ret_tr (if !NULL) must be free'd by the caller.
 */          
void
P7Maxmodelmaker(MSA *msa, char **dsq, float maxgap,
		struct p7prior_s *prior, 
		float *null, float null_p1, float mpri, 
		struct plan7_s **ret_hmm, struct p7trace_s  ***ret_tr)
{
  int     idx;			/* counter for seqs                */
  int     i, j;			/* positions in alignment          */
  int     x;			/* counter for syms or transitions */
  float **matc;                 /* count vectors: [1..alen][0..19] */ 
  float   cij[8], tij[8];	/* count and score transit vectors */
  float   matp[MAXABET];        /* match emission vector           */
  float   insp[MAXABET];	/* insert score vector             */
  float   insc[MAXABET];	/* insert count vector             */
  float  *sc;                   /* DP scores [0,1..alen,alen+1]    */
  int    *tbck;                 /* traceback ptrs for sc           */
  int    *matassign;            /* match assignments [1..alen]     */ 
  int    *insopt;		/* number of inserted chars [0..nseq-1] */
  int     first, last;		/* positions of first and last cols [1..alen] */
  float   bm1, bm2;		/* estimates for start,internal b->m t's  */
  int     est_M;		/* estimate for the size of the model */
  float   t_me;			/* estimate for an internal M->E transition */
  float   new, bestsc;		/* new score, best score so far     */
  int     code;			/* optimization: return code from build_cij() */
  int     ngap;			/* gap count in a column                      */
  float   wgtsum;		/* sum of weights; do not assume it is nseq */

  /* Allocations
   */
  matc      = (float **) MallocOrDie (sizeof(float *) * (msa->alen+1));
  sc        = (float *)  MallocOrDie (sizeof(float)   * (msa->alen+2));
  tbck      = (int *)    MallocOrDie (sizeof(int)     * (msa->alen+2));
  matassign = (int *)    MallocOrDie (sizeof(int)     * (msa->alen+1));
  insopt    = (int *)    MallocOrDie (sizeof(int)     * msa->nseq);    
  for (i = 0; i < msa->alen; i++) {
    matc[i+1] = (float *) MallocOrDie (Alphabet_size * sizeof(float));
    FSet(matc[i+1], Alphabet_size, 0.);
  }

  /* Precalculations
   */
  for (i = 0; i < msa->alen; i++) 
    for (idx = 0; idx < msa->nseq; idx++) 
      if (!isgap(msa->aseq[idx][i])) 
	P7CountSymbol(matc[i+1], SymbolIndex(msa->aseq[idx][i]), msa->wgt[idx]);
  mpri = sreLOG2(mpri);
  
  FCopy(insp, prior->i[0], Alphabet_size);
  FNorm(insp, Alphabet_size);
  wgtsum = FSum(msa->wgt, msa->nseq);
  for (x = 0; x < Alphabet_size; x++)
    insp[x] = sreLOG2(insp[x] / null[x]);
  
  /* Estimate the relevant special transitions.
   */
  est_M = estimate_model_length(msa);
  t_me  = 0.5 / (float) (est_M-1);
  bm1   = 0.5;
  bm2   = 0.5 / (float) (est_M-1);
  bm1   = sreLOG2(bm1 / null_p1);
  bm2   = sreLOG2(bm2 / null_p1);

  /* Estimate the position of the last match-assigned column
   * by counting gap frequencies.
   */ 
  maxgap = 0.5;
  for (last = msa->alen; last >= 1; last--) {
    ngap = 0;
    for (idx = 0; idx < msa->nseq; idx++)
      if (isgap(msa->aseq[idx][last-1])) ngap++;
    if ((float) ngap / (float) msa->nseq <= maxgap)
      break;
  }

  /* Initialization
   */
  sc[last]   = 0.;
  tbck[last] = 0;

  /* Set ME gaps to '_'
   */
  for (idx = 0; idx < msa->nseq; idx++) 
    for (i = last; i > 0 && isgap(msa->aseq[idx][i-1]); i--)
      msa->aseq[idx][i-1] = '_';

  /* Main recursion moves from right to left.
   */
  for (i = last-1; i > 0; i--) {
				/* Calculate match emission scores for i  */
    FCopy(matp, matc[i], Alphabet_size);
    P7PriorifyEmissionVector(matp, prior, prior->mnum, prior->mq, prior->m, NULL); 
    for (x = 0; x < Alphabet_size; x++)
      matp[x] = sreLOG2(matp[x] / null[x]);

				/* Initialize insert counters to zero */
    FSet(insc, Alphabet_size, 0.);
    for (idx = 0; idx < msa->nseq; idx++) insopt[idx] = 0;

    sc[i] = -FLT_MAX; 
    for (j = i+1; j <= last; j++) {
			/* build transition matrix for column pair i,j */
      code = build_cij(msa->aseq, msa->nseq, insopt, i, j, msa->wgt, cij);
      if (code == -1) break;	/* no j to our right can work for us */
      if (code == 1) {
	FCopy(tij, cij, 7);
	P7PriorifyTransitionVector(tij, prior, prior->tq); 
	FNorm(tij, 3);
	tij[TMM] = sreLOG2(tij[TMM] / null_p1); 
	tij[TMI] = sreLOG2(tij[TMI] / null_p1); 
	tij[TMD] = sreLOG2(tij[TMD]); 
	tij[TIM] = sreLOG2(tij[TIM] / null_p1); 
	tij[TII] = sreLOG2(tij[TII] / null_p1); 
	tij[TDM] = sreLOG2(tij[TDM] / null_p1); 
	tij[TDD] = sreLOG2(tij[TDD]); 
  				/* calculate the score of using this j. */
	new = sc[j] +  FDot(tij, cij, 7) + FDot(insp, insc, Alphabet_size);

	SQD_DPRINTF2(("%3d %3d new=%6.2f scj=%6.2f m=%6.2f i=%6.2f t=%6.2f\n",
	       i, j, new, sc[j], FDot(matp, matc[i], Alphabet_size), 
	       FDot(insp, insc, Alphabet_size), FDot(tij, cij, 7)));

				/* keep it if it's better */
	if (new > sc[i]) {
	  sc[i]   = new;
	  tbck[i] = j;
	} 
      }
				/* bump insc, insopt insert symbol counters */
      FAdd(insc, matc[j], Alphabet_size);
      for (idx = 0; idx < msa->nseq; idx++)
	if (!isgap(msa->aseq[idx][j-1])) insopt[idx]++;
    }
				/* add in constant contributions for col i */
				/* note ad hoc scaling of mpri by wgtsum (us. nseq)*/
    sc[i] += FDot(matp, matc[i], Alphabet_size) + mpri * wgtsum;
  } /* end loop over start positions i */

  /* Termination: place the begin state.
   * log odds score for S->N->B is all zero except for NB transition, which
   * is a constant. So we only have to evaluate BM transitions.
   */
  bestsc = -FLT_MAX;
  first  = 0;
  for (i = 1; i <= last; i++) {
    new = sc[i]; 
    for (idx = 0; idx < msa->nseq; idx++) {
      if (isgap(msa->aseq[idx][i-1])) 
	new += bm2;		/* internal B->M transition */
      else
	new += bm1;		/* B->M1 transition         */
    }
    if (new > bestsc) {
      bestsc = new;
      first  = i;
    }
  }

  /* Traceback
   */
  matassign[0] = 0;
  for (i = 1; i <= msa->alen; i++) matassign[i] = ASSIGN_INSERT; 
  for (i = first; i != 0; i = tbck[i]) {
    matassign[i] &= ~ASSIGN_INSERT;
    matassign[i] |= ASSIGN_MATCH;
  }

  /* Hand matassign off for remainder of model construction
   */
  /*  print_matassign(matassign, ainfo->alen); */
  matassign2hmm(msa, dsq, matassign, ret_hmm, ret_tr);

  /* Clean up.
   */
  for (i = 1; i <= msa->alen; i++) free(matc[i]);
  free(matc);
  free(sc);
  free(tbck);
  free(matassign);
  free(insopt);
}


/* Function: build_cij()
 * 
 * Purpose:  Construct a counts vector for transitions between
 *           column i and column j in a multiple alignment.
 *
 *           '_' gap characters indicate "external" gaps which
 *           are to be dealt with by B->M and M->E transitions. 
 *           These characters must be placed by a preprocessor.
 *
 *           insopt is an "insert optimization" -- an incrementor
 *           which keeps track of the number of insert symbols
 *           between i and j.
 *       
 * Args:     aseqs  - multiple alignment. [0.nseq-1][0.alen-1]
 *           nseq   - number of seqs in aseqs
 *           insopt - number of inserts per seq between i/j [0.nseq-1]
 *           i      - i column [1.alen], off by one from aseqs
 *           j      - j column [1.alen], off by one from aseqs
 *           wgt    - per-seq weights [0.nseq-1]
 *           cij    - transition count vectors [0..7]
 *           
 * Return:   -1 if an illegal transition was seen for this i/j assignment *and*
 *              we are guaranteed that any j to the right will also
 *              have illegal transitions.
 *           0  if an illegal transition was seen, but a j further to the
 *              right may work.  
 *           1 if all transitions were legal.
 */          
static int 
build_cij(char **aseqs, int nseq, int *insopt, int i, int j,
          float *wgt, float *cij)
{
  int idx;			/* counter for seqs */

  i--;				/* make i,j relative to aseqs [0..alen-1] */
  j--;
  FSet(cij, 8, 0.);		/* zero cij */
  for (idx = 0; idx < nseq; idx++) {
    if (insopt[idx] > 0) {
      if (isgap(aseqs[idx][i])) return -1; /* D->I prohibited. */
      if (isgap(aseqs[idx][j])) return 0;  /* I->D prohibited. */
      cij[TMI] += wgt[idx];
      cij[TII] += (insopt[idx]-1) * wgt[idx];
      cij[TIM] += wgt[idx];
    } else {
      if (!isgap(aseqs[idx][i])) {
	if (aseqs[idx][j] == '_')      ; /* YO! what to do with trailer? */
	else if (isgap(aseqs[idx][j])) cij[TMD] += wgt[idx];
	else                           cij[TMM] += wgt[idx];
      } else {			/* ignores B->E possibility */
	if (aseqs[idx][j] == '_')      continue;
	else if (isgap(aseqs[idx][j])) cij[TDD] += wgt[idx];
	else                           cij[TDM] += wgt[idx];
      }
    }
  }
  return 1;
}


/* Function: estimate_model_length()
 * 
 * Purpose:  Return a decent guess about the length of the model,
 *           based on the lengths of the sequences.
 *           
 *           Algorithm is dumb: use weighted average length.
 *           
 *           Don't assume that weights sum to nseq!
 */
static int
estimate_model_length(MSA *msa)
{
  int   idx;
  float total = 0.;
  float wgtsum = 0.;

  for (idx = 0; idx < msa->nseq; idx++)
    {
      total  += msa->wgt[idx] * DealignedLength(msa->aseq[idx]);
      wgtsum += msa->wgt[idx];
    }
    
  return (int) (total / wgtsum);
}
 

/* Function: matassign2hmm()
 * 
 * Purpose:  Given an assignment of alignment columns to match vs.
 *           insert, finish the final part of the model construction 
 *           calculation that is constant between model construction
 *           algorithms.
 *           
 * Args:     msa       - multiple sequence alignment
 *           dsq       - digitized unaligned aseq's
 *           matassign - 1..alen bit flags for column assignments
 *           ret_hmm   - RETURN: counts-form HMM
 *           ret_tr    - RETURN: array of tracebacks for aseq's
 *                         
 * Return:   (void)
 *           ret_hmm and ret_tr alloc'ed here for the calling
 *           modelmaker function.
 */
static void
matassign2hmm(MSA *msa, char **dsq, int *matassign, 
	      struct plan7_s **ret_hmm, struct p7trace_s ***ret_tr)
{
  struct plan7_s    *hmm;       /* RETURN: new hmm                     */
  struct p7trace_s **tr;        /* fake tracebacks for each seq        */
  int      M;                   /* length of new model in match states */
  int      idx;                 /* counter over sequences              */
  int      apos;                /* counter for aligned columns         */

				/* how many match states in the HMM? */
  M = 0;
  for (apos = 1; apos <= msa->alen; apos++) {
    if (matassign[apos] & ASSIGN_MATCH) 
      M++;
  }

  if (M == 0) 
    Die("No conserved consensus columns found; aborting construction!\n\
This is an unusual situation. Reexamine your sequence alignment. It is\n\
probably unusually full of gaps, or lots of sequence fragments. You may be\n\
able to force HMMER to model it; see the --fast (and --gapmax), or --hand\n\
options to hmmbuild.");

				/* delimit N-terminal tail */
  for (apos=1; matassign[apos] & ASSIGN_INSERT && apos <= msa->alen; apos++)
    matassign[apos] |= EXTERNAL_INSERT_N;
  if (apos <= msa->alen) matassign[apos] |= FIRST_MATCH;

				/* delimit C-terminal tail */
  for (apos=msa->alen; matassign[apos] & ASSIGN_INSERT && apos > 0; apos--)
    matassign[apos] |= EXTERNAL_INSERT_C;
  if (apos > 0) matassign[apos] |= LAST_MATCH;

  /* print_matassign(matassign, msa->alen);  */

				/* make fake tracebacks for each seq */
  fake_tracebacks(msa->aseq, msa->nseq, msa->alen, matassign, &tr);
				/* build model from tracebacks */
  hmm = AllocPlan7(M);
  ZeroPlan7(hmm);
  for (idx = 0; idx < msa->nseq; idx++) {
    /* P7PrintTrace(stdout, tr[idx], NULL, NULL);   */
    P7TraceCount(hmm, dsq[idx], msa->wgt[idx], tr[idx]);
  }
				/* annotate new model */
  annotate_model(hmm, matassign, msa);

  /* Set #=RF line of alignment to reflect our assignment
   * of match, delete. matassign is valid from 1..alen and is off
   * by one from msa->rf.
   */
  if (msa->rf != NULL) free(msa->rf);
  msa->rf = (char *) MallocOrDie (sizeof(char) * (msa->alen + 1));
  for (apos = 0; apos < msa->alen; apos++)
    msa->rf[apos] = matassign[apos+1] & ASSIGN_MATCH ? 'x' : '.';
  msa->rf[msa->alen] = '\0';

				/* Cleanup and return. */
  if (ret_tr != NULL) *ret_tr = tr;
  else   { for (idx = 0; idx < msa->nseq; idx++) P7FreeTrace(tr[idx]); free(tr); }
  if (ret_hmm != NULL) *ret_hmm = hmm; else FreePlan7(hmm);
  return;
}
  


/* Function: fake_tracebacks()
 * 
 * Purpose:  From a consensus assignment of columns to MAT/INS, construct fake
 *           tracebacks for each individual sequence.
 *           
 * Note:     Fragment tolerant by default. Internal entries are 
 *           B->M_x, instead of B->D1->D2->...->M_x; analogously
 *           for internal exits. 
 *           
 * Args:     aseqs     - alignment [0..nseq-1][0..alen-1]
 *           nseq      - number of seqs in alignment
 *           alen      - length of alignment in columns
 *           matassign - assignment of column; [1..alen] (off one from aseqs)
 *           ret_tr    - RETURN: array of tracebacks
 *           
 * Return:   (void)
 *           ret_tr is alloc'ed here. Caller must free.
 */          
static void
fake_tracebacks(char **aseq, int nseq, int alen, int *matassign,
		struct p7trace_s ***ret_tr)
{
  struct p7trace_s **tr;
  int  idx;                     /* counter over sequences          */
  int  i;                       /* position in raw sequence (1..L) */
  int  k;                       /* position in HMM                 */
  int  apos;                    /* position in alignment columns   */
  int  tpos;			/* position in traceback           */

  tr = (struct p7trace_s **) MallocOrDie (sizeof(struct p7trace_s *) * nseq);
  
  for (idx = 0; idx < nseq; idx++)
    {
      P7AllocTrace(alen+6, &tr[idx]); /* allow room for S,N,B,E,C,T */
      
				/* all traces start with S state... */
      tr[idx]->statetype[0] = STS;
      tr[idx]->nodeidx[0]   = 0;
      tr[idx]->pos[0]       = 0;
				/* ...and transit to N state; N-term tail
				   is emitted on N->N transitions */
      tr[idx]->statetype[1] = STN;
      tr[idx]->nodeidx[1]   = 0;
      tr[idx]->pos[1]       = 0;
      
      i = 1;
      k = 0;
      tpos = 2;
      for (apos = 0; apos < alen; apos++)
        {
	  tr[idx]->statetype[tpos] = STBOGUS; /* bogus, deliberately, to debug */

	  if (matassign[apos+1] & FIRST_MATCH)
	    {			/* BEGIN */
	      tr[idx]->statetype[tpos] = STB;
	      tr[idx]->nodeidx[tpos]   = 0;
	      tr[idx]->pos[tpos]       = 0;
	      tpos++;
	    }

	  if (matassign[apos+1] & ASSIGN_MATCH && ! isgap(aseq[idx][apos]))
	    {			/* MATCH */
	      k++;		/* move to next model pos */
	      tr[idx]->statetype[tpos] = STM;
	      tr[idx]->nodeidx[tpos]   = k;
	      tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
	    }	      
          else if (matassign[apos+1] & ASSIGN_MATCH)
            {                   /* DELETE */
		/* being careful about S/W transitions; no B->D transitions */
	      k++;		/* *always* move on model when ASSIGN_MATCH */
	      if (tr[idx]->statetype[tpos-1] != STB)
		{
		  tr[idx]->statetype[tpos] = STD;
		  tr[idx]->nodeidx[tpos]   = k;
		  tr[idx]->pos[tpos]       = 0;
		  tpos++;
		}
            }
          else if (matassign[apos+1] & EXTERNAL_INSERT_N &&
		   ! isgap(aseq[idx][apos]))
            {                   /* N-TERMINAL TAIL */
              tr[idx]->statetype[tpos] = STN;
              tr[idx]->nodeidx[tpos]   = 0;
              tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
            }
	  else if (matassign[apos+1] & EXTERNAL_INSERT_C &&
		   ! isgap(aseq[idx][apos]))
	    {			/* C-TERMINAL TAIL */
	      tr[idx]->statetype[tpos] = STC;
              tr[idx]->nodeidx[tpos]   = 0;
              tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
	    }
	  else if (! isgap(aseq[idx][apos]))
	    {			/* INSERT */
	      tr[idx]->statetype[tpos] = STI;
              tr[idx]->nodeidx[tpos]   = k;
              tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
	    }

	  if (matassign[apos+1] & LAST_MATCH)
	    {			/* END */
	      /* be careful about S/W transitions; may need to roll
	       * back over some D's because there's no D->E transition
	       */
	      while (tr[idx]->statetype[tpos-1] == STD) 
		tpos--;
	      tr[idx]->statetype[tpos] = STE;
	      tr[idx]->nodeidx[tpos]   = 0;
	      tr[idx]->pos[tpos]       = 0;
	      tpos++;
				/* and then transit E->C;
				   alignments that use J are undefined;
				   C-term tail is emitted on C->C transitions */
	      tr[idx]->statetype[tpos] = STC;
	      tr[idx]->nodeidx[tpos]   = 0;
	      tr[idx]->pos[tpos]       = 0;
	      tpos++;
	    }
        }
                                /* all traces end with T state */
      tr[idx]->statetype[tpos] = STT;
      tr[idx]->nodeidx[tpos]   = 0;
      tr[idx]->pos[tpos]       = 0;
      tr[idx]->tlen = ++tpos;
				/* deal with DI, ID transitions */
				/* k == M here */
      trace_doctor(tr[idx], k, NULL, NULL);

    }    /* end for sequence # idx */

  *ret_tr = tr;
  return;
}

/* Function: trace_doctor()
 * 
 * Purpose:  Plan 7 disallows D->I and I->D "chatter" transitions.
 *           However, these transitions may be implied by many
 *           alignments for hand- or heuristic- built HMMs.
 *           trace_doctor() collapses I->D or D->I into a
 *           single M position in the trace. 
 *           Similarly, B->I and I->E transitions may be implied
 *           by an alignment.
 *           
 *           trace_doctor does not examine any scores when it does
 *           this. In ambiguous situations (D->I->D) the symbol
 *           will be pulled arbitrarily to the left, regardless
 *           of whether that's the best column to put it in or not.
 *           
 * Args:     tr      - trace to doctor
 *           M       - length of model that traces are for 
 *           ret_ndi - number of DI transitions doctored
 *           ret_nid - number of ID transitions doctored
 * 
 * Return:   (void)
 *           tr is modified
 */               
static void
trace_doctor(struct p7trace_s *tr, int mlen, int *ret_ndi, int *ret_nid)
{
  int opos;			/* position in old trace                 */
  int npos;			/* position in new trace (<= opos)       */
  int ndi, nid;			/* number of DI, ID transitions doctored */

				/* overwrite the trace from left to right */
  ndi  = nid  = 0;
  opos = npos = 0;
  while (opos < tr->tlen) {
      /* fix implied D->I transitions; D transforms to M, I pulled in */
    if (tr->statetype[opos] == STD && tr->statetype[opos+1] == STI) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos]; /* D transforms to M      */
      tr->pos[npos]       = tr->pos[opos+1];   /* insert char moves back */
      opos += 2;
      npos += 1;
      ndi++;
    } /* fix implied I->D transitions; D transforms to M, I is pushed in */
    else if (tr->statetype[opos]== STI && tr->statetype[opos+1]== STD) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos+1];/* D transforms to M    */
      tr->pos[npos]       = tr->pos[opos];      /* insert char moves up */
      opos += 2;
      npos += 1;
      nid++; 
    } /* fix implied B->I transitions; pull I back to its M */
    else if (tr->statetype[opos]== STI && tr->statetype[opos-1]== STB) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos]; /* offending I transforms to M */
      tr->pos[npos]       = tr->pos[opos];
      opos++;
      npos++;
    } /* fix implied I->E transitions; push I to next M */
    else if (tr->statetype[opos]== STI && tr->statetype[opos+1]== STE) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos]+1;/* offending I transforms to M */
      tr->pos[npos]       = tr->pos[opos];
      opos++;
      npos++;
    } /* rare: N-N-B-E becomes N-B-M_1-E (swap B,N) */
    else if (tr->statetype[opos]==STB && tr->statetype[opos+1]==STE
	     && tr->statetype[opos-1]==STN && tr->pos[opos-1] > 0) {   
      tr->statetype[npos]   = STM;
      tr->nodeidx[npos]     = 1;
      tr->pos[npos]         = tr->pos[opos-1];
      tr->statetype[npos-1] = STB;
      tr->nodeidx[npos-1]   = 0;
      tr->pos[npos-1]       = 0;
      opos++;
      npos++;
    } /* rare: B-E-C-C-x becomes B-M_M-E-C-x (swap E,C) */
    else if (tr->statetype[opos]==STE && tr->statetype[opos-1]==STB
	     && tr->statetype[opos+1]==STC 
	     && tr->statetype[opos+2]==STC) { 
      tr->statetype[npos]   = STM;
      tr->nodeidx[npos]     = mlen;
      tr->pos[npos]         = tr->pos[opos+2];
      tr->statetype[npos+1] = STE;
      tr->nodeidx[npos+1]   = 0;
      tr->pos[npos+1]       = 0;
      tr->statetype[npos+2] = STC; /* first C must be a nonemitter  */
      tr->nodeidx[npos+2]   = 0;
      tr->pos[npos+2]       = 0;
      opos+=3;
      npos+=3;
    } /* everything else is just copied */
    else {
      tr->statetype[npos] = tr->statetype[opos];
      tr->nodeidx[npos]   = tr->nodeidx[opos];
      tr->pos[npos]       = tr->pos[opos];
      opos++;
      npos++;
    }
  }
  tr->tlen = npos;

  if (ret_ndi != NULL) *ret_ndi = ndi;
  if (ret_nid != NULL) *ret_nid = nid;
  return;
}


/* Function: annotate_model()
 * 
 * Purpose:  Add rf, cs optional annotation to a new model.
 * 
 * Args:     hmm       - new model
 *           matassign - which alignment columns are MAT; [1..alen]
 *           msa       - alignment, including annotation to transfer
 *           
 * Return:   (void)
 */
static void
annotate_model(struct plan7_s *hmm, int *matassign, MSA *msa)
{                      
  int   apos;			/* position in matassign, 1.alen  */
  int   k;			/* position in model, 1.M         */
  char *pri;			/* X-PRM, X-PRI, X-PRT annotation */

  /* Transfer reference coord annotation from alignment,
   * if available
   */
  if (msa->rf != NULL) {
    hmm->rf[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH) /* ainfo is off by one from HMM */
	hmm->rf[k++] = (msa->rf[apos-1] == ' ') ? '.' : msa->rf[apos-1];
    hmm->rf[k] = '\0';
    hmm->flags |= PLAN7_RF;
  }

  /* Transfer consensus structure annotation from alignment, 
   * if available
   */
  if (msa->ss_cons != NULL) {
    hmm->cs[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH)
	hmm->cs[k++] = (msa->ss_cons[apos-1] == ' ') ? '.' : msa->ss_cons[apos-1];
    hmm->cs[k] = '\0';
    hmm->flags |= PLAN7_CS;
  }

  /* Transfer surface accessibility annotation from alignment,
   * if available
   */
  if (msa->sa_cons != NULL) {
    hmm->ca[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH)
	hmm->ca[k++] = (msa->sa_cons[apos-1] == ' ') ? '.' : msa->sa_cons[apos-1];
    hmm->ca[k] = '\0';
    hmm->flags |= PLAN7_CA;
  }

  /* Store the alignment map
   */
  for (apos = k = 1; apos <= msa->alen; apos++)
    if (matassign[apos] & ASSIGN_MATCH)
      hmm->map[k++] = apos;
  hmm->flags |= PLAN7_MAP;

  /* Translate and transfer X-PRM annotation. 
   * 0-9,[a-zA-Z] are legal; translate as 0-9,10-35 into hmm->mpri.
   * Any other char is translated as -1, and this will be interpreted
   * as a flag that means "unknown", e.g. use the normal mixture Dirichlet
   * procedure for this column.
   */
  if ((pri = MSAGetGC(msa, "X-PRM")) != NULL)
    {
      hmm->mpri = MallocOrDie(sizeof(int) * (hmm->M+1));
      for (apos = k = 1; apos <= msa->alen; apos++)
	if (matassign[apos] & ASSIGN_MATCH)
	  {
	    if      (isdigit((int) pri[apos-1])) hmm->mpri[k] = pri[apos-1] - '0';
	    else if (islower((int) pri[apos-1])) hmm->mpri[k] = pri[apos-1] - 'a' + 10;
	    else if (isupper((int) pri[apos-1])) hmm->mpri[k] = pri[apos-1] - 'A' + 10;
	    else hmm->mpri[k] = -1;
	    k++;
	  }
    }
  /* And again for X-PRI annotation on insert priors:
   */
  if ((pri = MSAGetGC(msa, "X-PRI")) != NULL)
    {
      hmm->ipri = MallocOrDie(sizeof(int) * (hmm->M+1));
      for (apos = k = 1; apos <= msa->alen; apos++)
	if (matassign[apos] & ASSIGN_MATCH)
	  {
	    if      (isdigit((int) pri[apos-1])) hmm->ipri[k] = pri[apos-1] - '0';
	    else if (islower((int) pri[apos-1])) hmm->ipri[k] = pri[apos-1] - 'a' + 10;
	    else if (isupper((int) pri[apos-1])) hmm->ipri[k] = pri[apos-1] - 'A' + 10;
	    else hmm->ipri[k] = -1;
	    k++;
	  }
    }
  /* And one last time for X-PRT annotation on transition priors:
   */
  if ((pri = MSAGetGC(msa, "X-PRT")) != NULL)
    {
      hmm->tpri = MallocOrDie(sizeof(int) * (hmm->M+1));
      for (apos = k = 1; apos <= msa->alen; apos++)
	if (matassign[apos] & ASSIGN_MATCH)
	  {
	    if      (isdigit((int) pri[apos-1])) hmm->tpri[k] = pri[apos-1] - '0';
	    else if (islower((int) pri[apos-1])) hmm->tpri[k] = pri[apos-1] - 'a' + 10;
	    else if (isupper((int) pri[apos-1])) hmm->tpri[k] = pri[apos-1] - 'A' + 10;
	    else hmm->tpri[k] = -1;
	    k++;
	  }
    }

}

static void
print_matassign(int *matassign, int alen)
{
  int apos;

  for (apos = 0; apos <= alen; apos++) {
    printf("%3d %c %c %c\n", 
	   apos,
	   (matassign[apos] & ASSIGN_MATCH) ? 'x':' ',
	   (matassign[apos] & FIRST_MATCH || matassign[apos] & LAST_MATCH) ? '<' : ' ',
	   (matassign[apos] & EXTERNAL_INSERT_N ||
	    matassign[apos] & EXTERNAL_INSERT_C) ? '|':' ');
  }
}
