/* masks.c
 * Sequence masking routines; corrections for biased composition
 * target sequences. 
 * 
 * The Claverie/States XNU code is not used by default because I 
 * consider X'ing out sequence to be too black/white and too 
 * aggressive, but it's available as an option.
 * 
 * The Wooton/Federhen SEG code was studied, but deemed too 
 * nonportable to include; it would've suffered the same drawback
 * as XNU.
 * 
 * The TraceScoreCorrection() code is the default.
 * 
 * SRE, Tue Nov 18 10:12:28 1997
 * SVN $Id: masks.c 1382 2005-05-03 22:25:36Z eddy $
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "squid.h"

#include "plan7.h"
#include "structs.h"
#include "funcs.h"

/* The PAM120 score matrix, in HMMER's AMINO_ALPHABET alphabetic order 
 */
static int xpam120[23][23] = {
  { 3, -3,  0,  0, -4,  1, -3, -1, -2, -3, -2, -1,  1, -1, -3,  1,  1,  0, -7, -4,  1,  0,  0 },
  {-3,  9, -7, -7, -6, -4, -4, -3, -7, -7, -6, -5, -4, -7, -4,  0, -3, -3, -8, -1, -4, -6,  0 },
  { 0, -7,  5,  3, -7,  0,  0, -3, -1, -5, -4,  2, -3,  1, -3,  0, -1, -3, -8, -5,  5,  3,  0 },
  { 0, -7,  3,  5, -7, -1, -1, -3, -1, -4, -3,  1, -2,  2, -3, -1, -2, -3, -8, -5,  3,  5,  0 },
  {-4, -6, -7, -7,  8, -5, -3,  0, -7,  0, -1, -4, -5, -6, -5, -3, -4, -3, -1,  4, -4, -5,  0 },
  { 1, -4,  0, -1, -5,  5, -4, -4, -3, -5, -4,  0, -2, -3, -4,  1, -1, -2, -8, -6,  1, -1,  0 },
  {-3, -4,  0, -1, -3, -4,  7, -4, -2, -3, -4,  2, -1,  3,  1, -2, -3, -3, -3, -1,  2,  2,  0 },
  {-1, -3, -3, -3,  0, -4, -4,  6, -3,  1,  1, -2, -3, -3, -2, -2,  0,  3, -6, -2, -2, -2,  0 },
  {-2, -7, -1, -1, -7, -3, -2, -3,  5, -4,  0,  1, -2,  0,  2, -1, -1, -4, -5, -5,  1,  0,  0 },
  {-3, -7, -5, -4,  0, -5, -3,  1, -4,  5,  3, -4, -3, -2, -4, -4, -3,  1, -3, -2, -3, -2,  0 },
  {-2, -6, -4, -3, -1, -4, -4,  1,  0,  3,  8, -3, -3, -1, -1, -2, -1,  1, -6, -4, -3, -1,  0 },
  {-1, -5,  2,  1, -4,  0,  2, -2,  1, -4, -3,  4, -2,  0, -1,  1,  0, -3, -4, -2,  4,  1,  0 },
  { 1, -4, -3, -2, -5, -2, -1, -3, -2, -3, -3, -2,  6,  0, -1,  1, -1, -2, -7, -6, -1,  0,  0 },
  {-1, -7,  1,  2, -6, -3,  3, -3,  0, -2, -1,  0,  0,  6,  1, -2, -2, -3, -6, -5,  1,  5,  0 },
  {-3, -4, -3, -3, -5, -4,  1, -2,  2, -4, -1, -1, -1,  1,  6, -1, -2, -3,  1, -5, -1,  0,  0 },
  { 1,  0,  0, -1, -3,  1, -2, -2, -1, -4, -2,  1,  1, -2, -1,  3,  2, -2, -2, -3,  1,  0,  0 },
  { 1, -3, -1, -2, -4, -1, -3,  0, -1, -3, -1,  0, -1, -2, -2,  2,  4,  0, -6, -3,  1, -1,  0 },
  { 0, -3, -3, -3, -3, -2, -3,  3, -4,  1,  1, -3, -2, -3, -3, -2,  0,  5, -8, -3, -2, -2,  0 },
  {-7, -8, -8, -8, -1, -8, -3, -6, -5, -3, -6, -4, -7, -6,  1, -2, -6, -8, 12, -2, -5, -6,  0 },
  {-4, -1, -5, -5,  4, -6, -1, -2, -5, -2, -4, -2, -6, -5, -5, -3, -3, -3, -2,  8, -2, -4,  0 },
  { 1, -4,  5,  3, -4,  1,  2, -2,  1, -3, -3,  4, -1,  1, -1,  1,  1, -2, -5, -2,  6,  4,  0 },
  { 0, -6,  3,  5, -5, -1,  2, -2,  0, -2, -1,  1,  0,  5,  0,  0, -1, -2, -6, -4,  4,  6,  0 },
  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
};


/* Function: XNU()
 * Date:     18 Nov 1997 [StL]
 * 
 * Purpose:  x-out of repetitive sequence. XNU tends to be
 *           good at x'ing out short period tandem repeats.
 *           
 * Note:     Apply /only/ to protein sequence.            
 * 
 * Args:     dsq: 1..len digitized sequence
 *           len: length of dsq
 *           
 * Return:   number of characters x'ed out.
 */            
int
XNU(unsigned char *dsq, int len)
{
  int    i,k,off,sum,beg,end,top;
  int    topcut,fallcut;
  double s0;
  int    noff = 4;		/* maximum search offset */
  int    mcut = 1;
  double pcut = 0.01;
  int   *hit;
  double lambda = 0.346574;
  double K      = 0.2;
  double H      = 0.664;
  int    xnum   = 0;

  if (len == 0) return 0;

  hit = MallocOrDie(sizeof(int) * (len+1));
  for (i=1; i<=len; i++) hit[i]=0;

  /*
  ** Determine the score cutoff so that pcut will be the fraction
  ** of random sequence eliminated assuming lambda, K, and H are
  ** characteristic of the database as a whole
  */
  s0 = - log( pcut*H / (noff*K) ) / lambda;
  if (s0>0) topcut = floor(s0 + log(s0)/lambda + 0.5);
  else topcut = 0;
  fallcut = (int)log(K/0.001)/lambda;

  for (off=mcut; off<=noff; off++) {
    sum=top=0;
    beg=off;
    end=0;

    for (i=off+1; i<=len; i++) {
      sum += xpam120[dsq[i]][dsq[i-off]];
      if (sum>top) {
	top=sum;
	end=i;
      }
      if (top>=topcut && top-sum>fallcut) {
	for (k=beg; k<=end; k++) 
	  hit[k] = hit[k-off] = 1;
	sum=top=0;
	beg=end=i+1;
      } else if (top-sum>fallcut) {
	sum=top=0;
	beg=end=i+1;
      }
      if (sum<0) {
	beg=end=i+1;
	sum=top=0;
      }
    }
    if (top>=topcut) {
      for (k=beg; k<=end; k++) 
	hit[k] = hit[k-off] = 1;
    }
  }
  
  /* Now mask off detected repeats
   */
  for (i=1; i<=len; i++) 
    if (hit[i]) { xnum++; dsq[i] = Alphabet_iupac-1;} /* e.g. 'X' */

  free(hit);
  return xnum;
}


/* Function: TraceScoreCorrection()
 * Date:     Sun Dec 21 12:05:47 1997 [StL]
 * 
 * Purpose:  Calculate a correction (in integer log_2 odds) to be
 *           applied to a sequence, using a second null model, 
 *           based on a traceback. M/I emissions are corrected;
 *           C/N/J are not -- as if the nonmatching part and 
 *           matching part were each generated by the best null model.
 *           The null model is constructed /post hoc/ as the
 *           average over all the M,I distributions used by the trace.
 *           
 * Return:   the log_2-odds score correction.          
 */
float
TraceScoreCorrection(struct plan7_s *hmm, struct p7trace_s *tr, unsigned char *dsq)
{
  float p[MAXABET];		/* null model distribution */
  int   sc[MAXCODE];		/* null model scores       */
  int   x;
  int   tpos;
  int   score;

  /* Rarely, the alignment was totally impossible, and tr is NULL.
   */
  if (tr == NULL) return 0.0;

  /* Set up model: average over the emission distributions of
   * all M, I states that appear in the trace. Ad hoc? Sure, you betcha. 
   */
  FSet(p, Alphabet_size, 0.0);
  for (tpos = 0; tpos < tr->tlen; tpos++)
     if      (tr->statetype[tpos] == STM) 
       FAdd(p, hmm->mat[tr->nodeidx[tpos]], Alphabet_size);
     else if (tr->statetype[tpos] == STI) 
       FAdd(p, hmm->ins[tr->nodeidx[tpos]], Alphabet_size);
  FNorm(p, Alphabet_size);

  for (x = 0; x < Alphabet_size; x++)
    sc[x] = Prob2Score(p[x], hmm->null[x]);
				/* could avoid this chunk if we knew
				   we didn't need any degenerate char scores */
  for (x = Alphabet_size; x < Alphabet_iupac; x++)
    sc[x] = DegenerateSymbolScore(p, hmm->null, x);
					       

  /* Score all the M,I state emissions that appear in the trace.
   */
   score = 0;
   for (tpos = 0; tpos < tr->tlen; tpos++)
     if (tr->statetype[tpos] == STM || tr->statetype[tpos] == STI)
       score += sc[dsq[tr->pos[tpos]]];

   /* Apply an ad hoc 8 bit fudge factor penalty;
    * interpreted as a prior, saying that the second null model is 
    * 1/2^8 (1/256) as likely as the standard null model
    */
   score -= 8 * INTSCALE;	

   /* Return the correction to the bit score.
    */
   return Scorify(ILogsum(0, score));	
}


/* THE FOLLOWING CODE IS IN DEVELOPMENT.
 * it is commented out of the current release deliberately.
 * If you activate it, I'm not responsible for the consequences.
 */
#if MICHAEL_JORDAN_BUYS_THE_PACERS
/* Function: NewTraceScoreCorrection()
 * Date:     Wed Feb 17 14:32:45 1999 [StL]
 * 
 * Purpose:  Calculate a correction (in integer log_2 odds) to be
 *           applied to a sequence, using a second null model, 
 *           based on sequence endpoints. M/I emissions are corrected;
 *           C/N/J are not -- as if the nonmatching part and 
 *           matching part were each generated by the best null model.
 *           Each null model is constructed /post hoc/ from the
 *           sequence composition of each matching domain (e.g.
 *           a null2 model is constructed for each domain in a 
 *           multihit trace).
 *           
 *           Constraints on the construction of this function include:
 *            1) Paracel hardware can't deal with trace-dependent
 *               null2 models. Original implementation of 
 *               TraceScoreCorrection() was dependent on traceback
 *               and could not be reproduced on GeneMatcher.
 *               GeneMatcher may be able to deal w/ sequence endpoint
 *               dependent rescoring, though.
 *               Although this function looks like it's trace-
 *               dependent (because it's being passed a p7trace_s
 *               structure), it's really not; only the sequence
 *               endpoints are being used.
 *               
 *            2) It is desirable that for multihit traces,
 *               per-domain scores sum to the per-sequence score.
 *               Otherwise people see this as a "bug" (cf. 
 *               bug #2, David Kerk, NRC). HMMER calculates the
 *               per-domain scores by going through a separate
 *               TraceScore() call for each one and separately
 *               correcting them with TraceScoreCorrection(),
 *               so we have to do each domain in a full trace
 *               by a similar mechanism -- even if this means that
 *               we're adopting a very dubiously post hoc
 *               null model.
 *           
 * Return:   the log_2-odds score correction.          
 */
float
NewTraceScoreCorrection(struct plan7_s *hmm, struct p7trace_s *tr, unsigned char *dsq)
{
  float ct[MAXABET];		/* counts of observed residues            */
  float p[MAXABET];		/* null2 model distribution (also counts) */
  float sc[MAXCODE];		/* null2 model scores (as floats not int) */

  int   x;
  int   tpos;
  int   score;			/* tmp score for real HMM, integer logodds */
  float hmmscore;		/* score for real HMM for this domain */
  float null2score;		/* score for null2 model for this domain */


  float totscore;		/* overall score for trace */
  float maxscore;		/* best score so far for single domain */
  int   in_domain;		/* flag for whether we're counting this domain */
  unsigned char sym;		/* digitized symbol in dsq */
  int   ndom;			/* number of domains counted towards score */

  int   nsym;			/* number of symbols in this alignment */

  totscore  = 0.;
  maxscore  = -FLT_MAX;
  in_domain = FALSE;
  ndom      = 0;
  for (tpos = 0; tpos < tr->tlen; tpos++)
    {
				/* detect start of domain; start at N or J */
      if (tpos < tr->tlen-1 && tr->statetype[tpos+1] == STB)
	{
	  FCopy(ct, hmm->null, Alphabet_size);  /* simple Dirichlet prior */
	  score      = 0;
	  null2score = 0.;
	  nsym       = 0;
	  in_domain  = TRUE;
	}
		/* Count stuff in domain starting with N->B or J->B transition */
      if (in_domain) {
	sym = dsq[tr->pos[tpos]];

				/* count emitted symbols in domain */
	if (tr->statetype[tpos] == STM || tr->statetype[tpos] == STI)
	  {
	    P7CountSymbol(ct, sym, 1.0);
	    nsym++;
	  }

				/* score emitted symbols in domain towards HMM */
	if (tr->statetype[tpos] == STM) 
	  score += hmm->msc[sym][tr->nodeidx[tpos]];
	else if (tr->statetype[tpos] == STI) 
	  score += hmm->isc[sym][tr->nodeidx[tpos]];
				/* score transitions in domain towards HMM */
	score += TransitionScoreLookup(hmm, 
				       tr->statetype[tpos], tr->nodeidx[tpos],
				       tr->statetype[tpos+1], tr->nodeidx[tpos+1]);
      }

      
      if (tr->statetype[tpos] == STE) /* done w/ a domain; calc its score */
	{
				      /* convert counts to null2 prob distribution */
	  FCopy(p, ct, Alphabet_size);
	  FNorm(p, Alphabet_size); 
			              /* Convert probs to log-odds_e scores */
				      /* p can't be zero, because of prior  */
	  for (x = 0; x < Alphabet_size; x++)
	    sc[x] = log(p[x] / hmm->null[x]);
				      /* null2 score = counts \dot scores */
	  null2score = FDot(ct, sc, Alphabet_size);
	  
	  printf("NSYM = %d   NULL2 = %.1f\n", nsym, null2score);

	  /* Apply an ad hoc 12 bit fudge factor penalty, per domain.
	   * Interpreted probabilistically, saying that there's about
           * a 1/256 probability to transition into the second null model.
	   */
	  null2score  -= 12.;
	  
	  /* Now correct score1 using the null2 score.
	   * If it's still > 0, add it to accumulated score.
	   */
	  hmmscore  = Scorify(score);
	  hmmscore -= 1.44269504 * LogSum(0, null2score);
	  if (hmmscore > 0.) { totscore += hmmscore; ndom++; }
	  if (hmmscore > maxscore) maxscore = hmmscore;

	  in_domain = FALSE;
	}
    }

  /* Single domain special case.
   */
  if (ndom == 0) totscore = maxscore;

   /* Return the correction to the bit score
    */
   return (P7TraceScore(hmm, dsq, tr) - totscore);
}
#endif /*0*/


float
SantaCruzCorrection(struct plan7_s *hmm, struct p7trace_s *tr, unsigned char *dsq)
{
  return 0.0;			/* UNFINISHED CODE */
}

/************************************************************
 * @LICENSE@
 ************************************************************/

