/************************************************************
 * @LICENSE@
 ************************************************************/

/* masks.c
 * SRE, Tue Nov 18 10:12:28 1997
 * 
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
 * RCS $Id$
 */

#include <stdio.h>
#include <math.h>

#include "squid.h"
#include "config.h"
#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

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
XNU(char *dsq, int len)
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
      sum += xpam120[(int) dsq[i]][(int) dsq[i-off]];
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


#if 0
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
TraceScoreCorrection(struct plan7_s *hmm, struct p7trace_s *tr, char *dsq)
{
  float p[MAXABET];		/* null model distribution */
  int   sc[MAXCODE];		/* null model scores       */
  int   x;
  int   tpos;
  int   score;

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
       score += sc[(int) dsq[tr->pos[tpos]]];

   /* Apply an ad hoc 8 bit fudge factor penalty;
    * interpreted as a prior, saying that the second null model is 
    * 1/2^8 (1/256) as likely as the standard null model
    */
   score -= 8 * INTSCALE;	

   /* Return the correction to the bit score.
    */
   return Scorify(ILogsum(0, score));	
}
#endif


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
TraceScoreCorrection(struct plan7_s *hmm, struct p7trace_s *tr, char *dsq)
{
  float ct[MAXABET];		/* counts of observed residues            */
  float p[MAXABET];		/* null2 model distribution (also counts) */
  float sc[MAXCODE];		/* null2 model scores (as floats not int) */

  int   x;
  int   tpos;
  int   score;

  score = 0;
  for (tpos = 0; tpos < tr->tlen; tpos++)
    {

      if (tr->statetype[tpos] == STB)         /* clear the counts       */
	FCopy(ct, hmm->null, Alphabet_size);  /* simple Dirichlet prior */

				      /* count emitted symbols in domain */
      if (tr->statetype[tpos] == STM || tr->statetype[tpos] == STI)
	P7CountSymbol(ct, dsq[tr->pos[tpos]], 1.0);
      
      if (tr->statetype[tpos] == STE) /* done w/ a domain; calc its score */
	{
				      /* convert counts to null2 prob distribution */
	  FCopy(p, ct, Alphabet_size);
	  FNorm(p, Alphabet_size); 
			              /* Convert probs to log-odds_e scores */
				      /* p can't be zero, because of prior  */
	  for (x = 0; x < Alphabet_size; x++)
	    sc[x] = log(p[x] / hmm->null[x]);
				      /* score = counts \dot scores */
	  score += FDot(ct, sc, Alphabet_size);
	  
	  /* Apply an ad hoc 8 bit fudge factor penalty, per domain.
	   * Interpreted probabilistically, saying that there's about
           * a 1/256 probability to transition into the second null model.
	   */
	  score -= 5.5452;	/* our scores are in nats right now, not bits */
	}
    }

   /* Return the correction to the bit score, in bits (hence the 1.44 conversion)
    */
   return (1.44269504 * LogSum(0, score));
}
