/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* masks.c
 * SRE, Tue Nov 18 10:12:28 1997
 * 
 * Sequence masking routines: XNU, SEG.
 * 
 * RCS $Id$
 */

#include <stdio.h>
#include <math.h>

#include "squid.h"
#include "config.h"
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
