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


void 
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

  if (len == 0) return;

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
    if (hit[i]) dsq[i] = Alphabet_iupac-1; /* e.g. 'X' */

  free(hit);
  return;
}



/* Function: MiniSEG()
 * Date:     SRE, Tue Dec 16 14:04:58 1997 [StL]
 * 
 * Purpose:  A two-state HMM with a random amino acid distribution in
 *           state R and and a biased one in state X. Viterbi parse 
 *           the sequence, then replace X out all residues
 *           assigned to state X.
 *           
 * Args:     dsq  - 1..len digitized sequence
 *           len  - length of sequence
 *           
 * Returns:  number of x's introduced.          
 */
int
MiniSEG(char *dsq, int len, float *ret_max)
{
  float xp[20] = { 0.0681, 0.0120, 0.0623, 0.0651, 0.0313, 0.0902, 0.0241, 0.0371, 0.0687, 0.0676,
		   0.0143, 0.0548, 0.0647, 0.0415, 0.0551, 0.0926, 0.0623, 0.0505, 0.0102, 0.0269
                 };
  float p1 = 0.9999;		/* R->R transition prob. */
  float p2 = 0.95;		/* X->X transition prob. */
  int re[23];
  int xe[23];
  int trx,trr,txr,txx;
  int x,i;
  int osx, nsx;			/* old S_x, new S_x: two rows of DP matrix */
  int osr, nsr;			/* old S_r, new S_r  */
  char *tbr;                    /* traceback pointers for R: 0 = from R, 1 = from X */
  char *tbx;                    /* traceback ptrs for X: 0 = from R, 1 = from X */
  int  state;
  int  score;			/* running score   */
  int  max;			/* maximum score   */
  int  xnum;			/* number x'ed out */
  
  /* Set up the miniseg two state HMM, as integer log probabilities.
   * Should be generalized: parameters elsewhere (and shared with prior.c),
   * and capable of doing other alphabets besides protein. miniseg
   * model could be a structure.
   */
  for (x = 0; x < 20; x++)
    {
      re[x] = Prob2Score(aafq[x], 1.);
      xe[x] = Prob2Score(xp[x], 1.);
    }
  re[20] = Prob2Score(aafq[2] + aafq[11], 1.); /* B = D,N */
  re[21] = Prob2Score(aafq[3] + aafq[13], 1.); /* Z = E,Q */
  re[22] = 0;			               /* X */
  xe[20] = Prob2Score(xp[2] + xp[11], 1.);     /* B = D,N */
  xe[21] = Prob2Score(xp[3] + xp[13], 1.);     /* Z = E,Q */
  xe[22] = 0;	                               /* X */
  trr    = Prob2Score(p1,    1.);
  trx    = Prob2Score(1.-p1, 1.);
  txx    = Prob2Score(p2,    1.);
  txr    = Prob2Score(1.-p2, 1.);
  
  /* Allocate for dp
   */
  tbr = MallocOrDie(sizeof(char) * (len+1));
  tbx = MallocOrDie(sizeof(char) * (len+1));

  /* Initialize at position 0, in R state. 
   * in either state
   */
  osx    = -999999;
  osr    = 0;
  tbr[0] = tbx[0] = -1;
  score  = 0;
  max    = 0;

  /* Dynamic programming, Viterbi fill
   */

  for (i = 1; i <= len; i++) 
    {
      nsx = osr + txr; tbx[i] = 0;
      if (osx + txx > nsx) { nsx = osx + txx; tbx[i] = 1; }
      nsx += xe[(int) dsq[i]];

      nsr = osr + trr; tbr[i] = 0;
      if (osx + txr > nsr) {nsr = osx + txr; tbr[i] = 1; }
      nsr += re[(int) dsq[i]];

      score += xe[(int) dsq[i]] - re[(int) dsq[i]];
      score = MAX(score, 0);
      if (score > max) max = score;

      osr = nsr; osx = nsx;
    }
      
  /* Traceback
   */
  state = (osx + txr > osr) ? 1 : 0;	         /* initialize */
  xnum = 0;
  for (i = len; i >= 1; i--)
    {
      if (state == 1) { xnum++; dsq[i] = Alphabet_iupac-1;} /* X out */
      state = (state == 0) ? tbr[i] : tbx[i];
    }

  free(tbr);
  free(tbx);
  *ret_max = Scorify(max);
  return xnum;
}
