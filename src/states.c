/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1997 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* states.c
 * 
 * alloc, free, and initialization of state structures
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "squid.h"
#include "structs.h"
#include "funcs.h"

struct hmm_struc *
AllocHMM(int M)               		/* length of model to make */
{
  struct hmm_struc *hmm;        /* RETURN: blank HMM */

  hmm        = (struct hmm_struc *)   MallocOrDie (sizeof(struct hmm_struc));
  hmm->ins   = (struct basic_state *) MallocOrDie (sizeof(struct basic_state) * (M+2));
  hmm->del   = (struct basic_state *) MallocOrDie (sizeof(struct basic_state) * (M+2));
  hmm->mat   = (struct basic_state *) MallocOrDie (sizeof(struct basic_state) * (M+2));
  hmm->ref   = (char *)  MallocOrDie ((M+2) * sizeof(char));
  hmm->cs    = (char *)  MallocOrDie ((M+2) * sizeof(char));
  hmm->xray  = (float *) MallocOrDie ((M+2) * sizeof(float) * NINPUTS);
  hmm->M     = M;
  hmm->name  = Strdup("unnamed"); /* name is not optional. */

  hmm->flags = 0;
  ZeroHMM(hmm);
  return hmm;
}

/* Function: ZeroHMM()
 * 
 * Purpose:  Zero emission and transition counts in an HMM.
 */
void
ZeroHMM(struct hmm_struc *hmm)
{
  int k, ts, idx;

  for (k = 0; k <= hmm->M+1; k++)
    {
      for (ts = 0; ts < 3; ts++)
	{
	  hmm->mat[k].t[ts] = 0.0;
	  hmm->ins[k].t[ts] = 0.0;
	  hmm->del[k].t[ts] = 0.0;
	}
      for (idx = 0; idx < Alphabet_size; idx++)
	{
	  hmm->mat[k].p[idx]   = 0.0;
	  hmm->ins[k].p[idx]   = 0.0;
	  hmm->del[k].p[idx]   = 0.0;
	}
    }
}


/* Function: LogifyHMM()
 * 
 * Purpose:  Convert a probability-form HMM to log probabilities.
 *           Best to do this on a modifiable copy of an HMM.
 */
void
LogifyHMM(struct hmm_struc *hmm)
{
  int k, ts, idx;

  for (k = 0; k <= hmm->M+1; k++)
    {
      for (ts = 0; ts < 3; ts++)
	{
	  hmm->mat[k].t[ts] = sreLOG2(hmm->mat[k].t[ts]);
	  hmm->ins[k].t[ts] = sreLOG2(hmm->ins[k].t[ts]);
	  hmm->del[k].t[ts] = sreLOG2(hmm->del[k].t[ts]);
	}
      for (idx = 0; idx < Alphabet_size; idx++)
	{
	  hmm->mat[k].p[idx]   = sreLOG2(hmm->mat[k].p[idx]);
	  hmm->ins[k].p[idx]   = sreLOG2(hmm->ins[k].p[idx]);
	}
    }
}

/* Function: LogoddsifyHMM()
 * 
 * Convert a probability form HMM to log odds scores.
 * Best to do this on a modifiable copy of an HMM.
 */
void
LogoddsifyHMM(struct hmm_struc *hmm)
{
  int k, ts, x;

  for (k = 0; k <= hmm->M+1; k++)
    {
      for (ts = 0; ts < 3; ts++)
	{
	  hmm->mat[k].t[ts] = sreLOG2(hmm->mat[k].t[ts]);
	  hmm->ins[k].t[ts] = sreLOG2(hmm->ins[k].t[ts]);
	  hmm->del[k].t[ts] = sreLOG2(hmm->del[k].t[ts]);
	}
      for (x = 0; x < Alphabet_size; x++)
	{
	  hmm->mat[k].p[x]   = sreLOG2(hmm->mat[k].p[x]) - sreLOG2(hmm->null[x]);
	  hmm->ins[k].p[x]   = sreLOG2(hmm->ins[k].p[x]) - sreLOG2(hmm->null[x]);
	}
    }
}


/* Function: WriteFlatPriorHMM()
 * 
 * Purpose:  Fill an HMM with expected probabilities according
 *           to a given prior. Used to construct "flat" initial
 *           models for hmmt.
 */
int
WriteFlatPriorHMM(struct hmm_struc *hmm, struct prior_s *prior)
{
  int    k;			/* counter across model                */
  int    q;			/* counter over mixtures               */
  int    x;			/* counter over symbols or transitions */
  float  malpha;		/* alpha for mixture                   */
  float  ialpha;		/* alpha for insert mixture            */
  float  dalpha;		/* alpha for delete mixture            */

  for (k = 0; k <= hmm->M; k++)
    {
				/* xray info for structure prior */
      if (prior->strategy == PRI_STRUCT)
	{
	  hmm->xray[k*NINPUTS + XRAY_bias] = 1.0;
	  hmm->xray[k*NINPUTS + XRAY_E]    = 0.0;
	  hmm->xray[k*NINPUTS + XRAY_H]    = 0.0;
	  hmm->xray[k*NINPUTS + XRAY_SA]   = 0.0;
	}
				/* match symbol emissions */
      for (x = 0; x < Alphabet_size; x++) 
	hmm->mat[k].p[x] = 0.0;
      if (k > 0)
	for (q = 0; q < prior->mnum; q++)
	  {
	    if (prior->strategy == PRI_STRUCT) 
	      prior->mq[q] = 1.0 / prior->mnum;
	    malpha = 0.0;
	    for (x = 0; x < Alphabet_size; x++)
	      malpha += prior->mat[q][x];
	    for (x = 0; x < Alphabet_size; x++)
	      hmm->mat[k].p[x] += prior->mq[q] * prior->mat[q][x] / malpha;
	  }
				/* insert emissions */
      for (x = 0; x < Alphabet_size; x++) 
	hmm->ins[k].p[x] = 0.0;
      for (q = 0; q < prior->inum; q++)
	{
	  if (prior->strategy == PRI_STRUCT) 
	    prior->iq[q] = 1.0 / prior->inum;
	  ialpha = 0.0;
	  for (x = 0; x < Alphabet_size; x++)
	    ialpha += prior->ins[q][x];	      
	  for (x = 0; x < Alphabet_size; x++)
	    hmm->ins[k].p[x] += prior->iq[q] * prior->ins[q][x] / ialpha;
	}

				/* state transitions  */
      for (x = 0; x < 3; x++)
	hmm->mat[k].t[x] = hmm->ins[k].t[x] = hmm->del[k].t[x] = 0.0;
      for (q = 0; q < prior->tnum; q++)
	{
	  if (prior->strategy == PRI_STRUCT) 
	    prior->tq[q] = 1.0 / prior->tnum;
	  malpha = ialpha = dalpha = 0.0;
	  for (x = 0; x < 3; x++)
	    {
	      malpha += prior->tm[q][x];
	      ialpha += prior->ti[q][x];
	      dalpha += prior->td[q][x];
	    }
	  for (x = 0; x < 3; x++)
	    {
	      hmm->mat[k].t[x] += prior->tq[q] * prior->tm[q][x] / malpha;
	      hmm->ins[k].t[x] += prior->tq[q] * prior->ti[q][x] / ialpha;
	      if (k > 0) hmm->del[k].t[x] += prior->tq[q] * prior->td[q][x] / dalpha;
	    }
	}
    }
				/* the final state never transits to d+1 */
  hmm->mat[hmm->M].t[DELETE] = 0.0;
  hmm->ins[hmm->M].t[DELETE] = 0.0;
  hmm->del[hmm->M].t[DELETE] = 0.0;
  Renormalize(hmm);
  return 1;
}


/* Function: HMMDup()
 * 
 * Purpose:  Create a duplicate copy of an HMM.
 * 
 * Return:   Pointer to the duplicate. 
 *           Caller is responsible for free'ing the duplicate.
 */
struct hmm_struc *
HMMDup(struct hmm_struc *hmm)
{
  struct hmm_struc *newhmm;

  if ((newhmm = AllocHMM(hmm->M)) == NULL)
    Die("AllocHMM() failed");
  HMMCopy(newhmm, hmm);
  return newhmm;
}


/* Function: HMMCopy()
 * 
 * Purpose:  Make a copy of hmm2 in hmm1.
 * 
 * Return:   (void)
 *           Caller promises that hmm1 and hmm2 have identical architectures.
 */
void
HMMCopy(struct hmm_struc *hmm1, struct hmm_struc *hmm2)
{
  int k, x, ts;

  hmm1->flags = hmm2->flags;
  if (hmm1->name != NULL) free(hmm1->name); 
  hmm1->name = Strdup(hmm2->name);

  if (hmm2->flags & HMM_REF)  strcpy(hmm1->ref, hmm2->ref);
  if (hmm2->flags & HMM_CS)   strcpy(hmm1->cs,  hmm2->cs);
  if (hmm2->flags & HMM_XRAY) 
    memcpy(hmm1->xray, hmm2->xray, NINPUTS * (hmm2->M+2) * sizeof(float));
  memcpy(hmm1->null, hmm2->null, sizeof(float) * Alphabet_size);

  for (k = 0; k <= hmm2->M+1; k++)
    {
			/* copy transition T's */
      for (ts = 0; ts < 3; ts++)
	{
	  hmm1->mat[k].t[ts] = hmm2->mat[k].t[ts];
	  hmm1->ins[k].t[ts] = hmm2->ins[k].t[ts];
	  hmm1->del[k].t[ts] = hmm2->del[k].t[ts];
	}
				/* copy symbol P tables */
      for (x = 0; x < Alphabet_size; x++) 
	{
	  hmm1->mat[k].p[x]   = hmm2->mat[k].p[x];
	  hmm1->ins[k].p[x]   = hmm2->ins[k].p[x];
	}
    }
  return;
}


int
FreeHMM(struct hmm_struc *hmm)
{
  if (hmm == NULL) return 0;
  free(hmm->ref);
  free(hmm->cs);
  free(hmm->xray);
  free(hmm->name);
  if (hmm->mat != NULL)  free (hmm->mat);
  if (hmm->ins != NULL)  free (hmm->ins);
  if (hmm->del != NULL)  free (hmm->del);
  free(hmm);
  return 1;
}


struct shmm_s *
AllocSearchHMM(int M)
{
  struct shmm_s *shmm;
  int            x;

  if ((shmm = (struct shmm_s *) malloc (sizeof(struct shmm_s))) == NULL)
    Die("malloc failed");
  for (x = 0; x < 26; x++)
    if ((shmm->m_emit[x] = (int *) calloc (M+1, sizeof(int))) == NULL ||
	(shmm->i_emit[x] = (int *) calloc (M+1, sizeof(int))) == NULL)
      Die("malloc failed");
  if ((shmm->t   = (int *)  malloc (sizeof(int)  * (9*(M+1)))) == NULL ||
      (shmm->ref = (char *) malloc (sizeof(char) * (M+2))) == NULL ||
      (shmm->cs  = (char *) malloc (sizeof(char) * (M+2))) == NULL)
    Die("malloc failed");
  shmm->flags = 0;
  shmm->name = Strdup("nameless");
  shmm->M = M;
  return shmm;
}

void
FreeSearchHMM(struct shmm_s *shmm)
{
  int x;

  for (x = 0; x < 26; x++)
    {
      free(shmm->m_emit[x]);
      free(shmm->i_emit[x]);
    }
  free(shmm->t);
  free(shmm->ref);
  free(shmm->cs);
  free(shmm->name);
  free(shmm);
}


/* Function: CountSymbol()
 * 
 * Purpose:  Given an observed symbol, and a number of counts to
 *           distribute (typically just 1.0), bump the appropriate counter(s).
 * 
 *           This is completely trivial only so long as the symbols
 *           always come from the expected alphabet; since we also
 *           have to deal with degenerate symbols for both nucleic
 *           acid and protein languages, we make a function to deal
 *           with this.
 *
 * Args:     sym      - observed symbol, e.g. `A' or `X'
 *           wt       - number of counts to distribute (e.g. 1.0)
 *           counters - array of 4 or 20 counters to increment
 *
 * Return:   Returns 1 on success and bumps the necessary counters.
 *           Returns 0 on failure and bumps each counter evenly, as
 *           if it saw a completely ambiguous symbol; this lets
 *           the caller silently accept garbage symbols, if it cares to.
 */
int
CountSymbol(char sym, float wt, float *counters)
{
  char *sptr;                   /* pointer into symbol in hmm->alphabet         */
  int   status;			/* RETURN: status; did we recognize the symbol? */
  char  symidx;			/* index of symbol in Alphabet_iupac */

  if ((sptr = strchr(Alphabet,sym)) != NULL) 
    {
      symidx = (char) (sptr - Alphabet);
      status = 1;
    }
  else 
    {
      symidx = (char) (Alphabet_iupac - 1);
      Warn("unrecognized character %c in CountSymbol()\n", sym);
      status = 0;
    }
  P7CountSymbol(counters, symidx, wt);
  return status;
}


/* Function: HMMDistance()
 * 
 * Purpose:  Test two models for how different they are, using
 *           a simple squared difference measure on all homologous
 *           parameters. They must have the same architecture:
 *           i.e. check that newhmm->M == oldhmm->M before calling.
 *           
 * Args:     newhmm  - new HMM, probability form
 *           oldhmm  - old HMM, probability form
 *           
 * Return:   distance.
 */
float
HMMDistance(struct hmm_struc *newhmm, struct hmm_struc *oldhmm)
{
  int    k,x, ts;
  float  distance = 0.0;

  for (k = 0; k <= newhmm->M; k++)
    {
				/* state transition distances */
      if (k > 0)
	{			
	  for (ts = 0; ts < 3; ts++)
	    distance += SQR( 100. * (newhmm->del[k].t[ts] - oldhmm->del[k].t[ts]));
	}
      for (ts = 0; ts < 3; ts++)
	distance += SQR( 100. * (newhmm->mat[k].t[ts] - oldhmm->mat[k].t[ts]));
      for (ts = 0; ts < 3; ts++)
	distance += SQR( 100. * (newhmm->ins[k].t[ts] - oldhmm->ins[k].t[ts]));
      
				/* symbol emission distances */
      if (k > 0)
	for (x = 0; x < Alphabet_size; x++)
	  distance += SQR( 100. * (newhmm->mat[k].p[x] - oldhmm->mat[k].p[x]));
      for (x = 0; x < Alphabet_size; x++)
	  distance += SQR( 100. * (newhmm->ins[k].p[x] - oldhmm->ins[k].p[x]));
    }
  distance = sqrt(distance) / newhmm->M;
  return distance;
}




/* Function: Renormalize()
 * 
 * Normalize all P distributions so they sum to 1.
 * P distributions that are all 0, or contain negative
 * probabilities, are left untouched.
 * 
 * Returns 1 on success, or 0 on failure.
 */
void
Renormalize(struct hmm_struc *hmm)
{
  int    k;			/* counter for states                  */

  for (k = 0; k <= hmm->M ; k++)
    {
				/* match state transition frequencies */
      FNorm(hmm->mat[k].t, 3);
      FNorm(hmm->ins[k].t, 3);
      if (k > 0) FNorm(hmm->del[k].t, 3);

      if (k > 0) FNorm(hmm->mat[k].p, Alphabet_size);
      FNorm(hmm->ins[k].p, Alphabet_size);
    }
}

