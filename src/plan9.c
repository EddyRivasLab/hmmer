/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* plan9.c
 * SRE, Wed Apr  8 07:35:30 1998
 *
 * alloc, free, and initialization of old Plan9 (HMMER 1.x) functions.
 * Rescued from the wreckage of HMMER 1.9m code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "squid.h"
#include "config.h"
#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


struct plan9_s *
P9AllocHMM(int M)               		/* length of model to make */
{
  struct plan9_s *hmm;        /* RETURN: blank HMM */

  hmm        = (struct plan9_s *)     MallocOrDie (sizeof(struct plan9_s));
  hmm->ins   = (struct basic_state *) MallocOrDie (sizeof(struct basic_state) * (M+2));
  hmm->del   = (struct basic_state *) MallocOrDie (sizeof(struct basic_state) * (M+2));
  hmm->mat   = (struct basic_state *) MallocOrDie (sizeof(struct basic_state) * (M+2));
  hmm->ref   = (char *)  MallocOrDie ((M+2) * sizeof(char));
  hmm->cs    = (char *)  MallocOrDie ((M+2) * sizeof(char));
  hmm->xray  = (float *) MallocOrDie ((M+2) * sizeof(float) * NINPUTS);
  hmm->M     = M;
  hmm->name  = Strdup("unnamed"); /* name is not optional. */

  hmm->flags = 0;
  P9ZeroHMM(hmm);
  return hmm;
}
int
P9FreeHMM(struct plan9_s *hmm)
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


/* Function: P9ZeroHMM()
 * 
 * Purpose:  Zero emission and transition counts in an HMM.
 */
void
P9ZeroHMM(struct plan9_s *hmm)
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





/* Function: P9Renormalize()
 * 
 * Normalize all P distributions so they sum to 1.
 * P distributions that are all 0, or contain negative
 * probabilities, are left untouched.
 * 
 * Returns 1 on success, or 0 on failure.
 */
void
P9Renormalize(struct plan9_s *hmm)
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

/* Function: P9DefaultNullModel()
 * 
 * Purpose:  Set up a default random sequence model, using
 *           global aafq[]'s for protein or 0.25 for nucleic 
 *           acid. randomseq is alloc'ed in caller. Alphabet information
 *           must already be known.
 */
void
P9DefaultNullModel(float *null)
{
  int x;
  if (Alphabet_type == hmmAMINO)
    for (x = 0; x < Alphabet_size; x++)
      null[x] = aafq[x];
  else if (Alphabet_type == hmmNUCLEIC)
    for (x = 0; x < Alphabet_size; x++)
      null[x] = 0.25;
  else
    Die("No support for non-protein, non-nucleic acid alphabets.");
}
