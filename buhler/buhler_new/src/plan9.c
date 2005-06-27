/* plan9.c
 *
 * alloc, free, and initialization of old Plan9 (HMMER 1.x) functions.
 * Rescued from the wreckage of HMMER 1.9m code.
 *
 * SRE, Wed Apr  8 07:35:30 1998
 * SVN $Id$
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "squid.h"

#include "plan7.h"
#include "structs.h"
#include "funcs.h"


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


/*****************************************************************
 * 
 * Plan9/Plan7 interface
 * 
 * Very important code during the evolutionary takeover by Plan7 --
 * convert between Krogh/Haussler and Plan7 models.
 *****************************************************************/

/* Function: Plan9toPlan7()
 * 
 * Purpose:  Convert an old HMM into Plan7. Configures it in
 *           ls mode.
 *           
 * Args:     hmm       - old ugly plan9 style HMM
 *           ret_plan7 - new wonderful Plan7 HMM
 *           
 * Return:   (void)    
 *           Plan7 HMM is allocated here. Free w/ FreePlan7().
 */               
void
Plan9toPlan7(struct plan9_s *hmm, struct plan7_s **ret_plan7)
{
  struct plan7_s *plan7;
  int k, x;

  plan7 = AllocPlan7(hmm->M);
  
  for (k = 1; k < hmm->M; k++)
    {
      plan7->t[k][TMM] = hmm->mat[k].t[MATCH];
      plan7->t[k][TMD] = hmm->mat[k].t[DELETE];
      plan7->t[k][TMI] = hmm->mat[k].t[INSERT];
      plan7->t[k][TDM] = hmm->del[k].t[MATCH];
      plan7->t[k][TDD] = hmm->del[k].t[DELETE];
      plan7->t[k][TIM] = hmm->ins[k].t[MATCH];
      plan7->t[k][TII] = hmm->ins[k].t[INSERT];
    }

  for (k = 1; k <= hmm->M; k++)
    for (x = 0; x < Alphabet_size; x++)
      plan7->mat[k][x] = hmm->mat[k].p[x];

  for (k = 1; k < hmm->M; k++)
    for (x = 0; x < Alphabet_size; x++)
      plan7->ins[k][x] = hmm->ins[k].p[x];

  plan7->tbd1 = hmm->mat[0].t[DELETE] / (hmm->mat[0].t[DELETE] + hmm->mat[0].t[MATCH]);
  
		/* We have to make up the null transition p1; use default */
  P7DefaultNullModel(plan7->null, &(plan7->p1));
  for (x = 0; x < Alphabet_size; x++)
    plan7->null[x] = hmm->null[x];
      
  if (hmm->name != NULL) 
    Plan7SetName(plan7, hmm->name);
  if (hmm->flags & HMM_REF) {
    strcpy(plan7->rf, hmm->ref);
    plan7->flags |= PLAN7_RF;
  }
  if (hmm->flags & HMM_CS) {
    strcpy(plan7->cs, hmm->cs);
    plan7->flags |= PLAN7_CS;
  }

  Plan7LSConfig(plan7);		/* configure specials for ls-style alignment */
  Plan7Renormalize(plan7);	/* mainly to correct for missing ID and DI */
  plan7->flags |= PLAN7_HASPROB;	/* probabilities are valid */
  plan7->flags &= ~PLAN7_HASBITS;	/* scores are not valid    */
  *ret_plan7 = plan7;
}



/************************************************************
 * @LICENSE@
 ************************************************************/

