/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1997 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/


/* plan7.c
 * SRE, Sat Nov 16 14:19:56 1996
 * 
 * Support for Plan 7 HMM data structure, plan7_s.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "funcs.h"
#include "config.h"
#include "structs.h"
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Functions: AllocPlan7(), AllocPlan7Shell(), AllocPlan7Body(), FreePlan7()
 * 
 * Purpose:   Allocate or free a Plan7 HMM structure.
 *            Can either allocate all at one (AllocPlan7()) or
 *            in two steps (AllocPlan7Shell(), AllocPlan7Body()).
 *            The two step method is used in hmmio.c where we start
 *            parsing the header of an HMM file but don't 
 *            see the size of the model 'til partway thru the header.
 */
struct plan7_s *
AllocPlan7(int M) 
{
  struct plan7_s *hmm;
  int k, x;

  hmm    = (struct plan7_s *) MallocOrDie (sizeof(struct plan7_s));
  hmm->M = M;

  hmm->name   = NULL;
  hmm->desc   = NULL;
  hmm->rf     = (char *) MallocOrDie ((M+2) * sizeof(char));
  hmm->cs     = (char *) MallocOrDie ((M+2) * sizeof(char));
  hmm->comlog= NULL; 
  hmm->nseq   = 0;
  hmm->ctime  = NULL;

  hmm->t      = (float **) MallocOrDie (M     *           sizeof(float *));
  hmm->tsc    = (int **)   MallocOrDie (M     *           sizeof(int *));
  hmm->mat    = (float **) MallocOrDie ((M+1) *           sizeof(float *));
  hmm->ins    = (float **) MallocOrDie (M     *           sizeof(float *));
  hmm->msc    = (int **)   MallocOrDie (MAXCODE   *       sizeof(int *));
  hmm->isc    = (int **)   MallocOrDie (MAXCODE   *       sizeof(int *)); 
  hmm->t[0]   = (float *)  MallocOrDie ((7*M)     *       sizeof(float));
  hmm->tsc[0] = (int *)    MallocOrDie ((7*M)     *       sizeof(int));
  hmm->mat[0] = (float *)  MallocOrDie ((MAXABET*(M+1)) * sizeof(float));
  hmm->ins[0] = (float *)  MallocOrDie ((MAXABET*M) *     sizeof(float));
  hmm->msc[0] = (int *)    MallocOrDie ((MAXCODE*(M+1)) * sizeof(int));
  hmm->isc[0] = (int *)    MallocOrDie ((MAXCODE*M) *     sizeof(int));

  /* note allocation strategy for important 2D arrays -- trying
   * to keep locality as much as possible, cache efficiency etc.
   */
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * MAXABET;
    if (k < M) {
      hmm->ins[k] = hmm->ins[0] + k * MAXABET;
      hmm->t[k]   = hmm->t[0]   + k * 7;
      hmm->tsc[k] = hmm->tsc[0] + k * 7;
    }
  }
  for (x = 1; x < MAXCODE; x++) {
    hmm->msc[x] = hmm->msc[0] + x * (M+1);
    hmm->isc[x] = hmm->isc[0] + x * M;
  }
  /* tsc[0] is used as a boundary condition sometimes [Viterbi()],
   * so set to -inf always.
   */
  for (x = 0; x < 7; x++)
    hmm->tsc[0][x] = -INFTY;

  hmm->begin  = (float *) MallocOrDie  ((M+1) * sizeof(float));
  hmm->bsc    = (int *)   MallocOrDie  ((M+1) * sizeof(int));
  hmm->end    = (float *) MallocOrDie  ((M+1) * sizeof(float));
  hmm->esc    = (int *)   MallocOrDie  ((M+1) * sizeof(int));

				/* DNA translation is not enabled by default */
  hmm->dnam   = NULL;
  hmm->dnai   = NULL;
  hmm->dna2   = -INFTY;
  hmm->dna4   = -INFTY;
			/* statistical parameters set to innocuous empty values */
  hmm->mu     = 0.; 
  hmm->lambda = 0.;
  
  hmm->flags = 0;
  return hmm;
}  
struct plan7_s *
AllocPlan7Shell(void) 
{
  struct plan7_s *hmm;

  hmm    = (struct plan7_s *) MallocOrDie (sizeof(struct plan7_s));
  hmm->M = 0;

  hmm->name   = NULL;
  hmm->desc   = NULL;
  hmm->rf     = NULL;
  hmm->cs     = NULL;
  hmm->comlog = NULL; 
  hmm->nseq   = 0;
  hmm->ctime  = NULL;

  hmm->t      = NULL;
  hmm->tsc    = NULL;
  hmm->mat    = NULL;
  hmm->ins    = NULL;
  hmm->msc    = NULL;
  hmm->isc    = NULL;

  hmm->begin  = NULL;
  hmm->bsc    = NULL;
  hmm->end    = NULL;
  hmm->esc    = NULL;
				/* DNA translation is not enabled by default */
  hmm->dnam   = NULL;
  hmm->dnai   = NULL;
  hmm->dna2   = -INFTY;
  hmm->dna4   = -INFTY;
			/* statistical parameters set to innocuous empty values */
  hmm->mu     = 0.; 
  hmm->lambda = 0.;
  
  hmm->flags = 0;
  return hmm;
}  
void
AllocPlan7Body(struct plan7_s *hmm, int M) 
{
  int k, x;

  hmm->M = M;

  hmm->rf     = (char *)   MallocOrDie ((M+2) * sizeof(char));
  hmm->cs     = (char *)   MallocOrDie ((M+2) * sizeof(char));

  hmm->t      = (float **) MallocOrDie (M     *           sizeof(float *));
  hmm->tsc    = (int **)   MallocOrDie (M     *           sizeof(int *));
  hmm->mat    = (float **) MallocOrDie ((M+1) *           sizeof(float *));
  hmm->ins    = (float **) MallocOrDie (M     *           sizeof(float *));
  hmm->msc    = (int **)   MallocOrDie (MAXCODE   *       sizeof(int *));
  hmm->isc    = (int **)   MallocOrDie (MAXCODE   *       sizeof(int *)); 
  hmm->t[0]   = (float *)  MallocOrDie ((7*M)     *       sizeof(float));
  hmm->tsc[0] = (int *)    MallocOrDie ((7*M)     *       sizeof(int));
  hmm->mat[0] = (float *)  MallocOrDie ((MAXABET*(M+1)) * sizeof(float));
  hmm->ins[0] = (float *)  MallocOrDie ((MAXABET*M) *     sizeof(float));
  hmm->msc[0] = (int *)    MallocOrDie ((MAXCODE*(M+1)) * sizeof(int));
  hmm->isc[0] = (int *)    MallocOrDie ((MAXCODE*M) *     sizeof(int));

  /* note allocation strategy for important 2D arrays -- trying
   * to keep locality as much as possible, cache efficiency etc.
   */
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * MAXABET;
    if (k < M) {
      hmm->ins[k] = hmm->ins[0] + k * MAXABET;
      hmm->t[k]   = hmm->t[0]   + k * 7;
      hmm->tsc[k] = hmm->tsc[0] + k * 7;
    }
  }
  for (x = 1; x < MAXCODE; x++) {
    hmm->msc[x] = hmm->msc[0] + x * (M+1);
    hmm->isc[x] = hmm->isc[0] + x * M;
  }
  /* tsc[0] is used as a boundary condition sometimes [Viterbi()],
   * so set to -inf always.
   */
  for (x = 0; x < 7; x++)
    hmm->tsc[0][x] = -INFTY;

  hmm->begin  = (float *) MallocOrDie  ((M+1) * sizeof(float));
  hmm->bsc    = (int *)   MallocOrDie  ((M+1) * sizeof(int));
  hmm->end    = (float *) MallocOrDie  ((M+1) * sizeof(float));
  hmm->esc    = (int *)   MallocOrDie  ((M+1) * sizeof(int));

  return;
}  


void
FreePlan7(struct plan7_s *hmm)
{
  if (hmm->name    != NULL) free(hmm->name);
  if (hmm->desc    != NULL) free(hmm->desc);
  if (hmm->rf      != NULL) free(hmm->rf);
  if (hmm->cs      != NULL) free(hmm->cs);
  if (hmm->comlog  != NULL) free(hmm->comlog);
  if (hmm->ctime   != NULL) free(hmm->ctime);
  if (hmm->bsc     != NULL) free(hmm->bsc);
  if (hmm->begin   != NULL) free(hmm->begin);
  if (hmm->esc     != NULL) free(hmm->esc);
  if (hmm->end     != NULL) free(hmm->end);
  if (hmm->msc     != NULL) free(hmm->msc[0]);
  if (hmm->mat     != NULL) free(hmm->mat[0]);
  if (hmm->isc     != NULL) free(hmm->isc[0]);
  if (hmm->ins     != NULL) free(hmm->ins[0]);
  if (hmm->tsc     != NULL) free(hmm->tsc[0]);
  if (hmm->t       != NULL) free(hmm->t[0]);
  if (hmm->msc     != NULL) free(hmm->msc);
  if (hmm->mat     != NULL) free(hmm->mat);
  if (hmm->isc     != NULL) free(hmm->isc);
  if (hmm->ins     != NULL) free(hmm->ins);
  if (hmm->tsc     != NULL) free(hmm->tsc);
  if (hmm->t       != NULL) free(hmm->t);
  if (hmm->dnam    != NULL) free(hmm->dnam);
  if (hmm->dnai    != NULL) free(hmm->dnai);
  free(hmm);
}

/* Function: ZeroPlan7()
 * 
 * Purpose:  Zeros the counts/probabilities fields in a model.  
 *           Leaves null model untouched. 
 */
void
ZeroPlan7(struct plan7_s *hmm)
{
  int k;
  for (k = 1; k < hmm->M; k++)
    {
      FSet(hmm->t[k], 7, 0.);
      FSet(hmm->mat[k], Alphabet_size, 0.);
      FSet(hmm->ins[k], Alphabet_size, 0.);
    }
  FSet(hmm->mat[hmm->M], Alphabet_size, 0.);
  hmm->tbd1 = 0.;
  FSet(hmm->begin+1, hmm->M, 0.);
  FSet(hmm->end+1, hmm->M, 0.);
  for (k = 0; k < 4; k++)
    FSet(hmm->xt[k], 2, 0.);
  hmm->flags &= ~PLAN7_HASBITS;	/* invalidates scores */
  hmm->flags &= ~PLAN7_HASPROB;	/* invalidates probabilities */
}


/* Function: Plan7SetName()
 * 
 * Purpose:  Change the name of a Plan7 HMM. Convenience function.
 *      
 * Note:     Trailing whitespace and \n's are chopped.     
 */
void
Plan7SetName(struct plan7_s *hmm, char *name)
{
  if (hmm->name != NULL) free(hmm->name);
  hmm->name = Strdup(name);
  StringChop(hmm->name);
}
/* Function: Plan7SetDescription()
 * 
 * Purpose:  Change the description line of a Plan7 HMM. Convenience function.
 * 
 * Note:     Trailing whitespace and \n's are chopped.
 */
void
Plan7SetDescription(struct plan7_s *hmm, char *desc)
{
  if (hmm->desc != NULL) free(hmm->desc);
  hmm->desc = Strdup(desc);
  StringChop(hmm->desc);
}

/* Function: Plan7ComlogAppend()
 * Date:     SRE, Wed Oct 29 09:57:30 1997 [TWA 721 over Greenland] 
 * 
 * Purpose:  Concatenate command line options and append to the
 *           command line log.
 */
void
Plan7ComlogAppend(struct plan7_s *hmm, int argc, char **argv)
{
  int len;
  int i;

  /* figure out length of command line, w/ spaces and \n */
  len = argc;
  for (i = 0; i < argc; i++)
    len += strlen(argv[i]);

  /* allocate */
  if (hmm->comlog != NULL)
    {
      len += strlen(hmm->comlog);
      hmm->comlog = ReallocOrDie(hmm->comlog, sizeof(char)* (len+1));
    }
  else
    {
      hmm->comlog = MallocOrDie(sizeof(char)* (len+1));
      *(hmm->comlog) = '\0'; /* need this to make strcat work */
    }

  /* append */
  strcat(hmm->comlog, "\n");
  for (i = 0; i < argc; i++)
    {
      strcat(hmm->comlog, argv[i]);
      if (i < argc-1) strcat(hmm->comlog, " ");
    }
}

/* Function: Plan7SetCtime()
 * Date:     SRE, Wed Oct 29 11:53:19 1997 [TWA 721 over the Atlantic]
 * 
 * Purpose:  Set the ctime field in a new HMM to the current time.
 */
void
Plan7SetCtime(struct plan7_s *hmm)
{
  time_t date = time(NULL);
  if (hmm->ctime != NULL) free(hmm->ctime);
  hmm->ctime = Strdup(ctime(&date));
  StringChop(hmm->ctime);
}


/* Function: Plan7SetNullModel()
 * 
 * Purpose:  Set the null model section of an HMM.
 *           Convenience function.
 */
void
Plan7SetNullModel(struct plan7_s *hmm, float null[MAXABET], float p1)
{
  int x;
  for (x = 0; x < Alphabet_size; x++)
    hmm->null[x] = null[x];
  hmm->p1 = p1;
}


/* Function: P7Logoddsify()
 * 
 * Purpose:  Take an HMM with valid probabilities, and
 *           fill in the integer log-odds score section of the model.
 *           
 *    Notes on log-odds scores:
 *         type of parameter       probability        score
 *         -----------------       -----------        ------
 *         any emission             p_x           log_2 p_x/null_x
 *             N,J,C /assume/ p_x = null_x so /always/ score zero.  
 *         transition to emitters   t_x           log_2 t_x/p1
 *            (M,I; N,C; J)
 *             NN and CC loops are often equal to p1, so usu. score zero.
 *         C->T transition          t_x            log_2 t_x/p2 
 *            often zero, usu. C->T = p2. 
 *         all other transitions    t_x           log_2 t_x
 *             (no null model counterpart, so null prob is 1)    
 *             
 *    Notes on entry/exit scores, B->M and M->E:
 *         The probability form model includes delete states 1 and M. 
 *         these states are removed from a search form model to
 *         prevent B->D...D->E->J->B mute cycles, which would complicate
 *         dynamic programming algorithms. The data-independent
 *         S/W B->M and M->E transitions are folded together with
 *         data-dependent B->D...D->M and M->D...D->E paths.
 *         
 *         This process is referred to in the code as "wing folding"
 *         or "wing retraction"... the analogy is to a swept-wing
 *         fighter in landing vs. high speed flight configuration.
 *         
 *    Note on Viterbi vs. forward flag:     
 *         Wing retraction must take forward vs. Viterbi
 *         into account. If forward, sum two paths; if Viterbi, take
 *         max. I tried to slide this by as a sum, without
 *         the flag, but Alex detected it as a bug, because you can
 *         then find cases where the Viterbi score doesn't match
 *         the P7TraceScore().
 *             
 * Args:      hmm          - the hmm to calculate scores in.
 *            viterbi_mode - TRUE to fold wings in Viterbi configuration.
 *                  
 * Return:    (void)
 *            hmm scores are filled in.
 */  
void
P7Logoddsify(struct plan7_s *hmm, int viterbi_mode)
{
  int k;			/* counter for model position */
  int x;			/* counter for symbols        */
  float accum;
  float tbm, tme;

  if (hmm->flags & PLAN7_HASBITS) return;

  /* Symbol emission scores
   */
  for (k = 1; k <= hmm->M; k++) 
    {
				/* match/insert emissions in main model */
      for (x = 0; x < Alphabet_size; x++) 
	{
	  hmm->msc[x][k] = Prob2Score(hmm->mat[k][x], hmm->null[x]);
	  if (k < hmm->M) 
	    hmm->isc[x][k] =  Prob2Score(hmm->ins[k][x], hmm->null[x]); 
	}
				/* degenerate match/insert emissions */
      for (x = Alphabet_size; x < Alphabet_iupac; x++) 
	{
	  hmm->msc[x][k] = DegenerateSymbolScore(hmm->mat[k], hmm->null, x);
	  if (k < hmm->M)
	    hmm->isc[x][k] = DegenerateSymbolScore(hmm->ins[k], hmm->null, x);
	}
    }

  /* State transitions.
   * 
   * A note on "folding" of D_1 and D_M.
   * These two delete states are folded out of search form models
   * in order to prevent null cycles in the dynamic programming
   * algorithms (see code below). However, we use their log transitions
   * when we save the model! So the following log transition probs
   * are used *only* in save files, *never* in search algorithms:
   *    log (tbd1), D1 -> M1, D1 -> D2
   *    Mm-1 -> Dm, Dm-1 -> Dm
   *    
   * In a search algorithm, these have to be interpreted as -INFTY    
   * because their contributions are folded into bsc[] and esc[]
   * entry/exit scores. They can't be set to -INFTY here because
   * we need them in save files.
   */
  for (k = 1; k < hmm->M; k++)
    {
      hmm->tsc[k][TMM] = Prob2Score(hmm->t[k][TMM], hmm->p1);
      hmm->tsc[k][TMI] = Prob2Score(hmm->t[k][TMI], hmm->p1);
      hmm->tsc[k][TMD] = Prob2Score(hmm->t[k][TMD], 1.0);
      hmm->tsc[k][TIM] = Prob2Score(hmm->t[k][TIM], hmm->p1);
      hmm->tsc[k][TII] = Prob2Score(hmm->t[k][TII], hmm->p1);
      hmm->tsc[k][TDM] = Prob2Score(hmm->t[k][TDM], hmm->p1);
      hmm->tsc[k][TDD] = Prob2Score(hmm->t[k][TDD], 1.0);
    }

  /* B->M entry transitions. Note how D_1 is folded out.
   * M1 is just B->M1
   * M2 is sum (or max) of B->M2 and B->D1->M2
   * M_k is sum (or max) of B->M_k and B->D1...D_k-1->M_k
   */
  accum = hmm->tbd1;
  for (k = 1; k <= hmm->M; k++)
    {
      tbm = hmm->begin[k];	/* B->M_k part */
      if (k > 1) {		/* B->D1...D_k-1->M_k part we get from accum*/
	if (viterbi_mode) tbm =  MAX(tbm, accum * hmm->t[k-1][TDM]);
	else              tbm += accum * hmm->t[k-1][TDM];
	accum *= hmm->t[k-1][TDD];
      }
      hmm->bsc[k] = Prob2Score(tbm, hmm->p1); /* entries are to emitters */
    }

  /* M->E exit transitions. Note how D_M is folded out.
   * M_M is 1 by definition
   * M_M-1 is sum of M_M-1->E and M_M-1->D_M->E, where D_M->E is 1 by definition
   * M_k is sum of M_k->E and M_k->D_k+1...D_M->E
   */
  hmm->esc[hmm->M] = 0;
  accum = 1.;
  for (k = hmm->M-1; k >= 1; k--)
    {
      if (viterbi_mode) tme = MAX(hmm->end[k], accum * hmm->t[k][TMD]);
      else              tme = hmm->end[k] + accum * hmm->t[k][TMD];

      accum      *= hmm->t[k][TDD];
      hmm->esc[k] = Prob2Score(tme, 1.); /* exits are to non-emitter */
    }

				/* special transitions */
  hmm->xsc[XTN][LOOP] = Prob2Score(hmm->xt[XTN][LOOP], hmm->p1);
  hmm->xsc[XTN][MOVE] = Prob2Score(hmm->xt[XTN][MOVE], 1.0);
  hmm->xsc[XTE][LOOP] = Prob2Score(hmm->xt[XTE][LOOP], 1.0);
  hmm->xsc[XTE][MOVE] = Prob2Score(hmm->xt[XTE][MOVE], 1.0);
  hmm->xsc[XTC][LOOP] = Prob2Score(hmm->xt[XTC][LOOP], hmm->p1);
  hmm->xsc[XTC][MOVE] = Prob2Score(hmm->xt[XTC][MOVE], 1.-hmm->p1);
  hmm->xsc[XTJ][LOOP] = Prob2Score(hmm->xt[XTJ][LOOP], hmm->p1);
  hmm->xsc[XTJ][MOVE] = Prob2Score(hmm->xt[XTJ][MOVE], 1.0);

  hmm->flags |= PLAN7_HASBITS;	/* raise the log-odds ready flag */
}



/* Function: Plan7Renormalize()
 * 
 * Purpose:  Take an HMM in counts form, and renormalize
 *           all of its probability vectors. Also enforces
 *           Plan7 restrictions on nonexistent transitions.
 *           
 * Args:     hmm - the model to renormalize.
 *                 
 * Return:   (void)
 *           hmm is changed.
 */                          
void
Plan7Renormalize(struct plan7_s *hmm)
{
  int   k;			/* counter for model position */
  int   st;			/* counter for special states */
  float d;			/* denominator */

				/* match emissions */
  for (k = 1; k <= hmm->M; k++) 
    FNorm(hmm->mat[k], Alphabet_size);
				/* insert emissions */
  for (k = 1; k < hmm->M; k++)
    FNorm(hmm->ins[k], Alphabet_size);
				/* begin transitions */
  d = FSum(hmm->begin+1, hmm->M) + hmm->tbd1;
  FScale(hmm->begin+1, hmm->M, 1./d);
  hmm->tbd1 /= d;
				/* main model transitions */
  for (k = 1; k < hmm->M; k++)
    {
      d = FSum(hmm->t[k], 3) + hmm->end[k]; 
      FScale(hmm->t[k], 3, 1./d);
      hmm->end[k] /= d;

      FNorm(hmm->t[k]+3, 2);	/* insert */
      FNorm(hmm->t[k]+5, 2);	/* delete */
    }
				/* null model emissions */
  FNorm(hmm->null, Alphabet_size);
				/* special transitions  */
  for (st = 0; st < 4; st++)
    FNorm(hmm->xt[st], 2);
				/* enforce nonexistent transitions */
				/* (is this necessary?) */
  hmm->t[0][TDM] = hmm->t[0][TDD] = 0.0;

  hmm->flags &= ~PLAN7_HASBITS;	/* clear the log-odds ready flag */
  hmm->flags |= PLAN7_HASPROB;	/* set the probabilities OK flag */
}
  

/* Function: Plan7NakedConfig()
 * 
 * Purpose:  Set the alignment-independent, algorithm-dependent parameters
 *           of a Plan7 model so that no special states (N,C,J) emit anything:
 *           one simple, full global pass through the model.
 * 
 * Args:     hmm - the plan7 model
 *                 
 * Return:   (void)
 *           The HMM is modified; algorithm dependent parameters are set.
 *           Previous scores are invalidated if they existed.
 */
void
Plan7NakedConfig(struct plan7_s *hmm)                           
{
  hmm->xt[XTN][MOVE] = 1.;	      /* disallow N-terminal tail */
  hmm->xt[XTN][LOOP] = 0.;
  hmm->xt[XTE][MOVE] = 1.;	      /* only 1 domain/sequence ("global" alignment) */
  hmm->xt[XTE][LOOP] = 0.;
  hmm->xt[XTC][MOVE] = 1.;	      /* disallow C-terminal tail */
  hmm->xt[XTC][LOOP] = 0.;
  hmm->xt[XTJ][MOVE] = 0.;	      /* J state unused */
  hmm->xt[XTJ][LOOP] = 1.;
  FSet(hmm->begin+2, hmm->M-1, 0.);   /* disallow internal entries. */
  hmm->begin[1]    = 1. - hmm->tbd1;
  FSet(hmm->end+1,   hmm->M-1, 0.);   /* disallow internal exits. */
  hmm->end[hmm->M] = 1.;
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}
   
/* Function: Plan7GlobalConfig()
 * 
 * Purpose:  Set the alignment-independent, algorithm-dependent parameters
 *           of a Plan7 model to global (Needleman/Wunsch) configuration.
 * 
 *           Like a non-looping hmmls, since we actually allow flanking
 *           N and C terminal sequence. 
 *           
 * Args:     hmm - the plan7 model
 *                 
 * Return:   (void)
 *           The HMM is modified; algorithm dependent parameters are set.
 *           Previous scores are invalidated if they existed.
 */
void
Plan7GlobalConfig(struct plan7_s *hmm)                           
{
  hmm->xt[XTN][MOVE] = 1. - hmm->p1;  /* allow N-terminal tail */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 1.;	      /* only 1 domain/sequence ("global" alignment) */
  hmm->xt[XTE][LOOP] = 0.;
  hmm->xt[XTC][MOVE] = 1. - hmm->p1;  /* allow C-terminal tail */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 0.;	      /* J state unused */
  hmm->xt[XTJ][LOOP] = 1.;
  FSet(hmm->begin+2, hmm->M-1, 0.);   /* disallow internal entries. */
  hmm->begin[1]    = 1. - hmm->tbd1;
  FSet(hmm->end+1,   hmm->M-1, 0.);   /* disallow internal exits. */
  hmm->end[hmm->M] = 1.;
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}
   
/* Function: Plan7LSConfig()
 * 
 * Purpose:  Set the alignment independent parameters of a Plan7 model
 *           to hmmls (global in HMM, local in sequence) configuration.
 *           
 * Args:     hmm  - the plan7 model
 *                 
 * Return:   (void);
 *           the HMM probabilities are modified.
 */
void
Plan7LSConfig(struct plan7_s *hmm)
{
  hmm->xt[XTN][MOVE] = 1.-hmm->p1;    /* allow N-terminal tail */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 0.5;	     /* expectation 2 domains/seq */
  hmm->xt[XTE][LOOP] = 0.5;
  hmm->xt[XTC][MOVE] = 1.-hmm->p1;    /* allow C-terminal tail */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 1.-hmm->p1;   /* allow J junction state */
  hmm->xt[XTJ][LOOP] = hmm->p1;
  FSet(hmm->begin+2, hmm->M-1, 0.);  /* start at M1/D1 */
  hmm->begin[1]    = 1. - hmm->tbd1;
  FSet(hmm->end+1,   hmm->M-1, 0.);  /* end at M_m/D_m */
  hmm->end[hmm->M] = 1.;
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}  
                             

/* Function: Plan7SWConfig()
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a Plan7 model to hmmsw (Smith/Waterman) configuration.
 *           
 * Notes:    entry/exit is asymmetric because of the left/right
 *           nature of the HMM/profile. Entry probability is distributed
 *           simply by assigning p_x = pentry / (M-1) to M-1 
 *           internal match states. However, the same approach doesn't
 *           lead to a flat distribution over exit points. Exit p's
 *           must be corrected for the probability of a previous exit
 *           from the model. Requiring a flat distribution over exit
 *           points leads to an easily solved piece of algebra, giving:
 *                      p_1 = pexit / (M-1)
 *                      p_x = p_1 / (1 - (x-1) p_1)
 *           
 * Args:     hmm    - the Plan7 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over M-1 match states
 *           pexit  - probability of an internal exit somewhere; 
 *                    will be distributed over M-1 match states.
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
Plan7SWConfig(struct plan7_s *hmm, float pentry, float pexit)
{
  float basep;			/* p1 for exits: the base p */
  int   k;			/* counter over states      */
  float d;

  /* Configure special states.
   */
  hmm->xt[XTN][MOVE] = 1-hmm->p1;    /* allow N-terminal tail */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 1.;	     /* disallow jump state   */
  hmm->xt[XTE][LOOP] = 0.;
  hmm->xt[XTC][MOVE] = 1-hmm->p1;    /* allow C-terminal tail */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 1.;           /* J is unused */
  hmm->xt[XTJ][LOOP] = 0.;

  /* Configure entry.
   */
  hmm->begin[1] = (1. - pentry) * (1. - hmm->tbd1);
  FSet(hmm->begin+2, hmm->M-1, (pentry * (1.- hmm->tbd1)) / (float)(hmm->M-1));
  
  /* Configure exit.
   */
  hmm->end[hmm->M] = 1.0;
  basep = pexit / (float) (hmm->M-1);
  for (k = 1; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) (k-1));
  
  /* Renormalize transitions to deal with new exits.
   */
  for (k = 1; k < hmm->M; k++)
    {
      d = FSum(hmm->t[k], 3);
      FScale(hmm->t[k], 3, 1./(d + d*hmm->end[k]));
    }

  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}

/* Function: Plan7FSConfig()
 * Date:     SRE, Fri Jan  2 15:34:40 1998 [StL]
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a Plan7 model to hmmfs (multihit Smith/Waterman) configuration.
 *           
 *           See comments on Plan7SWConfig() for explanation of
 *           how pentry and pexit are used.
 *           
 * Args:     hmm    - the Plan7 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over M-1 match states
 *           pexit  - probability of an internal exit somewhere; 
 *                    will be distributed over M-1 match states.
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
Plan7FSConfig(struct plan7_s *hmm, float pentry, float pexit)
{
  float basep;			/* p1 for exits: the base p */
  int   k;			/* counter over states      */
  float d;

  /* Configure special states.
   */
  hmm->xt[XTN][MOVE] = 1-hmm->p1;    /* allow N-terminal tail     */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 0.5;	     /* allow loops / multihits   */
  hmm->xt[XTE][LOOP] = 0.5;
  hmm->xt[XTC][MOVE] = 1-hmm->p1;    /* allow C-terminal tail     */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 1.-hmm->p1;   /* allow J junction between domains */
  hmm->xt[XTJ][LOOP] = hmm->p1;

  /* Configure entry.
   */
  hmm->begin[1] = (1. - pentry) * (1. - hmm->tbd1);
  FSet(hmm->begin+2, hmm->M-1, (pentry * (1.-hmm->tbd1)) / (float)(hmm->M-1));
  
  /* Configure exit.
   */
  hmm->end[hmm->M] = 1.0;
  basep = pexit / (float) (hmm->M-1);
  for (k = 1; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) (k-1));
  
  /* Renormalize transitions to deal with new exits.
   */
  for (k = 1; k < hmm->M; k++)
    {
      d = FSum(hmm->t[k], 3);
      FScale(hmm->t[k], 3, 1./(d + d*hmm->end[k]));
    }

  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}




/* Function: Plan7ESTConfig()
 * 
 * Purpose:  Configure a Plan7 model for EST Smith/Waterman
 *           analysis.
 *           
 *           OUTDATED; DO NOT USE WITHOUT RECHECKING
 *           
 * Args:     hmm        - hmm to configure.
 *           aacode     - 0..63 vector mapping genetic code to amino acids
 *           estmodel   - 20x64 translation matrix, w/ codon bias and substitution error
 *           dna2       - probability of a -1 frameshift in a triplet
 *           dna4       - probability of a +1 frameshift in a triplet     
 */ 
void
Plan7ESTConfig(struct plan7_s *hmm, int *aacode, float **estmodel, 
	       float dna2, float dna4)
{
  int k;
  int x;
  float p;
  float *tripnull;		/* UNFINISHED!!! */

				/* configure specials */
  hmm->xt[XTN][MOVE] = 1./351.;
  hmm->xt[XTN][LOOP] = 350./351.;
  hmm->xt[XTE][MOVE] = 1.;
  hmm->xt[XTE][LOOP] = 0.;
  hmm->xt[XTC][MOVE] = 1./351.;
  hmm->xt[XTC][LOOP] = 350./351.;
  hmm->xt[XTJ][MOVE] = 1.;
  hmm->xt[XTJ][LOOP] = 0.;
				/* configure entry/exit */
  hmm->begin[1] = 0.5;
  FSet(hmm->begin+2, hmm->M-1, 0.5 / ((float)hmm->M - 1.));
  hmm->end[hmm->M] = 1.;
  FSet(hmm->end, hmm->M-1, 0.5 / ((float)hmm->M - 1.));

				/* configure dna triplet/frameshift emissions */
  for (k = 1; k <= hmm->M; k++)
    {
				/* translate aa to triplet probabilities */
      for (x = 0; x < 64; x++) {
	p =  hmm->mat[k][aacode[x]] * estmodel[aacode[x]][x] * (1.-dna2-dna4);
	hmm->dnam[x][k] = Prob2Score(p, tripnull[x]);

	p = hmm->ins[k][aacode[x]] * estmodel[aacode[x]][x] * (1.-dna2-dna4);
	hmm->dnai[x][k] = Prob2Score(p, tripnull[x]);
      }
      hmm->dnam[64][k] = 0;	/* ambiguous codons score 0 (danger?) */
      hmm->dna2 = Prob2Score(dna2, 1.);
      hmm->dna4 = Prob2Score(dna4, 1.);
    }
}
	  
/* Function: PrintPlan7Stats()
 * 
 * Purpose:  Given a newly constructed HMM and the tracebacks
 *           of the sequences it was trained on, print out all
 *           the interesting information at the end of hmmb 
 *           and hmmt runs that convinces the user we actually
 *           did something.
 *           
 * Args:     fp   - where to send the output (stdout, usually)
 *           hmm  - the new HMM, probability form
 *           dsq  - digitized training seqs
 *           nseq - number of dsq's
 *           tr   - array of tracebacks for dsq
 *                  
 * Return:   (void)
 */
void
PrintPlan7Stats(FILE *fp, struct plan7_s *hmm, char **dsq, int nseq,
		struct p7trace_s **tr)
{
  int   idx;			/* counter for sequences                */
  float score;			/* an individual trace score            */
  float total, best, worst;	/* for the avg. and range of the scores */
  float sqsum, stddev;		/* for the std. deviation of the scores */

  P7Logoddsify(hmm, TRUE);	/* make sure model scores are ready */

				/* find individual trace scores */
  score = P7TraceScore(hmm, dsq[0], tr[0]);
  total = best = worst = score;
  sqsum = score * score;
  for (idx = 1; idx < nseq; idx++) {
    /* P7PrintTrace(stdout, tr[idx], hmm, dsq[idx]); */
    score  = P7TraceScore(hmm, dsq[idx], tr[idx]);
    total += score;
    sqsum += score * score;
    if (score > best)  best = score;
    if (score < worst) worst = score;
  }
  if (nseq > 1) {
    stddev = (sqsum - (total * total / (float) nseq)) / ((float) nseq - 1.);
    stddev = (stddev > 0) ? sqrt(stddev) : 0.0;
  } else stddev = 0.0;
				/* print out stuff. */
  fprintf(fp, "Average score:  %10.2f bits\n", total / (float) nseq);
  fprintf(fp, "Minimum score:  %10.2f bits\n", worst);
  fprintf(fp, "Maximum score:  %10.2f bits\n", best);
  fprintf(fp, "Std. deviation: %10.2f bits\n", stddev);
}

/* Function: DegenerateSymbolScore()
 * 
 * Purpose:  Given a sequence character x and an hmm emission probability
 *           vector, calculate the log-odds (base 2) score of
 *           the symbol.
 *          
 *           Easy if x is in the emission alphabet, but not so easy
 *           is x is a degenerate symbol. The "correct" Bayesian
 *           philosophy is to calculate score(X) by summing over
 *           p(x) for all x in the degenerate symbol X to get P(X),
 *           doing the same sum over the prior to get F(X), and
 *           doing log_2 (P(X)/F(X)). This gives an X a zero score,
 *           for instance.
 *           
 *           Though this is correct in a formal Bayesian sense --
 *           we have no information on the sequence, so we can't
 *           say if it's random or model, so it scores zero --
 *           it sucks, big time, for scoring biological sequences.
 *           Sequences with lots of X's score near zero, while
 *           real sequences have average scores that are negative --
 *           so the X-laden sequences appear to be lifted out
 *           of the noise of a full histogram of a database search.
 *           Correct or not, this is highly undesirable.
 *           
 *           So therefore we calculated the expected score of
 *           the degenerate symbol by summing over all x in X:
 *                 e_x log_2 (p(x)/f(x))
 *           where the expectation of x, e_x, is calculated from
 *           the random model.
 *
 *           Empirically, this works; it also has a wooly hand-waving
 *           probabilistic justification that I'm happy enough about.
 *           
 * Args:     p      - probabilities of normal symbols
 *           null   - null emission model
 *           ambig  - index of the degenerate character in Alphabet[]
 *                    
 * Return:   the integer log odds score of x given the emission
 *           vector and the null model, scaled up by INTSCALE.              
 */
int 
DegenerateSymbolScore(float *p, float *null, int ambig)
{
  int x;
  float numer = 0.;
  float denom = 0.;

  for (x = 0; x < Alphabet_size; x++) {
    if (Degenerate[ambig][x]) {
      numer += null[x] * sreLOG2(p[x] / null[x]);
      denom += null[x];
    }
  }
  return (int) (INTSCALE * numer / denom);
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


