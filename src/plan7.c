/* plan7.c
 * Support for Plan 7 HMM data structure, plan7_s.
 * 
 * SRE, Sat Nov 16 14:19:56 1996
 * SVN $Id$
 */

#include "config.h"		/* must be included first */
#include "squidconf.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "squid.h"

#include "plan7.h"		/* the model structure */
#include "funcs.h"
#include "structs.h"


static char *score2ascii(int sc);

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

  hmm = AllocPlan7Shell();
  AllocPlan7Body(hmm, M);
  return hmm;
}  
struct plan7_s *
AllocPlan7Shell(void) 
{
  struct plan7_s *hmm;

  hmm    = (struct plan7_s *) MallocOrDie (sizeof(struct plan7_s));

  hmm->M    = 0;
  hmm->t    = NULL;
  hmm->mat  = NULL;
  hmm->ins  = NULL;

  hmm->begin  = NULL;
  hmm->end    = NULL;

  hmm->mode    = P7_NO_MODE;
  hmm->tsc     = hmm->msc     = hmm->isc     = NULL;
  hmm->tsc_mem = hmm->msc_mem = hmm->msc_mem = NULL;
  hmm->bsc = hmm->bsc_mem = NULL;
  hmm->esc = hmm->esc_mem = NULL;

  hmm->name     = NULL;
  hmm->acc      = NULL;
  hmm->desc     = NULL;
  hmm->rf       = NULL;
  hmm->cs       = NULL;
  hmm->ca       = NULL;
  hmm->comlog   = NULL; 
  hmm->nseq     = 0;
  hmm->ctime    = NULL;
  hmm->map      = NULL;
  hmm->checksum = 0;

  hmm->tpri = NULL;
  hmm->mpri = NULL;
  hmm->ipri = NULL;

  hmm->ga1 = hmm->ga2 = 0.0;
  hmm->tc1 = hmm->tc2 = 0.0;
  hmm->nc1 = hmm->nc2 = 0.0;
				/* DNA translation is not enabled by default */
  hmm->dnam   = NULL;
  hmm->dnai   = NULL;
  hmm->dna2   = -INFTY;
  hmm->dna4   = -INFTY;
			/* statistical parameters set to innocuous empty values */
  hmm->mu     = 0.; 
  hmm->lambda = 0.;
  hmm->kappa  = 0.;
  hmm->sigma  = 0.;
  hmm->Lbase  = 0;
  
  hmm->flags = 0;
  return hmm;
}  
#ifndef ALTIVEC  /* in Altivec port, this func is replaced in fast_algorithms.c */
void
AllocPlan7Body(struct plan7_s *hmm, int M) 
{
  int k, x;

  hmm->M = M;
  hmm->t      = MallocOrDie (M     *           sizeof(float *));
  hmm->mat    = MallocOrDie ((M+1) *           sizeof(float *));
  hmm->ins    = MallocOrDie (M     *           sizeof(float *));
  hmm->t[0]   = MallocOrDie ((7*M)     *       sizeof(float));
  hmm->mat[0] = MallocOrDie ((MAXABET*(M+1)) * sizeof(float));
  hmm->ins[0] = MallocOrDie ((MAXABET*M) *     sizeof(float));
  /* note allocation strategy for important 2D arrays -- trying
   * to keep locality as much as possible, cache efficiency etc.
   */
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * MAXABET;
    if (k < M) {
      hmm->ins[k] = hmm->ins[0] + k * MAXABET;
      hmm->t[k]   = hmm->t[0]   + k * 7;
    }
  }

  hmm->begin  = MallocOrDie  ((M+1) * sizeof(float));
  hmm->end    = MallocOrDie  ((M+1) * sizeof(float));

  hmm->tsc     = MallocOrDie (7     *           sizeof(int *));
  hmm->msc     = MallocOrDie (MAXCODE   *       sizeof(int *));
  hmm->isc     = MallocOrDie (MAXCODE   *       sizeof(int *)); 
  hmm->tsc_mem = MallocOrDie ((7*M)     *       sizeof(int));
  hmm->msc_mem = MallocOrDie ((MAXCODE*(M+1)) * sizeof(int));
  hmm->isc_mem = MallocOrDie ((MAXCODE*M) *     sizeof(int));
  hmm->tsc[0]  = hmm->tsc_mem;
  hmm->msc[0]  = hmm->msc_mem;
  hmm->isc[0]  = hmm->isc_mem;

  for (x = 1; x < MAXCODE; x++) {
    hmm->msc[x] = hmm->msc[0] + x * (M+1);
    hmm->isc[x] = hmm->isc[0] + x * M;
  }
  for (x = 0; x < 7; x++)
    hmm->tsc[x] = hmm->tsc[0] + x * M;

  /* tsc[x][0] is used as a boundary condition sometimes [Viterbi()],
   * so init (and keep) at -inf always.
   */
  for (x = 0; x < 7; x++)
    hmm->tsc[x][0] = -INFTY;

  hmm->bsc_mem  = MallocOrDie  ((M+1) * sizeof(int));
  hmm->esc_mem  = MallocOrDie  ((M+1) * sizeof(int));
  hmm->bsc = hmm->bsc_mem;
  hmm->esc = hmm->esc_mem;

  hmm->rf     = MallocOrDie ((M+2) * sizeof(char));
  hmm->cs     = MallocOrDie ((M+2) * sizeof(char));
  hmm->ca     = MallocOrDie ((M+2) * sizeof(char));
  hmm->map    = MallocOrDie ((M+1) * sizeof(int));

  return;
}  
#endif /* not the ALTIVEC version */

void
FreePlan7(struct plan7_s *hmm)
{
  if (hmm->mat     != NULL) free(hmm->mat[0]);
  if (hmm->ins     != NULL) free(hmm->ins[0]);
  if (hmm->t       != NULL) free(hmm->t[0]);
  if (hmm->mat     != NULL) free(hmm->mat);
  if (hmm->ins     != NULL) free(hmm->ins);
  if (hmm->t       != NULL) free(hmm->t);
  if (hmm->begin   != NULL) free(hmm->begin);
  if (hmm->end     != NULL) free(hmm->end);
  if (hmm->bsc_mem != NULL) free(hmm->bsc_mem);
  if (hmm->esc_mem != NULL) free(hmm->esc_mem);
  if (hmm->msc_mem != NULL) free(hmm->msc_mem);
  if (hmm->isc_mem != NULL) free(hmm->isc_mem);
  if (hmm->tsc_mem != NULL) free(hmm->tsc_mem);
  if (hmm->msc     != NULL) free(hmm->msc);
  if (hmm->isc     != NULL) free(hmm->isc);
  if (hmm->tsc     != NULL) free(hmm->tsc);
  if (hmm->name    != NULL) free(hmm->name);
  if (hmm->acc     != NULL) free(hmm->acc);
  if (hmm->desc    != NULL) free(hmm->desc);
  if (hmm->rf      != NULL) free(hmm->rf);
  if (hmm->cs      != NULL) free(hmm->cs);
  if (hmm->ca      != NULL) free(hmm->ca);
  if (hmm->comlog  != NULL) free(hmm->comlog);
  if (hmm->ctime   != NULL) free(hmm->ctime);
  if (hmm->map     != NULL) free(hmm->map);
  if (hmm->tpri    != NULL) free(hmm->tpri);
  if (hmm->mpri    != NULL) free(hmm->mpri);
  if (hmm->ipri    != NULL) free(hmm->ipri);
  if (hmm->dnam    != NULL) free(hmm->dnam);
  if (hmm->dnai    != NULL) free(hmm->dnai);
  free(hmm);
}

/* Function: ZeroPlan7()
 * 
 * Purpose:  Zeros the counts/probabilities fields in a model (both
 *           core and configured form).  
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
  hmm->tbm1 = 0.;
  FSet(hmm->begin+1, hmm->M, 0.);
  FSet(hmm->end+1, hmm->M, 0.);
  for (k = 0; k < 4; k++)
    FSet(hmm->xt[k], 2, 0.);
  
  hmm->kappa  = 0;
  hmm->sigma  = 0;
  hmm->mode   = P7_NO_MODE;
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
/* Function: Plan7SetAccession()
 * 
 * Purpose:  Change the accession number of a Plan7 HMM. Convenience function.
 *      
 * Note:     Trailing whitespace and \n's are chopped.     
 */
void
Plan7SetAccession(struct plan7_s *hmm, char *acc)
{
  if (hmm->acc != NULL) free(hmm->acc);
  hmm->acc = Strdup(acc);
  StringChop(hmm->acc);
  hmm->flags |= PLAN7_ACC;
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
  hmm->flags |= PLAN7_DESC;
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

/* Function:  Plan7Rescale()
 * Incept:    Steve Johnson, 3 May 2004
 *            eweights code incorp: SRE, Thu May 20 10:34:03 2004 [St. Louis]
 *
 * Purpose:   Scale a counts-based HMM by some factor, for
 *            adjusting counts to a new effective sequence number.
 *
 * Args:      hmm        - counts based HMM.
 *            scale      - scaling factor (e.g. eff_nseq/nseq); 1.0= no scaling.
 *
 * Returns:   (void)
 */
void 
Plan7Rescale(struct plan7_s *hmm, float scale)
{
  int k;
  int st;

  /* emissions and transitions in the main model.
   * Note that match states are 1..M, insert states are 1..M-1,
   * and only nodes 1..M-1 have a valid array of transitions.
   */
  for(k = 1; k <= hmm->M; k++) 
    FScale(hmm->mat[k], Alphabet_size, scale);
  for(k = 1; k <  hmm->M; k++) 
    FScale(hmm->ins[k], Alphabet_size, scale);
  for(k = 1; k <  hmm->M; k++) 
    FScale(hmm->t[k],   7,             scale);

  /* begin, end transitions; only valid [1..M] */
  FScale(hmm->begin+1, hmm->M, scale);
  FScale(hmm->end+1,   hmm->M, scale);
  
  /* B->D1 transition */
  hmm->tbd1 *= scale;

  /* special transitions */
  for (st = 0; st < 4; st++)
    FScale(hmm->xt[st], 2, scale);

  return;
}


/* Function: Plan7Renormalize()
 * 
 * Purpose:  Take an HMM in counts form, and renormalize
 *           all of its probability vectors. Also enforces
 *           Plan7 restrictions on nonexistent transitions.
 *           
 *           Note: only the core probability model is renormalized.
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
  float d;			/* denominator */

				/* match emissions */
  for (k = 1; k <= hmm->M; k++) 
    FNorm(hmm->mat[k], Alphabet_size);
				/* insert emissions */
  for (k = 1; k < hmm->M; k++)
    FNorm(hmm->ins[k], Alphabet_size);
                                /* tbd1,tbm1 */
  d = hmm->tbd1 + hmm->tbm1;
  hmm->tbm1 /= d;
  hmm->tbd1 /= d;
				/* main model transitions */
  for (k = 1; k < hmm->M; k++)
    {
      FNorm(hmm->t[k],   3);	/* match  */
      FNorm(hmm->t[k]+3, 2);	/* insert */
      FNorm(hmm->t[k]+5, 2);	/* delete */
    }
				/* enforce nonexistent transitions */
				/* (is this necessary?) */
  hmm->t[0][TDM] = hmm->t[0][TDD] = 0.0;

  hmm->flags &= ~PLAN7_HASBITS;	/* clear the log-odds ready flag */
  hmm->flags |= PLAN7_HASPROB;	/* set the probabilities OK flag */
}
  


#ifdef SRE_REMOVED
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
#endif /*SRE_REMOVED*/
	  

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

/* Function:  Plan7_DumpScores()
 * Incept:    SRE, Fri May  6 08:37:17 2005 [St. Louis]
 *
 * Purpose:   Debugging: print log-odds scores of a configured plan7
 *            model <hmm> to a stream <fp>, in roughly the same format
 *            as a save file.  
 */
void
Plan7_DumpScores(FILE *fp, struct plan7_s *hmm)
{
  int k;			/* counter for nodes */
  int x;			/* counter for symbols */
  int ts;			/* counter for state transitions */
  
  /* score2ascii() uses static buffer, and
   * can't be called twice in the same line; be careful 
   */
  fprintf(fp, "N: %6s ", score2ascii(hmm->xsc[XTN][MOVE]));
  fprintf(fp, "%6s\n",   score2ascii(hmm->xsc[XTN][LOOP]));

  for (k = 1; k <= hmm->M; k++)
    {
				/* Line 1: k, match emissions */
      fprintf(fp, " %5d ", k);
      for (x = 0; x < Alphabet_size; x++) 
        fprintf(fp, "%6s ", score2ascii(hmm->msc[x][k]));
      fputs("\n", fp);
				/* Line 2: insert emissions */
      fprintf(fp, "       ");
      for (x = 0; x < Alphabet_size; x++) 
	fprintf(fp, "%6s ", (k < hmm->M) ? score2ascii(hmm->isc[x][k]) : "*");
      fputs("\n", fp);
				/* Line 3: transition probs; begin, end */
      fprintf(fp, "       ");
      for (ts = 0; ts < 7; ts++)
	fprintf(fp, "%6s ", (k < hmm->M) ? score2ascii(hmm->tsc[ts][k]) : "*"); 
      fprintf(fp, "%6s ", score2ascii(hmm->bsc[k]));
      fprintf(fp, "%6s ", score2ascii(hmm->esc[k]));
      fputs("\n", fp);
    }
  fprintf(fp, "E: %6s ", score2ascii(hmm->xsc[XTE][MOVE]));
  fprintf(fp, "%6s\n",   score2ascii(hmm->xsc[XTE][LOOP])); 

  fprintf(fp, "J: %6s ", score2ascii(hmm->xsc[XTJ][MOVE]));
  fprintf(fp, "%6s\n",   score2ascii(hmm->xsc[XTJ][LOOP])); 

  fprintf(fp, "C: %6s ", score2ascii(hmm->xsc[XTC][MOVE]));
  fprintf(fp, "%6s\n",   score2ascii(hmm->xsc[XTC][LOOP])); 

  fputs("//\n", fp);
}

static char *
score2ascii(int sc)
{
  static char buf[8];
  if (sc == -INFTY) return "*";
  sprintf(buf, "%6d", sc);
  return buf;
}
    
/* Function:  Plan7_DumpCounts()
 * Incept:    SRE, Sun May  8 08:34:55 2005 [St. Louis]
 *
 * Purpose:   Debugging: given an HMM that contains weighted counts,
 *            (not normalized yet, but may or may not contain
 *            Dirichlet pseudocounts), dump those counts out as a save
 *            file.  
 */
void
Plan7_DumpCounts(FILE *fp, struct plan7_s *hmm)
{
  int k;			/* counter for nodes */
  int x;			/* counter for symbols */
  int ts;			/* counter for state transitions */
  
  fprintf(fp, "tbm1,tbd1: %5.1f %5.1f\n", hmm->tbm1, hmm->tbd1);
  for (k = 1; k <= hmm->M; k++)
    {
				/* Line 1: k, match emissions */
      fprintf(fp, " %5d ", k);
      for (x = 0; x < Alphabet_size; x++) 
        fprintf(fp, "%5.1f ", hmm->mat[k][x]);
      fputs("\n", fp);
				/* Line 2: insert emissions */
      fprintf(fp, "       ");
      for (x = 0; x < Alphabet_size; x++) 
	fprintf(fp, "%5.1f ", (k < hmm->M) ? hmm->ins[k][x] : 0);
      fputs("\n", fp);
				/* Line 3: transition probs */
      fprintf(fp, "       ");
      for (ts = 0; ts < 7; ts++)
	fprintf(fp, "%5.1f ", (k < hmm->M) ? hmm->t[k][ts] : 0); 
      fputs("\n", fp);
    }
  fputs("//\n", fp);
}



/************************************************************
 * @LICENSE@
 ************************************************************/


