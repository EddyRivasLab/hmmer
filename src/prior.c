/* prior.c
 * SRE, Mon Nov 18 15:44:08 1996
 * 
 * Support for Dirichlet prior data structure, p7prior_s.
 */

#include "config.h"
#include "structs.h"
#include "funcs.h" 
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static struct p7prior_s *default_amino_prior(void);
static struct p7prior_s *default_nucleic_prior(void);

/* Function: P7AllocPrior(), P7FreePrior()
 * 
 * Purpose:  Allocation and free'ing of a prior structure.
 *           Very simple, but might get more complex someday.
 */
struct p7prior_s *
P7AllocPrior(void)
{ return (struct p7prior_s *) MallocOrDie (sizeof(struct p7prior_s)); }
void
P7FreePrior(struct p7prior_s *pri)
{ free(pri); }


/* Function: P7LaplacePrior()
 * 
 * Purpose:  Create a Laplace plus-one prior. (single component Dirichlets). 
 *           Global alphabet info is assumed to have been set already.
 *
 * Args:     (void)
 *
 * Return:   prior. Allocated here; call FreePrior() to free it.
 */ 
struct p7prior_s *
P7LaplacePrior(void)
{
  struct p7prior_s *pri;
  
  pri = P7AllocPrior();
  pri->strategy = PRI_DCHLET;

  pri->tnum     = 1;
  pri->tq[0]    = 1.;
  FSet(pri->t[0], 8, 1.); 
  
  pri->mnum  = 1;
  pri->mq[0] = 1.;
  FSet(pri->m[0], Alphabet_size, 1.);

  pri->inum  = 1;
  pri->iq[0] = 1.;
  FSet(pri->i[0], Alphabet_size, 1.);

  return pri;
}

/* Function: P7DefaultPrior()
 * 
 * Purpose:  Set up a somewhat more realistic single component
 *           Dirichlet prior than Laplace.
 */ 
struct p7prior_s *
P7DefaultPrior(void)
{
  switch (Alphabet_type) {
  case hmmAMINO:     return default_amino_prior();
  case hmmNUCLEIC:   return default_nucleic_prior();
  case hmmNOTSETYET: Die("Can't set prior; alphabet type not set yet");
  }
  /*NOTREACHED*/
  return NULL;
}

/* Function: P7ReadPrior()
 * 
 * Purpose:  Input a prior from disk file.
 */
struct p7prior_s *
P7ReadPrior(char *prifile) 
{
  FILE             *fp;
  struct p7prior_s *pri;
  char             *sptr;
  int               q, x;

  if ((fp = fopen(prifile, "r")) == NULL)
    Die("Failed to open HMMER prior file %s\n", prifile);
  pri = P7AllocPrior();

  /* First entry is the strategy: 
   * Only standard Dirichlet prior (simple or mixture) is supported in Plan7 so far
   */
  sptr = Getword(fp, sqdARG_STRING);
  s2upper(sptr);
  if      (strcmp(sptr, "DIRICHLET") == 0) pri->strategy = PRI_DCHLET;
  else Die("No such prior strategy %s; failed to parse file %s", sptr, prifile);

  /* Second entry is the alphabet type:
   * Amino or Nucleic
   */
  sptr = Getword(fp, sqdARG_STRING);
  s2upper(sptr);
  if (strcmp(sptr, "AMINO") == 0)
    { 
      if (Alphabet_type != hmmAMINO)
	Die("HMM and/or sequences are DNA/RNA; can't use protein prior %s", prifile);
    }
  else if (strcmp(sptr, "NUCLEIC") == 0)
    {
      if (Alphabet_type != hmmNUCLEIC)
	Die("HMM and/or sequences are protein; can't use DNA/RNA prior %s", prifile);
    }
  else 
    Die("Alphabet \"%s\" in prior file %s isn't valid.", sptr, prifile);

  /* State transition priors:
   * # of mixtures.
   * then for each mixture:
   *    prior P(q)
   *    Dirichlet terms for Tmm, Tmi, Tmd, Tim, Tii, Tid, Tdm, Tdi, Tdd
   */
  pri->tnum = atoi(Getword(fp, sqdARG_INT));
  if (pri->tnum < 0)
    Die("%d is bad; need at least one state transition mixture component", pri->tnum);
  if (pri->tnum > MAXDCHLET)
    Die("%d is bad, too many transition components (MAXDCHLET = %d)\n", MAXDCHLET);
  for (q = 0; q < pri->tnum; q++)
    {
      pri->tq[q]    = (float) atof(Getword(fp, sqdARG_FLOAT));
      for (x = 0; x < 7; x++) 
	pri->t[q][x] = (float) atof(Getword(fp, sqdARG_FLOAT));
    }

  /* Match emission priors:
   * # of mixtures.
   * then for each mixture:
   *    prior P(q)
   *    Dirichlet terms for Alphabet_size symbols in Alphabet
   */
  pri->mnum = atoi(Getword(fp, sqdARG_INT));
  if (pri->mnum < 0)
    Die("%d is bad; need at least one match emission mixture component", pri->mnum);
  if (pri->mnum > MAXDCHLET)
    Die("%d is bad; too many match components (MAXDCHLET = %d)\n", pri->mnum, MAXDCHLET);

  for (q = 0; q < pri->mnum; q++)
    {
      pri->mq[q] = (float) atof(Getword(fp, sqdARG_FLOAT));
      for (x = 0; x < Alphabet_size; x++) 
	pri->m[q][x] = (float) atof(Getword(fp, sqdARG_FLOAT));
    }
  
  /* Insert emission priors:
   * # of mixtures.
   * then for each mixture component:
   *    prior P(q)
   *    Dirichlet terms for Alphabet_size symbols in Alphabet
   */
  pri->inum = atoi(Getword(fp, sqdARG_INT));
  if (pri->inum < 0)
    Die("%d is bad; need at least one insert emission mixture component", pri->inum);
  if (pri->inum > MAXDCHLET)
    Die("%d is bad; too many insert components (MAXDCHLET = %d)\n", pri->inum,  MAXDCHLET);
  for (q = 0; q < pri->inum; q++)
    {
      pri->iq[q]  = (float) atof(Getword(fp, sqdARG_FLOAT));
      for (x = 0; x < Alphabet_size; x++) 
	pri->i[q][x] = (float) atof(Getword(fp, sqdARG_FLOAT));
    }

  fclose(fp);
  return pri;
}


/* Function: PAMPrior()
 * 
 * Purpose:  Produces an ad hoc "Dirichlet mixture" prior for
 *           match emissions, using a PAM matrix. 
 *           
 *           Side effect notice: PAMPrior() replaces the match
 *           emission section of an existing Dirichlet prior,
 *           which is /expected/ to be a simple one-component 
 *           kind of prior. The insert emissions /must/ be a
 *           one-component prior (because of details in how 
 *           PriorifyEmissionVector() is done). However, 
 *           the transitions /could/ be a mixture Dirichlet prior 
 *           without causing problems. In other words, the
 *           -p and -P options of hmmb can coexist, but there
 *           may be conflicts. PAMPrior() checks for these,
 *           so there's no serious problem, except that the
 *           error message from PAMPrior() might be confusing to
 *           a user. 
 */
void
PAMPrior(char *pamfile, struct p7prior_s *pri, float wt)
{
  FILE  *fp;
  int  **pam;
  float  scale;
  int    xi, xj;
  int    idx1, idx2;

  if (Alphabet_type != hmmAMINO)
    Die("PAM prior is only valid for protein sequences");
  if (pri->strategy != PRI_DCHLET)
    Die("PAM prior may only be applied over an existing Dirichlet prior");
  if (pri->inum != 1)
    Die("PAM prior requires that the insert emissions be a single Dirichlet");
  if (MAXDCHLET < 20)
    Die("Whoa, code is misconfigured; MAXDCHLET must be >= 20 for PAM prior");
  if ((fp = fopen(pamfile, "r")) == NULL &&
      (fp = EnvFileOpen(pamfile, "BLASTMAT")) == NULL)
    Die("Failed to open PAM scoring matrix file %s", pamfile);
  if (! ParsePAMFile(fp, &pam, &scale))
    Die("Failed to parse PAM scoring matrix file %s", pamfile);
  fclose(fp);

  pri->strategy = PRI_PAM;
  pri->mnum     = 20;
  
  /* Convert PAM entries back to conditional prob's P(xj | xi),
   * which we'll use as "pseudocounts" weighted by wt.
   */
  for (xi = 0; xi < Alphabet_size; xi++)
    for (xj = 0; xj < Alphabet_size; xj++)
      {
        idx1 = Alphabet[xi] - 'A';
        idx2 = Alphabet[xj] - 'A';
        pri->m[xi][xj] = aafq[xj] * exp((float) pam[idx1][idx2] * scale);
      }
  
  /* Normalize so that rows add up to wt.
   * i.e. Sum(xj) mat[xi][xj] = wt for every row xi
   */
  for (xi = 0; xi < Alphabet_size; xi++)
    {
      pri->mq[xi] = 1. / Alphabet_size;
      FNorm(pri->m[xi], Alphabet_size);
      FScale(pri->m[xi], Alphabet_size, wt);
    }

  Free2DArray(pam,27);
}


/* Function: P7DefaultNullModel()
 * 
 * Purpose:  Set up a default random sequence model, using
 *           global aafq[]'s for protein or 1/Alphabet_size for anything
 *           else. randomseq is alloc'ed in caller. Alphabet information
 *           must already be known.
 */
void
P7DefaultNullModel(float *null, float *ret_p1)
{
  int x;
  if (Alphabet_type == hmmAMINO) {
    for (x = 0; x < Alphabet_size; x++)
      null[x] = aafq[x];
    *ret_p1 = 350./351.;	/* rationale: approx avg protein length. */
  } else {
    for (x = 0; x < Alphabet_size; x++)
      null[x] = 1.0 / (float) Alphabet_size;
    *ret_p1 = 1000./1001.;	/* rationale: approx inter-Alu distance. */
  }
}

void
P7ReadNullModel(char *rndfile, float *null, float *ret_p1)
{
  FILE *fp;
  char *s;
  int   x;
  int   type = 0; 

  if ((fp = fopen(rndfile, "r")) == NULL)
    Die("Failed to open null model file %s\n", rndfile);
  if ((s = Getword(fp, sqdARG_STRING)) == NULL) goto FAILURE;
  s2upper(s);
  if      (strcmp(s, "NUCLEIC") == 0) type = hmmNUCLEIC;
  else if (strcmp(s, "AMINO")   == 0) type = hmmAMINO;
  else    goto FAILURE;
				/* check/set alphabet type */
  if (Alphabet_type == 0) 
    SetAlphabet(type);
  else if (Alphabet_type != type)
    Die("Alphabet type conflict; null model in %s is inappropriate\n", rndfile);
				/* parse the file */
  for (x = 0; x < Alphabet_size; x++) {
    if ((s = Getword(fp, sqdARG_FLOAT)) == NULL) goto FAILURE;
    null[x] = atof(s);
  }
  if ((s = Getword(fp, sqdARG_FLOAT)) == NULL) goto FAILURE;
  *ret_p1 = atof(s);

  fclose(fp);
  return;

FAILURE:
  fclose(fp);
  Die("%s is not in HMMER null model file format\n");
}


/* Function: P7PriorifyHMM()
 * 
 * Purpose:  Add pseudocounts to an HMM using Dirichlet priors,
 *           and renormalize the HMM.
 * 
 * Args:     hmm -- the HMM to add counts to (counts form)
 *           pri -- the Dirichlet prior to use
 *           
 * Return:   (void)
 *           HMM returns in probability form.
 */          
void
P7PriorifyHMM(struct plan7_s *hmm, struct p7prior_s *pri)
{
  int k;			/* counter for model position   */
  float d;			/* a denominator */

  /* Model-dependent transitions are handled simply; Laplace.
   */
  FSet(hmm->begin+2, hmm->M-1, 0.);     /* wipe internal BM entries */
  FSet(hmm->end+1, hmm->M-1, 0.);	/* wipe internal ME exits   */
  d = hmm->tbd1 + hmm->begin[1] + 2.;
  hmm->tbd1        = (hmm->tbd1 + 1.)/ d;
  hmm->begin[1]    = (hmm->begin[1] + 1.)/ d;
  hmm->end[hmm->M] = 1.0;

  /* Main model transitions and emissions
   */
  for (k = 1; k < hmm->M; k++)
    {
      P7PriorifyTransitionVector(hmm->t[k], pri);
      P7PriorifyEmissionVector(hmm->mat[k], pri, pri->mnum, pri->mq, pri->m, NULL);
      P7PriorifyEmissionVector(hmm->ins[k], pri, pri->inum, pri->iq, pri->i, NULL);
    }
  P7PriorifyEmissionVector(hmm->mat[hmm->M], pri, pri->mnum, pri->mq, pri->m, NULL);

  Plan7Renormalize(hmm);
}


/* Function: P7PriorifyEmissionVector()
 * 
 * Purpose:  Add prior pseudocounts to an observed 
 *           emission count vector and renormalize. 
 *
 *           Can return the posterior mixture probabilities
 *           P(q | counts) if ret_mix[MAXDCHLET] is passed.
 *           Else, pass NULL.  
 * 
 * Args:     vec     - the 4 or 20-long vector of counts to modify
 *           pri     - prior data structure
 *           num     - pri->mnum or pri->inum; # of mixtures
 *           eq      - pri->mq or pri->iq; prior mixture probabilities
 *           e       - pri->i or pri->m; Dirichlet components          
 *           ret_mix - filled with posterior mixture probabilities, or NULL
 *                   
 * Return:   (void)
 *           The counts in vec are changed and normalized to probabilities.
 */                  
void
P7PriorifyEmissionVector(float *vec, struct p7prior_s *pri, 
		       int num, float eq[MAXDCHLET], float e[MAXDCHLET][MAXABET],
		       float *ret_mix)
{
  int   x;                      /* counter over vec                     */
  int   q;                      /* counter over mixtures                */
  float mix[MAXDCHLET];         /* posterior distribution over mixtures */
  float totc;                   /* total counts                         */
  float tota;                   /* total alpha terms                    */
  float xi;                     /* X_i term, Sjolander eq. 41           */

  /* Calculate mix[], which is the posterior probability
   * P(q | n) of mixture component q given the count vector n
   * (side effect note: note that an insert vector in a PAM prior
   * is passed with num = 1, bypassing pam prior code; this means
   * that inserts cannot be mixture Dirichlets...)
   */
  mix[0] = 1.0;
  if (pri->strategy == PRI_DCHLET && num > 1) 
    {
      for (q = 0; q < num; q++) 
	{
	  mix[q] =  eq[q] > 0.0 ? log(eq[q]) : -999.;
	  mix[q] += Logp_cvec(vec, Alphabet_size, e[q]);
	}
      LogNorm(mix, num);      /* now mix[q] is P(component_q | n) */
    }
  else if (pri->strategy == PRI_PAM && num > 1) 
    {		/* pam prior uses aa frequencies as `P(q|n)' */
      for (q = 0; q < Alphabet_size; q++) 
	mix[q] = vec[q];
      FNorm(mix, Alphabet_size);
    }

  /* Convert the counts to probabilities, following Sjolander (1996) 
   */
  totc = FSum(vec, Alphabet_size);
  for (x = 0; x < Alphabet_size; x++) {
    xi = 0.0;
    for (q = 0; q < num; q++) {
      tota = FSum(e[q], Alphabet_size);
      xi += mix[q] * (vec[x] + e[q][x]) / (totc + tota);
    }
    vec[x] = xi;
  }
  FNorm(vec, Alphabet_size);

  if (ret_mix != NULL)
    for (q = 0; q < num; q++)
      ret_mix[q] = mix[q];
}



/* Function: P7PriorifyTransitionVector()
 * 
 * Purpose:  Add prior pseudocounts to transition vector,
 *           which contains three different probability vectors
 *           for m, d, and i. 
 *           
 * Args:     t     - state transitions, counts: 3 for M, 2 for I, 2 for D.   
 *           prior - Dirichlet prior information
 *           
 * Return:   (void)
 *           t is changed, and renormalized -- comes back as
 *           probability vectors.
 */          
void
P7PriorifyTransitionVector(float *t, struct p7prior_s *prior)
{
  int   ts;
  int   q;
  float mix[MAXDCHLET];
  float totm, totd, toti;       /* total counts in three transition vecs */
  float xi;                     /* Sjolander's X_i term */

  mix[0] = 1.0;			/* default is simple one component */
  if ((prior->strategy == PRI_DCHLET || prior->strategy == PRI_PAM) && prior->mnum > 1)
    {
      for (q = 0; q < prior->tnum; q++)
        {
          mix[q] =  prior->tq[q] > 0.0 ? log(prior->tq[q]) : -999.;
          mix[q] += Logp_cvec(t,   3, prior->t[q]);   /* 3 match  */
          mix[q] += Logp_cvec(t+3, 2, prior->t[q]+3); /* 2 insert */
	  mix[q] += Logp_cvec(t+5, 2, prior->t[q]+5); /* 2 delete */
        }
      LogNorm(mix, prior->tnum); /* mix[q] is now P(q | counts) */
    }
				/* precalc some denominators */
  totm = FSum(t,3);		
  toti = t[TIM] + t[TII];
  totd = t[TDM] + t[TDD];

  for (ts = 0; ts < 7; ts++)  
    {
      xi = 0.0;
      for (q = 0; q < prior->tnum; q++)
        {
	  switch (ts) {
	  case TMM: case TMI: case TMD: 
	    xi += mix[q] * (t[ts] + prior->t[q][ts]) / 
	      (totm + FSum(prior->t[q], 3)); 
	    break;
	  case TIM: case TII: 
	    xi += mix[q] * (t[ts] + prior->t[q][ts]) / 
	      (toti + prior->t[q][TIM] + prior->t[q][TII]);
	    break;
	  case TDM: case TDD: 
	    xi += mix[q] * (t[ts] + prior->t[q][ts]) / 
	      (totd + prior->t[q][TDM] + prior->t[q][TDD]);
	    break;
	  }
        }
      t[ts] = xi;
    }
  FNorm(t,   3);		/* match  */
  FNorm(t+3, 2);		/* insert */
  FNorm(t+5, 2);		/* delete */
}


/* Function: default_amino_prior()
 * 
 * Purpose:  Set the default protein prior.
 */
struct p7prior_s *
default_amino_prior(void)
{
  struct p7prior_s *pri;
  int x;

  pri = P7AllocPrior();
  pri->strategy = PRI_DCHLET;

  /* Transition priors are subjective, but borrowed from GJM's estimations
   * on Pfam
   */
  pri->tnum     = 1;
  pri->tq[0]    = 1.0;
  pri->t[0][TMM]   = 0.7939;
  pri->t[0][TMI]   = 0.0278;
  pri->t[0][TMD]   = 0.0135;
  pri->t[0][TIM]   = 0.1551;
  pri->t[0][TII]   = 0.1331;
  pri->t[0][TDM]   = 0.9002;
  pri->t[0][TDD]   = 0.5630;
  
  /* Subjective match emission priors: Swissprot 34 frequencies * 3.
   */
  pri->mnum  = 1;
  pri->mq[0] = 1.0;
  for (x = 0; x < 20; x++)
    pri->m[0][x] = aafq[x] * 3.;
  
  /* These insert emission priors are subjective. Observed frequencies
   * were obtained from PFAM 1.0, 10 Nov 96; 
   *      see ~/projects/plan7/InsertStatistics.
   * Inserts are slightly biased towards polar residues and away from
   * hydrophobic residues.
   */
  pri->inum  = 1;
  pri->iq[0] = 1.;
  pri->i[0][0]  = 681.;         /* A */
  pri->i[0][1]  = 120.;         /* C */
  pri->i[0][2]  = 623.;         /* D */
  pri->i[0][3]  = 651.;         /* E */
  pri->i[0][4]  = 313.;         /* F */
  pri->i[0][5]  = 902.;         /* G */
  pri->i[0][6]  = 241.;         /* H */
  pri->i[0][7]  = 371.;         /* I */
  pri->i[0][8]  = 687.;         /* K */
  pri->i[0][9]  = 676.;         /* L */
  pri->i[0][10] = 143.;         /* M */
  pri->i[0][11] = 548.;         /* N */
  pri->i[0][12] = 647.;         /* P */
  pri->i[0][13] = 415.;         /* Q */
  pri->i[0][14] = 551.;         /* R */
  pri->i[0][15] = 926.;         /* S */
  pri->i[0][16] = 623.;         /* T */
  pri->i[0][17] = 505.;         /* V */
  pri->i[0][18] = 102.;         /* W */
  pri->i[0][19] = 269.;         /* Y */

  return pri;
}


/* Function: default_nucleic_prior()
 * 
 * Purpose:  Set the default DNA prior. (for now, a Laplace)
 */
struct p7prior_s *
default_nucleic_prior(void)
{
  struct p7prior_s *pri;

  pri = P7AllocPrior();
  pri->strategy = PRI_DCHLET;

  pri->tnum     = 1;
  pri->tq[0]    = 1.;
  FSet(pri->t[0], 7, 1.); 
  
  pri->mnum  = 1;
  pri->mq[0] = 1.;
  FSet(pri->m[0], Alphabet_size, 1.);

  pri->inum  = 1;
  pri->iq[0] = 1.;
  FSet(pri->i[0], Alphabet_size, 1.);

  return pri;
}

