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
  char  *blastpamfile;            /* BLAST looks in aa/ subdirectory of BLASTMAT */
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

  blastpamfile = FileConcat("aa", pamfile);

  if ((fp = fopen(pamfile, "r")) == NULL &&
      (fp = EnvFileOpen(pamfile, "BLASTMAT")) == NULL &&
      (fp = EnvFileOpen(blastpamfile, "BLASTMAT")) == NULL)
    Die("Failed to open PAM scoring matrix file %s", pamfile);
  if (! ParsePAMFile(fp, &pam, &scale))
    Die("Failed to parse PAM scoring matrix file %s", pamfile);
  fclose(fp);
  free(blastpamfile);

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
  Die("%s is not in HMMER null model file format", rndfile);
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
  int q, x;
				/* default match mixture coefficients */
  static float defmq[9] = {
    0.178091, 0.056591, 0.0960191, 0.0781233, 0.0834977, 
    0.0904123, 0.114468, 0.0682132, 0.234585 };

				/* default match mixture Dirichlet components */
  static float defm[9][20] = {
    { 0.270671, 0.039848, 0.017576, 0.016415, 0.014268, 
      0.131916, 0.012391, 0.022599, 0.020358, 0.030727, 
      0.015315, 0.048298, 0.053803, 0.020662, 0.023612,
      0.216147, 0.147226, 0.065438, 0.003758, 0.009621 },
    { 0.021465, 0.010300, 0.011741, 0.010883, 0.385651, 
      0.016416, 0.076196, 0.035329, 0.013921, 0.093517, 
      0.022034, 0.028593, 0.013086, 0.023011, 0.018866, 
      0.029156, 0.018153, 0.036100, 0.071770, 0.419641 },
    { 0.561459, 0.045448, 0.438366, 0.764167, 0.087364,
      0.259114, 0.214940, 0.145928, 0.762204, 0.247320,
      0.118662, 0.441564, 0.174822, 0.530840, 0.465529, 
      0.583402, 0.445586, 0.227050, 0.029510, 0.121090 },
    { 0.070143, 0.011140, 0.019479, 0.094657, 0.013162, 
      0.048038, 0.077000, 0.032939, 0.576639, 0.072293, 
      0.028240, 0.080372, 0.037661, 0.185037, 0.506783, 
      0.073732, 0.071587, 0.042532, 0.011254, 0.028723 },
    { 0.041103, 0.014794, 0.005610, 0.010216, 0.153602, 
      0.007797, 0.007175, 0.299635, 0.010849, 0.999446, 
      0.210189, 0.006127, 0.013021, 0.019798, 0.014509, 
      0.012049, 0.035799, 0.180085, 0.012744, 0.026466 },
    { 0.115607, 0.037381, 0.012414, 0.018179, 0.051778, 
      0.017255, 0.004911, 0.796882, 0.017074, 0.285858, 
      0.075811, 0.014548, 0.015092, 0.011382, 0.012696, 
      0.027535, 0.088333, 0.944340, 0.004373, 0.016741 },
    { 0.093461, 0.004737, 0.387252, 0.347841, 0.010822, 
      0.105877, 0.049776, 0.014963, 0.094276, 0.027761, 
      0.010040, 0.187869, 0.050018, 0.110039, 0.038668, 
      0.119471, 0.065802, 0.025430, 0.003215, 0.018742 },
    { 0.452171, 0.114613, 0.062460, 0.115702, 0.284246,
      0.140204, 0.100358, 0.550230, 0.143995, 0.700649, 
      0.276580, 0.118569, 0.097470, 0.126673, 0.143634, 
      0.278983, 0.358482, 0.661750, 0.061533, 0.199373 },
    { 0.005193, 0.004039, 0.006722, 0.006121, 0.003468, 
      0.016931, 0.003647, 0.002184, 0.005019, 0.005990, 
      0.001473, 0.004158, 0.009055, 0.003630, 0.006583, 
      0.003172, 0.003690, 0.002967, 0.002772, 0.002686 },
  };

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
  
  /* Match emission priors are a mixture Dirichlet,
   * from Kimmen Sjolander (Blocks9)
   */
  pri->mnum  = 9;
  for (q = 0; q < pri->mnum; q++) 
    {
      pri->mq[q] = defmq[q];
      for (x = 0; x < 20; x++)
	pri->m[q][x] = defm[q][x];
    }

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
 * Purpose:  Set the default DNA prior. (for now, almost a Laplace)
 */
struct p7prior_s *
default_nucleic_prior(void)
{
  struct p7prior_s *pri;

  pri = P7AllocPrior();
  pri->strategy = PRI_DCHLET;

  /* The use of the Pfam-trained amino acid transition priors
   * here is TOTALLY bogus. But it works better than a straight
   * Laplace, esp. for Maxmodelmaker(). For example, a Laplace 
   * prior builds M=1 models for a single sequence GAATTC (at
   * one time an open "bug").
   */
  pri->tnum        = 1;
  pri->tq[0]       = 1.;
  pri->t[0][TMM]   = 0.7939;
  pri->t[0][TMI]   = 0.0278;
  pri->t[0][TMD]   = 0.0135;
  pri->t[0][TIM]   = 0.1551;
  pri->t[0][TII]   = 0.1331;
  pri->t[0][TDM]   = 0.9002;
  pri->t[0][TDD]   = 0.5630;
  
  pri->mnum  = 1;
  pri->mq[0] = 1.;
  FSet(pri->m[0], Alphabet_size, 1.);

  pri->inum  = 1;
  pri->iq[0] = 1.;
  FSet(pri->i[0], Alphabet_size, 1.);

  return pri;
}

