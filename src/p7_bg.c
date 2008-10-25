/* P7_BG: the null (background) model
 * 
 * Contents:
 *     1. P7_BG object: allocation, initialization, destruction.
 *     2. Standard iid null model ("null1")
 *     3. Filter null model. 
 * 
 * SRE, Fri Jan 12 13:31:26 2007 [Janelia] [Ravel, Bolero]
 * SVN $Id$
 */

#include "p7_config.h"		/* must be included first */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_hmm.h"
#include "esl_random.h"		/* needed for now in SetFilter() routines; will disappear */
#include "esl_vectorops.h"

#include "hmmer.h"


/*****************************************************************
 * 1. The P7_BG object: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  p7_bg_Create()
 * Incept:    SRE, Fri Jan 12 13:32:51 2007 [Janelia]
 *
 * Purpose:   Allocate a <P7_BG> object for digital alphabet <abc>,
 *            initializes it to appropriate default values, and
 *            returns a pointer to it.
 *            
 *            For protein models, default iid background frequencies
 *            are set (by <p7_AminoFrequencies()>) to average
 *            SwissProt residue composition. For DNA, RNA and other
 *            alphabets, default frequencies are set to a uniform
 *            distribution.
 *            
 *            The model composition <bg->mcomp[]> is not initialized
 *            here; neither is the filter null model <bg->fhmm>.  To
 *            use the filter null model, caller will want to
 *            initialize these fields by calling
 *            <p7_bg_SetFilterByHMM()>.
 *
 * Throws:    <NULL> on allocation failure.
 *
 * Xref:      STL11/125.
 */
P7_BG *
p7_bg_Create(const ESL_ALPHABET *abc)
{
  P7_BG *bg = NULL;
  int    status;

  ESL_ALLOC(bg, sizeof(P7_BG));
  bg->f     = NULL;
  bg->mcomp = NULL;
  bg->fhmm  = NULL;

  ESL_ALLOC(bg->f,     sizeof(float) * abc->K);
  ESL_ALLOC(bg->mcomp, sizeof(float) * abc->K);
  if ((bg->fhmm = esl_hmm_Create(abc, 2)) == NULL) goto ERROR;

  if       (abc->type == eslAMINO)
    {
      if (p7_AminoFrequencies(bg->f) != eslOK) goto ERROR;
    }
  else
    esl_vec_FSet(bg->f, abc->K, 1. / (float) abc->K);

  bg->p1  = 350./351.;
  bg->abc = abc;
  return bg;

 ERROR:
  p7_bg_Destroy(bg);
  return NULL;
}


/* Function:  p7_bg_CreateUniform()
 * Synopsis:  Creates background model with uniform freqs.
 * Incept:    SRE, Sat Jun 30 10:25:27 2007 [Janelia]
 *
 * Purpose:   Creates a background model for alphabet <abc>
 *            with uniform residue frequencies.
 */
P7_BG *
p7_bg_CreateUniform(const ESL_ALPHABET *abc)
{
  P7_BG *bg = NULL;
  int    status;

  ESL_ALLOC(bg, sizeof(P7_BG));
  bg->f     = NULL;
  bg->mcomp = NULL;
  bg->fhmm  = NULL;

  ESL_ALLOC(bg->f,     sizeof(float) * abc->K);
  ESL_ALLOC(bg->mcomp, sizeof(float) * abc->K);
  if ((bg->fhmm = esl_hmm_Create(abc, 2)) == NULL) goto ERROR;

  esl_vec_FSet(bg->f, abc->K, 1. / (float) abc->K);
  bg->p1  = 350./351.;
  bg->abc = (ESL_ALPHABET *) abc; /* safe: we're just keeping a reference */
  return bg;

 ERROR:
  p7_bg_Destroy(bg);
  return NULL;
}


/* Function:  p7_bg_Dump()
 * Synopsis:  Outputs <P7_BG> object as text, for diagnostics.
 * Incept:    SRE, Fri May 25 08:07:11 2007 [Janelia]
 *
 * Purpose:   Given a null model <bg>, dump it as text to stream <fp>.
 */
int
p7_bg_Dump(FILE *ofp, const P7_BG *bg)
{
  esl_vec_FDump(ofp, bg->f, bg->abc->K, bg->abc->sym);
  return eslOK;
}



/* Function:  p7_bg_Destroy()
 * Incept:    SRE, Fri Jan 12 14:04:30 2007 [Janelia]
 *
 * Purpose:   Frees a <P7_BG> object.
 *
 * Returns:   (void)
 *
 * Xref:      STL11/125.
 */
void
p7_bg_Destroy(P7_BG *bg)
{
  if (bg != NULL) {
    if (bg->f     != NULL) free(bg->f);
    if (bg->mcomp != NULL) free(bg->mcomp);
    if (bg->fhmm  != NULL) esl_hmm_Destroy(bg->fhmm);
    free(bg);
  }
  return;
}


/* Function:  p7_bg_SetLength()
 * Synopsis:  Set the null model length distribution.
 * Incept:    SRE, Thu Aug 28 08:44:22 2008 [Janelia]
 *
 * Purpose:   Sets the geometric null model length 
 *            distribution in <bg> to a mean of <L> residues.
 */
int
p7_bg_SetLength(P7_BG *bg, int L)
{
  bg->p1 = (float) L / (float) (L+1);
  
  bg->fhmm->t[0][0] = bg->p1;
  bg->fhmm->t[0][1] = 1.0f - bg->p1;

  return eslOK;
}




/*****************************************************************
 * 2. Standard iid null model ("null1")
 *****************************************************************/

/* Function:  p7_bg_NullOne()
 * Incept:    SRE, Mon Apr 23 08:13:26 2007 [Janelia]
 *
 * Purpose:   Calculate the null1 lod score, for sequence <dsq>
 *            of length <L> "aligned" to the base null model <bg>. 
 * 
 * Note:      Because the residue composition in null1 <bg> is the
 *            same as the background used to calculate residue
 *            scores in profiles and null models, all we have to
 *            do here is score null model transitions.
 */
int
p7_bg_NullOne(const P7_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc)
{
  *ret_sc = (float) L * log(bg->p1) + log(1.-bg->p1);
  return eslOK;
}





/*****************************************************************
 * 3. Filter null model
 *****************************************************************/

static int create_filter_hmm(P7_BG *bg, int M);
int
p7_bg_SetFilterByHMM(P7_BG *bg, P7_HMM *hmm)
{
  float *mocc = NULL;
  float *iocc = NULL;
  int    status;
  int    k;
  float  norm = 0.0;

  ESL_ALLOC(mocc, sizeof(float) * (hmm->M+1));
  ESL_ALLOC(iocc, sizeof(float) * (hmm->M+1));

  p7_hmm_CalculateOccupancy(hmm, mocc, iocc);

  esl_vec_FSet(bg->mcomp, bg->abc->K, 0.0);
  esl_vec_FAddScaled(bg->mcomp, hmm->ins[0], iocc[0], hmm->abc->K);
  for (k = 1; k <= hmm->M; k++)
    {
      esl_vec_FAddScaled(bg->mcomp, hmm->mat[k], mocc[k], hmm->abc->K);
      esl_vec_FAddScaled(bg->mcomp, hmm->ins[k], iocc[k], hmm->abc->K);
    }
  norm  = esl_vec_FSum(mocc, hmm->M+1);
  norm += esl_vec_FSum(iocc, hmm->M+1);
  esl_vec_FScale(bg->mcomp, hmm->abc->K, 1.0 / norm);

  create_filter_hmm(bg, hmm->M);
  free(mocc);
  free(iocc);
  return eslOK;

 ERROR:
  if (mocc != NULL) free(mocc);
  if (iocc != NULL) free(iocc);
  return status;
}


#if 0
/* Function:  p7_bg_SetFilterByHMM()
 * Synopsis:  Parameterizes the filter null model from an HMM.
 * Incept:    SRE, Thu Aug 28 08:48:35 2008 [Janelia]
 *
 * Purpose:   Given a particular model <hmm>, parameterize a two-state
 *            model-dependent filter null model in <bg> in which one
 *            state is the usual iid background model, the second
 *            state emits the average residue composition of sequences
 *            emitted by <hmm>, and the transition probabilities are
 *            set to arbitrary hardcoded values that happen to work
 *            empirically.
 */
int
p7_bg_SetFilterByHMM(P7_BG *bg, P7_HMM *hmm)
{
  ESL_RANDOMNESS  *rng  = esl_randomness_CreateTimeseeded();     /* source of randomness     FIXME: see below  */
  ESL_SQ          *sq   = esl_sq_CreateDigital(bg->abc);         /* FIXME: don't need to do this by sampling; can do it analytically */
  int              nseq = 1000;
  int              pos;
  int              status;

  esl_vec_FSet(bg->mcomp, bg->abc->K, 1.0); /* a +1 prior */

  while (nseq--)
    {
      if ((status = p7_CoreEmit(rng, hmm, sq, NULL)) != eslOK) goto ERROR;

      for (pos = 1; pos <= sq->n; pos++)
	bg->mcomp[sq->dsq[pos]] += 1.0;

      esl_sq_Reuse(sq);
    }
  
  esl_vec_FNorm(bg->mcomp, bg->abc->K);

  create_filter_hmm(bg, hmm->M);

  esl_sq_Destroy(sq);
  esl_randomness_Destroy(rng);
  return eslOK;

 ERROR:
  esl_sq_Destroy(sq);
  esl_randomness_Destroy(rng);
  return status;
}
#endif


/* Function:  p7_bg_FilterScore()
 * Synopsis:  Calculates the filter null model score.
 * Incept:    SRE, Thu Aug 28 08:52:53 2008 [Janelia]
 *
 * Purpose:   Calculates the filter null model <bg> score for sequence
 *            <dsq> of length <L>, and return it in 
 *            <*ret_sc>.
 *            
 *            The score is calculated as an HMM Forward score using
 *            the two-state filter null model. It is a log-odds ratio,
 *            relative to the iid background frequencies, in nats:
 *            same as main model Forward scores.
 *
 *            The filter null model has no length distribution of its
 *            own; the same geometric length distribution (controlled
 *            by <bg->p1>) that the null1 model uses is imposed.
 */
int
p7_bg_FilterScore(P7_BG *bg, ESL_DSQ *dsq, int L, float *ret_sc)
{
  ESL_HMX *hmx = esl_hmx_Create(L, bg->fhmm->M); /* optimization target: this can be a 2-row matrix, and it can be stored in <bg>. */
  float nullsc;		                  	 /* (or it could be passed in as an arg, but for sure it shouldn't be alloc'ed here */
  
  esl_hmm_Forward(dsq, L, bg->fhmm, hmx, &nullsc);

  /* impose the length distribution */
  *ret_sc = nullsc + (float) L * logf(bg->p1) + logf(1.-bg->p1);
  esl_hmx_Destroy(hmx);
  return eslOK;
}


/* create_filter_hmm()
 * assumes bg->mcomp has just been set; 
 * now create the two-state HMM;
 * parameters are arbitrary hardwired guesses
 */
static int
create_filter_hmm(P7_BG *bg, int M)
{
  float L0 = 400.0;		/* mean length in state 0 */
  float L1 = (float) M / 8.0; 	/* mean length in state 1 */

  /* State 0 is the normal iid model. */
  bg->fhmm->t[0][0] =   L0 / (L0+1.0f);
  bg->fhmm->t[0][1] = 1.0f / (L0+1.0f);
  bg->fhmm->t[0][2] = 1.0f;          	/* 1.0 transition to E means we'll set length distribution externally. */
  esl_vec_FCopy(bg->f, bg->abc->K, bg->fhmm->e[0]);

  /* State 1 is the potentially biased model composition. */
  bg->fhmm->t[1][0] = 1.0f / (L1+1.0f);
  bg->fhmm->t[1][1] =   L1 / (L1+1.0f);
  bg->fhmm->t[1][2] = 1.0f;         	/* 1.0 transition to E means we'll set length distribution externally. */
  esl_vec_FCopy(bg->mcomp, bg->abc->K, bg->fhmm->e[1]);

  bg->fhmm->pi[0] = 0.999;
  bg->fhmm->pi[1] = 0.001;

  esl_hmm_Configure(bg->fhmm, bg->f);
  return eslOK;
}


/*****************************************************************
 * x. Benchmark
 *****************************************************************/
#ifdef p7BG_BENCHMARK
/*
   gcc -O2 -Wall -msse2 -std=gnu99 -o p7_bg_benchmark -I. -L. -I../easel -L../easel -Dp7BG_BENCHMARK p7_bg.c -lhmmer -leasel -lm
   ./p7_bg_benchmark <hmmfile>
 */ 
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",      0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0",     NULL,      NULL,    NULL, "length of random target seqs",              0 },
  { "-N",        eslARG_INT,    "100", NULL, "n>0",     NULL,      NULL,    NULL, "number of random target seqs",              0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark timing for calculating null model scores";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  int             i;
 
  /* Read one HMM from <hmmfile> */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  bg = p7_bg_Create(abc);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    p7_bg_SetFilterByHMM(bg, hmm);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7BG_BENCHMARK*/



/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7BG_EXAMPLE
/*
   gcc -O2 -Wall -msse2 -std=gnu99 -o p7_bg_example -I. -L. -I../easel -L../easel -Dp7BG_EXAMPLE p7_bg.c -lhmmer -leasel -lm
   ./p7_bg_example <hmmfile> <seqfile>
 */ 
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",      0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of calculating null model scores";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  ESL_SQ         *sq      = NULL;
  float           nullsc, filtersc, H;
  int             status;
 
  /* Read one HMM from <hmmfile> */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  /* Open <seqfile> for reading */
  status = esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  sq = esl_sq_CreateDigital(abc);
  bg = p7_bg_Create(abc);
  p7_bg_SetFilterByHMM(bg, hmm);

  H = esl_vec_FEntropy(bg->f, bg->abc->K);
  printf("bg iid H = %.4f\n", H);

  H = esl_vec_FEntropy(bg->mcomp, bg->abc->K);
  printf("modelcomp H = %.4f\n", H);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_bg_SetLength(bg, sq->n);

      p7_bg_NullOne    (bg, sq->dsq, sq->n, &nullsc);
      p7_bg_FilterScore(bg, sq->dsq, sq->n, &filtersc);

      printf("%-20s %5d %8.5f %8.5f %8.5f\n", sq->name, (int) sq->n, nullsc, filtersc, filtersc-nullsc);

      esl_sq_Reuse(sq);
    }
  if (status != eslEOF) 
    esl_fatal("Parse failed, line %ld, file %s:\n%s", (long) sqfp->linenumber, sqfp->filename, sqfp->errbuf);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7BG_EXAMPLE*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/


  
