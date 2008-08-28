/* Support for the background model, P7_BG
 * 
 * Contents:
 *     1. The P7_BG object: allocation, initialization, destruction.
 *     2. Simple null model scores
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
p7_bg_Dump(FILE *ofp, P7_BG *bg)
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
    if (bg->f != NULL) free(bg->f);
    free(bg);
  }
  return;
}


int
p7_bg_SetLength(P7_BG *bg, int L)
{
  bg->p1 = (float) L / (float) (L+1);
  return eslOK;
}



/*****************************************************************
 * 2. Simple null model scores
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


/* bg->mcomp has been set; create the two-state HMM */
static int
create_filter_hmm(P7_BG *bg)
{
  float L0 = 400.0;		/* mean length in state 0 */
  float L1 = 400.0;		/* mean length in state 1 */

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

  bg->fhmm->pi[0] = 0.99;
  bg->fhmm->pi[1] = 0.01;

  esl_hmm_Configure(bg->fhmm, bg->f);
  return eslOK;
}

int
p7_bg_SetFilterByHMM(P7_BG *bg, P7_HMM *hmm)
{
  ESL_RANDOMNESS  *rng  = esl_randomness_CreateTimeseeded();     /* source of randomness                    */
  ESL_SQ          *sq   = esl_sq_CreateDigital(bg->abc);
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

  create_filter_hmm(bg);

  esl_sq_Destroy(sq);
  esl_randomness_Destroy(rng);
  return eslOK;

 ERROR:
  esl_sq_Destroy(sq);
  esl_randomness_Destroy(rng);
  return status;
}


int
p7_bg_FilterScore(P7_BG *bg, ESL_DSQ *dsq, int L, float *ret_sc)
{
  ESL_HMX *hmx = esl_hmx_Create(L, bg->fhmm->M); /* optimization target: this can be a 2-row matrix, and it can be stored in <bg>. */
  float nullsc;
  int   i;
  
  esl_hmm_Forward(dsq, L, bg->fhmm, hmx, &nullsc);

  /* impose the length distribution */
  *ret_sc = nullsc + (float) L * logf(bg->p1) + logf(1.-bg->p1);
  esl_hmx_Destroy(hmx);
  return eslOK;
}
