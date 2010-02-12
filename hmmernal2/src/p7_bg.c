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
  bg->f = NULL;

  ESL_ALLOC(bg->f, sizeof(float) * abc->K);
  if       (abc->type == eslAMINO)
    {
      if (p7_AminoFrequencies(bg->f) != eslOK) goto ERROR;
    }
  else
    esl_vec_FSet(bg->f, abc->K, 1. / (float) abc->K);

  bg->p1  = 350./351.;
  bg->abc = (ESL_ALPHABET *) abc; /* safe: we're just keeping a reference */
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
  bg->f = NULL;
  ESL_ALLOC(bg->f, sizeof(float) * abc->K);

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

