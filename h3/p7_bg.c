/* Support for the background model, P7_BG
 * 
 * Contents:
 *     1. The P7_BG object: allocation, initialization, destruction.
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
p7_bg_Create(ESL_ALPHABET *abc)
{
  P7_BG *bg = NULL;
  int    status;

  ESL_ALLOC(bg, sizeof(P7_BG));
  bg->f = NULL;

  ESL_ALLOC(bg->f, sizeof(float) * abc->K);
  if       (abc->type == eslAMINO && p7_AminoFrequencies(bg->f) != eslOK) goto ERROR;
  else     esl_vec_FSet(bg->f, abc->K, 1. / (float) abc->K);

  bg->p1  = 350./351.;
  bg->abc = abc;
  return bg;

 ERROR:
  p7_bg_Destroy(bg);
  return NULL;
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
