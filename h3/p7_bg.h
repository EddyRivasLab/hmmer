/* The background model structure.
 * 
 * SRE, Fri Jan 12 13:26:30 2007 [Janelia]
 * SVN $Id$
 */
#ifdef  P7_BGH_INCLUDED
#define P7_BGH_INCLUDED

#include "p7_config.h"

#include "esl_alphabet.h"

/* P7_BG
 * The null hypothesis, sequence background model.
 */
typedef struct {
  ESL_ALPHABET *abc;		/* reference to alphabet in use       */
  
  float  p1;			/* null model's self-loop probability */
  float *f;			/* residue frequencies [0..K-1] */
} P7_BG;

extern P7_BG *p7_bg_Create(ESL_ALPHABET *abc);
extern void   p7_bg_Destroy(P7_BG *bg);

#endif /*P7_BGH_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
