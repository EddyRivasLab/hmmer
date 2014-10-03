#ifndef p7PRIOR_INCLUDED
#define p7PRIOR_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dirichlet.h"

#include "base/p7_hmm.h"

typedef struct p7_prior_s {
  ESL_MIXDCHLET *tm;		/*  match transitions */
  ESL_MIXDCHLET *ti;		/* insert transitions */
  ESL_MIXDCHLET *td;		/* delete transitions */
  ESL_MIXDCHLET *em;		/*  match emissions   */
  ESL_MIXDCHLET *ei;		/* insert emissions   */
} P7_PRIOR;

extern P7_PRIOR  *p7_prior_CreateAmino(void);
extern P7_PRIOR  *p7_prior_CreateNucleic(void);
extern P7_PRIOR  *p7_prior_CreateLaplace(const ESL_ALPHABET *abc);
extern void       p7_prior_Destroy(P7_PRIOR *pri);

extern int        p7_ParameterEstimation(P7_HMM *hmm, const P7_PRIOR *pri);

#endif /*p7PRIOR_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

