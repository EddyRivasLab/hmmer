#ifndef h4PRIOR_INCLUDED
#define h4PRIOR_INCLUDED

#include "h4_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dirichlet.h"

typedef struct {
  ESL_MIXDCHLET *tm;		/*  match transitions */
  ESL_MIXDCHLET *ti;		/* insert transitions */
  ESL_MIXDCHLET *td;		/* delete transitions */
  ESL_MIXDCHLET *em;		/*  match emissions   */
} H4_PRIOR;

extern H4_PRIOR *h4_prior_CreateAmino(void);
extern H4_PRIOR *h4_prior_CreateNucleic(void);
extern H4_PRIOR *h4_prior_CreateLaplace(const ESL_ALPHABET *abc);
extern void      h4_prior_Destroy(H4_PRIOR *pri);


#endif // h4PRIOR_INCLUDED
