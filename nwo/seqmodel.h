#ifndef h4SEQMODEL_INCLUDED
#define h4SEQMODEL_INCLUDED

#include "h4_config.h"

#include "easel.h"
#include "esl_alphabet.h"

#include "h4_profile.h"

extern int h4_seqmodel(const ESL_DSQ *dsq, int L, const ESL_ALPHABET *abc, H4_PROFILE **ret_hmm, char *errbuf);

#endif //h4SEQMODEL_INCLUDED
