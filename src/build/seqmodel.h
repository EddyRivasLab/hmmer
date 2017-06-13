#ifndef p7SEQMODEL_INCLUDED
#define p7SEQMODEL_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"

#include "base/p7_hmm.h"

extern int p7_Seqmodel(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int M, char *name,
		       ESL_DMATRIX *P, float *f, double popen, double pextend,
		       P7_HMM **ret_hmm);

#endif /*p7SEQMODEL_INCLUDED*/
