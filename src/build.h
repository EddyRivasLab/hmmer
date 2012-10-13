#ifndef P7_BUILD_INCLUDED
#define P7_BUILD_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_msa.h"

#include "base/p7_hmm.h"
#include "base/p7_trace.h"

#include "p7_builder.h"

extern int p7_Handmodelmaker(ESL_MSA *msa,                P7_BUILDER *bld, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);
extern int p7_Fastmodelmaker(ESL_MSA *msa, float symfrac, P7_BUILDER *bld, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);

#endif /*P7_BUILD_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
