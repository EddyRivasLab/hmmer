#ifndef p7HMMPGMD2MSA_INCLUDED
#define p7HMMPGMD2MSA_INCLUDED

#include "p7_config.h"

#include "esl_msa.h"

#include "base/p7_hmm.h"

extern int hmmpgmd2msa(void *data, P7_HMM *hmm, ESL_SQ *qsq,  int *incl, int incl_size, int *excl, int excl_size, ESL_MSA **ret_msa);

#endif /*p7HMMPGMD2MSA_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
