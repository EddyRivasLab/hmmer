/* The all-encompassing include file for HMMER.
 * 
 * SRE, Wed Jan  3 13:46:42 2007 [Janelia] [Philip Glass, The Fog of War]
 * SVN $Id$
 */
#ifndef P7_HMMERH_INCLUDED
#define P7_HMMERH_INCLUDED

/* Structures
 */
#include "p7_hmm.h"
#include "p7_trace.h"

/* build.c
 */
extern int p7_Handmodelmaker(ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);
extern int p7_Fastmodelmaker(ESL_MSA *msa, float symfrac, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);



#endif /*P7_HMMERH_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
