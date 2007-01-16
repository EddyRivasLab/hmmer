/* The all-encompassing include file for HMMER.
 * 
 * SRE, Wed Jan  3 13:46:42 2007 [Janelia] [Philip Glass, The Fog of War]
 * SVN $Id$
 */
#ifndef P7_HMMERH_INCLUDED
#define P7_HMMERH_INCLUDED

/* Structures
 */
#include "easel.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_sqio.h"

#include "p7_hmm.h"
#include "p7_hmmfile.h"
#include "p7_trace.h"

/* build.c
 */
extern int p7_Handmodelmaker(ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);
extern int p7_Fastmodelmaker(ESL_MSA *msa, float symfrac, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);

/* emit.c
 */
extern int p7_CoreEmit(ESL_RANDOMNESS *r, P7_HMM *hmm, ESL_SQ *sq, P7_TRACE *tr);

/* errors.c
 */
extern void p7_Die(char *format, ...);

/* hmmer.c
 */
extern int   p7_Prob2Score(float p, float null);
extern int   p7_LL2Score(float ll, float null);
extern float p7_Score2Prob(int sc, float null);
extern int   p7_AminoFrequencies(float *f);



#endif /*P7_HMMERH_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
