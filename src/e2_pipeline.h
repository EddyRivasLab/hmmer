/* e2_pipeline
 *
 *   
*/
#ifndef E2_PIPELINE_INCLUDED
#define E2_PIPELINE_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_getopts.h"	/* ESL_GETOPTS           */
#include "esl_random.h"

#include "e1_rate.h"
#include "e1_bg.h"
#include "e2_gmx.h"
#include "e2_profilesq.h"
#include "e2_trace.h"
#include "evohmmer.h"

#include "hmmer.h"

/*****************************************************************
 * 7. E2_PIPELINE: E2 pipeline
 *typedef struct e2_trace_s {
****************************************************************/
typedef struct e2_pipeline_s {
  E2_GMX       *gx1;
  E2_GMX       *gx2;

  char          errbuf[eslERRBUFSIZE];
} E2_PIPELINE;



/* e2_pipeline.c */
extern E2_PIPELINE *e2_pipeline_Create(ESL_GETOPTS *go, int L1_hint, int L2_hint, int M);
extern int          e2_pipeline_Reuse(E2_PIPELINE *pli);
extern void         e2_pipeline_Destroy(E2_PIPELINE *pli);
extern int          e2_Pipeline(ESL_RANDOMNESS *r, E2_PIPELINE *pli, const PSQ *sq1, const PSQ *sq2, const float *fres, 
				E1_RATE *R, P7_RATE *R7, E1_BG *bg, P7_BG *bg7,
				float time1, float time2, PSQ **ret_sqa, E2_TRACE **ret_tr, 
				float *ret_sc, float *ret_accsc, E2_ALI e2ali, 
				int mode, int decode, int do_viterbi, double tol, char *errbuf, int verbose);


#endif
