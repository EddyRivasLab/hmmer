#ifndef p7TOPHITS_OUTPUT_INCLUDED
#define p7TOPHITS_OUTPUT_INCLUDED

#include <p7_config.h>

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_keyhash.h"
#include "esl_msa.h"

#include "base/p7_tophits.h"
#include "search/p7_pipeline.h"

extern int p7_tophits_ComputeNhmmerEvalues(P7_TOPHITS *th, double N, int W);
extern int p7_tophits_RemoveDuplicates(P7_TOPHITS *th, int using_bit_cutoffs);
extern int p7_tophits_Threshold(P7_TOPHITS *th, P7_PIPELINE *pli);
extern int p7_tophits_CompareRanking(P7_TOPHITS *th, ESL_KEYHASH *kh, int *opt_nnew);
extern int p7_tophits_Targets(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw);
extern int p7_tophits_Domains(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw);
extern int p7_tophits_Alignment(const P7_TOPHITS *th, const ESL_ALPHABET *abc, 
				ESL_SQ **inc_sqarr, P7_TRACE **inc_trarr, int inc_n, int optflags,
				ESL_MSA **ret_msa);
extern int p7_tophits_AliScores(FILE *ofp, char *qname, P7_TOPHITS *th );


#endif /*p7TOPHITS_OUTPUT_INCLUDED*/
