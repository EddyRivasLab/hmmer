/* lsjfuncs.h
 * Declarations of external functions used in lsj_eweight.c
 * (Entropy-based sequence weighting)
 *
 * Steve Johnson
 * SVN $Id$
 */


#include "config.h"

#include "squid.h"
#include "msa.h"

#include "plan7.h"
#include "structs.h"


extern float Eweight(struct plan7_s *hmm,  struct p7prior_s *pri, 
		     float numb_seqs, float entwgt);
extern void ModelContent(float *ent1, float *ent2, int M);

/************************************************************
 * @LICENSE@
 ************************************************************/

