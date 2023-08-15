/* msatree - funtions to build a tree from an alignment
 *
 */
#ifndef E2_MSA_INCLUDED
#define E2_MSA_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_sq.h"
#include "esl_random.h"
#include "esl_tree.h"

#include "hmmer.h"

#include "e2.h"
#include "e1_bg.h"
#include "e1_rate.h"
#include "e2_pipeline.h"
#include "evohmmer.h"

struct e2_data {
  ESL_RANDOMNESS *r;
  PSQ            *sql;
  PSQ            *sqr;
  float          *frq;
  float           sc;
  float           accsc;
  double          timel;
  double          timer;
  E2_PIPELINE    *pli;
  E2_ALI          e2ali;
  E2_OPT          e2opt;
  int             mode;
  int             do_viterbi;
  E1_RATE        *R;
  P7_RATE        *R7;
  float           tsat;
  E1_BG          *bg;
  P7_BG          *bg7;
  ESL_MSA        *e2msa;
  E2_TRACE       *tr;
  int             it;
  float           tol;
  char           *errbuf;
  int             verbose;
};

extern int e2_msa(ESL_RANDOMNESS *r, E1_RATE *R, P7_RATE *R7, int n, ESL_SQ **seq, ESL_MSA *msa, float *msafrq, ESL_TREE *T, ESL_MSA **ret_omsa, float *ret_sc, E2_PIPELINE *pli, 
		  E1_BG *bg, P7_BG *bg7, E2_ALI e2ali, E2_OPT optimize, int mode, int do_viterbi, double tol, char *errbuf, int verbose);
extern int e2_Optimize(ESL_RANDOMNESS *r, E2_PIPELINE *pli, PSQ *sql, PSQ *sqr, float *msafrq, E1_RATE *R, P7_RATE *R7, E1_BG *bg, P7_BG *bg7,
		       double *ret_timel, double *ret_timer, PSQ **ret_sqa, E2_TRACE **ret_tr, float *ret_sc, E2_ALI e2ali, E2_OPT e2optimize, 
		       int mode, int do_viterbi, double tol, char *errbuf, int verbose);

#endif /*E2_MSA_INCLUDED*/


/************************************************************
 * @LICENSE@
 ************************************************************/
