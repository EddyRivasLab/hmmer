/*  e2_tree 
 *
 * ER, Wed Mar  5 21:34:31 EST 2014 [Janelia] 
 * SVN $Id:$
 */

#include "p7_config.h"

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include <math.h>
#include <float.h>
	
#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "e1_bg.h"
#include "e1_rate.h"
#include "e2_msa.h"
#include "e2_pipeline.h"
#include "e2_trace.h"
#include "e2tracealign.h"
#include "msatree.h"
#include "evohmmer.h"

int
e2_tree_UPGMA(ESL_TREE **ret_T, int n, ESL_SQ **seq, ESL_MSA *msa, float *frq, ESL_RANDOMNESS *r, E2_PIPELINE *pli,
	      E1_RATE *R, P7_RATE *R7, E1_BG *bg, P7_BG *bg7, E2_ALI e2ali, 
	      int mode, int do_viterbi, float fixtime, float tinit, double tol, char *errbuf, int verbose)
{
  ESL_TREE     *T    = NULL;
  ESL_DMATRIX  *D    = NULL;
  PSQ         **psq  = NULL;
  ESL_SQ       *sq   = NULL;
  double        dis;
  float         sc;
  double        timel, timer;
  int           N = (msa)? msa->nseq : n;
  int           i, j;
  int           status;

  
  ESL_ALLOC(psq, sizeof(PSQ *) * N);
  for (i = 0; i < N; i++) {

    if (msa) status = esl_sq_FetchFromMSA(msa, i, &sq); /* extract the seqs from the msa */
    else sq = seq[i];
    
    switch(e2ali) {
    case E2:
      psq[i] = psq_CreateFrom(sq->name, sq->desc, sq->acc, sq->abc, sq->dsq, sq->n);	
      break;
    case E2HMMER:
      psq[i] = psq_CreateFrom(sq->name, sq->desc, sq->acc, sq->abc, sq->dsq, sq->n);	
      break;
    case E2F:
      if (msa)
	psq[i] = psq_CreateFrom(sq->name, sq->desc, sq->acc, sq->abc, msa->ax[i], msa->alen);
      else goto ERROR;
      break;
    case E2FHMMER:
      psq[i] = psq_CreateFrom(sq->name, sq->desc, sq->acc, sq->abc, sq->dsq, sq->n);
      break;
    default:
      ESL_XFAIL(eslFAIL, errbuf, "unknown alicase\n"); goto ERROR;
    }	
    
    if (msa) { esl_sq_Destroy(sq); sq = NULL; }
  }	
  
  /* distances */
  D = esl_dmatrix_Create(N, N);
  esl_dmatrix_Set(D,    0.0);

  for (i = 0; i < N; i++) 
    for (j = i+1; j < N; j++)
      {
	if (fixtime > 0.0) {
	  timel = timer = fixtime;
	}
	else {
	  timel = timer = (tinit > 0)? tinit : 0.9;
	  status = e2_Optimize(r, pli, psq[i], psq[j], frq, R, R7, bg, bg7, &timel, &timer, NULL, NULL, &sc, e2ali, OPTTIME,
			       mode, do_viterbi, tol, errbuf, verbose);
	  }
	dis = (timel+timer);
      
	D->mx[i][j] = D->mx[j][i] = dis;
      }
  if (verbose) esl_dmatrix_Dump(stdout, D, NULL, NULL);
  if (esl_tree_UPGMA(D, &T) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to create tree\n");
  if (verbose) Tree_Dump(stdout, T, "the guide Tree");
  
#if 0
  /* root the Tree */
  if (Tree_InterLeafMaxDistRooted(T, NULL, errbuf, verbose) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to root the tree\n"); 
  if (verbose) Tree_Dump(stdout, T, "rooted guide Tree");
#endif

  esl_dmatrix_Destroy(D); 
  if (psq) {
    for (i = 0; i < N; i++) 
      psq_Destroy(psq[i]);
    free(psq);
  }

  *ret_T     = T;
  return eslOK;

 ERROR:
  if (T) esl_tree_Destroy(T);
  ret_T    = NULL;
  if (D) esl_dmatrix_Destroy(D);
  if (sq) esl_sq_Destroy(sq);
  if (psq) {
    for (i = 0; i < N; i++) psq_Destroy(psq[i]);
    free(psq);
  }
  return status;
}
