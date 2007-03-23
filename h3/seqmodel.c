/* Creating profile HMMs from single sequences.
 * 
 * Contents:
 *   1. Routines in the exposed API.
 * 
 * SRE, Fri Mar 23 07:54:02 2007 [Janelia] [Decembrists, Picaresque]
 * SVN $Id$
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

/*****************************************************************
 * 1. Routines in the exposed API.
 *****************************************************************/


int
p7_Seqmodel(ESL_ALPHABET *abc, ESL_DSQ *dsq, int M, ESL_DMATRIX *P, 
	    float *f, double tmi, double tii, double tmd, double tdd,
	    P7_HMM **ret_hmm)
{
  int     status;
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[HMM created from a query sequence]";
  int     k;

  hmm = p7_hmm_Create(M, abc);
  if (hmm == NULL) { status = eslEMEM; goto ERROR; }
  
  for (k = 0; k <= M; k++)
    {
      /* Use P matrix to obtain match emissions. */
      if (k > 0) esl_vec_D2F(P->mx[(int) dsq[k]], abc->K, hmm->mat[k]);

      /* Set inserts to background for now. This will be improved. */
      esl_vec_FCopy(f, abc->K, hmm->ins[k]);

      hmm->t[k][p7_TMM] = 1.0 - tmi - tmd;
      hmm->t[k][p7_TMI] = tmi;
      hmm->t[k][p7_TMD] = tmd;
      hmm->t[k][p7_TIM] = 1.0 - tii;
      hmm->t[k][p7_TII] = tii;
      hmm->t[k][p7_TDM] = 1.0 - tdd;
      hmm->t[k][p7_TDD] = tdd;
    }

  /* Deal w/ special stuff at node 0, M, overwriting a little of what we
   * just did. 
   */
  hmm->t[M][p7_TMM] = 1.0 - tmi;
  hmm->t[M][p7_TMD] = 0.;
  hmm->t[M][p7_TDM] = 1.0;
  hmm->t[M][p7_TDD] = 0.;
  
  /* Add mandatory annotation
   */
  p7_hmm_SetName(hmm, "query-based-hmm");
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  hmm->nseq     = 0;
  p7_hmm_SetCtime(hmm);
  hmm->checksum = 0;

  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
