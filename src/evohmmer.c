/* evohmmer - funtions to evolve parameters of an HMM
 * 
 * Contents:
 *   1. Miscellaneous functions for evoH3
 *   2. Unit tests
 *   3. Test driver
 *   4. License and copyright 
 *
 * ER, Tue Sep 27 13:19:30 2011 [Janelia] 
 * SVN $Id:$
 */

#include "p7_config.h"

#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_rootfinder.h"
#include "esl_scorematrix.h"
#include "esl_vectorops.h"
#include "hmmer.h"

#include "e2.h"
#include "e1_rate.h"
#include "e1_model.h"
#include "evohmmer.h"
#include "ratematrix.h"

static int    er_entropy_target_f(double weight, void *params, double *ret_fx);
static double er_MeanEntropy(float **vec, int M, int K);
static double er_MeanRelativeEntropy(float **vec, int M, float *f, int K);

/*****************************************************************
 * 1. Miscellaneous functions for evoH3
 *****************************************************************/

/* Function:  p7_RateCreate()
 * Synopsis:  Allocate a p7_RATE structure
 * Incept:    ER, Tue Sep 27 14:29:47 2011  [Janelia]
 *
 * Purpose:  
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success.
 * 
 * Note:      
 *
 */            
P7_RATE *
p7_RateCreate(const ESL_ALPHABET *abc, int M, EMRATE *emR, ESL_SCOREMATRIX *S, EVOM evomodel, float fixtime, float betainf, float tol)
{ 
  P7_RATE *R = NULL;
  int      status;

  ESL_ALLOC(R, sizeof(P7_RATE));

  R->nrate    = 1;
  R->nbern    = 0;
  R->M        = M;
  R->abc_r    = abc;
  R->betainf  = betainf;
  R->tol      = tol;
  R->fixtime  = fixtime;
  R->evomodel = evomodel;
  if (!(R->evomodel == AIF || R->evomodel == AGAX)) {
    free(R);
    return NULL;
  }
  R->emR = emR;
  R->S   = S;

  // bellow here will be allocated with p7_RateAllocate()
  R->name  = NULL;
  R->acc   = NULL;
  R->desc  = NULL;

  R->pzero = NULL;
  R->pstar = NULL;
  R->pinfy = NULL;
  R->ins   = NULL;

  R->dtval = NULL;
  R->Pdt   = NULL;
  R->pdt   = NULL;

  /* the rates for transition probabilities */
  R->e1R  = NULL;

  R->allocated = FALSE;
  R->done      = FALSE;
  
  return R;

 ERROR:
  return NULL;
}

extern int 
p7_RateAllocate(P7_RATE *R)
{
  int   M = R->M;
  int   ndt;
  int   ndtl, ndtr, ndts;
  float dtpre, dtpos;
  float dtmin, dtmax;
  int   isfixtime = (R->fixtime >= 0.0)? TRUE : FALSE;
  int   m;
  int   t;
  int   status;

  if (R == NULL) return eslFAIL;
  if (R->allocated) return eslOK;
  
  if (R->M == 0) return eslFAIL;
  if (R->abc_r == NULL) return eslFAIL;

  
  /* store the match emissions form the training set */
  ESL_ALLOC(R->pzero,    sizeof(float *) * (M+1)         );  R->pzero[0] = NULL;
  ESL_ALLOC(R->pstar,    sizeof(float *) * (M+1)         );  R->pstar[0] = NULL;
  ESL_ALLOC(R->pinfy,    sizeof(float *) * (M+1)         );  R->pinfy[0] = NULL;
  ESL_ALLOC(R->ins,      sizeof(float *) * (M+1)         );  R->ins[0]   = NULL;
  ESL_ALLOC(R->pzero[0], sizeof(float  ) * (M+1) * R->abc_r->K);
  ESL_ALLOC(R->pstar[0], sizeof(float  ) * (M+1) * R->abc_r->K);
  ESL_ALLOC(R->pinfy[0], sizeof(float  ) * (M+1) * R->abc_r->K);
  ESL_ALLOC(R->ins[0],   sizeof(float  ) * (M+1) * R->abc_r->K);
  for (m = 1; m <= M; m ++) R->pzero[m] = R->pzero[0] + R->abc_r->K * m;
  for (m = 1; m <= M; m ++) R->pstar[m] = R->pstar[0] + R->abc_r->K * m;
  for (m = 1; m <= M; m ++) R->pinfy[m] = R->pinfy[0] + R->abc_r->K * m;
  for (m = 1; m <= M; m ++) R->ins[m]   = R->ins[0]   + R->abc_r->K * m;
  
  ndt    = (isfixtime)? 1 : NDT;
  ndts   = (isfixtime)? 0 : NDTS;
  ndtl   = (isfixtime)? 0 : NDTL;
  ndtr   = (isfixtime)? 0 : NDTR;
  dtmin  = (isfixtime)? R->fixtime : DTMIN;
  dtmax  = (isfixtime)? R->fixtime : DTMAX;
  dtpre  = (isfixtime)? R->fixtime : DTPRE;
  dtpos  = (isfixtime)? R->fixtime : DTPOS;
  R->ndt = ndt;
  
  if (ndt > 0) {
    ESL_ALLOC(R->dtval, sizeof(float) * (ndt));      
    R->dtval[0] = dtmin;
    for (t = 1; t <= ndtl; t ++) 
      R->dtval[t] = R->dtval[t-1] + (float)(dtpre-dtmin)/(float)ndtl;
    for (t = ndtl+1; t < ndtl+1+ndts; t ++) 
      R->dtval[t] = R->dtval[t-1] + (float)(1.0-dtpre)/(float)(ndts-1);
    for (t = ndtl+1+ndts; t < ndtl+1+2*ndts; t ++) 
      R->dtval[t] = R->dtval[t-1] + (float)(dtpos-1.0)/(float)(ndts-1);
    for (t = ndtl+1+2*ndts; t < ndt; t ++) 
      R->dtval[t] = R->dtval[t-1] + (float)(dtmax-dtpos)/(float)ndtr;
  }
  
  ESL_ALLOC(R->Pdt, sizeof(ESL_DMATRIX **) * (ndt));  
  ESL_ALLOC(R->pdt, sizeof(float       **) * (ndt));  
  for (t = 0; t < ndt; t ++) {
    ESL_ALLOC(R->Pdt[t], sizeof(ESL_DMATRIX *) * (M+1) );  R->Pdt[t][0] = NULL;
    ESL_ALLOC(R->pdt[t], sizeof(float       *) * (M+1) );  R->pdt[t][0] = NULL;
    
    for (m = 0; m <= M; m ++) 
      R->Pdt[t][m] = esl_dmatrix_Create(R->abc_r->K, R->abc_r->K);
    
    ESL_ALLOC(R->pdt[t][0],   sizeof(float   ) * (M+1) * R->abc_r->K);
    for (m = 1; m <= M; m ++) R->pdt[t][m] = R->pdt[t][0] + R->abc_r->K * m;
    
  } 

  /* the rates for transition probabilities */
  ESL_ALLOC(R->e1R, sizeof(E1_RATE *) * (M+1));
  for (m = 0; m <= M; m ++) {   
     R->e1R[m] = e1_rate_Create(R->abc_r, R->evomodel);
   }

  R->allocated = TRUE;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_RateCompare()
 * Synopsis:  Compare two p7_RATE structures
 * Incept:    ER, Tue Sep 27 14:29:47 2011  [Janelia]
 *
 * Purpose:  
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success.
 * 
 * Note:      
 *
 */            
double 
p7_RateCompare(P7_RATE *R1, P7_RATE *R2, double tol)
{ 
  double dist = 0.0;
  int    m;

  if (strcmp(R1->name, R2->name)   != 0)  return -1.0;

  if (R1->evomodel    != R2->evomodel)    return -1.0;
  if (R1->M           != R2->M)           return -1.0;
  if (R1->abc_r->type != R2->abc_r->type) return -1.0;
  if (R1->fixtime     != R2->fixtime)     return -1.0;
  if (R1->ndt         != R2->ndt)         return -1.0;
  if (R1->betainf     != R2->betainf)     return -1.0;
  if (R1->tol         != R2->tol)         return -1.0;

  for (m = 0; m <= R1->M; m ++) {
    if (esl_vec_FCompare(R1->pzero[m], R2->pzero[m], R1->abc_r->K, tol) != eslOK) return -1.0;
    if (esl_vec_FCompare(R1->pstar[m], R2->pstar[m], R1->abc_r->K, tol) != eslOK) return -1.0;
    if (esl_vec_FCompare(R1->pinfy[m], R2->pinfy[m], R1->abc_r->K, tol) != eslOK) return -1.0;
    if (esl_vec_FCompare(R1->ins[m],   R2->ins[m],   R1->abc_r->K, tol) != eslOK) return -1.0;
  }
  
  for (m = 0; m <= R1->M; m ++) 
    dist += e1_rate_Compare(R1->e1R[m], R2->e1R[m], tol);

  return dist;
}

/* Function:  p7_RateCopy()
 * Synopsis:  Copy a p7_RATE structure
 * Incept:    ER, Tue Sep 27 14:29:47 2011  [Janelia]
 *
 * Purpose:  
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success.
 * 
 * Note:      
 *
 */            
int
p7_RateCopy(P7_RATE *R, P7_RATE *Rcopy)
{ 
  int m;
  int t;
  int status;

  Rcopy->evomodel  = R->evomodel;
  Rcopy->nrate     = R->nrate;
  Rcopy->nbern     = R->nbern;
  Rcopy->M         = R->M;
  Rcopy->abc_r     = R->abc_r;
  Rcopy->fixtime   = R->fixtime;
  Rcopy->ndt       = R->ndt;
  Rcopy->betainf   = R->betainf;
  Rcopy->tol       = R->tol;
  Rcopy->emR       = R->emR;
  Rcopy->S         = R->S;
  Rcopy->allocated = R->allocated;
  Rcopy->done      = R->done;
  if (!R->allocated || !R->done) return eslOK; // do not go further is nothing is allocated  yet in the rate
  
  if ((status = esl_strdup(R->name,   -1, &(Rcopy->name)))   != eslOK) goto ERROR;
  if ((status = esl_strdup(R->acc,    -1, &(Rcopy->acc)))    != eslOK) goto ERROR;
  if ((status = esl_strdup(R->desc,   -1, &(Rcopy->desc)))   != eslOK) goto ERROR;

  for (t = 0; t < R->ndt; t++)
    Rcopy->dtval[t] = R->dtval[t];

  for (m = 0; m <= R->M; m ++) {
    esl_vec_FCopy(R->pzero[m], R->abc_r->K, Rcopy->pzero[m]);
    esl_vec_FCopy(R->pstar[m], R->abc_r->K, Rcopy->pstar[m]);
    esl_vec_FCopy(R->pinfy[m], R->abc_r->K, Rcopy->pinfy[m]);
    esl_vec_FCopy(R->ins[m],   R->abc_r->K, Rcopy->ins[m]);
    }

  for (t = 0; t < R->ndt; t++) 
    for (m = 0; m <= R->M; m ++) {
      esl_dmatrix_Copy(R->Pdt[t][m], Rcopy->Pdt[t][m]);
      esl_vec_FCopy(R->pdt[t][m], R->abc_r->K, Rcopy->pdt[t][m]);
    }

  for (m = 0; m <= R->M; m ++) 
      e1_rate_Copy(R->e1R[m], Rcopy->e1R[m]);
  
 ERROR:
  return status;
}

/* Function:  p7_RateClone(()
 * Synopsis:  clones an hmm rate
 * Incept:    ER, Wed Jan 25 13:41:23 EST 2012 [Janelia]
 *
 * Purpose:  
 *            
 * Args:      
 *
 * Returns:   pointer to cloned rate
 *
 */    
P7_RATE *        
p7_RateClone(P7_RATE *R)
{ 
  P7_RATE *Rclone = NULL;

  if (R == NULL) return Rclone;

  Rclone = p7_RateCreate(R->abc_r, R->M, R->emR, R->S, R->evomodel, R->fixtime, R->betainf, R->tol);
  p7_RateAllocate(Rclone);
  p7_RateCopy(R, Rclone);

  return Rclone;
}

/* Function:  p7_RateDestroy()
 * Synopsis:  Free a p7_RATE structure
 * Incept:    ER, Tue Sep 27 14:29:47 2011  [Janelia]
 *
 * Purpose:  
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success.
 * 
 * Note:      
 *
 */            
void
p7_RateDestroy(P7_RATE *R)
{ 
  int t;
  int m;

  if (R == NULL) return;

  if (R->name) free(R->name);
  if (R->acc)  free(R->acc);
  if (R->desc) free(R->desc);

  if (R->pzero    != NULL) {
    if (R->pzero[0] != NULL) free(R->pzero[0]);
    free(R->pzero);
  }
  if (R->pstar    != NULL) {
    if (R->pstar[0] != NULL) free(R->pstar[0]);
    free(R->pstar);
  }
  if (R->pinfy    != NULL) { 
    if (R->pinfy[0] != NULL) free(R->pinfy[0]);
    free(R->pinfy);
  }
  if (R->ins      != NULL) {
    if (R->ins[0]   != NULL) free(R->ins[0]);
    free(R->ins);
  }
  if (R->dtval    != NULL) free(R->dtval);
  
  if (R->Pdt != NULL) {
    for (t = 0; t < R->ndt; t ++) {
      for (m = 0; m <= R->M; m ++) 
	if (R->Pdt[t][m] != NULL) esl_dmatrix_Destroy(R->Pdt[t][m]);
      if (R->Pdt[t]    != NULL) free(R->Pdt[t]);
    }
    free(R->Pdt);
  }

  if (R->pdt != NULL) {
    for (t = 0; t < R->ndt; t ++) {
      if (R->pdt[t] != NULL) {
	if (R->pdt[t][0] != NULL) free(R->pdt[t][0]);
	free(R->pdt[t]);
      }
    }
    free(R->pdt);
  }

  if (R->e1R) {
    for (m = 0; m <= R->M; m ++) 
      e1_rate_Destroy(R->e1R[m]);
    free(R->e1R);
  }

  free(R); 
}

/* Function:  p7_RateDump()
 * Synopsis:  Dump a p7_RATE structure
 * Incept:    ER, Tue Sep 27 14:29:47 2011  [Janelia]
 *
 * Purpose:  
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success.
 * 
 * Note:      
 *
 */            
void
p7_RateDump(FILE *fp, P7_RATE *R)
{ 
  int   m;

  for (m = 0; m <= R->M; m ++) 
    e1_rate_Dump(fp, R->e1R[m]);

  if (R->pzero) {
    /* dump the pzero match probs */
    for (m = 0; m <= R->M; m ++) {
      esl_vec_FDump(fp, R->pzero[m], R->abc_r->K, NULL);
    }
  }
  if (R->pstar) {
    /* dump the pstar match probs */
    for (m = 0; m <= R->M; m ++) {
      esl_vec_FDump(fp, R->pstar[m], R->abc_r->K, NULL);
    }
  }
  if (R->pinfy) {
    /* dump the pinfy match probs */
    for (m = 0; m <= R->M; m ++) {
      esl_vec_FDump(fp, R->pinfy[m], R->abc_r->K, NULL);
    }
  }
  if (R->ins) {
    /* dump the ins match probs */
    for (m = 0; m <= R->M; m ++) {
      esl_vec_FDump(fp, R->ins[m], R->abc_r->K, NULL);
    }
  }
}

/* Function:  p7_RateValidate()
 * Synopsis:  Validate a p7_RATE structure
 * Incept:    ER, Thu Sep 29 11:10:59 2011  [Janelia]
 *
 * Purpose:  
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success.
 * 
 * Note:      
 *
 */            
int
p7_RateValidate(P7_RATE *R, char *errbuf)
{ 
  int m;
  int t;
  int status;

  if (R == NULL) return eslOK;

  if (R->evomodel == UN) ESL_XFAIL(eslFAIL, errbuf, "Unknown evomodel for this rate");
  if (R->M < 0) ESL_XFAIL(eslFAIL, errbuf, "Rate has bad dimensions");
  
  if (R->pzero) {
    for (m = 0; m <= R->M; m ++) {
      if (esl_vec_FValidate(R->pzero[m], R->abc_r->K, R->tol, errbuf) != eslOK) {
	esl_vec_FDump(stdout, R->pzero[m], R->abc_r->K, NULL);
	ESL_XFAIL(eslFAIL, errbuf, "match zero probabilities at state %d/%d  did not validate", m, R->M);
      }
    }
  }

  if (R->pstar) {
    for (m = 0; m <= R->M; m ++) {
      if (esl_vec_FValidate(R->pstar[m], R->abc_r->K, R->tol, errbuf) != eslOK)  {
	esl_vec_FDump(stdout, R->pstar[m], R->abc_r->K, NULL);
	ESL_XFAIL(eslFAIL, errbuf, "match star probabilities at state %d/%d did not validate", m, R->M);
      }
    }
  }

  if (R->pinfy) {
    for (m = 0; m <= R->M; m ++) {
      if (esl_vec_FValidate(R->pinfy[m], R->abc_r->K, R->tol, errbuf) != eslOK) {  
	esl_vec_FDump(stdout, R->pinfy[m], R->abc_r->K, NULL);
	ESL_XFAIL(eslFAIL, errbuf, "match infy probabilities at state %d/%d  did not validate", m, R->M);
      }
    }
  }
  if (R->ins) {
    for (m = 0; m <= R->M; m ++) {
      if (esl_vec_FValidate(R->ins[m], R->abc_r->K, R->tol, errbuf) != eslOK) {  
	esl_vec_FDump(stdout, R->ins[m], R->abc_r->K, NULL);
	ESL_XFAIL(eslFAIL, errbuf, "match infy probabilities at state %d/%d  did not validate", m, R->M);
      }
    }
  }
   
  for (t = 0; t < R->ndt; t ++) {
    if (R->pdt[t]) {
      for (m = 0; m <= R->M; m ++) {
	if (esl_vec_FValidate(R->pdt[t][m], R->abc_r->K, R->tol, errbuf) != eslOK) {
	  esl_vec_FDump(stdout, R->pdt[t][m], R->abc_r->K, NULL);
	  ESL_XFAIL(eslFAIL, errbuf, "dt for t=%d probabilities at state %d/%d  did not validate", t, m, R->M);
	}
      }
    }
  } 

  if (R->e1R) {
    for (m = 0; m <= R->M; m ++) 
      if ((status = e1_rate_Validate(R->e1R[m], R->tol, errbuf)) != eslOK) goto ERROR;  
  }
  
  return eslOK;

  ERROR:
    return status;
}

/* Function:  p7_CalculateRate()
 * Synopsis:  
 * Incept:    ER, Tue Sep 27 14:29:47 2011  [Janelia]
 *
 * Purpose:  
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success.
 * 
 * Note:      
 *
 */                       
extern int 
p7_RateCalculate(const P7_HMM *hmm, const P7_BG *bg, P7_RATE *R, char *errbuf, int verbose)
{ 
  enum emevol_e  emevol;
  double         time;
  int            M = hmm->M;
  int            ndt;
  int            ndtl, ndtr, ndts;
  float          dtpre, dtpos;
  float          dtmin, dtmax;
  int            isfixtime = (R->fixtime >= 0.0)? TRUE : FALSE;
  int            t;
  int            m;
  int            i, j;
  int            status;

  if (R->done) return eslOK;
  
  // check if we need to allocate 
  status = p7_RateAllocate(R);
  if (status != eslOK) goto ERROR;
  
  if ((status = esl_strdup(hmm->name,   -1, &(R->name)))   != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->acc,    -1, &(R->acc)))    != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->desc,   -1, &(R->desc)))   != eslOK) goto ERROR;

  /* the emissions rate matrix, if provided */
  if      (R->S && R->emR) emevol = emBYRATE; /* if both given, use the rate, it is faster */
  else if (R->S)           emevol = emBYSCMX;
  else if (R->emR)         emevol = emBYRATE;
  else                     emevol = emNONE;

  switch(emevol) {
  case emBYSCMX: 
    printf("not integrated yet\n"); status = eslFAIL; goto ERROR;
    break;
  case emBYRATE: 
    for (m = 0; m <= R->M; m ++) 
      ratematrix_emrate_Copy(R->emR, R->e1R[m]->em);
   break;
  case emNONE: 
    R->e1R = NULL;
    break;
  default: printf("not such type of emission evolution\n"); status = eslFAIL; goto ERROR;
  }

  switch(emevol) {
  case emBYSCMX: 
    printf("not integrated yet\n"); status = eslFAIL; goto ERROR;
    break;
  case emBYRATE:
    for (m = 0; m <= R->M; m ++) 
      esl_vec_FCopy(hmm->mat[m], R->abc_r->K, R->pstar[m]);

    /* increase the information content of the match emission to pzero */
    status = p7_CalculatePzero(R, hmm, errbuf, verbose);
    if (status != eslOK) goto ERROR;

    /* assign the given match emission to pstar */
    for (m = 0; m <= R->M; m ++) 
      esl_vec_FCopy(hmm->mat[m], R->abc_r->K, R->pstar[m]);
 
    /* calculate pinfy */
    status = p7_CalculatePinfy(R, hmm, bg, errbuf, verbose);
   if (status != eslOK) goto ERROR;

    for (m = 0; m <= R->M; m ++) {
      ratematrix_FRateFromExchange(R->e1R[m]->em->E, R->pstar[m], R->e1R[m]->em->Qstar);
      ratematrix_FRateFromExchange(R->e1R[m]->em->E, R->pinfy[m], R->e1R[m]->em->Qinfy);
    }    
    
    /* assign the given ins emission to inst */
    for (m = 0; m <= R->M; m ++) 
      esl_vec_FCopy(hmm->ins[m], R->abc_r->K, R->ins[m]);
    
    for (t = 0; t < R->ndt; t ++) {
      if (R->dtval[t] < 1.0) time = R->dtval[t] / (1.0 - R->dtval[t]);
      else                   time = R->dtval[t] - 1.0;

      for (m = 0; m <= R->M; m ++) {
 	if (R->dtval[t] < 1.0)
	  status = ratematrix_CalculateConditionalsFromRate(time, R->e1R[m]->em->Qstar, R->Pdt[t][m], R->tol, errbuf, verbose);
	else 
	  status = ratematrix_CalculateConditionalsFromRate(time, R->e1R[m]->em->Qinfy, R->Pdt[t][m], R->tol, errbuf, verbose);
	
	esl_vec_FSet(R->pdt[t][m], hmm->abc->K, 0.0);	
	for (i = 0; i < hmm->abc->K; i ++) 
	  for (j = 0; j < hmm->abc->K; j ++) 
	    R->pdt[t][m][i] += (R->dtval[t] < 1.0)? R->Pdt[t][m]->mx[j][i] * R->pzero[m][j] : R->Pdt[t][m]->mx[j][i] * R->pstar[m][j];	

	esl_vec_FNorm(R->pdt[t][m], hmm->abc->K);
     }
    }
    
 #if 0
    double  *psat_star = NULL;
    double  *psat_infy = NULL;
    for (m = 0; m <= R->M; m ++) {
      printf("\nm=%d\n", m);
      esl_vec_FDump(stdout, R->pzero[m], R->abc_r->K, NULL);
      esl_vec_FDump(stdout, R->pstar[m], R->abc_r->K, NULL);
      esl_vec_FDump(stdout, R->pinfy[m], R->abc_r->K, NULL);   
      esl_vec_FDump(stdout, R->ins[m],   R->abc_r->K, NULL);   
 
      for (t = 0; t < R->ndt; t ++) {
	esl_vec_FDump(stdout, R->pdt[t][m], R->abc_r->K, NULL);
      }
#if 0
      ratematrix_SaturationTime(R->e1R[m].em->Qstar, NULL, &psat_star, tol, errbuf, verbose);
      ratematrix_SaturationTime(R->e1R[m].em->Qinfy, NULL, &psat_infy, tol, errbuf, verbose);
      esl_vec_FDump(stdout, psat_star, R->abc_r->K, NULL);
      esl_vec_FDump(stdout, psat_infy, R->abc_r->K, NULL);
      free(psat_star); psat_star = NULL;
      free(psat_infy); psat_infy = NULL;
#endif
    }
#endif

   break;
  case emNONE: 
    R->e1R = NULL;
    break;
  default: printf("not such type of emission evolution\n"); status = eslFAIL; goto ERROR;
  }

  if (p7_RateTransitions(hmm, R, errbuf, verbose) != eslOK) { status = eslFAIL; goto ERROR; }
  if (p7_RateValidate(R, errbuf) != eslOK) { status = eslFAIL; goto ERROR; }

#if 0
  printf("^^ testing the rate\n");
  status = p7_RateTest(hmm, R, bg, errbuf, verbose);
  if (status != eslOK) { printf("p7_RateTest() failed\n"); goto ERROR; }
#endif

  R->done = TRUE;
  return eslOK;  

 ERROR:
  return status;
}


extern int 
p7_RateConstruct(const P7_HMM *hmm, const P7_BG *bg, HMMRATE *hmmrate, P7_RATE **ret_R, char *errbuf, int verbose)
{ 
  P7_RATE       *R = NULL;
  int            status;

  R = p7_RateCreate(hmm->abc, hmm->M, hmmrate->emR, hmmrate->S, hmmrate->evomodel, hmmrate->fixtime, hmmrate->betainf, hmmrate->tol);
  if (R == NULL) { status = eslFAIL; goto ERROR; }
		   
  status = p7_RateCalculate(hmm, bg, R, errbuf, verbose);
  if (status != eslOK) goto ERROR;
    
  if (0&&hmmrate->statfp) {
    fprintf(hmmrate->statfp, "%*s time\tME\tMRE\n", (int)strlen(hmm->name), "#anchor hmm");
    fprintf(hmmrate->statfp, "%s %4.4f\t%2.4f\t%2.4f \n",
	    hmm->name, 0.0,         er_MeanEntropy(R->pzero, R->M, R->abc_r->K), er_MeanRelativeEntropy(R->pzero, R->M, bg->f, R->abc_r->K));
    fprintf(hmmrate->statfp, "%s %4.4f\t%2.4f\t%2.4f \n",
	    hmm->name, 1.0,         p7_MeanMatchEntropy(hmm),                    p7_MeanMatchRelativeEntropy(hmm, bg));
    fprintf(hmmrate->statfp, "%s %4.4f\t%2.4f\t%2.4f \n",
	    hmm->name, eslINFINITY, er_MeanEntropy(R->pinfy, R->M, R->abc_r->K), er_MeanRelativeEntropy(R->pinfy, R->M, bg->f, R->abc_r->K));
  }

  *ret_R = R;
  return eslOK;  

 ERROR:
  if (R) p7_RateDestroy(R);
  return status;
}


int
p7_RateTransitions(const P7_HMM *hmm, P7_RATE *R, char *errbuf, int verbose)
{
  struct rateparam_s  rateparam;	
  E1_RATE            *e1R;
  double              betainfk;   /* if beta > betainf, need to modify it state by state */
  double              gammaM;
  double              gammaD;
  double              betaM;
  double              eta;
  int                 m;
  int                 status;
  
  if (R->evomodel != AGAX) ESL_XFAIL(eslFAIL, errbuf, "not a plan7 evomodel");

  for (m = 0; m <= R->M; m ++) {

    e1R = R->e1R[m];
    if (m > 0 && hmm->t[m][p7H_DM] == 0.0) 
      ESL_XFAIL(eslFAIL, errbuf, "Model evolved to infinit time at position k=%d. Cannot extract rates.\n", m);
    
    /* Assign fundamental probabilities from the HMM transitions */
    betainfk = R->betainf;
    betaM    = hmm->t[m][p7H_MI]; 
    eta      = hmm->t[m][p7H_II];
    gammaD   = hmm->t[m][p7H_DD];

    if (betaM >= betainfk) {
      betainfk = (betaM + 0.05 < 1.0)? betaM + 0.05 : 1.0;
    }
  
    gammaM = (betaM < 1.0)? hmm->t[m][p7H_MD] / (1.0 - betaM) : 0.0;
 
    rateparam.muAM =  e1_rate_CalculateAncestralDeletionRates(gammaM, 1.0);
    rateparam.muAD =  e1_rate_CalculateAncestralDeletionRates(gammaD, 1.0);
    rateparam.muAI =  0.0;
    rateparam.muEM = -1.0;
    rateparam.muED =  0.0;
    rateparam.muI  =  0.0;
    rateparam.ldEM = -1.0;
    rateparam.ldED =  0.0;
    rateparam.ldI  = -1.0;
    rateparam.sI   =  eta;
    rateparam.sD   =  0.0;
    rateparam.vI   =  0.0;
    rateparam.vD   =  0.0;
    rateparam.rM   = -1.0;
    rateparam.rD   = -1.0;
    rateparam.rI   = -1.0;
    
    status = e1_rate_CalculateInsertRates(R->evomodel, &rateparam, betainfk, betainfk, eta, betaM, 0.0, R->tol, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    status = e1_rate_AssignTransitionsFromRates(e1R, rateparam, errbuf, verbose);
    if (status != eslOK) goto ERROR;
    
#if 0
    if (m == 59) {
      printf("\nhmm[k=%d] | MM %f MI %f MD %f IM %f II %f DM %f DD %f \n", m, hmm->t[m][p7H_MM], hmm->t[m][p7H_MI], hmm->t[m][p7H_MD], 
	     hmm->t[m][p7H_IM], hmm->t[m][p7H_II], hmm->t[m][p7H_DM], hmm->t[m][p7H_DD]);
      printf("%s-RATES:state=%d\tmuAM=%.4f\tmuAD=%.4f\tldEM=%.4f\tmuEM=%.4f\tldED=%.4f\tmuED=%.4f\tsI=%.4f\tbetainfk=%.4f\n", 
	     e1_rate_EvomodelType(R->evomodel), m, 
	     e1R->muA[e1R_S], e1R->muA[e1R_D], 
	     e1R->ldE[e1R_S], e1R->muE[e1R_S], 
	     e1R->ldE[e1R_D], e1R->muE[e1R_D], e1R->sI, betainfk);
    }
#endif
  }

  return eslOK;
  
 ERROR:
  return status;
}

/* Function:  p7_EvolveFromRate()
 * Synopsis:  
 * Incept:    ER, Tue Sep 27 14:29:47 2011  [Janelia]
 *
 * Purpose:  
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success.
 * 
 * Note:      
 *
 */            
int 
p7_EvolveFromRate(FILE *statfp, P7_HMM *hmm, const P7_RATE *R, const P7_BG *bg, double time, char *errbuf, int verbose)
 {
  struct e1_params  p;				
  E1_RATE          *e1R;
  double            gammaM, gammaD, gammaI;
  double            betaM, betaD;
  double            eta;
  double            tM, tD, tI;
  double            itM, itD, itI;
  double            fsmall = 1e-8;
  int               m; 
  int               k;
  int               t;
  int               dtselect;
  int               i;
  int               status;

  if (time < 0.0) ESL_XFAIL(eslFAIL, errbuf, "Not a valid time (%f)\n", time);
  if (R->M != hmm->M) ESL_XFAIL(eslFAIL, errbuf, "Rate dim (%d) does not correspond to HMM dim (%d)\n", R->M, hmm->M);
  
  if (hmm->name) free(hmm->name); hmm->name = NULL;
  if (hmm->acc)  free(hmm->acc);  hmm->acc  = NULL;
  if (hmm->desc) free(hmm->desc); hmm->desc = NULL;
  if ((status = esl_strdup(R->name,   -1, &(hmm->name)))   != eslOK) goto ERROR;
  if ((status = esl_strdup(R->acc,    -1, &(hmm->acc)))    != eslOK) goto ERROR;
  if ((status = esl_strdup(R->desc,   -1, &(hmm->desc)))   != eslOK) goto ERROR;

  if (time < 1.0) {
    for (t = 0; t < R->ndt; t++) {
      if (R->dtval[t] >= time) { dtselect = t; break; }
    }
    if (t == R->ndt) dtselect = R->ndt;	
  }
  else {
    for (t = R->ndt-1; t >= 0; t--) {
      if (R->dtval[t] <= time) { dtselect = t; break; }
    }
    if (t < 0) dtselect = R->ndt;	
  }
  
  for (k = 0; k <= hmm->M; k++)
    {
      /* pmat */ 
      for (i = 0; i < hmm->abc->K; i ++) 
	hmm->mat[k][i] = (dtselect < R->ndt)? R->pdt[dtselect][k][i] : R->pstar[k][i];	    
      /* pins */ 
      esl_vec_FCopy(R->ins[k], R->abc_r->K, hmm->ins[k]);
    }
  
  p.tol    = R->tol;          
  p.errbuf = errbuf;          
  p.fsmall = fsmall;    
  
  for (m = 0; m <= R->M; m ++) {
    
    e1R = R->e1R[m];
 
    /* Calculate the fundamental probabilites */
    if (e1R->muE[e1R_S] == 0.0 && e1R->muE[e1R_D] == 0.0 && e1R->ldE[e1R_S] == 0.0 && e1R->ldE[e1R_D] == 0.0) {
      /* not a valide case, return a trivial no-evolution solution */
      gammaM = 0.0;
      gammaD = 0.0;
      gammaI = 0.0;
      betaM  = 0.0;
      betaD  = 0.0;
      eta    = 0.0;
    }
    else {
      eta = e1R->sI;
      // hack to have t_II = 0 for time = 0
      if (time < 1.0) {
	eta *= (exp(time) - 1.0) / (exp(1.0) - 1.0);
      }

      if (time < eslINFINITY) {
	gammaM = 1.0 - exp(-e1R->muA[e1R_S] * time);
	gammaD = 1.0 - exp(-e1R->muA[e1R_D] * time);
 	gammaI = 1.0 - exp(-e1R->muA[e1R_I] * time);
      }
      else {
	gammaM = 1-fsmall;
	gammaD = 1-fsmall;
	gammaI = 1-fsmall;
      }
            
      p.time           = time;
      p.rateparam.muEM = e1R->muE[e1R_D];         
      p.rateparam.muED = 0.0;         
      p.rateparam.muI  = e1R->muE[e1R_I];         
      p.rateparam.ldEM = e1R->ldE[e1R_S];         
      p.rateparam.ldED = 0.0;         
      p.rateparam.ldI  = 0.0;
      if (R->evomodel == AIF) p.rateparam.rI = e1R->rI;          
      if (R->evomodel == AGA) p.rateparam.sI = e1R->sI;          
      p.special_case   = FALSE;
    
      status = e1_model_AG_BetaFunc(&p, &betaM, &betaD);
      if (status != eslOK) ESL_XFAIL(eslFAIL, p.errbuf, "AGA_betaM failed for m=%d", m);
     }
  
    /* Assign fundamental probabilities to the HMM transitions (avoid absolute zeros) */
    hmm->t[m][p7H_MI] = (betaM       > fsmall)? betaM                  : fsmall;
    hmm->t[m][p7H_II] = (eta         > fsmall)? eta                    : fsmall;
    
    hmm->t[m][p7H_MD] = (1.0 - betaM > fsmall)? (1.0 - betaM) * gammaM : fsmall;
    hmm->t[m][p7H_DD] = (1.0 - betaD > fsmall)? (1.0 - betaD) * gammaD : fsmall;
    
    hmm->t[m][p7H_MM] = ((1.0 - betaM) * (1.0 - gammaM) > fsmall)? (1.0 - betaM) * (1.0 - gammaM) : fsmall; 
    hmm->t[m][p7H_DM] = ((1.0 - betaD) * (1.0 - gammaD) > fsmall)? (1.0 - betaD) * (1.0 - gammaD) : fsmall; 
    hmm->t[m][p7H_IM] = ((1.0 - eta)   * (1.0 - gammaI) > fsmall)? (1.0 - eta)   * (1.0 - gammaI) : fsmall; 

    /* last node is a special: TMD should be zero, TDM should be one  */
    if (m == hmm->M && hmm->t[m][p7H_MD] > 0.0)   hmm->t[m][p7H_MD] = 0.0;
    if (m == hmm->M && hmm->t[m][p7H_DD] > 0.0) { hmm->t[m][p7H_DD] = 0.0; hmm->t[m][p7H_DM] = 1.0; }

    /* normalize (this takes care of D->I = 0 and I->D = 0) */
    tM = hmm->t[m][p7H_MM] + hmm->t[m][p7H_MD] + hmm->t[m][p7H_MI];
    itM = (tM > 0.)? 1./tM : 0.0;
    hmm->t[m][p7H_MM] *= itM;
    hmm->t[m][p7H_MD] *= itM;
    hmm->t[m][p7H_MI] *= itM;
    tD = hmm->t[m][p7H_DM] + hmm->t[m][p7H_DD];
    itD = (tD > 0.)? 1./tD : 0.0;
    hmm->t[m][p7H_DM] *= itD;
    hmm->t[m][p7H_DD] *= itD;
    tI = hmm->t[m][p7H_IM] + hmm->t[m][p7H_II];
    itI = (tI > 0.)? 1./tI : 0.0;
    hmm->t[m][p7H_IM] *= itI;
    hmm->t[m][p7H_II] *= itI;

 
#if 0
    if (1||verbose) {
      printf("k=%d t=%f | MM %f MI %f MD %f IM %f II %f DM %f DD %f \n", m, time, hmm->t[m][p7H_MM], hmm->t[m][p7H_MI], hmm->t[m][p7H_MD], 
	     hmm->t[m][p7H_IM], hmm->t[m][p7H_II], hmm->t[m][p7H_DM], hmm->t[m][p7H_DD]);
    }
#endif
  }

  if (0&&statfp) {
    fprintf(statfp, "%*s time\tME\tMRE\n", (int)strlen(hmm->name), "#ehmm");
    fprintf(statfp, "%s %4.4f\t%2.4f\t%2.4f \n", hmm->name, time, p7_MeanMatchEntropy(hmm), p7_MeanMatchRelativeEntropy(hmm, bg));
  }
  
#if 0
  if (1||verbose) 
    p7_hmm_Dump(stdout, hmm);
#endif

  if ((status = p7_hmm_SetConsensus(hmm, NULL)) != eslOK) return status;
  if ((status = p7_hmm_Validate(hmm, errbuf, 0.0001)) != eslOK) { printf("p7_hmm_Validate %s\n", errbuf); return status; }
  
 return eslOK;
 
 ERROR:
 return status;
}

/* Function:  p7_Evolve()
 * Synopsis:  given a parameterized HMM, change the parameters
 *            to a different time t'/t.
 * Incept:    ER, Tue Sep 27 14:29:47 2011  [Janelia] 
 *
 * Purpose:  
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success.
 * 
 * Note:      
 *
 */            
int
p7_Evolve(P7_HMM *hmm, const P7_BG *bg, double time, HMMRATE *hmmrate, char *errbuf, int verbose)
{ 
  P7_RATE *R = NULL;
  int      status;
  
  if ((status = p7_RateConstruct(hmm, bg, hmmrate, &R, errbuf, verbose))               != eslOK) goto ERROR; 
  if ((status = p7_EvolveFromRate(hmmrate->statfp, hmm, R, bg, time, errbuf, verbose)) != eslOK) goto ERROR; 
  
  if (R != NULL) p7_RateDestroy(R);
  return eslOK;
  
 ERROR:
  return status;
}

int
p7_RateTest(const P7_HMM *h1, P7_RATE *R, const P7_BG *bg, char *errbuf, int verbose)
{ 
  P7_HMM  *h2   = p7_hmm_Clone(h1);
  double   time = 1.0;
  int      k;
  int      a;
  int      t;
  int      status;
  
  p7_EvolveFromRate(NULL, h2, R, bg, time, errbuf, verbose);

  for (k = 0; k <= h1->M; k++)	/* (it's safe to include 0 here.) */
    {
      if ((status = esl_vec_FCompare(h1->mat[k], h2->mat[k], h1->abc->K, R->tol)) != eslOK) {
	for (a = 0; a < h1->abc->K; a ++) printf("mat[%d] %f %f diff=%f %f\n", k, h1->mat[k][a], h2->mat[k][a],
						 fabs(h1->mat[k][a]-h2->mat[k][a]), 2.0*fabs(h1->mat[k][a]-h2->mat[k][a])/fabs(h1->mat[k][a]+h2->mat[k][a]));
	goto ERROR;
      }
      
      if ((status = esl_vec_FCompare(h1->ins[k], h2->ins[k], h1->abc->K, R->tol)) != eslOK) {
	for (a = 0; a < h1->abc->K; a ++) printf("ins[%d] %f %f diff=%f %f\n", k, h1->ins[k][a], h2->ins[k][a],
						 fabs(h1->ins[k][a]-h2->ins[k][a]), 2.0*fabs(h1->ins[k][a]-h2->ins[k][a])/fabs(h1->ins[k][a]+h2->ins[k][a]));
	goto ERROR;
      }
      
    if ((status = esl_vec_FCompare(h1->t[k],   h2->t[k],   7,          R->tol)) != eslOK) {
      for (t = 0; t < 7; t ++) printf("k = %d/%d t[%d] %f %f diff=%f %f\n", k, h1->M, t, h1->t[k][t], h2->t[k][t],
				      fabs(h1->t[k][t]-h2->t[k][t]), 2.0*fabs(h1->t[k][t]-h2->t[k][t])/fabs(h1->t[k][t]+h2->t[k][t]));
	goto ERROR;
      }
      
    }
       
  p7_hmm_Destroy(h2);
  return eslOK;
  
 ERROR:
  if (h2) p7_hmm_Destroy(h2);
  return status;
}

/* Function:  p7_RestoreHMMmat()
 * Synopsis:  
 * Incept:    ER, Mon Mar 26 13:08:58 EDT 2012  [Janelia]
 *
 * Purpose:  Restore the match emissions to their 
 *           original values. This is used in the pipeline
 *           before the MSV filter.
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success.
 * 
 * Note:      
 *
 */            
int 
p7_RestoreHMMmat(P7_HMM *hmm, P7_RATE *R, char *errbuf, int verbose)
{
  int k;

   for (k = 0; k <= hmm->M; k++) esl_vec_FCopy(R->pstar[k], hmm->abc->K, hmm->mat[k]);

  return eslOK;
}

int 
p7_CalculatePzero(P7_RATE *R, const P7_HMM *hmm, char *errbuf, int verbose)
{
  double     etarget;
  double     cut;
  int        m;
  int        status;
 
  esl_vec_FCopy(hmm->mat[0], R->abc_r->K, R->pzero[0]);
  for (m = 1; m <= R->M; m ++) {   
    esl_vec_FCopy(hmm->mat[m], R->abc_r->K, R->pzero[m]);
    esl_vec_FNorm(R->pzero[m], hmm->abc->K);       
  }

#if 1
  etarget = p7_MeanMatchEntropy(hmm)/8.0;
  status = er_EntropyWeight(R->pzero, R->M, R->abc_r->K, hmm->mat, etarget, &cut);
  if (status != eslOK) goto ERROR;
#endif

  if (verbose) {
    printf("mean entropy pzero %f  pmat %f\n", 
	   er_MeanEntropy(R->pzero, R->M, R->abc_r->K), p7_MeanMatchEntropy(hmm));
  }

  if (verbose) {
    for (m = 1; m <= R->M; m ++) {   
      printf("\nm %d\n", m);
      esl_vec_FDump(stdout, R->pzero[m], R->abc_r->K, NULL);
      esl_vec_FDump(stdout, hmm->mat[m], R->abc_r->K, NULL);
    } 
  }
  
   return eslOK;
  
 ERROR:
 return status;
}

int 
p7_CalculatePinfy(P7_RATE *R, const P7_HMM *hmm, const P7_BG *bg, char *errbuf, int verbose)
{
    P7_HMM    *h2 = NULL;
    double     weight = 0.5;
    int        m;
    int        status;
    
#if 0
    if ((h2 = p7_hmm_Clone(hmm)) == NULL) { status = eslEMEM; goto ERROR; }
    
    /* scale the match transitions */
    p7_hmm_Scale(h2, weight);
    p7_ParameterEstimation(h2, NULL);
    
    for (m = 0; m <= R->M; m ++) {   
      esl_vec_FAdd(h2->mat[m], bg->f, R->abc_r->K);
      esl_vec_FNorm(h2->mat[m], h2->abc->K);      
      esl_vec_FCopy(h2->mat[m], R->abc_r->K, R->pinfy[m]);
    } 
      
    if (verbose) {
      printf("relative entropy pinfy %f pmat %f weigth %f\n", 
	     p7_MeanMatchRelativeEntropy(h2, bg), p7_MeanMatchRelativeEntropy(hmm, bg), weight);
    }
#endif
    
    /* assign the match emission of bg to infty */
    for (m = 0; m <= R->M; m ++) 
      esl_vec_FCopy(bg->f, R->abc_r->K, R->pinfy[m]);

    if (verbose) {
      printf("mean entropy pinfty %f  pmat %f\n", 
	     er_MeanEntropy(R->pinfy, R->M, R->abc_r->K), p7_MeanMatchEntropy(hmm));
    }

#if 0
    if (verbose) {
      for (m = 1; m <= R->M; m ++) {   
	printf("\nm %d\n", m);
	esl_vec_FDump(stdout, hmm->mat[m], R->abc_r->K, NULL);
	esl_vec_FDump(stdout, R->pinfy[m], R->abc_r->K, NULL);
      } 
    }    
#endif

    if (h2)  p7_hmm_Destroy(h2);
    return eslOK;
    
 ERROR:
    if (h2)  p7_hmm_Destroy(h2);
    return status;
}

/* Function:  p7_EntropyWeight()
 * Incept:    SRE, Fri May  4 15:32:59 2007 [Janelia]
 *
 * Purpose:   Use the "entropy weighting" algorithm to determine
 *            what effective sequence number we should use, and 
 *            return it in <ret_Neff>. 
 *            
 *            Caller provides a count-based <hmm>, and the
 *            Dirichlet prior <pri> that's to be used to parameterize
 *            models; neither of these will be modified. 
 *            Caller also provides the relative entropy
 *            target in bits in <etarget>. 
 *            
 *            <ret_Neff> will range from 0 to the true number of
 *            sequences counted into the model, <hmm->nseq>.
 *
 * Returns:   <eslOK> on success. 
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
er_EntropyWeight(float **prob, int M, int K, float **pref, double etarget, double *ret_cut)
{
  ESL_ROOTFINDER         *R = NULL;
  struct entropy_param_s  p;
  double                  cut;
  double                  maxcut;
  double                  c;
  double                  fx;
  double                  tol = 0.01;
  int                     m;
  int                     status;


  cut = 0.00;
  maxcut = esl_vec_FMax(prob[1], K);
  for (m = 2; m <= M; m++) {
    c = (double)esl_vec_FMax(prob[m], K);
    maxcut = (c>maxcut)? c : maxcut;
  }

  /* Store parameters in the structure we'll pass to the rootfinder
   */
  p.pref    = pref;
  p.prob    = prob;
  p.M       = M;
  p.K       = K;
  p.etarget = etarget;
  
  if ((status = er_entropy_target_f(cut, &p, &fx)) != eslOK) goto ERROR;

  if (fabs(fx) > tol)
    {
      if ((R = esl_rootfinder_Create(er_entropy_target_f, &p)) == NULL) {status = eslEMEM; goto ERROR;}
      esl_rootfinder_SetAbsoluteTolerance(R, tol); /* getting Neff to ~2 sig digits is fine */
      if ((status = esl_root_Bisection(R, 0., maxcut, &cut)) != eslOK) goto ERROR;

      esl_rootfinder_Destroy(R);
    }
  
  *ret_cut = cut;
  return eslOK;

 ERROR:
  *ret_cut = 0.0;
  return status;
}

/* Evaluate fx = mean entropy - etarget, which we want to be = 0,
 * for effective sequence number <x>.
 */
static int
er_entropy_target_f(double cut, void *params, double *ret_fx)
{
  struct entropy_param_s *p = (struct entropy_param_s *) params;
  double                  rcut;
  int                     m;
  int                     k;

  for (m = 1; m <= p->M; m++) {
    esl_vec_FCopy(p->pref[m], p->K, p->prob[m]);
    rcut = ESL_MIN(cut, (double)esl_vec_FMax(p->prob[m], p->K));
    for (k = 0; k < p->K; k++) {
      p->prob[m][k] = (p->prob[m][k] >= rcut)? p->prob[m][k] : 0.0;
    }
        
    esl_vec_FNorm(p->prob[m], p->K);
  }
  
  *ret_fx = er_MeanEntropy(p->prob, p->M, p->K) - p->etarget;
  return eslOK;
}


static double
er_MeanEntropy(float **vec, int M, int K)
{
  int    m;
  double H = 0.;

  for (m = 1; m <= M; m++)
    H += esl_vec_FEntropy(vec[m], K);
  H /= (double) M;
  return H;
}

static double
er_MeanRelativeEntropy(float **vec, int M, float *f, int K)
{
  int    m;
  double KL = 0.;

  for (m = 1; m <= M; m++)
    KL += esl_vec_FRelEntropy(vec[m], f, K);
  KL /= (double) M;
  return KL;
}

/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef evoHMMER_TESTDRIVE

static int
utest_hmmevolve(FILE *hmmfp, P7_HMM *hmm, double etaz, double etainf, int N, double tol, char *errbuf, int verbose)
{
  char    *msg = "evoHMMER unit test failed";
  P7_RATE *R = NULL;
  double   ttz = 0.0;
  double   tt;
  int      x;
  int      status;

  if (p7_CalculateRate(hmm, NULL, NULL, &R, etaz, etainf, tol, errbuf, verbose) != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
  p7_RateDump(stdout, R);

  for (x = 0; x <= 20*N; x ++) {
    tt = (x<=N)? ttz + x*(1.0-ttz)/N : ttz + x*(1.0-ttz)/N;
    printf("\nTIME = %f\n", tt);
    if (p7_EvolveFromRate(hmm, R, tt, tol, errbuf, verbose)  != eslOK) { printf("%s\n", errbuf); esl_fatal(msg); }
    if ((status = p7_hmmfile_WriteASCII(hmmfp,  -1, hmm))         != eslOK) ESL_FAIL(status, errbuf, "HMM save failed");
  }

  if (R != NULL) p7_RateDestroy(R);

  return eslOK;
}

#endif /*evoHMMER_TESTDRIVE*/

  

/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef evoHMMER_TESTDRIVE

/* gcc -o evohmmer_utest -g -Wall -I../easel -L../easel -I. -L. -DevoHMMER_TESTDRIVE evohmmer.c ratematrix.c -lhmmer -leasel -lm
 * ./evohmmer_utest ../tutorial/fn3.hmm ../tutorial/evofn3.hmm  
 */
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type       default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",         eslARG_NONE,   FALSE, NULL, NULL,   NULL, NULL, NULL, "show brief help on version and usage",              0 },
  { "--etaz",      eslARG_REAL,  "0.01", NULL, "x>=0", NULL, NULL, NULL, "eta at time zero",                                 0 },
  { "--etainfty", eslARG_REAL,  "0.99", NULL, "x>=0", NULL, NULL, NULL, "eta at time infinity",                              0 },
  { "-v",         eslARG_NONE,   FALSE, NULL, NULL,   NULL, NULL, NULL, "be verbose",                                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile_in> <evohmmfile_out>";
static char banner[] = "test driver for evohmmer.c";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go  = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char          errbuf[eslERRBUFSIZE];
  char         *hmmfile;
  char         *evohmmfile;
  FILE         *hmmfp = NULL;                                                     /* HMM output file handle  */
  ESL_ALPHABET *abc = NULL;                                                       /* digital alphabet        */
  P7_HMMFILE   *hfp = NULL;                                                       /* open input HMM file     */
  P7_HMM       *hmm = NULL;
  double        etaz,;
  double        etainf;
  double        tol = 0.0001;
  int           N = 100;
  int           status = eslOK;
  int           hstatus = eslOK;
  int           verbose;

  if (esl_opt_ArgNumber(go) != 2)                    { puts("Incorrect number of command line arguments");          exit(1); }
  if ((hmmfile    = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <hmmfile> argument on command line");    exit(1); }
  if ((evohmmfile = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <evohmmfile> argument on command line"); exit(1); }

 /* Open the query profile HMM file */
  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);  

 /* Open evoHMM output file */
  hmmfp = fopen(evohmmfile, "w");
  if (hmmfp == NULL) p7_Fail("Failed to open HMM file %s for writing", evohmmfile);

  /* Options */
  etaz    = esl_opt_GetReal(go, "--etaz");
  etainf  = esl_opt_GetReal(go, "--etainfty");
  verbose = esl_opt_GetBoolean(go, "-v");
  fprintf(stdout, "# query HMM file:           %s\n", hmmfile);
  fprintf(stdout, "# eta at time zero:       = %g\n", etaz);
  fprintf(stdout, "# eta at time infinity:   = %g\n", etainf);

  /* read the HMM */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);

  utest_hmmevolve(hmmfp, hmm, etaz, etainf, N, tol, errbuf, verbose);

  esl_getopts_Destroy(go);
  p7_hmmfile_Close(hfp);
  if (hmmfp) fclose(hmmfp);
  return 0;
}
#endif /*evoHMMER_TESTDRIVE*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
