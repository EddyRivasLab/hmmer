/* E2's pipeline
 *  
 * Contents:
 *   1. E2_PIPELINE: allocation, initialization, destruction
 *   2. Pipeline API
 *   3. Example 1: search mode (in a sequence db)
 *   4. Example 2: scan mode (in an HMM db)
 *   5. Copyright and license information
 * 
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

#include "easel.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_vectorops.h"


#include "hmmer.h"

#include "e2.h"
#include "e1_model.h"
#include "e2_generic_fwdback.h"
#include "e2_generic_decoding.h"
#include "e2_generic_optacc.h"
#include "e2_generic_viterbi.h"
#include "e2_generic_vtrace.h"
#include "e2f_generic_fwdback.h"
#include "e2f_generic_decoding.h"
#include "e2f_generic_optacc.h"
#include "e2fhmmer_generic_fwdback.h"
#include "e2fhmmer_generic_decoding.h"
#include "e2fhmmer_generic_optacc.h"
#include "modelconfig.h"
#include "e2_profile.h"
#include "e2hmmer_profile.h"
#include "e2_pipeline.h"
#include "e2_msa.h"

static int hmm_align(P7_HMM *evo7l, P7_HMM *evo7r, PSQ **ret_sql, PSQ **ret_sqr, char *errbuf, int verbose);

/*****************************************************************
 * 1. The E2_PIPELINE object: allocation, initialization, destruction.
 *****************************************************************/
/* Function:  e2_pipeline_Create()
 * Synopsis:  Create a new E2 pipeline.
 *
 * Purpose:   
 *            
 * Returns:   ptr to new <E2_PIPELINE> object on success. Caller frees this
 *            with <e2_pipeline_Destroy()>.
 *
 * Throws:    <NULL> on allocation failure.
 */
E2_PIPELINE *
e2_pipeline_Create(ESL_GETOPTS *go, int L1_hint, int L2_hint, int M)
{
  E2_PIPELINE *pli  = NULL;
  int          status;

  ESL_ALLOC(pli, sizeof(E2_PIPELINE));

  if ((pli->gx1 = e2_gmx_Create(M, L1_hint, L2_hint)) == NULL) goto ERROR;
  if ((pli->gx2 = e2_gmx_Create(M, L1_hint, L2_hint)) == NULL) goto ERROR;
 
  return pli;

 ERROR:
  e2_pipeline_Destroy(pli);
  return NULL;

}

/* Function:  e2_pipeline_Reuse()
 * Synopsis:  Reuse a pipeline for next target.
 *
 * Purpose:   Reuse <pli> for next target sequence (search mode)
 *            or model (scan mode). 
 *            
 *            May eventually need to distinguish from reusing pipeline
 *            for next query, but we're not really focused on multiquery
 *            use of hmmscan/hmmsearch/phmmer for the moment.
 */
int
e2_pipeline_Reuse(E2_PIPELINE *pli)
{
  e2_gmx_Reuse(pli->gx1);
  e2_gmx_Reuse(pli->gx2);
 return eslOK;
}



/* Function:  e2_pipeline_Destroy()
 * Synopsis:  Free a <E2_PIPELINE> object.
 *
 * Purpose:   Free a <E2_PIPELINE> object.
 */
void
e2_pipeline_Destroy(E2_PIPELINE *pli)
{
  if (pli == NULL) return;  
  e2_gmx_Destroy(pli->gx1);
  e2_gmx_Destroy(pli->gx2);
  free(pli);
}
/*---------------- end, E2_PIPELINE object ----------------------*/


/*****************************************************************
 * 2. The pipeline API.
 *****************************************************************/
/* Function:  e2_Pipeline()
 * Synopsis:  E2's pipeline.
 *
 * Purpose:   Run E2's pipeline
 *            
 * Returns:   <eslOK> on success.
 *            
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      
 */
int
e2_Pipeline(ESL_RANDOMNESS *r, E2_PIPELINE *pli, const PSQ *psql, const PSQ *psqr, const float *frq, E1_RATE *R, P7_RATE *R7, E1_BG *bg, P7_BG *bg7, 
	    float timel, float timer, PSQ **ret_sqa, E2_TRACE **ret_tr, float *ret_sc, 
	    float *ret_accsc, E2_ALI e2ali, int mode, int decode, int do_viterbi, double tol, char *errbuf, int verbose)
{
  E1_MODEL        *evol  = NULL;
  E1_MODEL        *evor  = NULL;
  P7_HMM          *evo7l = NULL;
  P7_HMM          *evo7r = NULL;
  E2_TRACE        *tr    = NULL;
  E2_PROFILE      *gm    = NULL;
  E2HMMER_PROFILE *gm7   = NULL;
  ESL_ALPHABET    *abc   = (ESL_ALPHABET *)psql->abc;
  PSQ             *sql   = NULL;
  PSQ             *sqr   = NULL;
  PSQ             *sqa   = NULL;
  float           *ancf;
  float           *insf;
  float            vsc;        // the viterbi score
  float            fsc, bsc;
  float            accscore;
  float            L;         /* average length of the two sequences */
  int              kmax;
  int              status;

  if (psql->abc->type != psqr->abc->type) {status = eslFAIL; goto ERROR; }

  /* Allocations */
  sql = psq_Clone(psql);
  sqr = psq_Clone(psqr);
  tr  = e2_trace_CreateWithPP();
  p7_FLogsumInit();
  
  /* Evolve the models */
  L = 0.5 * (sql->n + sqr->n); 
  if (e2ali == E2 || e2ali == E2F) {
    
    ancf = (float *)frq; // residue distribution for ancestral sequence
    insf = (float *)frq; // residue distribution for insertions
    
    evol = e1_model_Create(R, (timel==0.0)? 1e-5:timel, ancf, insf, mode, L, abc, tol, errbuf, FALSE); if (evol == NULL) {status = eslFAIL; goto ERROR; }
    evor = e1_model_Create(R, (timer==0.0)? 1e-5:timer, ancf, insf, mode, L, abc, tol, errbuf, FALSE); if (evor == NULL) {status = eslFAIL; goto ERROR; } 

    if (verbose) {
      e1_model_DumpTransitions(stdout, evol);
      e1_model_DumpTransitions(stdout, evor);
    }
    
    /* Configure a profile per HMM state*/
    gm = e2_profile_Create(abc);
    e2_ProfileConfig(evol, evor, ancf, gm, L, e2ali, verbose);
    
    /* grow matrices if needeed */
    e2_gmx_GrowTo(pli->gx1, 1, sql->n, sqr->n);/* expand matrices if needed */
    e2_gmx_GrowTo(pli->gx2, 1, sql->n, sqr->n);/* expand matrices if needed */
  }
  else  {
    /* build evohmms */
    evo7l = p7_hmm_Create(R7->M, abc);
    evo7r = p7_hmm_Create(R7->M, abc);
    
    status = p7_EvolveFromRate(stdout, evo7l, R7, bg7, timel, 0.0001, errbuf, verbose); if (status != eslOK) goto ERROR;
    status = p7_EvolveFromRate(stdout, evo7r, R7, bg7, timer, 0.0001, errbuf, verbose); if (status != eslOK) goto ERROR;
    
    if (e2ali == E2FHMMER) {
      /* Align each sequence to the corresponding hmm.
       * Use the two aligned sequences to proceed with the e2fhmmer algorithm 
       */
      status = hmm_align(evo7l, evo7r, &sql, &sqr, errbuf, verbose); if (status != eslOK) goto ERROR;
    }
    
    /* Configure a profile per HMM state */
    gm7 = e2hmmer_profile_Create(R7->M, abc);
    status = e2hmmer_ProfileConfig(R7, timel, timer, evo7l, evo7r, bg7, gm7, L, e2ali, e2_GLOBAL, verbose); if (status != eslOK) goto ERROR;
    
    /* grow matrices if needeed */
    e2_gmx_GrowTo(pli->gx1, R7->M, sql->n, sqr->n);/* expand matrices if needed */
    e2_gmx_GrowTo(pli->gx2, R7->M, sql->n, sqr->n);/* expand matrices if needed */
   }


  /* Run Forward, Backward; do OA fill and trace */
    switch(e2ali) {
    case E2:
      if (do_viterbi) {
	status = e2_GViterbi(sql, sqr, gm, pli->gx1, &vsc);
	if (status != eslOK) goto ERROR;
      }
      else {
	status = e2_GForward (sql, sqr, gm, pli->gx1, &fsc);
	if (status != eslOK) goto ERROR;
	status = e2_GBackward(sql, sqr, gm, pli->gx2, &bsc);
	if (status != eslOK) goto ERROR;
	status = e2_GDecoding(gm, pli->gx1, pli->gx2, pli->gx2);                                    /* <gx2> is now the posterior decoding matrix */
      }
      if (status != eslOK) goto ERROR;
      if (decode) {
	if (do_viterbi) {
	  status = e2_GTrace(sql, sqr, gm, pli->gx1, tr);
	  if (status != eslOK) goto ERROR;
	}
	else {
	  status = e2_GOptimalAccuracy(gm, pli->gx2, pli->gx1, &accscore);	                  /* <gx1> is now the OA matrix */
	  if (status != eslOK) goto ERROR;
	  status = e2_GOATrace(r, evol, evor, ancf, gm, pli->gx2, pli->gx1, tr, sql, sqr, &sqa);  
	  if (status != eslOK) goto ERROR;
	}
      }
      break;
    case E2F:
      if (sql->n != sqr->n) { 
	printf("e2_pipeline(): e2fix=TRUE but sqs are not aligned. sqln %" PRId64 "  sqrn %" PRId64 " \n", sql->n, sqr->n); 
	return eslFAIL;
      }
      
      if (do_viterbi)  ESL_FAIL(eslFAIL, errbuf, "Viterbi not implemented");
      status = e2f_GForward (sql, sqr, gm, pli->gx1, &fsc); 
      if (status != eslOK) goto ERROR;
      status = e2f_GBackward(sql, sqr, gm, pli->gx2, &bsc);
      if (status != eslOK) goto ERROR;
      status = e2f_GDecoding(gm, pli->gx1, pli->gx2, pli->gx2, errbuf);                    /* <gx2> is now the posterior decoding matrix */
      if (status != eslOK) goto ERROR;
      if (decode) {
	status = e2f_GOptimalAccuracy(sql, sqr, gm, pli->gx2, pli->gx1, &accscore);	 /* <gx1> is now the OA matrix */
	if (status != eslOK) goto ERROR;
	status = e2f_GOATrace(r, evol, evor, ancf, gm, pli->gx2, pli->gx1, tr, sql, sqr, &sqa);  
	if (status != eslOK) goto ERROR;
      }
      break;
    case E2HMMER:
      printf("TODO: E2HMMER alignment method not implemented yet\n"); status = eslFAIL; goto ERROR;
      break;
    case E2FHMMER:
      if (sql->n != sqr->n) { printf("e2_pipeline(): e2fix=TRUE but sqs are not aligned.\n"); return eslFAIL; }
      if (do_viterbi)  ESL_FAIL(eslFAIL, errbuf, "Viterbi not implemented");
      
      status = e2fhmmer_GForward (sql, sqr, gm7, pli->gx1, &fsc); 
      if (status != eslOK) goto ERROR;
      status = e2fhmmer_GBackward(sql, sqr, gm7, pli->gx2, &bsc);
      if (status != eslOK) goto ERROR;
      status = e2fhmmer_GDecoding(gm7, pli->gx1, pli->gx2, fsc, pli->gx2);                         /* <gx2> is now the posterior decoding matrix */
      if (status != eslOK) goto ERROR;
      if (decode) {
	status = e2fhmmer_GOptimalAccuracy(gm7, pli->gx2, pli->gx1, sql, sqr, &accscore, &kmax);   /* <gx1> is now the OA matrix */
	if (status != eslOK) goto ERROR;
	status = e2fhmmer_GOATrace(r, R7, timel, timer, bg7, gm7, kmax, pli->gx2, pli->gx1, tr, sql, sqr, &sqa);  
	if (status != eslOK) goto ERROR;
      }
      break;
    default:
      ESL_FAIL(eslFAIL, errbuf, "unknown alignment method");
      break;
    }
    if (verbose) {
      if (do_viterbi)  printf("vit = %.4f nats\n", vsc);
      else           { printf("fwd = %.4f nats\n", fsc);
                       printf("bck = %.4f nats\n", bsc);}
    }
#if 0
    if (fabs(fsc-bsc) > 0.01*0.5*(sql->n+sqr->n)) ESL_XFAIL(eslFAIL, errbuf, "fwd=%f and bck=%f  do not agree", fsc, bsc);
#endif

    if (verbose && decode) {
      printf("ancestral sequence length %d\n", (int)sqa->n);
      printf("acc = %.4f (%.2f%%)\n", accscore, (sqa->n > 0)? accscore * 100. / (float) sqa->n : accscore * 100);
    }
    if (verbose) e2_trace_Dump(stdout, tr, gm, gm7, sql, sqr, sqa);
    
    if (evol) e1_model_Destroy(evol);
    if (evor) e1_model_Destroy(evor);
    if (gm)   e2_profile_Destroy(gm);
    if (sql)  psq_Destroy(sql);
    if (sqr)  psq_Destroy(sqr);
   
    if (evo7l) p7_hmm_Destroy(evo7l);
    if (evo7r) p7_hmm_Destroy(evo7r);
    if (gm7)   e2hmmer_profile_Destroy(gm7);
    
    if (ret_sqa)   *ret_sqa   = sqa; else psq_Destroy(sqa);
    if (ret_tr)    *ret_tr    = tr;  else e2_trace_Destroy(tr);
    if (ret_sc)    *ret_sc    = fsc;
    if (ret_accsc) *ret_accsc = accscore;
    
    return eslOK;
    
 ERROR:
    if (evol)  e1_model_Destroy(evol);
    if (evor)  e1_model_Destroy(evor);
    if (gm)    e2_profile_Destroy(gm);
    if (evo7l) p7_hmm_Destroy(evo7l);
    if (evo7r) p7_hmm_Destroy(evo7r);
    if (gm7)   e2hmmer_profile_Destroy(gm7);
    if (tr)    e2_trace_Destroy(tr);
    if (sqa)   psq_Destroy(sqa);
    return status;
}




/*------------------- end, pipeline API -------------------------*/


static int
hmm_align(P7_HMM *evo7l, P7_HMM *evo7r, PSQ **ret_sql, PSQ **ret_sqr, char *errbuf, int verbose)
{
  P7_TRACE     **tr      = NULL;	/* array of tracebacks             */
  PSQ           *sql     = *ret_sql;
  PSQ           *sqr     = *ret_sqr;
  PSQ           *newl    = NULL;
  PSQ           *newr    = NULL;
  char          *asql    = NULL;
  char          *asqr    = NULL;
  ESL_SQ       **sq      = NULL;	/* array of sequences              */
  ESL_MSA       *msa     = NULL;	/* resulting multiple alignment    */
  ESL_ALPHABET  *abc     =  (ESL_ALPHABET *)evo7l->abc;
  int            msaopts = 0;        	/* flags to p7_tracealign_Seqs()   */
  int            totseq  = 2;
  int            M       = evo7l->M;
  int            idx;
  int            status;

  msaopts |= p7_ALL_CONSENSUS_COLS; /* default as of 3.1 */

  ESL_ALLOC(tr, sizeof(P7_TRACE *) * totseq);
  for (idx = 0; idx < totseq; idx++)
    tr[idx] = p7_trace_CreateWithPP();

  psq_ConvertToAseq(sql, &asql);
  psq_ConvertToAseq(sqr, &asqr);
  if (0&&verbose) {
    printf("asql %s\n",asql);
    printf("asqr %s\n", asqr);
  }

  ESL_ALLOC(sq, sizeof(ESL_SQ *) * (totseq));
  sq[0] = esl_sq_CreateFrom(evo7l->name, asql, evo7l->desc, evo7l->acc, NULL);
  sq[1] = esl_sq_CreateFrom(evo7r->name, asqr, evo7r->desc, evo7r->acc, NULL);
  esl_sq_Digitize(abc, sq[0]);
  esl_sq_Digitize(abc, sq[1]);

  p7_tracealign_ComputeTraces(evo7l, sq, 0, 1, tr);
  p7_tracealign_ComputeTraces(evo7r, sq, 1, 1, tr);

  p7_tracealign_Seqs(sq, tr, totseq, M, msaopts, evo7l, &msa);
  esl_msa_Digitize(abc, msa, NULL);

  if (verbose) esl_msafile_Write(stdout, msa, eslMSAFILE_STOCKHOLM);
  
  newl = psq_CreateFrom(NULL, NULL, NULL, abc, msa->ax[0], msa->alen);	
  newr = psq_CreateFrom(NULL, NULL, NULL, abc, msa->ax[1], msa->alen);	
  
  if (0&&verbose) {
    psq_ConvertToAseq (newl, &asql);
    psq_ConvertToAseq (newr, &asqr);
    printf("asql %s\n", asql);
    printf("asqr %s\n", asqr);
    if (asql) free(asql); asql = NULL;
    if (asqr) free(asqr); asqr = NULL;
   }
  
  *ret_sql = newl; psq_Destroy(sql);
  *ret_sqr = newr; psq_Destroy(sqr);

  for (idx = 0; idx < totseq; idx++) esl_sq_Destroy(sq[idx]);    
  for (idx = 0; idx < totseq; idx++) p7_trace_Destroy(tr[idx]); 
  free(asql);
  free(asqr); 
  free(sq);
  free(tr); 
  esl_msa_Destroy(msa);

  return eslOK;

 ERROR:
  for (idx = 0; idx <= totseq; idx++) if (sq[idx]) esl_sq_Destroy(sq[idx]);    
  for (idx = 0; idx <  totseq; idx++) if (tr[idx]) p7_trace_Destroy(tr[idx]); 
  if (asql) free(asql);
  if (asqr) free(asqr);
  if (sq) free(sq);
  if (tr) free(tr);
  if (msa) esl_msa_Destroy(msa);
  return status;
}


/*****************************************************************
 * @LICENSE@
 *
 * SVN $URL: $
 * SVN $Id: $
 *****************************************************************/
