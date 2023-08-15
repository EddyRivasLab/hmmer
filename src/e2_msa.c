/*  e2_msa 
 *
 * ER, Sat Jul  6 09:11:17 EDT 2013 [Janelia] 
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
#include "e2_train.h"
#include "e2_tree.h"
#include "e2tracealign.h"
#include "evohmmer.h"
#include "msatree.h"
#include "minimize.h"

static void   e2_bracket_define_direction   (double *u, long np, struct e2_data *data);
static void   e2hmm_bracket_define_direction(double *u, long np, struct e2_data *data);
static void   e2_gd_define_stepsize         (double *u, long np, struct e2_data *data);
static void   e2hmm_gd_define_stepsize      (double *u, long np, struct e2_data *data);
static void   e2_pack_paramvector           (double *p, long np, struct e2_data *data);
static void   e2_unpack_paramvector         (double *p, long np, struct e2_data *data);
static void   e2hmm_pack_paramvector        (double *p, long np, struct e2_data *data);
static void   e2hmm_unpack_paramvector      (double *p, long np, struct e2_data *data);
static double e2_func                       (double *p, long np, void *dptr);

int
e2_msa(ESL_RANDOMNESS *r, E1_RATE *R, P7_RATE *R7, int n, ESL_SQ **seq, ESL_MSA *msa, float *msafrq,
       ESL_TREE *T, ESL_MSA **ret_msa, float *ret_sc, E2_PIPELINE *pli, 
       E1_BG *bg, P7_BG *bg7, E2_ALI e2ali, E2_OPT e2optimize, int mode, int do_viterbi, double tol, char *errbuf, int verbose)
{
  E1_RATE     *Rv   = NULL;	                            /* node index stack */
  P7_RATE     *R7v  = NULL;	                            /* node index stack */
  ESL_STACK   *vs   = NULL;	                            /* node index stack */
  E2_TRACE   **tr   = NULL;
  ESL_SQ      *sq   = NULL;
  PSQ        **sqn  = NULL;                                  /* a profile sequence for each internal node */
  PSQ         *sql;                                          /* convenience pointer to sq in left branch */
  PSQ         *sqr;                                          /* convenience pointer to sq in left branch */
  float        sc;
  float        totsc = 0.;
  int          nnodes;
  int          v;
  int          which;
  int          status;

  if (e2ali == E2NONE)  ESL_XFAIL(eslFAIL, errbuf, "undefined alignment method");
  if (T     == NULL)    ESL_XFAIL(eslFAIL, errbuf, "missing tree");

  /* allocate the profile sequence for internal nodes */
  nnodes = (T->N > 1)? T->N-1 : T->N;
  ESL_ALLOC(tr,  sizeof(E2_TRACE *) * nnodes);
  ESL_ALLOC(sqn, sizeof(PSQ      *) * (nnodes+T->N));
  for (v = 0; v < nnodes;      v ++) tr[v]  = NULL;
  for (v = 0; v < nnodes+T->N; v ++) sqn[v] = NULL;
 
  if (e2ali == E2 || e2ali == E2F) Rv  = e1_rate_Create(R->em->abc_r, R->evomodel);
  else                             R7v = p7_RateCreate(R7->M, R7->abc_r, R7->evomodel, -1.0, -1.0);       

  /* make a copy of the rate, but use the last one optimized, 
   * when you move to the next node
   */
  if (R)  e1_rate_Copy(R,  Rv);
  if (R7) p7_RateCopy(R7, R7v);

  /* PostOrder trasversal */
  if ((vs = esl_stack_ICreate())   == NULL) { status = eslFAIL; goto ERROR; }
  if (esl_stack_IPush(vs, nnodes-1) != eslOK) { status = eslFAIL; goto ERROR; }
  while (esl_stack_IPop(vs, &v) == eslOK)
    {
      if (T->left[v] <= 0) { /* dealign seq and convert to a psq */
	which = -T->left[v];

	if (msa) status = esl_sq_FetchFromMSA(msa, which, &sq); /* extract the seqs from the msa */
	else sq = seq[which];
	if (status != eslOK) { printf("esl_sq_FetchFromMSA() failed\n");  goto ERROR; }	
  
	switch(e2ali) {
	case E2:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, sq->desc, sq->acc, sq->abc, sq->dsq, sq->n);	
	  break;
	case E2HMMER:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, sq->desc, sq->acc, sq->abc, sq->dsq, sq->n);	
	  break;
	case E2F:
	  if (msa)
	    sqn[nnodes+which] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, msa->ax[which], msa->alen);
	  break;
	case E2FHMMER:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, sq->desc, sq->acc, sq->abc, sq->dsq, sq->n);
	  break;
	default:
	  status = eslFAIL; printf("unknown alicase\n"); goto ERROR;
	}	
	sql = sqn[nnodes+which];
	if (msa) { esl_sq_Destroy(sq); sq = NULL; }
     }
      else sql = sqn[T->left[v]];
        
      if (T->right[v] <= 0) { /* dealign seq and convert to a psq */
	which = -T->right[v];
	if (msa) status = esl_sq_FetchFromMSA(msa, which, &sq); /* extract the seqs from the msa */
	else sq = seq[which];
	if (status != eslOK) { printf("esl_sq_FetchFromMSA() failed\n");  goto ERROR; }
	
	switch(e2ali) {
	case E2:
	if (e2ali == E2 || e2ali == E2HMMER) 
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, sq->desc, sq->acc, sq->abc, sq->dsq, sq->n);
	  break;
	case E2HMMER:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, sq->desc, sq->acc, sq->abc, sq->dsq, sq->n);
	  break;
	case E2F:
	  if (msa) sqn[nnodes+which] = psq_CreateFrom(sq->name, msa->desc, msa->acc, msa->abc, msa->ax[which], msa->alen);
	  break;
	case E2FHMMER:
	  sqn[nnodes+which] = psq_CreateFrom(sq->name, sq->desc, sq->acc, sq->abc, sq->dsq, sq->n);
	  break;
	default:
	  status = eslFAIL; printf("unknown alicase\n"); goto ERROR;
	}
	
	sqr = sqn[nnodes+which];
	if (msa) { esl_sq_Destroy(sq); sq = NULL; }
      }
      else sqr = sqn[T->right[v]];
      
      if (sql != NULL && sqr != NULL) { /* ready to go: find ancestral profile sq running the e2 algorithm */
	if (verbose) 
	  printf("\nNODE %d parent %d | l:%s %d (%f,len=%d) r:%s %d (%f,len=%d)\n", 
		 v, T->parent[v], sql->name, T->left[v], T->ld[v], (int)sql->n, sqr->name, T->right[v], T->rd[v], (int)sqr->n);
	
	/* do the optimization for (sql,sqr) */
	if (R)  e1_rate_Copy(R,  Rv);
	if (R7) p7_RateCopy(R7, R7v);

	status = e2_Optimize(r, pli, sql, sqr, msafrq, Rv, R7v, bg, bg7, &(T->ld[v]), &(T->rd[v]), &(sqn[v]), &(tr[v]),
			     &sc, e2ali, e2optimize, mode, do_viterbi, tol, errbuf, verbose);
	if (status != eslOK) goto ERROR;
	
	/* push parent into stack unless already at the root */
	if (v > 0 && esl_stack_IPush(vs, T->parent[v]) != eslOK) { status = eslFAIL; goto ERROR; }; 
      }
      else if (sql == NULL) { /* not ready: push left child  into stack */	
	if (esl_stack_IPush(vs, T->left[v])   != eslOK) { status = eslFAIL; goto ERROR; };
      }
      else if (sqr == NULL) { /* not ready: push right child into stack */	
  	if (esl_stack_IPush(vs, T->right[v])  != eslOK) { status = eslFAIL; goto ERROR; };
      }
      
      totsc += sc;
    }

  /* new msa */
  if (e2_Tracealign((msa)?msa->name:seq[0]->name, T, sqn, tr, ret_msa, e2ali, errbuf, verbose) != eslOK) { status = eslFAIL; goto ERROR; };
  if (ret_sc) *ret_sc = totsc;
  
  /* clean up */
  if (tr) {
    for (v = 0; v < nnodes;      v ++) e2_trace_Destroy(tr[v]);
    free(tr);
  }
  if (sqn) {
    for (v = 0; v < nnodes+T->N; v ++) psq_Destroy(sqn[v]);
    free(sqn);
  }
  esl_stack_Destroy(vs);
  if (Rv) e1_rate_Destroy(Rv);
  if (R7v) p7_RateDestroy(R7v);
  return eslOK;
  
 ERROR:
  if (tr) {
    for (v = 0; v < nnodes;      v ++) e2_trace_Destroy(tr[v]);
    free(tr);
  }
  if (sqn) {
    for (v = 0; v < nnodes+T->N; v ++) psq_Destroy(sqn[v]);
    free(sqn);
  }
  if (vs) esl_stack_Destroy(vs);
  if (Rv) e1_rate_Destroy(Rv); Rv = NULL;
  if (R7v) p7_RateDestroy(R7v); R7v = NULL;
  if (T) esl_tree_Destroy(T); T = NULL;
  return status;
}

int
e2_Optimize(ESL_RANDOMNESS *r, E2_PIPELINE *pli, PSQ *sql, PSQ *sqr, float *frq, E1_RATE *R, P7_RATE *R7, E1_BG *bg, P7_BG *bg7,
	    double *ret_timel, double *ret_timer, PSQ **ret_sqa, E2_TRACE **ret_tr, float *ret_sc, E2_ALI e2ali, E2_OPT e2optimize,  
	    int mode, int do_viterbi, double tol, char *errbuf, int verbose)
{
  struct e2_data   data;
  double          *p;	               /* parameter vector                  */
  double          *u;                  /* max initial step size vector      */
  double          *wrk; 	       /* 4 tmp vectors of length nbranches */
  double           fx;
  double           firststep;
  int              nvariables;
  int              status;
  

  /* Copy shared info into the "data" structure
   */
  data.r           = r;
  data.sql         = sql;
  data.sqr         = sqr;
  data.frq         = frq;
  data.timel       = *ret_timel;
  data.timer       = *ret_timer;
  data.pli         = pli;
  data.e2ali       = e2ali;
  data.e2opt       = e2optimize;
  data.mode        = mode;
  data.do_viterbi  = do_viterbi;
  data.R           = R;
  data.R7          = R7;
  data.bg          = bg;
  data.bg7         = bg7;
  data.errbuf      = errbuf;
  data.it          = 0;
  data.tol         = tol;
  data.errbuf      = errbuf;
  data.verbose     = verbose;

  if (data.e2opt == OPTNONE) 
    return e2_Pipeline(data.r, data.pli, data.sql, data.sqr, data.frq, data.R, data.R7, data.bg, data.bg7, data.timel, data.timer,
		       ret_sqa, ret_tr, ret_sc, NULL, data.e2ali, data.mode, TRUE, data.do_viterbi, data.tol, data.errbuf, data.verbose);

  /* optimize */
  if (e2ali == E2 || e2ali == E2F)  
    {
      if      (data.e2opt == OPTTIME) { nvariables = 1;                         } // timel=timer 
      else if (data.e2opt == OPTPARA) { nvariables = R->nrate + R->nbern;       } // params,
      else if (data.e2opt == OPTBOTH) { nvariables = R->nrate + R->nbern + 1;   } // params, timel=timer
      else                            { status = eslFAIL; goto ERROR; }     
    }
  else  
    {
      nvariables = 2; /* timel, timer */
    }

  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (nvariables+1));
  ESL_ALLOC(u,   sizeof(double) * (nvariables+1));
  ESL_ALLOC(wrk, sizeof(double) * (nvariables+1) * 4);

  if (data.e2ali == E2 || data.e2ali == E2F) e2_pack_paramvector   (p, (long)nvariables, &data);
  else                                       e2hmm_pack_paramvector(p, (long)nvariables, &data);

  /* pass problem to the optimizer
   */
#if 0
  // Define the step size vector u.
  if (e2ali == E2 || e2ali == E2F)  e2_gd_define_stepsize   (u, nvariables, &data);
  else                              e2hmm_gd_define_stepsize(u, nvariables, &data);


  if ((status = esl_min_ConjugateGradientDescent(p, u, nvariables, 
						 &e2_func,
						 NULL, 
						 (void *) (&data), 
						 tol, wrk, &fx))  != eslOK)
    esl_fatal("e2_msa(): bad conjugate gradient descent");
#else
  // Define the optimization direction
  if (e2ali == E2 || e2ali == E2F)  e2_bracket_define_direction   (u, nvariables, &data);
  else                              e2hmm_bracket_define_direction(u, nvariables, &data);
  firststep = 1e-0;

  if ((status = min_Bracket(p, u, (long)nvariables, firststep,
			    &e2_func,
			    (void *) (&data), 
			    tol, wrk, &fx))  != eslOK)
    esl_fatal("e2_msa(): bad bracket minimization");
#endif
  
  if (data.e2ali == E2 || data.e2ali == E2F) e2_unpack_paramvector(p, (double)nvariables, &data);
  else                                       e2hmm_unpack_paramvector(p, (double)nvariables, &data);

  if (verbose) {
    if (data.e2ali == E2 || data.e2ali == E2F) {
    printf("  optimal [%d]: rI %f rM %f rD %f | ldE %f muE %f ldI %f muI %f muA %f| time %f | sc %f | model %s\n\n",  
	   data.it, data.R->rI, data.R->rM, data.R->rD, 
	   data.R->ldE[e1R_S], data.R->muE[e1R_S], data.R->ldE[e1R_I], data.R->muE[e1R_I], data.R->muA[e1R_S], 
	   data.timel, fx, e1_rate_EvomodelType(data.R->evomodel));
   }
    else {
     printf(" optimal [%d]: timel %f timer %f | sc %f \n", 
	    data.it, data.timel, data.timer, fx);
    }
  }

  *ret_timel =  data.timel;
  *ret_timer =  data.timer;
  *ret_sc    = -fx;

  /* clean up */
  free(u);
  free(p);
  free(wrk);
  return eslOK;
  
 ERROR:
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  if (wrk != NULL) free(wrk);
  return status;  
}


/* -- internal functions -- */

static void
e2_bracket_define_direction(double *u, long np, struct e2_data *data)
{
  int x = 0;
  int y;

   if (data->e2opt == OPTTIME || 
       data->e2opt == OPTBOTH   ) 
    {
       u[x++] = 0.25;
    }
  for (y = x; y < np; y++) u[y] = 0.25;
  u[np] = 0.25;
}

static void
e2hmm_bracket_define_direction(double *u, long np, struct e2_data *data)
{
  int x = 0;
  int y;

  if (data->e2opt == OPTTIME || 
      data->e2opt == OPTBOTH   ) 
    {
      u[x++] = 0.1;
    }
  for (y = x; y < np; y++) u[y] = ((y-x)< data->R->nrate)? 0.1 : 0.1;
  u[np] = 10.;
}


static void
e2_gd_define_stepsize(double *u, long np, struct e2_data *data)
{
  int x = 0;
  int y;

  if (data->e2opt == OPTTIME || 
      data->e2opt == OPTBOTH   ) 
    {
      u[x++] = 0.1;
    }
  for (y = x; y < np; y++) u[y] = ((y-x)< data->R->nrate)? 0.1 : 100.1;
  u[np] = 0.1;
}

static void
e2hmm_gd_define_stepsize(double *u, long np, struct e2_data *data)
{
  int x = 0;
  int y;

  if (data->e2opt == OPTTIME || 
      data->e2opt == OPTBOTH   ) 
    {
      u[x++] = 0.1;
    }
  for (y = x; y < np; y++) u[y] = ((y-x)< data->R->nrate)? 0.1 : 100.1;
  u[np] = 0.1;
}

static void
e2_pack_paramvector(double *p, long np, struct e2_data *data)
{
  double             time;
  int                x = 0; 
  int                status;

  if (data->e2opt == OPTTIME || data->e2opt == OPTBOTH) {
    time = 0.5 * (data->timel + data->timer);
    p[x++] = time;
  }
 
  if (data->e2opt == OPTPARA || data->e2opt == OPTBOTH) {
    status = e2_transitions_pack_paramvector(&x, p, np, data->R, data->errbuf, data->verbose);
    if (status != eslOK) { printf("error at e2_pack_paramvector()\n%s\n", data->errbuf); exit(1); }
  }    
} 

static void
e2hmm_pack_paramvector(double *p, long np, struct e2_data *data)
{
  int x = 0;

  p[x++] = data->timel;
  p[x++] = data->timer;
} 

static void 
e2_unpack_paramvector(double *p, long np, struct e2_data *data) 
{ 
  double time;
  double tmax = 1.9;
  int    x = 0;
  int    status;
  
  if (data->e2opt == OPTTIME || data->e2opt == OPTBOTH) {
    time = fabs(p[x++]);
    if (time > tmax) time = tmax;
    data->timel = data->timer = time;
  }
  
  if (data->e2opt == OPTPARA || data->e2opt == OPTBOTH) {
    status = e2_transitions_unpack_paramvector(&x, p, np, data->R, data->errbuf, data->verbose);
    if (status != eslOK) { printf("error at e2_unpack_paramvector()\n"); exit(1); }
  }
}

static void 
e2hmm_unpack_paramvector(double *p, long np, struct e2_data *data) 
{ 
  float  tmax = 500.0;
  int    x = 0;
  
  data->timel = fabs(p[x++]); 
  data->timer = fabs(p[x++]); 

  if (data->timel > tmax) data->timel = tmax;
  if (data->timer > tmax) data->timer = tmax;
}

static double 
e2_func(double *p, long np, void *dptr)
{
  struct e2_data *data = (struct e2_data *) dptr;
  int    status;

  if (data->e2ali == E2 || data->e2ali == E2F) e2_unpack_paramvector(p, (double)np, data);
  else                                         e2hmm_unpack_paramvector(p, (double)np, data);

  data->it ++;
  status = e2_Pipeline(data->r, data->pli, data->sql, data->sqr, data->frq, data->R, 
		       data->R7, data->bg, data->bg7, data->timel, data->timer, 
		       NULL, NULL, &(data->sc), &(data->accsc), data->e2ali, data->mode, FALSE, data->do_viterbi, data->tol, data->errbuf, data->verbose);
  if (status != eslOK) { printf("error at e2_func()\n%s\n", data->errbuf); exit(1); }
  if (data->sc <= -eslINFINITY) return 20000.;  

#if 0
  if (data->e2ali == E2 || data->e2ali == E2F)
    printf("  it %d: rI %f rM %f rD %f | ldE %f muE %f ldI %f muI %f muA %f | timel %f timer %f | e2_func %f\n",  
	   data->it, data->R->rI, data->R->rM, data->R->rD, 
	   data->R->ldE[e1R_S], data->R->muE[e1R_S], data->R->ldE[e1R_I], data->R->muE[e1R_I], data->R->muA[e1R_S], 
	   data->timel, data->timer, data->sc);
  else 
    printf(" it: %d timel %f timer %f e2_func %f\n\n", data->it, data->timel, data->timer, data->sc);
#endif
  
  return -data->sc;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
