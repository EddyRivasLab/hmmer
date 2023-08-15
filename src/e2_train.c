/*  e2_train
 *
 * ER, Sat Mar 22 12:30:04 EDT 2014 [Janelia] 
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
#include "fastsp.h"
#include "msamanip.h"
#include "msatree.h"
#include "minimize.h"

static void   e2_train_bracket_define_direction(double *u, long np);
static void   e2_train_gd_define_stepsize      (double *u, long np, struct e2_train_data *data);
static void   e2_train_NM_define_stepsize      (double *u, long np, struct e2_train_data *data);
static void   e2_train_pack_paramvector        (double *p, long np, struct e2_train_data *data);
static void   e2_train_unpack_paramvector      (double *p, long np, struct e2_train_data *data);
static double e2_train_func                    (double *p, long np, void *dptr);

int
e2_train(ESL_RANDOMNESS *r, int nmsa, ESL_TREE **Tlist, ESL_MSA **msalist, float **msafrq, E1_RATE *R, E2_PIPELINE *pli, E1_BG *bg, int mode, int do_viterbi, int amoeba, 
	 double tol, char *errbuf, int verbose)
{
  struct e2_train_data   data;
  double                *p;	         /* parameter vector                  */
  double                *u;              /* max initial step size vector      */
  double                *wrk; 	         /* 4 tmp vectors of length nbranches */
  ESL_DMATRIX           *wrks = NULL;
  double                 firststep;
  double                 fx;
  int                    nvariables;
  int                    status;
  
  /* Copy shared info into the "data" structure
   */
  data.r           = r;
  data.pli         = pli;
  data.e2ali       = E2F;  // training: alignment is fixed
  data.mode        = mode;
  data.do_viterbi  = do_viterbi;
  data.R           = R;
  data.bg          = bg;
  data.nmsa        = nmsa;
  data.Tlist       = Tlist;
  data.msalist     = msalist;
  data.msafrq      = msafrq;
  data.it          = 0;
  data.errbuf      = errbuf;
  data.tol         = tol;
  data.errbuf      = errbuf;
  data.verbose     = verbose;
  
  nvariables = R->nrate + R->nbern;
  
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (nvariables+1));
  ESL_ALLOC(u,   sizeof(double) * (nvariables+1));
  ESL_ALLOC(wrk, sizeof(double) * (nvariables+1) * 4);
  
  e2_train_pack_paramvector(p, (long)nvariables, &data);
  
  /* pass problem to the optimizer
   */
 #if 0
  e2_train_bracket_define_direction(u, nvariables); // define the optimization direction
  firststep = 1e-0;
  
  if ((status = min_Bracket(p, u, (long)nvariables, firststep,
			    &e2_train_func,
			    (void *) (&data), 
			    tol, wrk, &fx))  != eslOK)
    esl_fatal("e2_train(): bad bracket minimization");
#endif

  if (amoeba) {
    wrks = esl_dmatrix_Create(nvariables+1, nvariables);
    e2_train_gd_define_stepsize(u, nvariables, &data); // Define the step size vector u.
    
    if ((status = min_NelderMead(p, u, (long)nvariables, 
				 &e2_train_func,
				 (void *) (&data), 
				 tol, wrk, wrks, &fx))  != eslOK)
      esl_fatal("e2_train(): bad NelderMead optimiziation");
    
  }
  else {
    e2_train_NM_define_stepsize(u, nvariables, &data); // Define the step size vector u.
    
    if ((status = esl_min_ConjugateGradientDescent(p, u, nvariables, 
						   &e2_train_func,
						   NULL, 
						   (void *) (&data), 
						   tol, wrk, &fx))  != eslOK)
      esl_fatal("e2_train(): bad conjugate gradient descent");
  }
  
  e2_train_unpack_paramvector(p, (double)nvariables, &data);
  
  if (1||verbose) {
    printf("OPT [%d]: rI %f rM %f rD %f | sI %f | ldEM %f muEM %f ldED %f muED %f ldI %f muI %f muAM %f muAD %f muAI %f | sc %f | model %s\n\n",  
	   data.it, data.R->rI, data.R->rM, data.R->rD, data.R->sI, 
	   data.R->ldE[e1R_S], data.R->muE[e1R_S], data.R->ldE[e1R_D], data.R->muE[e1R_D], 
	   data.R->ldE[e1R_I], data.R->muE[e1R_I], 
	   data.R->muA[e1R_S], data.R->muA[e1R_D], data.R->muA[e1R_I], 
	   fx, e1_rate_EvomodelType(data.R->evomodel));
  }
  
  /* clean up */
  free(u);
  free(p);
  free(wrk);
  if (wrks) esl_dmatrix_Destroy(wrks);
  return eslOK;
  
 ERROR:
  if (p)   free(p);
  if (u)   free(u);
  if (wrk) free(wrk);
  if (wrks) esl_dmatrix_Destroy(wrks);
  return status;  
}

int 
e2_transitions_pack_paramvector(int *ret_x, double *p, long np, E1_RATE *R, char *errbuf, int verbose)
{
  struct rateparam_s rateparam;
  int                rparam;
  double             inc     = 1e-1;
  double             ldmin   = 0.01;
  double             ldmax   = 0.69;
  double             mumin   = 0.01;
  double             mumax   = 0.7;
  double             rMmin   = 0.80;
  double             rMmax   = 0.9999;
  double             rDImin  = 0.30;
  double             rDImax  = 0.95;
  double             sImin  = 0.30;
  double             sImax  = 0.95;
  int                x       = *ret_x; 
  int                status  = eslOK;

  rateparam.muAM = R->muA[e1R_S];
  rateparam.muAD = R->muA[e1R_D];
  rateparam.muAI = R->muA[e1R_I];
  rateparam.ldEM = R->ldE[e1R_S];
  rateparam.muEM = R->muE[e1R_S];
  rateparam.ldED = R->ldE[e1R_D];
  rateparam.muED = R->muE[e1R_D];
  rateparam.ldI  = R->ldE[e1R_I];
  rateparam.muI  = R->muE[e1R_I];
  rateparam.rI   = R->rI;
  rateparam.rM   = R->rM;
  rateparam.rD   = R->rD;
  rateparam.sI   = R->sI;
  rateparam.sD   = R->sD;
  rateparam.vI   = R->vI;
  rateparam.vD   = R->vD;

  if (rateparam.ldI  > -1.0 && rateparam.ldI  < ldmin)  rateparam.ldI  = ldmin  + inc;
  if (rateparam.ldEM > -1.0 && rateparam.ldEM < ldmin)  rateparam.ldEM = ldmin  + inc;
  if (rateparam.ldED > -1.0 && rateparam.ldED < ldmin)  rateparam.ldED = ldmin  + inc;
  if (rateparam.muI  > -1.0 && rateparam.muI  < mumin)  rateparam.muI  = mumin  + inc;
  if (rateparam.muEM > -1.0 && rateparam.muEM < mumin)  rateparam.muEM = mumin  + inc;
  if (rateparam.muED > -1.0 && rateparam.muED < mumin)  rateparam.muED = mumin  + inc;
  if (rateparam.muAM > -1.0 && rateparam.muAM < mumin)  rateparam.muAM = mumin  + inc;
  if (rateparam.muAD > -1.0 && rateparam.muAD < mumin)  rateparam.muAD = mumin  + inc;
  if (rateparam.muAI > -1.0 && rateparam.muAI < mumin)  rateparam.muAI = mumin  + inc;
  if (rateparam.rM   > -1.0 && rateparam.rM   < rMmin)  rateparam.rM   = rMmin  + inc;
  if (rateparam.rD   > -1.0 && rateparam.rD   < rDImin) rateparam.rD   = rDImin + inc;
  if (rateparam.rI   > -1.0 && rateparam.rI   < rDImin) rateparam.rI   = rDImin + inc;
  if (rateparam.sI   > -1.0 && rateparam.sI   < sImin)  rateparam.sI   = sImin  + inc;

  if (rateparam.ldI  > ldmax)  rateparam.ldI  = ldmax  - inc;
  if (rateparam.ldEM > ldmax)  rateparam.ldEM = ldmax  - inc;
  if (rateparam.ldED > ldmax)  rateparam.ldED = ldmax  - inc;
  if (rateparam.muI  > mumax)  rateparam.muI  = mumax  - inc;
  if (rateparam.muEM > mumax)  rateparam.muEM = mumax  - inc;
  if (rateparam.muED > mumax)  rateparam.muED = mumax  - inc;
  if (rateparam.muAM > mumax)  rateparam.muAM = mumax  - inc;
  if (rateparam.muAD > mumax)  rateparam.muAD = mumax  - inc;
  if (rateparam.muAI > mumax)  rateparam.muAI = mumax  - inc;
  if (rateparam.rM   > rMmax)  rateparam.rM   = rMmax  - inc;
  if (rateparam.rD   > rDImax) rateparam.rD   = rDImax - inc;
  if (rateparam.rI   > rDImax) rateparam.rI   = rDImax - inc;
  if (rateparam.sI   > sImax)  rateparam.sI   = sImax  - inc;

  switch(R->evomodel) {
  case AALI:
    p[x++] = RATE2PARAM(rateparam.ldI);  
    p[x++] = RATE2PARAM(rateparam.muI);      
    p[x++] = RATE2PARAM(rateparam.muAM); 
    p[x++] = RATE2PARAM(rateparam.muAD); 
    p[x++] = RATE2PARAM(rateparam.muAI); 
    break;
  case LI:
    p[x++] = RATE2PARAM(rateparam.ldI);  
    p[x++] = RATE2PARAM(rateparam.muI);      
    p[x++] = RATE2PARAM(rateparam.muAM); 
    break;
  case LR:
    p[x++] = RATE2PARAM(rateparam.ldI);           
    p[x++] = RATE2PARAM(rateparam.muAM); 
    break;
  case AFG:
    p[x++] = RATE2PARAM(rateparam.ldI);      
    p[x++] = RATE2PARAM(rateparam.muI);      
    p[x++] = RATE2PARAM(rateparam.muAM); 
    
    p[x++] = BERN2PARAMLOG(rateparam.rI);
    p[x++] = BERN2PARAMLOG(rateparam.rM);
    p[x++] = BERN2PARAMLOG(rateparam.rD);
    break;
  case AFGX:
    p[x++] = RATE2PARAM(rateparam.ldI);      
    p[x++] = RATE2PARAM(rateparam.muI);      
    p[x++] = RATE2PARAM(rateparam.muAM); 
    p[x++] = RATE2PARAM(rateparam.muAD); 
    p[x++] = RATE2PARAM(rateparam.muAI); 
    
    p[x++] = BERN2PARAMLOG(rateparam.rI);
    p[x++] = BERN2PARAMLOG(rateparam.rM);
    p[x++] = BERN2PARAMLOG(rateparam.rD);
    break;
  case AFGR:
    p[x++] = RATE2PARAM(rateparam.ldI);          
    p[x++] = RATE2PARAM(rateparam.muAM); 
     
    p[x++] = BERN2PARAMLOG(rateparam.rI);
    p[x++] = BERN2PARAMLOG(rateparam.rM);
    break;
  case AFR:
    p[x++] = RATE2PARAM(rateparam.ldI);            
    p[x++] = RATE2PARAM(rateparam.muAM); 
    
    p[x++] = BERN2PARAMLOG(rateparam.rI);
    break;
  case AIF:
    p[x++] = RATE2PARAM(rateparam.ldI);      
    p[x++] = RATE2PARAM(rateparam.muI);      
    p[x++] = RATE2PARAM(rateparam.muAM); 
 
    p[x++] = BERN2PARAMLOG(rateparam.rI);
    break;
  case AIFX:
    p[x++] = RATE2PARAM(rateparam.ldI);      
    p[x++] = RATE2PARAM(rateparam.muI);      
    p[x++] = RATE2PARAM(rateparam.muAM); 
    p[x++] = RATE2PARAM(rateparam.muAD); 
    p[x++] = RATE2PARAM(rateparam.muAI); 
 
    p[x++] = BERN2PARAMLOG(rateparam.rI);
    break;
  case GG:
    printf("model does not have analytical solution\n"); exit(1);
   break;
  case AG:
    printf("model does not have analytical solution\n"); exit(1);
   break;
  case AGA:
    p[x++] = RATE2PARAM(rateparam.ldEM);      
    p[x++] = RATE2PARAM(rateparam.muEM);      
    p[x++] = RATE2PARAM(rateparam.muAM); 

    p[x++] = BERN2PARAMLOG(rateparam.sI);
    break;
  case AGAX:
    p[x++] = RATE2PARAM(rateparam.ldEM);      
    p[x++] = RATE2PARAM(rateparam.muEM);      
    p[x++] = RATE2PARAM(rateparam.muAM); 
    p[x++] = RATE2PARAM(rateparam.muAD); 
    p[x++] = RATE2PARAM(rateparam.muAI); 

    p[x++] = BERN2PARAMLOG(rateparam.sI);
    break;
  case TKF91: /* pack muI-ldI, since we have the constrain muI > ldI */
    p[x++] = RATE2PARAM(rateparam.ldI);      
    p[x++] = RATE2PARAM(rateparam.muI - rateparam.ldI);      
    break;
  case TKF92:
    p[x++] = RATE2PARAM(rateparam.ldI);      
    p[x++] = RATE2PARAM(rateparam.muI - rateparam.ldI);      
     
    p[x++] = BERN2PARAMLOG(rateparam.rI);
    break;
  case FID:
    p[x++] = RATE2PARAM(rateparam.ldI);      
     
    p[x++] = BERN2PARAMLOG(rateparam.rI);
    break;
  case GTKF92:
    p[x++] = RATE2PARAM(rateparam.ldI);      
    p[x++] = RATE2PARAM(rateparam.muI - rateparam.ldI);      
     
    p[x++] = BERN2PARAMLOG(rateparam.rI);
    p[x++] = BERN2PARAMLOG(rateparam.rM);
    p[x++] = BERN2PARAMLOG(rateparam.rD);
    break;
  case GRTKF92:
    p[x++] = RATE2PARAM(rateparam.ldI);      
    p[x++] = RATE2PARAM(rateparam.muI - rateparam.ldI);      
    
    p[x++] = BERN2PARAMLOG(rateparam.rI);
    p[x++] = BERN2PARAMLOG(rateparam.rM);
    break;
  case ITKF92:
    p[x++] = RATE2PARAM(rateparam.ldI);      
    p[x++] = RATE2PARAM(rateparam.muI - rateparam.ldI);      
    
    p[x++] = BERN2PARAMLOG(rateparam.rI);
    break;
  default:
    ESL_XFAIL(eslFAIL, errbuf, "could not find this evomodel");
    break;
    }
  
  rparam = x - 1;
  if (rparam > R->nrate + R->nbern) ESL_XFAIL(eslFAIL, errbuf, "too many parameters? is %d should be %d \n", rparam, R->nrate + R->nbern); 

  *ret_x = x;

  return status;

 ERROR:
  return status;
} 


int
e2_transitions_unpack_paramvector(int *ret_x, double *p, long np, E1_RATE *R, char *errbuf, int verbose) 
{ 
  struct rateparam_s rateparam;  
  double             small   = 1e-4;
  double             inc     = 1e-4;
  double             ldmin   = 0.01;
  double             ldmax   = 0.69;
  double             mumin   = 0.01+inc;
  double             mumax   = 0.7;
  double             rMmin   = 0.80;
  double             rMmax   = 0.9999;
  double             rDImin  = 0.30;
  double             rDImax  = 0.95;
  double             sImin   = 0.30;
  double             sImax   = 0.95;
  int                rparam;
  int                x       = *ret_x;
  int                status  = eslOK;
  
  rateparam.ldEM = -1.0;
  rateparam.ldED = -1.0;
  rateparam.muEM = -1.0;
  rateparam.muED = -1.0;
  rateparam.ldI  = -1.0;
  rateparam.muI  = -1.0;
  rateparam.muAM = -1.0;
  rateparam.muAD = -1.0;
  rateparam.muAI = -1.0;
  rateparam.rI   = -1.0;
  rateparam.rM   = -1.0;
  rateparam.rD   = -1.0;
  rateparam.sI   = -1.0;
  rateparam.sD   = -1.0;
  rateparam.vI   = -1.0;
  rateparam.vD   = -1.0;
  
  switch(R->evomodel) {
  case AALI:
    rateparam.ldI  = PARAM2RATE(p[x]); x++;
    rateparam.muI  = PARAM2RATE(p[x]); x++; 
    rateparam.muAM = PARAM2RATE(p[x]); x++; 
    rateparam.muAD = PARAM2RATE(p[x]); x++; 
    rateparam.muAI = PARAM2RATE(p[x]); x++; 
    break;
  case LI:
    rateparam.ldI  = PARAM2RATE(p[x]); x++;
    rateparam.muI  = PARAM2RATE(p[x]); x++; 
    rateparam.muAM = PARAM2RATE(p[x]); x++; 
    break;
  case LR:
    rateparam.ldI  = PARAM2RATE(p[x]); x++; 
    rateparam.muAM = PARAM2RATE(p[x]); x++;
    rateparam.muI  = rateparam.ldI + rateparam.muAM; 
    break;
  case AFG:
    rateparam.ldI  = PARAM2RATE(p[x]); x++;
    rateparam.muI  = PARAM2RATE(p[x]); x++;
    rateparam.muAM = PARAM2RATE(p[x]); x++;
    
    rateparam.rI   = PARAM2BERNLOG(p[x]); x++;
    rateparam.rM   = PARAM2BERNLOG(p[x]); x++;
    rateparam.rD   = PARAM2BERNLOG(p[x]); x++;
    break;
  case AFGX:
    rateparam.ldI  = PARAM2RATE(p[x]); x++;
    rateparam.muI  = PARAM2RATE(p[x]); x++;
    rateparam.muAM = PARAM2RATE(p[x]); x++;
    rateparam.muAD = PARAM2RATE(p[x]); x++;
    rateparam.muAI = PARAM2RATE(p[x]); x++;
    
    rateparam.rI   = PARAM2BERNLOG(p[x]); x++;
    rateparam.rM   = PARAM2BERNLOG(p[x]); x++;
    rateparam.rD   = PARAM2BERNLOG(p[x]); x++;
    break;
  case AFGR:
    rateparam.ldI  = PARAM2RATE(p[x]); x++; 
    rateparam.muAM = PARAM2RATE(p[x]); x++; 
    rateparam.muI  = rateparam.ldI + rateparam.muAM; 
    
    rateparam.rI   = PARAM2BERNLOG(p[x]); x++;
    rateparam.rM   = PARAM2BERNLOG(p[x]); x++;
    break;
  case AFR:
    rateparam.ldI  = PARAM2RATE(p[x]); x++;
    rateparam.muAM = PARAM2RATE(p[x]); x++; 
    rateparam.muI  = rateparam.ldI + rateparam.muAM; 
    
    rateparam.rI = PARAM2BERNLOG(p[x]); x++;
    break;
  case AIF:
    rateparam.ldI  = PARAM2RATE(p[x]); x++;
    rateparam.muI  = PARAM2RATE(p[x]); x++;
    rateparam.muAM = PARAM2RATE(p[x]); x++;
    
    rateparam.rI = PARAM2BERNLOG(p[x]); x++;
    break;
  case AIFX:
    rateparam.ldI  = PARAM2RATE(p[x]); x++;
    rateparam.muI  = PARAM2RATE(p[x]); x++;
    rateparam.muAM = PARAM2RATE(p[x]); x++;
    rateparam.muAD = PARAM2RATE(p[x]); x++;
    rateparam.muAI = PARAM2RATE(p[x]); x++;
    
    rateparam.rI = PARAM2BERNLOG(p[x]); x++;
    break;
  case GG:
    printf("model does not have analytical solution\n"); exit(1);
    break;
  case AG:
    printf("model does not have analytical solution\n"); exit(1);
   break;
  case AGA:
    rateparam.ldEM = PARAM2RATE(p[x]); x++; 
    rateparam.muEM = PARAM2RATE(p[x]); x++; 
    rateparam.muAM = PARAM2RATE(p[x]); x++; 

    rateparam.sI   = PARAM2BERNLOG(p[x]); x++;
    break;
   case AGAX:
    rateparam.ldEM = PARAM2RATE(p[x]); x++; 
    rateparam.muEM = PARAM2RATE(p[x]); x++; 
    rateparam.muAM = PARAM2RATE(p[x]); x++; 
    rateparam.muAD = PARAM2RATE(p[x]); x++; 
    rateparam.muAI = PARAM2RATE(p[x]); x++; 

    rateparam.sI   = PARAM2BERNLOG(p[x]); x++;
    break;
  case TKF91:
    rateparam.ldI  = PARAM2RATE(p[x]);  x++;
    rateparam.muI  = PARAM2RATE(p[x]) + rateparam.ldI + small; x++; 
    break;
  case TKF92:
    rateparam.ldI  = PARAM2RATE(p[x]); x++;
    rateparam.muI  = PARAM2RATE(p[x]) + rateparam.ldI + small; x++; 
    
    rateparam.rI = PARAM2BERNLOG(p[x]); x++;
    break;
  case FID:
     rateparam.ldI = PARAM2RATE(p[x]); x++;
     rateparam.muI = rateparam.ldI;
    
    rateparam.rI = PARAM2BERNLOG(p[x]); x++;
    break;
  case GTKF92:
    rateparam.ldI = PARAM2RATE(p[x]); x++;
    rateparam.muI = PARAM2RATE(p[x]) + rateparam.ldI + small; x++; 
    
    rateparam.rI  = PARAM2BERNLOG(p[x]); x++;
    rateparam.rM  = PARAM2BERNLOG(p[x]); x++;
    rateparam.rD  = PARAM2BERNLOG(p[x]); x++;
    break;
  case GRTKF92:
    rateparam.ldI = PARAM2RATE(p[x]); x++;
    rateparam.muI = PARAM2RATE(p[x]) + rateparam.ldI + small; x++; 
    
    rateparam.rI = PARAM2BERNLOG(p[x]); x++;
    rateparam.rM = PARAM2BERNLOG(p[x]); x++;
    break;
  case ITKF92:
    rateparam.ldI = PARAM2RATE(p[x]); x++;
    rateparam.muI = PARAM2RATE(p[x]) + rateparam.ldI + small; x++; 
    
    rateparam.rI = PARAM2BERNLOG(p[x]); x++;
    break;
  default:
    ESL_XFAIL(eslFAIL, errbuf, "could not find this evomodel");
    break;
  }

  if (rateparam.ldI  > -1.0 && rateparam.ldI  < ldmin)  rateparam.ldI  = ldmin  + inc;
  if (rateparam.ldEM > -1.0 && rateparam.ldEM < ldmin)  rateparam.ldEM = ldmin  + inc;
  if (rateparam.ldED > -1.0 && rateparam.ldED < ldmin)  rateparam.ldED = ldmin  + inc;
  if (rateparam.muI  > -1.0 && rateparam.muI  < mumin)  rateparam.muI  = mumin  + inc;
  if (rateparam.muEM > -1.0 && rateparam.muEM < mumin)  rateparam.muEM = mumin  + inc;
  if (rateparam.muED > -1.0 && rateparam.muED < mumin)  rateparam.muED = mumin  + inc;
  if (rateparam.muAM > -1.0 && rateparam.muAM < mumin)  rateparam.muAM = mumin  + inc;
  if (rateparam.muAD > -1.0 && rateparam.muAD < mumin)  rateparam.muAD = mumin  + inc;
  if (rateparam.muAI > -1.0 && rateparam.muAI < mumin)  rateparam.muAI = mumin  + inc;
  if (rateparam.rM   > -1.0 && rateparam.rM   < rMmin)  rateparam.rM   = rMmin  + inc;
  if (rateparam.rD   > -1.0 && rateparam.rD   < rDImin) rateparam.rD   = rDImin + inc;
  if (rateparam.rI   > -1.0 && rateparam.rI   < rDImin) rateparam.rI   = rDImin + inc;
  if (rateparam.sI   > -1.0 && rateparam.sI   < sImin)  rateparam.sI   = sImin  + inc;

  if (rateparam.ldI  > ldmax)  rateparam.ldI  = ldmax  - inc;
  if (rateparam.ldEM > ldmax)  rateparam.ldEM = ldmax  - inc;
  if (rateparam.ldED > ldmax)  rateparam.ldED = ldmax  - inc;
  if (rateparam.muI  > mumax)  rateparam.muI  = mumax  - inc;
  if (rateparam.muEM > mumax)  rateparam.muEM = mumax  - inc;
  if (rateparam.muED > mumax)  rateparam.muED = mumax  - inc;
  if (rateparam.muAM > mumax)  rateparam.muAM = mumax  - inc;
  if (rateparam.muAD > mumax)  rateparam.muAD = mumax  - inc;
  if (rateparam.muAI > mumax)  rateparam.muAI = mumax  - inc;
  if (rateparam.rM   > rMmax)  rateparam.rM   = rMmax  - inc;
  if (rateparam.rD   > rDImax) rateparam.rD   = rDImax - inc;
  if (rateparam.rI   > rDImax) rateparam.rI   = rDImax - inc;
  if (rateparam.sI   > sImax)  rateparam.sI   = sImax  - inc;

  rparam = x - 1;
  if (rparam > R->nrate + R->nbern) ESL_XFAIL(eslFAIL, errbuf, "too many parameters? is %d should be %d \n", rparam, R->nrate + R->nbern); 
  
  status = e1_rate_AssignTransitionsFromRates(R, rateparam, errbuf, verbose);

  *ret_x = x;

  return status;

 ERROR:
  return status;
}


/* -- internal functions -- */

static void
e2_train_bracket_define_direction(double *u, long np)
{
  int x;

  for (x = 0; x < np; x++) u[x] = 1.0;
  u[np] = 1.0;
}


static void
e2_train_gd_define_stepsize(double *u, long np, struct e2_train_data *data)
{
  int x;

  for (x = 0; x < np; x++) u[x] = (x < data->R->nrate)? 0.9 : 1e+4;
  u[np] = 1.0;
}

static void
e2_train_NM_define_stepsize(double *u, long np, struct e2_train_data *data)
{
  int x;

  for (x = 0; x < np; x++) u[x] = (x < data->R->nrate)? 0.9 : 1e+4;
  u[np] = 1.0;
}

static void
e2_train_pack_paramvector(double *p, long np, struct e2_train_data *data)
{
  int  x = 0; 
  int  status;

  status = e2_transitions_pack_paramvector(&x, p, np, data->R, data->errbuf, data->verbose);
  if (status != eslOK) { printf("error at train_pack_paramvector()\n%s\n", data->errbuf); exit(1); }
}
static void
e2_train_unpack_paramvector(double *p, long np, struct e2_train_data *data)
{
  int  x = 0; 
  int  status;

  status = e2_transitions_unpack_paramvector(&x, p, np, data->R, data->errbuf, data->verbose);
  if (status != eslOK) { printf("error at train_unpack_paramvector().\n%s\n", data->errbuf); exit(1); }
}
     
static double 
e2_train_func(double *p, long np, void *dptr)
{
  struct e2_train_data *data = (struct e2_train_data *) dptr;
  float     tinit = -1.0;
  float     sc;
  int       jump = 10;
  int       n;
  int       status;

  e2_train_unpack_paramvector(p, (double)np, data);
  
  data->sc = 0.0;
  for (n = 0; n < data->nmsa; n ++) {
    
    if (data->it % jump == jump-1) {
      tinit = esl_tree_er_AverageBL(data->Tlist[n]);
      if (data->Tlist[n]) esl_tree_Destroy(data->Tlist[n]); data->Tlist[n] = NULL;
      e2_tree_UPGMA(&data->Tlist[n], 0, NULL, data->msalist[n], data->msafrq[n], data->r, data->pli,
		    data->R, NULL, data->bg, NULL, data->e2ali, data->mode, data->do_viterbi, -1.0, tinit, 
		    data->tol, data->errbuf, data->verbose);
    }
    
    status = e2_msa(data->r, data->R, NULL, 0, NULL, data->msalist[n], data->msafrq[n], data->Tlist[n], NULL,
		    &sc, data->pli, data->bg, NULL, data->e2ali, OPTNONE, 
		    data->mode, data->do_viterbi, data->tol, data->errbuf, data->verbose);
    if (status != eslOK) { printf("error at e2_train_func()\n%s\n", data->errbuf); exit(1); }
    
    data->sc += sc;
  } 
 
  data->it ++;

#if 1
  printf("IT %d: rI %f rM %f rD %f | sI %f | ldEM %f muEM %f ldED %f muED %f ldI %f muI %f muAM %f muAD %f muAI %f | sc %f\n", data->it, 
	 data->R->rI, data->R->rM, data->R->rD, data->R->sI, 
	 data->R->ldE[e1R_S], data->R->muE[e1R_S], data->R->ldE[e1R_D], data->R->muE[e1R_D], 
	 data->R->ldE[e1R_I], data->R->muE[e1R_I], 
	 data->R->muA[e1R_S], data->R->muA[e1R_D], data->R->muA[e1R_I], 
	 data->sc);
#endif
  
  return -data->sc;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
