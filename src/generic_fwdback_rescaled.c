/*
 *  gcc -o generic_fwdback_rescaled_example -mssse3 -I. -I../easel -L../easel -L. generic_fwdback_rescaled.c  -lhmmer -leasel
 */

#include "p7_config.h"

#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

extern int p7_GForwardOdds(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx, float *opt_sc);
extern int p7_profile_CopyInfoFromHMM(P7_PROFILE *gm, const P7_HMM *hmm);
extern int p7_profile_ReconfigLengthInOdds(P7_PROFILE *gm, int L);
extern int p7_profile_ConfigInOdds(const P7_HMM *hmm, const P7_BG *bg, P7_PROFILE *gm, int L, int mode, float *ret_ddscale);

#define STYLES     "--fs,--sw,--ls,--s"	               /* Exclusive choice for alignment mode     */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "--fs",      eslARG_NONE,"default",NULL, NULL, STYLES,  NULL, NULL, "multihit local alignment",                         0 },
  { "--sw",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "unihit local alignment",                           0 },
  { "--ls",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "multihit glocal alignment",                        0 },
  { "--s",       eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "unihit glocal alignment",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "testbed for Farrar DD-scaledl Forward, generic implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *fwd     = NULL;
  P7_GMX         *bck     = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fsc, bsc;
  float           nullsc;
  float           ddscale;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Open the sqfile */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);

  /* Now reconfig the models however we were asked to */
  if      (esl_opt_GetBoolean(go, "--fs"))  p7_profile_ConfigInOdds(hmm, bg, gm, sq->n, p7_LOCAL,    &ddscale);
  else if (esl_opt_GetBoolean(go, "--sw"))  p7_profile_ConfigInOdds(hmm, bg, gm, sq->n, p7_UNILOCAL, &ddscale);
  else if (esl_opt_GetBoolean(go, "--ls"))  p7_profile_ConfigInOdds(hmm, bg, gm, sq->n, p7_GLOCAL,   &ddscale);
  else if (esl_opt_GetBoolean(go, "--s"))   p7_profile_ConfigInOdds(hmm, bg, gm, sq->n, p7_UNIGLOCAL,&ddscale);
  
  /* Allocate matrices */
  fwd = p7_gmx_Create(gm->M, sq->n);
  bck = p7_gmx_Create(gm->M, sq->n);

  printf("%-15s   %-10s %-10s   %-10s %-10s\n", "# seq name",      "fwd (raw)",   "bck (raw) ",  "fwd (bits)",  "bck (bits)");
  printf("%15s   %10s %10s   %10s %10s\n",      "#--------------", "----------",  "----------",  "----------",  "----------");

  while ( (status = esl_sqio_Read(sqfp, sq)) != eslEOF)
    {
      if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
      else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

      /* Resize the DP matrices if necessary */
      p7_gmx_GrowTo(fwd, gm->M, sq->n);
      p7_gmx_GrowTo(bck, gm->M, sq->n);

      /* Set the profile and null model's target length models */
      p7_bg_SetLength(bg, sq->n);
      p7_profile_ReconfigLengthInOdds(gm, sq->n);

      /* Run Forward */
      p7_GForwardOdds (sq->dsq, sq->n, gm, fwd, &fsc);
      fsc += ddscale;

      /* Those scores are partial log-odds likelihoods in nats.
       * Subtract off the rest of the null model, convert to bits.
       */
      p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

      printf("%-15s   %10.4f %10.4f   %10.4f %10.4f\n", 
	     sq->name, 
	     fsc, bsc, 
	     (fsc - nullsc) / eslCONST_LOG2, (bsc - nullsc) / eslCONST_LOG2);

      p7_gmx_Reuse(fwd);
      p7_gmx_Reuse(bck);
      esl_sq_Reuse(sq);
    }

  /* Cleanup */
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}

int
p7_GForwardOdds(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx, float *opt_sc)
{
  float const *tsc = gm->tsc;
  float      **dp  = gx->dp;
  float       *xmx = gx->xmx;
  int          M   = gm->M;
  int          i,k;
  float        esc  = p7_profile_IsLocal(gm) ? 1.0f : 0.0f;
  float        totscale = 0.0f;
  float        rescale;

  /* Initialization of the zero row */
  XMX(0,p7G_N) = 1.0f;		              /* S->N, p=1 */
  XMX(0,p7G_B) = gm->xsc[p7P_N][p7P_MOVE];    /* S->N->B, 1*t_NB */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = 0.0f; /* need seq to get here */
  for (k = 0; k <= M; k++) 
    MMX(0,k) = IMX(0,k) = DMX(0,k) = 0.0f;
  
  /* Recursion
   * boundary conditions: tsc[0][] = impossible for all eight transitions (no node 0)
   *                      D_1 wastefully calculated (doesn't exist)
   */
  for (i = 1; i <= L; i++)
    {
      float const *rsc = gm->rsc[dsq[i]];
      float sc;

      MMX(i,0) = IMX(i,0) = DMX(i,0) = 0.0f;
      XMX(i, p7G_E) = 0.0f;
      
      for (k = 1; k < M; k++)
	{
	  /* match state */
	  sc =   MMX(i-1,k-1)   * TSC(p7P_MM,k-1) 
	       + IMX(i-1,k-1)   * TSC(p7P_IM,k-1)
	       + XMX(i-1,p7G_B) * TSC(p7P_BM,k-1)
	       + DMX(i-1,k-1)   * TSC(p7P_DM,k-1);
	  MMX(i,k) = sc * MSC(k);

	  /* insert state */
	  sc =   MMX(i-1,k) * TSC(p7P_MI,k)
	       + IMX(i-1,k) * TSC(p7P_II,k);
	  IMX(i,k) = sc * ISC(k);

	  /* delete state */
	  DMX(i,k) =   MMX(i,k-1) * TSC(p7P_MD,k-1)
	             + DMX(i,k-1) * TSC(p7P_DD,k-1);
	  
	  /* E state update */
	  XMX(i,p7G_E) +=   MMX(i,k) * esc
                          + DMX(i,k) * esc;
	}
      /* unrolled match state M_M */
      sc =   MMX(i-1,M-1)   * TSC(p7P_MM,M-1) 
	   + IMX(i-1,M-1)   * TSC(p7P_IM,M-1)
	   + XMX(i-1,p7G_B) * TSC(p7P_BM,M-1)
	   + DMX(i-1,M-1)   * TSC(p7P_DM,M-1);
      MMX(i,M) = sc * MSC(M);
      IMX(i,M) = 0.0f;

      /* unrolled delete state D_M */
      DMX(i,M) =   MMX(i,M-1) * TSC(p7P_MD,M-1)
	         + DMX(i,M-1) * TSC(p7P_DD,M-1);

      /* unrolled E state update */
      XMX(i,p7G_E) += MMX(i,M) + DMX(i,M);

      /* J state */
      XMX(i,p7G_J) =   XMX(i-1,p7G_J) * gm->xsc[p7P_J][p7P_LOOP]
                     + XMX(i,  p7G_E) * gm->xsc[p7P_E][p7P_LOOP];
      /* C state */
      XMX(i,p7G_C) =   XMX(i-1,p7G_C) * gm->xsc[p7P_C][p7P_LOOP]
		     + XMX(i,  p7G_E) * gm->xsc[p7P_E][p7P_MOVE];
      /* N state */
      XMX(i,p7G_N) =   XMX(i-1,p7G_N) * gm->xsc[p7P_N][p7P_LOOP];

      /* B state */
      XMX(i,p7G_B) =   XMX(i,  p7G_N) * gm->xsc[p7P_N][p7P_MOVE]
		     + XMX(i,  p7G_J) * gm->xsc[p7P_J][p7P_MOVE];

      /* sparse rescaling */
      if (XMX(i,p7G_E) > 1.0e4)
	{
	  rescale   = 1.0 / XMX(i,p7G_E);
	  totscale += log(XMX(i,p7G_E));

	  XMX(i,p7G_N) *= rescale;
	  XMX(i,p7G_B) *= rescale;
	  XMX(i,p7G_E)  = 1.0;
	  XMX(i,p7G_J) *= rescale;
	  XMX(i,p7G_C) *= rescale;
	  for (k = 1; k <= M; k++)
	    {
	      MMX(i,k) *= rescale;
	      DMX(i,k) *= rescale;
	      IMX(i,k) *= rescale;
	    }
	}
    }

  if (opt_sc != NULL) *opt_sc = log(XMX(L,p7G_C) * gm->xsc[p7P_C][p7P_MOVE]) + totscale;
  gx->M = M;
  gx->L = L;
  return eslOK;
}


/* this is a copy of a block of modelconfig.c::p7_ProfileConfig() 
 * ...and could replace it
 */
int
p7_profile_CopyInfoFromHMM(P7_PROFILE *gm, const P7_HMM *hmm)
{
  int z;
  int status;

  /* Contract checks */
  if (gm->abc->type != hmm->abc->type) ESL_XEXCEPTION(eslEINVAL, "HMM and profile alphabet don't match");
  if (hmm->M > gm->allocM)             ESL_XEXCEPTION(eslEINVAL, "profile too small to hold HMM");
  if (! (hmm->flags & p7H_CONS))       ESL_XEXCEPTION(eslEINVAL, "HMM must have a consensus to transfer to the profile");

  /* Copy some pointer references and other info across from HMM  */
  gm->M                = hmm->M;
  gm->max_length       = hmm->max_length;
  gm->mode             = p7_NO_MODE;
  gm->roff             = -1;
  gm->eoff             = -1;
  gm->offs[p7_MOFFSET] = -1;
  gm->offs[p7_FOFFSET] = -1;
  gm->offs[p7_POFFSET] = -1;
  if (gm->name != NULL) free(gm->name);
  if (gm->acc  != NULL) free(gm->acc);
  if (gm->desc != NULL) free(gm->desc);
  if ((status = esl_strdup(hmm->name,   -1, &(gm->name))) != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->acc,    -1, &(gm->acc)))  != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->desc,   -1, &(gm->desc))) != eslOK) goto ERROR;
  if (hmm->flags & p7H_RF)   strcpy(gm->rf,        hmm->rf);
  if (hmm->flags & p7H_CONS) strcpy(gm->consensus, hmm->consensus); /* must be present, actually, so the flag test is just for symmetry w/ other optional HMM fields */
  if (hmm->flags & p7H_CS)   strcpy(gm->cs,        hmm->cs);
  for (z = 0; z < p7_NEVPARAM; z++) gm->evparam[z] = hmm->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) gm->cutoff[z]  = hmm->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) gm->compo[z]   = hmm->compo[z];

  return eslOK;

 ERROR:
  return status;
}


int
p7_profile_ReconfigLengthInOdds(P7_PROFILE *gm, int L)
{
  float ploop, pmove;
  
  /* Configure N,J,C transitions so they bear L/(2+nj) of the total
   * unannotated sequence length L. 
   */
  pmove = (2.0f + gm->nj) / ((float) L + 2.0f + gm->nj); /* 2/(L+2) for sw; 3/(L+3) for fs */
  ploop = 1.0f - pmove;
  gm->xsc[p7P_N][p7P_LOOP] =  gm->xsc[p7P_C][p7P_LOOP] = gm->xsc[p7P_J][p7P_LOOP] = ploop;
  gm->xsc[p7P_N][p7P_MOVE] =  gm->xsc[p7P_C][p7P_MOVE] = gm->xsc[p7P_J][p7P_MOVE] = pmove;
  gm->L = L;
  return eslOK;
}


/* like p7_ProfileConfig(), but in odds ratios rather than log 
 * In local mode, transitions/emissions are odds ratios
 * In glocal mode, they are scaled odds ratios, divided by t_k(DD),
 * using Farrar's trick.
 */
int
p7_profile_ConfigInOdds(const P7_HMM *hmm, const P7_BG *bg, P7_PROFILE *gm, int L, int mode, float *ret_ddscale)
{
  int    k,x;
  float *occ = NULL;
  float  ddscale;
  float *tp, *rp;
  float  sc[p7_MAXCODE];
  float  Z;
  int    status;

  if ( (status = p7_profile_CopyInfoFromHMM(gm, hmm)) != eslOK) return status;
  gm->mode = mode;

  /* Entry scores. Recall k-1,k off-by-one storage issue here. 
   * p7P_TSC(gm, k-1, p7P_BM) is the t_BMk transition 
   */
  if (p7_profile_IsLocal(gm))
    {
      /* Local mode entry:  occ[k] /( \sum_i occ[i] * (M-i+1))
       * (Reduces to uniform 2/(M(M+1)) for occupancies of 1.0)  */
      Z = 0.;
      ESL_ALLOC(occ, sizeof(float) * (hmm->M+1));

      if ((status = p7_hmm_CalculateOccupancy(hmm, occ, NULL)) != eslOK) goto ERROR;
      for (k = 1; k <= hmm->M; k++) 
	Z += occ[k] * (float) (hmm->M-k+1);
      for (k = 1; k <= hmm->M; k++) 
	p7P_TSC(gm, k-1, p7P_BM) = occ[k] / Z; /* note off-by-one: entry at Mk stored as [k-1][BM] */
      free(occ);
    }
  else	/* glocal modes: left wing retraction; must be DD-scaled, else odds ratios will underflow */
    {	
      p7P_TSC(gm, 0, p7P_BM) = 1.0 - hmm->t[0][p7H_MD];
      for (k = 1; k < hmm->M; k++) 
	p7P_TSC(gm, k, p7P_BM) =  hmm->t[0][p7H_MD] * hmm->t[k][p7H_DM] / hmm->t[k][p7H_DD];
    }

  /* E state loop/move probabilities: nonzero for MOVE allows loops/multihits
   * N,C,J transitions are set later by target length model config
   */
  if (p7_profile_IsMultihit(gm)) {
    gm->xsc[p7P_E][p7P_MOVE] = 0.5f;
    gm->xsc[p7P_E][p7P_LOOP] = 0.5f;
    gm->nj                   = 1.0f;
  } else {
    gm->xsc[p7P_E][p7P_MOVE] = 1.0f;   
    gm->xsc[p7P_E][p7P_LOOP] = 0.0f;
    gm->nj                   = 0.0f;
  }
  
  /* main profile transition scores */
  if (p7_profile_IsLocal(gm))
    {
      for (k = 1; k < gm->M; k++) 
	{
	  tp = gm->tsc + k * p7P_NTRANS;
	  tp[p7P_MM] = hmm->t[k][p7H_MM];
	  tp[p7P_MI] = hmm->t[k][p7H_MI];
	  tp[p7P_MD] = hmm->t[k][p7H_MD];
	  tp[p7P_IM] = hmm->t[k][p7H_IM];
	  tp[p7P_II] = hmm->t[k][p7H_II];
	  tp[p7P_DM] = hmm->t[k][p7H_DM];
	  tp[p7P_DD] = hmm->t[k][p7H_DD];
	}
    }
  else  /* glocal mode used DD-scaling */
    {
      for (k = 1; k < gm->M; k++)
	{
	  tp      = gm->tsc + k * p7P_NTRANS;
	  ddscale = 1.0f/hmm->t[k][p7H_DD];
      
	  tp[p7P_MM] = hmm->t[k][p7H_MM] * ddscale;
	  tp[p7P_MI] = hmm->t[k][p7H_MI];
	  tp[p7P_MD] = hmm->t[k][p7H_MD] * ddscale;
	  tp[p7P_IM] = hmm->t[k][p7H_IM] * ddscale;
	  tp[p7P_II] = hmm->t[k][p7H_II];
	  tp[p7P_DM] = hmm->t[k][p7H_DM] * ddscale;
	  tp[p7P_DD] = 1.0;
	}
    }
				 
  /* compute log(\prod_k=1..M-1 t_k(DD), total rescale factor */
  ddscale = 0.0f;
  if (! p7_profile_IsLocal(gm)) /* glocal mode: nonzero ddscale  */
    {
      for (k = 1; k < gm->M; k++)
	ddscale += log(hmm->t[k][p7H_DD]);
    }

  /* match residue scores */
  /* we still compute this in log space, then exp() it back,
   * because degenerate residue scores are avg log odds.
   */
  sc[hmm->abc->K]     = -eslINFINITY; /* gap character */
  sc[hmm->abc->Kp-2]  = -eslINFINITY; /* nonresidue character */
  sc[hmm->abc->Kp-1]  = -eslINFINITY; /* missing data character */
  for (k = 1; k <= hmm->M; k++) {
    for (x = 0; x < hmm->abc->K; x++) 
      sc[x] = log(hmm->mat[k][x] / bg->f[x]);
    esl_abc_FExpectScVec(hmm->abc, sc, bg->f); 
    for (x = 0; x < hmm->abc->Kp; x++) {
      rp = gm->rsc[x] + k * p7P_NR;
      rp[p7P_MSC] = exp(sc[x]);
    }
  }

  

  /* Remaining specials, [NCJ][MOVE | LOOP] are set by ReconfigLengthInOdds() */
  gm->L = 0;			/* force ReconfigLengthInOdds to reconfig */
  if ((status = p7_profile_ReconfigLengthInOdds(gm, L)) != eslOK) goto ERROR;

  *ret_ddscale = ddscale;
  return eslOK;

 ERROR:
  if (occ) free(occ);
  return status;
}




