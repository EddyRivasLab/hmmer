/* The Forwards filter:
 *   - in striped SIMD vectors;
 *   - in checkpointed memory, O(M sqrt(L)); 
 *   - in probability space, with sparse rescaling;
 *   - restricted to multihit local alignment mode only (numeric range issues);
 *   - closely tied to the implementation of the Backwards filter in bckfilter.c.
 *
 * In the acceleration pipeline:
 * SSVFilter --> MSVFilter --> VitFilter --> FwdFilter --> BckFilter
 *                                           ^^^^^^^^^
 *                                        (you are here)
 * Contents:
 *    x. 
 *    x. 
 *    x. Copyright and license information.
 */
#include "p7_config.h"

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_sse.h"
#include "p7_filtermx.h"

static inline float forward_row(int Q, const P7_OPROFILE *om, const __m128 *rp, const __m128 *dpp, __m128 *dpc);

/* Function:  p7_ForwardFilter()
 * Synopsis:  Checkpointed, SIMD striped vector Forward calculation.
 *
 * Purpose:   Calculate the Forward algorithm for target sequence <dsq>
 *            of <L> residues aligned to query profile <om>, using the
 *            checkpointed DP matrix <ox> provided by the caller. Upon
 *            successful return, <ox> contains the filled Forward
 *            matrix, and <*opt_sc> optionally contains the raw Forward
 *            score in nats.
 *            
 * Args:      dsq    - digital target sequence, 1..L
 *            L      - length of dsq, residues
 *            om     - optimized profile (multihit local)
 *            ox     - checkpointed DP matrix
 *            opt_sc - optRETURN: raw Forward score (nats)
 *
 * Returns:   <eslOK> on success, <ox> contains the checkpointed
 *            Forward matrix calculation, and <*opt_sc> optionally
 *            has the raw Forward score in nats.
 *
 * Throws:    (no abnormal error conditions)
 * 
 * Xref:      For layout of checkpointed <ox> see exegesis in p7_filtermx.h.
 */
int
p7_ForwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *opt_sc)
{
  int           Q     = P7F_NQF(om->M);                  /* segment length; # of MDI vectors on each row */
  __m128       *dpp   = (__m128 *) ox->dpf[ox->R0-1];    /* dpp points at prev row. start on dpf[2]; rows 0,1 are Backwards */
  __m128       *dpc   = NULL;		                 /* dpc points at current row */
  float        *xc    = (float *) (dpp + Q*p7F_NSCELLS); /* specials E,N,JJ,J,B,CC,C,SCALE */
  const __m128  zerov = _mm_setzero_ps();		 
  float         totsc = 0.0f;
  int     q;			/* counter over vectors 0..Q-1                        */
  int     i;			/* counter over residues/rows 1..L                    */
  int     b;			/* counter down through checkpointed blocks, Rb+Rc..1 */
  int     w;			/* counter down through rows in a checkpointed block  */

  /* Initialization of the zero row, including specials */
  for (q = 0; q < p7F_NSCELLS*Q; q++) dpp[q] = zerov;
  xc[p7F_N]     = 1.;
  xc[p7F_B]     = om->xf[p7O_N][p7O_MOVE]; 
  xc[p7F_E]     = xc[p7F_JJ] = xc[p7F_J]  = xc[p7F_CC] = xc[p7F_C]  = 0.;			
  xc[p7F_SCALE] = 1.;			   

  /* Phase one: the "a" region: all rows in this region are saved */
  for (i = 1; i <= ox->La; i++)
    {
      dpc = (__m128 *) ox->dpf[ox->R0+ox->R]; ox->R++; /* idiomatic for "get next save/checkpoint row" */
      totsc += forward_row(Q, om, om->rfv[dsq[i]], dpp, dpc);
      dpp = dpc;	    	           /* current row becomes prev row */
    }

  /* Phase two: "b" and "c" regions: partially and fully checkpointed */
  for (b = ox->Rb + ox->Rc, w = (ox->Rb ? ox->Lb : ox->Rc+1); i <= L; i++)
    {
      if (! (--w)) { 		                   /* we're on the last row in segment: this row is saved    */
	dpc = (__m128 *) ox->dpf[ox->R0+ox->R]; ox->R++;  /* idiomatic for "get next save/checkpoint row"    */
	w = b;  			           /* next segment has this many rows, ending in a saved row */
	b--;					   /* decrement segment number counter; last segment is r=1  */
      } else dpc = (__m128 *) ox->dpf[i%2];        /* idiomatic for "get next tmp row", 0/1; i%2 makes sure dpp != dpc */
      
      totsc += forward_row(Q, om, om->rfv[dsq[i]], dpp, dpc);
      dpp = dpc;
    }

  ox->M = om->M;
  ox->L = L;
  xc    = (float *) (dpc + Q*p7F_NSCELLS);

  ESL_DASSERT1( (ox->R == ox->Ra+ox->Rb+ox->Rc) );
  ESL_DASSERT1( ( (! isnan(xc[p7F_C])) && (! isinf(xc[p7F_C]))) );
  
  if (opt_sc) *opt_sc = totsc + logf(xc[p7F_C] * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


static inline float
forward_row(int Q, const P7_OPROFILE *om, const __m128 *rp, const __m128 *dpp, __m128 *dpc)
{
  const    __m128 zerov = _mm_setzero_ps();
  const    __m128 *tp   = om->tfv;
  const    float  *xp   = (float *) (dpp + Q * p7F_NSCELLS);
  float           *xc   = (float *) (dpc + Q * p7F_NSCELLS);
  __m128          dcv   = _mm_setzero_ps();
  __m128          xEv   = _mm_setzero_ps();
  __m128          xBv   = _mm_set1_ps(xp[p7F_B]);
  __m128 mpv, dpv, ipv;
  __m128 sv;
  int    q;
  int    j;

  mpv = esl_sse_rightshift_ps(P7F_MQ(dpp, Q-1), zerov); 
  ipv = esl_sse_rightshift_ps(P7F_IQ(dpp, Q-1), zerov); 
  dpv = esl_sse_rightshift_ps(P7F_DQ(dpp, Q-1), zerov); 

  /* DP recursion for main states, all but the D->D path */
  for (q = 0; q < Q; q++)
    {
      /* Calculate M(i,q); hold it in tmp var <sv> */
      sv     =                _mm_mul_ps(xBv, *tp);  tp++; /* B->Mk    */
      sv     = _mm_add_ps(sv, _mm_mul_ps(mpv, *tp)); tp++; /* Mk-1->Mk */
      sv     = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++; /* Ik-1->Mk */
      sv     = _mm_add_ps(sv, _mm_mul_ps(dpv, *tp)); tp++; /* Dk-1->Dk */
      sv     = _mm_mul_ps(sv, *rp);                  rp++; /* e_Mk(x_i)*/
      xEv    = _mm_add_ps(xEv, sv);			   /* Mk->E    */

      /* Advance on previous row, picking up M,D,I values. */
      mpv = *dpp++;
      dpv = *dpp++;
      ipv = *dpp++;

      /* Delayed store of M,D */
      P7F_MQ(dpc, q) = sv;		
      P7F_DQ(dpc, q) = dcv;

      /* Partial calculation of *next* D(i,q+1); M->D only; delay storage, hold in dcv */
      dcv    = _mm_mul_ps(sv, *tp); tp++;

      /* Calculate and store I(i,q) */
      sv             =                _mm_mul_ps(mpv, *tp);  tp++;
      P7F_IQ(dpc, q) = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++;
    }

  /* Now the DD paths. We would rather not serialize them but 
   * in an accurate Forward calculation, we have few options.
   * dcv has carried through from end of q loop above; store it 
   * in first pass, we add M->D and D->D path into DMX.
   */ 
  /* We're almost certainly're obligated to do at least one complete 
   * DD path to be sure: 
   */
  dcv            = esl_sse_rightshift_ps(dcv, zerov);
  P7F_DQ(dpc, 0) = zerov;
  tp             = om->tfv + 7*Q;	/* set tp to start of the DD's */
  for (q = 0; q < Q; q++) 
    {
      P7F_DQ(dpc,q) = _mm_add_ps(dcv, P7F_DQ(dpc,q));	
      dcv           = _mm_mul_ps(P7F_DQ(dpc,q), *tp); tp++; /* extend DMO(q), so we include M->D and D->D paths */
    }
  
  /* now. on small models, it seems best (empirically) to just go
   * ahead and serialize. on large models, we can do a bit better,
   * by testing for when dcv (DD path) accrued to DMO(q) is below
   * machine epsilon for all q, in which case we know DMO(q) are all
   * at their final values. The tradeoff point is (empirically) somewhere around M=100,
   * at least on my desktop. We don't worry about the conditional here;
   * it's outside any inner loops.
   */
  if (om->M < 100)
    {			/* Fully serialized version */
      for (j = 1; j < 4; j++)
	{
	  dcv = esl_sse_rightshift_ps(dcv, zerov);
	  tp  = om->tfv + 7*Q;	/* reset tp to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    { /* note, extend dcv, not DMO(q); only adding DD paths now */
	      P7F_DQ(dpc,q) = _mm_add_ps(dcv, P7F_DQ(dpc,q));	
	      dcv           = _mm_mul_ps(dcv, *tp);   tp++; 
	    }	    
	}
    } 
  else
    {			/* Slightly parallelized version, but which incurs some overhead */
      for (j = 1; j < 4; j++)
	{
	  register __m128 cv = zerov;	/* keeps track of whether any DD's change DMO(q) */

	  dcv = esl_sse_rightshift_ps(dcv, zerov);
	  tp  = om->tfv + 7*Q;	/* set tp to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
	      sv            = _mm_add_ps(dcv, P7F_DQ(dpc,q));	
	      cv            = _mm_or_ps(cv, _mm_cmpgt_ps(sv, P7F_DQ(dpc,q))); 
	      P7F_DQ(dpc,q) = sv;	                                    /* store new DMO(q) */
	      dcv           = _mm_mul_ps(dcv, *tp);   tp++;                 /* note, extend dcv, not DMO(q) */
	    }	    
	  if (! _mm_movemask_ps(cv)) break; /* DD's didn't change any DMO(q)? Then done, break out. */
	}
    }
  
  /* Add D's to xEv */
  for (q = 0; q < Q; q++) xEv = _mm_add_ps(P7F_DQ(dpc,q), xEv);

  /* Specials, in order: E N JJ J B CC C */
  esl_sse_hsum_ps(xEv, &xc[p7F_E]);  
  xc[p7F_N]  =                                       xp[p7F_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7F_JJ] =                                       xp[p7F_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7F_J]  = xc[p7F_JJ]                          + xp[p7F_E] * om->xf[p7O_E][p7O_LOOP];
  xc[p7F_B]  = xc[p7F_N] * om->xf[p7O_N][p7O_MOVE] + xc[p7F_J] * om->xf[p7O_J][p7O_MOVE];
  xc[p7F_CC] =                                       xp[p7F_C] * om->xf[p7O_C][p7O_LOOP];
  xc[p7F_C]  = xc[p7F_CC]                          + xc[p7F_E] * om->xf[p7O_E][p7O_MOVE];
  
  /* Sparse rescaling. xE above threshold? Then trigger a rescaling event.            */
  if (xc[p7F_E] > 1.0e4)	/* that's a little less than e^10, ~10% of our dynamic range */
    {
      xc[p7F_N]  /= xc[p7F_E];
      xc[p7F_JJ] /= xc[p7F_E];
      xc[p7F_J]  /= xc[p7F_E];
      xc[p7F_B]  /= xc[p7F_E];
      xc[p7F_CC] /= xc[p7F_E];
      xc[p7F_C]  /= xc[p7F_E];
      xEv = _mm_set1_ps(1.0 / xc[p7F_E]);

      for (q = 0; q < Q; q++)
	{
	  P7F_MQ(dpc,q) = _mm_mul_ps(P7F_MQ(dpc,q), xEv);
	  P7F_DQ(dpc,q) = _mm_mul_ps(P7F_DQ(dpc,q), xEv);
	  P7F_IQ(dpc,q) = _mm_mul_ps(P7F_IQ(dpc,q), xEv);
	}

      xc[p7F_SCALE] = xc[p7F_E];
      xc[p7F_E]     = 1.0f;
    }
  else xc[p7F_SCALE] = 1.0f;

  return logf(xc[p7F_SCALE]);
}

/*****************************************************************
 * x. Benchmark
 *****************************************************************/
/* Difference between gcc -g vs. icc -O3 is large! */

#ifdef p7FWDFILTER_BENCHMARK
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_sse.h"
#include "p7_filtermx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,   "2000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for checkpointed ForwardFilter()";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_FILTERMX    *ox      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);

  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, L);

  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  ox  = p7_filtermx_Create(om->M, L, ESL_MBYTES(32));

  /* Baseline time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_ForwardFilter(dsq, L, om, ox, &sc);
      p7_filtermx_Reuse(ox);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_filtermx_Destroy(ox);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7FWDFILTER_BENCHMARK*/
/*-------------------- end, benchmark ---------------------------*/


/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7FWDFILTER_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "impl_sse.h"
#include "p7_filtermx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example driver, ForwardFilter()";

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
  P7_OPROFILE    *om      = NULL;
  P7_GMX         *gx      = NULL;
  P7_FILTERMX    *ox      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fraw, nullsc, fsc;
  float           gfraw, gfsc;
  float           gmem, cmem;
  double          P, gP;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  /* Open sequence file for reading */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);               
  gm = p7_profile_Create(hmm->M, abc); 
  p7_profile_Config(gm, hmm, bg);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  /* p7_oprofile_Dump(stdout, om);  */

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_bg_SetLength(bg,            sq->n);

      if (!ox) ox = p7_filtermx_Create(gm->M, sq->n, ESL_MBYTES(32));  
      else          p7_filtermx_GrowTo(ox, om->M, sq->n); 

      if (!gx) gx = p7_gmx_Create     (gm->M, sq->n);
      else          p7_gmx_GrowTo     (gx, gm->M, sq->n); 

      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);
    
      p7_ForwardFilter(sq->dsq, sq->n, om, ox, &fraw);
      p7_GForward     (sq->dsq, sq->n, gm, gx, &gfraw);

      fsc  =  (fraw-nullsc) / eslCONST_LOG2;
      gfsc = (gfraw-nullsc) / eslCONST_LOG2;
      P  = esl_exp_surv(fsc,   om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
      gP = esl_exp_surv(gfsc,  gm->evparam[p7_FTAU],  gm->evparam[p7_FLAMBDA]);

      gmem = (float) p7_gmx_Sizeof(gx) / 1000000.;
      cmem = (float) p7_filtermx_Sizeof(ox) / 1000000.;

      printf("%-30s\t%-20s\t%9.2g\t%6.1f\t%9.2g\t%6.1f\t%6.2fM\t%6.2fM\n", sq->name, hmm->name, P, fsc, gP, gfsc, gmem, cmem);

      esl_sq_Reuse(sq);
      p7_gmx_Reuse(gx);
      p7_filtermx_Reuse(ox);
    }

  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_filtermx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7FWDFILTER_EXAMPLE*/
/*---------------------- end, example ---------------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
                                          
