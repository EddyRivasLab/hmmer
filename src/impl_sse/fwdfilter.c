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
 *    1. Forward and Backward API calls.
 *    2. Internal functions: inlined fwd/bck/decoding recursions.
 *    x. Copyright and license information.
 */
#include "p7_config.h"

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "p7_gbands.h"

#include "impl_sse.h"
#include "p7_filtermx.h"


static inline float forward_row(int Q, const P7_OPROFILE *om, const __m128 *rp, const __m128 *dpp, __m128 *dpc);
static inline void  backward_row_main  (const P7_OPROFILE *om, const __m128 *rp, const __m128 *dpp, __m128 *dpc, int Q, float scalefactor);
static inline void  backward_row_L     (const P7_OPROFILE *om,                                      __m128 *dpc, int Q, float scalefactor);
static inline float backward_row_zero  (const P7_OPROFILE *om, const __m128 *rp, const __m128 *dpp, __m128 *dpc, int Q);
static inline void  backward_row_finish(const P7_OPROFILE *om, __m128 xEv, __m128 *dpc, int Q);
static inline void  backward_row_rescale(float *xc, __m128 *dpc, int Q, float scalefactor);
static inline void  posterior_decode_row(int rowi, const __m128 *fwd, const __m128 *bck, int Q, P7_GBANDS *bnd);

static inline void  alt_backward_row_main(const P7_OPROFILE *om, const __m128 *rp, __m128 *dpp, __m128 *dpc, int Q, float scalefactor);
static inline void  alt_backward_row_L(const P7_OPROFILE *om,  __m128 *dpc, int Q, float scalefactor);
static inline float alt_backward_row_zero(const P7_OPROFILE *om, const __m128 *rp, __m128 *dpp, __m128 *dpc, int Q);
static inline void  alt_backward_row_finish(const P7_OPROFILE *om, __m128 *dpc, int Q, __m128 dcv);
static inline void  alt_backward_row_rescale(float *xc, __m128 *dpc, int Q, float scalefactor);

/*****************************************************************
 * 1. Forward and Backward API calls
 *****************************************************************/

/* Function:  p7_ForwardFilter()
 * Synopsis:  Checkpointed striped vector Forward calculation.
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

#if p7_DEBUGGING
  if (ox->debugging) p7_filtermx_DumpFBHeader(ox)
#endif

  /* Initialization of the zero row, including specials */
  for (q = 0; q < p7F_NSCELLS*Q; q++) dpp[q] = zerov;
  xc[p7F_N]     = 1.;
  xc[p7F_B]     = om->xf[p7O_N][p7O_MOVE]; 
  xc[p7F_E]     = xc[p7F_JJ] = xc[p7F_J]  = xc[p7F_CC] = xc[p7F_C]  = 0.;			
  xc[p7F_SCALE] = 1.;			   
#if p7_DEBUGGING
  if (ox->debugging) p7_filtermx_DumpFBRow(ox, 0, dpp, "fwd O");
#endif

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


/* Function:  p7_BackwardFilter()
 * Synopsis:  
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
p7_BackwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, P7_GBANDS *bnd, float *opt_sc)
{
  int Q = P7F_NQF(om->M);
  __m128 *fwd;
  __m128 *bck;
  __m128 *dpp;
  __m128 *rp;
  float  *xf;
  float   totsc = 0.;
  int     i, b, w, i2;
  
  /* Row L is a special case for Backwards; so in checkpointed Backward,
   * we have to special-case the last block, rows L-1 and L
   */
  i = L;
  ox->R--;
  fwd = (__m128 *) ox->dpf[ox->R0 + ox->R];   /* pop row for fwd[L] off the checkpointed stack */
  xf  = (float *) (fwd + Q*p7F_NSCELLS);
  bck = (__m128 *) ox->dpf[i%2];	      /* get tmp space for bck[L]                      */
  alt_backward_row_L(om, bck, Q, xf[p7F_SCALE]);  /* calculate bck[L] row                          */
  posterior_decode_row(i, fwd, bck, Q, bnd);
  totsc += logf(xf[p7F_SCALE]);
  i--;
  dpp = bck;


  /* If there's any checkpointing, there's an L-1 row to fill now. */
  if (ox->Rb+ox->Rc > 0)
    {
      /* Compute fwd[L-1] from last checkpoint, which we know is fwd[L-2] */
      dpp = (__m128 *) ox->dpf[ox->R0+ox->R-1];  /* fwd[L-2] values, already known        */
      fwd = (__m128 *) ox->dpf[ox->R0+ox->R];    /* get free row memory from top of stack */
      rp  = om->rfv[dsq[i]];			 /* emission scores on row i, for fwd     */
      forward_row(Q, om, rp, dpp, fwd);          /* calculate fwd[L-1]                    */

      /* Compute bck[L-1] from bck[L]. */
      xf  = (float *) (fwd + Q*p7F_NSCELLS);
      dpp = (__m128 *) ox->dpf[(i+1)%2]; 
      bck = (__m128 *) ox->dpf[i%2];             /* get space for bck[L-1]                */
      rp  = om->rfv[dsq[i+1]];			 /* emission scores on row i+1, for bck   */
      alt_backward_row_main(om, rp, dpp, bck, Q, xf[p7F_SCALE]);

      /* And decode. */
      posterior_decode_row(i, fwd, bck, Q, bnd);
      totsc += logf(xf[p7F_SCALE]);
      dpp = bck;
      i--;			/* i is now L-2 if there's checkpointing; else it's L-1 */
    }

  /* Main loop for checkpointed regions (b,c) */
   for (b = 2; b <= ox->Rb+ox->Rc; b++)
    {				/* i=L-2 as we enter here, and <dpp> is on bck[L-1] */
      w = (b <= ox->Rc ? b+1 : ox->Lb);

      /* We know current row i (r=R0+R-1) ends a block and is checkpointed in fwd. */
      ox->R--;
      fwd = (__m128 *) ox->dpf[ox->R0+ox->R];      /* pop checkpointed forward row off "stack" */
      xf  = (float *) (fwd + Q*p7F_NSCELLS);

      /* Calculate bck[i]; <dpp> is already bck[i+1] */
      bck = (__m128 *) ox->dpf[i%2];	    /* get available tmp memory for row     */
      rp  = om->rfv[dsq[i+1]];	            /* emission scores for row i+1, for bck */
      alt_backward_row_main(om, rp, dpp, bck, Q, xf[p7F_SCALE]);
      totsc += logf(xf[p7F_SCALE]);

      /* And decode checkpointed row i. */
      posterior_decode_row(i, fwd, bck, Q, bnd);
      
      /* The rest of the rows in the block weren't checkpointed.
       * Compute Forwards from last checkpoint ...
       */
      dpp = (__m128 *) ox->dpf[ox->R0+ox->R-1];       /* get last Fwd checkpoint. */
      for (i2 = i-w+1; i2 <= i-1; i2++)
	{
	  fwd = (__m128 *) ox->dpf[ox->R0+ox->R]; ox->R++;  /* push new forward row on "stack"     */
	  rp  = om->rfv[dsq[i2]];		            /* emission scores for row i2, for fwd */
	  forward_row(Q, om, rp, dpp, fwd);
	  dpp = fwd;	  
	}

      /* ... and compute Backwards over the block we just calculated, while decoding. */
      dpp = bck;
      for (i2 = i-1; i2 >= i-w+1; i2--)
	{
	  ox->R--;
	  fwd = (__m128 *) ox->dpf[ox->R0+ox->R]; /* pop just-calculated forward row i2 off "stack" */
	  xf  = (float *) (fwd + Q*p7F_NSCELLS);
	  bck = (__m128 *) ox->dpf[i2%2];	  /* get available for calculating bck[i2]          */
	  rp  = om->rfv[dsq[i2+1]];               /* emission scores for row i2+1, for bck          */
	  alt_backward_row_main(om, rp, dpp, bck, Q, xf[p7F_SCALE]);
	  posterior_decode_row(i2, fwd, bck, Q, bnd);
	  totsc += logf(xf[p7F_SCALE]);
	  dpp = bck;
	}
      i -= w;
    }
   /* now i=La as we leave the checkpointed regions; or i=L-1 if there was no checkpointing */
   

   /* The uncheckpointed "a" region */
   for (; i >= 1; i--)
     {
       ox->R--; 
       fwd = (__m128 *) ox->dpf[ox->R0+ox->R]; /* pop off calculated row fwd[i]           */
       xf  = (float *) (fwd + Q*p7F_NSCELLS);
       bck = (__m128 *) ox->dpf[i%2];	       /* get open space for bck[i]               */
       rp  = om->rfv[dsq[i+1]];		       /* emission scores for row i+1, for bck[i] */
       alt_backward_row_main(om, rp, dpp, bck, Q, xf[p7F_SCALE]);
       posterior_decode_row(i, fwd, bck, Q, bnd);
       totsc += logf(xf[p7F_SCALE]);
       dpp = bck;
    }

   /* To get the backward score, we need to complete row 0 too, where
    * only the N and B states are reachable: B from Mk, and N from B.
    * Seeing the backward score match the forward score is useful in
    * debugging. But to do posterior decoding, which is all we need in 
    * production code, we can stop at row 1.
    */
   bck = (__m128 *) ox->dpf[0];
   rp  = om->rfv[dsq[1]];
   totsc += alt_backward_row_zero(om, rp, dpp, bck, Q);

   if (opt_sc) *opt_sc = totsc;
   return eslOK;
}
/*----------- end forward/backward API calls --------------------*/



/*****************************************************************
 * 2. Internal functions: inlined recursions
 *****************************************************************/

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


static inline void
alt_backward_row_main(const P7_OPROFILE *om, const __m128 *rp, __m128 *dpp, __m128 *dpc, int Q, float scalefactor)
{
  float       * const xc = (float *) (dpc + Q * p7F_NSCELLS); /* E N JJ J B CC C SCALE */
  const float * const xp = (float *) (dpp + Q * p7F_NSCELLS);
  const __m128 *tp, *tpdd;
  const __m128  zerov = _mm_set1_ps(0.0f);
  __m128        xBv   = zerov;
  __m128       *dp;
  __m128        xEv;
  __m128        dcv, mcv, ipv, mpv;
  int           q;
  __m128 tmmv, timv, tdmv;				   /* copies of transition prob quads; a leftshift is needed as boundary cond */
 
  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. */
  dp  = dpp;
  for (q = 0; q < Q; q++)
    {
      *dp = _mm_mul_ps(*dp, *rp); rp++;
      xBv = _mm_add_ps(xBv, *dp); dp+= p7F_NSCELLS;
    }

  /* Specials. Dependencies dictate partial order C,CC,B < N,J,JJ < E */
  xc[p7F_C] = xc[p7F_CC] = xp[p7F_C] * om->xf[p7O_C][p7O_LOOP];
  esl_sse_hsum_ps(xBv, &(xc[p7F_B]));
  xc[p7F_J] = xc[p7F_JJ] = xc[p7F_B] * om->xf[p7O_J][p7O_MOVE] + xp[p7F_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7F_N]              = xc[p7F_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7F_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7F_E]              = xc[p7F_C] * om->xf[p7O_E][p7O_MOVE] + xc[p7F_J] * om->xf[p7O_E][p7O_LOOP];

  /* Initialize for the row calculation */
  mpv  = esl_sse_leftshift_ps(*dpp,       zerov); /* [1 5 9 13] -> [5 9 13 x], M(i+1,k+1) * e(M_k+1, x_{i+1}) */
  tmmv = esl_sse_leftshift_ps(om->tfv[1], zerov);
  timv = esl_sse_leftshift_ps(om->tfv[2], zerov);
  tdmv = esl_sse_leftshift_ps(om->tfv[3], zerov);
  xEv  = _mm_set1_ps(xc[p7F_E]);
  tp   = om->tfv + 7*Q - 1;
  tpdd = tp + Q;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--)
    {
      ipv                 = P7F_IQ(dpp, q);
      P7F_IQ(dpc,q)       = _mm_add_ps( _mm_mul_ps(ipv, *tp),   _mm_mul_ps(mpv, timv)); tp--;   /* II,IM; I is done         */
      mcv                 = _mm_add_ps( _mm_mul_ps(ipv, *tp),   _mm_mul_ps(mpv, tmmv)); tp-=2;  /* MI,MM; ME,MD remain      */
      dcv                 = _mm_add_ps( _mm_mul_ps(dcv, *tpdd), _mm_mul_ps(mpv, tdmv)); tpdd--; /* DM and one segment of DD */

      P7F_DQ(dpc,q) = dcv = _mm_add_ps( xEv, dcv);
      P7F_MQ(dpc,q)       = _mm_add_ps( xEv, mcv);

      mpv  = P7F_MQ(dpp, q);
      tdmv = *tp; tp--;
      timv = *tp; tp--;
      tmmv = *tp; tp-=2;
    }
  
  alt_backward_row_finish(om, dpc, Q, dcv);
  alt_backward_row_rescale(xc, dpc, Q, scalefactor);
}

static inline void
alt_backward_row_finish(const P7_OPROFILE *om, __m128 *dpc, int Q, __m128 dcv)
{
  const __m128 zerov = _mm_setzero_ps();
  const __m128 *tp;
  __m128 *dp;
  int           j,q;
  
  if (om->M < 100)
    {
      for (j = 1; j < 4; j++)
	{
	  dcv = esl_sse_leftshift_ps(dcv, zerov);  
	  tp  = om->tfv + 8*Q - 1;	/* <*tp> now the [4 8 12 x] TDD quad */
	  dp  = dpc + Q*p7F_NSCELLS - 2; /* init to point at D(i,q) vector  */
	  for (q = Q-1; q >= 0; q--)
	    {
	      dcv = _mm_mul_ps(dcv, *tp); tp--;
	      *dp = _mm_add_ps(*dp, dcv); dp -= p7F_NSCELLS;
	    }
	}
    }
  else
    {
      __m128 sv;
      __m128 cv;			/* keeps track of whether any DD addition changes DQ(q) value */
      for (j = 1; j < 4; j++)
	{
	  dcv = esl_sse_leftshift_ps(dcv, zerov);
	  tp  = om->tfv + 8*Q - 1;	/* <*tp> now the [4 8 12 x] TDD quad */
	  dp  = dpc + Q*p7F_NSCELLS - 2;
	  cv  = zerov;
	  for (q = Q-1; q >= 0; q--)
	    { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
	      dcv  = _mm_mul_ps(dcv, *tp); tp--;
	      sv   = _mm_add_ps(*dp, dcv);
	      cv   = _mm_or_ps(cv, _mm_cmpgt_ps(sv, *dp)); /* if DD path changed DQ(dpc,q), cv knows it now */
	      *dp  = sv; 
	      dp  -= p7F_NSCELLS;
	    }
	  if (! _mm_movemask_ps(cv)) break; /* if no DD path changed DQ(q) in this segment, then done, no more segments needed */
	}
    }

  /* M->D paths */
  dcv = esl_sse_leftshift_ps(P7F_DQ(dpc, 0), zerov);
  tp  = om->tfv + 7*Q - 3;	 /* <*tp> is now the [4 8 12 x] tMk->Dk+1 quad */
  dp  = dpc + (Q-1)*p7F_NSCELLS; /* <*dp> is now on M(i,q) vector */
  for (q = Q-1; q >= 0; q--)
    {
      *dp  = _mm_add_ps(*dp, _mm_mul_ps(dcv, *tp)); tp -= 7; 
      dcv  = *(dp+1);                               dp -= p7F_NSCELLS;
    }
}

static inline void
alt_backward_row_L(const P7_OPROFILE *om,  __m128 *dpc, int Q, float scalefactor)
{
  const __m128  zerov = _mm_setzero_ps();
  float        *xc    = (float *) (dpc + Q * p7F_NSCELLS);
  const __m128 *tpdd  = om->tfv + 8*Q - 1;
  __m128       *dp;
  __m128       xEv, dcv;
  int          q;

  /* Backwards from T <- C,CC <- E;  all other specials unreachable, impossible on row L.
   * specials are stored in order E N JJ J B CC C.  
   */
  xc[p7F_C] = xc[p7F_CC] = om->xf[p7O_C][p7O_MOVE];
  xc[p7F_B] = xc[p7F_J] = xc[p7F_JJ] = xc[p7F_N] = 0.0;
  xc[p7F_E] = xc[p7F_C] * om->xf[p7O_E][p7O_MOVE];

  xEv = _mm_set1_ps(xc[p7F_E]);
  dp = dpc;
  for (q = Q-1; q >= 0; q--) 
    {
      *dp++       = xEv; 	                                      /* M */
      *dp++ = dcv = _mm_add_ps(xEv, _mm_mul_ps(dcv, *tpdd)); tpdd--;  /* D */
      *dp++       = zerov;		                              /* I */
    }
  alt_backward_row_finish(om, dpc, Q, dcv);
  alt_backward_row_rescale(xc, dpc, Q, scalefactor);
}


static inline float
alt_backward_row_zero(const P7_OPROFILE *om, const __m128 *rp, __m128 *dpp, __m128 *dpc, int Q)
{
  const __m128 zerov = _mm_setzero_ps();
  float        *xc   = (float *) (dpc + Q * p7F_NSCELLS); /* special states on current row i      */
  float        *xp   = (float *) (dpp + Q * p7F_NSCELLS); /* special states on "previous" row i+1 */
  __m128       *dp;
  __m128        xBv  = zerov;
  int           q;

  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. */
  dp  = dpp;
  for (q = 0; q < Q; q++)
    {
      *dp = _mm_mul_ps(*dp, *rp); rp++;
      xBv = _mm_add_ps(xBv, *dp); dp+= p7F_NSCELLS;
    }

  xc[p7F_C] = xc[p7F_CC] = 0.0f;
  esl_sse_hsum_ps(xBv, &(xc[p7F_B]));
  xc[p7F_J] = xc[p7F_JJ] = 0.0f;
  xc[p7F_N]              = xc[p7F_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7F_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7F_E]              = 0.0f;
  
  /* not needed in production code: */
  for (q = 0; q < Q; q++)
    P7F_MQ(dpc, q) = P7F_IQ(dpc, q) = P7F_DQ(dpc, q) = zerov;

  return logf(xc[p7F_N]);
}


static inline void
alt_backward_row_rescale(float *xc, __m128 *dpc, int Q, float scalefactor)
{
  if (scalefactor > 1.0f)
    {
      __m128  sv = _mm_set1_ps(1.0 / scalefactor);
      __m128 *dp = dpc;
      int     q;

      xc[p7F_E]  /= scalefactor;
      xc[p7F_N]  /= scalefactor;
      xc[p7F_JJ] /= scalefactor;
      xc[p7F_J]  /= scalefactor;
      xc[p7F_B]  /= scalefactor;
      xc[p7F_CC] /= scalefactor;
      xc[p7F_C]  /= scalefactor;

      for (q = 0; q < Q; q++) 
	{
	  *dp = _mm_mul_ps(*dp, sv); dp++; /* M */
	  *dp = _mm_mul_ps(*dp, sv); dp++; /* D */
	  *dp = _mm_mul_ps(*dp, sv); dp++; /* I */
	}
    }
  xc[p7F_SCALE] = scalefactor;
}


static inline void
backward_row_main(const P7_OPROFILE *om, const __m128 *rp, const __m128 *dpp, __m128 *dpc, int Q, float scalefactor)
{
  const __m128  zerov = _mm_setzero_ps();
  float        *xc    = (float *) (dpc + Q * p7F_NSCELLS); /* special states on current row i      */
  float        *xp    = (float *) (dpp + Q * p7F_NSCELLS); /* special states on "previous" row i+1 */
  const __m128 *tp;					   /* <tp> will step backwards thru transition prob quads */
  __m128 mcv;					           /* M value for current quad */
  __m128 mpv, ipv;					   /* MI values for quads on i+1 row we're connecting to */
  __m128 tmmv, timv, tdmv;				   /* copies of transition prob quads; a leftshift is needed as boundary cond */
  __m128 xBv, xEv;					   /* splatted B, E values */
  int    q;
  
  mpv = _mm_mul_ps(P7F_MQ(dpp,0), rp[0]); /* precalc M(i+1,k+1) * e(M_k+1, x_{i+1}) */
  mpv = esl_sse_leftshift_ps(mpv, zerov); /* [1 5 9 13] -> [5 9 13 x], M(i+1,k+1) * e(M_k+1, x_{i+1}) */

  tmmv = esl_sse_leftshift_ps(om->tfv[1], zerov);
  timv = esl_sse_leftshift_ps(om->tfv[2], zerov);
  tdmv = esl_sse_leftshift_ps(om->tfv[3], zerov);

  rp  += Q-1;
  xBv = zerov;
  tp  = om->tfv + 7*Q - 1; /* <*tp> is now the [4 8 12 x] TII transition quad, the last one  */
  for (q = Q-1; q >= 0; q--)     /* backwards stride */
    {
      ipv = P7F_IQ(dpp,q); /* assumes emission odds ratio of 1.0; i+1's IMO(q) now free */
      P7F_IQ(dpc,q) = _mm_add_ps(_mm_mul_ps(ipv, *tp), _mm_mul_ps(mpv, timv));   tp--;
      P7F_DQ(dpc,q) =                                  _mm_mul_ps(mpv, tdmv); 
      mcv           = _mm_add_ps(_mm_mul_ps(ipv, *tp), _mm_mul_ps(mpv, tmmv));   tp-= 2;
	  
      mpv           = _mm_mul_ps(P7F_MQ(dpp,q), *rp);  rp--;  /* obtain mpv for next q. i+1's MMO(q) is freed  */
      P7F_MQ(dpc,q) = mcv;

      tdmv = *tp--;
      timv = *tp--;
      tmmv = *tp--;

      xBv = _mm_add_ps(xBv, _mm_mul_ps(mpv, *tp)); tp--;
    }

  xc[p7F_C] = xc[p7F_CC] = xp[p7F_C] * om->xf[p7O_C][p7O_LOOP];
  esl_sse_hsum_ps(xBv, &(xc[p7F_B])); /* collect the B values together: horizontal sum in xBv vector  */
  xc[p7F_J] = xc[p7F_JJ] = xc[p7F_B] * om->xf[p7O_J][p7O_MOVE] + xp[p7F_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7F_N]              = xc[p7F_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7F_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7F_E]              = xc[p7F_C] * om->xf[p7O_E][p7O_MOVE] + xc[p7F_J] * om->xf[p7O_E][p7O_LOOP];
  
  xEv = _mm_set1_ps(xc[p7F_E]);
  backward_row_finish(om, xEv, dpc, Q);
  backward_row_rescale(xc, dpc, Q, scalefactor);
}
  

/* backward_row_L()
 * 
 * The last row in the backwards matrix is a special case,
 * with no contribution from any next row. We don't explicitly
 * represent the terminal T state or the C->T transition, so
 * it's easier to break out this special case than to try 
 * to initialize an L+1 row and use the normal backwards
 * recursion.
 */
static inline void
backward_row_L(const P7_OPROFILE *om,  __m128 *dpc, int Q, float scalefactor)
{
  const __m128 zerov = _mm_setzero_ps();
  float       *xc    = (float *) (dpc + Q * p7F_NSCELLS);
  __m128       xEv;
  int          q;

  for (q = 0; q < Q; q++) 
    P7F_MQ(dpc, q) = P7F_DQ(dpc, q) = P7F_IQ(dpc, q) = zerov; 

  /* Backwards from T <- C,CC <- E;  all other specials unreachable, impossible on row L.
   * specials are stored in order E N JJ J B CC C.  
   */
  xc[p7F_C] = xc[p7F_CC] = om->xf[p7O_C][p7O_MOVE];
  xc[p7F_B] = xc[p7F_J] = xc[p7F_JJ] = xc[p7F_N] = 0.0;
  xc[p7F_E] = xc[p7F_C] * om->xf[p7O_E][p7O_MOVE];
  xEv = _mm_set1_ps(xc[p7F_E]);

  backward_row_finish(om, xEv, dpc, Q);
  backward_row_rescale(xc, dpc, Q, scalefactor);
}
  
static inline float
backward_row_zero(const P7_OPROFILE *om, const __m128 *rp, const __m128 *dpp, __m128 *dpc, int Q)
{
  const __m128 zerov = _mm_setzero_ps();
  float        *xc   = (float *) (dpc + Q * p7F_NSCELLS); /* special states on current row i      */
  float        *xp   = (float *) (dpp + Q * p7F_NSCELLS); /* special states on "previous" row i+1 */
  __m128        xBv  = zerov;
  const __m128 *tp   = om->tfv;
  __m128        mpv;
  int           q;

  for (q = 0; q < Q; q++)
    {
      P7F_MQ(dpc, q) = P7F_IQ(dpc, q) = P7F_DQ(dpc, q) = zerov;

      mpv = _mm_mul_ps( P7F_MQ(dpp, q), *rp); rp++;
      mpv = _mm_mul_ps( mpv,            *tp); tp+=7;
      xBv = _mm_add_ps( xBv, mpv);
    }

   xc[p7F_C] = xc[p7F_CC] = 0.0f;
   esl_sse_hsum_ps(xBv, &(xc[p7F_B]));
   xc[p7F_J] = xc[p7F_JJ] = 0.0f;
   xc[p7F_N]              = xc[p7F_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7F_N] * om->xf[p7O_N][p7O_LOOP];
   xc[p7F_E]              = 0.0f;

   return logf(xc[p7F_N]);
}
  
/* backward_row_finish()
 * 
 * Given: current backward row <dpc> of Q vectors that already
 * includes any contributions from next row i+1; 
 * and given xEv vector containing splatted E value on current row;
 * finish the backward calculation on the row, which means:
 *   adding {MD}k->E paths;
 *   adding DD paths (serialized, fully or partially);
 *   and adding the Mk->Dk+1 paths.
 */
static inline void
backward_row_finish(const P7_OPROFILE *om, __m128 xEv, __m128 *dpc, int Q)
{
  const __m128 zerov = _mm_setzero_ps();
  const __m128 *tp;
  __m128 dpv, dcv;
  int    j,q;

  /* First DD pass adds {MD}k->E path, and one segment of Dk->Dk+1 paths   */
  tp  = om->tfv + 8*Q - 1;	           /* <*tp> now the [4 8 12 x] TDD quad */
  dpv = _mm_add_ps(P7F_DQ(dpc,0), xEv);	   /* <dpv>: add Dk->E piece to D[1 5 9 13] vector */
  dpv = esl_sse_leftshift_ps(dpv, zerov);  /* [1 5 9 13] -> [5 9 13 x] */
  for (q = Q-1; q >= 0; q--)
    {
      dcv           = _mm_mul_ps(dpv, *tp); tp--;
      P7F_DQ(dpc,q) = _mm_add_ps(P7F_DQ(dpc,q), _mm_add_ps(dcv, xEv));
      dpv           = P7F_DQ(dpc,q);
      P7F_MQ(dpc,q) = _mm_add_ps(P7F_MQ(dpc,q), xEv);
    }
  
  /* DD paths have to be at least partially serialized, as in Forward;
   * see comments there.  In short, we have two options, a fully
   * serialized option, and a partially serialized option that has
   * more overhead. They trade off at around M~100 in query length,
   * empirically.
   */
  if (om->M < 100)
    {				/* full serialization */
      for (j = 1; j < 4; j++)	/* three passes: we've already done 1 segment, we need 4 total */
	{
	  dcv = esl_sse_leftshift_ps(dcv, zerov);  
	  tp  = om->tfv + 8*Q - 1;	/* <*tp> now the [4 8 12 x] TDD quad */
	  for (q = Q-1; q >= 0; q--)
	    {
	      dcv           = _mm_mul_ps(dcv, *tp); tp--;
	      P7F_DQ(dpc,q) = _mm_add_ps(P7F_DQ(dpc,q), dcv);
	    }
	}
    }
  else 
    {  /* Slightly parallelized version.  */
      __m128 sv;
      __m128 cv;			/* keeps track of whether any DD addition changes DQ(q) value */
      for (j = 1; j < 4; j++)
	{
	  dcv = esl_sse_leftshift_ps(dcv, zerov);
	  tp  = om->tfv + 8*Q - 1;	/* <*tp> now the [4 8 12 x] TDD quad */
	  cv  = zerov;
	  for (q = Q-1; q >= 0; q--)
	    { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
	      dcv           = _mm_mul_ps(dcv, *tp); tp--;
	      sv            = _mm_add_ps(P7F_DQ(dpc,q), dcv);
	      cv            = _mm_or_ps(cv, _mm_cmpgt_ps(sv, P7F_DQ(dpc,q))); /* if DD path changed DQ(dpc,q), cv knows it now */
	      P7F_DQ(dpc,q) = sv;
	    }
	  if (! _mm_movemask_ps(cv)) break; /* if no DD path changed DQ(q) in this segment, then done, no more segments needed */
	}
    }
 
  /* Finally, we can add the Mk<-Dk+1 path */   
  tp  = om->tfv + 7*Q - 3;	                        /* <*tp> now the [4 8 12 x] Mk->Dk+1 quad    */
  dcv = esl_sse_leftshift_ps(P7F_DQ(dpc, 0), zerov);	/* leftshift: [1 5 9 13] -> [5 9 13 x]       */
  for (q = Q-1; q >= 0; q--)
    {
      P7F_MQ(dpc,q) = _mm_add_ps(P7F_MQ(dpc,q), _mm_mul_ps(dcv, *tp)); tp -= 7;
      dcv           = P7F_DQ(dpc,q);
    }
  /* Calculation of the backwards row is now complete, ready for rescaling */
}
  

/* backward_row_rescale()
 * Rescale the entire row (dpc[0..Q-1] vectors, xc floats) using
 * the scalefactor that Forward set, which caller provides.
 */
static inline void
backward_row_rescale(float *xc, __m128 *dpc, int Q, float scalefactor)
{
  __m128 sv;
  int    q;

  if (scalefactor > 1.0f)
    {
      xc[p7F_E]  /= scalefactor;
      xc[p7F_N]  /= scalefactor;
      xc[p7F_JJ] /= scalefactor;
      xc[p7F_J]  /= scalefactor;
      xc[p7F_B]  /= scalefactor;
      xc[p7F_CC] /= scalefactor;
      xc[p7F_C]  /= scalefactor;

      sv = _mm_set1_ps(1.0 / scalefactor);
      for (q = 0; q < Q; q++) {
	P7F_MQ(dpc,q) = _mm_mul_ps(P7F_MQ(dpc,q), sv);
	P7F_DQ(dpc,q) = _mm_mul_ps(P7F_DQ(dpc,q), sv);
	P7F_IQ(dpc,q) = _mm_mul_ps(P7F_IQ(dpc,q), sv);
      }
    }
  xc[p7F_SCALE] = scalefactor;
}


static inline void
posterior_decode_row(int rowi, const __m128 *fwd, const __m128 *bck, int Q, P7_GBANDS *bnd)
{
  const  float *xf = (float *) (fwd + Q*p7F_NSCELLS);
  const  float *xb = (float *) (bck + Q*p7F_NSCELLS);
  __m128 threshv   = _mm_set1_ps(0.02);
  __m128 onev      = _mm_set1_ps(1.0);
  __m128 kv        = _mm_set_ps((float) 3*Q, (float) 2*Q, (float) Q, 0.0f); /* iv is going to contain current i coord in each vector cell */
  __m128 minkv     = _mm_set1_ps( (float) Q*4.0+1);
  __m128 maxkv     = _mm_set1_ps( 0.0f);
  __m128 mask;
  __m128 pv;
  int    q;
  float  ktmp;
  int    ka, kb;
  
  if (xf[p7F_N] * xb[p7F_N] + xf[p7F_JJ] * xb[p7F_JJ] + xf[p7F_CC] * xb[p7F_CC] >= 0.9) return;

  /* ka = min k that satisfies threshold posterior prob; kb = max k.
   * Because of striping, identifying ka..kb is slightly nonobvious.
   * What we do is collect min_k, max_k in each vector cell; these are
   * the min/max in each segment (1..Q, Q+1..2Q, 2Q+1..3Q, 3Q+1..4Q).
   * These are in <minkv>, <maxkv>. Finally we do a horizontal min/max
   * on these vectors, and cast the result back to the integer coords
   * we want.
   * 
   * kv vector contains the current i indices for each cell; we
   * increment it by one (<onev>) each loop.
   *
   * All this coord computation is done as floats, because SSE2
   * doesn't give us the integer _max/min we want. A float has exact
   * range up to 2^24: 16.77M, more than enough for our design limit
   * of L<=100K.
   */
  for (q = 0; q < Q; q++)
    {
      pv   =                _mm_mul_ps(P7F_MQ(fwd, q), P7F_MQ(bck, q));
      pv   = _mm_add_ps(pv, _mm_mul_ps(P7F_IQ(fwd, q), P7F_IQ(bck, q)));
      mask = _mm_cmpge_ps(pv, threshv);

      kv    = _mm_add_ps(kv, onev);
      minkv = _mm_min_ps( minkv, _mm_and_ps(kv, mask));
      maxkv = _mm_max_ps( maxkv, _mm_and_ps(kv, mask));
    }
  esl_sse_hmax_ps(maxkv, &ktmp); kb = (int) ktmp;
  esl_sse_hmin_ps(minkv, &ktmp); ka = (int) ktmp;

  if (ka) p7_gbands_Prepend(bnd, rowi, ka, kb);
}
/*------------------ end, inlined recursions -------------------*/





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
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Forward",                         0 },
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
  P7_GBANDS      *bnd     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, bench_time, Mcs;

  impl_Init();

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
  bnd = p7_gbands_Create();

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
      if (! esl_opt_GetBoolean(go, "-F")) 
	p7_BackwardFilter(dsq, L, om, ox, bnd, NULL);

      p7_filtermx_Reuse(ox);
      p7_gbands_Reuse(bnd);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_filtermx_Destroy(ox);
  p7_gbands_Destroy(bnd);
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
  P7_GBANDS      *bnd     = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fraw, braw, nullsc, fsc, bsc;
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

  bnd = p7_gbands_Create();

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
    
      p7_ForwardFilter (sq->dsq, sq->n, om, ox, &fraw);
      p7_BackwardFilter(sq->dsq, sq->n, om, ox, bnd, &braw);

      p7_GForward     (sq->dsq, sq->n, gm, gx, &gfraw);

      fsc  =  (fraw-nullsc) / eslCONST_LOG2;
      bsc  =  (braw-nullsc) / eslCONST_LOG2;
      gfsc = (gfraw-nullsc) / eslCONST_LOG2;
      P  = esl_exp_surv(fsc,   om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
      gP = esl_exp_surv(gfsc,  gm->evparam[p7_FTAU],  gm->evparam[p7_FLAMBDA]);

      gmem = (float) p7_gmx_Sizeof(gx) / 1000000.;
      cmem = (float) p7_filtermx_Sizeof(ox) / 1000000.;

      printf("%-30s\t%-20s\t%9.2g\t%6.1f\t%6.1f\t%9.2g\t%6.1f\t%6.2fM\t%6.2fM\n", sq->name, hmm->name, P, fsc, bsc, gP, gfsc, gmem, cmem);

      esl_sq_Reuse(sq);
      p7_gmx_Reuse(gx);
      p7_gbands_Reuse(bnd);
      p7_filtermx_Reuse(ox);
    }

  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_gbands_Destroy(bnd);
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
                                          
