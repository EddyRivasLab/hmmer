/* Forward/Backward filters; x86 SSE version
 * 
 * See fbfilter.md for notes.
 *
 * This file is conditionally compiled when eslENABLE_SSE4 is defined.
 *
 * Contents:
 *    1. Forward and Backward filter: SSE implementations.
 *    2. Internal functions: inlined recursions.
 *    3. Internal debugging tools.
 */
#include <h4_config.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "h4_checkptmx.h"
#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_refmx.h"
#include "h4_sparsemask.h"

#include "fbfilter.h"

#ifdef eslENABLE_SSE4
#include <x86intrin.h>
#include "esl_sse.h"


/* Vectorized DP recursions are tediously lengthy. For some semblance
 * of clarity, we break them out into one-page-ish chunks, using
 * static inlined functions.
 */
static inline float forward_row_sse      (ESL_DSQ xi, const H4_PROFILE *hmm, const H4_MODE *mo, const __m128 *dpp, __m128 *dpc, int Q);
static inline void  backward_row_main_sse(ESL_DSQ xi, const H4_PROFILE *hmm, const H4_MODE *mo,       __m128 *dpp, __m128 *dpc, int Q, float scalefactor);
static inline void  backward_row_L_sse   (            const H4_PROFILE *hmm, const H4_MODE *mo,                    __m128 *dpc, int Q, float scalefactor);
static inline void  backward_row_finish_sse(const H4_PROFILE *hmm, __m128 *dpc, int Q, __m128 dcv);
static inline void  backward_row_rescale_sse(float *xc, __m128 *dpc, int Q, float scalefactor);
static inline int   posterior_decode_row_sse(H4_CHECKPTMX *cpx, int rowi, H4_SPARSEMASK *sm, float sm_thresh, float overall_sc);


#if eslDEBUGLEVEL > 0
static inline float backward_row_zero_sse(ESL_DSQ x1, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx);
static        void  save_debug_row_pp_sse(H4_CHECKPTMX *cpx, __m128 *dpc, int i);
static        void  save_debug_row_fb_sse(H4_CHECKPTMX *cpx, H4_REFMX *gx, __m128 *dpc, int i, float totscale);
#endif // eslDEBUGLEVEL > 0




/* Function:  h4_fwdfilter_sse()
 * Synopsis:  SSE implementation of the Forward filter
 * Incept:    SRE, Mon 22 Jul 2019 (Albany NY)
 *
 * Purpose:   Forward filter, SSE version. Compare digital sequence
 *            <dsq> of length <L> to profile <hmm> using alignment
 *            mode <mo>. Use or reuse <cpx> as the DP matrix;
 *            reallocate and reinitialize <cpx> as needed (caller
 *            doesn't need to worry about this itself).
 *
 *            In general <mo> should be in multihit local mode and
 *            configured by the caller for length <L> by a call to
 *            <h4_mode_SetLength(mo, L)>. However, alternatively, <mo>
 *            can be in dual glocal/local mode or even glocal mode if
 *            that's convenient; the FB filter ignores the local vs
 *            glocal configuration and always treats <mo> as if it's
 *            in local-only mode. Multihit local alignment mode is a
 *            requirement of the numerics of sparse rescaling, to
 *            avoid underflow [Eddy11].
 *
 *            Return the raw Forward score, in bits, in <*opt_sc>. The
 *            checkpointed DP matrix in <cpx> is filled, and ready for
 *            a call to <h4_bckfilter_sse()>.
 */
int
h4_fwdfilter_sse(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc)
{
  __m128       *dpp   = NULL;                     // dpp=prev row. start on dpf[2]; rows 0,1=Backwards
  __m128       *dpc   = NULL;                     // dpc points at current row
  const __m128  zerov = _mm_setzero_ps();     
  float        *xc    = NULL;                     // specials E,N,JJ,J,B,CC,C,SCALE    
  float         totsc = 0.0f;                     // accumulates Forward score in bits
  int Q;                                          // segment length; # of MID vectors on each row
  int q;                                          // counter over vectors 0..Q-1
  int i;                                          // counter over residues/rows 1..L 
  int b;                                          // counter down through checkpointed blocks, Rb+Rc..1
  int w;                                          // counter down through rows in a checkpointed block
#if eslDEBUGLEVEL > 0
  int do_dumping = (cpx->dfp? TRUE : FALSE);
#endif

  /* First make sure <cpx> is allocated big enough.
   * (DO NOT set any ptrs into the matrix until after this potential reallocation!)
   * Then set the size of the problem in <cpx> immediately, not later,
   * because debugging dumps need this information, for example.
   * (Ditto for any debugging copy of the fwd mx).
   */
  h4_checkptmx_Reinit(cpx, hmm->M, L);
  cpx->M   = hmm->M;       
  cpx->L   = L;
  cpx->Vf  = h4_VWIDTH_SSE / sizeof(float);
  cpx->Q   = Q = H4_Q(hmm->M, cpx->Vf);           // Q, just for shorthand avoidance of cpx->Q everywhere below.
  dpp     = (__m128 *) cpx->dpf[cpx->R0-1];       // dpp=prev row. start on dpf[2]; rows 0,1=Backwards
  xc      = (float *)  (dpp + Q*h4C_NSCELLS);     // specials E,N,JJ,J,B,CC,C,SCALE 
 
#if eslDEBUGLEVEL > 0
  if (do_dumping) h4_checkptmx_DumpFBHeader(cpx);
  if (cpx->fwd) { cpx->fwd->M = hmm->M; cpx->fwd->L = L; cpx->fwd->type = h4R_FORWARD; }
#endif

  /* Initialization of the zero row, including specials */
  for (q = 0; q < h4C_NSCELLS*Q; q++) dpp[q] = zerov;
  xc[h4C_N]     = 1.;
  xc[h4C_B]     = mo->xf[h4_N][h4_MOVE]; 
  xc[h4C_E]     = xc[h4C_JJ] = xc[h4C_J]  = xc[h4C_CC] = xc[h4C_C]  = 0.;                       
  xc[h4C_SCALE] = 1.;

#if eslDEBUGLEVEL > 0
  if (do_dumping) h4_checkptmx_DumpFBRow(cpx, 0, (float *) dpp, "f1 O"); 
  if (cpx->fwd)   save_debug_row_fb_sse(cpx, cpx->fwd, dpp, 0, totsc); 
#endif

  /* Phase one: the "a" region: all rows in this region are saved */
  for (i = 1; i <= cpx->La; i++)
    {   
      dpc = (__m128 *) cpx->dpf[cpx->R0+cpx->R]; cpx->R++;    // idiomatic for "get next save/checkpoint row" 

      totsc += forward_row_sse(dsq[i], hmm, mo, dpp, dpc, Q);
      dpp = dpc;         // current row becomes prev row

#if eslDEBUGLEVEL > 0
      if (do_dumping) h4_checkptmx_DumpFBRow(cpx, i, (float *) dpc, "f1 O"); 
      if (cpx->fwd)   save_debug_row_fb_sse(cpx, cpx->fwd, dpc, i, totsc); 
#endif
    }

  /* Phase two: "b" and "c" regions: partially and fully checkpointed */
  for (b = cpx->Rb + cpx->Rc, w = (cpx->Rb ? cpx->Lb : cpx->Rc+1); i <= L; i++)
    {
      /* this section, deciding whether to get a checkpointed vs. tmp_check
       * row, is why we set <fwd>, <dpp> here, rather than having the
       * inlined forward_row() set them.
       */
      if (! (--w)) {                               // we're on the last row in segment: this row is saved 
        dpc = (__m128 *) cpx->dpf[cpx->R0+cpx->R]; cpx->R++;  // idiomatic for "get next save/checkpoint row" 
        w = b;                                     // next segment has this many rows, ending in a saved row
        b--;                                       // decrement segment number counter; last segment is r=1
      } else      
        dpc = (__m128 *) cpx->dpf[i%2];            // idiomatic for "next tmp row", 0/1; i%2 makes sure dpp != dpc 

      totsc += forward_row_sse(dsq[i], hmm, mo, dpp, dpc, Q);
      dpp = dpc;

#if eslDEBUGLEVEL > 0
      if (do_dumping) h4_checkptmx_DumpFBRow(cpx, i, (float *) dpc, w ? "f1 X" : "f1 O"); 
      if (cpx->fwd)   save_debug_row_fb_sse(cpx, cpx->fwd, dpc, i, totsc); 
#endif
    }

  xc  = (float *) (dpc + Q * h4C_NSCELLS);
  ESL_DASSERT1( (cpx->R == cpx->Ra+cpx->Rb+cpx->Rc) );
  ESL_DASSERT1( ( (! isnan(xc[h4C_C])) && (! isinf(xc[h4C_C]))) );
  if (opt_sc) *opt_sc = totsc + esl_log2f(xc[h4C_C] * mo->xf[h4_C][h4_MOVE]) - mo->nullsc;
  return eslOK;
}



/* Function:  h4_bckfilter_sse()
 * Synopsis:  Backward and decoding steps of the F/B filter.
 * Incept:    SRE, Tue 23 Jul 2019
 *
 * Purpose:   Backwards and posterior decoding, resulting in a sparse
 *            mask. Compare profile <hmm> to sequence <dsq> of length
 *            <L>.  Use DP matrix <cpx>, which already contains the
 *            results of a Forward filter calculation.  Caller
 *            provides an allocated sparse mask structure <sm> of any
 *            size; this structure will be reinitialized and
 *            reallocated as needed. <sm_thresh> is a posterior
 *            decoding probability threshold for including a DP cell
 *            in the sparse mask; if the sum of decodings for M/D/I at
 *            k,i is >= <sm_thresh>, cell <i,k> is included.
 */
int
h4_bckfilter_sse(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh)
{
  int     Q = cpx->Q;     // short/clarity
  float  *xf;
  __m128 *fwd;
  __m128 *bck;
  __m128 *dpp;
  float   Tvalue;
  int     i, b, w, i2;
  int     status;
#if eslDEBUGLEVEL > 0
  int do_dumping = (cpx->dfp? TRUE : FALSE);
#endif

  h4_sparsemask_Reinit(sm, hmm->M, L);

#if eslDEBUGLEVEL > 0
  if (cpx->bck) { cpx->bck->M = hmm->M; cpx->bck->L = L; cpx->bck->type = h4R_BACKWARD; }
  if (cpx->pp)  { cpx->pp->M  = hmm->M; cpx->pp->L  = L; cpx->pp->type  = h4R_DECODING; }
#endif

  /* Row L is a special case for Backwards; so in checkpointed Backward,
   * we have to special-case the last block, rows L-1 and L
   */  
  i = L;
  cpx->R--;
  fwd = (__m128 *) cpx->dpf[cpx->R0 + cpx->R];        // pop row for fwd[L] off the checkpointed stack
  xf  = (float *) (fwd + Q*h4C_NSCELLS);              // set <xf> to the specials on fwd[L] row
  Tvalue = xf[h4C_C] * mo->xf[h4_C][h4_MOVE];         // i.e. scaled fwd[L] val at T state = scaled overall score
  bck = (__m128 *) cpx->dpf[i%2];                     // get tmp space for bck[L]
  backward_row_L_sse(hmm, mo, bck, Q, xf[h4C_SCALE]); // calculate bck[L] row

#if eslDEBUGLEVEL > 0
  cpx->bcksc = esl_log2f(xf[h4C_SCALE]);
  if (do_dumping) { 
    h4_checkptmx_DumpFBRow(cpx, L, (float *) fwd, "f2 O"); 
    h4_checkptmx_DumpFBRow(cpx, L, (float *) bck, "bck");  
  }
  if (cpx->bck) save_debug_row_fb_sse(cpx, cpx->bck, bck, L, cpx->bcksc); 
#endif

  if ( (status = posterior_decode_row_sse(cpx, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
  i--;
  dpp = bck;

  /* If there's any checkpointing, there's an L-1 row to fill now. */
  if (cpx->Rb+cpx->Rc > 0)
    {
      /* Compute fwd[L-1] from last checkpoint, which we know is fwd[L-2] */
      dpp = (__m128 *) cpx->dpf[cpx->R0+cpx->R-1];  // fwd[L-2] values, already known
      fwd = (__m128 *) cpx->dpf[cpx->R0+cpx->R];    // get free row memory from top of stack
      forward_row_sse(dsq[i], hmm, mo, dpp, fwd, Q);    // calculate fwd[L-1]
#if eslDEBUGLEVEL > 0
      if (do_dumping) h4_checkptmx_DumpFBRow(cpx, i, (float *) fwd, "f2 X");
#endif

      /* Compute bck[L-1] from bck[L]. */
      xf  = (float *) (fwd + Q*h4C_NSCELLS);
      dpp = (__m128 *) cpx->dpf[(i+1)%2]; 
      bck = (__m128 *) cpx->dpf[i%2];             // get space for bck[L-1]
      backward_row_main_sse(dsq[i+1], hmm, mo, dpp, bck, Q, xf[h4C_SCALE]);
#if eslDEBUGLEVEL > 0
      cpx->bcksc += esl_log2f(xf[h4C_SCALE]);
      if (do_dumping) h4_checkptmx_DumpFBRow(cpx, i, (float *) bck, "bck");
      if (cpx->bck)   save_debug_row_fb_sse(cpx, cpx->bck, bck, i, cpx->bcksc); 
#endif
      /* And decode. */
      if ( (status = posterior_decode_row_sse(cpx, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
      dpp = bck;
      i--;  // i is now L-2 if there's checkpointing; else it's L-1
    }

  /* Main loop for checkpointed regions (b,c) */
  for (b = 2; b <= cpx->Rb+cpx->Rc; b++)
    { // i=L-2 as we enter here, and <dpp> is on bck[L-1] 
      w = (b <= cpx->Rc ? b+1 : cpx->Lb);

      /* We know current row i (r=R0+R-1) ends a block and is checkpointed in fwd. */
      cpx->R--;
      fwd = (__m128 *) cpx->dpf[cpx->R0+cpx->R];        // pop checkpointed forward row off "stack" 
      xf  = (float *)  (fwd + Q * h4C_NSCELLS);

      /* Calculate bck[i]; <dpp> is already bck[i+1] */
      bck = (__m128 *) cpx->dpf[i%2];        // get available tmp memory for row
      backward_row_main_sse(dsq[i+1], hmm, mo, dpp, bck, Q, xf[h4C_SCALE]);
#if eslDEBUGLEVEL > 0
      cpx->bcksc += esl_log2f(xf[h4C_SCALE]);
      if (do_dumping) { 
        h4_checkptmx_DumpFBRow(cpx, i, (float *) fwd, "f2 O");
        h4_checkptmx_DumpFBRow(cpx, i, (float *) bck, "bck"); 
      }
      if (cpx->bck) save_debug_row_fb_sse(cpx, cpx->bck, bck, i, cpx->bcksc); 
#endif
      /* And decode checkpointed row i. */
      if ( (status = posterior_decode_row_sse(cpx, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
      
      /* The rest of the rows in the block weren't checkpointed.
       * Compute Forwards from last checkpoint ...
       */
      dpp = (__m128 *) cpx->dpf[cpx->R0+cpx->R-1];  // get last Fwd checkpoint. 
      for (i2 = i-w+1; i2 <= i-1; i2++)
        {
          fwd = (__m128 *) cpx->dpf[cpx->R0+cpx->R]; cpx->R++;  // push new forward row on "stack"
          forward_row_sse(dsq[i2], hmm, mo, dpp, fwd, Q);
          dpp = fwd;      
        }

      /* ... and compute Backwards over the block we just calculated, while decoding. */
      dpp = bck;
      for (i2 = i-1; i2 >= i-w+1; i2--)
        {
          cpx->R--;
          fwd = (__m128 *) cpx->dpf[cpx->R0+cpx->R]; // pop just-calculated forward row i2 off "stack"
          xf  = (float *)  (fwd + Q * h4C_NSCELLS);
          bck = (__m128 *) cpx->dpf[i2%2];           // get available for calculating bck[i2] 
          backward_row_main_sse(dsq[i2+1], hmm, mo, dpp, bck, Q, xf[h4C_SCALE]);
#if eslDEBUGLEVEL > 0
          cpx->bcksc += esl_log2f(xf[h4C_SCALE]);
          if (do_dumping) { 
            h4_checkptmx_DumpFBRow(cpx, i2, (float *) fwd, "f2 X"); 
            h4_checkptmx_DumpFBRow(cpx, i2, (float *) bck, "bck"); 
          }
          if (cpx->bck) save_debug_row_fb_sse(cpx, cpx->bck, bck, i2, cpx->bcksc); 
#endif

          if ((status = posterior_decode_row_sse(cpx, i2, sm, sm_thresh, Tvalue)) != eslOK) return status;
          dpp = bck;
        }
      i -= w;
    }
  /* now i=La as we leave the checkpointed regions; or i=L-1 if there was no checkpointing */

  /* the uncheckpointed "a" region */
  for (; i >= 1; i--)
    {
      cpx->R--; 
      fwd = (__m128 *) cpx->dpf[cpx->R0+cpx->R]; // pop off calculated row fwd[i]
      xf  = (float *)  (fwd + Q * h4C_NSCELLS);
      bck = (__m128 *) cpx->dpf[i%2];          // get open space for bck[i]
      backward_row_main_sse(dsq[i+1], hmm, mo, dpp, bck, Q, xf[h4C_SCALE]);

#if eslDEBUGLEVEL > 0
      cpx->bcksc += esl_log2f(xf[h4C_SCALE]);
      if (do_dumping) { 
        h4_checkptmx_DumpFBRow(cpx, i, (float *) fwd, "f2 O"); 
        h4_checkptmx_DumpFBRow(cpx, i, (float *) bck, "bck"); 
      }
      if (cpx->bck) save_debug_row_fb_sse(cpx, cpx->bck, bck, i, cpx->bcksc); 
#endif
      if ((status = posterior_decode_row_sse(cpx, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
      dpp = bck;
    }

#if eslDEBUGLEVEL > 0
  /* To get the backward score, we need to complete row 0 too, where
   * only the N and B states are reachable: B from Mk, and N from B.
   * Seeing the backward score match the forward score is useful in
   * debugging. But to do posterior decoding, which is all we need in 
   * production code, we can stop at row 1.
   */
  float xN;
  cpx->R--;
  fwd = (__m128 *) cpx->dpf[cpx->R0];
  bck = (__m128 *) cpx->dpf[i%2];              
  xN = backward_row_zero_sse(dsq[1], hmm, mo, cpx); 
  if (do_dumping) { 
    h4_checkptmx_DumpFBRow(cpx, 0, (float *) fwd, "f2 O"); 
    h4_checkptmx_DumpFBRow(cpx, 0, (float *) bck, "bck"); 
  }
  if (cpx->bck) save_debug_row_fb_sse(cpx, cpx->bck, bck, 0, cpx->bcksc); 
  if ((status = posterior_decode_row_sse(cpx, 0, sm, sm_thresh, Tvalue)) != eslOK) return status;
  cpx->bcksc += xN;
  cpx->bcksc -= mo->nullsc;
#endif

  h4_sparsemask_Finish(sm);
  return eslOK;
}
/*----------- end forward/backward API calls --------------------*/



/*****************************************************************
 * 2. Internal functions: inlined recursions
 *****************************************************************/

/* forward_row_sse()
 * Computes one row of the checkpointed Forward matrix.
 *
 * <xi>   dsq[i]; residue on this row.
 * <hmm>  query model
 * <mo>   alignment mode
 * <dpp>  ptr to previous row; may be a tmp row or a checkpointed row.
 * <dpc>  ptr to current row; ditto.
 * <Q>    number of vectors in <dpp> and <dpc>; cpx->Q.
 * 
 * Upon return, <dpc> contains scaled Forward values.
 * Returns log2(scalevalue), for the caller to accumulate
 * as it accumulates the total Forward score.
 */
static inline float
forward_row_sse(ESL_DSQ xi, const H4_PROFILE *hmm, const H4_MODE *mo, const __m128 *dpp, __m128 *dpc, int Q)
{
  const __m128 *rp   = (__m128 *) hmm->rfv[xi];
  const __m128 zerov = _mm_setzero_ps();
  const __m128 *tp   = (__m128 *) hmm->tfv[0];              // this steps thru all transitions except DD's
  const __m128 *itsc = (__m128 *) hmm->tfv[Q];              // this steps thru DD transitions
  const float  *xp   = (float *) (dpp + Q * h4C_NSCELLS);   // specials on previous row
  float        *xc   = (float *) (dpc + Q * h4C_NSCELLS);   //      ... and current row
  __m128       xEv   = zerov;
  __m128       xBv   = _mm_set1_ps(xp[h4C_B]);
  __m128 mpv, ipv, dpv;
  __m128 mcv, icv, dcv;
  __m128 sv;
  register __m128 cv;
  int    q;
  int    z;

  mpv = esl_sse_rightshiftz_float(H4C_MQ(dpp, Q-1));  // "rightshiftz" versions shift zero on, not -inf
  ipv = esl_sse_rightshiftz_float(H4C_IQ(dpp, Q-1));
  dpv = esl_sse_rightshiftz_float(H4C_DQ(dpp, Q-1));
  dcv = zerov;

  /* DP recursion for main states, all but the D->D path */
  for (q = 0; q < Q; q++)
    {
      /* Calculate M(i,q); hold it in tmp var <sv> */
      mcv =                 _mm_mul_ps(xBv, *tp);  tp++; // B->Mk    
      mcv = _mm_add_ps(mcv, _mm_mul_ps(mpv, *tp)); tp++; // Mk-1->Mk 
      mcv = _mm_add_ps(mcv, _mm_mul_ps(ipv, *tp)); tp++; // Ik-1->Mk 
      mcv = _mm_add_ps(mcv, _mm_mul_ps(dpv, *tp)); tp++; // Dk-1->Mk 
      mcv = _mm_mul_ps(mcv, *rp);                  rp++; // e_Mk(x_i)
      xEv = _mm_add_ps(mcv, xEv);                        // Mk->E    

      /* Advance on previous row, picking up M,I,D values. */
      mpv = *dpp++;
      ipv = *dpp++;
      dpv = *dpp++;

      /* Delayed store of M,D */
      H4C_MQ(dpc, q) = mcv;              
      H4C_DQ(dpc, q) = dcv;

      /* Calculate and store I(i,q) */
      icv                  =                 _mm_mul_ps(mpv, *tp);  tp++;  // M->I
      icv                  = _mm_add_ps(icv, _mm_mul_ps(ipv, *tp)); tp++;  // I->I
      icv = H4C_IQ(dpc, q) = _mm_add_ps(icv, _mm_mul_ps(dpv, *tp)); tp++;  // D->I, and store icv.

      /* Calculation of *next* D(i,q+1); delay storage, hold in dcv */
      dcv    =                 _mm_mul_ps(dcv, itsc[q]);
      dcv    = _mm_add_ps(dcv, _mm_mul_ps(mcv, *tp));     tp++;
      dcv    = _mm_add_ps(dcv, _mm_mul_ps(icv, *tp));     tp++;
    }

  /* Now the equivalent of Farrar's "lazy F" step. DD paths, alas,
   * don't parallelize well. So far we've only done DD paths within 1
   * segment. Now, alas, we need to extend them. I've tried three
   * ways:
   * 
   *   1. just go ahead and fully serialize, run V-1 loops
   *   2. stop when extending path mass drops below machine eps.
   *   3. stop when that mass is < a specified relative eps.
   *
   * Empirically, (2) seems to run fastest, on a wide range of
   * different model lengths. In the past (in H3) (1) was faster for
   * small models, so I'm leaving that code here pseudogenized for
   * future exaptation.  (3) uses an expensive division, and I wasn't
   * able to see any benefit.
   *
   * H3 kept the _mm_movemask_ps() test outside the inner loop, on
   * the theory that the conditional branch would be expensive, but 
   * upon testing, it's better to put it inside.
   */
  for (z = 1; z < 4; z++) // Vf = 4 for SSE
    {
      cv  = zerov;
      dcv = esl_sse_rightshiftz_float(dcv);   
      for (q = 0; q < Q; q++)
	{
	  sv = _mm_add_ps(dcv, H4C_DQ(dpc,q));                  // sv is tmp var for D(i,q)
	  cv = _mm_or_ps(cv, _mm_cmpgt_ps(sv, H4C_DQ(dpc, q))); // testing whether any sv element > D(i,q) element, without conditional branch
	  H4C_DQ(dpc,q) = sv;                                   // store sv as the updated D(i,q)
	  dcv = _mm_mul_ps(dcv, itsc[q]);                       // and extend DD paths one more.
	  if (! _mm_movemask_ps(cv)) break;                     // if we didn't change any D(i,q), break out.
	}
      if (! _mm_movemask_ps(cv)) break;                         // if we broke out because we didn't change D(i,q), break again.
    }
 
  // pseudogenized alternative version, fully serialized
  // for (z = 1; z < 4; z++)   
  //   {
  //     dcv  = esl_sse_rightshiftz_float(dcv);
  //     for (q = 0; q < Q; q++)
  // 	  {
  // 	    H4C_DQ(dpc,q) = _mm_add_ps(dcv, H4C_DQ(dpc, q));
  // 	    dcv           = _mm_mul_ps(dcv, itsc[q]);         
  // 	  }
  //   }

  /* Add Dk's to xEv */
  for (q = 0; q < Q; q++) xEv = _mm_add_ps(H4C_DQ(dpc,q), xEv);

  /* Specials, in order: E N JJ J B CC C */
  esl_sse_hsum_ps(xEv, &(xc[h4C_E]));  
  xc[h4C_N]  =                                     xp[h4C_N] * mo->xf[h4_N][h4_LOOP];
  xc[h4C_JJ] =                                     xp[h4C_J] * mo->xf[h4_J][h4_LOOP];
  xc[h4C_J]  = xc[h4C_JJ]                        + xc[h4C_E] * mo->xf[h4_E][h4_LOOP];
  xc[h4C_B]  = xc[h4C_N] * mo->xf[h4_N][h4_MOVE] + xc[h4C_J] * mo->xf[h4_J][h4_MOVE];
  xc[h4C_CC] =                                     xp[h4C_C] * mo->xf[h4_C][h4_LOOP];
  xc[h4C_C]  = xc[h4C_CC]                        + xc[h4C_E] * mo->xf[h4_E][h4_MOVE];

  /* Sparse rescaling. xE above threshold? Then trigger a rescaling event.             */
  if (xc[h4C_E] > 1.0e4)  /* that's a little less than e^10, ~10% of our dynamic range */
    {
      xc[h4C_N]  /= xc[h4C_E];
      xc[h4C_JJ] /= xc[h4C_E];
      xc[h4C_J]  /= xc[h4C_E];
      xc[h4C_B]  /= xc[h4C_E];
      xc[h4C_CC] /= xc[h4C_E];
      xc[h4C_C]  /= xc[h4C_E];
      xEv = _mm_set1_ps(1.0 / xc[h4C_E]);  // this is re-using xEv is a tmp var to scale the row

      for (q = 0; q < Q; q++)
        {
          H4C_MQ(dpc,q) = _mm_mul_ps(H4C_MQ(dpc,q), xEv);
          H4C_DQ(dpc,q) = _mm_mul_ps(H4C_DQ(dpc,q), xEv);
          H4C_IQ(dpc,q) = _mm_mul_ps(H4C_IQ(dpc,q), xEv);
        }

      xc[h4C_SCALE] = xc[h4C_E];
      xc[h4C_E]     = 1.0f;
    }
  else xc[h4C_SCALE] = 1.0f;

  return esl_log2f(xc[h4C_SCALE]);
}


/* backward_row_main_sse()
 * 
 * Backward calculation for rows 1..L-1. Rows 0, L are special cases.
 * 
 * <xi>        = residue dsq[i+1] on NEXT row
 * <hmm>       = query model
 * <mo>        = alignment mode
 * <dpp>       = 'previous' row i+1
 * <dpc>       = current row i
 * Q           = number of vectors in row; cpx->Q
 * scalefactor = scalefactor that Forward set for this row.
 *
 * Upon return, <dpc> contains backwards values. 
 *
 * Also, <dpp> is modified as a side effect, and is no longer
 * a valid backward row; caller should consider it destroyed. Its M
 * values now include the emission score contribution. This is
 * done for efficiency; better to add them in once and be done with
 * it.
 */
static inline void
backward_row_main_sse(ESL_DSQ xi, const H4_PROFILE *hmm, const H4_MODE *mo, __m128 *dpp, __m128 *dpc, int Q, float scalefactor)
{
  const __m128 *rp       = (__m128 *) hmm->rfv[xi];            // emission scores on row i+1, for bck; xi = dsq[i+1]  
  float       * const xc = (float *) (dpc + Q * h4C_NSCELLS);  // E N JJ J B CC C SCALE 
  const float * const xp = (float *) (dpp + Q * h4C_NSCELLS);
  const __m128  zerov    = _mm_setzero_ps();
  const __m128 *tpdd     = (__m128 *) hmm->tfv[Q];
  const __m128 *tp;
  __m128        xBv, xEv;
  __m128       *dp;
  __m128        dcv, mcv, ipv, mpv;
  __m128        tmmv, timv, tdmv;    /* copies of transition prob quads; a leftshift is needed as boundary cond */
  int           q;
 
  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. 
   * This invalidates <dpp> as a backwards row; its values are now
   * intermediates in the calculation of the current <dpc> row.
   */
  dp  = dpp;
  xBv = zerov;
  tp  = (__m128 *) hmm->tfv[0];          /* on first transition vector */
  for (q = 0; q < Q; q++)
    {
      *dp = _mm_mul_ps(*dp, *rp); rp++;
      xBv = _mm_add_ps(xBv, _mm_mul_ps(*dp, *tp)); dp += h4C_NSCELLS; tp += h4_NVT-1;
    }

  /* Specials. Dependencies dictate partial order C,CC,B < N,J,JJ < E */
  xc[h4C_C] = xc[h4C_CC] = xp[h4C_C] * mo->xf[h4_C][h4_LOOP];
  esl_sse_hsum_ps(xBv, &(xc[h4C_B]));
  xc[h4C_J] = xc[h4C_JJ] = xc[h4C_B] * mo->xf[h4_J][h4_MOVE] + xp[h4C_J] * mo->xf[h4_J][h4_LOOP];
  xc[h4C_N]              = xc[h4C_B] * mo->xf[h4_N][h4_MOVE] + xp[h4C_N] * mo->xf[h4_N][h4_LOOP];
  xc[h4C_E]              = xc[h4C_C] * mo->xf[h4_E][h4_MOVE] + xc[h4C_J] * mo->xf[h4_E][h4_LOOP];

  /* Initialize for the row calculation */
  mpv  = esl_sse_leftshiftz_float(*dpp); /* [1 5 9 13] -> [5 9 13 x], M(i+1,k+1) * e(M_k+1, x_{i+1}) */
  tmmv = esl_sse_leftshiftz_float( ((__m128 *) hmm->tfv[0])[h4_VMM]);
  timv = esl_sse_leftshiftz_float( ((__m128 *) hmm->tfv[0])[h4_VIM]);
  tdmv = esl_sse_leftshiftz_float( ((__m128 *) hmm->tfv[0])[h4_VDM]);
  xEv  = _mm_set1_ps(xc[h4C_E]);
  tp   = tpdd;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--)
    {
      ipv = H4C_IQ(dpp, q);

      // D<-{MIE} and one segment of DDs  (remaining DD passes are in _finish_sse).
      tp -= 3;  // tp skips DD ID MD, and is now on DI
      dcv =                 _mm_add_ps( xEv, _mm_mul_ps(dcv, tpdd[q]));
      dcv =                 _mm_add_ps( dcv, _mm_mul_ps(ipv, *tp));      
      dcv = H4C_DQ(dpc,q) = _mm_add_ps( dcv, _mm_mul_ps(mpv, tdmv));
      tp--;     // now on II

      // I<-{MI}; I<-D part delayed until after DD passes
      H4C_IQ(dpc,q) =  _mm_add_ps( _mm_mul_ps(ipv, *tp),  _mm_mul_ps(mpv, timv));
      tp--;     // now on MI

      // M<-{MIE}; M<-D part delayed until after DD passes
      mcv           = _mm_add_ps( xEv, _mm_mul_ps(mpv, tmmv));
      H4C_MQ(dpc,q) = _mm_add_ps( mcv, _mm_mul_ps(ipv, *tp));
      tp--;     // now on DM
      
      // pickups for next loop iteration
      mpv  = H4C_MQ(dpp, q);
      tdmv = *tp; tp--;
      timv = *tp; tp--;
      tmmv = *tp; tp--;
    }
  backward_row_finish_sse(hmm, dpc, Q, dcv);
  backward_row_rescale_sse(xc, dpc, Q, scalefactor);
}


/* backward_row_L_sse()
 * 
 * Backward calculation for row L; 
 * a special case because the matrix has no 'previous' row L+1.
 * Otherwise identical to backward_row_main().
 */
static inline void
backward_row_L_sse(const H4_PROFILE *hmm, const H4_MODE *mo, __m128 *dpc, int Q, float scalefactor)
{
  const __m128  zerov = _mm_setzero_ps();
  float        *xc    = (float *) (dpc + Q * h4C_NSCELLS);
  const __m128 *tpdd  = (__m128 *) hmm->tfv[Q];
  __m128       *dp;
  __m128       xEv, dcv;
  int          q;

  /* Backwards from T <- C,CC <- E;  all other specials unreachable, impossible on row L.
   * specials are stored in order E N JJ J B CC C.  
   */
  xc[h4C_C] = xc[h4C_CC] = mo->xf[h4_C][h4_MOVE];
  xc[h4C_B] = xc[h4C_J] = xc[h4C_JJ] = xc[h4C_N] = 0.0;
  xc[h4C_E] = xc[h4C_C] * mo->xf[h4_E][h4_MOVE];

  xEv  = _mm_set1_ps(xc[h4C_E]);
  dp   = dpc + Q*h4C_NSCELLS -1;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--) 
    {
      *dp-- = dcv = _mm_add_ps(xEv, _mm_mul_ps(dcv, tpdd[q]));  // D
      *dp--       = zerov;                                      // I
      *dp--       = xEv;                                        // M
    }
  backward_row_finish_sse(hmm, dpc, Q, dcv);
  backward_row_rescale_sse(xc, dpc, Q, scalefactor);
}


/* backward_row_finish_sse()
 * 
 * Finishing DD/MD/ID path contributions along the row.  Both L and
 * main (1..L-1) recursions share this.  First segment of DD
 * contributions has been calculated, and up to three remain; then
 * once D values are known, we can add the M<-D and I<-D
 * contributions.
 * 
 * Caller gets <dcv> from doing the first DD pass.  <dcv> contains D
 * scores for first vector, [1 5 9 13] in my standard example. We need
 * it as an initialization condition; we leftshift it to [5 9 13 *]
 * and continue propagating DD path backwards.
 * 
 * The DD path calculation is the most annoying part of 
 * the otherwise gorgeous Farrar algorithm; it has to be
 * at least partially serialized, if not fully.
 * 
 * <hmm>  = query model
 * <dpc> = current row, partially calculated
 * Q     = width of rows in # vectors; cpx->Q
 * dcv   = the first D vector [1 5 9 13] from caller's earlier calculation
 */
static inline void
backward_row_finish_sse(const H4_PROFILE *hmm, __m128 *dpc, int Q, __m128 dcv)
{
  const __m128 zerov = _mm_setzero_ps();
  const __m128 *tp;
  __m128        dpv;
  __m128       *dp;
  __m128        sv;
  __m128        cv;  // keeps track of whether any DD addition changes DQ(q) value 
  int           j,q;
  
  for (j = 1; j < 4; j++) // 4 = vector width in floats
    {
      dcv = esl_sse_leftshiftz_float(dcv);
      tp  = (__m128 *) hmm->tfv[Q] + Q - 1;   // position on last DD transition vector
      dp  = dpc + Q*h4C_NSCELLS - 1;          // position on last D cell on row
      cv  = zerov;
      for (q = Q-1; q >= 0; q--)
	{ // using cmpgt below tests if DD changed any D(q) without conditional branch (i.e. no if) 
	  dcv  = _mm_mul_ps(dcv, *tp); tp--;
	  sv   = _mm_add_ps(*dp, dcv);
	  cv   = _mm_or_ps(cv, _mm_cmpgt_ps(sv, *dp)); // if DD path changed DQ(dpc,q), cv bits know it now 
	  *dp  = sv; 
	  dp  -= h4C_NSCELLS;  // hop to q-1 D state
	  if (! _mm_movemask_ps(cv)) break; 
	}
      if (! _mm_movemask_ps(cv)) break; 
    }

#if 0  // fully serialized version. 
  for (j = 1; j < 4; j++)
    {
      dcv = esl_sse_leftshiftz_float(dcv);     /* [1 5 9 13] => [5 9 13 *]          */
      tp  = ((__m128 *) hmm->tfv[Q]) + Q - 1;  /* <*tp> now the [4 8 12 x] TDD quad */
      dp  = dpc + Q*h4C_NSCELLS - 1;           /* init to point at D(i,q) vector    */
      for (q = Q-1; q >= 0; q--)
	{
	  dcv = _mm_mul_ps(dcv, *tp); tp--;
	  *dp = _mm_add_ps(*dp, dcv); dp -= h4C_NSCELLS;
	}
    }
#endif


  /* Now we can pull the M<-D and I<-D path contributions
   * these couldn't be added to M/I until we'd finished calculating D values on row.
   */
  dpv = esl_sse_leftshiftz_float(H4C_DQ(dpc, 0));
  tp  = (__m128 *) hmm->tfv[Q] - 1; // positioned on last ID transition
  dp  = dpc + Q*h4C_NSCELLS - 1;    // positioned on last I state
  for (q = Q-1; q >= 0; q--)
    {
      dcv = *dp--;
      *dp = _mm_add_ps(*dp, _mm_mul_ps(dpv, *tp)); dp--;  tp--;           // dp moves to M state; tp to MD transition
      *dp = _mm_add_ps(*dp, _mm_mul_ps(dpv, *tp)); dp--;  tp -= h4_NVT-2; // dp moves to prev D state; tp to prev ID transition
      dpv = dcv;
    }
}


/* backward_row_rescale_sse()
 * 
 * Sparse rescaling, using the scalefactor that Forward set and 
 * Backward shares.
 * 
 * xc          - ptr to specials on row (floats)
 * dpc         - ptr to vectors for row      
 * Q           - number of vectors; cpx->Q
 * scalefactor - scalefactor that Forward set for this row
 *
 * Upon return, values in current row <dpc> have been rescaled.
 */
static inline void
backward_row_rescale_sse(float *xc, __m128 *dpc, int Q, float scalefactor)
{
  if (scalefactor > 1.0f)
    {
      __m128  sv = _mm_set1_ps(1.0 / scalefactor);
      __m128 *dp = dpc;
      int     q;

      xc[h4C_E]  /= scalefactor;
      xc[h4C_N]  /= scalefactor;
      xc[h4C_JJ] /= scalefactor;
      xc[h4C_J]  /= scalefactor;
      xc[h4C_B]  /= scalefactor;
      xc[h4C_CC] /= scalefactor;
      xc[h4C_C]  /= scalefactor;

      for (q = 0; q < Q; q++) 
        {
          *dp = _mm_mul_ps(*dp, sv); dp++; // M
          *dp = _mm_mul_ps(*dp, sv); dp++; // I
          *dp = _mm_mul_ps(*dp, sv); dp++; // D
        }
    }
  xc[h4C_SCALE] = scalefactor;
}

/* posterior_decode_row_sse()
 *
 * In production code, we don't have to save results of posterior
 * decoding; we immediately use them to calculate sparse mask on the row.
 *
 * In debugging code, we save pp's, overwriting the forward row as we
 * go; when all pp's on the row have been calculated, if there's a
 * cpx->pp matrix, we copy the pp values there.
 * 
 * Understanding how sparse rescaling interacts with posterior decoding:
 * Let F_i, B_i be the correct (unscaled) forward/backward terms for a
 * given cell on row i that we would've calculated, had we had
 * infinite numerical range. Let T be the overall Forward or Backward
 * score of the sequence comparison; we call it T because we'll use
 * the Forward score in the terminal T state.
 * 
 * Let f_i, b_i, t be the forward/backward terms we actually calculate
 * using sparse rescaling. Let s_i be the scale terms on every row
 * i=1..L. (s_0 = 1). Sparse rescaling means that:
 *     f_i = F_i / \sum_j=1^i s_j 
 *     b_i = B_j / \sum_j=i^L s_j 
 *     t   = T   / \sum_j=1^L s_j 
 * Observe that the s_i term is present in both f_i and f_b.
 * 
 * The posterior decoding we want to calculate is:
 *     (F_i * B_i) / T
 * Substituting and cancelling scale terms gives:
 *     (f_i * b_i) * (s_i / t)
 * So in the code below: caller passes us <overall_sc> (i.e. t), and
 * we set a multiplicative "scaleterm" of (s_i / t).
 * 
 * Note that because the debugging version of posterior decoding
 * code overwrites the forward row, if the caller wants to 
 * dump, copy, or test anything on the forward row, it must do it
 * BEFORE calling posterior_decode_row().
 * 
 * The threshold for determining 'significant' posterior alignment
 * probability is passed as <sm_thresh> (typically
 * <h4SPARSIFY_THRESH>, in <h4_config.h.in>). We take advantage of the
 * fact that if N/C/J emission postprobs exceed <1.0 - sm_thresh>, we
 * don't even need to look at the vectorized row; no cell can exceed
 * threshold.
 *
 * Can throw <eslEINVAL> on bad code (something awry in data structure initialization)
 *           <eslEMEM> on allocation failure in sparsemask
 */
static inline int
posterior_decode_row_sse(H4_CHECKPTMX *cpx, int rowi, H4_SPARSEMASK *sm, float sm_thresh, float overall_sc)
{
  int             Q        = cpx->Q;
  __m128        *fwd       = (__m128 *) cpx->dpf[cpx->R0 + cpx->R]; // a calculated fwd row R has been popped off
  const  __m128 *bck       = (__m128 *) cpx->dpf[rowi%2];
  float         *xf        = (float *) (fwd + Q * h4C_NSCELLS);
  const  float  *xb        = (float *) (bck + Q * h4C_NSCELLS);
  const __m128   threshv   = _mm_set1_ps(sm_thresh);
  float          scaleterm = xf[h4C_SCALE] / overall_sc;          // see comments above, on how rescaling affects posterior decoding equations 
  const __m128   cv        = _mm_set1_ps(scaleterm);
  __m128 mask;
  __m128 pv;
  float  pnonhomology;
  int    maskbits;              // xxxx 4-bit mask for which cells 0..3 have passed threshold (if any) 
  int    q,r;
  int    status;

  /* test to see if *any* cells can meet threshold, before wasting
   * time looking at them all.
   * 
   * Useful side effect: row 0 automatically fails this test (all pp
   * in S->N->B), so posterior_decode_row() can be called on row 0 (in
   * debugging code, we need to decode and store row zero specials),
   * without triggering contract check failures in p7_sparsemask_* API
   * functions that are checking for i=1..L.
   * 
   * This code block MUST NOT CHANGE the contents of fwd, bck vectors;
   * in debugging, we will use them again to recalculate and store the
   * decoding.
   */
  pnonhomology = (xf[h4C_N] * xb[h4C_N] + xf[h4C_JJ] * xb[h4C_JJ] + xf[h4C_CC] * xb[h4C_CC]) * scaleterm;
  if (pnonhomology <= 1.0f - sm_thresh)
    {
      if ((status = h4_sparsemask_StartRow(sm, rowi)) != eslOK) return status;
      for (q = Q-1; q >= 0; q--)                   // reverse, because SPARSEMASK is entirely in reversed order 
        {
          pv       =                _mm_mul_ps(H4C_MQ(fwd, q), H4C_MQ(bck, q));
          pv       = _mm_add_ps(pv, _mm_mul_ps(H4C_IQ(fwd, q), H4C_IQ(bck, q)));
          pv       = _mm_add_ps(pv, _mm_mul_ps(H4C_DQ(fwd, q), H4C_DQ(bck, q)));
          pv       = _mm_mul_ps(pv, cv);           // pv is now the posterior probability of elements q,r=0..3 
          mask     = _mm_cmpge_ps(pv, threshv);    // mask now has all 0's in elems r that failed thresh; all 1's for r that passed 
          maskbits = _mm_movemask_ps(mask);        // maskbits is now something like 0100: 1's indicate which cell passed. 
          
          for (r = 0; r < 4; r++)  // 4 = width of vectors in floats
            if ( maskbits & (1<<r)) 
              if ((status = h4_sparsemask_Add(sm, q, r)) != eslOK) return status;
        }
      if ((status = h4_sparsemask_FinishRow(sm)) != eslOK) return status;
    }

#if eslDEBUGLEVEL > 0
  xf[h4C_E]  = xf[h4C_E]  * xb[h4C_E]  * scaleterm;
  xf[h4C_N]  = (rowi == 0 ? 1.0f : xf[h4C_N]  * xb[h4C_N]  * scaleterm);
  xf[h4C_JJ] = xf[h4C_JJ] * xb[h4C_JJ] * scaleterm;
  xf[h4C_J]  = xf[h4C_J]  * xb[h4C_J]  * scaleterm;
  xf[h4C_B]  = xf[h4C_B]  * xb[h4C_B]  * scaleterm;
  xf[h4C_CC] = xf[h4C_CC] * xb[h4C_CC] * scaleterm;
  xf[h4C_C]  = xf[h4C_C]  * xb[h4C_C]  * scaleterm;

  for (q = 0; q < Q; q++)       
    {
      H4C_MQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(H4C_MQ(fwd, q), H4C_MQ(bck, q)));
      H4C_IQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(H4C_IQ(fwd, q), H4C_IQ(bck, q)));
      H4C_DQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(H4C_DQ(fwd, q), H4C_DQ(bck, q)));
    }

  if (cpx->pp)  save_debug_row_pp_sse(cpx, fwd, rowi);
#endif
  return eslOK;
}
/*------------------ end, inlined recursions -------------------*/




/*****************************************************************
 * 3. Debugging tools
 *****************************************************************/
#if eslDEBUGLEVEL > 0
/* backward_row_zero_sse()
 * 
 * Slightly peculiar but true: in production code we don't
 * need to calculate backward row 0, because we're only doing
 * posterior decoding on residues 1..L.  We only need row 0
 * if we need a complete Backward calculation -- which happens
 * when we're in debugging mode, and we're going to test the
 * Backwards score or compare the Backwards matrix to the
 * reference implementation.
 * 
 * x1    - residue dsq[1] on row 1
 * hmm - query model
 * cpx  - DP matrix; dpf[0] is the current Backward row
 * 
 * Upon return, cpx->dpf[0] (backward row 0) has been calculated.
 * Returns the log of the value in the N cell at row 0. When this 
 * is added to the sum of the logs of all the scalefactors (which
 * the caller has been accumulating in <cpx->bcksc>), then
 * the <cpx->bcksc> value is finished and equal to the Backwards
 * score in bits.
 */
static inline float
backward_row_zero_sse(ESL_DSQ x1, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx)
{
  int          Q     = cpx->Q;
  __m128       *dpc  = (__m128 *) cpx->dpf[0];
  __m128       *dpp  = (__m128 *) cpx->dpf[1];
  const __m128 *rp   = (__m128 *) hmm->rfv[x1];
  const __m128 zerov = _mm_setzero_ps();
  float        *xc   = (float *) (dpc + Q * h4C_NSCELLS); // special states on current row i      
  float        *xp   = (float *) (dpp + Q * h4C_NSCELLS); // special states on "previous" row i+1 
  __m128        xBv  = zerov;
  __m128       *dp;
  __m128       *tp;
  int           q;

  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. */
  dp  = dpp;
  tp  = (__m128 *) hmm->tfv[0];
  for (q = 0; q < Q; q++)
    {
      *dp = _mm_mul_ps(*dp, *rp); rp++;
      xBv = _mm_add_ps(xBv, _mm_mul_ps(*dp, *tp)); dp+= h4C_NSCELLS; tp += h4_NVT-1;
    }

  /* Only B,N,E will decode to nonzero posterior probability; C,J,E
   * can't be reached.  However, Backwards recursion gives C,J,E
   * nonzero values, because it's the Forward term that makes them
   * impossible. Though it's tempting to set these values to 0.0 (since
   * we know they're impossible), just do the calculation anyway;
   * by convention, all our DP algorithms (sparse, reference, and here)
   * do this, and we have unit tests comparing the values.
   */
  xc[h4C_C] = xc[h4C_CC] = xp[h4C_C] * mo->xf[h4_C][h4_LOOP];
  esl_sse_hsum_ps(xBv, &(xc[h4C_B]));
  xc[h4C_J] = xc[h4C_JJ] = xc[h4C_B] * mo->xf[h4_J][h4_MOVE] + xp[h4C_J] * mo->xf[h4_J][h4_LOOP];
  xc[h4C_N]              = xc[h4C_B] * mo->xf[h4_N][h4_MOVE] + xp[h4C_N] * mo->xf[h4_N][h4_LOOP];
  xc[h4C_E]              = xc[h4C_C] * mo->xf[h4_E][h4_MOVE] + xc[h4C_J] * mo->xf[h4_E][h4_LOOP];
  
  /* not needed in production code: */
  for (q = 0; q < Q; q++)
    H4C_MQ(dpc, q) = H4C_IQ(dpc, q) = H4C_DQ(dpc, q) = zerov;

  return log2f(xc[h4C_N]);
}


static void
save_debug_row_pp_sse(H4_CHECKPTMX *cpx, __m128 *dpc, int i)
{
  union { __m128 v; float x[4]; } u;  // 4 = number of floats per SSE vector
  int      Q  = cpx->Q;
  float  *xc  = (float *) (dpc + Q*h4C_NSCELLS);
  int     q,k,z,y;

  if (! cpx->pp) return;
  
  H4R_XMX(cpx->pp,i,h4R_E)  = xc[h4C_E];
  H4R_XMX(cpx->pp,i,h4R_N)  = xc[h4C_N];
  H4R_XMX(cpx->pp,i,h4R_J)  = xc[h4C_J];
  H4R_XMX(cpx->pp,i,h4R_B)  = xc[h4C_B];
  H4R_XMX(cpx->pp,i,h4R_L)  = xc[h4C_B]; // all mass in local path 
  H4R_XMX(cpx->pp,i,h4R_G)  = 0.0;       // ... none in glocal     
  H4R_XMX(cpx->pp,i,h4R_C)  = xc[h4C_C];
  H4R_XMX(cpx->pp,i,h4R_JJ) = xc[h4C_JJ];
  H4R_XMX(cpx->pp,i,h4R_CC) = xc[h4C_CC];
  
  /* in a posterior decoding matrix (prob-space), all k=0 cells are -inf by convention */
  for (y = 0; y < h4R_NSCELLS; y++) H4R_MX(cpx->pp,i,0,y) = -eslINFINITY;

  /* ... all mass is on local path for the filter, so all glocal cells are 0.0 */
  for (k =1; k <= cpx->M; k++)
    {
      H4R_MX(cpx->pp,i,k,h4R_MG) = 0.0f;
      H4R_MX(cpx->pp,i,k,h4R_IG) = 0.0f;
      H4R_MX(cpx->pp,i,k,h4R_DG) = 0.0f;
    }

  /* now the transfer from checkptmx decoding to local cells of refmx */
  for (q = 0; q < Q; q++)
    { // 4 = number of floats per SSE vector
      u.v = H4C_MQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= cpx->M) H4R_MX(cpx->pp,i,k,h4R_ML) = u.x[z]; }
      u.v = H4C_IQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= cpx->M) H4R_MX(cpx->pp,i,k,h4R_IL) = u.x[z]; }
      u.v = H4C_DQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= cpx->M) H4R_MX(cpx->pp,i,k,h4R_DL) = u.x[z]; }
    }
}


/* save_debug_row_fb_sse()
 * 
 * Debugging only. Transfer fwd or bck values (sparse scaled,
 * prob space) from a vectorized row, to appropriate row of <cpx->fwd>
 * or <cpx->bck> (log space, inclusive of partial sum of scalefactors);
 * <cpx->fwd> and <cpx->bck> should be identical (within numerical error
 * tolerance) to a reference implementation Forward/Backward in log
 * space.
 */
static void
save_debug_row_fb_sse(H4_CHECKPTMX *cpx, H4_REFMX *gx, __m128 *dpc, int i, float totscale)
{
  union { __m128 v; float x[4]; } u;  // 4 = number of floats per SSE vector
  int      Q  = cpx->Q;
  float  *xc  = (float *) (dpc + Q*h4C_NSCELLS);
  int     q,k,z;

  if (! gx) return;
  
  H4R_XMX(gx,i,h4R_E)  = esl_log2f(xc[h4C_E]) + totscale;
  H4R_XMX(gx,i,h4R_N)  = esl_log2f(xc[h4C_N]) + totscale;
  H4R_XMX(gx,i,h4R_J)  = esl_log2f(xc[h4C_J]) + totscale;
  H4R_XMX(gx,i,h4R_B)  = esl_log2f(xc[h4C_B]) + totscale;
  H4R_XMX(gx,i,h4R_L)  = H4R_XMX(gx,i,h4R_B);         // filter is local-mode. all mass assigned to local path 
  H4R_XMX(gx,i,h4R_G)  = -eslINFINITY;                // ... and no mass assigned to glocal path               
  H4R_XMX(gx,i,h4R_C)  = esl_log2f(xc[h4C_C]) + totscale;
  H4R_XMX(gx,i,h4R_JJ) = -eslINFINITY;                // JJ only saved in decoding, not fwd/bck 
  H4R_XMX(gx,i,h4R_CC) = -eslINFINITY;                // ... CC, ditto                          
  
  /* in H4_REFMX, all k=0 cells are initialized to -eslINFINITY;
   * set all glocal cells to -eslINFINITY too: */
  for (k =1; k <= cpx->M; k++)
    {
      H4R_MX(gx,i,k,h4R_MG) = -eslINFINITY;
      H4R_MX(gx,i,k,h4R_IG) = -eslINFINITY;
      H4R_MX(gx,i,k,h4R_DG) = -eslINFINITY;
    }

  /* now the transfer from checkptmx (scaled prob-space) to local cells of refmx (log-space): */
  for (q = 0; q < Q; q++)
    { // 4 = number of floats per SSE vector
      u.v = H4C_MQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= cpx->M) H4R_MX(gx,i,k,h4R_ML) = esl_log2f(u.x[z]) + totscale; }
      u.v = H4C_IQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= cpx->M) H4R_MX(gx,i,k,h4R_IL) = esl_log2f(u.x[z]) + totscale; }
      u.v = H4C_DQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= cpx->M) H4R_MX(gx,i,k,h4R_DL) = esl_log2f(u.x[z]) + totscale; }
    }
}
#endif // eslDEBUGLEVEL > 0




#else // ! eslENABLE_SSE4
/* provide callable functions even when we're `./configure --disable-sse` */
int
h4_fwdfilter_sse(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(hmm); ESL_UNUSED(mo); ESL_UNUSED(cpx); ESL_UNUSED(opt_sc); 
  esl_fatal("SSE support was not enabled at compile time. Can't use h4_fwdfilter_sse().");
  return eslFAIL; // NOTREACHED
}
int
h4_bckfilter_sse(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(hmm); ESL_UNUSED(mo); ESL_UNUSED(cpx); ESL_UNUSED(sm); ESL_UNUSED(sm_thresh); 
  esl_fatal("SSE support was not enabled at compile time. Can't use h4_bckfilter_sse().");
  return eslFAIL; // NOTREACHED
}
#endif // eslENABLE_SSE4 or not



