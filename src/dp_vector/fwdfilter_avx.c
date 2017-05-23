/* Forwards/Backwards filters: x86 AVX vector implementation.
 *
 * This file is conditionally compiled, when eslENABLE_AVX is defined.
 *
 * Contents:
 *    1. Forward and Backward filter: AVX implementations.
 *    2. Internal functions: inlined recursions.
 *    3. Internal debugging tools.
 *
 * See fwdfilter.md for notes, and fwdfilter_sse.c, the "primary"
 * vector implementation, for more fully commented code.
 */
#include "p7_config.h"
#ifdef eslENABLE_AVX

#include <x86intrin.h>

#include "easel.h"
#include "esl_avx.h"
#include "esl_vectorops.h"

#include "dp_reference/p7_refmx.h"
#include "dp_sparse/p7_sparsemx.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/fwdfilter.h"

static inline float forward_row_avx         (ESL_DSQ xi, const P7_OPROFILE *om, const __m256 *dpp, __m256 *dpc, int Q);
static inline void  backward_row_main_avx   (ESL_DSQ xi, const P7_OPROFILE *om,       __m256 *dpp, __m256 *dpc, int Q, float scalefactor);
static inline void  backward_row_L_avx      (            const P7_OPROFILE *om,                    __m256 *dpc, int Q, float scalefactor);
static inline void  backward_row_finish_avx (            const P7_OPROFILE *om,                    __m256 *dpc, int Q, __m256 dcv);
static inline void  backward_row_rescale_avx(float *xc, __m256 *dpc, int Q, float scalefactor);
static inline int   posterior_decode_row_avx(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc);

#if eslDEBUGLEVEL > 0
static inline float backward_row_zero_avx(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox);
static        void  save_debug_row_pp_avx(P7_CHECKPTMX *ox,               __m256 *dpc, int i);
static        void  save_debug_row_fb_avx(P7_CHECKPTMX *ox, P7_REFMX *gx, __m256 *dpc, int i, float totscale);
#endif


/*****************************************************************
 * 1. Forward and Backward filter: AVX implementations
 *****************************************************************/

/* Function:  p7_ForwardFilter_avx()
 * See:       fwdfilter.c::p7_ForwardFilter()
 */
int
p7_ForwardFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc)
{
  __m256       *dpp   = NULL;                    // dpp=prev row. start on dpf[2]; rows 0,1=Backwards
  __m256       *dpc   = NULL;                    // dpc points at current row
  const __m256  zerov = _mm256_setzero_ps();     
  float        *xc    = NULL;                    // specials E,N,JJ,J,B,CC,C,SCALE
  float         totsc = 0.0f;                    // accumulates Forward score in nats
  int     Q;                                     // segment length; # of MDI vectors on each row
  int     q;                                     // counter over vectors 0..Q-1
  int     i;                                     // counter over residues/rows 1..L
  int     b;                                     // counter down through checkpointed blocks, Rb+Rc..1
  int     w;                                     // counter down through rows in a checkpointed block

  /* First make sure <ox> is allocated big enough.
   * (DO NOT set any ptrs into the matrix until after this potential reallocation!)
   * Then set the size of the problem in <ox> immediately, not later,
   * because debugging dumps need this information, for example.
   * (Ditto for any debugging copy of the fwd mx).
   */
  p7_checkptmx_Reinit(ox, om->M, L);
  ox->M  = om->M;       
  ox->L  = L;
  ox->Vf = p7_VWIDTH_AVX / sizeof(float); 
  ox->Q  = Q = P7_Q(om->M, ox->Vf);
  dpp    =  (__m256 *) ox->dpf[ox->R0-1];   
  xc     =  (float *) (dpp + Q*p7C_NSCELLS);

#if eslDEBUGLEVEL > 0
  ox->dump_flags |= p7_SHOW_LOG;                     // also sets for Backward dumps, since <ox> shared 
  if (ox->do_dumping) p7_checkptmx_DumpFBHeader(ox);
  if (ox->fwd) { ox->fwd->M = om->M; ox->fwd->L = L; ox->fwd->type = p7R_FORWARD; }
#endif

  /* Initialization of the zero row, including specials */
  for (q = 0; q < p7C_NSCELLS*Q; q++) dpp[q] = zerov;
  xc[p7C_N]     = 1.;
  xc[p7C_B]     = om->xf[p7O_N][p7O_MOVE]; 
  xc[p7C_E]     = xc[p7C_JJ] = xc[p7C_J]  = xc[p7C_CC] = xc[p7C_C]  = 0.;     
  xc[p7C_SCALE] = 1.;

#if eslDEBUGLEVEL > 0
  if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, 0, (float *) dpp, "f1 O"); 
  if (ox->fwd)        save_debug_row_fb_avx(ox, ox->fwd, dpp, 0, totsc); 
#endif

  /* Phase one: the "a" region: all rows in this region are saved */
  for (i = 1; i <= ox->La; i++)
    {
      dpc = (__m256 *) ox->dpf[ox->R0+ox->R]; ox->R++;  // idiomatic for "get next save/checkpoint row"

      totsc += forward_row_avx(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;   // current row becomes prev row 

#if eslDEBUGLEVEL > 0
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, (float *) dpc, "f1 O"); 
      if (ox->fwd)        save_debug_row_fb_avx(ox, ox->fwd, dpc, i, totsc); 
#endif
    }


  /* Phase two: "b" and "c" regions: partially and fully checkpointed */
  for (b = ox->Rb + ox->Rc, w = (ox->Rb ? ox->Lb : ox->Rc+1); i <= L; i++)
    {
      /* this section, deciding whether to get a checkpointed vs. tmp_check
       * row, is why we set <fwd>, <dpp> here, rather than having the
       * inlined forward_row() set them.
       */
      if (! (--w)) {                               // we're on the last row in segment: this row is saved 
        dpc = (__m256 *) ox->dpf[ox->R0+ox->R]; ox->R++;  // idiomatic for "get next save/checkpoint row" 
        w = b;                                     // next segment has this many rows, ending in a saved row
        b--;                                       // decrement segment number counter; last segment is r=1
      } else 
        dpc = (__m256 *) ox->dpf[i%2];             // idiomatic for "next tmp row", 0/1; i%2 makes sure dpp != dpc

      totsc += forward_row_avx(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;

#if eslDEBUGLEVEL > 0
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, (float *) dpc, w ? "f1 X" : "f1 O"); 
      if (ox->fwd)        save_debug_row_fb_avx(ox, ox->fwd, dpc, i, totsc); 
#endif
    }

  xc     = (float *) (dpc + Q*p7C_NSCELLS);

  ESL_DASSERT1( (ox->R == ox->Ra+ox->Rb+ox->Rc) );
  ESL_DASSERT1( ( (! isnan(xc[p7C_C])) && (! isinf(xc[p7C_C]))) );

  if (opt_sc) *opt_sc = totsc + logf(xc[p7C_C] * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}




/* Function:  p7_BackwardFilter_avx()
 * See:       fwdfilter.c::p7_BackwardFilter()
 */
int
p7_BackwardFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh)
{
  int     Q    = ox->Q;
  float  *xf;
  __m256 *fwd;
  __m256 *bck;
  __m256 *dpp;
  float   Tvalue;
  int     i, b, w, i2;
  int     status;

  p7_sparsemask_Reinit(sm, om->M, L);

#if eslDEBUGLEVEL > 0
  if (ox->bck) { ox->bck->M = om->M; ox->bck->L = L; ox->bck->type = p7R_BACKWARD; }
  if (ox->pp)  { ox->pp->M  = om->M; ox->pp->L  = L; ox->pp->type  = p7R_DECODING; }
#endif

  /* Row L is a special case for Backwards; so in checkpointed Backward,
   * we have to special-case the last block, rows L-1 and L
   */  
  i = L;
  ox->R--;
  fwd = (__m256 *) ox->dpf[ox->R0 + ox->R];      // pop row for fwd[L] off the checkpointed stack
  xf  = (float *) (fwd + Q*p7C_NSCELLS);
  Tvalue = xf[p7C_C] * om->xf[p7O_C][p7O_MOVE];  // i.e. scaled fwd[L] val at T state = scaled overall score
  bck = (__m256 *) ox->dpf[i%2];                 // get tmp space for bck[L]
  backward_row_L_avx(om, bck, Q, xf[p7C_SCALE]); // calculate bck[L] row

#if eslDEBUGLEVEL > 0
  ox->bcksc = logf(xf[p7C_SCALE]);
  if (ox->do_dumping) {
    p7_checkptmx_DumpFBRow(ox, L, (float *) fwd, "f2 O"); 
    p7_checkptmx_DumpFBRow(ox, L, (float *) bck, "bck");  
  }
  if (ox->bck) save_debug_row_fb_avx(ox, ox->bck, bck, L, ox->bcksc); 
#endif

  if ( (status = posterior_decode_row_avx(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
  i--;
  dpp = bck;

  /* If there's any checkpointing, there's an L-1 row to fill now. */
  if (ox->Rb+ox->Rc > 0)
    {
      /* Compute fwd[L-1] from last checkpoint, which we know is fwd[L-2] */
      dpp = (__m256 *) ox->dpf[ox->R0+ox->R-1];  // fwd[L-2] values, already known
      fwd = (__m256 *) ox->dpf[ox->R0+ox->R];    // get free row memory from top of stack
      forward_row_avx(dsq[i], om, dpp, fwd, Q);  // calculate fwd[L-1]
#if eslDEBUGLEVEL > 0
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, (float *) fwd, "f2 X");
#endif

      /* Compute bck[L-1] from bck[L]. */
      xf  = (float *) (fwd + Q*p7C_NSCELLS);
      dpp = (__m256 *) ox->dpf[(i+1)%2]; 
      bck = (__m256 *) ox->dpf[i%2];   
      backward_row_main_avx(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#if eslDEBUGLEVEL > 0
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, (float *) bck, "bck");
      if (ox->bck)        save_debug_row_fb_avx(ox, ox->bck, bck, i, ox->bcksc); 
#endif
      /* And decode. */
      if ( (status = posterior_decode_row_avx(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
      dpp = bck;
      i--;  // i is now L-2 if there's checkpointing; else it's L-1 
    }

  /* Main loop for checkpointed regions (b,c) */
  for (b = 2; b <= ox->Rb+ox->Rc; b++)
    { // i=L-2 as we enter here, and <dpp> is on bck[L-1] 
      w = (b <= ox->Rc ? b+1 : ox->Lb);

      /* We know current row i (r=R0+R-1) ends a block and is checkpointed in fwd. */
      ox->R--;
      fwd = (__m256 *) ox->dpf[ox->R0+ox->R];        // pop checkpointed forward row off "stack"
      xf  = (float *) (fwd + Q*p7C_NSCELLS);

      /* Calculate bck[i]; <dpp> is already bck[i+1] */
      bck = (__m256 *) ox->dpf[i%2];                 // get available tmp memory for row
      backward_row_main_avx(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#if eslDEBUGLEVEL > 0
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) { 
        p7_checkptmx_DumpFBRow(ox, i, (float *) fwd, "f2 O"); 
        p7_checkptmx_DumpFBRow(ox, i, (float *) bck, "bck"); 
      }
      if (ox->bck) save_debug_row_fb_avx(ox, ox->bck, bck, i, ox->bcksc); 
#endif
      /* And decode checkpointed row i. */
      if ( (status = posterior_decode_row_avx(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
      
      /* The rest of the rows in the block weren't checkpointed.
       * Compute Forwards from last checkpoint ...
       */
      dpp = (__m256 *) ox->dpf[ox->R0+ox->R-1];             // get last Fwd checkpoint.
      for (i2 = i-w+1; i2 <= i-1; i2++)
        {
          fwd = (__m256 *) ox->dpf[ox->R0+ox->R]; ox->R++;  // push new forward row on "stack"
          forward_row_avx(dsq[i2], om, dpp, fwd, Q);
          dpp = fwd;    
        }

      /* ... and compute Backwards over the block we just calculated, while decoding. */
      dpp = bck;
      for (i2 = i-1; i2 >= i-w+1; i2--)
        {
          ox->R--;
          fwd = (__m256 *) ox->dpf[ox->R0+ox->R]; // pop just-calculated forward row i2 off "stack"
          xf  = (float *) (fwd + Q*p7C_NSCELLS);
          bck = (__m256 *) ox->dpf[i2%2];         // get available for calculating bck[i2]
          backward_row_main_avx(dsq[i2+1], om, dpp, bck, Q, xf[p7C_SCALE]);

#if eslDEBUGLEVEL > 0
          ox->bcksc += logf(xf[p7C_SCALE]);
          if (ox->do_dumping) { 
            p7_checkptmx_DumpFBRow(ox, i2, (float *) fwd, "f2 X"); 
            p7_checkptmx_DumpFBRow(ox, i2, (float *) bck, "bck"); 
          }
          if (ox->bck) save_debug_row_fb_avx(ox, ox->bck, bck, i2, ox->bcksc); 
#endif
          if ((status = posterior_decode_row_avx(ox, i2, sm, sm_thresh, Tvalue)) != eslOK) return status;
          dpp = bck;
        }
      i -= w;
    }
  /* now i=La as we leave the checkpointed regions; or i=L-1 if there was no checkpointing */


  /* The uncheckpointed "a" region */
  for (; i >= 1; i--)
    {
      ox->R--; 
      fwd = (__m256 *) ox->dpf[ox->R0+ox->R]; // pop off calculated row fwd[i] 
      xf  = (float *) (fwd + Q*p7C_NSCELLS);
      bck = (__m256 *) ox->dpf[i%2];          // get open space for bck[i]
      backward_row_main_avx(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);

#if eslDEBUGLEVEL > 0
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) { 
        p7_checkptmx_DumpFBRow(ox, i, (float *) fwd, "f2 O"); 
        p7_checkptmx_DumpFBRow(ox, i, (float *) bck, "bck"); 
      }
      if (ox->bck) save_debug_row_fb_avx(ox, ox->bck, bck, i, ox->bcksc); 
#endif
 
      if ((status = posterior_decode_row_avx(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
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
  ox->R--;
  fwd = (__m256 *) ox->dpf[ox->R0];
  bck = (__m256 *) ox->dpf[i%2];              
  xN = backward_row_zero_avx(dsq[1], om, ox); 
  if (ox->do_dumping) { 
    p7_checkptmx_DumpFBRow(ox, 0, (float *) fwd, "f2 O"); 
    p7_checkptmx_DumpFBRow(ox, 0, (float *) bck, "bck"); 
  }
  if (ox->bck) save_debug_row_fb_avx(ox, ox->bck, bck, 0, ox->bcksc); 
  if ((status = posterior_decode_row_avx(ox, 0, sm, sm_thresh, Tvalue)) != eslOK) return status;
  ox->bcksc += xN;
#endif

  p7_sparsemask_Finish(sm);
  return eslOK;
}
/*----------- end forward/backward API calls --------------------*/




/*****************************************************************
 * 2. Internal functions: inlined recursions
 *****************************************************************/

static inline float
forward_row_avx(ESL_DSQ xi, const P7_OPROFILE *om, const __m256 *dpp, __m256 *dpc, int Q)
{
  const    __m256 *rp   = (__m256 *) om->rfv[xi];
  const    __m256 zerov = _mm256_setzero_ps();
  const    __m256 *tp   = (__m256 *) om->tfv;
  const    float  *xp   = (float *) (dpp + Q * p7C_NSCELLS);
  float           *xc   = (float *) (dpc + Q * p7C_NSCELLS);
  __m256          dcv   = _mm256_setzero_ps();
  __m256          xEv   = _mm256_setzero_ps();
  __m256          xBv   = _mm256_set1_ps(xp[p7C_B]);
  __m256 mpv, dpv, ipv;
  __m256 sv;
  int    q;
  int    j;

  mpv = esl_avx_leftshift_ps(P7C_MQ(dpp, Q-1)); 
  ipv = esl_avx_leftshift_ps(P7C_IQ(dpp, Q-1)); 
  dpv = esl_avx_leftshift_ps(P7C_DQ(dpp, Q-1)); 

  /* DP recursion for main states, all but the D->D path */
  for (q = 0; q < Q; q++)
    {
      /* Calculate M(i,q); hold it in tmp var <sv> */
      sv   =                   _mm256_mul_ps(xBv, *tp);  tp++; /* B->Mk    */
      sv   = _mm256_add_ps(sv, _mm256_mul_ps(mpv, *tp)); tp++; /* Mk-1->Mk */
      sv   = _mm256_add_ps(sv, _mm256_mul_ps(ipv, *tp)); tp++; /* Ik-1->Mk */
      sv   = _mm256_add_ps(sv, _mm256_mul_ps(dpv, *tp)); tp++; /* Dk-1->Dk */
      sv   = _mm256_mul_ps(sv, *rp);                     rp++; /* e_Mk(x_i)*/
      xEv  = _mm256_add_ps(xEv, sv);                           /* Mk->E    */

      /* Advance on previous row, picking up M,D,I values. */
      mpv = *dpp++;
      dpv = *dpp++;
      ipv = *dpp++;

      /* Delayed store of M,D */
      P7C_MQ(dpc, q) = sv;    
      P7C_DQ(dpc, q) = dcv;

      /* Partial calculation of *next* D(i,q+1); M->D only; delay storage, hold in dcv */
      dcv    = _mm256_mul_ps(sv, *tp); tp++;

      /* Calculate and store I(i,q) */
      sv             =                   _mm256_mul_ps(mpv, *tp);  tp++;
      P7C_IQ(dpc, q) = _mm256_add_ps(sv, _mm256_mul_ps(ipv, *tp)); tp++;
    }

  /* Now the DD paths. We would rather not serialize them but 
   * in an accurate Forward calculation, we have few options.
   * dcv has carried through from end of q loop above; store it 
   * in first pass, we add M->D and D->D path into DMX.
   */ 
  /* We're almost certainly're obligated to do at least one complete 
   * DD path to be sure: 
   */
  dcv            = esl_avx_leftshift_ps(dcv);  
  // Function naming issue: The esl_sse_(right/left)shift_ps functions 
  // describe the direction of logical shift of the vectors.  The esl_avx_leftshift functions that I wrote for
  // other filters talk about the direction of the shift instruction, which is influenced by the little-endian 
  // nature of x86.  For consistency with myself, I'm sticking to function names that match the direction of the 
  // shift instruction.
  P7C_DQ(dpc, 0) = zerov;
  tp             = ((__m256 *) om->tfv) + 7*Q; /* set tp to start of the DD's */
  for (q = 0; q < Q; q++) 
    {
      P7C_DQ(dpc,q) = _mm256_add_ps(dcv, P7C_DQ(dpc,q)); 
      dcv           = _mm256_mul_ps(P7C_DQ(dpc,q), *tp); tp++; /* extend DMO(q), so we include M->D and D->D paths */
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
    {  // Fully serialized version 
      for (j = 1; j < 8; j++)  // 8 = vector width, in floats
        { 
          dcv = esl_avx_leftshift_ps(dcv);
          tp  = ((__m256 *) om->tfv) + 7*Q;  /* reset tp to start of the DD's */
          for (q = 0; q < Q; q++) 
            { /* note, extend dcv, not DMO(q); only adding DD paths now */
              P7C_DQ(dpc,q) = _mm256_add_ps(dcv, P7C_DQ(dpc,q)); 
              dcv           = _mm256_mul_ps(dcv, *tp);   tp++; 
            }     
        }
    } 
  else
    {  // Slightly parallelized version, but which incurs some overhead 
      for (j = 1; j < 8; j++) // 8 = vector width, in floats
        {
          register __m256 cv = zerov; /* keeps track of whether any DD's change DMO(q) */

          dcv = esl_avx_leftshift_ps(dcv);
          tp  = ((__m256 *) om->tfv) + 7*Q;  /* set tp to start of the DD's */
          for (q = 0; q < Q; q++) 
            { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
              sv            = _mm256_add_ps(dcv, P7C_DQ(dpc,q)); 
              cv            = _mm256_or_ps(cv, _mm256_cmp_ps(sv, P7C_DQ(dpc,q), 14));   // 14 = code for compare greater than
              P7C_DQ(dpc,q) = sv;                                                       // store new DMO(q) 
              dcv           = _mm256_mul_ps(dcv, *tp);   tp++;                          // note, extend dcv, not DMO(q) 
            }     
          if (! _mm256_movemask_ps(cv)) break; /* DD's didn't change any DMO(q)? Then done, break out. */
        }
    }
  
  /* Add Dk's to xEv */
  for (q = 0; q < Q; q++) xEv = _mm256_add_ps(P7C_DQ(dpc,q), xEv);

  /* Specials, in order: E N JJ J B CC C */
  esl_avx_hsum_ps(xEv, &xc[p7C_E]);  
  xc[p7C_N]  =                                       xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_JJ] =                                       xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_J]  = xc[p7C_JJ]                          + xc[p7C_E] * om->xf[p7O_E][p7O_LOOP];
  xc[p7C_B]  = xc[p7C_N] * om->xf[p7O_N][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_J][p7O_MOVE];
  xc[p7C_CC] =                                       xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  xc[p7C_C]  = xc[p7C_CC]                          + xc[p7C_E] * om->xf[p7O_E][p7O_MOVE];

  /* Sparse rescaling. xE above threshold? Then trigger a rescaling event.             */
  if (xc[p7C_E] > 1.0e4)  /* that's a little less than e^10, ~10% of our dynamic range */
    {
      xc[p7C_N]  /= xc[p7C_E];
      xc[p7C_JJ] /= xc[p7C_E];
      xc[p7C_J]  /= xc[p7C_E];
      xc[p7C_B]  /= xc[p7C_E];
      xc[p7C_CC] /= xc[p7C_E];
      xc[p7C_C]  /= xc[p7C_E];
      xEv = _mm256_set1_ps(1.0 / xc[p7C_E]);

      for (q = 0; q < Q; q++)
        {
          P7C_MQ(dpc,q) = _mm256_mul_ps(P7C_MQ(dpc,q), xEv);
          P7C_DQ(dpc,q) = _mm256_mul_ps(P7C_DQ(dpc,q), xEv);
          P7C_IQ(dpc,q) = _mm256_mul_ps(P7C_IQ(dpc,q), xEv);
        }

      xc[p7C_SCALE] = xc[p7C_E];
      xc[p7C_E]     = 1.0f;
    }
  else xc[p7C_SCALE] = 1.0f;

  return logf(xc[p7C_SCALE]);
}


static inline void
backward_row_main_avx(ESL_DSQ xi, const P7_OPROFILE *om, __m256 *dpp, __m256 *dpc, int Q, float scalefactor)
{
  const __m256 *rp       = (__m256 *) om->rfv[xi];            // emission scores on row i+1, for bck; xi = dsq[i+1]
  float       * const xc = (float *) (dpc + Q * p7C_NSCELLS); // E N JJ J B CC C SCALE 
  const float * const xp = (float *) (dpp + Q * p7C_NSCELLS);
  const __m256 *tp, *tpdd;
  const __m256  zerov = _mm256_setzero_ps();
  __m256        xBv;
  __m256       *dp;
  __m256        xEv;
  __m256        dcv, mcv, ipv, mpv;
  int           q;
  __m256 tmmv, timv, tdmv;           /* copies of transition prob quads; a leftshift is needed as boundary cond */
 
  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. 
   * This invalidates <dpp> as a backwards row; its values are now
   * intermediates in the calculation of the current <dpc> row.
   */
  dp  = dpp;
  xBv = zerov;
  tp  = (__m256 *) om->tfv;    /* on first transition vector */
  for (q = 0; q < Q; q++)
    {
      *dp = _mm256_mul_ps(*dp, *rp); rp++;
      xBv = _mm256_add_ps(xBv, _mm256_mul_ps(*dp, *tp)); dp+= p7C_NSCELLS; tp += 7;
    }

  /* Specials. Dependencies dictate partial order C,CC,B < N,J,JJ < E */
  xc[p7C_C] = xc[p7C_CC] = xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  esl_avx_hsum_ps(xBv, &(xc[p7C_B]));
  xc[p7C_J] = xc[p7C_JJ] = xc[p7C_B] * om->xf[p7O_J][p7O_MOVE] + xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_N]              = xc[p7C_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_E]              = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_E][p7O_LOOP];

  /* Initialize for the row calculation */
  mpv  = esl_avx_rightshift_ps(*dpp); /* [1 5 9 13] -> [5 9 13 x], M(i+1,k+1) * e(M_k+1, x_{i+1}) */
  tmmv = esl_avx_rightshift_ps( ((__m256 *) om->tfv)[1]);
  timv = esl_avx_rightshift_ps( ((__m256 *) om->tfv)[2]);
  tdmv = esl_avx_rightshift_ps( ((__m256 *) om->tfv)[3]);
  xEv  = _mm256_set1_ps(xc[p7C_E]);
  tp   = ((__m256 *) om->tfv) + 7*Q - 1;
  tpdd = tp + Q;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--)
    {
      ipv                 = P7C_IQ(dpp, q);
      P7C_IQ(dpc,q)       = _mm256_add_ps( _mm256_mul_ps(ipv, *tp),   _mm256_mul_ps(mpv, timv)); tp--;   /* II,IM; I is done         */
      mcv                 = _mm256_add_ps( _mm256_mul_ps(ipv, *tp),   _mm256_mul_ps(mpv, tmmv)); tp-=2;  /* MI,MM; ME,MD remain      */
      dcv                 = _mm256_add_ps( _mm256_mul_ps(dcv, *tpdd), _mm256_mul_ps(mpv, tdmv)); tpdd--; /* DM and one segment of DD */

      P7C_DQ(dpc,q) = dcv = _mm256_add_ps( xEv, dcv);
      P7C_MQ(dpc,q)       = _mm256_add_ps( xEv, mcv);

      mpv  = P7C_MQ(dpp, q);
      tdmv = *tp; tp--;
      timv = *tp; tp--;
      tmmv = *tp; tp-=2;
    }
  backward_row_finish_avx(om, dpc, Q, dcv);
  backward_row_rescale_avx(xc, dpc, Q, scalefactor);
}


static inline void
backward_row_L_avx(const P7_OPROFILE *om,  __m256 *dpc, int Q, float scalefactor)
{
  const __m256  zerov = _mm256_setzero_ps();
  float        *xc    = (float *) (dpc + Q * p7C_NSCELLS);
  const __m256 *tpdd;
  __m256       *dp;
  __m256       xEv, dcv;
  int          q;

  /* Backwards from T <- C,CC <- E;  all other specials unreachable, impossible on row L.
   * specials are stored in order E N JJ J B CC C.  
   */
  xc[p7C_C] = xc[p7C_CC] = om->xf[p7O_C][p7O_MOVE];
  xc[p7C_B] = xc[p7C_J] = xc[p7C_JJ] = xc[p7C_N] = 0.0;
  xc[p7C_E] = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE];

  xEv  = _mm256_set1_ps(xc[p7C_E]);
  dp   = dpc + Q*p7C_NSCELLS - 1;
  tpdd = ((__m256 *) om->tfv) + 8*Q - 1;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--) 
    {
      *dp--       = zerov;                                                  /* I */
      *dp-- = dcv = _mm256_add_ps(xEv, _mm256_mul_ps(dcv, *tpdd)); tpdd--;  /* D */
      *dp--       = xEv;                                                    /* M */
    }
  backward_row_finish_avx(om, dpc, Q, dcv);
  backward_row_rescale_avx(xc, dpc, Q, scalefactor);
}



static inline void
backward_row_finish_avx(const P7_OPROFILE *om, __m256 *dpc, int Q, __m256 dcv)
{
  const __m256 zerov = _mm256_setzero_ps();
  const __m256 *tp;
  __m256       *dp;
  int           j,q;
  
  /* See notes on forward calculation: 
   * we have two options, either full serialization
   * or on long models, it becomes worthwhile to
   * check that all propagating DD contributions have
   * become negligble.
   */
  if (om->M < 100)
    { /* Full serialization */
      for (j = 1; j < 8; j++)
        {
          dcv = esl_avx_rightshift_ps(dcv);     /* [1 5 9 13] => [5 9 13 *]          */
          tp  = ((__m256 *) om->tfv) + 8*Q - 1; /* <*tp> now the [4 8 12 x] TDD quad */
          dp  = dpc + Q*p7C_NSCELLS - 2;        /* init to point at D(i,q) vector    */
          for (q = Q-1; q >= 0; q--)
            {
              dcv = _mm256_mul_ps(dcv, *tp); tp--;
              *dp = _mm256_add_ps(*dp, dcv); dp -= p7C_NSCELLS;
            }
        }
    }
  else
    { /* With check for early convergence */
      __m256 sv;
      __m256 cv;  /* keeps track of whether any DD addition changes DQ(q) value */
      for (j = 1; j < 8; j++) // 8 = vector width in floats
        {
          dcv = esl_avx_rightshift_ps(dcv);
          tp  = ((__m256 *) om->tfv) + 8*Q - 1;  
          dp  = dpc + Q*p7C_NSCELLS - 2;
          cv  = zerov;
          for (q = Q-1; q >= 0; q--)
            { /* using cmpgt below tests if DD changed any DMO(q) without conditional branch (i.e. no if) */
              dcv  = _mm256_mul_ps(dcv, *tp); tp--;
              sv   = _mm256_add_ps(*dp, dcv);
              cv   = _mm256_or_ps(cv, _mm256_cmp_ps(sv, *dp, 14)); /* if DD path changed DQ(dpc,q), cv bits know it now */
              *dp  = sv; 
              dp  -= p7C_NSCELLS;
            }
          if (! _mm256_movemask_ps(cv)) break; /* if no DD path changed DQ(q) in this segment, then done, no more segments needed */
        }
    }

  /* Finally, M->D path contribution
   * these couldn't be added to M until we'd finished calculating D values on row.
   */
  dcv = esl_avx_rightshift_ps(P7C_DQ(dpc, 0));
  tp  = ((__m256 *) om->tfv) + 7*Q - 3;   
  dp  = dpc + (Q-1)*p7C_NSCELLS; 
  for (q = Q-1; q >= 0; q--)
    {
      *dp  = _mm256_add_ps(*dp, _mm256_mul_ps(dcv, *tp)); tp -= 7; 
      dcv  = *(dp+1);                               dp -= p7C_NSCELLS;
    }
}



static inline void
backward_row_rescale_avx(float *xc, __m256 *dpc, int Q, float scalefactor)
{
  if (scalefactor > 1.0f)
    {
      __m256  sv = _mm256_set1_ps(1.0 / scalefactor);
      __m256 *dp = dpc;
      int     q;

      xc[p7C_E]  /= scalefactor;
      xc[p7C_N]  /= scalefactor;
      xc[p7C_JJ] /= scalefactor;
      xc[p7C_J]  /= scalefactor;
      xc[p7C_B]  /= scalefactor;
      xc[p7C_CC] /= scalefactor;
      xc[p7C_C]  /= scalefactor;

      for (q = 0; q < Q; q++) 
        {
          *dp = _mm256_mul_ps(*dp, sv); dp++; /* M */
          *dp = _mm256_mul_ps(*dp, sv); dp++; /* D */
          *dp = _mm256_mul_ps(*dp, sv); dp++; /* I */
        }
    }
  xc[p7C_SCALE] = scalefactor;
}


static inline int
posterior_decode_row_avx(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc)
{
  int             Q        = ox->Q;
  __m256        *fwd       = (__m256 *) ox->dpf[ox->R0 + ox->R]; 
  const  __m256 *bck       = (__m256 *) ox->dpf[rowi%2];
  float         *xf        = (float *) (fwd + Q*p7C_NSCELLS);
  const  float  *xb        = (float *) (bck + Q*p7C_NSCELLS);
  const __m256   threshv   = _mm256_set1_ps(sm_thresh);
  float          scaleterm = xf[p7C_SCALE] / overall_sc; 
  const __m256   cv        = _mm256_set1_ps(scaleterm);
  float  pnonhomology;
  __m256 mask;
  int    maskbits;
  __m256 pv;
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
   * This code block MAY NOT CHANGE the contents of fwd, bck vectors;
   * in debugging, we will use them again to recalculate and store the
   * decoding.
   */
  // Need to figure out how to put the sparsemask stuff back when we want to just run AVX 
  pnonhomology = (xf[p7C_N] * xb[p7C_N] + xf[p7C_JJ] * xb[p7C_JJ] + xf[p7C_CC] * xb[p7C_CC]) * scaleterm;
  if (pnonhomology <= 1.0f - sm_thresh)
    {  
      if ((status = p7_sparsemask_StartRow(sm, rowi)) != eslOK) return status;
      for (q = Q-1; q >= 0; q--)             // reverse, because SPARSEMASK is entirely in reversed order 
        {
          pv       =                   _mm256_mul_ps(P7C_MQ(fwd, q), P7C_MQ(bck, q));
          pv       = _mm256_add_ps(pv, _mm256_mul_ps(P7C_IQ(fwd, q), P7C_IQ(bck, q)));
          pv       = _mm256_add_ps(pv, _mm256_mul_ps(P7C_DQ(fwd, q), P7C_DQ(bck, q)));
          pv       = _mm256_mul_ps(pv, cv);           // pv is now the posterior probability of elements q,r=0..3 
          mask     = _mm256_cmp_ps(pv, threshv, 13);  // 13 is magic number for >= 
          // mask now has all 0's in elems r that failed thresh; all 1's for r that passed 
          maskbits = _mm256_movemask_ps(mask);    // maskbits is now something like 0100: 1's indicate which cell passed. 
  
          for (r = 0; r < 8; r++)  // 8 = width of vector in floats
            if ( maskbits & (1<<r)) 
              if ((status = p7_sparsemask_Add(sm, q, r)) != eslOK) return status;
        }
      if ((status = p7_sparsemask_FinishRow(sm)) != eslOK) return status;
    }

#if eslDEBUGLEVEL > 0
  xf[p7C_E]  = xf[p7C_E]  * xb[p7C_E]  * scaleterm;
  xf[p7C_N]  = (rowi == 0 ? 1.0f : xf[p7C_N]  * xb[p7C_N]  * scaleterm);
  xf[p7C_JJ] = xf[p7C_JJ] * xb[p7C_JJ] * scaleterm;
  xf[p7C_J]  = xf[p7C_J]  * xb[p7C_J]  * scaleterm;
  xf[p7C_B]  = xf[p7C_B]  * xb[p7C_B]  * scaleterm;
  xf[p7C_CC] = xf[p7C_CC] * xb[p7C_CC] * scaleterm;
  xf[p7C_C]  = xf[p7C_C]  * xb[p7C_C]  * scaleterm;

  for (q = 0; q < Q; q++) 
    {
      P7C_MQ(fwd, q) = _mm256_mul_ps(cv, _mm256_mul_ps(P7C_MQ(fwd, q), P7C_MQ(bck, q)));
      P7C_DQ(fwd, q) = _mm256_mul_ps(cv, _mm256_mul_ps(P7C_DQ(fwd, q), P7C_DQ(bck, q)));
      P7C_IQ(fwd, q) = _mm256_mul_ps(cv, _mm256_mul_ps(P7C_IQ(fwd, q), P7C_IQ(bck, q)));
    }

  if (ox->pp)  save_debug_row_pp_avx(ox, fwd, rowi);
#endif
  return eslOK;
}
/*------------------ end, inlined recursions -------------------*/


/*****************************************************************
 * 3. More internal functions: debugging tools
 *****************************************************************/
#if eslDEBUGLEVEL > 0

static inline float
backward_row_zero_avx(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox)
{
  int          Q     = ox->Q;
  __m256       *dpc  = (__m256 *) ox->dpf[0];
  __m256       *dpp  = (__m256 *) ox->dpf[1];
  const __m256 *rp   = (__m256 *) om->rfv[x1];
  const __m256 zerov = _mm256_setzero_ps();
  float        *xc   = (float *) (dpc + Q * p7C_NSCELLS); /* special states on current row i  */
  float        *xp   = (float *) (dpp + Q * p7C_NSCELLS); /* special states on "previous" row i+1 */
  __m256       *dp;
  __m256       *tp;
  __m256        xBv  = zerov;
  int           q;

  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. */
  dp  = dpp;
  tp  = (__m256 *) om->tfv;
  for (q = 0; q < Q; q++)
    {
      *dp = _mm256_mul_ps(*dp, *rp); rp++;
      xBv = _mm256_add_ps(xBv, _mm256_mul_ps(*dp, *tp)); dp+= p7C_NSCELLS; tp += 7;
    }

  /* Only B,N,E will decode to nonzero posterior probability; C,J,E
   * can't be reached.  However, Backwards recursion gives C,J,E
   * nonzero values, because it's the Forward term that makes them
   * impossible. Though it's tempting to set these values to 0.0 (since
   * we know they're impossible), just do the calculation anyway;
   * by convention, all our DP algorithms (sparse, reference, and here)
   * do this, and we have unit tests comparing the values.
   */
  xc[p7C_C] = xc[p7C_CC] = xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  esl_avx_hsum_ps(xBv, &(xc[p7C_B]));
  xc[p7C_J] = xc[p7C_JJ] = xc[p7C_B] * om->xf[p7O_J][p7O_MOVE] + xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_N]              = xc[p7C_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_E]              = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_E][p7O_LOOP];
  
  /* not needed in production code: */
  for (q = 0; q < Q; q++)
    P7C_MQ(dpc, q) = P7C_IQ(dpc, q) = P7C_DQ(dpc, q) = zerov;

  return logf(xc[p7C_N]);
}

static void
save_debug_row_pp_avx(P7_CHECKPTMX *ox, __m256 *dpc, int i)
{
  union { __m256 v; float x[8]; } u;  // 8 = number of floats per AVX vector
  int      Q  = ox->Q;
  float  *xc  = (float *) (dpc + Q*p7C_NSCELLS);
  int     q,k,z,s;

  if (! ox->pp) return;
  
  P7R_XMX(ox->pp,i,p7R_E)  = xc[p7C_E];
  P7R_XMX(ox->pp,i,p7R_N)  = xc[p7C_N];
  P7R_XMX(ox->pp,i,p7R_J)  = xc[p7C_J];
  P7R_XMX(ox->pp,i,p7R_B)  = xc[p7C_B];
  P7R_XMX(ox->pp,i,p7R_L)  = xc[p7C_B]; /* all mass in local path */
  P7R_XMX(ox->pp,i,p7R_G)  = 0.0;       /* ... none in glocal     */
  P7R_XMX(ox->pp,i,p7R_C)  = xc[p7C_C];
  P7R_XMX(ox->pp,i,p7R_JJ) = xc[p7C_JJ];
  P7R_XMX(ox->pp,i,p7R_CC) = xc[p7C_CC];
  
  /* in a posterior decoding matrix (prob-space), all k=0 cells are 0.0 */
  for (s = 0; s < p7R_NSCELLS; s++) P7R_MX(ox->pp,i,0,s) = 0.0f;

  /* ... all mass is on local path for the filter, so all glocal cells are 0.0 */
  for (k =1; k <= ox->M; k++)
    {
      P7R_MX(ox->pp,i,k,p7R_MG) = 0.0f;
      P7R_MX(ox->pp,i,k,p7R_IG) = 0.0f;
      P7R_MX(ox->pp,i,k,p7R_DG) = 0.0f;
    }

  /* now the transfer from checkptmx decoding to local cells of refmx */
  for (q = 0; q < Q; q++)
    { // 8 = number of floats per AVX vector
      u.v = P7C_MQ(dpc, q); for (z = 0; z < 8; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_ML) = u.x[z]; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < 8; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_DL) = u.x[z]; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < 8; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_IL) = u.x[z]; }
    }
}



static void
save_debug_row_fb_avx(P7_CHECKPTMX *ox, P7_REFMX *gx, __m256 *dpc, int i, float totscale)
{
  union { __m256 v; float x[8]; } u;  // 8 = number of floats per AVX vector
  int      Q  = ox->Q;
  float  *xc  = (float *) (dpc + Q*p7C_NSCELLS);
  int     q,k,z;

  if (! gx) return;
  
  P7R_XMX(gx,i,p7R_E)  = logf(xc[p7C_E]) + totscale;
  P7R_XMX(gx,i,p7R_N)  = logf(xc[p7C_N]) + totscale;
  P7R_XMX(gx,i,p7R_J)  = logf(xc[p7C_J]) + totscale;
  P7R_XMX(gx,i,p7R_B)  = logf(xc[p7C_B]) + totscale;
  P7R_XMX(gx,i,p7R_L)  = P7R_XMX(gx,i,p7R_B);         /* filter is local-mode. all mass assigned to local path */
  P7R_XMX(gx,i,p7R_G)  = -eslINFINITY;                /* ... and no mass assigned to glocal path               */
  P7R_XMX(gx,i,p7R_C)  = logf(xc[p7C_C]) + totscale;
  P7R_XMX(gx,i,p7R_JJ) = -eslINFINITY;                /* JJ only saved in decoding, not fwd/bck */
  P7R_XMX(gx,i,p7R_CC) = -eslINFINITY;                /* ... CC, ditto                          */
  
  /* in P7_REFMX, all k=0 cells are initialized to -eslINFINITY;
   * set all glocal cells to -eslINFINITY too: */
  for (k =1; k <= ox->M; k++)
    {
      P7R_MX(gx,i,k,p7R_MG) = -eslINFINITY;
      P7R_MX(gx,i,k,p7R_IG) = -eslINFINITY;
      P7R_MX(gx,i,k,p7R_DG) = -eslINFINITY;
    }

  /* now the transfer from checkptmx (scaled prob-space) to local cells of refmx (log-space): */
  for (q = 0; q < Q; q++)
    { // 8 = number of floats per AVX vector
      u.v = P7C_MQ(dpc, q); for (z = 0; z < 8; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_ML) = logf(u.x[z]) + totscale; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < 8; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_DL) = logf(u.x[z]) + totscale; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < 8; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_IL) = logf(u.x[z]) + totscale; }
    }
}

#endif // eslDEBUGLEVEL > 0
/*---------------------------- end, debugging tools ---------------------------*/

#else // ! eslENABLE_AVX

/* Standard compiler-pleasing mantra for an #ifdef'd-out, empty code file. */
void p7_fwdfilter_avx_silence_hack(void) { return; }
#if defined p7FWDFILTER_AVX_TESTDRIVE || p7FWDFILTER_AVX_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_AVX or not
