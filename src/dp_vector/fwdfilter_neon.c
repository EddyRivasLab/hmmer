/* Forwards/Backwards filters: ARM NEON vector implementation.
 * Original author: Tyler Camp (U. Texas, Austin)
 *
 * See fwdfilter.md for notes.
 *
 * Refer to fwdfilter_sse.c for more fully commented code. That's our
 * "primary" vector implementation. Other vector implementations
 * follow it with less commentary.
 *
 * This file is conditionally compiled when eslENABLE_NEON is defined.
 *
 * Contents:
 *    1. Forward and Backward filter: NEON implementations.
 *    2. Internal functions: inlined recursions.
 *    3. Internal debugging tools.
 */
#include "p7_config.h"
#ifdef p7_ENABLE_NEON

#include <arm_neon.h>

#include "easel.h"
#include "esl_neon.h"
#include "esl_vectorops.h"

#include "dp_reference/p7_refmx.h"
#include "dp_sparse/p7_sparsemx.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/fwdfilter.h"

static inline float forward_row_neon         (ESL_DSQ xi, const P7_OPROFILE *om, const esl_neon_128f_t *dpp, esl_neon_128f_t *dpc, int Q);
static inline void  backward_row_main_neon   (ESL_DSQ xi, const P7_OPROFILE *om,       esl_neon_128f_t *dpp, esl_neon_128f_t *dpc, int Q, float scalefactor);
static inline void  backward_row_L_neon      (            const P7_OPROFILE *om,                             esl_neon_128f_t *dpc, int Q, float scalefactor);
static inline void  backward_row_finish_neon (            const P7_OPROFILE *om,                             esl_neon_128f_t *dpc, int Q, esl_neon_128f_t dcv);
static inline void  backward_row_rescale_neon(float *xc, esl_neon_128f_t *dpc, int Q, float scalefactor);
static inline int   posterior_decode_row_neon(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc);

#if eslDEBUGLEVEL > 0
static inline float backward_row_zero_neon(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox);
static        void  save_debug_row_pp_neon(P7_CHECKPTMX *ox,               esl_neon_128f_t *dpc, int i);
static        void  save_debug_row_fb_neon(P7_CHECKPTMX *ox, P7_REFMX *gx, esl_neon_128f_t *dpc, int i, float totscale);
#endif


/*****************************************************************
 * 1. Forward and Backward filter: NEON implementations
 *****************************************************************/

/* Function:  p7_ForwardFilter_neon()
 * See:       fwdfilter.c::p7_ForwardFilter()
 */
int
p7_ForwardFilter_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc)
{
  int              Q     = P7_Q(om->M, 4);    /* segment length; # of MDI vectors on each row      */
  esl_neon_128f_t *dpp   = NULL;              /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  esl_neon_128f_t *dpc   = NULL;              /* dpc points at current row                         */
  float           *xc    = NULL;              /* specials E,N,JJ,J,B,CC,C,SCALE                    */
  float            totsc = 0.0f;              /* accumulates Forward score in nats                 */
  esl_neon_128f_t  place;
  place.f32x4 = vdupq_n_f32(0);         
  const esl_neon_128f_t zerov = place; 
  int     q ;                                 /* counter over vectors 0..Q-1                        */
  int     i;                                  /* counter over residues/rows 1..L                    */
  int     b;                                  /* counter down through checkpointed blocks, Rb+Rc..1 */
  int     w;                                  /* counter down through rows in a checkpointed block  */

  /* First make sure <ox> is allocated big enough.
   * (DO NOT set any ptrs into the matrix until after this potential reallocation!)
   * Then set the size of the problem in <ox> immediately, not later,
   * because debugging dumps need this information, for example.
   * (Ditto for any debugging copy of the fwd mx).
   */
  p7_checkptmx_Reinit(ox, om->M, L);
  ox->M  = om->M;       
  ox->L  = L;
  ox->V  = p7_VWIDTH_NEON / sizeof(float);
  ox->Qf = Q = P7_Q(om->M, ox->V);
  dpp    =  (esl_neon_128f_t *) ox->dpf[ox->R0-1]; 
  xc     =  (float *) (dpp + Q*p7C_NSCELLS);       

#if eslDEBUGLEVEL > 0
  ox->dump_flags |= p7_SHOW_LOG;                     /* also sets for Backward dumps, since <ox> shared */
  if (ox->do_dumping)  p7_checkptmx_DumpFBHeader(ox);
  if (ox->fwd)  { ox->fwd->M = om->M; ox->fwd->L = L; ox->fwd->type = p7R_FORWARD; }
#endif

  /* Initialization of the zero row, including specials */
  for (q = 0; q < p7C_NSCELLS*Q; q++) dpp[q] = zerov;
  xc[p7C_N]     = 1.;
  xc[p7C_B]     = om->xf[p7O_N][p7O_MOVE]; 
  xc[p7C_E]     = xc[p7C_JJ] = xc[p7C_J]  = xc[p7C_CC] = xc[p7C_C]  = 0.;                       
  xc[p7C_SCALE] = 1.;                      

#if eslDEBUGLEVEL > 0
  if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, 0, (float *) dpp, "f1 O"); 
  if (ox->fwd)        save_debug_row_fb_neon(ox, ox->fwd, dpp, 0, totsc); 
#endif

  /* Phase one: the "a" region: all rows in this region are saved */
  for (i = 1; i <= ox->La; i++)
    {
      dpc = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R]; ox->R++;    // idiomatic for "get next save/checkpoint row" 

      totsc += forward_row_neon(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;    // current row becomes prev row

#if eslDEBUGLEVEL > 0
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, (float *) dpc, "f1 O"); 
      if (ox->fwd)        save_debug_row_fb_neon(ox, ox->fwd, dpc, i, totsc); 
#endif
    }

  /* Phase two: "b" and "c" regions: partially and fully checkpointed */
  for (b = ox->Rb + ox->Rc, w = (ox->Rb ? ox->Lb : ox->Rc+1); i <= L; i++)
    {
      /* this section, deciding whether to get a checkpointed vs. tmp
       * row, is why we set <fwd>, <dpp> here, rather than having the
       * inlined forward_row_neon() set them.
       */
      if (! (--w)) {                                  // we're on the last row in segment: this row is saved
        dpc = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R]; ox->R++;  // idiomatic for "get next save/checkpoint row"
        w = b;                                        // next segment has this many rows, ending in a saved row
        b--;                                          // decrement segment number counter; last segment is r=1
      } else dpc = (esl_neon_128f_t *) ox->dpf[i%2];  // idiomatic for "next tmp row", 0/1; i%2 makes sure dpp != dpc
      
      totsc += forward_row_neon(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;

#if eslDEBUGLEVEL > 0
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, (float *) dpc, w ? "f1 X" : "f1 O"); 
      if (ox->fwd)        save_debug_row_fb_neon(ox, ox->fwd, dpc, i, totsc); 
#endif
    }

  xc     = (float *) (dpc + Q*p7C_NSCELLS);

  ESL_DASSERT1( (ox->R == ox->Ra+ox->Rb+ox->Rc) );
  ESL_DASSERT1( ( (! isnan(xc[p7C_C])) && (! isinf(xc[p7C_C]))) );

  if (opt_sc) *opt_sc = totsc + logf(xc[p7C_C] * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


/* Function:  p7_BackwardFilter()
 * See:       fwdfilter.c::p7_BackwardFilter()
 */
int
p7_BackwardFilter_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh)
{
  int              Q = ox->Q;
  esl_neon_128f_t *fwd;
  esl_neon_128f_t *bck;
  esl_neon_128f_t *dpp;
  float           *xf;
  float            Tvalue;
  int              i, b, w, i2;
  int              status;

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
  fwd = (esl_neon_128f_t *) ox->dpf[ox->R0 + ox->R]; // pop row for fwd[L] off the checkpointed stack 
  xf  = (float *) (fwd + Q*p7C_NSCELLS);
  Tvalue = xf[p7C_C] * om->xf[p7O_C][p7O_MOVE];      // i.e. scaled fwd[L] val at T state = scaled overall score
  bck = (esl_neon_128f_t *) ox->dpf[i%2];            // get tmp space for bck[L]
  backward_row_L_neon(om, bck, Q, xf[p7C_SCALE]);    // calculate bck[L] row

#if eslDEBUGLEVEL > 0
  ox->bcksc = logf(xf[p7C_SCALE]);
  if (ox->do_dumping) { 
    p7_checkptmx_DumpFBRow(ox, L, (float *) fwd, "f2 O"); 
    p7_checkptmx_DumpFBRow(ox, L, (float *) bck, "bck");  
  }
  if (ox->bck) save_debug_row_fb_neon(ox, ox->bck, bck, L, ox->bcksc); 
#endif

  if ( (status = posterior_decode_row_neon(ox, i, sm, sm_thresh, Tvalue)) != eslOK)  return status;
  i--;
  dpp = bck;

  /* If there's any checkpointing, there's an L-1 row to fill now. */
  if (ox->Rb+ox->Rc > 0)
    {
      /* Compute fwd[L-1] from last checkpoint, which we know is fwd[L-2] */
      dpp = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R-1];  /* fwd[L-2] values, already known        */
      fwd = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R];    /* get free row memory from top of stack */
      forward_row_neon(dsq[i], om, dpp, fwd, Q);          /* calculate fwd[L-1]                    */
#if eslDEBUGLEVEL > 0
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, (float *) fwd, "f2 X");
#endif

      /* Compute bck[L-1] from bck[L]. */
      xf  = (float *) (fwd + Q*p7C_NSCELLS);
      dpp = (esl_neon_128f_t *) ox->dpf[(i+1)%2]; 
      bck = (esl_neon_128f_t *) ox->dpf[i%2];
      backward_row_main_neon(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#if eslDEBUGLEVEL > 0
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, bck, "bck");
      if (ox->bck)        save_debug_row_fb_neon(ox, ox->bck, bck, i, ox->bcksc); 
#endif

      /* And decode. */
      if ( (status = posterior_decode_row_neon(ox, i, sm, sm_thresh, Tvalue)) != eslOK)  return status;
      dpp = bck;
      i--;  // i is now L-2 if there's checkpointing; else it's L-1 
    }

  /* Main loop for checkpointed regions (b,c) */
  for (b = 2; b <= ox->Rb+ox->Rc; b++)
    { // i=L-2 as we enter here, and <dpp> is on bck[L-1] 
      w = (b <= ox->Rc ? b+1 : ox->Lb);

      /* We know current row i (r=R0+R-1) ends a block and is checkpointed in fwd. */
      ox->R--;
      fwd = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R];  // pop checkpointed forward row off "stack" 
      xf  = (float *) (fwd + Q*p7C_NSCELLS);
      Tvalue = xf[p7C_C] * om->xf[p7O_C][p7O_MOVE];     // i.e. scaled fwd[L] val at T state = scaled overall score

      /* Calculate bck[i]; <dpp> is already bck[i+1] */
      bck = (esl_neon_128f_t *) ox->dpf[i%2];           // get available tmp memory for row 
      backward_row_main_neon(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);

#if eslDEBUGLEVEL > 0
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) { 
        p7_checkptmx_DumpFBRow(ox, i, (float *) fwd, "f2 O");    
        p7_checkptmx_DumpFBRow(ox, i, (float *) bck, "bck"); 
      }
      if (ox->bck)        save_debug_row_fb_neon(ox, ox->bck, bck, i, ox->bcksc); 
#endif

      /* And decode checkpointed row i. */
      if ( (status = posterior_decode_row_neon(ox, i, sm, sm_thresh, Tvalue)) != eslOK)  return status;

      /* The rest of the rows in the block weren't checkpointed.
       * Compute Forwards from last checkpoint ...
       */
      dpp = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R-1];       
      for (i2 = i-w+1; i2 <= i-1; i2++)
        {
          fwd = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R]; ox->R++;  /* push new forward row on "stack"     */
          forward_row_neon(dsq[i2], om, dpp, fwd, Q);
          dpp = fwd;      
        }

      /* ... and compute Backwards over the block we just calculated, while decoding. */
      dpp = bck;
      for (i2 = i-1; i2 >= i-w+1; i2--)
        {
          ox->R--;
          fwd = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R];
          xf  = (float *) (fwd + Q*p7C_NSCELLS);
          bck = (esl_neon_128f_t *) ox->dpf[i2%2];        
          backward_row_main_neon(dsq[i2+1], om, dpp, bck, Q, xf[p7C_SCALE]);

#if eslDEBUGLEVEL > 0
          ox->bcksc += logf(xf[p7C_SCALE]);
          if (ox->do_dumping) { 
            p7_checkptmx_DumpFBRow(ox, i2, (float *) fwd, "f2 X"); 
            p7_checkptmx_DumpFBRow(ox, i2, (float *) bck, "bck"); 
          }
          if (ox->bck)        save_debug_row_fb_neon(ox, ox->bck, bck, i2, ox->bcksc); 
#endif

          if ((status = posterior_decode_row_neon(ox, i2, sm, sm_thresh, Tvalue)) != eslOK)  return status;
          dpp = bck;
        }
      i -= w;
    }
   /* now i=La as we leave the checkpointed regions; or i=L-1 if there was no checkpointing */
   

   /* The uncheckpointed "a" region */
   for (; i >= 1; i--)
     {
       ox->R--; 
       fwd = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R]; /* pop off calculated row fwd[i]           */
       xf  = (float *) (fwd + Q*p7C_NSCELLS);
       bck = (esl_neon_128f_t *) ox->dpf[i%2];         /* get open space for bck[i]               */
       backward_row_main_neon(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);

#if eslDEBUGLEVEL > 0
       ox->bcksc += logf(xf[p7C_SCALE]);
       if (ox->do_dumping) { 
         p7_checkptmx_DumpFBRow(ox, i, (float *) fwd, "f2 O");
         p7_checkptmx_DumpFBRow(ox, i, (float *) bck, "bck"); 
       }
       if (ox->bck)        save_debug_row_fb_neon(ox, ox->bck, bck, i, ox->bcksc); 
#endif

       if ((status = posterior_decode_row_neon(ox, i, sm, sm_thresh, Tvalue)) != eslOK)  return status;
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
   fwd = (esl_neon_128f_t *) ox->dpf[ox->R0];
   bck = (esl_neon_128f_t *) ox->dpf[i%2];             
   xN = backward_row_zero_neon(dsq[1], om, ox); 
   if (ox->do_dumping) { 
     p7_checkptmx_DumpFBRow(ox, 0, (float *) fwd, "f2 O");
     p7_checkptmx_DumpFBRow(ox, 0, (float *) bck, "bck"); 
   }
   if (ox->bck)        save_debug_row_fb_neon(ox, ox->bck, bck, 0, ox->bcksc); 
   if ((status = posterior_decode_row_neon(ox, 0, sm, sm_thresh, Tvalue)) != eslOK)  return status;
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
forward_row_neon(ESL_DSQ xi, const P7_OPROFILE *om, const esl_neon_128f_t *dpp, esl_neon_128f_t *dpc, int Q)
{
  const  esl_neon_128f_t *rp   = (esl_neon_128f_t *) om->rfv[xi];
  esl_neon_128f_t         place;
  place.f32x4 = vdupq_n_f32(0);
  const    esl_neon_128f_t  zerov = place;
  const    esl_neon_128f_t *tp    = (esl_neon_128f_t *) om->tfv;
  const    float           *xp    = (float *) (dpp + Q * p7C_NSCELLS);
  float                    *xc    = (float *) (dpc + Q * p7C_NSCELLS);
  esl_neon_128f_t           dcv;
  esl_neon_128f_t           xEv;
  esl_neon_128f_t           xBv;
  esl_neon_128f_t mpv, dpv, ipv;
  esl_neon_128f_t sv;
  int    q;
  int    j;

  dcv.f32x4 = vdupq_n_f32(0); 
  xEv.f32x4 = vdupq_n_f32(0);  
  xBv.f32x4 = vdupq_n_f32(xp[p7C_B]);   

  mpv = esl_neon_rightshift_float(P7C_MQ(dpp, Q-1), zerov); 
  ipv = esl_neon_rightshift_float(P7C_IQ(dpp, Q-1), zerov); 
  dpv = esl_neon_rightshift_float(P7C_DQ(dpp, Q-1), zerov); 

  /* DP recursion for main states, all but the D->D path */
  for (q = 0; q < Q; q++)
    {
      /* Calculate M(i,q); hold it in tmp var <sv> */
      sv.f32x4     =                vmulq_f32(xBv.f32x4, (*tp).f32x4);       tp++; /* B->Mk    */
      sv.f32x4     = vaddq_f32(sv.f32x4, vmulq_f32(mpv.f32x4, (*tp).f32x4)); tp++; /* Mk-1->Mk */
      sv.f32x4     = vaddq_f32(sv.f32x4, vmulq_f32(ipv.f32x4, (*tp).f32x4)); tp++; /* Ik-1->Mk */
      sv.f32x4     = vaddq_f32(sv.f32x4, vmulq_f32(dpv.f32x4, (*tp).f32x4)); tp++; /* Dk-1->Dk */
      sv.f32x4     = vmulq_f32(sv.f32x4, (*rp).f32x4);                       rp++; /* e_Mk(x_i)*/
      xEv.f32x4    = vaddq_f32(xEv.f32x4, sv.f32x4);                               /* Mk->E    */

      /* Advance on previous row, picking up M,D,I values. */
      mpv = *dpp++;
      dpv = *dpp++;
      ipv = *dpp++;

      /* Delayed store of M,D */
      P7C_MQ(dpc, q) = sv;              
      P7C_DQ(dpc, q) = dcv;

      /* Partial calculation of *next* D(i,q+1); M->D only; delay storage, hold in dcv */
      dcv.f32x4    = vmulq_f32(sv.f32x4, (*tp).f32x4); tp++;

      /* Calculate and store I(i,q) */
      sv.f32x4             =                     vmulq_f32(mpv.f32x4, (*tp).f32x4);  tp++;
      P7C_IQ(dpc, q).f32x4 = vaddq_f32(sv.f32x4, vmulq_f32(ipv.f32x4, (*tp).f32x4)); tp++;
    }

  /* Now the DD paths. We would rather not serialize them but 
   * in an accurate Forward calculation, we have few options.
   * dcv has carried through from end of q loop above; store it 
   * in first pass, we add M->D and D->D path into DMX.
   */ 
  /* We're almost certainly're obligated to do at least one complete 
   * DD path to be sure: 
   */
  dcv            = esl_neon_rightshift_float(dcv, zerov);
  P7C_DQ(dpc, 0) = zerov;
  tp             = ((esl_neon_128f_t *) om->tfv) + 7*Q;       /* set tp to start of the DD's */
  for (q = 0; q < Q; q++) 
    {
      P7C_DQ(dpc,q).f32x4 = vaddq_f32(dcv.f32x4, P7C_DQ(dpc,q).f32x4);  
      dcv.f32x4           = vmulq_f3l2(P7C_DQ(dpc,q).f32x4, (*tp).f32x4); tp++; /* extend DMO(q), so we include M->D and D->D paths */
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
    {                   /* Fully serialized version */
      for (j = 1; j < 4; j++)
        {
          dcv = esl_neon_rightshift_float(dcv, zerov);
          tp  = ((esl_neon_128f_t *) om->tfv) + 7*Q;  /* reset tp to start of the DD's */
          for (q = 0; q < Q; q++) 
            { /* note, extend dcv, not DMO(q); only adding DD paths now */
              P7C_DQ(dpc,q).f32x4 = vaddq_f32(dcv.f32x4, P7C_DQ(dpc,q).f32x4);  
              dcv.f32x4           = vmulq_f32(dcv.f32x4, (*tp).f32x4);   tp++; 
            }       
        }
    } 
  else
    {                   /* Slightly parallelized version, but which incurs some overhead */
      for (j = 1; j < 4; j++)  // 4 = vector width, in floats
        {
          register esl_neon_128f_t cv = zerov;  /* keeps track of whether any DD's change DMO(q) */

          dcv = esl_neon_rightshift_float(dcv, zerov);
          tp  = ((esl_neon_128f_t *) om->tfv) + 7*Q;  /* set tp to start of the DD's */
          for (q = 0; q < Q; q++) 
            { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
              sv.f32x4            = vaddq_f32(dcv.f32x4, P7C_DQ(dpc,q).f32x4);
                  esl_neon_128i_t tmp;
                  tmp.u32x4 = vreinterpretq_u32_f32(cv.f32x4);  
              tmp.u32x4 = vorrq_u32(vreinterpretq_u32_f32(cv.f32x4), vcgtq_f32(sv.f32x4, P7C_DQ(dpc,q).f32x4)); 
              cv.f32x4 = vreinterpretq_f32_u32(tmp.u32x4); 
              P7C_DQ(dpc,q) = sv;                                          // store new DMO(q) 
              dcv.f32x4     = vmulq_f32(dcv.f32x4, (*tp).f32x4);   tp++;   // note, extend dcv, not DMO(q) 
            }       
          union { esl_neon_128f_t v; uint32_t i[4]; } s;
          s.v = cv;
          int mask = ((s.i[0] & 0x1) | (s.i[1] & 0x2) | (s.i[2] & 0x4) | (s.i[3] & 0x8));
          if (!mask) break; /* DD's didn't change any DMO(q)? Then done, break out. */
        }
    }
  
  /* Add Dk's to xEv */
  for (q = 0; q < Q; q++) xEv.f32x4 = vaddq_f32(P7C_DQ(dpc,q).f32x4, xEv.f32x4);

  /* Specials, in order: E N JJ J B CC C */
  esl_neon_hsum_float(xEv, &xc[p7C_E]);  
  xc[p7C_N]  =                                       xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_JJ] =                                       xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_J]  = xc[p7C_JJ]                          + xc[p7C_E] * om->xf[p7O_E][p7O_LOOP];
  xc[p7C_B]  = xc[p7C_N] * om->xf[p7O_N][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_J][p7O_MOVE];
  xc[p7C_CC] =                                       xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  xc[p7C_C]  = xc[p7C_CC]                          + xc[p7C_E] * om->xf[p7O_E][p7O_MOVE];

  /* Sparse rescaling. xE above threshold? Then trigger a rescaling event.            */
  if (xc[p7C_E] > 1.0e4)        /* that's a little less than e^10, ~10% of our dynamic range */
    {
      xc[p7C_N]  /= xc[p7C_E];
      xc[p7C_JJ] /= xc[p7C_E];
      xc[p7C_J]  /= xc[p7C_E];
      xc[p7C_B]  /= xc[p7C_E];
      xc[p7C_CC] /= xc[p7C_E];
      xc[p7C_C]  /= xc[p7C_E];
      xEv.f32x4 = vdupq_n_f32(1.0 / xc[p7C_E]);

      for (q = 0; q < Q; q++)
        {
          P7C_MQ(dpc,q).f32x4 = vmulq_f32(P7C_MQ(dpc,q).f32x4, xEv.f32x4);
          P7C_DQ(dpc,q).f32x4 = vmulq_f32(P7C_DQ(dpc,q).f32x4, xEv.f32x4);
          P7C_IQ(dpc,q).f32x4 = vmulq_f32(P7C_IQ(dpc,q).f32x4, xEv.f32x4);
        }

      xc[p7C_SCALE] = xc[p7C_E];
      xc[p7C_E]     = 1.0f;
    }
  else xc[p7C_SCALE] = 1.0f;

  return logf(xc[p7C_SCALE]);
}


static inline void
backward_row_main_neon(ESL_DSQ xi, const P7_OPROFILE *om, esl_neon_128f_t *dpp, esl_neon_128f_t *dpc, int Q, float scalefactor)
{
  const esl_neon_128f_t *rp       = (esl_neon_128f_t *) om->rfv[xi];  /* emission scores on row i+1, for bck; xi = dsq[i+1]  */
  float       * const xc = (float *) (dpc + Q * p7C_NSCELLS);         /* E N JJ J B CC C SCALE */
  const float * const xp = (float *) (dpp + Q * p7C_NSCELLS);
  const esl_neon_128f_t *tp, *tpdd;
  esl_neon_128f_t place;
  place.f32x4 = vdupq_n_f32(0.0f);
  const esl_neon_128f_t  zerov = place;
  esl_neon_128f_t        xBv;
  esl_neon_128f_t       *dp;
  esl_neon_128f_t        xEv;
  esl_neon_128f_t        dcv, mcv, ipv, mpv;
  int           q;
  esl_neon_128f_t tmmv, timv, tdmv;  
 
  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. 
   * This invalidates <dpp> as a backwards row; its values are now
   * intermediates in the calculation of the current <dpc> row.
   */
  dp  = dpp;
  xBv = zerov;
  tp  = (esl_neon_128f_t *) om->tfv;                /* on first transition vector */
  for (q = 0; q < Q; q++)
    {
      (*dp).f32x4 = vmulq_f32((*dp).f32x4, (*rp).f32x4); rp++;
      xBv.f32x4 = vaddq_f32(xBv.f32x4, vmulq_f32((*dp).f32x4, (*tp).f32x4)); dp+= p7C_NSCELLS; tp += 7;
    }

  /* Specials. Dependencies dictate partial order C,CC,B < N,J,JJ < E */
  xc[p7C_C] = xc[p7C_CC] = xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  esl_neon_hsum_float(xBv, &(xc[p7C_B]));
  xc[p7C_J] = xc[p7C_JJ] = xc[p7C_B] * om->xf[p7O_J][p7O_MOVE] + xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_N]              = xc[p7C_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_E]              = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_E][p7O_LOOP];

  /* Initialize for the row calculation */
  mpv  = esl_neon_leftshift_float(*dpp,       zerov); /* [1 5 9 13] -> [5 9 13 x], M(i+1,k+1) * e(M_k+1, x_{i+1}) */
  tmmv = esl_neon_leftshift_float((esl_neon_128f_t) om->tfv[1], zerov);
  timv = esl_neon_leftshift_float((esl_neon_128f_t) om->tfv[2], zerov);
  tdmv = esl_neon_leftshift_float((esl_neon_128f_t) om->tfv[3], zerov);
  xEv.f32x4  = vdupq_n_f32(xc[p7C_E]);
  tp   = ((esl_neon_128f_t *) om->tfv) + 7*Q - 1;
  tpdd = tp + Q;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--)
    {
      ipv                 = P7C_IQ(dpp, q);
      P7C_IQ(dpc,q).f32x4 = vaddq_f32( vmulq_f32(ipv.f32x4, (*tp).f32x4),   vmulq_f32(mpv.f32x4, timv.f32x4)); tp--; /* II,IM; I is done */
      mcv.f32x4 = vaddq_f32( vmulq_f32(ipv.f32x4, (*tp).f32x4),   vmulq_f32(mpv.f32x4, tmmv.f32x4)); tp-=2;  /* MI,MM; ME,MD remain      */
      dcv.f32x4 = vaddq_f32( vmulq_f32(dcv.f32x4, (*tpdd).f32x4), vmulq_f32(mpv.f32x4, tdmv.f32x4)); tpdd--; /* DM and one segment of DD */

      P7C_DQ(dpc,q).f32x4 = dcv.f32x4 = vaddq_f32( xEv.f32x4, dcv.f32x4);
      P7C_MQ(dpc,q).f32x4       = vaddq_f32( xEv.f32x4, mcv.f32x4);

      mpv  = P7C_MQ(dpp, q);
      tdmv = *tp; tp--;
      timv = *tp; tp--;
      tmmv = *tp; tp-=2;
    }
  backward_row_finish_neon(om, dpc, Q, dcv);
  backward_row_rescale_neon(xc, dpc, Q, scalefactor);
}


static inline void
backward_row_L_neon(const P7_OPROFILE *om,  esl_neon_128f_t *dpc, int Q, float scalefactor)
{
  esl_neon_128f_t tmp;
  tmp.f32x4 = vdupq_n_f32(0);
  const esl_neon_128f_t  zerov = tmp;
  float        *xc    = (float *) (dpc + Q * p7C_NSCELLS);
  const esl_neon_128f_t *tpdd;
  esl_neon_128f_t       *dp;
  esl_neon_128f_t       xEv, dcv;
  int                   q;

  /* Backwards from T <- C,CC <- E;  all other specials unreachable, impossible on row L.
   * specials are stored in order E N JJ J B CC C.  
   */
  xc[p7C_C] = xc[p7C_CC] = om->xf[p7O_C][p7O_MOVE];
  xc[p7C_B] = xc[p7C_J] = xc[p7C_JJ] = xc[p7C_N] = 0.0;
  xc[p7C_E] = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE];

  xEv.f32x4  = vdupq_n_f32(xc[p7C_E]);
  dp   = dpc + Q*p7C_NSCELLS - 1;
  tpdd = ((esl_neon_128f_t *) om->tfv) + 8*Q - 1;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--) 
    {
      *dp--      = zerov;                                                             /* I */
      dcv.f32x4 = vaddq_f32(xEv.f32x4, vmulq_f32(dcv.f32x4, (*tpdd).f32x4)); tpdd--;  /* D */
      *dp--     = dcv; 
      *dp--     = xEv;                                                                /* M */
    }
  backward_row_finish_neon(om, dpc, Q, dcv);
  backward_row_rescale_neon(xc, dpc, Q, scalefactor);
}


static inline void
backward_row_finish_neon(const P7_OPROFILE *om, esl_neon_128f_t *dpc, int Q, esl_neon_128f_t dcv)
{
  esl_neon_128f_t tmp;
  tmp.f32x4 = vdupq_n_f32(0);
  const esl_neon_128f_t zerov = tmp;
  const esl_neon_128f_t *tp;
  esl_neon_128f_t       *dp;
  int                    j,q;
  
  /* See notes on forward calculation: 
   * we have two options, either full serialization
   * or on long models, it becomes worthwhile to
   * check that all propagating DD contributions have
   * become negligble.
   */
  if (om->M < 100)
    { /* Full serialization */
      for (j = 1; j < 4; j++)
        {
          dcv = esl_neon_leftshift_float(dcv, zerov);     /* [1 5 9 13] => [5 9 13 *]          */
          tp  = ((esl_neon_128f_t *) om->tfv) + 8*Q - 1;  /* <*tp> now the [4 8 12 x] TDD quad */
          dp  = dpc + Q*p7C_NSCELLS - 2;                  /* init to point at D(i,q) vector    */
          for (q = Q-1; q >= 0; q--)
            {
              dcv.f32x4 = vmulq_f32(dcv.f32x4, (*tp).f32x4); tp--;
              (*dp).f32x4 = vaddq_f32((*dp).f32x4, dcv.f32x4); dp -= p7C_NSCELLS;
            }
        }
    }
  else
    { /* With check for early convergence */
      esl_neon_128f_t sv;
      esl_neon_128i_t cv;  /* keeps track of whether any DD addition changes DQ(q) value */
      for (j = 1; j < 4; j++)
        {
          dcv = esl_neon_leftshift_float(dcv, zerov);
          tp  = ((esl_neon_128f_t *) om->tfv) + 8*Q - 1;      
          dp  = dpc + Q*p7C_NSCELLS - 2;
          cv.u32x4  = vreinterpretq_u32_f32(zerov.f32x4);
          for (q = Q-1; q >= 0; q--)
            { /* using cmpgt below tests if DD changed any DMO(q) without conditional branch (i.e. no if) */
              dcv.f32x4  = vmulq_f32(dcv.f32x4, (*tp).f32x4); tp--;
              sv.f32x4   = vaddq_f32((*dp).f32x4, dcv.f32x4);
              cv.u32x4   = vorrq_u32(cv.u32x4, vcgtq_f32(sv.f32x4, (*dp).f32x4)); /* if DD path changed DQ(dpc,q), cv bits know it now */
              *dp  = sv; 
              dp  -= p7C_NSCELLS;
            }
          union { esl_neon_128i_t v; uint32_t i[4]; }s;
          s.v = cv;
          int mask = ((s.i[0] & 0x1) | (s.i[1] & 0x2) | (s.i[2] & 0x4) | (s.i[3] & 0x8));
          if (!mask) break;  /* if no DD path changed DQ(q) in this segment, then done, no more segments needed */
        }
    }

  /* Finally, M->D path contribution
   * these couldn't be added to M until we'd finished calculating D values on row.
   */
  dcv = esl_neon_leftshift_float(P7C_DQ(dpc, 0), zerov);
  tp  = ((esl_neon_128f_t *) om->tfv) + 7*Q - 3;       
  dp  = dpc + (Q-1)*p7C_NSCELLS; 
  for (q = Q-1; q >= 0; q--)
    {
      (*dp).f32x4  = vaddq_f32((*dp).f32x4, vmulq_f32(dcv.f32x4, (*tp).f32x4)); tp -= 7; 
      dcv  = *(dp+1);                                                           dp -= p7C_NSCELLS;
    }
}


static inline void
backward_row_rescale_neon(float *xc, esl_neon_128f_t *dpc, int Q, float scalefactor)
{
  if (scalefactor > 1.0f)
    {
      esl_neon_128f_t  sv;
          sv.f32x4 = vdupq_n_f32(1.0 / scalefactor);
      esl_neon_128f_t *dp = dpc;
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
          (*dp).f32x4 = vmulq_f32((*dp).f32x4, sv.f32x4); dp++; /* M */
          (*dp).f32x4 = vmulq_f32((*dp).f32x4, sv.f32x4); dp++; /* D */
          (*dp).f32x4 = vmulq_f32((*dp).f32x4, sv.f32x4); dp++; /* I */
        }
    }
  xc[p7C_SCALE] = scalefactor;
}


static inline int
posterior_decode_row_neon(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc)
{
  int                     Q         = ox->Q;
  esl_neon_128f_t        *fwd       = (esl_neon_128f_t *) ox->dpf[ox->R0 + ox->R];
  const  esl_neon_128f_t *bck       = (esl_neon_128f_t *) ox->dpf[rowi%2];
  float                  *xf        = (float *) (fwd + Q*p7C_NSCELLS);
  const  float           *xb        = (float *) (bck + Q*p7C_NSCELLS);
  esl_neon_128f_t         tmp; 
  tmp.f32x4   = vdupq_n_f32(sm_thresh);
  const esl_neon_128f_t   threshv   = tmp;
  float                   scaleterm = xf[p7C_SCALE] / overall_sc; 
  esl_neon_128f_t         place;
  place.f32x4 = vdupq_n_f32(scaleterm);  
  const esl_neon_128f_t   cv = place;
  float                   pnonhomology;
  esl_neon_128i_t         mask;
  int                     maskbits;
  esl_neon_128f_t         pv;
  int                     q,r;
  int                     status;

  /* test to see if *any* cells can meet threshold, before wasting
   * time looking at them all.
   * 
   * Useful side effect: row 0 automatically fails this test (all pp
   * in S->N->B), so posterior_decode_row_neon() can be called on row 0 (in
   * debugging code, we need to decode and store row zero specials),
   * without triggering contract check failures in p7_sparsemask_* API
   * functions that are checking for i=1..L.
   * 
   * This code block MAY NOT CHANGE the contents of fwd, bck vectors;
   * in debugging, we will use them again to recalculate and store the
   * decoding.
   */
  pnonhomology = (xf[p7C_N] * xb[p7C_N] + xf[p7C_JJ] * xb[p7C_JJ] + xf[p7C_CC] * xb[p7C_CC]) * scaleterm;
  if (pnonhomology <= 1.0f - sm_thresh)
    {
      if ((status = p7_sparsemask_StartRow(sm, rowi)) != eslOK) return status;
      for (q = Q-1; q >= 0; q--)                   // reverse, because SPARSEMASK is entirely in reversed order 
        {
          pv.f32x4       =                     vmulq_f32(P7C_MQ(fwd, q).f32x4, P7C_MQ(bck, q).f32x4);
          pv.f32x4       = vaddq_f32(pv.f32x4, vmulq_f32(P7C_IQ(fwd, q).f32x4, P7C_IQ(bck, q).f32x4));
          pv.f32x4       = vaddq_f32(pv.f32x4, vmulq_f32(P7C_DQ(fwd, q).f32x4, P7C_DQ(bck, q).f32x4));
          pv.f32x4       = vmulq_f32(pv.f32x4, cv.f32x4);         // pv is now the posterior probability of elements q,r=0..3 
          mask.u32x4     = vcgeq_f32(pv.f32x4, threshv.f32x4);    // mask now has all 0's in elems r that failed thresh; all 1's for r that passed 
          union { esl_neon_128i_t v; uint32_t i[4]; }s;
          s.v = mask;
          maskbits = ((s.i[0] & 0x1) | (s.i[1] & 0x2) | (s.i[2] & 0x4) | (s.i[3] & 0x8));
          /* maskbits is now something like 0100: 1's indicate which cell passed. */ 
          
          for (r = 0; r < 4; r++)   // 4 = number of floats per NEON vector
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
      P7C_MQ(fwd, q).f32x4 = vmulq_f32(cv.f32x4, vmulq_f32(P7C_MQ(fwd, q).f32x4, P7C_MQ(bck, q).f32x4));
      P7C_DQ(fwd, q).f32x4 = vmulq_f32(cv.f32x4, vmulq_f32(P7C_DQ(fwd, q).f32x4, P7C_DQ(bck, q).f32x4));
      P7C_IQ(fwd, q).f32x4 = vmulq_f32(cv.f32x4, vmulq_f32(P7C_IQ(fwd, q).f32x4, P7C_IQ(bck, q).f32x4));
    }

  if (ox->pp)  save_debug_row_pp_neon(ox, fwd, rowi);
#endif
  return eslOK;
}
/*------------------ end, inlined recursions -------------------*/





/*****************************************************************
 * 3. More internal functions: debugging tools
 *****************************************************************/
#if eslDEBUGLEVEL > 0

static inline float
backward_row_zero_neon(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox)
{
  int                    Q    = ox->Q;
  esl_neon_128f_t       *dpc  = (esl_neon_128f_t *) ox->dpf[0];
  esl_neon_128f_t       *dpp  = (esl_neon_128f_t *) ox->dpf[1];
  const esl_neon_128f_t *rp   = (esl_neon_128f_t *) om->rfv[x1];
  esl_neon_128f_t        place;
  place.f32x4 = vdupq_n_f32(0); 
  const esl_neon_128f_t zerov = place;
  float                 *xc   = (float *) (dpc + Q * p7C_NSCELLS); /* special states on current row i  */
  float                 *xp   = (float *) (dpp + Q * p7C_NSCELLS); /* special states on "previous" row i+1 */
  esl_neon_128f_t       *dp;
  esl_neon_128f_t       *tp;
  esl_neon_128f_t        xBv  = zerov;
  int                    q;

  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. */
  dp  = dpp;
  tp  = (esl_neon_128f_t *) om->tfv;
  for (q = 0; q < Q; q++)
    {
      (*dp).f32x4 = vmulq_f32((*dp).f32x4, (*rp).f32x4); rp++;
      xBv.f32x4 = vaddq_f32(xBv.f32x4, vmulq_f32((*dp).f32x4, (*tp).f32x4)); dp+= p7C_NSCELLS; tp += 7;
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
  esl_neon_hsum_float(xBv, &(xc[p7C_B]));
  xc[p7C_J] = xc[p7C_JJ] = xc[p7C_B] * om->xf[p7O_J][p7O_MOVE] + xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_N]              = xc[p7C_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_E]              = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_E][p7O_LOOP];
  
  /* not needed in production code: */
  for (q = 0; q < Q; q++)
    P7C_MQ(dpc, q) = P7C_IQ(dpc, q) = P7C_DQ(dpc, q) = zerov;

  return logf(xc[p7C_N]);
}


void
save_debug_row_pp_neon(P7_CHECKPTMX *ox, esl_neon_128f_t *dpc, int i)
{
  union { esl_neon_128f_t v; float x[4]; } u;  // 4 = number of floats per NEON vector
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
    { // 4 = number of floats per NEON vector
      u.v = P7C_MQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_ML) = u.x[z]; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_DL) = u.x[z]; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_IL) = u.x[z]; }
    }
}

void
save_debug_row_fb_neon(P7_CHECKPTMX *ox, P7_REFMX *gx, esl_neon_128f_t *dpc, int i, float totscale)
{
  union { esl_neon_128f_t v; float x[4]; } u;  // 4 = number of floats in NEON vector
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
    { // 4 = number of floats per NEON vector
      u.v = P7C_MQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_ML) = logf(u.x[z]) + totscale; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_DL) = logf(u.x[z]) + totscale; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_IL) = logf(u.x[z]) + totscale; }
    }
}
#endif // eslDEBUGLEVEL > 0
/*---------------- end, debugging tools -------------------------*/


                                         
#else // ! eslENABLE_NEON

/* Standard compiler-pleasing mantra for an #ifdef'd-out, empty code file. */
void p7_fwdfilter_neon_silence_hack(void) { return; }
#if defined p7FWDFILTER_NEON_TESTDRIVE || p7FWDFILTER_NEON_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_NEON or not
