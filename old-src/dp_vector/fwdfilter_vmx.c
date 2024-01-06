/* Forwards/Backwards filters: x86 SSE vector implementation.
 * 
 * See fwdfilter.md for notes.
 *
 * This file is conditionally compiled when eslENABLE_VMX is defined.
 *
 * Contents:
 *    1. Forward and Backward filter: SSE implementations.
 *    2. Internal functions: inlined recursions.
 *    3. Internal debugging tools.
 */
#include <p7_config.h>

#include <float.h>
#include "easel.h"
#include "esl_vectorops.h"

#include "dp_reference/p7_refmx.h"
#include "dp_sparse/p7_sparsemx.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/fwdfilter.h"

#ifdef eslENABLE_VMX
#include <altivec.h>
#include "esl_vmx.h"

/* Vectorized DP recursions are tediously lengthy. For some semblance
 * of clarity, we break them out into one-page-ish chunks, using
 * static inlined functions.
 */
static inline float forward_row_vmx        (ESL_DSQ xi, const P7_OPROFILE *om, const vector float *dpp, vector float *dpc, int Q);
static inline void  backward_row_main_vmx  (ESL_DSQ xi, const P7_OPROFILE *om,       vector float *dpp, vector float *dpc, int Q, float scalefactor);
static inline void  backward_row_L_vmx     (            const P7_OPROFILE *om,                    vector float *dpc, int Q, float scalefactor);
static inline void  backward_row_finish_vmx(            const P7_OPROFILE *om,                    vector float *dpc, int Q, vector float dcv);
static inline void  backward_row_rescale_vmx(float *xc, vector float *dpc, int Q, float scalefactor);
static inline int   posterior_decode_row_vmx(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc);

#if eslDEBUGLEVEL > 0
static inline float backward_row_zero_vmx(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox);
static        void  save_debug_row_pp_vmx(P7_CHECKPTMX *ox,               vector float *dpc, int i);
static        void  save_debug_row_fb_vmx(P7_CHECKPTMX *ox, P7_REFMX *gx, vector float *dpc, int i, float totscale);
#endif

/*****************************************************************
 * 1. Forward and Backward filter: VMX implementations.
 *****************************************************************/

/* Function:  p7_ForwardFilter_vmx()
 * See:       fwdfilter.c::p7_ForwardFilter()
 */
int
p7_ForwardFilter_vmx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc)
{
  vector float       *dpp   = NULL;                 // dpp=prev row. start on dpf[2]; rows 0,1=Backwards
  vector float       *dpc   = NULL;                 // dpc points at current row
  const vector float  zerov = {0.0, 0.0, 0.0, 0.0};     
  float        *xc    = NULL;                 // specials E,N,JJ,J,B,CC,C,SCALE    
  float         totsc = 0.0f;                 // accumulates Forward score in nats 
  int     Q;                                  // segment length; # of MDI vectors on each row
  int     q;                                  // counter over vectors 0..Q-1
  int     i;                                  // counter over residues/rows 1..L 
  int     b;                                  // counter down through checkpointed blocks, Rb+Rc..1
  int     w;                                  // counter down through rows in a checkpointed block

  /* First make sure <ox> is allocated big enough.
   * (DO NOT set any ptrs into the matrix until after this potential reallocation!)
   * Then set the size of the problem in <ox> immediately, not later,
   * because debugging dumps need this information, for example.
   * (Ditto for any debugging copy of the fwd mx).
   */
  p7_checkptmx_Reinit(ox, om->M, L);
  ox->M   = om->M;       
  ox->L   = L;
  ox->Vf  = p7_VWIDTH_VMX / sizeof(float);
  ox->Q   = Q = P7_Q(om->M, ox->Vf);             // Q, just for shorthand avoidance of ox->Q everywhere below.
  dpp     = (vector float *) ox->dpf[ox->R0-1];        // dpp=prev row. start on dpf[2]; rows 0,1=Backwards
  xc      = (float *)  (dpp + Q*p7C_NSCELLS);    // specials E,N,JJ,J,B,CC,C,SCALE 
 
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
  if (ox->fwd)        save_debug_row_fb_vmx(ox, ox->fwd, dpp, 0, totsc); 
#endif

  /* Phase one: the "a" region: all rows in this region are saved */
  for (i = 1; i <= ox->La; i++)
    {   
      dpc = (vector float *) ox->dpf[ox->R0+ox->R]; ox->R++;    // idiomatic for "get next save/checkpoint row" 

      totsc += forward_row_vmx(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;         // current row becomes prev row

#if eslDEBUGLEVEL > 0
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, (float *) dpc, "f1 O"); 
      if (ox->fwd)        save_debug_row_fb_vmx(ox, ox->fwd, dpc, i, totsc); 
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
        dpc = (vector float *) ox->dpf[ox->R0+ox->R]; ox->R++;  // idiomatic for "get next save/checkpoint row" 
        w = b;                                     // next segment has this many rows, ending in a saved row
        b--;                                       // decrement segment number counter; last segment is r=1
      } else      
        dpc = (vector float *) ox->dpf[i%2];            // idiomatic for "next tmp row", 0/1; i%2 makes sure dpp != dpc 

      totsc += forward_row_vmx(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;

#if eslDEBUGLEVEL > 0
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, (float *) dpc, w ? "f1 X" : "f1 O"); 
      if (ox->fwd)        save_debug_row_fb_vmx(ox, ox->fwd, dpc, i, totsc); 
#endif
    }

  xc  = (float *) (dpc + Q * p7C_NSCELLS);

  ESL_DASSERT1( (ox->R == ox->Ra+ox->Rb+ox->Rc) );
  ESL_DASSERT1( ( (! isnan(xc[p7C_C])) && (! isinf(xc[p7C_C]))) );

  if (opt_sc) *opt_sc = totsc + logf(xc[p7C_C] * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}



/* Function:  p7_BackwardFilter_vmx()
 * See:       fwdfilter.c::p7_BackwardFilter()
 */
int
p7_BackwardFilter_vmx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh)
{
  int     Q    = ox->Q;     // short/clarity
  float  *xf;
  vector float *fwd;
  vector float *bck;
  vector float *dpp;
  float   Tvalue;
  int     i, b, w, i2;
  int     status;

  p7_sparsemask_Reinit(sm, om->M, L);

  /* Contract checks */
  //ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  //ESL_DASSERT1(( om->L    == L ));        /*  ... and it's easy to forget to set <om> that way           */
  //ESL_DASSERT1(( om->nj   == 1.0f ));     /*  ... hence the check                                        */
  /*  ... which you can disable, if you're playing w/ config     */
  // SRE: disabled 21 Jul 14 because unit tests in sparse_asc_fwdback need to create
  // sparse mask w/ checkpointed F/B w/ unilocal L=0 models

#if eslDEBUGLEVEL > 0
  /* Debugging instrumentations. */
  if (ox->bck) { ox->bck->M = om->M; ox->bck->L = L; ox->bck->type = p7R_BACKWARD; }
  if (ox->pp)  { ox->pp->M  = om->M; ox->pp->L  = L; ox->pp->type  = p7R_DECODING; }
#endif

  /* Row L is a special case for Backwards; so in checkpointed Backward,
   * we have to special-case the last block, rows L-1 and L
   */  
  i = L;
  ox->R--;
  fwd = (vector float *) ox->dpf[ox->R0 + ox->R];      // pop row for fwd[L] off the checkpointed stack
  xf  = (float *) (fwd + Q*p7C_NSCELLS);
  Tvalue = xf[p7C_C] * om->xf[p7O_C][p7O_MOVE];  // i.e. scaled fwd[L] val at T state = scaled overall score
  bck = (vector float *) ox->dpf[i%2];                 // get tmp space for bck[L]
  backward_row_L_vmx(om, bck, Q, xf[p7C_SCALE]); // calculate bck[L] row

#if eslDEBUGLEVEL > 0
  ox->bcksc = logf(xf[p7C_SCALE]);
  if (ox->do_dumping) { 
    p7_checkptmx_DumpFBRow(ox, L, (float *) fwd, "f2 O"); 
    if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, L, (float *) bck, "bck");  
  }
  if (ox->bck) save_debug_row_fb_vmx(ox, ox->bck, bck, L, ox->bcksc); 
#endif

  if ( (status = posterior_decode_row_vmx(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
  i--;
  dpp = bck;

  /* If there's any checkpointing, there's an L-1 row to fill now. */
  if (ox->Rb+ox->Rc > 0)
    {
      /* Compute fwd[L-1] from last checkpoint, which we know is fwd[L-2] */
      dpp = (vector float *) ox->dpf[ox->R0+ox->R-1];  // fwd[L-2] values, already known
      fwd = (vector float *) ox->dpf[ox->R0+ox->R];    // get free row memory from top of stack
      forward_row_vmx(dsq[i], om, dpp, fwd, Q);  // calculate fwd[L-1]
#if eslDEBUGLEVEL > 0
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, (float *) fwd, "f2 X");
#endif

      /* Compute bck[L-1] from bck[L]. */
      xf  = (float *) (fwd + Q*p7C_NSCELLS);
      dpp = (vector float *) ox->dpf[(i+1)%2]; 
      bck = (vector float *) ox->dpf[i%2];             // get space for bck[L-1]
      backward_row_main_vmx(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#if eslDEBUGLEVEL > 0
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, (float *) bck, "bck");
      if (ox->bck)        save_debug_row_fb_vmx(ox, ox->bck, bck, i, ox->bcksc); 
#endif
      /* And decode. */
      if ( (status = posterior_decode_row_vmx(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
      dpp = bck;
      i--;  // i is now L-2 if there's checkpointing; else it's L-1
    }

  /* Main loop for checkpointed regions (b,c) */
  for (b = 2; b <= ox->Rb+ox->Rc; b++)
    { // i=L-2 as we enter here, and <dpp> is on bck[L-1] 
      w = (b <= ox->Rc ? b+1 : ox->Lb);

      /* We know current row i (r=R0+R-1) ends a block and is checkpointed in fwd. */
      ox->R--;
      fwd = (vector float *) ox->dpf[ox->R0+ox->R];        // pop checkpointed forward row off "stack" 
      xf  = (float *)  (fwd + Q * p7C_NSCELLS);

      /* Calculate bck[i]; <dpp> is already bck[i+1] */
      bck = (vector float *) ox->dpf[i%2];        // get available tmp memory for row
      backward_row_main_vmx(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#if eslDEBUGLEVEL > 0
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) { 
        p7_checkptmx_DumpFBRow(ox, i, (float *) fwd, "f2 O");
        p7_checkptmx_DumpFBRow(ox, i, (float *) bck, "bck"); 
      }
      if (ox->bck) save_debug_row_fb_vmx(ox, ox->bck, bck, i, ox->bcksc); 
#endif
      /* And decode checkpointed row i. */
      if ( (status = posterior_decode_row_vmx(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
      
      /* The rest of the rows in the block weren't checkpointed.
       * Compute Forwards from last checkpoint ...
       */
      dpp = (vector float *) ox->dpf[ox->R0+ox->R-1];  // get last Fwd checkpoint. 
      for (i2 = i-w+1; i2 <= i-1; i2++)
        {
          fwd = (vector float *) ox->dpf[ox->R0+ox->R]; ox->R++;  // push new forward row on "stack"
          forward_row_vmx(dsq[i2], om, dpp, fwd, Q);
          dpp = fwd;      
        }

      /* ... and compute Backwards over the block we just calculated, while decoding. */
      dpp = bck;
      for (i2 = i-1; i2 >= i-w+1; i2--)
        {
          ox->R--;
          fwd = (vector float *) ox->dpf[ox->R0+ox->R]; // pop just-calculated forward row i2 off "stack"
          xf  = (float *)  (fwd + Q * p7C_NSCELLS);
          bck = (vector float *) ox->dpf[i2%2];         // get available for calculating bck[i2] 
          backward_row_main_vmx(dsq[i2+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#if eslDEBUGLEVEL > 0
          ox->bcksc += logf(xf[p7C_SCALE]);
          if (ox->do_dumping) { 
            p7_checkptmx_DumpFBRow(ox, i2, (float *) fwd, "f2 X"); 
            p7_checkptmx_DumpFBRow(ox, i2, (float *) bck, "bck"); 
          }
          if (ox->bck) save_debug_row_fb_vmx(ox, ox->bck, bck, i2, ox->bcksc); 
#endif

          if ((status = posterior_decode_row_vmx(ox, i2, sm, sm_thresh, Tvalue)) != eslOK) return status;
          dpp = bck;
        }
      i -= w;
    }
  /* now i=La as we leave the checkpointed regions; or i=L-1 if there was no checkpointing */

  /* the uncheckpointed "a" region */
  for (; i >= 1; i--)
    {
      ox->R--; 
      fwd = (vector float *) ox->dpf[ox->R0+ox->R]; // pop off calculated row fwd[i]
      xf  = (float *)  (fwd + Q * p7C_NSCELLS);
      bck = (vector float *) ox->dpf[i%2];          // get open space for bck[i]
      backward_row_main_vmx(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);

#if eslDEBUGLEVEL > 0
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) { 
        p7_checkptmx_DumpFBRow(ox, i, (float *) fwd, "f2 O"); 
        p7_checkptmx_DumpFBRow(ox, i, (float *) bck, "bck"); 
      }
      if (ox->bck) save_debug_row_fb_vmx(ox, ox->bck, bck, i, ox->bcksc); 
#endif
      if ((status = posterior_decode_row_vmx(ox, i, sm, sm_thresh, Tvalue)) != eslOK) return status;
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
  fwd = (vector float *) ox->dpf[ox->R0];
  bck = (vector float *) ox->dpf[i%2];              
  xN = backward_row_zero_vmx(dsq[1], om, ox); 
  if (ox->do_dumping) { 
    p7_checkptmx_DumpFBRow(ox, 0, (float *) fwd, "f2 O"); 
    p7_checkptmx_DumpFBRow(ox, 0, (float *) bck, "bck"); 
  }
  if (ox->bck) save_debug_row_fb_vmx(ox, ox->bck, bck, 0, ox->bcksc); 
  if ((status = posterior_decode_row_vmx(ox, 0, sm, sm_thresh, Tvalue)) != eslOK) return status;
  ox->bcksc += xN;
#endif

  p7_sparsemask_Finish(sm);
  return eslOK;
}
/*----------- end forward/backward API calls --------------------*/




/*****************************************************************
 * 2. Internal functions: inlined recursions
 *****************************************************************/

/* forward_row_vmx()
 * 
 * <xi>   dsq[i]; residue on this row.
 * <om>   query model
 * <dpp>  ptr to previous row; may be a tmp row or a checkpointed row.
 * <dpc>  ptr to current row; ditto.
 * Q      number of vectors in <dpp>, <dpc>; ox->Q.
 * 
 * Upon return, <dpc> contains scaled Forward values.
 * Returns log(scalevalue), for the caller to accumulate
 * as it accumulates the total Forward score.
 */
static inline float
forward_row_vmx(ESL_DSQ xi, const P7_OPROFILE *om, const vector float *dpp, vector float *dpc, int Q)
{
  const    vector float *rp   = (vector float *) om->rfv[xi];
  const    vector float zerov = {0.0,  0.0, 0.0, 0.0};
  const    vector float *tp   = (vector float *) om->tfv;
  const    float  *xp   = (float *) (dpp + Q * p7C_NSCELLS);
  float           *xc   = (float *) (dpc + Q * p7C_NSCELLS);
  vector float          dcv   = zerov;
  vector float          xEv   = zerov;
  vector float          xBv   = vec_splats(xp[p7C_B]);
  vector float mpv, dpv, ipv;
  vector float sv;
  int    q;
  int    j;

  mpv = esl_vmx_rightshiftz_float(P7C_MQ(dpp, Q-1));
  ipv = esl_vmx_rightshiftz_float(P7C_IQ(dpp, Q-1));
  dpv = esl_vmx_rightshiftz_float(P7C_DQ(dpp, Q-1));

  /* DP recursion for main states, all but the D->D path */
  for (q = 0; q < Q; q++)
    {
      /* Calculate M(i,q); hold it in tmp var <sv> */
      // Altivec only provides multiply-add, not just multiply.  Pass vector of zeroes to third argument of madd to
      // get multiply behavior
      sv     =                vec_madd(xBv, *tp, zerov);  tp++; /* B->Mk    */
      sv     = vec_madd(mpv, *tp, sv); tp++; /* Mk-1->Mk */
      sv     = vec_madd(ipv, *tp, sv); tp++; /* Ik-1->Mk */
      sv     = vec_madd(dpv, *tp, sv); tp++; /* Dk-1->Dk */
      sv     = vec_madd(sv, *rp, zerov);                  rp++; /* e_Mk(x_i)*/
      xEv    = vec_add(xEv, sv);                        /* Mk->E    */

      /* Advance on previous row, picking up M,D,I values. */
      mpv = *dpp++;
      dpv = *dpp++;
      ipv = *dpp++;

      /* Delayed store of M,D */
      P7C_MQ(dpc, q) = sv;              
      P7C_DQ(dpc, q) = dcv;

      /* Partial calculation of *next* D(i,q+1); M->D only; delay storage, hold in dcv */
      dcv    = vec_madd(sv, *tp, zerov); tp++;

      /* Calculate and store I(i,q) */
      sv             =                vec_madd(mpv, *tp, zerov);  tp++;
      P7C_IQ(dpc, q) = vec_madd(ipv, *tp, sv); tp++;
    }

  /* Now the DD paths. We would rather not serialize them but 
   * in an accurate Forward calculation, we have few options.
   * dcv has carried through from end of q loop above; store it 
   * in first pass, we add M->D and D->D path into DMX.
   */ 
  /* We're almost certainly obligated to do at least one complete 
   * DD path to be sure: 
   */
  dcv            = esl_vmx_rightshiftz_float(dcv);
  P7C_DQ(dpc, 0) = zerov;
  tp             = ((vector float *) om->tfv) + 7*Q;       /* set tp to start of the DD's */
  for (q = 0; q < Q; q++) 
    {
      P7C_DQ(dpc,q) = vec_add(dcv, P7C_DQ(dpc,q));   
      dcv           = vec_madd(P7C_DQ(dpc,q), *tp, zerov); tp++; /* extend DMO(q), so we include M->D and D->D paths */
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
    {                         // Fully serialized version 
      for (j = 1; j < 4; j++) // 4 = vector width, in floats
        {
          dcv = esl_vmx_rightshiftz_float(dcv);
          tp  = ((vector float *) om->tfv) + 7*Q;  /* reset tp to start of the DD's */
          for (q = 0; q < Q; q++) 
            { /* note, extend dcv, not DMO(q); only adding DD paths now */
              P7C_DQ(dpc,q) = vec_add(dcv, P7C_DQ(dpc,q));   
              dcv           = vec_madd(dcv, *tp, zerov);   tp++; 
            }       
        }
    } 
  else
    {                   /* Slightly parallelized version, but which incurs some overhead */
      for (j = 1; j < 4; j++)
        {
          register vector float cv = zerov;   /* keeps track of whether any DD's change DMO(q) */

          dcv = esl_vmx_rightshiftz_float(dcv);
          tp  = ((vector float *) om->tfv) + 7*Q;  /* set tp to start of the DD's */
          for (q = 0; q < Q; q++) 
            { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
              sv            = vec_add(dcv, P7C_DQ(dpc,q));         
             // This line gets moved into the if at the end of the loop because of the vec_any_gt function
              cv            = vec_or(cv, vec_cmpgt(sv, P7C_DQ(dpc,q))); 
              P7C_DQ(dpc,q) = sv;                                           /* store new DMO(q) */
              dcv           = vec_madd(dcv, *tp, zerov);   tp++;                 /* note, extend dcv, not DMO(q) */
            }       
          if (! vec_any_gt(cv, zerov)) break; /* DD's didn't change any DMO(q)? Then done, break out. */
        }
    }
  
  /* Add Dk's to xEv */
  for (q = 0; q < Q; q++) xEv = vec_add(P7C_DQ(dpc,q), xEv);

  /* Specials, in order: E N JJ J B CC C */
  xc[p7C_E]=esl_vmx_hsum_float(xEv);
  //esl_vmx_hsum_float(xEv, &xc[p7C_E]);  
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
      xEv = vec_splats((float) 1.0 / xc[p7C_E]);

      for (q = 0; q < Q; q++)
        {
          P7C_MQ(dpc,q) = vec_madd(P7C_MQ(dpc,q), xEv, zerov);
          P7C_DQ(dpc,q) = vec_madd(P7C_DQ(dpc,q), xEv, zerov);
          P7C_IQ(dpc,q) = vec_madd(P7C_IQ(dpc,q), xEv, zerov);
        }

      xc[p7C_SCALE] = xc[p7C_E];
      xc[p7C_E]     = 1.0f;
    }
  else xc[p7C_SCALE] = 1.0f;

  return logf(xc[p7C_SCALE]);
}


/* backward_row_main_vmx()
 * 
 * Backward calculation for rows 1..L-1. Rows 0, L are special cases.
 * 
 * <xi>        = residue dsq[i+1] on NEXT row
 * <om>        = query model
 * <dpp>       = 'previous' row i+1
 * <dpc>       = current row i
 * Q           = number of vectors in row; ox->Q
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
backward_row_main_vmx(ESL_DSQ xi, const P7_OPROFILE *om, vector float *dpp, vector float *dpc, int Q, float scalefactor)
{
  const vector float *rp       = (vector float *) om->rfv[xi];            // emission scores on row i+1, for bck; xi = dsq[i+1]  
  float       * const xc = (float *) (dpc + Q * p7C_NSCELLS); // E N JJ J B CC C SCALE 
  const float * const xp = (float *) (dpp + Q * p7C_NSCELLS);
  const vector float *tp, *tpdd;
  const vector float  zerov = {0.0, 0.0, 0.0, 0.0};
  vector float        xBv;
  vector float       *dp;
  vector float        xEv;
  vector float        dcv, mcv, ipv, mpv;
  vector float        tmmv, timv, tdmv;    /* copies of transition prob quads; a leftshift is needed as boundary cond */
  int           q;
 
  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. 
   * This invalidates <dpp> as a backwards row; its values are now
   * intermediates in the calculation of the current <dpc> row.
   */
  dp  = dpp;
  xBv = zerov;
  tp  = (vector float *) om->tfv;                /* on first transition vector */
  for (q = 0; q < Q; q++)
    {
      *dp = vec_madd(*dp, *rp, zerov); rp++;
      xBv = vec_madd(*dp, *tp, xBv); dp+= p7C_NSCELLS; tp += 7;
    }

  /* Specials. Dependencies dictate partial order C,CC,B < N,J,JJ < E */
  xc[p7C_C] = xc[p7C_CC] = xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  xc[p7C_B] = esl_vmx_hsum_float(xBv);
  //esl_vmx_hsum_float(xBv, &(xc[p7C_B]));
  xc[p7C_J] = xc[p7C_JJ] = xc[p7C_B] * om->xf[p7O_J][p7O_MOVE] + xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_N]              = xc[p7C_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_E]              = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_E][p7O_LOOP];

  /* Initialize for the row calculation */
  mpv  = esl_vmx_leftshiftz_float(*dpp); /* [1 5 9 13] -> [5 9 13 x], M(i+1,k+1) * e(M_k+1, x_{i+1}) */
  tmmv = esl_vmx_leftshiftz_float( ((vector float *) om->tfv)[1]);
  timv = esl_vmx_leftshiftz_float( ((vector float *) om->tfv)[2]);
  tdmv = esl_vmx_leftshiftz_float( ((vector float *) om->tfv)[3]);
  xEv  = vec_splats(xc[p7C_E]);
  tp   = ((vector float *) om->tfv) + 7*Q - 1;
  tpdd = tp + Q;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--)
    {
      ipv                 = P7C_IQ(dpp, q);
      P7C_IQ(dpc,q)       = vec_madd(ipv, *tp, vec_madd(mpv, timv, zerov)); tp--;   /* II,IM; I is done         */
      mcv                 = vec_madd(ipv, *tp, vec_madd(mpv, tmmv, zerov)); tp-=2;  /* MI,MM; ME,MD remain      */
      dcv                 = vec_madd(dcv, *tpdd, vec_madd(mpv, tdmv, zerov)); tpdd--; /* DM and one segment of DD */

      P7C_DQ(dpc,q) = dcv = vec_add( xEv, dcv);
      P7C_MQ(dpc,q)       = vec_add( xEv, mcv);

      mpv  = P7C_MQ(dpp, q);
      tdmv = *tp; tp--;
      timv = *tp; tp--;
      tmmv = *tp; tp-=2;
    }
  backward_row_finish_vmx(om, dpc, Q, dcv);
  backward_row_rescale_vmx(xc, dpc, Q, scalefactor);
}


/* backward_row_L_vmx()
 * 
 * Backward calculation for row L; 
 * a special case because the matrix <ox> has no 'previous' row L+1.
 * Otherwise identical to backward_row_main().
 */
static inline void
backward_row_L_vmx(const P7_OPROFILE *om,  vector float *dpc, int Q, float scalefactor)
{
  const vector float  zerov = {0.0, 0.0, 0.0, 0.0};
  float        *xc    = (float *) (dpc + Q * p7C_NSCELLS);
  const vector float *tpdd;
  vector float       *dp;
  vector float       xEv, dcv;
  int          q;

  /* Backwards from T <- C,CC <- E;  all other specials unreachable, impossible on row L.
   * specials are stored in order E N JJ J B CC C.  
   */
  xc[p7C_C] = xc[p7C_CC] = om->xf[p7O_C][p7O_MOVE];
  xc[p7C_B] = xc[p7C_J] = xc[p7C_JJ] = xc[p7C_N] = 0.0;
  xc[p7C_E] = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE];

  xEv  = vec_splats(xc[p7C_E]);
  dp   = dpc + Q*p7C_NSCELLS - 1;
  tpdd = ((vector float *) om->tfv) + 8*Q - 1;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--) 
    {
      *dp--       = zerov;                                            /* I */
      *dp-- = dcv = vec_madd(dcv, *tpdd, xEv); tpdd--;  /* D */
      *dp--       = xEv;                                              /* M */
    }
  backward_row_finish_vmx(om, dpc, Q, dcv);
  backward_row_rescale_vmx(xc, dpc, Q, scalefactor);
}



/* backward_row_finish_vmx()
 * 
 * Finishing DD and MD path contributions along the row.
 * Both L and main (1..L-1) recursions share this.
 * First segment of DD contributions has been calculated,
 * and up to three remain; then once D values are known,
 * we can add the M->D contribution to M.
 * 
 * Caller gets <dcv> from doing the first DD pass. 
 * <dcv> contains D scores for first vector, [1 5 9 13] in 
 * the example. We need it as an initialization condition;
 * we leftshift it to [5 9 13 *] and continue propagating
 * DD path.
 * 
 * The DD path calculation is the most annoying part of 
 * the otherwise gorgeous Farrar algorithm; it has to be
 * at least partially serialized, if not fully.
 * 
 * <om>  = query model
 * <dpc> = current row, partially calculated
 * Q     = width of rows in # vectors; ox->Q
 * dcv   = the first D vector [1 5 9 13] from caller's earlier calculation
 */
static inline void
backward_row_finish_vmx(const P7_OPROFILE *om, vector float *dpc, int Q, vector float dcv)
{
  const vector float zerov = {0.0, 0.0, 0.0, 0.0};
  const vector float *tp;
  vector float       *dp;
  int           j,q;
  
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
          dcv = esl_vmx_leftshiftz_float(dcv);     /* [1 5 9 13] => [5 9 13 *]          */
          tp  = ((vector float *) om->tfv) + 8*Q - 1;   /* <*tp> now the [4 8 12 x] TDD quad */
          dp  = dpc + Q*p7C_NSCELLS - 2;          /* init to point at D(i,q) vector    */
          for (q = Q-1; q >= 0; q--)
            {
              dcv = vec_madd(dcv, *tp, zerov); tp--;
              *dp = vec_add(*dp, dcv); dp -= p7C_NSCELLS;
            }
        }
    }
  else
    { /* With check for early convergence */
      vector float sv;
      vector float cv;  /* keeps track of whether any DD addition changes DQ(q) value */
      for (j = 1; j < 4; j++) // 4 = vector width in floats
        {
          dcv = esl_vmx_leftshiftz_float(dcv);
          tp  = ((vector float *) om->tfv) + 8*Q - 1;      
          dp  = dpc + Q*p7C_NSCELLS - 2;
          cv  = zerov;
          for (q = Q-1; q >= 0; q--)
            { /* using cmpgt below tests if DD changed any DMO(q) without conditional branch (i.e. no if) */
              dcv  = vec_madd(dcv, *tp, zerov); tp--;
              sv   = vec_add(*dp, dcv);
              cv   = vec_or(cv, vec_cmpgt(sv, *dp)); /* if DD path changed DQ(dpc,q), cv bits know it now */
              *dp  = sv; 
              dp  -= p7C_NSCELLS;
            }
          if (! vec_any_gt(cv, zerov)) break; /* if no DD path changed DQ(q) in this segment, then done, no more segments needed */
        }
    }

  /* Finally, M->D path contribution
   * these couldn't be added to M until we'd finished calculating D values on row.
   */
  dcv = esl_vmx_leftshiftz_float(P7C_DQ(dpc, 0));
  tp  = ((vector float *) om->tfv) + 7*Q - 3;       
  dp  = dpc + (Q-1)*p7C_NSCELLS; 
  for (q = Q-1; q >= 0; q--)
    {
      *dp  = vec_madd(dcv, *tp, *dp); tp -= 7; 
      dcv  = *(dp+1);                               dp -= p7C_NSCELLS;
    }
}



/* backward_row_rescale_vmx()
 * 
 * Sparse rescaling, using the scalefactor that Forward set and 
 * Backward shares.
 * 
 * xc          - ptr to specials on row (floats)
 * dpc         - ptr to vectors for row      
 * Q           - number of vectors; ox->Q
 * scalefactor - scalefactor that Forward set for this row
 *
 * Upon return, values in current row <dpc> have been rescaled.
 */
static inline void
backward_row_rescale_vmx(float *xc, vector float *dpc, int Q, float scalefactor)
{
  if (scalefactor > 1.0f)
    {
      vector float zerov = {0.0, 0.0, 0.0, 0.0};
      vector float sv = vec_splats((float) 1.0 / scalefactor);
      vector float *dp = dpc;
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
          *dp = vec_madd(*dp, sv, zerov); dp++; /* M */
          *dp = vec_madd(*dp, sv, zerov); dp++; /* D */
          *dp = vec_madd(*dp, sv, zerov); dp++; /* I */
        }
    }
  xc[p7C_SCALE] = scalefactor;
}


/* posterior_decode_row_vmx()
 *
 * In production code, we don't have to save results of posterior
 * decoding; we immediately use them to calculate sparse mask on the row.
 *
 * In debugging code, we save pp's, overwriting the forward row as we
 * go; when all pp's on the row have been calculated, if there's a
 * ox->pp matrix, we copy the pp values there.
 * 
 * Understanding how sparse rescaling interacts with posterior decoding:
 * Let F_i, B_i be the correct (unscaled) forward/backward terms for a
 * given cell on row i that we would've calculated, had we had
 * infinite numerical range. Let T be the overall Forward or Backward
 * score of the sequence comparison; we call it T, because we'll use
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
 * <p7_SPARSIFY_THRESH>, in <p7_config.h.in>). We take
 * advantage of the fact that if N/C/J emission postprobs exceed
 * <1.0 - sm_thresh>, we don't even need to look at the vectorized
 * row; no cell can exceed threshold.
 *
 * Can throw <eslEINVAL> on bad code (something awry in data structure initialization)
 *           <eslEMEM> on allocation failure in sparsemask
 */
static inline int
posterior_decode_row_vmx(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc)
{
  int             Q        = ox->Q;
  vector float        *fwd       = (vector float *) ox->dpf[ox->R0 + ox->R]; /* a calculated fwd row R has been popped off */
  const  vector float *bck       = (vector float *) ox->dpf[rowi%2];
  float         *xf        = ox->dpf[ox->R0 + ox->R] + ox->Vf * Q * p7C_NSCELLS;
  const  float  *xb        = ox->dpf[rowi%2]         + ox->Vf * Q * p7C_NSCELLS;
  const vector float   threshv   = vec_splats(sm_thresh);
  float          scaleterm = xf[p7C_SCALE] / overall_sc; /* see comments above, on how rescaling affects posterior decoding equations */
  const vector float   cv        = vec_splats(scaleterm);
  float  pnonhomology;
  vector bool mask;
  vector float pv;
  int    q;
  int    status;
  vector float zerov = {0.0, 0.0, 0.0, 0.0};

  // Extracting a conditional from a vector is very hard in VMX, so instead we define these four selection vectors
  // that we'll use in conjunction with vec_any_gt to determine if the corresponding element of a vector is = 1
  vector bool sel1 = {1, 0, 0, 0}; 
  vector bool sel2 = {0, 1, 0, 0};
  vector bool sel3 = {0, 0, 1, 0}; 
  vector bool sel4 = {0, 0, 0, 1};
  vector bool zeroes = {0, 0, 0, 0}; 
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
  pnonhomology = (xf[p7C_N] * xb[p7C_N] + xf[p7C_JJ] * xb[p7C_JJ] + xf[p7C_CC] * xb[p7C_CC]) * scaleterm;
  if (pnonhomology <= 1.0f - sm_thresh)
    {
      if ((status = p7_sparsemask_StartRow(sm, rowi)) != eslOK) return status;
      for (q = Q-1; q >= 0; q--)                   // reverse, because SPARSEMASK is entirely in reversed order 
        {
          pv       =                vec_madd(P7C_MQ(fwd, q), P7C_MQ(bck, q), zerov);
          pv       = vec_madd(P7C_IQ(fwd, q), P7C_IQ(bck, q), pv);
          pv       = vec_madd(P7C_DQ(fwd, q), P7C_DQ(bck, q), pv);
          pv       = vec_madd(pv, cv, zerov);           // pv is now the posterior probability of elements q,r=0..3 
          mask     = vec_cmpge(pv, threshv);    // mask now has all 0's in elems r that failed thresh; all 1's for r that passed 

          // Painful unrolled check for whether each element of mask is != 0
          if(vec_any_ne(vec_and(mask, sel1), zeroes)){
            if ((status = p7_sparsemask_Add(sm, q, 0)) != eslOK) return status;
          }
          if(vec_any_ne(vec_and(mask, sel2), zeroes)){
            if ((status = p7_sparsemask_Add(sm, q, 1)) != eslOK) return status;
          }
          if(vec_any_ne(vec_and(mask, sel3), zeroes)){
            if ((status = p7_sparsemask_Add(sm, q, 2)) != eslOK) return status;
          }
          if(vec_any_ne(vec_and(mask, sel4), zeroes)){
            if ((status = p7_sparsemask_Add(sm, q, 3)) != eslOK) return status;
          }

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
      P7C_MQ(fwd, q) = vec_madd(cv, vec_madd(P7C_MQ(fwd, q), P7C_MQ(bck, q), zerov), zerov);
      P7C_DQ(fwd, q) = vec_madd(cv, vec_madd(P7C_DQ(fwd, q), P7C_DQ(bck, q), zerov), zerov);
      P7C_IQ(fwd, q) = vec_madd(cv, vec_madd(P7C_IQ(fwd, q), P7C_IQ(bck, q), zerov), zerov);
    }

  if (ox->pp)  save_debug_row_pp_vmx(ox, fwd, rowi);
#endif
  return eslOK;
}
/*------------------ end, inlined recursions -------------------*/



/*****************************************************************
 * 3. More internal functions: debugging tools
 *****************************************************************/
#if eslDEBUGLEVEL > 0

/* backward_row_zero_vmx()
 * 
 * Slightly peculiar but true: in production code we don't
 * need to calculate backward row 0, because we're only doing
 * posterior decoding on residues 1..L.  We only need row 0
 * if we need a complete Backward calculation -- which happens
 * when we're in debugging mode, and we're going to test the
 * Backwards score or compare the Backwards matrix to the
 * reference implementation.
 * 
 * x1  - residue dsq[1] on row 1
 * om  - query model
 * ox  - DP matrix; dpf[0] is the current Backward row
 * 
 * Upon return, ox->dpf[0] (backward row 0) has been calculated.
 * Returns the log of the value in the N cell at row 0. When this 
 * is added to the sum of the logs of all the scalefactors (which
 * the caller has been accumulating in <ox->bcksc>), then
 * the <ox->bcksc> value is finished and equal to the Backwards
 * raw score in nats.
 */
static inline float
backward_row_zero_vmx(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox)
{
  int          Q     = ox->Q;
  vector float       *dpc  = (vector float *) ox->dpf[0];
  vector float       *dpp  = (vector float *) ox->dpf[1];
  const vector float *rp   = (vector float *) om->rfv[x1];
  const vector float zerov = {0.0, 0.0, 0.0, 0.0};
  float        *xc   = (float *) (dpc + Q * p7C_NSCELLS); /* special states on current row i      */
  float        *xp   = (float *) (dpp + Q * p7C_NSCELLS); /* special states on "previous" row i+1 */
  vector float       *dp;
  vector float       *tp;
  vector float        xBv  = zerov;
  int           q;

  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. */
  dp  = dpp;
  tp  = (vector float *) om->tfv;
  for (q = 0; q < Q; q++)
    {
      *dp = vec_madd(*dp, *rp, zerov); rp++;
      xBv = vec_madd(*dp, *tp, xBv); dp+= p7C_NSCELLS; tp += 7;
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
  xc[p7C_B] = esl_vmx_hsum_float(xBv);
  xc[p7C_J] = xc[p7C_JJ] = xc[p7C_B] * om->xf[p7O_J][p7O_MOVE] + xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_N]              = xc[p7C_B] * om->xf[p7O_N][p7O_MOVE] + xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_E]              = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_E][p7O_LOOP];
  
  /* not needed in production code: */
  for (q = 0; q < Q; q++)
    P7C_MQ(dpc, q) = P7C_IQ(dpc, q) = P7C_DQ(dpc, q) = zerov;

  return logf(xc[p7C_N]);
}

static void
save_debug_row_pp_vmx(P7_CHECKPTMX *ox, vector float *dpc, int i)
{
  union { vector float v; float x[4]; } u;  // 4 = number of floats per SSE vector
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
    { // 4 = number of floats per SSE vector
      u.v = P7C_MQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_ML) = u.x[z]; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_DL) = u.x[z]; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_IL) = u.x[z]; }
    }
}


/* save_debug_row_fb_vmx()
 * 
 * Debugging only. Transfer posterior decoding values (sparse scaled,
 * prob space) from a vectorized row, to appropriate row of <ox->fwd>
 * or <ox->bck> (log space, inclusive of partial sum of scalefactors);
 * <ox->fwd> and <ox->bck> should be identical (within numerical error
 * tolerance) to a reference implementation Forward/Backward in log
 * space.
 */
static void
save_debug_row_fb_vmx(P7_CHECKPTMX *ox, P7_REFMX *gx, vector float *dpc, int i, float totscale)
{
  union { vector float v; float x[4]; } u;  // 4 = number of floats per SSE vector
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
    { // 4 = number of floats per SSE vector
      u.v = P7C_MQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_ML) = logf(u.x[z]) + totscale; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_DL) = logf(u.x[z]) + totscale; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < 4; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_IL) = logf(u.x[z]) + totscale; }
    }

}

#endif // eslDEBUGLEVEL > 0
/*---------------------------- end, debugging tools ---------------------------*/

#else // ! eslENABLE_VMX
/* provide callable functions even when we're `./configure --disable-vmx` */
int
p7_ForwardFilter_vmx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(om); ESL_UNUSED(ox); ESL_UNUSED(opt_sc); 
  esl_fatal("Altivec/VMX support was not enabled at compile time. Can't use p7_ForwardFilter_vmx().");
  return eslFAIL; // NOTREACHED
}
int
p7_BackwardFilter_vmx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(om); ESL_UNUSED(ox); ESL_UNUSED(sm); ESL_UNUSED(sm_thresh); 
  esl_fatal("Altivec/VMX support was not enabled at compile time. Can't use p7_BackwardFilter_vmx().");
  return eslFAIL; // NOTREACHED
}
/* Standard compiler-pleasing mantra for an #ifdef'd-out, empty code file. */
void p7_fwdfilter_vmx_silence_hack(void) { return; }
#if defined p7FWDFILTER_VMX_TESTDRIVE || p7FWDFILTER_VMX_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_VMX or not
