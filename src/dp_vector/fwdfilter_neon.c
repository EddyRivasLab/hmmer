/* Forwards/Backwards filters: ARM NEON version.
 *
 * Ported from Intel/SSE: Tyler Camp (University of Texas, Austin)
 * See Intel/SSE version for general notes.
 *
 * Contents:
 *    1. ForwardFilter() and BackwardFilter() API.
 *    2. Internal functions: inlined recursions.
 *    3. Debugging tools.
 *    4. Stats driver (memory requirements)
 *    5. Benchmark driver.
 *    6. Unit tests.
 *    7. Test driver.
 *    8. Example.
 */
#include "p7_config.h"
#if p7_CPU_ARCH == arm || p7_CPU_ARCH == arm64
#include <arm_neon.h>
#endif
#include "easel.h"
#include "esl_vectorops.h"
#include "esl_neon.h"

#include "dp_reference/p7_refmx.h"
#include "dp_sparse/p7_sparsemx.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/fwdfilter.h"


/* vectorized DP recursions are tediously lengthy, so for some
 * semblance of clarity, they're broken out into one-page-ish
 * chunks, using static inlined functions.
 */
#ifdef HAVE_NEON
static inline float forward_row      (ESL_DSQ xi, const P7_OPROFILE *om, const esl_neon_128f_t *dpp, esl_neon_128f_t *dpc, int Q);
static inline void  backward_row_main(ESL_DSQ xi, const P7_OPROFILE *om,       esl_neon_128f_t *dpp, esl_neon_128f_t *dpc, int Q, float scalefactor);
static inline void  backward_row_L   (            const P7_OPROFILE *om,                             esl_neon_128f_t *dpc, int Q, float scalefactor);
static inline void  backward_row_finish(          const P7_OPROFILE *om,                             esl_neon_128f_t *dpc, int Q, esl_neon_128f_t dcv);
static inline void  backward_row_rescale(float *xc, esl_neon_128f_t *dpc, int Q, float scalefactor);
static inline int   posterior_decode_row(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc);
#endif
#ifdef p7_DEBUGGING
static inline float backward_row_zero(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox);
static        void  save_debug_row_pp(P7_CHECKPTMX *ox,               debug_print *dpc, int i);
static        void  save_debug_row_fb(P7_CHECKPTMX *ox, P7_REFMX *gx, debug_print *dpc, int i, float totscale);
#endif

/*****************************************************************
 * 1. Forward and Backward API calls
 *****************************************************************/

/* Function:  p7_ForwardFilter()
 * Synopsis:  Checkpointed striped vector Forward calculation, producing score.
 *
 * Purpose:   Calculate the Forward algorithm for target sequence <dsq>
 *            of <L> residues aligned to query profile <om>, using the
 *            checkpointed DP matrix <ox> provided by the caller. Upon
 *            successful return, <ox> contains the filled Forward
 *            matrix, and <*opt_sc> optionally contains the raw Forward
 *            score in nats.
 *            
 *            <ox> will be reallocated, if needed, for the MxL problem;
 *            caller does not need to call <p7_checkptmx_GrowTo()> itself.
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
 * Throws:    <eslEMEM> on reallocation error.
 * 
 * Xref:      For layout of checkpointed <ox> see exegesis in p7_checkptmx.h.
 */
int
p7_ForwardFilter_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc)
{
  #ifdef HAVE_NEON
  int              Q     = P7_NVF(om->M);   /* segment length; # of MDI vectors on each row      */
  esl_neon_128f_t *dpp   = NULL;            /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  esl_neon_128f_t *dpc   = NULL;	    /* dpc points at current row                         */
  float           *xc    = NULL;            /* specials E,N,JJ,J,B,CC,C,SCALE                    */
  float            totsc = 0.0f;	    /* accumulates Forward score in nats                 */
  esl_neon_128f_t  place;
  place.f32x4 = vdupq_n_f32(0);		
  const esl_neon_128f_t zerov = place; 
  int     q ;			            /* counter over vectors 0..Q-1                        */
  int     i;			            /* counter over residues/rows 1..L                    */
  int     b;			            /* counter down through checkpointed blocks, Rb+Rc..1 */
  int     w;			            /* counter down through rows in a checkpointed block  */


  /* Contract checks */
  //  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  //  ESL_DASSERT1(( om->L    == L ));	      /*  ... and it's easy to forget to set <om> that way           */
  //  ESL_DASSERT1(( om->nj   == 1.0f ));     /*  ... hence the check                                        */
                                              /*  ... which you can disable, if you're playing w/ config     */
  // SRE: disabled 21 Jul 14 because unit tests in sparse_asc_fwdback need to create
  // sparse mask w/ checkpointed F/B w/ unilocal L=0 models

  /* Make sure <ox> is allocated big enough.
   * DO NOT set any ptrs into the matrix until after this potential reallocation!
   */
  p7_checkptmx_GrowTo(ox, om->M, L);
  dpp =  (esl_neon_128f_t *) ox->dpf[ox->R0-1];    /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  xc  =  (float *) (dpp + Q*p7C_NSCELLS);          /* specials E,N,JJ,J,B,CC,C,SCALE    */

  /* Set the size of the problem in <ox> now, not later
   * Debugging dumps need this information, for example
   * (Ditto for any debugging copy of the fwd mx)
   */
  ox->M  = om->M;	
  ox->L  = L;
  ox->Qf = Q;
#ifdef p7_DEBUGGING
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
#ifdef p7_DEBUGGING
  if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, 0, dpp, "f1 O"); 
  if (ox->fwd)        save_debug_row_fb(ox, ox->fwd, dpp, 0, totsc); 
#endif

  /* Phase one: the "a" region: all rows in this region are saved */
  for (i = 1; i <= ox->La; i++)
    {
      dpc = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R]; ox->R++;    /* idiomatic for "get next save/checkpoint row" */

      totsc += forward_row(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;	    	                          /* current row becomes prev row */
#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, dpc, "f1 O"); 
      if (ox->fwd)        save_debug_row_fb(ox, ox->fwd, dpc, i, totsc); 
#endif
    }

  /* Phase two: "b" and "c" regions: partially and fully checkpointed */
  for (b = ox->Rb + ox->Rc, w = (ox->Rb ? ox->Lb : ox->Rc+1); i <= L; i++)
    {
      /* this section, deciding whether to get a checkpointed vs. tmp
       * row, is why we set <fwd>, <dpp> here, rather than having the
       * inlined forward_row() set them.
       */
      if (! (--w)) { 		                   /* we're on the last row in segment: this row is saved    */
	dpc = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R]; ox->R++;  /* idiomatic for "get next save/checkpoint row"    */
	w = b;  			           /* next segment has this many rows, ending in a saved row */
	b--;					   /* decrement segment number counter; last segment is r=1  */
      } else dpc = (esl_neon_128f_t *) ox->dpf[i%2];        /* idiomatic for "next tmp row", 0/1; i%2 makes sure dpp != dpc */
      
      totsc += forward_row(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;
#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, dpc, w ? "f1 X" : "f1 O"); 
      if (ox->fwd)        save_debug_row_fb(ox, ox->fwd, dpc, i, totsc); 
#endif
    }

  xc     = (float *) (dpc + Q*p7C_NSCELLS);

  ESL_DASSERT1( (ox->R == ox->Ra+ox->Rb+ox->Rc) );
  ESL_DASSERT1( ( (! isnan(xc[p7C_C])) && (! isinf(xc[p7C_C]))) );

  if (opt_sc) *opt_sc = totsc + logf(xc[p7C_C] * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
  #endif /* HAVE_NEON */
  #ifndef HAVE_NEON
    return eslNORESULT;
  #endif
}


/* Function:  p7_BackwardFilter()
 * Synopsis:  Checkpointed striped vector Backward calculation, producing sparse mask.
 *
 * Purpose:   Given a target sequence <dsq> of length <L>, a query model
 *            <om>, and a DP matrix <ox> resulting from a successful
 *            call to <p7_ForwardFilter()>; calculate the Backward and
 *            posterior decoding algorithms. On each row <i=1..L>, use
 *            posterior decoding to determine which <i,k> cells pass a
 *            significance threshold for inclusion in the sparse DP
 *            mask. Store those sparse cells in <sm>, which the caller
 *            allocates (or reuses) and provides. <sm> will be resized
 *            here as needed; caller does not need to call 
 *            <p7_sparsemask_Reinit()> on it.
 *            
 *            
 * Args:      dsq       - digital target sequence, 1..L
 *            L         - length of dsq, residues
 *            om        - optimized profile (multihit local)
 *            ox        - checkpointed DP matrix, ForwardFilter already run
 *            sm        - allocated P7_SPARSEMASK structure to hold sparse DP mask
 *            sm_thresh - Threshold for determining 'significant' posterior
 *                        alignment probability, passed in to posterior_decode_row()
 *
 * Throws:    <eslEINVAL> if something's awry with a data structure's internals.
 *            <eslEMEM> on allocation failure.
 */
int
p7_BackwardFilter_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh)
{
  #ifdef HAVE_NEON
  int              Q = ox->Qf;
  esl_neon_128f_t *fwd;
  esl_neon_128f_t *bck;
  esl_neon_128f_t *dpp;
  float           *xf;
  float            Tvalue;
  int              i, b, w, i2;
  int              status;

  p7_sparsemask_Reinit(sm, om->M, L);

  /* Contract checks */
  //ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  //ESL_DASSERT1(( om->L    == L ));	    /*  ... and it's easy to forget to set <om> that way           */
  //ESL_DASSERT1(( om->nj   == 1.0f ));	    /*  ... hence the check                                        */
                                            /*  ... which you can disable, if you're playing w/ config     */
  // SRE: disabled 21 Jul 14 because unit tests in sparse_asc_fwdback need to create
  // sparse mask w/ checkpointed F/B w/ unilocal L=0 models

#ifdef p7_DEBUGGING
  /* Debugging instrumentations. */
   if (ox->bck) { ox->bck->M = om->M; ox->bck->L = L; ox->bck->type = p7R_BACKWARD; }
   if (ox->pp)  { ox->pp->M  = om->M; ox->pp->L  = L; ox->pp->type  = p7R_DECODING; }
#endif

  /* Row L is a special case for Backwards; so in checkpointed Backward,
   * we have to special-case the last block, rows L-1 and L
   */
  i = L;
  ox->R--;
  fwd = (esl_neon_128f_t *) ox->dpf[ox->R0 + ox->R];     /* pop row for fwd[L] off the checkpointed stack */
  xf  = (float *) (fwd + Q*p7C_NSCELLS);
  Tvalue = xf[p7C_C] * om->xf[p7O_C][p7O_MOVE];          /* i.e. scaled fwd[L] val at T state = scaled overall score */
  bck = (esl_neon_128f_t *) ox->dpf[i%2];	         /* get tmp space for bck[L]                                 */
  backward_row_L(om, bck, Q, xf[p7C_SCALE]);             /* calculate bck[L] row                                     */
#ifdef p7_DEBUGGING
  ox->bcksc = logf(xf[p7C_SCALE]);
  if (ox->do_dumping) { p7_checkptmx_DumpFBRow(ox, L, fwd, "f2 O"); if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, L, bck, "bck");  }
  if (ox->bck)          save_debug_row_fb(ox, ox->bck, bck, L, ox->bcksc); 
#endif
  if ( (status = posterior_decode_row(ox, i, sm, sm_thresh, Tvalue)) != eslOK)  return status;
  i--;
  dpp = bck;


  /* If there's any checkpointing, there's an L-1 row to fill now. */
  if (ox->Rb+ox->Rc > 0)
    {
      /* Compute fwd[L-1] from last checkpoint, which we know is fwd[L-2] */
      dpp = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R-1];  /* fwd[L-2] values, already known        */
      fwd = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R];    /* get free row memory from top of stack */
      forward_row(dsq[i], om, dpp, fwd, Q);               /* calculate fwd[L-1]                    */
#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, fwd, "f2 X");
#endif

      /* Compute bck[L-1] from bck[L]. */
      xf  = (float *) (fwd + Q*p7C_NSCELLS);
      dpp = (esl_neon_128f_t *) ox->dpf[(i+1)%2]; 
      bck = (esl_neon_128f_t *) ox->dpf[i%2];             /* get space for bck[L-1]                */
      backward_row_main(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#ifdef p7_DEBUGGING
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, bck, "bck");
      if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i, ox->bcksc); 
#endif
      /* And decode. */
      if ( (status = posterior_decode_row(ox, i, sm, sm_thresh, Tvalue)) != eslOK)  return status;
	   dpp = bck;
      i--;			/* i is now L-2 if there's checkpointing; else it's L-1 */
    }

  /* Main loop for checkpointed regions (b,c) */
   for (b = 2; b <= ox->Rb+ox->Rc; b++)
    {				/* i=L-2 as we enter here, and <dpp> is on bck[L-1] */
      w = (b <= ox->Rc ? b+1 : ox->Lb);

      /* We know current row i (r=R0+R-1) ends a block and is checkpointed in fwd. */
      ox->R--;
      fwd = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R];      /* pop checkpointed forward row off "stack" */
      xf  = (float *) (fwd + Q*p7C_NSCELLS);

      /* Calculate bck[i]; <dpp> is already bck[i+1] */
      bck = (esl_neon_128f_t *) ox->dpf[i%2];	    /* get available tmp memory for row     */
      backward_row_main(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#ifdef p7_DEBUGGING
      ox->bcksc += logf(xf[p7C_SCALE]);
      if (ox->do_dumping) { p7_checkptmx_DumpFBRow(ox, i, fwd, "f2 O");	p7_checkptmx_DumpFBRow(ox, i, bck, "bck"); }
      if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i, ox->bcksc); 
#endif
      /* And decode checkpointed row i. */
      if ( (status = posterior_decode_row(ox, i, sm, sm_thresh, Tvalue)) != eslOK)  return status;
      /* The rest of the rows in the block weren't checkpointed.
       * Compute Forwards from last checkpoint ...
       */
      dpp = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R-1];       /* get last Fwd checkpoint. */
      for (i2 = i-w+1; i2 <= i-1; i2++)
	{
	  fwd = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R]; ox->R++;  /* push new forward row on "stack"     */
	  forward_row(dsq[i2], om, dpp, fwd, Q);
	  dpp = fwd;	  
	}

      /* ... and compute Backwards over the block we just calculated, while decoding. */
      dpp = bck;
      for (i2 = i-1; i2 >= i-w+1; i2--)
	{
	  ox->R--;
	  fwd = (esl_neon_128f_t *) ox->dpf[ox->R0+ox->R]; /* pop just-calculated forward row i2 off "stack" */
	  xf  = (float *) (fwd + Q*p7C_NSCELLS);
	  bck = (esl_neon_128f_t *) ox->dpf[i2%2];	  /* get available for calculating bck[i2]          */
	  backward_row_main(dsq[i2+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#ifdef p7_DEBUGGING
	  ox->bcksc += logf(xf[p7C_SCALE]);
	  if (ox->do_dumping) { p7_checkptmx_DumpFBRow(ox, i2, fwd, "f2 X"); p7_checkptmx_DumpFBRow(ox, i2, bck, "bck"); }
	  if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i2, ox->bcksc); 
#endif
	  if ((status = posterior_decode_row(ox, i2, sm, sm_thresh, Tvalue)) != eslOK)  return status;
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
       bck = (esl_neon_128f_t *) ox->dpf[i%2];	       /* get open space for bck[i]               */
       backward_row_main(dsq[i+1], om, dpp, bck, Q, xf[p7C_SCALE]);
#ifdef p7_DEBUGGING
       ox->bcksc += logf(xf[p7C_SCALE]);
       if (ox->do_dumping) { p7_checkptmx_DumpFBRow(ox, i, fwd, "f2 O"); p7_checkptmx_DumpFBRow(ox, i, bck, "bck"); }
       if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i, ox->bcksc); 
#endif
       if ((status = posterior_decode_row(ox, i, sm, sm_thresh, Tvalue)) != eslOK)  return status;
       dpp = bck;
     }

#ifdef p7_DEBUGGING
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
   xN = backward_row_zero(dsq[1], om, ox); 
   if (ox->do_dumping) { p7_checkptmx_DumpFBRow(ox, 0, fwd, "f2 O"); p7_checkptmx_DumpFBRow(ox, 0, bck, "bck"); }
   if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, 0, ox->bcksc); 
   if ((status = posterior_decode_row(ox, 0, sm, sm_thresh, Tvalue)) != eslOK)  return status;
   ox->bcksc += xN;
#endif

   p7_sparsemask_Finish(sm);
   return eslOK;
  #endif /* HAVE_NEON */
  #ifndef HAVE_NEON
    return eslNORESULT;
  #endif
}
/*----------- end forward/backward API calls --------------------*/



/*****************************************************************
 * 2. Internal functions: inlined recursions
 *****************************************************************/

/* forward_row()
 * 
 * <xi>   dsq[i]; residue on this row.
 * <om>   query model
 * <dpp>  ptr to previous row; may be a tmp row or a checkpointed row.
 * <dpc>  ptr to current row; ditto.
 * Q      number of vectors in <dpp>, <dpc>; ox->Qf.
 * 
 * Upon return, <dpc> contains scaled Forward values.
 * Returns log(scalevalue), for the caller to accumulate
 * as it accumulates the total Forward score.
 */
#ifdef HAVE_NEON
static inline float
forward_row_neon(ESL_DSQ xi, const P7_OPROFILE *om, const esl_neon_128f_t *dpp, esl_neon_128f_t *dpc, int Q)
{
  const  esl_neon_128f_t *rp   = om->rfv[xi];
  esl_neon_128f_t         place;
  place.f32x4 = vdupq_n_f32(0);

  const    esl_neon_128f_t  zerov = place;
  const    esl_neon_128f_t *tp    = om->tfv;
  const    float           *xp    = (float *) (dpp + Q * p7C_NSCELLS);
  float                    *xc    = (float *) (dpc + Q * p7C_NSCELLS);
  esl_neon_128f_t           dcv;
     dcv.f32x4 = vdupq_n_f32(0); 
  esl_neon_128f_t           xEv;
     xEv.f32x4 = vdupq_n_f32(0);  
  esl_neon_128f_t           xBv;
     xBv.f32x4 = vdupq_n_f32(xp[p7C_B]);   
  esl_neon_128f_t mpv, dpv, ipv;
  esl_neon_128f_t sv;
  int    q;
  int    j;

  mpv = esl_neon_rightshift_float(P7C_MQ(dpp, Q-1), zerov); 
  //union { esl_neon_128f_t v; float f[4]; }ok;
  //ok.v = mpv;
  // printf("before: %f %f %f %f\n", ok.f[0], ok.f[1], ok.f[2], ok.f[3]);
 
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
      xEv.f32x4    = vaddq_f32(xEv.f32x4, sv.f32x4);			           /* Mk->E    */

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
  //union { esl_neon_128f_t v; float f[4]; }ok;
  //printf("before: %f %f %f %f\n", ok.f[0], ok.f[1], ok.f[2], ok.f[3]);
  dcv            = esl_neon_rightshift_float(dcv, zerov);
  // printf("after: %f %f %f %f\n", ok.f[0], ok.f[1], ok.f[2], ok.f[3]);
  P7C_DQ(dpc, 0) = zerov;
  tp             = om->tfv + 7*Q;	/* set tp to start of the DD's */
  for (q = 0; q < Q; q++) 
    {
      P7C_DQ(dpc,q).f32x4 = vaddq_f32(dcv.f32x4, P7C_DQ(dpc,q).f32x4);	
      dcv.f32x4           = vmulq_f32(P7C_DQ(dpc,q).f32x4, (*tp).f32x4); tp++; /* extend DMO(q), so we include M->D and D->D paths */
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
	  dcv = esl_neon_rightshift_float(dcv, zerov);
	  tp  = om->tfv + 7*Q;	/* reset tp to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    { /* note, extend dcv, not DMO(q); only adding DD paths now */
	      P7C_DQ(dpc,q).f32x4 = vaddq_f32(dcv.f32x4, P7C_DQ(dpc,q).f32x4);	
	      dcv.f32x4           = vmulq_f32(dcv.f32x4, (*tp).f32x4);   tp++; 
	    }	    
	}
    } 
  else
    {			/* Slightly parallelized version, but which incurs some overhead */
      for (j = 1; j < 4; j++)
	{
	  register esl_neon_128f_t cv = zerov;	/* keeps track of whether any DD's change DMO(q) */

	  dcv = esl_neon_rightshift_float(dcv, zerov);
	  tp  = om->tfv + 7*Q;	/* set tp to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
	      sv.f32x4            = vaddq_f32(dcv.f32x4, P7C_DQ(dpc,q).f32x4);
		  esl_neon_128i_t tmp;
		  tmp.u32x4 = vreinterpretq_u32_f32(cv.f32x4);	
	      tmp.u32x4 = vorrq_u32(vreinterpretq_u32_f32(cv.f32x4), vcgtq_f32(sv.f32x4, P7C_DQ(dpc,q).f32x4)); 
	      cv.f32x4 = vreinterpretq_f32_u32(tmp.u32x4); 
		  P7C_DQ(dpc,q) = sv;	                                    /* store new DMO(q) */
	      dcv.f32x4           = vmulq_f32(dcv.f32x4, (*tp).f32x4);   tp++;                 /* note, extend dcv, not DMO(q) */
	    }	    
	  union { esl_neon_128f_t v; uint32_t i[4]; }s;
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
  if (xc[p7C_E] > 1.0e4)	/* that's a little less than e^10, ~10% of our dynamic range */
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
#endif

/* backward_row_main()
 * 
 * Backward calculation for rows 1..L-1. Rows 0, L are special cases.
 * 
 * <xi>        = residue dsq[i+1] on NEXT row
 * <om>        = query model
 * <dpp>       = 'previous' row i+1
 * <dpc>       = current row i
 * Q           = number of vectors in row; ox->Qf
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

#ifdef HAVE_NEON
static inline void
backward_row_main_neon(ESL_DSQ xi, const P7_OPROFILE *om, esl_neon_128f_t *dpp, esl_neon_128f_t *dpc, int Q, float scalefactor)
{
  const esl_neon_128f_t *rp       = om->rfv[xi];			      /* emission scores on row i+1, for bck; xi = dsq[i+1]  */
  float       * const xc = (float *) (dpc + Q * p7C_NSCELLS); /* E N JJ J B CC C SCALE */
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
  esl_neon_128f_t tmmv, timv, tdmv;				   /* copies of transition prob quads; a leftshift is needed as boundary cond */
 
  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. 
   * This invalidates <dpp> as a backwards row; its values are now
   * intermediates in the calculation of the current <dpc> row.
   */
  dp  = dpp;
  xBv = zerov;
  tp  = om->tfv;		/* on first transition vector */
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
  tmmv = esl_neon_leftshift_float(om->tfv[1], zerov);
  timv = esl_neon_leftshift_float(om->tfv[2], zerov);
  tdmv = esl_neon_leftshift_float(om->tfv[3], zerov);
  xEv.f32x4  = vdupq_n_f32(xc[p7C_E]);
  tp   = om->tfv + 7*Q - 1;
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
  backward_row_finish(om, dpc, Q, dcv);
  backward_row_rescale(xc, dpc, Q, scalefactor);
}
#endif

/* backward_row_L()
 * 
 * Backward calculation for row L; 
 * a special case because the matrix <ox> has no 'previous' row L+1.
 * Otherwise identical to backward_row_main().
 */

#ifdef HAVE_NEON
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
  int          q;

  /* Backwards from T <- C,CC <- E;  all other specials unreachable, impossible on row L.
   * specials are stored in order E N JJ J B CC C.  
   */
  xc[p7C_C] = xc[p7C_CC] = om->xf[p7O_C][p7O_MOVE];
  xc[p7C_B] = xc[p7C_J] = xc[p7C_JJ] = xc[p7C_N] = 0.0;
  xc[p7C_E] = xc[p7C_C] * om->xf[p7O_E][p7O_MOVE];

  xEv.f32x4  = vdupq_n_f32(xc[p7C_E]);
  dp   = dpc + Q*p7C_NSCELLS - 1;
  tpdd = om->tfv + 8*Q - 1;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--) 
    {
      *dp--      = zerov;		                                              /* I */
      dcv.f32x4 = vaddq_f32(xEv.f32x4, vmulq_f32(dcv.f32x4, (*tpdd).f32x4)); tpdd--;  /* D */
      *dp--     = dcv; 
      *dp--     = xEv; 	                                                              /* M */
    }
  backward_row_finish(om, dpc, Q, dcv);
  backward_row_rescale(xc, dpc, Q, scalefactor);
}
#endif


/* backward_row_finish()
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
 * Q     = width of rows in # vectors; ox->Qf
 * dcv   = the first D vector [1 5 9 13] from caller's earlier calculation
 */

#ifdef HAVE_NEON
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
	  dcv = esl_neon_leftshift_float(dcv, zerov); /* [1 5 9 13] => [5 9 13 *]          */
	  tp  = om->tfv + 8*Q - 1;	          /* <*tp> now the [4 8 12 x] TDD quad */
	  dp  = dpc + Q*p7C_NSCELLS - 2;          /* init to point at D(i,q) vector    */
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
	  tp  = om->tfv + 8*Q - 1;	
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
  tp  = om->tfv + 7*Q - 3;	 
  dp  = dpc + (Q-1)*p7C_NSCELLS; 
  for (q = Q-1; q >= 0; q--)
    {
      (*dp).f32x4  = vaddq_f32((*dp).f32x4, vmulq_f32(dcv.f32x4, (*tp).f32x4)); tp -= 7; 
      dcv  = *(dp+1);                                                           dp -= p7C_NSCELLS;
    }
}
#endif

/* backward_row_rescale()
 * 
 * Sparse rescaling, using the scalefactor that Forward set and 
 * Backward shares.
 * 
 * xc          - ptr to specials on row (floats)
 * dpc         - ptr to vectors for row      
 * Q           - number of vectors; ox->Qf
 * scalefactor - scalefactor that Forward set for this row
 *
 * Upon return, values in current row <dpc> have been rescaled.
 */
#ifdef HAVE_NEON
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
#endif

/* posterior_decode_row()
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

#ifdef HAVE_NEON
static inline int
posterior_decode_row_neon(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc)
{
  int             Q        = ox->Qf;
  esl_neon_128f_t        *fwd       = (esl_neon_128f_t *) ox->dpf[ox->R0 + ox->R]; /* a calculated fwd row R has been popped off */
  const  esl_neon_128f_t *bck       = (esl_neon_128f_t *) ox->dpf[rowi%2];
  float         *xf        = (float *) (fwd + Q*p7C_NSCELLS);
  const  float  *xb        = (float *) (bck + Q*p7C_NSCELLS);

  esl_neon_128f_t tmp; 
     tmp.f32x4   = vdupq_n_f32(sm_thresh);
  const esl_neon_128f_t   threshv = tmp;
  float          scaleterm = xf[p7C_SCALE] / overall_sc; /* see comments above, on how rescaling affects posterior decoding equations */
  esl_neon_128f_t place;
     place.f32x4 = vdupq_n_f32(scaleterm);  
  const esl_neon_128f_t   cv = place;
  float  pnonhomology;
  esl_neon_128i_t mask;
  int    maskbits;		/* xxxx 4-bit mask for which cells 0..3 have passed threshold (if any) */
  esl_neon_128f_t pv;
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
  pnonhomology = (xf[p7C_N] * xb[p7C_N] + xf[p7C_JJ] * xb[p7C_JJ] + xf[p7C_CC] * xb[p7C_CC]) * scaleterm;
  if (pnonhomology <= 1.0f - sm_thresh)
    {
      if ((status = p7_sparsemask_StartRow(sm, rowi)) != eslOK) return status;
      for (q = Q-1; q >= 0; q--)	           // reverse, because SPARSEMASK is entirely in reversed order 
	{
	  pv.f32x4       =                     vmulq_f32(P7C_MQ(fwd, q).f32x4, P7C_MQ(bck, q).f32x4);
	  pv.f32x4       = vaddq_f32(pv.f32x4, vmulq_f32(P7C_IQ(fwd, q).f32x4, P7C_IQ(bck, q).f32x4));
	  pv.f32x4       = vaddq_f32(pv.f32x4, vmulq_f32(P7C_DQ(fwd, q).f32x4, P7C_DQ(bck, q).f32x4));
	  pv.f32x4       = vmulq_f32(pv.f32x4, cv.f32x4);           // pv is now the posterior probability of elements q,r=0..3 
	  mask.u32x4     = vcgeq_f32(pv.f32x4, threshv.f32x4);    // mask now has all 0's in elems r that failed thresh; all 1's for r that passed 
	  union { esl_neon_128i_t v; uint32_t i[4]; }s;
	  s.v = mask;
	  maskbits = ((s.i[0] & 0x1) | (s.i[1] & 0x2) | (s.i[2] & 0x4) | (s.i[3] & 0x8));
	  /* maskbits is now something like 0100: 1's indicate which cell passed. */ 
	  
	  for (r = 0; r < p7_VNF; r++) 
	    if ( maskbits & (1<<r)) 
	      if ((status = p7_sparsemask_Add(sm, q, r)) != eslOK) return status;
	}
      if ((status = p7_sparsemask_FinishRow(sm)) != eslOK) return status;
    }

#ifdef p7_DEBUGGING
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

  if (ox->pp)  save_debug_row_pp(ox, fwd, rowi);
#endif
  return eslOK;
}
#endif
/*------------------ end, inlined recursions -------------------*/





/*****************************************************************
 *  Debugging functions 
 *****************************************************************/
/* backward_row_zero()
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
backward_row_zero_neon(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox)
{
#ifdef HAVE_NEON
#ifdef p7_DEBUGGING
  int                    Q    = ox->Qf;
  esl_neon_128f_t       *dpc  = (esl_neon_128f_t *) ox->dpf[0];
  esl_neon_128f_t       *dpp  = (esl_neon_128f_t *) ox->dpf[1];
  const esl_neon_128f_t *rp   = om->rfv[x1];
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
  tp  = om->tfv;
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
  #endif /* p7_DEBUGGING */
  #endif /* HAVE_NEON */
  #ifndef HAVE_NEON
    return 0.0;
  #endif
}

/* save_debug_row_pp()
 * 
 * Debugging only. Transfer posterior decoding values from a
 * vectorized row to appropriate row of <ox->pp>, as probabilities.
 *
 * Zero wherever there's no emission:
 *   all values on row 0 
 *   all D states
 *   all E,B states
 * Zero in row positions that can't be reached in some state:
 *   all insert states on row 1
 *   J(1),C(1)
 *   N(L), J(L) on row L
 * Zero wherever there's no state:
 *   D1, IM states
 */

static void
save_debug_row_pp_neon(P7_CHECKPTMX *ox, debug_print *dpc, int i)
{
  #ifdef HAVE_NEON
  #ifdef p7_DEBUGGING
  union { esl_neon_128f_t v; float x[p7_VNF]; } u;
  int      Q  = ox->Qf;
  float  *xc  = (float *) (dpc + Q*p7C_NSCELLS);
  int     q,k,z,s;

  if (! ox->pp) return;
  
  P7R_XMX(ox->pp,i,p7R_E)  = xc[p7C_E];
  P7R_XMX(ox->pp,i,p7R_N)  = xc[p7C_N];
  P7R_XMX(ox->pp,i,p7R_J)  = xc[p7C_J];
  P7R_XMX(ox->pp,i,p7R_B)  = xc[p7C_B];
  P7R_XMX(ox->pp,i,p7R_L)  = xc[p7C_B];	/* all mass in local path */
  P7R_XMX(ox->pp,i,p7R_G)  = 0.0;	/* ... none in glocal     */
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
    {
      u.v = P7C_MQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_ML) = u.x[z]; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_DL) = u.x[z]; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_IL) = u.x[z]; }
    }
  #endif /* p7_DEBUGGING */
  #endif /* HAVE_NEON */
}
/* save_debug_row_fb()
 * 
 * Debugging only. Transfer posterior decoding values (sparse scaled,
 * prob space) from a vectorized row, to appropriate row of <ox->fwd>
 * or <ox->bck> (log space, inclusive of partial sum of scalefactors);
 * <ox->fwd> and <ox->bck> should be identical (within numerical error
 * tolerance) to a reference implementation Forward/Backward in log
 * space.
 */

static void
save_debug_row_fb_neon(P7_CHECKPTMX *ox, P7_REFMX *gx, debug_print *dpc, int i, float totscale)
{
  #ifdef HAVE_NEON
  #ifdef p7_DEBUGGING
  union { esl_neon_128f_t v; float x[p7_VNF]; } u;
  int      Q  = ox->Qf;
  float  *xc  = (float *) (dpc + Q*p7C_NSCELLS);
  int     q,k,z;

  if (! gx) return;
  
  P7R_XMX(gx,i,p7R_E)  = logf(xc[p7C_E]) + totscale;
  P7R_XMX(gx,i,p7R_N)  = logf(xc[p7C_N]) + totscale;
  P7R_XMX(gx,i,p7R_J)  = logf(xc[p7C_J]) + totscale;
  P7R_XMX(gx,i,p7R_B)  = logf(xc[p7C_B]) + totscale;
  P7R_XMX(gx,i,p7R_L)  = P7R_XMX(gx,i,p7R_B);         /* filter is local-mode. all mass assigned to local path */
  P7R_XMX(gx,i,p7R_G)  = -eslINFINITY;		      /* ... and no mass assigned to glocal path               */
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
    {
      u.v = P7C_MQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_ML) = logf(u.x[z]) + totscale; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_DL) = logf(u.x[z]) + totscale; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_IL) = logf(u.x[z]) + totscale; }
    }
  #endif /* HAVE_NEON */
  #endif /* p7_DEBUGGING */
}

/*---------------- end, debugging tools -------------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
                                         
