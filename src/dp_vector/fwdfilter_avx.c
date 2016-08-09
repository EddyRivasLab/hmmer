/* Forwards/Backwards filters
 *
 * This is code in the acceleration pipeline:
 * SSVFilter -> MSVFilter -> VitFilter -> ForwardFilter -> BackwardFilter
 *                                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 *                                                (you are here)
 *
 * Striped SIMD vector implementation [Farrar07], using checkpointing 
 * to guarantee O(M sqrt L) memory [Grice97,TarnasHughey98,Newberg08],
 * in probability space using sparse rescaling [Eddy11].
 *
 * Called in two pieces. p7_ForwardFilter() returns the Forward score
 * in nats. Caller chooses whether or not to proceed with Backward and
 * posterior decoding. If caller proceeds, p7_BackwardFilter() does
 * Backward and posterior decoding. Based on the posterior decoding
 * probabilities on each row i, it determines which cells are to be
 * added to the sparse DP mask for subsequent local/glocal
 * reprocessing.  The sparse DP mask is returned in a <P7_SPARSEMASK>
 * structure.
 * 
 * Any cell (i,k) with total posterior probability (i.e., summed over
 * M,D,I) >= sm_thresh is marked and included in the sparse mask.
 * Cells needed for glocal entry/exit delete paths (G->DDDD->Mk,
 * Mk->DDDD->E) are not marked, because the filter only uses local
 * alignment, not glocal. The default for sm_thresh is
 * <p7_SPARSIFY_THRESH>, in <p7_config.h.in>; currently set
 * to 0.01.
 *
 * ForwardFilter() and BackwardFilter() are a dependent pair, sharing
 * the same DP matrix object. They must be called sequentially,
 * because BackwardFilter() is doing posterior decoding on the fly,
 * and for this it needs both forward and backward scores; it uses
 * checkpointed Forward rows that have been stored by the preceding
 * ForwardFilter() call.
 * 
 * Prob-space implementations using sparse rescaling require multihit
 * local alignment mode for numeric range reasons. Unihit or glocal
 * will result in errors due to underflow.
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
 *    9. Notes
 *       a. On debugging and testing methods.
 *       b. On running time, in theory and in practice.
 *       c. Copyright and license information.
 */
#include "p7_config.h"

#if p7_CPU_ARCH == x86
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#endif
#include "easel.h"

#ifdef HAVE_AVX2
#include <immintrin.h>
#include "esl_avx.h"
#endif

#include "esl_vectorops.h"

#include "dp_reference/p7_refmx.h"
#include "dp_sparse/p7_sparsemx.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/fwdfilter.h"


/* vectorized DP recursions are tediously lengthy, so for some
 * semblance of clarity, they're broken out into one-page-ish
 * chunks, using static inlined functions.
 */
#ifdef HAVE_AVX2
static inline float
forward_row_avx(ESL_DSQ xi, const P7_OPROFILE *om, const __m256 *dpp, __m256 *dpc, int Q);
static inline void  backward_row_main_avx(ESL_DSQ xi, const P7_OPROFILE *om,       __m256 *dpp, __m256 *dpc, int Q, float scalefactor);
static inline void  backward_row_L_avx  (            const P7_OPROFILE *om,                    __m256 *dpc, int Q, float scalefactor);
static inline void  backward_row_finish_avx(          const P7_OPROFILE *om,                    __m256 *dpc, int Q, __m256 dcv);
static inline void  backward_row_rescale_avx(float *xc, __m256 *dpc, int Q, float scalefactor);
static inline int   posterior_decode_row_avx(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc);
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
p7_ForwardFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc)
{
//  printf("Calling ForwardFilter\n");
  #ifdef HAVE_AVX2
  int           Q_AVX     = P7_NVF_AVX(om->M);                   /* segment length; # of MDI vectors on each row      */
  __m256       *dpp_AVX   = NULL;                            /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  __m256       *dpc_AVX   = NULL;                    /* dpc points at current row         */
  const __m256  zerov_AVX = _mm256_setzero_ps();     
  float        *xc_AVX    = NULL;                            /* specials E,N,JJ,J,B,CC,C,SCALE    */
  float         totsc_AVX = 0.0f;                          /* accumulates Forward score in nats */
  int     q_AVX;      /* counter over vectors 0..Q-1                        */
   ox->Qf_AVX = Q_AVX;

  int     i;			/* counter over residues/rows 1..L                    */
  int     b;			/* counter down through checkpointed blocks, Rb+Rc..1 */
  int     w;			/* counter down through rows in a checkpointed block  */

  /* Contract checks */
  //  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  //  ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  //  ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  // SRE: disabled 21 Jul 14 because unit tests in sparse_asc_fwdback need to create
  // sparse mask w/ checkpointed F/B w/ unilocal L=0 models

  /* Make sure <ox> is allocated big enough.
   * DO NOT set any ptrs into the matrix until after this potential reallocation!
   */
  p7_checkptmx_GrowTo(ox, om->M, L);
 
  dpp_AVX =  (__m256 *) ox->dpf_AVX[ox->R0-1];    /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  xc_AVX  =  (float *) (dpp_AVX + Q_AVX*p7C_NSCELLS); /* specials E,N,JJ,J,B,CC,C,SCALE    */

  /* Set the size of the problem in <ox> now, not later
   * Debugging dumps need this information, for example
   * (Ditto for any debugging copy of the fwd mx)
   */
  ox->M  = om->M;	
  ox->L  = L;
 
#ifdef p7_DEBUGGING
  ox->dump_flags |= p7_SHOW_LOG;                     /* also sets for Backward dumps, since <ox> shared */
  if (ox->do_dumping) p7_checkptmx_DumpFBHeader(ox);
  if (ox->fwd) { ox->fwd->M = om->M; ox->fwd->L = L; ox->fwd->type = p7R_FORWARD; }
#endif

  /* Initialization of the zero row, including specials */

  for (q_AVX = 0; q_AVX < p7C_NSCELLS*Q_AVX; q_AVX++) dpp_AVX[q_AVX] = zerov_AVX;
  xc_AVX[p7C_N]     = 1.;
  xc_AVX[p7C_B]     = om->xf[p7O_N][p7O_MOVE]; 
  xc_AVX[p7C_E]     = xc_AVX[p7C_JJ] = xc_AVX[p7C_J]  = xc_AVX[p7C_CC] = xc_AVX[p7C_C]  = 0.;     
  xc_AVX[p7C_SCALE] = 1.;


#ifdef p7_DEBUGGING
  if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, 0, dpp, "f1 O"); 
  if (ox->fwd)        save_debug_row_fb(ox, ox->fwd, dpp, 0, totsc); 
#endif

  /* Phase one: the "a" region: all rows in this region are saved */
  for (i = 1; i <= ox->La; i++)
    {

      dpc_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX]; ox->R_AVX++;    /* idiomatic for "get next save/checkpoint row" */

      totsc_AVX += forward_row_avx(dsq[i], om, dpp_AVX, dpc_AVX, Q_AVX);
      dpp_AVX = dpc_AVX;   /* current row becomes prev row */

#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, dpc, "f1 O"); 
      if (ox->fwd)        save_debug_row_fb(ox, ox->fwd, dpc, i, totsc); 
#endif
    }

  /* Phase two: "b" and "c" regions: partially and fully checkpointed */
  for (b = ox->Rb + ox->Rc, w = (ox->Rb ? ox->Lb : ox->Rc+1); i <= L; i++)
    {
      /* this section, deciding whether to get a checkpointed vs. tmp_check
       * row, is why we set <fwd>, <dpp> here, rather than having the
       * inlined forward_row() set them.
       */
      if (! (--w)) { 		                   /* we're on the last row in segment: this row is saved    */
    
  dpc_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX]; ox->R_AVX++;  /* idiomatic for "get next save/checkpoint row"    */


	w = b;  			           /* next segment has this many rows, ending in a saved row */
	b--;					   /* decrement segment number counter; last segment is r=1  */
      } else{

       dpc_AVX = (__m256 *) ox->dpf_AVX[i%2];        /* idiomatic for "next tmp row", 0/1; i%2 makes sure dpp != dpc */

      }


      totsc_AVX += forward_row_avx(dsq[i], om, dpp_AVX, dpc_AVX, Q_AVX);
      dpp_AVX = dpc_AVX;

#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_checkptmx_DumpFBRow(ox, i, dpc, w ? "f1 X" : "f1 O"); 
      if (ox->fwd)        save_debug_row_fb(ox, ox->fwd, dpc, i, totsc); 
#endif

    }


  xc_AVX     = (float *) (dpc_AVX + Q_AVX*p7C_NSCELLS);

  ESL_DASSERT1( (ox->R_AVX == ox->Ra+ox->Rb+ox->Rc) );
  ESL_DASSERT1( ( (! isnan(xc_AVX[p7C_C])) && (! isinf(xc_AVX[p7C_C]))) );

  if (opt_sc) *opt_sc = totsc_AVX + logf(xc_AVX[p7C_C] * om->xf[p7O_C][p7O_MOVE]);

  return eslOK;
  #endif //HAVE_AVX2
  #ifndef HAVE_AVX2
  return eslENORESULT;  // stub so there's still a function here if we don't support AVX
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
p7_BackwardFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh)
{
//printf("Calling Backward Filter\n");

#ifdef HAVE_AVX2
    float  *xf_AVX;
  int Q_AVX = ox->Qf_AVX;
  __m256 *fwd_AVX;
  __m256 *bck_AVX;
  __m256 *dpp_AVX;
  int i_AVX, i;

  float   Tvalue, Tvalue_AVX, Tvalue_AVX_512;
  int     b, w, i2;
  int     status;

  p7_sparsemask_Reinit(sm, om->M, L);

  /* Contract checks */
  //ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  //ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  //ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
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
 
   i_AVX = L;
  ox->R_AVX--;
  fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0 + ox->R_AVX];      /* pop row for fwd[L] off the checkpointed stack */
  xf_AVX  = (float *) (fwd_AVX + Q_AVX*p7C_NSCELLS);
  Tvalue_AVX = xf_AVX[p7C_C] * om->xf[p7O_C][p7O_MOVE];  /* i.e. scaled fwd[L] val at T state = scaled overall score */
  bck_AVX = (__m256 *) ox->dpf_AVX[i_AVX%2];           /* get tmp space for bck[L]                                 */
  backward_row_L_avx(om, bck_AVX, Q_AVX, xf_AVX[p7C_SCALE]);     /* calculate bck[L] row                                     */

  if ( (status = posterior_decode_row_avx(ox, i_AVX, sm, sm_thresh, Tvalue_AVX)) != eslOK) return status;
  i_AVX--;
  dpp_AVX = bck_AVX;


  /* If there's any checkpointing, there's an L-1 row to fill now. */
  if (ox->Rb+ox->Rc > 0)
    {
  //    printf("L-1 row\n");

      /* Compute fwd[L-1] from last checkpoint, which we know is fwd[L-2] */
      dpp_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX-1];  /* fwd[L-2] values, already known        */
      fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX];    /* get free row memory from top of stack */
      forward_row_avx(dsq[i_AVX], om, dpp_AVX, fwd_AVX, Q_AVX);      /* calculate fwd[L-1]                    */


      /* Compute bck[L-1] from bck[L]. */
      xf_AVX  = (float *) (fwd_AVX + Q_AVX*p7C_NSCELLS);
      dpp_AVX = (__m256 *) ox->dpf_AVX[(i_AVX+1)%2]; 
      bck_AVX = (__m256 *) ox->dpf_AVX[i_AVX%2];             /* get space for bck[L-1]                */
      backward_row_main_avx(dsq[i_AVX+1], om, dpp_AVX, bck_AVX, Q_AVX, xf_AVX[p7C_SCALE]);

      /* And decode. */
      if ( (status = posterior_decode_row_avx(ox, i_AVX, sm, sm_thresh, Tvalue_AVX)) != eslOK) return status;
      dpp_AVX = bck_AVX;
      i_AVX--;      /* i is now L-2 if there's checkpointing; else it's L-1 */
    }

  /* Main loop for checkpointed regions (b,c) */
   for (b = 2; b <= ox->Rb+ox->Rc; b++)
    {				/* i=L-2 as we enter here, and <dpp> is on bck[L-1] */
      w = (b <= ox->Rc ? b+1 : ox->Lb);
//printf("Checkpointed region\n");

      /* We know current row i (r=R0+R-1) ends a block and is checkpointed in fwd. */
      ox->R_AVX--;
      fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX];      /* pop checkpointed forward row off "stack" */
      xf_AVX  = (float *) (fwd_AVX + Q_AVX*p7C_NSCELLS);
      Tvalue_AVX = xf_AVX[p7C_C] * om->xf[p7O_C][p7O_MOVE];  /* i.e. scaled fwd[L] val at T state = scaled overall score */
      /* Calculate bck[i]; <dpp> is already bck[i+1] */
      bck_AVX = (__m256 *) ox->dpf_AVX[i_AVX%2];      /* get available tmp memory for row     */
      backward_row_main_avx(dsq[i_AVX+1], om, dpp_AVX, bck_AVX, Q_AVX, xf_AVX[p7C_SCALE]);

      /* And decode checkpointed row i. */
      if ( (status = posterior_decode_row_avx(ox, i_AVX, sm, sm_thresh, Tvalue_AVX)) != eslOK) return status;
      
      /* The rest of the rows in the block weren't checkpointed.
       * Compute Forwards from last checkpoint ...
       */
      dpp_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX-1];       /* get last Fwd checkpoint. */
      for (i2 = i_AVX-w+1; i2 <= i_AVX-1; i2++)
  {
    fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX]; ox->R_AVX++;  /* push new forward row on "stack"     */
    forward_row_avx(dsq[i2], om, dpp_AVX, fwd_AVX, Q_AVX);
    dpp_AVX = fwd_AVX;    
  }

      /* ... and compute Backwards over the block we just calculated, while decoding. */
      dpp_AVX = bck_AVX;
      for (i2 = i_AVX-1; i2 >= i_AVX-w+1; i2--)
  {
    ox->R_AVX--;
    fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX]; /* pop just-calculated forward row i2 off "stack" */
    xf_AVX  = (float *) (fwd_AVX + Q_AVX*p7C_NSCELLS);
    bck_AVX = (__m256 *) ox->dpf_AVX[i2%2];   /* get available for calculating bck[i2]          */
    backward_row_main_avx(dsq[i2+1], om, dpp_AVX, bck_AVX, Q_AVX, xf_AVX[p7C_SCALE]);

    if ((status = posterior_decode_row_avx(ox, i2, sm, sm_thresh, Tvalue_AVX)) != eslOK) return status;
    dpp_AVX = bck_AVX;
  }
      i_AVX -= w;

    }
   /* now i=La as we leave the checkpointed regions; or i=L-1 if there was no checkpointing */
   /* The uncheckpointed "a" region */
   for (i = i_AVX; i >= 1; i--)
     {
    //  printf("uncheckpointed region\n");
     // Use regular i, not i_AVX here because it's the loop iterator variable
       ox->R_AVX--; 
       fwd_AVX = (__m256 *) ox->dpf_AVX[ox->R0+ox->R_AVX]; /* pop off calculated row fwd[i]           */
       xf_AVX  = (float *) (fwd_AVX + Q_AVX*p7C_NSCELLS);
       bck_AVX = (__m256 *) ox->dpf_AVX[i%2];        /* get open space for bck[i]               */
       backward_row_main_avx(dsq[i+1], om, dpp_AVX, bck_AVX, Q_AVX, xf_AVX[p7C_SCALE]);
 
       if ((status = posterior_decode_row_avx(ox, i, sm, sm_thresh, Tvalue_AVX)) != eslOK) return status;
       dpp_AVX = bck_AVX;

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
   fwd = (__m128 *) ox->dpf[ox->R0];
   bck = (__m128 *) ox->dpf[i%2];	       
   xN = backward_row_zero(dsq[1], om, ox); 
   if (ox->do_dumping) { p7_checkptmx_DumpFBRow(ox, 0, fwd, "f2 O"); p7_checkptmx_DumpFBRow(ox, 0, bck, "bck"); }
   if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, 0, ox->bcksc); 
   if ((status = posterior_decode_row(ox, 0, sm, sm_thresh, Tvalue)) != eslOK) return status;
   ox->bcksc += xN;
#endif

   p7_sparsemask_Finish_avx(sm);
   return eslOK;
#endif
#ifndef HAVE_AVX2
   return eslENORESULT;
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

#ifdef HAVE_AVX2
static inline float
forward_row_avx(ESL_DSQ xi, const P7_OPROFILE *om, const __m256 *dpp, __m256 *dpc, int Q)
{
  const    __m256 *rp   = om->rfv_AVX[xi];
  const    __m256 zerov = _mm256_setzero_ps();
  const    __m256 *tp   = om->tfv_AVX;
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
printf("Q= %d\n", Q);
  /* DP recursion for main states, all but the D->D path */
  for (q = 0; q < Q; q++)
    {
      /* Calculate M(i,q); hold it in tmp var <sv> */
      sv     =                _mm256_mul_ps(xBv, *tp);  tp++; /* B->Mk    */
      sv     = _mm256_add_ps(sv, _mm256_mul_ps(mpv, *tp)); tp++; /* Mk-1->Mk */
      sv     = _mm256_add_ps(sv, _mm256_mul_ps(ipv, *tp)); tp++; /* Ik-1->Mk */
      sv     = _mm256_add_ps(sv, _mm256_mul_ps(dpv, *tp)); tp++; /* Dk-1->Dk */
      sv     = _mm256_mul_ps(sv, *rp);                  rp++; /* e_Mk(x_i)*/
      xEv    = _mm256_add_ps(xEv, sv);        /* Mk->E    */

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
      sv             =                _mm256_mul_ps(mpv, *tp);  tp++;
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
  dcv            = esl_avx_leftshift_ps(dcv);  // Function naming issue: The esl_sse_(right/left)shift_ps functions 
  // describe the direction of logical shift of the vectors.  The esl_avx_leftshift functions that I wrote for
  // other filters talk about the direction of the shift instruction, which is influenced by the little-endian 
  // nature of x86.  For consistency with myself, I'm sticking to function names that match the direction of the 
  // shift instruction
  P7C_DQ(dpc, 0) = zerov;
  tp             = om->tfv_AVX + 7*Q; /* set tp to start of the DD's */
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
    {     /* Fully serialized version */
      for (j = 1; j < 8 /*was 4 */; j++)  /* What is this number?  Should it be vector length? Guess yes because of the 
      rightshift inside the loop*/
  { 
    dcv = esl_avx_leftshift_ps(dcv);
    tp  = om->tfv_AVX + 7*Q;  /* reset tp to start of the DD's */
    for (q = 0; q < Q; q++) 
      { /* note, extend dcv, not DMO(q); only adding DD paths now */
        P7C_DQ(dpc,q) = _mm256_add_ps(dcv, P7C_DQ(dpc,q)); 
        dcv           = _mm256_mul_ps(dcv, *tp);   tp++; 
      }     
  }
    } 
  else
    {     /* Slightly parallelized version, but which incurs some overhead */
      for (j = 1; j < 8; j++) /* ditto for this number, see previous j loop */
  {
    register __m256 cv = zerov; /* keeps track of whether any DD's change DMO(q) */

    dcv = esl_avx_leftshift_ps(dcv);
    tp  = om->tfv_AVX + 7*Q;  /* set tp to start of the DD's */
    for (q = 0; q < Q; q++) 
      { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
        sv            = _mm256_add_ps(dcv, P7C_DQ(dpc,q)); 
        cv            = _mm256_or_ps(cv, _mm256_cmp_ps(sv, P7C_DQ(dpc,q), 14));   // 14 = code for compare greater than
        P7C_DQ(dpc,q) = sv;                                     /* store new DMO(q) */
        dcv           = _mm256_mul_ps(dcv, *tp);   tp++;                 /* note, extend dcv, not DMO(q) */
      }     
    if (! _mm256_movemask_ps(cv)) break; /* DD's didn't change any DMO(q)? Then done, break out. */
  }
    }
  
  /* Add Dk's to xEv */
  for (q = 0; q < Q; q++) xEv = _mm256_add_ps(P7C_DQ(dpc,q), xEv);

  /* Specials, in order: E N JJ J B CC C */
  /* xc, xp are pointers into vectors provided as inputs, so should be no problem with running multiple
  ISAs at same time */

  esl_avx_hsum_ps(xEv, &xc[p7C_E]);  
  xc[p7C_N]  =                                       xp[p7C_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7C_JJ] =                                       xp[p7C_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7C_J]  = xc[p7C_JJ]                          + xc[p7C_E] * om->xf[p7O_E][p7O_LOOP];
  xc[p7C_B]  = xc[p7C_N] * om->xf[p7O_N][p7O_MOVE] + xc[p7C_J] * om->xf[p7O_J][p7O_MOVE];
  xc[p7C_CC] =                                       xp[p7C_C] * om->xf[p7O_C][p7O_LOOP];
  xc[p7C_C]  = xc[p7C_CC]                          + xc[p7C_E] * om->xf[p7O_E][p7O_MOVE];

  /* Sparse rescaling. xE above threshold? Then trigger a rescaling event.            */
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

#ifdef HAVE_AVX2
static inline void
backward_row_main_avx(ESL_DSQ xi, const P7_OPROFILE *om, __m256 *dpp, __m256 *dpc, int Q, float scalefactor)
{
  const __m256 *rp       = om->rfv_AVX[xi];           /* emission scores on row i+1, for bck; xi = dsq[i+1]  */
  float       * const xc = (float *) (dpc + Q * p7C_NSCELLS); /* E N JJ J B CC C SCALE */
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
  tp  = om->tfv_AVX;    /* on first transition vector */
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
  mpv  = esl_avx_rightshift_ps(*dpp      ); /* [1 5 9 13] -> [5 9 13 x], M(i+1,k+1) * e(M_k+1, x_{i+1}) */
  tmmv = esl_avx_rightshift_ps(om->tfv_AVX[1]);
  timv = esl_avx_rightshift_ps(om->tfv_AVX[2]);
  tdmv = esl_avx_rightshift_ps(om->tfv_AVX[3]);
  xEv  = _mm256_set1_ps(xc[p7C_E]);
  tp   = om->tfv_AVX + 7*Q - 1;
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
#endif

/* backward_row_L()
 * 
 * Backward calculation for row L; 
 * a special case because the matrix <ox> has no 'previous' row L+1.
 * Otherwise identical to backward_row_main().
 */
#ifdef HAVE_AVX2
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
  tpdd = om->tfv_AVX + 8*Q - 1;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--) 
    {
      *dp--       = zerov;                                  /* I */
      *dp-- = dcv = _mm256_add_ps(xEv, _mm256_mul_ps(dcv, *tpdd)); tpdd--;  /* D */
      *dp--       = xEv;                                        /* M */
    }
  backward_row_finish_avx(om, dpc, Q, dcv);
  backward_row_rescale_avx(xc, dpc, Q, scalefactor);
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

#ifdef HAVE_AVX2
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
    dcv = esl_avx_rightshift_ps(dcv); /* [1 5 9 13] => [5 9 13 *]          */
    tp  = om->tfv_AVX + 8*Q - 1;            /* <*tp> now the [4 8 12 x] TDD quad */
    dp  = dpc + Q*p7C_NSCELLS - 2;          /* init to point at D(i,q) vector    */
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
      for (j = 1; j < 8; j++)
  {
    dcv = esl_avx_rightshift_ps(dcv);
    tp  = om->tfv_AVX + 8*Q - 1;  
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
  tp  = om->tfv_AVX + 7*Q - 3;   
  dp  = dpc + (Q-1)*p7C_NSCELLS; 
  for (q = Q-1; q >= 0; q--)
    {
      *dp  = _mm256_add_ps(*dp, _mm256_mul_ps(dcv, *tp)); tp -= 7; 
      dcv  = *(dp+1);                               dp -= p7C_NSCELLS;
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

#ifdef HAVE_AVX2
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

#ifdef HAVE_AVX2
static inline int
posterior_decode_row_avx(P7_CHECKPTMX *ox, int rowi, P7_SPARSEMASK *sm, float sm_thresh, float overall_sc)
{
  int             Q        = ox->Qf_AVX;
  __m256        *fwd       = (__m256 *) ox->dpf_AVX[ox->R0 + ox->R_AVX]; /* a calculated fwd row R has been popped off */
  const  __m256 *bck       = (__m256 *) ox->dpf_AVX[rowi%2];
  float         *xf        = (float *) (fwd + Q*p7C_NSCELLS);
  const  float  *xb        = (float *) (bck + Q*p7C_NSCELLS);
  const __m256   threshv   = _mm256_set1_ps(sm_thresh);
  float          scaleterm = xf[p7C_SCALE] / overall_sc; /* see comments above, on how rescaling affects posterior decoding equations */
  const __m256   cv        = _mm256_set1_ps(scaleterm);
  float  pnonhomology;
  __m256 mask;
  int    maskbits;    /* xxxx 4-bit mask for which cells 0..3 have passed threshold (if any) */
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
      if ((status = p7_sparsemask_StartRow_avx(sm, rowi)) != eslOK) return status;
      for (q = Q-1; q >= 0; q--)             // reverse, because SPARSEMASK is entirely in reversed order 
  {
    pv       =                _mm256_mul_ps(P7C_MQ(fwd, q), P7C_MQ(bck, q));
    pv       = _mm256_add_ps(pv, _mm256_mul_ps(P7C_IQ(fwd, q), P7C_IQ(bck, q)));
    pv       = _mm256_add_ps(pv, _mm256_mul_ps(P7C_DQ(fwd, q), P7C_DQ(bck, q)));
    pv       = _mm256_mul_ps(pv, cv);           // pv is now the posterior probability of elements q,r=0..3 
    mask     = _mm256_cmp_ps(pv, threshv, 13);  // 13 is magic number for >= 
     // mask now has all 0's in elems r that failed thresh; all 1's for r that passed 
    maskbits = _mm256_movemask_ps(mask);    // maskbits is now something like 0100: 1's indicate which cell passed. 
  
    for (r = 0; r < p7_VNF_AVX; r++) 
      if ( maskbits & (1<<r)) 
         if ((status = p7_sparsemask_Add_avx(sm, q, r)) != eslOK) return status;
  }
      if ((status = p7_sparsemask_FinishRow_avx(sm)) != eslOK) return status;
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
      P7C_MQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7C_MQ(fwd, q), P7C_MQ(bck, q)));
      P7C_DQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7C_DQ(fwd, q), P7C_DQ(bck, q)));
      P7C_IQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7C_IQ(fwd, q), P7C_IQ(bck, q)));
    }

  if (ox->pp)  save_debug_row_pp(ox, fwd, rowi);
#endif
  return eslOK;
}
#endif

/*------------------ end, inlined recursions -------------------*/
/*------------------ Debugging Functions ---------------------*/
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
backward_row_zero_avx(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox)
{
#ifdef HAVE_AVX2 
#ifdef p7_DBUGGING
  int          Q     = ox->Qf_AVX;
  __m256       *dpc  = (__m256 *) ox->dpf_AVX[0];
  __m256       *dpp  = (__m256 *) ox->dpf_AVX[1];
  const __m256 *rp   = om->rfv_AVX[x1];
  const __m256 zerov = _mm256_setzero_ps();
  float        *xc   = (float *) (dpc + Q * p7C_NSCELLS); /* special states on current row i  */
  float        *xp   = (float *) (dpp + Q * p7C_NSCELLS); /* special states on "previous" row i+1 */
  __m256       *dp;
  __m256       *tp;
  __m256        xBv  = zerov;
  int           q;

  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. */
  dp  = dpp;
  tp  = om->tfv;
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
 #endif //p7_DEBUGGING 
#endif // HAVE_AVX2
#ifndef HAVE_AVX2
  return 0.0
#endif
}

static void
save_debug_row_pp_avx(P7_CHECKPTMX *ox, debug_print *dpc, int i)
{
#ifdef HAVE_AVX2  
 #ifdef p7_DEBUGGING 
  union { __m256 v; float x[p7_VNF_AVX]; } u;
  int      Q  = ox->Qf_AVX;
  float  *xc  = (float *) (dpc + Q*p7C_NSCELLS);
  int     q,k,z,s;

  if (! ox->pp) return;
  
  P7R_XMX(ox->pp,i,p7R_E)  = xc[p7C_E];
  P7R_XMX(ox->pp,i,p7R_N)  = xc[p7C_N];
  P7R_XMX(ox->pp,i,p7R_J)  = xc[p7C_J];
  P7R_XMX(ox->pp,i,p7R_B)  = xc[p7C_B];
  P7R_XMX(ox->pp,i,p7R_L)  = xc[p7C_B]; /* all mass in local path */
  P7R_XMX(ox->pp,i,p7R_G)  = 0.0; /* ... none in glocal     */
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
      u.v = P7C_MQ(dpc, q); for (z = 0; z < p7_VNF_AVX; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_ML) = u.x[z]; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < p7_VNF_AVX; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_DL) = u.x[z]; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < p7_VNF_AVX; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(ox->pp,i,k,p7R_IL) = u.x[z]; }
    }
  #endif //p7_DEBUGGING  
 #endif   // HAVE_AVX2

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
save_debug_row_fb_avx(P7_CHECKPTMX *ox, P7_REFMX *gx, debug_print *dpc, int i, float totscale)
{
 #ifdef HAVE_AVX2 
 #ifdef p7_DEBUGGING 
  union { __m256 v; float x[p7_VNF_AVX]; } u;
  int      Q  = ox->Qf_AVX;
  float  *xc  = (float *) (dpc + Q*p7C_NSCELLS);
  int     q,k,z;

  if (! gx) return;
  
  P7R_XMX(gx,i,p7R_E)  = logf(xc[p7C_E]) + totscale;
  P7R_XMX(gx,i,p7R_N)  = logf(xc[p7C_N]) + totscale;
  P7R_XMX(gx,i,p7R_J)  = logf(xc[p7C_J]) + totscale;
  P7R_XMX(gx,i,p7R_B)  = logf(xc[p7C_B]) + totscale;
  P7R_XMX(gx,i,p7R_L)  = P7R_XMX(gx,i,p7R_B);         /* filter is local-mode. all mass assigned to local path */
  P7R_XMX(gx,i,p7R_G)  = -eslINFINITY;          /* ... and no mass assigned to glocal path               */
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
      u.v = P7C_MQ(dpc, q); for (z = 0; z < p7_VNF_AVX; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_ML) = logf(u.x[z]) + totscale; }
      u.v = P7C_DQ(dpc, q); for (z = 0; z < p7_VNF_AVX; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_DL) = logf(u.x[z]) + totscale; }
      u.v = P7C_IQ(dpc, q); for (z = 0; z < p7_VNF_AVX; z++) { k = q+Q*z+1; if (k <= ox->M) P7R_MX(gx,i,k,p7R_IL) = logf(u.x[z]) + totscale; }
    }
 #endif   //p7_DEBUGGING
 #endif // HAVE_AVX2   
}
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/




                                          
