/* Forwards/Backwards filters.
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
 * in nats. Caller then chooses whether or not to proceed with Backward 
 * and posterior decoding. p7_BackwardFilter() then does Backward and 
 * posterior decoding, and based on the posterior decoding probabilities 
 * on each row i, it calculates a band ka..kb in which nonnegligible 
 * alignment mass appears to lie, returning these bands in a <P7_GBANDS> 
 * structure.
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
 *    4. Benchmark driver.
 *    5. Unit tests.
 *    6. Test driver.
 *    7. Example.
 *    8. Notes
 *       a. On debugging and testing methods.
 *       b. On running time, in theory and in practice.
 *       c. Copyright and license information.
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


/* vectorized DP recursions are tediously lengthy, so for some
 * semblance of clarity, they're broken out into one-page-ish
 * chunks, using static inlined functions.
 */
static inline float forward_row      (ESL_DSQ xi, const P7_OPROFILE *om, const __m128 *dpp, __m128 *dpc, int Q);
static inline void  backward_row_main(ESL_DSQ xi, const P7_OPROFILE *om,       __m128 *dpp, __m128 *dpc, int Q, float scalefactor);
static inline void  backward_row_L   (            const P7_OPROFILE *om,                    __m128 *dpc, int Q, float scalefactor);
static inline void  backward_row_finish(          const P7_OPROFILE *om,                    __m128 *dpc, int Q, __m128 dcv);
static inline void  backward_row_rescale(float *xc, __m128 *dpc, int Q, float scalefactor);
static inline void  posterior_decode_row(P7_FILTERMX *ox, int rowi, P7_GBANDS *bnd, float overall_sc);

#ifdef p7_DEBUGGING
static inline float backward_row_zero(ESL_DSQ x1, const P7_OPROFILE *om, P7_FILTERMX *ox);
static        void  save_debug_row_pp(P7_FILTERMX *ox,             __m128 *dpc, int i);
static        void  save_debug_row_fb(P7_FILTERMX *ox, P7_GMX *gx, __m128 *dpc, int i, float totscale);
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
  int           Q     = P7F_NQF(om->M);                  /* segment length; # of MDI vectors on each row      */
  __m128       *dpp   = (__m128 *) ox->dpf[ox->R0-1];    /* dpp=prev row. start on dpf[2]; rows 0,1=Backwards */
  __m128       *dpc   = NULL;		                 /* dpc points at current row         */
  float        *xc    = (float *) (dpp + Q*p7F_NSCELLS); /* specials E,N,JJ,J,B,CC,C,SCALE    */
  float         totsc = 0.0f;	                         /* accumulates Forward score in nats */
  const __m128  zerov = _mm_setzero_ps();		 
  int     q;			/* counter over vectors 0..Q-1                        */
  int     i;			/* counter over residues/rows 1..L                    */
  int     b;			/* counter down through checkpointed blocks, Rb+Rc..1 */
  int     w;			/* counter down through rows in a checkpointed block  */

  /* Set the size of the problem in <ox> now, not later
   * Debugging dumps need this information, for example
   */
  ox->M  = om->M;	
  ox->L  = L;
  ox->Qf = Q;
#ifdef p7_DEBUGGING
  ox->dump_flags |= p7_SHOW_LOG;                     /* also sets for Backward dumps, since <ox> shared */
  if (ox->do_dumping) p7_filtermx_DumpFBHeader(ox);
#endif

  /* Initialization of the zero row, including specials */
  for (q = 0; q < p7F_NSCELLS*Q; q++) dpp[q] = zerov;
  xc[p7F_N]     = 1.;
  xc[p7F_B]     = om->xf[p7O_N][p7O_MOVE]; 
  xc[p7F_E]     = xc[p7F_JJ] = xc[p7F_J]  = xc[p7F_CC] = xc[p7F_C]  = 0.;			
  xc[p7F_SCALE] = 1.;			   
#ifdef p7_DEBUGGING
  if (ox->do_dumping) p7_filtermx_DumpFBRow(ox, 0, dpp, "f1 O"); 
  if (ox->fwd)        save_debug_row_fb(ox, ox->fwd, dpp, 0, totsc); 
#endif

  /* Phase one: the "a" region: all rows in this region are saved */
  for (i = 1; i <= ox->La; i++)
    {
      dpc = (__m128 *) ox->dpf[ox->R0+ox->R]; ox->R++;    /* idiomatic for "get next save/checkpoint row" */
      totsc += forward_row(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;	    	                          /* current row becomes prev row */
#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_filtermx_DumpFBRow(ox, i, dpc, "f1 O"); 
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
	dpc = (__m128 *) ox->dpf[ox->R0+ox->R]; ox->R++;  /* idiomatic for "get next save/checkpoint row"    */
	w = b;  			           /* next segment has this many rows, ending in a saved row */
	b--;					   /* decrement segment number counter; last segment is r=1  */
      } else dpc = (__m128 *) ox->dpf[i%2];        /* idiomatic for "next tmp row", 0/1; i%2 makes sure dpp != dpc */
      
      totsc += forward_row(dsq[i], om, dpp, dpc, Q);
      dpp = dpc;
#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_filtermx_DumpFBRow(ox, i, dpc, w ? "f1 X" : "f1 O"); 
      if (ox->fwd)        save_debug_row_fb(ox, ox->fwd, dpc, i, totsc); 
#endif
    }

  xc     = (float *) (dpc + Q*p7F_NSCELLS);

  ESL_DASSERT1( (ox->R == ox->Ra+ox->Rb+ox->Rc) );
  ESL_DASSERT1( ( (! isnan(xc[p7F_C])) && (! isinf(xc[p7F_C]))) );
#ifdef p7_DEBUGGING
  if (ox->fwd) { ox->fwd->M = om->M; ox->fwd->L = L; }
#endif

  if (opt_sc) *opt_sc = totsc + logf(xc[p7F_C] * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


/* Function:  p7_BackwardFilter()
 * Synopsis:  Checkpointed striped vector Backward calculation, producing bands.
 *
 * Purpose:   Given a target sequence <dsq> of length <L>, a query model
 *            <om>, and a DP matrix <ox> resulting from a successful
 *            call to <p7_ForwardFilter()>. Calculate the Backward and
 *            posterior decoding algorithms. On each row <i=1..L>, use
 *            posterior decoding to determine a band <ka..kb> (or no
 *            band at all) in which significant posterior alignment
 *            probability falls. The bands are stored in the <bnd> structure,
 *            which the caller allocates (or reuses) and provides.
 *            
 *            Currently the rule for determining 'significant'
 *            posterior alignment probability is hardcoded. If the
 *            probability of nonhomology (N/C/J emissions) is 
 *            0.9 (<p7_BANDS_THRESH1>) or more, the row is skipped
 *            altogether. If the sum of M+I emission in any k,i cell
 *            is 0.02 (<p7_BANDS_THRESH2>) or more, the cell is
 *            considered to have 'significant' posterior probability.
 *            The minimum and maximum significant k positions on 
 *            each row determine the band <ka..kb>.
 *            
 * Args:      dsq    - digital target sequence, 1..L
 *            L      - length of dsq, residues
 *            om     - optimized profile (multihit local)
 *            ox     - checkpointed DP matrix, ForwardFilter already run  
 *            bnd    - allocated P7_GBANDS structure to hold posterior bands
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_BackwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, P7_GBANDS *bnd)
{
  int Q = ox->Qf;
  __m128 *fwd;
  __m128 *bck;
  __m128 *dpp;
  float  *xf;
  float   Tvalue;
  int     i, b, w, i2;
  
  /* Row L is a special case for Backwards; so in checkpointed Backward,
   * we have to special-case the last block, rows L-1 and L
   */
  i = L;
  ox->R--;
  fwd = (__m128 *) ox->dpf[ox->R0 + ox->R];      /* pop row for fwd[L] off the checkpointed stack */
  xf  = (float *) (fwd + Q*p7F_NSCELLS);
  Tvalue = xf[p7F_C] * om->xf[p7O_C][p7O_MOVE];  /* i.e. scaled fwd[L] val at T state = scaled overall score */
  bck = (__m128 *) ox->dpf[i%2];	         /* get tmp space for bck[L]                                 */
  backward_row_L(om, bck, Q, xf[p7F_SCALE]);     /* calculate bck[L] row                                     */
#ifdef p7_DEBUGGING
  ox->bcksc = logf(xf[p7F_SCALE]);
  if (ox->do_dumping) { p7_filtermx_DumpFBRow(ox, L, fwd, "f2 O"); if (ox->do_dumping) p7_filtermx_DumpFBRow(ox, L, bck, "bck");  }
  if (ox->bck)          save_debug_row_fb(ox, ox->bck, bck, L, ox->bcksc); 
#endif
  posterior_decode_row(ox, i, bnd, Tvalue);
  i--;
  dpp = bck;


  /* If there's any checkpointing, there's an L-1 row to fill now. */
  if (ox->Rb+ox->Rc > 0)
    {
      /* Compute fwd[L-1] from last checkpoint, which we know is fwd[L-2] */
      dpp = (__m128 *) ox->dpf[ox->R0+ox->R-1];  /* fwd[L-2] values, already known        */
      fwd = (__m128 *) ox->dpf[ox->R0+ox->R];    /* get free row memory from top of stack */
      forward_row(dsq[i], om, dpp, fwd, Q);      /* calculate fwd[L-1]                    */
#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_filtermx_DumpFBRow(ox, i, fwd, "f2 X");
#endif

      /* Compute bck[L-1] from bck[L]. */
      xf  = (float *) (fwd + Q*p7F_NSCELLS);
      dpp = (__m128 *) ox->dpf[(i+1)%2]; 
      bck = (__m128 *) ox->dpf[i%2];             /* get space for bck[L-1]                */
      backward_row_main(dsq[i+1], om, dpp, bck, Q, xf[p7F_SCALE]);
#ifdef p7_DEBUGGING
      ox->bcksc += logf(xf[p7F_SCALE]);
      if (ox->do_dumping) p7_filtermx_DumpFBRow(ox, i, bck, "bck");
      if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i, ox->bcksc); 
#endif
      /* And decode. */
      posterior_decode_row(ox, i, bnd, Tvalue);
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
      backward_row_main(dsq[i+1], om, dpp, bck, Q, xf[p7F_SCALE]);
#ifdef p7_DEBUGGING
      ox->bcksc += logf(xf[p7F_SCALE]);
      if (ox->do_dumping) { p7_filtermx_DumpFBRow(ox, i, fwd, "f2 O");	p7_filtermx_DumpFBRow(ox, i, bck, "bck"); }
      if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i, ox->bcksc); 
#endif
      /* And decode checkpointed row i. */
      posterior_decode_row(ox, i, bnd, Tvalue);
      
      /* The rest of the rows in the block weren't checkpointed.
       * Compute Forwards from last checkpoint ...
       */
      dpp = (__m128 *) ox->dpf[ox->R0+ox->R-1];       /* get last Fwd checkpoint. */
      for (i2 = i-w+1; i2 <= i-1; i2++)
	{
	  fwd = (__m128 *) ox->dpf[ox->R0+ox->R]; ox->R++;  /* push new forward row on "stack"     */
	  forward_row(dsq[i2], om, dpp, fwd, Q);
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
	  backward_row_main(dsq[i2+1], om, dpp, bck, Q, xf[p7F_SCALE]);
#ifdef p7_DEBUGGING
	  ox->bcksc += logf(xf[p7F_SCALE]);
	  if (ox->do_dumping) { p7_filtermx_DumpFBRow(ox, i2, fwd, "f2 X"); p7_filtermx_DumpFBRow(ox, i2, bck, "bck"); }
	  if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i2, ox->bcksc); 
#endif
	  posterior_decode_row(ox, i2, bnd, Tvalue);
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
       backward_row_main(dsq[i+1], om, dpp, bck, Q, xf[p7F_SCALE]);
#ifdef p7_DEBUGGING
       ox->bcksc += logf(xf[p7F_SCALE]);
       if (ox->do_dumping) { p7_filtermx_DumpFBRow(ox, i, fwd, "f2 O"); p7_filtermx_DumpFBRow(ox, i, bck, "bck"); }
       if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, i, ox->bcksc); 
#endif
       posterior_decode_row(ox, i, bnd, Tvalue);
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
   fwd = (__m128 *) ox->dpf[ox->R0];
   bck = (__m128 *) ox->dpf[i%2];	       
   xN = backward_row_zero(dsq[1], om, ox); 
   if (ox->do_dumping) { p7_filtermx_DumpFBRow(ox, 0, fwd, "f2 O"); p7_filtermx_DumpFBRow(ox, 0, bck, "bck"); }
   if (ox->bck)        save_debug_row_fb(ox, ox->bck, bck, 0, ox->bcksc); 
   posterior_decode_row(ox, 0, bnd, Tvalue);
   ox->bcksc += xN;

   if (ox->bck) { ox->bck->M = om->M; ox->bck->L = L; }
   if (ox->pp)  { ox->pp->M  = om->M; ox->pp->L  = L; }
#endif

   bnd->L = L;
   bnd->M = om->M;
   p7_gbands_Reverse(bnd);
   return eslOK;
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
static inline float
forward_row(ESL_DSQ xi, const P7_OPROFILE *om, const __m128 *dpp, __m128 *dpc, int Q)
{
  const    __m128 *rp   = om->rfv[xi];
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
  
  /* Add Dk's to xEv */
  for (q = 0; q < Q; q++) xEv = _mm_add_ps(P7F_DQ(dpc,q), xEv);

  /* Specials, in order: E N JJ J B CC C */
  esl_sse_hsum_ps(xEv, &xc[p7F_E]);  
  xc[p7F_N]  =                                       xp[p7F_N] * om->xf[p7O_N][p7O_LOOP];
  xc[p7F_JJ] =                                       xp[p7F_J] * om->xf[p7O_J][p7O_LOOP];
  xc[p7F_J]  = xc[p7F_JJ]                          + xc[p7F_E] * om->xf[p7O_E][p7O_LOOP];
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
static inline void
backward_row_main(ESL_DSQ xi, const P7_OPROFILE *om, __m128 *dpp, __m128 *dpc, int Q, float scalefactor)
{
  const __m128 *rp       = om->rfv[xi];			      /* emission scores on row i+1, for bck; xi = dsq[i+1]  */
  float       * const xc = (float *) (dpc + Q * p7F_NSCELLS); /* E N JJ J B CC C SCALE */
  const float * const xp = (float *) (dpp + Q * p7F_NSCELLS);
  const __m128 *tp, *tpdd;
  const __m128  zerov = _mm_set1_ps(0.0f);
  __m128        xBv;
  __m128       *dp;
  __m128        xEv;
  __m128        dcv, mcv, ipv, mpv;
  int           q;
  __m128 tmmv, timv, tdmv;				   /* copies of transition prob quads; a leftshift is needed as boundary cond */
 
  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. 
   * This invalidates <dpp> as a backwards row; its values are now
   * intermediates in the calculation of the current <dpc> row.
   */
  dp  = dpp;
  xBv = zerov;
  tp  = om->tfv;		/* on first transition vector */
  for (q = 0; q < Q; q++)
    {
      *dp = _mm_mul_ps(*dp, *rp); rp++;
      xBv = _mm_add_ps(xBv, _mm_mul_ps(*dp, *tp)); dp+= p7F_NSCELLS; tp += 7;
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
  backward_row_finish(om, dpc, Q, dcv);
  backward_row_rescale(xc, dpc, Q, scalefactor);
}


/* backward_row_L()
 * 
 * Backward calculation for row L; 
 * a special case because the matrix <ox> has no 'previous' row L+1.
 * Otherwise identical to backward_row_main().
 */
static inline void
backward_row_L(const P7_OPROFILE *om,  __m128 *dpc, int Q, float scalefactor)
{
  const __m128  zerov = _mm_setzero_ps();
  float        *xc    = (float *) (dpc + Q * p7F_NSCELLS);
  const __m128 *tpdd;
  __m128       *dp;
  __m128       xEv, dcv;
  int          q;

  /* Backwards from T <- C,CC <- E;  all other specials unreachable, impossible on row L.
   * specials are stored in order E N JJ J B CC C.  
   */
  xc[p7F_C] = xc[p7F_CC] = om->xf[p7O_C][p7O_MOVE];
  xc[p7F_B] = xc[p7F_J] = xc[p7F_JJ] = xc[p7F_N] = 0.0;
  xc[p7F_E] = xc[p7F_C] * om->xf[p7O_E][p7O_MOVE];

  xEv  = _mm_set1_ps(xc[p7F_E]);
  dp   = dpc + Q*p7F_NSCELLS - 1;
  tpdd = om->tfv + 8*Q - 1;
  dcv  = zerov;
  for (q = Q-1; q >= 0; q--) 
    {
      *dp--       = zerov;		                              /* I */
      *dp-- = dcv = _mm_add_ps(xEv, _mm_mul_ps(dcv, *tpdd)); tpdd--;  /* D */
      *dp--       = xEv; 	                                      /* M */
    }
  backward_row_finish(om, dpc, Q, dcv);
  backward_row_rescale(xc, dpc, Q, scalefactor);
}



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
static inline void
backward_row_finish(const P7_OPROFILE *om, __m128 *dpc, int Q, __m128 dcv)
{
  const __m128 zerov = _mm_setzero_ps();
  const __m128 *tp;
  __m128       *dp;
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
	  dcv = esl_sse_leftshift_ps(dcv, zerov); /* [1 5 9 13] => [5 9 13 *]          */
	  tp  = om->tfv + 8*Q - 1;	          /* <*tp> now the [4 8 12 x] TDD quad */
	  dp  = dpc + Q*p7F_NSCELLS - 2;          /* init to point at D(i,q) vector    */
	  for (q = Q-1; q >= 0; q--)
	    {
	      dcv = _mm_mul_ps(dcv, *tp); tp--;
	      *dp = _mm_add_ps(*dp, dcv); dp -= p7F_NSCELLS;
	    }
	}
    }
  else
    { /* With check for early convergence */
      __m128 sv;
      __m128 cv;  /* keeps track of whether any DD addition changes DQ(q) value */
      for (j = 1; j < 4; j++)
	{
	  dcv = esl_sse_leftshift_ps(dcv, zerov);
	  tp  = om->tfv + 8*Q - 1;	
	  dp  = dpc + Q*p7F_NSCELLS - 2;
	  cv  = zerov;
	  for (q = Q-1; q >= 0; q--)
	    { /* using cmpgt below tests if DD changed any DMO(q) without conditional branch (i.e. no if) */
	      dcv  = _mm_mul_ps(dcv, *tp); tp--;
	      sv   = _mm_add_ps(*dp, dcv);
	      cv   = _mm_or_ps(cv, _mm_cmpgt_ps(sv, *dp)); /* if DD path changed DQ(dpc,q), cv bits know it now */
	      *dp  = sv; 
	      dp  -= p7F_NSCELLS;
	    }
	  if (! _mm_movemask_ps(cv)) break; /* if no DD path changed DQ(q) in this segment, then done, no more segments needed */
	}
    }

  /* Finally, M->D path contribution
   * these couldn't be added to M until we'd finished calculating D values on row.
   */
  dcv = esl_sse_leftshift_ps(P7F_DQ(dpc, 0), zerov);
  tp  = om->tfv + 7*Q - 3;	 
  dp  = dpc + (Q-1)*p7F_NSCELLS; 
  for (q = Q-1; q >= 0; q--)
    {
      *dp  = _mm_add_ps(*dp, _mm_mul_ps(dcv, *tp)); tp -= 7; 
      dcv  = *(dp+1);                               dp -= p7F_NSCELLS;
    }
}

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
static inline void
backward_row_rescale(float *xc, __m128 *dpc, int Q, float scalefactor)
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


/* posterior_decode_row()
 *
 * In production code, we don't have to save results of posterior
 * decoding; we immediately use them to calculate bands on the row.
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
 */
#define p7_BANDS_THRESH1  0.9
#define p7_BANDS_THRESH2  0.02
static inline void
posterior_decode_row(P7_FILTERMX *ox, int rowi, P7_GBANDS *bnd, float overall_sc)
{
  int             Q  = ox->Qf;
  __m128        *fwd = (__m128 *) ox->dpf[ox->R0 + ox->R]; /* a calculated fwd row R has been popped off */
  const  __m128 *bck = (__m128 *) ox->dpf[rowi%2];
  float         *xf  = (float *) (fwd + Q*p7F_NSCELLS);
  const  float  *xb  = (float *) (bck + Q*p7F_NSCELLS);
  __m128 threshv     = _mm_set1_ps(p7_BANDS_THRESH2);
  __m128 onev        = _mm_set1_ps(1.0);
  __m128 kv          = _mm_set_ps((float) 3*Q, (float) 2*Q, (float) Q, 0.0f); /* iv is going to contain current i coord in each vector cell */
  __m128 minkv       = _mm_set1_ps( (float) Q*4.0+1);
  __m128 maxkv       = _mm_set1_ps( 0.0f);
  float  scaleterm   = xf[p7F_SCALE] / overall_sc; /* see comments above, on how rescaling affects posterior decoding equations */
  float  pnonhomology;
  __m128 mask;
  __m128 pv;
  __m128 cv;
  int    q;
  float  ktmp;
  int    ka, kb;
  
#ifdef p7_DEBUGGING
  /* in debugging code, we store the pp's in the fwd row's space, and defer the pnonhomology 
   * threshold test until after we've calculated and saved the entire row 
   */
  xf[p7F_N]               = (rowi == 0 ? 0.0f : xf[p7F_N]  * xb[p7F_N]  * scaleterm);
  xf[p7F_J] = xf[p7F_JJ]  = xf[p7F_JJ] * xb[p7F_JJ] * scaleterm;
  xf[p7F_C] = xf[p7F_CC]  = xf[p7F_CC] * xb[p7F_CC] * scaleterm;
  xf[p7F_E] = xf[p7F_B] = 0.0f;
  pnonhomology = xf[p7F_N] + xf[p7F_J] + xf[p7F_C];
#else
  /* in production code, we don't need to store, and we may immediately threshold on pnonhomology */
  pnonhomology = (xf[p7F_N] * xb[p7F_N] + xf[p7F_JJ] * xb[p7F_JJ] + xf[p7F_CC] * xb[p7F_CC]) * scaleterm;
  if (pnonhomology >= p7_BANDS_THRESH1) return;
#endif
  
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
  cv = _mm_set1_ps(scaleterm);
  for (q = 0; q < Q; q++)
    {
      pv   =                _mm_mul_ps(P7F_MQ(fwd, q), P7F_MQ(bck, q));
      pv   = _mm_add_ps(pv, _mm_mul_ps(P7F_IQ(fwd, q), P7F_IQ(bck, q)));
      pv   = _mm_add_ps(pv, _mm_mul_ps(P7F_DQ(fwd, q), P7F_DQ(bck, q)));
      pv   = _mm_mul_ps(pv, cv);
      mask = _mm_cmpge_ps(pv, threshv);

#ifdef p7_DEBUGGING
      P7F_MQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7F_MQ(fwd, q), P7F_MQ(bck, q)));
      P7F_DQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7F_DQ(fwd, q), P7F_DQ(bck, q)));
      P7F_IQ(fwd, q) = _mm_mul_ps(cv, _mm_mul_ps(P7F_IQ(fwd, q), P7F_IQ(bck, q)));
#endif

      kv    = _mm_add_ps(kv, onev);
      minkv = _mm_min_ps( minkv, esl_sse_select_ps(minkv, kv, mask));
      maxkv = _mm_max_ps( maxkv, esl_sse_select_ps(maxkv, kv, mask));
    }
  esl_sse_hmax_ps(maxkv, &ktmp); kb = (int) ktmp;
  esl_sse_hmin_ps(minkv, &ktmp); ka = (int) ktmp;

#ifdef p7_DEBUGGING
  if (ox->pp)  save_debug_row_pp(ox, fwd, rowi);
  if (pnonhomology >= p7_BANDS_THRESH1) return;
#endif

  if (kb) p7_gbands_Prepend(bnd, rowi, ka, kb);
}
/*------------------ end, inlined recursions -------------------*/





/*****************************************************************
 * 3. Debugging tools
 *****************************************************************/
#ifdef p7_DEBUGGING

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
backward_row_zero(ESL_DSQ x1, const P7_OPROFILE *om, P7_FILTERMX *ox)
{
  int          Q     = ox->Qf;
  __m128       *dpc  = (__m128 *) ox->dpf[0];
  __m128       *dpp  = (__m128 *) ox->dpf[1];
  const __m128 *rp   = om->rfv[x1];
  const __m128 zerov = _mm_setzero_ps();
  float        *xc   = (float *) (dpc + Q * p7F_NSCELLS); /* special states on current row i  */
  float        *xp   = (float *) (dpp + Q * p7F_NSCELLS); /* special states on "previous" row i+1 */
  __m128       *dp;
  __m128       *tp;
  __m128        xBv  = zerov;
  int           q;

  /* On "previous" row i+1: include emission prob, and sum to get xBv, xB. */
  dp  = dpp;
  tp  = om->tfv;
  for (q = 0; q < Q; q++)
    {
      *dp = _mm_mul_ps(*dp, *rp); rp++;
      xBv = _mm_add_ps(xBv, _mm_mul_ps(*dp, *tp)); dp+= p7F_NSCELLS; tp += 7;
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
save_debug_row_pp(P7_FILTERMX *ox, __m128 *dpc, int i)
{
  union { __m128 v; float x[p7_VNF]; } u;
  int      Q  = ox->Qf;
  float  *xc  = (float *) (dpc + Q*p7F_NSCELLS);
  float **dp;
  float  *xmx;
  int     q,k,z;

  if (! ox->pp) return;
  dp  = ox->pp->dp;    	/* sets up {MDI}MX() macros in <pp> */
  xmx = ox->pp->xmx;	/* sets up XMX() macro in <pp>      */
  
  XMX(i,p7G_E) = xc[p7F_E];
  XMX(i,p7G_N) = xc[p7F_N];
  XMX(i,p7G_J) = xc[p7F_J];
  XMX(i,p7G_B) = xc[p7F_B];
  XMX(i,p7G_C) = xc[p7F_C];
  
  MMX(i,0) = DMX(i,0) = IMX(i,0) = 0.0;
  for (q = 0; q < Q; q++)
    {
      u.v = P7F_MQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) MMX(i,k) = u.x[z]; }
      u.v = P7F_DQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) DMX(i,k) = u.x[z]; }
      u.v = P7F_IQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) IMX(i,k) = u.x[z]; }
    }
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
save_debug_row_fb(P7_FILTERMX *ox, P7_GMX *gx, __m128 *dpc, int i, float totscale)
{
  union { __m128 v; float x[p7_VNF]; } u;
  int      Q  = ox->Qf;
  float  *xc  = (float *) (dpc + Q*p7F_NSCELLS);
  float **dp;
  float  *xmx;
  int     q,k,z;

  if (! gx) return;
  dp  = gx->dp; 	/* sets up {MDI}MX() macros */
  xmx = gx->xmx;	/* sets up XMX() macro      */
  
  XMX(i,p7G_E) = logf(xc[p7F_E]) + totscale;
  XMX(i,p7G_N) = logf(xc[p7F_N]) + totscale;
  XMX(i,p7G_J) = logf(xc[p7F_J]) + totscale;
  XMX(i,p7G_B) = logf(xc[p7F_B]) + totscale;
  XMX(i,p7G_C) = logf(xc[p7F_C]) + totscale;
  
  MMX(i,0) = DMX(i,0) = IMX(i,0) = -eslINFINITY;
  for (q = 0; q < Q; q++)
    {
      u.v = P7F_MQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) MMX(i,k) = logf(u.x[z]) + totscale; }
      u.v = P7F_DQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) DMX(i,k) = logf(u.x[z]) + totscale; }
      u.v = P7F_IQ(dpc, q); for (z = 0; z < p7_VNF; z++) { k = q+Q*z+1; if (k <= ox->M) IMX(i,k) = logf(u.x[z]) + totscale; }
    }
}
#endif
/*---------------- end, debugging tools -------------------------*/


/*****************************************************************
 * 4. Benchmark
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
  bnd = p7_gbands_Create(om->M, L);

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
	p7_BackwardFilter(dsq, L, om, ox, bnd);

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
 * 5. Unit tests
 *****************************************************************/
#ifdef p7FWDFILTER_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

#include "esl_sqio.h"

/* Compare scores of Forward, Backward to those from the reference
 * implementation.
 * 
 */
static void
utest_scores(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, int64_t ramlimit)
{
  char        *msg    = "fbfilter scores unit test failed";
  P7_HMM      *hmm    = NULL;
  P7_PROFILE  *gm     = NULL;
  P7_OPROFILE *om     = NULL;
  ESL_DSQ     *dsqmem = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_DSQ     *dsq    = NULL;
  int          tL     = 0;
  ESL_SQ      *sq     = esl_sq_CreateDigital(abc);
  P7_FILTERMX *ox     = p7_filtermx_Create(M, L, ramlimit);
  P7_GMX      *fwd    = p7_gmx_Create   (M, L);
  P7_GMX      *bck    = p7_gmx_Create   (M, L);
  P7_GMX      *pp     = p7_gmx_Create   (M, L);
  P7_GBANDS   *bnd    = p7_gbands_Create(M, L);
  float        tol2   = ( p7_logsum_IsSlowExact() ? 0.001  : 0.1);   /* absolute agreement of reference (log-space) and vector (prob-space) depends on whether we're using LUT-based logsum() */
  float fsc1, fsc2;
  float bsc2;
#ifdef p7_DEBUGGING
  float        bsc1;
  float        tol1   = 0.0001;	                                     /* forward and backward scores from same implementation type should agree with high tolerance */
  float        ptol   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01);  /* posterior decoding values differ by no more than this */
#endif

#ifdef p7_DEBUGGING
  /* We set the debugging tools to record full pp, fwd, bck matrices
   * for comparison to reference implementation: 
   */
  ox->pp  = p7_gmx_Create(M,L);	
  ox->fwd = p7_gmx_Create(M,L);	
  ox->bck = p7_gmx_Create(M,L);	
#endif

  if ( p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om) != eslOK) esl_fatal(msg);

  while (N--)
    {
      /* A mix of generated (homologous) and random (nonhomologous) sequences */
      if (esl_rnd_Roll(r, 2)) 
	{
	  esl_rsq_xfIID(r, bg->f, abc->K, L, dsqmem);  
	  dsq = dsqmem;  
	  tL = L;     
	}
      else  
	{
	  do {
	    esl_sq_Reuse(sq);
	    p7_ProfileEmit(r, hmm, gm, bg, sq, NULL);
	  } while (sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
	  dsq = sq->dsq; 
	  tL = sq->n; 
	}
	
      if ( p7_gmx_GrowTo     (fwd, M, tL) != eslOK) esl_fatal(msg);
      if ( p7_gmx_GrowTo     (bck, M, tL) != eslOK) esl_fatal(msg);
      if ( p7_gmx_GrowTo     (pp,  M, tL) != eslOK) esl_fatal(msg);
      if ( p7_filtermx_GrowTo(ox,  M, tL) != eslOK) esl_fatal(msg);

      p7_ForwardFilter (dsq, tL, om, ox, &fsc1);
      p7_BackwardFilter(dsq, tL, om, ox,  bnd);

      p7_GForward (dsq, tL, gm, fwd,  &fsc2);
      p7_GBackward(dsq, tL, gm, bck,  &bsc2);
      p7_GDecoding(         gm, fwd, bck, pp);

#ifdef p7_DEBUGGING
      /* vector Forward and Backward scores should agree with high tolerance.
       * Backward score is only available in debugging mode 
       */
      bsc1 = ox->bcksc;
      if (fabs(fsc1-bsc1) > tol1)    esl_fatal(msg);
#endif

      /* reference scores should agree with tolerance depending on whether logsum compiled to use LUT or not */
      if (fabs(fsc2-bsc2) > tol2)    esl_fatal(msg);

      /* Reference and vector implementations should agree depending on logsum */
      if (fabs(fsc1-fsc2) > tol2) esl_fatal(msg);

#ifdef p7_DEBUGGING
      /* Compare all DP cell values to reference implementation,
       * in fwd, bck, and pp matrices.
       */
      if (p7_gmx_Compare(pp,  ox->pp,  ptol) != eslOK) esl_fatal(msg);
      if (p7_gmx_Compare(fwd, ox->fwd, tol2) != eslOK) esl_fatal(msg);
      if (p7_gmx_Compare(bck, ox->bck, tol2) != eslOK) esl_fatal(msg);
#endif
      esl_sq_Reuse(sq);
      p7_gmx_Reuse(fwd);
      p7_gmx_Reuse(bck);
      p7_gmx_Reuse(pp);
      p7_filtermx_Reuse(ox);
      p7_gbands_Reuse(bnd);
    }

  free(dsqmem);
  esl_sq_Destroy(sq);
  p7_gbands_Destroy(bnd);
  p7_filtermx_Destroy(ox);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_gmx_Destroy(pp);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7FWDFILTER_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/


/*****************************************************************
 * 6. Test driver
 *****************************************************************/
#ifdef p7FWDFILTER_TESTDRIVE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for checkpointed vector SSE Forward, Backward implementations";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  p7_FLogsumInit();
  impl_Init();

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(32));   /* normal sized models              */
  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(0));    /* zero memory: force checkpointing */
  utest_scores(r, abc, bg, 1, L, 10, ESL_MBYTES(32));   /* size 1 models                    */
  utest_scores(r, abc, bg, M, 1, 10, ESL_MBYTES(32));   /* size 1 sequences                 */

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(32));  
  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(0));   
  utest_scores(r, abc, bg, 1, L, 10, ESL_MBYTES(32));  
  utest_scores(r, abc, bg, M, 1, 10, ESL_MBYTES(32));  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7FWDFILTER_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/


/*****************************************************************
 * 7. Example
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
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "output in one line awkable format",                0 },
#ifdef p7_DEBUGGING
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump vector DP matrices for examination (verbose)",0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump recorded forward matrix",                     0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump recorded backward matrix",                    0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump recorded posterior prob matrix",              0 },
#endif
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
  float           fraw, nullsc, fsc, bsc;
  float           gfraw, gbraw, gfsc, gbsc;
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

  /* Initially allocate matrices for a default sequence size, 500; 
   * we resize as needed for each individual target seq 
   */
  ox  = p7_filtermx_Create(gm->M, 500, ESL_MBYTES(32));  
  gx  = p7_gmx_Create     (gm->M, 500);
  bnd = p7_gbands_Create(gm->M, 500);
#ifdef p7_DEBUGGING
  /* When the p7_DEBUGGING flag is up, <ox> matrix has the ability to
   * record generic, complete <fwd>, <bck>, and <pp> matrices, for
   * comparison to reference implementation, even when checkpointing.
   */
  if (esl_opt_GetBoolean(go, "-D")) p7_filtermx_SetDumpMode(ox, stdout, TRUE);
  if (esl_opt_GetBoolean(go, "-F")) ox->fwd = p7_gmx_Create(gm->M, 100);
  if (esl_opt_GetBoolean(go, "-B")) ox->bck = p7_gmx_Create(gm->M, 100);
  if (esl_opt_GetBoolean(go, "-P")) ox->pp  = p7_gmx_Create(gm->M, 100);
#endif

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_bg_SetLength(bg,            sq->n);

      p7_filtermx_GrowTo(ox,  om->M, sq->n); 
      p7_gmx_GrowTo     (gx,  gm->M, sq->n); 
      p7_gbands_Reinit  (bnd, gm->M, sq->n); 

      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);
    
      p7_ForwardFilter (sq->dsq, sq->n, om, ox, &fraw);
      p7_BackwardFilter(sq->dsq, sq->n, om, ox, bnd);

      p7_GForward     (sq->dsq, sq->n, gm, gx, &gfraw);
      p7_GBackward    (sq->dsq, sq->n, gm, gx, &gbraw);

      bsc = 0.0;		/* Backward score only available in debugging mode */
#ifdef p7_DEBUGGING
      if (esl_opt_GetBoolean(go, "-F")) p7_gmx_Dump(stdout, ox->fwd, p7_DEFAULT);
      if (esl_opt_GetBoolean(go, "-B")) p7_gmx_Dump(stdout, ox->bck, p7_DEFAULT);
      if (esl_opt_GetBoolean(go, "-P")) p7_gmx_Dump(stdout, ox->pp, p7_DEFAULT);
      bsc  =  (ox->bcksc-nullsc) / eslCONST_LOG2;
#endif

      fsc  =  (fraw-nullsc) / eslCONST_LOG2;
      gfsc = (gfraw-nullsc) / eslCONST_LOG2;
      gbsc = (gbraw-nullsc) / eslCONST_LOG2;
      P  = esl_exp_surv(fsc,   om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
      gP = esl_exp_surv(gfsc,  gm->evparam[p7_FTAU],  gm->evparam[p7_FLAMBDA]);

      gmem = (float) p7_gmx_Sizeof(gx) / 1000000.;
      cmem = (float) p7_filtermx_Sizeof(ox) / 1000000.;

      p7_gbands_Dump(stdout, bnd);	  

      if (esl_opt_GetBoolean(go, "-1")) 
	printf("%-30s\t%-20s\t%9.2g\t%7.4f\t%7.4f\t%9.2g\t%6.1f\t%6.2fM\t%6.2fM\n", sq->name, hmm->name, P, fsc, bsc, gP, gfsc, gmem, cmem);
      else
	{
	  
	  printf("target sequence:      %s\n",        sq->name);
	  printf("fwd filter raw score: %.4f nats\n", fraw);
#ifdef p7_DEBUGGING
	  printf("bck filter raw score: %.4f nats\n", ox->bcksc);
#endif
	  printf("null score:           %.2f nats\n", nullsc);
	  printf("per-seq score:        %.2f bits\n", fsc);
	  printf("P-value:              %g\n",        P);
	  printf("GForward raw score:   %.2f nats\n", gfraw);
	  printf("GBackward raw score:  %.2f nats\n", gbraw);
	  printf("GForward seq score:   %.2f bits\n", gfsc);
	  printf("GForward P-value:     %g\n",        gP);
	  printf("RAM, f/b filter:      %6.2fM\n",    cmem);
	  printf("RAM, generic:         %6.2fM\n",    gmem);
	}

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
 * 9. Notes
 *****************************************************************/

/* [a.]  On debugging and testing methods.
 * 
 *    When compiled with <p7_DEBUGGING> (specifically, when
 *    fwdfilter.c and p7_filtermx.[ch] are thus compiled), the
 *    <P7_FILTERMX> structure is augmented with additional fields for
 *    debugging dumps and unit test comparisons.
 *    
 *   %% Dumping vector matrices for examination
 *      The values from the vectorized <P7_FILTERMX> can be dumped
 *      during DP calculations by calling <p7_filtermx_SetDumpMode()>
 *      on the object. Dumping has to happen during DP, not after,
 *      because of the way checkpointing discards rows as it goes (and
 *      for posterior decoding, the implementation never stores a row
 *      at all). Dumped rows are prefixed by a tag "f1 O", "f1 X", "f2
 *      O", "f2 X", or "bck", indicating backward (bck), first pass
 *      Forward (f1), second pass Forward (f2), checkpointed rows (O)
 *      that get saved and recalled, and discarded rows (X) that get
 *      recalculated in the second Forward pass.
 *      
 *      This capability is most useful for examining small DP matrices
 *      by hand; see fwdfilter_example, -D option.
 *      
 *   %% Saving matrices for comparison to reference
 *      With the <p7_DEBUGGING> compile flag, a caller may
 *      additionally provide an allocated <P7_GMX> to the
 *      <P7_FILTERMX>, to enable storage of all DP matrix values in a
 *      form suitable for a <p7_gmx_Compare()> call against <P7_GMX>
 *      DP matrices calculated by the reference implementation.
 *      Caller does something like <ox->fwd = p7_gmx_Create(M,L)> to
 *      save a Forward matrix, and/or analogously for <ox->bck> and/or
 *      <ox->pp> for Backward and Decoding.
 *      
 *      This capability is most useful for unit tests and automated
 *      comparison to the reference implementation. See utest_scores().
 *      
 *   %% High-precision comparison to reference
 *      Normally the reference implementation uses a table-driven
 *      log-sum-exp approximation (see p7_logsum.c), in order to do
 *      stable numerical calculations in log space. This introduces
 *      nonnegligible numerical error into DP calculations, so
 *      comparisons between a probability space vector implementation
 *      and the reference implementation must allow a large amount of
 *      numeric slop. At a cost of about 20x in speed, if the
 *      p7_LOGSUM_SLOWEXACT flag is compiled in, the p7_FLogsum()
 *      function uses the (more) exact calculation, allowing DP cells
 *      values to be compared more stringently.
 *      
 *      This capability is reflected in unit tests that set tolerances
 *      for floating-point comparison, after checking the flag with a
 *      <p7_logsum_IsSlowExact()> call. See utest_scores() for
 *      example.
 */

/* [b.] Running time, in theory and in practice.
 *
 *    Checkpointing requires more time in theory, but in practice you
 *    probably won't notice. The checkpointing method requires
 *    recalculation of Forward rows that weren't checkpointed, meaning
 *    up to two Forwards passes (plus one Backwards and one posterior
 *    decoding pass) over the DP matrix, rather than one each. First,
 *    using the Forward score as a filter minimizes this penalty:
 *    relatively few sequences pass on to the BackwardFilter
 *    call. Second, memory management in the checkpointed P7_FILTERMX
 *    structure uses "partial checkpointing" to minimize the use of
 *    checkpointing; for most comparisons, of all but the longest
 *    query/target combinations, DP calculations will fit in the
 *    available memory, and checkpointing is not invoked. And finally,
 *    the implementation here has been further optimized, such that
 *    it's actually slightly faster (with checkpointing) than the
 *    original HMMER3.0 implementation of Forward and Backward without
 *    checkpointing.
 * 
 */

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/




                                          
