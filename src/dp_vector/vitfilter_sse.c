/* Viterbi filter implementation; SSE version.
 * 
 * This is a SIMD vectorized, striped, interleaved, one-row, reduced
 * precision (epi16) implementation of the Viterbi algorithm.
 * 
 * It calculates a close approximation of the Viterbi score, in
 * limited precision (signed words: 16 bits) and range. It may overflow on
 * high scoring sequences, but this indicates that the sequence is a
 * high-scoring hit worth examining more closely anyway.  It will not
 * underflow, in local alignment mode.
 * 
 * Contents:
 *   1. Viterbi filter implementation.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#if p7_CPU_ARCH == intel
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#endif /* intel arch */

#include "easel.h"
#include "esl_sse.h"

#include "esl_gumbel.h"

#include "base/p7_hmmwindow.h"
#include "search/p7_pipeline.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"


/*****************************************************************
 * 1. Viterbi filter implementation.
 *****************************************************************/

/* Function:  p7_ViterbiFilter()
 * Synopsis:  Calculates Viterbi score, vewy vewy fast, in limited precision.
 *
 * Purpose:   Calculates an approximation of the Viterbi score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and a preallocated one-row DP matrix <ox>. Return the 
 *            estimated Viterbi score (in nats) in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow. 
 *            
 *            <ox> will be resized if needed. It's fine if it was just
 *            <_Reuse()'d> from a previous, smaller profile comparison.
 *            
 *            The model must be in a local alignment mode; other modes
 *            cannot provide the necessary guarantee of no underflow.
 *            
 *            This is a striped SIMD Viterbi implementation using Intel
 *            SSE/SSE2 integer intrinsics \citep{Farrar07}, in reduced
 *            precision (signed words, 16 bits).
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: Viterbi score (in nats)          
 *
 * Returns:   <eslOK> on success;
 *            <eslERANGE> if the score overflows; in this case
 *            <*ret_sc> is <eslINFINITY>, and the sequence can 
 *            be treated as a high-scoring hit.
 *            <ox> may be reallocated.
 *
 * Xref:      [Farrar07] for ideas behind striped SIMD DP.
 *            J2/46-47 for layout of HMMER's striped SIMD DP.
 *            J2/50 for single row DP.
 *            J2/60 for reduced precision (epu8)
 *            J2/65 for initial benchmarking
 *            J2/66 for precision maximization
 *            J4/138-140 for reimplementation in 16-bit precision
 *            J9/110-111 for reimplementation with P7_FILTERMX, memory share w/ checkpointed DP matrix
 *            J10/101 for separating P7_FILTERMX from P7_CHECKPTMX again: don't share these
 */
int
p7_ViterbiFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
#ifdef HAVE_SSE2
//printf("Starting ViterbiFilter\n");
  int i;        /* counter over sequence positions 1..L                      */
  
  register __m128i mpv, dpv, ipv;  /* previous row values                                       */
  register __m128i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128i dcv;		   /* delayed storage of D(i,q+1)                               */
  register __m128i xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i Dmaxv;          /* keeps track of maximum D cell on row                      */
  int16_t  xE, xB, xC, xJ, xN;	   /* special states' scores                                    */
  int16_t  Dmax;		   /* maximum D cell score on row                               */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = P7_NVW(om->M);    /* segment length: # of vectors                              */
  __m128i *dp;
  __m128i *rsc;			   /* will point at om->ru[x] for residue x[i]                  */
  __m128i *tsc;			   /* will point into (and step thru) om->tu                    */
  __m128i negInfv;
  int     status;

  /* Contract checks */
  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  /* note however that ViterbiFilter numerics are only guaranteed for local alignment, not glocal        */


  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of Vit filter matrix failed");
  dp = ox->dp;           /* enables MMXf(), IMXf(), DMXf() access macros. Must be set AFTER the GrowTo, because ox->dp may get moved */

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_VITFILTER;

  /* -infinity is -32768 */
  negInfv = _mm_set1_epi16(-32768);
  negInfv = _mm_srli_si128(negInfv, 14);  /* negInfv = 16-byte vector, 14 0 bytes + 2-byte value=-32768, for an OR operation. */

  /* Initialization. In unsigned arithmetic, -infinity is -32768
   */
  for (q = 0; q < Q; q++)
    MMXf(q) = IMXf(q) = DMXf(q) = _mm_set1_epi16(-32768);
  xN   = om->base_w;
  xB   = xN + om->xw[p7O_N][p7O_MOVE];
  xJ   = -32768;
  xC   = -32768;
  xE   = -32768;

#ifdef p7_DEBUGGING
  if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, 0, xE, 0, xJ, xB, xC); /* first 0 is <rowi>: do header. second 0 is xN: always 0 here. */
#endif

  for (i = 1; i <= L; i++)
    {
// check each iteration         
      rsc   = om->rwv[dsq[i]];
      tsc   = om->twv;
      dcv   = _mm_set1_epi16(-32768);      /* "-infinity" */
      xEv   = _mm_set1_epi16(-32768);     
      Dmaxv = _mm_set1_epi16(-32768);     
      xBv   = _mm_set1_epi16(xB);

      /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically; replace it with -32768.
       */
      mpv = MMXf(Q-1);  mpv = _mm_slli_si128(mpv, 2);  mpv = _mm_or_si128(mpv, negInfv);
      dpv = DMXf(Q-1);  dpv = _mm_slli_si128(dpv, 2);  dpv = _mm_or_si128(dpv, negInfv);
      ipv = IMXf(Q-1);  ipv = _mm_slli_si128(ipv, 2);  ipv = _mm_or_si128(ipv, negInfv);

      for (q = 0; q < Q; q++)
      {
    
        /* Calculate new MMXf(i,q); don't store it yet, hold it in sv. */
        sv   =                    _mm_adds_epi16(xBv, *tsc);  tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(mpv, *tsc)); tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(ipv, *tsc)); tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(dpv, *tsc)); tsc++;
        sv   = _mm_adds_epi16(sv, *rsc);                      rsc++;
        xEv  = _mm_max_epi16(xEv, sv);

        /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
         * {MDI}MX(q) is then the current, not the prev row
         */
        mpv = MMXf(q);
        dpv = DMXf(q);
        ipv = IMXf(q);

        /* Do the delayed stores of {MD}(i,q) now that memory is usable */
        MMXf(q) = sv;
        DMXf(q) = dcv;

        /* Calculate the next D(i,q+1) partially: M->D only;
               * delay storage, holding it in dcv
         */
        dcv   = _mm_adds_epi16(sv, *tsc);  tsc++;
        Dmaxv = _mm_max_epi16(dcv, Dmaxv);

        /* Calculate and store I(i,q) */
        sv     =                    _mm_adds_epi16(mpv, *tsc);  tsc++;
        IMXf(q)= _mm_max_epi16 (sv, _mm_adds_epi16(ipv, *tsc)); tsc++;
      }

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      xE = esl_sse_hmax_epi16(xEv);
      if (xE >= 32767) { *ret_sc = eslINFINITY; return eslERANGE; }	/* immediately detect overflow */
      xN = xN + om->xw[p7O_N][p7O_LOOP];
      xC = ESL_MAX(xC + om->xw[p7O_C][p7O_LOOP], xE + om->xw[p7O_E][p7O_MOVE]);
      xJ = ESL_MAX(xJ + om->xw[p7O_J][p7O_LOOP], xE + om->xw[p7O_E][p7O_LOOP]);
      xB = ESL_MAX(xJ + om->xw[p7O_J][p7O_MOVE], xN + om->xw[p7O_N][p7O_MOVE]);
      /* and now xB will carry over into next i, and xC carries over after i=L */

      /* Finally the "lazy F" loop (sensu [Farrar07]). We can often
       * prove that we don't need to evaluate any D->D paths at all.
       *
       * The observation is that if we can show that on the next row,
       * B->M(i+1,k) paths always dominate M->D->...->D->M(i+1,k) paths
       * for all k, then we don't need any D->D calculations.
       * 
       * The test condition is:
       *      max_k D(i,k) + max_k ( TDD(k-2) + TDM(k-1) - TBM(k) ) < xB(i)
       * So:
       *   max_k (TDD(k-2) + TDM(k-1) - TBM(k)) is precalc'ed in om->dd_bound;
       *   max_k D(i,k) is why we tracked Dmaxv;
       *   xB(i) was just calculated above.
       */
      Dmax = esl_sse_hmax_epi16(Dmaxv);
      if (Dmax + om->ddbound_w > xB) 
	{
   // printf("vitfilter SSE doing DD calc, Dmax = %i\n", Dmax);
	  /* Now we're obligated to do at least one complete DD path to be sure. */
	  /* dcv has carried through from end of q loop above */
	  dcv = _mm_slli_si128(dcv, 2); 
	  dcv = _mm_or_si128(dcv, negInfv);
	  tsc = om->twv + 7*Q;	/* set tsc to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    {
	      DMXf(q) = _mm_max_epi16(dcv, DMXf(q));	
	      dcv     = _mm_adds_epi16(DMXf(q), *tsc); tsc++;
	    }

	  /* We may have to do up to three more passes; the check
	   * is for whether crossing a segment boundary can improve
	   * our score. 
	   */
	  do {
	    dcv = _mm_slli_si128(dcv, 2);
	    dcv = _mm_or_si128(dcv, negInfv);
	    tsc = om->twv + 7*Q;	/* set tsc to start of the DD's */
	    for (q = 0; q < Q; q++) 
	      {
		if (! esl_sse_any_gt_epi16(dcv, DMXf(q))) break;
		DMXf(q) = _mm_max_epi16(dcv, DMXf(q));	
		dcv     = _mm_adds_epi16(DMXf(q), *tsc);   tsc++;
	      }	    
	  } while (q == Q);
	}
      else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
	{
 //   printf("vitfilter SSE skipping DD calc, Dmax = %i\n", Dmax);
	  dcv = _mm_slli_si128(dcv, 2);
	  DMXf(0) = _mm_or_si128(dcv, negInfv);
	}

#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, i, xE, 0, xJ, xB, xC);   
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T */
  if (xC > -32768)
    {
      *ret_sc = (float) xC + (float) om->xw[p7O_C][p7O_MOVE] - (float) om->base_w;
      /* *ret_sc += L * om->ncj_roundoff;  see J4/150 for rationale: superceded by -3.0nat approximation*/
      *ret_sc /= om->scale_w;
      *ret_sc -= 3.0; /* the NN/CC/JJ=0,-3nat approximation: see J5/36. That's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ contrib */
    }
  else  *ret_sc = -eslINFINITY;
#endif // HAVE_SSE2   Leave return so that there's some function to link if we aren't compiling SSE.
  return eslOK;
}
/*---------------- end, p7_ViterbiFilter_sse() ----------------------*/

/* Function:  p7_ViterbiFilter_longtarget()
 * Synopsis:  Finds windows within potentially long sequence blocks with Viterbi
 *            scores above threshold (vewy vewy fast, in limited precision)
 *
 * Purpose:   Calculates an approximation of the Viterbi score for regions
 *            of sequence <dsq>, using optimized profile <om>, and a pre-
 *            allocated one-row DP matrix <ox>, and captures the positions
 *            at which such regions exceed the score required to be
 *            significant in the eyes of the calling function (usually
 *            p=0.001).
 *
 *            The resulting landmarks are converted to subsequence
 *            windows by the calling function
 *
 *            The model must be in a local alignment mode; other modes
 *            cannot provide the necessary guarantee of no underflow.
 *
 *            This is a striped SIMD Viterbi implementation using Intel
 *            SSE/SSE2 integer intrinsics \citep{Farrar07}, in reduced
 *            precision (signed words, 16 bits).
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues
 *            om      - optimized profile
 *            ox      - DP matrix
 *            filtersc   - null or bias correction, required for translating a P-value threshold into a score threshold
 *            P          - p-value below which a region is captured as being above threshold
 *            windowlist - RETURN: preallocated array of hit windows (start and end of diagonal) for the above-threshold areas
 *
 * Returns:   <eslOK> on success;
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small, or if
 *            profile isn't in a local alignment mode. (Must be in local
 *            alignment mode because that's what helps us guarantee
 *            limited dynamic range.)
 *
 * Xref:      See p7_ViterbiFilter()
 */
int
p7_ViterbiFilter_longtarget_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox,
                            float filtersc, double P, P7_HMM_WINDOWLIST *windowlist)
{
  #ifdef p7_build_SSE // this funciton has not been converted to AVX
  register __m128i mpv, dpv, ipv; /* previous row values                                       */
  register __m128i sv;            /* temp storage of 1 curr row value in progress              */
  register __m128i dcv;           /* delayed storage of D(i,q+1)                               */
  register __m128i xEv;           /* E state: keeps max for Mk->E as we go                     */
  register __m128i xBv;           /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i Dmaxv;         /* keeps track of maximum D cell on row                      */
  int16_t  xE, xB, xC, xJ, xN;    /* special states' scores                                    */
  int16_t  Dmax;                  /* maximum D cell score on row                               */
  int i;                          /* counter over sequence positions 1..L                      */
  int q;                          /* counter over vectors 0..nq-1                              */
  int Q        = P7_NVW(om->M);   /* segment length: # of vectors                              */
  __m128i *dp  = ox->dp;          /* using {MDI}MXf(q) macro requires initialization of <dp>   */
  __m128i *rsc;                   /* will point at om->ru[x] for residue x[i]                  */
  __m128i *tsc;                   /* will point into (and step thru) om->tu                    */
  __m128i negInfv;
  int16_t sc_thresh;
  float   invP;
  int z;
  union { __m128i v; int16_t i[8]; } tmp;
  int status;

  windowlist->count = 0;

  /* Contract checks / argument validation */
  ESL_DASSERT1( (om->mode == p7_LOCAL || om->mode == p7_UNILOCAL) ); /* VF numerics only work for local alignment */

  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of Vit filter matrix failed");

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_VITFILTER;

/*
 *  In p7_ViterbiFilter, converting from a scaled int Viterbi score
 *  S (aka xE the score getting to state E) to a probability
 *  goes like this:
 *    vsc =  S + om->xw[p7O_E][p7O_MOVE] + om->xw[p7O_C][p7O_MOVE] - om->base_w
 *    ret_sc /= om->scale_w;
 *    vsc -= 3.0;
 *    P  = esl_gumbel_surv((vfsc - filtersc) / eslCONST_LOG2  ,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
 *  and we're computing the threshold vsc, so invert it:
 *    (vsc - filtersc) /  eslCONST_LOG2 = esl_gumbel_invsurv( P, om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA])
 *    vsc = filtersc + eslCONST_LOG2 * esl_gumbel_invsurv( P, om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA])
 *    vsc += 3.0
 *    vsc *= om->scale_w
 *    S = vsc - (float)om->xw[p7O_E][p7O_MOVE] - (float)om->xw[p7O_C][p7O_MOVE] + (float)om->base_w
 */
  invP = esl_gumbel_invsurv(P, om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
  sc_thresh =   (int) ceil ( ( (filtersc + (eslCONST_LOG2 * invP) + 3.0) * om->scale_w )
                - (float)om->xw[p7O_E][p7O_MOVE] - (float)om->xw[p7O_C][p7O_MOVE] + (float)om->base_w );

  /* -infinity is -32768 */
  negInfv = _mm_set1_epi16(-32768);
  negInfv = _mm_srli_si128(negInfv, 14);  /* negInfv = 16-byte vector, 14 0 bytes + 2-byte value=-32768, for an OR operation. */

  /* Initialization. In unsigned arithmetic, -infinity is -32768
   */
  for (q = 0; q < Q; q++)
    MMXf(q) = IMXf(q) = DMXf(q) = _mm_set1_epi16(-32768);
  xN   = om->base_w;
  xB   = xN + om->xw[p7O_N][p7O_MOVE];
  xJ   = -32768;
  xC   = -32768;
  xE   = -32768;

#if p7_DEBUGGING
  if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, 0, xE, 0, xJ, xB, xC); /* first 0 is <rowi>: do header. second 0 is xN: always 0 here. */
#endif


  for (i = 1; i <= L; i++)
  {
      rsc   = om->rwv[dsq[i]];
      tsc   = om->twv;
      dcv   = _mm_set1_epi16(-32768);      /* "-infinity" */
      xEv   = _mm_set1_epi16(-32768);
      Dmaxv = _mm_set1_epi16(-32768);
      xBv   = _mm_set1_epi16(xB);

      /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12.
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically; replace it with -32768.
       */
      mpv = MMXf(Q-1);  mpv = _mm_slli_si128(mpv, 2);  mpv = _mm_or_si128(mpv, negInfv);
      dpv = DMXf(Q-1);  dpv = _mm_slli_si128(dpv, 2);  dpv = _mm_or_si128(dpv, negInfv);
      ipv = IMXf(Q-1);  ipv = _mm_slli_si128(ipv, 2);  ipv = _mm_or_si128(ipv, negInfv);

      for (q = 0; q < Q; q++)
      {
        /* Calculate new MMXf(i,q); don't store it yet, hold it in sv. */
        sv   =                    _mm_adds_epi16(xBv, *tsc);  tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(mpv, *tsc)); tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(ipv, *tsc)); tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(dpv, *tsc)); tsc++;
        sv   = _mm_adds_epi16(sv, *rsc);                      rsc++;
        xEv  = _mm_max_epi16(xEv, sv);

        /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
         * {MDI}MX(q) is then the current, not the prev row
         */
        mpv = MMXf(q);
        dpv = DMXf(q);
        ipv = IMXf(q);

        /* Do the delayed stores of {MD}(i,q) now that memory is usable */
        MMXf(q) = sv;
        DMXf(q) = dcv;

        /* Calculate the next D(i,q+1) partially: M->D only;
               * delay storage, holding it in dcv
         */
        dcv   = _mm_adds_epi16(sv, *tsc);  tsc++;
        Dmaxv = _mm_max_epi16(dcv, Dmaxv);

        /* Calculate and store I(i,q) */
        sv     =                    _mm_adds_epi16(mpv, *tsc);  tsc++;
        IMXf(q)= _mm_max_epi16 (sv, _mm_adds_epi16(ipv, *tsc)); tsc++;
      }

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      xE = esl_sse_hmax_epi16(xEv);

      if (xE >= sc_thresh) {
        //hit score threshold. Add a window to the list, then reset scores.

        /* Unpack and unstripe, then find the position responsible for the hit */
        for (q = 0; q < Q; q++) {
          tmp.v = MMXf(q);
          for (z = 0; z < 8; z++)  { // unstripe
            if ( tmp.i[z] == xE && (q+Q*z+1) <= om->M) {
              // (q+Q*z+1) is the model position k at which the xE score is found
              p7_hmmwindow_new(windowlist, 0, i, 0, (q+Q*z+1), 1, 0.0, p7_NOCOMPLEMENT );
            }
          }
          MMXf(q) = IMXf(q) = DMXf(q) = _mm_set1_epi16(-32768); //reset score to start search for next vit window.
        }

      } else {


        xN = xN + om->xw[p7O_N][p7O_LOOP];
        xC = ESL_MAX(xC + om->xw[p7O_C][p7O_LOOP], xE + om->xw[p7O_E][p7O_MOVE]);
        xJ = ESL_MAX(xJ + om->xw[p7O_J][p7O_LOOP], xE + om->xw[p7O_E][p7O_LOOP]);
        xB = ESL_MAX(xJ + om->xw[p7O_J][p7O_MOVE], xN + om->xw[p7O_N][p7O_MOVE]);
        /* and now xB will carry over into next i, and xC carries over after i=L */

        /* Finally the "lazy F" loop (sensu [Farrar07]). We can often
         * prove that we don't need to evaluate any D->D paths at all.
         *
         * The observation is that if we can show that on the next row,
         * B->M(i+1,k) paths always dominate M->D->...->D->M(i+1,k) paths
         * for all k, then we don't need any D->D calculations.
         *
         * The test condition is:
         *      max_k D(i,k) + max_k ( TDD(k-2) + TDM(k-1) - TBM(k) ) < xB(i)
         * So:
         *   max_k (TDD(k-2) + TDM(k-1) - TBM(k)) is precalc'ed in om->dd_bound;
         *   max_k D(i,k) is why we tracked Dmaxv;
         *   xB(i) was just calculated above.
         */
        Dmax = esl_sse_hmax_epi16(Dmaxv);
        if (Dmax + om->ddbound_w > xB)
        {
          /* Now we're obligated to do at least one complete DD path to be sure. */
          /* dcv has carried through from end of q loop above */
          dcv = _mm_slli_si128(dcv, 2);
          dcv = _mm_or_si128(dcv, negInfv);
          tsc = om->twv + 7*Q;  /* set tsc to start of the DD's */
          for (q = 0; q < Q; q++)
          {
            DMXf(q) = _mm_max_epi16(dcv, DMXf(q));
            dcv     = _mm_adds_epi16(DMXf(q), *tsc); tsc++;
          }

          /* We may have to do up to three more passes; the check
           * is for whether crossing a segment boundary can improve
           * our score.
           */
          do {
            dcv = _mm_slli_si128(dcv, 2);
            dcv = _mm_or_si128(dcv, negInfv);
            tsc = om->twv + 7*Q;  /* set tsc to start of the DD's */
            for (q = 0; q < Q; q++)
            {
              if (! esl_sse_any_gt_epi16(dcv, DMXf(q))) break;
              DMXf(q) = _mm_max_epi16(dcv, DMXf(q));
              dcv     = _mm_adds_epi16(DMXf(q), *tsc);   tsc++;
            }
          } while (q == Q);
        }
        else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
        {
          dcv = _mm_slli_si128(dcv, 2);
          DMXf(0) = _mm_or_si128(dcv, negInfv);
        }
      }
#if p7_DEBUGGING
      if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, i, xE, 0, xJ, xB, xC);
#endif
  } /* end loop over sequence residues 1..L */
#endif /* p7_build_SSE */

  return eslOK;

}
/*---------------- end, p7_ViterbiFilter_longtarget() ----------------------*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

