/* The MSV filter implementation; SSE version.
 * 
 * A "filter" is a one-row, O(M), DP implementation that calculates an
 * approximated nat score (i.e. in limited precision - here, uchar)
 * and may have limited numeric range. It will return <eslERANGE> if
 * its numeric range is exceeded, in which case the caller will have
 * to obtain the score by another (probably slower) method.
 * 
 * Contents:
 *   1. p7_MSVFilter() implementation
 *   2. Benchmark driver
 *   3. Unit tests
 *   4. Test driver
 *   5. Example
 *   6. Copyright and license information
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"
#include "esl_gumbel.h"

#include "base/p7_hmmwindow.h"
#include "fm/fm.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"
#include "dp_vector/ssvfilter.h"
#include "dp_vector/msvfilter.h"


/*****************************************************************
 * 1. The p7_MSVFilter() DP implementation.
 *****************************************************************/
 
/* Function:  p7_MSVFilter()
 * Synopsis:  Calculates MSV score, vewy vewy fast, in limited precision.
 *
 * Purpose:   Calculates an approximation of the MSV score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and the one-row DP matrix <ox>. Return the 
 *            estimated MSV score (in nats) in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow.
 *            
 *            <ox> will be resized if needed. It's fine if it was
 *            just <_Reuse()'d> from a previous, smaller profile.
 *            
 *            The model may be in any mode, because only its match
 *            emission scores will be used. The MSV filter inherently
 *            assumes a multihit local mode, and uses its own special
 *            state transition scores, not the scores in the profile.
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - filter DP matrix (one row)
 *            ret_sc  - RETURN: MSV score (in nats)          
 *                      
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the score overflows the limited range; in
 *            this case, this is a high-scoring hit.
 *            <ox> may have been resized.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_MSVFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
  register __m128i mpv;            /* previous row values                                       */
  register __m128i xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128i biasv;	   /* emission bias in a vector                                 */
  uint8_t  xJ;                     /* special states' scores                                    */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = P7_NVB(om->M);    /* segment length: # of vectors                              */
  __m128i *dp;   	           /* the dp row memory                                         */
  __m128i *rsc;			   /* will point at om->rbv[x] for residue x[i]                 */
  __m128i xJv;                     /* vector for states score                                   */
  __m128i tjbmv;                   /* vector for cost of moving {JN}->B->M                      */
  __m128i tecv;                    /* vector for E->C  cost                                     */
  __m128i basev;                   /* offset for scores                                         */
  __m128i ceilingv;                /* saturated simd value used to test for overflow            */
  __m128i tempv;                   /* work vector                                               */
  int     cmp;
  int     status;

  /* Try highly optimized Knudsen SSV filter first. 
   * Note that SSV doesn't use any main memory (from <ox>) at all! 
   */
  if (( status = p7_SSVFilter(dsq, L, om, ret_sc)) != eslENORESULT) return status;

  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of MSV filter matrix failed");
  dp = ox->dp;			/* This MUST be set AFTER the GrowTo() call. */

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_MSVFILTER;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.  */
  biasv = _mm_set1_epi8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */
  for (q = 0; q < Q; q++) dp[q] = _mm_setzero_si128();
  xJ   = 0;

  /* saturate simd register for overflow test */
  ceilingv = _mm_cmpeq_epi8(biasv, biasv);
  basev    = _mm_set1_epi8((int8_t) om->base_b);
  tjbmv    = _mm_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  tecv     = _mm_set1_epi8((int8_t) om->tec_b);
  xJv      = _mm_subs_epu8(biasv, biasv);
  xBv      = _mm_subs_epu8(basev, tjbmv);

#ifdef p7_DEBUGGING
  if (ox->do_dumping)
    {
      uint8_t xB;
      xB = _mm_extract_epi16(xBv, 0);
      xJ = _mm_extract_epi16(xJv, 0);
      p7_filtermx_DumpMFRow(ox, 0, 0, 0, xJ, xB, xJ);
    }
#endif


  for (i = 1; i <= L; i++)
    {
      rsc = om->rbv[dsq[i]];
      xEv = _mm_setzero_si128();      

      /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically, which is our -infinity.
       */
      mpv = _mm_slli_si128(dp[Q-1], 1);   
      for (q = 0; q < Q; q++)
      {
        /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
        sv   = _mm_max_epu8(mpv, xBv);
        sv   = _mm_adds_epu8(sv, biasv);
        sv   = _mm_subs_epu8(sv, *rsc);   rsc++;
        xEv  = _mm_max_epu8(xEv, sv);

        mpv   = dp[q];   	  /* Load {MDI}(i-1,q) into mpv */
        dp[q] = sv;       	  /* Do delayed store of M(i,q) now that memory is usable */
      }

      /* test for the overflow condition */
      tempv = _mm_adds_epu8(xEv, biasv);
      tempv = _mm_cmpeq_epi8(tempv, ceilingv);
      cmp   = _mm_movemask_epi8(tempv);

      /* Now the "special" states, which start from Mk->E (->C, ->J->B)
       * Use shuffles instead of shifts so when the last max has completed,
       * the last four elements of the simd register will contain the
       * max value.  Then the last shuffle will broadcast the max value
       * to all simd elements.
       */
      tempv = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(2, 3, 0, 1));
      xEv   = _mm_max_epu8(xEv, tempv);
      tempv = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(0, 1, 2, 3));
      xEv   = _mm_max_epu8(xEv, tempv);
      tempv = _mm_shufflelo_epi16(xEv, _MM_SHUFFLE(2, 3, 0, 1));
      xEv   = _mm_max_epu8(xEv, tempv);
      tempv = _mm_srli_si128(xEv, 1);
      xEv   = _mm_max_epu8(xEv, tempv);
      xEv   = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(0, 0, 0, 0));

      /* immediately detect overflow */
      if (cmp != 0x0000) { *ret_sc = eslINFINITY; return eslERANGE; }

      xEv = _mm_subs_epu8(xEv, tecv);
      xJv = _mm_max_epu8(xJv,xEv);
      xBv = _mm_max_epu8(basev, xJv);
      xBv = _mm_subs_epu8(xBv, tjbmv);
	  
#ifdef p7_DEBUGGING
      if (ox->do_dumping)
	{
	  uint8_t xB, xE;
	  xB = _mm_extract_epi16(xBv, 0);
	  xE = _mm_extract_epi16(xEv, 0);
	  xJ = _mm_extract_epi16(xJv, 0);
	  p7_filtermx_DumpMFRow(ox, i, xE, 0, xJ, xB, xJ);
	}
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  xJ = (uint8_t) _mm_extract_epi16(xJv, 0);
  *ret_sc = ((float) (xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
  return eslOK;
}
/*------------------ end, p7_MSVFilter() ------------------------*/



/* Function:  p7_SSVFilter_longtarget()
 * Synopsis:  Finds windows with SSV scores above some threshold (vewy vewy fast, in limited precision)
 *
 * Purpose:   Calculates an approximation of the MSV score for regions of
 *            sequence <dsq> of length <L> residues, using optimized profile
 *            <om>, and a preallocated one-row DP matrix <ox>, and captures
 *            the positions at which such regions exceed the score required
 *            to be significant in the eyes of the calling function (usually
 *            p=0.02).  Note that this variant performs only SSV computations,
 *            never passing through the J state. See comments in
 *            p7_MSVFilter_longtarget() for details
 *
 *            Rather than simply capturing positions at which a score threshold
 *            is reached, this function establishes windows around those
 *            high-scoring positions. p7_MSVFilter_longtarget then merges
 *            overlapping windows.
 *
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues
 *            om      - optimized profile
 *            ox      - DP matrix
 *            msvdata    - compact representation of substitution scores, for backtracking diagonals
 *            sc_thresh  - a precomputed threshold required for a window to be considered "high-scoring"
 *            sc_threshv - a precomputed vector containing a vector of sc_thresh values
 *            starts  - RETURN: array of start positions for windows surrounding above-threshold areas
 *            ends    - RETURN: array of end positions for windows surrounding above-threshold areas
 *			  hit_cnt - RETURN: count of entries in the above two arrays
 *
 *
 * Note:      We misuse the matrix <ox> here, using only a third of the
 *            first dp row, accessing it as <dp[0..Q-1]> rather than
 *            in triplets via <{MDI}MX(q)> macros, since we only need
 *            to store M state values. We know that if <ox> was big
 *            enough for normal DP calculations, it must be big enough
 *            to hold the MSVFilter calculation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
static int
p7_SSVFilter_longtarget(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_FILTERMX *ox, const P7_SCOREDATA *msvdata,
                     uint8_t sc_thresh, __m128i sc_threshv, P7_HMM_WINDOWLIST *windowlist)
{

  register __m128i mpv;            /* previous row values                                       */
  register __m128i xEv;		   /* E state: keeps max for Mk->E for a single iteration       */
  register __m128i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128i biasv;	   /* emission bias in a vector                                 */
  uint8_t  xJ;                     /* special states' scores                                    */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = P7_NVB(om->M);    /* segment length: # of vectors                              */
  __m128i *dp  = ox->dp;           /* one DP row of memory                                      */
  __m128i *rsc;			   /* will point at om->rbv[x] for residue x[i]                 */
  __m128i tecv;                    /* vector for E->C  cost                                     */
  __m128i tjbmv;                   /* vector for J->B move cost + B->M move costs               */
  __m128i basev;                   /* offset for scores                                         */
  __m128i ceilingv;                /* saturated simd value used to test for overflow            */
  __m128i tempv;                   /* work vector                                               */
  int cmp;
  int k;
  int n;
  int end;
  int rem_sc;
  int start;
  int target_end;
  int target_start;
  int max_end;
  int max_sc;
  int sc;
  int pos_since_max;
  float ret_sc;
  union { __m128i v; uint8_t b[16]; } u;
  int   status;

  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of SSV filter matrix failed");

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_SSVFILTER;

  /* Initialization. In offset unsigned  arithmetic, -infinity is 0, and 0 is om->base.
   */
  biasv = _mm_set1_epi8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */
  ceilingv = _mm_cmpeq_epi8(biasv, biasv);
  for (q = 0; q < Q; q++) dp[q] = _mm_setzero_si128();
  xJ   = 0;

  basev = _mm_set1_epi8((int8_t) om->base_b);
  tecv = _mm_set1_epi8((int8_t) om->tec_b);
  tjbmv = _mm_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);

  xBv = _mm_subs_epu8(basev, tjbmv);

  for (i = 1; i <= L; i++) {
    rsc = om->rbv[dsq[i]];
    xEv = _mm_setzero_si128();

	  /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12.
	   * Because ia32 is littlendian, this means a left bit shift.
	   * Zeros shift on automatically, which is our -infinity.
	   */
	  mpv = _mm_slli_si128(dp[Q-1], 1);
	  for (q = 0; q < Q; q++) {
		  /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
		  sv   = _mm_max_epu8(mpv, xBv);
		  sv   = _mm_adds_epu8(sv, biasv);
		  sv   = _mm_subs_epu8(sv, *rsc);   rsc++;
		  xEv  = _mm_max_epu8(xEv, sv);

		  mpv   = dp[q];   	  /* Load {MDI}(i-1,q) into mpv */
		  dp[q] = sv;       	  /* Do delayed store of M(i,q) now that memory is usable */
	  }

	  /* test if the pthresh significance threshold has been reached;
	   * note: don't use _mm_cmpgt_epi8, because it's a signed comparison, which won't work on uint8s */
	  tempv = _mm_adds_epu8(xEv, sc_threshv);
	  tempv = _mm_cmpeq_epi8(tempv, ceilingv);
	  cmp = _mm_movemask_epi8(tempv);

	  if (cmp != 0) {  //hit pthresh, so add position to list and reset values

	    //figure out which model state hit threshold
	    end = -1;
	    rem_sc = -1;
	    for (q = 0; q < Q; q++) {  /// Unpack and unstripe, so we can find the state that exceeded pthresh
          u.v = dp[q];
          for (k = 0; k < 16; k++) { // unstripe
            //(q+Q*k+1) is the model position k at which the xE score is found
            if (u.b[k] >= sc_thresh && u.b[k] > rem_sc && (q+Q*k+1) <= om->M) {
              end = (q+Q*k+1);
              rem_sc = u.b[k];
            }
          }
          dp[q] = _mm_set1_epi8(0); // while we're here ... this will cause values to get reset to xB in next dp iteration
	    }

	    //recover the diagonal that hit threshold
	    start = end;
	    target_end = target_start = i;
	    sc = rem_sc;
	    while (rem_sc > om->base_b - om->tjb_b - om->tbm_b) {
	      rem_sc -= om->bias_b -  msvdata->msv_scores[start*om->abc->Kp + dsq[target_start]];
	      --start;
	      --target_start;
	    }
	    start++;
	    target_start++;


	    //extend diagonal further with single diagonal extension
	    k = end+1;
	    n = target_end+1;
	    max_end = target_end;
	    max_sc = sc;
	    pos_since_max = 0;
	    while (k<om->M && n<=L) {
	      sc += om->bias_b -  msvdata->msv_scores[k*om->abc->Kp + dsq[n]];

	      if (sc >= max_sc) {
	        max_sc = sc;
	        max_end = n;
	        pos_since_max=0;
	      } else {
	        pos_since_max++;
	        if (pos_since_max == 5)
	          break;
	      }
	      k++;
	      n++;
	    }

	    end  +=  (max_end - target_end);
	    k    +=  (max_end - target_end);
      target_end = max_end;

      ret_sc = ((float) (max_sc - om->tjb_b) - (float) om->base_b);
      ret_sc /= om->scale_b;
      ret_sc -= 3.0; // that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ

      p7_hmmwindow_new(windowlist, 0, target_start, k, end, end-start+1 , ret_sc, fm_nocomplement);

      i = target_end; // skip forward
	  }


  } /* end loop over sequence residues 1..L */

  return eslOK;

}
/*------------------ end, p7_SSVFilter_longtarget() ------------------------*/



/* Function:  p7_MSVFilter_longtarget()
 * Synopsis:  Finds windows with MSV scores above some threshold (vewy vewy fast, in limited precision)
 *
 * Purpose:   Calculates an approximation of the MSV score for regions of
 *            sequence <dsq>, using optimized profile <om>, and a pre-
 *            allocated one-row DP matrix <ox>, and captures the positions
 *            at which such regions exceed the score required to be
 *            significant in the eyes of the calling function (usually
 *            p=0.02).  Note that this will typically simply call
 *            p7_SSVFilter_longtarget(), as the score threshold is nearly
 *            always lower than the score required for a pass through the
 *            J-state to yield a positive score.
 *
 *
 *
 * Args:      dsq        - digital target sequence, 1..L
 *            L          - length of dsq in residues
 *            om         - optimized profile
 *            ox         - DP matrix
 *            msvdata    - compact representation of substitution scores, for backtracking diagonals
 *            bg         - the background model, required for translating a P-value threshold into a score threshold
 *            P          - p-value below which a region is captured as being above threshold
 *            windowlist - RETURN: preallocated array of hit windows (start and end of diagonal) for the above-threshold areas
 *
 * Note:      We misuse the matrix <ox> here, using only a third of the
 *            first dp row, accessing it as <dp[0..Q-1]> rather than
 *            in triplets via <{MDI}MX(q)> macros, since we only need
 *            to store M state values. We know that if <ox> was big
 *            enough for normal DP calculations, it must be big enough
 *            to hold the MSVFilter calculation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_MSVFilter_longtarget(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_FILTERMX *ox, const P7_SCOREDATA *msvdata, P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist)
{
#ifdef ALLOW_MSV
  register __m128i mpv;            /* previous row values                                       */
  register __m128i xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128i biasv;	   /* emission bias in a vector                                 */
  uint8_t  xJ;                     /* special states' scores                                    */
  int i;			   /* counter over sequence positions 1..L                      */
//  int j;
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = P7_NVB(om->M);    /* segment length: # of vectors                              */
  __m128i *dp  = ox->dp;           /* DP row                                                    */
  __m128i *rsc;		           /* will point at om->rbv[x] for residue x[i]                 */
  __m128i xEv_prev;      	   /* E state: keeps the max seen up to the current position i, for window overlap test.  Note that this will always just be xJv + tecv    */
  __m128i xJv;                     /* vector for states score                                   */
  __m128i tecv;                    /* vector for E->C  cost                                     */
  __m128i tjbmv;                   /* vector for [JN]->B->M move cost                           */
  __m128i basev;                   /* offset for scores                                         */
  __m128i jthreshv;                /* pushes value to saturation if it's above pthresh          */
  __m128i ceilingv;                /* saturated simd value used to test for overflow            */
  __m128i tempv;                   /* work vector                                               */
#endif

  /*
   * Computing the score required to let P meet the F1 prob threshold
   * In original code, converting from a scaled int MSV
   * score S (the score getting to state E) to a probability goes like this:
   *  usc =  S - om->tec_b - om->tjb_b - om->base_b;
   *  usc /= om->scale_b;
   *  usc -= 3.0;
   *  P = f ( (usc - nullsc) / eslCONST_LOG2 , mu, lambda)
   * and we're computing the threshold usc, so reverse it:
   *  (usc - nullsc) /  eslCONST_LOG2 = inv_f( P, mu, lambda)
   *  usc = nullsc + eslCONST_LOG2 * inv_f( P, mu, lambda)
   *  usc += 3
   *  usc *= om->scale_b
   *  S = usc + om->tec_b + om->tjb_b + om->base_b
   *
   *  Here, I compute threshold with length model based on max_length.  Doesn't
   *  matter much - in any case, both the bg and om models will change with roughly
   *  1 bit for each doubling of the length model, so they offset.
   */
  float nullsc;
  __m128i sc_threshv;
  uint8_t sc_thresh;

  float invP = esl_gumbel_invsurv(P, om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);

  p7_bg_SetLength(bg, om->max_length);
  p7_oprofile_ReconfigMSVLength(om, om->max_length);
  p7_bg_NullOne  (bg, dsq, om->max_length, &nullsc);

  sc_thresh = (int) ceil( ( ( nullsc  + (invP * eslCONST_LOG2) + 3.0 )  * om->scale_b ) + om->base_b +  om->tec_b  + om->tjb_b );
  sc_threshv = _mm_set1_epi8((int8_t) 255 - sc_thresh);

#ifdef ALLOW_MSV
  //see ~wheelert/notebook/2012/0326-nhmmer-msv-ssv/00NOTES  - notes on 'Mar 31'
  uint8_t jthresh = om->base_b + om->tec_b + 1;
  jthreshv = _mm_set1_epi8((int8_t) 255 - (int)jthresh);


  if (force_ssv || (sc_thresh < jthresh)) {
#endif

    p7_SSVFilter_longtarget(dsq, L, om, ox, msvdata, (uint8_t)sc_thresh, sc_threshv, windowlist);

#ifdef ALLOW_MSV
  } else {

	  /*e.g. if base=190, tec=3, tjb=22 then a score of 217 would be required for a
	   * second ssv-hit to improve the score of an earlier one. Usually,
	   * sc_thresh < base. So only bother with msv if sc thresh is uncommonly
	   *  high. This will require mu of something in the neighborhood of -6.4 bits
	   *  (see TJW notes Mon Jun  7 10:19:34 EDT 2010 for details)
	   */
	  int hit_pthresh;
	  int e_unchanged;

	  //used to track neighboring MSV peaks
	  int first_peak = -1;
	  int last_peak = -1;
	  int hit_jthresh;
	  __m128i *dp  = ox->dp; 

	  /* Resize the filter mx as needed */
	  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of MSV filter matrix failed");
	  
	  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
	  ox->M    = om->M;
	  ox->type = p7F_MSVFILTER;

	  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.
	   */
	  biasv = _mm_set1_epi8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */
	  for (q = 0; q < Q; q++) dp[q] = _mm_setzero_si128();
	  xJ   = 0;


	  /* saturate simd register for overflow test */
	  ceilingv = _mm_cmpeq_epi8(biasv, biasv);

	  basev = _mm_set1_epi8((int8_t) om->base_b);

	  tecv = _mm_set1_epi8((int8_t) om->tec_b);
	  tjbmv = _mm_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);

	  xJv = _mm_subs_epu8(biasv, biasv);
	  xBv = _mm_subs_epu8(basev, tjbmv);

	  xEv_prev = _mm_setzero_si128();


	#ifdef p7_DEBUGGING
	  if (ox->do_dumping)
		{
		  uint8_t xB;
		  xB = _mm_extract_epi16(xBv, 0);
		  xJ = _mm_extract_epi16(xJv, 0);
		  p7_filtermx_DumpMFRow(ox, 0, 0, 0, xJ, xB, xJ);
		}
	#endif


	  for (i = 1; i <= L; i++)
	  {

		  rsc = om->rbv[dsq[i]];
		  xEv = _mm_setzero_si128();

		  /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12.
		   * Because ia32 is littlendian, this means a left bit shift.
		   * Zeros shift on automatically, which is our -infinity.
		   */
		  mpv = _mm_slli_si128(dp[Q-1], 1);
		  for (q = 0; q < Q; q++)
		  {
			  /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
			  sv   = _mm_max_epu8(mpv, xBv);
			  sv   = _mm_adds_epu8(sv, biasv);
			  sv   = _mm_subs_epu8(sv, *rsc);   rsc++;
			  xEv  = _mm_max_epu8(xEv, sv);

			  mpv   = dp[q];   	  /* Load {MDI}(i-1,q) into mpv */
			  dp[q] = sv;       	  /* Do delayed store of M(i,q) now that memory is usable */
		  }

		  if (last_peak == -1) { // haven't seen an above-jthresh peak recently, check if we've hit jthresh now
			  tempv = _mm_adds_epu8(xEv, jthreshv);
			  tempv = _mm_cmpeq_epi8(tempv, ceilingv);
			  hit_jthresh = _mm_movemask_epi8(tempv);
		  } else {
			  hit_jthresh = 0xffff; // already passed some above-jthresh peak in the recent past
		  }


		  /* Now the "special" states, which start from Mk->E (->C, ->J->B)
		   * Use shuffles instead of shifts so when the last max has completed,
		   * the last four elements of the simd register will contain the
		   * max value.  Then the last shuffle will broadcast the max value
		   * to all simd elements.
		   *
		   * This is placed here so it can be run while movemask above is blocking
		   */
		  tempv = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(2, 3, 0, 1));
		  xEv = _mm_max_epu8(xEv, tempv);
		  tempv = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(0, 1, 2, 3));
		  xEv = _mm_max_epu8(xEv, tempv);
		  tempv = _mm_shufflelo_epi16(xEv, _MM_SHUFFLE(2, 3, 0, 1));
		  xEv = _mm_max_epu8(xEv, tempv);
		  tempv = _mm_srli_si128(xEv, 1);
		  xEv = _mm_max_epu8(xEv, tempv);
		  xEv = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(0, 0, 0, 0));


		  if (hit_jthresh != 0x0000) {
			  //test if xE has increased this iteration.  If we go a long way without increasing, then we'll give up on the current window
			  tempv = _mm_max_epu8(xEv, xEv_prev);
			  tempv = _mm_cmpeq_epi8(tempv, xEv_prev); //if these are equal, then xEv didn't go up
			  e_unchanged = _mm_movemask_epi8(tempv);


			  if (e_unchanged == 0x0000) { // there was an increase. Establish peaks
				  xEv_prev = xEv;
				  if (first_peak == -1)
					  first_peak = i;
				  last_peak = i;

				  //check if we've hit pthresh.  If so, reset values and add this region to the list of hits
				  tempv = _mm_adds_epu8(xEv, sc_threshv);
				  tempv = _mm_cmpeq_epi8(tempv, ceilingv);
				  hit_pthresh = _mm_movemask_epi8(tempv);
				  if (hit_pthresh != 0)
				  {

			      int start = first_peak - om->max_length + 1;
			      int length = i-first_peak+1 + 2 * om->max_length;
			      fm_newWindow(windowlist, 0, start, -1, -1, length , -1, fm_nocomplement );

					  // got one peak.  Start scanning for the next one, as though starting from scratch
					  first_peak = last_peak = -1;

					  //reset xE, xJ, and xB
					  xJv = _mm_setzero_si128();  //_mm_subs_epu8(biasv, biasv);
					  xBv = _mm_subs_epu8(basev, tjbmv);
					  for (q = 0; q < Q; q++)
							dp[q] = _mm_set1_epi8(0); // this will cause values to get reset to xB in next iteration

				  }
			  } else if (i > last_peak + om->max_length) {
				  //been a long time since the last above-jthresh peak was seen;
				  // give up looking for a J-state-induced extension that will drive the score above sc_thresh

				  first_peak = last_peak = -1;

				  //Reset values as if the previous peak hadn't existed, retaining the growth of the current peak
				  //See end of TW notes feb 25, 2010
				  tempv = _mm_subs_epu8(xEv_prev, basev);
				  tempv = _mm_subs_epu8(tempv, tecv);
				  for (q = 0; q < Q; q++)
					  dp[q] = _mm_subs_epu8(dp[q], tempv);

				  //reset xE, xJ, and xB
				  xEv_prev = xEv = xJv = _mm_setzero_si128();
				  xBv = _mm_subs_epu8(basev, tjbmv);

			  } //else ... nothing to do, just extend another step

		  }


		  xEv = _mm_subs_epu8(xEv, tecv);
		  xJv = _mm_max_epu8(xJv,xEv);

		  xBv = _mm_max_epu8(basev, xJv);
		  xBv = _mm_subs_epu8(xBv, tjbmv);

	  } /* end loop over sequence residues 1..L */
  }
#endif

  return eslOK;


}
/*------------------ end, p7_MSVFilter_longtarget() ------------------------*/



/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
/* The benchmark driver an additional non-benchmarking options
 * to facilitate small-scale (by-eye) comparison of MSV scores against
 * other implementations, for debugging purposes.
 * The -x option compares against an emulation that should give
 * exactly the same scores. The emulation is achieved by jiggering the
 * fp scores in a generic profile to disallow gaps, have the same
 * rounding and precision as the uint8_t's MSVFilter() is using, and
 * to make the same post-hoc corrections for the NN, CC, JJ
 * contributions to the final nat score; under these contrived
 * circumstances, p7_ReferenceViterbi() gives the same scores as
 * p7_MSVFilter().
 * 
 * For using -x, you probably also want to limit the number of
 * generated target sequences, using -N10 or -N100 for example.
 */
#ifdef p7MSVFILTER_BENCHMARK
/* 
   ./msvfilter_benchmark <hmmfile>            runs benchmark 
   ./msvfilter_benchmark -N100 -x <hmmfile>   compare scores to exact emulation
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "equate scores to trusted implementation (debug)",  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for MSVFilter() implementation";

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
  P7_REFMX       *gx      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc1, sc2;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, L);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  if (esl_opt_GetBoolean(go, "-x")) p7_profile_SameAsMF(om, gm);

  ox = p7_filtermx_Create(gm->M);
  gx = p7_refmx_Create(gm->M, L);

  /* Get a baseline time: how long it takes just to generate the sequences */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_MSVFilter    (dsq, L, om, ox, &sc1);   

      /* -x option: compare generic to fast score in a way that should give exactly the same result */
      if (esl_opt_GetBoolean(go, "-x"))
	{
	  p7_ReferenceViterbi(dsq, L, gm, gx, NULL, &sc2); 
	  sc2 /= om->scale_b;
	  if (om->mode == p7_UNILOCAL)   sc2 -= 2.0; /* that's ~ L \log \frac{L}{L+2}, for our NN,CC,JJ */
	  else if (om->mode == p7_LOCAL) sc2 -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
	  printf("%.4f %.4f\n", sc1, sc2);  
	}
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_filtermx_Destroy(ox);
  p7_refmx_Destroy(gx);
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
#endif /*p7MSVFILTER_BENCHMARK*/
/*------------------ end, benchmark driver ----------------------*/




/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7MSVFILTER_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_viterbi.h"

/* utest_comparison()
 * 
 * MSV is tested against reference Viterbi, after configuring a
 * generic profile such that its scores should match MSV scores,
 * using p7_profile_SameAsMF().
 *
 * We sample a random model of length <M>, and score <N> random test
 * sequences of length <L>. 
 * 
 * Because MSV scores can overflow, wedon't sample high-scoring
 * homologs for this test. Indeed, here we assume that we don't
 * accidentally generate a high-scoring random sequence that overflows
 * MSVFilter()'s limited range -- if we do, the score comparison will
 * fail.
 */
static void
utest_comparison(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_FILTERMX *ox  = p7_filtermx_Create(M);
  P7_REFMX    *gx  = p7_refmx_Create(M, L);
  float sc1, sc2;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_profile_SameAsMF(om, gm);
#if 0
  p7_oprofile_Dump(stdout, om);                   /* dumps the optimized profile */
  p7_filtermx_SetDumpMode(ox, stdout, TRUE);      /* makes the fast DP algorithms dump their matrices */
#endif

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_MSVFilter        (dsq, L, om, ox,       &sc1);
      p7_ReferenceViterbi (dsq, L, gm, gx, NULL, &sc2);
#if 0
      p7_refmx_Dump(stdout, gx);   /* dumps a generic DP matrix */
#endif

      sc2 = sc2 / om->scale_b - 3.0f;
      if (fabs(sc1-sc2) > 0.001) esl_fatal("msv filter unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);

      p7_filtermx_Reuse(ox);
      p7_refmx_Reuse(gx);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_filtermx_Destroy(ox);
  p7_refmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7MSVFILTER_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7MSVFILTER_TESTDRIVE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the SSE MSVFilter() implementation";

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

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("MSVFilter() tests, DNA\n");
  utest_comparison(r, abc, bg, M, L, N);   /* normal sized models */
  utest_comparison(r, abc, bg, 1, L, 10);  /* size 1 models       */
  utest_comparison(r, abc, bg, M, 1, 10);  /* size 1 sequences    */

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("MSVFilter() tests, protein\n");
  utest_comparison(r, abc, bg, M, L, N);   
  utest_comparison(r, abc, bg, 1, L, 10);  
  utest_comparison(r, abc, bg, M, 1, 10);  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);

  fprintf(stderr, "#  status = ok\n");
  return eslOK;
}
#endif /*p7MSVFILTER_TESTDRIVE*/



/*****************************************************************
 * 5. Example
 *****************************************************************/

#ifdef p7MSVFILTER_EXAMPLE
/* A minimal example.
   Also useful for debugging on small HMMs and sequences.

   ./msvfilter_example <hmmfile> <seqfile>
 */ 
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "output in one line awkable format",                0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "output in profmark format",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of MSV filter algorithm";

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
  P7_FILTERMX    *ox      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           msvraw, nullsc, msvscore;
  double          P;
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
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, sq->n);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  ox = p7_filtermx_Create(gm->M);

  /* Useful to place and compile in for debugging: 
     p7_oprofile_Dump(stdout, om);                  dumps the optimized profile
     p7_filtermx_SetDumpMode(ox, stdout, TRUE);     makes the fast DP algorithms dump their matrices
     p7_oprofile_SameMSV(om, gm);
  */

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_bg_SetLength(bg,            sq->n);

      p7_filtermx_GrowTo(ox, om->M) ;

      p7_MSVFilter   (sq->dsq, sq->n, om, ox, &msvraw);  
      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);
      msvscore = (msvraw - nullsc) / eslCONST_LOG2;
      P        = esl_gumbel_surv(msvscore,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);

      if (esl_opt_GetBoolean(go, "-1"))
	{
	  printf("%-30s  %-20s  %9.2g  %7.2f\n", sq->name, hmm->name, P, msvscore);
	}
      else if (esl_opt_GetBoolean(go, "-P"))
	{ /* output suitable for direct use in profmark benchmark postprocessors: */
	  printf("%g  %.2f  %s  %s\n", P, msvscore, sq->name, hmm->name);
	}
      else
	{
	  printf("target sequence:         %s\n",        sq->name);
	  printf("msv filter raw score:    %.2f nats\n", msvraw);
	  printf("null score:              %.2f nats\n", nullsc);
	  printf("per-seq score:           %.2f bits\n", msvscore);
	  printf("P-value:                 %g\n",        P);
	  puts("\n");
	}
      
      p7_filtermx_Reuse(ox);
      esl_sq_Reuse(sq);
    }

  /* cleanup */
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_filtermx_Destroy(ox);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7MSVFILTER_EXAMPLE*/
/*---------------------- end, example ---------------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/



