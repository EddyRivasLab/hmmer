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
 * 
 * SRE, Sun Nov 25 11:26:48 2007 [Casa de Gatos]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"
#include "esl_gumbel.h"

#include "hmmer.h"
#include "impl_sse.h"

/*****************************************************************
 * 1. The p7_MSVFilter() DP implementation.
 *****************************************************************/
 
/* Function:  p7_MSVFilter()
 * Synopsis:  Calculates MSV score, vewy vewy fast, in limited precision.
 * Incept:    SRE, Wed Dec 26 15:12:25 2007 [Janelia]
 *
 * Purpose:   Calculates an approximation of the MSV score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and a preallocated one-row DP matrix <ox>. Return the 
 *            estimated MSV score (in nats) in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow.
 *            
 *            The model may be in any mode, because only its match
 *            emission scores will be used. The MSV filter inherently
 *            assumes a multihit local mode, and uses its own special
 *            state transition scores, not the scores in the profile.
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: MSV score (in nats)          
 *                      
 * Note:      We misuse the matrix <ox> here, using only a third of the
 *            first dp row, accessing it as <dp[0..Q-1]> rather than
 *            in triplets via <{MDI}MX(q)> macros, since we only need
 *            to store M state values. We know that if <ox> was big
 *            enough for normal DP calculations, it must be big enough
 *            to hold the MSVFilter calculation.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the score overflows the limited range; in
 *            this case, this is a high-scoring hit.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_MSVFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m128i mpv;            /* previous row values                                       */
  register __m128i xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128i biasv;	   /* emission bias in a vector                                 */
  uint8_t  xJ;                     /* special states' scores                                    */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQB(om->M);   /* segment length: # of vectors                              */
  __m128i *dp  = ox->dpb[0];	   /* we're going to use dp[0][0..q..Q-1], not {MDI}MX(q) macros*/
  __m128i *rsc;			   /* will point at om->rbv[x] for residue x[i]                 */

  __m128i xJv;                     /* vector for states score                                   */
  __m128i tjbmv;                   /* vector for cost of moving from either J or N through B to an M state */
  __m128i tecv;                    /* vector for E->C  cost                                     */
  __m128i basev;                   /* offset for scores                                         */
  __m128i ceilingv;                /* saturateed simd value used to test for overflow           */
  __m128i tempv;                   /* work vector                                               */

  int cmp;
  int status = eslOK;

  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ16)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M   = om->M;

  /* Try highly optimized ssv filter first */
  status = p7_SSVFilter(dsq, L, om, ret_sc);
  if (status != eslENORESULT) return status;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.
   */
  biasv = _mm_set1_epi8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */
  for (q = 0; q < Q; q++) dp[q] = _mm_setzero_si128();
  xJ   = 0;

  /* saturate simd register for overflow test */
  ceilingv = _mm_cmpeq_epi8(biasv, biasv);
  basev = _mm_set1_epi8((int8_t) om->base_b);

  tjbmv = _mm_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  tecv = _mm_set1_epi8((int8_t) om->tec_b);

  xJv = _mm_subs_epu8(biasv, biasv);
  xBv = _mm_subs_epu8(basev, tjbmv);

#if p7_DEBUGGING
  if (ox->debugging)
    { 
      uint8_t xB;
      xB = _mm_extract_epi16(xBv, 0);
      xJ = _mm_extract_epi16(xJv, 0);
      p7_omx_DumpMFRow(ox, 0, 0, 0, xJ, xB, xJ);
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
      cmp = _mm_movemask_epi8(tempv);

      /* Now the "special" states, which start from Mk->E (->C, ->J->B)
       * Use shuffles instead of shifts so when the last max has completed,
       * the last four elements of the simd register will contain the
       * max value.  Then the last shuffle will broadcast the max value
       * to all simd elements.
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

      /* immediately detect overflow */
      if (cmp != 0x0000)
	{ 
	  *ret_sc = eslINFINITY; 
	  return eslERANGE; 
	}

      xEv = _mm_subs_epu8(xEv, tecv);
      xJv = _mm_max_epu8(xJv,xEv);
      
      xBv = _mm_max_epu8(basev, xJv);
      xBv = _mm_subs_epu8(xBv, tjbmv);
	  
#if p7_DEBUGGING
      if (ox->debugging)
	{
	  uint8_t xB, xE;
	  xB = _mm_extract_epi16(xBv, 0);
	  xE = _mm_extract_epi16(xEv, 0);
	  xJ = _mm_extract_epi16(xJv, 0);
	  p7_omx_DumpMFRow(ox, i, xE, 0, xJ, xB, xJ);   
	}
#endif
    } /* end loop over sequence residues 1..L */

  xJ = (uint8_t) _mm_extract_epi16(xJv, 0);

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  *ret_sc = ((float) (xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */

  return eslOK;
}
/*------------------ end, p7_MSVFilter() ------------------------*/



/* Function:  p7_SSVFilter_longtarget()
 * Synopsis:  Finds windows with SSV scores above some threshold (vewy vewy fast, in limited precision)
 * Incept:    TJW, Mon Feb 22 16:27:25 2010 [Janelia] - based on p7_MSVFilter
 *
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
 *            sc_threshv - a precomputed vector containing the score threshold required for a window to be considered "high-scoring"
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
int
p7_SSVFilter_longtarget(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const FM_HMMDATA *hmmdata, uint8_t sc_thresh, __m128i sc_threshv, FM_WINDOWLIST *windowlist)
{

  register __m128i mpv;            /* previous row values                                       */
  register __m128i xEv;		   /* E state: keeps max for Mk->E for a single iteration       */
  register __m128i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128i biasv;	   /* emission bias in a vector                                 */
  uint8_t  xJ;                     /* special states' scores                                    */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQB(om->M);   /* segment length: # of vectors                              */
  __m128i *dp  = ox->dpb[0];	   /* we're going to use dp[0][0..q..Q-1], not {MDI}MX(q) macros*/
  __m128i *rsc;			   /* will point at om->rbv[x] for residue x[i]                 */
  __m128i tecv;                    /* vector for E->C  cost                                     */
  __m128i tjbmv;                   /* vector for J->B move cost + B->M move costs               */
  __m128i basev;                   /* offset for scores                                         */
  __m128i ceilingv;                /* saturated simd value used to test for overflow           */
  __m128i tempv;                   /* work vector                                               */
  int cmp;
  int status;
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
  uint8_t *scores = NULL;

  ESL_ALLOC(scores, (1+16*Q) * sizeof(uint8_t) );

  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ16)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M   = om->M;


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
          for (k = 0; k < 16; k++) scores[q+Q*k+1] = u.b[k]; // unstripe
          dp[q] = _mm_set1_epi8(0); // while we're here ... this will cause values to get reset to xB in next dp iteration
	    }
	    for (k = 1; k <= om->M; k++) {
          if (scores[k] >= sc_thresh && scores[k] > rem_sc) {
            end = k;
            rem_sc = scores[k];
          }
	    }
	    //recover the diagonal that hit threshold
	    start = end;
	    target_end = target_start = i;
	    sc = rem_sc;
	    while (rem_sc > om->base_b - om->tjb_b - om->tbm_b) {
	      //printf("%d (%d)\n", rem_sc, om->bias_b -  hmmdata->s.scores_b[start][dsq[target_start]] );
	      rem_sc -= om->bias_b -  hmmdata->s.scores_b[start][dsq[target_start]];
	      --start;
	      --target_start;
	      //if ( start == 0 || target_start==0)    break;
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
	      sc += om->bias_b -  hmmdata->s.scores_b[start][dsq[n]];
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

      fm_newWindow(windowlist, 0, target_start, k, end, end-start+1 , ret_sc, fm_nocomplement );

      i = target_end; // skip forward
	  }


  } /* end loop over sequence residues 1..L */


  free(scores);
  return eslOK;

ERROR:
  if (scores != NULL) free(scores);
  return eslEMEM;
}
/*------------------ end, p7_SSVFilter_longtarget() ------------------------*/



/* Function:  p7_MSVFilter_longtarget()
 * Synopsis:  Finds windows with MSV scores above some threshold (vewy vewy fast, in limited precision)
 * Incept:    TJW, Mon Feb 22 16:27:25 2010 [Janelia] - based on p7_MSVFilter
 *
 * Purpose:   Calculates an approximation of the MSV score for regions of
 *            sequence <dsq> of length <L> residues, using optimized profile
 *            <om>, and a preallocated one-row DP matrix <ox>, and captures
 *            the positions at which such regions exceed the score required
 *            to be significant in the eyes of the calling function (usually
 *            p=0.02).  Note that this will typically simply call
 *            p7_SSVFilter_longtarget(), as the score threshold is nearly
 *            always lower than the score required for a pass through the
 *            J-state to yield a positive score.
 *
 *            Rather than simply capturing positions at which a score threshold
 *            is reached, this function (and SSVFilter_longtarget) establishes
 *            windows around those high-scoring positions, and merges overlapping
 *            windows.
 *
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues
 *            om      - optimized profile
 *            ox      - DP matrix
 *            bg      - the background model, required for translating a P-value threshold into a score threshold
 *            P       - p-value below which a region is captured as being above threshold
 *            starts  - RETURN: array of start positions for windows surrounding above-threshold areas
 *            ends    - RETURN: array of end positions for windows surrounding above-threshold areas
 *			  hit_cnt - RETURN: count of entries in the above two arrays
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
p7_MSVFilter_longtarget(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const FM_HMMDATA *hmmdata, P7_BG *bg, double P, FM_WINDOWLIST *windowlist, int do_biasfilter)
{


  register __m128i mpv;        /* previous row values                                       */
  register __m128i xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128i biasv;	   /* emission bias in a vector                                 */
  uint8_t  xJ;                 /* special states' scores                                    */
  int i;			           /* counter over sequence positions 1..L                      */
  int j;
  int q;			           /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQB(om->M);   /* segment length: # of vectors                              */
  __m128i *dp  = ox->dpb[0];	   /* we're going to use dp[0][0..q..Q-1], not {MDI}MX(q) macros*/
  __m128i *rsc;			        /* will point at om->rbv[x] for residue x[i]                 */

  __m128i xEv_prev;            	 /* E state: keeps the max seen up to the current position i, for window overlap test.  Note that this will always just be xJv + tecv    */
  __m128i xJv;                     /* vector for states score                                   */
  __m128i tecv;                    /* vector for E->C  cost                                     */
  __m128i tjbmv;                    /* vector for [JN]->B->M move cost                                  */
  __m128i basev;                   /* offset for scores                                         */
  __m128i sc_threshv;               /* pushes value to saturation if it's above pthresh  */
  __m128i jthreshv;                /* pushes value to saturation if it's above pthresh  */
  __m128i ceilingv;                /* saturateed simd value used to test for overflow           */
  __m128i tempv;                   /* work vector                                               */

  FM_WINDOW *prev_window;
  FM_WINDOW *curr_window;
  int window_start;
  int window_end;

  float nullsc;
  float bias_sc;
  float biasP;
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
   *  Here, I compute threshold with length model based on max_length.  Usually, the
   *  length of a window returned by this scan will be 2*max_length-1 or longer.  Doesn't
   *  really matter - in any case, both the bg and om models will change with roughly
   *  1 bit for each doubling of the length model, so they offset.
   */
  float invP = esl_gumbel_invsurv(P, om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  p7_bg_SetLength(bg, om->max_length);
  p7_oprofile_ReconfigMSVLength(om, om->max_length);
  p7_bg_NullOne  (bg, dsq, om->max_length, &nullsc);

  uint8_t sc_thresh = (int) ceil( ( ( nullsc  + (invP * eslCONST_LOG2) + 3.0 )  * om->scale_b ) + om->base_b +  om->tec_b  + om->tjb_b );
  sc_threshv = _mm_set1_epi8((int8_t) 255 - sc_thresh);

  uint8_t jthresh = om->base_b + om->tjb_b + om->tec_b + 1;  //+1 because the score of a pass through the model must be positive to contribute to MSV score
  jthreshv = _mm_set1_epi8((int8_t) 255 - (int)jthresh);


  if ( sc_thresh < jthresh) {

    p7_SSVFilter_longtarget(dsq, L, om, ox, hmmdata, (uint8_t)sc_thresh, sc_threshv, windowlist);

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


	  /* Check that the DP matrix is ok for us. */
	  if (Q > ox->allocQ16)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
	  ox->M   = om->M;

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


	#if p7_DEBUGGING
	  if (ox->debugging)
		{
		  uint8_t xB;
		  xB = _mm_extract_epi16(xBv, 0);
		  xJ = _mm_extract_epi16(xJv, 0);
		  p7_omx_DumpMFRow(ox, 0, 0, 0, xJ, xB, xJ);
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


  //filter for biased composition here
  if ( windowlist->count > 0 ) {

    if (do_biasfilter) {
      j = 0;
      for (i=0; i<windowlist->count; i++) {
        curr_window = windowlist->windows+i;

        p7_bg_FilterScore(bg, dsq+curr_window->n, curr_window->length, &bias_sc);
        bias_sc = (curr_window->score - bias_sc) / eslCONST_LOG2;
        biasP = esl_gumbel_surv(bias_sc,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);

        if (biasP <= P ) { // keep it
          windowlist->windows[j] = windowlist->windows[i];
          j++;
        }
      }
      windowlist->count = j;
    }

    //widen window

    for (i=0; i<windowlist->count; i++) {
      curr_window = windowlist->windows+i;

//      printf ("%d -> %d  ||| ", curr_window->n, curr_window->n + curr_window->length-1);

//      window_start = ESL_MAX( 1,   curr_window->n + curr_window->length - om->max_length) ;
//      window_end   = ESL_MIN( L ,  curr_window->n + curr_window->length + om->max_length - 2);

      window_start = ESL_MAX( 1,   curr_window->n - om->max_length * hmmdata->prefix_lengths[curr_window->k - curr_window->length + 1] + 1 - 10)  ;
      window_end   = ESL_MIN( L ,  curr_window->n + curr_window->length + om->max_length * hmmdata->suffix_lengths[curr_window->k] - 1 + 10)  ;

//      window_start = ESL_MAX( 1,   ( curr_window->n + curr_window->length )  - (om->max_length * hmmdata->prefix_lengths[curr_window->k ]) )  ;
//      window_end   = ESL_MIN( L ,  ( curr_window->n )                        + (om->max_length * hmmdata->suffix_lengths[curr_window->k - curr_window->length + 1]))   ;



//printf ("%d -> %d\n", window_start, window_end);
      curr_window->n = window_start;
      curr_window->length = window_end - window_start + 1;
    }


    // merge overlapping windows, compressing list in place.
    int new_hit_cnt = 0;
    for (i=1; i<windowlist->count; i++) {
      prev_window = windowlist->windows+new_hit_cnt;
      curr_window = windowlist->windows+i;
      if (curr_window->n <= prev_window->n + prev_window->length ) {
        //merge windows
        if (  curr_window->n + curr_window->length >  prev_window->n + prev_window->length )
          prev_window->length = curr_window->n + curr_window->length - prev_window->n;

      } else {
        new_hit_cnt++;
        windowlist->windows[new_hit_cnt] = windowlist->windows[i];
      }
    }
    windowlist->count = new_hit_cnt+1;

	  if ( windowlist->windows[0].n  <  1)
	    windowlist->windows[0].n =  1;

	  if ( windowlist->windows[windowlist->count].n + windowlist->windows[windowlist->count].length - 1  >  L)
	    windowlist->windows[windowlist->count].length =  L - windowlist->windows[windowlist->count].n + 1;

  }
  return eslOK;


}
/*------------------ end, p7_MSVFilter_longtarget() ------------------------*/



/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
/* The benchmark driver has some additional non-benchmarking options
 * to facilitate small-scale (by-eye) comparison of MSV scores against
 * other implementations, for debugging purposes.
 * 
 * The -c option compares against p7_GMSV() scores. This allows
 * measuring the error inherent in the SSE implementation's reduced
 * precision (p7_MSVFilter() runs in uint8_t; p7_GMSV() uses floats).
 * 
 * The -x option compares against an emulation that should give
 * exactly the same scores. The emulation is achieved by jiggering the
 * fp scores in a generic profile to disallow gaps, have the same
 * rounding and precision as the uint8_t's MSVFilter() is using, and
 * to make the same post-hoc corrections for the NN, CC, JJ
 * contributions to the final nat score; under these contrived
 * circumstances, p7_GViterbi() gives the same scores as
 * p7_MSVFilter().
 * 
 * For using either -c or -x, you probably also want to limit the
 * number of generated target sequences, using -N10 or -N100 for
 * example.
 */
#ifdef p7MSVFILTER_BENCHMARK
/* 
   gcc -o msvfilter-benchmark -std=gnu99 -g -Wall -msse2 -I.. -L.. -I../../easel -L../../easel -Dp7MSVFILTER_BENCHMARK msvfilter.c -lhmmer -leasel -lm 
   icc -o msvfilter-benchmark -O3 -static -I.. -L.. -I../../easel -L../../easel -Dp7MSVFILTER_BENCHMARK msvfilter.c -lhmmer -leasel -lm 

   ./benchmark-msvfilter <hmmfile>            runs benchmark 
   ./benchmark-msvfilter -N100 -c <hmmfile>   compare scores to generic impl
   ./benchmark-msvfilter -N100 -x <hmmfile>   compare scores to exact emulation
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-x", "compare scores to generic implementation (debug)", 0 }, 
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-c", "equate scores to trusted implementation (debug)",  0 },
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
  P7_OMX         *ox      = NULL;
  P7_GMX         *gx      = NULL;
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
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  if (esl_opt_GetBoolean(go, "-x")) p7_profile_SameAsMF(om, gm);

  ox = p7_omx_Create(gm->M, 0, 0);
  gx = p7_gmx_Create(gm->M, L);

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

      /* -c option: compare generic to fast score */
      if (esl_opt_GetBoolean(go, "-c")) 
	{
	  p7_GMSV    (dsq, L, gm, gx, 2.0, &sc2); 
	  printf("%.4f %.4f\n", sc1, sc2);  
	}

      /* -x option: compare generic to fast score in a way that should give exactly the same result */
      if (esl_opt_GetBoolean(go, "-x"))
	{
	  p7_GViterbi(dsq, L, gm, gx, &sc2); 
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
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
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

/* 
 * We can check that scores are identical (within machine error) to
 * scores of generic DP with scores rounded the same way.  Do this for
 * a random model of length <M>, for <N> test sequences of length <L>.
 * 
 * We assume that we don't accidentally generate a high-scoring random
 * sequence that overflows MSVFilter()'s limited range.
 * 
 */
static void
utest_msv_filter(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *ox  = p7_omx_Create(M, 0, 0);
  P7_GMX      *gx  = p7_gmx_Create(M, L);
  float sc1, sc2;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_profile_SameAsMF(om, gm);
#if 0
  p7_oprofile_Dump(stdout, om);              /* dumps the optimized profile */
  p7_omx_SetDumpMode(stdout, ox, TRUE);      /* makes the fast DP algorithms dump their matrices */
#endif

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_MSVFilter(dsq, L, om, ox, &sc1);
      p7_GViterbi (dsq, L, gm, gx, &sc2);
#if 0
      p7_gmx_Dump(stdout, gx, p7_DEFAULT);   /* dumps a generic DP matrix */
#endif

      sc2 = sc2 / om->scale_b - 3.0f;
      if (fabs(sc1-sc2) > 0.001) esl_fatal("msv filter unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7MSVFILTER_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7MSVFILTER_TESTDRIVE
/* 
   gcc -g -Wall -msse2 -std=gnu99 -I.. -L.. -I../../easel -L../../easel -o msvfilter_utest -Dp7MSVFILTER_TESTDRIVE msvfilter.c -lhmmer -leasel -lm
   ./msvfilter_utest
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "impl_sse.h"

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

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("MSVFilter() tests, DNA\n");
  utest_msv_filter(r, abc, bg, M, L, N);   /* normal sized models */
  utest_msv_filter(r, abc, bg, 1, L, 10);  /* size 1 models       */
  utest_msv_filter(r, abc, bg, M, 1, 10);  /* size 1 sequences    */

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("MSVFilter() tests, protein\n");
  utest_msv_filter(r, abc, bg, M, L, N);   
  utest_msv_filter(r, abc, bg, 1, L, 10);  
  utest_msv_filter(r, abc, bg, M, 1, 10);  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*VITFILTER_TESTDRIVE*/



/*****************************************************************
 * 5. Example
 *****************************************************************/

#ifdef p7MSVFILTER_EXAMPLE
/* A minimal example.
   Also useful for debugging on small HMMs and sequences.

   gcc -g -Wall -msse2 -std=gnu99 -I.. -L.. -I../../easel -L../../easel -o msvfilter_example -Dp7MSVFILTER_EXAMPLE msvfilter.c -lhmmer -leasel -lm
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
#include "impl_sse.h"

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
  P7_OMX         *ox      = NULL;
  P7_GMX         *gx      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           msvraw, nullsc, msvscore;
  float           graw, gscore;
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
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  /* allocate DP matrices, both a generic and an optimized one */
  ox = p7_omx_Create(gm->M, 0, 0); /* one row version */
  gx = p7_gmx_Create(gm->M, sq->n);

  /* Useful to place and compile in for debugging: 
     p7_oprofile_Dump(stdout, om);              dumps the optimized profile
     p7_omx_SetDumpMode(stdout, ox, TRUE);      makes the fast DP algorithms dump their matrices
     p7_gmx_Dump(stdout, gx, p7_DEFAULT);       dumps a generic DP matrix
     p7_oprofile_SameMSV(om, gm);
  */
  //p7_oprofile_Dump(stdout, om);
  //p7_omx_SetDumpMode(stdout, ox, TRUE);    

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_ReconfigLength(gm,          sq->n);
      p7_bg_SetLength(bg,            sq->n);
      p7_omx_GrowTo(ox, om->M, 0,    sq->n); 
      p7_gmx_GrowTo(gx, gm->M,       sq->n); 

      p7_MSVFilter   (sq->dsq, sq->n, om, ox, &msvraw);  
      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);
      msvscore = (msvraw - nullsc) / eslCONST_LOG2;
      P        = esl_gumbel_surv(msvscore,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);

      p7_GMSV(sq->dsq, sq->n, gm, gx, 2.0, &graw);
      gscore   = (graw - nullsc) / eslCONST_LOG2;
      gP       = esl_gumbel_surv(gscore,  gm->evparam[p7_MMU],  gm->evparam[p7_MLAMBDA]);

      if (esl_opt_GetBoolean(go, "-1"))
	{
	  printf("%-30s\t%-20s\t%9.2g\t%7.2f\t%9.2g\t%7.2f\n", sq->name, hmm->name, P, msvscore, gP, gscore);
	}
      else if (esl_opt_GetBoolean(go, "-P"))
	{ /* output suitable for direct use in profmark benchmark postprocessors: */
	  printf("%g\t%.2f\t%s\t%s\n", P, msvscore, sq->name, hmm->name);
	}
      else
	{
	  printf("target sequence:      %s\n",        sq->name);
	  printf("msv filter raw score: %.2f nats\n", msvraw);
	  printf("null score:           %.2f nats\n", nullsc);
	  printf("per-seq score:        %.2f bits\n", msvscore);
	  printf("P-value:              %g\n",        P);
	  printf("GMSV raw score:       %.2f nats\n", graw);
	  printf("GSMV per-seq score:   %.2f bits\n", gscore);
	  printf("GSMV P-value:         %g\n",        gP);
	}
      
      esl_sq_Reuse(sq);
    }

  /* cleanup */
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_omx_Destroy(ox);
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
#endif /*p7MSVFILTER_EXAMPLE*/
/*---------------------- end, example ---------------------------*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/



