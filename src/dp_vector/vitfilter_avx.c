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

#include "easel.h"
#include "esl_sse.h"

#if p7_CPU_ARCH == intel
#include <immintrin.h>
#ifdef HAVE_AVX2
  #include "esl_avx.h"
#endif
#endif /* intel arch */

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
p7_ViterbiFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
#ifdef HAVE_AVX2
//printf("Starting ViterbiFilter\n");
  int i;        /* counter over sequence positions 1..L                      */
  
  register __m256i mpv_AVX, dpv_AVX, ipv_AVX; /*  prvious row values                                       */
  register __m256i sv_AVX;       /* temp storage of 1 curr row value in progress              */
  register __m256i dcv_AVX;      /* delayed storage of D(i,q+1)                               */
  register __m256i xEv_AVX;      /* E state: keeps max for Mk->E as we go                     */
  register __m256i xBv_AVX;      /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m256i Dmaxv_AVX;          /* keeps track of maximum D cell on row                      */
  int16_t  xE_AVX, xB_AVX, xC_AVX, xJ_AVX, xN_AVX;     /* special states' scores                                    */
  int16_t  Dmax_AVX;       /* maximum D cell score on row                               */
  int q_AVX;         /* counter over vectors 0..nq-1                              */
  int Q_AVX        = P7_NVW_AVX(om->M);    /* segment length: # of vectors                              */
  __m256i *dp_AVX;
  __m256i *rsc_AVX;        /* will point at om->ru[x] for residue x[i]                  */
  __m256i *tsc_AVX;        /* will point into (and step thru) om->tu                    */
  __m256i negInfv_AVX;

  int     status;

  /* Contract checks */
  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  /* note however that ViterbiFilter numerics are only guaranteed for local alignment, not glocal        */


  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of Vit filter matrix failed");

  dp_AVX = ox->dp_AVX;           /* enables MMX_AVXf(), IMX_AVXf(), DMX_AVXf() access macros. Must be set AFTER the GrowTo, because ox->dp may get moved */

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_VITFILTER;

  /* -infinity is -32768 */
  // Want -32768 in the low 16 bits of negInfv_AVX
  negInfv_AVX = _mm256_setzero_si256();
  negInfv_AVX = _mm256_insert_epi16(negInfv_AVX, -32768, 0);
  /* Initialization. In unsigned arithmetic, -infinity is -32768
   */
  for (q_AVX = 0; q_AVX < Q_AVX; q_AVX++)
    MMX_AVXf(q_AVX) = IMX_AVXf(q_AVX) = DMX_AVXf(q_AVX) = _mm256_set1_epi16(-32768);
  xN_AVX   = om->base_w;
  xB_AVX   = xN_AVX + om->xw[p7O_N][p7O_MOVE];
  xJ_AVX   = -32768;
  xC_AVX   = -32768;
  xE_AVX   = -32768;

#ifdef p7_DEBUGGING
  if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, 0, xE, 0, xJ, xB, xC); /* first 0 is <rowi>: do header. second 0 is xN: always 0 here. */
#endif

  for (i = 1; i <= L; i++)
    {
      rsc_AVX   = om->rwv_AVX[dsq[i]];
      tsc_AVX   = om->twv_AVX;
      dcv_AVX   = _mm256_set1_epi16(-32768);      /* "-infinity" */
      xEv_AVX   = _mm256_set1_epi16(-32768);     
      Dmaxv_AVX = _mm256_set1_epi16(-32768);     
      xBv_AVX   = _mm256_set1_epi16(xB_AVX);

      /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically; replace it with -32768.
       */
    
      mpv_AVX = MMX_AVXf(Q_AVX-1);  
      mpv_AVX = esl_avx_leftshift_two(mpv_AVX);
      mpv_AVX = _mm256_or_si256(mpv_AVX, negInfv_AVX);
      
      dpv_AVX = DMX_AVXf(Q_AVX-1);  
       dpv_AVX = esl_avx_leftshift_two(dpv_AVX);
       dpv_AVX = _mm256_or_si256(dpv_AVX, negInfv_AVX);
      
      ipv_AVX = IMX_AVXf(Q_AVX-1);  
      ipv_AVX = esl_avx_leftshift_two(ipv_AVX);
       ipv_AVX = _mm256_or_si256(ipv_AVX, negInfv_AVX);

      for (q_AVX = 0; q_AVX < Q_AVX; q_AVX++)
      {
    
        /* Calculate new MMXf(i,q); don't store it yet, hold it in sv. */
        sv_AVX   =                    _mm256_adds_epi16(xBv_AVX, *tsc_AVX);  tsc_AVX++;
        sv_AVX   = _mm256_max_epi16 (sv_AVX, _mm256_adds_epi16(mpv_AVX, *tsc_AVX)); tsc_AVX++;
        sv_AVX   = _mm256_max_epi16 (sv_AVX, _mm256_adds_epi16(ipv_AVX, *tsc_AVX)); tsc_AVX++;
        sv_AVX   = _mm256_max_epi16 (sv_AVX, _mm256_adds_epi16(dpv_AVX, *tsc_AVX)); tsc_AVX++;
        sv_AVX   = _mm256_adds_epi16(sv_AVX, *rsc_AVX);                      rsc_AVX++;
        xEv_AVX  = _mm256_max_epi16(xEv_AVX, sv_AVX);

        /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
         * {MDI}MX(q) is then the current, not the prev row
         */
        mpv_AVX = MMX_AVXf(q_AVX);
        dpv_AVX = DMX_AVXf(q_AVX);
        ipv_AVX = IMX_AVXf(q_AVX);

        /* Do the delayed stores of {MD}(i,q) now that memory is usable */
        MMX_AVXf(q_AVX) = sv_AVX;
        DMX_AVXf(q_AVX) = dcv_AVX;

        /* Calculate the next D(i,q+1) partially: M->D only;
               * delay storage, holding it in dcv
         */
        dcv_AVX   = _mm256_adds_epi16(sv_AVX, *tsc_AVX);  tsc_AVX++;
        Dmaxv_AVX = _mm256_max_epi16(dcv_AVX, Dmaxv_AVX);

        /* Calculate and store I(i,q) */
        sv_AVX     =                    _mm256_adds_epi16(mpv_AVX, *tsc_AVX);  tsc_AVX++;
        IMX_AVXf(q_AVX)= _mm256_max_epi16 (sv_AVX, _mm256_adds_epi16(ipv_AVX, *tsc_AVX)); tsc_AVX++;
      }

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      //Now, find the horizontal max of the 16-bit values in xEv_AVX. 
      xE_AVX = esl_avx_hmax_epi16(xEv_AVX);

      if (xE_AVX >= 32767) {*ret_sc = eslINFINITY; return eslERANGE; } /* immediately detect overflow */
      xN_AVX = xN_AVX + om->xw[p7O_N][p7O_LOOP];
      xC_AVX = ESL_MAX(xC_AVX + om->xw[p7O_C][p7O_LOOP], xE_AVX + om->xw[p7O_E][p7O_MOVE]);
      xJ_AVX = ESL_MAX(xJ_AVX + om->xw[p7O_J][p7O_LOOP], xE_AVX + om->xw[p7O_E][p7O_LOOP]);
      xB_AVX = ESL_MAX(xJ_AVX + om->xw[p7O_J][p7O_MOVE], xN_AVX + om->xw[p7O_N][p7O_MOVE]);
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
      
      // compute Dmax_AVX = horizontal 16-bit max of Dmaxv_AVX
      Dmax_AVX = esl_avx_hmax_epi16(Dmaxv_AVX);

      if (Dmax_AVX + om->ddbound_w > xB_AVX) 
  {

  
//printf("vitfilter AVX doing DD calc, Dmax_AVX = %i\n", Dmax_AVX);
    /* Now we're obligated to do at least one complete DD path to be sure. */
    /* dcv has carried through from end of q loop above */
      dcv_AVX = esl_avx_leftshift_two(dcv_AVX);
 
    dcv_AVX = _mm256_or_si256(dcv_AVX, negInfv_AVX);
    tsc_AVX = om->twv_AVX + 7*Q_AVX;  /* set tsc to start of the DD's */
    for (q_AVX = 0; q_AVX < Q_AVX; q_AVX++) 
      {
        DMX_AVXf(q_AVX) = _mm256_max_epi16(dcv_AVX, DMX_AVXf(q_AVX));  
        dcv_AVX     = _mm256_adds_epi16(DMX_AVXf(q_AVX), *tsc_AVX); tsc_AVX++;
      }

    /* We may have to do up to three more passes; the check
     * is for whether crossing a segment boundary can improve
     * our score. 
     */
    do {
         dcv_AVX = esl_avx_leftshift_two(dcv_AVX);
 
      dcv_AVX = _mm256_or_si256(dcv_AVX, negInfv_AVX);
      tsc_AVX = om->twv_AVX + 7*Q_AVX;  /* set tsc to start of the DD's */
      for (q_AVX = 0; q_AVX < Q_AVX; q_AVX++) 
        {
            __m256i check1 = _mm256_cmpgt_epi16(dcv_AVX, DMX_AVXf(q_AVX));
            int cont = _mm256_movemask_epi8(check1); // Continue if any of the 16-bit field s in dcv are greater than 
            // the one we compared to
            if (cont == 0) break;
            DMX_AVXf(q_AVX) = _mm256_max_epi16(dcv_AVX, DMX_AVXf(q_AVX));  
            dcv_AVX   = _mm256_adds_epi16(DMX_AVXf(q_AVX), *tsc_AVX); tsc_AVX++;
         }     
    } while (q_AVX == Q_AVX);
  }
      else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
  {
//printf("vitfilter AVX skipping DD calc Dmax_AVX = %i\n", Dmax_AVX);
  dcv_AVX = esl_avx_leftshift_two(dcv_AVX);
   
    DMX_AVXf(0) = _mm256_or_si256(dcv_AVX, negInfv_AVX);
  }

#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, i, xE, 0, xJ, xB, xC);   
#endif
    } /* end loop over sequence residues 1..L */


  /* finally C->T */
  if (xC_AVX > -32768)
    {
      *ret_sc = (float) xC_AVX + (float) om->xw[p7O_C][p7O_MOVE] - (float) om->base_w;
      /* *ret_sc += L * om->ncj_roundoff;  see J4/150 for rationale: superceded by -3.0nat approximation*/
      *ret_sc /= om->scale_w;
      *ret_sc -= 3.0; /* the NN/CC/JJ=0,-3nat approximation: see J5/36. That's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ contrib */
    }
  else  *ret_sc = -eslINFINITY;
#endif // HAVE_AVX2. Leave the return statement so that there's a function to link if we don't have AVX2 support
  return eslOK;
}


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

