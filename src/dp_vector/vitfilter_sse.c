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
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_gumbel.h"

#include "search/p7_pipeline.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"

#ifdef eslENABLE_SSE4
#include <x86intrin.h>
#include "esl_sse.h"




/*****************************************************************
 * 1. Viterbi filter implementation.
 *****************************************************************/

/* Function:  p7_ViterbiFilter_sse()
 * Synopsis:  Calculates Viterbi score, vewy vewy fast, in limited precision.
 * See:       vitfilter.c::p7_ViterbiFilter()
 */
int
p7_ViterbiFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
  int i;                           /* counter over sequence positions 1..L                      */
  register __m128i mpv, dpv, ipv;  /* previous row values                                       */
  register __m128i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128i dcv;		   /* delayed storage of D(i,q+1)                               */
  register __m128i xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i Dmaxv;          /* keeps track of maximum D cell on row                      */
  int16_t  xE, xB, xC, xJ, xN;	   /* special states' scores                                    */
  int16_t  Dmax;		   /* maximum D cell score on row                               */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q = P7_Q(om->M, p7_VWIDTH_SSE / sizeof(int16_t));  /* segment length: # of vectors        */
  __m128i *dp;
  __m128i *rsc;			   /* will point at om->ru[x] for residue x[i]                  */
  __m128i *tsc;			   /* will point into (and step thru) om->tu                    */
  __m128i neginfmask = _mm_insert_epi16( _mm_setzero_si128(), -32768, 0);
  int     status;

  /* Contract checks */
  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  /* note however that ViterbiFilter numerics are only guaranteed for local alignment, not glocal        */

  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_Reinit(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of Vit filter matrix failed");
  dp = (__m128i *) ox->dp;       /* enables MMXf(), IMXf(), DMXf() access macros. Must be set AFTER the Reinit, because ox->dp may get moved */

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->Vw   = p7_VWIDTH_SSE / sizeof(int16_t);
  ox->type = p7F_VITFILTER;
  ESL_DASSERT1(( ox->Vw == om->V / sizeof(int16_t)));

  /* Initialization. In int16_t, our -infinity is -32768  */
  for (q = 0; q < Q; q++)
    MMXf(q) = IMXf(q) = DMXf(q) = _mm_set1_epi16(-32768);
  xN   = om->base_w;
  xB   = xN + om->xw[p7O_N][p7O_MOVE];
  xJ   = -32768;
  xC   = -32768;
  xE   = -32768;

#if eslDEBUGLEVEL > 0
  if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, 0, xE, 0, xJ, xB, xC); /* first 0 is <rowi>: do header. second 0 is xN: always 0 here. */
#endif

  for (i = 1; i <= L; i++)
    {
      rsc   = (__m128i *) om->rwv[dsq[i]];
      tsc   = (__m128i *) om->twv;
      dcv   = _mm_set1_epi16(-32768);      /* "-infinity" */
      xEv   = _mm_set1_epi16(-32768);     
      Dmaxv = _mm_set1_epi16(-32768);     
      xBv   = _mm_set1_epi16(xB);

      /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12. */
      mpv = esl_sse_rightshift_int16(MMXf(Q-1), neginfmask);
      dpv = esl_sse_rightshift_int16(DMXf(Q-1), neginfmask);
      ipv = esl_sse_rightshift_int16(IMXf(Q-1), neginfmask);

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
	  /* Now we're obligated to do at least one complete DD path to be sure. */
	  /* dcv has carried through from end of q loop above */
          dcv = esl_sse_rightshift_int16(dcv, neginfmask);
	  tsc = (__m128i *) om->twv + 7*Q;	/* set tsc to start of the DD's */
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
            dcv = esl_sse_rightshift_int16(dcv, neginfmask);
	    tsc = (__m128i *) om->twv + 7*Q;	/* set tsc to start of the DD's */
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
	  DMXf(0) = esl_sse_rightshift_int16(dcv, neginfmask);
	}

#if eslDEBUGLEVEL > 0
      if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, i, xE, 0, xJ, xB, xC);   
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T */
  if (xC > -32768)
    {
      *ret_sc = (float) xC + (float) om->xw[p7O_C][p7O_MOVE] - (float) om->base_w;
      *ret_sc /= om->scale_w;
      *ret_sc -= 3.0; /* the NN/CC/JJ=0,-3nat approximation: see J5/36. That's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ contrib */
    }
  else  *ret_sc = -eslINFINITY;
  return eslOK;
}

#else // ! eslENABLE_SSE4
/* provide a callable function even when we're `./configure --disable-sse` */
int
p7_ViterbiFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(om); ESL_UNUSED(ox); ESL_UNUSED(ret_sc);  // shut up, compiler, I know what I'm doing
  esl_fatal("SSE support was not enabled at compile time. Can't use p7_ViterbiFilter_sse().");
  return eslFAIL; // NOTREACHED
}
#if defined p7VITFILTER_SSE_TESTDRIVE || p7VITFILTER_SSE_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_SSE4 or not




