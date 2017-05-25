/* Viterbi filter implementation, ARM NEON version.
 * See vitfilter_sse.c for more detailed commentary.
 */
#include "p7_config.h"
#ifdef eslENABLE_NEON

#include <stdio.h>
#include <math.h>

#include <arm_neon.h>

#include "easel.h"
#include "esl_gumbel.h"
#include "esl_neon.h"

#include "search/p7_pipeline.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"

/*****************************************************************
 * 1. Viterbi filter implementation.
 *****************************************************************/

/* Function:  p7_ViterbiFilter_neon()
 * Synopsis:  Calculates Viterbi score, vewy vewy fast, in limited precision.
 * See:       vitfilter.c::p7_ViterbiFilter()
 */
int
p7_ViterbiFilter_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
  register esl_neon_128i_t mpv, dpv, ipv;  /* previous row values                                       */
  register esl_neon_128i_t sv;		   /* temp storage of 1 curr row value in progress              */
  register esl_neon_128i_t dcv;		   /* delayed storage of D(i,q+1)                               */
  register esl_neon_128i_t xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register esl_neon_128i_t xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register esl_neon_128i_t Dmaxv;          /* keeps track of maximum D cell on row                      */
  int16_t  xE, xB, xC, xJ, xN;	           /* special states' scores                                    */
  int16_t  Dmax;		           /* maximum D cell score on row                               */
  int i;			           /* counter over sequence positions 1..L                      */
  int q;			           /* counter over vectors 0..nq-1                              */
  int Q  = P7_Q(om->M, p7_VWIDTH_NEON / sizeof(int16_t));  /* segment length: # of vectors              */
  esl_neon_128i_t *dp;
  esl_neon_128i_t *rsc;			   /* will point at om->ru[x] for residue x[i]                  */
  esl_neon_128i_t *tsc;			   /* will point into (and step thru) om->tu                    */
  esl_neon_128i_t  negInfv;
  int              status;

  /* Contract checks */
  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  /* note however that ViterbiFilter numerics are only guaranteed for local alignment, not glocal        */


  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_Reinit(ox, om->M)) != eslOK) ESL_EXCEPTION(status, "Reallocation of Vit filter matrix failed");
  dp = (esl_neon_128i_t *) ox->dp;   /* enables MMXf(), IMXf(), DMXf() access macros. Must be set AFTER the Reinit(), because ox->dp may get moved */

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->Vw   = p7_VWIDTH_NEON / sizeof(int16_t);
  ox->type = p7F_VITFILTER;
  ESL_DASSERT1(( ox->Vw = om->V / sizeof(int16_t)));

  /* -infinity is -32768 */
  negInfv.s16x8 = vdupq_n_s16(0);
  negInfv.s16x8 = vsetq_lane_s16(-32768, negInfv.s16x8, 7);  /* negInfv = 16-byte vector, 14 0 bytes + 2-byte value=-32768, for an OR operation. */

  /* Initialization. In unsigned arithmetic, -infinity is -32768
   */
  for (q = 0; q < Q; q++)
    MMXf(q).s16x8 = IMXf(q).s16x8 = DMXf(q).s16x8 =vdupq_n_s16(-32768);
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
      rsc   = (esl_neon_128i_t *) om->rwv[dsq[i]];
      tsc   = (esl_neon_128i_t *) om->twv;
      dcv.s16x8   = vdupq_n_s16(-32768);      /* "-infinity" */
      xEv.s16x8   = vdupq_n_s16(-32768);     
      Dmaxv.s16x8 = vdupq_n_s16(-32768);     
      xBv.s16x8   = vdupq_n_s16(xB);

      /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically; replace it with -32768.
       */
      mpv = MMXf(Q-1);  mpv.s16x8 = vextq_s16(negInfv.s16x8, mpv.s16x8, 7);
      dpv = DMXf(Q-1);  dpv.s16x8 = vextq_s16(negInfv.s16x8, dpv.s16x8, 7);
      ipv = IMXf(Q-1);  ipv.s16x8 = vextq_s16(negInfv.s16x8, ipv.s16x8, 7);

      for (q = 0; q < Q; q++)
      {
        /* Calculate new MMXf(i,q); don't store it yet, hold it in sv. */
        sv.s16x8   =                    vqaddq_s16(xBv.s16x8, (*tsc).s16x8);  tsc++;
        sv.s16x8   = vmaxq_s16(sv.s16x8, vqaddq_s16(mpv.s16x8, (*tsc).s16x8)); tsc++;
        sv.s16x8   = vmaxq_s16(sv.s16x8, vqaddq_s16(ipv.s16x8, (*tsc).s16x8)); tsc++;
        sv.s16x8   = vmaxq_s16(sv.s16x8, vqaddq_s16(dpv.s16x8, (*tsc).s16x8)); tsc++;
        sv.s16x8   = vqaddq_s16(sv.s16x8, (*rsc).s16x8);                      rsc++;
        xEv.s16x8  = vmaxq_s16(xEv.s16x8, sv.s16x8);

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
   
     	dcv.s16x8   = vqaddq_s16(sv.s16x8, (*tsc).s16x8);  tsc++;
        Dmaxv.s16x8 = vmaxq_s16(dcv.s16x8, Dmaxv.s16x8);

        /* Calculate and store I(i,q) */
        sv.s16x8      =                      vqaddq_s16(mpv.s16x8, (*tsc).s16x8);  tsc++;
        IMXf(q).s16x8 = vmaxq_s16 (sv.s16x8, vqaddq_s16(ipv.s16x8, (*tsc).s16x8)); tsc++;
      }

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      xE = esl_neon_hmax_s16(xEv);
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
      Dmax = esl_neon_hmax_s16(Dmaxv);
      if (Dmax + om->ddbound_w > xB) 
	{
	  /* Now we're obligated to do at least one complete DD path to be sure. */
	  /* dcv has carried through from end of q loop above */
	  dcv.s16x8 = vextq_s16(negInfv.s16x8, dcv.s16x8, 7); 
	  tsc = (esl_neon_128i_t *) om->twv + 7*Q;	/* set tsc to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    {
	      DMXf(q).s16x8 = vmaxq_s16(dcv.s16x8, DMXf(q).s16x8);	
	      dcv.s16x8     = vqaddq_s16(DMXf(q).s16x8, (*tsc).s16x8); tsc++;
	    }

	  /* We may have to do up to three more passes; the check
	   * is for whether crossing a segment boundary can improve
	   * our score. 
	   */
	  do {
	    dcv.s16x8 = vextq_s16(negInfv.s16x8, dcv.s16x8, 7); 
	    tsc = (esl_neon_128i_t *) om->twv + 7*Q;	/* set tsc to start of the DD's */
	    for (q = 0; q < Q; q++) 
	      {
		if (! esl_neon_any_gt_s16(dcv, DMXf(q))) break;
		DMXf(q).s16x8 = vmaxq_s16(dcv.s16x8, DMXf(q).s16x8);	
		dcv.s16x8     = vqaddq_s16(DMXf(q).s16x8, (*tsc).s16x8);   tsc++;
	      }	    
	  } while (q == Q);
	}
      else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
	{
	  
	  DMXf(0).s16x8 = vextq_s16(negInfv.s16x8, dcv.s16x8, 7);
	}
	  
#if eslDEBUGLEVEL > 0
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
  return eslOK;
}


#else // ! eslENABLE_NEON

/* Standard compiler-pleasing mantra for an #ifdef'd-out, empty code file. */
void p7_vitfilter_neon_silence_hack(void) { return; }
#if defined p7VITFILTER_NEON_TESTDRIVE || p7VITFILTER_NEON_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_NEON or not
