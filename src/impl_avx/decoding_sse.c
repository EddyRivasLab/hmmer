/* Posterior decoding algorithms; SSE versions.
 * 
 * Contents:
 *   1. Posterior decoding algorithms.
 */

#include <p7_config.h>

#include <stdio.h>
#include <math.h>

#include <x86intrin.h>

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_avx.h"

int
p7_Decoding_sse(const P7_OPROFILE *om, const P7_OMX *oxf, P7_OMX *oxb, P7_OMX *pp)
{
  __m128 *ppv;
  __m128 *fv;
  __m128 *bv;
  __m128  totrv;
  int    L  = oxf->L;
  int    M  = om->M;
  int    Q  = p7O_NQF(M);	
  int    i,q;
  float  scaleproduct = 1.0 / oxb->xmx[p7X_N];

  pp->M = M;
  pp->L = L;
  pp->last_written_by = sse;
  ppv = pp->dpf[0];
  for (q = 0; q < Q; q++) {
    *ppv = _mm_setzero_ps(); ppv++;
    *ppv = _mm_setzero_ps(); ppv++;
    *ppv = _mm_setzero_ps(); ppv++;
  }
  pp->xmx[p7X_E] = 0.0;
  pp->xmx[p7X_N] = 0.0;
  pp->xmx[p7X_J] = 0.0;
  pp->xmx[p7X_C] = 0.0;
  pp->xmx[p7X_B] = 0.0;

  for (i = 1; i <= L; i++)
    {
      ppv   =  pp->dpf[i];
      fv    = oxf->dpf[i];
      bv    = oxb->dpf[i];
      totrv = _mm_set1_ps(scaleproduct * oxf->xmx[i*p7X_NXCELLS+p7X_SCALE]);

      for (q = 0; q < Q; q++)
	{
	  /* M */
	  *ppv = _mm_mul_ps(*fv,  *bv);
	  *ppv = _mm_mul_ps(*ppv,  totrv);
	  ppv++;  fv++;  bv++;

	  /* D */
	  *ppv = _mm_setzero_ps();
	  ppv++;  fv++;  bv++;

	  /* I */
	  *ppv = _mm_mul_ps(*fv,  *bv);
	  *ppv = _mm_mul_ps(*ppv,  totrv);
	  ppv++;  fv++;  bv++;
	}
      pp->xmx[i*p7X_NXCELLS+p7X_E] = 0.0;
      pp->xmx[i*p7X_NXCELLS+p7X_N] = oxf->xmx[(i-1)*p7X_NXCELLS+p7X_N] * oxb->xmx[i*p7X_NXCELLS+p7X_N] * om->xf[p7O_N][p7O_LOOP] * scaleproduct;
      pp->xmx[i*p7X_NXCELLS+p7X_J] = oxf->xmx[(i-1)*p7X_NXCELLS+p7X_J] * oxb->xmx[i*p7X_NXCELLS+p7X_J] * om->xf[p7O_J][p7O_LOOP] * scaleproduct;
      pp->xmx[i*p7X_NXCELLS+p7X_C] = oxf->xmx[(i-1)*p7X_NXCELLS+p7X_C] * oxb->xmx[i*p7X_NXCELLS+p7X_C] * om->xf[p7O_C][p7O_LOOP] * scaleproduct;
      pp->xmx[i*p7X_NXCELLS+p7X_B] = 0.0;

      if (oxb->has_own_scales) scaleproduct *= oxf->xmx[i*p7X_NXCELLS+p7X_SCALE] /  oxb->xmx[i*p7X_NXCELLS+p7X_SCALE];
    }

  if (isinf(scaleproduct)) return eslERANGE;
  else                     return eslOK;
}


/*------------------ end, posterior decoding --------------------*/

