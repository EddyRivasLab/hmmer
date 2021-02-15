/* Routines for using two H4_REFMX's for ASC UP/DOWN matrices.
 *
 * There's no special data structure for reference ASC dynamic programming. Instead,
 * we use two H4_REFMX structures for ASC UP and DOWN matrices.
 *
 * Debugging and development routines specific for ASC DP matrices are gathered here.
 * Some unit tests in <reference_asc.c> create contrived comparisons where the
 * ensemble of valid paths in anchor set constrained DP is the same as in standard
 * unconstrained DP, so ASC decoding and standard decoding matrices should be equal,
 * cell by cell. The routines here are used to do those comparisons.
 *
 * Because these routines are used in unit tests, they call ESL_FAIL() with a failure
 * message when a comparison fails. This pattern allows you to set a debugging
 * breakpoint in esl_fail() to catch exactly where a failure happens.
 *
 * We need to marginalize over ASC UP and DOWN matrices before comparing to a
 * standard DP cell. This need dictates a different access pattern in the matrices.
 * See h4_anchorset.md for further explanation of the access pattern.
 *
 * Comparisons are done with an absolute floating point tolerance, not relative;
 * see h4_refmx_Compare() for additional explanation.
 */
#include "h4_config.h"

#include <math.h>

#include "h4_anchorset.h"
#include "h4_refmx.h"

#include "logsum.h"

/* Function:  h4_ascmx_fb_compare_std()
 * Synopsis:  Compare ASC Fwd or Bck matricee to standard, for equality
 * Incept:    SRE, Sat 13 Feb 2021
 *
 * Purpose:   Compare ASC Forward or Backward UP/DOWN matrices <ascu>/<ascd>
 *            (constrained by anchorset <anch>) for equality to a standard Forward or
 *            Backward matrix <std>. If the absolute difference between any
 *            corresponding values exceeds <abstol>, return <eslFAIL>.
 */
int
h4_ascmx_fb_compare_std(const H4_REFMX *std, const H4_REFMX *ascu, const H4_REFMX *ascd, const H4_ANCHORSET *anch, float abstol)
{
  char  msg[] = "comparison of ASC DP F|B matrix to standard DP F|B matrix failed";
  int   M     = std->M;
  int   L     = std->L;
  float rtol  = 0.0;
  float val;
  int   d,i,k,s;

  /* Contract checks */
  if (ascu->M != M || ascd->M != M) ESL_FAIL(eslFAIL, NULL, msg);
  if (ascu->L != L || ascd->L != L) ESL_FAIL(eslFAIL, NULL, msg);

  if      (std->type == h4R_FORWARD)  { if (ascu->type != h4R_ASC_FWD_UP || ascd->type != h4R_ASC_FWD_DOWN) ESL_FAIL(eslFAIL, NULL, msg); }
  else if (std->type == h4R_BACKWARD) { if (ascu->type != h4R_ASC_BCK_UP || ascd->type != h4R_ASC_BCK_DOWN) ESL_FAIL(eslFAIL, NULL, msg); }
  else    ESL_FAIL(eslFAIL, NULL, msg);

  for (d = 1, i = 0; i <= L; i++)
    {
      if (i == anch->a[d].i0) d++;  // d = idx of current UP matrix; d-1 is current DOWN
      for (k = 1; k <= M; k++)      //   ... alternatively, you can think d = idx of next anchor we'll reach (1..D; D+1 when we're past final anchor)
	for (s = 0; s < h4R_NSCELLS; s++)
	  {
            /* It's possible for an i,k to be valid in both DOWN and UP matrices, in which
             * case some paths are in each; we have to sum to get what's in the
             * unconstrained matrix.  Thus we traverse the afu/afd matrices in an unusual
             * pattern, sweeping out all i,k and testing which ones are valid in afu and
             * afd.
             */
	    val = -eslINFINITY;
	    if (k >= anch->a[d-1].k0) val = h4_logsum(val, H4R_MX(ascd,i,k,s)); // DOWN. sentinel k0(0)   = M+1, so no k gets evaluated for DOWN(d-1=0)
	    if (k <  anch->a[d].k0)   val = h4_logsum(val, H4R_MX(ascu,i,k,s)); // UP    sentinel k0(D+1) = 0,   so no k gets evaluated for UP(d=D+1)
	    
            /* Standard DP Backwards matrix can reach cells that ASC can't, through local exits.
             * So as a special case, ignore Backward cells where ASC val is -inf.
             */
            if (std->type != h4R_BACKWARD || val != -eslINFINITY) {
              if (esl_FCompareNew(H4R_MX(std,i,k,s), val, rtol, abstol) != eslOK) ESL_FAIL(eslFAIL, NULL, msg);
            }
	  }

      for (s = 0; s < h4R_NXCELLS; s++)
	{
	  /* special exception: ignore L state in this test, for Backwards. Reference
	   * backward can have valid prob for {MI}(k<k0) on the anchor row, and via M->I,
	   * reference can reach an Mk on any row i that ASC can't reach, in addition to
	   * Mk that ASC can reach; thus L can legitimately differ between ASC and
	   * reference Backward. [H5/26]
	   */
	  if (std->type == h4R_BACKWARD && s == h4R_L) continue;

	  if (esl_FCompareNew(H4R_XMX(std,i,s), H4R_XMX(ascd,i,s), rtol, abstol) != eslOK) ESL_FAIL(eslFAIL, NULL, msg);
	}
    }
  return eslOK;
}




/* Function:  h4_ascmx_pp_compare_std()
 * Synopsis:  Compare ASC Decoding to standard decoding, for equality
 * Incept:    SRE, Fri 12 Feb 2021
 *
 * Purpose:   Like <h4_ascmx_fb_compare_std()> but for posterior decoding, not
 *            forward/backward matrices.
 */
int
h4_ascmx_pp_compare_std(const H4_REFMX *rxd, const H4_REFMX *apu, const H4_REFMX *apd, const H4_ANCHORSET *anch, float abstol)
{
  char  msg[] = "comparison of ASC decoding to standard decoding matrices failed";
  int   M     = rxd->M;
  int   L     = rxd->L;
  float rtol  = 0.0;      // all floating point comparisons are absolute
  int   d, i, k, s;
  float ascval;

  if (apu->M != M || apd->M != M) ESL_FAIL(eslFAIL, NULL, msg);
  if (apu->L != L || apd->L != L) ESL_FAIL(eslFAIL, NULL, msg);
  if (rxd->type != h4R_DECODING || apu->type != h4R_ASC_DECODE_UP || apd->type != h4R_ASC_DECODE_DOWN) ESL_FAIL(eslFAIL, NULL, msg);

  /* We use an access pattern that lets us check UP and DOWN sectors
   * for the same i,k,s,d cell at the same time, because of the
   * UP/DOWN marginalization.
   */
  for (d = 1, i = 0; i <= L; i++)
    {      
      if (i == anch->a[d].i0) d++;         // d = index of current UP matrix. d-1 = DOWN matrix;  d will be D+1 for final chunk
      for (k = 1; k <= M; k++)             //   ... other way to think about it: d = index of next anchor we will reach
	for (s = 0; s < h4R_NSCELLS; s++)
	  {
	    ascval = 0.0;
	    if (k >= anch->a[d-1].k0) ascval += H4R_MX(apd,i,k,s);   // sentinel k0(0)   = M+1, so no k gets evaluated for DOWN(d-1=0) 
	    if (k <  anch->a[d].k0)   ascval += H4R_MX(apu,i,k,s);   // sentinel k0(D+1) = 0,   so no k gets evaluated for UP(d=D+1)
	    
	    if (esl_FCompareNew(H4R_MX(rxd,i,k,s), ascval, rtol, abstol) != eslOK) ESL_FAIL(eslFAIL, NULL, msg);
	  }
      for (s = 0; s < h4R_NXCELLS; s++)
	if (esl_FCompareNew(H4R_XMX(rxd,i,s), H4R_XMX(apd,i,s), rtol, abstol) != eslOK)  ESL_FAIL(eslFAIL, NULL, msg);
    }
  return eslOK;
}

/* Function:  h4_ascmx_pp_compare_path()
 * Synopsis:  Compare a single-pathed test's path to a ASC decoding matrix
 * Incept:    SRE, Sat 13 Feb 2021
 *
 * Purpose:   For a single pathed test, the cells in a posterior decoding matrix
 *            have values 1.0 or 0.0; specifically, 1.0 wherever the single path
 *            visits a cell. Test this, given the single path <pi> and
 *            ASC posterior decoding matrix <apu/apd>. Return <eslFAIL> if it's
 *            not so; <eslOK> if it is.
 */
int
h4_ascmx_pp_compare_path(const H4_PATH *pi, const H4_REFMX *apu, const H4_REFMX *apd, const H4_ANCHORSET *anch, float abstol)
{
  H4_REFMX *rxd = h4_refmx_Create(apu->M, apu->L);
  int       status;

  rxd->M    = apu->M;
  rxd->L    = apu->L;
  rxd->type = h4R_DECODING;

  if ((status = h4_refmx_SetValues(rxd, 0.0)) != eslOK) goto ERROR;
  if ((status = h4_refmx_CountPath(pi, rxd))  != eslOK) goto ERROR;
  if ((status = h4_ascmx_pp_compare_std(rxd, apu, apd, anch, abstol)) != eslOK) goto ERROR;

  h4_refmx_Destroy(rxd);
  return eslOK;

 ERROR:
  h4_refmx_Destroy(rxd);
  return status;
}


/* Function:  h4_ascmx_compare_asc()
 * Synopsis:  Compare ASC matrices for equality
 * Incept:    SRE, Sat 13 Feb 2021
 *
 * Purpose:   Compare ASC UP/DOWN matrices <au1>/<ad1> to another pair <au2>/<ad2> for
 *            equality. If the absolute difference between any corresponding values
 *            exceeds <abstol>, return <eslFAIL>.
 */
int
h4_ascmx_compare_asc(const H4_REFMX *au1, const H4_REFMX *ad1, const H4_REFMX *au2, const H4_REFMX *ad2, const H4_ANCHORSET *anch, float abstol)
{
  char  msg[] = "comparison of two reference ASC DP matrices failed";
  int   M     = au1->M;
  int   L     = au1->L;
  float rtol  = 0.0;     // all floating point comparisons are absolute
  int   d, i, k, s;

  if (ad1->M != M || au2->M != M || ad2->M != M) ESL_FAIL(eslFAIL, NULL, msg);
  if (ad1->L != L || au2->L != L || ad2->L != L) ESL_FAIL(eslFAIL, NULL, msg);

  if      (au1->type == h4R_ASC_FWD_UP)    { if (ad1->type != h4R_ASC_FWD_DOWN    || au2->type != h4R_ASC_FWD_UP    || ad2->type != h4R_ASC_FWD_DOWN)     ESL_FAIL(eslFAIL, NULL, msg); }
  else if (au1->type == h4R_ASC_BCK_UP)    { if (ad1->type != h4R_ASC_BCK_DOWN    || au2->type != h4R_ASC_BCK_UP    || ad2->type != h4R_ASC_BCK_DOWN)     ESL_FAIL(eslFAIL, NULL, msg); }
  else if (au1->type == h4R_ASC_DECODE_UP) { if (ad1->type != h4R_ASC_DECODE_DOWN || au2->type != h4R_ASC_DECODE_UP || ad2->type != h4R_ASC_DECODE_DOWN)  ESL_FAIL(eslFAIL, NULL, msg); }
  else     ESL_FAIL(eslFAIL, NULL, msg);

  d = 1;
  for (i = 0; i <= L; i++)
    {
      if (i == anch->a[d].i0) d++;  

      /* DOWN row, if one exists for this i */
      for (k = anch->a[d-1].k0; k <= M; k++)   // sentinel k0(0) = M+1, so no k gets evaluated at d=1 for nonexistent DOWN(0)
	for (s = 0; s < h4R_NSCELLS; s++)   //   ... i.e., first leg has only an UP(1) matrix.
	  if (esl_FCompareNew(H4R_MX(ad1,i,k,s), H4R_MX(ad2,i,k,s), rtol, abstol) != eslOK) ESL_FAIL(eslFAIL, NULL, msg);

      /* specials exist for all rows, stored in DOWN matrix */
      for (s = 0; s < h4R_NXCELLS; s++)
	if (esl_FCompareNew( H4R_XMX(ad1,i,s), H4R_XMX(ad2,i,s), rtol, abstol) != eslOK) ESL_FAIL(eslFAIL, NULL, msg);
      
      /* UP row, if one exists for this i */
      for (k = 1; k < anch->a[d].k0; k++)   // sentinel k0(D+1) = 0, so no k gets evaluated at d=D+1 for nonexistent UP(D+1)
	for (s = 0; s < h4R_NSCELLS; s++)
	  if (esl_FCompareNew( H4R_MX(au1,i,k,s), H4R_MX(au2,i,k,s), rtol, abstol) != eslOK) ESL_FAIL(eslFAIL, NULL, msg);
    }
  return eslOK;
}



/* Function:  h4_ascmx_maxdiff_std()
 * Synopsis:  Find the max abs difference between cells of ASC vs. standard DP matrix.
 * Incept:    SRE, Sat 13 Feb 2021
 */
int
h4_ascmx_maxdiff_std(const H4_REFMX *std, const H4_REFMX *ascu, const H4_REFMX *ascd, const H4_ANCHORSET *anch, float *ret_maxdiff)
{
  int   d,i,k,s;
  float ascval;
  float diff;
  float maxdiff = 0.;

  d = 1;
  for (i = 0; i <= std->L; i++)
    {
      if (i == anch->a[d].i0) d++;
      for (k = 1; k <= std->M; k++)
	for (s = 0; s < h4R_NSCELLS; s++)
	  {
	    ascval = 0.0;	
	    if (k >= anch->a[d-1].k0) ascval += H4R_MX(ascd,i,k,s); 
	    if (k <  anch->a[d].k0)   ascval += H4R_MX(ascu,i,k,s); 

	    diff = ascval - H4R_MX(std,i,k,s);  
	    if (fabs(diff) > fabs(maxdiff)) maxdiff = diff;
	  }
      for (s = 0; s < h4R_NXCELLS; s++)
	{
	  diff = H4R_XMX(ascd,i,s) - H4R_XMX(std,i,s);
	  if (fabs(diff) > fabs(maxdiff)) maxdiff = diff;
	}
    }
  *ret_maxdiff = maxdiff;
  return eslOK;
}


/* Function:  h4_ascmx_pp_maxdiff_path()
 * Synopsis:  Find max absolute diff from expected 0.0/1.0 in ASC decoding of singlepath test
 * Incept:    SRE, Sat 13 Feb 2021
 */
int
h4_ascmx_pp_maxdiff_path(const H4_PATH *pi, const H4_REFMX *apu, const H4_REFMX *apd, const H4_ANCHORSET *anch, float *ret_maxdiff)
{
  H4_REFMX *rxd = h4_refmx_Create(apu->M, apu->L);
  int       status;

  if ((status = h4_refmx_SetValues(rxd, 0.0)) != eslOK) goto ERROR;
  if ((status = h4_refmx_CountPath(pi, rxd))  != eslOK) goto ERROR;
  if ((status = h4_ascmx_maxdiff_std(rxd, apu, apd, anch, ret_maxdiff)) != eslOK) goto ERROR;

  h4_refmx_Destroy(rxd);
  return eslOK;

 ERROR:
  h4_refmx_Destroy(rxd);
  return status;
}


/* Function:  h4_ascmx_pp_max_overage()
 * Synopsis:  Find maximum excess over 1.0 in ASC decoding matrix
 * Incept:    SRE, Sat 13 Feb 2021
 */
int
h4_ascmx_pp_max_overage(const H4_REFMX *apu, const H4_REFMX *apd, const H4_ANCHORSET *anch, float *ret_max_overage)
{
  int d, i, k, s;
  float max = 0.;
  
  d = 1;
  for (i = 0; i <= apu->L; i++)
    {
      if (i == anch->a[d].i0) d++;

      for (k = anch->a[d-1].k0; k <= apu->M; k++)
	for (s = 0; s < h4R_NSCELLS; s++)
	  if (H4R_MX(apd,i,k,s) > max) max = H4R_MX(apd,i,k,s);

      for (s = 0; s < h4R_NXCELLS; s++)
	if (H4R_XMX(apd,i,s) > max) max = H4R_XMX(apd,i,s);

      for (k = 1; k < anch->a[d].k0; k++)
	for (s = 0; s < h4R_NSCELLS; s++)
	  if (H4R_MX(apu,i,k,s) > max) max = H4R_MX(apu,i,k,s);
    }
  
  *ret_max_overage = (max > 1.0 ? max-1.0 : 0.0);
  return eslOK;
}


/* Function:  h4_ascmx_pp_Validate()
 * Synopsis:  Validate ASC decoding UP/DOWN matrices
 * Incept:    SRE, Sat 13 Feb 2021
 *
 * Purpose:   Testing here is less thorough than <h4_refmx_Validate()>. Here we're just
 *            looking for obviously bad values that are not probabilities: values <
 *            0.0 or > 1.0 + <abstol>.
 */
int
h4_ascmx_pp_Validate(const H4_REFMX *apu, const H4_REFMX *apd, const H4_ANCHORSET *anch, float abstol, char *errbuf)
{
  int L = apu->L;
  int M = apu->M;
  int d,i,k,s;
  float val;
  
  if (apd->L != L || apd->M != M)       ESL_FAIL(eslFAIL, errbuf, "dimensions of UP, DOWN matrices don't match");
  if (apu->type != h4R_ASC_DECODE_UP)   ESL_FAIL(eslFAIL, errbuf, "apu is not an ASC UP matrix");
  if (apd->type != h4R_ASC_DECODE_DOWN) ESL_FAIL(eslFAIL, errbuf, "apd is not an ASC DOWN matrix");

  d = 1;
  for (i = 0; i <= L; i++)
    {
      if (i == anch->a[d].i0) d++;

      /* DOWN sector for domain d, if one exists for this i */
      for (k = anch->a[d-1].k0; k <= M; k++)
	for (s = 0; s < h4R_NSCELLS; s++)
	  {
	    val = H4R_MX(apd,i,k,s);
	    if (isnan(val) || val < 0. || val > 1.0+abstol)
	      ESL_FAIL(eslFAIL, errbuf, "bad DOWN value %f at i=%d, k=%d, %s", val, i, k, h4_refmx_DecodeState(s));
	  }	

      /* specials exist on all rows, stored in <apd> */
      for (s = 0; s < h4R_NXCELLS; s++)
	{
	  val = H4R_XMX(apd,i,s);
	  if (isnan(val) || val < 0. || val > 1.0+abstol) 
	    ESL_FAIL(eslFAIL, errbuf, "bad special value at i=%d, %s", i, h4_refmx_DecodeSpecial(s));
	}	

      /* UP sector for domain d */
      for (k = 1; k < anch->a[d].k0; k++)
	for (s = 0; s < h4R_NSCELLS; s++)
	  {
	    val = H4R_MX(apu,i,k,s);
	    if (isnan(val) || val < 0. || val > 1.0+abstol)
	      ESL_FAIL(eslFAIL, errbuf, "bad UP value %f at i=%d, k=%d, %s", val, i, k, h4_refmx_DecodeState(s));
	    }	
    }
  return eslOK;
}


