#include "hmmer.h"
#include "impl_avx.h"

//descriptions of functions can be found in p7_omx.c
P7_OMX *
p7_omx_Create_avx(int allocM, int allocL, int allocXL)
{
  P7_OMX  *ox     = NULL;
  int      i;
  int      status;

  ESL_ALLOC(ox, sizeof(P7_OMX));
  ox->dp_mem = NULL;
  ox->dpb    = NULL;
  ox->dpw    = NULL;
  ox->dpf    = NULL;
  ox->xmx    = NULL;
  ox->x_mem  = NULL;
  ox->last_written_by = none;
  /* DP matrix will be allocated for allocL+1 rows 0,1..L; allocQ4*p7X_NSCELLS columns */
  ox->allocR   = allocL+1;
  ox->validR   = ox->allocR;
  ox->allocQ4_avx  = p7O_NQF_AVX(allocM);
  ox->allocQ8_avx  = p7O_NQW_AVX(allocM);
  ox->allocQ16_avx = p7O_NQB_AVX(allocM);
  ox->ncells   = (int64_t) ox->allocR * (int64_t) ox->allocQ4_avx * 8;      /* # of DP cells allocated, where 1 cell contains MDI */

  ESL_ALLOC(ox->dp_mem_avx, sizeof(__m256) * (int64_t) ox->allocR * (int64_t) ox->allocQ4_avx * p7X_NSCELLS + 32); 
 
  ESL_ALLOC(ox->dpb_avx,    sizeof(__m256i *) * ox->allocR);
  ESL_ALLOC(ox->dpw_avx,    sizeof(__m256i *) * ox->allocR);
  ESL_ALLOC(ox->dpf_avx,    sizeof(__m256  *) * ox->allocR);

  ox->dpb_avx[0] = (__m256i *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));
  ox->dpw_avx[0] = (__m256i *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));
  ox->dpf_avx[0] = (__m256  *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));

  for (i = 1; i <= allocL; i++) {
    ox->dpf_avx[i] = ox->dpf_avx[0] + (int64_t) i * (int64_t) ox->allocQ4_avx  * p7X_NSCELLS;
    ox->dpw_avx[i] = ox->dpw_avx[0] + (int64_t) i * (int64_t) ox->allocQ8_avx  * p7X_NSCELLS;
    ox->dpb_avx[i] = ox->dpb_avx[0] + (int64_t) i * (int64_t) ox->allocQ16_avx;
  }

  ox->allocXR = allocXL+1;
  ESL_ALLOC(ox->x_mem,  sizeof(float) * ox->allocXR * p7X_NXCELLS + 15);  //pad these out for 512-bit vectors
  ox->xmx = (float *) ( ( (unsigned long int) ((char *) ox->x_mem  + 15) & (~0xf)));

  ox->M              = 0;
  ox->L              = 0;
  ox->totscale       = 0.0;
  ox->has_own_scales = TRUE;	/* most matrices are Forward, control their own scale factors */
#if eslDEBUGLEVEL > 0
  ox->debugging = FALSE;
  ox->dfp       = NULL;
#endif
  return ox;

 ERROR:
  p7_omx_Destroy(ox);
  return NULL;
}

P7_OMX *
p7_omx_Create_test_sse_avx(int allocM, int allocL, int allocXL)
{
  P7_OMX  *ox     = NULL;
  int      i;
  int      status;

  ESL_ALLOC(ox, sizeof(P7_OMX));
  ox->dp_mem = NULL;
  ox->dpb    = NULL;
  ox->dpw    = NULL;
  ox->dpf    = NULL;
  ox->xmx    = NULL;
  ox->x_mem  = NULL;
  ox->last_written_by = none;
  /* DP matrix will be allocated for allocL+1 rows 0,1..L; allocQ4*p7X_NSCELLS columns */
  ox->allocR   = allocL+1;
  ox->validR   = ox->allocR;
  ox->allocQ4  = p7O_NQF(allocM);
  ox->allocQ8  = p7O_NQW(allocM);
  ox->allocQ16 = p7O_NQB(allocM);
  ox->allocQ4_avx  = p7O_NQF_AVX(allocM);
  ox->allocQ8_avx  = p7O_NQW_AVX(allocM);
  ox->allocQ16_avx = p7O_NQB_AVX(allocM);
  ox->ncells   = (int64_t) ox->allocR * (int64_t) ox->allocQ4 * 4;      /* # of DP cells allocated, where 1 cell contains MDI */

  ESL_ALLOC(ox->dp_mem, sizeof(__m128) * (int64_t) ox->allocR * (int64_t) ox->allocQ4 * p7X_NSCELLS + 15);  /* floats always dominate; +15 for alignment */
  ESL_ALLOC(ox->dp_mem_avx, sizeof(__m256) * (int64_t) ox->allocR * (int64_t) ox->allocQ4_avx * p7X_NSCELLS + 32); 
  ESL_ALLOC(ox->dpb,    sizeof(__m128i *) * ox->allocR);
  ESL_ALLOC(ox->dpw,    sizeof(__m128i *) * ox->allocR);
  ESL_ALLOC(ox->dpf,    sizeof(__m128  *) * ox->allocR);
  ESL_ALLOC(ox->dpb_avx,    sizeof(__m256i *) * ox->allocR);
  ESL_ALLOC(ox->dpw_avx,    sizeof(__m256i *) * ox->allocR);
  ESL_ALLOC(ox->dpf_avx,    sizeof(__m256  *) * ox->allocR);

  ox->dpb[0] = (__m128i *) ( ( (unsigned long int) ((char *) ox->dp_mem + 15) & (~0xf)));
  ox->dpw[0] = (__m128i *) ( ( (unsigned long int) ((char *) ox->dp_mem + 15) & (~0xf)));
  ox->dpf[0] = (__m128  *) ( ( (unsigned long int) ((char *) ox->dp_mem + 15) & (~0xf)));
  ox->dpb_avx[0] = (__m256i *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));
  ox->dpw_avx[0] = (__m256i *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));
  ox->dpf_avx[0] = (__m256  *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));

  for (i = 1; i <= allocL; i++) {
    ox->dpf[i] = ox->dpf[0] + (int64_t) i * (int64_t) ox->allocQ4  * p7X_NSCELLS;
    ox->dpw[i] = ox->dpw[0] + (int64_t) i * (int64_t) ox->allocQ8  * p7X_NSCELLS;
    ox->dpb[i] = ox->dpb[0] + (int64_t) i * (int64_t) ox->allocQ16;
    ox->dpf_avx[i] = ox->dpf_avx[0] + (int64_t) i * (int64_t) ox->allocQ4_avx  * p7X_NSCELLS;
    ox->dpw_avx[i] = ox->dpw_avx[0] + (int64_t) i * (int64_t) ox->allocQ8_avx  * p7X_NSCELLS;
    ox->dpb_avx[i] = ox->dpb_avx[0] + (int64_t) i * (int64_t) ox->allocQ16_avx;
  }

  ox->allocXR = allocXL+1;
  ESL_ALLOC(ox->x_mem,  sizeof(float) * ox->allocXR * p7X_NXCELLS + 15);  //pad these out for 512-bit vectors
  ox->xmx = (float *) ( ( (unsigned long int) ((char *) ox->x_mem  + 15) & (~0xf)));

  ox->M              = 0;
  ox->L              = 0;
  ox->totscale       = 0.0;
  ox->has_own_scales = TRUE;	/* most matrices are Forward, control their own scale factors */
#if eslDEBUGLEVEL > 0
  ox->debugging = FALSE;
  ox->dfp       = NULL;
#endif
  return ox;

 ERROR:
  p7_omx_Destroy(ox);
  return NULL;
}

int
p7_omx_GrowTo_avx(P7_OMX *ox, int allocM, int allocL, int allocXL)
{
  void   *p;

  int     nqf_avx    = p7O_NQF_AVX(allocM);	   /* segment length; total # of striped vectors for uchar */
  int     nqw_avx    = p7O_NQW_AVX(allocM);	   /* segment length; total # of striped vectors for float */
  int     nqb_avx    = p7O_NQB_AVX(allocM);	   /* segment length; total # of striped vectors for float */

  int64_t ncells = (int64_t) (allocL+1) * (int64_t) nqf_avx * 8;
  int     reset_row_pointers = FALSE;
  int     i;
  int     status;
 
  /* If all possible dimensions are already satisfied, the matrix is fine */
  if (ox->allocQ4_avx*8 >= allocM && ox->validR > allocL && ox->allocXR >= allocXL+1) return eslOK;

  /* If the main matrix is too small in cells, reallocate it; 
   * and we'll need to realign/reset the row pointers later.
   */
  if (ncells > ox->ncells)
    {
      ESL_RALLOC(ox->dp_mem_avx, p, sizeof(__m256) * (int64_t) (allocL+1) * (int64_t) nqf_avx * p7X_NSCELLS + 31);
      ox->ncells = ncells;
      reset_row_pointers = TRUE;
    }

  /* If the X beams are too small, reallocate them. */
  if (allocXL+1 >= ox->allocXR)
    {
      ESL_RALLOC(ox->x_mem, p,  sizeof(float) * (allocXL+1) * p7X_NXCELLS + 15); 
      ox->allocXR = allocXL+1;
      ox->xmx     = (float *) ( ( (unsigned long int) ((char *) ox->x_mem  + 15) & (~0xf)));
    }

  /* If there aren't enough rows, reallocate the row pointers; we'll
   * realign and reset them later.
   */
  if (allocL >= ox->allocR)
    {
      ESL_RALLOC(ox->dpb_avx, p, sizeof(__m256i *) * (allocL+1));
      ESL_RALLOC(ox->dpw_avx, p, sizeof(__m256i *) * (allocL+1));
      ESL_RALLOC(ox->dpf_avx, p, sizeof(__m256  *) * (allocL+1));
      ox->allocR         = allocL+1;
      reset_row_pointers = TRUE;
    }

  /* must we widen the rows? */
  if (allocM > ox->allocQ4*4)
    reset_row_pointers = TRUE;

  /* must we set some more valid row pointers? */
  if (allocL >= ox->validR)
    reset_row_pointers = TRUE;

  /* now reset the row pointers, if needed */
  if (reset_row_pointers)
    {
      ox->dpb_avx[0] = (__m256i *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));
      ox->dpw_avx[0] = (__m256i *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));
      ox->dpf_avx[0] = (__m256  *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));
      ox->validR = ESL_MIN( ox->ncells / (nqf_avx * 8), ox->allocR);
      for (i = 1; i < ox->validR; i++)
	{
	  ox->dpb_avx[i] = ox->dpb_avx[0] + (int64_t) i * (int64_t) nqb_avx;
	  ox->dpw_avx[i] = ox->dpw_avx[0] + (int64_t) i * (int64_t) nqw_avx * p7X_NSCELLS;
	  ox->dpf_avx[i] = ox->dpf_avx[0] + (int64_t) i * (int64_t) nqf_avx * p7X_NSCELLS;
	}

      ox->allocQ4_avx  = nqf_avx;
      ox->allocQ8_avx  = nqw_avx;
      ox->allocQ16_avx = nqb_avx;

    }
  
  ox->M = 0;
  ox->L = 0;
  return eslOK;

 ERROR:
  return status;
}  

int
p7_omx_GrowTo_test_sse_avx(P7_OMX *ox, int allocM, int allocL, int allocXL)
{
  void   *p;
  int     nqf    = p7O_NQF(allocM);	   /* segment length; total # of striped vectors for uchar */
  int     nqw    = p7O_NQW(allocM);	   /* segment length; total # of striped vectors for float */
  int     nqb    = p7O_NQB(allocM);	   /* segment length; total # of striped vectors for float */
  int     nqf_avx    = p7O_NQF_AVX(allocM);	   /* segment length; total # of striped vectors for uchar */
  int     nqw_avx    = p7O_NQW_AVX(allocM);	   /* segment length; total # of striped vectors for float */
  int     nqb_avx    = p7O_NQB_AVX(allocM);	   /* segment length; total # of striped vectors for float */

  int64_t ncells = (int64_t) (allocL+1) * (int64_t) nqf_avx * 8;
  int     reset_row_pointers = FALSE;
  int     i;
  int     status;
 
  /* If all possible dimensions are already satisfied, the matrix is fine */
  if (ox->allocQ4_avx*8 >= allocM && ox->validR > allocL && ox->allocXR >= allocXL+1) return eslOK;

  /* If the main matrix is too small in cells, reallocate it; 
   * and we'll need to realign/reset the row pointers later.
   */
  if (ncells > ox->ncells)
    {
      ESL_RALLOC(ox->dp_mem, p, sizeof(__m128) * (int64_t) (allocL+1) * (int64_t) nqf * p7X_NSCELLS + 15);
      ESL_RALLOC(ox->dp_mem_avx, p, sizeof(__m256) * (int64_t) (allocL+1) * (int64_t) nqf_avx * p7X_NSCELLS + 31);
      ox->ncells = ncells;
      reset_row_pointers = TRUE;
    }

  /* If the X beams are too small, reallocate them. */
  if (allocXL+1 >= ox->allocXR)
    {
      ESL_RALLOC(ox->x_mem, p,  sizeof(float) * (allocXL+1) * p7X_NXCELLS + 15); 
      ox->allocXR = allocXL+1;
      ox->xmx     = (float *) ( ( (unsigned long int) ((char *) ox->x_mem  + 15) & (~0xf)));
    }

  /* If there aren't enough rows, reallocate the row pointers; we'll
   * realign and reset them later.
   */
  if (allocL >= ox->allocR)
    {
      ESL_RALLOC(ox->dpb, p, sizeof(__m128i *) * (allocL+1));
      ESL_RALLOC(ox->dpw, p, sizeof(__m128i *) * (allocL+1));
      ESL_RALLOC(ox->dpf, p, sizeof(__m128  *) * (allocL+1));
      ESL_RALLOC(ox->dpb_avx, p, sizeof(__m256i *) * (allocL+1));
      ESL_RALLOC(ox->dpw_avx, p, sizeof(__m256i *) * (allocL+1));
      ESL_RALLOC(ox->dpf_avx, p, sizeof(__m256  *) * (allocL+1));
      ox->allocR         = allocL+1;
      reset_row_pointers = TRUE;
    }

  /* must we widen the rows? */
  if (allocM > ox->allocQ4*4)
    reset_row_pointers = TRUE;

  /* must we set some more valid row pointers? */
  if (allocL >= ox->validR)
    reset_row_pointers = TRUE;

  /* now reset the row pointers, if needed */
  if (reset_row_pointers)
    {
      ox->dpb[0] = (__m128i *) ( ( (unsigned long int) ((char *) ox->dp_mem + 15) & (~0xf)));
      ox->dpw[0] = (__m128i *) ( ( (unsigned long int) ((char *) ox->dp_mem + 15) & (~0xf)));
      ox->dpf[0] = (__m128  *) ( ( (unsigned long int) ((char *) ox->dp_mem + 15) & (~0xf)));
      ox->dpb_avx[0] = (__m256i *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));
      ox->dpw_avx[0] = (__m256i *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));
      ox->dpf_avx[0] = (__m256  *) ( ( (unsigned long int) ((char *) ox->dp_mem_avx + 31) & (~0x1f)));
      ox->validR = ESL_MIN( ox->ncells / (nqf * 4), ox->allocR);
      for (i = 1; i < ox->validR; i++)
	{
	  ox->dpb[i] = ox->dpb[0] + (int64_t) i * (int64_t) nqb;
	  ox->dpw[i] = ox->dpw[0] + (int64_t) i * (int64_t) nqw * p7X_NSCELLS;
	  ox->dpf[i] = ox->dpf[0] + (int64_t) i * (int64_t) nqf * p7X_NSCELLS;
	  ox->dpb_avx[i] = ox->dpb_avx[0] + (int64_t) i * (int64_t) nqb_avx;
	  ox->dpw_avx[i] = ox->dpw_avx[0] + (int64_t) i * (int64_t) nqw_avx * p7X_NSCELLS;
	  ox->dpf_avx[i] = ox->dpf_avx[0] + (int64_t) i * (int64_t) nqf_avx * p7X_NSCELLS;
	}

      ox->allocQ4  = nqf;
      ox->allocQ8  = nqw;
      ox->allocQ16 = nqb;
      ox->allocQ4_avx  = nqf_avx;
      ox->allocQ8_avx  = nqw_avx;
      ox->allocQ16_avx = nqb_avx;
    }
  
  ox->M = 0;
  ox->L = 0;
  return eslOK;

 ERROR:
  return status;
}  


void
p7_omx_Destroy_avx(P7_OMX *ox)
{
  if (ox == NULL) return;
  if (ox->x_mem   != NULL) free(ox->x_mem);
  if (ox->dp_mem_avx  != NULL) free(ox->dp_mem_avx);
  if (ox->dpf_avx     != NULL) free(ox->dpf_avx);
  if (ox->dpw_avx     != NULL) free(ox->dpw_avx);
  if (ox->dpb_avx     != NULL) free(ox->dpb_avx);
  free(ox);
  return;
}

int
p7_omx_FDeconvert_avx(P7_OMX *ox, P7_GMX *gx)
{
  int Q = p7O_NQF_AVX(ox->M);
  int i, q, r, k;
  union { __m256 v; float p[8]; } u;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx; 			    

  for (i = 0; i <= ox->L; i++)
    {
      MMX(i,0) = DMX(i,0) = IMX(i,0) = -eslINFINITY;
      for (q = 0; q < Q; q++)
	{
	  u.v = MMO(ox->dpf_avx[i],q);  for (r = 0; r < 8; r++) { k = (Q*r)+q+1; if (k <= ox->M) MMX(i, (Q*r)+q+1) = u.p[r]; }
	  u.v = DMO(ox->dpf_avx[i],q);  for (r = 0; r < 8; r++) { k = (Q*r)+q+1; if (k <= ox->M) DMX(i, (Q*r)+q+1) = u.p[r]; }
	  u.v = IMO(ox->dpf_avx[i],q);  for (r = 0; r < 8; r++) { k = (Q*r)+q+1; if (k <= ox->M) IMX(i, (Q*r)+q+1) = u.p[r]; }
	}
      XMX(i,p7G_E) = ox->xmx[i*p7X_NXCELLS+p7X_E];
      XMX(i,p7G_N) = ox->xmx[i*p7X_NXCELLS+p7X_N];
      XMX(i,p7G_J) = ox->xmx[i*p7X_NXCELLS+p7X_J];
      XMX(i,p7G_B) = ox->xmx[i*p7X_NXCELLS+p7X_B];
      XMX(i,p7G_C) = ox->xmx[i*p7X_NXCELLS+p7X_C];
    }
  gx->L = ox->L;
  gx->M = ox->M;
  return eslOK;
}


int
p7_omx_DumpMFRow_avx(P7_OMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC)
{
  __m256i *dp = ox->dpb_avx[0];	
  int      M  = ox->M;
  int      Q  = p7O_NQB_AVX(M);
  uint8_t *v  = NULL;		/* array of unstriped scores  */
  int      q,z,k;
  union { __m256i v; uint8_t i[32]; } tmp;
  int      status;

  ESL_ALLOC(v, sizeof(unsigned char) * ((Q*32)+1));
  v[0] = 0;

  /* Header (if we're on the 0th row)  */
  if (rowi == 0)
    {
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M;  k++) fprintf(ox->dfp, "%3d ", k);
      fprintf(ox->dfp, "%3s %3s %3s %3s %3s\n", "E", "N", "J", "B", "C");
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M+5;  k++) fprintf(ox->dfp, "%3s ", "---");
      fprintf(ox->dfp, "\n");
    }

  /* Unpack and unstripe, then print M's. */
  for (q = 0; q < Q; q++) {
    tmp.v = dp[q];
    for (z = 0; z < 32; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d M ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", v[k]);

  /* The specials */
  fprintf(ox->dfp, "%3d %3d %3d %3d %3d\n", xE, xN, xJ, xB, xC);

  /* I's are all 0's; print just to facilitate comparison. */
  fprintf(ox->dfp, "%4d I ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", 0);
  fprintf(ox->dfp, "\n");

  /* D's are all 0's too */
  fprintf(ox->dfp, "%4d D ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", 0);
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

ERROR:
  free(v);
  return status;
}

int
p7_omx_DumpVFRow_avx(P7_OMX *ox, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC)
{
  __m256i *dp = ox->dpw_avx[0];	/* must set <dp> before using {MDI}MX macros */
  int      M  = ox->M;
  int      Q  = p7O_NQW_AVX(M);
  int16_t *v  = NULL;		/* array of unstriped, uninterleaved scores  */
  int      q,z,k;
  union { __m256i v; int16_t i[16]; } tmp;
  int      status;

  ESL_ALLOC(v, sizeof(int16_t) * ((Q*16)+1));
  v[0] = 0;

  /* Header (if we're on the 0th row)
   */
  if (rowi == 0)
    {
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M;  k++) fprintf(ox->dfp, "%6d ", k);
      fprintf(ox->dfp, "%6s %6s %6s %6s %6s\n", "E", "N", "J", "B", "C");
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M+5;  k++) fprintf(ox->dfp, "%6s ", "------");
      fprintf(ox->dfp, "\n");
    }

  /* Unpack and unstripe, then print M's. */
  for (q = 0; q < Q; q++) {
    tmp.v = MMXo(q);
    for (z = 0; z < 16; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d M ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%6d ", v[k]);

  /* The specials */
  fprintf(ox->dfp, "%6d %6d %6d %6d %6d\n", xE, xN, xJ, xB, xC);

  /* Unpack and unstripe, then print I's. */
  for (q = 0; q < Q; q++) {
    tmp.v = IMXo(q);
    for (z = 0; z < 16; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d I ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%6d ", v[k]);
  fprintf(ox->dfp, "\n");

  /* Unpack, unstripe, then print D's. */
  for (q = 0; q < Q; q++) {
    tmp.v = DMXo(q);
    for (z = 0; z < 16; z++) v[q+Q*z+1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d D ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%6d ", v[k]);
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

ERROR:
  free(v);
  return status;

}

int
p7_omx_DumpFBRow_avx(P7_OMX *ox, int logify, int rowi, int width, int precision, float xE, float xN, float xJ, float xB, float xC)
{
  __m256 *dp;
  int      M  = ox->M;
  int      Q  = p7O_NQF_AVX(M);
  float   *v  = NULL;		/* array of uninterleaved, unstriped scores  */
  int      q,z,k;
  union { __m256 v; float x[8]; } tmp;
  int      status;

  dp = (ox->allocR == 1) ? ox->dpf_avx[0] : ox->dpf_avx[rowi];	  /* must set <dp> before using {MDI}MX macros */

  ESL_ALLOC(v, sizeof(float) * ((Q*8)+1));
  v[0] = 0.;

  if (rowi == 0)
    {
      fprintf(ox->dfp, "      ");
      for (k = 0; k <= M;  k++) fprintf(ox->dfp, "%*d ", width, k);
      fprintf(ox->dfp, "%*s %*s %*s %*s %*s\n", width, "E", width, "N", width, "J", width, "B", width, "C");
      fprintf(ox->dfp, "      ");
      for (k = 0; k <= M+5;  k++) fprintf(ox->dfp, "%*s ", width, "--------");
      fprintf(ox->dfp, "\n");
    }

  /* Unpack, unstripe, then print M's. */
  for (q = 0; q < Q; q++) {
    tmp.v = MMXo(q);
    for (z = 0; z < 8; z++) v[q+Q*z+1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d M ", rowi);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0. ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);

 /* The specials */
  if (logify) fprintf(ox->dfp, "%*.*f %*.*f %*.*f %*.*f %*.*f\n",
		      width, precision, xE == 0. ? -eslINFINITY : log(xE),
		      width, precision, xN == 0. ? -eslINFINITY : log(xN),
		      width, precision, xJ == 0. ? -eslINFINITY : log(xJ),
		      width, precision, xB == 0. ? -eslINFINITY : log(xB), 
		      width, precision, xC == 0. ? -eslINFINITY : log(xC));
  else        fprintf(ox->dfp, "%*.*f %*.*f %*.*f %*.*f %*.*f\n",
		      width, precision, xE,   width, precision, xN, width, precision, xJ, 
		      width, precision, xB,   width, precision, xC);

  /* Unpack, unstripe, then print I's. */
  for (q = 0; q < Q; q++) {
    tmp.v = IMXo(q);
    for (z = 0; z < 8; z++) v[q+Q*z+1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d I ", rowi);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0. ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);
  fprintf(ox->dfp, "\n");

  /* Unpack, unstripe, then print D's. */
  for (q = 0; q < Q; q++) {
    tmp.v = DMXo(q);
    for (z = 0; z < 8; z++) v[q+Q*z+1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d D ", rowi);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0. ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

ERROR:
  free(v);
  return status;
}
