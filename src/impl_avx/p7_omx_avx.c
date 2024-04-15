#include "hmmer.h"
#include "impl_avx.h"

int
p7_omx_FDeconvert_avx(P7_OMX *ox, P7_GMX *gx)
{
  int Q = p7O_NQF_AVX(ox->M);
  int i, q, r, k;
  union { __m256 v; float p[4]; } u;
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