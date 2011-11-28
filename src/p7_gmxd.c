

P7_GMXD *
p7_gmxd_Create(int M, int L)
{
  P7_GMXD *gxd = NULL;
  int      r;
  int      status;


  ESL_ALLOC(gxd, sizeof(P7_GMXD));
  gxd->dp_mem = NULL;
  gxd->dp     = NULL;

  gxd->allocR = L+1;
  gxd->allocW = (M+1) * p7GD_NSCELLS + p7GD_NXCELLS;
  gxd->allocN = (int64_t) gxd->allocR * (int64_t) gxd->allocW;

  ESL_ALLOC(gxd->dp_mem, sizeof(float  ) * gxd->allocN);
  ESL_ALLOC(gxd->dp,     sizeof(float *) * gxd->allocR);
  for (r = 0; r < gxd->allocR; r++)
    gxd->dp[r] = gxd->dp_mem + (r * gxd->allocW);

  gxd->validR = gxd->allocR;
  gxd->M      = 0;
  gxd->L      = 0;
  return gxd;

 ERROR:
  if (gxd) p7_gmxd_Destroy(gxd);
  return NULL;
}

int
p7_gmxd_GrowTo(P7_GMXD *gxd, int M, int L)
{
  int      W        = (M+1) * p7GD_NSCELLS + p7GD_NXCELLS;
  int      R        = L+1;
  uint64_t N        = (int64_t) R * (int64_t) W;
  int      do_reset = FALSE;
  int      r;
  int      status;

  /* are we already big enough? */
  if (W <= gxd->allocW && R <= gxd->validR) return eslOK;

  /* must we reallocate the matrix cells? */
  if (N > gxd->allocN)
    {
      ESL_REALLOC(gxd->dp_mem, sizeof(float) * N);
      gxd->allocN = N;
      do_reset    = TRUE;
    }
  
  /* must we reallocate the row pointers? */
  if (R > gxd->allocR)
    {
      ESL_REALLOC(gxd->dp, sizeof(float *) * R);
      gxd->allocR = R;
      do_reset    = TRUE;
    }

  /* must we widen the rows? */
  if (W > gxd->allocW) do_reset = TRUE;

  /* must we set some more valid row pointers? */
  if (R > gxd->validR) do_reset = TRUE;

  /* resize rows, reset valid row pointers */
  if (do_reset)
    {
      gxd->allocW = W;
      gxd->validR = ESL_MIN(gxd->allocR, (int) ( gxd->allocN / (uint64_t) gxd->allocW));
      for (r = 0; r < gxd->validR; r++)
	gxd->dp[r] = gxd->dp_mem + (r * gxd->allocW);
    }

  gxd->M = 0;
  gxd->L = 0;
  return eslOK;

 ERROR:
  return status;
}
	

void
p7_gmxd_Destroy(P7_GMXD *gxd)
{
  if (!gxd) return;
  if (gxd->dp_mem) free(gxd->dp_mem);
  if (gxd->dp)     free(gxd->dp);
  free(gxd);
}
