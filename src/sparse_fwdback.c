
#define P7_SXV(val, k, s) (val + k*p7B_NSCELLS + s)

int
p7_SparseForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_SPARSEMX *sx, float *opt_sc)
{
  P7_SPARSEMASK *sm   = sx->sm;
  int i;
  float const *tsc    = gm->tsc;		 /* sets up TSC() macro, access to profile's transitions */
  float const *rsc;				 /* will be set up for MSC(), ISC() macros for residue scores */
  int   nskip;

#ifdef P7_DEBUGGING  
  if (L != sx->L) ESL_EXCEPTION("L, sx->L disagree: sparse matrix wasn't allocated or reinitialized for this sequence");
#endif

  for (i = 1; i <= L; i++)
    {
      if (! sm->nin[i]) { nskip++; continue; }   /* skip rows that have no included cells */




      

      rsc = gm->rsc[dsq[i]];	/* now MSC(k), ISC(k) residue score macros work */

      inp = sx->kin[i-1];
      inc = sx->kin[i];

      dpc = sx->val + ( i    * sx->W * p7R_NSCELLS); /* dpc (current row values) will step through all sparse cells/values, contiguously */
      dpp = sx->val + ((i-1) * sx->W * p7R_NSCELLS); /* dpp (prev row values) will need to hop around in sparse cells on prev row  */
      
      for (kp = 0, kc = 0, k=-1; inc[kc] != p7_SX_SENTINEL; kc++)
	{
	  if (inc[kc] != k+1) { dlc = dgc = -eslINFINITY; } /* inc[kc] >= 1 (including SENTINEL) so k=-1 init suffices to make this test fail, and init dlc/dgc to -inf*/
	  k = inc[kc];
	  
	  while (inp[kp] < k-1)  kp++;
	  cp = (inp[kp] == k-1) ? dpp[kp*p7R_NSCELLS] : sx->val; /* sx->val, the first cell, is all -inf */

	  *dpc++ = mlc = MSC(k) + p7_FLogsum( p7_FLogsum( cp[p7R_ML] + TSC(p7P_MM, k-1),
							  cp[p7R_IL] + TSC(p7P_IM, k-1)),
					      p7_FLogsum( cp[p7R_DL] + TSC(p7P_DM, k-1),
							         xL  + TSC(p7P_LM, k-1)));

	  *dpc++ = mgc = MSC(k) + p7_FLogsum( p7_FLogsum( cp[p7R_MG] + TSC(p7P_MM, k-1),
							  cp[p7R_IG] + TSC(p7P_IM, k-1)),
					      p7_FLogsum( cp[p7R_DG] + TSC(p7P_DM, k-1),
							         xG  + TSC(p7P_GM, k-1)));

	  while (inp[kp] < k) kp++;
	  cp = (inp[kp] == k) ? dpp[kp*p7R_NSCELLS] : sx->val; /* sx->val, the first cell, is all -inf */

	  *dpc++ = ISC(k) + p7_FLogsum( cp[p7R_ML] + TSC(p7P_MI,k),  cp[p7R_IL] + TSC(p7P_II, k));
	  *dpc++ = ISC(k) + p7_FLogsum( cp[p7R_MG] + TSC(p7P_MI,k),  cp[p7R_IG] + TSC(p7P_II, k));
	  
	  xE = p7_FLogsum(xE, p7_FLogsum(mlc, dlc));

	  *dpc++ = dlc;
	  *dpc++ = dgc;
	  dlc = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
	  dgc = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	}
      



    }



}
