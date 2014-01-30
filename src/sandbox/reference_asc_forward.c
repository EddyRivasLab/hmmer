
int
p7_ReferenceASCForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
		       P7_COORD2 *anch, int D,
		       P7_REFMX *mxu, P7_REFMX *mxd, float *opt_sc)
{

  int  M = gm->M;		/* for a bit more clarity */
  int d, i, k, s, iu;
  float dlv, dgv;

  /* Initialize specials on row 0 for the d=0 UP matrix sector */
  dpc = mxu->dp[0] + (M+1) * p7R_NSCELLS;
  dpc[p7R_E]      = -eslINFINITY;
  dpc[p7R_N]      = 0.0;
  dpc[p7R_J]      = -eslINFINITY;                                  
  dpc[p7R_B]      = gm->xsc[p7P_N][p7P_MOVE];                      
  dpc[p7R_L] = xL = gm->xsc[p7P_N][p7P_MOVE] + gm->xsc[p7P_B][0];  
  dpc[p7R_G] = xG = gm->xsc[p7P_N][p7P_MOVE] + gm->xsc[p7P_B][1];  
  dpc[p7R_C]      = -eslINFINITY;                                  
  dpc[p7R_JJ]     = -eslINFINITY;                             
  dpc[p7R_CC]     = -eslINFINITY;                             
  
  /* Initialize row 1 of the d=0 UP matrix sector */
  rsc = gm->rsc[dsq[i]] + p7P_NR; /* position at k=1 */
  tsc = gm->tsc;		  /* position at k=0 */
  dpc = mxu->dp[1];
  for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; /* dpc now on k=1 */
  dlv = dgv = -eslINFINITY;

  for (k = 1; k <= anch[0].n2; k++)
    {
      /* match states MLk, MGk are reachable only by B->{LG}->Mk entry*/
      mlv = *dpc++ = *rsc +  xL  + *(tsc + p7P_LM);
      mgv = *dpc++ = *rsc +  xG  + *(tsc + p7P_GM);

      tsc += p7P_NTRANS;    /* tsc advances to transitions in states k     */

      /* Insert states are not reachable on this initial row */
      *dpc++ = -eslINFINITY;
      *dpc++ = -eslINFINITY;
      rsc+=2;		/* rsc advances to next match state emission for next k */

      /* E state is not reachable from UP matrix sector */

      /* Delete state, deferred storage trick */
      *dpc++ = dlv;
      *dpc++ = dgv;
      dlv = p7_FLogsum( mlv + *(tsc + p7P_MD), dlv + *(tsc + p7P_DD));
      dgv = p7_FLogsum( mgv + *(tsc + p7P_MD), dgv + *(tsc + p7P_DD));
    }      
  dpc = mxu->dp[1] + (M+1) * p7R_NSCELLS;
  dpc[p7R_E]      = -eslINFINITY;
  dpc[p7R_N]      = gm->xsc[p7P_N][p7P_LOOP];
  dpc[p7R_J]      = -eslINFINITY;                                  
  dpc[p7R_B]      = dpc[p7R_N] + gm->xsc[p7P_N][p7P_MOVE]
  dpc[p7R_L] = xL = dpc[p7R_B] + gm->xsc[p7P_B][0];  
  dpc[p7R_G] = xG = dpc[p7R_B] + gm->xsc[p7P_B][1];  
  dpc[p7R_C]      = -eslINFINITY;                                  
  dpc[p7R_JJ]     = -eslINFINITY;                             
  dpc[p7R_CC]     = -eslINFINITY;                             
  


  /* More initialization; the rest of the d=0 UP matrix sector */
  for (i = 2; i <= anch[0].n1; i++)
    {
      rsc = gm->rsc[dsq[i]] + p7P_NR;	
      tsc = gm->tsc;			

      dpp = rmx->dp[i-1];               
      dpc = rmx->dp[i];                 
      for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;
      dlv = dgv = -eslINFINITY;

      for (k = 1; k <= anch[0].n2; k++)
	{
	  mlv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_ML) + *(tsc + p7P_MM),
						       *(dpp+p7R_IL) + *(tsc + p7P_IM)),
					    p7_FLogsum(*(dpp+p7R_DL) + *(tsc + p7P_DM),
						       xL            + *(tsc + p7P_LM)));

	  mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_MG) + *(tsc + p7P_MM),
						       *(dpp+p7R_IG) + *(tsc + p7P_IM)),
 					    p7_FLogsum(*(dpp+p7R_DG) + *(tsc + p7P_DM),
						       xG            + *(tsc + p7P_GM)));
	  rsc++;              
	  tsc += p7P_NTRANS;  
	  dpp += p7R_NSCELLS;

	  *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_ML) + *(tsc + p7P_MI), *(dpp + p7R_IL) + *(tsc + p7P_II));
	  *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_MG) + *(tsc + p7P_MI), *(dpp + p7R_IG) + *(tsc + p7P_II));
	  rsc++;		
	  /* E state is not reachable from UP sector */
	  *dpc++ = dlv;
	  *dpc++ = dgv;
	  dlv = p7_FLogsum( mlv + *(tsc + p7P_MD), dlv + *(tsc + p7P_DD));
	  dgv = p7_FLogsum( mgv + *(tsc + p7P_MD), dgv + *(tsc + p7P_DD));
	}

      

  

  for (d = 0; d < D-1; d++)
    {

      /* "Down" matrix: matrix section downstream of anchor d */
      
      /* We have already calculated the anchor(d) i,k cell  */


      /* Initialize the rest of the topmost row on the section, anch[d].i */
      i = anch[d].n1;

      dpc  = mxd->dp[i];
      dpc += p7R_NSCELLS * (anch[d].n2 + 1);

      dlv = mlv + (*tsc + p7P_MD); /* mlv = from the anchor. */
      dgv = mgv + (*tsc + p7P_MD);
      

      for (k = anch[d].n2 + 1; k < M; k++)
	{
	  *dpc++ = -eslINFINITY;
	  *dpc++ = -eslINFINITY;
	  *dpc++ = -eslINFINITY;
	  *dpc++ = -eslINFINITY;
	  *dpc++ = dlv;
	  *dpc++ = dgv;
	  tsc   += p7P_NTRANS;
	  
	  dlv    = dlv + (*tsc + p7P_DD);
	  dgv    = dgv + (*tsc + p7P_DD);
	}
      


    }
  



}
