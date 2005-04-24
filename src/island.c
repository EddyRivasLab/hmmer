/* island.c
 * SRE, Thu Apr 21 09:57:21 2005
 * SVN $Id$
 * 
 * Implementation of the "island method" for estimating extreme value
 * distribution parameters. 
 * 
 * Currently experimental. The island method does not appear to work
 * well on glocal alignment score distributions, which are not EVD,
 * although the direct method works acceptably well. The island method
 * does look somewhat promising for low-resolution estimation of
 * Gumbel parameters for fully local (sw or fs) alignment, but these
 * are not HMMER default styles at present, so island estimation is
 * getting shelved for the time being. (SRE, Fri Apr 22 15:08:46 2005;
 * xref STL9/p75).
 * 
 * References: 
 *  [Altschul01] SF Altschul, R Bundschuh, R Olsen, T Hwa. The estimation of
 *               statistical parameters for local alignment score
 *               distributions. Nucleic Acids Research 29:351-361, 2001.
 *               
 *  [Olsen99]    R Olsen, R Bundschuh, T Hwa. Rapid assessment of extremal
 *               statistics for gapped local alignment. Proc Int Conf
 *               Intell Syst Mol Biol 1999; 211-222.            
 *               
 */

#include "config.h"
#include "squidconf.h"

#include "structs.h"
#include "funcs.h"
#include "squid.h"

/* Function: IslandViterbi()
 * Date:     SRE, Thu Apr 21 11:23:36 2005 [St. Louis]
 *
 * Purpose:  The Olsen/Bundschuh/Hwa island algorithm, extended to Plan
 *           7 profile HMM algorithm.
 *           
 *           Every S->B->M entry is counted as an independent island
 *           beginning; every M->E->C^n->T is counted as an
 *           independent island ending. But an island may still
 *           contain multiple hits to the model, passing through the J
 *           state. 
 *           
 *           An unsorted list of island scores is allocated here, and
 *           returned via <*ret_isle_sc>; it contains <*ret_inum>
 *           scores. This array can be quite long (10^5 islands?).
 *           I haven't yet tried to estimate or control its size.
 *           
 *           Also returns the island lengths, defined (following
 *           Altschul) as the average of the aligned consensus model
 *           length M (matches + deletes) and the aligned sequence
 *           length L (residues): (M states + D states + residues) /
 *           2. The caller may want to use this info for edge length
 *           corrections. Island lengths are calculated recursively in
 *           a bit of a tricky way, using a path-info-propagating DP
 *           trick, adding 2 per aligned match, 1 per delete, 1 per
 *           insert.
 *           
 *           Island length might be a little deceptive for glocal
 *           alignment; because of wing folding, it does *not* count
 *           the D's in leading and trailing D->D paths in the model,
 *           because wing folding makes these get counted as internal
 *           M entries and exits. This should be dealt with in the
 *           future if island lengths need to be accurate for glocal
 *           or global alignment styles. 
 *           
 *           Adapted from P7ParsingViterbi().
 *
 * Args:     dsq          - sequence in digitized form
 *           L            - length of dsq
 *           hmm          - the model (log odds scores ready)
 *           ret_isle_sc  - RETURN: unsorted list of island scores
 *           ret_isle_len - RETURN: corresponding island lengths
 *
 * Returns:  Optimal score, log P(S|M)/P(S|R), as a bit score.
 * 
 * Xref:     STL9/75.
 */
float
IslandViterbi(unsigned char *dsq, int L, struct plan7_s *hmm, 
	      int **ret_isle_sc, int **ret_isle_len, int *ret_inum)
{
  struct dpmatrix_s *mx;        /* 2 rows of score matrix */
  struct dpmatrix_s *inmx;      /* propagated island lengths */
  struct dpmatrix_s *itag;	/* propagated island number tags (~ I(r,z))*/
  int  **xmx, **mmx, **dmx, **imx; /* convenience ptrs to score matrix   */  
  int  **xn, **mn, **dn, **in;     /* convenience ptrs to island lengths */
  int  **xi, **mi, **di, **ii;     /* convenience ptrs to island tags    */
  
  int   *isc;		/* isc[i] = max sc of island i (~ \sigma(i) in Hwa) */
  int   *ilen;		/* length of island i */
  int    ialloc;	/* current allocation size for isc */
  int    inum;		/* current # of islands, in isc[0..inum-1] */
  
  int    sc;			/* integer score of optimal alignment  */
  int    i,k,tpos;		/* index for seq, model, trace position */
  int    cur, prv;		/* indices for rolling dp matrix */
  int    curralloc;		/* size of allocation for tr */

  /* Alloc a DP matrix, two rows, linear memory O(M).
   * Alloc the same to hold propagated island lengths, and
   * propagated island tag #.
   * Allocate for the max island scores, and lengths.
   */
  mx    = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  inmx  = AllocPlan7Matrix(2, hmm->M, &xn,  &mn,  &in,  &dn);
  itag  = AllocPlan7Matrix(2, hmm->M, &xi,  &mi,  &ii,  &di);

  ialloc = L;			/* just an initial guess. */
  inum   = 0;
  isc  = MallocOrDie(sizeof(int) * L);
  ilen = MallocOrDie(sizeof(int) * L);

  /* Initialization of the zero row.
   */
  xmx[0][XMN] = 0;		                     /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */
  xi[0][XMB]  = -1;

  /* Recursion. Done as a pull. 
   * Rolling index trick. Trace ptr propagation trick.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   *    
   * Notes on start position and island tag propagation.
   *  - In the path B->E, we propagate the i that B was aligned to 
   *    in the optimal alignment, via mtr, dtr, and itr. 
   *  - When we reach an E, we record the i of the B it started from in etr.
   *  - In a looping path E->J...->B or terminal path E->C...->T, we propagate
   *    the i that E was aligned to in the optimal alignment via xtr[][XMC]
   *    and xtr[][XMJ].
   *  - When we enter B, we record the i of the best previous E, or 0 if there
   *    isn't one, in btr.
   */
  for (i = 1; i <= L; i++) {
    cur = i % 2;
    prv = !cur;
    
    mmx[cur][0] = imx[cur][0] = dmx[cur][0] = -INFTY;

    for (k = 1; k <= hmm->M; k++) {
				/* match state */
      mmx[cur][k] = -INFTY;
      if ((sc = mmx[prv][k-1] + hmm->tsc[TMM][k-1]) > -INFTY) { 
	mmx[cur][k] = sc; 
	mn[cur][k]  = mn[prv][k-1] + 2; 
	mi[cur][k]  = mi[prv][k-1];  
      }
      if ((sc = imx[prv][k-1] + hmm->tsc[TIM][k-1]) > mmx[cur][k]) {
	mmx[cur][k] = sc;
	mn[cur][k]  = in[prv][k-1] + 2;
	mi[cur][k]  = ii[prv][k-1];
      }
      if ((sc = dmx[prv][k-1] + hmm->tsc[TDM][k-1]) > mmx[cur][k]) {
	mmx[cur][k] = sc;
	mn[cur][k]  = dn[prv][k-1] + 2;
	mi[cur][k]  = di[prv][k-1];
      }
      if ((sc = xmx[prv][XMB] + hmm->bsc[k]) > mmx[cur][k]) {
	mmx[cur][k] = sc; 

	if (xi[prv][XMB] == -1) {
	  mi[cur][k]  = inum;
	  mn[cur][k]  = 2;     /* init a new island if B->M is best path. */

	  isc[inum] = -INFTY;
	  inum++;			
	  if (inum == ialloc) { /* reallocation needed? */
	    ialloc *= 2;
	    isc  = ReallocOrDie(isc,  sizeof(int) * ialloc);
	    ilen = ReallocOrDie(ilen, sizeof(int) * ialloc);
	  }
	} else {
	  mi[cur][k] = xi[prv][XMB]; /* propagating a loop */
	  mn[cur][k] = xn[prv][XMB];
	}
      }

      if (hmm->msc[dsq[i]][k] != -INFTY)
	mmx[cur][k] += hmm->msc[dsq[i]][k];
      else
	mmx[cur][k] = -INFTY;

				/* delete state */
      dmx[cur][k] = -INFTY;
      if ((sc = mmx[cur][k-1] + hmm->tsc[TMD][k-1]) > -INFTY) {
	dmx[cur][k] = sc; 
	dn[cur][k]  = mn[cur][k-1] + 1;	
	di[cur][k]  = mi[cur][k-1];
      }
      if ((sc = dmx[cur][k-1] + hmm->tsc[TDD][k-1]) > dmx[cur][k]) {
	dmx[cur][k] = sc; 
	dn[cur][k]  = dn[cur][k-1] + 1;
	di[cur][k]  = di[cur][k-1];
      }

				/* insert state */
      if (k < hmm->M) {
	imx[cur][k] = -INFTY;
	if ((sc = mmx[prv][k] + hmm->tsc[TMI][k]) > -INFTY) {
	  imx[cur][k] = sc; 
	  in[cur][k]  = mn[prv][k] + 1; 
	  ii[cur][k]  = mi[prv][k];
	}
	if ((sc = imx[prv][k] + hmm->tsc[TII][k]) > imx[cur][k]) {
	  imx[cur][k] = sc; 
	  in[cur][k]  = in[prv][k] + 1;
	  ii[cur][k]  = ii[prv][k];
	}
	if (hmm->isc[dsq[i]][k] != -INFTY)
	  imx[cur][k] += hmm->isc[dsq[i]][k];
	else
	  imx[cur][k] = -INFTY;
      }
    }

    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */
				/* N state */
    xmx[cur][XMN] = -INFTY;
    if ((sc = xmx[prv][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
      xmx[cur][XMN] = sc;
				/* E state */
    xmx[cur][XME] = -INFTY;
    for (k = 1; k <= hmm->M; k++)
      {
	if ((sc =  mmx[cur][k] + hmm->esc[k]) > xmx[cur][XME]) {
	  xmx[cur][XME] = sc;
	  xi[cur][XME]  = mi[cur][k];
	  xn[cur][XME]  = mn[cur][k];
	}
	/* calculate what island sc would be, if we ended it here. */
	sc = sc 
	  + hmm->xsc[XTE][MOVE] 
	  + (L-i)* hmm->xsc[XTC][LOOP] 
	  + hmm->xsc[XTC][MOVE];
	if (sc > isc[mi[cur][k]]) {
	  isc[mi[cur][k]]  = sc;
	  ilen[mi[cur][k]] = mn[cur][k] / 2;
	}
      }
				/* J state */
    xmx[cur][XMJ] = -INFTY;
    if ((sc = xmx[prv][XMJ] + hmm->xsc[XTJ][LOOP]) > -INFTY) {
      xmx[cur][XMJ] = sc; 
      xn[cur][XMJ]  = xn[prv][XMJ]; 
      xi[cur][XMJ]  = xi[prv][XMJ];
    }
    if ((sc = xmx[cur][XME]   + hmm->xsc[XTE][LOOP]) > xmx[cur][XMJ]) {
      xmx[cur][XMJ] = sc;
      xi[cur][XMJ]  = xi[cur][XME];
      xn[cur][XMJ]  = xn[cur][XME]; 
    }
				/* B state */
    xmx[cur][XMB] = -INFTY;
    if ((sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY) {
      xmx[cur][XMB] = sc; 
      xi[cur][XMB]  = -1;       /* if coming from N, we are islandless */
      xn[cur][XMB]  = 0;	
    }
    if ((sc = xmx[cur][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[cur][XMB]) {
      xmx[cur][XMB] = sc; 
      xi[cur][XMB]  = xi[cur][XMJ]; /* if from J, then propagate island tag */
      xn[cur][XMB]  = xn[cur][XMJ];
    }
				/* C state */
    xmx[cur][XMC] = -INFTY;
    if ((sc = xmx[prv][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      xmx[cur][XMC] = sc; 
    if ((sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > xmx[cur][XMC])
      xmx[cur][XMC] = sc; 
  }
				/* T state (not stored) */
  sc = xmx[cur][XMC] + hmm->xsc[XTC][MOVE];
  /* sc is the overall optimal score. */

  FreePlan7Matrix(mx);
  FreePlan7Matrix(inmx);
  FreePlan7Matrix(itag);

  *ret_isle_sc  = isc;
  *ret_isle_len = ilen;
  *ret_inum     = inum;
  return Scorify(sc);
}





/************************************************************
 * @LICENSE@
 ************************************************************/
