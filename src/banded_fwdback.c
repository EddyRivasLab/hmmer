/* Banded Forward/Backward, posterior decoding, and optimal accuracy alignment,
 * with a dual-mode glocal/local model
 */

#include "easel.h"

#include "hmmer.h"
#include "p7_bandmx.h"

int 
p7_BandedForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_BANDMX *bmx, float *opt_sc)
{
  int         *bnd_ip = bmx->bnd->imem;          /* ptr to current ia, ib segment band in bmx->bnd  */
  int         *bnd_kp = bmx->bnd->kmem;		 /* ptr to current ka, kb row band in bmx->bnd      */
  float       *dpc    = bmx->dp1;	         /* ptr to current DP matrix cell; <dp1>=Fwd matrix */
  float       *xpc    = bmx->xmx1;		 /* ptr to current special cell   */
  float const *tsc    = gm->tsc;		 /* sets up TSC() macro, access to profile's transitions */
  float const *rsc;				 /* will be set up for MSC(), ISC() macros for residue scores */
  float       *dpp;	                  	 /* ptr to previous DP matrix cell */
  float       *last_dpc;			 /* used to reinitialize dpp after each row        */
  int          ia, ib;				 /* current segment band is rows ia..ib            */
  int          last_ib;				 /* intersegment interval is last_ib+1..ia-1       */
  int          kac, kbc;			 /* current row band is kac..kbc                   */
  int          kap, kbp;			 /* previous row band is kap..kbp                  */
  float        xE, xN, xJ, xB, xL, xG, xC;       /* tmp scores on special states. only stored when in row bands */
  float        mlp, mgp, ilp, igp, dlp, dgp;     /* M,I,D cell values from previous row i-1     */
  float        dlc, dgc;			 /* precalculated D(i,k+1) value on current row */
  float        mlc, mgc;			 /* temporary score calculations M(i,k)         */
  int          g, i, k;				 /* indices running over segments, residues (rows) x_i, model positions (cols) k  */
  
  xN      = 0.0f;
  xJ      = -eslINFINITY;
  xC      = -eslINFINITY;
  last_ib = 0;

  for (g = 0; g < bmx->bnd->nseg; g++)
    {
      ia = *bnd_ip++;
      ib = *bnd_ip++;

      /* kap,kbp initialization for i=ia:
       *  left overhang dpp advance must always eval to 0, 
       *  {m,i,d}vp initialization must always eval to -eslINFINITY.
       */
      kap = kbp = gm->M+1;   
      dpp = dpc;		/* re-initialize dpp */
      
      /* re-initialization: specials for previous row ia-1 just outside banded segment:  */
      xE  = -eslINFINITY;
      xN  = xN + (ia - last_ib - 1) * gm->xsc[p7P_N][p7P_LOOP];
      xJ  = xJ + (ia - last_ib - 1) * gm->xsc[p7P_J][p7P_LOOP];
      xB  = p7_FLogsum( xN + gm->xsc[p7P_N][p7P_MOVE], xJ + gm->xsc[p7P_J][p7P_MOVE]);
      xL  = xB + gm->xsc[p7P_B][0]; /* B->L */
      xG  = xB + gm->xsc[p7P_B][1]; /* B->G */
      xC  = xC + (ia - last_ib - 1) * gm->xsc[p7P_C][p7P_LOOP];

      for (i = ia; i <= ib; i++)
	{
	  rsc       = gm->rsc[dsq[i]];   /* sets up MSC(k), ISC(k) residue scores for this row i */
	  dlc = dgc = -eslINFINITY;
	  xE        = -eslINFINITY;
	  last_dpc  = dpc;

	  kac      = *bnd_kp++;         /* current row's band is cells k=kac..kbc  */
	  kbc      = *bnd_kp++; 

	  /* dpp must advance by any left overhang of previous row; but no more than the entire row */
	  dpp += (kac-1 > kap ? ESL_MIN(kac-kap-1, kbp-kap+1) : 0);

	  if (kac > kap && kac-1 <= kbp) { mlp = *dpp++; mgp = *dpp++; ilp = *dpp++; igp = *dpp++; dlp = *dpp++; dgp = *dpp++;       }
	  else                           { mlp =         mgp =         ilp =         igp =         dlp =         dgp = -eslINFINITY; }

	  for (k = kac; k <= kbc; k++)
	    {
	      *dpc++ = mlc = MSC(k) + p7_FLogsum( p7_FLogsum(mlp + TSC(p7P_MM, k-1), 
							     ilp + TSC(p7P_IM, k-1)),
						  p7_FLogsum(dlp + TSC(p7P_DM, k-1),
							     xL  + TSC(p7P_LM, k-1)));

	      *dpc++ = mgc = MSC(k) + p7_FLogsum( p7_FLogsum(mgp + TSC(p7P_MM, k-1), 
							     igp + TSC(p7P_IM, k-1)),
						  p7_FLogsum(dgp + TSC(p7P_DM, k-1),
							     xG  + TSC(p7P_GM, k-1))); /* wing retracted glocal entry */

	      /* This "if" seems unavoidable. */
	      if (k >= kap && k <= kbp) { mlp = *dpp++; mgp = *dpp++; ilp = *dpp++; igp = *dpp++; dlp = *dpp++; dgp = *dpp++;       }
	      else                      { mlp =         mgp =         ilp =         igp =         dlp =         dgp = -eslINFINITY; }

	      *dpc++ = ISC(k) + p7_FLogsum( mlp + TSC(p7P_MI, k), ilp + TSC(p7P_II, k)); /* IL */     // correctly sets -inf at k=M boundary, because TSC() is -inf
	      *dpc++ = ISC(k) + p7_FLogsum( mgp + TSC(p7P_MI, k), igp + TSC(p7P_II, k)); /* IG */     // ditto

	      xE = p7_FLogsum( p7_FLogsum(mlc, dlc),
			       p7_FLogsum(xE,  mgc + TSC(p7P_MGE, k))); 
			       
	      /* Delayed store of Dk; advance calculation of next D_k+1 */
	      *dpc++ = dlc;
	      *dpc++ = dgc;
	      dlc    = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));	          // dlc=-inf at k=M boundary, doesn't matter, won't be used.
	      dgc    = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));	          // ditto 
	    }

	  *xpc++ = xE;
	  *xpc++ = xN = xN + gm->xsc[p7P_N][p7P_LOOP];
	  *xpc++ = xJ = p7_FLogsum( xJ + gm->xsc[p7P_J][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_MOVE]);
	  *xpc++ = xB = p7_FLogsum( xJ + gm->xsc[p7P_J][p7P_MOVE],  xN + gm->xsc[p7P_N][p7P_MOVE]);
	  *xpc++ = xL = xB + gm->xsc[p7P_B][0]; 
	  *xpc++ = xG = xB + gm->xsc[p7P_B][1]; 
	  *xpc++ = xC = p7_FLogsum( xE + gm->xsc[p7P_E][p7P_MOVE],  xC + gm->xsc[p7P_C][p7P_LOOP]);

	  dpp = last_dpc;	/* this skips any right overhang on the previous row, so dpp advances (if necessary) to start of curr row */
	  kap = kac;
	  kbp = kbc;
	}
      last_ib = ib;
    }

  /* last_ib+1..L is outside any band segment, so it can only run through xC. */
  if (opt_sc != NULL) *opt_sc = xC + (L-last_ib) *  gm->xsc[p7P_C][p7P_LOOP] + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}

int
p7_BandedBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_BANDMX *bmx, float *opt_sc)
{
  int   *bnd_ip;		/* ptr into <bnd>'s array of segments ia..ib coords                   */
  int   *bnd_kp;		/* ptr into <bnd>'s array of <ka..kb> band extent on each banded row  */
  float *dpc;			/* ptr into current row's memory; starts at end, decrements; kbc..kac */
  float *dpn;			/* ptr into next row's memory; decrements from end; kbn..kan          */
  float *rsc;			/* ptr into emission scores on current row i                          */
  float *rsn;			/* ptr into emission scores on next row i+1                           */
  float *xc;			/* ptr into current row's special states in <bmx>                     */
  int    kac, kbc;		/* band extent on current row i is kac..kbc                           */
  int    kan, kbn;		/* band extent on next row, i+1, is kan..kbn                          */
  int    k, ia, ib, i, g;	/* indices used for traversing banded rows, segments in <bnd>         */
  int    last_ia;		/* lets us calculate gap size between segments ..ib] x x x [last_ia.. */
  float *last_dpc;              /* as we finish each row i, lets us set <dpn> where <dpc> started     */
  float *tsc;			/* ptr into <gm>'s transition scores, enabling TSC() macro            */
  float  xC,xG,xL,xB,xJ,xN,xE;	/* tmp vars for special cells on this row i                           */
  float  dlc,dgc;		/* tmp vars for D vals on current row i                               */
  float  mln,mgn,iln,ign;	/* tmp vars for M+e, I+e on next row i+1                              */
  float  mlc,mgc;		/* tmp vars for M(i,k) on current row i                               */

  tsc = gm->tsc;

  /* In Backwards, various pointers get set to the _end_ of arrays, which we traverse backwards            */
  bnd_ip  = bmx->bnd->imem + (           2 * bmx->bnd->nseg - 1);   /* set to _last_ array elem            */
  bnd_kp  = bmx->bnd->kmem + (p7_GBANDS_NK * bmx->bnd->nrow - 1);   /* ditto                               */
  dpc     = bmx->dp2  + (p7B_NSCELLS * bmx->bnd->ncell - 1);        /* Backwards scores go to <dp2> matrix */
  xc      = bmx->xmx2 + (p7B_NXCELLS * bmx->bnd->nrow - 1);         /* last elem in xmx2                   */
  last_ia = L+1;		/* assure a boundary condition (last_ia - L - 1) = 0, when first ib==L     */

  /* Traversing banded rows takes two loops: <g> over segments ia..ib, then <i> over ib..ia rows */
  for (g = bmx->bnd->nseg-1; g >= 0; g--)
    {
      ib = *bnd_ip--;		/* note ib picked up before ia, because we're accessing in reverse */
      ia = *bnd_ip--;

      /* kan..kbn initialization, where we have no next row i+1:
       *   right overhang dpn advance must always eval to 0;
       *   {m,i,d}{lg}n initialization must evaluate to -eslINFINITY;
       *   kbn ? test must eval to FALSE if there's no next row 
       */
      kan = kbn = 0;
      dpn = last_dpc = dpc;

      /* reinitialization of the specials, for next row ib+1 just outside
       * the banded segment we're about to start. The number of residues
       * in the intersegment gap is (last_ia - ib - 1). They must be accounted
       * for by that many C->C or J->J transitions.
       *     ib] ib+1 ... ... ... [last_ia ... ... ]
       *           <----------------- xC
       *             x   x   x   x   : last_ia - (ib+1) 
       * note the boundary condition, if ib==L, (last_ia - ib - 1) must be 0; this
       * dictates initialization of last_ia to L+1.
       */
      xC = xC + (last_ia - ib - 1) * gm->xsc[p7P_C][p7P_LOOP];
      xJ = xJ + (last_ia - ib - 1) * gm->xsc[p7P_J][p7P_LOOP];
      xG = xL = xB = xN = xE = -eslINFINITY;

      /* xE was reachable from xJ, xC on prev row ib+1, but this will
       * have no effect when we calculate first row i=ib below, so may
       * as well set it to -eslINFINITY instead of wasting calculation
       * on it.
       */
      for (i = ib; i >= ia; i--)
	{
	  dlc = dgc = -eslINFINITY;    /* initialize: D(i,kb+1) is out of band, so -inf */

	  kbc = *bnd_kp--;	       /* pick up kac..kbc band on row i; b before a    */
	  kac = *bnd_kp--;

	  /* <rsc>, <rsn> will always point to the M emission score, and will always decrement by 2 (p7P_NR).
	   * When we want the I emission score, we will access *(rsc+1), *(rsc+p7P_I)
	   */
	  rsc =        gm->rsc[dsq[i]]   + p7P_NR * kbc;;              /* <rsc> now on last e(Mk,i) on current row, and ready to decrement  */
	  rsn = (kbn ? gm->rsc[dsq[i+1]] + p7P_NR * (kbc+1) : NULL);   /* <rsn> now on last e(Mk+1,i+1) on next row, and ready to decrement */

	  /* dpn initialization removes any excess right overhang of previous row,
           * but no more than the entire previous row.
	   *    ... kbc]    current row
           *        kbn]                   : do nothing, ok, dpp decrement will get delayed
           *        ... kbn]               : do nothing, this is just right
           *        ...  ooo kbn]          : -(kbn-kbc-1), and dpp points at [ooo]
           *             xxx [kan ... kbn] : -(kbn-kan+1), and dpp points at [kbc] on cur row; dpp decrement won't happen
	   */
	  dpn -= (kbn > kbc+1 ?  : p7B_NSCELLS * ESL_MIN(kbn-kbc-1, kbn-kan+1), 0);

	  /* if k+1 exists on next row, pick up mln, mgn, and move dpp back one supercell.
	   * i.e.:    kbc]      but not   kbc]
           *          ... kbn]                 xxx [kan ...]
	   * -4 = p7G_NSCELLS-P7G_MG+1;  -5 = p7G_NSCELLS-P7G_ML+1
	   */
	  if (kbn>kbc && kbc+1 >= kan) { mgn = *(dpn-4) + *rsn; mln = *(dpn-5) + *rsn; dpn -= p7G_NSCELLS;  rsn -= p7P_NR; } 
	  else                         { mgn = mln = -eslINFINITY; }
	  
	  /* ditto for k, now without moving dpn; dpn stays on k=kbc in prv row
	   * -2 = p7G_NSCELLS-p7G_IG+1; -3 = p7G_NSCELLS-p7G_IL+1 
	   */
	  if (kbn >= kbc && kbc >= kan) { ign = *(dpn-2) + *(rsn+p7P_I); iln = *(dpn-3) + *(rsn+p7P_I); } /* <rsn> stays on last e(Mk,i+1) - we'll need to pick up M score */
	  else                          { ign = iln = -eslINFINITY;       }

	  *xc-- = xC = xC + gm->xsc[p7P_C][p7P_LOOP];    // 
	  *xc-- = xG;					 // delayed store; was calculated on prv row i+1
	  *xc-- = xL;					 // delayed store, ditto
	  *xc-- = xB = p7_FLogsum( xL + gm->xsc[p7P_B][0],
				   xG + gm->xsc[p7P_B][1]);
	  *xc-- = xJ = p7_FLogsum( xJ + gm->xsc[p7P_J][p7P_LOOP],
				   xB + gm->xsc[p7P_J][p7P_MOVE]);
	  *xc-- = xN = p7_FLogsum( xN + gm->xsc[p7P_N][p7P_LOOP],
				   xB + gm->xsc[p7P_N][p7P_MOVE]);
	  *xc-- = xE = p7_FLogsum( xJ + gm->xsc[p7P_E][p7P_LOOP],
				   xC + gm->xsc[p7P_E][p7P_MOVE]);

	  for (k = kbc; k >= kac; k--)
	    {
	      /* M(i,k) calculations need to use dgc,dlc before we
	       * change them, while they're still D(i,k+1). But we
	       * can't store them just yet; <dpc> is pointing at
	       * DG state right now. 
	       */
	      mgc = p7_FLogsum( p7_FLogsum(mgn + TSC(p7P_MM,  k),    // Mg(i,k) =   Mg(i+1,k+1) * t(k,MM) * e_M(k+1, x_i+1)      | mgn = log[ Mg(i+1,k+1) * e_M(k+1,x_i+1)]
					   ign + TSC(p7P_MI,  k)),   //           + Ig(i+1,k)   * t(k,MI) * e_I(k, x_i+1)        | ign = log[ Ig(i+1,k)   * e_I(k,  x_i+1)]
				p7_FLogsum(dgc + TSC(p7P_MD,  k),    //           + Dg(i,  k+1) * t(k,DD)                        | dgc = Dg(i,k+1), wrapped around from prev loop iteration
					   xE  + TSC(p7P_MGE, k)));  //           + E(i)        * \prod_j=k^M-1 t(j,DD)          | t(k,MGE) = right wing retraction, log[ \prod_j=k^M-1 t(j,DD) ]
	      mlc = p7_FLogsum( p7_FLogsum(mln + TSC(p7P_MM,  k),    // Ml(i,k) is essentially the same recursion...
					   iln + TSC(p7P_MI,  k)),   //
				p7_FLogsum(dlc + TSC(p7P_MD,  k),    //
					   xE));                     // ... with the difference that the Mk->E transition prob in local alignment is 1.0
	      
	      /* Accumulate xG, xL as we sweep over the current row; 
	       * they get used to initialize the next row (i-1) 
	       */
	      xG = p7_FLogsum( xG, mgc + *rsc + TSC(p7P_GM, k-1));   // t(k-1,p7P_GM) = left wing retraction, entry at Mk, stored off-by-one at k-1 
	      xL = p7_FLogsum( xL, mlc + *rsc + TSC(p7P_LM, k-1));   // t(k-1,p7P_LM) = uniform local entry at Mk, off-by-one storage at k-1 
	      rsc -= p7P_NR;                                         // skip the e_I(k,x_i) insert emission; here we're concerned only with {LG}->Mk entries

	      /* D(i,k) calculations; 
	       * we can start storing results in <dpc> and decrementing;
	       * and we're done with dgc,dlc, so can update them for
               * time thru the k loop.
               * Dg has no exit to E; ...Dm->E paths have been wrapped into the right wing retraction in t(k,MGE).
	       */
	      *dpc-- = dgc = p7_FLogsum(            mgn  + TSC(p7P_DM, k),  // Dg(i,k) =  Mg(i+1,k+1) * t(k,DM) * e_M(k+1, x_i+1)
					            dgc  + TSC(p7P_DD, k)); //          + Dg(i,  k+1) * t(k,DD)  
	      *dpc-- = dlc = p7_FLogsum( p7_FLogsum(mln  + TSC(p7P_DM, k),  // Dl(i,k) is essentially the same recursion...
						    dlc  + TSC(p7P_DD, k)), //        
                 				    xE);                    // ... except it also includes Dk->E local exit transition, always 1.0

	      /* I(i,k) calculations */
	      *dpc-- = p7_FLogsum( mgn + TSC(p7P_IM, k),    // Ig(i,k) = Mg(i+1,k+1) * t(k,IM) * e_M(k+1, x_i+1)
				   ign + TSC(p7P_II, k));   //         + Ig(i+1,k)   * t(k,II) * e_I(k,   x_i+1)
	      *dpc-- = p7_FLogsum( mln + TSC(p7P_IM, k),    // Il(i,k) is essentially the same
				   iln + TSC(p7P_II, k));

	      /* Two-step pickup:
               * get mgn,mln from current k on next row i+1; decrement <dpn>; then pick up ign,iln.
	       * Include the emission scores; leave <rsn> on current k, where it will become k+1 as we loop around.
	       * Thus, as we roll around and reenter the k loop, after decrementing k:
               *     mgn will be M(i+1,k+1) + e(Mk+1, xi+1), and
	       *     ign will be I(i+1,k)   + e(Ik, xi+1)
	       */
	      if (kbn > kbc && kbc+1 >= kan) { mgn = *(dpn-4) + *rsn; mln = *(dpn-5) + *rsn; dpn -= p7G_NSCELLS; rsn -= p7P_NR; } 
	      else                           { mgn = mln = -eslINFINITY; }
	      if (kbn >= kbc && kbc >= kan)  { ign = *(dpn-2) + *(rsn+p7P_I); iln = *(dpn-3) + *(rsn+p7P_I); }
	      else                           { ign = iln = -eslINFINITY;       }
	      
	      /* Store mgc, mlc. */
	      *dpc-- = mgc;
	      *dpc-- = mlc;
	    }

	  kbn = kbc;
	  kan = kac;
	  dpn = last_dpc;	/* skips any remaining right overhang on row i+1 */
	}
      last_ia = ia;
    }

  /* 1..last_ia-1 is outside any band segment; (last_ia-1) residues that must go thru N; S->N transition prob = 1.0 */
  if (opt_sc) *opt_sc = xN + (last_ia-1) * gm->xsc[p7P_N][p7P_LOOP];
  return eslOK;
}

int 
p7_BandedDecoding(const P7_PROFILE *gm, P7_BANDMX *bmx)
{
  return eslOK;
}


int 
p7_BandedAlignment(const P7_PROFILE *gm, P7_BANDMX *bmx, P7_TRACE *tr)
{
  return eslOK;

}
