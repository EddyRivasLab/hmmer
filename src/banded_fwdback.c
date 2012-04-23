/* Banded Forward/Backward, posterior decoding, and optimal accuracy alignment,
 * with a dual-mode glocal/local model.
 * 
 * Contents:
 *    1. Banded Forward
 *    2. Banded Backward
 *    3. Banded posterior decoding
 *    4. Banded alignment (MEG, gamma-centroid)
 *    5. Banded alignment traceback 
 *    x. Example
 *    x. Copyright and license.
 *    
 * See also:
 *    p7_bandmx.[ch] : P7_BANDMX banded DP matrix structure used here.   
 */

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_bandmx.h"

static int traceback(const P7_PROFILE *gm, const P7_BANDMX *bmd, const P7_BANDMX *bma, P7_TRACE *tr);

/*****************************************************************
 * 1. Banded Forward
 *****************************************************************/

/* Function:  p7_BandedForward()
 * Synopsis:  Forward algorithm, in banded DP
 *
 * Purpose:   Compute the Forward algorithm for comparing query profile <gm>
 *            to digitized target sequence <dsq> of length <L>, using the
 *            caller-allocated banded DP matrix space provided in <bmf>.
 *            
 *            Upon return, <bmf> contains the computed Forward matrix,
 *            and <opt_sc> optionally contains the Forward raw lod
 *            score in nats.
 *
 * Args:      dsq    - target sequence, digital, 1..L
 *            L      - length of <dsq> in residues
 *            gm     - query profile
 *            bmf    - RESULT:    Forward matrix (caller allocates the space)
 *            opt_sc - optRETURN: Forward raw lod score, nats
 *
 * Returns:   <eslOK> on success, and <*opt_sc> is the Forward score,
 *            and <bmf> contains the computed Forward matrix.
 *
 * Throws:    (no abnormal error conditions)
 */
int 
p7_BandedForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_BANDMX *bmf, float *opt_sc)
{
  int         *bnd_ip = bmf->bnd->imem;          /* ptr to current ia, ib segment band in bmx->bnd  */
  int         *bnd_kp = bmf->bnd->kmem;		 /* ptr to current ka, kb row band in bmx->bnd      */
  float       *dpc    = bmf->dp;	         /* ptr to current DP matrix cell */
  float       *xpc    = bmf->xmx;		 /* ptr to current special cell   */
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

  for (g = 0; g < bmf->bnd->nseg; g++)
    {
      ia = *bnd_ip++;
      ib = *bnd_ip++;

      /* kap,kbp initialization for i=ia:
       *  left overhang dpp advance must always eval to 0, 
       *  {m,i,d}vp initialization must always eval to -eslINFINITY.
       */
      kap = kbp = gm->M+1;   
      dpp = dpc;		/* re-initialize dpp */
      
      /* re-initialization: specials for previous row ia-1 just outside banded segment; these are stored!       */
      /* in general, ia > last_ib+1 (else ia would have merged w/ prv segment). Exception: ia==1 at the first segment, where last_ib==0 as a boundary condition */
      *xpc++ = xE  = -eslINFINITY;
      *xpc++ = xN  = xN + ( ia == last_ib+1 ? 0.0f : (ia - last_ib - 1) * gm->xsc[p7P_N][p7P_LOOP]); /* watch out for 0*-inf special case */
      *xpc++ = xJ  = xJ + ( ia == last_ib+1 ? 0.0f : (ia - last_ib - 1) * gm->xsc[p7P_J][p7P_LOOP]);
      *xpc++ = xB  = p7_FLogsum( xN + gm->xsc[p7P_N][p7P_MOVE], xJ + gm->xsc[p7P_J][p7P_MOVE]);
      *xpc++ = xL  = xB + gm->xsc[p7P_B][0]; /* B->L */
      *xpc++ = xG  = xB + gm->xsc[p7P_B][1]; /* B->G */
      *xpc++ = xC  = xC + ( ia == last_ib+1 ? 0.0f : (ia - last_ib - 1) * gm->xsc[p7P_C][p7P_LOOP]);
      *xpc++       = -eslINFINITY; /* JJ: this space only used in a Decoding matrix. */
      *xpc++       = -eslINFINITY; /* CC: this space only used in a Decoding matrix. */

      for (i = ia; i <= ib; i++)
	{
	  rsc       = gm->rsc[dsq[i]];   /* sets up MSC(k), ISC(k) residue scores for this row i */
	  dlc = dgc = -eslINFINITY;
	  xE        = -eslINFINITY;
	  last_dpc  = dpc;

	  kac      = *bnd_kp++;         /* current row's band is cells k=kac..kbc  */
	  kbc      = *bnd_kp++; 

	  /* dpp must advance by any left overhang of previous row; but no more than the entire row */
	  dpp += (kac-1 > kap ? ESL_MIN(kac-kap-1, kbp-kap+1) * p7B_NSCELLS : 0);

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

	      *dpc++ = ISC(k) + p7_FLogsum( mlp + TSC(p7P_MI, k), ilp + TSC(p7P_II, k)); /* IL */     // correctly sets -inf at k=M boundary, because TSC(II,M) is set to -inf
	      *dpc++ = ISC(k) + p7_FLogsum( mgp + TSC(p7P_MI, k), igp + TSC(p7P_II, k)); /* IG */     // ditto

	      xE = p7_FLogsum( xE, p7_FLogsum(mlc, dlc)); /* local exit paths */
			       
	      /* Delayed store of Dk; advance calculation of next D_k+1 */
	      *dpc++ = dlc;
	      *dpc++ = dgc;
	      dlc    = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));	          // dlc=-inf at k=M boundary, doesn't matter, won't be used.
	      dgc    = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));	          // ditto 
	    }

	  *xpc++ = xE = p7_FLogsum( xE, dgc + TSC(p7P_DGE, kbc));   // glocal exit path(s) added on, from D(kbc+1)->E
	  *xpc++ = xN = xN + gm->xsc[p7P_N][p7P_LOOP];
	  *xpc++ = xJ = p7_FLogsum( xJ + gm->xsc[p7P_J][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_LOOP]);
	  *xpc++ = xB = p7_FLogsum( xJ + gm->xsc[p7P_J][p7P_MOVE],  xN + gm->xsc[p7P_N][p7P_MOVE]);
	  *xpc++ = xL = xB + gm->xsc[p7P_B][0]; 
	  *xpc++ = xG = xB + gm->xsc[p7P_B][1]; 
	  *xpc++ = xC = p7_FLogsum( xE + gm->xsc[p7P_E][p7P_MOVE],  xC + gm->xsc[p7P_C][p7P_LOOP]);
	  *xpc++      = -eslINFINITY; /* JJ: this space only used in a Decoding matrix. */
	  *xpc++      = -eslINFINITY; /* CC: this space only used in a Decoding matrix. */

	  dpp = last_dpc;	/* this skips any right overhang on the previous row, so dpp advances (if necessary) to start of curr row */
	  kap = kac;
	  kbp = kbc;
	}
      last_ib = ib;
    }

  /* last_ib+1..L is outside any band segment, so it can only run through xC. */
  if (opt_sc != NULL) *opt_sc = xC + ( L-last_ib ? (L-last_ib) *  gm->xsc[p7P_C][p7P_LOOP] : 0.0f) + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}

/*****************************************************************
 * 2. Banded Backward
 *****************************************************************/

/* Function:  p7_BandedBackward()
 * Synopsis:  Backward algorithm, in banded DP
 *
 * Purpose:   Compute the Backward algorithm for comparing query profile <gm>
 *            to digitized target sequence <dsq> of length <L>, using the
 *            caller-allocated banded DP matrix space provided in <bmb>.
 *            
 *            Upon return, <bmb> contains the computed Backward matrix,
 *            and <opt_sc> optionally contains the Backward raw lod
 *            score in nats.
 *
 * Args:      dsq    - target sequence, digital, 1..L
 *            L      - length of <dsq> in residues
 *            gm     - query profile
 *            bmf    - RESULT:    Backward matrix (caller allocates the space)
 *            opt_sc - optRETURN: Backward raw lod score, nats
 *
 * Returns:   <eslOK> on success, and <*opt_sc> is the Backward score,
 *            and <bmf> contains the computed Backward matrix.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_BandedBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_BANDMX *bmb, float *opt_sc)
{
  int   *bnd_ip;		/* ptr into <bnd>'s array of segments ia..ib coords                   */
  int   *bnd_kp;		/* ptr into <bnd>'s array of <ka..kb> band extent on each banded row  */
  float *dpc;			/* ptr into current row's memory; starts at end, decrements; kbc..kac */
  float *dpn;			/* ptr into next row's memory; decrements from end; kbn..kan          */
  float *rsc;			/* ptr into emission scores on current row i                          */
  float *rsn;			/* ptr into emission scores on next row i+1                           */
  float *xc;			/* ptr into current row's special states in <bmb>                     */
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
  bnd_ip  = bmb->bnd->imem + (           2 * bmb->bnd->nseg - 1);             /* set to _last_ array elem            */
  bnd_kp  = bmb->bnd->kmem + (p7_GBANDS_NK * bmb->bnd->nrow - 1);             /* ditto                               */
  dpc     = bmb->dp  + (p7B_NSCELLS * bmb->bnd->ncell - 1);                   /* Backwards scores go to <dp> matrix  */
  xc      = bmb->xmx + (p7B_NXCELLS * (bmb->bnd->nrow+bmb->bnd->nseg) - 1);   /* last elem in xmx                    */

  /* Initialization condition for xC/xN/xJ follows from a convention
   * that as we start row i, they already include all contribution
   * from row i+1... as we start, that's L+1.
   */
  xC      = gm->xsc[p7P_C][p7P_MOVE];
  xJ      = -eslINFINITY;
  xN      = -eslINFINITY;
  last_ia = L+1;		/* assure a boundary condition (last_ia - L - 1) = 0, when first ib==L     */

  /* Traversing banded rows takes two loops: <g> over segments ia..ib, then <i> over ib..ia rows */
  for (g = bmb->bnd->nseg-1; g >= 0; g--)
    {
      ib = *bnd_ip--;		/* note ib picked up before ia, because we're accessing in reverse */
      ia = *bnd_ip--;

      /* kan..kbn initialization, where we have no next row i+1:
       *   right overhang dpn advance must always eval to 0;
       *   {m,i,d}{lg}n initialization must evaluate to -eslINFINITY;
       *   kbn ? test must eval to FALSE if there's no next row 
       */
      kan = kbn = 0;
      dpn = dpc;

      /* reinitialization of the specials, for this row i=ib as we
       * start a banded segment. We last left xN/xJ/xC set for row
       * (last_ia-1). The number of residues in the intersegment
       * gap is (last_ia - ib - 1). Because we defined segments to be
       * separated by at least 1 row, (last_ia-ib-1) is at least 1,
       * except when we're at ib=L where it's 0. We have to account
       * for (last_ia-ib-1) transitions (CC,NN,JJ) to reset the
       * xC/xN/xJ values for row ib. To avoid a pathological case
       * where ib=L, last_ia-ib-1, and a loop penalty happens to be
       * -inf, we have to avoid multiplying 0*-inf (=NaN), so we 
       * wrap in a (last_ia-ib-1 > 0) test.
       */
      if (last_ia-ib-1) {
	xC += (last_ia - ib - 1) * gm->xsc[p7P_C][p7P_LOOP];
	xJ += (last_ia - ib - 1) * gm->xsc[p7P_J][p7P_LOOP];
	xN += (last_ia - ib - 1) * gm->xsc[p7P_N][p7P_LOOP];
      }
      xG = xL = xB = xE = -eslINFINITY;

      for (i = ib; i >= ia; i--)
	{
	  kbc      = *bnd_kp--;	       /* pick up kac..kbc band on row i; b before a    */
	  kac      = *bnd_kp--;
	  last_dpc = dpc;

	  /* <rsc>, <rsn> will always point to the M emission score, and will always decrement by 2 (p7P_NR).
	   * When we want the I emission score, we will access *(rsc+1), *(rsc+p7P_I)
	   */
	  rsc =        gm->rsc[dsq[i]]   + p7P_NR * kbc;                          /* <rsc> now on last e(Mk,i) on current row, and ready to decrement  */
	  rsn = (kbn ? gm->rsc[dsq[i+1]] + p7P_NR * ESL_MIN(kbn,kbc+1) : NULL);   /* <rsn> now on last e(Mk+1,i+1) (or e(Mk,i+1)) on next row, and ready to decrement */

	  /* dpn, rsn initialization now removes any excess right overhang of previous row,
           * but no more than the entire previous row.
	   *    ... kbc]    current row
           *        kbn]                   : do nothing, ok, dpp decrement will get delayed
           *        ... kbn]               : do nothing, this is just right
           *        ...  ooo kbn]          : -(kbn-kbc-1), and dpp points at [ooo]
           *             xxx [kan ... kbn] : -(kbn-kan+1), and dpp points at [kbc] on cur row; dpp decrement won't happen
	   */
	  dpn -= (kbn > kbc+1 ? p7B_NSCELLS * ESL_MIN(kbn-kbc-1, kbn-kan+1) : 0);

	  /* if k+1 exists on next row, pick up mln, mgn, and move dpp back one supercell.
	   * i.e.:    kbc]      but not   kbc]
           *          ... kbn]                 xxx [kan ...]
	   * -4 = p7B_NSCELLS-P7G_MG+1;  -5 = p7B_NSCELLS-P7G_ML+1
	   */
	  if (kbn>kbc && kbc+1 >= kan) { mgn = *(dpn-4) + *rsn; mln = *(dpn-5) + *rsn; dpn -= p7B_NSCELLS;  rsn -= p7P_NR; } 
	  else                         { mgn = mln = -eslINFINITY; }
	  
	  /* xC,xJ,xN,xL,xG already contain their contribution from next row i+1 */
	  *xc-- = -eslINFINITY;	/* CC only stored in a Decoding matrix. */
	  *xc-- = -eslINFINITY; /* JJ only stored in a Decoding matrix. */
	  *xc-- = xC;		
	  *xc-- = xG;					 // delayed store; was calculated on prv row i+1
	  *xc-- = xL;					 // delayed store, ditto
	  *xc-- = xB = p7_FLogsum( xL + gm->xsc[p7P_B][0],
				   xG + gm->xsc[p7P_B][1]);
	  *xc-- = xJ = p7_FLogsum( xJ,
				   xB + gm->xsc[p7P_J][p7P_MOVE]);
	  *xc-- = xN = p7_FLogsum( xN,
				   xB + gm->xsc[p7P_N][p7P_MOVE]);
	  *xc-- = xE = p7_FLogsum( xJ + gm->xsc[p7P_E][p7P_LOOP],
				   xC + gm->xsc[p7P_E][p7P_MOVE]);

	  xC += gm->xsc[p7P_C][p7P_LOOP];
	  xJ += gm->xsc[p7P_J][p7P_LOOP];
	  xN += gm->xsc[p7P_N][p7P_LOOP];
	  xG  = -eslINFINITY;
	  xL  = -eslINFINITY;

	  dlc = -eslINFINITY;             /* initialize: DL(i,kb+1) is out of band, so -inf */
	  dgc = xE + TSC(p7P_DGE,kbc);    /*      whereas the DG(kbc) is reachable on a glocal path. This sets dgc to D(i,k+1), Dk+1->...->Dm->E wing retracted prob. At kbc=M, TSC(DGE,M) = 0.  */
	  for (k = kbc; k >= kac; k--)
	    {
	      /* Pick up iln, ign from (k,i+1) 
	       * -2 = p7B_NSCELLS-p7G_IG+1; -3 = p7B_NSCELLS-p7G_IL+1 
               * dpn and rsn both stay on k cell for now
	       */
	      if (k >= kan && k <= kbn) { ign = *(dpn-2) + *(rsn+p7P_I); iln = *(dpn-3) + *(rsn+p7P_I); } /* <rsn> stays on last e(Mk,i+1) - we'll need to pick up M score */
	      else                      { ign = iln = -eslINFINITY;       }

	      /* M(i,k) calculations need to use dgc,dlc before we
	       * change them, while they're still D(i,k+1). But we
	       * can't store them just yet; <dpc> is pointing at
	       * DG state right now. 
	       */
	      mgc = p7_FLogsum( p7_FLogsum(mgn + TSC(p7P_MM, k),                      // Mg(i,k) =   Mg(i+1,k+1) * t(k,MM) * e_M(k+1, x_i+1)      | mgn = log[ Mg(i+1,k+1) * e_M(k+1,x_i+1)]
					   ign + TSC(p7P_MI, k)),                     //           + Ig(i+1,k)   * t(k,MI) * e_I(k, x_i+1)        | ign = log[ Ig(i+1,k)   * e_I(k,  x_i+1)]
				           dgc + TSC(p7P_MD, k));                     //           + Dg(i,  k+1) * t(k,DD)                        | dgc = Dg(i,k+1), wrapped around from prev loop iteration. This has tricky boundary conditions at kbc=k, and kbc=k=M
	      mlc = p7_FLogsum( p7_FLogsum(mln + TSC(p7P_MM, k),                      // Ml(i,k) is essentially the same recursion...
					   iln + TSC(p7P_MI, k)),                     //
				p7_FLogsum(dlc + TSC(p7P_MD, k),                      //
					   xE));                                      // ... with the difference that there is also a Mk->E transition prob in local alignment, always 1.0
	      
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
					            dgc  + TSC(p7P_DD, k)); //          + Dg(i,  k+1) * t(k,DD)                    | tricky boundary condition: at k=kbc, dgc was initialized to Dg(i,k+1) using TSC(DGE,kbc)
	      *dpc-- = dlc = p7_FLogsum( p7_FLogsum(mln  + TSC(p7P_DM, k),  // Dl(i,k) is essentially the same recursion...
						    dlc  + TSC(p7P_DD, k)), //        
                 				    xE);                    // ... except it also includes Dk->E local exit transition, always 1.0
	      /* I(i,k) calculations */
	      *dpc-- = p7_FLogsum( mgn + TSC(p7P_IM, k),    // Ig(i,k) = Mg(i+1,k+1) * t(k,IM) * e_M(k+1, x_i+1)
				   ign + TSC(p7P_II, k));   //         + Ig(i+1,k)   * t(k,II) * e_I(k,   x_i+1)
	      *dpc-- = p7_FLogsum( mln + TSC(p7P_IM, k),    // Il(i,k) is essentially the same
				   iln + TSC(p7P_II, k));

	      /* Now pick up mgn,mln from current k on next row i+1; decrement <dpn> and <rsn> by one k unit.
	       * Thus, as we roll around and reenter the k loop, after decrementing k:
               *     mgn will be M(i+1,k+1) + e(Mk+1, xi+1), and
	       *     ign will be I(i+1,k)   + e(Ik, xi+1)
	       */
	      if (k   >= kan && k   <= kbn) { mgn = *(dpn-4) + *rsn; mln = *(dpn-5) + *rsn; dpn -= p7B_NSCELLS; rsn -= p7P_NR; } 
	      else                          { mgn = mln = -eslINFINITY; }
	      
	      /* Store mgc, mlc. */
	      *dpc-- = mgc;
	      *dpc-- = mlc;
	    }
	  kbn = kbc;
	  kan = kac;
	  dpn = last_dpc;	/* skips any remaining right overhang on row i+1 */
	}
      
      /* As we terminate a segment, we need to store specials for row ia-1.
       * Some of these (C,G,L) were already calculated as a lookahead, above.
       * Others (B,J,N,E) still need a calculation.
       * xB, xJ, xN need to be set here, to be ready for init of next segment.
       */
      *xc-- = -eslINFINITY;	/* CC only stored in a Decoding matrix. */
      *xc-- = -eslINFINITY; /* JJ only stored in a Decoding matrix. */
      *xc-- = xC;
      *xc-- = xG;
      *xc-- = xL;
      *xc-- = xB = p7_FLogsum( xL + gm->xsc[p7P_B][0],        xG + gm->xsc[p7P_B][1]);
      *xc-- = xJ = p7_FLogsum( xJ,       	              xB + gm->xsc[p7P_J][p7P_MOVE]);
      *xc-- = xN = p7_FLogsum( xN,		              xB + gm->xsc[p7P_N][p7P_MOVE]);
      *xc-- = xE = p7_FLogsum( xJ + gm->xsc[p7P_E][p7P_LOOP], xC + gm->xsc[p7P_E][p7P_MOVE]);
      last_ia = ia;
    }

  if (opt_sc) *opt_sc = xN + ((last_ia-1) ? (last_ia-1) * gm->xsc[p7P_N][p7P_LOOP] : 0.0f);
  return eslOK;
}


/*****************************************************************
 * 3. Banded posterior decoding.
 *****************************************************************/

/* Function:  p7_BandedDecoding()
 * Synopsis:  Posterior decoding, in banded DP.
 *
 * Purpose:   Given a banded Forward and Backward matrices <bmf> and
 *            <bmb>, for a comparison of query profile <gm> to a target
 *            sequence that received the overall raw lod score <totsc>,
 *            compute and store a posterior decoding matrix in
 *            <bmd> (which caller has allocated).
 *
 * Args:      gm    - query profile
 *            totsc - Forward (or Backward) raw lod score
 *            bmf   - computed Forward matrix
 *            bmb   - computed Backward matrix
 *            bmd   - RESULT: decoding matrix, space provided by caller
 *
 * Returns:   <eslOK> on success, and <bmd> contains the computed
 *            posterior decoding matrix.
 *
 * Throws:    (no abnormal error conditions)
 */
int 
p7_BandedDecoding(const P7_PROFILE *gm, float totsc, const P7_BANDMX *bmf, P7_BANDMX *bmb, P7_BANDMX *bmd)
{
  const int   *bnd_ip = bmf->bnd->imem;          /* ptr to current ia, ib segment band in bmx->bnd  */
  const int   *bnd_kp = bmf->bnd->kmem;		 /* ptr to current ka, kb row band in bmx->bnd      */
  const float *dpf    = bmf->dp;
  float       *dpb    = bmb->dp;
  float       *dpd    = bmd->dp;
  const float *xpf    = bmf->xmx;
  float       *xpb    = bmb->xmx;
  float       *xpd    = bmd->xmx;
  float        xN, xJ, xC;
  float        norm;
  int          ia,ib;
  int          last_ib;
  int          ka,kb;
  int          g,i,k,x;

  xN      = 0.0f; 
  xJ = xC = -eslINFINITY; 
  last_ib = 0;

  for (g = 0; g < bmf->bnd->nseg; g++)
    {
      ia = *bnd_ip++;
      ib = *bnd_ip++;

      /* specials on row ia-1: [ E N J B L G C JJ CC] 
       * 
       * On unbanded row, we know that N/J/C scores could only have
       * come from prv N/C/J, not from E->{CJ} for example; so xpf+xpb
       * works as a special case here. 
       */
      norm        = 0.0f;
      xpd[p7B_E]  = 0.0f;                                
      xpd[p7B_N]  = expf(xpf[p7B_N] + xpb[p7B_N] - totsc);  xN = xpf[p7B_N]; norm += xpd[p7B_N];
      xpd[p7B_J]  = expf(xpf[p7B_J] + xpb[p7B_J] - totsc);  
      xpd[p7B_B]  = expf(xpf[p7B_B] + xpb[p7B_B] - totsc);
      xpd[p7B_L]  = expf(xpf[p7B_L] + xpb[p7B_L] - totsc);
      xpd[p7B_G]  = expf(xpf[p7B_G] + xpb[p7B_G] - totsc);
      xpd[p7B_C]  = expf(xpf[p7B_C] + xpb[p7B_C] - totsc);
      xpd[p7B_JJ] = xpd[p7B_J];                             xJ = xpf[p7B_J]; norm += xpd[p7B_JJ];
      xpd[p7B_CC] = xpd[p7B_C];                             xC = xpf[p7B_C]; norm += xpd[p7B_CC];

      norm = 1.0f/norm;
      for (x = 0; x < p7B_NXCELLS; x++)	xpd[x] *= norm;

      xpf += p7B_NXCELLS;
      xpb += p7B_NXCELLS;
      xpd += p7B_NXCELLS;

      for (i = ia; i <= ib; i++)
	{
	  ka = *bnd_kp++;         /* current row's band is cells k=kac..kbc  */
	  kb = *bnd_kp++; 

	  norm = 0.0;
	  for (k = ka; k <= kb; k++)
	    { 
	      dpd[p7B_ML] = expf(dpf[p7B_ML] + dpb[p7B_ML] - totsc); norm += dpd[p7B_ML];
	      dpd[p7B_MG] = expf(dpf[p7B_MG] + dpb[p7B_MG] - totsc); norm += dpd[p7B_MG];
	      dpd[p7B_IL] = expf(dpf[p7B_IL] + dpb[p7B_IL] - totsc); norm += dpd[p7B_IL];
	      dpd[p7B_IG] = expf(dpf[p7B_IG] + dpb[p7B_IG] - totsc); norm += dpd[p7B_IG];
	      dpd[p7B_DL] = expf(dpf[p7B_DL] + dpb[p7B_DL] - totsc);                       // nonemitters don't count toward normalization
	      dpd[p7B_DG] = expf(dpf[p7B_DG] + dpb[p7B_DG] - totsc);                       // the normalization term should be ~totsc, except for local numerical error in DP matrices

	      dpd += p7B_NSCELLS; 
	      dpf += p7B_NSCELLS; 
	      dpb += p7B_NSCELLS;
	    }

	  xpd[p7B_E]  = expf(xpf[p7B_E] + xpb[p7B_E] - totsc);                                  
	  xpd[p7B_N]  = expf(xpf[p7B_N] + xpb[p7B_N] - totsc); xN = xpf[p7B_N]; norm += xpd[p7B_N];
	  xpd[p7B_J]  = expf(xpf[p7B_J] + xpb[p7B_J] - totsc);  
	  xpd[p7B_B]  = expf(xpf[p7B_B] + xpb[p7B_B] - totsc);                                  
	  xpd[p7B_L]  = expf(xpf[p7B_L] + xpb[p7B_L] - totsc);                                  
	  xpd[p7B_G]  = expf(xpf[p7B_G] + xpb[p7B_G] - totsc);                                  
	  xpd[p7B_C]  = expf(xpf[p7B_C] + xpb[p7B_C] - totsc);                                  
	  xpd[p7B_JJ] = expf(  xJ       + xpb[p7B_J] + gm->xsc[p7P_J][p7P_LOOP] - totsc); xJ = xpf[p7B_J]; norm += xpd[p7B_JJ];
	  xpd[p7B_CC] = expf(  xC       + xpb[p7B_C] + gm->xsc[p7P_C][p7P_LOOP] - totsc); xC = xpf[p7B_C]; norm += xpd[p7B_CC];

	  norm = 1.0f/norm;
	  dpd -= (kb-ka+1)*p7B_NSCELLS;
	  for (k = ka; k <= kb; k++) 
	    {
	      for (x = p7B_ML; x <= p7B_IG; x++) dpd[x] *= norm;
	      dpd += p7B_NSCELLS;
	    }
	  for (x = 0; x < p7B_NXCELLS; x++) xpd[x] *= norm;
	  
	  xpd += p7B_NXCELLS;
	  xpf += p7B_NXCELLS;
	  xpb += p7B_NXCELLS;
	}
    }
  return eslOK;
}


/*****************************************************************
 * 4. Banded alignment (MEG = maximum expected gain estimator; "gamma-centroid")
 *****************************************************************/

/* Function:  p7_BandedAlign()
 * Synopsis:  Gamma-centroid alignment, in banded DP
 *
 * Purpose:   Compute a gamma-centroid alignment of the query
 *            profile <gm> to a target sequence, given the posterior
 *            decoding matrix <bmd> for that query/target comparison,
 *            using the banded DP matrix space in <bma> for storage.
 *            Return the alignment traceback in <tr>, space preallocated
 *            and initialized by the caller and grown here as needed.
 *            Also optionally return the gain score,
 *            the total sum of pp - (1-(1+gamma)) for each state 
 *            in the state path.
 *
 *            <gamma> is the parameter of gamma-centroid alignment.
 *            Higher <gamma> increases alignment sensitivity; lower
 *            <gamma> increases specificity.  Given posterior
 *            probabilities pp(i,x) for every state x at every target
 *            position i in the DP matrix, the algorithm finds a state
 *            path consistent with the model that maximizes the sum of
 *            pp(i,x) - 1/(1+gamma).  Thus states with posterior
 *            probabilities < 1/(1+gamma) will be penalized, and those
 *            > 1/(1+gamma) are rewarded.
 *
 * Args:      gm       - query profile (we need to check consistency of alignment against nonzero transitions)
 *            gamma    - gamma-centroid parameter
 *            bmd      - banded posterior probability matrix, previously calculated by caller
 *            bma      - RESULT: filled alignment banded DP matrix; caller provides the allocated space
 *            tr       - RESULT: alignment traceback for the entire target sequence.
 *            opt_gain - optRETURN: gain score
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      [HamadaAsai11] for gamma-centroid estimation
 *            SRE:J9/137 for notes on extending Hamada/Asai gamma-centroid
 *                       from simple residue ij alignment to full state path
 *            [Kall05] for more on why the delta function on transitions is needed
 */
int 
p7_BandedAlign(const P7_PROFILE *gm, float gamma, const P7_BANDMX *bmd, P7_BANDMX *bma, P7_TRACE *tr, float *opt_gain)
{
  int         *bnd_ip = bmd->bnd->imem;          /* ptr to current ia, ib segment band in bmx->bnd       */
  int         *bnd_kp = bmd->bnd->kmem;		 /* ptr to current ka, kb row band in bmx->bnd           */
  float       *ppp    = bmd->dp;		 /* ptr into posterior probs (our 'scores'): current banded row, main states  */
  float       *ppx    = bmd->xmx;		 /* ptr into posterior probs: current row, special states */
  float       *dpc    = bma->dp;		 /* ptr into current banded row of MEA DP matrix recursion, main states  */
  float       *xc     = bma->xmx;                /* ptr into row of MEA DP matrix, special states  */
  float       *dpp    = NULL;
  float const *tsc    = gm->tsc;		 /* sets up TSC() macro, access to profile's transitions */
  float        xE, xN, xJ, xL, xG, xC;           /* tmp scores on special states. only stored when in row bands */
  float        mlp, mgp, ilp, igp, dlp, dgp;     /* M,I,D cell values from previous row i-1     */
  float        dlc,dgc;
  float        mlc,mgc;
  int          ia,ib;		/* banded segment of rows ia..ib   1<=ia<=ib<=L */
  int          kac,kbc;		/* band kac..kbc on current row i  1<=ka<=kb<=M */
  int          kap,kbp;		/* band kap..kbp on previous row i-1            */
  int          last_ib;
  int          g, i, k;
  float        gammaterm = -1.0f/(1.0f + gamma);
  
  xN      = 2.0 + 2*gammaterm;	/* S->N always have pp=1.0 at start of any trace. */
  xJ      = -eslINFINITY;
  xC      = -eslINFINITY;
  last_ib = 0;

  for (g = 0; g < bmd->bnd->nseg; g++)
    {
      ia = *bnd_ip++;
      ib = *bnd_ip++;

      /* kap,kbp initialization for i=ia:
       *  left overhang dpp advance must always eval to 0, and
       *  {m,i,d}vp initialization must always eval to -eslINFIyNITY.
       */
      kap = kbp = gm->M+1;   
      dpp = dpc;

      /* Reinitialize (and store) specials on row ia-1 [E N J B L G C]*/
      /* The unbanded residues (last_ib+1)..(ia-1) must all be emitted
       * by N, C, or J.  [(ia-1)-(last_ib+1)+1 of them, which means
       * (ia-last_ib-1) of them.] Only N->N, C->C, J->J paths are
       * possible.  Thus we know posterior probability ppx for
       * N(ia-1), J(ia-1), C(ia-1) is also the post prob for all
       * N(last_ib+1..ia-1), etc.  Hence the
       * ppx[p7B_N]*(ia-last_ib-1) term.
       *
       * (ia-last_ib-1) can be 0 for first segment with ia=1,
       * last_ib=0 where ia-1 overlaps our i=0 initialization case;
       * thus must guard against -inf * 0 = NaN; hence the test on
       * (ia-last_ib-1 == 0)
       */
      xc[p7B_E]      = -eslINFINITY;	                                                 
      xc[p7B_N] = xN = P7_DELTAT( xN + ( (ia-last_ib-1 == 0) ? 0.0f : (ppx[p7B_N]+gammaterm) * (ia-last_ib-1) ), gm->xsc[p7P_N][p7P_LOOP]);
      xc[p7B_J] = xJ = P7_DELTAT( xJ + ( (ia-last_ib-1 == 0) ? 0.0f : (ppx[p7B_J]+gammaterm) * (ia-last_ib-1) ), gm->xsc[p7P_J][p7P_LOOP]);
      xc[p7B_B]      = ppx[p7B_B] + gammaterm + ESL_MAX( P7_DELTAT(xc[p7B_N], gm->xsc[p7P_N][p7P_MOVE]), P7_DELTAT(xc[p7B_J], gm->xsc[p7P_J][p7P_MOVE]));
      xc[p7B_L] = xL = ppx[p7B_L] + gammaterm + P7_DELTAT(xc[p7B_B], gm->xsc[p7P_B][0]);                                  
      xc[p7B_G] = xG = ppx[p7B_G] + gammaterm + P7_DELTAT(xc[p7B_B], gm->xsc[p7P_B][1]);                                  
      xc[p7B_C] = xC = P7_DELTAT( xC + ( (ia-last_ib-1 == 0) ? 0.0f : (ppx[p7B_C]+gammaterm) * (ia-last_ib-1) ), gm->xsc[p7P_C][p7P_LOOP]);
      xc[p7B_JJ]     = -eslINFINITY;
      xc[p7B_CC]     = -eslINFINITY;

      xc  += p7B_NXCELLS;
      ppx += p7B_NXCELLS;
 
      for (i = ia; i <= ib; i++)
	{
	  kac       = *bnd_kp++;         /* current row's band is cells k=kac..kbc  */
	  kbc       = *bnd_kp++; 
	  dlc = dgc = -eslINFINITY;	  
	  xE        = -eslINFINITY;

	  /* dpp must advance by any left overhang of previous row; but no more than the entire row */
	  dpp += (kac-1 > kap ? ESL_MIN(kac-kap-1, kbp-kap+1) * p7B_NSCELLS : 0);

	  /* pick up scores we need from prev row -- or -eslINFINITY if that cell isn't in the band */
	  if (kac > kap && kac-1 <= kbp) { mlp = *dpp++; mgp = *dpp++; ilp = *dpp++; igp = *dpp++; dlp = *dpp++; dgp = *dpp++;       }
	  else                           { mlp =         mgp =         ilp =         igp =         dlp =         dgp = -eslINFINITY; }

	  for (k = kac; k <= kbc; k++)
	    {
	      mlc = *dpc++ = (*ppp++) + gammaterm + ESL_MAX( ESL_MAX( P7_DELTAT(mlp, TSC(p7P_MM, k-1)),
								      P7_DELTAT(ilp, TSC(p7P_IM, k-1))),
							     ESL_MAX( P7_DELTAT(dlp, TSC(p7P_DM, k-1)),
								      P7_DELTAT( xL, TSC(p7P_LM, k-1))));

	      mgc = *dpc++ = (*ppp++) + gammaterm + ESL_MAX( ESL_MAX( P7_DELTAT(mgp, TSC(p7P_MM, k-1)),
								      P7_DELTAT(igp, TSC(p7P_IM, k-1))),
							     ESL_MAX( P7_DELTAT(dgp, TSC(p7P_DM, k-1)),
								      P7_DELTAT( xG, TSC(p7P_GM, k-1))));

	      /* This "if" seems unavoidable, as we pick up vals from i-1,k on prev row */
	      if (k >= kap && k <= kbp) { mlp = *dpp++; mgp = *dpp++; ilp = *dpp++; igp = *dpp++; dlp = *dpp++; dgp = *dpp++;       }
	      else                      { mlp =         mgp =         ilp =         igp =         dlp =         dgp = -eslINFINITY; }
	      
	      /* IL/IG states */
	      *dpc++ = (*ppp++) + gammaterm + ESL_MAX( P7_DELTAT(mlp, TSC(p7P_MI, k)), P7_DELTAT(ilp, TSC(p7P_II, k))); /* IL */
	      *dpc++ = (*ppp++) + gammaterm + ESL_MAX( P7_DELTAT(mgp, TSC(p7P_MI, k)), P7_DELTAT(igp, TSC(p7P_II, k))); /* IG */

	      /* E state update with {ML,DL}->E local exits */
	      xE = ESL_MAX(xE, ESL_MAX(mlc, dlc));
	      
	      /* DL/DG states delayed storage trick */
	      *dpc++ = dlc = dlc + (*ppp++) + gammaterm;
	      *dpc++ = dgc = dgc + (*ppp++) + gammaterm;
	      dlc = ESL_MAX( P7_DELTAT(mlc, TSC(p7P_MD, k)), P7_DELTAT(dlc, TSC(p7P_DD, k))); 
	      dgc = ESL_MAX( P7_DELTAT(mgc, TSC(p7P_MD, k)), P7_DELTAT(dgc, TSC(p7P_DD, k))); 
	    }

	  xc[p7B_E]      =  ppx[p7B_E] + gammaterm + ESL_MAX( xE, P7_DELTAT(dgc, TSC(p7P_DGE, kbc))); /* dgc includes Mkb-> exit as Mkb->Dkb+1 */
	  xc[p7B_N] = xN =  ppx[p7B_N] + gammaterm +          P7_DELTAT(       xN, gm->xsc[p7P_N][p7P_LOOP]);
	  xc[p7B_J] = xJ =  ppx[p7B_J] + gammaterm + ESL_MAX( P7_DELTAT(       xJ, gm->xsc[p7P_J][p7P_LOOP]),  P7_DELTAT(xc[p7B_E], gm->xsc[p7P_E][p7P_LOOP]));
	  xc[p7B_B]      =  ppx[p7B_B] + gammaterm + ESL_MAX( P7_DELTAT(xc[p7B_J], gm->xsc[p7P_J][p7P_MOVE]),  P7_DELTAT(xc[p7B_N], gm->xsc[p7P_N][p7P_MOVE]));
	  xc[p7B_L] = xL =  ppx[p7B_L] + gammaterm +          P7_DELTAT(xc[p7B_B], gm->xsc[p7P_B][0]);
	  xc[p7B_G] = xG =  ppx[p7B_G] + gammaterm +          P7_DELTAT(xc[p7B_B], gm->xsc[p7P_B][1]);
	  xc[p7B_C] = xC =  ppx[p7B_C] + gammaterm + ESL_MAX( P7_DELTAT(       xC, gm->xsc[p7P_C][p7P_LOOP]),  P7_DELTAT(xc[p7B_E], gm->xsc[p7P_E][p7P_MOVE]));
	  xc[p7B_JJ]     =  -eslINFINITY;
	  xc[p7B_CC]     =  -eslINFINITY;

	  xc  += p7B_NXCELLS;
	  ppx += p7B_NXCELLS;
	  dpp = dpc - (kbc-kac+1)*p7B_NSCELLS;	/* skips any right overhang on the previous row: so dpp advances (if necessary) to start of curr row, then i rolls around */
	  kap = kac;
	  kbp = kbc;
	}
      last_ib = ib;
    }

  if (opt_gain) *opt_gain = xC + 1.0f + gammaterm; /* C->T, and T has pp=1.0 by construction */
  return (tr ? traceback(gm, bmd, bma, tr) : eslOK);
}

/*****************************************************************
 * 5. Banded MGE alignment traceback
 *****************************************************************/

static inline int
select_ml(const P7_PROFILE *gm, int k, const float *dpk, const float *xp)
{
  int   state[4] = { p7T_ML, p7T_IL, p7T_DL, p7T_L };
  float path[4];

  path[0] = P7_DELTAT( dpk[p7B_ML], P7P_TSC(gm, k-1, p7P_MM));
  path[1] = P7_DELTAT( dpk[p7B_IL], P7P_TSC(gm, k-1, p7P_IM));
  path[2] = P7_DELTAT( dpk[p7B_DL], P7P_TSC(gm, k-1, p7P_DM));
  path[3] = P7_DELTAT(   xp[p7B_L], P7P_TSC(gm, k-1, p7P_LM));
  return state[esl_vec_FArgMax(path, 4)];
}
static inline int
select_mg(const P7_PROFILE *gm, int k, const float *dpk, const float *xp)
{
  int   state[4] = { p7T_MG, p7T_IG, p7T_DG, p7T_G };
  float path[4];

  path[0] = P7_DELTAT( dpk[p7B_MG], P7P_TSC(gm, k-1, p7P_MM));
  path[1] = P7_DELTAT( dpk[p7B_IG], P7P_TSC(gm, k-1, p7P_IM));
  path[2] = P7_DELTAT( dpk[p7B_DG], P7P_TSC(gm, k-1, p7P_DM));
  path[3] = P7_DELTAT(   xp[p7B_G], P7P_TSC(gm, k-1, p7P_GM));
  return state[esl_vec_FArgMax(path, 4)];
}
static inline int
select_il(const P7_PROFILE *gm, int k, const float *dpk)
{
  float path[2];
  path[0] = P7_DELTAT( dpk[p7B_ML], P7P_TSC(gm, k, p7P_MI));
  path[1] = P7_DELTAT( dpk[p7B_IL], P7P_TSC(gm, k, p7P_II));
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_IL);
}
static inline int
select_ig(const P7_PROFILE *gm, int k, const float *dpk)
{
  float path[2];
  path[0] = P7_DELTAT( dpk[p7B_MG], P7P_TSC(gm, k, p7P_MI));
  path[1] = P7_DELTAT( dpk[p7B_IG], P7P_TSC(gm, k, p7P_II));
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_IG);
}
static inline int
select_dl(const P7_PROFILE *gm, int k, const float *dpk)
{
  float path[2];
  path[0] = P7_DELTAT( dpk[p7B_ML], P7P_TSC(gm, k-1, p7P_MD));
  path[1] = P7_DELTAT( dpk[p7B_DL], P7P_TSC(gm, k-1, p7P_DD));
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_DL);
}
static inline int
select_dg(const P7_PROFILE *gm, int k, const float *dpk)
{
  float path[2];
  path[0] = P7_DELTAT( dpk[p7B_MG], P7P_TSC(gm, k-1, p7P_MD));
  path[1] = P7_DELTAT( dpk[p7B_DG], P7P_TSC(gm, k-1, p7P_DD));
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_DG);
}
static inline int
select_e(const P7_PROFILE *gm, const float *dp, int ka, int kb, int *ret_k)
{
  float max  = -eslINFINITY;
  int   smax = -1;
  int   kmax = -1;
  int   k;

  for (k = ka; k < kb; k++)
    {
      if (dp[p7B_ML] >= max) { max = dp[p7B_ML]; smax = p7T_ML; kmax = k; }
      if (dp[p7B_DL] >= max) { max = dp[p7B_DL]; smax = p7T_DL; kmax = k; }
      dp+= p7B_NSCELLS;
    }
  if (dp[p7B_ML] >= max) { max = dp[p7B_ML]; smax = p7T_ML; kmax = k; }
  if (dp[p7B_DL] >= max) { max = dp[p7B_DL]; smax = p7T_DL; kmax = k; }

  /* Glocal Mk,Dk->E: Mkb + t(MD,kb->kb+1) + t(Dkb+1->E) wing retraction */
  /* remember DGE is stored off by one; TSC(gm,kb,DGE) is t(Dkb+1->E) wing retraction */
  /* to work on boundary condition kb=M, requires TSC(gm,M,DGE) = TSC(gm,M,MD) = TSC(gm,M,DD) = 0.0 */
  /* to work on boundary condition kb=M-1, requires TSC(gm,M-1,DGE) = 0.0 */
  /* and those boundary conditions are enforced: see modelconfig.c */
  if ( P7_DELTAT(dp[p7B_MG], P7P_TSC(gm, kb, p7P_MD) + P7P_TSC(gm, kb, p7P_DGE)) >= max) { max = dp[p7B_MG]; smax = p7T_MG; kmax = kb; }
  if ( P7_DELTAT(dp[p7B_DG], P7P_TSC(gm, kb, p7P_DD) + P7P_TSC(gm, kb, p7P_DGE)) >= max) { max = dp[p7B_DG]; smax = p7T_DG; kmax = kb; }

  *ret_k = kmax;
  return   smax;
}
static inline int
select_j(const P7_PROFILE *gm, const float *xc)
{
  float path[2];
  path[0] = P7_DELTAT(*(xc-p7B_NXCELLS+p7B_J), gm->xsc[p7P_J][p7P_LOOP]); /* i.e. xp[p7B_J] on prv row i-1. */
  path[1] = P7_DELTAT( xc[p7B_E],               gm->xsc[p7P_E][p7P_LOOP]);
  return ( (path[0] > path[1]) ? p7T_J : p7T_E);
}
static inline int
select_c(const P7_PROFILE *gm, const float *xc)
{
  float path[2];
  path[0] = P7_DELTAT(*(xc-p7B_NXCELLS+p7B_C), gm->xsc[p7P_C][p7P_LOOP]); /* i.e. xp[p7B_C] on prv row i-1. */
  path[1] = P7_DELTAT(xc[p7B_E],               gm->xsc[p7P_E][p7P_MOVE]);
  return ( (path[0] > path[1]) ? p7T_C : p7T_E);
}
static inline int
select_b(const P7_PROFILE *gm, const float *xc)
{
  float path[2];
  path[0] = P7_DELTAT(xc[p7B_J], gm->xsc[p7P_J][p7P_MOVE]);
  path[1] = P7_DELTAT(xc[p7B_N], gm->xsc[p7P_N][p7P_MOVE]);
  return ( (path[0] > path[1]) ? p7T_J : p7T_N);
}

static int
traceback(const P7_PROFILE *gm, const P7_BANDMX *bmd, const P7_BANDMX *bma, P7_TRACE *tr)
{
  int   *bnd_ip;		/* ptr into <bnd>'s array of segments ia..ib coords                   */
  int   *bnd_kp;		/* ptr into <bnd>'s array of <ka..kb> band extent on each banded row  */
  const float *dp;		/* points to next banded row we will trace to (i or i-1, depending)   */
  const float *xc;		/* points to specials on current row i                                */
  const float *ppp;		/* points to main states, row i, in banded post prob matrix           */
  const float *ppx;		/* points to specials, row i, in posterior probability matrix         */
  int ia,ib;
  int ka,kb;
  int g,i,k,k2;
  int scur, snxt;
  float ppv;
  int   last_ia;
  int   decrement_i;
  int   do_NCJ_emit;
  int   decrement_dp;
  int   status;
  

  /* In Backwards, various pointers get set to the _end_ of arrays, which we traverse backwards            */
  bnd_ip  = bma->bnd->imem + (           2 * bma->bnd->nseg - 1);             /* set to _last_ array elem            */
  bnd_kp  = bma->bnd->kmem + (p7_GBANDS_NK * bma->bnd->nrow - 1);             /* ditto                               */

  dp  = bma->dp  + (p7B_NSCELLS * bma->bnd->ncell);                          /* initialize to just off edge (can't dereference yet!) */
  ppp = bmd->dp  + (p7B_NSCELLS * bma->bnd->ncell);                          /* initialize to just off edge (can't dereference yet!) */
  xc  = bma->xmx + (p7B_NXCELLS * (bma->bnd->nrow+bma->bnd->nseg));          /* ditto: just off edge in xmx;  no dereference yet    */
  ppx = bmd->xmx + (p7B_NXCELLS * (bma->bnd->nrow+bma->bnd->nseg));          /* ditto for ppx */
  
  last_ia = bma->bnd->L+1;
  i       = bma->bnd->L;
  k       = 0;
  if ((status = p7_trace_AppendWithPP(tr, p7T_T, 0, 0, 0.0)) != eslOK) return status;
  if ((status = p7_trace_AppendWithPP(tr, p7T_C, 0, 0, 0.0)) != eslOK) return status; 
  scur = p7T_C;
  ppv  = 1.0f;

  for (g = bmd->bnd->nseg-1; g >= 0; g--)
    {
      ib = *bnd_ip--;		/* note ib picked up before ia, because we're accessing in reverse */
      ia = *bnd_ip--;

      /* intersegment residues ib+1..last_ia-1 were emitted from <sprv> (NCJ) with posterior prob <ppv> that we rolled around from prev segment */
      for (i = last_ia-1; i > ib; i--)
	if ((status = p7_trace_AppendWithPP(tr, scur, 0, i, ppv)) != eslOK) return status;
      /* that loop left us at i=ib, scur={NCJ}, k=undefined */

      /* As we enter the segment: dp,ppp pointing to last ia; xc,ppx to last_ia-1 */
      kb   = *bnd_kp--;
      ka   = *bnd_kp--;	
      dp  -= (kb-ka+1)*p7B_NSCELLS;
      ppp -= (kb-ka+1)*p7B_NSCELLS;
      xc  -= p7B_NXCELLS;	      /* xc reinitialized to row ib                 */
      ppx -= p7B_NXCELLS;	      /* ditto ppx                                  */

      /* When we trace a segment back to {NCJ} on ia-1, we're done with that segment. */
      while (i >= ia || (scur != p7T_N && scur != p7T_J && scur != p7T_C))
	{
	  decrement_dp = FALSE;
	  decrement_i  = FALSE;
	  do_NCJ_emit  = FALSE;

	  /* dp always points to the next main banded row we'll be tracing to;
	   * so dp moves up a row when scur is an M or I or on an NN/CC/JJ emit.
	   * We set <decrement_dp> for the M/I case.
	   *
	   * xc (and ppx), in contrast, always point to current row i, and they
	   * move whenever i decrements: that's whenever *sprv* is M or I,
	   * or on an NN/CC/JJ emit. We set <decrement_i> for the M/I case.
	   *
	   * When we emit from an NN/CC/JJ on transition, these decrements have 
	   * to be deferred until after we've picked up the posterior probability
	   * of the *downstream* N/C/J. We set do_NCJ_emit for this case, which
	   * implies both decrement_dp and decrement_i, deferred.
	   */
	  switch (scur) {
	  case p7T_ML:  snxt = ( (i-1>=ia && k-1>=ka && k-1<=kb) ? select_ml(gm, k, dp+(k-ka-1)*p7B_NSCELLS, xc-p7B_NXCELLS) : p7T_L);     k--; decrement_i  = TRUE; break;
	  case p7T_MG:  snxt = ( (i-1>=ia && k-1>=ka && k-1<=kb) ? select_mg(gm, k, dp+(k-ka-1)*p7B_NSCELLS, xc-p7B_NXCELLS) : p7T_G);     k--; decrement_i  = TRUE; break;
	  case p7T_IL:  snxt = ( (i-1>=ia &&   k>=ka &&   k<=kb) ? select_il(gm, k, dp+(k-ka)  *p7B_NSCELLS)                 : p7T_BOGUS);      decrement_i  = TRUE; break;
	  case p7T_IG:  snxt = ( (i-1>=ia &&   k>=ka &&   k<=kb) ? select_ig(gm, k, dp+(k-ka)  *p7B_NSCELLS)                 : p7T_BOGUS);      decrement_i  = TRUE; break;
	  case p7T_DL:  snxt = (                       (k-1>=ka) ? select_dl(gm, k, dp+(k-ka-1)*p7B_NSCELLS)                 : p7T_BOGUS); k--;                      break;                         
	  case p7T_DG:  snxt = (                       (k-1>=ka) ? select_dg(gm, k, dp+(k-ka-1)*p7B_NSCELLS)                 : p7T_BOGUS); k--;                      break;

 	  case p7T_N:   snxt = p7T_N;             do_NCJ_emit = TRUE;                           break;
	  case p7T_J:   snxt = select_j(gm, xc);  do_NCJ_emit = (snxt==p7T_J ? TRUE : FALSE);   break;
	  case p7T_C:   snxt = select_c(gm, xc);  do_NCJ_emit = (snxt==p7T_C ? TRUE : FALSE);   break;

	  case p7T_E:   snxt = select_e (gm, dp, ka, kb, &k); break;
	  case p7T_B:   snxt = select_b (gm, xc);             break;                             
	  case p7T_L:   snxt = p7T_B;                         break;
	  case p7T_G:   snxt = p7T_B;                         break;
	  default:      ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in traceback");
	  }
	  if (snxt == p7T_ML || snxt == p7T_MG || snxt == p7T_IL || snxt == p7T_IG) decrement_dp = TRUE;

	  /* Pick up posterior probability annotation before moving
	   * <dp> or <xc> where snxt has connected to. We annotate
	   * only residue emission ppv on alignments, not the
	   * fully decoded state ppv; moreover, for M/I, we marginalize over
	   * glocal/local.
	   */
	  switch (snxt) {
	  case p7T_ML: case p7T_MG: ppv = ppp[(k-ka)*p7B_NSCELLS+p7B_ML] + ppp[(k-ka)*p7B_NSCELLS+p7B_MG]; break;  // marginalized over glocal/local
	  case p7T_IL: case p7T_IG: ppv = ppp[(k-ka)*p7B_NSCELLS+p7B_IL] + ppp[(k-ka)*p7B_NSCELLS+p7B_IG]; break;
	  case p7T_DL: case p7T_DG: ppv = 0.0f;                               break;
	  case p7T_E:               ppv = 0.0f;                               break;
	  case p7T_N:               ppv = (do_NCJ_emit ? ppx[p7B_N]  : 0.0f); break;
	  case p7T_J:               ppv = (do_NCJ_emit ? ppx[p7B_JJ] : 0.0f); break; /* note, pick up the *residue emission* pp for annotation, not the occupancy pp */
	  case p7T_B:               ppv = 0.0f;                               break;
	  case p7T_L:               ppv = 0.0f;                               break;
	  case p7T_G:               ppv = 0.0f;                               break;
	  case p7T_C:               ppv = (do_NCJ_emit ? ppx[p7B_CC] : 0.0f); break; /* ditto */
	  default:     ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in traceback");
	  }

	  /* Now we can decrement dp/ppp, if we're supposed to. Do this before decrementing i */
	  if ((i>ia+1 && decrement_dp) || (i>ia && do_NCJ_emit))
	    {
	      kb   = *bnd_kp--;
	      ka   = *bnd_kp--;	
	      dp  -= (kb-ka+1)*p7B_NSCELLS;
	      ppp -= (kb-ka+1)*p7B_NSCELLS;
	    }

	  /* Now we decrement i - before we attach snxt - if scur on row i is on an M/I state; snxt is on row i-1 */
	  if (decrement_i)
	    {
	      xc  -= p7B_NXCELLS; /* xc might reach ia-1 here, boundary condition for the segment; we need to trace to an NCJ on that boundary before calling segment done */
	      ppx -= p7B_NXCELLS;
	      i--;
	    }

	  /* Glocal B->G->Mk left wing retraction entry: unfold it */
	  if (snxt == p7T_G) {
	    while (k>1) {
	      if ( (status = p7_trace_AppendWithPP(tr, p7T_DG, k-1, i, 0.0)) != eslOK) return status;
	      k--;
	    }
	  }	  
	  /* Glocal Mk->E right wing retraction: off band edge kb, Mkb->Dkb+1->E or Dkb->Dkb+1->E */
	  if (scur == p7T_E && (snxt == p7T_MG || snxt == p7T_DG))
	    {
	      for (k2 = gm->M; k2 > k; k2--)
		if ( (status = p7_trace_AppendWithPP(tr, p7T_DG, k2, i, 0.0)) != eslOK) return status;
	    }
	  /* at last, append the traceback state itself */
	  if ( (status = p7_trace_AppendWithPP(tr, snxt, k, i, ppv)) != eslOK) return status; /* remember, on a NN/CC/JJ, we're doing a deferred emit of the *previous* i */

	  /* Deferred i decrement for NN/CC/JJ */
	  if (do_NCJ_emit)
	    {
	      xc  -= p7B_NXCELLS; /* xc might reach ia-1 here, boundary condition for the segment; we need to trace to an NCJ on that boundary before calling segment done */
	      ppx -= p7B_NXCELLS;
	      i--;
	    }

	  /* and around we go. */
	  scur = snxt;
	} /* end loop over segment ia..ib.*/
      /* as we finish the segment: i  = ia-1, xc is on specials for ia-1, dp is still on ia; snxt = NCJ, which will emit the entire intersegment with the post prob that's at this boundary row ia-1 */
      last_ia = ia;
      switch (snxt) {
      case p7T_N: ppv = ppx[p7B_N];  break;
      case p7T_C: ppv = ppx[p7B_CC]; break;
      case p7T_J: ppv = ppx[p7B_JJ]; break;
      default:    ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in traceback");
      }
    }
  /* as we finish the last segment: i=ia-1 (which could be 0), xc is on ia-1, sprv=N 
   * residues 1..ia-1 are emitted on N with post prob 1.0
   */
  for (; i >= 1; i--)
    if ( (status = p7_trace_AppendWithPP(tr, p7T_N, 0, i, 1.0f)) != eslOK) return status; 
  if ( (status = p7_trace_AppendWithPP(tr, p7T_S, 0, i, 0.0f))   != eslOK) return status; 
  
  tr->M = bma->bnd->M;
  tr->L = bma->bnd->L;
  return p7_trace_Reverse(tr);
}




/*****************************************************************
 * x. Benchmark driver
 *****************************************************************/ 
#ifdef p7BANDED_FWDBACK_BENCHMARK
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "p7_gbands.h"
#include "p7_bandmx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,   "2000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Backward",                        0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Forward",                         0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for banded dual local/glocal Forward/Backward implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GBANDS      *bnd     = NULL;
  P7_BANDMX      *bmf     = NULL;
  P7_BANDMX      *bmb     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, bench_time, Mcs;

  /* Initialize log-sum calculator */
  impl_Init();
  p7_FLogsumInit();

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);

  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, L);

  bnd = p7_gbands_Create(hmm->M, L);
  p7_gbands_SetFull(bnd);
  
  bmf = p7_bandmx_Create(bnd);
  bmb = p7_bandmx_Create(bnd);

  /* Baseline time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_bandmx_Reinit(bmf, bnd);
      p7_bandmx_Reinit(bmb, bnd);

      if (! esl_opt_GetBoolean(go, "-B"))  p7_BandedForward  (dsq, L, gm, bmf, &sc);
      if (! esl_opt_GetBoolean(go, "-F"))  p7_BandedBackward (dsq, L, gm, bmb, NULL);

      p7_bandmx_Reuse(bmf);
      p7_bandmx_Reuse(bmb);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_bandmx_Destroy(bmb);
  p7_bandmx_Destroy(bmf);
  p7_gbands_Destroy(bnd);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7BANDED_FWDBACK_BENCHMARK*/
/*------------------ end, benchmark -----------------------------*/


/*****************************************************************
 * x. Unit tests
 *****************************************************************/
#ifdef p7BANDED_FWDBACK_TESTDRIVE
#include "esl_random.h"
#include "esl_sq.h"

static void
utest_single_paths(ESL_GETOPTS *go)
{
  char            msg[] = "single path unit test failed";
  int             M     = esl_opt_GetInteger(go, "-M");
  ESL_ALPHABET   *abc   = esl_alphabet_Create(eslCOINS);
  ESL_RANDOMNESS *rng   = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  P7_HMM         *hmm   = NULL;
  P7_BG          *bg    = p7_bg_Create(abc);
  P7_PROFILE     *gm    = p7_profile_Create(M, abc);
  ESL_SQ         *sq    = esl_sq_CreateDigital(abc);
  P7_TRACE       *tr    = p7_trace_Create();
  P7_GBANDS      *bnd   = NULL;
  P7_BANDMX      *bmf   = p7_bandmx_Create(NULL);	  
  P7_BANDMX      *bmb   = p7_bandmx_Create(NULL);	  
  float           tsc,fsc,bsc;

  /* Create a profile that has only a single possible path thru it */
  p7_hmm_SampleSinglePathed(rng, M, abc, &hmm);
  p7_profile_ConfigUniglocal(gm, hmm, bg, 0); 

  /* Sample that sequence */
  p7_ProfileEmit(rng, hmm, gm, bg, sq, tr);

  /* Initialize bands, banded matrix */
  bnd = p7_gbands_Create(M, sq->n);
  p7_gbands_SetFull(bnd);

  p7_bandmx_Reinit(bmf, bnd);
  p7_bandmx_Reinit(bmb, bnd);

  /* Since all the probability mass is in a single path,
   * the trace score is equal to the forward and backward scores.
   */
  p7_trace_Score(tr, sq->dsq, gm, &tsc);
  p7_BandedForward (sq->dsq, sq->n, gm, bmf, &fsc);
  p7_BandedBackward(sq->dsq, sq->n, gm, bmb, &bsc);

  //printf("n   = %" PRId64 "\n", sq->n);
  //printf("tsc = %.2f nats\n",   tsc);
  //printf("fsc = %.2f nats\n",   fsc);
  //printf("bsc = %.2f nats\n",   bsc);


  if ( esl_FCompareAbs(fsc, tsc, 0.0001) != eslOK) esl_fatal(msg);
  if ( esl_FCompareAbs(bsc, tsc, 0.0001) != eslOK) esl_fatal(msg);

  p7_bandmx_Destroy(bmb);
  p7_bandmx_Destroy(bmf);
  p7_gbands_Destroy(bnd);
  p7_trace_Destroy(tr);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_randomness_Destroy(rng);
  esl_alphabet_Destroy(abc);
}

#endif /*p7BANDED_FWDBACK_TESTDRIVE*/

/*---------------- end, unit tests ------------------------------*/

/*****************************************************************
 * x. Test driver
 *****************************************************************/
#ifdef p7BANDED_FWDBACK_TESTDRIVE

#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",        eslARG_INT,     "10", NULL, NULL,  NULL,  NULL, NULL, "set length of sampled models <n>",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for banded Forward/Backward dual-mode implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);

  p7_FLogsumInit();

  utest_single_paths(go);

  esl_getopts_Destroy(go);
  return 0;
}


#endif /*p7BANDED_FWDBACK_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/

/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7BANDED_FWDBACK_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_regexp.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_gbands.h"
#include "p7_bandmx.h"

static int parse_coord_string(const char *cstring, int *ret_start, int *ret_end);

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",              0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "1 line per sequence, tabular summary output",       0 },
  { "-i",        eslARG_STRING, FALSE, NULL, NULL,   NULL,  NULL, NULL, "when dumping, restrict dump to rows <i1>..<i2>",    0 },
  { "-k",        eslARG_STRING, FALSE, NULL, NULL,   NULL,  NULL, NULL, "when dumping, restrict dump to columns <k1>..<k2>", 0 },
  { "-G",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump DP bands for examination",                     0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Forward DP matrix for examination",            0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Backward DP matrix for examination",           0 },
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump posterior Decoding matrix for examination",    0 },
  { "-A",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump alignment DP matrix for examination",          0 },
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump alignment Traceback for examination",          0 },
#ifdef p7_DEBUGGING
  { "--fF",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Forward filter matrix for examination",      0 },
  { "--fB",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Backward filter matrix for examination",     0 },
  { "--fD",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Decoding filter matrix for examination",     0 },
  { "--fT",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump superannotated alignment Traceback",         0 },
#endif
  { "--fs",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "config model for multihit local mode only",       0 },
  { "--sw",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "config model for unihit local mode only",         0 },
  { "--ls",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "config model for multihit glocal mode only",      0 },
  { "--s",       eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "config model for unihit glocal mode only",        0 },
  { "--full",    eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "wide open bands",                                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of Forward/Backward, banded dual local/glocal implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_FILTERMX    *ox      = NULL;
  P7_GBANDS      *bnd     = NULL;
  P7_BANDMX      *bmf     = NULL;
  P7_BANDMX      *bmb     = NULL;
  P7_BANDMX      *bmd     = NULL;
  P7_BANDMX      *bma     = NULL;
  P7_TRACE       *tr      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fsc, bsc;
  float           nullsc;
  float           gain;
  int             istart, iend, kstart, kend;
  char           *cstring;
  int             status;

  /* Initialize log-sum calculator */
  impl_Init();
  p7_FLogsumInit();

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Determine coords of dump windows */
  istart = iend = 0;
  kstart = kend = 0;
  if ( esl_opt_IsOn(go, "-i")) {
    cstring = esl_opt_GetString(go, "-i");
    parse_coord_string(cstring, &istart, &iend);
  }
  if ( esl_opt_IsOn(go, "-k")) {
    cstring = esl_opt_GetString(go, "-k");
    parse_coord_string(cstring, &kstart, &kend);
  }

  /* Open sequence database */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);

  /* Now reconfig the models however we were asked to */
  if      (esl_opt_GetBoolean(go, "--fs"))  p7_profile_ConfigLocal    (gm, hmm, bg, 400);
  else if (esl_opt_GetBoolean(go, "--sw"))  p7_profile_ConfigUnilocal (gm, hmm, bg, 400);
  else if (esl_opt_GetBoolean(go, "--ls"))  p7_profile_ConfigGlocal   (gm, hmm, bg, 400);
  else if (esl_opt_GetBoolean(go, "--s"))   p7_profile_ConfigUniglocal(gm, hmm, bg, 400);
  else                                      p7_profile_Config         (gm, hmm, bg);

  /* Rearrange into an optimized profile */
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  /* Allocate bands, matrices */
  bnd    = p7_gbands_Create(gm->M, 400); /* L=400 is a dummy; Reinit() will fix */
  bmf    = p7_bandmx_Create(NULL);	
  bmb    = p7_bandmx_Create(NULL);	
  bmd    = p7_bandmx_Create(NULL);	
  bma    = p7_bandmx_Create(NULL);	
  ox     = p7_filtermx_Create(om->M, 400, ESL_MBYTES(32));
  tr     = p7_trace_CreateWithPP();

#ifdef p7_DEBUGGING
  /* Under debugging mode only, we can also dump internals of the ForwardFilter/BackwardFilter vector calculations */
  if (esl_opt_GetBoolean(go, "--fF")) ox->fwd = p7_refmx_Create(gm->M, 100);
  if (esl_opt_GetBoolean(go, "--fB")) ox->bck = p7_refmx_Create(gm->M, 100);
  if (esl_opt_GetBoolean(go, "--fD") || 
      esl_opt_GetBoolean(go, "--fT")) ox->pp  = p7_refmx_Create(gm->M, 100);
#endif

  if (esl_opt_GetBoolean(go, "-1")) {
    printf("%-30s   %-10s %-10s   %-10s %-10s\n", "# seq name",      "fwd (raw)",   "bck (raw) ",  "fwd (bits)",  "bck (bits)");
    printf("%-30s   %10s %10s   %10s %10s\n",     "#--------------", "----------",  "----------",  "----------",  "----------");
  }

  while ( (status = esl_sqio_Read(sqfp, sq)) != eslEOF)
    {
      if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
      else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

      /* Set the profile and null model's target length models */
      p7_bg_SetLength           (bg, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_oprofile_ReconfigLength(om, sq->n);

      /* Determine bands */
      p7_gbands_Reinit(bnd, gm->M, sq->n);
      p7_filtermx_GrowTo(ox, om->M, sq->n);
      p7_ForwardFilter (sq->dsq, sq->n, om, ox, &fsc);
      p7_BackwardFilter(sq->dsq, sq->n, om, ox, bnd);

      /* --full overrides the bands; but we still need to run the fwdfilter, for refmx_Dumps and trace_DumpSuper below */
      if (esl_opt_GetBoolean(go, "--full")) p7_gbands_SetFull(bnd);

#ifdef p7_DEBUGGING
      if (esl_opt_GetBoolean(go, "--fF")) p7_refmx_Dump(stdout, ox->fwd);
      if (esl_opt_GetBoolean(go, "--fB")) p7_refmx_Dump(stdout, ox->bck);
      if (esl_opt_GetBoolean(go, "--fD")) p7_refmx_Dump(stdout, ox->pp);
#endif
      if (esl_opt_GetBoolean(go, "-G")) p7_gbands_Dump(stdout, bnd);

      /* Resize banded DP matrix if necessary */
      p7_bandmx_Reinit(bmf, bnd);
      p7_bandmx_Reinit(bmb, bnd);
      p7_bandmx_Reinit(bmd, bnd);
      p7_bandmx_Reinit(bma, bnd);

      /* Run banded Forward, Backward */
      p7_BandedForward (sq->dsq, sq->n, gm, bmf, &fsc);
      p7_BandedBackward(sq->dsq, sq->n, gm, bmb, &bsc);
      p7_BandedDecoding(gm, fsc, bmf, bmb, bmd);
      p7_BandedAlign   (gm, /*gamma=*/1.0, bmd, bma, tr, &gain);

      if (esl_opt_GetBoolean(go, "-F"))   p7_bandmx_DumpWindow(stdout, bmf, (istart ? istart: 0), (iend ? iend: sq->n), (kstart? kstart : 0), (kend? kend:gm->M));
      if (esl_opt_GetBoolean(go, "-B"))   p7_bandmx_DumpWindow(stdout, bmb, (istart ? istart: 0), (iend ? iend: sq->n), (kstart? kstart : 0), (kend? kend:gm->M));
      if (esl_opt_GetBoolean(go, "-D"))   p7_bandmx_DumpWindow(stdout, bmd, (istart ? istart: 0), (iend ? iend: sq->n), (kstart? kstart : 0), (kend? kend:gm->M));
      if (esl_opt_GetBoolean(go, "-A"))   p7_bandmx_DumpWindow(stdout, bma, (istart ? istart: 0), (iend ? iend: sq->n), (kstart? kstart : 0), (kend? kend:gm->M));
      if (esl_opt_GetBoolean(go, "-T"))   p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);
#ifdef p7_DEBUGGING
      if (esl_opt_GetBoolean(go, "--fT")) p7_trace_DumpSuper    (stdout, tr, gm, sq->dsq, /*gamma=*/1.0f, ox->pp, bmd);
#endif

      /* Those scores are partial log-odds likelihoods in nats.
       * Subtract off the rest of the null model, convert to bits.
       */
      p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

      if (esl_opt_GetBoolean(go, "-1")) 
	{
	  printf("%-30s   %10.4f %10.4f   %10.4f %10.4f  %10.2f  %10.3fM  %10.3fM\n", 
		 sq->name, 
		 fsc, bsc, 
		 (fsc - nullsc) / eslCONST_LOG2, (bsc - nullsc) / eslCONST_LOG2,
		 gain,
		 (double) (p7_filtermx_Sizeof(ox) / 1000000),
		 (double) (p7_bandmx_Sizeof(bmf)  / 1000000));
	}
      else
	{
	  printf("target sequence:      %s\n",         sq->name);
	  printf("fwd raw score:        %.4f nats\n",  fsc);
	  printf("bck raw score:        %.4f nats\n",  bsc);
	  printf("null score:           %.2f nats\n",  nullsc);
	  printf("per-seq score:        %.2f bits\n",  (fsc - nullsc) / eslCONST_LOG2);
	  printf("gain score:           %.2f \n",      gain);
	  printf("RAM usage, filter:    %.3fM\n",      (double) (p7_filtermx_Sizeof(ox) / 1000000));
	  printf("RAM usage, bands:     %.3fM (x4)\n", (double) (p7_bandmx_Sizeof(bmf)  / 1000000));
	}

      p7_trace_Reuse(tr);
      p7_filtermx_Reuse(ox);
      p7_bandmx_Reuse(bmf);
      p7_bandmx_Reuse(bmb);
      p7_gbands_Reuse(bnd);
      esl_sq_Reuse(sq);
    }

  /* Cleanup */
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_gbands_Destroy(bnd);
  p7_bandmx_Destroy(bma);
  p7_bandmx_Destroy(bmd);
  p7_bandmx_Destroy(bmb);
  p7_bandmx_Destroy(bmf);
  p7_filtermx_Destroy(ox);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}

static int
parse_coord_string(const char *cstring, int *ret_start, int *ret_end)
{
  ESL_REGEXP *re = esl_regexp_Create();
  char        tok1[32];
  char        tok2[32];

  if (esl_regexp_Match(re, "^(\\d+)\\D+(\\d*)$", cstring) != eslOK) esl_fatal("-c takes arg of subseq coords <from>..<to>; %s not recognized", cstring);
  if (esl_regexp_SubmatchCopy(re, 1, tok1, 32)            != eslOK) esl_fatal("Failed to find <from> coord in %s", cstring);
  if (esl_regexp_SubmatchCopy(re, 2, tok2, 32)            != eslOK) esl_fatal("Failed to find <to> coord in %s",   cstring);
  
  *ret_start = atol(tok1);
  *ret_end   = (tok2[0] == '\0') ? 0 : atol(tok2);
  
  esl_regexp_Destroy(re);
  return eslOK;
}


#endif /*p7BANDED_FWDBACK_EXAMPLE*/

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
