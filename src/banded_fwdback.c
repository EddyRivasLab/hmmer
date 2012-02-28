/* Banded Forward/Backward, posterior decoding, and optimal accuracy alignment,
 * with a dual-mode glocal/local model
 * 
 * Contents:
 *    x. Banded Forward
 *    x. Banded Backward
 *    x. Example
 *    x. Copyright and license.
 */

#include "easel.h"

#include "hmmer.h"
#include "p7_bandmx.h"


/*****************************************************************
 * 1. Banded Forward
 *****************************************************************/

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
      xN  = xN + ( ia == last_ib+1 ? 0.0f : (ia - last_ib - 1) * gm->xsc[p7P_N][p7P_LOOP]); /* watch out for 0*-inf special case */
      xJ  = xJ + ( ia == last_ib+1 ? 0.0f : (ia - last_ib - 1) * gm->xsc[p7P_J][p7P_LOOP]);
      xB  = p7_FLogsum( xN + gm->xsc[p7P_N][p7P_MOVE], xJ + gm->xsc[p7P_J][p7P_MOVE]);
      xL  = xB + gm->xsc[p7P_B][0]; /* B->L */
      xG  = xB + gm->xsc[p7P_B][1]; /* B->G */
      xC  = xC + ( ia == last_ib+1 ? 0.0f : (ia - last_ib - 1) * gm->xsc[p7P_C][p7P_LOOP]);

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
  if (opt_sc != NULL) *opt_sc = xC + ( L-last_ib ? (L-last_ib) *  gm->xsc[p7P_C][p7P_LOOP] : 0.0f) + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}

/*****************************************************************
 * 2. Banded Backward
 *****************************************************************/

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

  /* Initialization condition for xC/xN/xJ follows from a convention
   * that as we start row i, they already include all contribution
   * from row i+1... as we start, that's L+1.
   */
  xC      = gm->xsc[p7P_C][p7P_MOVE];
  xJ      = -eslINFINITY;
  xN      = -eslINFINITY;
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
      
      /* As we terminate a segment, we need to leave xN/xC/xJ calculated for row ia-1 */
      xB = p7_FLogsum( xL + gm->xsc[p7P_B][0], xG + gm->xsc[p7P_B][1]);
      xJ = p7_FLogsum( xJ,       	       xB + gm->xsc[p7P_J][p7P_MOVE]);
      xN = p7_FLogsum( xN,		       xB + gm->xsc[p7P_N][p7P_MOVE]);
      last_ia = ia;
    }

  if (opt_sc) *opt_sc = xN + ((last_ia-1) ? (last_ia-1) * gm->xsc[p7P_N][p7P_LOOP] : 0.0f);
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
  P7_BANDMX      *bmx     = NULL;
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

  bnd = p7_gbands_Create();
  p7_gbands_SetFull(bnd, gm->M, L);
  
  bmx = p7_bandmx_Create(bnd);

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

      p7_bandmx_GrowTo(bmx, bnd);

      if (! esl_opt_GetBoolean(go, "-B"))  p7_BandedForward  (dsq, L, gm, bmx, &sc);
      if (! esl_opt_GetBoolean(go, "-F"))  p7_BandedBackward (dsq, L, gm, bmx, NULL);

      p7_bandmx_Reuse(bmx);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_bandmx_Destroy(bmx);
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
  P7_GBANDS      *bnd   = p7_gbands_Create();
  P7_BANDMX      *bmx   = p7_bandmx_Create(NULL);	  
  float           tsc,fsc,bsc;

  /* Create a profile that has only a single possible path thru it */
  p7_hmm_SampleSinglePathed(rng, M, abc, &hmm);
  p7_profile_ConfigUniglocal(gm, hmm, bg, 0); 

  /* Sample that sequence */
  p7_ProfileEmit(rng, hmm, gm, bg, sq, tr);

  /* Initialize bands, banded matrix */
  p7_gbands_SetFull(bnd, M, sq->n);
  p7_bandmx_GrowTo(bmx, bnd);

  /* Since all the probability mass is in a single path,
   * the trace score is equal to the forward and backward scores.
   */
  p7_trace_Score(tr, sq->dsq, gm, &tsc);
  p7_BandedForward (sq->dsq, sq->n, gm, bmx, &fsc);
  p7_BandedBackward(sq->dsq, sq->n, gm, bmx, &bsc);

  printf("n   = %" PRId64 "\n", sq->n);
  printf("tsc = %.2f nats\n",   tsc);
  printf("fsc = %.2f nats\n",   fsc);
  printf("bsc = %.2f nats\n",   bsc);


  if ( esl_FCompareAbs(fsc, tsc, 0.0001) != eslOK) esl_fatal(msg);
  if ( esl_FCompareAbs(bsc, tsc, 0.0001) != eslOK) esl_fatal(msg);

  p7_bandmx_Destroy(bmx);
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
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_gbands.h"
#include "p7_bandmx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",            0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "1 line per sequence, tabular summary output",     0 },
  { "-A",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump OA alignment DP matrix for examination",     0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Backward DP matrix for examination",         0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Forward DP matrix for examination",          0 },
  { "-G",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump DP bands for examination",                   0 },
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump posterior Decoding matrix for examination",  0 },
#ifdef p7_DEBUGGING
  { "--fF",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Forward filter matrix for examination",      0 },
  { "--fB",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Backward filter matrix for examination",     0 },
  { "--fD",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Decoding filter matrix for examination",     0 },
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
  P7_BANDMX      *bmx     = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fsc, bsc;
  float           nullsc;
  int             status;

  /* Initialize log-sum calculator */
  impl_Init();
  p7_FLogsumInit();

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
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
  bnd = p7_gbands_Create();
  bmx = p7_bandmx_Create(NULL);	
  ox  = p7_filtermx_Create(om->M, 400, ESL_MBYTES(32));

#ifdef p7_DEBUGGING
  /* Under debugging mode only, we can also dump internals of the ForwardFilter/BackwardFilter vector calculations */
  if (esl_opt_GetBoolean(go, "--fF")) ox->fwd = p7_gmx_Create(gm->M, 100);
  if (esl_opt_GetBoolean(go, "--fB")) ox->bck = p7_gmx_Create(gm->M, 100);
  if (esl_opt_GetBoolean(go, "--fD")) ox->pp  = p7_gmx_Create(gm->M, 100);
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
      if (esl_opt_GetBoolean(go, "--full"))
	p7_gbands_SetFull(bnd, gm->M, sq->n);
      else
	{
	  p7_filtermx_GrowTo(ox, om->M, sq->n);
	  p7_ForwardFilter (sq->dsq, sq->n, om, ox, &fsc);
	  p7_BackwardFilter(sq->dsq, sq->n, om, ox, bnd);
	}

#ifdef p7_DEBUGGING
      if (esl_opt_GetBoolean(go, "--fF")) p7_gmx_Dump(stdout, ox->fwd, p7_DEFAULT);
      if (esl_opt_GetBoolean(go, "--fB")) p7_gmx_Dump(stdout, ox->bck, p7_DEFAULT);
      if (esl_opt_GetBoolean(go, "--fD")) p7_gmx_Dump(stdout, ox->pp,  p7_DEFAULT);
#endif

      if (esl_opt_GetBoolean(go, "-G")) p7_gbands_Dump(stdout, bnd);

      /* Resize banded DP matrix if necessary */
      p7_bandmx_GrowTo(bmx, bnd);

      /* Run banded Forward, Backward */
      p7_BandedForward (sq->dsq, sq->n, gm, bmx, &fsc);
      p7_BandedBackward(sq->dsq, sq->n, gm, bmx, &bsc);

      if (esl_opt_GetBoolean(go, "-F")) p7_bandmx_Dump(stdout, bmx, p7B_FORWARD);
      if (esl_opt_GetBoolean(go, "-B")) p7_bandmx_Dump(stdout, bmx, p7B_BACKWARD);

      /* Those scores are partial log-odds likelihoods in nats.
       * Subtract off the rest of the null model, convert to bits.
       */
      p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

      if (esl_opt_GetBoolean(go, "-1")) 
	{
	  printf("%-30s   %10.4f %10.4f   %10.4f %10.4f\n", 
		 sq->name, 
		 fsc, bsc, 
		 (fsc - nullsc) / eslCONST_LOG2, (bsc - nullsc) / eslCONST_LOG2);
	}
      else
	{
	  printf("target sequence:      %s\n",        sq->name);
	  printf("fwd raw score:        %.4f nats\n", fsc);
	  printf("bck raw score:        %.4f nats\n", bsc);
	  printf("null score:           %.2f nats\n", nullsc);
	  printf("per-seq score:        %.2f bits\n", (fsc - nullsc) / eslCONST_LOG2);
	  printf("RAM usage in DP:      %.3fM\n",     (double) (p7_bandmx_Sizeof(bmx) / 1000000));
	}

      p7_filtermx_Reuse(ox);
      p7_bandmx_Reuse(bmx);
      p7_gbands_Reuse(bnd);
      esl_sq_Reuse(sq);
    }

  /* Cleanup */
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_gbands_Destroy(bnd);
  p7_bandmx_Destroy(bmx);
  p7_filtermx_Destroy(ox);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7BANDED_FWDBACK_EXAMPLE*/

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
