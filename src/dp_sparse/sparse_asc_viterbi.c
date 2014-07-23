/* 
 */

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_coords2.h"

#include "misc/logsum.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/p7_spascmx.h"

int
p7_sparse_asc_Viterbi(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D, 
		      const P7_SPARSEMASK *sm, P7_SPARSEMX *asv, float *opt_sc)
{
  float const *tsc = gm->tsc;	 // sets up TSC() macro, access to profile's transitions      
  float const *rsc = NULL;	 // will be set up for MSC(), ISC() macros for residue scores for each row i
  int    g         = 1;	         // index of next or current segment. When we enter it, and while we're in it, in_seg = TRUE. s=1..S, with sentinels.
  int    in_seg    = FALSE;      //   ... this bumps to TRUE when we see ia(g), first row of segment; to FALSE on ib(g), last row.
  int    d         = 1;          // idx of next domain anchor we reach. Row i may be in sector UP(d), DOWN(d-1).
  int    ndown     = 0;          //   ... this bumps to 1 when i reaches an anchor; row # of DOWN sector; back to 0 at seg end
  float *dpc       = NULL;       // all indexing is relative! This ptr walks through the sparse ASC matrix main cells.
  float *xc        = NULL;       //   ... and this one, through the specials.
  float *dpp       = NULL;       // ptr for stepping thru DP cells on a previous row. Gets initialized by...
  float *last_down = NULL;       //   ... remembering the <dpc> start of each DOWN...
  float *last_up   = NULL;       //   ... or UP sector as we compute them.
  float  mlc, mgc;               // tmp variables for computing MLk, MGk score
  float  dlc, dgc;               //   ... and DLk, DGk.
  float  xE,xN,xJ,xB,xL,xG,xC;   // current special state values on this row; specials are not always stored
  int i;                         // index for dsq[i] seq position/rows 0..L
  int k;                         // index of model positions 1..M
  int z;                         // index for sparse cell indices k[i][z] along a row
  int y;			 //   ... and for previous row, k[i-1][y].
  int s;                         // index of individual states
  int ngap;
  int status;

  /* Contract checks, argument validation */
  ESL_DASSERT1(( sm->L == L     ));
  ESL_DASSERT1(( sm->M == gm->M ));

  /* Assure that <asv> is allocated large enough. */
  if (( status = p7_spascmx_Reinit(asv, sm, anch, D)) != eslOK) return status;
  asv->type = p7S_ASC_VITERBI;
  
  dpc = asv->dp;   // you can't initialize <dpc> until AFTER the _Reinit() call above, because Reinit might move the memory.
  xc  = asv->xmx;  //   ... ditto for <xc>.
  xN  = 0.;
  for (i = 0; i <= L; i++)
    {
      /* Mechanics of traversing a sparse ASC matrix, row by row */
      if      (i == anch[d].i0) { ndown = 1; d++; }  // when i reaches next anchor; bump d to next domain index, and DOWN sector is active...
      else if (sm->n[i] == 0)   { ndown = 0;      }  //  ... until when we reach end of segment, when DOWN becomes inactive again.
      else if (ndown)           { ndown++;        }  // counting ndown lets us determine when we're in 1st, 2nd row; useful for boundary conditions on UP.

      if      (i >  sm->seg[g].ib)  { in_seg  = FALSE; g++; }  // g bumps, to start expecting to see start of segment <g> next.
      else if (i == sm->seg[g].ia)  { in_seg  = TRUE;       }  // g might be S+1, but this is safe because of sentinel seg[S+1].ia=ib=L+2      

      xE = -eslINFINITY;

      /*****************************************************************
       * Computing row i when it is in a DOWN sector
       *****************************************************************/
      if (ndown == 1)                                // topmost row of a DOWN sector is a special case for initialization.
	{
	  rsc       = gm->rsc[dsq[i]];	             // now MSC(k), ISC(k) residue score macros work 
	  dpp       = last_up;                       // yes, UP: this initialization of i0,k0 in DOWN is the connecting thread from UP to DOWN.
	  last_down = dpc;
	  y = z     = 0;

	  while ( sm->k[i][z] < anch[d-1].k0) z++;   // skip z to the anchor cell for row i. (Which MUST exist; we don't even need to check z<n[i])
	  k = sm->k[i][z];

	  mlc = xL + TSC(p7P_LM, k-1);
	  mgc = xG + TSC(p7P_GM, k-1);
	  while (dpp && y < sm->n[i-1] && sm->k[i-1][y] <  k-1) { y++; dpp += p7S_NSCELLS; }  // dpp is an UP row. i-1,k-1 may not exist (UP sector can be empty), hence the check for non-NULL dpp
	  if    (dpp && y < sm->n[i-1] && sm->k[i-1][y] == k-1) {                             //   ... but usually, there's an UP sector with a i-1,k-1 cell 
	    mlc = ESL_MAX( ESL_MAX( dpp[p7R_ML] + TSC(p7P_MM, k-1),
				    dpp[p7R_IL] + TSC(p7P_IM, k-1)),
			   ESL_MAX( dpp[p7R_DL] + TSC(p7P_DM, k-1),
				    mlc));
	    mgc = ESL_MAX( ESL_MAX( dpp[p7R_MG] + TSC(p7P_MM, k-1),
				    dpp[p7R_IG] + TSC(p7P_IM, k-1)),
			   ESL_MAX( dpp[p7R_DG] + TSC(p7P_DM, k-1),
				    mgc));
	  }
	  *dpc++ = mlc = MSC(k) + mlc;
	  *dpc++ = mgc = MSC(k) + mgc;
	  *dpc++ = -eslINFINITY;       // IL(i0,k0)
	  *dpc++ = -eslINFINITY;       // IG(i0,k0)
	  xE     = ESL_MAX(xE, mlc);               // dlc is -inf, not included in sum
	  *dpc++ = -eslINFINITY;       // DL(i0,k0)
	  *dpc++ = -eslINFINITY;       // DG(i0,k0)
	  if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) {                     // is there a (i,k+1) cell to our right? 
	    dlc =  mlc + TSC(p7P_MD, k);                                    //  ... then ->Dk+1 is precalculated as usual.
	    dgc =  mgc + TSC(p7P_MD, k);                                    // Since we know we're the first cell (the anchor) we initialize dlc, dgc here.
	  } else {                                                          // If not, we MUST add {MD}l->Dk+1..E glocal exit path, even from internal sparse cells - not just last cell! 
	    xE  = ESL_MAX( xE, mgc + TSC(p7P_DGE, k) + TSC(p7P_MD, k));
	    dlc = dgc = -eslINFINITY;
	  }
	  
	  for (z = z+1; z < sm->n[i]; z++)                                  // Now the rest of the anchor row for DOWN. Can only be a ->DDD->E path.
	    {
	      k = sm->k[i][z];
	      
	      *dpc++ = -eslINFINITY;
	      *dpc++ = -eslINFINITY;
	      *dpc++ = -eslINFINITY;
	      *dpc++ = -eslINFINITY;
	      xE     = ESL_MAX(xE, dlc);
	      *dpc++ = dlc;
	      *dpc++ = dgc;
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) {               // see comment above for ->Dk+1 precalculation.
		dlc = dlc + TSC(p7P_DD, k);
		dgc = dgc + TSC(p7P_DD, k);
	      } else {   
		xE  = ESL_MAX( xE, dgc + TSC(p7P_DGE, k) + TSC(p7P_DD, k)); 
		dlc = dgc = -eslINFINITY;
	      }
	    }
	}

      else if (ndown)                                             // Main recursion for DOWN rows. 
	{
	  rsc = gm->rsc[dsq[i]];	                          // now MSC(k), ISC(k) residue score macros work 
	  y = z = 0; 
	  dlc = dgc = -eslINFINITY;

	  dpp        = last_down;
	  last_down  = dpc;

	  while (z < sm->n[i] && sm->k[i][z] < anch[d-1].k0) z++;   // skip sparse cells that aren't in DOWN sector (which is k0..M) 
	  for (; z < sm->n[i]; z++)                                 // then calculate the rest of the row.
	    {                                  
	      k = sm->k[i][z];  // for notational convenience

	      /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	      mlc = mgc = -eslINFINITY;
	      while (y < sm->n[i-1] && sm->k[i-1][y] < anch[d-1].k0) y++;                           // skip cells on prev row that aren't in DOWN at all
	      while (y < sm->n[i-1] && sm->k[i-1][y] < k-1)        { y++; dpp += p7S_NSCELLS; }     // skip cells that exist in sparse ASC matrix, but aren't (i-1,k-1)
	      if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) {
		mlc = ESL_MAX( ESL_MAX( dpp[p7R_ML] + TSC(p7P_MM, k-1),
					dpp[p7R_IL] + TSC(p7P_IM, k-1)),
			                dpp[p7R_DL] + TSC(p7P_DM, k-1));
		mgc = ESL_MAX( ESL_MAX( dpp[p7R_MG] + TSC(p7P_MM, k-1),
					dpp[p7R_IG] + TSC(p7P_IM, k-1)),
			                dpp[p7R_DG] + TSC(p7P_DM, k-1));
	      }
	      *dpc++ = mlc = MSC(k) + mlc;
	      *dpc++ = mgc = MSC(k) + mgc;

	      /* Try to find cell i-1,k; then compute I(i,k) from it */
	      if (y < sm->n[i-1] && sm->k[i-1][y] < k) { y++; dpp += p7S_NSCELLS; }
	      if (y < sm->n[i-1] && sm->k[i-1][y] == k) {
		*dpc++ = ESL_MAX( dpp[p7R_ML] + TSC(p7P_MI,k),  dpp[p7R_IL] + TSC(p7P_II, k)); // +ISC(k) if we weren't enforcing it to zero
		*dpc++ = ESL_MAX( dpp[p7R_MG] + TSC(p7P_MI,k),  dpp[p7R_IG] + TSC(p7P_II, k)); // ditto
	      } else {
		*dpc++ = -eslINFINITY;
		*dpc++ = -eslINFINITY;
	      }

	      /* Local exit paths */
	      xE = ESL_MAX(xE, ESL_MAX(mlc, dlc));

	      /* Delayed store of Dk. */
	      *dpc++ = dlc;
	      *dpc++ = dgc;

	      /* Advance calculation of next D_k+1 */
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) { 
		dlc = ESL_MAX( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
		dgc = ESL_MAX( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	      } else {
		xE  = ESL_MAX( xE, TSC(p7P_DGE, k) + ESL_MAX( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k))); 
		dlc = dgc = -eslINFINITY;
	      }
	    } // end loop over z: we've finished sparse row i for DOWN sector.
	} // end of <ndown> block, calculation of DOWN row i
	  

      


      /*****************************************************************
       * Computing row i when it's in an UP sector 
       *****************************************************************/

      if (in_seg && anch[d].i0 <= sm->seg[g].ib)                             // d may be D+1 here: if so, sentinel makes the the comparison to seg[g].ib fail
	{
	  if (ndown == 1)       // The first row of an UP matrix for subsequent domains in a segment is the anchor row i0.
	    {                   // All sparse cells in this row are unreachable, initialized to -inf. They only get used for G->DDDD->Mk,i+1 decoding, wing unfolding.
	      last_up = dpc;
	      for (z = 0; z < sm->n[i] && sm->k[i][z] < anch[d].k0; z++)
		for (s = 0; s < p7S_NSCELLS; s++)
		  *dpc++ = -eslINFINITY;
	    }

	  else if (i == sm->seg[g].ia)         // The first row of UP(d) when d is the first domain in a segment
	    {                                  // is the first row of the segment. Cells on this row can be reached
	      rsc       = gm->rsc[dsq[i]];     // by {GL}->Mk entries, followed by {Dk,Mk}->Dk+1 delete transitions.
	      dpp       = NULL;                // The (ndown==1) code must precede this, to deal with a case of an anchor on first row of a segment, which means a nonexistent UP sector and immediate init in DOWN
	      last_up   = dpc;
	      dlc = dgc = -eslINFINITY;

	      for (z=0, y=0; z < sm->n[i] && sm->k[i][z] < anch[d].k0; z++)  // for all sparse cells in UP sector on this row
		{
		  k = sm->k[i][z];

		  *dpc++ = mlc = xL  + TSC(p7P_LM, k-1) + MSC(k);
		  *dpc++ = mgc = xG  + TSC(p7P_GM, k-1) + MSC(k);
		  *dpc++       = -eslINFINITY;
		  *dpc++       = -eslINFINITY;
		  *dpc++       = dlc;
		  *dpc++       = dgc;
		  // Advance calculation of next D_k+1, if it's there
		  if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) {  // is there an (i,k+1) cell to our right?
		    dlc = ESL_MAX( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
		    dgc = ESL_MAX( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
		  } else 
		    dlc = dgc = -eslINFINITY;
		}
	    }

	  else                      // all other cases: we're within UP(d), not on its first row; 
	    {                       // so i-1 row exists in UP(d) and sparse mask; we may look at dp[i-1] and sm->k[i-1]
	      rsc = gm->rsc[dsq[i]];	                   // now MSC(k), ISC(k) residue score macros work 
	      dpp       = last_up;
	      last_up   = dpc;
	      dlc = dgc = -eslINFINITY;

	      for (z=0, y=0; z < sm->n[i] && sm->k[i][z] < anch[d].k0; z++)  // for all sparse cells in UP sector on this row
		{
		  k = sm->k[i][z];

		  // Try to find cell i-1,k-1. Then compute M(i,k) from it.
		  mlc = xL  + TSC(p7P_LM, k-1);
		  mgc = xG  + TSC(p7P_GM, k-1);
		  while (y < sm->n[i-1] && sm->k[i-1][y] < k-1)  { y++; dpp+= p7S_NSCELLS; }
		  if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) {
		    mlc = ESL_MAX( ESL_MAX( dpp[p7R_ML] + TSC(p7P_MM, k-1),
					    dpp[p7R_IL] + TSC(p7P_IM, k-1)),
				   ESL_MAX( dpp[p7R_DL] + TSC(p7P_DM, k-1),
					    mlc));
		    mgc = ESL_MAX( ESL_MAX( dpp[p7R_MG] + TSC(p7P_MM, k-1),
					    dpp[p7R_IG] + TSC(p7P_IM, k-1)),
				   ESL_MAX( dpp[p7R_DG] + TSC(p7P_DM, k-1),
					    mgc));
		  }
		  *dpc++ = mlc = MSC(k) + mlc;
		  *dpc++ = mgc = MSC(k) + mgc;

		  // Try to find cell i-1, k. Then compute I(i,k) from it.
		  if (y < sm->n[i-1] && sm->k[i-1][y] < k) { y++; dpp += p7S_NSCELLS; }
		  if (y < sm->n[i-1] && sm->k[i-1][y] == k) {
		    *dpc++ = ESL_MAX( dpp[p7R_ML] + TSC(p7P_MI,k),  dpp[p7R_IL] + TSC(p7P_II, k)); // +ISC(k) if we weren't enforcing it to zero
		    *dpc++ = ESL_MAX( dpp[p7R_MG] + TSC(p7P_MI,k),  dpp[p7R_IG] + TSC(p7P_II, k)); // ditto
		  } else {
		    *dpc++ = -eslINFINITY;
		    *dpc++ = -eslINFINITY;
		  }

		  // Delayed store of Dk
		  *dpc++ = dlc;
		  *dpc++ = dgc;

		  // Advance calculation of next D_k+1, if it's there
		  if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) {  // is there an (i,k+1) cell to our right?
		    dlc = ESL_MAX( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
		    dgc = ESL_MAX( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
		  } else 
		    dlc = dgc = -eslINFINITY;
	      
		} // End of loop over z, sparse cells on row i
	    } // End of dealing with all but the first row in an UP sector
	} // End of UP sector.

      /*****************************************************************
       * Specials. Always calculated; sometimes stored.
       *****************************************************************/

      if (i == sm->seg[g].ia-1)
	{
	  ngap =  (g == 1 ? i : sm->seg[g].ia - sm->seg[g-1].ib - 1);

	  if      (d == 1)             xN =  xN + (ngap == 0 ? 0. : ngap * gm->xsc[p7P_N][p7P_LOOP]);
	  else                         xN = -eslINFINITY;

	  if      (d == 1 || d == D+1) xJ = -eslINFINITY;
	  else                         xJ = ESL_MAX(xJ + (float) ngap * gm->xsc[p7P_J][p7P_LOOP], 
						    xE +                gm->xsc[p7P_E][p7P_LOOP]);

	  if      (d == D+1)           xC = ESL_MAX(xC + (float) ngap * gm->xsc[p7P_C][p7P_LOOP], 
						    xE +                gm->xsc[p7P_E][p7P_MOVE]);
	  else                         xC = -eslINFINITY;

	  xB = (d == 1 ? xN + gm->xsc[p7P_N][p7P_MOVE] : xJ + gm->xsc[p7P_J][p7P_MOVE]);
	  xL = xB + gm->xsc[p7P_B][0]; /* B->L */
	  xG = xB + gm->xsc[p7P_B][1]; /* B->G */
	}
      else if (i >= sm->seg[g].ia && i <= sm->seg[g].ib)
	{
	  xN = (d == 1 ? xN + gm->xsc[p7P_N][p7P_LOOP] : -eslINFINITY);

	  if      (d == 1 || d == D+1) xJ = -eslINFINITY;
	  else if (ndown == 1)         xJ = xE + gm->xsc[p7P_E][p7P_LOOP];  // don't propagate J->J across DOWN sector boundaries
	  else                         xJ = ESL_MAX(xJ + gm->xsc[p7P_J][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_LOOP]);

	  xB = (d == 1 ? xN + gm->xsc[p7P_N][p7P_MOVE] : xJ + gm->xsc[p7P_J][p7P_MOVE]);
	  xL = xB + gm->xsc[p7P_B][0]; /* B->L */
	  xG = xB + gm->xsc[p7P_B][1]; /* B->G */
	  xC = (d == D+1 ? ESL_MAX(xE + gm->xsc[p7P_E][p7P_MOVE],  xC + gm->xsc[p7P_C][p7P_LOOP]) : -eslINFINITY);
	}
      
      /* ASC only stores specials on ia(g)-1..ib(g) for segments g that contain >= 1 anchor. simple to say, a little awkward to test for: */
      if ( (i >= sm->seg[g].ia-1 && anch[d].i0 <= sm->seg[g].ib) || ndown)   // d is next anchor we'll see. i0(d) <= ib until we're in last DOWN sector in seg.
	{                                                                    // in last DOWN sector in seg, i0(d) is in a later seg... hence we need the <ndown> test too here
	  *xc++ = xE;
	  *xc++ = xN;
	  *xc++ = xJ;
	  *xc++ = xB;
	  *xc++ = xL;
	  *xc++ = xG;
	  *xc++ = xC;
	  *xc++ = -eslINFINITY; // JJ; only valid in a decoding mx
	  *xc++ = -eslINFINITY; // CC; ditto
	}
    } // end of i=0..L. Sparse ASC DP matrix complete. Mission accomplished.

  ngap = L - sm->seg[sm->S].ib;
  xC   = (ngap == 0 ? xC : xC + (float) ngap * gm->xsc[p7P_C][p7P_LOOP]);
  if (opt_sc) *opt_sc = xC + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}

/*****************************************************************
 * x. Example
 *****************************************************************/ 
#ifdef p7SPARSE_ASC_VITERBI_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "include all cells in sparse mx",                   0 },
  { "-s",         eslARG_INT,    "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
 
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of sparse ASC Viterbi";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_CHECKPTMX   *cx      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  P7_TRACE       *vtr     = p7_trace_Create();
  P7_SPARSEMX    *sxf     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxd     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *asv     = p7_sparsemx_Create(NULL);
  P7_ANCHORS     *anch    = p7_anchors_Create();
  float           fsc, vsc, asc;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
 
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Get a sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create (hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config   (gm, hmm, bg);
  p7_oprofile_Convert (gm, om);

  p7_bg_SetLength     (bg, sq->n);
  p7_profile_SetLength(gm, sq->n);
  p7_oprofile_ReconfigLength(om, sq->n);

  /* We need a sparse mask <sm>.
   * To get it, run checkpointed Fwd/Bck/Decoding
   */
  cx = p7_checkptmx_Create(hmm->M, sq->n, ESL_MBYTES(32));
  sm = p7_sparsemask_Create(gm->M, sq->n);
  if (esl_opt_GetBoolean(go, "-a")) 
    p7_sparsemask_AddAll(sm);
  else {
    p7_ForwardFilter (sq->dsq, sq->n, om, cx, /*fsc=*/NULL);
    p7_BackwardFilter(sq->dsq, sq->n, om, cx, sm, p7_SPARSEMASK_THRESH_DEFAULT);
  }

  /* We need an anchor set <anch>.
   * To get it, get the sparse Viterbi trace, and select an anchor.
   * To select the anchor, we also have to decode.
   */
  p7_SparseViterbi (sq->dsq, sq->n, gm, sm, sxf, vtr, &vsc);
  p7_SparseForward (sq->dsq, sq->n, gm, sm, sxf,      &fsc);
  p7_SparseBackward(sq->dsq, sq->n, gm, sm, sxd,      &fsc);
  p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxd, sxd);
  p7_sparse_anchors_SetFromTrace(sxd, vtr, anch);

  /* Finally...
   * Run sparse ASC Viterbi.
   */
  p7_sparse_asc_Viterbi(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asv, &asc);
  
  printf("Sparse Viterbi score     = %.2f nats\n", vsc);
  printf("Sparse ASC Viterbi score = %.2f nats\n", asc);

  p7_anchors_Destroy(anch);
  p7_sparsemx_Destroy(asv);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemx_Destroy(sxf);
  p7_trace_Destroy(vtr);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(cx);
  p7_profile_Destroy(gm);  p7_oprofile_Destroy(om);
  esl_sqfile_Close(sqfp);  esl_sq_Destroy(sq);
  p7_hmmfile_Close(hfp);   p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPARSE_ASC_VITERBI_EXAMPLE*/
/*---------------------- end, example ---------------------------*/
