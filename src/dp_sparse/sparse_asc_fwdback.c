/* Production implementation of anchor set constrained (ASC) Forward
 * and Backward, using sparse dynamic programming.
 * 
 * Contents:
 *    1. Sparse ASC Forward
 *    2. Sparse ASC Backward
 *    3. Unit tests
 *    4. Test driver
 *    5. Example
 *    6. Notes
 *    7. Copyright and license information
 */

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_coords2.h"

#include "misc/logsum.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/p7_spascmx.h"


/*****************************************************************
 * 1. Sparse ASC Forward 
 *****************************************************************/

/* Function:  
 * Synopsis:  
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:  
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      
 */
int
p7_sparse_asc_Forward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_COORD2 *anch, int D, 
		      const P7_SPARSEMASK *sm, P7_SPARSEMX *asf, float *opt_sc)
{
  float const *tsc = gm->tsc;	 // sets up TSC() macro, access to profile's transitions      
  float const *rsc = NULL;	 // will be set up for MSC(), ISC() macros for residue scores for each row i
  int    g         = 1;	         // index of next or current segment. When we enter it, and while we're in it, in_seg = TRUE. s=1..S, with sentinels.
  int    in_seg    = FALSE;      //   ... this bumps to TRUE when we see ia(g), first row of segment; to FALSE on ib(g), last row.
  int    d         = 0;          // how many domain anchors we've reached so far. Row i may be in sector UP(d), DOWN(d-1).
  int    ndown     = 0;          //   ... this bumps to 1 when i reaches an anchor; row # of DOWN sector; back to 0 at seg end
  float *dpc       = asf->dp;    // all indexing is relative! This ptr walks through the sparse ASC matrix main cells.
  float *xc        = asf->xmx;   //   ... and this one, through the specials.
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
  int status;

  /* Contract checks, argument validation */
  ESL_DASSERT1(( sm->L == L     ));
  ESL_DASSERT1(( sm->M == gm->M ));

  /* Assure that <asf> is allocated large enough. */
  if (( status = p7_spascmx_Reinit(asf, sm, anch, D)) != eslOK) return status;
  asf->type = p7S_ASC_FWD;

  for (i = 0; i <= L; i++)
    {
      /* Mechanics of traversing a sparse ASC matrix, row by row */
      if      (d < D && i == anch[d].n1) { ndown = 1; d++; }  // when i reaches next anchor; bump d to next domain index, and DOWN sector is active...
      else if (sm->n[i] == 0)            { ndown = 0;      }  //  ... until when we reach end of segment, when DOWN becomes inactive again.
      else if (ndown)                    { ndown++;        }  // counting ndown lets us determine when we're in 1st, 2nd row; useful for boundary conditions on UP.

      if      (i >  sm->seg[g].ib)  { in_seg  = FALSE; g++; }  // g bumps, to start expecting to see start of segment <g> next.
      else if (i == sm->seg[g].ia)  { in_seg  = TRUE;       }  // g might be S+1, but this is safe because of sentinel seg[S+1].ia=ib=L+2      

      xE = -eslINFINITY;

      /*****************************************************************
       * Computing row i when it is in a DOWN sector
       *****************************************************************/
      if (ndown == 1)                                    // topmost row of a DOWN sector is a special case for initialization.
	{
	  rsc       = gm->rsc[dsq[i]];	                   // now MSC(k), ISC(k) residue score macros work 
	  dpp       = last_up;                             // yes, UP: this initialization of i0,k0 in DOWN is the connecting thread from UP to DOWN.
	  last_down = dpc;
	  y = z     = 0;

	  while ( sm->k[i][z]   <  anch[d-1].n2) z++;                                  // skip z to the anchor cell for row i. (Which MUST exist; we don't even need to check z<n[i])
	  k = sm->k[i][z];

	  mlc = xL + TSC(p7P_LM, k-1);
	  mgc = xG + TSC(p7P_GM, k-1);
	  while (y < sm->n[i-1] && sm->k[i-1][y] <  k-1) { y++; dpp += p7S_NSCELLS; }  // dpp is an UP row. i-1,k-1 may not exist (UP sector can be empty)
	  if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) {                             //   ... but usually, there's an UP sector with a i-1,k-1 cell 
	    mlc = p7_FLogsum( p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MM, k-1),
					  dpp[p7R_IL] + TSC(p7P_IM, k-1)),
			      p7_FLogsum( dpp[p7R_DL] + TSC(p7P_DM, k-1),
					  mlc));
	    mgc = p7_FLogsum( p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MM, k-1),
					  dpp[p7R_IG] + TSC(p7P_IM, k-1)),
			      p7_FLogsum( dpp[p7R_DG] + TSC(p7P_DM, k-1),
					  mgc));
	  }
	  *dpc++ = mlc = MSC(k) + mlc;
	  *dpc++ = mgc = MSC(k) + mgc;
	  *dpc++ = -eslINFINITY;       // IL(i0,k0)
	  *dpc++ = -eslINFINITY;       // IG(i0,k0)
	  xE     = p7_FLogsum(xE, mlc);
	  *dpc++ = -eslINFINITY;       // DL(i0,k0)
	  *dpc++ = -eslINFINITY;       // DG(i0,k0)
	  if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) {                     // is there a (i,k+1) cell to our right? 
	    dlc =  mlc + TSC(p7P_MD, k);                                    //  ... then ->Dk+1 is precalculated as usual.
	    dgc =  mgc + TSC(p7P_MD, k);                                    // Since we know we're the first cell (the anchor) we initialize dlc, dgc here.
	  } else {                                                          // If not, we MUST add {MD}l->Dk+1..E glocal exit path, even from internal sparse cells - not just last cell! 
	    xE  = p7_FLogsum( xE, dgc + TSC(p7P_DGE, k) + TSC(p7P_MD, k));  // this line must come BEFORE you set dgc to -inf.
	    dlc = dgc = -eslINFINITY;
	  }
	  
	  for (z = z+1; z < sm->n[i]; z++)                                  // Now the rest of the anchor row for DOWN. Can only be a ->DDD->E path.
	    {
	      k = sm->k[i][z];
	      
	      *dpc++ = -eslINFINITY;
	      *dpc++ = -eslINFINITY;
	      *dpc++ = -eslINFINITY;
	      *dpc++ = -eslINFINITY;
	      xE     = p7_FLogsum(xE, dlc);
	      *dpc++ = dlc;
	      *dpc++ = dgc;
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) {               // see comment above for ->Dk+1 precalculation.
		dlc = dlc + TSC(p7P_DD, k);
		dgc = dgc + TSC(p7P_DD, k);
	      } else {   
		xE  = p7_FLogsum( xE, dgc + TSC(p7P_DGE, k) + TSC(p7P_DD, k)); 
		dlc = dgc = -eslINFINITY;
	      }
	    }
	}

      else if (ndown)                                             // Main recursion for DOWN rows. 
	{
	  rsc = gm->rsc[dsq[i]];	                            // now MSC(k), ISC(k) residue score macros work 
	  y = z = 0; 
	  dlc = dgc = -eslINFINITY;

	  dpp        = last_down;
	  last_down  = dpc;

	  while (z < sm->n[i] && sm->k[i][z] < anch[d-1].n2) z++;   // skip sparse cells that aren't in DOWN sector (which is k0..M) 
	  for (; z < sm->n[i]; z++)                                 // then calculate the rest of the row.
	    {                                                       //   note: d-1 is safe, because ndown can only be TRUE for d>=1
	      k = sm->k[i][z];  // for notational convenience

	      /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	      mlc = mgc = -eslINFINITY;
	      while (y < sm->n[i-1] && sm->k[i-1][y] < anch[d-1].n2) y++;                           // skip cells on prev row that aren't in DOWN at all
	      while (y < sm->n[i-1] && sm->k[i-1][y] < k-1)        { y++; dpp += p7S_NSCELLS; }     // skip cells that exist in sparse ASC matrix, but aren't (i-1,k-1)
	      if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) {
		mlc = p7_FLogsum( p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MM, k-1),
					      dpp[p7R_IL] + TSC(p7P_IM, k-1)),
				              dpp[p7R_DL] + TSC(p7P_DM, k-1));
		mgc = p7_FLogsum( p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MM, k-1),
					      dpp[p7R_IG] + TSC(p7P_IM, k-1)),
				              dpp[p7R_DG] + TSC(p7P_DM, k-1));
	      }
	      *dpc++ = mlc = MSC(k) + mlc;
	      *dpc++ = mgc = MSC(k) + mgc;

	      /* Try to find cell i-1,k; then compute I(i,k) from it */
	      if (y < sm->n[i-1] && sm->k[i-1][y] < k) { y++; dpp += p7S_NSCELLS; }
	      if (y < sm->n[i-1] && sm->k[i-1][y] == k) {
		*dpc++ = p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MI,k),  dpp[p7R_IL] + TSC(p7P_II, k)); // +ISC(k) if we weren't enforcing it to zero
		*dpc++ = p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MI,k),  dpp[p7R_IG] + TSC(p7P_II, k)); // ditto
	      } else {
		*dpc++ = -eslINFINITY;
		*dpc++ = -eslINFINITY;
	      }

	      /* Local exit paths */
	      xE = p7_FLogsum(xE, p7_FLogsum(mlc, dlc));

	      /* Delayed store of Dk. */
	      *dpc++ = dlc;
	      *dpc++ = dgc;

	      /* Advance calculation of next D_k+1 */
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) { 
		dlc = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
		dgc = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	      } else {
		xE  = p7_FLogsum( xE, TSC(p7P_DGE, k) + p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k))); 
		dlc = dgc = -eslINFINITY;
	      }
	    } // end loop over z: we've finished sparse row i for DOWN sector.
	} // end of <ndown> block, calculation of DOWN row i
	  

      


      /*****************************************************************
       * Computing row i when it's in an UP sector 
       *****************************************************************/

      if (i == sm->seg[g].ia || (ndown == 2 && d < D))                   // top row of UP. No previous row, so you can't enter from M,I; only by {LG}->Mk entry, followed by ->DDD paths.
	{
	  rsc       = gm->rsc[dsq[i]];
	  dpp       = NULL;
	  last_up   = dpc;
	  dlc = dgc = -eslINFINITY;

	  for (z=0, y=0; z < sm->n[i] && sm->k[i][z] < anch[d].n2; z++)  // for all sparse cells in UP sector on this row
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
		dlc = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
		dgc = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	      } else 
		dlc = dgc = -eslINFINITY;
	    }
	}

      else if (in_seg && d < D && anch[d].n1 <= sm->seg[g].ib)      // if i is in UP sector [d]. d may be =D here, hence the guard d<D before anch[d] access
	{
	  rsc = gm->rsc[dsq[i]];	                            // now MSC(k), ISC(k) residue score macros work 
	  dpp       = last_up;
	  last_up   = dpc;
	  dlc = dgc = -eslINFINITY;

	  for (z=0, y=0; z < sm->n[i] && sm->k[i][z] < anch[d].n2; z++)  // for all sparse cells in UP sector on this row
	    {
	      k = sm->k[i][z];

	      // Try to find cell i-1,k-1. Then compute M(i,k) from it.
	      mlc = xL  + TSC(p7P_LM, k-1);
	      mgc = xG  + TSC(p7P_GM, k-1);
	      while (y < sm->n[i-1] && sm->k[i-1][y] < k-1)  { y++; dpp+= p7S_NSCELLS; }
	      if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) {
		mlc = p7_FLogsum( p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MM, k-1),
					      dpp[p7R_IL] + TSC(p7P_IM, k-1)),
				  p7_FLogsum( dpp[p7R_DL] + TSC(p7P_DM, k-1),
					      mlc));
		mgc = p7_FLogsum( p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MM, k-1),
					      dpp[p7R_IG] + TSC(p7P_IM, k-1)),
				  p7_FLogsum( dpp[p7R_DG] + TSC(p7P_DM, k-1),
					      mgc));
	      }
	      *dpc++ = mlc = MSC(k) + mlc;
	      *dpc++ = mgc = MSC(k) + mgc;

	      // Try to find cell i-1, k. Then compute I(i,k) from it.
	      if (y < sm->n[i-1] && sm->k[i-1][y] < k) { y++; dpp += p7S_NSCELLS; }
	      if (y < sm->n[i-1] && sm->k[i-1][y] == k) {
		*dpc++ = p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MI,k),  dpp[p7R_IL] + TSC(p7P_II, k)); // +ISC(k) if we weren't enforcing it to zero
		*dpc++ = p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MI,k),  dpp[p7R_IG] + TSC(p7P_II, k)); // ditto
	      } else {
		*dpc++ = -eslINFINITY;
		*dpc++ = -eslINFINITY;
	      }

	      // Delayed store of Dk
	      *dpc++ = dlc;
	      *dpc++ = dgc;

	      // Advance calculation of next D_k+1, if it's there
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) {  // is there an (i,k+1) cell to our right?
		dlc = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
		dgc = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	      } else 
		dlc = dgc = -eslINFINITY;
	      
	    } // End of loop over z, sparse cells on row i


	} // End of row i in UP sector


      /*****************************************************************
       * Specials. Always calculated; sometimes stored.
       *****************************************************************/

      xN = (i == 0 ? 0. : (d == 0 ? xN + gm->xsc[p7P_N][p7P_LOOP] : -eslINFINITY));
      xJ = (d > 0 && d < D  ? p7_FLogsum(xJ + gm->xsc[p7P_J][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_LOOP]) : -eslINFINITY);
      xB = (d == 0 ? xN + gm->xsc[p7P_N][p7P_MOVE] : xJ + gm->xsc[p7P_J][p7P_MOVE]);
      xL = xB + gm->xsc[p7P_B][0]; /* B->L */
      xG = xB + gm->xsc[p7P_B][1]; /* B->G */
      xC = (d == D ? p7_FLogsum(xE + gm->xsc[p7P_E][p7P_MOVE],  xC + gm->xsc[p7P_C][p7P_LOOP]) : -eslINFINITY);
      
      if ( (i >= sm->seg[g].ia-1 && d < D && anch[d].n1 <= sm->seg[g].ib) || ndown)
	{
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

  if (opt_sc) *opt_sc = xC + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}
/*--------------- end, sparse ASC Forward -----------------------*/



/*****************************************************************
 * 2. Sparse ASC Backward
 *****************************************************************/

int
p7_sparse_asc_Backward(void)
{

  return eslOK;
}


/*-------------- end, sparse ASC Backward -----------------------*/


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7SPARSE_ASC_FWDBACK_TESTDRIVE


/* "compare_reference" unit test.
 * 
 * When we include all i,k supercells in the sparse mask, then sparse
 * DP calculations give the same results as reference calculations,
 * even at the level of individual DP cells.
 * 
 * Sample a random profile of length <M>.  Generate <N> sequences from
 * that profile, using a length model of <L> during sampling.  For
 * each sampled sequence, Make a sparse mask that contains all i,k
 * cells. Make an anchor set from the generating trace (the anchor set
 * just needs to be reasonable, not optimal).
 * 
 * Then:
 *   1. Reference and sparse Forward scores must be identical (within 
 *      numerical tolerance).
 *   2. Cells of reference and sparse DP matrix have identical values
 *      (within numerical tolerance).
 *   3. Sparse Fwd matrix structure passes Validate().
 */
static void
utest_compare_reference(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int L, int N)
{
  char           msg[] = "sparse_asc_fwdback :: compare_reference unit test failed";
  P7_BG         *bg    = p7_bg_Create(abc);
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = p7_profile_Create(M, abc);
  P7_SPARSEMASK *sm    = p7_sparsemask_Create(M, L);
  P7_ANCHORS    *anch  = p7_anchors_Create();
  P7_REFMX      *afu   = p7_refmx_Create(100,100);
  P7_REFMX      *afd   = p7_refmx_Create(100,100);
  P7_SPARSEMX   *asf   = p7_sparsemx_Create(NULL);
  float          sc1, sc2;
  int            idx;

  /* Sample a profile. Config as usual, multihit dual-mode local/glocal. */
  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);

  for (idx = 0; idx < N; idx++)
    {
      /* Generate (sample) a sequence from the profile */
      if ( p7_profile_SetLength(gm, L)  != eslOK) esl_fatal(msg);   /* config to generate mean length of L (length was probably reset by last emitted seq) */
      do {
	esl_sq_Reuse(sq);
	p7_ProfileEmit(rng, hmm, gm, bg, sq, gtr);
      } while (sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
      if ( p7_profile_SetLength(gm, sq->n) != eslOK) esl_fatal(msg);

      /* Mark all cells in sparse mask */
      if ( p7_sparsemask_Reinit(sm, M, sq->n) != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_AddAll(sm)           != eslOK) esl_fatal(msg);

      /* Use generating trace to create a plausible anchor set */
      if ( p7_anchors_SampleFromTrace(rng, gtr, anch) != eslOK) esl_fatal(msg);

      /* reference ASC forward calculation */
      if ( p7_ReferenceASCForward(sq->dsq, sq->n, gm, anch->a, anch->D, afu, afd, &sc1) != eslOK) esl_fatal(msg);
      
      /* sparse ASC forward calculation */
      if ( p7_sparse_asc_Forward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asf, &sc2)   != eslOK) esl_fatal(msg);
      


      p7_sparsemask_Reuse(sm);
      p7_anchors_Reuse(anch);
      p7_trace_Reuse(gtr);
      p7_refmx_Reuse(afu);
      p7_refmx_Reuse(afd);
      p7_sparsemx_Reuse(asf);
      esl_sq_Reuse(sq);
    }

  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(afd);
  p7_anchors_Destroy(anch);
  p7_sparsemask_Destroy(sm);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}


#endif /*p7SPARSE_ASC_FWDBACK_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/



/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7SPARSE_ASC_FWDBACK_TESTDRIVE


#endif /*p7SPARSE_ASC_FWDBACK_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/


/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7SPARSE_ASC_FWDBACK_EXAMPLE

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
static char usage[]  = "[-options] <hmmfile> <seqfile> <ndom> [<i0> <k0>]...";
static char banner[] = "example of ASC Forward production (sparse) implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, -1, argc, argv, banner, usage);
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
  P7_TRACE       *vtr     = NULL;
  P7_REFMX       *rxf     = NULL;
  P7_REFMX       *rxd     = NULL;
  P7_REFMX       *afu     = NULL;
  P7_REFMX       *afd     = NULL;
  P7_SPARSEMX    *asf     = NULL;
  P7_COORDS2     *anch    = p7_coords2_Create(0,0);
  float          *wrk     = NULL;
  P7_COORDS2_HASH *hashtbl = p7_coords2_hash_Create(0,0,0);
  float           fsc, vsc, asc, asc_sparse;
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
   * To get it, run the reference prototype code;
   * (we don't have the MPAS algorithm in its sparse production form yet)
   */
  vtr = p7_trace_Create();
  rxf = p7_refmx_Create(gm->M, sq->n);
  rxd = p7_refmx_Create(gm->M, sq->n);
  afu = p7_refmx_Create(gm->M, sq->n);
  afd = p7_refmx_Create(gm->M, sq->n);

  p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf, vtr, &vsc);
  p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc);   
  p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd,      NULL);   
  p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd);   

  p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, vtr, &wrk, hashtbl,
		       afu, afd, anch, &asc, NULL, NULL);


  //p7_refmx_Dump(stdout, afu);
  //p7_refmx_Dump(stdout, afd);

  /* Finally...
   * Run sparse ASC Forward.
   */
  asf = p7_sparsemx_Create(sm);
  p7_sparse_asc_Forward(sq->dsq, sq->n, gm, anch->arr, anch->n, sm, asf, &asc_sparse);

  

  //p7_spascmx_Dump(stdout, asf, anch->arr, anch->n);

  printf("Reference ASC fwd score = %.2f nats\n", asc);
  printf("Sparse ASC fwd score    = %.2f nats\n", asc_sparse);

  p7_coords2_hash_Destroy(hashtbl);
  if (wrk) free(wrk);
  p7_sparsemx_Destroy(asf);
  p7_coords2_Destroy(anch);
  p7_trace_Destroy(vtr);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(cx);
  p7_refmx_Destroy(afu);   p7_refmx_Destroy(afd);
  p7_refmx_Destroy(rxf);   p7_refmx_Destroy(rxd);
  p7_profile_Destroy(gm);  p7_oprofile_Destroy(om);
  esl_sqfile_Close(sqfp);  esl_sq_Destroy(sq);
  p7_hmmfile_Close(hfp);   p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPARSE_ASC_FWDBACK_EXAMPLE*/
/*---------------------- end, example ---------------------------*/


/*****************************************************************
 * 6. Notes
 *****************************************************************
 *
 * 
 *****************************************************************/






/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
