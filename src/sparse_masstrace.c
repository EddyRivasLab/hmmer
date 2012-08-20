/* Mass trace algorithm, sparse DP version:
 * finding bounds of a domain envelope.
 */

#include "p7_config.h"

#include "easel.h"

#include "hmmer.h"
#include "p7_sparsemx.h"
#include "sparse_envscore.h"

int
p7_sparse_masstrace_Up(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *fwd, P7_SPARSEMX *mass, P7_TRACE *tr, int z, float massthresh, int *ret_iae, int *ret_kae)
{
  const P7_SPARSEMASK *sm = fwd->sm;
  int   st0 = tr->st[z];	/* anchor point's state (p7T_{MDI}) */
  int   k0  = tr->k[z];		/* anchor point's k position in profile (1..M) */
  int   i0  = tr->i[z];		/* anchor point's i positon in sequence (1..L) */
  int   iae = i0;		/* RETURN: envelope start coord in sequence (1..i0) */
  int   kae = k0;		/* RETURN: envelope start coord in profile (1..k0) */
  float *rhoc;			/* ptr into current row i rho (trace mass) cells */
  float *rhon;			/* ptr into next row i+1 rho (trace mass) cells */
  const float *dpc;		/* ptr into current row i Forward DP cells */
  const float *dpn;		/* ptr into next row i+1 Forward DP cells */
  float *last_rhoc;
  const float *last_dpc;
  const float *rsn;		/* residue score vector for next row i+1. Enables MSN(k) notation macro. */
  const float *tsc = gm->tsc;   /* transition score vector. Enables TSC(k) notation macro. */
  int   i, k, v,w;
  float  rowmass;
  float *kmass           = NULL;
  int    iae_proven_done = FALSE;
  int    kae_proven_done = FALSE;
  int    status;

  /* we should avoid this tmp alloc. find space for kmass in sparsemx. */
  ESL_ALLOC(kmass, sizeof(float) * (sm->M+1));
  for (k = 0; k<= sm->M; k++) 
    kmass[k] = 0.;

  dpc  = fwd->dp  + p7S_NSCELLS*((sm->k[i0] + sm->n[i0] - 1) - sm->kmem);      // unwarranted pointer-arithmetic chumminess with sparsemask, sparsemx: assumes that dp[] and kmem[] are EXACTLY 1:1 with each other
  rhoc = mass->dp + p7S_NSCELLS*((sm->k[i0] + sm->n[i0] - 1) - sm->kmem);      // sm->k[i+1]-1 doesn't work because of i=L edge case (no k[L+1]) 

  /* ok, so, the following is insane.
   * note that local, glocal paths don't cross; 
   * if st0 is {MDI}L, only {MDI}L can be reached; analogous for {MDI}G.
   * so we only need to compute the L or G interleaved half of the sparse matrix, depending on st0
   * so...
   * the supercells are [ ML MG IL IG DL DG ]
   * dpc, rhoc currently point at the final supercell on row i: specifically at its ML cell
   * so... 
   * if st0 is {MDI}G glocal, offset ptrs by +1 cell; now our ptrs are on MG 
   * because of relative indexing in dp (everything done by stepping exactly in supercell units),
   * this single offset suffices to keep dpc,dpn,rhoc,rhon on L vs. G in all the rest of the code here
   * we access everything with L indices, but with the offset, these will exactly be the G cells
   * seriously, trust me.
   */
  if (st0 == p7T_MG || st0 == p7T_IG || st0 == p7T_DG)
    { dpc += 1; rhoc += 1; }
  last_dpc  = dpc;
  last_rhoc = rhoc;

  /* special case the first row (i0): initialize on i0,k0, then pull delete path to the left.  */
  /* first, skip to find k0, bravely assuming that k0 MUST be in the sparse list for i0 */
  for (v = sm->n[i0]-1; sm->k[i0][v] != k0; v--)
    { dpc  -= p7S_NSCELLS; rhoc -= p7S_NSCELLS;  }

  /* now v is the index of k0 in row i0's sparse cell list;
   *     *dpc is the sparse DP supercell \alpha(i0,k0) for X = {ML MG IL IG DL DG}
   *     *rho is \rho(i0,k0) supercell
   */
  switch (st0) {
  case p7T_ML: case p7T_MG: rhoc[p7S_ML] = 1.; rhoc[p7S_IL] = 0.;  rhoc[p7S_DL] = 0.;  break;
  case p7T_IL: case p7T_IG: rhoc[p7S_ML] = 0.; rhoc[p7S_IL] = 1.;  rhoc[p7S_DL] = 0.;  break;
  case p7T_DL: case p7T_DG: rhoc[p7S_ML] = 0.; rhoc[p7S_IL] = 0.;  rhoc[p7S_DL] = 1.;  break;
  default:     ESL_EXCEPTION(eslEINCONCEIVABLE, "you know it is");
  }
  kmass[k0]  = 1.0;
  dpc       -= p7S_NSCELLS;
  rhoc      -= p7S_NSCELLS;
  
  /* now pull to the left. If we didn't start on a D, or
   * if we don't have contiguous supercells on the row, this
   * is all an expensive way of zeroing the row: but it's
   * clearer this way than a lot of special case branching,
   * especially since it's obvious where the same ->D case is
   * in the main recursion later.
   */ 
  for (v = v-1; v >= 0; v--) 
    {
      k = sm->k[i0][v];
      if (sm->k[i0][v+1] == k+1) {
	rhoc[p7S_ML] = rhoc[p7S_DL+p7S_NSCELLS] * exp( dpc[p7S_ML] + TSC(p7P_MD, k) - dpc[p7S_DL+p7S_NSCELLS]);
	rhoc[p7S_IL] = 0.0f;
	rhoc[p7S_DL] = rhoc[p7S_DL+p7S_NSCELLS] * exp( dpc[p7S_DL] + TSC(p7P_DD, k) - dpc[p7S_DL+p7S_NSCELLS]);
      } else { rhoc[p7S_ML] = rhoc[p7S_IL] = rhoc[p7S_DL] = 0.0f; }

      kmass[k] += rhoc[p7S_ML] + rhoc[p7S_DL];
      if (kmass[k] >= massthresh) kae = k; 
      dpc      -= p7S_NSCELLS;
      rhoc     -= p7S_NSCELLS;
    }

  
  /* The main recursion.  */
  for (i = i0-1; i >= 1; i--)
    {
      dpn     = last_dpc;
      rhon    = last_rhoc;
      rsn     = gm->rsc[dsq[i+1]];  // MSN() notation macro now valid
      rowmass = 0.;

      last_dpc  = dpc;
      last_rhoc = rhoc;

      w = sm->n[i+1]-1;
      v = sm->n[i] - 1;
      while (v >= 0 && sm->k[i][v] > k0) { v--; dpc -= p7S_NSCELLS; rhoc -= p7S_NSCELLS; }
      if    (v < 0) break;	/* no cells on row at all? trace mass can't flow back any further then; we're done for sure. */
      for (; v >= 0; v--)	/* for all sparse k on row, such that k <= k0. if no such k exist, this code doesn't execute, and mass flow is done */
	{
	  k = sm->k[i][v];

	  /* Try to find the M(i+1,k+1) cell on row i+1. If it exists: apportion its mass to our current cells {MID}ik */
	  while (w >= 0 && sm->k[i+1][w]  > k+1) { w--; dpn -= p7S_NSCELLS; rhon -= p7S_NSCELLS; }
	  if    (w >= 0 && sm->k[i+1][w] == k+1 && k < k0) {  // note k<k0 test; for k=k0, Mk+1 was never initialized in rhon[]
	    rhoc[p7S_ML]  = rhon[p7S_ML] * exp( dpc[p7S_ML] + TSC(p7P_MM, k) + MSN(k+1) - dpn[p7S_ML]);
	    rhoc[p7S_IL]  = rhon[p7S_ML] * exp( dpc[p7S_IL] + TSC(p7P_IM, k) + MSN(k+1) - dpn[p7S_ML]);
	    rhoc[p7S_DL]  = rhon[p7S_ML] * exp( dpc[p7S_DL] + TSC(p7P_DM, k) + MSN(k+1) - dpn[p7S_ML]);
	  } else { rhoc[p7S_ML] = rhoc[p7S_IL] = rhoc[p7S_DL] = 0.0f; }

	  /* Try to find I(i+1,k) cell on row i+1; if exists, apportion its mass */
	  while (w >= 0 && sm->k[i+1][w]  > k) { w--; dpn -= p7S_NSCELLS; rhon -= p7S_NSCELLS; }
	  if    (w >= 0 && sm->k[i+1][w] == k) { 
	    rhoc[p7S_ML] += rhon[p7S_IL] * exp( dpc[p7S_ML] + TSC(p7P_MI, k) - dpn[p7S_IL]); // insert scores ISN(k) assumed to be zero
	    rhoc[p7S_IL] += rhon[p7S_IL] * exp( dpc[p7S_IL] + TSC(p7P_II, k) - dpn[p7S_IL]);
	  }

	  /* If v+1 is k+1 ... then v,v+1 are contiguous supercells ... and D(i,k+1) exists, so apportion its mass. */
	  if (v < sm->n[i]-1 && sm->k[i][v+1] == k+1 && k < k0) { // k<k0 test: don't look for Dk0+1
	    rhoc[p7S_ML] += rhoc[p7S_DL+p7S_NSCELLS] * exp( dpc[p7S_ML] + TSC(p7P_MD, k) - dpc[p7S_DL+p7S_NSCELLS]);
	    rhoc[p7S_DL] += rhoc[p7S_DL+p7S_NSCELLS] * exp( dpc[p7S_DL] + TSC(p7P_DD, k) - dpc[p7S_DL+p7S_NSCELLS]);
	  }

	  kmass[k] += rhoc[p7S_ML] + rhoc[p7S_DL]; /* kmass[k] is a lower bound on how much probability mass is flowing leftwards thru this column  */
	  if (k < kae && kmass[k] >= massthresh) kae = k; 
	  if (kae == 1 || kmass[k] + rowmass < massthresh) kae_proven_done = TRUE; /* kmass[k] + rowmass is the upper bound on what can flow leftward thru k */
	  rowmass  += rhoc[p7S_ML] + rhoc[p7S_IL]; /* how much total probability mass is flowing upwards through this row  */

	  rhoc -= p7S_NSCELLS;
	  dpc  -= p7S_NSCELLS;
	}

      if (rowmass < massthresh)  iae_proven_done = TRUE; else iae = i;
      if (iae_proven_done && kae_proven_done) break;
    }
  
  free(kmass);
  mass->type = p7S_MASSTRACE;
  *ret_iae   = iae;
  *ret_kae   = kae;
  return eslOK;

 ERROR:
  if (kmass) free(kmass);
  return status;
}


int
p7_sparse_masstrace_Down(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *fwd, P7_SPARSEMX *mass, P7_TRACE *tr, int z, float massthresh, int *ret_ibe, int *ret_kbe)
{
  const P7_SPARSEMASK *sm = fwd->sm; 
  int          st0 = tr->st[z];	/* anchor point's state (p7T_{MDI}) */
  int          k0  = tr->k[z];	/* anchor point's k position in profile (1..M) */
  int          i0  = tr->i[z];	/* anchor point's i positon in sequence (1..L) */
  int          ibe = i0;	/* RETURN: envelope end coord in sequence (i0..L) */
  int          kbe = k0;	/* RETURN: envelope end coord in profile (10..M) */
  float       *rhoc;  		/* ptr that steps through sparse mass supercells on current row i */
  float       *rhop;		/* ptr that steps through sparse mass supercells on previous row i-1 */
  float       *last_rhoc;	/* used to set the next rhop */
  const float *dpc;		/* ptr that steps through sparse Forward supercells on current row i */
  const float *dpp;		/* ptr that steps through sparse Forward supercells on previous row i-1 */
  const float *last_dpc;	/* used to set the next dpp */
  const float *rsc;		/* residue score vector on current row. Enables MSC(k) notation macro */
  const float *tsc   = gm->tsc; /* transition score vector. Enables TSC(k) notation macro */
  float        rowmass;
  float       *kmass = NULL;
  int ibe_proven_done = FALSE;
  int kbe_proven_done = FALSE;
  int v;			/* v index steps through list of sparse supercells k on current row i    */
  int w;			/* w index steps through list of sparse supercells k on previous row i-1 */
  int k;			/* index in profile positions 1..M */
  int i;			/* index in sequence positions 1..L */
  int status;

  /* we should avoid this tmp alloc. find space for kmass in sparsemx. */
  ESL_ALLOC(kmass, sizeof(float) * (sm->M+1));
  for (k = 0; k<= sm->M; k++) 
    kmass[k] = 0.;

  dpc  = fwd->dp  + p7S_NSCELLS*(sm->k[i0] - sm->kmem);  // unwarranted ptr-arithmetic chumminess with sparse matrix layout
  rhoc = mass->dp + p7S_NSCELLS*(sm->k[i0] - sm->kmem);
  
  /* See comment in _Down() about the following magic: */
  if (st0 == p7T_MG || st0 == p7T_IG || st0 == p7T_DG)
    { dpc += 1; rhoc += 1; }

  last_dpc  = dpc;  // these two lines must be placed after the magic above, of course
  last_rhoc = rhoc;

  /* Start initializing the first row, i0. 
   * Up() is responsible for k<k0. We are (Down() is) responsible for k>k0.
   * In both Up() and Down(), (s0,i0,k0) anchor itself is set to 1.0, and
   * it doesn't hurt to redo that here (just in case Up() isn't called too).
   */
  /* first, skip to k0. We trust that it exists in sparse cell list: caller said so. */
  for (v = 0; sm->k[i0][v] != k0; v++) { dpc += p7S_NSCELLS; rhoc += p7S_NSCELLS; }

  /* now v is the index of k0 in row i0's sparse cell list, and rhoc/dpc are on the
   * anchor cell
   */
  switch (st0) {
  case p7T_ML: case p7T_MG: rhoc[p7S_ML] = 1.; rhoc[p7S_IL] = 0.;  rhoc[p7S_DL] = 0.;  break;
  case p7T_IL: case p7T_IG: rhoc[p7S_ML] = 0.; rhoc[p7S_IL] = 1.;  rhoc[p7S_DL] = 0.;  break;
  case p7T_DL: case p7T_DG: rhoc[p7S_ML] = 0.; rhoc[p7S_IL] = 0.;  rhoc[p7S_DL] = 1.;  break;
  default:     ESL_EXCEPTION(eslEINCONCEIVABLE, "you know it is");
  }
  kmass[k0]  = 1.0;
  dpc       += p7S_NSCELLS;
  rhoc      += p7S_NSCELLS;

  /* now pull to the right, on the delete path.
   * if we started on an I, this is just an expensive way to zero remaining sparse cells k>k0 
   */
  for (v = v+1; v < sm->n[i0]; v++) 
    {
      k = sm->k[i0][v];
      rhoc[p7S_ML] = rhoc[p7S_IL] = 0.;
      if (sm->k[i0][v-1] == k-1) 
	{
	  rhoc[p7S_DL] = 
	    rhoc[p7S_ML-p7S_NSCELLS] * exp ( dpc[p7S_DL] + TSC(p7P_MD, k-1) - dpc[p7S_ML-p7S_NSCELLS]) +  // yes, those array indices are negative;
	    rhoc[p7S_DL-p7S_NSCELLS] * exp ( dpc[p7S_DL] + TSC(p7P_DD, k-1) - dpc[p7S_DL-p7S_NSCELLS]);   // but that's valid C,  seriously.
	}
      else rhoc[p7S_DL] = 0.;

      kmass[k] += rhoc[p7S_ML] + rhoc[p7S_DL];
      if (kmass[k] >= massthresh) kbe = k; 
      dpc  += p7S_NSCELLS;
      rhoc += p7S_NSCELLS;
    }

  /* The main recursion */
  for (i = i0 + 1; i <= L; i++)
    {
      dpp     = last_dpc;
      rhop    = last_rhoc;
      rsc     = gm->rsc[dsq[i]];  // MSC(k) notation macro now valid
      rowmass = 0.;
      
      last_dpc  = dpc;
      last_rhoc = rhoc;

      v = w = 0;
      while (v <  sm->n[i] && sm->k[i][v] < k0) { v++; dpc += p7S_NSCELLS; rhoc += p7S_NSCELLS; }  // for Down pass, k >= k0
      if    (v == sm->n[i]) break;  // no cells on row i at all? trace mass can't flow any more. break completely out (of i loop), we're done.
      for (; v <  sm->n[i]; v++)
	{
	  k = sm->k[i][v];

	  /* Try to find supercell {MID}(i-1,k-1). If it exists, each cell apportions some of its mass to M(i,k) */
	  while (w < sm->n[i-1] && sm->k[i-1][w]  < k-1) { w++; dpp += p7S_NSCELLS; rhop += p7S_NSCELLS; }
	  if    (w < sm->n[i-1] && sm->k[i-1][w] == k-1 && k > k0) // k > k0 test is there to prevent looking at k0-1 supercell
	    {
	      rhoc[p7S_ML] = 
		rhop[p7S_ML] * exp( dpc[p7S_ML] + TSC(p7P_MM, k-1) + MSC(k) - dpp[p7S_ML]) +
		rhop[p7S_IL] * exp( dpc[p7S_ML] + TSC(p7P_IM, k-1) + MSC(k) - dpp[p7S_IL]) +
		rhop[p7S_DL] * exp( dpc[p7S_ML] + TSC(p7P_DM, k-1) + MSC(k) - dpp[p7S_DL]);
	    }
	  else rhoc[p7S_ML] = 0.;

	  /* Try to find supercell {MID}(i-1,k), which is either w or w+1 if it exists; if so, its MI cells apportion mass to I(i,k) */
	  while (w < sm->n[i-1] && sm->k[i-1][w]  < k) { w++; dpp += p7S_NSCELLS; rhop += p7S_NSCELLS; }
	  if    (w < sm->n[i-1] && sm->k[i-1][w] == k && k < gm->M)  // k=M check because Im doesn't exist, and we need to avoid a -inf - -inf = nan
	    {
	      rhoc[p7S_IL] = 
		rhop[p7S_ML] * exp( dpc[p7S_IL] + TSC(p7P_MI, k) - dpp[p7S_ML]) +   // here we're assuming ISC(k)=0
		rhop[p7S_IL] * exp( dpc[p7S_IL] + TSC(p7P_II, k) - dpp[p7S_IL]);    // ditto

	    }
	  else rhoc[p7S_IL] = 0.;

	  /* If v-1 is k-1, then v-1,v are contiguous, and supercell i,k-1) exists; if so, its MD cells appportion mass to D(i,k)  */
	  if (v > 0 && sm->k[i][v-1] == k-1 && k > k0) 
	    {
	      rhoc[p7S_DL] = 
		rhoc[p7S_ML-p7S_NSCELLS] * exp ( dpc[p7S_DL] + TSC(p7P_MD, k-1) - dpc[p7S_ML-p7S_NSCELLS]) +  // yes, those array indices are negative;
		rhoc[p7S_DL-p7S_NSCELLS] * exp ( dpc[p7S_DL] + TSC(p7P_DD, k-1) - dpc[p7S_DL-p7S_NSCELLS]);   // but that's valid C,  seriously.
	    }
	  else rhoc[p7S_DL] = 0.;

	  kmass[k] += rhoc[p7S_ML] + rhoc[p7S_DL];       // lower bound on how much mass is flowing right, through column k (partial sum, rows i0..i). Don't count I.
	  if (k > kbe && kmass[k] > massthresh) kbe = k; // update kbe envelope end bound: rightmost k that satisfies threshold
	  if (kbe == gm->M || kmass[k] + rowmass < massthresh) kbe_proven_done = TRUE;  // *upper* bound on mass flowing right, by adding total mass that's still to the left of k 
	  rowmass  += rhoc[p7S_ML] + rhoc[p7S_IL]; // how much mass is still flowing down, through this row i. Don't count D.

	  rhoc += p7S_NSCELLS;  // advance rhoc, dpc ptrs by one supercell;
	  dpc  += p7S_NSCELLS;  // we will figure out its k index when the v loop rolls around now...
	}

      if (rowmass < massthresh) ibe_proven_done = TRUE; else ibe = i;
      if (ibe_proven_done && kbe_proven_done) break;
    } // end loop over i=i0..L

  free(kmass);
  mass->type = p7S_MASSTRACE;
  *ret_ibe   = ibe;
  *ret_kbe   = kbe;
  return eslOK;

 ERROR:
  if (kmass) free(kmass);
  return status;
}


/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7SPARSE_MASSTRACE_EXAMPLE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_sparsemx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",              0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of Forward/Backward, sparse dual implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_HMM         *hmm     = NULL;
  ESL_SQ         *sq      = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  P7_SPARSEMX    *sxv     = NULL;
  P7_SPARSEMX    *sxf     = NULL;
  P7_SPARSEMX    *sxb     = NULL;
  P7_SPARSEMX    *sxd     = NULL;
  P7_SPARSEMX    *sxm     = NULL;
  P7_TRACE       *tr      = p7_trace_CreateWithPP();
  float           vsc, fsc;
  int             iae,ibe;
  int             kae,kbe;
  float           envsc, envsc2;
  int             d;
  int             status;

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
 
  /* Read in one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);
  esl_sqfile_Close(sqfp);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);

  /* Allocate bands, matrices */
  sm  = p7_sparsemask_Create(gm->M, sq->n);
  p7_sparsemask_AddAll(sm);
  sxv = p7_sparsemx_Create(sm);
  sxf = p7_sparsemx_Create(sm);
  sxb = p7_sparsemx_Create(sm);
  sxd = p7_sparsemx_Create(sm);
  sxm = p7_sparsemx_Create(sm);

  /* Set the profile and null model's target length models */
  p7_bg_SetLength           (bg, sq->n);
  p7_profile_SetLength      (gm, sq->n);

  /* Sparse DP calculations */
  p7_SparseViterbi (sq->dsq, sq->n, gm, sxv, tr, &vsc);
  p7_SparseForward (sq->dsq, sq->n, gm, sxf,     &fsc);
  p7_SparseBackward(sq->dsq, sq->n, gm, sxb,     NULL);
  p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxb, sxd);
  p7_sparsemx_TracePostprobs(sxd, tr);
  p7_trace_Index(tr);

  //p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);
  
  for (d = 0; d < tr->ndom; d++)
    {
      p7_sparsemx_Reinit(sxm, sm);

      p7_sparse_masstrace_Up  (sq->dsq, sq->n, gm, sxf, sxm, tr, tr->anch[d], 0.1, &iae, &kae);
      p7_sparse_masstrace_Down(sq->dsq, sq->n, gm, sxb, sxm, tr, tr->anch[d], 0.1, &ibe, &kbe);
      p7_sparsemx_Reuse(sxm);
      
      p7_sparsemx_ApproxEnvScore(gm, sxf, iae, ibe, &envsc);

      p7_sparsemx_Reinit(sxm, sm);
      p7_SparseEnvScore(sq->dsq, sq->n, gm, iae, ibe, kae, kbe, sxm, &envsc2);

      printf("# domain %3d  iali: %d..%d [%daa]  ienv: %d..%d [%daa]  kali: %d..%d [%daa]  kenv: %d..%d [%daa]  envsc: %.2f  envsc2: %.2f\n",
	     d,
	     tr->sqfrom[d],  tr->sqto[d],  tr->sqto[d]-tr->sqfrom[d]+1, iae, ibe, ibe-iae+1,
	     tr->hmmfrom[d], tr->hmmto[d], tr->hmmto[d]-tr->hmmfrom[d]+1, kae, kbe, kbe-kae+1,
	     envsc, envsc2);

      p7_sparsemx_Reuse(sxm);
    }
  
  /* Cleanup */
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_sparsemx_Destroy(sxv);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemx_Destroy(sxm);
  p7_sparsemask_Destroy(sm);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
/*------------------ end, example driver ------------------------*/
#endif /*p7SPARSE_MASSTRACE_EXAMPLE*/



/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

/* 
 *  The probability mass that flows through state X at position i,k is calculated recursively:
 *  
 *  \rho^X(i,k) = \sum_Y \rho^Y(i_Y,k_Y) * 
 *                          exp \left(  \alpha^X(i,k) + \tau_k(XY) + \epsilon(Y_k_Y, x_i_Y) - \alpha^Y(i_Y, k_Y)} \right)
 *  where the notation is as follows:
 *    X = current state
 *    Y = state that X can transition to (MID only; exclude E)
 *    i_Y,k_Y    = condensed notation for next i,k. i_M,k_M = i+1,k+1; i_I,k_I = i+1,k; i_D,k_D = i,k+1.
 *    \alpha()   = Forward matrix, log space in nats
 *    \tau_k(XY) = log transition prob, for example log t(Mk->Mk+1) is \tau_k(MM)
 *    \epsilon() = emission scores, for example e(Mk+1, x_i+1)
 *    
 *  A little useful guidance:
 *    
 *    this thing: \alpha^X(i,k) + \tau_k(XY) + \epsilon(Y_k_Y, x_i_Y)
 *    is the XY component (edge) of the Forward recursion into state Y, accounting for the X->Y state transition.
 *    
 *    this thing: \alpha^Y(i_Y, k_Y)} 
 *    is the sum of all such edges into Y: i.e., the Forward recursion.
 *    
 *    if these were in probability space, the ratio of one edge into Y
 *    over all edges into Y gives us the fraction of posterior
 *    probability mass that flows back that edge: call this the _edge
 *    weight_. The sum of all edge weights into Y is 1. Since
 *    everything's in null-scaled log space, the ratio shows up as a
 *    subtraction and an exp().
 *    
 *    For a Y=Mk match state, one such incoming edge weight is the
 *    B->Mk edge; mass that flows back that edge disappears from our
 *    recursion.  Essentially, it's this lossage that we're tracking
 *    back: how far do we have to trace back in i,k before enough mass
 *    is lost from the main state alignment ensemble that the mass
 *    we're still tracing back drops below the <massthresh> threshold?
 *    
 *    In a stochastic traceback, we would trace back from Y
 *    stochastically, according to its edge weight distribution
 *    (inclusive of the B->Mk edges). Notationally, we implement a
 *    stochastic traceback as a _push_ from a current state Y to a
 *    next state X, whereas the mass traceback is implemented as a
 *    _pull_ to a current state X from connected states Y.
 */


