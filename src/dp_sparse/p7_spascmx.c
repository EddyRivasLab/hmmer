/* The p7_spascmx module does not provide an independent data
 * structure, but rather a set of routines for using P7_SPARSEMX and
 * P7_SPARSEMASK for sparse anchor set contrained DP calculations.
 * 
 * To be included in a sparse ASC matrix, a main supercell i,k must
 * satisfy all of the following three conditions:
 * 
 * 1. It's in the sparse mask. 
 * 2. It's in an ASC UP or DOWN sector.
 * 3. It is "connected to" its anchor.
 * 
 * In more detail:
 *    
 * 1. Row i is in the sparse mask segment g if  seg[g].ia <= i <= seg[g].ib;
 *    Cell (i,k) is in the mask if for some z,  sm->k[i][z] == k
 *    
 * 2. Suppose d is the index of the anchor on or after i:
 *       d = argmin_d anch[d].i0 >= i
 *    Then (i,k) is in UP(d) if k < anch[d].k0;
 *    and (i,k) is in DOWN(d-1) if k >= anch[d-1].k0
 *    (It is possible to be in both UP and DOWN.)
 *    
 * 3. In principle, we only need to calculate (i,k,s) cell if
 *    there is a valid DP path that connects it to "its anchor":
 *    where "its anchor" means anch[d].i0,k0 for a cell in UP(d)
 *    (the anchor down and right from i,k), or anch[d-1].i0,k0 
 *    for a cell in DOWN(d-1) (anchor up/left from i,k). However,
 *    calculating whether an individual i,k,s cell is connected
 *    to its anchor appears to require a DP calculation all of its
 *    own, so we don't do that.
 *    
 *    Instead, we use a less stringent criterion for "connected to its
 *    anchor", based only on whether rows i,i0 are connected, rather
 *    than a detailed path through cells. For cell (i,k),
 *    row i and its anchor i0 must be in the same segment g:
 *          seg[g].ia <= (i, i0) <= seg[g].ib
 *          
 * Another way to think about why we have the 3rd criterion for
 * connectedness, instead of just taking the intersection of the
 * sparse mask and the ASC UP/DOWN sectors: Consider a segment that
 * contains 0 anchors. Here we don't have to do anything (not even
 * storing specials), because no path through any cells in this
 * segment can pass thru an anchor. Consider a segment that contains 1
 * anchor. Now we'll just have an UP and a DOWN sector, with no rows
 * that have both UP and DOWN cells. The only time we can have a row
 * with both UP and DOWN cells in it is when there's 2 or more anchors
 * in the same sparse segment.
 * 
 *          
 * For a simple example of traversing a sparse ASC matrix, see
 * p7_spascmx_MinSizeof().
 * 
 * 
 * Contents:
 *    1. Using P7_SPARSEMX for ASC calculations.
 * 
 */

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"

#include "base/p7_coords2.h"
#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/p7_spascmx.h"


/*****************************************************************
 * 1. Using P7_SPARSEMX for ASC calculations
 ***************************************************************** 


/* Function:  p7_spascmx_Reinit()
 * Synopsis:  Reinitialize, reallocate sparse ASC matrix for new DP problem.
 *
 * Purpose:   Reinitialize and, if necessary, reallocate an existing sparse
 *            matrix <asx> for a new ASC DP calculation that will be
 *            constrained both by the sparse mask <sm> and the 
 *            anchor array <anch>.
 *            
 *            <asx> keeps an internal pointer to <sm>, so the caller
 *            must not modify <sm> while <asx> remains in use.
 *
 * Args:      asx   : sparse ASC matrix to reinitialize
 *            sm    : sparse mask that will constrain the new ASC DP calculation
 *            anch  : anchor array (1..D, with sentinels) that constrains new ASC DP calc
 *            D     : number of domains defined by anchors in <anch>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error. Now <asx> is in an undefined state,
 *            and can only be <_Destroy>'ed.
 */
int
p7_spascmx_Reinit(P7_SPARSEMX *asx, const P7_SPARSEMASK *sm, const P7_ANCHOR *anch, int D)
{
  int64_t dalloc_req;    // denominated in i,k supercells, each w/ p7_NSCELLS; number of i,k main supercells stored by a sparse ASC DP calculation
  int     xalloc_req;    // denominated in i supercells, each w/ p7_NXCELLS:   number of rows that have specials stored
  int     status;

  p7_spascmx_MinSizeof(sm, anch, D, &dalloc_req, &xalloc_req);

  if (dalloc_req > asx->dalloc) {
    ESL_REALLOC(asx->dp, sizeof(float) * p7S_NSCELLS * dalloc_req);
    asx->dalloc = dalloc_req;
  }
  if (xalloc_req > asx->xalloc) {
    ESL_REALLOC(asx->xmx, sizeof(float) * p7S_NXCELLS * xalloc_req);
    asx->xalloc = xalloc_req;
  }
  asx->sm   = sm;
  asx->type = p7S_UNSET;   // DP routines themselves set this.
  return eslOK;
  
 ERROR: 
  return status;
}


/* Function:  p7_spascmx_MinSizeof()
 * Synopsis:  Calculate minimum allocation size needed for sparse ASC DP calculation.
 *
 * Purpose:   For a sparse ASC DP calculation constrained by sparse mask
 *            <sm> and anchor array <anch> for <1..D> domains,
 *            calculate minimum required allocation size. Return the
 *            total size in bytes. You might use this for collecting
 *            statistics on how much memory is required by sparse ASC
 *            DP calculations.
 *            
 *            Optionally, return in <opt_dalloc> the number of (i,k)
 *            main supercells that need to be stored, and in
 *            <opt_xalloc> the number of i rows for which specials
 *            need to be stored. <p7_spascmx_Reinit()> uses these
 *            numbers when it reallocates a matrix for a new DP problem.
 *            
 *            This routine also makes a good example of how to
 *            traverse a sparse ASC matrix.
 *
 * Args:      sm         : sparse mask that constrains the DP calc
 *            anch       : anchor set array 1..D (with sentinels 0,D+1) 
 *            D          : number of domains/anchors in <anch>
 *            opt_dalloc : optRETURN: number of main i,k supercells that need to be stored
 *            opt_xalloc : optRETURN: number of rows for which specials need to be stored
 *
 * Returns:   Minimum allocation size required for the complete <P7_SPARSEMX>
 *            structure, in bytes.
 *
 * Throws:    (no abnormal error conditions)
 */
size_t
p7_spascmx_MinSizeof(const P7_SPARSEMASK *sm, const P7_COORD2 *anch, int D, int64_t *opt_dalloc, int *opt_xalloc)
{
  size_t  n       = sizeof(P7_SPARSEMX);
  int     g       = 1;		// index of next or current segment. When we enter it, and while we're in it, in_seg = TRUE.
  int     in_seg  = FALSE;      //   ... this bumps to TRUE when we see ia(g), first row of segment; to FALSE on ib(g), last row.
  int     d       = 1;          // index of next anchor we will reach; thus a current cell may be in UP(d) or DOWN(d-1) sector.
  int     ndown   = 0;          //   ... this bumps to 1 when i reaches an anchor, then counts rows in the DOWN sector, then goes to 0 when we leave seg g. Using counter allows knowing when we're on top DOWN row.
  int64_t dalloc  = 0;          // number of supercells in main DP matrix
  int     xalloc  = 0;          // number of special supercells 
  int     i,z;

  for (i = 0; i <= sm->L; i++)
    {
      if      (i == anch[d].i0)     { ndown = 1;  d++;     }    // when i reaches next anchor; bump d to next domain index, and DOWN sector is active...
      else if (sm->n[i] == 0)       { ndown = 0;           }    //  ... until when we reach end of segment, when DOWN becomes inactive again.
      else                          { ndown++;             }    // counting ndown lets us easily test if we're on the special top row. Not used in this routine; here for pedagogy, since this is a traversal example

      if      (i >  sm->seg[g].ib)  { in_seg = FALSE; g++; }    // g bumps, to start expecting to see start of segment <g> next.
      else if (i == sm->seg[g].ia)  { in_seg = TRUE;       }    // g might be S+1, but this is safe because of sentinel seg[S+1].ia=ib=L+2      

      if (in_down)                                              // if i is in a DOWN sector:
	{
	  for (z  = 0; z < sm->n[i]; z++)
	    if (sm->k[i][z] >= anch[d-1].k0) break;             // d-1 is safe here; you can't be in_down with d=0.
	  dalloc += (sm->n[i] - z);                             // z is now on first cell in DOWN row; remainder of line is valid
	}

      if ( (i >= sm->seg[g].ia-1 && anch[d].i0 <= sm->seg[g].ib) ||  // if we need to store specials for row i because of UP xB...
	   (in_down))                                                //   .. or because of DOWN xE ...
	xalloc++;
      
      if (in_seg && anch[d].i0 <= sm->seg[g].ib)               // if i is in an UP sector:
	{
	  for (z = 0; z < sm->n[i]; z++)
	    if (sm->k[i][z] >= anch[d].k0) break;   	       
	  dalloc += z;                                        // z is now +1 past the last sparse cell on the UP row
	}
    }

  n += dalloc * sizeof(float) * p7S_NSCELLS;
  n += xalloc * sizeof(float) * p7S_NXCELLS;
  if (opt_dalloc) *opt_dalloc  = dalloc;
  if (opt_xalloc) *opt_xalloc  = xalloc;
  return n;
}


/* useful side effect: pass -1 for k0, and it will print the whole line
 *
  */

static void
dump_up_header(FILE *fp, int k1, int k2)
{
  int  width     = 9;
  int k;

  fprintf(fp, "\n# UP component(s) of sparse matrix\n");
  fprintf(fp, "       ");
  for (k = k1; k <= k2;         k++) fprintf(fp, "%*d ", width, k);
  fprintf(fp, "\n");

  fprintf(fp, "       ");
  for (k = k1; k <= k2;        k++) fprintf(fp, "%*.*s ", width, width, "----------");
  fprintf(fp, "\n");
}

static void
dump_down_header(FILE *fp, int k1, int k2)
{
  int  width     = 9;
  int k,x;

  fprintf(fp, "\n# DOWN component(s) of sparse matrix\n");
  fprintf(fp, "       ");
  for (k = k1; k <= k2;         k++) fprintf(fp, "%*d ", width, k);
  for (x = 0;  x < p7S_NXCELLS; x++) fprintf(fp, "%*s ", width, p7_sparsemx_DecodeSpecial(x));
  fprintf(fp, "\n");

  fprintf(fp, "       ");
  for (k = k1; k <= k2;        k++) fprintf(fp, "%*.*s ", width, width, "----------");
  for (x = 0; x < p7S_NXCELLS; x++) fprintf(fp, "%*.*s ", width, width, "----------");
  fprintf(fp, "\n");
}


static void 
dump_up_row(FILE *fp, int i, const P7_SPARSEMASK *sm, const float *dpc, int k0, int k1, int k2, int s)
{
  int  width     = 9;
  int  precision = 4;
  int  k,z;


  fprintf(fp, "%3d %2s ", i, p7_sparsemx_DecodeState(s));
  for (z = 0, k = k1; k <= k2 && k < k0; k++) {   
    while (z < sm->n[i] && sm->k[i][z]  < k) z++;
    if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(fp, "%*.*f ", width, precision, *(dpc + z*p7R_NSCELLS + s));
    else                                     fprintf(fp, "%*s ",   width, "......");
  }
  for ( ; k <= k2; k++) {
    while (z < sm->n[i] && sm->k[i][z]  < k) z++;
    if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(fp, "%*s ", width, "...xx.");
    else                                     fprintf(fp, "%*s ", width, "......");
  }
  fprintf(fp, "\n");
}

static void
dump_down_row(FILE *fp, int i, const P7_SPARSEMASK *sm, const float *dpc, int k0, int k1, int k2, int s)
{
  int width     = 9;
  int precision = 4;
  int k,z0,z;

  fprintf(fp, "%3d %2s ", i, p7_sparsemx_DecodeState(s));
  for (z = 0, k = k1; k <= k2 && k < k0; k++) {
    while (z < sm->n[i] && sm->k[i][z]  < k) z++;
    if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(fp, "%*s ", width, "...xx.");
    else                                     fprintf(fp, "%*s ", width, "......");
  }
  z0 = z + 1;
  for ( ; k <= k2; k++) {
    while (z < sm->n[i] && sm->k[i][z]  < k) z++;
    if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(fp, "%*.*f ", width, precision, *(dpc + (z-z0)*p7R_NSCELLS + s));
    else                                     fprintf(fp, "%*s ",   width, "......");
  }
}


int
p7_spascmx_Dump(FILE *fp, const P7_SPARSEMX *asx, const P7_COORD2 *anch, int D)
{
  const P7_SPARSEMASK *sm  = asx->sm;
  const float         *dpc;
  const float         *xc;
  int   width     = 9;
  int   precision = 4;
  int   i1        = 0;
  int   i2        = sm->L;
  int   k1        = 0;
  int   k2        = sm->M;
  int   in_seg, in_up, in_down, in_x;
  int   i,g,d,z,s,k0;

  /* First pass: print the UP sectors, skip DOWN; do nothing to xc yet  */
  dump_up_header(fp, k1, k2);
  d        = 0;
  in_seg   = FALSE;
  in_up    = FALSE;
  in_down  = FALSE;
  g        = 1;
  dpc      = asx->dp;
  for (i = 0; i <= sm->L; i++)
    {
      if      (d < D && i == anch[d].n1) { in_down = TRUE;  d++;  }  
      else if (sm->n[i] == 0)            { in_down = FALSE;       }  

      if      (i >  sm->seg[g].ib)  { in_seg = FALSE; g++; }
      else if (i == sm->seg[g].ia)  { in_seg = TRUE; }      

      if (in_seg && d < D && anch[d].n1 <= sm->seg[g].ib)  { in_up = TRUE;  k0 = anch[d].n2; }
      else                                                 { in_up = FALSE; k0 = -1; }

      if (in_down)
	{
	  for (z  = 0; z < sm->n[i]; z++)
	    if (sm->k[i][z] >= anch[d-1].n2) break;  // d-1 safe because you can't be in DOWN sector with d=0.
	  dpc += (sm->n[i] - z) * p7S_NSCELLS;       // skip past DOWN rows in the matrix in this pass
	}

      if (i >= i1 && i <= i2) { 
	for (s = 0; s < p7S_NSCELLS; s++) dump_up_row(fp, i, sm, dpc, k0, k1, k2, s);
	fprintf(fp, "\n");
      }
	  
      if (in_up) {
	for (z = 0; z < sm->n[i]; z++)
	  if (sm->k[i][z] >= anch[d].n2) break;    
	dpc += z * p7S_NSCELLS;
      }
    }
  
  /* Second pass: print DOWN sectors and xc rows */
  dump_down_header(fp, k1, k2);
  d        = 0;
  in_seg   = FALSE;
  in_up    = FALSE;
  in_down  = FALSE;
  g        = 1;
  dpc      = asx->dp;
  xc       = asx->xmx;
  for (i = 0; i <= sm->L; i++)
    {
      if      (d < D && i == anch[d].n1) { in_down = TRUE;  d++;  }  
      else if (sm->n[i] == 0)            { in_down = FALSE;       }  

      if      (i >  sm->seg[g].ib)  { in_seg = FALSE; g++; }
      else if (i == sm->seg[g].ia)  { in_seg = TRUE; }      

      k0 = (in_down ? anch[d-1].n2 : sm->M+1);

      if ( (i >= sm->seg[g].ia-1 && d < D && anch[d].n1 <= sm->seg[g].ib) ||  in_down) in_x = TRUE;
      else in_x = FALSE;

      if (i >= i1 && i <= i2) {
	dump_down_row(fp, i, sm, dpc, k0, k1, k2, p7S_ML);
	if (in_x) for (s = 0; s < p7S_NXCELLS; s++) fprintf(fp, "%*.*f ", width, precision, xc[s]);
	else      for (s = 0; s < p7S_NXCELLS; s++) fprintf(fp, "%*s ", width, ".....");
	fprintf(fp, "\n");
	for (s = 1; s < p7S_NSCELLS; s++) {
	  dump_down_row(fp, i, sm, dpc, k0, k1, k2, s);
	  fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
      }

      if (in_down) {
	for (z  = 0; z < sm->n[i]; z++)
	  if (sm->k[i][z] >= anch[d-1].n2) break;  // d-1 safe because you can't be in DOWN sector with d=0.
	dpc += (sm->n[i] - z) * p7S_NSCELLS;       // skip past DOWN rows in the matrix in this pass
      }

      if (in_x) 
	xc += p7S_NXCELLS;

      if (in_seg && d < D && anch[d].n1 <= sm->seg[g].ib) {
	for (z = 0; z < sm->n[i]; z++)
	  if (sm->k[i][z] >= anch[d].n2) break;    
	dpc += z * p7S_NSCELLS;
      }
    }

  return eslOK;
}


int
p7_spascmx_CompareReference(const P7_SPARSEMX *sx, const P7_REFMX *rxu, const P7_REFMX *rxd, float tol)
{
  const P7_SPARSEMASK *sm = sx->sm;
  int i;                              // index over rows in DP matrices, 0.1..L
  int g      = 1;                     // idx of next or current segment, 1..S, with sentinels. When we enter it, & while we're in it, in_seg = TRUE
  int in_seg = FALSE;                 //  ... => TRUE when starting ia(g), 1st row of segment; => FALSE when we pass ib(g).
  int d      = 1;                     // idx of next domain anchor we will see. Row i may be in sector UP(d), DOWN(d-1). 
  int ndown  = 0;                     // row # of DOWN sector, 1..; becomes 1 when i reaches anchor.
  int z;                              // index over sparse cell list for a row

  for (i = 0; i <= sm->L; i++)
    {
      if      (i == anch[d].n1)    { ndown = 1; d++;      }
      else if (sm->n[i] == 0)      { ndown = 0;           }
      else if (ndown)              { ndown++;             }

      if      (i > sm->seg[g].ib)  { in_seg = FALSE; g++; }
      else if (i == sm->seg[g].ia) { in_seg = TRUE;       }

      if (in_down)
	{
	  for (z = 0; z < sm->n[i]; z++)              // Skip ahead in sparse cell list to first z in DOWN sector (anch[d-1].k0).
	    if (sm->k[i][z] >= anch[d-1].n2) break;

	}




}



/*****************************************************************
 * x. Statistics collection driver.
 *****************************************************************/
#ifdef p7SPASCMX_STATS
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type           default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",          eslARG_NONE,   FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",                   0 },
  { "-s",          eslARG_INT,      "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "statistics collection on sparse ASC DP matrices";

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
  P7_COORDS2     *anch    = p7_coords2_Create(0,0);
  P7_REFMX       *rxf     = NULL;
  P7_REFMX       *rxd     = NULL;
  P7_REFMX       *afu     = NULL;
  P7_REFMX       *afd     = NULL;
  P7_TRACE       *vtr     = NULL;
  P7_FILTERMX    *fx      = p7_filtermx_Create(100);
  P7_CHECKPTMX   *cx      = p7_checkptmx_Create(100, 100, ESL_MBYTES(p7_RAMLIMIT));
  P7_SPARSEMASK  *sm      = p7_sparsemask_Create(100, 100);
  float          *wrk     = NULL;
  P7_COORDS2_HASH *hashtbl = p7_coords2_hash_Create(0,0,0);
  float           fsc, vsc, asc;
  int             dalloc, xalloc, spascmxsize;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_oprofile_Convert(gm, om);

  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Allocate DP matrices and tracebacks */
  vtr = p7_trace_Create();
  rxf = p7_refmx_Create(gm->M, 100);
  rxd = p7_refmx_Create(gm->M, 100);
  afu = p7_refmx_Create(gm->M, 100);
  afd = p7_refmx_Create(gm->M, 100);

  /* For each target sequence... */
  while (( status = esl_sqio_Read(sqfp, sq)) == eslOK) 
    {
      /* Set the profile and null model's target length models */
      p7_bg_SetLength           (bg, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_oprofile_ReconfigLength(om, sq->n);

      if (( status = p7_pipeline_AccelerationFilter(sq->dsq, sq->n, om, bg, fx, cx, sm)) == eslOK)
	{
	  /* First pass analysis */
	  p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf, vtr, &vsc);
	  p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc);
	  p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd, NULL);   
	  p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd);   

	  /* Find most probable anchor set */
	  p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, vtr, &wrk, hashtbl,
			       afu, afd, anch, &asc, NULL, NULL);

      
	  printf("%-15s  %-30s  ", gm->name, sq->name);

	  /* Reference matrix size: cells, special rows, total in bytes */
	  printf("%10d %10d %10d ",
		 (int) (sq->n+1) * (gm->M+1),
		 (int) (sq->n+1),
		 (int) p7_refmx_MinSizeof(gm->M, sq->n));
		 
	  /* Sparse matrix size: cells, special rows, total in bytes */
	  printf("%10d %10d %10d ",
		 (int) sm->ncells,
		 (int) (sm->nrow + sm->S),
		 (int) p7_sparsemx_MinSizeof(sm));

	  /* Sparse ASC matrix size: UP cells, DOWN cells, total cells, special rows, total in bytes */
	  spascmxsize = p7_spascmx_MinSizeof(sm, anch->arr, anch->n, &dalloc, &xalloc);
	  printf("%10d %10d %10d\n",
		 (int) dalloc,
		 (int) xalloc,
		 (int) spascmxsize);

	  p7_trace_Reuse(vtr);
	  p7_refmx_Reuse(rxf);   p7_refmx_Reuse(rxd);
	  p7_refmx_Reuse(afu);   p7_refmx_Reuse(afd);
	  p7_coords2_hash_Reuse(hashtbl);
	  p7_coords2_Reuse(anch);
	}

      esl_sq_Reuse(sq);
      p7_filtermx_Reuse(fx);
      p7_checkptmx_Reuse(cx);
      p7_sparsemask_Reuse(sm);
    }
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslEOF)     p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  esl_sqfile_Close(sqfp);
  p7_filtermx_Destroy(fx);
  p7_checkptmx_Destroy(cx);
  p7_sparsemask_Destroy(sm);
  p7_coords2_hash_Destroy(hashtbl);
  if (wrk) free(wrk);
  p7_trace_Destroy(vtr);
  p7_refmx_Destroy(afd);  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(rxd);  p7_refmx_Destroy(rxf);
  p7_coords2_Destroy(anch);
  esl_sq_Destroy(sq);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPASCMX_STATS*/
