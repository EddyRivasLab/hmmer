#include "p7_config.h"

#include <stdlib.h>

#include "easel.h"

#include "base/p7_coords2.h"
#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/p7_spascmx.h"


int
p7_spascmx_Reinit(P7_SPARSEMX *asx, const P7_SPARSEMASK *sm, const P7_COORD2 *anch, int D)
{
  int64_t dalloc_req;
  int     xalloc_req;
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
  asx->type = p7S_UNSET;
  return eslOK;
  
 ERROR: 
  return status;
}



size_t
p7_spascmx_MinSizeof(const P7_SPARSEMASK *sm, const P7_COORD2 *anch, int D, int64_t *opt_dalloc, int *opt_xalloc)
{
  size_t  n       = sizeof(P7_SPARSEMX);
  int     g       = 1;		// index of next or current segment. When we enter it, and while we're in it, in_seg = TRUE.
  int     in_seg  = FALSE;      //   ... this bumps to TRUE when we see ia(g), first row of segment; to FALSE on ib(g), last row.
  int     d       = 0;          // how many domain anchors we've reached so far. Row i may be in sector UP(d), DOWN(d-1).
  int     in_down = FALSE;      //   ... this bumps to TRUE when i reaches an anchor; back to FALSE when we end a segment. 
  int64_t dalloc  = 0;          // number of supercells in main DP matrix
  int     xalloc  = 0;          // number of special supercells 
  int     i,z;
  
  /* A tad more exegesis...  You might think we're just going to take
   * the intersection of the sparse mask and the UP/DOWN sectors of
   * the ASC matrix. But there's an additional trick. We only need to
   * evaluate UP sector d's rows i when they're in the same segment as
   * the anchor(d); likewise, for DOWN sector d's rows i. 
   * 
   * For a segment that contains 0 anchors, we don't have to do
   * anything, not even storing specials. For a segment that contains
   * 1 anchor, we'll just have an UP and a DOWN sector, with no rows
   * that have both UP and DOWN segments. The only time we have
   * a row with both UP and DOWN cells in it is when there's 2 or
   * more anchors in the same sparse segment.
   */

  for (i = 0; i <= sm->L; i++)
    {
      if      (d < D && i == anch[d].n1) { in_down = TRUE;  d++; }   // when i reaches next anchor; bump d to next domain index, and DOWN sector is active...
      else if (sm->n[i] == 0)            { in_down = FALSE;      }   //  ... until when we reach end of segment, when DOWN becomes inactive again.

      if      (i >  sm->seg[g].ib)  { in_seg = FALSE; g++; }         // g bumps, to start expecting to see start of segment <g> next.
      else if (i == sm->seg[g].ia)  { in_seg = TRUE; }               // g might be S+1, but this is safe because of sentinel seg[S+1].ia=ib=L+2      

      if (in_down)                                                   // if i is in a DOWN sector:
	{
	  for (z  = 0; z < sm->n[i]; z++)
	    if (sm->k[i][z] >= anch[d-1].n2) break;                  // d-1 is safe here; you can't be in_down with d=0.
	  dalloc += (sm->n[i] - z);                                  // z is now on first cell in DOWN row; remainder of line is valid
	}

      if ( (i >= sm->seg[g].ia-1 && d < D && anch[d].n1 <= sm->seg[g].ib) ||  // if we need to store specials for row i because of UP xB...
	   (in_down))                                                //   .. or because of DOWN xE ...
	xalloc++;
      
      if (in_seg && d < D && anch[d].n1 <= sm->seg[g].ib)            // if i is in an UP sector:
	{
	  for (z = 0; z < sm->n[i]; z++)
	    if (sm->k[i][z] >= anch[d].n2) break;   	            
	  dalloc += z;                                               // z is now +1 past the last sparse cell on the UP row
	}
    }

  n      += dalloc * sizeof(float) * p7S_NSCELLS;
  n      += xalloc * sizeof(float) * p7S_NXCELLS;
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
