/* Created an animated gif figure, showing how sampled suboptimal
 * traces fill in to become the decoding matrix.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_regexp.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "sandbox/reference_mpl_fwd.h"

/* pixel in a ppm file = 0..255 R,G,B */
typedef struct {
  uint8_t r;
  uint8_t g;
  uint8_t b;
} RGBPIXEL;

/* data for a .ppm file representation of a 2D matrix */
typedef struct {
  int        cellsize;		/* in pixels: each cell is cellsize x cellsize        */
  int        nrowpix;		/* actual image dimensions in pixels: rows...         */
  int        ncolpix;		/*   ... and columns.                                 */
  RGBPIXEL **ppmdata;		/* ppmdata[0] is allocated as a vector.               */
} MXIMAGE;

/* 10-class RGB heatmap scheme:
 * colorbrewer2.org 9-class Reds + zero class as light blue
 */
RGBPIXEL heatmap_colors[10] = {
  { 198, 219, 239 }, 
  { 255, 245, 240 },
  { 254, 224, 210 },
  { 252, 187, 161 },
  { 252, 146, 114 },
  { 251, 106,  74 },
  { 239,  59,  44 },
  { 203,  24,  29 },
  { 165,  15,  21 },
  { 103,   0,  13 },
};

static ESL_OPTIONS options[] = {
  /* name           type           default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",          eslARG_NONE,    NULL,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",                    0 },
  { "-i",          eslARG_STRING,  NULL,  NULL, NULL,   NULL,  NULL, NULL, "when plotting, restrict plot to rows <i1>..<i2>",         0 },
  { "-k",          eslARG_STRING,  NULL,  NULL, NULL,   NULL,  NULL, NULL, "when plotting, restrict plot to columns <k1>..<k2>",      0 },
  { "-N",          eslARG_INT,      "10", NULL, NULL,   NULL,  NULL, NULL, "number of frames/samples to produce",                     0 },
  { "--seed",      eslARG_INT,      "0",  NULL, NULL,   NULL,  NULL, NULL, "set the random number seed to <n>",                       0 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile> <outpfx>";
static char banner[] = "animation of sampled traces becoming a decoding matrix";

/* mximage routines defined later */
static MXIMAGE *mximage_Create(int nrow, int ncol, int cellsize);
static void     mximage_Destroy(MXIMAGE *mxi);
static int      mximage_SetAll(MXIMAGE *mxi, RGBPIXEL color);
static int      mximage_SetCell(MXIMAGE *mxi, int i, int j, RGBPIXEL color);
static int      mximage_Write(FILE *ofp, MXIMAGE *mxi);
static int      mximage_plot_refmx(MXIMAGE *mxi, P7_REFMX *pp, int ia, int ib, int ka, int kb);
static int      mximage_plot_trace(MXIMAGE *mxi, P7_TRACE *tr, int ia, int ib, int ka, int kb, RGBPIXEL color);

/* these functions are duplicated from decoding-study */
static int parse_coord_string(const char *cstring, int *ret_start, int *ret_end);

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go         = p7_CreateDefaultApp(options, 3, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng        = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  char           *hmmfile    = esl_opt_GetArg(go, 1);
  char           *seqfile    = esl_opt_GetArg(go, 2);
  char           *outpfx     = esl_opt_GetArg(go, 3);
  int             N          = esl_opt_GetInteger(go, "-N");
  P7_HMMFILE     *hfp        = NULL;
  ESL_ALPHABET   *abc        = NULL;
  ESL_SQ         *sq         = NULL;
  int             format     = eslSQFILE_UNKNOWN;
  ESL_SQFILE     *sqfp       = NULL;
  P7_BG          *bg         = NULL;
  P7_HMM         *hmm        = NULL;
  P7_PROFILE     *gm         = NULL;           /* profile in H4's standard dual-mode local/glocal */
  P7_REFMX       *vit        = p7_refmx_Create(100, 100);
  P7_REFMX       *fwd        = p7_refmx_Create(100, 100);
  P7_REFMX       *bck        = p7_refmx_Create(100, 100);
  P7_REFMX       *pp         = p7_refmx_Create(100, 100);
  P7_REFMX       *mpl        = p7_refmx_Create(100, 100);
  P7_REFMX       *tct        = NULL; /* counts accumulated from sampled traces */
  P7_REFMX       *ppa        = NULL; /* normalized counts from sampled traces */
  P7_TRACE       *tr         = p7_trace_Create();
  P7_COORDS2     *dom        = p7_coords2_Create(16, 64);
  FILE           *outfp      = NULL;
  char           *outfile    = NULL;
  float          *wrk        = NULL; /* tmp space needed by stochastic trace */
  MXIMAGE        *mxi        = NULL;
  RGBPIXEL        blackpixel = { 0, 0, 0 };
  int             cellsz     = 2;
  float           fsc,   vsc,   bsc,   mplsc;
  int             s;
  float           nullsc;
  int             ia,ib,ka,kb;	/* optional bounds of a plot window */
  int             status;

  /* Determine coords of dump plot, if any */
  ia = ib = ka = kb = 0;
  if ( esl_opt_IsOn(go, "-i")) parse_coord_string(esl_opt_GetString(go, "-i"), &ia, &ib);
  if ( esl_opt_IsOn(go, "-k")) parse_coord_string(esl_opt_GetString(go, "-k"), &ka, &kb);

  /* Read in one HMM. Set alphabet to whatever the HMM's alphabet is. */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  /* Configure profile from HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config     (gm, hmm, bg);
  
  /* Open sequence file, read one sequence from it */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));
  else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);  
  esl_sqfile_Close(sqfp);

  /* Set length models */
  p7_bg_SetLength     (bg,  sq->n);
  p7_profile_SetLength(gm,  sq->n);

  p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

  p7_ReferenceViterbi (sq->dsq, sq->n, gm, vit, tr,  &vsc);
  p7_ReferenceForward (sq->dsq, sq->n, gm, fwd,      &fsc);   
  p7_ReferenceBackward(sq->dsq, sq->n, gm, bck,      &bsc);
  p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, bck, pp);
  p7_trace_Index(tr);
  p7_coords2_SetFromTrace(dom, tr);
  p7_ReferenceMPLForward(sq->dsq, sq->n, gm, dom->arr, dom->n, mpl, &mplsc);

  /* set up the movie image structure */
  /* The ppm image has the seq running horizontally, profile running vertically */
  mxi = mximage_Create((ka && kb) ? kb-ka+1 : gm->M, (ia && ib) ? ib-ia+1 : sq->n, cellsz);

  /* Frame 0 of the movie is just the Viterbi parse */
  /* you have to make filenames as .000001.ppm, not just .1.ppm, or unix cmdline won't keep them in numeric order */
  if (esl_sprintf(&outfile, "%s.%06d.ppm", outpfx, 0) != eslOK) p7_Fail("allocation failed for outfile");
  if ((outfp = fopen(outfile, "w"))                   == NULL)  p7_Fail("failed to open %s for writing", outfile);  

  mximage_SetAll(mxi, heatmap_colors[0]);
  mximage_plot_trace(mxi, tr, (ia ? ia : 1), (ib ? ib : sq->n), (ka ? ka : 1), (kb ? kb : gm->M), blackpixel);
  mximage_Write(outfp, mxi);

  fclose(outfp);
  free(outfile);

  p7_trace_Reuse(tr);
  p7_coords2_Reuse(dom);
  p7_refmx_Reuse(mpl);

  /* Frames 1..N of the movie are suboptimal traces, accumulating
   * toward a pp matrix
   */

  tct = p7_refmx_Create(gm->M, sq->n); /* counts matrix */
  ppa = p7_refmx_Create(gm->M, sq->n); /* normalized counts - pp approx by sampling */
  p7_refmx_SetType  (tct, gm->M, sq->n, p7R_DECODING);
  p7_refmx_SetType  (ppa, gm->M, sq->n, p7R_DECODING);
  p7_refmx_SetValues(tct, 0.0);
  p7_refmx_SetValues(ppa, 0.0);

  for (s = 1; s <= N; s++)
    {
      if (esl_sprintf(&outfile, "%s.%06d.ppm", outpfx, s) != eslOK) p7_Fail("allocation failed for outfile");
      if ((outfp = fopen(outfile, "w"))                  == NULL) p7_Fail("failed to open %s for writing", outfile);

      p7_reference_trace_Stochastic(rng, &wrk, gm, fwd, tr);	  

      p7_trace_Index(tr);
      p7_coords2_SetFromTrace(dom, tr);
      p7_trace_Score(tr, sq->dsq, gm, &vsc);
      p7_ReferenceMPLForward(sq->dsq, sq->n, gm, dom->arr, dom->n, mpl, &mplsc);

      p7_refmx_CountTrace(tr, tct);
      p7_refmx_Copy(tct, ppa);
      p7_refmx_Rescale(ppa, 1./(float)s);

      mximage_plot_refmx(mxi, ppa, (ia ? ia : 1), (ib ? ib : sq->n), (ka ? ka : 1), (kb ? kb : gm->M));
      mximage_plot_trace(mxi, tr,  (ia ? ia : 1), (ib ? ib : sq->n), (ka ? ka : 1), (kb ? kb : gm->M), blackpixel);
      
      mximage_Write(outfp, mxi);

      fclose(outfp);
      free(outfile);
      p7_trace_Reuse(tr);
      p7_coords2_Reuse(dom);
      p7_refmx_Reuse(mpl);
    }
  
  /* Frame N+1 of the movie is the actual pp matrix */
  if (esl_sprintf(&outfile, "%s.%06d.ppm", outpfx, N+1) != eslOK) p7_Fail("allocation failed for outfile");
  if ((outfp = fopen(outfile, "w"))                    == NULL) p7_Fail("failed to open %s for writing", outfile);

  mximage_plot_refmx(mxi, pp, (ia ? ia : 1), (ib ? ib : sq->n), (ka ? ka : 1), (kb ? kb : gm->M));

  mximage_Write(outfp, mxi);

  fclose(outfp);
  free(outfile);

  mximage_Destroy(mxi);
  p7_refmx_Destroy(tct);
  p7_refmx_Destroy(ppa);

  p7_trace_Destroy(tr);
  p7_refmx_Destroy(mpl);
  p7_refmx_Destroy(pp);
  p7_refmx_Destroy(vit);
  p7_refmx_Destroy(bck);
  p7_refmx_Destroy(fwd);
  free(wrk);
  p7_coords2_Destroy(dom);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
}



/*****************************************************************
 * mximage routines 
 *****************************************************************/

/* Create a ppm image for a matrix of dimensions <nrow>*<ncol>,
 * with each matrix cell a square <cellsize> pixels large.
 * The image, then, is <nrow*cellsize> pixels high by 
 * <ncol*cellsize> pixels wide.
 */

static MXIMAGE *
mximage_Create(int nrow, int ncol, int cellsize)
{
  MXIMAGE *mxi = NULL;
  int      i;
  int      status;

  ESL_ALLOC(mxi, sizeof(MXIMAGE));
  mxi->ppmdata = NULL;

  mxi->cellsize = cellsize;
  mxi->nrowpix  = nrow * cellsize;
  mxi->ncolpix  = ncol * cellsize;

  ESL_ALLOC(mxi->ppmdata, sizeof(RGBPIXEL *) * mxi->nrowpix);
  mxi->ppmdata[0] = NULL;

  ESL_ALLOC(mxi->ppmdata[0], sizeof(RGBPIXEL) * mxi->nrowpix * mxi->ncolpix);
  for (i = 1; i < mxi->nrowpix; i++)
    mxi->ppmdata[i] = mxi->ppmdata[0] + (i * mxi->ncolpix);
  
  return mxi;

 ERROR:
  return NULL;
}

static void
mximage_Destroy(MXIMAGE *mxi)
{
  if (mxi) {
    if (mxi->ppmdata) {
      if (mxi->ppmdata[0]) free(mxi->ppmdata[0]);
      free(mxi->ppmdata);
    }
    free(mxi);
  }
  return;
}

int
mximage_SetAll(MXIMAGE *mxi, RGBPIXEL color)
{
  int i;
  for (i = 0; i < mxi->nrowpix * mxi->ncolpix; i++)
    mxi->ppmdata[0][i] = color;
  return eslOK;
}

int
mximage_SetCell(MXIMAGE *mxi, int i, int j, RGBPIXEL color)
{
  int imgi,imgj;
  for (imgi = i*mxi->cellsize; imgi < (i+1)*mxi->cellsize; imgi++)
    for (imgj = j*mxi->cellsize; imgj < (j+1)*mxi->cellsize; imgj++)
      mxi->ppmdata[imgi][imgj] = color;
  return eslOK;
}

int
mximage_Write(FILE *ofp, MXIMAGE *mxi)
{
  fprintf(ofp, "P6\n");
  fprintf(ofp, "%d %d\n", mxi->ncolpix, mxi->nrowpix);
  fprintf(ofp, "255\n");
  fwrite(mxi->ppmdata[0], sizeof(RGBPIXEL), mxi->nrowpix * mxi->ncolpix, ofp);
  return eslOK;
}

/* transposes matrix, so image is
 * rows[0..M-1], cols[0..L-1], with sequence
 * running horizontally left-right, and model running
 * vertically top-bottom.
 */
int
mximage_plot_refmx(MXIMAGE *mxi, P7_REFMX *pp, int ia, int ib, int ka, int kb)
{
  int   i,k;
  int   bin;
  float val;

  for (i = ia; i <= ib; i++)
    for (k = ka; k <= kb; k++)
      {
	val = 
	  pp->dp[i][k * p7R_NSCELLS + p7R_ML] + pp->dp[i][k * p7R_NSCELLS + p7R_MG] + 
	  pp->dp[i][k * p7R_NSCELLS + p7R_IL] + pp->dp[i][k * p7R_NSCELLS + p7R_IG] + 
	  pp->dp[i][k * p7R_NSCELLS + p7R_DL] + pp->dp[i][k * p7R_NSCELLS + p7R_DG];

	bin = val / 0.1;
	if (bin > 9) bin = 9;

	mximage_SetCell(mxi, k-ka, i-ia, heatmap_colors[bin]); /* note k <=> i transpose */
      }
  return eslOK;
}
  
int    
mximage_plot_trace(MXIMAGE *mxi, P7_TRACE *tr, int ia, int ib, int ka, int kb, RGBPIXEL color)
{
  int i,k,z;
  
  for (z = 0; z < tr->N; z++)
    {
      i = (tr->i[z] ? tr->i[z] : i); /* last i we emitted: coord we'll use for D states */
      k = tr->k[z];

      if (   p7_trace_IsMain(tr->st[z])   /* any {MDI}{LG} */
	     && ( i >= ia && i <= ib) 
	     && ( k >= ka && k <= kb))
	mximage_SetCell(mxi, k-ka, i-ia, color);
    }
  return eslOK;
}

/*****************************************************************
 * Functions duplicated from decoding-study
 *****************************************************************/

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

#if 0
static int
test_output(void)
{
  MXIMAGE *mxi = mximage_Create(10, 20, 1);
  FILE    *ofp = fopen("test.ppm", "w");
  int      i;

  mximage_SetAll(mxi, heatmap_colors[0]);

  for (i = 0; i < 10; i++)
    mximage_SetCell(mxi, i,    i, heatmap_colors[9]);
  mximage_Write(ofp, mxi);
		 
  mximage_Destroy(mxi);
  fclose(ofp);
  return eslOK;
}
#endif

/* 
 * 
 * Notes on what we're going to do.
 * 
 * 1. Generate each frame in ppm binary format
 *    P6\n
 *    3 2\n
 *    255\n
 *    255   0   0     0 255   0     0   0 255
 *    255 255   0   255 255 255     0   0   0
 *    
 * 2. We need the heatmap scheme in RGB instead of CMYK
 *    colorbrewer2.org, 9-class red, plus blue for zero class
 *    r = 222  255 254 252 252 251 239 203 165 103   
 *    g = 235  245 224 187 146 106  59  24  15   0
 *    b = 247  240 210 161 114  74  44  29  21  13
 *    
 *  alt., for blue: try 198,219,239;
 *                   or 158,202,225
 *                   or 107,174,214
 *                   
 * http://jupiter.ethz.ch/~pjt/makingMovies.html
 * 
 * 3. ImageMagick to convert each frame
 *       convert -quality 100 image123.ppm image123.{jpg,gif}
 *       
 * 4. One option for converting frames to animation: 
 *       convert -delay 6 -quality 95 test*ppm movie.mpg      
 *       convert -delay 20 test*ppm movie.gif
 *       
 *       
 */


/*
0. Set image to background (blue)
1. Write a DP matrix heat map to an image 
2. Write a trace to an image, in some color (black)
3. Save image to file in ppm format
4. matrix_image
   write one cell (in specified color)

  5. Add trace to a count-based pp matrix.  p7_refmx_CountTrace()
  6. Copy count matrix to pp - need to write p7_refmx_Copy()
  7. Renormalize pp matrix. p7_refmx_Rescale()
*/
