/* Compare a profile to a sequence, with the ability to output a
 * variety of data (statistics, visualizations) about what's going 
 * on with posterior decoding, compared to alignment.
 * 
 * One of two vehicles for asking questions about the importance of using
 * ensemble calculations instead of optimal alignment. See also:
 * decoding-hunt.
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

static ESL_OPTIONS options[] = {
  /* name           type           default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",          eslARG_NONE,   FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",                    0 },
  { "-i",          eslARG_STRING, FALSE,  NULL, NULL,   NULL,  NULL, NULL, "when plotting, restrict plot to rows <i1>..<i2>",         0 },
  { "-k",          eslARG_STRING, FALSE,  NULL, NULL,   NULL,  NULL, NULL, "when plotting, restrict plot to columns <k1>..<k2>",      0 },
  { "--diplot",    eslARG_OUTFILE, NULL,  NULL, NULL,   NULL,  NULL, NULL, "output domain inference plot to <f> (xmgrace format)",    0 },
  { "--heatmap",   eslARG_OUTFILE, NULL,  NULL, NULL,   NULL,  NULL, NULL, "output heat map of decoding matrix to <f> (PostScript)",  0 },
  { "--local",     eslARG_NONE,    NULL,  NULL, NULL,   NULL,  NULL, NULL, "configure the model for local-only alignment",            0 },
  { "--seed",      eslARG_INT,      "0",  NULL, NULL,   NULL,  NULL, NULL, "set the random number seed to <n>",                       0 },  
  { "--splot",     eslARG_OUTFILE, NULL,  NULL, NULL,   NULL,  NULL, NULL, "plot suboptimal path, domain inference stype, to <f>",    0 },
  { "--smap",      eslARG_OUTFILE, NULL,  NULL, NULL,   NULL,  NULL, NULL, "plot suboptimal path in matrix heat map style to <f>",    0 },
  { "--vplot",     eslARG_OUTFILE, NULL,  NULL, NULL,   NULL,  NULL, NULL, "plot Viterbi path in domain inference plot style to <f>", 0 },
  { "--vmap",      eslARG_OUTFILE, NULL,  NULL, NULL,   NULL,  NULL, NULL, "plot Viterbi path in matrix heat map style to <f>",       0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "testing importance of ensemble calculations";

static int parse_coord_string(const char *cstring, int *ret_start, int *ret_end);

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go         = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng        = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  char           *hmmfile    = esl_opt_GetArg(go, 1);
  P7_HMMFILE     *hfp        = NULL;
  ESL_ALPHABET   *abc        = NULL;
  char           *seqfile    = esl_opt_GetArg(go, 2);
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
  P7_TRACE       *tr         = p7_trace_Create();
  P7_COORDS2     *dom        = p7_coords2_Create(16, 64);
  char           *difile     = esl_opt_GetString(go, "--diplot");
  FILE           *difp       = NULL;
  char           *heatfile   = esl_opt_GetString(go, "--heatmap");
  FILE           *heatfp     = NULL;
  char           *vdifile    = esl_opt_GetString(go, "--vplot");
  FILE           *vdifp      = NULL;
  char           *vheatfile  = esl_opt_GetString(go, "--vmap");
  FILE           *vheatfp    = NULL;
  char           *sdifile    = esl_opt_GetString(go, "--splot");
  FILE           *sdifp      = NULL;
  char           *sheatfile  = esl_opt_GetString(go, "--smap");
  FILE           *sheatfp    = NULL;
  float           fsc,   vsc,   bsc,   mplsc;
  int             d;
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
  if (esl_opt_GetBoolean(go, "--local")) p7_profile_ConfigLocal(gm, hmm, bg, 100);
  else                                   p7_profile_Config     (gm, hmm, bg);
  
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

  /* Open other optional plot and data files for output */
  if (difile     && (difp     = fopen(difile,     "w")) == NULL) p7_Fail("Failed to open domain inference plot file %s for writing.", difile);
  if (heatfile   && (heatfp   = fopen(heatfile,   "w")) == NULL) p7_Fail("Failed to open heatmap plot file %s for writing.", heatfile);
  if (vdifile    && (vdifp    = fopen(vdifile,    "w")) == NULL) p7_Fail("Failed to open Viterbi domain inference plot file %s for writing.", vdifile);
  if (vheatfile  && (vheatfp  = fopen(vheatfile,  "w")) == NULL) p7_Fail("Failed to open Viterbi heatmap plot file %s for writing.", vheatfile);
  if (sdifile    && (sdifp    = fopen(sdifile,    "w")) == NULL) p7_Fail("Failed to open suboptimal domain inference plot file %s for writing.", sdifile);
  if (sheatfile  && (sheatfp  = fopen(sheatfile,  "w")) == NULL) p7_Fail("Failed to open suboptimal heatmap plot file %s for writing.", sheatfile);

  /* Set length models */
  p7_bg_SetLength     (bg,  sq->n);
  p7_profile_SetLength(gm,  sq->n);

  p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

  p7_ReferenceViterbi (sq->dsq, sq->n, gm, vit, tr,  &vsc);
  p7_ReferenceForward (sq->dsq, sq->n, gm, fwd,      &fsc);   
  p7_ReferenceBackward(sq->dsq, sq->n, gm, bck,      &bsc);
  p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, bck, pp);

  /* Create a coord2 list, crudely. We can do better later. */
  p7_trace_Index(tr);
  p7_coords2_SetFromTrace(dom, tr);

  p7_ReferenceMPLForward(sq->dsq, sq->n, gm, dom->arr, dom->n, mpl, &mplsc);
	     
  if (difp)    p7_refmx_PlotDomainInference(difp,   pp, (ia ? ia : 1), (ib ? ib : sq->n), tr);
  if (heatfp)  p7_refmx_PlotHeatMap        (heatfp, pp, (ia ? ia : 1), (ib ? ib : sq->n), (ka ? ka : 1), (kb ? kb : gm->M));

  if (vdifp)   p7_trace_PlotDomainInference(vdifp,   tr, (ia ? ia : 1), (ib ? ib : sq->n));
  if (vheatfp) p7_trace_PlotHeatMap        (vheatfp, tr, (ia ? ia : 1), (ib ? ib : sq->n), (ka ? ka : 1), (kb ? kb : gm->M));

  printf("Viterbi:\n");
  for (d = 0; d < dom->n; d++)
    printf("   domain %3d : %5d ... %5d\n", d+1, dom->arr[d].n1, dom->arr[d].n2);
  printf("tracesc:  %.2f\n", vsc);
  printf("    MPL:  %.2f\n", mplsc);
  printf("    fwd:  %.2f\n", fsc);
  printf("      %%:  %.2g\n", exp(mplsc-fsc));

  p7_trace_Reuse(tr);
  p7_coords2_Reuse(dom);
  p7_refmx_Reuse(mpl);

  /* One suboptimal alignment - sampled from fwd matrix */
  p7_reference_trace_Stochastic(rng, NULL, gm, fwd, tr);	  
  p7_trace_Index(tr);
  p7_coords2_SetFromTrace(dom, tr);

  p7_trace_Score(tr, sq->dsq, gm, &vsc);
  p7_ReferenceMPLForward(sq->dsq, sq->n, gm, dom->arr, dom->n, mpl, &mplsc);

  printf("suboptimal (rng seed %" PRIu32 ")\n", esl_randomness_GetSeed(rng));
  for (d = 0; d < dom->n; d++)
    printf("   domain %3d : %5d ... %5d\n", d+1, dom->arr[d].n1, dom->arr[d].n2);
  printf("tracesc:  %.2f\n", vsc);
  printf("    MPL:  %.2f\n", mplsc);
  printf("    fwd:  %.2f\n", fsc);
  printf("      %%:  %.2g\n", exp(mplsc-fsc));

  if (sdifp)   p7_trace_PlotDomainInference(sdifp,   tr, (ia ? ia : 1), (ib ? ib : sq->n));
  if (sheatfp) p7_trace_PlotHeatMap        (sheatfp, tr, (ia ? ia : 1), (ib ? ib : sq->n), (ka ? ka : 1), (kb ? kb : gm->M));

  p7_trace_Reuse(tr);
  p7_coords2_Reuse(dom);
  p7_refmx_Reuse(mpl);


  if (difp)    fclose(difp);
  if (heatfp)  fclose(heatfp);
  if (vdifp)   fclose(vdifp);
  if (vheatfp) fclose(vheatfp);
  if (sdifp)   fclose(sdifp);
  if (sheatfp) fclose(sheatfp);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(mpl);
  p7_refmx_Destroy(pp);
  p7_refmx_Destroy(vit);
  p7_refmx_Destroy(bck);
  p7_refmx_Destroy(fwd);

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




