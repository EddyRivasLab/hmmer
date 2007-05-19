/* main() for scoring a profile HMM against simulated random sequences
 * 
 * SRE, Thu Feb  8 14:45:23 2007 [UA 1130 from Vancouver]
 * SVN $Id$
 */

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_vectorops.h"
#include "esl_dmatrix.h"
#include "esl_ratematrix.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"

#include "hmmer.h"

#define ALGORITHMS "--fwd,--viterbi,--island"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",   1 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0",     NULL,      NULL,    NULL, "length of random target seqs",           1 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0",     NULL,      NULL,    NULL, "number of random target seqs",           1 },
  { "-2",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "configure profile in old HMMER2 style",  1 },
  { "-o",        eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,    NULL, "output score histogram to <f>, xmgrace xy format", 1 },
  { "--xfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,    NULL, "output bitscores as binary double-prec vector to <f>", 1},
  { "--ips",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,    NULL, "output .ps heatmap of i endpoints to <f>", 1 },
  { "--kps",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,    NULL, "output .ps heatmap of k endpoints to <f>", 1 },

  { "--viterbi", eslARG_NONE,  "TRUE", NULL, NULL, ALGORITHMS,     NULL,    NULL, "Score seqs with the Viterbi algorithm",  2},
  { "--fwd",     eslARG_NONE,   FALSE, NULL, NULL, ALGORITHMS,     NULL,    NULL, "Score seqs with the Forward algorithm",  2},
  { "--island",  eslARG_NONE,   FALSE, NULL, NULL, ALGORITHMS,     NULL,    NULL, "Score seqs with the island algorithm",   2},

  { "-m",        eslARG_INFILE,  NULL, NULL, NULL,      NULL,      NULL, "-s,-u,-M,-S", "input HMM from file <f>",                3 },
  { "-s",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "-m,-u,-S",    "sample a uniform-transition HMM",        3 },
  { "-u",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "-s,-m,-S",    "sample an ungapped HMM",                 3 },
  { "-S",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "-m,-s,-u",    "sample a sequence and make HMM from it", 3 },

  { "-M",        eslARG_INT,     "50", NULL, "n>0",     NULL,      NULL,    "-m", "length of a sampled query HMM or seq",         4 },
  { "-t",        eslARG_REAL,   "2.0", NULL, "x>0",     NULL,      "-S",    NULL, "branch length, for parameterizing seq-based query", 4 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "hmmsim [-options]";

int
main(int argc, char **argv)
{
  int status;			       /* status of a function call               */
  ESL_ALPHABET    *abc     = NULL;     /* sequence alphabet                       */
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                 */
  ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                    */
  P7_HMM          *hmm     = NULL;     /* HMM to emit from                        */
  P7_TRACE        *tr      = NULL;     /* sampled trace                           */
  P7_GMX          *mx      = NULL;     /* DP matrix                               */
  ESL_DSQ         *dsq     = NULL;     /* sampled digital sequence                */
  int              sc;		       /* a Viterbi score (internal: SILO)        */
  int              nseq;	       /* counter over sequences                  */
  char             errbuf[eslERRBUFSIZE];
  ESL_HISTOGRAM   *h       = NULL;
  ESL_DMATRIX     *imx     = NULL;
  ESL_DMATRIX     *kmx     = NULL;
  int              i1,i2,k1,k2;	       /* endpoints of a domain                   */
  int              which;	       /* counter over domains                    */
  int              ndom;	       /* total hits counted into imx, kmx        */
  double           expect;
  FILE            *fp;
  double           tot_ilen, tot_klen;

  char            *hmmfile;            /* file to read HMM(s) from                */
  int              do_swlike;	       /* TRUE to sample a uniform-transition HMM */
  int              do_ungapped;	       /* TRUE to sample an ungapped HMM          */
  int              do_seqsample;       /* TRUE to make model from random sequence */
  int              L;		       /* length of generated seqs                */
  int              M;		       /* length of a sampled HMM                 */
  int              N;		       /* number of seqs to generate              */
  double           t;		       /* branch length for seq query model       */
  int              do_h2;              /* TRUE to use H2 exit/entry configuration */
  enum { DO_VITERBI, DO_FORWARD, DO_ISLAND } algorithm_choice;
  char            *histfile;	       /* histogram output file                   */
  char            *ipsfile;	       /* i endpoint heatmap file                 */
  char            *kpsfile;	       /* k endpoint heatmap file                 */
  char            *xfile;	       /* data vector file                        */

  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    puts(usage);
    puts("\ngeneral options are:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 2 = indentation; 80=textwidth*/
    puts("\nalternative scoring algorithms :");
    esl_opt_DisplayHelp(stdout, go, 2, 2, 80); /* 2 = indentation; 80=textwidth*/
    puts("\nalternatives for random query model :");
    esl_opt_DisplayHelp(stdout, go, 3, 2, 80); /* 2 = indentation; 80=textwidth*/
    puts("\ncontrolling character of random query model :");
    esl_opt_DisplayHelp(stdout, go, 4, 2, 80); /* 2 = indentation; 80=textwidth*/
    return eslOK;
  }

  hmmfile       = esl_opt_GetString (go, "-m");
  do_swlike     = esl_opt_GetBoolean(go, "-s");
  do_ungapped   = esl_opt_GetBoolean(go, "-u");
  do_seqsample  = esl_opt_GetBoolean(go, "-S");
  L             = esl_opt_GetInteger(go, "-L");
  M             = esl_opt_GetInteger(go, "-M");
  N             = esl_opt_GetInteger(go, "-N");
  do_h2         = esl_opt_GetBoolean(go, "-2");
  t             = esl_opt_GetReal   (go, "-t");
  histfile      = esl_opt_GetString (go, "-o");
  ipsfile       = esl_opt_GetString (go, "--ips");
  kpsfile       = esl_opt_GetString (go, "--kps");
  xfile         = esl_opt_GetString (go, "--xfile");
  if      (esl_opt_GetBoolean(go, "--viterbi")) algorithm_choice = DO_VITERBI;
  else if (esl_opt_GetBoolean(go, "--fwd"))     algorithm_choice = DO_FORWARD;
  else if (esl_opt_GetBoolean(go, "--island"))  algorithm_choice = DO_ISLAND;
  
  if (esl_opt_ArgNumber(go) != 0) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return eslFAIL;
  }

  /*****************************************************************
   * Initializations, including opening and reading the HMM file 
   *****************************************************************/

  if ((r       = esl_randomness_CreateTimeseeded()) == NULL)  esl_fatal("Failed to create rng");

  if (algorithm_choice == DO_ISLAND) h = esl_histogram_Create(-50.5, 50.5, 1.0);
  else                               h = esl_histogram_CreateFull(-50.5, 50.5, 1.0);

  if (hmmfile != NULL)		/* Read in an HMM from file... */
    {
      P7_HMMFILE      *hfp     = NULL;    

      status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
      if (status == eslENOTFOUND) esl_fatal("Failed to open hmm file %s for reading.\n", hmmfile);
      else if (status != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", hmmfile);

      if ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslOK) {
	if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
	else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s", hmmfile);
	else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets", hmmfile);
	else                             esl_fatal("Unexpected error in reading HMMs");
      }
      M = hmm->M;
      p7_hmmfile_Close(hfp);
    }
  else if (do_seqsample)	/* or generate a single-sequence HMM from a random query sequence... */
    {
      double   pi[20];
      abc = esl_alphabet_Create(eslAMINO);    
      ESL_DMATRIX *Q     = esl_dmatrix_Create(abc->K,abc->K);
      ESL_DMATRIX *P     = esl_dmatrix_Create(abc->K,abc->K);
      ESL_DSQ     *query = malloc(sizeof(ESL_DSQ) * (M+2));
      P7_BG       *bg    = p7_bg_Create(abc);

      esl_vec_F2D(bg->f, 20, pi);
      esl_rmx_SetWAG(Q, pi);
      esl_dmx_Exp(Q, t, P);

      printf("expected score   of WAG at t=%f  is %f bits\n", t, esl_rmx_ExpectedScore  (P, pi));
      printf("relative entropy of WAG at t=%f  is %f bits\n", t, esl_rmx_RelativeEntropy(P, pi));

      esl_rnd_xfIID(r, bg->f, abc->K, M, query);
      p7_Seqmodel(abc, query, M, P, bg->f, 0.05, 0.5, 0.05, 0.2, &hmm); /* tmi, tii, tmd, tdd */
      
      esl_dmatrix_Destroy(Q);
      esl_dmatrix_Destroy(P);
      free(dsq);
      p7_bg_Destroy(bg);
    }
  else				/* or generate a simulated HMM. */
    {
      abc = esl_alphabet_Create(eslAMINO);    
      if      (do_ungapped) p7_hmm_SampleUngapped(r, M, abc, &hmm);
      else if (do_swlike)   p7_hmm_SampleUniform (r, M, abc, 0.05, 0.5, 0.05, 0.2, &hmm); /* tmi, tii, tmd, tdd */
      else                  p7_hmm_Sample        (r, M, abc, &hmm);
    }

  if ((status  = p7_trace_Create(256, &tr))         != eslOK) esl_fatal("Failed to allocate trace");
  if ((mx      = p7_gmx_Create(hmm->M, L))          == NULL)  esl_fatal("failed to create dp matrix");
  if ((hmm->bg = p7_bg_Create(abc))                 == NULL)  esl_fatal("failed to create null model");
  if ((hmm->gm = p7_profile_Create(hmm->M, abc))    == NULL)  esl_fatal("failed to create profile");
  if ((dsq     = malloc(sizeof(ESL_DSQ) * (L+2)))   == NULL)  esl_fatal("failed to create dsq");
  if ((imx     = esl_dmatrix_CreateUpper(L))        == NULL)  esl_fatal("failed to create imx");
  if ((kmx     = esl_dmatrix_CreateUpper(M))        == NULL)  esl_fatal("failed to create kmx");
  esl_dmatrix_SetZero(imx);
  esl_dmatrix_SetZero(kmx);


  /*****************************************************************
   * Choose how to configure the profile, and do it.
   *****************************************************************/
  if (do_h2) {
    if (p7_H2_ProfileConfig(hmm, hmm->gm, p7_UNILOCAL) != eslOK) esl_fatal("failed to configure profile");
  } else {
    if (p7_ProfileConfig(hmm, hmm->gm, p7_UNILOCAL)    != eslOK) esl_fatal("failed to configure profile");
    if (p7_ReconfigLength(hmm->gm,  L)                 != eslOK) esl_fatal("failed to reconfig profile L");
    if (p7_hmm_Validate    (hmm,     0.0001, NULL)     != eslOK) esl_fatal("whoops, HMM is bad!");
    if (p7_profile_Validate(hmm->gm, 0.0001)           != eslOK) esl_fatal("whoops, profile is bad!");
  }



  /*****************************************************************
   * Align to simulated sequences, collect statistics.
   *****************************************************************/

  tot_ilen = tot_klen = 0.;

  for (ndom = 0, nseq = 1; nseq <= N; nseq++) 
    {
      if (esl_rnd_xfIID(r, hmm->bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");


      if (algorithm_choice == DO_VITERBI) 
	{
	  if (p7_Viterbi(dsq, L, hmm->gm, mx, tr, &sc)     != eslOK) esl_fatal("viterbi failed");
	  if (p7_trace_Validate(tr, abc, dsq, errbuf)      != eslOK) esl_fatal("trace validation failed:\n%s", errbuf);
	  /*printf("score = %6.2f bits\n", p7_SILO2Bitscore(sc));*/
	  esl_histogram_Add(h, p7_SILO2Bitscore(sc));
	  
	  for (which = 1; p7_trace_GetDomainCoords(tr, which, &i1, &i2, &k1, &k2) != eslEOD; which++)
	    {
	      imx->mx[i1-1][i2-1] += 1.;
	      kmx->mx[k1-1][k2-1] += 1.;
	      ndom++;

	      tot_ilen += i2-i1+1;
	      tot_klen += k2-k1+1;
	      /* printf("   i: %d..%d  k: %d..%d\n", i1, i2, k1, k2); */
	    }
	}
      else if (algorithm_choice == DO_FORWARD) 
	{
	  if (p7_Forward(dsq, L, hmm->gm, mx, &sc)     != eslOK) esl_fatal("forward failed");
	  esl_histogram_Add(h, p7_SILO2Bitscore(sc));
	}
      else if (algorithm_choice == DO_ISLAND) 
	{ /* island needs mx for at least 4 rows */
	  if (p7_island_Viterbi(dsq, L, hmm->gm, mx, h)        != eslOK) esl_fatal("island algorithm failed");
	}
    }

  
 /*****************************************************************
  * Output to xmgrace file.
  *****************************************************************/

  if (histfile != NULL) 
    {
      FILE   *hfp;
      double *xv;
      int     n;
      double  mu, lambda;
      double  param[2];

      if ((hfp = fopen(histfile, "w")) == NULL)  esl_fatal("failed to open histogram file %s", histfile);

      if (algorithm_choice == DO_VITERBI)
	{
	  esl_histogram_GetData(h, &xv, &n);

	  esl_gumbel_FitComplete(xv, n, &mu, &lambda);
	  param[0] = mu;
	  param[1] = lambda; 
	  esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, &param);
	  esl_histogram_PlotSurvival(hfp, h);

	  fprintf(stderr, "ML mu,lambda = %f %f\n", mu, lambda);

	  lambda = 0.693;
	  esl_gumbel_FitCompleteLoc(xv, n, 0.693, &mu);
	  param[0] = mu;
	  param[1] = 0.693; 
	  esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, &param);
	  esl_histogram_PlotSurvival(hfp, h);

	  fprintf(stderr, "ML mu (lambda=0.693) = %f\n", mu);
	}

      else if (algorithm_choice == DO_FORWARD)
	{
	  esl_histogram_GetTailByMass(h, 0.01, &xv, &n, NULL);
	  esl_exp_FitComplete(xv, n, &mu, &lambda);

	  param[0] = mu;
	  param[1] = lambda;
	  esl_histogram_SetExpectedTail(h, mu, 0.02, &esl_exp_generic_cdf, &param);
	  esl_histogram_PlotSurvival(hfp, h);

	  param[0] = mu;
	  param[1] = 0.693; 
	  esl_histogram_SetExpectedTail(h, mu, 0.02, &esl_exp_generic_cdf, &param);
	  esl_histogram_PlotSurvival(hfp, h);

	  fprintf(stderr, "mu     = %f\n", mu);
	  fprintf(stderr, "lambda = %f\n", lambda);
	}
      else if (algorithm_choice == DO_ISLAND)
	{
	  double actual_mass;

	  esl_histogram_SetTailByMass(h, 0.1, &actual_mass);
	  esl_exp_FitCompleteBinned(h, &mu, &lambda);

	  param[0] = mu;
	  param[1] = lambda;
	  esl_histogram_SetExpectedTail(h, mu, actual_mass, &esl_exp_generic_cdf, &param);
	  esl_histogram_PlotSurvival(hfp, h);

	  param[0] = mu;
	  param[1] = 0.693; 
	  esl_histogram_SetExpectedTail(h, mu, actual_mass, &esl_exp_generic_cdf, &param);
	  esl_histogram_PlotSurvival(hfp, h);

	  fprintf(stderr, "mu     = %f\n", mu);
	  fprintf(stderr, "lambda = %f\n", lambda);
	}

      fclose(hfp);
    }

  /*****************************************************************
   * Optionally, output raw data
   *****************************************************************/
  
  if (xfile != NULL) 
    {
      double *xv;
      int     n;

      if (algorithm_choice == DO_ISLAND)      esl_fatal("Island stats don't keep raw data");
      if ((fp = fopen(xfile, "wb")) == NULL)  esl_fatal("Failed to open output raw data file %s", xfile);

      esl_histogram_GetData(h, &xv, &n);
      if (fwrite(xv, sizeof(double), n, fp) != n) esl_fatal("failed to save binary data to %s\n", xfile);
      fclose(fp);
    }


  /*****************************************************************
   * Convert i,k mx's to log-odds relative to uniform distribution.
   *****************************************************************/
  
  if (algorithm_choice == DO_VITERBI) 
    {
      expect = (double) ndom * 2. / (double) (L * (L+1));
      for (i1 = 0; i1 < L; i1++)
	for (i2 = i1; i2 < L; i2++)
	  if (expect == 0. && imx->mx[i1][i2] == 0.) imx->mx[i1][i2] = 0.;
	  else if (expect == 0.)             	 imx->mx[i1][i2] = eslINFINITY;
	  else if (imx->mx[i1][i2] == 0.)            imx->mx[i1][i2] = -eslINFINITY;
	  else                                       imx->mx[i1][i2] = log(imx->mx[i1][i2] / expect) / log(2.0);


      expect = (double) ndom * 2. / (double) (M * (M+1));
      for (k1 = 0; k1 < M; k1++)
	for (k2 = k1; k2 < M; k2++)
	  if (expect == 0. && kmx->mx[k1][k2] == 0.) kmx->mx[k1][k2] = 0.;
	  else if (expect == 0.)             	 kmx->mx[k1][k2] = eslINFINITY;
	  else if (kmx->mx[k1][k2] == 0.)            kmx->mx[k1][k2] = -eslINFINITY;
	  else                                       kmx->mx[k1][k2] = log(kmx->mx[k1][k2] / expect) / log(2.0);
   
      /* esl_dmatrix_Dump(stdout, kmx, NULL, NULL); */

      if (kpsfile != NULL) {
	if ((fp = fopen(kpsfile, "w")) == NULL) esl_fatal("Failed to open output postscript file %s", kpsfile);
	dmx_Visualize(fp, kmx, -4., 5.);
	fclose(fp);
      }
      if (ipsfile != NULL) {
	if ((fp = fopen(ipsfile, "w")) == NULL) esl_fatal("Failed to open output postscript file %s", ipsfile);
	dmx_Visualize(fp, imx, -4., 5.); 
	fclose(fp);
      }

      printf("i matrix values range from %f to %f\n", esl_dmx_Min(imx), esl_dmx_Max(imx));
      printf("k matrix values range from %f to %f\n", esl_dmx_Min(kmx), esl_dmx_Max(kmx));
      printf("average alignment length on model: %f\n", tot_klen / (double) ndom);
      printf("average alignment length on seq:   %f\n", tot_ilen / (double) ndom);
    }	 

  /*****************************************************************
   * Cleanup, exit.
   *****************************************************************/
  esl_histogram_Destroy(h);
  esl_dmatrix_Destroy(kmx);
  esl_dmatrix_Destroy(imx);
  free(dsq);
  p7_profile_Destroy(hmm->gm);  
  p7_bg_Destroy(hmm->bg);
  p7_gmx_Destroy(mx);
  p7_trace_Destroy(tr);
  esl_randomness_Destroy(r);

  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;


}
