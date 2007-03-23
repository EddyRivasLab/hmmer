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
#include "esl_vectorops.h"
#include "esl_dmatrix.h"
#include "esl_ratematrix.h"

#include "hmmer.h"

#define MODELOPTS "-m,-u,-s,-S"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",   1 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0",     NULL,      NULL,    NULL, "length of random target seqs",           1 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0",     NULL,      NULL,    NULL, "number of random target seqs",           1 },
  { "-2",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "configure profile in old HMMER2 style",  1 },
  { "--ips",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,    NULL, "output .ps heatmap of i endpoints to <f>", 1 },
  { "--kps",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,    NULL, "output .ps heatmap of k endpoints to <f>", 1 },

  { "-m",        eslARG_INFILE,  NULL, NULL, NULL,      NULL,      NULL, "-s,-u,-M,-S", "input HMM from file <f>",                2 },
  { "-s",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "-m,-u,-S",    "sample a uniform-transition HMM",        2 },
  { "-u",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "-s,-m,-S",    "sample an ungapped HMM",                 2 },
  { "-S",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "-m,-s,-u",    "sample a sequence and make HMM from it", 2 },

  { "-M",        eslARG_INT,     "50", NULL, "n>0",     NULL,      NULL,    "-m", "length of a sampled query HMM or seq",         3 },
  { "-t",        eslARG_REAL,   "1.0", NULL, "x>0",     NULL,      "-S",    NULL, "branch length, for parameterizing seq-based query", 3 },
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
  ESL_DMATRIX     *imx     = NULL;
  ESL_DMATRIX     *kmx     = NULL;
  int              i1,i2,k1,k2;	       /* endpoints of a domain                   */
  int              which;	       /* counter over domains                    */
  int              ndom;	       /* total hits counted into imx, kmx        */
  double           expect;
  FILE            *fp;

  char            *hmmfile;            /* file to read HMM(s) from                */
  int              do_swlike;	       /* TRUE to sample a uniform-transition HMM */
  int              do_ungapped;	       /* TRUE to sample an ungapped HMM          */
  int              do_seqsample;       /* TRUE to make model from random sequence */
  int              L;		       /* length of generated seqs                */
  int              M;		       /* length of a sampled HMM                 */
  int              N;		       /* number of seqs to generate              */
  double           t;		       /* branch length for seq query model       */
  int              do_h2;              /* TRUE to use H2 exit/entry configuration */
  char            *ipsfile;	       /* i endpoint heatmap file                 */
  char            *kpsfile;	       /* k endpoint heatmap file                 */

  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);
  if (esl_opt_IsSet(go, "-h")) {
    puts(usage);
    puts("\ngeneral options are:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 2 = indentation; 80=textwidth*/
    puts("\nalternatives for random query model :");
    esl_opt_DisplayHelp(stdout, go, 2, 2, 80); /* 2 = indentation; 80=textwidth*/
    puts("\ncontrolling character of random query model :");
    esl_opt_DisplayHelp(stdout, go, 3, 2, 80); /* 2 = indentation; 80=textwidth*/
    return eslOK;
  }

  esl_opt_GetStringOption (go, "-m",     &hmmfile);
  esl_opt_GetBooleanOption(go, "-s",     &do_swlike);
  esl_opt_GetBooleanOption(go, "-u",     &do_ungapped);
  esl_opt_GetBooleanOption(go, "-S",     &do_seqsample);
  esl_opt_GetIntegerOption(go, "-L",     &L);
  esl_opt_GetIntegerOption(go, "-M",     &M);
  esl_opt_GetIntegerOption(go, "-N",     &N);
  esl_opt_GetBooleanOption(go, "-2",     &do_h2);
  esl_opt_GetDoubleOption (go, "-t",     &t);
  esl_opt_GetStringOption (go, "--ips",  &ipsfile);
  esl_opt_GetStringOption (go, "--kps",  &kpsfile);

  if (esl_opt_ArgNumber(go) != 0) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return eslFAIL;
  }

  /*****************************************************************
   * Initializations, including opening and reading the HMM file 
   *****************************************************************/

  if ((r       = esl_randomness_CreateTimeseeded()) == NULL)  esl_fatal("Failed to create rng");

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
  else if (do_seqsample)
    {
      abc = esl_alphabet_Create(eslAMINO);    
      ESL_DMATRIX *Q     = esl_dmatrix_Create(abc->K,abc->K);
      ESL_DMATRIX *P     = esl_dmatrix_Create(abc->K,abc->K);
      ESL_DSQ     *query = malloc(sizeof(ESL_DSQ) * (M+2));
      P7_BG       *bg    = p7_bg_Create(abc);
  
      esl_rmx_SetWAG(Q, NULL);
      esl_dmx_Exp(Q, t, P);

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
  for (ndom = 0, nseq = 1; nseq <= N; nseq++) 
    {
      if (esl_rnd_xfIID(r, hmm->bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");
      if (p7_Viterbi(dsq, L, hmm->gm, mx, tr, &sc)     != eslOK) esl_fatal("viterbi failed");
      if (p7_trace_Validate(tr, abc, dsq, errbuf)      != eslOK) esl_fatal("trace validation failed:\n%s", errbuf);
      /*printf("score = %6.2f bits\n", p7_SILO2Bitscore(sc));*/

      for (which = 1; p7_trace_GetDomainCoords(tr, which, &i1, &i2, &k1, &k2) != eslEOD; which++)
	{
	  imx->mx[i1-1][i2-1] += 1.;
	  kmx->mx[k1-1][k2-1] += 1.;
	  ndom++;
	  /* printf("   i: %d..%d  k: %d..%d\n", i1, i2, k1, k2); */
	}
    }

  
  /*****************************************************************
   * Convert i,k mx's to log-odds relative to uniform distribution.
   *****************************************************************/
  
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

  /*****************************************************************
   * Cleanup, exit.
   *****************************************************************/
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
