/* hmmer statsim :: scoring profile HMM(s) against simulated sequences
 *
 * Testbed for statistical distributions of HMMER4 scores on synthetic sequences.
 */
#include <h4_config.h>

#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dsq.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_histogram.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_subcmd.h"

#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"

#include "calibrate.h"
#include "general.h"
#include "ssvfilter.h"

static ESL_OPTIONS statsim_options[] = {
   /* name          type        default  env   range  toggles reqs incomp  help                                             docgroup*/
  { "-h",          eslARG_NONE,  FALSE,  NULL,  NULL,  NULL,  NULL, NULL,  "show brief help on version and usage",             1 },
  { "-o",       eslARG_OUTFILE,   NULL,  NULL,  NULL,  NULL,  NULL, NULL,  "direct output to file <f>, not stdout",            1 },
  { "-s",           eslARG_INT,    "0",  NULL,"n>=0",  NULL,  NULL, NULL,  "set random number generator seed to <n>",          1 },
  { "-L",           eslARG_INT,  "100",  NULL, "n>0",  NULL,  NULL, NULL,  "length of random target seqs",                     1 },
  { "-N",           eslARG_INT, "1000",  NULL, "n>0",  NULL,  NULL, NULL,  "number of random target seqs",                     1 },

  { "--EsL",        eslARG_INT,  "200",  NULL, "n>0",  NULL,  NULL, NULL,  "length of sequences for SSV Gumbel mu fit",        1 },   
  { "--EsN",        eslARG_INT,  "200",  NULL, "n>0",  NULL,  NULL, NULL,  "number of sequences for SSV Gumbel mu fit",        1 },   

  { "--xyfile", eslARG_OUTFILE,   NULL,  NULL,  NULL,  NULL,  NULL, NULL,  "output P(S>x) survival plots to <f> in xy format", 1 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct statsim_cfg_s {
  uint32_t rng_seed;
  int      N;
  int      L;

  int      EsL;
  int      EsN;

  FILE   *ofp;          // tabular summary output
  FILE   *survfp;       // optional output of survival plots (xmgrace .xy format)
};

static void process_workunit(struct statsim_cfg_s *cfg, const H4_PROFILE *hmm, double *scores, double *ret_mu, double *ret_lambda);
static void output_results  (struct statsim_cfg_s *cfg, char *hmmname, double *scores, double emu, double elambda);


int
h4_cmd_statsim(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  struct statsim_cfg_s cfg;
  ESL_GETOPTS  *go      = esl_subcmd_CreateDefaultApp(topcmd, sub, statsim_options, argc, argv);
  char         *hmmfile = esl_opt_GetArg(go, 1);
  char         *outfile = NULL;
  H4_HMMFILE   *hfp     = NULL;
  H4_PROFILE   *hmm     = NULL;
  ESL_ALPHABET *abc     = NULL;
  double       *scores  = NULL;
  int           nhmm    = 0;
  double        emu, elambda;
  int           status;

  cfg.rng_seed = esl_opt_GetInteger(go, "-s");
  cfg.N        = esl_opt_GetInteger(go, "-N");
  cfg.L        = esl_opt_GetInteger(go, "-L");
  cfg.EsN      = esl_opt_GetInteger(go, "--EsN");
  cfg.EsL      = esl_opt_GetInteger(go, "--EsL");

  if (( outfile = esl_opt_GetString(go, "-o")) != NULL)
    {
      if (( cfg.ofp = fopen(outfile, "w")) == NULL)
        esl_fatal("Failed to open file %s for -o output", outfile);
    }
  else cfg.ofp = stdout;

  if (( outfile = esl_opt_GetString(go, "--xyfile")) != NULL)
    {
      if (( cfg.survfp = fopen(outfile, "w")) == NULL)
        esl_fatal("Failed to open file %s for --xyfile output", outfile);
    }
  else cfg.survfp = NULL;

  ESL_ALLOC(scores, sizeof(double) * cfg.N);

  status = h4_hmmfile_Open(hmmfile, NULL, &hfp);
  if (status != eslOK) esl_fatal("Error: failed to open %s for reading profile HMM(s)\n%s\n", strcmp(hmmfile, "-") == 0? "<stdin>" : hmmfile, hfp->errmsg);

  while ((status = h4_hmmfile_Read(hfp, &abc, &hmm)) == eslOK)
    {
      process_workunit(&cfg, hmm, scores, &emu, &elambda);
      output_results  (&cfg, hmm->name, scores, emu, elambda);
      nhmm++;
    }
  if      (status == eslEFORMAT)      esl_fatal("Parse failed, bad profile HMM file format in %s:\n   %s",  strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status == eslEOF && !nhmm) esl_fatal("Empty input? No profile HMM found in %s\n",                strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status != eslEOF)          esl_fatal("Unexpected error reading profile HMM from %s (code %d)\n", strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, status);

  status = eslOK;
 ERROR:
  free(scores);
  if (cfg.ofp != stdout) fclose(cfg.ofp);
  if (cfg.survfp)        fclose(cfg.survfp);
  h4_hmmfile_Close(hfp);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return status;
}



static void
process_workunit(struct statsim_cfg_s *cfg, const H4_PROFILE *hmm, double *scores, double *ret_mu, double *ret_lambda)
{
  ESL_RANDOMNESS *rng = esl_randomness_Create(cfg->rng_seed);
  ESL_DSQ        *dsq = NULL;
  H4_MODE        *mo  = NULL;
  int             i;
  float           sc;
  double          mu, lambda;
  int             status;

  if (!(hmm->flags & h4_HASVECS)) esl_fatal("profile not vectorized");

  if ((mo = h4_mode_Create()) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (cfg->L + 2));

  // SSV only uses mo->nullsc. It doesn't use any H4_MODE transitions.
  if ((status = h4_mode_SetLength(mo, cfg->L)) != eslOK) goto ERROR;

  h4_lambda(hmm, &lambda);
  h4_ssv_mu(rng, hmm, cfg->EsL, cfg->EsN, lambda, &mu);

  for (i = 0; i < cfg->N; i++)
    {
      esl_rsq_xfIID(rng, hmm->f, hmm->abc->K, cfg->L, dsq);
      h4_ssvfilter(dsq, cfg->L, hmm, &sc);   // can return eslERANGE if we overflowed SSV filter score range, and sc is only a lower bound
      scores[i] = sc - mo->nullsc;
    }

  *ret_mu     = mu;
  *ret_lambda = lambda;
 ERROR:
  h4_mode_Destroy(mo);
  esl_randomness_Destroy(rng);
  free(dsq);
}


static void
output_results(struct statsim_cfg_s *cfg, char *hmmname, double *scores, double emu, double elambda)
{
  ESL_HISTOGRAM *h = esl_histogram_CreateFull(-50., 50., 0.2);
  double         x10;               // 10th highest score
  double         mu, lambda, E10;   // mu, lambda from ML Gumbel fitting; and E-value of 10th highest score using them
  double         mufix,  E10fix;    //  ... using lambda fixed at log 2
  double         mufix2, E10fix2;   //  ... using H4-estimated lambda (with edge correction)
  double         E10e;              //  ... using H4 estimates for mu, lambda
  double         tailp = 1.0;       // Unused for now. Viterbi scores fit to complete distribution.
  int            i;
 
  for (i = 0; i < cfg->N; i++)
    esl_histogram_Add(h, scores[i]);

  esl_histogram_GetRank(h, 10, &x10);
  
  /* mu, lambda, E10 fields: ML Gumbel fit to the observed data */
  esl_gumbel_FitComplete(scores, cfg->N, &mu, &lambda);
  E10    = cfg->N * esl_gumbel_surv(x10, mu, lambda); 

  /* mufix, E10fix fields:   assume lambda = log2; fit an ML mu to the data */
  esl_gumbel_FitCompleteLoc(scores, cfg->N, 0.693147, &mufix);
  E10fix = cfg->N * esl_gumbel_surv(x10, mufix, 0.693147); 

  /* mufix2, E10fix2 fields: assume edge-corrected H3 lambda estimate; fit ML mu */
  esl_gumbel_FitCompleteLoc(scores, cfg->N, elambda, &mufix2);
  E10fix2 = cfg->N * esl_gumbel_surv(x10, mufix2, elambda); 
      
  /* emu, elambda, E10e:     use H4 estimates (emu, elambda) */
  E10e    = cfg->N * esl_gumbel_surv(x10, emu,   elambda); 
      
  fprintf(cfg->ofp, "%-20s  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", 
          hmmname, tailp, mu, lambda, E10, mufix, E10fix, mufix2, E10fix2, emu, elambda, E10e);

  
  if (cfg->survfp) {
    esl_histogram_PlotSurvival(cfg->survfp, h);
    esl_gumbel_Plot(cfg->survfp, mu,    lambda,   esl_gumbel_surv, h->xmin - 5., h->xmax + 5., 0.1);
    esl_gumbel_Plot(cfg->survfp, mufix, 0.693147, esl_gumbel_surv, h->xmin - 5., h->xmax + 5., 0.1);
  }

}
  

