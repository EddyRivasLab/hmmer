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

#include "h4_filtermx.h"
#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"

#include "calibrate.h"
#include "general.h"
#include "ssvfilter.h"
#include "vitfilter.h"
#include "fbfilter.h"

#define ALGORITHMS "--ssvfilter,--vitfilter,--fwdfilter,--vit,--fwd,--asc"

#define h4STATSIM_DO_SSVFILTER  1
#define h4STATSIM_DO_VITFILTER  2
#define h4STATSIM_DO_FWDFILTER  3
#define h4STATSIM_DO_VIT        4
#define h4STATSIM_DO_FWD        5 
#define h4STATSIM_DO_ASC        6

static ESL_OPTIONS statsim_options[] = {
   /* name          type        default  env   range  toggles reqs incomp  help                                             docgroup*/
  { "-h",          eslARG_NONE,  FALSE,  NULL,  NULL,  NULL,  NULL, NULL,  "show brief help on version and usage",             1 },
  { "-o",       eslARG_OUTFILE,   NULL,  NULL,  NULL,  NULL,  NULL, NULL,  "direct output to file <f>, not stdout",            1 },
  { "-s",           eslARG_INT,    "0",  NULL,"n>=0",  NULL,  NULL, NULL,  "set random number generator seed to <n>",          1 },
  { "-L",           eslARG_INT,  "100",  NULL, "n>0",  NULL,  NULL, NULL,  "length of random target seqs",                     1 },
  { "-N",           eslARG_INT, "1000",  NULL, "n>0",  NULL,  NULL, NULL,  "number of random target seqs",                     1 },

  /* Which type of score to collect - which algorithm to run */
  { "--ssvfilter", eslARG_NONE,  FALSE,  NULL,  NULL,ALGORITHMS,NULL, NULL, "collect SSV filter scores (local)",                               2 },
  { "--vitfilter", eslARG_NONE,  FALSE,  NULL,  NULL,ALGORITHMS,NULL, NULL, "collect Viterbi filter scores (local)",                           2 },
  { "--fwdfilter", eslARG_NONE,  FALSE,  NULL,  NULL,ALGORITHMS,NULL, NULL, "collect Forward filter scores (local)",                           2 },
  { "--vit",       eslARG_NONE,  FALSE,  NULL,  NULL,ALGORITHMS,NULL, NULL, "collect Viterbi scores (dual-mode glocal/local)",                 2 },
  { "--fwd",       eslARG_NONE,"default",NULL,  NULL,ALGORITHMS,NULL, NULL, "collect Forward scores (dual-mode glocal/local)",                 2 },
  { "--asc",       eslARG_NONE,  FALSE,  NULL,  NULL,ALGORITHMS,NULL, NULL, "collect ASC Forward scores (dual-mode glocal/local)",             2 },
  { "--reference", eslARG_NONE,  FALSE,  NULL,  NULL,  NULL,    NULL, NULL, "use reference implementation, not production sparse|vector code", 2 },
  
  /* Outputs */
  { "--xyfile", eslARG_OUTFILE,   NULL,  NULL,  NULL,  NULL,  NULL, NULL,  "output P(S>x) survival plots to <f> in xy format", 3 },

  /* Statistical fitting parameters */
  { "--EsL",        eslARG_INT,  "200",  NULL, "n>0",  NULL,  NULL, NULL,  "length of sequences for SSV Gumbel mu fit",           4 },   
  { "--EsN",        eslARG_INT,  "200",  NULL, "n>0",  NULL,  NULL, NULL,  "number of sequences for SSV Gumbel mu fit",           4 },   
  { "--EvL",        eslARG_INT,  "200",  NULL, "n>0",  NULL,  NULL, NULL,  "length of sequences for Vit filter Gumbel mu fit",    4 },   
  { "--EvN",        eslARG_INT,  "200",  NULL, "n>0",  NULL,  NULL, NULL,  "number of sequences for Vit filter Gumbel mu fit",    4 },   
  { "--EfL",        eslARG_INT,  "100",  NULL, "n>0",  NULL,  NULL, NULL,  "length of sequences for Fwd filter exp tail tau fit", 4 },   
  { "--EfN",        eslARG_INT,  "200",  NULL, "n>0",  NULL,  NULL, NULL,  "number of sequences for Fwd filter exp tail tau fit", 4 },   
  { "--Eft",       eslARG_REAL, "0.04",  NULL,"0<x<=1",NULL,  NULL, NULL,  "tail mass for Fwd filter exponential tail tau fit",   4 },   
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct statsim_cfg_s {
  uint32_t rng_seed;
  int      N;
  int      L;

  int      which_score_type;
  int      use_reference_impl;

  FILE   *ofp;          // tabular summary output
  FILE   *survfp;       // optional output of survival plots (xmgrace .xy format)

  int      EsL;
  int      EsN;
  int      EvL;
  int      EvN;
  int      EfL;
  int      EfN;
  int      Eft;
};

static ESL_GETOPTS *process_cmdline(const char *topcmd, const ESL_SUBCMD *sub, const ESL_OPTIONS *suboptions, int argc, char **argv);
static void         initialize_cfg(ESL_GETOPTS *go, struct statsim_cfg_s *cfg);

static void process_workunit(struct statsim_cfg_s *cfg, const H4_PROFILE *hmm, double *scores, double *ret_mu, double *ret_lambda);
static void output_results  (struct statsim_cfg_s *cfg, char *hmmname, double *scores, double emu, double elambda);


int
h4_cmd_statsim(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  struct statsim_cfg_s cfg;
  ESL_GETOPTS  *go      = process_cmdline(topcmd, sub, statsim_options, argc, argv);
  char         *hmmfile = esl_opt_GetArg(go, 1);
  H4_HMMFILE   *hfp     = NULL;
  H4_PROFILE   *hmm     = NULL;
  ESL_ALPHABET *abc     = NULL;
  double       *scores  = NULL;
  int           nhmm    = 0;
  double        emu, elambda;
  int           status;

  initialize_cfg(go, &cfg);

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

static ESL_GETOPTS *
process_cmdline(const char *topcmd, const ESL_SUBCMD *sub, const ESL_OPTIONS *suboptions, int argc, char **argv)
{
  ESL_GETOPTS *go        = esl_getopts_Create(suboptions);
  char        *lastslash = strrchr(topcmd, '/');

  if (lastslash) topcmd = lastslash+1;

  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK) 
    {
      esl_printf("Failed to parse command line: %s\n", go->errbuf);
      esl_printf("Usage:\n  %s %s %s\n", topcmd, sub->subcmd, sub->usage);
      esl_printf("\nTo see more help on available options, do `%s %s -h`\n\n", topcmd, sub->subcmd);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      esl_printf("%s %s : %s\n", topcmd, sub->subcmd, sub->description);
      esl_printf("\nUsage:\n  %s %s %s\n", topcmd, sub->subcmd, sub->usage);
      esl_printf("\nOptions:\n");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      esl_printf("\noptions for choosing which type of score to collect:\n");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      esl_printf("\noptions for outputs:\n");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      esl_printf("\ncustomizing parameters for statistical fitting:\n");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != sub->nargs) 
    {
      esl_printf("Incorrect number of command line arguments.\n");
      esl_printf("Usage:\n  %s %s %s\n", topcmd, sub->subcmd, sub->usage);
      esl_printf("\nTo see more help on available options, do `%s %s -h`\n\n", topcmd, sub->subcmd);
      exit(1);
    }
  return go;
}

static void
initialize_cfg(ESL_GETOPTS *go, struct statsim_cfg_s *cfg)
{
  char *outfile = NULL;
  
  cfg->rng_seed = esl_opt_GetInteger(go, "-s");
  cfg->N        = esl_opt_GetInteger(go, "-N");
  cfg->L        = esl_opt_GetInteger(go, "-L");

  if      (esl_opt_GetBoolean(go, "--ssvfilter"))  cfg->which_score_type = h4STATSIM_DO_SSVFILTER;
  else if (esl_opt_GetBoolean(go, "--vitfilter"))  cfg->which_score_type = h4STATSIM_DO_VITFILTER;
  else if (esl_opt_GetBoolean(go, "--fwdfilter"))  cfg->which_score_type = h4STATSIM_DO_FWDFILTER;
  else if (esl_opt_GetBoolean(go, "--vit"))        cfg->which_score_type = h4STATSIM_DO_VIT;
  else if (esl_opt_GetBoolean(go, "--fwd"))        cfg->which_score_type = h4STATSIM_DO_FWD;
  else if (esl_opt_GetBoolean(go, "--asc"))        cfg->which_score_type = h4STATSIM_DO_ASC;
  else esl_fatal("can't happen. no default score type set in options?");

  cfg->use_reference_impl = esl_opt_GetBoolean(go, "--reference");

  if (( outfile = esl_opt_GetString(go, "-o")) != NULL)
    {
      if (( cfg->ofp = fopen(outfile, "w")) == NULL)
        esl_fatal("Failed to open file %s for -o output", outfile);
    }
  else cfg->ofp = stdout;

  if (( outfile = esl_opt_GetString(go, "--xyfile")) != NULL)
    {
      if (( cfg->survfp = fopen(outfile, "w")) == NULL)
        esl_fatal("Failed to open file %s for --xyfile output", outfile);
    }
  else cfg->survfp = NULL;

  cfg->EsN = esl_opt_GetInteger(go, "--EsN");
  cfg->EsL = esl_opt_GetInteger(go, "--EsL");
  cfg->EvN = esl_opt_GetInteger(go, "--EvN");
  cfg->EvL = esl_opt_GetInteger(go, "--EvL");
  cfg->EfN = esl_opt_GetInteger(go, "--EfN");
  cfg->EfL = esl_opt_GetInteger(go, "--EfL");
  cfg->Eft = esl_opt_GetReal   (go, "--Eft");
}


static void
process_workunit(struct statsim_cfg_s *cfg, const H4_PROFILE *hmm, double *scores, double *ret_mu, double *ret_lambda)
{
  ESL_RANDOMNESS *rng = esl_randomness_Create(cfg->rng_seed);
  ESL_DSQ        *dsq = NULL;
  H4_MODE        *mo  = NULL;
  H4_FILTERMX    *fx  = NULL;
  int             i;
  float           sc;
  double          mu, lambda;
  int             status;

  if (!(hmm->flags & h4_HASVECS)) esl_fatal("profile not vectorized");

  if ((mo = h4_mode_Create()) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (cfg->L + 2));
  if ((status = h4_mode_SetLength(mo, cfg->L)) != eslOK) goto ERROR;  


  /* Allocate any DP matrix(s) we need */
  if (cfg->use_reference_impl)
    {
      esl_fatal("unimplemented");
    }
  else
    {
      if (cfg->which_score_type == h4STATSIM_DO_VITFILTER)
        fx = h4_filtermx_Create(hmm->M);
    }
  

  /* (Re)calibrate the model as directed */
  h4_lambda(hmm, &lambda);
  switch (cfg->which_score_type) {
  case h4STATSIM_DO_SSVFILTER: h4_ssv_mu(rng, hmm, cfg->EsL, cfg->EsN, lambda, &mu); break;
  case h4STATSIM_DO_VITFILTER: h4_vit_mu(rng, hmm, cfg->EvL, cfg->EvN, lambda, &mu); break;
  default: esl_fatal("unimplemented");
  }


  /* Collect scores on N random sequences */
  for (i = 0; i < cfg->N; i++)
    {
      esl_rsq_xfIID(rng, hmm->f, hmm->abc->K, cfg->L, dsq);

      if (cfg->use_reference_impl)
        {
          esl_fatal("unimplemented");
        }
      else
        {
          switch (cfg->which_score_type) {
          case h4STATSIM_DO_SSVFILTER: h4_ssvfilter(dsq, cfg->L, hmm, mo,     &sc); break;  // can return eslERANGE if we overflowed SSV filter score range, and sc is only a lower bound
          case h4STATSIM_DO_VITFILTER: h4_vitfilter(dsq, cfg->L, hmm, mo, fx, &sc); break;
          case h4STATSIM_DO_FWDFILTER: esl_fatal("unimplemented");
          case h4STATSIM_DO_VIT:       esl_fatal("unimplemented");
          case h4STATSIM_DO_FWD:       esl_fatal("unimplemented");
          case h4STATSIM_DO_ASC:       esl_fatal("unimplemented");
          default: esl_fatal("can't happen. no such score type?");
          }
        }
      scores[i] = sc;
    }

  // status == eslOK if all is well, and we flow through to cleanup
  *ret_mu     = mu;
  *ret_lambda = lambda;
 ERROR:
  h4_filtermx_Destroy(fx);
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
  

