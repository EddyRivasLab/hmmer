#include <h4_config.h>

#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_subcmd.h"

#include "h4_profile.h"
#include "h4_hmmfile.h"

#include "build.h"

#define ALPHOPTS "--amino,--dna,--rna"        // Exclusive options for alphabet choice      (default = autodetect)
#define WGTOPTS  "--wnone,--wgiven"           // Exclusive options for relative weighting   (default = PB weights)
#define EFFOPTS  "--enone,--eset"             // Exclusive options for effective seq number (default = entropy weighting)

static ESL_OPTIONS build_options[] = {
  /* name             type      default                         env  range      toggles reqs incomp  help                                                docgroup*/
  { "-h",           eslARG_NONE,  FALSE,                        NULL, NULL,      NULL, NULL, NULL,     "show brief help on version and usage",                     1 },
  /* assert the content of the input MSA file, instead of autoguessing */
  { "--informat",   eslARG_STRING, NULL,                        NULL, NULL,      NULL, NULL, NULL,     "specify the input MSA file is in format <s>",              2 }, 
  { "--dna",        eslARG_NONE,  FALSE,                        NULL, NULL,  ALPHOPTS, NULL, NULL,     "use DNA alphabet",                                         2 },
  { "--rna",        eslARG_NONE,  FALSE,                        NULL, NULL,  ALPHOPTS, NULL, NULL,     "use RNA alphabet",                                         2 },
  { "--amino",      eslARG_NONE,  FALSE,                        NULL, NULL,  ALPHOPTS, NULL, NULL,     "use protein alphabet",                                     2 },
  /* options controlling model architecture (consensus column determination) */
  { "--hand",       eslARG_NONE,  FALSE,                        NULL, NULL,      NULL, NULL, NULL,     "manual construction (requires reference annotation line)", 3 },
  { "--fragthresh", eslARG_REAL,  ESL_STR(h4BUILD_FRAGTHRESH),  NULL, "0<=x<=1", NULL, NULL, "--hand", "seq is fragment is aspan/alen < fragthresh",               3 },
  { "--symfrac",    eslARG_REAL,  ESL_STR(h4BUILD_SYMFRAC),     NULL, "0<=x<=1", NULL, NULL, "--hand", "col is consensus if nres/(nres+ngap) >= symfrac",          3 },
  { "--nsamp",      eslARG_INT,   ESL_STR(h4BUILD_NSAMP),       NULL, "n>=0",    NULL, NULL, "--hand", "number of seqs to sample for determining consensus cols",  3 },
  { "--maxfragdens",eslARG_REAL,  ESL_STR(h4BUILD_MAXFRAGDENS), NULL, "0<=x<=1", NULL, NULL, "--hand", "if sample has > maxdens*nseq frags, use all seqs",         3 },
  /* alternative relative weighting strategies */
  { "--wnone",      eslARG_NONE,  FALSE,                        NULL, NULL,   WGTOPTS, NULL, NULL,     "no relative weighting; set them all to 1",                 4 },
  { "--wgiven",     eslARG_NONE,  FALSE,                        NULL, NULL,   WGTOPTS, NULL, NULL,     "use weights given in input MSA file",                      4 },
  /* effective sequence number (absolute weighting) */
  { "--enone",      eslARG_NONE,  FALSE,                        NULL, NULL,   EFFOPTS, NULL, NULL,     "no effective seq # weighting; use nseq",                   5 },
  { "--eset",       eslARG_REAL,  FALSE,                        NULL, NULL,   EFFOPTS, NULL, NULL,     "set eff seq # for all models to <x>",                      5 },
  { "--etarg",      eslARG_REAL,  FALSE,                        NULL, "x>0",     NULL, NULL, EFFOPTS,  "set entropy weighting target to <x> bits",                 5 },
  { "--esigma",     eslARG_REAL,  FALSE,                        NULL, "x>0",     NULL, NULL, EFFOPTS,  "set sigma parameter for entropy weighting to <x>",         5 },
  /* other options */
  { "--seed",       eslARG_INT,   ESL_STR(h4BUILD_RNGSEED),     NULL, NULL,      NULL, NULL, NULL,     "set RNG seed to <n>; if 0, then seed quasirandomly",       6 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static ESL_GETOPTS *process_cmdline(const char *topcmd, const ESL_SUBCMD *sub, const ESL_OPTIONS *suboptions, int argc, char **argv);

int
h4_cmd_build(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS  *go      = process_cmdline(topcmd, sub, build_options, argc, argv);
  char         *msafile = esl_opt_GetArg(go, 1);
  char         *hmmfile = esl_opt_GetArg(go, 2);
  int           infmt   = eslMSAFILE_UNKNOWN;
  H4_BUILD_CONFIG *cfg  = NULL;
  ESL_ALPHABET *abc     = NULL;
  ESL_MSAFILE  *afp     = NULL;
  ESL_MSA      *msa     = NULL;
  H4_PROFILE   *hmm     = NULL;
  FILE         *ofp     = NULL;
  int           nali    = 0;
  int           status;

  if (esl_opt_IsOn(go, "--informat") &&
      (infmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));

  if      (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO); 
  
  /* Open MSA file before opening output <hmmfile>.
   * If user messed up order of arguments, we want to detect that before overwriting their data.
   */
  if (( status = esl_msafile_Open(&abc, msafile, NULL, infmt, NULL, &afp)) != eslOK)
    esl_msafile_OpenFailure(afp, status);

  if (strcmp(hmmfile, "-") == 0) ofp = stdout;
  else if (( ofp = fopen(hmmfile, "w")) == NULL)
    esl_fatal("couldn't open output HMM file %s for writing", hmmfile);

  cfg = h4_build_config_Create(abc);

  if      (esl_opt_GetBoolean(go, "--hand"))   cfg->arch_strategy = h4BUILD_ARCH_GIVEN;

  if      (esl_opt_GetBoolean(go, "--wnone"))  cfg->wgt_strategy  = h4BUILD_WGT_NONE;
  else if (esl_opt_GetBoolean(go, "--wgiven")) cfg->wgt_strategy  = h4BUILD_WGT_GIVEN;

  if      (esl_opt_GetBoolean(go, "--enone"))  cfg->effn_strategy = h4BUILD_EFFN_NONE;
  else if (esl_opt_IsOn      (go, "--eset")) { cfg->effn_strategy = h4BUILD_EFFN_GIVEN; cfg->effn_set = esl_opt_GetReal(go, "--eset"); }

  cfg->fragthresh    = esl_opt_GetReal   (go, "--fragthresh");
  cfg->symfrac       = esl_opt_GetReal   (go, "--symfrac");
  cfg->nsamp         = esl_opt_GetInteger(go, "--nsamp");
  cfg->maxfragdens   = esl_opt_GetReal   (go, "--maxfragdens");
  cfg->seed          = esl_opt_GetInteger(go, "--seed");

  if (esl_opt_IsOn(go, "--etarg"))  cfg->re_target = esl_opt_GetReal(go, "--etarg");
  if (esl_opt_IsOn(go, "--esigma")) cfg->re_sigma  = esl_opt_GetReal(go, "--esigma");
  
  while ((status = esl_msafile_Read(afp, &msa)) == eslOK)
    {
      nali++;

      if (!msa->name) esl_msa_SetName(msa, "test", -1);  // SRE: temporary

      h4_Build(cfg, msa, &hmm, NULL);
  
      h4_hmmfile_Write(ofp, hmm);

      esl_msa_Destroy(msa);
      h4_profile_Destroy(hmm);
    }
  if (nali == 0 || status != eslEOF) esl_msafile_ReadFailure(afp, status); /* a convenience, like esl_msafile_OpenFailure() */
  esl_msafile_Close(afp);
  
  if (ofp != stdout) fclose(ofp);
  esl_alphabet_Destroy(abc);
  h4_build_config_Destroy(cfg);
  esl_getopts_Destroy(go);
  return 0;
}


/* The build command has a multipart help page, so this is a
 * customized copy of esl_subcmd_CreateDefaultApp().
 */
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
      esl_printf("\noptions for asserting MSA file format, instead of autoguessing:\n");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      esl_printf("\noptions for determining consensus columns:\n");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      esl_printf("\noptions for relative sequence weighting:\n");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      esl_printf("\noptions for effective sequence number:\n");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      esl_printf("\nother options:\n");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
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


