#include <h4_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_subcmd.h"

extern int h4_cmd_build   (const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv);
extern int h4_cmd_kiteline(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv);
extern int h4_cmd_statsim (const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv);

static void top_usage(const char *topcmd);
static void top_help (const char *topcmd);

ESL_SUBCMD subcommands[] = {
  { h4_cmd_build,    "build",    2, "[-options] <msafile> <hmmfile>", "build profile(s) from multiple alignment(s)"       },
  { h4_cmd_kiteline, "kiteline", 2, "[-options] <hmmfile> <seqfile>", "prototype: search profile(s) against sequences"    },
  { h4_cmd_statsim,  "statsim",  1, "[-options] <hmmfile>",           "test statistical distributions of sequence scores" },
};

static ESL_OPTIONS top_options[] = {
  /* name         type          default  env  range tog's   reqs incomp  help                       docgroup*/
  { "-h",         eslARG_NONE,   FALSE, NULL, NULL,  NULL, NULL,   NULL, "show overall brief help summary", 1 },
  { "--version",  eslARG_NONE,   FALSE, NULL, NULL,  NULL, NULL,   NULL, "show version number",             1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};



int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(top_options);
  int ncmds =  sizeof(subcommands) / sizeof(ESL_SUBCMD);
  int idx;
 
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("Failed to parse command line: %s\n\n",  go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("Failed to parse command line: %s\n\n",  go->errbuf);
  
  if      (esl_opt_GetBoolean(go, "--version") == TRUE)  esl_printf("%s\n", HMMER_VERSION);
  else if (esl_opt_GetBoolean(go, "-h")        == TRUE)  top_help(argv[0]); 
  else if (argc - go->optind == 0)                       top_help(argv[0]); 
  else {
    for (idx = 0; idx < ncmds; idx++)
      if (strcmp(go->argv[go->optind], subcommands[idx].subcmd) == 0) break;

    if (idx == ncmds) top_usage(argv[0]); 
    else              subcommands[idx].func(argv[0], &subcommands[idx], argc-go->optind, argv+go->optind);
  }
  
  esl_getopts_Destroy(go);
  return 0;
}


static void
top_usage(const char *topcmd)
{
  char *lastslash = strrchr(topcmd, '/');
  if (lastslash) topcmd = lastslash+1;

  esl_printf("Usage:\n");
  esl_printf("  %s -h                : show overall brief help summary\n",     topcmd);
  esl_printf("  %s --version         : show version number\n",                 topcmd);
  esl_printf("  %s <cmd> -h          : show brief help for a HMMER command\n", topcmd);
  esl_printf("  %s <cmd> [<args>...] : run a HMMER command\n",                 topcmd);
}

static void
top_help(const char *topcmd)
{
  int ncmds =  sizeof(subcommands) / sizeof(ESL_SUBCMD);
  int i;

  esl_printf("HMMER: biological sequence analysis using profile hidden Markov models\n");
  esl_printf("Version %s (%s): %s\n\n", HMMER_VERSION, HMMER_DATE, HMMER_URL);
  top_usage(topcmd);
  esl_printf("\nAvailable commands:\n");
  for (i = 0; i < ncmds; i++)
    esl_printf("  %-12s %s\n", subcommands[i].subcmd, subcommands[i].description);
}
