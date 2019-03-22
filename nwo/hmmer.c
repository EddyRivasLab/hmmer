#include "h4_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_subcmd.h"

#include "cmd_build.h"

ESL_SUBCMD subcommands[] = {
  { h4_cmd_build, "build", 2, "[-options] <msafile> <hmmfile>", "build profile(s) from multiple alignment(s)" },
};

static ESL_OPTIONS top_options[] = {
  /* name         type          default  env  range tog's   reqs incomp  help                       docgroup*/
  { "-h",         eslARG_NONE,   FALSE, NULL, NULL,  NULL, NULL,   NULL, "show overall brief help summary", 1 },
  { "--version",  eslARG_NONE,   FALSE, NULL, NULL,  NULL, NULL,   NULL, "show version number",             1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

  

static int
top_usage(const char *topcmd)
{
  char *lastslash = strrchr(topcmd, '/');
  if (lastslash) topcmd = lastslash+1;

  if (printf("Usage:\n")                                                               < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if (printf("  %s -h                : show overall brief help summary\n",     topcmd) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if (printf("  %s --version         : show version number\n",                 topcmd) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if (printf("  %s <cmd> -h          : show brief help for a HMMER command\n", topcmd) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if (printf("  %s <cmd> [<args>...] : run a HMMER command\n",                 topcmd) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  return eslOK;
}

static int
top_help(const char *topcmd)
{
  int   ncmds     =  sizeof(subcommands) / sizeof(ESL_SUBCMD);
  int   i;
  int   status;

  if ( printf("HMMER: biological sequence analysis using profile hidden Markov models\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if ( printf("Version %s (%s): %s\n\n", HMMER_VERSION, HMMER_DATE, HMMER_URL)            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if (( status = top_usage(topcmd)) != eslOK) return status;
  if ( printf("\nAvailable commands:\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  for (i = 0; i < ncmds; i++)
    if ( printf("  %-12s %s\n", subcommands[i].subcmd, subcommands[i].description)        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");

  return eslOK;
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(top_options);
  int ncmds =  sizeof(subcommands) / sizeof(ESL_SUBCMD);
  int idx;
  int status;
 
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("Failed to parse command line: %s\n\n",  go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("Failed to parse command line: %s\n\n",  go->errbuf);
  
  if (esl_opt_GetBoolean(go, "--version") == TRUE) { printf("%s\n", HMMER_VERSION); status = eslOK; goto DONE; }
  if (esl_opt_GetBoolean(go, "-h")        == TRUE) { status = top_help(argv[0]);    goto DONE; }
  if (argc - go->optind == 0)                      { status = top_help(argv[0]);    goto DONE; }

  for (idx = 0; idx < ncmds; idx++)
    if (strcmp(go->argv[go->optind], subcommands[idx].subcmd) == 0) break;
  if (idx == ncmds) { status = top_usage(argv[0]); goto DONE; }

  status = subcommands[idx].func(argv[0], &subcommands[idx], argc-go->optind, argv+go->optind);
  
 DONE:
  esl_getopts_Destroy(go);
  return status;
}

