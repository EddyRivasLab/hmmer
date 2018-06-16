#include "h4_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"

extern int h4_cmd_build(int argc, char **argv);


static ESL_OPTIONS options[] = {
  /* name         type          default  env  range tog's   reqs incomp  help                       docgroup*/
  { "--help",     eslARG_NONE,   FALSE, NULL, NULL,  NULL, NULL,   NULL, "show overall brief help summary", 1 },
  { "--version",  eslARG_NONE,   FALSE, NULL, NULL,  NULL, NULL,   NULL, "show version number",             1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


static struct subcommands_s
  {
    char *cmdname;
    int (*func)(int argc, char **argv);
    char *helpline;
  }
  subcommands[] =
    {
      { "build", h4_cmd_build, "construct new profile HMM(s) from multiple sequence alignment(s)" },
    };
  

static int
usage(void)
{
  if (printf("Usage:\n")                                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if (printf("  hmmer --help            : show overall brief help summary\n")            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if (printf("  hmmer --version         : show version number\n")                        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if (printf("  hmmer <cmd> --help      : show brief help for a HMMER command\n")        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if (printf("  hmmer <cmd> [<args>...] : run a HMMER command\n")                        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  return eslOK;
}

static int
help(void)
{
  int ncmds =  sizeof(subcommands) / sizeof(struct subcommands_s);
  int i;
  int status;

  if ( printf("HMMER: biological sequence analysis using profile hidden Markov models\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if ( printf("Version %s (%s): %s\n\n", HMMER_VERSION, HMMER_DATE, HMMER_URL)            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  if (( status = usage()) != eslOK) return status;
  if ( printf("\nCommon commands:\n")                                                     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");
  for (i = 0; i < ncmds; i++)
    if ( printf("  %-12s %s\n", subcommands[i].cmdname, subcommands[i].helpline)          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "printf failed");

  return eslOK;
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int ncmds =  sizeof(subcommands) / sizeof(struct subcommands_s);
  int i;
  int status;
 
  if (esl_opt_ProcessEnvironment(go)         != eslOK) { if ((status = esl_printf("Failed to process environment: %s\n\n", go->errbuf)) != eslOK) goto DONE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) { if ((status = esl_printf("Failed to parse command line: %s\n\n",  go->errbuf)) != eslOK) goto DONE; }
  if (esl_opt_VerifyConfig(go)               != eslOK) { if ((status = esl_printf("Failed to parse command line: %s\n\n",  go->errbuf)) != eslOK) goto DONE; }

  if (esl_opt_GetBoolean(go, "--version") == TRUE) { status = esl_printf("%s\n", HMMER_VERSION); goto DONE; }
  if (esl_opt_GetBoolean(go, "--help")    == TRUE) { status = help();                            goto DONE; }
  if (argc - go->optind == 0)                      { status = help();                            goto DONE; }

  for (i = 0; i < ncmds; i++)
    {
      if (strcmp(go->argv[go->optind], subcommands[i].cmdname) == 0)
	{
	  status = subcommands[i].func(argc-go->optind, argv+go->optind);
	  goto DONE;
	}
    }

  if (i == ncmds) status = usage();
  
 DONE:
  esl_getopts_Destroy(go);
  return status;
}

