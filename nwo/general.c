#include "h4_config.h"

#include <stdlib.h>

#include "easel.h"
#include "esl_getopts.h"


/* Function:  h4_CreateDefaultApp()
 * Synopsis:  Initialize a small HMMER program such as a test driver or example.
 * Incept:    SRE, Fri 22 Jun 2018
 *
 * Purpose:   Based on <esl_getopts_CreateDefaultApp()> but specialized for HMMER.
 *            See documentation in <easel/esl_getopts.c>.
 * 
 *            The <options> list is assumed to contain '--help' (help) and '--version' (version).
 *
 * Args:      options : Easel options structures
 *            nargs   : expected number of command line arguments, aside from options
 *            argc    : number of command line arguments in <argv>
 *            argv    : command line arguments, inclusive of argv[0] program name
 *            banner  : short one-line description of program
 *            usage   : short one-line summary of command line usage
 *
 * Returns:   ptr to a new ESL_GETOPTS object.
 * 
 *            On command line parsing errors, print an error message
 *            to <stderr> and exit with abnormal (1) status.
 *            
 *            When the standard `--help` option is seen, print a help
 *            page to <stdout> and exit with normal (0) status.
 *            
 *            When the standard `--version` option is seen, print the
 *            HMMER version number and exit with normal (0) status.
 *            
 * Throws:    On exceptions and errors, calls esl_fatal() to exit with
 *            abnormal (1) status.
 */
ESL_GETOPTS *
h4_CreateDefaultApp(ESL_OPTIONS *options, int nargs, int argc, char **argv, char *banner, char *usage)
{
  ESL_GETOPTS *go       = NULL;
  char        *progname = NULL;

  //if (  h4_Init()                               != eslOK) goto ERROR;  // we'll probably need this someday
  if (  esl_FileTail(argv[0], FALSE, &progname) != eslOK) goto ERROR; // remove any leading directory path
  if (( go = esl_getopts_Create(options))       == NULL)  goto ERROR;

  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK) 
    {
      if ( esl_fprintf(stderr,   "Failed to parse command line: %s\n", go->errbuf)                    != eslOK) goto ERROR;
      if ( esl_fprintf(stderr, "\nUsage: %s %s\n", progname, usage)                                   != eslOK) goto ERROR;
      if ( esl_fprintf(stderr, "\nTo see more help on available options, do %s --help\n\n", progname) != eslOK) goto ERROR;
      free(progname);
      exit(1);
    }

  if (esl_opt_GetBoolean(go, "--version") == TRUE)
    {
      if ( esl_printf("%s\n", HMMER_VERSION) != eslOK) goto ERROR;
      exit(0);
    }

  if (esl_opt_GetBoolean(go, "--help") == TRUE) 
    {
      if (banner)
	{
	  if ( esl_printf("%s : %s\n", progname, banner)                               != eslOK) goto ERROR;
	  if ( esl_printf("HMMER %s (%s); %s\n", HMMER_VERSION, HMMER_DATE, HMMER_URL) != eslOK) goto ERROR;
	}

      if (usage)
	{
	  if ( esl_printf("\nUsage: %s %s\n", progname, usage) != eslOK) goto ERROR;
	}

      if ( esl_printf("\nOptions:\n")                != eslOK) goto ERROR;
      if ( esl_opt_DisplayHelp(stdout, go, 0, 2, 80) != eslOK) goto ERROR;

      free(progname);
      exit(0);
    }

  if (nargs != -1 && esl_opt_ArgNumber(go) != nargs) 
    {
      if ( esl_fprintf(stderr, "Incorrect number of command line arguments.")                         != eslOK) goto ERROR;
      if ( esl_fprintf(stderr, "\nUsage: %s %s\n", progname, usage)                                   != eslOK) goto ERROR;
      if ( esl_fprintf(stderr, "\nTo see more help on available options, do %s --help\n\n", progname) != eslOK) goto ERROR;
      free(progname);
      exit(1);
    }

  free(progname);
  return go;

 ERROR: // only reached if a nondefault, nonfatal exception handler happens to be in use
  esl_fatal("internal failure in h4_CreateDefaultApp()");
}
