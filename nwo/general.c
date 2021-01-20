#include "h4_config.h"

#include <stdlib.h>

#include "easel.h"
#include "esl_getopts.h"

#include "logsum.h"
#include "simdvec.h"

/* Function:  h4_Init()
 * Synopsis:  Initialize a new H4 process or thread.
 * Incept:    SRE, Fri 17 May 2019
 *
 * Purpose:   Initialize a new H4 process or thread:
 *
 *              - initialize lookup table for fast table-driven
 *                approximation of the log-sum-exp operation, used by
 *                non-SIMD Forward/Backward code.
 *
 *              - set processor flags to turn off denormalized
 *                floating point math; performance penalty is too 
 *                high.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
h4_Init(void)
{
  // Initialize the lookup table for fast approximate log-sum-exp operations (see logsum.md)
  h4_logsum_Init();

  // Set processor flags to turn off denormalized floating point math (see simdvec.md) 
  h4_simdvec_Init();

  return eslOK;
}



/* Function:  h4_CreateDefaultApp()
 * Synopsis:  Initialize a small HMMER program such as a test driver or example.
 * Incept:    SRE, Fri 22 Jun 2018
 *
 * Purpose:   Based on <esl_getopts_CreateDefaultApp()> but specialized for HMMER.
 *            See documentation in <easel/esl_getopts.c>.
 * 
 *            The <options> list is assumed to contain '-h' (help) and '--version' (version).
 *
 *            If <nargs> is >=0, the number of command line arguments is checked. If you
 *            don't want this checked (for example, if the program allows a variable number
 *            of arguments), set <nargs> to -1.
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
 *            When the standard `-h` option is seen, print a help
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

  if ( h4_Init() != eslOK) goto ERROR; 

  if ( esl_FileTail(argv[0], FALSE, &progname) != eslOK) goto ERROR; // remove any leading directory path
  if (( go = esl_getopts_Create(options))       == NULL)  goto ERROR;

  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK) 
    {
      if ( esl_fprintf(stderr,   "Failed to parse command line: %s\n", go->errbuf)                    != eslOK) goto ERROR;
      if ( esl_fprintf(stderr, "\nUsage: %s %s\n", progname, usage)                                   != eslOK) goto ERROR;
      if ( esl_fprintf(stderr, "\nTo see more help on available options, do %s -h\n\n", progname) != eslOK) goto ERROR;
      free(progname);
      exit(1);
    }

  if (esl_opt_GetBoolean(go, "--version") == TRUE)
    {
      if ( esl_printf("%s\n", HMMER_VERSION) != eslOK) goto ERROR;
      exit(0);
    }

  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
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
      if ( esl_fprintf(stderr, "Incorrect number of command line arguments.")                     != eslOK) goto ERROR;
      if ( esl_fprintf(stderr, "\nUsage: %s %s\n", progname, usage)                               != eslOK) goto ERROR;
      if ( esl_fprintf(stderr, "\nTo see more help on available options, do %s -h\n\n", progname) != eslOK) goto ERROR;
      free(progname);
      exit(1);
    }

  free(progname);
  return go;

 ERROR: // only reached if a nondefault, nonfatal exception handler happens to be in use
  esl_fatal("internal failure in h4_CreateDefaultApp()");
}



/* Function:  h4_AminoFrequencies()
 *
 * Purpose:   Fills a vector <f> with amino acid background frequencies,
 *            in [A..Y] alphabetic order, same order that Easel digital
 *            alphabet uses. Caller must provide <f> allocated for at
 *            least 20 floats.
 *            
 *            These were updated 4 Sept 2007, from Swiss-Prot 50.8,
 *            (Oct 2006), counting over 85956127 (86.0M) residues.
 *
 * Returns:   <eslOK> on success.
 */
int
h4_AminoFrequencies(float *f)
{
  f[0] = 0.0787945;	// A
  f[1] = 0.0151600;	// C
  f[2] = 0.0535222;	// D
  f[3] = 0.0668298;	// E
  f[4] = 0.0397062;	// F
  f[5] = 0.0695071;	// G
  f[6] = 0.0229198;	// H
  f[7] = 0.0590092;	// I
  f[8] = 0.0594422;	// K
  f[9] = 0.0963728;	// L
  f[10]= 0.0237718;	// M
  f[11]= 0.0414386;	// N
  f[12]= 0.0482904;	// P
  f[13]= 0.0395639;	// Q
  f[14]= 0.0540978;	// R
  f[15]= 0.0683364;	// S
  f[16]= 0.0540687;	// T
  f[17]= 0.0673417;	// V
  f[18]= 0.0114135;	// W
  f[19]= 0.0304133;	// Y
  return eslOK;
}
