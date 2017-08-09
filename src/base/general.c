/* A few general routines used throughout HMMER.
 * 
 * Contents:
 *   1. Miscellaneous functions for H3
 *   2. Error handling (Die, Fail)
 */
#include "p7_config.h"

#include <inttypes.h>
#include <math.h>
#include <float.h>
#include <syslog.h>
#include "misc/logsum.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include "dp_vector/simdvec.h"
#include "base/general.h"
#include "easel.h"

#include "esl_getopts.h"





/*****************************************************************
 * 1. Miscellaneous functions for H3
 *****************************************************************/

/* Function:  p7_Init()
 * Synopsis:  Initialize a new HMMER thread or process.
 *
 * Purpose:   Initialize a new HMMER thread or executable process.
 *            
 *            - initialize the lookup table for our fast table-driven
 *              approximation of log-sum-exp.
 *
 *            - set processor flags to turn off denormalized
 *              floating point math; performance penalty is too
 *              high.
 *
 * Args:      none.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_Init(void)
{
  /* We do a lot of log-sum-exp operations using a table-driven
   * approximation; initialize the table.
   */
  p7_FLogsumInit();

  /* Initialize any cpu flags we need to set for SIMD vector processing.
   * On Intel/AMD platforms, for example, we turn off denormalized fp.
   * This is implemented in dp_vector, so most of our code is independent
   * of whether we're SSE vs. another vector flavor (VMX, AVX, whatever).
   */
  p7_simdvec_Init();

  return eslOK;
}


/* Function:  p7_banner()
 * Synopsis:  print standard HMMER application output header
 * Incept:    SRE, Wed May 23 10:45:53 2007 [Janelia]
 *
 * Purpose:   Print the standard HMMER command line application banner
 *            to <fp>, constructing it from <progname> (the name of the
 *            program) and a short one-line description <banner>.
 *            For example, 
 *            <p7_banner(stdout, "hmmsim", "collect profile HMM score distributions");>
 *            might result in:
 *            
 *            \begin{cchunk}
 *            # hmmsim :: collect profile HMM score distributions
 *            # HMMER 4.0 (Sept 2014); hmmer.org
 *            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *            \end{cchunk}
 *              
 *            <progname> would typically be an application's
 *            <argv[0]>, rather than a fixed string. This allows the
 *            program to be renamed, or called under different names
 *            via symlinks. Any path in the <progname> is discarded;
 *            for instance, if <progname> is "/usr/local/bin/hmmsim",
 *            "hmmsim" is used as the program name.
 *            
 * Note:    
 *    Needs to pick up preprocessor #define's from p7_config.h,
 *    as set by ./configure:
 *            
 *    symbol          example
 *    ------          ----------------
 *    HMMER_VERSION   "4.0"
 *    HMMER_DATE      "Sept 2014"
 *    HMMER_URL       "hmmer.org"
 *
 * Returns:   (void)
 */
void
p7_banner(FILE *fp, char *progname, char *banner)
{
  char *appname = NULL;

  esl_FileTail(progname, FALSE, &appname);  // remove leading directory path

  fprintf(fp, "# %s :: %s\n", appname? appname : progname, banner);
  fprintf(fp, "# HMMER %s (%s); %s\n", HMMER_VERSION, HMMER_DATE, HMMER_URL);
  fprintf(fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  if (appname != NULL) free(appname);
  return;
}


/* Function:  p7_CreateDefaultApp()
 * Synopsis:  Initialize a small/simple/standard HMMER application
 *
 * Purpose:   Identical to <esl_getopts_CreateDefaultApp()>, but 
 *            specialized for HMMER. See documentation in 
 *            <easel/esl_getopts.c>. 
 *
 * Args:      options - array of <ESL_OPTIONS> structures for getopts
 *            nargs   - number of cmd line arguments expected (excl. of cmdname)
 *            argc    - <argc> from main()
 *            argv    - <argv> from main()
 *            banner  - optional one-line description of program (or NULL)
 *            usage   - optional one-line usage hint (or NULL)
 *
 * Returns:   ptr to new <ESL_GETOPTS> object.
 * 
 *            On command line errors, this routine prints an error
 *            message to <stderr> then calls <exit(1)> to halt
 *            execution with abnormal (1) status.
 *            
 *            If the standard <-h> option is seen, the routine prints
 *            the help page (using the data in the <options> structure),
 *            then calls <exit(0)> to exit with normal (0) status.
 *            
 * Xref:      J7/3
 * 
 * Note:      The only difference between this and esl_getopts_CreateDefaultApp()
 *            is to call p7_banner() instead of esl_banner(), to get HMMER
 *            versioning info into the header. There ought to be a better way
 *            (perhaps using PACKAGE_* define's instead of HMMER_* vs. EASEL_*
 *            define's in esl_banner(), thus removing the need for p7_banner).
 */
ESL_GETOPTS *
p7_CreateDefaultApp(ESL_OPTIONS *options, int nargs, int argc, char **argv, char *banner, char *usage)
{
  ESL_GETOPTS *go = NULL;

  p7_Init();

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK) 
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      if (usage != NULL) esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      if (banner != NULL) p7_banner(stdout, argv[0], banner);
      if (usage  != NULL) esl_usage (stdout, argv[0], usage);
      puts("\nOptions:");
      esl_opt_DisplayHelp(stdout, go, 0, 2, 80);
      exit(0);
    }
  if (nargs != -1 && esl_opt_ArgNumber(go) != nargs) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  return go;
}


/* Function:  p7_AminoFrequencies()
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
p7_AminoFrequencies(float *f)
{
  f[0] = 0.0787945;		/* A */
  f[1] = 0.0151600;		/* C */
  f[2] = 0.0535222;		/* D */
  f[3] = 0.0668298;		/* E */
  f[4] = 0.0397062;		/* F */
  f[5] = 0.0695071;		/* G */
  f[6] = 0.0229198;		/* H */
  f[7] = 0.0590092;		/* I */
  f[8] = 0.0594422;		/* K */
  f[9] = 0.0963728;		/* L */
  f[10]= 0.0237718;		/* M */
  f[11]= 0.0414386;		/* N */
  f[12]= 0.0482904;		/* P */
  f[13]= 0.0395639;		/* Q */
  f[14]= 0.0540978;		/* R */
  f[15]= 0.0683364;		/* S */
  f[16]= 0.0540687;		/* T */
  f[17]= 0.0673417;		/* V */
  f[18]= 0.0114135;		/* W */
  f[19]= 0.0304133;		/* Y */
  return eslOK;
}

/*****************************************************************
 * 2. Error handling.
 *****************************************************************/
/* 
 * HMMER's fatal error messages distinguish between user errors
 * ("failure", with p7_Fail()) and internal faults ("death", with
 * p7_Die()). For now, though, there is no difference between the two
 * functions. Someday we might have p7_Die() print a comforting
 * apology, or provide some help on how to report bugs to us;
 * p7_Fail() might provide some pointers on where to read more
 * documentation.
 */

/* Function:  p7_Die()
 * Synopsis:  Handle a fatal exception (something that's our fault)
 */
void
p7_Die(char *format, ...)
{
  va_list  argp;
#ifdef HAVE_MPI
  int      mpiflag;
#endif
                                /* format the error mesg */
  fprintf(stderr, "\nFATAL: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
#ifdef HAVE_MPI
  MPI_Initialized(&mpiflag);
  if (mpiflag) MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  exit(1);
}

/* Function:  p7_Fail()
 * Synopsis:  Handle a user error (something that's the user's fault).
 */
void
p7_Fail(char *format, ...)
{
  va_list  argp;

  /* Check whether we are running as a daemon so we can do the right thing about logging instead of printing errors */
  int parent_pid;
  parent_pid = getppid();


  if(parent_pid != 1){ // We are not running as a daemon, so just print the error message
                                /* format the error mesg */
    fprintf(stderr, "\nError: ");
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fflush(stderr);
    exit(1);
  }
  else{ // I am running as a daemon, so log the error to /sys/log
    vsyslog(LOG_ERR, format, argp);
    exit(1);
  }

}




