/* weeviterbi_test.c
 * Wed Mar  4 17:30:39 1998
 * 
 * Test driver for Myers/Miller/Hirschberg linear memory Viterbi tracebacks.
 * 
 * RCS $Id$
 */

#include <stdio.h>
#include <time.h>
#include <math.h>

#include "structs.h"
#include "funcs.h"
#include "globals.h"
#include "squid.h"

static char banner[] = "\
weeviterbi_test : testing of Plan7 Myers/Miller/Hirschberg Viterbi traceback code";

static char usage[] = "\
Usage: testdriver [-options]\n\
  Available options are:\n\
  -h              : help; display this usage info\n\
  -v              : be verbose\n\
";

static char experts[] = "\
  --hmm <f>       : use HMM in file <f>\n\
  --seq <f>       : use seq(s) in file <f>\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-h",       TRUE,  sqdARG_NONE },
  { "-v",       TRUE,  sqdARG_NONE },
  { "--hmm",    FALSE, sqdARG_STRING },
  { "--seq",    FALSE, sqdARG_STRING },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char    *hmmfile;	        /* file to read HMM(s) from                */
  HMMFILE *hmmfp;               /* opened hmmfile for reading              */
  char    *seqfile;             /* file to read target sequence(s) from    */ 
  SQFILE   *sqfp;               /* opened seqfile for reading              */
  char     *seq;		/* target sequence                         */
  SQINFO    sqinfo;	        /* optional info for seq                   */
  char     *dsq;		/* digitized target sequence               */
  struct plan7_s  *hmm;         /* HMM to search with                      */ 
  struct p7trace_s  *t1;	/* standard Viterbi traceback              */
  struct p7trace_s  *t2;	/* WeeViterbi traceback                    */
  int       nseq;
  float     sc1,sc2;		/* scores from Viterbi, WeeViterbi         */

  int be_verbose;

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  be_verbose = FALSE;
  hmmfile    = "weeviterbi_test.hmm";
  seqfile    = "weeviterbi_test.seq";

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-v")       == 0) be_verbose = TRUE;
    else if (strcmp(optname, "--hmm")    == 0) hmmfile    = optarg;
    else if (strcmp(optname, "--seq")    == 0) seqfile    = optarg;
    else if (strcmp(optname, "-h")       == 0) {
      Banner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(0);
    }
  }
  if (argc - optind != 0)
    Die("Incorrect number of arguments.\n%s\n", usage);

  /*********************************************** 
   * Open test sequence file
   ***********************************************/

  if ((sqfp = SeqfileOpen(seqfile, SQFILE_UNKNOWN, "BLASTDB")) == NULL)
    Die("Failed to open sequence database file %s\n%s\n", seqfile, usage);

  /*********************************************** 
   * Open HMM file 
   * Read a single HMM from it. (Config HMM, if necessary).
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, NULL)) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  if (!HMMFileRead(hmmfp, &hmm)) 
    Die("Failed to read any HMMs from %s\n", hmmfile);
  if (hmm == NULL) 
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  P7Logoddsify(hmm, TRUE);

  /*********************************************** 
   * Search HMM against each sequence
   ***********************************************/

  nseq = 0;
  while (ReadSeq(sqfp, sqfp->format, &seq, &sqinfo)) 
    {
      nseq++;
      dsq = DigitizeSequence(seq, sqinfo.len);

      sc1 = P7Viterbi(dsq, sqinfo.len, hmm, &t1);
      sc2 = P7WeeViterbi(dsq, sqinfo.len, hmm, &t2);

      if (be_verbose)
	{
	  printf("test sequence %d: %s %s\n",
		 nseq, sqinfo.name, 
		 sqinfo.flags & SQINFO_DESC ? sqinfo.desc : "");
	  printf("** P7Viterbi trace:\n");
	  P7PrintTrace(stdout, t1, hmm, dsq); 
	  printf("** P7WeeViterbi trace:\n");
	  P7PrintTrace(stdout, t2, hmm, dsq); 
	}

      if (! TraceVerify(t1, hmm->M, sqinfo.len))
	Die("Trace verify failed on Viterbi for seq #%d, %s\n", nseq, sqinfo.name);
      if (! TraceVerify(t2, hmm->M, sqinfo.len))
	Die("Trace verify failed on WeeViterbi for seq #%d, %s\n", nseq, sqinfo.name);
      if (sc1 != sc2)
	Die("Scores for the two Viterbi implementations are unequal (%.1f,%.1f)", sc1, sc2);
      if (! TraceCompare(t1, t2))
	Die("WeeViterbi() trace is not identical to Viterbi() trace");

      FreeSequence(seq, &sqinfo); 
      P7FreeTrace(t1);
      P7FreeTrace(t2);
      free(dsq);
    }

  FreePlan7(hmm);
  SeqfileClose(sqfp);
  HMMFileClose(hmmfp);

  return EXIT_SUCCESS;
}
