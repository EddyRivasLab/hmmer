/* trace_test.c
 * Test driver for Viterbi tracebacks.
 * 
 * Mon Feb  2 07:57:47 1998
 * CVS $Id$
 */


#include "config.h"

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "squid.h"

#include "plan7.h"
#include "structs.h"
#include "funcs.h"
#include "globals.h"


static char banner[] = "\
trace_test : testing of Plan7 Viterbi traceback code";

static char usage[] = "\
Usage: trace_test [-options]\n\
  Available options are:\n\
  -h              : help; display this usage info\n\
  -v              : be verbose\n\
";

static char experts[] = "\
  --hmm <f>       : use HMM in file <f>\n\
  --seq <f>       : use seq(s) in file <f>\n\
  --small         : run P7SmallViterbi()\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-h",       TRUE,  sqdARG_NONE },
  { "-v",       TRUE,  sqdARG_NONE },
  { "--hmm",    FALSE, sqdARG_STRING },
  { "--seq",    FALSE, sqdARG_STRING },
  { "--small",  FALSE, sqdARG_NONE },
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
  unsigned char   *dsq;		/* digitized target sequence               */
  struct plan7_s  *hmm;         /* HMM to search with                      */ 
  cust_dpmatrix_s *mx;          /* reusable, growable DP matrix            */
  struct p7trace_s  *tr;	/* traceback                               */
  int       nseq;
  float     sc;

  int be_verbose;
  int do_small;			/* TRUE to invoke P7SmallViterbi */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  be_verbose = FALSE;
  hmmfile    = "trace_test.hmm";
  seqfile    = "trace_test.seq";
  do_small   = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-v")       == 0) be_verbose = TRUE;
    else if (strcmp(optname, "--hmm")    == 0) hmmfile    = optarg;
    else if (strcmp(optname, "--seq")    == 0) seqfile    = optarg;
    else if (strcmp(optname, "--small")  == 0) do_small   = TRUE;
    else if (strcmp(optname, "-h")       == 0) {
      HMMERBanner(stdout, banner);
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
   * Read a single HMM from it.
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, NULL)) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  if (!HMMFileRead(hmmfp, &hmm)) 
    Die("Failed to read any HMMs from %s\n", hmmfile);
  if (hmm == NULL) 
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);

  /*********************************************** 
   * Search HMM against each sequence
   ***********************************************/

  nseq = 0;
  mx = CreateDPMatrix(1, hmm->M, 25, 0);
  while (ReadSeq(sqfp, sqfp->format, &seq, &sqinfo)) 
    {
      nseq++;
      dsq = DigitizeSequence(seq, sqinfo.len);

      if (do_small) sc = P7SmallViterbi(dsq, sqinfo.len, hmm, mx, &tr);
      else          sc = Viterbi(dsq, sqinfo.len, hmm, mx, &tr);

      if (be_verbose)
	{
	  printf("test sequence %d: score %.1f : %s %s\n",
		 nseq, sc, sqinfo.name, 
		 sqinfo.flags & SQINFO_DESC ? sqinfo.desc : "");
	  P7PrintTrace(stdout, tr, hmm, dsq); 
	}

      if (! TraceVerify(tr, hmm->M, sqinfo.len))
	Die("Trace verify failed on seq #%d, %s\n", nseq, sqinfo.name);

      FreeSequence(seq, &sqinfo); 
      P7FreeTrace(tr);
      free(dsq);
    }

  FreeDPMatrix(mx);
  FreePlan7(hmm);
  HMMFileClose(hmmfp);
  SeqfileClose(sqfp);

  return EXIT_SUCCESS;
}
