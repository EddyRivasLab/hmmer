/* trace_test.c
 * Mon Feb  2 07:57:47 1998
 * cp trace_test.c ../src/testdriver.c; cd ../src; make testdriver
 * 
 * Test driver for Viterbi tracebacks.
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

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static char banner[] = "\
trace_test : testing of Plan7 Viterbi traceback code";

static char usage[] = "\
Usage: testdriver [-options]\n\
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
  int       format;	        /* format of seqfile                       */
  char     *seq;		/* target sequence                         */
  SQINFO    sqinfo;	        /* optional info for seq                   */
  char     *dsq;		/* digitized target sequence               */
  struct plan7_s  *hmm;         /* HMM to search with                      */ 
  struct p7trace_s  *tr;	/* traceback                               */
  int       nseq;
  float     sc;

  int be_verbose;
  int do_small;			/* TRUE to invoke P7SmallViterbi */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif
  
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

  if (! SeqfileFormat(seqfile, &format, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: 
      Die("Sequence file %s could not be opened for reading", seqfile); break;
    case SQERR_FORMAT: 
    default:           
      Die("Failed to determine format of sequence file %s", seqfile);
    }
  if ((sqfp = SeqfileOpen(seqfile, format, "BLASTDB")) == NULL)
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
  while (ReadSeq(sqfp, format, &seq, &sqinfo)) 
    {
      nseq++;
      dsq = DigitizeSequence(seq, sqinfo.len);

      if (do_small) sc = P7SmallViterbi(dsq, sqinfo.len, hmm, &tr);
      else          sc = P7Viterbi(dsq, sqinfo.len, hmm, &tr);

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

  FreePlan7(hmm);
  HMMFileClose(hmmfp);
  SeqfileClose(sqfp);

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) Die("trace_test failed memory test");
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  
  return EXIT_SUCCESS;
}
