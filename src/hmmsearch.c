/* hmmsearch.c
 * SRE, Tue Jan  7 17:19:20 1997
 *
 * main() for HMM global alignment searching.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "structs.h"		/* data structures, macros, #define's   */
#include "config.h"		/* compile-time configuration constants */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "squid.h"		/* general sequence analysis library    */

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static char banner[] = "hmms7 - HMM search: global alignment";

static char usage[]  = "\
Usage: hmms [-options] <hmmfile> <sequence file or database>\n\
  Available options are:\n\
   -h        : help; print brief help on version and usage\n\
   -t <x>    : report only hits above a score of <x> [default 0]\n\
   --forward : use the full Forward() algorithm instead of Viterbi\n\
\n";


static struct opt_s OPTIONS[] = {
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-t",        TRUE,  sqdARG_FLOAT },
  { "--forward", FALSE, sqdARG_NONE},
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv) 
{
  char            *hmmfile;	/* file to read HMM(s) from                */
  HMMFILE *hmmfp;       /* opened hmmfile for reading              */
  char    *seqfile;     /* file to read target sequence(s) from    */ 
  SQFILE   *sqfp;        /* opened seqfile for reading              */
  int       format;	/* format of seqfile                       */
  char     *seq;		/* target sequence                         */
  SQINFO    sqinfo;	/* optional info for seq                   */
  char     *dsq;		/* digitized target sequence               */
  float     sc;		/* score of an HMM search                  */
  int       validfit;         /* TRUE if histogram gave us an EVD fit     */
  struct plan7_s  *hmm;         /* HMM to search with                      */ 
  struct histogram_s *histogram;/* histogram of all scores              */
  struct dpmatrix_s *mx;	/* DP matrix after alignment */
  struct p7trace_s  *tr;	/* traceback */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */
  float cutoff;			/* report scores above this                 */
  int   do_forward;		/* TRUE to use Forward() not Viterbi()      */
  int   working;

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  cutoff     = 0.;
  do_forward = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-t")        == 0) cutoff     = atof(optarg);
    else if (strcmp(optname, "--forward") == 0) do_forward = TRUE;
    else if (strcmp(optname, "-h") == 0) {
      Banner(stdout, banner);
      puts(usage);
      exit(0);
    }
  }
  if (argc - optind != 2)
    Die("Incorrect number of arguments.\n%s\n", usage);

  hmmfile = argv[optind++];
  seqfile = argv[optind++]; 

  /*********************************************** 
   * Open sequence database (might be in BLASTDB or current directory)
   ***********************************************/

  if (! SeqfileFormat(seqfile, &format, "BLASTDB"))
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
   * Open HMM file (might be in HMMERDB or current directory)
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  
  /*********************************************** 
   * Allocate histogram structure
   ***********************************************/

  histogram = AllocHistogram(-200, 200, 100);

  /*********************************************** 
   * Show the banner
   ***********************************************/

  Banner(stdout, banner);
  printf(   "HMM file:                 %s\n", hmmfile);
  printf(   "Sequence database:        %s\n", seqfile);
  printf(   "Search strategy:          global alignment\n");
  printf(   "Cutoff at score:          %.2f\n", cutoff);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  /*********************************************** 
   * Search each HMM against each sequence
   ***********************************************/

  while (HMMFileRead(hmmfp, &hmm)) {
    if (hmm == NULL) 
      Die("HMM file %s may be corrupt or in incorrect format; parse failed");
    Plan7LSConfig(hmm); 
    /*Plan7SWConfig(hmm, 0.001, 0.001); */
    Plan7Logoddsify(hmm);

    working = 1;
#pragma parallel local(seq, sqinfo, dsq, sc)
#pragma shared(sqfp, format, hmm, cutoff, histogram, do_forward, working)
{
    while (working) {
#pragma critical
      if (! ReadSeq(sqfp, format, &seq, &sqinfo)) { working = 0; break; }

      s2upper(seq);
      dsq = DigitizeSequence(seq, sqinfo.len);
      
      if (do_forward)
	sc = Plan7Forward(dsq, sqinfo.len, hmm, NULL);
      else
	{
	  sc = Plan7Viterbi(dsq, sqinfo.len, hmm, &mx);
	  P7ViterbiTrace(hmm, dsq, sqinfo.len, mx, &tr);
	  /*	  P7PrintTrace(stdout, tr, hmm, dsq); */
	  P7FreeTrace(tr);
	  FreePlan7Matrix(mx);
	}

#pragma critical
      if (sc > cutoff)
	printf("%-7.1f %12s %s\n", sc, sqinfo.name,
	       sqinfo.flags & SQINFO_DESC ? sqinfo.desc : "");

#pragma critical
      AddToHistogram(histogram, sc);

      free(dsq);
      FreeSequence(seq, &sqinfo); 
    }
}

    FreePlan7(hmm);
    SeqfileRewind(sqfp);

  }


  
  /**************************************************
   * Do extreme value distribution black magic to find P, E values
   **************************************************/

  validfit = ExtremeValueFitHistogram(histogram, 20.);
  PrintASCIIHistogram(stdout, histogram);

  /*********************************************** 
   * Clean-up and exit.
   ***********************************************/

  FreeHistogram(histogram);
  HMMFileClose(hmmfp);
  SeqfileClose(sqfp);
  SqdClean();

#ifdef MEMDEBUG
  current_size = malloc_size(&histid2);
  if (current_size != orig_size) malloc_list(2, histid1, histid2);
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  return 0;
}


