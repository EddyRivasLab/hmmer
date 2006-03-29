/* parsingviterbi_test.c
 * Wed Mar  4 15:07:37 1998
 * cp trace_test.c ../src/testdriver.c; cd ../src; make testdriver
 * 
 * Test driver for P7ParsingViterbi(); alignment in linear memory.
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
parsingviterbi_test : testing of Plan7 linear memory alignment code";

static char usage[] = "\
Usage: parsingviterbi_test [-options]\n\
  Available options are:\n\
  -h              : help; display this usage info\n\
  -v              : be verbose\n\
";

static char experts[] = "\
\n";

static struct opt_s OPTIONS[] = {
  { "-h",       TRUE,  sqdARG_NONE },
  { "-v",       TRUE,  sqdARG_NONE },
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
  struct p7trace_s  *tr1;	/* traceback from P7Viterbi()              */
  struct p7trace_s  *tr2;	/* traceback from P7ParsingViterbi()       */
  int       nseq;
  float     sc1, sc2;		/* scores from Viterbi, ParsingViterbi()   */

  struct p7trace_s **tarr;	/* array of decomposed Viterbi traces      */
  int    ntr;			/* number of traces                        */
  int    i1,i2,k1,k2;		/* starts, stops in seq, model for Viterbi */
  int    idx;			/* index of a decomposed trace             */

  int be_verbose;

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

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-v")       == 0) be_verbose = TRUE;
    else if (strcmp(optname, "-h")       == 0) {
      Banner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(0);
    }
  }
  if (argc - optind != 0)
    Die("Incorrect number of arguments.\n%s\n", usage);

  hmmfile = "fn3.hmm";
  seqfile = "titin.fa";

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
   * Search HMM against each sequence, using both
   * normal Viterbi and P7ParsingViterbi.
   ***********************************************/

  nseq = 0;
  while (ReadSeq(sqfp, format, &seq, &sqinfo)) 
    {
      nseq++;
      dsq = DigitizeSequence(seq, sqinfo.len);

      sc1  = P7Viterbi(dsq, sqinfo.len, hmm, &tr1);
      sc2  = P7ParsingViterbi(dsq, sqinfo.len, hmm, &tr2);

      if (be_verbose)
	{
	  printf("test sequence %d: %s %s\n",
		 nseq, sqinfo.name, 
		 sqinfo.flags & SQINFO_DESC ? sqinfo.desc : "");
	  for (idx = 0; idx < tr2->tlen; idx++)
	    printf("%1s  %d\n", Statetype(tr2->statetype[idx]), tr2->pos[idx]);
	}

      if (sc1 != sc2)
	Die("Scores for the two Viterbi implementations are unequal (%d,%d)", sc1, sc2);

      TraceDecompose(tr1, &tarr, &ntr);
      if (ntr == 0)
	Die("ntr == 0 can't happen");
      if (ntr != (tr2->tlen/2) -1)
	Die("# of domains for the two Viterbi implementations are unequal (%d, %d)",
	    ntr, (tr2->tlen/2) -1); 
      
      for (idx = 0; idx < ntr; idx++)
	{
	  TraceSimpleBounds(tarr[idx], &i1, &i2, &k1, &k2);
	  
	  if (i1 != tr2->pos[idx*2 + 1] + 1)
	    Die("Start positions %d and %d disagree for domain %d\n", 
		i1, tr2->pos[idx*2 + 1] + 1, idx);
	  if (i2 != tr2->pos[idx*2 + 2])
	    Die("End positions %d and %d disagree for domain %d\n", 
		i2, tr2->pos[idx*2 + 2], idx);
	}


      for (idx = 0; idx < ntr; idx++)
	P7FreeTrace(tarr[idx]);
      free(tarr);
      FreeSequence(seq, &sqinfo); 
      P7FreeTrace(tr1);
      P7FreeTrace(tr2);
      free(dsq);
    }

  FreePlan7(hmm);
  HMMFileClose(hmmfp);
  SeqfileClose(sqfp);

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) Die("evd_test failed memory test");
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  
  return EXIT_SUCCESS;
}
