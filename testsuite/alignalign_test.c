/* alignalign_test.c
 * Sun Jul  5 13:42:41 1998
 * 
 * Test driver for P7ViterbiAlignAlignment().
 * 
 * The test is to 
 *    1) read an alignment and a corresponding HMM
 *    2) align the alignment to the HMM to get a master trace
 *    3) map the alignment to the HMM to get another master trace
 *    4) Test that the two master traces are identical; if not, fail.
 *         This doesn't have to be true always, but it's true for the
 *         fn3 test example.
 *    5) Get imposed traces for each sequence
 *    6) Viterbi align individual seqs to the model;
 *           compare the imposed trace with the Viterbi trace;
 *    7) If an excessive number of individual traces differ from
 *          those imposed by master, fail.
 * 
 * RCS $Id$
 */

#include <stdio.h>

#include "structs.h"
#include "funcs.h"
#include "globals.h"
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static char banner[] = "\
alignalign_test : testing of P7ViterbiAlignAlignment() code";

static char usage[] = "\
Usage: alignalign_test [-options]\n\
  Available options are:\n\
  -h              : help; display this usage info\n\
  -o <f>          : output alignment to <f> for comparing diffs\n\
  -v              : be verbose\n\
";

static char experts[] = "\
  --ali <f>       : read alignment from <f>\n\
  --hmm <f>       : read HMM from <f>\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-h",       TRUE,  sqdARG_NONE },
  { "-o",       TRUE,  sqdARG_NONE },
  { "-v",       TRUE,  sqdARG_NONE },
  { "--ali",    FALSE, sqdARG_STRING },
  { "--hmm",    FALSE, sqdARG_STRING },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char    *hmmfile;	        /* file to read HMM(s) from                */
  HMMFILE *hmmfp;               /* opened hmmfile for reading              */
  struct plan7_s  *hmm;         /* HMM to search with                      */ 
  char    *afile;               /* file to read alignment from             */
  int      format;              /* format determined for afile             */
  char   **aseq;                /* aligned sequences                       */
  char   **rseq;                /* raw, dealigned aseq                     */
  AINFO    ainfo;               /* alignment information                   */
  char     *dsq;		/* digitized target sequence               */
  struct p7trace_s  *mtr;	/* master traceback from alignment         */
  struct p7trace_s  *maptr;     /* master traceback from mapping           */
  struct p7trace_s **tr;        /* individual tracebacks imposed by mtr    */
  struct p7trace_s **itr;       /* individual trace from P7Viterbi()       */
  int       idx;		/* counter for seqs                        */
  int       ndiff;		/* number of differing traces              */
  
  char *ofile;
  FILE *ofp;
  struct p7trace_s **otr;       /* traces for output alignment             */
  char **odsq;                  /* digitized sequences for output          */
  SQINFO *osqinfo;              /* sqinfo for output                       */
  float  *owgt;                 /* wgt for output                          */
  int     oi;			/* counter for seqs in output              */
  char  **oaseq;
  AINFO   oainfo;

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
  hmmfile    = "fn3.hmm";
  afile      = "fn3.seed";
  ofile      = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-o")       == 0) ofile      = optarg; 
    else if (strcmp(optname, "-v")       == 0) be_verbose = TRUE;
    else if (strcmp(optname, "--ali")    == 0) afile      = optarg;
    else if (strcmp(optname, "--hmm")    == 0) hmmfile    = optarg;
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
   * Open test alignment file
   ***********************************************/

  if (! SeqfileFormat(afile, &format, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: 
      Die("Alignment file %s could not be opened for reading", afile); break;
    case SQERR_FORMAT: 
    default:           
      Die("Failed to determine format of alignment file %s", afile);
    }

  if (! ReadAlignment(afile, format, &aseq, &ainfo))
    Die("Failed to read aligned sequence file %s", afile);
  for (idx = 0; idx < ainfo.nseq; idx++)
    s2upper(aseq[idx]);
  DealignAseqs(aseq, ainfo.nseq, &rseq);


  /*********************************************** 
   * Open HMM file 
   * Read a single HMM from it. (Config HMM, if necessary).
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, NULL)) == NULL)
    Die("Failed to open HMM file %s\n", hmmfile);
  if (!HMMFileRead(hmmfp, &hmm)) 
    Die("Failed to read any HMMs from %s\n", hmmfile);
  if (hmm == NULL) 
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  P7Logoddsify(hmm, TRUE);

  if (! (hmm->flags & PLAN7_MAP))
    Die("HMM in %s has no map", hmmfile);
  if (GCGMultchecksum(aseq, ainfo.nseq) != hmm->checksum)
    Die("Checksum for alignment in %s does not match that in HMM (%d != %d)", 
	afile, GCGMultchecksum(aseq, ainfo.nseq), hmm->checksum);

  /* Allocations for output alignment
   */
  otr  = MallocOrDie(sizeof(struct p7trace_s *) * ainfo.nseq * 2);
  odsq = MallocOrDie(sizeof(char *) * ainfo.nseq * 2);
  osqinfo = MallocOrDie(sizeof(SQINFO) * ainfo.nseq * 2);
  owgt = MallocOrDie(sizeof(float) * ainfo.nseq*2);
  FSet(owgt, ainfo.nseq*2, 1.0);

  /*********************************************** 
   * Align alignment to HMM, get master trace
   ***********************************************/

  mtr   = P7ViterbiAlignAlignment(aseq, &ainfo, hmm);
  maptr = MasterTraceFromMap(hmm->map, hmm->M, ainfo.alen);
  itr = MallocOrDie(sizeof(struct p7trace_s *) * ainfo.nseq);

				/* verify master traces */
  if (! TraceVerify(mtr, hmm->M, ainfo.alen))
    Die("Trace verify on P7ViterbiAlignAlignment() result failed\n");
  if (! TraceVerify(maptr, hmm->M, ainfo.alen))
    Die("Trace verify on MasterTraceFromMap() result failed\n");
  if (! TraceCompare(mtr, maptr))
    Die("Master traces differ for alignment versus map\n");

				/* impose master trace on individuals */
  ImposeMasterTrace(aseq, ainfo.nseq, mtr, &tr);

				/* align individuals, compare traces */
  oi = 0;
  ndiff = 0;
  for (idx = 0; idx < ainfo.nseq; idx++)
    {
      dsq = DigitizeSequence(rseq[idx], ainfo.sqinfo[idx].len);
      P7Viterbi(dsq, ainfo.sqinfo[idx].len, hmm, &(itr[idx]));
      
      /* The trace given by the alignment automatically goes in
       * the output alignment
       */
      odsq[oi] = DigitizeSequence(rseq[idx], ainfo.sqinfo[idx].len);
      SeqinfoCopy(&(osqinfo[oi]), &(ainfo.sqinfo[idx]));
      otr[oi] = tr[idx];
      oi++;

      if (! TraceCompare(itr[idx], tr[idx]))
	{
				/* put indiv Viterbi trace in output for comparison */
	  odsq[oi] = DigitizeSequence(rseq[idx], ainfo.sqinfo[idx].len);
	  SeqinfoCopy(&(osqinfo[oi]), &(ainfo.sqinfo[idx]));
	  otr[oi] = itr[idx];
	  oi++;
	  ndiff++;
	}
      
      free(dsq);
    }

  /* Output the alignment
   */
  P7Traces2Alignment(odsq, osqinfo, owgt, oi, hmm->M, otr, FALSE,
		     &oaseq, &oainfo);
  if (ofile != NULL) 
    {
      if ((ofp = fopen(ofile,"w")) == NULL)
	Die("alignalign_test: failed to open output file %s", ofile);
      WriteSELEX(stdout, oaseq, &oainfo, 50);
      fclose(ofp);
    }

  /* Determine success/failure. Leave as flag in ndiff.
   */
  if (ndiff > ainfo.nseq / 2) {
    if (be_verbose) printf("alignalign: Test FAILED; %d/%d differ\n", ndiff, ainfo.nseq);
    ndiff = EXIT_FAILURE;
  } else {
    if (be_verbose) printf("alignalign: Test passed; %d/%d differ\n", ndiff, ainfo.nseq);
    ndiff = EXIT_SUCCESS;
  }

  /* Cleanup.
   */
  P7FreeTrace(mtr);
  P7FreeTrace(maptr);
  for (idx = 0; idx < ainfo.nseq; idx++)
    {
      P7FreeTrace(tr[idx]);
      P7FreeTrace(itr[idx]);
    }
  for (idx = 0; idx < oi; idx++)
    free(odsq[idx]);
  free(tr);
  free(itr);
  free(otr);
  free(odsq);
  free(owgt);
  free(osqinfo);
  Free2DArray(rseq, ainfo.nseq);
  FreeAlignment(aseq, &ainfo);
  FreeAlignment(oaseq, &oainfo);
  FreePlan7(hmm);
  SqdClean();

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) Die("trace_test failed memory test");
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  
  return ndiff;
}
