/* viterbi_exercise.c
 * SRE, Mon Mar  9 07:55:47 1998 [St. Louis]
 * 
 * Exercise the various Viterbi algorithms, big and small.
 * 
 * CVS $Id$
 */


#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "structs.h"
#include "funcs.h"
#include "globals.h"
#include "squid.h"

static char banner[] = "\
viterbi_exercise : testing of Plan7 Viterbi code";

static char usage[] = "\
Usage: testdriver [-options]\n\
  Available options are:\n\
  -h              : help; display this usage info\n\
  -v              : be verbose\n\
";

static char experts[] = "\
  --hmm <f>       : use HMM in file <f>\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-h",       TRUE,  sqdARG_NONE },
  { "-v",       TRUE,  sqdARG_NONE },
  { "--hmm",    FALSE, sqdARG_STRING },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char    *hmmfile;	        /* file to read HMM(s) from                */
  HMMFILE *hmmfp;               /* opened hmmfile for reading              */
  struct plan7_s  *hmm;         /* the HMM to search with                  */ 
  char    *dsq;			/* digitized target sequence               */
  char    *seq;
  SQINFO   sqinfo;
  int      L;			/* length of dsq                           */
  struct dpmatrix_s *mx;        /* growable, reusable DP matrix            */
  struct p7trace_s  *tr1;	/* traceback                               */
  struct p7trace_s  *tr2;	/* another traceback                       */
  int       nseq;
  float     sc1, sc2;		/* scores                                  */
  int       config;
  int       i;

  int be_verbose;

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */

  
  /*********************************************** 
   * Parse command line
   ***********************************************/

  be_verbose = FALSE;
  hmmfile    = "fn3.hmm";
  nseq       = 100;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-v")       == 0) be_verbose = TRUE;
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
   * Open HMM file 
   * Read a single HMM from it. 
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, NULL)) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  if (!HMMFileRead(hmmfp, &hmm)) 
    Die("Failed to read any HMMs from %s\n", hmmfile);
  if (hmm == NULL) 
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  Plan7Renormalize(hmm);

  /*********************************************** 
   * We cycle through different model configurations.
   * For each configuration, we repeat 100 times:
   *    - generate a sequence
   *    - score it by Viterbi and by SmallViterbi  
   *    - make sure they give OK and identical results 
   ***********************************************/

  for (config = 1; config <= 5; config++)
    {
      switch (config) {
      case 1: Plan7NakedConfig(hmm);            break;
      case 2: Plan7GlobalConfig(hmm);           break;
      case 3: Plan7LSConfig(hmm);               break;
      case 4: Plan7FSConfig(hmm, 0.5, 0.5);     break;
      case 5: Plan7SWConfig(hmm, 0.5, 0.5);     break;
      default: Die("never happens");
      }
      P7Logoddsify(hmm, TRUE);

      
      mx = CreatePlan7Matrix(1, hmm->M, 25, 0);
      for (i = 0; i < nseq; i++)
	{
	  EmitSequence(hmm, &dsq, &L, NULL);
	  sprintf(sqinfo.name, "seq%d", i+1);
	  sqinfo.len   = L;
	  sqinfo.flags = SQINFO_NAME | SQINFO_LEN;

	  sc1 = P7Viterbi(dsq, L, hmm, mx, &tr1);
	  sc2 = P7SmallViterbi(dsq, L, hmm, mx, &tr2);

	  if (be_verbose)
	    {
	      printf("Viterbi score: %.1f   SmallViterbi: %.1f\n", sc1, sc2);
	      P7PrintTrace(stdout, tr1, hmm, dsq);
	      P7PrintTrace(stdout, tr2, hmm, dsq);
	      
	      seq = DedigitizeSequence(dsq, L);
	      WriteSeq(stdout, SQFILE_FASTA, seq, &sqinfo);
	      free(seq);
	    }

	  if (sc1 != sc2) 
	    Die("Different scores from normal/small Viterbi");

	  if (fabs(sc1 - P7TraceScore(hmm, dsq, tr1)) > 0.1)
	    Die("P7Viterbi score doesn't match its TraceScore");
	  if (fabs(sc2 - P7TraceScore(hmm, dsq, tr2)) > 0.1)
	    Die("P7SmallViterbi score doesn't match its TraceScore");

	  if (! TraceVerify(tr1, hmm->M, L))
	    Die("TraceVerify() failed for a P7Viterbi trace");
	  if (! TraceVerify(tr2, hmm->M, L))
	    Die("TraceVerify() failed for a P7SmallViterbi trace");

	  if (tr1->tlen != tr2->tlen)
	    Die("Trace lengths differ for normal/small Viterbi");
	  if (! TraceCompare(tr1, tr2))
	    Die("Different traces from normal/small Viterbi");

	  P7FreeTrace(tr1);
	  P7FreeTrace(tr2);
	  free(dsq);
	}
    }

  FreePlan7Matrix(mx);
  return EXIT_SUCCESS;
}
