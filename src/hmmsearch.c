/* hmmsearch.c
 * SRE, Tue Jan  7 17:19:20 1997
 *
 * main() for HMM global alignment searching.
 * RCS $Header$
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

static char banner[] = "hmmsearch - search a sequence database with a profile HMM";

static char usage[]  = "\
Usage: hmmsearch [-options] <hmmfile> <sequence file or database>\n\
  Available options are:\n\
   -h        : help; print brief help on version and usage\n\
   --forward : use the full Forward() algorithm instead of Viterbi\n\
\n";


static struct opt_s OPTIONS[] = {
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "--forward", FALSE, sqdARG_NONE},
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


/* Function: record_by_domain()
 * Date:     SRE, Tue Nov  4 11:25:14 1997 [St. Louis]
 * 
 * Purpose:  Decompose a trace, and register scores, P-values, alignments,
 *           etc. for individual domains in a hitlist. 
 *           
 * Return:   (void)          
 */
static void
record_by_domain(struct tophit_s *h, 
		 struct plan7_s *hmm, char *dsq, SQINFO *sqinfo,
		 struct p7trace_s *tr, float domT)
{
  struct p7trace_s **tarr;      /* array of per-domain traces */
  struct fancyali_s *ali;       /* alignment of a domain      */ 
  int ntr;			/* number of domain traces    */
  int idx;			/* index for traces           */
  int k1, k2;			/* start, stop coord in model */
  int i1, i2;			/* start, stop in sequence    */
  float  score;
  double pvalue;

  TraceDecompose(tr, &tarr, &ntr);
  if (ntr == 0) return;		/* "can't happen" */

  for (idx = 0; idx < ntr; idx++)
    {
      /* Get the score and bounds of the match.
       */
      score  = P7TraceScore(hmm, dsq, tarr[idx]);
      pvalue = PValue(hmm, score); 
      TraceSimpleBounds(tarr[idx], &i1, &i2, &k1, &k2);
      
      /* Record the match. 
       * since we don't yet necessarily know our EVD statistics, we can't
       * screen by E-value.
       * Note that we store the name, desc of the /sequence/.
       * Also, we record the negative pvalue as a sort key.
       */
      if (score > domT)
	{
	  ali = CreateFancyAli(tarr[idx], hmm, dsq, sqinfo->name);
	  RegisterHit(h, -1.*pvalue, pvalue, score,
		      sqinfo->name, 
		      sqinfo->flags & SQINFO_DESC ? sqinfo->desc : NULL, 
		      i1,i2, sqinfo->len, 
		      k1,k2, hmm->M, 
		      idx+1, ntr,
		      ali);
	}
    }

  for (idx = 0; idx < ntr; idx++)
    P7FreeTrace(tarr[idx]);
  free(tarr);
  return;
}


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
  float     sc;	        	/* score of an HMM search                  */
  int       validfit;           /* TRUE if histogram gave us an EVD fit    */
  int       i; 
  struct plan7_s  *hmm;         /* HMM to search with                      */ 
  struct histogram_s *histogram;/* histogram of all scores                 */
  struct dpmatrix_s *mx;	/* DP matrix after alignment               */
  struct p7trace_s  *tr;	/* traceback                               */
  struct tophit_s   *ghit;      /* list of top hits for whole sequences    */
  struct tophit_s   *dhit;	/* list of top hits for domains            */

  double  pvalue;		/* pvalue of an HMM score                  */
  double  evalue;		/* evalue of an HMM score                  */
  int     sqfrom, sqto;		/* coordinates in sequence                 */
  int     hmmfrom, hmmto;	/* coordinate in HMM                       */
  char   *name, *desc;          /* hit sequence name and description       */
  int     sqlen;		/* length of seq that was hit              */
  int     nseq;			/* number of sequences searched            */
  int     domidx;		/* number of this domain                   */
  int     ndom;			/* total # of domains in this seq          */
  int     namewidth;		/* max width of sequence name              */

  float  globT;			/* T parameter: reporting threshold in bits */
  double globE;			/* E parameter: reporting thresh in E-value */
  int    globH;			/* H parameter: save at most top H hits     */
  int    globA;			/* A parameter: save at most top A hits     */

  float  domT;			/* T parameter for individual domains */
  double domE;			/* E parameter for individual domains */
  int    domH;			/* H parameter for individual domains */
  int    domA;			/* A parameter for individual domains */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */
  int   do_forward;		/* TRUE to use Forward() not Viterbi()      */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  do_forward  = FALSE;
  globT       = -999999;
  globE       = 10.0;
  globH       = 10000;
  globA       = 100;
  domT        = -999999;
  domE        = 0.01;
  domH        = 10000;
  domA        = 100;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "--forward") == 0) do_forward = TRUE;
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
   * Open HMM file (might be in HMMERDB or current directory).
   * Read a single HMM from it. (Config HMM, if necessary).
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  if (!HMMFileRead(hmmfp, &hmm)) 
    Die("Failed to read any HMMs from %s\n", hmmfile);
  if (hmm == NULL) 
    Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);
  Plan7Logoddsify(hmm);
  
  /*********************************************** 
   * Show the banner
   ***********************************************/

  Banner(stdout, banner);
  printf(   "HMM file:                 %s [%s]\n", hmmfile, hmm->name);
  printf(   "Sequence database:        %s\n", seqfile);
  printf(   "Search strategy:          global alignment\n");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  /*********************************************** 
   * Search HMM against each sequence
   ***********************************************/

  /* set up structures for storing output
     */
  histogram = AllocHistogram(-200, 200, 100);	 /* keeps full histogram    */
  ghit      = AllocTophits(globH, globA);      /* keeps global seq scores */
  dhit      = AllocTophits(domH, domA);        /* keeps domain scores     */

  nseq = 0;
  while (ReadSeq(sqfp, format, &seq, &sqinfo)) 
    {
      dsq = DigitizeSequence(seq, sqinfo.len);
      
      /* 1. Score the sequence. */
      if (do_forward) sc = Plan7Forward(dsq, sqinfo.len, hmm, NULL);
      else            sc = Plan7Viterbi(dsq, sqinfo.len, hmm, &mx);

      /* 2. Recover a trace.    */
      if (do_forward) Plan7Viterbi(dsq, sqinfo.len, hmm, &mx);
      P7ViterbiTrace(hmm, dsq, sqinfo.len, mx, &tr);

      /* 2. Store score/pvalue for global alignment. */
      if (sc > globT) 
	{
	  pvalue = PValue(hmm, sc);
	  RegisterHit(ghit, -1.*pvalue, pvalue, sc,
		      sqinfo.name, 
		      sqinfo.flags & SQINFO_DESC ? sqinfo.desc : NULL, 
		      0,0,0,                	/* seq positions  */
		      0,0,0,	                /* HMM positions  */
		      0, TraceDomainNumber(tr), /* domain info    */
		      NULL);	                /* alignment info */
	}
      AddToHistogram(histogram, sc);

      /* 3. Parse domains out of the trace */
      if (sc > domT) 
	{
	  record_by_domain(dhit, hmm, dsq, &sqinfo, tr, domT); 
	}

      P7FreeTrace(tr);
      FreePlan7Matrix(mx);
      free(dsq);
      FreeSequence(seq, &sqinfo); 
      nseq++;
    }

  /* We're done searching an HMM over the whole sequence database.
     * Fit the histogram now if we use it for E-values. The "20"
     * is a hint to the fitting program that we expect scores over
     * 20 to be true positives.
     */
  validfit = ExtremeValueFitHistogram(histogram, 20.);

  /* Now format and report our output 
   */

  /* 1. Report overall sequence hits (sorted on E-value) */
  printf("\nQuery HMM:  %s  %s\n", hmm->name, hmm->desc != NULL ? hmm->desc : "");
  FullSortTophits(ghit);
  namewidth = TophitsMaxName(ghit, globH);
  printf("\nScores for complete sequences (score includes all domains):\n");
  printf("%-*s %-*s %7s %10s %3s\n", namewidth, "Sequence", 52-namewidth, "Description", "Score", "E-value", " N ");
  printf("%-*s %-*s %7s %10s %3s\n", namewidth, "--------", 52-namewidth, "-----------", "-----", "-------", "---");
  for (i = 0; i < ghit->pos; i++)
    {
      GetRankedHit(ghit, i, 
		   &pvalue, &sc, &name, &desc,
		   NULL, NULL, NULL,               /* sequence positions */
		   NULL, NULL, NULL,               /* HMM positions      */
		   NULL, &ndom,	                   /* domain info        */
		   NULL);	                   /* alignment info     */
      evalue = pvalue * (double) nseq;
      if (evalue <= globE && sc >= globT) 
	printf("%-*s %-*.*s %7.1f %10.2g %3d\n", 
	       namewidth, name, 
	       52-namewidth, 52-namewidth, desc != NULL ? desc : "",
	       sc, evalue, ndom);
      else if (evalue > globE)
	{
	  if (i > 0) printf("          [no more scores below E = %.2g]\n", globE);
	  break;
	}
      else if (sc < globT)
	{
	  if (i > 0) printf("          [no more scores above T = %.1f]\n", globT);
	  break;
	}
    }
  if (i == 0)
    printf("          [no hits above E = %.2g, T = %.1f]\n", globE, globT);
  if (i == globH)
    printf("          [output cut off at H = %d top hits]\n", globH);

  /* 2. Report domain hits (also sorted on E-value) */
  FullSortTophits(dhit);
  namewidth = TophitsMaxName(dhit, domH);
  printf("\nParsed for domains:\n");
  printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
	 namewidth, "Sequence", "Domain ", "seq-f", "seq-t", "hmm-f", "hmm-t", "score", "E-value");
  printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
	 namewidth, "--------", "-------", "-----", "-----", "-----", "-----", "-----", "-------");
      
  for (i = 0; i < dhit->pos; i++)
    {
      GetRankedHit(dhit, i, 
		   &pvalue, &sc, &name, NULL,
		   &sqfrom, &sqto, &sqlen,            /* seq position info  */
		   &hmmfrom, &hmmto, NULL,            /* HMM position info  */
		   &domidx, &ndom,                    /* domain info        */
		   NULL);	                      /* alignment info     */
      evalue = pvalue * (double) nseq;
      if (evalue <= domE && sc >= domT) 
	printf("%-*s %3d/%-3d %5d %5d %c%c %5d %5d %c%c %7.1f %8.2g\n",
	       namewidth, name, 
	       domidx, ndom,
	       sqfrom, sqto, 
	       sqfrom == 1 ? '[' : '.', sqto == sqlen ? ']' : '.',
	       hmmfrom, hmmto,
	       hmmfrom == 1 ? '[':'.', hmmto == hmm->M ? ']' : '.',
	       sc, pvalue);
      else if (evalue > domE) {
	if (i > 0) printf("          [no more scores below E = %.2g]\n", domE);
	break;
      }
      else if (sc < domT) {
	if (i > 0) printf("          [no more scores above T = %.1f]\n", domT);
	break;
      }
    }
  if (i == 0)    printf("          [no hits above E = %.2g, T = %.1f]\n", domE, domT);
  if (i == domH) printf("          [output cut off at H = %d top hits]\n", domH);


  /* Alignment output would go here [UNFINISHED]
   */

  /* Histogram output */
  printf("\nHistogram of all scores:\n");
  PrintASCIIHistogram(stdout, histogram);

  /*********************************************** 
   * Clean-up and exit.
   ***********************************************/

  printf("//\n");
  FreePlan7(hmm);
  SeqfileRewind(sqfp);
  FreeTophits(ghit);
  FreeTophits(dhit);
  FreeHistogram(histogram);
  HMMFileClose(hmmfp);
  SeqfileClose(sqfp);
  SqdClean();

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) malloc_list(2, histid1, histid2);
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  return 0;
}


