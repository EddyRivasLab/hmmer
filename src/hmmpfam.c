/* hmmpfam.c
 * SRE, Mon Aug 25 17:03:14 1997: Denver 
 *
 * main() for production Pfam searching.
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

static void display_by_domain(struct plan7_s *hmm, char *dsq, SQINFO *sqinfo,
			      struct p7trace_s *tr);
static void record_by_domain(struct tophit_s *h, 
			     struct plan7_s *hmm, char *dsq, SQINFO *sqinfo,
			     struct p7trace_s *tr, double domE, float domT);

static char banner[] = "hmmpfam - search a single seq against HMM database";

static char usage[]  = "\
Usage: hmms [-options] <hmm database> <sequence file>\n\
  Available options are:\n\
   -h        : help; print brief help on version and usage\n\
   -n        : nucleic acid models/sequence (default protein)\n\
   --noxnu   : Turn off XNU filtering\n\
\n";


static struct opt_s OPTIONS[] = {
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-n",        TRUE,  sqdARG_NONE },
  { "--noxnu",   FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv) 
{
  char            *hmmfile;	/* file to read HMMs from                  */
  HMMFILE         *hmmfp;       /* opened hmmfile for reading              */
  char            *seqfile;     /* file to read target sequence from       */ 
  SQFILE          *sqfp;        /* opened seqfile for reading              */
  int              format;	/* format of seqfile                       */
  char              *seq;	/* target sequence                         */
  SQINFO             sqinfo;	/* optional info for seq                   */
  char              *dsq;	/* digitized target sequence               */
  struct plan7_s    *hmm;       /* current HMM to search with              */ 
  struct dpmatrix_s *mx;	/* DP matrix after alignment               */
  struct p7trace_s  *tr;	/* traceback of alignment                  */
  struct fancyali_s *ali;	/* an alignment for display                */

  float            score;	/* log-odds score in bits                  */
  double           pvalue;	/* pvalue of an HMM score                  */
  double           evalue;	/* evalue of an HMM score                  */
  int     sqfrom, sqto;		/* coordinates in sequence                 */
  int     hmmfrom, hmmto;	/* coordinate in HMM                       */
  char   *name, *desc;          /* hit HMM name and description            */
  int     hmmlen;		/* length of HMM hit                       */
  int     nhmm;			/* number of HMMs searched                 */


  struct tophit_s *ghit;        /* list of top hits and alignments for seq  */
  struct tophit_s *dhit;	/* list of top hits/alignments for domains  */
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
  int   do_nucleic;		/* TRUE to do DNA/RNA instead of protein    */
  int   do_xnu;			/* TRUE to do XNU filtering                 */
  
  int   i;

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  do_nucleic  = FALSE;
  do_xnu      = TRUE;
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
    if      (strcmp(optname, "-n")      == 0) { do_nucleic = TRUE; }
    else if (strcmp(optname, "--noxnu") == 0) { do_xnu     = FALSE;}
    else if (strcmp(optname, "-h")      == 0) {
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
   * Open sequence database (must be in curr directory);
   * get target sequence.
   ***********************************************/

  if (do_nucleic) SetAlphabet(hmmNUCLEIC);
  else            SetAlphabet(hmmAMINO);

  if (! SeqfileFormat(seqfile, &format, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: 
      Die("Sequence file %s could not be opened for reading", seqfile); 
      break;
    case SQERR_FORMAT: 
    default:           
      Die("Failed to determine format of sequence file %s", seqfile);
    }
  if ((sqfp = SeqfileOpen(seqfile, format, NULL)) == NULL)
    Die("Failed to open sequence file %s\n%s\n", seqfile, usage);

  /*********************************************** 
   * Open HMM database (might be in HMMERDB or current directory)
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("Failed to open HMM database %s\n%s", hmmfile, usage);

  /*********************************************** 
   * Show the banner
   ***********************************************/

  Banner(stdout, banner);
  printf(   "HMM file:                 %s\n", hmmfile);
  printf(   "Sequence database:        %s\n", seqfile);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  /*********************************************** 
   * Search each HMM against each sequence
   ***********************************************/

  while (ReadSeq(sqfp, format, &seq, &sqinfo)) 
    {
      ghit = AllocTophits(globH, globA); /* keeps full seq scores */
      dhit = AllocTophits(domH, domA);   /* keeps domain scores   */

      /* 1. Search sequence against every HMM.
       */
      dsq = DigitizeSequence(seq, sqinfo.len);
      if (do_xnu) XNU(dsq, sqinfo.len);
      
      nhmm = 0;
      while (HMMFileRead(hmmfp, &hmm)) {
	if (hmm == NULL) 
	  Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);
	Plan7LSConfig(hmm); 
	Plan7Logoddsify(hmm);

	/* 1a. Do the alignment (Viterbi), recover a trace
	 */
	score = Plan7Viterbi(dsq, sqinfo.len, hmm, &mx);
	P7ViterbiTrace(hmm, dsq, sqinfo.len, mx, &tr);

	/* 1b. Store scores/pvalue for each HMM aligned to this sequence, overall
	 */
	pvalue = PValue(hmm, score);
	if (pvalue < globE && score >= globT) 
	  { 
	    RegisterHit(ghit, score, pvalue, score,
			hmm->name, hmm->desc, 
			0,0,0,
			0,0,0,
			0,0,
			NULL);
	  }

	/* 1c. Store scores/evalue/position/alignment for each HMM
	 *    aligned to domains of this sequence; UNFINISHED
	 */
	record_by_domain(dhit, hmm, dsq, &sqinfo, tr, domE, domT);

	P7FreeTrace(tr);
	FreePlan7Matrix(mx);
	FreePlan7(hmm);
	nhmm++;
      }

      /* 2. Report the overall sequence hits
       */
      printf("Query:  %s  %s\n", sqinfo.name, 
	     sqinfo.flags & SQINFO_DESC ? sqinfo.desc : "");
      FullSortTophits(ghit);
				/* temporary */
      printf("  Scores for complete sequence (score includes all domains):\n");
      printf("  %-12s %-40s %7s %10s\n", "HMM", "Description", "Score", "E-value");
      printf("  %-12s %-40s %7s %10s\n", "---", "-----------", "-----", "-------");
      for (i = 0; i < ghit->pos; i++)
	{
	  GetRankedHit(ghit, i, 
		       &pvalue, &score, &name, &desc,
		       NULL, NULL, NULL, 
		       NULL, NULL, NULL, 
		       NULL, NULL,
		       NULL);
	  evalue = pvalue * (double) nhmm;
	  if (evalue <= globE && score >= globT) 
	    printf("  %-12s %-40s %7.1f %10.2g\n", 
		   name, desc != NULL ? desc : "",
		   score, evalue);
	}

      /* 3. Report domain hits (sorted on sqto coordinate)
       */
      FullSortTophits(dhit);
      printf("\n  Parsed for domains:\n");
      printf("  %-12s %6s %6s [%5s] %6s %6s [%5s] %7s %10s\n",
	     "HMM", "seq-f", "seq-t", "sqlen", "hmm-f", "hmm-t", "hmm-M", "score", "E-value");
      printf("  %-12s %6s %6s [%5s] %6s %6s [%5s] %7s %10s\n",
	     "---", "-----", "-----", "-----", "-----", "-----", "-----", "-----", "-------");
      
      for (i = 0; i < dhit->pos; i++)
	{
	  GetRankedHit(dhit, i, 
		       &pvalue, &score, &name, NULL,
		       &sqfrom, &sqto, NULL, 
		       &hmmfrom, &hmmto, &hmmlen, 
		       NULL, NULL,
		       NULL);
	  evalue = pvalue * (double) nhmm;
	  if (evalue <= domE && score >= domT) 
	  printf("  %-12s %6d %6d [%5d] %6d %6d [%5d] %7.1f %10.2g\n",
		 name, sqfrom, sqto, sqinfo.len, hmmfrom, hmmto, hmmlen,
		 score, pvalue);
	}
      printf("//\n");

      FreeSequence(seq, &sqinfo); 
      free(dsq);
      HMMFileRewind(hmmfp);
      FreeTophits(ghit);
      FreeTophits(dhit);
    }

  /*********************************************** 
   * Clean-up and exit.
   ***********************************************/
  SeqfileClose(sqfp);
  HMMFileClose(hmmfp);
  SqdClean();

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) malloc_list(2, histid1, histid2);
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  return 0;
}


/* Function: display_by_domain()
 * Date:     Sat Aug 30 11:21:35 1997 (Denver CO)  
 * 
 * Purpose:  Display scores, P-values, alignments for 
 *           individual domains.
 *           
 * Return:   (void)          
 */
static void
display_by_domain(struct plan7_s *hmm, char *dsq, SQINFO *sqinfo,
		  struct p7trace_s *tr)
{
  struct p7trace_s **tarr;      /* array of per-domain traces */
  struct fancyali_s *ali;	/* alignment for display      */
  int ntr;			/* number of domain traces    */
  int idx;			/* index for traces           */
  float sc;			/* score of a trace           */
  int k1, k2;			/* start, stop coord in model */
  int i1, i2;			/* start, stop in sequence    */

  TraceDecompose(tr, &tarr, &ntr);
  if (ntr == 0) return;		/* "can't happen" */

  for (idx = 0; idx < ntr; idx++)
    {
      /* Get the score of the match.
       */
      sc = P7TraceScore(hmm, dsq, tarr[idx]);
      /*P7PrintTrace(stdout, tarr[idx], hmm, dsq); */
      TraceSimpleBounds(tarr[idx], &i1, &i2, &k1, &k2);

      /* Here's where we'd evaluate the P-value of the score,
       * but NOT FINISHED
       */
      
      /* Print out the match
       */
      printf("   Domain %4d: %-7.1f %6d %6d [%d] %6d %6d [%d]\n",
	     idx+1, sc, i1, i2, sqinfo->len, k1, k2, hmm->M);
    }

  for (idx = 0; idx < ntr; idx++)
    {
      ali = CreateFancyAli(tarr[idx], hmm, dsq, sqinfo->name);
      PrintFancyAli(stdout, ali);
      P7FreeTrace(tarr[idx]);
      FreeFancyAli(ali);
    }
  free(tarr);
  return;
}


/* Function: record_by_domain()
 * Date:     SRE, Tue Oct 28 14:20:56 1997 [Cambridge UK]
 * 
 * Purpose:  Decompose a trace, and register scores, P-values, alignments,
 *           etc. for individual domains in a hitlist. 
 *           
 * Return:   (void)          
 */
static void
record_by_domain(struct tophit_s *h, 
		 struct plan7_s *hmm, char *dsq, SQINFO *sqinfo,
		 struct p7trace_s *tr, double domE, float domT)
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
      
      /* Record the match
       */
      if (pvalue <= domE && score > domT)
	{
	  ali = CreateFancyAli(tarr[idx], hmm, dsq, sqinfo->name);
	  RegisterHit(h, -1.*(double)i2, pvalue, score,
		      hmm->name, hmm->desc, 
		      i1,i2, sqinfo->len, 
		      k1,k2, hmm->M, 
		      0,0,
		      ali);
	}
    }

  for (idx = 0; idx < ntr; idx++)
    P7FreeTrace(tarr[idx]);
  free(tarr);
  return;
}
