/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998,
 * Sean R. Eddy and Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmmsearch.c
 * SRE, Tue Jan  7 17:19:20 1997
 *
 * Search a sequence database with a profile HMM.
 * RCS $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

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
   -A <n>    : sets alignment output limit to <n> best domain alignments\n\
   -E <x>    : sets E value cutoff (globE) to <x>\n\
   -T <x>    : sets T bit threshold (globT) to <x>\n\
";

static char experts[] = "\
   --domE <x>: sets domain Eval cutoff (2nd threshold) to <x>\n\
   --domT <x>: sets domain T bit thresh (2nd threshold) to <x>\n\
   --forward : use the full Forward() algorithm instead of Viterbi\n\
   --null2   : turn OFF the post hoc second null model\n\
   --xnu     : turn ON XNU filtering of target protein sequences\n\
";

static struct opt_s OPTIONS[] = {
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-A",        TRUE,  sqdARG_INT  },  
  { "-E",        TRUE,  sqdARG_FLOAT},  
  { "-T",        TRUE,  sqdARG_FLOAT},  
  { "--domE",    FALSE, sqdARG_FLOAT},
  { "--domT",    FALSE, sqdARG_FLOAT},
  { "--forward", FALSE, sqdARG_NONE },
  { "--null2",   FALSE, sqdARG_NONE },
  { "--xnu",     FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static void record_domains(struct tophit_s *h, 
			   struct plan7_s *hmm, char *dsq, SQINFO *sqinfo,
			   struct p7trace_s *tr, double whole_pval, float whole_sc,
			   int do_null2);

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
  int       i; 
  struct plan7_s  *hmm;         /* HMM to search with                      */ 
  struct histogram_s *histogram;/* histogram of all scores                 */
  struct dpmatrix_s *mx;	/* DP matrix after alignment               */
  struct p7trace_s  *tr;	/* traceback                               */
  struct fancyali_s *ali;       /* displayed alignment info                */ 
  struct tophit_s   *ghit;      /* list of top hits for whole sequences    */
  struct tophit_s   *dhit;	/* list of top hits for domains            */

  float     sc;	        	/* score of an HMM search                  */
  double  pvalue;		/* pvalue of an HMM score                  */
  double  evalue;		/* evalue of an HMM score                  */
  double  motherp;		/* pvalue of a whole seq HMM score         */
  float   mothersc;		/* score of a whole seq parent of domain   */
  int     sqfrom, sqto;		/* coordinates in sequence                 */
  int     hmmfrom, hmmto;	/* coordinate in HMM                       */
  char   *name, *desc;          /* hit sequence name and description       */
  int     sqlen;		/* length of seq that was hit              */
  int     nseq;			/* number of sequences searched            */
  int     domidx;		/* number of this domain                   */
  int     ndom;			/* total # of domains in this seq          */
  int     namewidth;		/* max width of sequence name              */

  int    Alimit;		/* A parameter limiting output alignments   */
  float  globT;			/* T parameter: keep only hits > globT bits */
  double globE;			/* E parameter: keep hits < globE E-value   */
  float  domT;			/* T parameter for individual domains       */
  double domE;			/* E parameter for individual domains       */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */
  int   do_null2;		/* TRUE to adjust scores with null model #2 */
  int   do_forward;		/* TRUE to use Forward() not Viterbi()      */
  int   do_xnu;			/* TRUE to filter sequences thru XNU        */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  do_forward  = FALSE;
  do_null2    = TRUE;
  do_xnu      = FALSE;

  Alimit      = INT_MAX;	/* no limit on alignment output     */
  globE       = 10.0;		/* use a reasonable Eval threshold; */
  globT       = -FLT_MAX;	/*   but no bit threshold,          */
  domT        = -FLT_MAX;	/*   no domain bit threshold,       */
  domE        = FLT_MAX;        /*   and no domain Eval threshold.  */

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-A") == 0)        Alimit     = atoi(optarg);  
    else if (strcmp(optname, "-E") == 0)        globE      = atof(optarg);
    else if (strcmp(optname, "-T") == 0)        globT      = atof(optarg);
    else if (strcmp(optname, "--domE")    == 0) domE       = atof(optarg);
    else if (strcmp(optname, "--domT")    == 0) domT       = atof(optarg);
    else if (strcmp(optname, "--forward") == 0) do_forward = TRUE;
    else if (strcmp(optname, "--null2")   == 0) do_null2   = FALSE;
    else if (strcmp(optname, "--xnu")     == 0) do_xnu     = TRUE;
    else if (strcmp(optname, "-h") == 0) {
      Banner(stdout, banner);
      puts(usage);
      puts(experts);
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
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  Plan7Logoddsify(hmm);
  
  /*********************************************** 
   * Show the banner
   ***********************************************/

  Banner(stdout, banner);
  printf(   "HMM file:                 %s [%s]\n", hmmfile, hmm->name);
  printf(   "Sequence database:        %s\n", seqfile);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  /*********************************************** 
   * Search HMM against each sequence
   ***********************************************/

  /* set up structures for storing output
     */
  histogram = AllocHistogram(-200, 200, 100);  /* keeps full histogram */
  ghit      = AllocTophits(200);         /* per-seq hits: 200=lumpsize */
  dhit      = AllocTophits(200);         /* domain hits:  200=lumpsize */

  nseq = 0;
  while (ReadSeq(sqfp, format, &seq, &sqinfo)) 
    {
      nseq++;
      dsq = DigitizeSequence(seq, sqinfo.len);

      if (do_xnu) XNU(dsq, sqinfo.len);
      
      /* 1. Score the sequence. */
      if (do_forward) sc  = Plan7Forward(dsq, sqinfo.len, hmm, NULL);
      else            sc  = Plan7Viterbi(dsq, sqinfo.len, hmm, &mx);

      /* 2. Recover a trace.    */
      if (do_forward) Plan7Viterbi(dsq, sqinfo.len, hmm, &mx);
      P7ViterbiTrace(hmm, dsq, sqinfo.len, mx, &tr);

      if (do_null2)  sc -= TraceScoreCorrection(hmm, tr, dsq);

#if DEBUGLEVEL >= DEBUG_LOTS
      P7PrintTrace(stdout, tr, hmm, dsq); 
#endif

      /* 2. Store score/pvalue for global alignment; will sort on score. 
       *    Keep all domains in a significant sequence hit.
       *    We can only make a lower bound estimate of E-value since
       *    we don't know the final value of nseq yet. 
       */
      pvalue = PValue(hmm, sc);
      if (sc > globT && pvalue * (double) nseq < globE) 
	{
	  RegisterHit(ghit, sc, pvalue, sc, 
		      0., 0.,	                /* no mother seq */
		      sqinfo.name, 
		      sqinfo.flags & SQINFO_DESC ? sqinfo.desc : NULL, 
		      0,0,0,                	/* seq positions  */
		      0,0,0,	                /* HMM positions  */
		      0, TraceDomainNumber(tr), /* domain info    */
		      NULL);	                /* alignment info */

	  record_domains(dhit, hmm, dsq, &sqinfo, tr, pvalue, sc, do_null2); 
	}
      AddToHistogram(histogram, sc);

      FreeSequence(seq, &sqinfo); 
      FreePlan7Matrix(mx);
      P7FreeTrace(tr);
      free(dsq);
    }

  /* We're done searching an HMM over the whole sequence database.
   * Set the theoretical EVD curve in our histogram using 
   * calibration in the HMM, if available. 
   */
  if (hmm->flags & PLAN7_STATS)
    ExtremeValueSetHistogram(histogram, hmm->mu, hmm->lambda, 
			     histogram->lowscore, histogram->highscore, 
			     1.0, 0);

  /* Now format and report our output 
   */

  /* 1. Report overall sequence hits (sorted on E-value) */
  printf("\nQuery HMM:  %s  %s\n", 
	 hmm->name, hmm->desc != NULL ? hmm->desc : "");
  if (hmm->flags & PLAN7_STATS)
    printf("  [HMM has been calibrated; E-values are empirical estimates]\n");
  else
    printf("  [No calibration for HMM; E-values are upper bounds]\n");

  FullSortTophits(ghit);
  namewidth = MAX(8, TophitsMaxName(ghit));

  printf("\nScores for complete sequences (score includes all domains):\n");
  printf("%-*s %-*s %7s %10s %3s\n", namewidth, "Sequence", 52-namewidth, "Description", "Score", "E-value", " N ");
  printf("%-*s %-*s %7s %10s %3s\n", namewidth, "--------", 52-namewidth, "-----------", "-----", "-------", "---");
  for (i = 0; i < ghit->num; i++)
    {
      GetRankedHit(ghit, i, 
		   &pvalue, &sc, NULL, NULL,
		   &name, &desc,
		   NULL, NULL, NULL,               /* sequence positions */
		   NULL, NULL, NULL,               /* HMM positions      */
		   NULL, &ndom,	                   /* domain info        */
		   NULL);	                   /* alignment info     */
      evalue = pvalue * (double) nseq;
      if (evalue < globE && sc > globT) 
	printf("%-*s %-*.*s %7.1f %10.2g %3d\n", 
	       namewidth, name, 
	       52-namewidth, 52-namewidth, desc != NULL ? desc : "",
	       sc, evalue, ndom);
      else if (evalue >= globE)
	{
	  if (i > 0) printf("\t[no more scores below E threshold]\n");
	  break;
	}
      else if (sc <= globT)
	{
	  if (i > 0) printf("\t[no more scores above T threshold]");
	  break;
	}
    }
  if (i == 0) printf("\t[no hits above thresholds]\n");




  /* 2. Report domain hits (also sorted on E-value) */
  FullSortTophits(dhit);
  namewidth = MAX(8, TophitsMaxName(dhit));

  printf("\nParsed for domains:\n");
  printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
	 namewidth, "Sequence", "Domain ", "seq-f", "seq-t", "hmm-f", "hmm-t", "score", "E-value");
  printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
	 namewidth, "--------", "-------", "-----", "-----", "-----", "-----", "-----", "-------");
      
  for (i = 0; i < dhit->num; i++)
    {
      GetRankedHit(dhit, i, 
		   &pvalue, &sc, &motherp, &mothersc,
		   &name, NULL,
		   &sqfrom, &sqto, &sqlen,            /* seq position info  */
		   &hmmfrom, &hmmto, NULL,            /* HMM position info  */
		   &domidx, &ndom,                    /* domain info        */
		   NULL);	                      /* alignment info     */
      evalue = pvalue * (double) nseq;

      if (motherp * (double) nseq >= globE || mothersc <= globT) 
	continue;
      else if (evalue < domE && sc > domT)
	printf("%-*s %3d/%-3d %5d %5d %c%c %5d %5d %c%c %7.1f %8.2g\n",
	       namewidth, name, 
	       domidx, ndom,
	       sqfrom, sqto, 
	       sqfrom == 1 ? '[' : '.', sqto == sqlen ? ']' : '.',
	       hmmfrom, hmmto,
	       hmmfrom == 1 ? '[':'.', hmmto == hmm->M ? ']' : '.',
	       sc, evalue);
      else if (evalue >= domE) {
	if (i > 0) printf("\t[no more scores below domE threshold]\n");
	break;
      }
      else if (sc <= domT) {
	if (i > 0) printf("\t[no more scores above domT threshold]\n");
	break;
      }
    }
  if (i == 0) printf("\t[no hits above thresholds\n");



  /* 3. Alignment output, also by domain.
   *    dhits is already sorted and namewidth is set, from above code.
   *    Number of displayed alignments is limited by Alimit parameter;
   *    also by domE (evalue threshold), domT (score theshold).
   */
  if (Alimit != 0)
    {
      printf("\nAlignments of top-scoring domains:\n");
      for (i = 0; i < dhit->num; i++)
	{
	  if (i == Alimit) break; /* limit to Alimit output alignments */
	  GetRankedHit(dhit, i, 
		       &pvalue, &sc, &motherp, &mothersc,
		       &name, NULL,
		       &sqfrom, &sqto, &sqlen,            /* seq position info  */
		       &hmmfrom, &hmmto, NULL,            /* HMM position info  */
		       &domidx, &ndom,                    /* domain info        */
		       &ali);	                      /* alignment info     */
	  evalue = pvalue * (double) nseq;

	  if (motherp * (double) nseq >= globE || mothersc <= globT) 
	    continue;
	  else if (evalue < domE && sc > domT) 
	    {
	      printf("%s: domain %d of %d, from %d to %d: score %.1f, E = %.2g\n", 
		     name, domidx, ndom, sqfrom, sqto, sc, evalue);
	      PrintFancyAli(stdout, ali);
	    }
	  else if (evalue >= domE) {
	    if (i > 0) printf("\t[no more alignments below domE threshold\n");
	    break;
	  }
	  else if (sc <= domT) {
	    if (i > 0) printf("\t[no more alignments above domT threshold\n");
	    break;
	  }
	}
      if (i == 0)      printf("\t[no hits above thresholds\n");
      if (i == Alimit) printf("\t[output cut off at A = %d top alignments]\n", Alimit);
    }

  /* 4. Histogram output */
  printf("\nHistogram of all scores:\n");
  PrintASCIIHistogram(stdout, histogram);

  /* 5. Tophits summaries, while developing...
   */
  printf("\nWhole sequence top hits:\n");
  TophitsReport(ghit, globE, nseq);
  printf("\nDomain top hits:\n");
  TophitsReport(dhit, domE, nseq);

  /*********************************************** 
   * Clean-up and exit.
   ***********************************************/

  FreeHistogram(histogram);
  HMMFileClose(hmmfp);
  SeqfileClose(sqfp);
  FreeTophits(ghit);
  FreeTophits(dhit);
  FreePlan7(hmm);
  SqdClean();

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) malloc_list(2, histid1, histid2);
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  return 0;
}


/* Function: record_domains()
 * Date:     SRE, Tue Nov  4 11:25:14 1997 [St. Louis]
 * 
 * Purpose:  Decompose a trace, and register scores, P-values, alignments,
 *           etc. for individual domains in a hitlist. 
 *           
 * Args:     hmm    - the HMM structure
 *           dsq    - digitized sequence 1..L
 *           sqinfo - contains name of sequence, etc.
 *           tr     - traceback of the whole sequence aligned to HMM
 *           whole_pval - P-value of complete alignment
 *           whole_sc   - score of complete alignment (bits)
 *           do_null2   - TRUE to use post hoc null model correction 
 *           
 * Return:   (void)          
 */
static void
record_domains(struct tophit_s *h, 
	       struct plan7_s *hmm, char *dsq, SQINFO *sqinfo,
	       struct p7trace_s *tr, double whole_pval, float whole_sc,
	       int do_null2)
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
  if (ntr == 0) Die("TraceDecompose() screwup"); /* "can't happen" (!) */

  for (idx = 0; idx < ntr; idx++)
    {
      /* Get the score and bounds of the match.
       */
      score  = P7TraceScore(hmm, dsq, tarr[idx]);
      if (do_null2) 
	score -= TraceScoreCorrection(hmm, tarr[idx], dsq);
      pvalue = PValue(hmm, score); 
      TraceSimpleBounds(tarr[idx], &i1, &i2, &k1, &k2);
      
#if DEBUGLEVEL >= DEBUG_LOTS
      P7PrintTrace(stdout, tarr[idx], hmm, dsq); 
#endif

      /* Record the match. Use score as the sort key.
       */
      ali = CreateFancyAli(tarr[idx], hmm, dsq, sqinfo->name);
      RegisterHit(h, score, pvalue, score, whole_pval, whole_sc,
		  sqinfo->name, 
		  sqinfo->flags & SQINFO_DESC ? sqinfo->desc : NULL, 
		  i1,i2, sqinfo->len, 
		  k1,k2, hmm->M, 
		  idx+1, ntr,
		  ali);
    }
  for (idx = 0; idx < ntr; idx++)
    P7FreeTrace(tarr[idx]);
  free(tarr);
  return;
}
