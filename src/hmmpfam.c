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

/* hmmpfam.c
 * SRE, Mon Aug 25 17:03:14 1997: Denver 
 *
 * Search a single sequence against an HMM database.
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

static void record_domains(struct tophit_s *h, 
			   struct plan7_s *hmm, char *dsq, SQINFO *sqinfo,
			   struct p7trace_s *tr, double whole_pval, float whole_sc,
			   int do_null2);

static char banner[] = "hmmpfam - search a single seq against HMM database";

static char usage[]  = "\
Usage: hmmpfam [-options] <hmm database> <sequence file>\n\
  Available options are:\n\
   -h        : help; print brief help on version and usage\n\
   -n        : nucleic acid models/sequence (default protein)\n\
   -A <n>    : sets alignment output limit to <n> best domain alignments\n\
   -E <x>    : sets E value cutoff (globE) to <x>\n\
   -T <x>    : sets T bit threshold (globT) to <x>\n\
   -Z <n>    : sets Z (# seqs) for E-value calculation\n\
";

static char experts[] = "\
   --domE <x>: sets domain Eval cutoff (2nd threshold) to <x>\n\
   --domT <x>: sets domain T bit thresh (2nd threshold) to <x>\n\
   --forward : use the full Forward() algorithm instead of Viterbi\n\
   --null2   : turn OFF the post hoc second null model\n\
   --xnu     : turn ON XNU filtering of query sequence\n\
\n";


static struct opt_s OPTIONS[] = {
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-n",        TRUE,  sqdARG_NONE },
  { "-A",        TRUE,  sqdARG_INT  },  
  { "-E",        TRUE,  sqdARG_FLOAT},  
  { "-T",        TRUE,  sqdARG_FLOAT},  
  { "-Z",        TRUE,  sqdARG_INT },
  { "--domE",    FALSE, sqdARG_FLOAT},
  { "--domT",    FALSE, sqdARG_FLOAT},
  { "--forward", FALSE, sqdARG_NONE },
  { "--null2",   FALSE, sqdARG_NONE },  
  { "--xnu",     FALSE, sqdARG_NONE },
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
  struct tophit_s   *ghit;      /* list of top hits and alignments for seq  */
  struct tophit_s   *dhit;	/* list of top hits/alignments for domains  */

  float   sc;			/* log-odds score in bits                  */
  double  pvalue;		/* pvalue of an HMM score                  */
  double  evalue;		/* evalue of an HMM score                  */
  double  motherp;		/* pvalue of a whole seq HMM score         */
  float   mothersc;		/* score of a whole seq parent of domain   */
  int     sqfrom, sqto;		/* coordinates in sequence                 */
  int     hmmfrom, hmmto;	/* coordinate in HMM                       */
  char   *name, *desc;          /* hit HMM name and description            */
  int     hmmlen;		/* length of HMM hit                       */
  int     nhmm;			/* number of HMMs searched                 */
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
  int   do_forward;		/* TRUE to use Forward() not Viterbi()      */
  int   do_nucleic;		/* TRUE to do DNA/RNA instead of protein    */
  int   do_null2;		/* TRUE to adjust scores with null model #2 */
  int   do_xnu;			/* TRUE to do XNU filtering                 */
  int   Z;			/* nseq to base E value calculation on      */
  
  int   i;

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  do_forward  = FALSE;
  do_nucleic  = FALSE;
  do_null2    = TRUE;
  do_xnu      = FALSE;
  Z           = 59021;		/* default: nseq in Swissprot34     */
  
  Alimit      = INT_MAX;	/* no limit on alignment output     */
  globE       = 10.0;		/* use a reasonable Eval threshold; */
  globT       = -FLT_MAX;	/*   but no bit threshold,          */
  domT        = -FLT_MAX;	/*   no domain bit threshold,       */
  domE        = FLT_MAX;        /*   and no domain Eval threshold.  */

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-n")        == 0) do_nucleic = TRUE; 
    else if (strcmp(optname, "-A")        == 0) Alimit     = atoi(optarg);  
    else if (strcmp(optname, "-E")        == 0) globE      = atof(optarg);
    else if (strcmp(optname, "-T")        == 0) globT      = atof(optarg);
    else if (strcmp(optname, "-Z")        == 0) Z          = atoi(optarg);
    else if (strcmp(optname, "--domE")    == 0) domE       = atof(optarg);
    else if (strcmp(optname, "--domT")    == 0) domT       = atof(optarg);
    else if (strcmp(optname, "--forward") == 0) do_forward = TRUE;
    else if (strcmp(optname, "--null2")   == 0) do_null2   = FALSE;
    else if (strcmp(optname, "--xnu")     == 0) do_xnu     = TRUE;
    else if (strcmp(optname, "-h")      == 0) {
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
  printf(   "Sequence file:            %s\n", seqfile);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  /*********************************************** 
   * Search each HMM against each sequence
   ***********************************************/

  while (ReadSeq(sqfp, format, &seq, &sqinfo)) 
    {
      ghit = AllocTophits(20);   /* keeps full seq scores */
      dhit = AllocTophits(20);   /* keeps domain scores   */

      /* 1. Search sequence against every HMM.
       */
      dsq = DigitizeSequence(seq, sqinfo.len);
      if (do_xnu) XNU(dsq, sqinfo.len);
      
      nhmm = 0;
      while (HMMFileRead(hmmfp, &hmm)) {
	if (hmm == NULL) 
	  Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);
	Plan7Logoddsify(hmm);

	/* 1a. Score sequence, do alignment (Viterbi), recover trace
	 */
	if (do_forward) sc  = Plan7Forward(dsq, sqinfo.len, hmm, NULL);
	else            sc  = Plan7Viterbi(dsq, sqinfo.len, hmm, &mx);

	if (do_forward) (void) Plan7Viterbi(dsq, sqinfo.len, hmm, &mx);
	P7ViterbiTrace(hmm, dsq, sqinfo.len, mx, &tr);

	if (do_null2)  sc -= TraceScoreCorrection(hmm, tr, dsq);

	/* 1b. Store scores/pvalue for each HMM aligned to this sequence, overall
	 */
	pvalue = PValue(hmm, sc);
	if (sc > globT && pvalue * (float) Z < globE) 
	  { 
	    RegisterHit(ghit, sc, pvalue, sc,
			0., 0.,	                  /* no mother seq */
			hmm->name, hmm->desc, 
			0,0,0,                	  /* seq positions  */
			0,0,0,	                  /* HMM positions  */
			0, TraceDomainNumber(tr), /* domain info    */
			NULL);	                  /* alignment info */
	    /* 1c. Store scores/evalue/position/alignment for each HMM
	     *    aligned to domains of this sequence; UNFINISHED
	     */
	    record_domains(dhit, hmm, dsq, &sqinfo, tr, pvalue, sc, do_null2);
	  }

	FreePlan7Matrix(mx);
	P7FreeTrace(tr);
	FreePlan7(hmm);
	nhmm++;
      }

      /* 2. (Done searching all HMMs for this query seq; start output)
       *    Report the overall sequence hits, sorted by significance.
       */
      printf("Query:  %s  %s\n", sqinfo.name, 
	     sqinfo.flags & SQINFO_DESC ? sqinfo.desc : "");
      FullSortTophits(ghit);
      namewidth =  MAX(8, TophitsMaxName(ghit));

      printf("\nScores for sequence family classification (score includes all domains):\n");
      printf("%-*s %-*s %7s %10s %3s\n", namewidth, "Sequence", 52-namewidth, "Description", "Score", "E-value", " N ");
      printf("%-*s %-*s %7s %10s %3s\n", namewidth, "--------", 52-namewidth, "-----------", "-----", "-------", "---");
      for (i = 0; i < ghit->num; i++)
	{
	  GetRankedHit(ghit, i, 
		       &pvalue, &sc, NULL, NULL,
		       &name, &desc,
		       NULL, NULL, NULL,           /* seq positions */
		       NULL, NULL, NULL,           /* HMM positions */
		       NULL, &ndom,                /* domain info   */
		       NULL);	                   /* alignment info*/

	  evalue = pvalue * (double) Z;
	  if (evalue < globE && sc >= globT) 
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

      /* 3. Report domain hits (sorted on sqto coordinate)
       */
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
		       &sqfrom, &sqto, NULL, 
		       &hmmfrom, &hmmto, &hmmlen, 
		       &domidx, &ndom,
		       NULL);
	  evalue = pvalue * (double) Z;
	  
	  if (motherp * (double) Z >= globE || mothersc <= globT) 
	    continue;
	  else if (evalue < domE && sc > domT)
	    printf("%-*s %3d/%-3d %5d %5d %c%c %5d %5d %c%c %7.1f %8.2g\n",
		   namewidth, name, 
		   domidx, ndom,
		   sqfrom, sqto, 
		   sqfrom == 1 ? '[' : '.', sqto == sqinfo.len ? ']' : '.',
		   hmmfrom, hmmto,
		   hmmfrom == 1 ? '[':'.',  hmmto == hmmlen ? ']' : '.',
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
			   &sqfrom, &sqto, NULL,              /* seq position info  */
			   &hmmfrom, &hmmto, &hmmlen,         /* HMM position info  */
			   &domidx, &ndom,                    /* domain info        */
			   &ali);	                      /* alignment info     */
	      evalue = pvalue * (double) Z;

	      if (motherp * (double) Z >= globE || mothersc <= globT) 
		continue;
	      else if (evalue < domE && sc > domT) 
		{
		  printf("%s: domain %d of %d, from %d to %d: score %.1f, E = %.2g\n", 
			 name, domidx, ndom, sqfrom, sqto, sc, evalue);
		  PrintFancyAli(stdout, ali);
		}
	      else if (evalue >= domE) {
		if (i > 0) printf("\t[no more alignments below E threshold\n");
		break;
	      }
	      else if (sc <= domT) {
		if (i > 0) printf("\t[no more alignments above T threshold\n");
		break;
	      }
	    }
	  if (i == 0)      printf("\t[no hits above thresholds\n");
	  if (i == Alimit) printf("\t[output cut off at A = %d top alignments]\n", Alimit);
	}


      printf("//\n");
      FreeSequence(seq, &sqinfo); 
      FreeTophits(ghit);
      FreeTophits(dhit);
      free(dsq);

      HMMFileRewind(hmmfp);
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



/* Function: record_domains()
 * Date:     SRE, Tue Oct 28 14:20:56 1997 [Cambridge UK]
 * 
 * Purpose:  Decompose a trace, and register scores, P-values, alignments,
 *           etc. for individual domains in a hitlist. Almost a duplicate
 *           of the eponymous function in hmmsearch, but we sort by
 *           start position of the domain. 
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
      
      /* Record the match
       */
      ali = CreateFancyAli(tarr[idx], hmm, dsq, sqinfo->name);
      RegisterHit(h, -1.*(double)i2,           /* sortkey= -(start) (so 1 comes first) */
		  pvalue, score, whole_pval, whole_sc,
		  hmm->name, hmm->desc, 
		  i1,i2, sqinfo->len, 
		  k1,k2, hmm->M, 
		  idx+1,ntr,
		  ali);
    }
  
  for (idx = 0; idx < ntr; idx++)
    P7FreeTrace(tarr[idx]);
  free(tarr);
  return;
}
