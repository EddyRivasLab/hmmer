/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmmalign.c
 * SRE, Thu Dec 18 16:05:29 1997 [St. Louis]
 * 
 * main() for aligning a set of sequences to an HMM.
 * RCS $Id$
 */ 

#include <stdio.h>
#include <stdlib.h>

#include "structs.h"		/* data structures, macros, #define's   */
#include "config.h"		/* compile-time configuration constants */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "squid.h"		/* general sequence analysis library    */

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static char banner[] = "hmmalign - align sequences to an HMM profile";

static char usage[]  = "\
Usage: hmmalign [-options] <hmm file> <sequence file>\n\
Available options are:\n\
   -h     : help; print brief help on version and usage\n\
   -m     : only print symbols aligned to match states\n\
   -o <f> : save alignment in file <f> in SELEX format\n\
   -q     : quiet - suppress verbose banner\n\
";

static char experts[] = "\
\n";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE   }, 
  { "-m", TRUE, sqdARG_NONE   } ,
  { "-o", TRUE, sqdARG_STRING },
  { "-q", TRUE, sqdARG_NONE   },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv) 
{
  char            *hmmfile;	/* file to read HMMs from                  */
  HMMFILE         *hmmfp;       /* opened hmmfile for reading              */
  char            *seqfile;     /* file to read target sequence from       */ 
  int              format;	/* format of seqfile                       */
  char           **rseq;        /* raw, unaligned sequences                */ 
  SQINFO          *sqinfo;      /* info associated with sequences          */
  char           **dsq;         /* digitized raw sequences                 */
  int              nseq;        /* number of sequences                     */  
  char           **aseq;        /* aligned sequences                       */
  AINFO            ainfo;       /* alignment information                   */
  float           *wgt;         /* per-sequence weights                    */
  int              i;
  struct plan7_s    *hmm;       /* HMM to align to                         */ 
  struct p7trace_s **tr;        /* traces for aligned sequences            */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */
  int   be_quiet;		/* TRUE to suppress verbose banner          */
  int   matchonly;		/* TRUE to show only match state syms       */
  char *outfile;                /* optional alignment output file           */
  FILE *ofp;                    /* handle on alignment output file          */

  
#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  matchonly = FALSE;
  outfile   = NULL;
  be_quiet  = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-m") == 0) matchonly= TRUE;
    else if (strcmp(optname, "-o") == 0) outfile  = optarg;
    else if (strcmp(optname, "-q") == 0) be_quiet = TRUE; 
    else if (strcmp(optname, "-h") == 0) 
      {
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
  * Open HMM file (might be in HMMERDB or current directory).
  * Read a single HMM from it.
  ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  if (!HMMFileRead(hmmfp, &hmm)) 
    Die("Failed to read any HMMs from %s\n", hmmfile);
  HMMFileClose(hmmfp);
  if (hmm == NULL) 
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  P7Logoddsify(hmm, TRUE);
  
  /*********************************************** 
   * Open sequence file in current directory.
   * Read all seqs from it.
   ***********************************************/

  if (! SeqfileFormat(seqfile, &format, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: 
      Die("Sequence file %s could not be opened for reading", seqfile); /*FALLTHRU*/
    case SQERR_FORMAT: 
    default:           
      Die("Failed to determine format of sequence file %s", seqfile);
    }
  if (! ReadMultipleRseqs(seqfile, format, &rseq, &sqinfo, &nseq))
    Die("Failed to read any sequences from file %s", seqfile);

  /*********************************************** 
   * Show the banner
   ***********************************************/

  if (! be_quiet) 
    {
      Banner(stdout, banner);
      printf(   "HMM file:             %s\n", hmmfile);
      printf(   "Sequence file:        %s\n", seqfile);
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
    }

  /*********************************************** 
   * Do the work
   ***********************************************/

  /* Allocations and initializations.
   */
  dsq = MallocOrDie(sizeof(char *) * nseq);
  tr  = MallocOrDie(sizeof(struct p7trace_s *) * nseq);
  wgt = MallocOrDie(sizeof(float) * nseq);
  FSet(wgt, nseq, 1.0);

  /* Align each sequence to the model, collect traces
   */
  for (i = 0; i < nseq; i++)
    {
      dsq[i] = DigitizeSequence(rseq[i], sqinfo[i].len);

      if (P7ViterbiSize(sqinfo[i].len, hmm->M) <= RAMLIMIT)
	(void) P7Viterbi(dsq[i], sqinfo[i].len, hmm, &(tr[i]));
      else
	(void) P7SmallViterbi(dsq[i], sqinfo[i].len, hmm, &(tr[i]));
    }

  /* Turn traces into a multiple alignment
   */ 
  P7Traces2Alignment(dsq, sqinfo, wgt, nseq, hmm->M, tr, matchonly,
		     &aseq, &ainfo);

  /*********************************************** 
   * Output the alignment
   ***********************************************/
  
  if (outfile != NULL && (ofp = fopen(outfile, "w")) != NULL)
    {
      WriteSELEX(ofp, aseq, &ainfo, 50);
      printf("Alignment saved in file %s\n", outfile);
      fclose(ofp);
    }
  else
    WriteSELEX(stdout, aseq, &ainfo, 50);

  /*********************************************** 
   * Cleanup and exit
   ***********************************************/
  
  for (i = 0; i < nseq; i++) 
    {
      P7FreeTrace(tr[i]);
      FreeSequence(rseq[i], &(sqinfo[i]));
      free(dsq[i]);
    }
  FreeAlignment(aseq, &ainfo);
  FreePlan7(hmm);
  free(sqinfo);
  free(dsq);
  free(wgt);
  free(tr);

  SqdClean();

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) malloc_list(2, histid1, histid2);
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  return 0;
}
