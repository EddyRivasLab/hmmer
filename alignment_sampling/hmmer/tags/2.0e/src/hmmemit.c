/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998, Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmmemit.c
 * SRE, Sun Mar  8 14:11:24 1998 [St. Louis]
 * 
 * main() for generating sequences from an HMM
 * RCS $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "structs.h"		/* data structures, macros, #define's   */
#include "config.h"		/* compile-time configuration constants */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "squid.h"		/* general sequence analysis library    */

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static char banner[] = "hmmemit - generate sequences from a profile HMM";

static char usage[]  = "\
Usage: hmmemit [-options] <hmm file>\n\
Available options are:\n\
   -a     : write generated sequences as an alignment, not FASTA\n\
   -h     : help; print brief help on version and usage\n\
   -n <n> : emit <n> sequences (default 10)\n\
   -o <f> : save sequences in file <f>\n\
   -q     : quiet - suppress verbose banner\n\
";

static char experts[] = "\
   --seed <n> : set random number seed to <n>\n\
";

static struct opt_s OPTIONS[] = {
  { "-a",        TRUE,  sqdARG_NONE },  
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-n",        TRUE,  sqdARG_INT},  
  { "-o",        TRUE,  sqdARG_STRING},
  { "-q",        TRUE,  sqdARG_NONE},  
  { "--seed",    FALSE, sqdARG_INT},
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv) 
{
  char            *hmmfile;	/* file to read HMMs from                  */
  FILE            *fp;          /* output file handle                      */
  HMMFILE         *hmmfp;       /* opened hmmfile for reading              */
  struct plan7_s  *hmm;         /* HMM to generate from                    */
  int              L;		/* length of a sequence                    */
  int              i;		/* counter over sequences                  */

  char            *ofile;       /* output sequence file                    */
  int              nseq;	/* number of seqs to sample                */
  int              seed;	/* random number generator seed            */
  int              be_quiet;	/* TRUE to silence header/footer           */
  int              do_alignment;/* TRUE to output in aligned format        */ 

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

  nseq         = 10;
  seed         = time ((time_t *) NULL);
  be_quiet     = FALSE;
  do_alignment = FALSE;  
  ofile        = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-a")     == 0) do_alignment = TRUE;
    else if (strcmp(optname, "-n")     == 0) nseq         = atoi(optarg); 
    else if (strcmp(optname, "-o")     == 0) ofile        = optarg;
    else if (strcmp(optname, "-q")     == 0) be_quiet     = TRUE;
    else if (strcmp(optname, "--seed") == 0) seed         = atoi(optarg);
    else if (strcmp(optname, "-h") == 0) 
      {
	Banner(stdout, banner);
	puts(usage);
	puts(experts);
	exit(0);
      }
  }
  if (argc - optind != 1)
    Die("Incorrect number of arguments.\n%s\n", usage);

  hmmfile = argv[optind++];

  sre_srandom(seed);

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

  /* Configure the HMM to shut off N,J,C emission: so we
   * do a simple single pass through the model.
   */
  Plan7NakedConfig(hmm);
  Plan7Renormalize(hmm);

  /*********************************************** 
   * Open the output file, or stdout
   ***********************************************/ 

   if (ofile == NULL) fp = stdout;
   else {
     if ((fp = fopen(ofile, "w")) == NULL)
       Die("Failed to open output file %s for writing", ofile);
   }
 
  /*********************************************** 
   * Show the options banner
   ***********************************************/

  if (! be_quiet) 
    {
      Banner(stdout, banner);
      printf("HMM file:             %s\n", hmmfile);
      printf("Random seed:          %d\n", seed);
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
    }

  /*********************************************** 
   * Do the work.
   * If we're generating an alignment, we have to collect
   * all our traces, then output. If we're generating unaligned
   * sequences, we can emit one at a time.
   ***********************************************/

  if (do_alignment)
    {
      struct p7trace_s **tr;        /* traces for aligned sequences            */
      char           **dsq;         /* digitized sequences                     */
      SQINFO          *sqinfo;      /* info about sequences (name/desc)        */
      char           **aseq;        /* sequence alignment                      */
      AINFO            ainfo;	    /* optional alignment info                 */
      float           *wgt;

      dsq    = MallocOrDie(sizeof(char *)             * nseq);
      tr     = MallocOrDie(sizeof(struct p7trace_s *) * nseq);
      sqinfo = MallocOrDie(sizeof(SQINFO)             * nseq);
      wgt    = MallocOrDie(sizeof(float)              * nseq);
      FSet(wgt, nseq, 1.0);

      for (i = 0; i < nseq; i++)
	{
	  EmitSequence(hmm, &(dsq[i]), &L, &(tr[i]));
	  sprintf(sqinfo[i].name, "seq%d", i+1);
	  sqinfo[i].len   = L;
	  sqinfo[i].flags = SQINFO_NAME | SQINFO_LEN;
	}

      P7Traces2Alignment(dsq, sqinfo, wgt, nseq, hmm->M, tr, FALSE, 
			 &aseq, &ainfo);

				/* Output the alignment */
      WriteSELEX(fp, aseq, &ainfo, 50);
      if (ofile != NULL && !be_quiet) printf("Alignment saved in file %s\n", ofile);

      /* Free memory
       */
      for (i = 0; i < nseq; i++) 
	{
	  P7FreeTrace(tr[i]);
	  free(dsq[i]);
	}
      FreeAlignment(aseq, &ainfo);
      free(sqinfo);
      free(dsq);
      free(wgt);
      free(tr);
    }
  else				/* unaligned sequence output */
    {
      struct p7trace_s *tr;         /* generated trace                        */
      char             *dsq;        /* digitized sequence                     */
      char             *seq;        /* alphabetic sequence                    */
      SQINFO            sqinfo;     /* info about sequence (name/len)         */

      for (i = 0; i < nseq; i++)
	{
	  EmitSequence(hmm, &dsq, &L, &tr);
	  sprintf(sqinfo.name, "seq%d", i+1);
	  sqinfo.len   = L;
	  sqinfo.flags = SQINFO_NAME | SQINFO_LEN;

	  seq = DedigitizeSequence(dsq, L);

	  WriteSeq(fp, kPearson, seq, &sqinfo);
	  
	  P7FreeTrace(tr);
	  free(dsq);
	  free(seq);
	}
    }

  if (ofile != NULL) fclose(fp);
  FreePlan7(hmm);
  SqdClean();

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) malloc_list(2, histid1, histid2);
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  return 0;
}

