/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1997 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmmcalibrate.c
 * SRE, Fri Oct 31 09:25:21 1997 [St. Louis]
 * 
 * Score an HMM against random sequence data sets;
 * set histogram fitting parameters.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>

#include "squid.h"
#include "config.h"
#include "structs.h"
#include "funcs.h"
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#include "globals.h"

static char banner[] = "hmmcalibrate -- calibrate HMM search statistics";

static char usage[] = "\
Usage: hmmcalibrate [-options] <hmmfile>\n\
Option        : Default : Description\n\
  -h          :       - : print short usage and version info, then exit\n\
  -l <x>      :     2.0 : set random sequence length multiplier to <x>\n\
  -n <n>      :   10000 : set number of sampled seqs to <n>\n\
\n\
  --histfile <f>        : save histograms to file <f>\n\
";

static struct opt_s OPTIONS[] = {
   { "-h", TRUE, sqdARG_NONE   },
   { "-l", TRUE, sqdARG_FLOAT  },
   { "-n", TRUE, sqdARG_INT    },
   { "--histfile", FALSE, sqdARG_STRING },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char    *hmmfile;             /* HMM file to open                */
  char    *tmpfile;             /* temporary calibrated HMM file   */
  HMMFILE *hmmfp;               /* opened hmm file pointer         */
  FILE    *outfp;               /* for writing HMM(s) into tmpfile */
  char    *mode;                /* write mode, "w" or "wb"         */
  struct plan7_s     *hmm;      /* the hidden Markov model         */
  struct histogram_s *hist;	/* score histogram                 */
  int     idx;			/* counter over sequences          */
  char   *seq;			/* a random sequence               */
  char   *dsq;			/* seq, digitized for alignment    */
  float   randomseq[MAXABET];	/* random sequence model           */
  float   p1;			/* random sequence model p1        */
  float   score;		/* score of an alignment           */
  float   max;			/* maximum score                   */
  int     sqlen;		/* length of sampled sequences     */

  int     nsample;		/* number of random seqs to sample */
  int     seed;			/* random number seed              */
  float   lenmult;		/* length multiplier               */
  char   *histfile;             /* histogram save file             */
  FILE   *hfp;                  /* open file pointer for histfile  */

  char *optname;		/* name of option found by Getopt() */
  char *optarg;			/* argument found by Getopt()       */
  int   optind;		        /* index in argv[]                  */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /***********************************************
   * Parse the command line
   ***********************************************/

  nsample  = 10000;
  lenmult  = 2.0;
  seed     = (int) time ((time_t *) NULL);
  histfile = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))
    {
      if      (strcmp(optname, "-l") == 0) lenmult = atof(optarg);
      else if (strcmp(optname, "-n") == 0) nsample = atoi(optarg); 
      else if (strcmp(optname, "--histfile") == 0) histfile = optarg;
      else if (strcmp(optname, "-h") == 0)
	{
	  Banner(stdout, banner);
	  puts(usage);
	  exit(0);
	}
    }

  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  hmmfile = argv[optind++];

  sre_srandom(seed);

  /***********************************************
   * Open our i/o file pointers, make sure all is well
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("failed to open HMM file %s for reading.", hmmfile);

				/* generate a tmp file name in hmmfile directory */
  tmpfile = MallocOrDie(strlen(hmmfile) + 5);
  strcpy(tmpfile, hmmfile);
  strcat(tmpfile, ".tmp");
  if (FileExists(tmpfile))
    Die("temporary file %s already exists, can't open for writing", tmpfile);
  if (hmmfp->is_binary) mode = "wb";
  else                  mode = "w"; 
  if ((outfp = fopen(tmpfile, mode)) == NULL)
    Die("temporary file %s couldn't be opened for writing", tmpfile); 

				/* histogram file */
  if (histfile != NULL)
    {
      if ((hfp = fopen(histfile, "w")) == NULL)
	Die("Failed to open histogram save file %s for writing\n", histfile);
    }

  /*********************************************** 
   * Show the banner
   ***********************************************/

  Banner(stdout, banner);
  printf("HMM file:                 %s\n", hmmfile);
  printf("Length multiplier:        %.1f\n", lenmult);
  printf("Number of samples:        %d\n", nsample);
  printf("random seed:              %d\n", seed);
  printf("histogram(s) saved to:    %s\n",
	 histfile != NULL ? histfile : "[not saved]");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");

  /***********************************************
   * Calibrate each model in turn
   ***********************************************/

  while (HMMFileRead(hmmfp, &hmm)) 
    {	
      if (hmm == NULL) 
	Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);
      Plan7Logoddsify(hmm);
				/* we could use the null model in the HMM? */
      P7DefaultNullModel(randomseq, &p1);
      
      hist = AllocHistogram(-200, 200, 100);

				/* choose length of random sequences */
      sqlen = lenmult * hmm->M;

      max = -FLT_MAX;
      for (idx = 0; idx < nsample; idx++)
	{
	  seq = RandomSequence(Alphabet, randomseq, Alphabet_size, sqlen);
	  dsq = DigitizeSequence(seq, sqlen);

	  score = Plan7Viterbi(dsq, sqlen, hmm, NULL);

	  AddToHistogram(hist, score);
	  if (score > max) max = score;

	  free(dsq);
	  free(seq);
	}

      ExtremeValueFitHistogram(hist, 200.);

      /* Set HMM EVD parameters 
       */
      hmm->mu      = hist->param[EVD_MU];
      hmm->lambda  = hist->param[EVD_LAMBDA];
      hmm->flags  |= PLAN7_STATS;

      /* Save HMM to tmpfile
       */
      if (hmmfp->is_binary) WriteBinHMM(outfp, hmm);
      else                  WriteAscHMM(outfp, hmm); 

      /* Output results
       */
      printf("HMM    : %s\n", hmm->name);
      printf("sqlen  : %12d\n", sqlen);
      printf("mu     : %12f\n", hmm->mu);
      printf("lambda : %12f\n", hmm->lambda);
      printf("max    : %12f\n", max);
      printf("//\n");

      if (histfile != NULL) 
	{
	  fprintf(hfp, "HMM: %s\n", hmm->name);
	  PrintASCIIHistogram(hfp, hist);
	  fprintf(hfp, "//\n");
	}

      FreePlan7(hmm);
      FreeHistogram(hist);
    }

  /* If all is well, remove original file and replace it
   * with the tmpfile. This code should be protected from
   * interrupt signals.
   */
  HMMFileClose(hmmfp);
  if (fclose(outfp)   != 0)          Die("fclose() failed on %s", tmpfile);
  if (remove(hmmfile) != 0)          Die("failed to remove %s",   hmmfile);
  if (rename(tmpfile, hmmfile) != 0) Die("failed to rename %s to %s", tmpfile, hmmfile); 

  /***********************************************
   * Exit
   ***********************************************/
  free(tmpfile);
  SqdClean();
#ifdef MEMDEBUG
  current_size = malloc_size(&histid2);
  if (current_size != orig_size)
    malloc_list(2, histid1, histid2);
  else
    fprintf(stderr, "[No memory leaks]\n");
#endif

  return 0;
}


