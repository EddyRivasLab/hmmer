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

/* hmmcalibrate.c
 * SRE, Fri Oct 31 09:25:21 1997 [St. Louis]
 * 
 * Score an HMM against random sequence data sets;
 * set histogram fitting parameters.
 * 
 * RCS $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
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
  -l <n>      :       - : fix random sequence length at <n>\n\
  -m <x>      :     350 : set random seq length mean at <x>\n\
  -n <n>      :   10000 : set number of sampled seqs to <n>\n\
  -s <x>      :     350 : set random seq length std. dev to <x>\n\
\n\
  --histfile <f> :    - : save histograms to file <f>\n\
";

static struct opt_s OPTIONS[] = {
   { "-h", TRUE, sqdARG_NONE  },
   { "-l", TRUE, sqdARG_INT   },
   { "-m", TRUE, sqdARG_FLOAT },
   { "-n", TRUE, sqdARG_INT   },
   { "-s", TRUE, sqdARG_FLOAT },   
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
  sigset_t blocksigs;		/* list of signals to protect from */

  int     nsample;		/* number of random seqs to sample */
  int     seed;			/* random number seed              */
  int     fixedlen;		/* fixed length, or 0 if unused    */
  float   lenmean;		/* mean of length distribution     */
  float   lensd;		/* std dev of length distribution  */
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
  fixedlen = 0;
  lenmean  = 325.;
  lensd    = 200.;
  seed     = (int) time ((time_t *) NULL);
  histfile = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))
    {
      if      (strcmp(optname, "-l") == 0) fixedlen = atoi(optarg);
      else if (strcmp(optname, "-m") == 0) lenmean  = atof(optarg); 
      else if (strcmp(optname, "-n") == 0) nsample  = atoi(optarg); 
      else if (strcmp(optname, "-s") == 0) lensd    = atof(optarg); 
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

  /* Generate calibrated HMM(s) in a tmp file in the current
   * directory. When we're finished, we delete the original
   * HMM file and rename() this one. That way, the worst
   * effect of a catastrophic failure should be that we
   * leave a tmp file lying around, but the original HMM
   * file remains uncorrupted. tmpnam() doesn't work here,
   * because it'll put the file in /tmp and we won't
   * necessarily be able to rename() it from there.
   */
  tmpfile = MallocOrDie(strlen(hmmfile) + 5);
  strcpy(tmpfile, hmmfile);
  strcat(tmpfile, ".xxx");	/* could be more inventive here... */
  if (FileExists(tmpfile))
    Die("temporary file %s already exists; please delete it first", tmpfile);
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
  if (fixedlen) 
    printf("Length fixed to:          %d\n", fixedlen);
  else {
    printf("Length distribution mean: %.0f\n", lenmean);
    printf("Length distribution s.d.: %.0f\n", lensd);
  }
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

      max = -FLT_MAX;
      for (idx = 0; idx < nsample; idx++)
	{
				/* choose length of random sequence */
	  if (fixedlen) sqlen = fixedlen;
	  else do sqlen = (int) Gaussrandom(lenmean, lensd); while (sqlen < 1);

	  seq = RandomSequence(Alphabet, randomseq, Alphabet_size, sqlen);
	  dsq = DigitizeSequence(seq, sqlen);

	  score = Plan7Viterbi(dsq, sqlen, hmm, NULL);

	  AddToHistogram(hist, score);
	  if (score > max) max = score;

	  free(dsq);
	  free(seq);
	}

      /* Fit an EVD to the observed histogram.
       * The TRUE left-censors and fits only the right slope of the histogram.
       * The 9999. is an arbitrary high number that means we won't trim outliers
       * on the right.
       */
      if (! ExtremeValueFitHistogram(hist, TRUE, 9999.))
	{ 
	  fclose(outfp); 
	  remove(tmpfile);
	  Die("Failed to fit the histogram; maybe you set -n too small?");
	}

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

  /* Now, carefully remove original file and replace it
   * with the tmpfile. Note the protection from signals;
   * we wouldn't want a user to ctrl-C just as we've deleted
   * their HMM but before the new one is moved.
   */
  HMMFileClose(hmmfp);
  if (fclose(outfp)   != 0) PANIC;
  if (sigemptyset(&blocksigs) != 0) PANIC;
  if (sigaddset(&blocksigs, SIGINT) != 0) PANIC;
  if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   PANIC;
  if (remove(hmmfile) != 0)                            PANIC;
  if (rename(tmpfile, hmmfile) != 0)                   PANIC;
  if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) PANIC;

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


