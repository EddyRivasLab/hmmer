/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmmconvert.c
 * SRE, Thu Oct 30 08:56:22 1997; St. Louis
 * 
 * main() for converting between HMM file formats, and
 * for converting HMMs to other software formats like GCG profiles.
 * 
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

static char banner[] = "hmmconvert - convert between profile HMM file formats";

static char usage[]  = "\
Usage: hmmconvert [-options] <old hmm file> <new hmm file>\n\
  Available options are:\n\
   -h        : help; print brief help on version and usage\n\
\n\
   -a        : convert to HMMER ASCII file (the default)\n\
   -b        : convert to HMMER binary file\n\
   -p        : convert to GCG Profile .prf format\n\
   -P        : convert to Compugen extended .eprf profile format\n\
\n\
   -A        : append mode; append to <new hmm file>\n\
   -F        : force mode; allow overwriting of existing files\n\
";

static char experts[] = "\
\n";


static struct opt_s OPTIONS[] = {
  { "-a",        TRUE,  sqdARG_NONE },
  { "-b",        TRUE,  sqdARG_NONE },
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-p",        TRUE,  sqdARG_NONE },
  { "-A",        TRUE,  sqdARG_NONE },
  { "-F",        TRUE,  sqdARG_NONE },
  { "-P",        TRUE,  sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int 
main(int argc, char **argv)
{
  char    *infile;              /* name of input HMM file                   */
  char    *outfile;             /* name of output HMM file                  */
  HMMFILE *infp;                /* input HMM file ptr                       */
  FILE    *outfp;               /* output HMM file ptr                      */
  char    *mode;                /* mode to open file in                     */
  struct plan7_s *hmm;          /* a profile HMM structure                  */
  int      nhmm;		/* number of HMMs converted                 */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */

  int   do_append;		/* TRUE to append to existing outfile       */
  int   do_force;		/* TRUE to allow overwriting */
  enum hmmfmt_e { P7ASCII, P7BINARY, GCGPROFILE, BICPROFILE } 
      outfmt;			/* output format */
      


#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/

  outfmt    = P7ASCII; 
  do_append = FALSE;
  do_force  = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-a") == 0) { outfmt    = P7ASCII;    }
    else if (strcmp(optname, "-b") == 0) { outfmt    = P7BINARY;   } 
    else if (strcmp(optname, "-p") == 0) { outfmt    = GCGPROFILE; } 
    else if (strcmp(optname, "-A") == 0) { do_append = TRUE;       } 
    else if (strcmp(optname, "-F") == 0) { do_force  = TRUE;       } 
    else if (strcmp(optname, "-P") == 0) { outfmt    = BICPROFILE; } 
    else if (strcmp(optname, "-h") == 0) {
      Banner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(0);
    }
  }
  if (argc - optind != 2)
    Die("Incorrect number of arguments.\n%s\n", usage);

  infile  = argv[optind++];
  outfile = argv[optind++]; 

  /*********************************************** 
   * Open input HMM database (might be in HMMERDB or current directory)
   ***********************************************/

  if ((infp = HMMFileOpen(infile, "HMMERDB")) == NULL)
    Die("Failed to open HMM database %s\n%s", infile, usage);

  /*********************************************** 
   * Open output HMM file
   ***********************************************/
  
  if (do_append) 
    {		/* If we're appending to a file, it needs to be Plan7 format */
      HMMFILE *test;

      if (FileExists(outfile)) { 
	test = HMMFileOpen(outfile, NULL);
	if (test == NULL) 
	  Die("%s not an HMM file; I refuse to append to it; using stdout instead",
	      outfile);
	HMMFileClose(test);
      }
      switch (outfmt) {
      case P7ASCII:    mode = "a";  break;
      case P7BINARY:   mode = "ab"; break;
      case GCGPROFILE: Die("You cannot append GCG profiles");
      case BICPROFILE: Die("You cannot append Compugen extended profiles");
      default:         Die("unexpected format");
      }
    }
  else
    {		/* else, we're writing a new file */
      if (! do_force && FileExists(outfile))
	Die("Output HMM file %s already exists. Please rename or delete it.", outfile); 
      switch (outfmt) {
      case P7ASCII:    mode = "w";  break;
      case P7BINARY:   mode = "wb"; break;
      case GCGPROFILE: mode = "w";  break;
      case BICPROFILE: mode = "w";  break;
      default:         Die("unexpected format");
      }
    }
  if ((outfp = fopen(outfile, mode)) == NULL) 
    Die("Failed to open output file %s for writing", outfile);

  /*********************************************** 
   * Show the banner
   ***********************************************/

  Banner(stdout, banner);
  printf(   "Input HMM file:           %s\n", infile);
  printf(   "Output HMM file:          %s\n", outfile);
  printf(   "Converting to:            ");
  switch (outfmt) {
  case P7ASCII:    puts("HMMER Plan7 ASCII");      break;
  case P7BINARY:   puts("HMMER Plan7 binary");     break;
  case GCGPROFILE: puts("GCG Profile .prf");       break;
  case BICPROFILE: puts("Compugen .eprf profile"); break;
  default:         Die("unexpected fault");
  }
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");

  /*********************************************** 
   * Do the conversion
   ***********************************************/

  nhmm = 0;
  while (HMMFileRead(infp, &hmm)) {
    if (hmm == NULL) 
      Die("HMM file %s may be corrupt or in incorrect format; parse failed", infile);

    switch(outfmt) {
    case P7ASCII:    WriteAscHMM(outfp, hmm);         break;
    case P7BINARY:   WriteBinHMM(outfp, hmm);         break;
    case GCGPROFILE: WriteProfile(outfp, hmm, FALSE); break;
    case BICPROFILE: WriteProfile(outfp, hmm, TRUE);  break;
    default:         Die("unexpected format");
    }

    printf(" - converted %s\n", hmm->name);
    FreePlan7(hmm);
    nhmm++;
  }
  printf("\n%d HMM(s) converted and written to %s\n", nhmm, outfile);

  /*********************************************** 
   * Clean-up and exit.
   ***********************************************/

  HMMFileClose(infp);
  fclose(outfp);
  SqdClean();
#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) malloc_list(2, histid1, histid2);
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  return EXIT_SUCCESS;
}
