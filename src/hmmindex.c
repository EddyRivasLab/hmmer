/************************************************************
 * @LICENSE@
 ************************************************************/

/* hmmindex.c
 * SRE, Wed Aug  5 11:05:03 1998 [St. Louis]
 * 
 * Create an SSI index file for an HMM database.
 * 
 * CVS $Id$
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"
#include "ssi.h"
#include "structs.h"
#include "funcs.h"
#include "globals.h"


static char banner[] = "hmmindex -- create SSI index for an HMM database";

static char usage[] = "\
Usage: hmmindex [-options] <hmmfile>\n\
Available options are:\n\
  -h             : print short usage and version info, then exit\n\
";

static char experts[] = "\
";

static struct opt_s OPTIONS[] = {
   { "-h", TRUE, sqdARG_NONE  },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char    *hmmfile;             /* HMM file to open                */
  SSIINDEX *ssi;                /* SSI index in memory             */
  char    *ssifile;             /* name of SSI index on disk       */
  HMMFILE *hmmfp;               /* opened hmm file pointer         */
  struct plan7_s     *hmm;      /* a hidden Markov model           */
  int     nhmm;		        /* counter over HMMs               */
  int     npri, nsec;		/* # of names, accessions          */
  int     fh;			/* file handle                     */
  int     status;		/* return status from SSI call     */

  char *optname;		/* name of option found by Getopt() */
  char *optarg;			/* argument found by Getopt()       */
  int   optind;		        /* index in argv[]                  */

  /***********************************************
   * Parse the command line
   ***********************************************/

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))
    {
      if (strcmp(optname, "-h")        == 0)
	{
	  HMMERBanner(stdout, banner);
	  puts(usage);
	  puts(experts);
	  exit(0);
	}
    }

  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  hmmfile = argv[optind++];

  /***********************************************
   * Open our input HMM file, make sure all is well with the output SSI filename
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, NULL)) == NULL)
    Die("failed to open HMM file %s for reading.", hmmfile);
  if (hmmfp->ssi != NULL)
    Die("SSI index already exists for %s.\nPlease delete it first.", hmmfile);
  
  ssifile = MallocOrDie(strlen(hmmfile) + 5);
  sprintf(ssifile, "%s%s", hmmfile, ".ssi");
  if (FileExists(ssifile))   /* shouldn't happen */
    Die("An SSI file %s already exists; please delete it first", ssifile);

  if ((ssi = SSICreateIndex(hmmfp->mode)) == NULL)
    Die("Failed to initialize the SSI index structure");
  if (SSIAddFileToIndex(ssi, hmmfile, hmmfp->is_binary, &fh) != 0)
    Die("SSIAddFileToIndex() failed");

  /*********************************************** 
   * Show the banner
   ***********************************************/

  HMMERBanner(stdout, banner);
  printf("HMM file:                 %s\n", hmmfile);
  if (hmmfp->mode == SSI_OFFSET_I64) 
    printf("Index file mode:          64-bit (large HMM file)\n");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");

  /***********************************************
   * Get offsets and names for every model; store in keylist
   ***********************************************/

  printf("Determining offsets for %s, please be patient...\n", hmmfile);

  nhmm = npri = nsec = 0;
  while (HMMFileRead(hmmfp, &hmm)) 
    {	
      if (hmm == NULL) 
	Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);

				/* record name of HMM as the primary retrieval key */
      status = SSIAddPrimaryKeyToIndex(ssi, hmm->name, fh, &(hmmfp->offset), NULL, 0);
      if (status != 0) Die("SSIAddPrimaryKeyToIndex() failed");
      npri++;

				/* record accession of HMM as a secondary retrieval key */
      if (hmm->flags & PLAN7_ACC) {
	status = SSIAddSecondaryKeyToIndex(ssi, hmm->acc, hmm->name);
	if (status != 0) Die("SSIAddSecondaryKeyToIndex() failed");
	nsec++;
      }

      nhmm++;
      FreePlan7(hmm);
    }
  HMMFileClose(hmmfp);

  /***********************************************
   * Output the SSI file
   ***********************************************/

  status = SSIWriteIndex(ssifile, ssi);
  if (status != 0) Die("SSIWriteIndex() failed");

  printf("Complete.\n");
  printf("HMM file:       %s\n", hmmfile);
  printf("SSI index:      %s\n", ssifile);
  printf("# of HMMS:      %d\n", nhmm);
  printf("HMM names:      %d\n", npri);
  printf("HMM accessions: %d\n", nsec);


  /***********************************************
   * Exit
   ***********************************************/

  free(ssifile);
  SSIFreeIndex(ssi);
  SqdClean();
  return 0;
}


