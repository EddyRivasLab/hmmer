/************************************************************
 * @LICENSE@
 ************************************************************/

/* hmmfetch.c
 * SRE, Wed Aug  5 14:26:51 1998 [St. Louis]
 * 
 * Recover a specific HMM file from an HMM database, using
 * a GSI index (created with hmmindex).
 * 
 * CVS $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"
#include "config.h"
#include "structs.h"
#include "funcs.h"
#include "version.h"

#include "globals.h"

static char banner[] = "hmmfetch -- retrieve specific HMM from an HMM database";

static char usage[] = "\
Usage: hmmfetch [-options] <hmmfile> <HMM name>\n\
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
  char    *key;			/* HMM name to retrieve            */
  HMMFILE *hmmfp;               /* opened hmm file pointer         */
  struct plan7_s *hmm;		/* a hidden Markov model           */

  char *optname;		/* name of option found by Getopt() */
  char *optarg;			/* argument found by Getopt()       */
  int   optind;		        /* index in argv[]                  */

  /***********************************************
   * Parse the command line
   ***********************************************/

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))
    {
      if (strcmp(optname, "-h") == 0)
	{
	  Banner(stdout, banner);
	  puts(usage);
	  puts(experts);
	  exit(0);
	}
    }

  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  hmmfile = argv[optind++];
  key     = argv[optind++];

  /***********************************************
   * Open HMM file, make sure GSI index exists
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("failed to open HMM file %s for reading.", hmmfile);
  if (hmmfp->gsi == NULL)
    Die("There is no GSI index for %s; you need to use hmmindex on it.", hmmfile);

  /***********************************************
   * find key in hmmfile; get HMM; show as ASCII
   ***********************************************/

  if (! HMMFilePositionByName(hmmfp, key))
    Die("No such hmm %s in HMM file %s\n", key, hmmfile);
  HMMFileRead(hmmfp, &hmm); 
  if (hmm == NULL) 
    Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);

  WriteAscHMM(stdout, hmm);

  FreePlan7(hmm);
  HMMFileClose(hmmfp);

  /***********************************************
   * Exit
   ***********************************************/

  SqdClean();
  return 0;
}


