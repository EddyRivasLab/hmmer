/* hmmfetch.c
 * Recover a specific HMM file from an HMM database, using
 * an SSI index (created with hmmindex).
 * 
 * SRE, Wed Aug  5 14:26:51 1998 [St. Louis]
 * SVN $Id$
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"

#include "plan7.h"
#include "structs.h"
#include "funcs.h"
#include "globals.h"

static char banner[] = "hmmfetch -- retrieve specific HMM from an HMM database";

static char usage[] = "\
Usage: hmmfetch [-options] <hmmfile> <HMM name>\n\
Available options are:\n\
  -h   : print short usage and version info, then exit\n\
  -n   : interpret <HMM name> instead as an HMM number (0..nhmm-1)\n\
";

static char experts[] = "\
";

static struct opt_s OPTIONS[] = {
   { "-h", TRUE, sqdARG_NONE  },
   { "-n", TRUE, sqdARG_NONE  },
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

  int   by_number;		/* fetch by number, not name        */
  int   nhmm;			/* hmm number */

  /***********************************************
   * Parse the command line
   ***********************************************/
  
  by_number = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))
    {
      if      (strcmp(optname, "-n") == 0) by_number = TRUE;
      else if (strcmp(optname, "-h") == 0)
	{
	  HMMERBanner(stdout, banner);
	  puts(usage);
	  puts(experts);
	  exit(0);
	}
    }

  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  hmmfile = argv[optind++];
  key     = argv[optind++];

  /***********************************************
   * Open HMM file, make sure SSI index exists
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("failed to open HMM file %s for reading.", hmmfile);
  if (hmmfp->ssi == NULL)
    Die("There is no SSI index for %s; you need to use hmmindex on it.", hmmfile);

  /***********************************************
   * find key in hmmfile; get HMM; show as ASCII
   ***********************************************/

  if (by_number) {
    if (! IsInt(key)) Die("%s does not appear to be a number.", key);
    nhmm = atoi(key);
    if (! HMMFilePositionByIndex(hmmfp, nhmm)) 
      Die("failed to position %s to HMM #%d", hmmfile, nhmm);
  } else {
    if (! HMMFilePositionByName(hmmfp, key))
      Die("No such hmm %s in HMM file %s\n", key, hmmfile);
  }

  if (! HMMFileRead(hmmfp, &hmm))
    Die("Unexpected end of HMM file");
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


/************************************************************
 * @LICENSE@
 ************************************************************/

