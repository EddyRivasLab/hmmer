/************************************************************
 * @LICENSE@
 ************************************************************/

/* hmmindex.c
 * SRE, Wed Aug  5 11:05:03 1998 [St. Louis]
 * 
 * Create a GSI index file for an HMM database.
 * 
 * RCS $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"
#include "config.h"
#include "structs.h"
#include "funcs.h"
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#include "globals.h"

static char banner[] = "hmmindex -- create GSI index for an HMM database";

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


struct gsikey_s {
  char        key[GSI_KEYSIZE];
  int         filenum;
  long        offset;
};
#define KEYBLOCK 100		/* malloc blocks of keys */

static int
gsikey_compare(const void *el1, const void *el2)
{
  struct gsikey_s *key1;
  struct gsikey_s *key2;
  key1 = (struct gsikey_s *) el1;
  key2 = (struct gsikey_s *) el2;

  return strcmp(key1->key, key2->key);
}

int
main(int argc, char **argv)
{
  char    *hmmfile;             /* HMM file to open                */
  char    *hmmtail;             /* HMMfile without directory path  */
  char    *gsifile;             /* GSI file to write               */
  HMMFILE *hmmfp;               /* opened hmm file pointer         */
  FILE    *outfp;               /* open gsifile for writing        */
  struct plan7_s     *hmm;      /* a hidden Markov model           */
  int     idx, nhmm;		/* counter over HMMs               */
  int     nkeys;		/* number of keys                  */
  long    offset;		/* offset in HMM file              */
  struct gsikey_s *keylist;     /* list of keys                    */
  char    fname[GSI_KEYSIZE];

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

  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  hmmfile = argv[optind++];

  /***********************************************
   * Open our i/o file pointers, make sure all is well
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, NULL)) == NULL)
    Die("failed to open HMM file %s for reading.", hmmfile);
  if (hmmfp->gsi != NULL)
    Die("GSI index already exists for %s. Please delete it first.", hmmfile);
  
  gsifile = MallocOrDie(strlen(hmmfile) + 5);
  strcpy(gsifile, hmmfile);
  strcat(gsifile, ".gsi");
  if (FileExists(gsifile))
    Die("GSI file %s already exists; please delete it first", gsifile);	/* shouldn't happen */
  if ((outfp = fopen(gsifile, "wb")) == NULL)
    Die("GSI file %s couldn't be opened for writing", gsifile); 

  keylist = MallocOrDie(sizeof(struct gsikey_s) * KEYBLOCK);

  /*********************************************** 
   * Show the banner
   ***********************************************/

  Banner(stdout, banner);
  printf("HMM file:                 %s\n", hmmfile);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");

  /***********************************************
   * Get offsets and names for every model; store in keylist
   ***********************************************/

  printf("Determining offsets for %s, please be patient...\n", hmmfile);
  offset = ftell(hmmfp->f);
  nhmm = nkeys = 0;
  while (HMMFileRead(hmmfp, &hmm)) 
    {	
      if (hmm == NULL) 
	Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);

      if (strlen(hmm->name) >= GSI_KEYSIZE )
	Warn("HMM %s name is too long to index in a GSI file; truncating it.", hmm->name);
				/* record name of HMM as a retrieval key */
      if (strlen(hmm->name) >= GSI_KEYSIZE) 
	Warn("HMM name %s is too long to be indexed: must be < %d\n", GSI_KEYSIZE);
      strncpy(keylist[nkeys].key, hmm->name, GSI_KEYSIZE-1);
      keylist[nkeys].key[GSI_KEYSIZE-1] = '\0';
      keylist[nkeys].filenum = 1;
      keylist[nkeys].offset  = offset;
      nkeys++;
      if (nkeys % KEYBLOCK == 0)	
	keylist = ReallocOrDie(keylist, sizeof(struct gsikey_s) * (nkeys + KEYBLOCK));

				/* record accession of HMM as a retrieval key */
      if (hmm->flags & PLAN7_ACC) {
	if (strlen(hmm->acc) >= GSI_KEYSIZE) 
	  Warn("HMM accession %s is too long to be indexed: must be < %d\n", GSI_KEYSIZE);
	strncpy(keylist[nkeys].key, hmm->acc, GSI_KEYSIZE-1);
	keylist[nkeys].key[GSI_KEYSIZE-1] = '\0';
	keylist[nkeys].filenum = 1;
	keylist[nkeys].offset  = offset;
	nkeys++;
	if (nkeys % KEYBLOCK == 0)	
	  keylist = ReallocOrDie(keylist, sizeof(struct gsikey_s) * (nkeys + KEYBLOCK));
      }

      offset = ftell(hmmfp->f);
      nhmm++;
      FreePlan7(hmm);
    }
  HMMFileClose(hmmfp);

  /***********************************************
   * Sort the keylist
   ***********************************************/

  printf("Sorting keys... \n");
  qsort((void *) keylist, nkeys, sizeof(struct gsikey_s), gsikey_compare);
  SQD_DPRINTF1(("(OK, done with qsort)\n"));

  /***********************************************
   * Output the GSI file
   ***********************************************/

  hmmtail = FileTail(hmmfile, FALSE);
  if (strlen(hmmtail) >= GSI_KEYSIZE)
    {
      Warn("HMM file name length is >= %d char. Truncating.", GSI_KEYSIZE);
      hmmtail[GSI_KEYSIZE-1] = '\0';
    }
  strcpy(fname, hmmtail);

  GSIWriteHeader(outfp, 1, nkeys);
  GSIWriteFileRecord(outfp, fname, 1, 0); /* this line is unused, so doesn't matter */
  for (idx = 0; idx < nkeys; idx++)
    GSIWriteKeyRecord(outfp, keylist[idx].key, keylist[idx].filenum, keylist[idx].offset);

  printf("Complete.\n");
  printf("GSI %s indexes %d keys (names+accessions) for %d HMMs in %s.\n", 
	 gsifile, nkeys, nhmm, hmmfile);

  /***********************************************
   * Exit
   ***********************************************/

  free(hmmtail);
  free(gsifile);
  free(keylist);
  if (fclose(outfp) != 0) PANIC;
  SqdClean();

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size)
    malloc_list(2, histid1, histid2);
  else
    fprintf(stderr, "[No memory leaks]\n");
#endif

  return 0;
}


