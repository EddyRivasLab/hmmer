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
  char    *pgsifile;            /* primary key GSI file to write   */
  char    *sgsifile;		/* secondary key GSI file to write */
  HMMFILE *hmmfp;               /* opened hmm file pointer         */
  FILE    *pfp;                 /* open gsifile for writing primary keys   */
  FILE    *sfp;                 /* open gsifile for writing secondary keys */
  struct plan7_s     *hmm;      /* a hidden Markov model           */
  int     idx, nhmm;		/* counter over HMMs               */
  int     npri;			/* number of primary keys          */
  int     nsec;			/* number of secondary keys        */
  long    offset;		/* offset in HMM file              */
  struct gsikey_s *pkeylist;    /* list of primary keys            */
  struct gsikey_s *skeylist;    /* list of primary keys            */
  char    fname[GSI_KEYSIZE];

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
    Die("Primary GSI index already exists for %s. Please delete it first.", hmmfile);
  if (hmmfp->sgsi != NULL)
    Die("Secondary GSI index already exists for %s. Please delete it first.", hmmfile);
  
  pgsifile = MallocOrDie(strlen(hmmfile) + 5);
  strcpy(pgsifile, hmmfile);
  strcat(pgsifile, ".gsi");
  if (FileExists(pgsifile))
    Die("Primary GSI file %s already exists; please delete it first", pgsifile);  /* shouldn't happen */
  if ((pfp = fopen(pgsifile, "wb")) == NULL)
    Die("Primary GSI file %s couldn't be opened for writing", pgsifile); 

  sgsifile = MallocOrDie(strlen(hmmfile) + 6);
  strcpy(sgsifile, hmmfile);
  strcat(sgsifile, ".sgsi");
  if (FileExists(sgsifile))
    Die("Secondary GSI file %s already exists; please delete it first", sgsifile); /* shouldn't happen */
  if ((sfp = fopen(sgsifile, "wb")) == NULL)
    Die("Secondary GSI file %s couldn't be opened for writing", sgsifile); 

  pkeylist = MallocOrDie(sizeof(struct gsikey_s) * KEYBLOCK);
  skeylist = MallocOrDie(sizeof(struct gsikey_s) * KEYBLOCK);

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
  nhmm = npri = nsec = 0;
  while (HMMFileRead(hmmfp, &hmm)) 
    {	
      if (hmm == NULL) 
	Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);

				/* record name of HMM as the primary retrieval key */
      if (strlen(hmm->name) >= GSI_KEYSIZE )
	Warn("HMM name %s is too long to be GSI indexed; truncating it.", hmm->name);
      strncpy(pkeylist[npri].key, hmm->name, GSI_KEYSIZE-1);
      pkeylist[npri].key[GSI_KEYSIZE-1] = '\0';
      pkeylist[npri].filenum = 1;
      pkeylist[npri].offset  = offset;
      npri++;
      if (npri % KEYBLOCK == 0)	
	pkeylist = ReallocOrDie(pkeylist, sizeof(struct gsikey_s) * (npri + KEYBLOCK));

				/* record accession of HMM as a secondary retrieval key */
      if (hmm->flags & PLAN7_ACC) {
	if (strlen(hmm->acc) >= GSI_KEYSIZE) 
	  Warn("HMM accession %s is too long to be GSI indexed; truncating it\n", hmm->acc);
	strncpy(skeylist[nsec].key, hmm->acc, GSI_KEYSIZE-1);
	skeylist[nsec].key[GSI_KEYSIZE-1] = '\0';
	skeylist[nsec].filenum = 1;
	skeylist[nsec].offset  = offset;
	nsec++;
	if (nsec % KEYBLOCK == 0)	
	  skeylist = ReallocOrDie(skeylist, sizeof(struct gsikey_s) * (nsec + KEYBLOCK));
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
  qsort((void *) pkeylist, npri, sizeof(struct gsikey_s), gsikey_compare);
  qsort((void *) skeylist, nsec, sizeof(struct gsikey_s), gsikey_compare);
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

  GSIWriteHeader(pfp, 1, npri);
  GSIWriteFileRecord(pfp, fname, 1, 0); /* this line is unused, so doesn't matter */
  for (idx = 0; idx < npri; idx++)
    GSIWriteKeyRecord(pfp, pkeylist[idx].key, pkeylist[idx].filenum, pkeylist[idx].offset);
  if (fclose(pfp) != 0) PANIC;

  GSIWriteHeader(sfp, 1, nsec);
  GSIWriteFileRecord(sfp, fname, 1, 0); /* this line is unused, so doesn't matter */
  for (idx = 0; idx < nsec; idx++)
    GSIWriteKeyRecord(sfp, skeylist[idx].key, skeylist[idx].filenum, skeylist[idx].offset);
  if (fclose(sfp) != 0) PANIC;

  printf("Complete.\n");
  printf("Created GSI indexes for %d names and %d accessions for %d HMMs in %s.\n", 
	 npri, nsec, nhmm, hmmfile);

  /***********************************************
   * Exit
   ***********************************************/

  free(hmmtail);
  free(pgsifile);
  free(sgsifile);
  free(pkeylist);
  free(skeylist);
  SqdClean();
  return 0;
}


