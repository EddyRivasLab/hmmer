/************************************************************
 * @LICENSE@
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
   --mapali <f>  : include alignment in file <f> using map in HMM\n\
   --withali <f> : include alignment to (fixed) alignment in file <f>\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE   }, 
  { "-m", TRUE, sqdARG_NONE   } ,
  { "-o", TRUE, sqdARG_STRING },
  { "-q", TRUE, sqdARG_NONE   },
  { "--mapali",  FALSE, sqdARG_STRING },
  { "--withali", FALSE, sqdARG_STRING },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static void include_alignment(char *seqfile, struct plan7_s *hmm, int do_mapped,
			      char ***rseq, char ***dsq, SQINFO **sqinfo, 
			      struct p7trace_s ***tr, int *nseq);

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
  char *withali;                /* name of additional alignment file to align */
  char *mapali;                 /* name of additional alignment file to map   */
  
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
  withali   = NULL;
  mapali    = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-m")        == 0) matchonly= TRUE;
    else if (strcmp(optname, "-o")        == 0) outfile  = optarg;
    else if (strcmp(optname, "-q")        == 0) be_quiet = TRUE; 
    else if (strcmp(optname, "--mapali")  == 0) mapali   = optarg;
    else if (strcmp(optname, "--withali") == 0) withali  = optarg;
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
  * 
  * Currently hmmalign disallows the J state and
  * only allows one domain per sequence. To preserve
  * the S/W entry information, the J state is explicitly
  * disallowed, rather than calling a Plan7*Config() function.
  * this is a workaround in 2.1 for the 2.0.x "yo!" bug.
  ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  if (!HMMFileRead(hmmfp, &hmm)) 
    Die("Failed to read any HMMs from %s\n", hmmfile);
  HMMFileClose(hmmfp);
  if (hmm == NULL) 
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  hmm->xt[XTE][MOVE] = 1.;	      /* only 1 domain/sequence ("global" alignment) */
  hmm->xt[XTE][LOOP] = 0.;
  P7Logoddsify(hmm, TRUE);
				/* do we have the map we might need? */
  if (mapali != NULL && ! (hmm->flags & PLAN7_MAP))
    Die("HMMER: HMM file %s has no map; you can't use --mapali.", hmmfile);

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

  /* Include an aligned alignment, if desired.
   */
  if (mapali != NULL)
    include_alignment(mapali, hmm, TRUE, &rseq, &dsq, &sqinfo, &tr, &nseq);
  if (withali != NULL) 
    include_alignment(withali, hmm, FALSE, &rseq, &dsq, &sqinfo, &tr, &nseq);

  /* Turn traces into a multiple alignment
   */ 
  wgt = MallocOrDie(sizeof(float) * nseq);
  FSet(wgt, nseq, 1.0);
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
  free(rseq);
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


/* Function: include_alignment()
 * Date:     SRE, Sun Jul  5 15:25:13 1998 [St. Louis]
 *
 * Purpose:  Given the name of a multiple alignment file,
 *           align that alignment to the HMM, and add traces
 *           to an existing array of traces. If do_mapped
 *           is TRUE, we use the HMM's map file. If not,
 *           we use P7ViterbiAlignAlignment().
 *
 * Args:     seqfile  - name of alignment file
 *           hmm      - model to align to
 *           do_mapped- TRUE if we're to use the HMM's alignment map
 *           rsq      - RETURN: array of rseqs to add to
 *           dsq      - RETURN: array of dsq to add to
 *           sqinfo   - RETURN: array of SQINFO to add to
 *           tr       - RETURN: array of traces to add to
 *           nseq     - RETURN: number of seqs           
 *
 * Returns:  new, realloc'ed arrays for rsq, dsq, sqinfo, tr; nseq is
 *           increased to nseq+ainfo.nseq.
 */
void
include_alignment(char *seqfile, struct plan7_s *hmm, int do_mapped,
		  char ***rsq, char ***dsq, SQINFO **sqinfo, 
		  struct p7trace_s ***tr, int *nseq)
{
  int format;			/* format of alignment file */
  char **aseq;			/* aligned seqs             */
  char **newdsq;
  char **newrseq;
  AINFO ainfo;			/* info that goes with aseq */
  int   idx;			/* counter over aseqs       */
  struct p7trace_s *master;     /* master trace             */
  struct p7trace_s **addtr;     /* individual traces for aseq */

  if (! SeqfileFormat(seqfile, &format, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: 
      Die("Alignment file %s could not be opened for reading", seqfile);
      /*FALLTHRU*/ /* a white lie to shut lint up */
    case SQERR_FORMAT: 
    default:           
      Die("Failed to determine format of alignment file %s", seqfile);
    }
				/* read the alignment from file */
  if (! ReadAlignment(seqfile, format, &aseq, &ainfo))
    Die("Failed to read aligned sequence file %s", seqfile);
  for (idx = 0; idx < ainfo.nseq; idx++)
    s2upper(aseq[idx]);
				/* Verify checksums before mapping */
  if (do_mapped && GCGMultchecksum(aseq, ainfo.nseq) != hmm->checksum)
    Die("The checksums for alignment file %s and the HMM alignment map don't match.", 
	seqfile);
				/* Get a master trace */
  if (do_mapped) master = MasterTraceFromMap(hmm->map, hmm->M, ainfo.alen);
  else           master = P7ViterbiAlignAlignment(aseq, &ainfo, hmm);

				/* convert to individual traces */
  ImposeMasterTrace(aseq, ainfo.nseq, master, &addtr);
				/* add those traces to existing ones */
  *tr = MergeTraceArrays(*tr, *nseq, addtr, ainfo.nseq);
  
				/* additional bookkeeping: add to dsq, sqinfo */
  *rsq = ReallocOrDie((*rsq), sizeof(char *) * (*nseq + ainfo.nseq));
  DealignAseqs(aseq, ainfo.nseq, &newrseq);
  for (idx = *nseq; idx < *nseq + ainfo.nseq; idx++)
    (*rsq)[idx] = newrseq[idx - (*nseq)];
  free(newrseq);

  *dsq = ReallocOrDie((*dsq), sizeof(char *) * (*nseq + ainfo.nseq));
  DigitizeAlignment(aseq, &ainfo, &newdsq);
  for (idx = *nseq; idx < *nseq + ainfo.nseq; idx++)
    (*dsq)[idx] = newdsq[idx - (*nseq)];
  free(newdsq);
				/* unnecessarily complex, but I can't be bothered... */
  *sqinfo = ReallocOrDie((*sqinfo), sizeof(SQINFO) * (*nseq + ainfo.nseq));
  for (idx = *nseq; idx < *nseq + ainfo.nseq; idx++)
    SeqinfoCopy(&((*sqinfo)[idx]), &(ainfo.sqinfo[idx - (*nseq)]));
  
  *nseq = *nseq + ainfo.nseq;

				/* Cleanup */
  P7FreeTrace(master);
  FreeAlignment(aseq, &ainfo);
				/* Return */
  return;
}



