/* Construct a training alignment/test sequences set from an MSA.
 * 
 * This procedure is used in constructing our internal RMARK and PMARK
 * benchmarks.
 * 
 * SRE, Thu Mar 27 11:05:38 2008 [Janelia]
 * SVN $Id$
 */

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_vectorops.h"


static char banner[] = "show summary statistics for a multiple sequence alignment file";
static char usage[]  = "[options] <msafile> <concatdb>\n\
The <msafile> must be in Stockholm format; it can be a multi-MSA file.";

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                                   docgroup */
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              0 },
  { "-1",       eslARG_REAL, "0.20", NULL,"0<x<=1.0",NULL,NULL,NULL,         "require all test seqs to have < x id to training",        0 },
  { "-2",       eslARG_REAL, "0.50", NULL,"0<x<=1.0",NULL,NULL,NULL,         "require all test seqs to have < x id to each other",      0 },
  { "-F",       eslARG_REAL, "0.70", NULL,"0<x<=1.0",NULL,NULL,NULL,         "filter out seqs <x*average length",                       0 },
  { "-N",       eslARG_INT,"75438310",NULL,NULL, NULL, NULL, NULL,           "specify # residues in <concatdb>",                        0 },
  { "--amino",  eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--dna,--rna",    "<msafile> contains protein alignments",                   0 },
  { "--dna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--rna",  "<msafile> contains DNA alignments",                       0 },
  { "--rna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--dna",  "<msafile> contains RNA alignments",                       0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

struct cfg_s {
  ESL_ALPHABET   *abc;          /* biological alphabet             */
  ESL_RANDOMNESS *r;            /* random number generator         */
  double          fragfrac;	/* seqs less than x*avg length are removed from alignment  */
  double          idthresh1;	/* fractional identity threshold for train/test split      */
  double          idthresh2;	/* fractional identity threshold for selecting test seqs   */
  FILE           *dbfp;		/* open concatdb file, source of nonhomologous seq frags   */
  int             N;		/* # of residues in concatdb                               */
};

static void
cmdline_failure(char *argv0, char *format, ...)
{
  va_list argp;
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  esl_usage(stdout, argv0, usage);
  printf("\nTo see more help on available options, do %s -h\n\n", argv0);
  exit(1);
}

static void
cmdline_help(char *argv0, ESL_GETOPTS *go) 
{
  esl_banner(stdout, argv0, banner);
  esl_usage (stdout, argv0, usage);
  puts("\n where general options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  exit(0);
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = NULL;	/* command line configuration      */
  struct cfg_s  cfg;     	/* application configuration       */
  char         *alifile = NULL;	/* alignment file name             */
  char         *dbfile  = NULL;	/* name of concatdb file           */
  int           fmt;		/* format code for alifile         */
  ESL_MSAFILE  *afp     = NULL;	/* open alignment file             */
  ESL_MSA      *origmsa = NULL;	/* one multiple sequence alignment */
  ESL_MSA      *msa     = NULL;	/* MSA after frags are removed     */
  int           status;		/* easel return code               */
  int           nali;		/* number of alignments read       */
  int           i;		/* counter over seqs               */
  int           rlen;		/* a raw (unaligned) seq length    */
  int           small, large;	/* smallest, largest sequence      */
  uint64_t      nres;		/* total # of residues in msa      */
  double        avgid;		/* average fractional pair id      */
  int    max_comparisons;       /* maximum # comparisons for avg id */
  
  /* Parse command line */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h"))                    cmdline_help(argv[0], go);
  if (esl_opt_ArgNumber(go)                  != 2)     cmdline_failure(argv[0], "Incorrect number of command line arguments\n");
  alifile = esl_opt_GetArg(go, 1);
  dbfile  = esl_opt_GetArg(go, 2);
  fmt     = eslMSAFILE_STOCKHOLM;

  /* Set up the configuration structure shared amongst functions here */
  cfg.r         = esl_randomness_Create(42);
  cfg.abc       = NULL;		          /* until we open the MSA file, below */
  cfg.fragfrac  = esl_opt_GetReal(go, "-F");
  cfg.idthresh1 = esl_opt_GetReal(go, "-1");
  cfg.idthresh2 = esl_opt_GetReal(go, "-2");
  cfg.dbfp      = fopen(dbfile, "rb");
  cfg.N         = esl_opt_GetInteger(go, "-N");

  if (cfg.dbfp == NULL) esl_fatal("Failed to open concatdb file %s for reading\n", dbfile);

  /* Open the MSA file */
  status = esl_msafile_Open(alifile, fmt, NULL, &afp);
  if      (status == eslENOTFOUND) esl_fatal("Alignment file %s doesn't exist or is not readable\n", alifile);
  else if (status == eslEFORMAT)   esl_fatal("Couldn't determine format of alignment %s\n", alifile);
  else if (status != eslOK)        esl_fatal("Alignment file open failed with error %d\n", status);

  if      (esl_opt_GetBoolean(go, "--amino"))   cfg.abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))     cfg.abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))     cfg.abc = esl_alphabet_Create(eslRNA);
  else {
    int type;
    status = esl_msafile_GuessAlphabet(afp, &type);
    if (status == eslEAMBIGUOUS)    esl_fatal("Failed to guess the bio alphabet used in %s.\nUse --dna, --rna, or --amino option to specify it.", alifile);
    else if (status == eslEFORMAT)  esl_fatal("Alignment file parse failed: %s\n", afp->errbuf);
    else if (status == eslENODATA)  esl_fatal("Alignment file %s is empty\n", alifile);
    else if (status != eslOK)       esl_fatal("Failed to read alignment file %s\n", alifile);
    abc = esl_alphabet_Create(type);
  }
  esl_msafile_SetDigital(afp, abc);

  while ((status = esl_msa_Read(afp, &origmsa)) == eslOK)
    {
      remove_fragments(&cfg, origmsa, &msa);
      separate_sets   (&cfg, msa, &trainmsa, &teststack);
      esl_stack_Shuffle(cfg.r, teststack);


      esl_msa_Destroy(trainmsa);
      esl_msa_Destroy(origmsa);
      esl_msa_Destroy(msa);
    }
  if      (status == eslEFORMAT)  esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", 
					    afp->linenumber, afp->fname, afp->errbuf, afp->buf);	
  else if (status != eslEOF)      esl_fatal("Alignment file read failed with error code %d\n", status);
  else if (nali   == 0)           esl_fatal("No alignments found in file %s\n", alifile);

  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  esl_msafile_Close(afp);
  esl_getopts_Destroy(go);
  return 0;
}
      
  
/* Step 1. Label all sequence fragments < fragfrac of average raw length */
static int
remove_fragments(struct cfg_s *cfg, ESL_MSA *msa, ESL_MSA **ret_filteredmsa)
{
  int     *useme    = NULL;
  double   len      = 0.0;
  int      status;

  for (i = 0; i < msa->nseq; i++) 
    len += esl_abc_dsqrlen(msa->abc, msa->ax[i]);
  len *= cfg->fragfrac / (double) msa->nseq;

  ESL_ALLOC(useme, sizeof(int) * msa->nseq);
  for (i = 0; i < msa->nseq; i++) 
    useme[i] = (esl_abc_dsqrlen(msa->abc, msa->ax[i]) < len) ? 0 : 1;

  return esl_msa_SequenceSubset(msa, useme, ret_filteredmsa);

 ERROR:
  if (useme != NULL) free(useme); 
  *ret_filteredmsa = NULL;
  return status;
}

/* Step 2. Extract the training set and test set.
 */
static int
separate_sets(struct cfg_s *cfg, ESL_MSA *msa, ESL_MSA *ret_trainmsa, ESL_STACK *ret_teststack)
{      
  ESL_MSA   *trainmsa  = NULL;
  ESL_MSA   *test_msa  = NULL;
  ESL_STACK *teststack = NULL;
  int *assignment = NULL;
  int *nin        = NULL;
  int *useme      = NULL;
  int  nc         = 0;
  int  c;
  int  ctrain;			/* index of the cluster that becomes the training alignment */
  int  ntrain;			/* number of seqs in the training alignment */
  int  nskip;
  int  i;
  int  status;

  if ((teststack = esl_stack_PCreate()) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(useme, sizeof(int) * msa->nseq);

  if ((status = esl_msacluster_SingleLinkage(msa, cfg->idthresh1, &assignment, &nin, &nc)) != eslOK) goto ERROR;
  ctrain = esl_vec_IArgMax(nin, nc);
  ntrain = esl_vec_IMax(nin, nc);

  for (i = 0; i < msa->nseq; i++) useme[i] = (!is_frag[i] && assignment[i] == ctrain) ? 1 : 0;
  if ((status = esl_msa_SequenceSubset(msa, useme, &trainmsa)) != eslOK) goto ERROR;

  free(nin);         nin        = NULL;
  free(assignment);  assignment = NULL;

  /* Put all the other sequences into an MSA of their own; from these, we'll
   * choose test sequences.
   */
  for (i = 0; i < msa->nseq; i++) useme[i] = (!is_frag[i] && assignment[i] != ctrain) ? 1 : 0;
  if ((status = esl_msa_SequenceSubset(msa, useme, &test_msa))                        != eslOK) goto ERROR;
  if ((status = esl_msacluster_SingleLinkage(msa, cfg->idthresh2, &assignment, &nin, &nc)) != eslOK) goto ERROR;

  for (c = 0; c < nc; c++)
    {
      nskip = esl_rnd_Choose(cfg->r, nin[c]); /* pick a random seq in this cluster to be the test. */
      for (i=0; nskip >= 0 && i < msa->nseq; i++)
	if (assignment[i] == c) nskip--;
      
      esl_sq_FetchFromMSA(msa, i, &sq);
      esl_stack_PPush(teststack, (void *) sq);
    }

  esl_msa_Destroy(test_msa);
  free(useme);
  free(nin);
  free(assignment);

  *ret_trainmsa  = msa;
  *ret_teststack = teststack;
  return eslOK;

 ERROR:
  if (useme      != NULL) free(useme);
  if (assignment != NULL) free(assignment);
  if (nin        != NULL) free(nin);
  esl_msa_Destroy(trainmsa); 
  esl_msa_Destroy(test_msa); 
  while (esl_stack_PPop(teststack, &sq) == eslOK) esl_sq_Destroy(sq);
  esl_stack_Destroy(teststack);
  *ret_trainmsa  = NULL;
  *ret_teststack = NULL;
  return status;
}


static int
synthesize_testseqs(struct cfg_s *cfg, ESL_STACK *teststack)
{
  ESL_SQ *domain1, *domain2;
  int     L;

  
  

  

  

}


static int
fetch_random_segment(struct cfg_s *cfg, int L, ESL_DSQ **ret_dsq)
{
  ESL_DSQ *dsq = NULL;
  int      i;
  int      status;

  if (L> cfg->N) esl_fatal("can't fetch a segment of length %d; database is only %d\n", L, cfg->N);

  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;

  i = esl_rnd_Choose(cfg->r, cfg->N-L);
  if (fseek(cfg->dbfp, i, SEEK_SET)                        != 0)     esl_fatal("fseek() failed");
  if (fread(dsq+1, sizeof(ESL_DSQ), L, cfg->dbfp)          != L)     esl_fatal("fread() failed");
  if (esl_rnd_XShuffleDP(cfg->r, dsq, L, cfg->abc->K, dsq) != eslOK) esl_fatal("esl's DP shuffle failed");
  *ret_dsq = dsq;
  return eslOK;

 ERROR:
  if (dsq != NULL) free(dsq);
  *ret_dsq = NULL;
  return status;
}
  

  
