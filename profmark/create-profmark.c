/* Construct a training alignment/test sequences set from an MSA.
 *
     gcc -g -Wall -I ../src -I ../easel -L ../src -L ../easel -o create-profmark create-profmark.c -lhmmer -leasel -lm
     ./create-profmark <basename> <msa Stockholm file> <FASTA db>
   For example:
     ./create-profmark pmark /misc/data0/databases/Pfam/Pfam-A.seed /misc/data0/databases/uniprot-7.0/uniprot_sprot.fasta
 *
 * This generates five output files:
 *   <basename>.tbl  - table summarizing the benchmark
 *   <basename>.msa  - MSA queries, stockholm format
 *   <basename>.fa   - sequence targets, fasta format
 *   <basename>.pos  - table summarizing positive test set
 *   <basename>.neg  - table summarizing negative test set                
 * 
 * SRE, Thu Mar 27 11:05:38 2008 [Janelia]
 * SVN $Id$
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_distance.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_vectorops.h"


static char banner[] = "construct a benchmark profile training/test set";
static char usage[]  = "[options] <basename> <msafile> <seqdb>\n";

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                                   docgroup */
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              0 },
  { "-1",       eslARG_REAL, "0.25", NULL,"0<x<=1.0",NULL,NULL,NULL,         "require all test seqs to have < x id to training",        0 },
  { "-2",       eslARG_REAL, "0.50", NULL,"0<x<=1.0",NULL,NULL,NULL,         "require all test seqs to have < x id to each other",      0 },
  { "-F",       eslARG_REAL, "0.70", NULL,"0<x<=1.0",NULL,NULL,NULL,         "filter out seqs <x*average length",                       0 },
  { "-N",       eslARG_INT,"200000", NULL, NULL, NULL, NULL, NULL,           "number of negative test seqs",                            0 },
  { "--minDPL", eslARG_INT,   "100", NULL, NULL, NULL, NULL, NULL,           "minimum segment length for DP shuffling",                 0 },
  { "--amino",  eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--dna,--rna",    "<msafile> contains protein alignments",                   0 },
  { "--dna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--rna",  "<msafile> contains DNA alignments",                       0 },
  { "--rna",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--amino,--dna",  "<msafile> contains RNA alignments",                       0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

struct testseq_s {
  int  L;			/* total length            */
  int  L1;			/* length of first spacer  */
  int  d1n;			/* length of first domain  */
  int  L2;			/* length of second spacer */
  int  d2n;			/* length of second domain */
  int  L3;			/* length of third spacer  */
};

struct cfg_s {
  ESL_ALPHABET   *abc;          /* biological alphabet             */
  ESL_RANDOMNESS *r;            /* random number generator         */
  double          fragfrac;	/* seqs less than x*avg length are removed from alignment  */
  double          idthresh1;	/* fractional identity threshold for train/test split      */
  double          idthresh2;	/* fractional identity threshold for selecting test seqs   */

  FILE           *out_msafp;	/* output: training MSAs  */
  FILE           *out_seqfp;	/* output: test sequences */
  FILE           *possummfp;    /* output: summary table of the positive test set */
  FILE           *negsummfp;	/* output: summary table of the negative test set */
  FILE           *tblfp;	/* output: summary table of the training set alignments */

  int            *db_lens;	/* array of sequence lengths from db  [0..db_nseq]         */
  int             db_maxL;	/* maximum seq length in db_lens                           */
  ESL_DSQ        *db_dsq;	/* concatenated digitized sequences from db [1..db_nres]   */
  int             db_nres;	/* # of residues in the db                                 */
  int             db_nseq;	/* # of sequences in the db                                */

  struct testseq_s *test_lens;	/* array of length info about positive test seqs */
  int             ntest;	/* number of positive test seqs                  */
};



static int process_dbfile      (struct cfg_s *cfg, char *dbfile, int dbfmt);
static int remove_fragments    (struct cfg_s *cfg, ESL_MSA *msa, ESL_MSA **ret_filteredmsa, int *ret_nfrags);
static int separate_sets       (struct cfg_s *cfg, ESL_MSA *msa, ESL_MSA **ret_trainmsa, ESL_STACK **ret_teststack);
static int synthesize_positives(ESL_GETOPTS *go, struct cfg_s *cfg, char *testname, ESL_STACK *teststack, int *ret_ntest);
static int synthesize_negatives(ESL_GETOPTS *go, struct cfg_s *cfg, int nneg);
static int set_random_segment  (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_DSQ *dsq, int L);

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
  esl_opt_DisplayHelp(stdout, go, 0, 2, 80);
  exit(0);
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = NULL;	/* command line configuration      */
  struct cfg_s  cfg;     	/* application configuration       */
  char         *basename= NULL;	/* base of the output file names   */
  char         *alifile = NULL;	/* alignment file name             */
  char         *dbfile  = NULL;	/* name of seq db file             */
  char          outfile[256];	/* name of an output file          */
  int           alifmt;		/* format code for alifile         */
  int           dbfmt;		/* format code for dbfile          */
  ESL_MSAFILE  *afp     = NULL;	/* open alignment file             */
  ESL_MSA      *origmsa = NULL;	/* one multiple sequence alignment */
  ESL_MSA      *msa     = NULL;	/* MSA after frags are removed     */
  ESL_MSA      *trainmsa= NULL;	/* training set, aligned           */
  ESL_STACK    *teststack=NULL; /* test set: stack of ESL_SQ ptrs  */
  int           status;		/* easel return code               */
  int           nfrags;		/* # of fragments removed          */
  int           ntestdom;       /* # of test domains               */
  int           ntest;		/* # of test sequences created     */
  int           nali;		/* number of alignments read       */
  double        avgid;
  
  
  /* Parse command line */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h"))                    cmdline_help(argv[0], go);
  if (esl_opt_ArgNumber(go)                  != 3)     cmdline_failure(argv[0], "Incorrect number of command line arguments\n");
  basename = esl_opt_GetArg(go, 1); 
  alifile  = esl_opt_GetArg(go, 2);
  dbfile   = esl_opt_GetArg(go, 3);
  alifmt   = eslMSAFILE_STOCKHOLM;
  dbfmt    = eslSQFILE_FASTA;

  /* Set up the configuration structure shared amongst functions here */
  cfg.r         = esl_randomness_Create(42);
  cfg.abc       = NULL;		          /* until we open the MSA file, below */
  cfg.fragfrac  = esl_opt_GetReal(go, "-F");
  cfg.idthresh1 = esl_opt_GetReal(go, "-1");
  cfg.idthresh2 = esl_opt_GetReal(go, "-2");
  cfg.test_lens = NULL;
  cfg.ntest     = 0;

  /* Open the output files */ 
  if (snprintf(outfile, 256, "%s.msa", basename) >= 256)  esl_fatal("Failed to construct output MSA file name");
  if ((cfg.out_msafp = fopen(outfile, "w"))      == NULL) esl_fatal("Failed to open MSA output file %s\n", outfile);
  if (snprintf(outfile, 256, "%s.fa",  basename) >= 256)  esl_fatal("Failed to construct output FASTA file name");
  if ((cfg.out_seqfp = fopen(outfile, "w"))      == NULL) esl_fatal("Failed to open FASTA output file %s\n", outfile);
  if (snprintf(outfile, 256, "%s.pos", basename) >= 256)  esl_fatal("Failed to construct pos test set summary file name");
  if ((cfg.possummfp = fopen(outfile, "w"))      == NULL) esl_fatal("Failed to open pos test set summary file %s\n", outfile);
  if (snprintf(outfile, 256, "%s.neg", basename) >= 256)  esl_fatal("Failed to construct neg test set summary file name");
  if ((cfg.negsummfp = fopen(outfile, "w"))      == NULL) esl_fatal("Failed to open neg test set summary file %s\n", outfile);
  if (snprintf(outfile, 256, "%s.tbl", basename) >= 256)  esl_fatal("Failed to construct benchmark table file name");
  if ((cfg.tblfp     = fopen(outfile, "w"))      == NULL) esl_fatal("Failed to open benchmark table file %s\n", outfile);

  /* Open the MSA file; determine alphabet */
  status = esl_msafile_Open(alifile, alifmt, NULL, &afp);
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
    cfg.abc = esl_alphabet_Create(type);
  }
  esl_msafile_SetDigital(afp, cfg.abc);

  /* Suck the whole dbfile into memory; make sure it's in the same alphabet */
  process_dbfile(&cfg, dbfile, dbfmt);

  /* Read and process MSAs one at a time  */
  while ((status = esl_msa_Read(afp, &origmsa)) == eslOK)
    {
      remove_fragments(&cfg, origmsa, &msa, &nfrags);
      separate_sets   (&cfg, msa, &trainmsa, &teststack);
      ntestdom = esl_stack_ObjectCount(teststack);

      if (ntestdom >= 2) 
	{
	  esl_stack_Shuffle(cfg.r, teststack);
	  synthesize_positives(go, &cfg, msa->name, teststack, &ntest);

	  esl_msa_MinimGaps(trainmsa, NULL);
	  esl_msa_Write(cfg.out_msafp, trainmsa, eslMSAFILE_STOCKHOLM);

	  esl_dst_XAverageId(cfg.abc, trainmsa->ax, trainmsa->nseq, 10000, &avgid);
	  fprintf(cfg.tblfp, "%-20s  %3.0f%% %6d %6d %6d %6d %6d %6d\n", msa->name, 100.*avgid, trainmsa->alen, msa->nseq, nfrags, trainmsa->nseq, ntestdom, ntest);
	  nali++;
	}

      esl_msa_Destroy(trainmsa);
      esl_msa_Destroy(origmsa);
      esl_msa_Destroy(msa);
    }
  if      (status == eslEFORMAT)  esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", 
					    afp->linenumber, afp->fname, afp->errbuf, afp->buf);	
  else if (status != eslEOF)      esl_fatal("Alignment file read failed with error code %d\n", status);
  else if (nali   == 0)           esl_fatal("No alignments found in file %s\n", alifile);

  if (nali > 0)
    synthesize_negatives(go, &cfg, esl_opt_GetInteger(go, "-N"));

  free(cfg.db_dsq);
  free(cfg.db_lens);
  fclose(cfg.out_msafp);
  fclose(cfg.out_seqfp);
  fclose(cfg.possummfp);
  fclose(cfg.negsummfp);
  fclose(cfg.tblfp);
  esl_randomness_Destroy(cfg.r);
  esl_alphabet_Destroy(cfg.abc);
  esl_msafile_Close(afp);
  esl_getopts_Destroy(go);
  return 0;
}
      
  
/* Suck the FASTA sequence database into memory;
 * upon return, cfg->db_dsq[1..cfg->db_nres] is the concatenated digitized database;
 * cfg->db_lens is an array [0..cfg->db_nseq-1] of the sequence lengths.
 */
static int
process_dbfile(struct cfg_s *cfg, char *dbfile, int dbfmt)
{
  ESL_SQFILE *dbfp  = NULL;
  ESL_SQ     *sq    = NULL;
  void       *p;
  int         i;
  int         status;
  

  /* Open the sequence file in digital mode */
  status = esl_sqfile_Open(dbfile, dbfmt, NULL, &dbfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file %s", dbfile);
  else if (status == eslEFORMAT)   esl_fatal("Format of seqfile %s unrecognized.", dbfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  sq  = esl_sq_CreateDigital(cfg->abc);

  cfg->db_dsq  = NULL;
  cfg->db_lens = NULL;
  cfg->db_nres = 0;
  cfg->db_nseq = 0;
  cfg->db_maxL = 0;

  /* Read each sequence; output digital to dbfp. */
  while ((status = esl_sqio_Read(dbfp, sq)) == eslOK)
    {
      ESL_RALLOC(cfg->db_dsq,  p, sizeof(ESL_DSQ) * (sq->n + cfg->db_nres)+2);
      ESL_RALLOC(cfg->db_lens, p, sizeof(int) * (cfg->db_nseq+1));

      memcpy(cfg->db_dsq+cfg->db_nres+1, sq->dsq+1, sizeof(ESL_DSQ) * sq->n);
      cfg->db_lens[cfg->db_nseq] = sq->n;
      
      if (sq->n > cfg->db_maxL) cfg->db_maxL = sq->n;

      cfg->db_nres += sq->n;
      cfg->db_nseq++;

      esl_sq_Reuse(sq);
    }
  if (status != eslEOF) esl_fatal("Something went wrong with reading the seq db");
  cfg->db_dsq[0] = cfg->db_dsq[cfg->db_nres+1] = eslDSQ_SENTINEL;

  /* Let's just check, shall we? */
  for (i = 1; i <= cfg->db_nres; i++)
    if (cfg->db_dsq[i] >= cfg->abc->Kp) esl_fatal("whoops, something's wrong in db string");

  esl_sqfile_Close(dbfp);
  esl_sq_Destroy(sq);
  return eslOK;

 ERROR:
  esl_fatal("Failure in process_dbfile");
  return status; /*silence compiler warnings*/
}


/* Step 1. Label all sequence fragments < fragfrac of average raw length */
static int
remove_fragments(struct cfg_s *cfg, ESL_MSA *msa, ESL_MSA **ret_filteredmsa, int *ret_nfrags)
{
  int     *useme    = NULL;
  double   len      = 0.0;
  int      i;
  int      nfrags;
  int      status;

  for (i = 0; i < msa->nseq; i++) 
    len += esl_abc_dsqrlen(msa->abc, msa->ax[i]);
  len *= cfg->fragfrac / (double) msa->nseq;

  ESL_ALLOC(useme, sizeof(int) * msa->nseq);
  for (nfrags = 0, i = 0; i < msa->nseq; i++) 
    useme[i] = (esl_abc_dsqrlen(msa->abc, msa->ax[i]) < len) ? 0 : 1;

  if ((status = esl_msa_SequenceSubset(msa, useme, ret_filteredmsa)) != eslOK) goto ERROR;
  *ret_nfrags = msa->nseq - esl_vec_ISum(useme, msa->nseq);

  free(useme);
  return eslOK;

 ERROR:
  if (useme != NULL) free(useme); 
  *ret_filteredmsa = NULL;
  return status;
}

/* Step 2. Extract the training set and test set.
 */
static int
separate_sets(struct cfg_s *cfg, ESL_MSA *msa, ESL_MSA **ret_trainmsa, ESL_STACK **ret_teststack)
{      
  ESL_MSA   *trainmsa  = NULL;
  ESL_MSA   *test_msa  = NULL;
  ESL_STACK *teststack = NULL;
  ESL_SQ    *sq        = NULL;
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

  for (i = 0; i < msa->nseq; i++) useme[i] = (assignment[i] == ctrain) ? 1 : 0;
  if ((status = esl_msa_SequenceSubset(msa, useme, &trainmsa)) != eslOK) goto ERROR;

  /* If all the seqs went into the training msa, none are left for testing; we're done here */
  if (trainmsa->nseq == msa->nseq) {
    free(useme);
    free(assignment);
    free(nin);
    *ret_trainmsa  = trainmsa;
    *ret_teststack = teststack;
    return eslOK;
  }

  /* Put all the other sequences into an MSA of their own; from these, we'll
   * choose test sequences.
   */
  for (i = 0; i < msa->nseq; i++) useme[i] = (assignment[i] != ctrain) ? 1 : 0;
  if ((status = esl_msa_SequenceSubset(msa, useme, &test_msa))                             != eslOK) goto ERROR;

  /* Cluster those test sequences. */
  free(nin);         nin        = NULL;
  free(assignment);  assignment = NULL;
  if ((status = esl_msacluster_SingleLinkage(test_msa, cfg->idthresh2, &assignment, &nin, &nc)) != eslOK) goto ERROR;
  for (c = 0; c < nc; c++)
    {
      nskip = esl_rnd_Roll(cfg->r, nin[c]); /* pick a random seq in this cluster to be the test. */
      for (i=0; i < test_msa->nseq; i++)
	if (assignment[i] == c) {
	  if (nskip == 0) {
	    esl_sq_FetchFromMSA(test_msa, i, &sq);
	    esl_stack_PPush(teststack, (void *) sq);
	    break;
	  } else nskip--;
	}
    }

  esl_msa_Destroy(test_msa);
  free(useme);
  free(nin);
  free(assignment);

  *ret_trainmsa  = trainmsa;
  *ret_teststack = teststack;
  return eslOK;

 ERROR:
  if (useme      != NULL) free(useme);
  if (assignment != NULL) free(assignment);
  if (nin        != NULL) free(nin);
  esl_msa_Destroy(trainmsa); 
  esl_msa_Destroy(test_msa); 
  while (esl_stack_PPop(teststack, (void **) &sq) == eslOK) esl_sq_Destroy(sq);
  esl_stack_Destroy(teststack);
  *ret_trainmsa  = NULL;
  *ret_teststack = NULL;
  return status;
}


/* Each test sequence will contain two domains.
 */
static int
synthesize_positives(ESL_GETOPTS *go, struct cfg_s *cfg, char *testname, ESL_STACK *teststack, int *ret_ntest)
{
  ESL_SQ *domain1, *domain2;
  ESL_SQ *sq;
  void   *p;
  int     L;			/* total length of synthetic test seq */
  int     d1n, d2n;		/* lengths of two domains             */
  int     L1,L2,L3;		/* lengths of three random regions    */
  int     i,j;
  int     ntest = 0;
  int     status;
  
  while (esl_stack_ObjectCount(teststack) >= 2)
    {
      ESL_RALLOC(cfg->test_lens, p, (cfg->ntest+1) * sizeof(struct testseq_s));

      /* Pop our two test domains off the stack */
      esl_stack_PPop(teststack, &p);  domain1 = p;
      esl_stack_PPop(teststack, &p);  domain2 = p;
      d1n = domain1->n;
      d2n = domain2->n;

      /* Select a random total sequence length */
      if (d1n+d2n > cfg->db_maxL) esl_fatal("can't construct test seq; no db seq >= %d residues\n", d1n+d2n);
      do {
	L = cfg->db_lens[esl_rnd_Roll(cfg->r, cfg->db_nseq)];
      } while (L <= d1n+d2n);

      /* Select random lengths of three flanking domains;
       * Imagine picking two "insert after" points i,j in sequence 1..L', for
       * L' = L-d1n-d2n (the total length of nonhomologous test seq)
       */
      do {
	i = esl_rnd_Roll(cfg->r, L - d1n - d2n + 1 ); /* i = 0..L' */
	j = esl_rnd_Roll(cfg->r, L - d1n - d2n + 1 ); /* j = 0..L' */
      } while (i > j);

      /* now 1           .. i         = random region 1 (if i==0, there's none); 
       *     i+1         .. i+d1n     = domain 1
       *     i+d1n+1     .. j+d1n     = random region 2 (if i==j, there's none);
       *     j+d1n+1     .. j+d1n+d2n = domain 2
       *     j+d1n+d2n+1 .. L         = random region 3 (if j == L-d1n-d2n, there's none);
       */
      L1 = i;			
      L2 = j-i;
      L3 = L - d1n - d2n - j;
      
      sq = esl_sq_CreateDigital(cfg->abc);
      esl_sq_GrowTo(sq, L);
      esl_sq_SetName(sq, "%s/%d-%d/%d-%d", testname, i+1, i+d1n, j+d1n+1, j+d1n+d2n);
      esl_sq_SetDesc(sq, "domains: %s %s", domain1->name, domain2->name);

      sq->dsq[0] = sq->dsq[L+1] = eslDSQ_SENTINEL;
      set_random_segment(go, cfg, sq->dsq+1,           L1);
      set_random_segment(go, cfg, sq->dsq+i+d1n+1,     L2);
      set_random_segment(go, cfg, sq->dsq+j+d1n+d2n+1, L3);
      memcpy(sq->dsq+i+1,     domain1->dsq+1, sizeof(ESL_DSQ) * d1n);
      memcpy(sq->dsq+j+d1n+1, domain2->dsq+1, sizeof(ESL_DSQ) * d2n);
      sq->n = L;

      cfg->test_lens[cfg->ntest].L   = L;
      cfg->test_lens[cfg->ntest].L1  = L1;
      cfg->test_lens[cfg->ntest].d1n = d1n;
      cfg->test_lens[cfg->ntest].L2  = L2;
      cfg->test_lens[cfg->ntest].d2n = d2n;
      cfg->test_lens[cfg->ntest].L3  = L3;
      cfg->ntest++;
      ntest++;

      esl_sqio_Write(cfg->out_seqfp, sq, eslSQFILE_FASTA);

      fprintf(cfg->possummfp, "%-35s %5d %24s %5d %5d %24s %5d %5d %5d %5d %5d\n", 
	      sq->name, sq->n,
	      domain1->name, i+1, i+d1n,
	      domain2->name, j+d1n+1, j+d1n+d2n,
	      L1, L2, L3);
      
      esl_sq_Destroy(domain1);
      esl_sq_Destroy(domain2);
      esl_sq_Destroy(sq);
    }

  *ret_ntest = ntest;
  return eslOK;

 ERROR:
  esl_fatal("Failure in synthesize_positives");
  return status;
}

static int
synthesize_negatives(ESL_GETOPTS *go, struct cfg_s *cfg, int nneg)
{
  ESL_SQ *sq = esl_sq_CreateDigital(cfg->abc);
  int     a;
  int     i;
  int     L1,L2,L3,d1n,d2n;

  for (i = 0; i < nneg; i++)
    {
      /* Select a random test seq, to use its same segments */
      a = esl_rnd_Roll(cfg->r, cfg->ntest);

      L1  = cfg->test_lens[a].L1;
      L2  = cfg->test_lens[a].L2;
      L3  = cfg->test_lens[a].L3;
      d1n = cfg->test_lens[a].d1n;
      d2n = cfg->test_lens[a].d2n;

      esl_sq_GrowTo(sq, cfg->test_lens[a].L);

      esl_sq_SetName(sq, "decoy%d", i+1);
      esl_sq_SetDesc(sq, "L=%d in segments: %d/%d/%d/%d/%d", cfg->test_lens[a].L, L1, d1n, L2, d2n, L3);

      sq->dsq[0] = sq->dsq[cfg->test_lens[a].L+1] = eslDSQ_SENTINEL;
      set_random_segment(go, cfg, sq->dsq+1,               L1);
      set_random_segment(go, cfg, sq->dsq+1+L1,            d1n);
      set_random_segment(go, cfg, sq->dsq+1+L1+d1n,        L2);
      set_random_segment(go, cfg, sq->dsq+1+L1+d1n+L2,     d2n);
      set_random_segment(go, cfg, sq->dsq+1+L1+d1n+L2+d2n, L3);
      sq->n = cfg->test_lens[a].L;
  
      esl_sqio_Write(cfg->out_seqfp, sq, eslSQFILE_FASTA);

      fprintf(cfg->negsummfp, "%-15s %5d %5d %5d %5d %5d %5d\n", 
	      sq->name, sq->n,
	      L1, d1n, L2, d2n, L3);

      esl_sq_Reuse(sq);
    }

  esl_sq_Destroy(sq);
  return eslOK;
}

/* Fetch in a random sequence of length <L> from the
 * the pre-digitized concatenated
 * sequence database, select a random subseq, then 
 * shuffle it by diresidue composition.
 */
static int
set_random_segment(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_DSQ *dsq, int L)
{
  char    *tmp_dsq = NULL;
  int      minDPL  = esl_opt_GetInteger(go, "--minDPL");
  int      i;
  int      status;

  if (L==0) return eslOK;
  if (L> cfg->db_nres) esl_fatal("can't fetch a segment of length %d; database is only %d\n", L, cfg->db_nres);

  ESL_ALLOC(tmp_dsq, sizeof(ESL_DSQ) * (L+2));
  tmp_dsq[0] = tmp_dsq[L+1] = eslDSQ_SENTINEL;

  i = esl_rnd_Roll(cfg->r, cfg->db_nres-L);
  memcpy(tmp_dsq+1, cfg->db_dsq+i+1, sizeof(ESL_DSQ) * L);

  if (L < minDPL) status = esl_rsq_XShuffle  (cfg->r, tmp_dsq, L, tmp_dsq);
  else            status = esl_rsq_XShuffleDP(cfg->r, tmp_dsq, L, cfg->abc->Kp, tmp_dsq);
  if (status != eslOK) esl_fatal("esl's shuffling failed");

  memcpy(dsq, tmp_dsq+1, sizeof(ESL_DSQ) * L);
  free(tmp_dsq);
  return eslOK;

 ERROR:
  if (tmp_dsq != NULL) free(tmp_dsq);
  return status;
}
  

  
