/* main() for profile HMM construction from a multiple sequence alignment
 * 
 * SRE, Wed Jan  3 11:03:47 2007 [Janelia] [The Chemical Brothers]
 * SVN $Id$
 */

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msaweight.h"
#include "esl_msacluster.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#define CONOPTS "--fast,--hand"                                /* Exclusive options for model construction                    */
#define EFFOPTS "--eent,--eclust,--eset,--enone"               /* Exclusive options for effective sequence number calculation */
#define WGTOPTS "--wgsc,--wblosum,--wpb,--wnone,--wgiven"      /* Exclusive options for relative weighting                    */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",     1 },
/* Alternate model construction strategies */
  { "--fast",    eslARG_NONE,"default",NULL, NULL,    CONOPTS,    NULL,     NULL, "assign cols w/ >= symfrac residues as consensus",       2 },
  { "--hand",    eslARG_NONE,   FALSE, NULL, NULL,    CONOPTS,    NULL,     NULL, "manual construction (requires reference annotation)",   2 },
  { "--symfrac", eslARG_REAL,   "0.5", NULL, "0<=x<=1", NULL,   "--fast",   NULL, "sets sym fraction controlling --fast construction",     2 },
/* Alternate relative sequence weighting strategies */
  /* --wme not implemented in HMMER3 yet */
  { "--wgsc",    eslARG_NONE,"default",NULL, NULL,    WGTOPTS,    NULL,      NULL, "Gerstein/Sonnhammer/Chothia tree weights",         3},
  { "--wblosum", eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff simple filter weights",                   3},
  { "--wpb",     eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff position-based weights",                  3},
  { "--wnone",   eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "don't do any relative weighting; set all to 1",    3},
  { "--wgiven",  eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "use weights as given in MSA file",                 3},
  { "--pbswitch",eslARG_INT,  "1000",  NULL,"n>0",       NULL,    NULL,      NULL, "set failover to efficient PB wgts at > <n> seqs",  3},
  { "--wid",     eslARG_REAL, "0.62",  NULL,"0<=x<=1",   NULL,"--wblosum",   NULL, "for --wblosum: set identity cutoff",               3},
/* Alternate effective sequence weighting strategies */
  { "--eent",    eslARG_NONE,"default",NULL, NULL,    EFFOPTS,    NULL,      NULL, "adjust eff seq # to achieve relative entropy target", 4},
  { "--eclust",  eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "eff seq # is # of single linkage clusters",           4},
  { "--enone",   eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "no effective seq # weighting: just use nseq",         4},
  { "--eset",    eslARG_REAL,   NULL,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "set eff seq # for all models to <x>",                 4},
  { "--ere",     eslARG_REAL,   NULL,  NULL,"x>0",       NULL, "--eent",     NULL, "for --eent: set target relative entropy to <x>",      4},
  { "--eX",      eslARG_REAL,  "6.0",  NULL,"x>0",       NULL, "--eent",  "--ere", "for --eent: set minimum total rel ent param to <x>",  4},
  { "--eid",     eslARG_REAL, "0.62",  NULL,"0<=x<=1",   NULL,"--eclust",    NULL, "for --eclust: set fractional identity cutoff to <x>", 4},
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char banner[] = "profile HMM construction from a multiple sequence alignment";
static char usage[]  = "[-options] <hmmfile output> <alignment file input>";

static int  set_relative_weights(ESL_GETOPTS *go, ESL_MSA *msa);
static void build_model(ESL_GETOPTS *go, ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);
static void set_model_name(P7_HMM *hmm, char *setname, char *msa_name, char *alifile, int nali);
static int  set_effective_seqnumber(const ESL_GETOPTS *go, const ESL_MSA *msa, P7_HMM *hmm, const P7_BG *bg, const P7_DPRIOR *prior);
static double default_target_relent(const ESL_ALPHABET *abc, int M, double eX);

int
main(int argc, char **argv)
{
  int              status;	/* status of a function call               */
  ESL_ALPHABET    *abc;		/* sequence alphabet                       */
  ESL_GETOPTS     *go;		/* command line processing                 */
  char            *alifile;     /* seqfile to read alignment from          */
  int              fmt;	        /* format of the alifile                   */
  ESL_MSAFILE     *afp;         /* open alignment file                     */
  ESL_MSA         *msa;         /* a multiple sequence alignment           */
  int              nali;	/* count number of alignments/HMMs         */
  char            *hmmfile;     /* file to write HMM to                    */
  FILE            *hmmfp;       /* HMM output file handle                  */
  P7_HMM          *hmm;         /* constructed HMM; written to hmmfile     */
  P7_BG		  *bg;		/* null model                              */
  P7_DPRIOR       *pri;		/* mixture Dirichlet prior for the HMM     */
  char             errbuf[eslERRBUFSIZE];

  /*****************************************************************
   * Parse the command line
   *****************************************************************/

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK) 
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\n  where options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\n  alternative model construction strategies:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      puts("\n  alternative relative sequence weighting strategies:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\n  alternate effective sequence weighting strategies:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 2) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if ((hmmfile = esl_opt_GetArg(go, 1)) == NULL) p7_Fail("%s\nUsage: %s\n", go->errbuf, usage);
  if ((alifile = esl_opt_GetArg(go, 2)) == NULL) p7_Fail("%s\nUsage: %s\n", go->errbuf, usage);
  fmt     = eslMSAFILE_UNKNOWN;  /* autodetect alignment format by default. */

  /*****************************************************************
   * Set up the alphabet and prior
   *****************************************************************/
  
  abc = esl_alphabet_Create(eslAMINO);
  pri = p7_dprior_CreateAmino();
  bg  = p7_bg_Create(abc);

  /*****************************************************************
   * Open the alignment file (it might have >1 alignment)
   *****************************************************************/

  status = esl_msafile_OpenDigital(abc, alifile, fmt, NULL, &afp); /* NULL= no database dir from the environment */
  if      (status == eslENOTFOUND) p7_Fail("Alignment file %s doesn't exist or isn't readable.\n",     alifile);
  else if (status == eslEFORMAT)   p7_Fail("Couldn't determine format of alignment file %s.\n",        alifile);
  else if (status != eslOK)        p7_Fail("Alignment file open unexpectedly failed with error %d.\n", status);


  /*****************************************************************
   * Open the HMM output file.
   *****************************************************************/

  hmmfp = fopen(hmmfile, "w");
  if (hmmfp == NULL) p7_Fail("Failed to open HMM file %s for writing", hmmfile);

  /*****************************************************************
   * Read alignments one at a time, build HMMs, and save them.
   *****************************************************************/

  nali = 0;
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;

      /* Print some stuff about what we're about to do.
       */
      if (msa->name != NULL) printf("Alignment:           %s\n",  msa->name);
      else                   printf("Alignment:           #%d\n", nali);
      printf                       ("Number of sequences: %d\n",  msa->nseq);
      printf                       ("Number of columns:   %d\n",  msa->alen);
      puts("");
      fflush(stdout);

      set_relative_weights(go, msa);                       /* msa->wgt vector gets set. */
      build_model(go, msa, &hmm, NULL);                    /* Build <hmm> containing weighted observed counts.  */
      hmm->bg = bg;
      set_model_name(hmm, NULL, msa->name, alifile, nali); /* hmm->name gets set */
      set_effective_seqnumber(go, msa, hmm, bg, pri);      /* rescale the total counts in the model */
      p7_ParameterEstimation(hmm, pri);                    /* apply the prior; counts -> probability parameters */

      if (p7_hmm_Validate(hmm, 0.0001, errbuf) != eslOK)
	p7_Fail("HMM validation failed:\n%s", errbuf);

      status = p7_hmmfile_Write(hmmfp, hmm); 
      if (status != eslOK) p7_Fail("Failed to write model to disk.");

      /* Print some stuff about what we've done.
       */
      printf("Built a model of %d nodes.\n", hmm->M);
      printf("Mean match relative entropy:  %.2f bits\n", p7_MeanMatchRelativeEntropy(hmm, bg));
      printf("Mean match information:       %.2f bits\n", p7_MeanMatchInfo(hmm, bg));

      p7_hmm_Destroy(hmm);
      esl_msa_Destroy(msa);
    }
  if (status == eslEFORMAT) 
    p7_Fail("\
Alignment file parse error, line %d of file %s:\n\
%s\n\
Offending line is:\n\
%s\n", afp->linenumber, afp->fname, afp->errbuf, afp->buf);
  else if (status != eslEOF)
    p7_Fail("Alignment file read unexpectedly failed with code %d\n", status);
      
  fclose(hmmfp);
  esl_msafile_Close(afp);
  esl_getopts_Destroy(go);
  p7_dprior_Destroy(pri);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  return 0;
}


static int
set_relative_weights(ESL_GETOPTS *go, ESL_MSA *msa)
{
  printf("%-40s ... ", "Relative sequence weighting"); 
  fflush(stdout);

  if      (esl_opt_GetBoolean(go, "--wnone"))
    esl_vec_DSet(msa->wgt, msa->nseq, 1.);
  else if (esl_opt_GetBoolean(go, "--wgiven"))
    ;
  else if (msa->nseq >= esl_opt_GetInteger(go, "--pbswitch") ||
	   esl_opt_GetBoolean(go, "--wpb"))
    esl_msaweight_PB(msa);
  else if (esl_opt_GetBoolean(go, "--wgsc"))
    esl_msaweight_GSC(msa);
  else if (esl_opt_GetBoolean(go, "--wblosum"))
    esl_msaweight_BLOSUM(msa, esl_opt_GetReal(go, "--wid"));

  printf("done.\n");
  return eslOK;
}

static void
build_model(ESL_GETOPTS *go, ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***ret_tr)
{
  int status;

  printf("%-40s ... ", "Constructing model architecture"); 
  fflush(stdout);

  if      (esl_opt_GetBoolean(go, "--fast")) status = p7_Fastmodelmaker(msa, esl_opt_GetReal(go, "--symfrac"), ret_hmm, ret_tr);
  else if (esl_opt_GetBoolean(go, "--hand")) status = p7_Handmodelmaker(msa, ret_hmm, ret_tr);

  if (status == eslOK) { printf("done.\n"); return; }

  printf("[failed]\n");
  if (status == eslENORESULT) 
    {
      if (esl_opt_GetBoolean(go, "--fast"))
	printf("Alignment %s has no consensus columns w/ > %d%% residues - can't build a model.\n", 
	       msa->name != NULL ? msa->name : "",
	       (int) (100 * esl_opt_GetReal(go, "--symfrac")));
      else
	printf("Alignment %s has no annotated consensus columns - can't build a model.\n", 
	       msa->name != NULL ? msa->name : "");
    }
  else if (status == eslEFORMAT)
    printf("Alignment %s has no reference annotation line\n", 
	   msa->name != NULL ? msa->name : "");      

  else if (status == eslEMEM)
    printf("Memory allocation failure in model construction.\n");

  else 
    printf("unknown internal error.\n");

  exit(1);
}






/* set_model_name()
 * Give the model a name.
 * We deal with this differently depending on whether
 * we're in an alignment database or a single alignment.
 * 
 * If a single alignment, priority is:
 *      1. Use -n <name> if set.
 *      2. Use msa->name (avail in Stockholm or SELEX formats only)
 *      3. If all else fails, use alignment file name without
 *         filename extension (e.g. "globins.slx" gets named "globins"
 *         
 * If a multiple MSA database (e.g. Stockholm/Pfam), 
 * only msa->name is applied. -n is not allowed.
 * if msa->name is unavailable, or -n was used,
 * a fatal error is thrown.
 * 
 * Because we can't tell whether we've got more than one
 * alignment 'til we're on the second one, these fatal errors
 * only happen after the first HMM has already been built.
 * Oh well.
 */
static void
set_model_name(P7_HMM *hmm, char *setname, char *msa_name, char *alifile, int nali)
{
  char *name;

  printf("%-40s ... ", "Set model name");
  fflush(stdout);

  if (nali == 1)		/* first (only?) HMM in file:  */
    {
      if      (setname  != NULL) esl_strdup(setname,  -1, &name);
      else if (msa_name != NULL) esl_strdup(msa_name, -1, &name);
      else                       esl_FileTail(alifile, TRUE, &name); /* TRUE=nosuffix */
    }
  else
    {
      if      (setname  != NULL) p7_Fail("Oops. Wait. You can't use -n with an alignment database.\n");
      else if (msa_name != NULL) esl_strdup(msa_name, -1, &name);
      else    p7_Fail("Oops. Wait. I need name annotation on each alignment.\n");
    }

  p7_hmm_SetName(hmm, name);
  free(name);
  printf("done. [%s]\n", hmm->name);
}


/* set_effective_seqnumber()
 * Incept:    SRE, Fri May 11 08:14:57 2007 [Janelia]
 *
 * <hmm> comes in with weighted observed counts, and it goes out with
 * those observed counts rescaled to sum to the "effective sequence
 * number". 
 *
 * <msa> is needed because we may need to see the sequences in order 
 * to determine effective seq #. (for --eclust)
 *
 * <prior> is needed because we may need to parameterize test models
 * looking for the right relative entropy. (for --eent, the default)
 */
static int
set_effective_seqnumber(const ESL_GETOPTS *go, const ESL_MSA *msa, P7_HMM *hmm, const P7_BG *bg, const P7_DPRIOR *prior)
{
  int    status = eslOK;
  double neff;

  printf("%-40s ... ", "Set effective sequence number");
  fflush(stdout);

  if      (esl_opt_GetBoolean(go, "--enone") == TRUE) 
    {
      neff = msa->nseq;
      printf("done. [--enone: neff=nseq=%d]\n", msa->nseq);
    }
  else if (! esl_opt_IsDefault(go, "--eset"))
    {
      neff = esl_opt_GetReal(go, "--eset");
      printf("done. [--eset: set to neff = %.2f]\n", neff);
    }
  else if (esl_opt_GetBoolean(go, "--eclust") == TRUE)
    {
      int nclust;
      status = esl_msacluster_SingleLinkage(msa, esl_opt_GetReal(go, "--eid"), NULL, &nclust);
      neff = (double) nclust;

      printf("done. [--eclust SLC at %.1f%%; neff = %.2f clusters]\n", 100. * esl_opt_GetReal(go, "--eid"), neff);
    }
  
  else if (esl_opt_GetBoolean(go, "--eent") == TRUE)
    {
      double etarget; 

      if (esl_opt_IsDefault(go, "--ere")) etarget = default_target_relent(hmm->abc, hmm->M, esl_opt_GetReal(go, "--eX"));
      else                                etarget = esl_opt_GetReal(go, "--ere");

      status = p7_EntropyWeight(hmm, bg, prior, etarget, &neff);
    
      printf("done. [etarget %.2f bits; neff %.2f]\n", etarget, neff);
    }
    
  if (status == eslOK) {
    hmm->eff_nseq = neff;
    p7_hmm_Scale(hmm, neff / (double) hmm->nseq);
    return eslOK;
  }

  printf("[failed]\n");
  if (status == eslEINVAL && esl_opt_GetBoolean(go, "--eclust"))
    printf("Alignment %s seems to be corrupt;\nat least one pairwise distance calculation failed.\n",
	   msa->name != NULL ? msa->name : "");      
  else if (status == eslEMEM)
    printf("Memory allocation failure.\n");
  else 
    printf("unknown internal error.\n");
  exit(1);
}


/* default_amino_target_relent()
 * Incept:    SRE, Fri May 25 15:14:16 2007 [Janelia]
 *
 * Purpose:   Implements a length-dependent calculation of the target rel entropy
 *            per position, attempting to ensure that the information content of
 *            the model is high enough to find local alignments; but don't set it
 *            below a hard alphabet-dependent limit (p7_ETARGET_AMINO, etc.). See J1/67 for
 *            notes.
 *            
 * Args:      M  - model length in nodes
 *            eX - X parameter: minimum total rel entropy target
 *
 * Xref:      J1/67.
 */
static double
default_target_relent(const ESL_ALPHABET *abc, int M, double eX)
{
  double etarget;

  etarget = 6.* (eX + log((double) ((M * (M+1)) / 2)) / log(2.))    / (double)(2*M + 4);

  switch (abc->type) {
  case eslAMINO:  if (etarget < p7_ETARGET_AMINO)  etarget = p7_ETARGET_AMINO; break;
  case eslDNA:    if (etarget < p7_ETARGET_DNA)    etarget = p7_ETARGET_DNA;   break;
  case eslRNA:    if (etarget < p7_ETARGET_DNA)    etarget = p7_ETARGET_DNA;   break;
  default:        if (etarget < p7_ETARGET_OTHER)  etarget = p7_ETARGET_OTHER; break;
  }
  return etarget;
}
