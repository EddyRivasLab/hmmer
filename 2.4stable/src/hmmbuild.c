/* hmmbuild.c
 * main() for HMM construction from an alignment.
 *
 * SRE, Mon Nov 18 12:41:29 1996
 * SVN $Id$
 */

#include "config.h"		/* compile-time configuration constants */
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "lsjfuncs.h"		/* Steve Johnson's additions            */

static char banner[] = "hmmbuild - build a hidden Markov model from an alignment";

static char usage[]  = "\
Usage: hmmbuild [-options] <hmmfile output> <alignment file>\n\
  Available options are:\n\
   -h     : help; print brief help on version and usage\n\
   -n <s> : name; name this (first) HMM <s>\n\
   -o <f> : re-save annotated alignment to <f>\n\
   -A     : append; append this HMM to <hmmfile>\n\
   -F     : force; allow overwriting of <hmmfile>\n\
\n\
  Alternative alignment types: (default: multihit full domain alignment; hmmls)\n\
   -f     : multihit local (hmmfs style)\n\
   -g     : single global alignment (hmms style, Needleman/Wunsch)\n\
   -s     : single local alignment (hmmsw style, Smith/Waterman)\n\
";

static char experts[] = "\
  Forcing an alphabet: (normally autodetected)\n\
   --amino   : override autodetection, assert that seqs are protein\n\
   --nucleic : override autodetection, assert that seqs are DNA/RNA\n\
\n\
  Alternative model construction strategies: (default: --fast)\n\
   --fast        : assign cols w/ >= symfrac residues to be consensus {default}\n\
   --hand        : manual construction (requires annotated alignment)\n\
   --symfrac <x> : for --fast: min fraction of syms in MAT col {0.50} [0<=x<=1]\n\
\n\
  Customization of null model and priors:\n\
   --null  <f>   : read null (random sequence) model from file <f>\n\
   --prior <f>   : read Dirichlet prior parameters from file <f>\n\
   --pam   <f>   : heuristic PAM-based prior, using BLAST PAM matrix in <f>\n\
   --pamwgt <x>  : for --pam: set weight on PAM-based prior to <x> {20.}[>=0]\n\
\n\
  Alternative effective sequence number strategies:\n\
   --effclust    : eff seq # is # of clusters by single linkage [default]\n\
   --effent      : adjust eff seq # to achieve entropy loss target\n\
   --effnone     : effective sequence number is just # of seqs\n\
   --effset <x>  : set effective sequence number to <x>\n\
   --eloss <x>   : for --effent: set target loss [defaults: fs=0.59; ls=1.30]\n\
   --eidlevel <x>: for --effclust: set identity cutoff to <x> {0.62}\n\
\n\
  Relative sequence weighting strategies: (default: GSC weights)\n\
   --wblosum     : Henikoff simple filter weights (see --idlevel)\n\
   --wgsc        : Gerstein/Sonnhammer/Chothia tree weights (default)\n\
   --wme         : maximum entropy (ME)\n\
   --wpb         : Henikoff position-based weights\n\
   --wvoronoi    : Sibbald/Argos Voronoi weights\n\
   --wnone       : don't do any relative weighting\n\
   --pbswitch <n>: set switch from GSC to position-based wgts at > n seqs\n\
   --widlevel <x>: for --wblosum: set identity cutoff {0.62}\n\
\n\
  Evolving information content of match emission probabilities (experimental)\n\
   --evolve       : evolve information content of models\n\
   --evolveic <x> : evolve match emission probabilities to average info content [>=0]\n\
   --matrix <f>   : use matrix in <f> rather than default WAG\n\
\n\
  Other expert options:\n\
   --binary      : save the model in binary format, not ASCII text\n\
   --cfile <f>   : save count vectors to <f>\n\
   --informat <s>: input alignment is in format <s>, not Stockholm\n\
   --verbose     : print boring information\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-f", TRUE, sqdARG_NONE },
  { "-g", TRUE, sqdARG_NONE },
  { "-h", TRUE, sqdARG_NONE },
  { "-n", TRUE, sqdARG_STRING},
  { "-o", TRUE, sqdARG_STRING},
  { "-s", TRUE, sqdARG_NONE },
  { "-A", TRUE, sqdARG_NONE },
  { "-F", TRUE, sqdARG_NONE },
  { "--amino",   FALSE, sqdARG_NONE  },
  { "--binary",  FALSE, sqdARG_NONE  },
  { "--cfile",   FALSE, sqdARG_STRING},
  { "--effnone", FALSE, sqdARG_NONE },
  { "--effclust",FALSE, sqdARG_NONE },
  { "--effent",  FALSE, sqdARG_NONE },
  { "--effset",  FALSE, sqdARG_FLOAT },
  { "--eidlevel",FALSE, sqdARG_FLOAT },
  { "--eloss",   FALSE, sqdARG_FLOAT },
  { "--evolve",  FALSE, sqdARG_NONE  },
  { "--evolveic", FALSE, sqdARG_FLOAT },
  { "--fast",    FALSE, sqdARG_NONE },
  { "--hand",    FALSE, sqdARG_NONE},
  { "--informat",FALSE, sqdARG_STRING },
  { "--matrix",  FALSE, sqdARG_STRING },
  { "--nucleic", FALSE, sqdARG_NONE },
  { "--null",    FALSE, sqdARG_STRING },
  { "--pam",     FALSE, sqdARG_STRING },
  { "--pamwgt",  FALSE, sqdARG_FLOAT },
  { "--pbswitch",FALSE, sqdARG_INT },
  { "--prior",   FALSE, sqdARG_STRING },
  { "--symfrac", FALSE, sqdARG_FLOAT },
  { "--verbose", FALSE, sqdARG_NONE  },
  { "--wgsc",    FALSE, sqdARG_NONE },
  { "--wblosum", FALSE, sqdARG_NONE },
  { "--widlevel",FALSE, sqdARG_FLOAT },
  { "--wme",     FALSE, sqdARG_NONE },
  { "--wnone",   FALSE, sqdARG_NONE },
  { "--wpb",     FALSE, sqdARG_NONE },
  { "--wvoronoi",FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))




enum wgtconfig_e { WGT_NONE, WGT_GSC, WGT_BLOSUM, WGT_PB, WGT_VORONOI, WGT_ME };
enum cconfig_e   { P7_FAST_CONSTRUCTION,  P7_HAND_CONSTRUCTION };
enum effconfig_e { EFF_NONE, EFF_USERSET, EFF_NCLUSTERS, EFF_ENTROPY };

/* Structure: p7config_s
 * Configuration options for a model under construction. We keep all
 * this stuff in a static structure in hmmbuild.c, to keep
 * it bundled out of the way, in the hope of making the main()
 * a little cleaner and more understandable.
 */
struct p7config_s {
  char *setname;		/* name of the HMM */

  /* The null model. */
  char *rndfile;                /* file to read random model from, or NULL for default */
  float randomseq[MAXABET];
  float p1;

  /* The priors. */
  char *prifile;		/* Dirichlet prior file to read          */
  char *pamfile;		/* PAM matrix file for heuristic prior   */
  float pamwgt;			/* weight on PAM for heuristic prior     */

  /* Fragment identification in alignments */
  float fragthresh;		/* rlen_i < fragthresh*mean rlen -> fragment  */

  /* Relative sequence weighting */
  enum wgtconfig_e w_strategy;	/* which relative weight strategy we're using */
  int   pbswitch;		/* if nseq >= this, failover to PB weights    */
  float widlevel;		/* --wblosum: frac id filtering level [0.62]  */

  /* Model construction */
  enum cconfig_e c_strategy;	/* which construction algorithm we're using    */
  float symfrac;		/* --fast: min frac of syms in MATCH col, 0..1 */
  int   symfrac_set;		/* T/F, was symfrac set at command line?       */
  
  /* Effective sequence number */
  enum effconfig_e eff_strategy;/* which eff sequence # calculation we're using */
  float eff_nseq;		/* effective sequence number                    */
  float eidlevel;		/* --effclust: frac id filter level [0.62]      */
  float eloss;		        /* --effent: target entropy loss                */
  int   eloss_set;		/* --effent: TRUE if eloss set on commandline   */
  float etarget;		/* --effent: target entropy (background - eloss)*/

  /* Phylogenetic extrapolation */
  int	evolve_ic;              /* TRUE to evolve to specified info content */
  float info;        		/* specified ave model info content      */
  char *matrixfile;	   	/* open file containing rate matrices    */
  
  /* Choice of algorithm mode */
  enum p7_algmode mode;

  /* Optional output files for hmmbuild
   */
  char *align_ofile;            /* name of output alignment file, or NULL */
  FILE *alignfp;
  char *cfile;			/* output file for count vectors, or NULL */
  FILE *cfp;

  /* Control over normal hmmbuild outputs
   */
  int   overwrite_protect;	/* TRUE to prevent overwriting HMM file  */
  int   verbose;		/* TRUE to show a lot of output          */
  int   do_append;		/* TRUE to append to hmmfile             */
  int   do_binary;		/* TRUE to write in binary format        */
};




static void  default_config(struct p7config_s *cfg);
static void  process_cmdline(int argc, char **argv, struct p7config_s *cfg, 
			     char **ret_hmmfile, char **ret_alifile, int *ret_format);
static void  verify_options(struct p7config_s *cfg, char *hmmfile);
static FILE *open_hmmfile(struct p7config_s *cfg, char *hmmfile);
static void  open_optional_outputfiles(struct p7config_s *cfg);
static void  print_config_header(struct p7config_s *cfg, char *alifile, int aliformat, char *hmmfile);
static void  set_alphabet(MSA *msa);
static void  set_relative_weights(MSA *msa, struct p7config_s *cfg);
static float set_entropy_target(int eloss_set, float eloss, int mode, float *randomseq);
static int   tag_candidate_seq_fragments(MSA *msa, float thresh, 
					 int *rlen, char *isfrag);
static void  construct_model(struct p7config_s *cfg, MSA *msa, unsigned char **dsq, char *isfrag, 
			     struct plan7_s **ret_hmm, struct p7trace_s ***ret_tr);
static void  set_effective_seqnumber(struct p7config_s *cfg, MSA *msa, struct plan7_s *hmm, struct p7prior_s *pri);
static void  save_countvectors(FILE *cfp, char *cfile, char *name, struct plan7_s *hmm);
static void  position_average_score(struct plan7_s *hmm, unsigned char **dsq, float *wgt,
				   int nseq, struct p7trace_s **tr, float *pernode,
				   float *ret_avg);
static float frag_trace_score(struct plan7_s *hmm, unsigned char *dsq, struct p7trace_s *tr,
			      float *pernode, float expected);
static void  maximum_entropy(struct plan7_s *hmm, unsigned char **dsq, MSA *msa,
			     float eff_nseq,
			     struct p7prior_s *prior, struct p7trace_s **tr);
static void  set_model_name(struct plan7_s *hmm, char *setname, char *msa_name, char *alifile, int nali);
static void  print_statistics(FILE *fp, struct plan7_s *hmm, unsigned char **dsq, int nseq,
			      struct p7trace_s **tr);
static void  save_hmmbuild_alignment(FILE *alignfp, MSA *msa, unsigned char **dsq, struct plan7_s *hmm,
				     struct p7trace_s **tr);
static void  print_all_scores(FILE *fp, struct plan7_s *hmm,
			      unsigned char **dsq, MSA *msa, struct p7trace_s **tr);





int
main(int argc, char **argv)
{
  struct p7config_s cfg;	/* holds a multitude of config options     */
  char            *hmmfile;     /* file to write HMM to                    */
  FILE            *hmmfp;       /* HMM output file handle                  */
  char            *alifile;     /* seqfile to read alignment from          */
  int              format;	/* format of the alifile                   */
  MSAFILE         *afp;         /* open alignment file                     */
  MSA             *msa;         /* a multiple sequence alignment           */
  unsigned char  **dsq;         /* digitized unaligned aseq's              */
  struct plan7_s  *hmm;         /* constructed HMM; written to hmmfile     */
  struct p7prior_s *pri;        /* Dirichlet priors to use                 */
  struct p7trace_s **tr;        /* fake tracebacks for aseq's              */
  int              idx;		/* counter for sequences                   */
  int              nali;	/* count number of alignments/HMMs         */
  int             *rlen;	/* raw (unaligned) seq lengths 0..nseq-1   */
  char            *isfrag;	/* TRUE/FALSE for candidate seq frags      */

  /***********************************************
   * Parse command line
   ***********************************************/

  format            = MSAFILE_UNKNOWN;  /* autodetect format by default.      */
  Alphabet_type     = hmmNOTSETYET;	/* initially unknown */
  default_config(&cfg);		        /* see buildconfig.h for this structure */
  process_cmdline(argc, argv, &cfg, &hmmfile, &alifile, &format);
  verify_options(&cfg, hmmfile);

  /***********************************************
   * Preliminaries: open our files for i/o; print header
   ***********************************************/

  if ((afp = MSAFileOpen(alifile, format, NULL)) == NULL)
    Die("Alignment file %s could not be opened for reading", alifile);
  hmmfp = open_hmmfile(&cfg, hmmfile);
  open_optional_outputfiles(&cfg);
  print_config_header(&cfg, alifile, afp->format, hmmfile);

  /***********************************************
   * Get alignment(s), build HMMs one at a time
   ***********************************************/

  nali = 0;
  while ((msa = MSAFileRead(afp)) != NULL)
    {
      /* Print some stuff about what we're about to do.
       */
      if (msa->name != NULL) printf("Alignment:           %s\n",  msa->name);
      else                   printf("Alignment:           #%d\n", nali+1);
      printf                       ("Number of sequences: %d\n",  msa->nseq);
      printf                       ("Number of columns:   %d\n",  msa->alen);
      puts("");
      fflush(stdout);

      /* Make alignment upper case, because some symbol counting
       * things are case-sensitive.
       */
      for (idx = 0; idx < msa->nseq; idx++)
	s2upper(msa->aseq[idx]);

      /* --- Post-alphabet initialization section ---
       * (This code must be delayed until after we've seen the
       * first alignment, because we have to see the alphabet type first.)
       */
      if (nali == 0)
	{
	  if (Alphabet_type == hmmNOTSETYET)  set_alphabet(msa);	

	  if (cfg.prifile == NULL)  pri = P7DefaultPrior();
	  else                      pri = P7ReadPrior(cfg.prifile);

	  if (cfg.pamfile != NULL)  PAMPrior(cfg.pamfile, pri, cfg.pamwgt);

	  if (cfg.rndfile == NULL)  P7DefaultNullModel(cfg.randomseq, &(cfg.p1));
	  else                      P7ReadNullModel(cfg.rndfile, cfg.randomseq, &(cfg.p1));

	  if (cfg.eff_strategy == EFF_ENTROPY) 
	    cfg.etarget = set_entropy_target(cfg.eloss_set, cfg.eloss, cfg.mode, cfg.randomseq);
	} /* -- this ends the post-alphabet initialization section -- */

      /* Prepare unaligned digitized sequences for internal use */
      DigitizeAlignment(msa, &dsq);

      /* Determine relative sequence weights (except ME, which comes later).
       * The weights are stored in msa->wgt.
       */
      set_relative_weights(msa, &cfg);

      /* Identify candidate seq frags, to inform the model
       * construction algorithms. Relative sequence weighting must precede this.
       * Upon return, rlen[i] is the unaligned length of seq i [0..nseq-1];
       * isfrag[i] is a 1/0 flag for whether we defined it as a frag or not.
       */
      rlen   = MallocOrDie(sizeof(int)  * msa->nseq);
      isfrag = MallocOrDie(sizeof(char) * msa->nseq);
      tag_candidate_seq_fragments(msa, cfg.fragthresh, rlen, isfrag);

      /* Build the model architecture, and collect counts in it.
       * Upon return, the core HMM contains counts.
       */
      construct_model(&cfg, msa, dsq, isfrag, &hmm, &tr);

      /* Determine "effective sequence number", and rescale counts. */
      set_effective_seqnumber(&cfg, msa, hmm, pri);

      /* Save the count vectors if asked. Used primarily for
       * making the data files for training priors.
       */
      if (cfg.cfile != NULL)
	save_countvectors(cfg.cfp, cfg.cfile, (msa->name != NULL ? msa->name : "-"), hmm);
      /* Plan7_DumpCounts(stdout, hmm); */

      /* Record the null model in the HMM;
       * add prior contributions in pseudocounts and renormalize.
       */
      printf("%-40s ... ", "Converting counts to probabilities"); fflush(stdout);
      Plan7SetNullModel(hmm, cfg.randomseq, cfg.p1);
      P7PriorifyHMM(hmm, pri);	/* this function also renormalizes. */
      printf("done.\n");

      /* if evolving information content: */
      if (cfg.evolve_ic)
      {
        /* if info wasn't defined by user, use default settings for gl or ll */
        if (cfg.info == 0.0) {
	  if (cfg.mode == P7_FS_MODE || cfg.mode == P7_SW_MODE) { cfg.info = 0.65; } /* ll default */
	  else { cfg.info = 1.4; } /* gl default */
	}
        AdjustAveInfoContent(hmm, cfg.info, cfg.matrixfile);  /* DJB */
      }

      /* Model configuration, temporary.  hmmbuild assumes that it's
       * given an alignment of single domains, and the alignment may
       * contain fragments. So, for the purpose of scoring the
       * sequences (or, optionally, MD/ME weighting), configure the
       * model into S/W mode. Later we'll configure the model
       * according to how the user wants to use it.
       */
      Plan7SWConfig(hmm, 0.5, 0.5);
      /* Plan7_DumpScores(stdout, hmm); */


      /* Maximum entropy is the one relative weighting strategy
       * that requires having a model architecture, and a model
       * scoring configuration.
       */
      if (cfg.w_strategy == WGT_ME) 
	maximum_entropy(hmm, dsq, msa, cfg.eff_nseq, pri, tr);

      /* Transfer other information from the alignment to
       * the HMM. This typically only works for Stockholm or SELEX format
       * alignments, so these things are conditional/optional.
       */
      if (msa->acc  != NULL) Plan7SetAccession(hmm,   msa->acc);
      if (msa->desc != NULL) Plan7SetDescription(hmm, msa->desc);

      if (msa->cutoff_is_set[MSA_CUTOFF_GA1] && msa->cutoff_is_set[MSA_CUTOFF_GA2])
	{ hmm->flags |= PLAN7_GA; hmm->ga1 = msa->cutoff[MSA_CUTOFF_GA1]; hmm->ga2 = msa->cutoff[MSA_CUTOFF_GA2]; }
      if (msa->cutoff_is_set[MSA_CUTOFF_TC1] && msa->cutoff_is_set[MSA_CUTOFF_TC2])
	{ hmm->flags |= PLAN7_TC; hmm->tc1 = msa->cutoff[MSA_CUTOFF_TC1]; hmm->tc2 = msa->cutoff[MSA_CUTOFF_TC2]; }
      if (msa->cutoff_is_set[MSA_CUTOFF_NC1] && msa->cutoff_is_set[MSA_CUTOFF_NC2])
	{ hmm->flags |= PLAN7_NC; hmm->nc1 = msa->cutoff[MSA_CUTOFF_NC1]; hmm->nc2 = msa->cutoff[MSA_CUTOFF_NC2]; }

      /* Record some other miscellaneous information in the HMM,
       * like its name, and how/when we built it.
       */
      set_model_name(hmm, cfg.setname, msa->name, alifile, nali);
      Plan7ComlogAppend(hmm, argc, argv);
      Plan7SetCtime(hmm);
      hmm->nseq = msa->nseq;
   
      /* Print information for the user
       */
      printf("\nConstructed a profile HMM (length %d)\n", hmm->M);
      print_statistics(stdout, hmm, dsq, msa->nseq, tr); 
      printf("\n");


      /* Configure the model for chosen algorithm.
       * For an ASCII text model format, this isn't actually necessary any more.
       */
      printf("%-40s ... ", "Finalizing model configuration");
      fflush(stdout);
      switch (cfg.mode) {
      case P7_S_MODE:   Plan7GlobalConfig(hmm);        break;
      case P7_SW_MODE:  Plan7SWConfig(hmm, 0.5, 0.5);  break;
      case P7_LS_MODE:  Plan7LSConfig(hmm);            break;
      case P7_FS_MODE:  Plan7FSConfig(hmm, 0.5, 0.5);  break;
      default:          Die("bogus configuration choice");
      }
      printf("done.\n");

      /* Save new HMM to disk: open a file for appending or writing.
       */
      printf("%-40s ... ", "Saving model to file");  fflush(stdout);
      if (cfg.do_binary) WriteBinHMM(hmmfp, hmm);
      else               WriteAscHMM(hmmfp, hmm);
      printf("done.\n");

      /* the annotated alignment may be resaved */
      if (cfg.alignfp != NULL) 
	save_hmmbuild_alignment(cfg.alignfp, msa, dsq, hmm, tr);

      /* Verbose output; show scores for each sequence
       */
      if (cfg.verbose)
	print_all_scores(stdout, hmm, dsq, msa, tr);

      /* Clean up before moving on to next alignment
       */
      for (idx = 0; idx < msa->nseq; idx++) P7FreeTrace(tr[idx]);
      free(tr);
      FreePlan7(hmm);
      Free2DArray((void **) dsq, msa->nseq); 
      free(rlen);
      free(isfrag);
      MSAFree(msa);
      fflush(hmmfp);
      if (cfg.cfp != NULL)     fflush(cfg.cfp);
      if (cfg.alignfp != NULL) fflush(cfg.alignfp);

      puts("//\n");
      nali++;
    }


  /* Clean up and exit
   */
  MSAFileClose(afp);
  fclose(hmmfp);
  if (cfg.cfp     != NULL) fclose(cfg.cfp);
  if (cfg.alignfp != NULL) fclose(cfg.alignfp);
  P7FreePrior(pri);
  SqdClean();
  return 0;
}



/* default_configuration()
 * 
 */
static void
default_config(struct p7config_s *cfg)
{
  cfg->setname      = NULL;
  cfg->rndfile      = NULL;
  FSet(cfg->randomseq, MAXABET, 0.); /* not set yet */
  cfg->p1           = 0.;	/* not set yet  */
  cfg->prifile      = NULL;
  cfg->pamfile      = NULL;
  cfg->pamwgt       = 20.;
  cfg->fragthresh   = 0.5;
  cfg->w_strategy   = WGT_GSC;
  cfg->pbswitch     = 1000;
  cfg->widlevel     = 0.62;
  cfg->c_strategy   = P7_FAST_CONSTRUCTION;
  cfg->symfrac      = 0.5;
  cfg->symfrac_set  = FALSE;
  cfg->eff_strategy = EFF_NCLUSTERS;
  cfg->eff_nseq     = 0.;	/* not set yet */
  cfg->eidlevel     = 0.62;
  cfg->eloss        = 0.0;	/* not set yet */
  cfg->eloss_set    = FALSE;
  cfg->etarget      = 0.0;	/* not set yet */
  cfg->evolve_ic    = FALSE;
  cfg->info         = 0.;	/* not set yet */
  cfg->matrixfile   = NULL;
  cfg->mode         = P7_LS_MODE;
  cfg->align_ofile  = NULL;
  cfg->cfile        = NULL;
  cfg->overwrite_protect = TRUE;
  cfg->verbose           = FALSE;
  cfg->do_append         = FALSE;
  cfg->do_binary         = FALSE;
}

static void
process_cmdline(int argc, char **argv, struct p7config_s *cfg, char **hmmfile, char **alifile, int *ret_format)
{
  char *optname;                /* name of option found by Getopt()      */
  char *optarg;                 /* argument found by Getopt()            */
  int   optind;                 /* index in argv[]                       */

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-f") == 0) cfg->mode              = P7_FS_MODE;
    else if (strcmp(optname, "-g") == 0) cfg->mode              = P7_S_MODE;
    else if (strcmp(optname, "-n") == 0) cfg->setname           = optarg;
    else if (strcmp(optname, "-o") == 0) cfg->align_ofile       = optarg;
    else if (strcmp(optname, "-s") == 0) cfg->mode              = P7_SW_MODE;
    else if (strcmp(optname, "-A") == 0) cfg->do_append         = TRUE;
    else if (strcmp(optname, "-F") == 0) cfg->overwrite_protect = FALSE;
    else if (strcmp(optname, "--amino")   == 0) SetAlphabet(hmmAMINO);
    else if (strcmp(optname, "--binary")  == 0) cfg->do_binary     = TRUE;
    else if (strcmp(optname, "--cfile")   == 0) cfg->cfile         = optarg;
    else if (strcmp(optname, "--effclust")== 0) cfg->eff_strategy  = EFF_NCLUSTERS;
    else if (strcmp(optname, "--effent")  == 0) cfg->eff_strategy  = EFF_ENTROPY;
    else if (strcmp(optname, "--effnone") == 0) cfg->eff_strategy  = EFF_NONE;
    else if (strcmp(optname, "--effset")  == 0) { cfg->eff_strategy= EFF_USERSET; cfg->eff_nseq = atof(optarg); }
    else if (strcmp(optname, "--eidlevel")== 0) cfg->eidlevel      = atof(optarg);
    else if (strcmp(optname, "--eloss")   == 0) { cfg->eloss       = atof(optarg); cfg->eloss_set  = TRUE; }
    else if (strcmp(optname, "--evolve")  == 0) cfg->evolve_ic     = TRUE;
    else if (strcmp(optname, "--evolveic")== 0) { cfg->evolve_ic   = TRUE; cfg->info          = atof(optarg); }
    else if (strcmp(optname, "--hand")    == 0) cfg->c_strategy    = P7_HAND_CONSTRUCTION;
    else if (strcmp(optname, "--matrix")  == 0) cfg->matrixfile    = optarg;
    else if (strcmp(optname, "--nucleic") == 0) SetAlphabet(hmmNUCLEIC);
    else if (strcmp(optname, "--null")    == 0) cfg->rndfile       = optarg;
    else if (strcmp(optname, "--pam")     == 0) cfg->pamfile       = optarg;
    else if (strcmp(optname, "--pamwgt")  == 0) cfg->pamwgt        = atof(optarg);
    else if (strcmp(optname, "--pbswitch")== 0) cfg->pbswitch      = atoi(optarg);
    else if (strcmp(optname, "--prior")   == 0) cfg->prifile       = optarg;
    else if (strcmp(optname, "--symfrac") == 0) { cfg->symfrac     = atof(optarg); cfg->symfrac_set = TRUE; }
    else if (strcmp(optname, "--verbose") == 0) cfg->verbose       = TRUE;
    else if (strcmp(optname, "--wgsc")    == 0) cfg->w_strategy    = WGT_GSC;
    else if (strcmp(optname, "--wblosum") == 0) cfg->w_strategy    = WGT_BLOSUM;
    else if (strcmp(optname, "--widlevel")== 0) cfg->widlevel      = atof(optarg);
    else if (strcmp(optname, "--wme")     == 0) cfg->w_strategy    = WGT_ME;
    else if (strcmp(optname, "--wpb")     == 0) cfg->w_strategy    = WGT_PB;
    else if (strcmp(optname, "--wnone")   == 0) cfg->w_strategy    = WGT_NONE;
    else if (strcmp(optname, "--wvoronoi")== 0) cfg->w_strategy    = WGT_VORONOI;
    else if (strcmp(optname, "--informat") == 0) {
      *ret_format = String2SeqfileFormat(optarg);
      if (*ret_format == MSAFILE_UNKNOWN)
	Die("unrecognized sequence file format \"%s\"", optarg);
      if (! IsAlignmentFormat(*ret_format))
	Die("%s is an unaligned format, can't read as an alignment", optarg);
    }
    else if (strcmp(optname, "-h") == 0) {
      HMMERBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }
  }
  if (argc - optind != 2)
    Die("Incorrect number of arguments.\n%s\n", usage);

  *hmmfile = argv[optind++];
  *alifile = argv[optind++];
}


/* verify_options()
 * 
 * Check that nothing silly has been done to the configuration.
 * More could be done here; eventually, we'll bring in Easel's
 * improved command line parser.
 */
static void
verify_options(struct p7config_s *cfg, char *hmmfile)
{
  if (cfg->symfrac < 0. || cfg->symfrac > 1.)
    Die("--symfrac must be a value from 0 to 1\n%s\n", usage);

  if (cfg->overwrite_protect && !cfg->do_append && FileExists(hmmfile))
    Die("HMM file %s already exists. Rename or delete it.", hmmfile);

  if (cfg->overwrite_protect && cfg->align_ofile != NULL && FileExists(cfg->align_ofile))
    Die("Alignment resave file %s exists. Rename or delete it.", cfg->align_ofile);

  if (cfg->symfrac_set && cfg->c_strategy  != P7_FAST_CONSTRUCTION)
    Die("using --symfrac only makes sense for default --fast construction strategy");
}

/* open_hmmfile()
 * 
 * Open the output hmmfile, return a writable FILE ptr.
 */
static FILE *
open_hmmfile(struct p7config_s *cfg, char *hmmfile)
{
  FILE *hmmfp;
  char  fpopts[3];   /* options to open with, e.g. "ab"  */

  if (cfg->do_append) strcpy(fpopts, "a");
  else                strcpy(fpopts, "w");
  if (cfg->do_binary) strcat(fpopts, "b");
  if ((hmmfp = fopen(hmmfile, fpopts)) == NULL)
    Die("Failed to open HMM file %s for %s\n", hmmfile,
	cfg->do_append ? "appending" : "writing");

  return hmmfp;
}

/* open_optional_outputfiles()
 * 
 * counts file, alignment resave file...
 */
static void
open_optional_outputfiles(struct p7config_s *cfg)
{
  				/* optional count vector save file */
  cfg->cfp = NULL;
  if (cfg->cfile != NULL)
    if ((cfg->cfp = fopen(cfg->cfile, "w")) == NULL)
      Die("Failed to open count vector file %s for writing\n", cfg->cfile);

				/* optional alignment resave file */
  cfg->alignfp = NULL;
  if (cfg->align_ofile != NULL)
    if ((cfg->alignfp = fopen(cfg->align_ofile, "w")) == NULL)
      Die("Failed to open alignment resave file %s for writing\n", cfg->align_ofile);
}


/* print_config_header()
 * 
 * Print out the configuration, at the start of an hmmbuild
 * output.
 */
static void
print_config_header(struct p7config_s *cfg, char *alifile, int aliformat, char *hmmfile)
{
  HMMERBanner(stdout, banner);
  printf("Alignment file:                    %s\n",
	 alifile);
  printf("File format:                       %s\n",
	 SeqfileFormat2String(aliformat));

  printf("Search algorithm configuration:    ");
  if      (cfg->mode == P7_S_MODE)   puts("Global alignment (hmms)");
  else if (cfg->mode == P7_SW_MODE)  puts("local Smith/Waterman (hmmsw)");
  else if (cfg->mode == P7_LS_MODE)  puts("glocal, multihit (hmmls)");
  else if (cfg->mode == P7_FS_MODE)  puts("local, multihit (hmmfs)");
  else Die("whoops. can't build in mode %d\n", cfg->mode);

  printf("Model construction strategy:       ");
  if (cfg->c_strategy == P7_HAND_CONSTRUCTION)    puts("Manual, from #=RF annotation");
  else if (cfg->c_strategy==P7_FAST_CONSTRUCTION) printf("Fast/ad hoc (symfrac %.2f)\n", cfg->symfrac);

  printf("Null model used:                   %s\n",
	 (cfg->rndfile == NULL) ? "(default)" : cfg->rndfile);

  printf("Prior used:                        %s\n",
	 (cfg->prifile == NULL) ? "(default)" : cfg->prifile);

  printf("Effective sequence # calculation:  ");
  if      (cfg->eff_strategy == EFF_NONE)      puts("none; using actual seq #");
  else if (cfg->eff_strategy == EFF_USERSET)   printf("set to %.2f\n", cfg->eff_nseq);
  else if (cfg->eff_strategy == EFF_NCLUSTERS) {
    puts("# single-linkage clusters");
    printf("  by single-linkage clustering at: %.2f identity\n", cfg->eidlevel);
  }
  else if (cfg->eff_strategy == EFF_ENTROPY) {
    puts("entropy targeting");
    if (cfg->eloss_set)
      printf("  mean target entropy loss:        %.2f bits\n", cfg->eloss);
    else
      printf("  mean target entropy loss:        (default)\n");
  }

  printf("Relative sequence weighting:       ");
  if      (cfg->w_strategy == WGT_NONE)   puts("none");
  else if (cfg->w_strategy == WGT_GSC)    puts("G/S/C tree weights");
  else if (cfg->w_strategy == WGT_BLOSUM) printf("BLOSUM filter at %.2f id\n", cfg->widlevel);
  else if (cfg->w_strategy == WGT_PB)     puts("Henikoff position-based");
  else if (cfg->w_strategy == WGT_VORONOI)puts("Sibbald/Argos Voronoi");
  else if (cfg->w_strategy == WGT_ME)     puts("Maximum entropy");

  printf("New HMM file:                      %s %s\n",
	 hmmfile, cfg->do_append? "[appending]" : "");
  if (cfg->cfile != NULL)
    printf("Count vectors saved to:            %s\n", cfg->cfile);
  if (cfg->align_ofile != NULL)
    printf("Annotated alignment(s) resaved to: %s\n", cfg->align_ofile);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
}



/* set_alphabet():
 * 
 * Called if alphabet wasn't set explicitly by --amino or --nucleic
 * on the command line.
 * Examine the multiple seq alignment <msa>, make a guess at the
 * alphabet, and set the globals accordingly.
 */
static void
set_alphabet(MSA *msa)
{
  printf("%-40s ... ", "Determining alphabet");
  fflush(stdout);
  DetermineAlphabet(msa->aseq, msa->nseq);
  if      (Alphabet_type == hmmNUCLEIC) puts("done. [DNA]");
  else if (Alphabet_type == hmmAMINO)   puts("done. [protein]");
  else                                  puts("done.");
}


/* set_entropy_target()
 * 
 * Called when we're using entropy-weighting to determine
 * effective sequence number.
 *   eloss_set: TRUE/FALSE, whether we set a target entropy loss
 *              w/ --eloss <x> on the command line;
 *   eloss:     if eloss_set was TRUE, what was the setting;
 *   mode:      alignment mode, example: P7_FS_MODE;
 *              affects choice of default eloss.
 *   randomseq: the null model, for calculating background
 *              sequence entropy that eloss is measured against.
 *              
 * Returns <etarget>, the target entropy of the model:
 *   H(background) - eloss.
 *                   
 * We have only tested default entropy targets for fs and ls mode, the
 * usual modes. We assume that -s behaves like fs mode, and -g behaves
 * like -s mode.  The default numbers are in config.h; they come from
 * LSJ's optimizations on the ASTRAL benchmark.
 */
static float
set_entropy_target(int eloss_set, float eloss, int mode, float *randomseq)
{
  float etarget;

  /* explicit setting? */
  if (eloss_set) return FEntropy(randomseq, Alphabet_size) - eloss;

  /* if not, protein defaults: */
  if (Alphabet_type == hmmAMINO) 
    {
      if   (mode == P7_FS_MODE || mode == P7_SW_MODE)
	etarget = FEntropy(randomseq, Alphabet_size) - ENTROPYLOSS_FS_AMINO_DEFAULT;
      else if (mode == P7_LS_MODE || mode == P7_S_MODE)
	etarget = FEntropy(randomseq, Alphabet_size) - ENTROPYLOSS_LS_AMINO_DEFAULT;
      else
	Die("I don't have a entropy loss default for alignment mode %d\n", mode);
      return etarget;
    }

  /* we don't have defaults available for DNA, or other alphabets for that matter */
  Die("\
--effent: entropy loss targeting:\n\
Default entropy loss targets are only available for protein models.\n\
To use --effent on DNA models (or other alphabets), you must set\n\
--eloss <x> explicitly, in addition to selecting --effent.\n");
  /*NOTREACHED*/
  return 0.;
}



/* set_relative_weights()
 * 
 * Given a multiple sequence alignment;
 * determine the relative sequence weights for it, and
 * set msa->wgt array accordingly.
 */
static void
set_relative_weights(MSA *msa, struct p7config_s *cfg)
{
  if (cfg->w_strategy == WGT_GSC     ||
      cfg->w_strategy == WGT_BLOSUM  ||
      cfg->w_strategy == WGT_VORONOI ||
      cfg->w_strategy == WGT_PB)
    {
      printf("%-40s ... ", "Weighting sequences heuristically");
      fflush(stdout);

      if (cfg->w_strategy != WGT_PB && msa->nseq >= cfg->pbswitch)
	{
	  printf("[big alignment! doing PB]... ");
	  PositionBasedWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
	}
      else if (cfg->w_strategy == WGT_GSC)
	GSCWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
      else if (cfg->w_strategy == WGT_BLOSUM)
	BlosumWeights(msa->aseq, msa->nseq, msa->alen, cfg->widlevel, msa->wgt);
      else if (cfg->w_strategy == WGT_PB)
	PositionBasedWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
      else if (cfg->w_strategy ==  WGT_VORONOI)
	VoronoiWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
    }
  printf("done.\n");
}



/* Function:  tag_candidate_seq_fragments()
 * Incept:    SRE, Fri May  6 13:32:20 2005 [St. Louis]
 *
 * Purpose:   Heuristic identification of candidate sequence
 *            "fragments": any individual seq with a length <
 *            <thresh> * mean length of all seqs is tagged
 *            as a candidate fragment. 
 *            
 *            The modelmakers then use this information. Any seq
 *            tagged as a candidate fragment that shows
 *            a B->D+->Mk internal entry path or a Mk->D+->E
 *            internal exit path has those transitions ignored,
 *            on the principle that it's just a local alignment,
 *            not a true deletion. Delete transitions on internal
 *            entry/exit paths are only counted on sequences that
 *            are not identified as candidate fragments.
 *            
 *            The mean length is a weighted average, taking relative
 *            seq weights into account.
 *            
 * Args:      msa    - the sequence alignment
 *            thresh - threshold fraction; rlen < thresh*mean rlen 
 *                      defines candidate fragment.
 *            rlen   - RETURN: raw (unaligned) lengths of seqs in residues,
 *                      [0..nseq-1]
 *                  old    (Caller allocates this for at least msa->nseq integers).
 *            isfrag - RETURN: TRUE/FALSE flags for whether seq i is a candidate
 *                      fragment or not; [0..nseq-1]
 *                      Caller allocates this for at least msa->nseq chars.
 *
 * Returns:   # of fragments that were flagged.
 */
static int
tag_candidate_seq_fragments(MSA *msa, float thresh, int *rlen, char *isfrag) 
{
  int i;
  float mean;
  float totwgt;
  int nfrags;

  printf("%-40s ... ", "Tagging putative sequence fragments");

  /* Calculate lengths of each seq, and mean length */
  mean   = 0.;
  totwgt = 0.;
  for (i = 0; i < msa->nseq; i++)
    {
      rlen[i] = DealignedLength(msa->aseq[i]);
      mean   += msa->wgt[i] * (float) rlen[i];
      totwgt += msa->wgt[i];
    }
  mean /= totwgt;

  /* any seq < threshold * mean is identified as a fragment */
  nfrags = 0;
  for (i = 0; i < msa->nseq; i++)
    if (rlen[i] < (int) (thresh * mean))
      { isfrag[i] = TRUE; nfrags++; }
    else
      isfrag[i] = FALSE;
    
  printf("done.\n");
  return nfrags;
}


/* construct_model()
 * 
 * Determine the architecture, collect counts; return 
 * a new HMM containing counts.
 */
static void
construct_model(struct p7config_s *cfg, MSA *msa, unsigned char **dsq, char *isfrag, 
		struct plan7_s **ret_hmm, struct p7trace_s ***ret_tr)
{
  int checksum;

  printf("%-40s ... ", "Constructing model architecture");
  fflush(stdout);

  checksum = GCGMultchecksum(msa->aseq, msa->nseq);
  if (cfg->c_strategy == P7_FAST_CONSTRUCTION)
    P7Fastmodelmaker(msa, dsq, isfrag, cfg->symfrac, ret_hmm, ret_tr);
  else if (cfg->c_strategy == P7_HAND_CONSTRUCTION)
    P7Handmodelmaker(msa, dsq, isfrag, ret_hmm, ret_tr);
  else 
    Die("no such construction strategy");
  (*ret_hmm)->checksum = checksum;
  printf("done.\n");
}


static void
set_effective_seqnumber(struct p7config_s *cfg, MSA *msa, struct plan7_s *hmm, struct p7prior_s *pri)
{
  if (cfg->eff_strategy == EFF_NONE)
    {
      printf("%-40s ... ", "Set effective seq # to just nseq");
      cfg->eff_nseq     = (float) msa->nseq;
    }

  else if (cfg->eff_strategy == EFF_USERSET)
    {
      printf("%-40s ... ", "Set effective seq #");
    }

  else if (cfg->eff_strategy == EFF_NCLUSTERS)
    {
      float *wgt;		/* dummy weights array to feed BlosumWeights*/
      wgt = MallocOrDie(sizeof(float) * msa->nseq);
      printf("%-40s ... ", "Determining eff seq # by clustering");
      fflush(stdout);
      BlosumWeights(msa->aseq, msa->nseq, msa->alen, cfg->eidlevel, wgt);
      cfg->eff_nseq = FSum(wgt, msa->nseq);
      free(wgt);
    }

  else if (cfg->eff_strategy == EFF_ENTROPY) 
    {
      printf("%-40s ... ", "Determining eff seq # by entropy target");
      fflush(stdout);
      cfg->eff_nseq = Eweight(hmm, pri, (float) msa->nseq, cfg->etarget);
    } 

  else
    Die("no effective seq # strategy: shouldn't happen");

  Plan7Rescale(hmm, cfg->eff_nseq / (float) msa->nseq);
  printf("done. [%.1f]\n", cfg->eff_nseq);
}


/* Function: print_all_scores()
 * 
 * Purpose:  For each training sequence, print its score under
 *           the final model.
 *           
 * Args:     fp   - where to print the output (usu. stdout)
 *           hmm  - newly constructed HMM, with prob's.
 *           dsq  - digitized unaligned training sequences.
 *           msa  - alignment and associated info 
 *           tr   - array of tracebacks
 *           
 * Return:   (void)                         
 */
static void
print_all_scores(FILE *fp, struct plan7_s *hmm,
		 unsigned char **dsq, MSA *msa, struct p7trace_s **tr)
{
  int idx;			/* counter for sequences */

  				/* make sure model scores are ready */
  P7Logoddsify(hmm, TRUE);
				/* header */
  fputs("**\n", fp);
  fputs("Individual training sequence scores:\n", fp);
				/* score for each sequence */
  for (idx = 0; idx < msa->nseq; idx++) 
    {
      fprintf(fp, "%7.2f %-12s %s\n", 
	      P7TraceScore(hmm, dsq[idx], tr[idx]),
	      msa->sqname[idx],
	      (MSAGetSeqDescription(msa,idx) != NULL) ? 
	       MSAGetSeqDescription(msa,idx) : "");
      P7PrintTrace(fp, tr[idx], hmm, dsq[idx]);
    }
  fputs("\n", fp);
}



/* Function: save_countvectors()
 * 
 * Purpose:  Save emission/transition count vectors to a file.
 *           Used for gathering the data on which to train a
 *           prior (e.g. mixture Dirichlet, etc.)
 *           
 *           The format of the file is one vector per line:
 *           M <f> <f> ...: 20 match emission counts in order AC..WY.
 *                          followed by two chars of CS, CA annotation.
 *           I <f> <f> ...: 20 insert emission counts in order AC..WY.
 *                          followed by two chars of CS, CA annotation.
 *           T <f> <f> ...: 7 transition counts in order TMM, TMI, TMD, 
 *                            TIM, TII, TDM, TDD. (see structs.h)
 *                            followed by four chars of structure
 *                            annotation: CS, CS of M+1; CA, CA of M+1. 
 *           
 * Args:     cfg    - config includes cfp, open counts file, and cfile, its name
 *           name   - name of alignment or HMM to associate with these vectors
 *           hmm    - counts-based HMM
 */
static void
save_countvectors(FILE *cfp, char *cfile, char *name, struct plan7_s *hmm)
{
  int k, x;

  printf("%-40s ... ", "Saving count vector file");
  fflush(stdout);
				/* match emission vectors */
  for (k = 1; k <= hmm->M; k++)
    {
      fputs("M ", cfp);
      for (x = 0; x < Alphabet_size; x++)
	fprintf(cfp, "%8.2f ", hmm->mat[k][x]);

      fprintf(cfp, "%15s %6d %6d ", name, hmm->map[k], k);
      if ((hmm->flags & PLAN7_CS) && hmm->flags & PLAN7_CA)
	fprintf(cfp, "%c %c", hmm->cs[k], hmm->ca[k]);
      else
	fputs("- -", cfp);
      fputs("\n", cfp);
    }
				/* insert emission vectors */
  for (k = 1; k < hmm->M; k++)
    {
      fputs("I ", cfp);
      for (x = 0; x < Alphabet_size; x++)
	fprintf(cfp, "%8.2f ", hmm->ins[k][x]);

      fprintf(cfp, "%15s %6d %6d ", name, hmm->map[k], k);
      if ((hmm->flags & PLAN7_CS) && hmm->flags & PLAN7_CA)
	fprintf(cfp, "%c %c", hmm->cs[k], hmm->ca[k]);
      else
	fputs("- -", cfp);

      fputs("\n", cfp);
    }
				/* transition vectors */
    for (k = 1; k < hmm->M; k++)
    {
      fputs("T ", cfp);

      for (x = 0; x < 7; x++)
	fprintf(cfp, "%8.2f ", hmm->t[k][x]);

      fprintf(cfp, "%15s %6d %6d ", name, hmm->map[k], k);
      if ((hmm->flags & PLAN7_CS) && hmm->flags & PLAN7_CA)
	fprintf(cfp, "%c %c %c %c", 
		hmm->cs[k], hmm->cs[k+1], 
		hmm->ca[k], hmm->ca[k+1]);
      else
	fputs("- -", cfp);
      fputs("\n", cfp);
    }
    printf("done. [%s]\n", cfile);
}


/* Function: position_average_score()
 * Date:     Wed Dec 31 09:36:35 1997 [StL]
 * 
 * Purpose:  Calculate scores from tracebacks, keeping them
 *           in a position specific array. The final array
 *           is normalized position-specifically too, according
 *           to how many sequences contributed data to this
 *           position. Used for compensating for sequence 
 *           fragments in ME and MD score optimization. 
 *           Very much ad hoc.
 *           
 *           Code related to (derived from) TraceScore().
 *           
 * Args:     hmm       - HMM structure, scores valid
 *           dsq       - digitized unaligned sequences
 *           wgt       - weights on the sequences
 *           nseq      - number of sequences
 *           tr        - array of nseq tracebacks that aligns each dsq to hmm
 *           pernode   - RETURN: [0]1..M array of position-specific avg scores
 *           ret_avg   - RETURN: overall average full-length, one-domain score
 *           
 * Return:   1 on success, 0 on failure.          
 *           pernode is malloc'ed [0]1..M by CALLER and filled here.
 */
static void
position_average_score(struct plan7_s    *hmm, 
		       unsigned char    **dsq, 
		       float             *wgt,
		       int                nseq,
		       struct p7trace_s **tr,
		       float             *pernode,
		       float             *ret_avg)
{
  unsigned char sym;
  int    pos;                   /* position in seq */
  int    tpos;                  /* position in trace/state sequence */
  float *counts;                /* counts at each position */
  float  avg;                   /* RETURN: average overall */
  int    k;                     /* counter for model position */
  int    idx;                   /* counter for sequence number */

  /* Allocations
   */
  counts = MallocOrDie ((hmm->M+1) * sizeof(float));
  FSet(pernode, hmm->M+1, 0.);
  FSet(counts,  hmm->M+1, 0.);

  /* Loop over traces, accumulate weighted scores per position
   */
  for (idx = 0; idx < nseq; idx++)
    for (tpos = 0; tpos < tr[idx]->tlen; tpos++)
      {
	pos = tr[idx]->pos[tpos];
	sym = dsq[idx][tr[idx]->pos[tpos]];
	k   = tr[idx]->nodeidx[tpos];

	/* Counts: how many times did we use this model position 1..M?
         * (weighted)
	 */
	if (tr[idx]->statetype[tpos] == STM || tr[idx]->statetype[tpos] == STD)
	  counts[k] += wgt[idx];

	/* Emission scores.
	 */
	if (tr[idx]->statetype[tpos] == STM) 
	  pernode[k] += wgt[idx] * Scorify(hmm->msc[sym][k]);
	else if (tr[idx]->statetype[tpos] == STI) 
	  pernode[k] += wgt[idx] * Scorify(hmm->isc[sym][k]);
	
	/* Transition scores.
	 */
	if (tr[idx]->statetype[tpos] == STM ||
	    tr[idx]->statetype[tpos] == STD ||
	    tr[idx]->statetype[tpos] == STI)
	  pernode[k] += wgt[idx] * 
	    Scorify(TransitionScoreLookup(hmm, tr[idx]->statetype[tpos], tr[idx]->nodeidx[tpos],
					  tr[idx]->statetype[tpos+1],tr[idx]->nodeidx[tpos+1]));
      }

  /* Divide accumulated scores by accumulated weighted counts
   */
  avg = 0.;
  for (k = 1; k <= hmm->M; k++)
    {
      pernode[k] /= counts[k];
      avg += pernode[k];
    }
  
  free(counts);
  *ret_avg = avg;
  return;
}


/* Function: frag_trace_score()
 * Date:     SRE, Wed Dec 31 10:03:47 1997 [StL]
 * 
 * Purpose:  Allow MD/ME optimization to be used for alignments
 *           that include fragments and multihits -- estimate a full-length
 *           per-domain score.
 *           
 *
 *           
 * Return:   "corrected" score.
 */
static float
frag_trace_score(struct plan7_s *hmm, unsigned char *dsq, struct p7trace_s *tr, 
                 float *pernode, float expected)
{
  float sc;			/* corrected score  */
  float fragexp;		/* expected score for a trace like this */
  int   tpos;			/* position in trace */

                                /* get uncorrected score */
  sc = P7TraceScore(hmm, dsq, tr);

                               /* calc expected score for trace like this */
  fragexp = 0.;
  for (tpos = 0; tpos < tr->tlen; tpos++)
    if (tr->statetype[tpos] == STM || tr->statetype[tpos] == STD)
      fragexp += pernode[tr->nodeidx[tpos]];

				/* correct for multihits */
  fragexp /= (float) TraceDomainNumber(tr);

                                /* extrapolate to full-length, one-hit score */
  sc = sc * expected / fragexp;
  return sc;
}


/* Function: maximum_entropy()
 * Date:     SRE, Fri Jan  2 10:56:00 1998 [StL]
 * 
 * Purpose:  Optimizes a model according to maximum entropy weighting.
 *           See Krogh and Mitchison (1995).
 *
 *           [Actually, we do minimum relative entropy, rather than
 *           maximum entropy. Same thing, though we refer to "ME"
 *           weights and models. The optimization is a steepest
 *           descents minimization of the relative entropy.]
 *           
 *           Expects to be called shortly after a Maxmodelmaker()
 *           or Handmodelmaker(), so that both a new model architecture
 *           (with MAP parameters) and fake tracebacks are available.
 *           
 *           Prints a summary of optimization progress to stdout.
 *           
 * Args:     hmm     - model. allocated, set with initial MAP parameters.
 *           dsq     - dealigned digitized seqs the model is based on
 *           ainfo   - extra info for aseqs
 *           nseq    - number of aseqs
 *           eff_nseq- effective sequence number; weights normalize up to this.
 *           prior   - prior distributions for parameterizing model
 *           tr      - array of fake traces for each sequence        
 *           
 * Return:   (void)
 *           hmm changed to an ME HMM
 *           ainfo changed, contains ME weights          
 */
static void
maximum_entropy(struct plan7_s *hmm, unsigned char **dsq, MSA *msa,
		float eff_nseq, struct p7prior_s *prior, struct p7trace_s **tr)
{
  float *wgt;                  /* current best set of ME weights   */
  float *new_wgt;              /* new set of ME weights to try     */
  float *sc;                    /* log-odds score of each sequence */
  float *grad;                  /* gradient                        */
  float  epsilon;               /* steepness of descent            */
  float  relative_entropy;      /* current best relative entropy   */
  float  new_entropy;           /* relative entropy at new weights */
  float  last_new_entropy;      /* last new_entropy we calc'ed     */
  float  use_epsilon;           /* current epsilon value in use    */
  int    idx;                   /* counter over sequences          */
  int    i1, i2;		/* counters for iterations         */

  float  converge_criterion;
  float  minw, maxw;            /* min, max weight                 */
  int    posw, highw;           /* number of positive weights      */
  float  mins, maxs, avgs;      /* min, max, avg score             */
  float *pernode;               /* expected score per node of HMM  */
  float  expscore;              /* expected score of complete HMM  */
  int    max_iter;		/* bulletproof against infinite loop bugs */

  epsilon  = 0.2;                /* works fine */
  max_iter = 666;
  
  /* Allocations
   */
  sc      = MallocOrDie (sizeof(float) * msa->nseq);
  wgt     = MallocOrDie (sizeof(float) * msa->nseq);
  new_wgt = MallocOrDie (sizeof(float) * msa->nseq);
  grad    = MallocOrDie (sizeof(float) * msa->nseq);
  pernode = MallocOrDie (sizeof(float) * (hmm->M+1));

  /* Initialization. Start with all weights == 1.0.
   * Find relative entropy and gradient.
   */
  Plan7SWConfig(hmm, 0.5, 0.5);
  P7Logoddsify(hmm, TRUE);

  FSet(wgt, msa->nseq, eff_nseq / (float) msa->nseq);
  position_average_score(hmm, dsq, wgt, msa->nseq, tr, pernode,&expscore);
  for (idx = 0; idx < msa->nseq; idx++) 
    sc[idx] = frag_trace_score(hmm, dsq[idx], tr[idx], pernode, expscore);
  relative_entropy = FSum(sc, msa->nseq) / (float) msa->nseq;
  for (idx = 0; idx < msa->nseq; idx++)
    grad[idx] = relative_entropy - sc[idx];

  printf("\n%-40s ...\n", "Maximum entropy weighting, iterative");
  printf("iter avg-sc min-sc max-sc min-wgt max-wgt +wgt ++wgt rel.ent convergence\n");
  printf("---- ------ ------ ------ ------- ------- ---- ----- ------- -----------\n");
  mins = maxs = avgs = sc[0];
  for (idx = 1; idx < msa->nseq; idx++)
    {
      if (sc[idx] < mins) mins = sc[idx];
      if (sc[idx] > maxs) maxs = sc[idx];
      avgs += sc[idx];
    }
  avgs /= (float) msa->nseq;
  printf("%4d %6.1f %6.1f %6.1f %7.2f %7.2f %4d %5d %7.2f %8s\n",
         0, avgs, mins, maxs, 1.0, 1.0, msa->nseq, 0, relative_entropy, "-");

  
  /* Steepest descents optimization;
   * iterate until relative entropy converges.
   */
  i1 = 0;
  while (++i1 < max_iter)
    {
      /* Gradient gives us a line of steepest descents.
       * (Roughly speaking, anyway. We actually have a constraint
       * that weights are nonnegative and normalized, and the
       * gradient doesn't take these into account.)
       * Look along this line, a distance of epsilon * gradient:
       * if new point is better, accept; if new point is worse,
       * move back along the line by half the distance and re-evaluate.
       */
      use_epsilon = epsilon;
      new_entropy = relative_entropy + 1.0;    /* just ensure new > old */

      i2 = 0; 
      while (new_entropy > relative_entropy && ++i2 < max_iter)
        {
          last_new_entropy = new_entropy;

                                /* find a new point in weight space */
          for (idx = 0; idx < msa->nseq; idx++)
            {
              new_wgt[idx] = wgt[idx] + use_epsilon * grad[idx];
              if (new_wgt[idx] < 0.) new_wgt[idx] = 0.0;
            }
          FNorm(new_wgt, msa->nseq);
          FScale(new_wgt, msa->nseq, eff_nseq);

                                /* Make new HMM using these weights */
          ZeroPlan7(hmm);
          for (idx = 0; idx < msa->nseq; idx++)
            P7TraceCount(hmm, dsq[idx], new_wgt[idx], tr[idx]);
          P7PriorifyHMM(hmm, prior);

  
                                /* Evaluate new point */
	  Plan7SWConfig(hmm, 0.5, 0.5);
	  P7Logoddsify(hmm, TRUE);
	  position_average_score(hmm, dsq, new_wgt, msa->nseq, tr, pernode, &expscore);
          for (idx = 0; idx < msa->nseq; idx++) 
	    sc[idx]      = frag_trace_score(hmm, dsq[idx], tr[idx], pernode, expscore);
	  new_entropy = FDot(sc, new_wgt, msa->nseq) / (float) msa->nseq;

          use_epsilon /= 2.0;
	  /* Failsafe: we're not converging. Set epsilon to zero,
	   * do one more round.
	   */
	  if (use_epsilon < 1e-6) use_epsilon = 0.0; 
	  if (use_epsilon == 0.0) break;
          
          /* Failsafe: avoid infinite loops. Sometimes the
             new entropy converges without ever being better 
             than the previous point, probably as a result
             of minor roundoff error. */
          if (last_new_entropy == new_entropy) break;
        }
      if (i2 == max_iter) printf("   -- exceeded maximum iterations; giving up --\n");

      /* Evaluate convergence before accepting the new weights;
       * then, accept the new point and evaluate the gradient there.
       */
      converge_criterion = fabs((relative_entropy-new_entropy)/relative_entropy);
      relative_entropy = new_entropy;
      FCopy(wgt, new_wgt, msa->nseq);
      for (idx = 0; idx < msa->nseq; idx++)
	grad[idx] = relative_entropy - sc[idx];

      /* Print some statistics about this iteration
       */
      mins = maxs = avgs = sc[0];
      minw = maxw = wgt[0];
      posw = (wgt[0] > 0.0) ? 1 : 0;
      highw = (wgt[0] > 1.0) ? 1 : 0;
      for (idx = 1; idx < msa->nseq; idx++)
        {
          if (sc[idx] < mins) mins = sc[idx];
          if (sc[idx] > maxs) maxs = sc[idx];
          if (wgt[idx] < minw) minw = wgt[idx];
          if (wgt[idx] > maxw) maxw = wgt[idx];
          if (wgt[idx] > 0.0)  posw++;
          if (wgt[idx] > 1.0)  highw++;
          avgs += sc[idx];
        }
      avgs /= (float) msa->nseq;
      printf("%4d %6.1f %6.1f %6.1f %7.2f %7.2f %4d %5d %7.2f %8.5f\n",
             i1, 
             avgs, mins, maxs, 
             minw, maxw, posw, highw,
             relative_entropy, converge_criterion);
      
      if (converge_criterion < 1e-5) break;
    }
  if (i1 == max_iter) printf("   -- exceeded maximum iterations; giving up --\n");

  /* Renormalize weights to sum to eff_nseq, and save.
   */
  FNorm(wgt, msa->nseq);
  FScale(wgt, msa->nseq, (float) eff_nseq);
  FCopy(msa->wgt, wgt, msa->nseq);
			/* Make final HMM using these adjusted weights */
  ZeroPlan7(hmm);
  for (idx = 0; idx < msa->nseq; idx++)
    P7TraceCount(hmm, dsq[idx], wgt[idx], tr[idx]);
  P7PriorifyHMM(hmm, prior);
                                
  /* Cleanup and return
   */
  free(pernode);
  free(new_wgt);
  free(grad);
  free(wgt);
  free(sc);
  printf("----------------------------------------------\n\n");
  return;
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
set_model_name(struct plan7_s *hmm, char *setname, char *msa_name, char *alifile, int nali)
{
  char *name;

  printf("%-40s ... ", "Set model name, record commandline");
  fflush(stdout);

  if (nali == 0)		/* first (only?) HMM in file:  */
    {
      if      (setname  != NULL) name = Strdup(setname);
      else if (msa_name != NULL) name = Strdup(msa_name);
      else                       name = FileTail(alifile, TRUE);
    }
  else
    {
      if (setname != NULL) 
	Die("Oops. Wait. You can't use -n with an alignment database.");
      else if (msa_name != NULL) name = Strdup(msa_name);
      else
	Die("Oops. Wait. I need name annotation on each alignment.\n");
    }
  Plan7SetName(hmm, name);
  free(name);
  printf("done. [%s]\n", hmm->name);
}

/* print_statistics()
 * 
 * Purpose:  Given a newly constructed HMM and the tracebacks
 *           of the sequences it was trained on, print out all
 *           the interesting information at the end of hmmbuild
 *           runs that convinces the user we actually
 *           did something.
 *           
 * Args:     fp   - where to send the output (stdout, usually)
 *           hmm  - the new HMM, probability form
 *           dsq  - digitized training seqs
 *           nseq - number of dsq's
 *           tr   - array of tracebacks for dsq
 */
static void
print_statistics(FILE *fp, struct plan7_s *hmm, unsigned char **dsq, int nseq,
		 struct p7trace_s **tr)
{
  int   idx;			/* counter for sequences                */
  float score;			/* an individual trace score            */
  float total, best, worst;	/* for the avg. and range of the scores */
  float sqsum, stddev;		/* for the std. deviation of the scores */

  P7Logoddsify(hmm, TRUE);
				/* find individual trace scores */
  score = P7TraceScore(hmm, dsq[0], tr[0]);
  total = best = worst = score;
  sqsum = score * score;
  for (idx = 1; idx < nseq; idx++) {
    /* P7PrintTrace(stdout, tr[idx], hmm, dsq[idx]); */
    score  = P7TraceScore(hmm, dsq[idx], tr[idx]);
    total += score;
    sqsum += score * score;
    if (score > best)  best = score;
    if (score < worst) worst = score;
  }
  if (nseq > 1) {
    stddev = (sqsum - (total * total / (float) nseq)) / ((float) nseq - 1.);
    stddev = (stddev > 0) ? sqrt(stddev) : 0.0;
  } else stddev = 0.0;
				/* print out stuff. */
  fprintf(fp, "Average score:  %10.2f bits\n", total / (float) nseq);
  fprintf(fp, "Minimum score:  %10.2f bits\n", worst);
  fprintf(fp, "Maximum score:  %10.2f bits\n", best);
  fprintf(fp, "Std. deviation: %10.2f bits\n", stddev);
}


/* save_hmmbuild_alignment()
 * 
 * hmmbuild modifies the input alignment in a number of ways:
 *   1. Implied D->I and I->D transitions are removed, by small
 *      alterations in the aligned residues;
 *   2. Weights have been determined (relative and absolute);
 *   3. Reference line is added, marking w/ x those columns that
 *      were called consensus (M/D).
 *
 * For caller's information, optionally output the alignment 
 * to a save file.
 */
static void
save_hmmbuild_alignment(FILE *alignfp, MSA *msa, unsigned char **dsq, struct plan7_s *hmm,
			struct p7trace_s **tr)
{
  MSA    *new_msa;
  SQINFO *sqinfo;
  int     idx;

  printf("%-40s ... ", "Saving annotated alignment");
  fflush(stdout);
  sqinfo  = MSAToSqinfo(msa);
  new_msa = P7Traces2Alignment(dsq, sqinfo, msa->wgt, msa->nseq, 
			       hmm->M, tr, FALSE);
  WriteStockholm(alignfp, new_msa);
  MSAFree(new_msa);
  for (idx = 0; idx < msa->nseq; idx++)
    FreeSequence(NULL, &(sqinfo[idx]));
  free(sqinfo);
  printf("done.\n");
}






/************************************************************
 * @LICENSE@
 ************************************************************/

