/************************************************************
 * @LICENSE@
 ************************************************************/

/* hmmbuild.c
 * SRE, Mon Nov 18 12:41:29 1996
 *
 * main() for HMM construction from an alignment.
 * CVS $Id$
 */

#include "config.h"		/* compile-time configuration constants */
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */

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
  Alternative search algorithm styles: (default: hmmls domain alignment)\n\
   -f     : multi-hit local (hmmfs style)\n\
   -g     : global alignment (hmms style, Needleman/Wunsch)\n\
   -s     : local alignment (hmmsw style, Smith/Waterman)\n\
";

static char experts[] = "\
  Alternative model construction strategies: (default: MAP)\n\
   --fast        : Krogh/Haussler fast heuristic construction (see --gapmax)\n\
   --hand        : manual construction (requires annotated alignment)\n\
\n\
  Expert customization of parameters and priors:\n\
   --null  <f>   : read null (random sequence) model from <f>\n\
   --pam   <f>   : heuristic PAM-based prior, using BLAST PAM matrix in <f>\n\
   --prior <f>   : read Dirichlet prior parameters from <f>\n\
\n\
  Alternative sequence weighting strategies: (default: GSC weights)\n\
   --wblosum     : Henikoff simple filter weights (see --idlevel)\n\
   --wgsc        : Gerstein/Sonnhammer/Chothia tree weights (default)\n\
   --wme         : maximum entropy (ME)\n\
   --wpb         : Henikoff position-based weights\n\
   --wvoronoi    : Sibbald/Argos Voronoi weights\n\
   --wnone       : don't do any weighting\n\
   --noeff       : don't use effective sequence number; just use nseq\n\
   --pbswitch <n>: set switch from GSC to position-based wgts at > n seqs\n\
\n\
  Forcing an alphabet: (normally autodetected)\n\
   --amino   : override autodetection, assert that seqs are protein\n\
   --nucleic : override autodetection, assert that seqs are DNA/RNA\n\
\n\
  Other expert options:\n\
   --archpri <x> : set architecture size prior to <x> {0.85} [0..1]\n\
   --binary      : save the model in binary format, not ASCII text\n\
   --cfile <f>   : save count vectors to <f>\n\
   --gapmax <x>  : max fraction of gaps in mat column {0.50} [0..1]\n\
   --idlevel <x> : set frac. id level used by eff. nseq and --wblosum {0.62}\n\
   --informat <s>: input alignment is in format <s>, not Stockholm\n\
   --pamwgt <x>  : set weight on PAM-based prior to <x> {20.}[>=0]\n\
   --swentry <x> : set S/W aggregate entry prob. to <x> {0.5}\n\
   --swexit <x>  : set S/W aggregate exit prob. to <x>  {0.5}\n\
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
  { "--archpri", FALSE, sqdARG_FLOAT }, 
  { "--binary",  FALSE, sqdARG_NONE  }, 
  { "--cfile",   FALSE, sqdARG_STRING},
  { "--fast",    FALSE, sqdARG_NONE},
  { "--gapmax",  FALSE, sqdARG_FLOAT },
  { "--hand",    FALSE, sqdARG_NONE},
  { "--idlevel", FALSE, sqdARG_FLOAT },
  { "--informat",FALSE, sqdARG_STRING },
  { "--noeff",   FALSE, sqdARG_NONE },
  { "--nucleic", FALSE, sqdARG_NONE },
  { "--null",    FALSE, sqdARG_STRING },
  { "--pam",     FALSE, sqdARG_STRING },
  { "--pamwgt",  FALSE, sqdARG_FLOAT },
  { "--pbswitch",FALSE, sqdARG_INT },
  { "--prior",   FALSE, sqdARG_STRING },
  { "--swentry", FALSE, sqdARG_FLOAT },
  { "--swexit",  FALSE, sqdARG_FLOAT },
  { "--verbose", FALSE, sqdARG_NONE  },
  { "--wgsc",    FALSE, sqdARG_NONE },
  { "--wblosum", FALSE, sqdARG_NONE },
  { "--wme",     FALSE, sqdARG_NONE },
  { "--wnone",   FALSE, sqdARG_NONE },
  { "--wpb",     FALSE, sqdARG_NONE },
  { "--wvoronoi",FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static void print_all_scores(FILE *fp, struct plan7_s *hmm, 
			     char **dsq, MSA *msa, struct p7trace_s **tr);
static void save_countvectors(FILE *cfp, char *name, struct plan7_s *hmm);
static void position_average_score(struct plan7_s *hmm, char **seq, float *wgt,
				   int nseq, struct p7trace_s **tr, float *pernode,
				   float *ret_avg);
static float frag_trace_score(struct plan7_s *hmm, char *dsq, struct p7trace_s *tr, 
			      float *pernode, float expected);
static void maximum_entropy(struct plan7_s *hmm, char **dsq, MSA *msa,
			    float eff_nseq, 
			    struct p7prior_s *prior, struct p7trace_s **tr);


int
main(int argc, char **argv) 
{
  char            *seqfile;     /* seqfile to read alignment from          */ 
  int              format;	/* format of seqfile                       */
  MSAFILE         *afp;         /* open alignment file                     */
  MSA             *msa;         /* a multiple sequence alignment           */
  char           **dsq;         /* digitized unaligned aseq's              */ 
  struct plan7_s  *hmm;         /* constructed HMM; written to hmmfile     */
  struct p7prior_s *pri;        /* Dirichlet priors to use                 */
  struct p7trace_s **tr;        /* fake tracebacks for aseq's              */ 
  char            *hmmfile;     /* file to write HMM to                    */
  FILE            *hmmfp;       /* HMM output file handle                  */
  char            *name;        /* name of the HMM                         */
  int              idx;		/* counter for sequences                   */
  float  randomseq[MAXABET];	/* null sequence model                     */
  float            p1;		/* null sequence model p1 transition       */
  int              nali;	/* count number of alignments/HMMs         */
  char             fpopts[3];   /* options to open a file with, e.g. "ab"  */
  int              checksum;	/* checksum of the alignment               */

  char *optname;                /* name of option found by Getopt()      */
  char *optarg;                 /* argument found by Getopt()            */
  int   optind;                 /* index in argv[]                       */

  enum p7_construction c_strategy; /* construction strategy choice        */
  enum p7_weight {		/* weighting strategy */
    WGT_NONE, WGT_GSC, WGT_BLOSUM, WGT_PB, WGT_VORONOI, WGT_ME} w_strategy;
  enum p7_config {              /* algorithm configuration strategy      */
    P7_BASE_CONFIG, P7_LS_CONFIG, P7_FS_CONFIG, P7_SW_CONFIG } cfg_strategy;
  float gapmax;			/* max frac gaps in mat col for -k       */
  int   overwrite_protect;	/* TRUE to prevent overwriting HMM file  */
  int   verbose;		/* TRUE to show a lot of output          */
  char *rndfile;		/* random sequence model file to read    */
  char *prifile;		/* Dirichlet prior file to read          */
  char *pamfile;		/* PAM matrix file for heuristic prior   */
  char *align_ofile;            /* name of output alignment file         */
  char *cfile;			/* output file for count vectors         */
  FILE *alignfp;                /* open filehandle for alignment resaves */
  FILE *cfp;                    /* open filehandle for count vector saves*/
  float archpri;		/* "architecture" prior on model size    */
  float pamwgt;			/* weight on PAM for heuristic prior     */
  int   do_append;		/* TRUE to append to hmmfile             */
  int   do_binary;		/* TRUE to write in binary format        */
  float blosumlevel;		/* BLOSUM frac id filtering level [0.62] */
  float swentry;		/* S/W aggregate entry probability       */
  float swexit;			/* S/W aggregate exit probability        */
  int   do_eff;			/* TRUE to set an effective seq number   */
  float eff_nseq;		/* effective sequence number             */
  int   pbswitch;		/* nseq >= this, switchover to PB weights*/
  char *setname;                /* NULL, or ptr to HMM name to set       */
  int   gapmax_set;		/* TRUE if gapmax was set on commandline */

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  format            = MSAFILE_UNKNOWN;        /* autodetect format by default. */
  c_strategy        = P7_MAP_CONSTRUCTION;
  w_strategy        = WGT_GSC;
  blosumlevel       = 0.62;
  cfg_strategy      = P7_LS_CONFIG;
  gapmax            = 0.5;
  overwrite_protect = TRUE;
  verbose           = FALSE;
  rndfile           = NULL;
  prifile           = NULL;
  pamfile           = NULL;
  align_ofile       = NULL;
  alignfp           = NULL;
  cfile             = NULL;
  cfp               = NULL;
  archpri           = 0.85;
  pamwgt            = 20.;
  Alphabet_type     = hmmNOTSETYET;	/* initially unknown */
  name              = NULL;
  do_append         = FALSE; 
  swentry           = 0.5;
  swexit            = 0.5;
  do_eff            = TRUE;
  do_binary         = FALSE;
  pbswitch          = 1000;
  setname           = NULL;
  gapmax_set        = FALSE;
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-f") == 0) cfg_strategy      = P7_FS_CONFIG;
    else if (strcmp(optname, "-g") == 0) cfg_strategy      = P7_BASE_CONFIG;
    else if (strcmp(optname, "-n") == 0) setname           = optarg; 
    else if (strcmp(optname, "-o") == 0) align_ofile       = optarg;
    else if (strcmp(optname, "-s") == 0) cfg_strategy      = P7_SW_CONFIG;
    else if (strcmp(optname, "-A") == 0) do_append         = TRUE; 
    else if (strcmp(optname, "-F") == 0) overwrite_protect = FALSE;
    else if (strcmp(optname, "--amino")   == 0) SetAlphabet(hmmAMINO);
    else if (strcmp(optname, "--archpri") == 0) archpri       = atof(optarg);
    else if (strcmp(optname, "--binary")  == 0) do_binary     = TRUE;
    else if (strcmp(optname, "--cfile")   == 0) cfile         = optarg;
    else if (strcmp(optname, "--fast")    == 0) c_strategy    = P7_FAST_CONSTRUCTION;
    else if (strcmp(optname, "--gapmax")  == 0) { gapmax      = atof(optarg); gapmax_set = TRUE; }
    else if (strcmp(optname, "--hand")    == 0) c_strategy    = P7_HAND_CONSTRUCTION;
    else if (strcmp(optname, "--idlevel") == 0) blosumlevel   = atof(optarg);
    else if (strcmp(optname, "--noeff")   == 0) do_eff        = FALSE;
    else if (strcmp(optname, "--nucleic") == 0) SetAlphabet(hmmNUCLEIC);
    else if (strcmp(optname, "--null")    == 0) rndfile       = optarg;
    else if (strcmp(optname, "--pam")     == 0) pamfile       = optarg;
    else if (strcmp(optname, "--pamwgt")  == 0) pamwgt        = atof(optarg);
    else if (strcmp(optname, "--pbswitch")== 0) pbswitch      = atoi(optarg);
    else if (strcmp(optname, "--prior")   == 0) prifile       = optarg;
    else if (strcmp(optname, "--swentry") == 0) swentry       = atof(optarg); 
    else if (strcmp(optname, "--swexit")  == 0) swexit        = atof(optarg); 
    else if (strcmp(optname, "--verbose") == 0) verbose       = TRUE;
    else if (strcmp(optname, "--wgsc")    == 0) w_strategy    = WGT_GSC;
    else if (strcmp(optname, "--wblosum") == 0) w_strategy    = WGT_BLOSUM; 
    else if (strcmp(optname, "--wme")     == 0) w_strategy    = WGT_ME;  
    else if (strcmp(optname, "--wpb")     == 0) w_strategy    = WGT_PB;  
    else if (strcmp(optname, "--wnone")   == 0) w_strategy    = WGT_NONE; 
    else if (strcmp(optname, "--wvoronoi")== 0) w_strategy    = WGT_VORONOI;
    else if (strcmp(optname, "--informat") == 0) {
      format = String2SeqfileFormat(optarg);
      if (format == MSAFILE_UNKNOWN) 
	Die("unrecognized sequence file format \"%s\"", optarg);
      if (! IsAlignmentFormat(format))
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

  hmmfile = argv[optind++];
  seqfile = argv[optind++]; 

  if (gapmax < 0. || gapmax > 1.) 
    Die("--gapmax must be a value from 0 to 1\n%s\n", usage);
  if (archpri < 0. || archpri > 1.)
    Die("--archpri must be a value from 0 to 1\n%s\n", usage);
  if (overwrite_protect && !do_append && FileExists(hmmfile))
    Die("HMM file %s already exists. Rename or delete it.", hmmfile); 
  if (overwrite_protect && align_ofile != NULL && FileExists(align_ofile))
    Die("Alignment resave file %s exists. Rename or delete it.", align_ofile); 
  if (gapmax_set && c_strategy  != P7_FAST_CONSTRUCTION)
    Die("using --gapmax only makes sense if you use --fast");

  /*********************************************** 
   * Preliminaries: open our files for i/o
   ***********************************************/

				/* Open the alignment */
  if ((afp = MSAFileOpen(seqfile, format, NULL)) == NULL)
    Die("Alignment file %s could not be opened for reading", seqfile);

				/* Open the HMM output file */
  if (do_append) strcpy(fpopts, "a");
  else           strcpy(fpopts, "w");
  if (do_binary) strcat(fpopts, "b");
  if ((hmmfp = fopen(hmmfile, fpopts)) == NULL)
    Die("Failed to open HMM file %s for %s\n", hmmfile, 
	do_append ? "appending" : "writing");

				/* Open the count vector save file */
  cfp = NULL;
  if (cfile != NULL)
    if ((cfp = fopen(cfile, "w")) == NULL)
      Die("Failed to open count vector file %s for writing\n", cfile);

				/* Open the alignment resave file */
  alignfp = NULL;
  if (align_ofile != NULL)
    if ((alignfp = fopen(align_ofile, "w")) == NULL)
      Die("Failed to open alignment resave file %s for writing\n", align_ofile);

  /*********************************************** 
   * Show the banner
   ***********************************************/

  HMMERBanner(stdout, banner);
  printf("Alignment file:                    %s\n", 
	 seqfile);
  printf("File format:                       %s\n", 
	 SeqfileFormat2String(afp->format));

  printf("Search algorithm configuration:    ");
  if      (cfg_strategy == P7_BASE_CONFIG)   puts("Global alignment (hmms)");
  else if (cfg_strategy == P7_SW_CONFIG)     {
    puts("Local  (hmmsw)");
    printf("S/W aggregate entry probability:   %.2f\n", swentry);
    printf("S/W aggregate exit probability:    %.2f\n", swexit);
  }
  else if (cfg_strategy == P7_LS_CONFIG)     puts("Multiple domain (hmmls)");
  else if (cfg_strategy == P7_FS_CONFIG)     {
    puts("Multiple local (hmmfs)");
    printf("S/W aggregate entry probability:   %.2f\n", swentry);
    printf("S/W aggregate exit probability:    %.2f\n", swexit);
  }

  printf("Model construction strategy:       ");
  if (c_strategy == P7_HAND_CONSTRUCTION)    puts("Manual, from #=RF annotation");
  else if (c_strategy==P7_FAST_CONSTRUCTION) printf("Fast/ad hoc (gapmax %.2f)\n", gapmax);
  else                                       printf("MAP (gapmax hint: %.2f)\n", gapmax);

  printf("Null model used:                   %s\n",
	 (rndfile == NULL) ? "(default)" : rndfile);

  printf("Prior used:                        %s\n",
	 (prifile == NULL) ? "(default)" : prifile);

  printf("Sequence weighting method:         ");
  if      (w_strategy == WGT_NONE)   puts("none");
  else if (w_strategy == WGT_GSC)    puts("G/S/C tree weights");
  else if (w_strategy == WGT_BLOSUM) printf("BLOSUM filter at %.2f id\n", blosumlevel);
  else if (w_strategy == WGT_PB)     puts("Henikoff position-based");
  else if (w_strategy == WGT_VORONOI)puts("Sibbald/Argos Voronoi");
  else if (w_strategy == WGT_ME)     puts("Maximum entropy");

  printf("New HMM file:                      %s %s\n",
	 hmmfile, do_append? "[appending]" : "");
  if (cfile != NULL)
    printf("Count vectors saved to:            %s\n", cfile);
  if (align_ofile != NULL)
    printf("Annotated alignment(s) resaved to: %s\n", align_ofile);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");


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

      /* Set up the alphabet globals:
       * either already set by --amino or --nucleic, or
       * we guess based on the first alignment we see
       */
      if (Alphabet_type == hmmNOTSETYET) 
	DetermineAlphabet(msa->aseq, msa->nseq);

      /* Do some initialization the first time through.
       * This code must be delayed until after we've seen the
       * first alignment, because we have to see the alphabet type first
       */
      if (nali == 0) 
	{
				/* Set up Dirichlet priors */
	  if (prifile == NULL)  pri = P7DefaultPrior();
	  else                  pri = P7ReadPrior(prifile);

	  if (pamfile != NULL)  PAMPrior(pamfile, pri, pamwgt);

				/* Set up the null/random seq model */
	  if (rndfile == NULL)  P7DefaultNullModel(randomseq, &p1);
	  else                  P7ReadNullModel(rndfile, randomseq, &p1);
	}

      /* Prepare unaligned digitized sequences for internal use 
       */
      DigitizeAlignment(msa, &dsq);
  
      /* In some respects we treat DNA more crudely for now;
       * for example, we can't do eff seq #, because it's
       * calibrated for protein.
       */
      if (Alphabet_type == hmmNUCLEIC) 
	do_eff = FALSE;	

      /* Determine "effective sequence number".
       * The BlosumWeights() routine is now an efficient O(N)
       * memory clustering algorithm that doesn't blow up on,
       * say, Pfam's GP120 alignment (13000+ sequences)
       */
      eff_nseq = (float) msa->nseq;
      if (do_eff)
	{
	  float *wgt;
	  printf("%-40s ... ", "Determining effective sequence number");
	  fflush(stdout);
				/* dummy weights array to feed BlosumWeights*/
	  wgt = MallocOrDie(sizeof(float) * msa->nseq);
	  BlosumWeights(msa->aseq, msa->nseq, msa->alen, blosumlevel, wgt);
	  eff_nseq = FSum(wgt, msa->nseq);

	  free(wgt);
	  printf("done. [%.0f]\n", eff_nseq);
	}


      /* Weight the sequences (optional),
       */
      if (w_strategy == WGT_GSC     || 
	  w_strategy == WGT_BLOSUM  || 
	  w_strategy == WGT_VORONOI ||
	  w_strategy == WGT_PB)
	{
	  printf("%-40s ... ", "Weighting sequences heuristically");
	  fflush(stdout);

	  if (w_strategy != WGT_PB && msa->nseq >= pbswitch)
	    {
	      printf("[big alignment! doing PB]... ");
	      PositionBasedWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
	    }
	  else if (w_strategy == WGT_GSC) 
	    GSCWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
	  else if (w_strategy == WGT_BLOSUM)
	    BlosumWeights(msa->aseq, msa->nseq, msa->alen, blosumlevel, msa->wgt);
	  else if (w_strategy == WGT_PB)
	    PositionBasedWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt);
	  else if (w_strategy ==  WGT_VORONOI)
	    VoronoiWeights(msa->aseq, msa->nseq, msa->alen, msa->wgt); 

	  printf("done.\n");
	}

      /* Set the effective sequence number (if do_eff is FALSE, eff_nseq 
       * was set to nseq).
       */
      FNorm(msa->wgt,  msa->nseq);
      FScale(msa->wgt, msa->nseq, eff_nseq);

      /* Build a model architecture.
       * If we're not doing MD or ME, that's all we need to do.
       * We get an allocated, counts-based HMM back.
       * 
       * Because the architecture algorithms are allowed to change
       * gap characters in the alignment, we have to calculate the
       * alignment checksum before we enter the algorithms.
       */
      printf("%-40s ... ", "Constructing model architecture");
      fflush(stdout);
      checksum = GCGMultchecksum(msa->aseq, msa->nseq);
      if (c_strategy == P7_FAST_CONSTRUCTION)
	P7Fastmodelmaker(msa, dsq, gapmax, &hmm, &tr);
      else if (c_strategy == P7_HAND_CONSTRUCTION)
	P7Handmodelmaker(msa, dsq, &hmm, &tr);
      else
	P7Maxmodelmaker(msa, dsq, gapmax, 
			pri, randomseq, p1, archpri, &hmm, &tr);
      hmm->checksum = checksum;
      printf("done.\n");

      /* Save the count vectors if asked. Used primarily for
       * making the data files for training priors.
       */
      if (cfile != NULL) 
	{
	  printf("%-40s ... ", "Saving count vector file");
	  fflush(stdout);
	  save_countvectors(cfp, 
			    (msa->name != NULL ? msa->name : "-"),
			    hmm); 
	  printf("done. [%s]\n", cfile);
	}

      /* Record the null model in the HMM;
       * add prior contributions in pseudocounts and renormalize.
       */
      printf("%-40s ... ", "Converting counts to probabilities");
      fflush(stdout);
      Plan7SetNullModel(hmm, randomseq, p1);
      P7PriorifyHMM(hmm, pri);
      printf("done.\n");

      /* Model configuration, temporary.
       * hmmbuild assumes that it's given an alignment of single domains,
       * and the alignment may contain fragments. So, for the purpose of
       * scoring the sequences (or, optionally, MD/ME weighting),
       * configure the model into hmmsw mode. Later we'll
       * configure the model according to how the user wants to
       * use it.
       */
      Plan7SWConfig(hmm, 0.5, 0.5);

      /* Do model-dependent "weighting" strategies.
       */
      if (w_strategy == WGT_ME)
	{
	  printf("\n%-40s ...\n", "Maximum entropy weighting, iterative");
	  maximum_entropy(hmm, dsq, msa, eff_nseq, pri, tr);
	  printf("----------------------------------------------\n\n");
	}

      /* Give the model a name.
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
      printf("%-40s ... ", "Setting model name, etc.");
      fflush(stdout);
      if (nali == 0)		/* first (only?) HMM in file:  */
	{
	  if      (setname != NULL)   name = Strdup(setname);
	  else if (msa->name != NULL) name = Strdup(msa->name);
	  else                        name = FileTail(seqfile, TRUE);
	}
      else
	{
	  if (setname != NULL) 
	    Die("Oops. Wait. You can't use -n with an alignment database.");
	  else if (msa->name != NULL) name = Strdup(msa->name);
	  else
	    Die("Oops. Wait. I need name annotation on each alignment.\n");
	}
      Plan7SetName(hmm, name);
      free(name);

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
       * like how/when we built it.
       */
      Plan7ComlogAppend(hmm, argc, argv);
      Plan7SetCtime(hmm);
      hmm->nseq = msa->nseq;
      printf("done. [%s]\n", hmm->name); 
   
      /* Print information for the user
       */
      printf("\nConstructed a profile HMM (length %d)\n", hmm->M);
      PrintPlan7Stats(stdout, hmm, dsq, msa->nseq, tr); 
      printf("\n");

      /* Configure the model for chosen algorithm
       */
      printf("%-40s ... ", "Finalizing model configuration");
      fflush(stdout);
      switch (cfg_strategy) {
      case P7_BASE_CONFIG:  Plan7GlobalConfig(hmm);              break;
      case P7_SW_CONFIG:    Plan7SWConfig(hmm, swentry, swexit); break;
      case P7_LS_CONFIG:    Plan7LSConfig(hmm);                  break;
      case P7_FS_CONFIG:    Plan7FSConfig(hmm, swentry, swexit); break;
      default:              Die("bogus configuration choice");
      }
      printf("done.\n");

      /* Save new HMM to disk: open a file for appending or writing.
       */
      printf("%-40s ... ", "Saving model to file");
      fflush(stdout);
      if (do_binary) WriteBinHMM(hmmfp, hmm);
      else           WriteAscHMM(hmmfp, hmm);
      printf("done.\n");

				/* the annotated alignment may be resaved */
      if (alignfp != NULL) 
	{
	  MSA    *new_msa;
	  SQINFO *sqinfo; 

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

      /* Verbose output; show scores for each sequence
       */
      if (verbose)
	print_all_scores(stdout, hmm, dsq, msa, tr);

      /* Clean up before moving on to next alignment
       */
      for (idx = 0; idx < msa->nseq; idx++) P7FreeTrace(tr[idx]);
      free(tr);
      FreePlan7(hmm);
      Free2DArray((void **) dsq, msa->nseq); 
      MSAFree(msa);
      fflush(hmmfp);
      if (cfp != NULL)     fflush(cfp);
      if (alignfp != NULL) fflush(alignfp);

      puts("//\n");
      nali++;
    }



  /* Clean up and exit
   */
  MSAFileClose(afp);
  fclose(hmmfp);
  if (cfp != NULL)     fclose(cfp);
  if (alignfp != NULL) fclose(alignfp);
  P7FreePrior(pri);
  SqdClean();
  return 0;
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
		 char **dsq, MSA *msa, struct p7trace_s **tr)
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
 * Args:     cfp    - open counts file 
 *           name   - name of alignment or HMM to associate with these vectors
 *           hmm    - counts-based HMM
 */
static void
save_countvectors(FILE *cfp, char *name, struct plan7_s *hmm)
{
  int k, x;
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
		       char             **dsq, 
		       float             *wgt,
		       int                nseq,
		       struct p7trace_s **tr,
		       float             *pernode,
		       float             *ret_avg)
{
  int    pos;                   /* position in seq */
  int    sym;
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
	sym = (int) dsq[idx][tr[idx]->pos[tpos]];
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
frag_trace_score(struct plan7_s *hmm, char *dsq, struct p7trace_s *tr, 
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
maximum_entropy(struct plan7_s *hmm, char **dsq, MSA *msa,
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

  FSet(wgt, msa->nseq, 1.0);
  position_average_score(hmm, dsq, wgt, msa->nseq, tr, pernode,&expscore);
  for (idx = 0; idx < msa->nseq; idx++) 
    sc[idx] = frag_trace_score(hmm, dsq[idx], tr[idx], pernode, expscore);
  relative_entropy = FSum(sc, msa->nseq) / (float) msa->nseq;
  for (idx = 0; idx < msa->nseq; idx++)
    grad[idx] = relative_entropy - sc[idx];

  
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
          FScale(new_wgt, msa->nseq, (float) msa->nseq);

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
  return;
}
