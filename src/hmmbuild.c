/* hmmbuild.c
 * SRE, Mon Nov 18 12:41:29 1996
 *
 * main() for HMM construction from an alignment.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "structs.h"		/* data structures, macros, #define's   */
#include "config.h"		/* compile-time configuration constants */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "squid.h"		/* general sequence analysis library    */

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static char banner[] = "hmmb - build a hidden Markov model from an alignment";

static char usage[]  = "\
Usage: hmmb [-options] <hmmfile output> <alignment file>\n\
  Available options are:\n\
   -h        : help; print brief help on version and usage\n\
   -A        : append; append this HMM to <hmmfile>\n\
   -b        : binary; write the HMM in binary format (default is ASCII)\n\
   -F        : force; allow overwriting of <hmmfile> or other output files\n\
   -n <s>    : name; name this HMM <s> (default is <hmmfile>)\n\
   -o <file> : re-save annotated alignment to <file>\n\
\n\
  Alternative model construction strategies: (default: MAP)\n\
   -k        : Krogh/Haussler fast heuristic construction (see --gapmax)\n\
   -m        : manual construction (requires SELEX file, #=RF annotation)\n\
\n\
  Alternative sequence weighting strategies: (default: none)\n\
   --wgsc       : Gerstein/Sonnhammer/Chothia tree weights\n\
   --wblosum <x>: Henikoff^2 BLOSUM filter weights, at <x> frac id (e.g. 0.62)\n\
   --wvoronoi   : Sibbald/Argos Voronoi weights\n\
\n\
  Alternative parameter optimization strategies: (default: MAP)\n\
   -d            : maximum discrimination (MD)\n\
   -e            : maximum relative entropy (MRE)\n\
   --star <file> : Star model (experimental)\n\
\n\
  Expert customization of parameters and priors:\n\
   -r <file> : read null (random sequence) model from <file>\n\
   -p <file> : read Dirichlet prior parameters from <file>\n\
   -P <file> : heuristic PAM-based prior, using BLAST PAM matrix in <file>\n\
\n\
  Forcing an alphabet: (normally autodetected)\n\
   --amino   : override autodetection, assert that seqs are protein\n\
   --nucleic : override autodetection, assert that seqs are DNA/RNA\n\
\n\
  Other expert options:\n\
   --archpri <x>: set architecture size prior to <x> {0.85} [0..1]\n\
   --cfile <file>: save count vectors to <file>\n\
   --gapmax <x> : max fraction of gaps in mat column {0.50} [0..1]\n\
   --pamwgt <x> : set weight on PAM-based prior to <x> {20.}[>=0]\n\
   --verbose    : print a lot of boring information\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-A", TRUE, sqdARG_NONE },
  { "-b", TRUE, sqdARG_NONE },
  { "-d", TRUE, sqdARG_NONE },
  { "-e", TRUE, sqdARG_NONE }, 
  { "-F", TRUE, sqdARG_NONE },
  { "-h", TRUE, sqdARG_NONE }, 
  { "-n", TRUE, sqdARG_STRING},  
  { "-k", TRUE, sqdARG_NONE },
  { "-m", TRUE, sqdARG_NONE },
  { "-o", TRUE, sqdARG_STRING},
  { "-p", TRUE, sqdARG_STRING},
  { "-r", TRUE, sqdARG_STRING},
  { "-P", TRUE, sqdARG_STRING}, 
  { "--amino",   FALSE, sqdARG_NONE  },
  { "--archpri", FALSE, sqdARG_FLOAT }, 
  { "--cfile",   FALSE, sqdARG_STRING},
  { "--gapmax",  FALSE, sqdARG_FLOAT },
  { "--nucleic", FALSE, sqdARG_NONE },
  { "--pamwgt",  FALSE, sqdARG_FLOAT },
  { "--star"  ,  FALSE, sqdARG_STRING },
  { "--verbose", FALSE, sqdARG_NONE  },
  { "--wgsc",    FALSE, sqdARG_NONE },
  { "--wblosum", FALSE, sqdARG_FLOAT },
  { "--wvoronoi",FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static void save_model(struct plan7_s *hmm, char *hmmfile, int do_append, int do_binary);
static void print_all_scores(FILE *fp, struct plan7_s *hmm, 
			     AINFO *ainfo, char **dsq, int nseq, 
			     struct p7trace_s **tr);
static void make_star_model(struct plan7_s *hmm, char *starfile, struct p7prior_s *pri); 
static void save_countvectors(char *cfile, struct plan7_s *hmm);


int
main(int argc, char **argv) 
{
  char            *seqfile;     /* seqfile to read alignment from          */ 
  int              format;	/* format of seqfile                       */
  char           **aseq;        /* aligned sequences. [0.nseq-1][0.alen-1] */
  AINFO            ainfo;	/* optional info attached to aseq's        */
  char           **dsq;         /* digitized unaligned aseq's              */ 
  struct plan7_s  *hmm;         /* constructed HMM; written to hmmfile     */
  struct p7prior_s *pri;        /* Dirichlet priors to use                 */
  struct p7trace_s **tr;        /* fake tracebacks for aseq's              */ 
  char            *hmmfile;     /* file to write HMM to                    */
  FILE            *fp;          /* OUTPUT file handle (misc.)              */
  char            *name;        /* name of the HMM                         */
  int              idx;		/* counter for sequences                   */
  float  randomseq[MAXABET];	/* null sequence model                     */
  float            p1;		/* null sequence model p1 transition       */

  char *optname;                /* name of option found by Getopt()      */
  char *optarg;                 /* argument found by Getopt()            */
  int   optind;                 /* index in argv[]                       */
  enum p7_construction c_strategy;/* construction strategy choice        */
  enum p7_param p_strategy;	/* parameterization strategy choice      */
  enum p7_weight {		/* weighting strategy */
    WGT_NONE, WGT_GSC, WGT_BLOSUM, WGT_VORONOI} w_strategy;

  float gapmax;			/* max frac gaps in mat col for -k       */
  int   overwrite_protect;	/* TRUE to prevent overwriting HMM file  */
  int   verbose;		/* TRUE to show a lot of output          */
  char *align_ofile;            /* name of output alignment file         */
  char *rndfile;		/* random sequence model file to read    */
  char *prifile;		/* Dirichlet prior file to read          */
  char *pamfile;		/* PAM matrix file for heuristic prior   */
  char *starfile;               /* Star matrix file (GJM format, experimental */
  char *cfile;			/* output file for count vectors         */
  float archpri;		/* "architecture" prior on model size    */
  float pamwgt;			/* weight on PAM for heuristic prior     */
  int   do_append;		/* TRUE to append to hmmfile             */
  int   do_binary;		/* TRUE to write in binary format        */
  float blosumlevel;		/* BLOSUM frac id filtering level, e.g. 0.62 */

#ifdef MEMDEBUG
  unsigned long histid1, histid2;
  size_t orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  c_strategy        = P7_MAP_CONSTRUCTION;
  p_strategy        = P7_MAP_PARAM;
  w_strategy        = WGT_NONE;
  gapmax            = 0.5;
  overwrite_protect = TRUE;
  verbose           = FALSE;
  align_ofile       = NULL;
  rndfile           = NULL;
  prifile           = NULL;
  pamfile           = NULL;
  starfile          = NULL; 
  cfile             = NULL;
  archpri           = 0.85;
  pamwgt            = 20.;
  Alphabet_type     = 0;	/* initially unknown */
  name              = NULL;
  do_append         = FALSE; 
  do_binary         = FALSE;
  
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-A") == 0) do_append         = TRUE; 
    else if (strcmp(optname, "-b") == 0) do_binary         = TRUE;
    else if (strcmp(optname, "-d") == 0) p_strategy        = P7_MD_PARAM;
    else if (strcmp(optname, "-e") == 0) p_strategy        = P7_MRE_PARAM;  
    else if (strcmp(optname, "-F") == 0) overwrite_protect = FALSE;
    else if (strcmp(optname, "-k") == 0) c_strategy        = P7_FAST_CONSTRUCTION;
    else if (strcmp(optname, "-m") == 0) c_strategy        = P7_HAND_CONSTRUCTION;
    else if (strcmp(optname, "-m") == 0) name              = Strdup(optarg); 
    else if (strcmp(optname, "-o") == 0) align_ofile       = optarg;
    else if (strcmp(optname, "-p") == 0) prifile           = optarg;
    else if (strcmp(optname, "-r") == 0) rndfile           = optarg;
    else if (strcmp(optname, "-P") == 0) pamfile           = optarg;
    else if (strcmp(optname, "--amino")   == 0) SetAlphabet(hmmAMINO);
    else if (strcmp(optname, "--archpri") == 0) archpri       = atof(optarg);
    else if (strcmp(optname, "--cfile")   == 0) cfile         = optarg;
    else if (strcmp(optname, "--gapmax")  == 0) gapmax        = atof(optarg);
    else if (strcmp(optname, "--nucleic") == 0) SetAlphabet(hmmNUCLEIC);
    else if (strcmp(optname, "--pamwgt")  == 0) pamwgt        = atof(optarg);
    else if (strcmp(optname, "--star")    == 0) starfile      = optarg; 
    else if (strcmp(optname, "--verbose") == 0) verbose       = TRUE;
    else if (strcmp(optname, "--wgsc")    == 0) w_strategy    = WGT_GSC;
    else if (strcmp(optname, "--wblosum") == 0) {w_strategy = WGT_BLOSUM; blosumlevel = atof(optarg);}
    else if (strcmp(optname, "--wvoronoi")== 0) w_strategy    = WGT_VORONOI;
    else if (strcmp(optname, "-h") == 0) {
      Banner(stdout, banner);
      puts(usage);
      exit(0);
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

  /*********************************************** 
   * Get sequence data
   ***********************************************/

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
				/* Set up the alphabet globals */
  if (Alphabet_type == 0) DetermineAlphabet(aseq, ainfo.nseq);

				/* Set up Dirichlet priors */
  if (prifile == NULL)  pri = P7DefaultPrior();
  else                  pri = P7ReadPrior(prifile);

  if (pamfile != NULL)  PAMPrior(pamfile, pri, pamwgt);

				/* Set up the null/random seq model */
  if (rndfile == NULL)  P7DefaultNullModel(randomseq, &p1);
  else                  P7ReadNullModel(rndfile, randomseq, &p1);

				/* Prepare sequences for internal use */
  DigitizeAlignment(aseq, &ainfo, &dsq);
  
  /*********************************************** 
   * Show the banner
   ***********************************************/

  Banner(stdout, banner);
  printf("Training alignment:                %s\n", seqfile);
  printf("Number of sequences:               %d\n", ainfo.nseq);

  printf("Model construction strategy:       ");
  if (c_strategy == P7_HAND_CONSTRUCTION)    puts("Manual, from #=RF annotation");
  else if (c_strategy==P7_FAST_CONSTRUCTION) printf("Fast/ad hoc (gapmax %.2f)\n", gapmax);
  else                                       printf("MAP (gapmax hint: %.2f)\n", gapmax);

  printf("Prior used:                        %s\n",
	 (prifile == NULL) ? "(default)" : prifile);

  printf("Prior strategy:                    ");
  if      (pri->strategy == PRI_DCHLET) puts("Dirichlet");
  else if (pri->strategy == PRI_PAM)    puts("PAM hack");
  if (pri->strategy == PRI_PAM)
    printf("PAM prior weight:                  %.1f\n", pamwgt);

  printf("Sequence weighting method:         ");
  if      (w_strategy == WGT_NONE)   puts("none");
  else if (w_strategy == WGT_GSC)    puts("GSC tree weights");
  else if (w_strategy == WGT_BLOSUM) printf("BLOSUM filter at %.2f id\n", blosumlevel);
  else if (w_strategy == WGT_VORONOI)puts("Sibbald/Argos Voronoi");

  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  
  /*********************************************** 
   * Build an HMM
   ***********************************************/

  /* Weight the sequences first, if asked.
   */
  if      (w_strategy == WGT_GSC)     GSCWeights(aseq, &ainfo);
  else if (w_strategy == WGT_BLOSUM)  BlosumWeights(aseq, &ainfo, blosumlevel);
  else if (w_strategy == WGT_VORONOI) VoronoiWeights(aseq, &ainfo); 

  /* Build a model architecture.
   * If we're not doing MD or MRE, that's all we need to do.
   * We get an allocated, counts-based HMM back.
   */
  if (c_strategy == P7_FAST_CONSTRUCTION)
    P7Fastmodelmaker(aseq, dsq, &ainfo, gapmax, &hmm, &tr);
  else if (c_strategy == P7_HAND_CONSTRUCTION)
    P7Handmodelmaker(aseq, dsq, &ainfo, &hmm, &tr);
  else
    P7Maxmodelmaker(aseq, dsq, &ainfo, gapmax, 
		    pri, randomseq, p1, archpri, &hmm, &tr);

  /* Save the count vectors if asked. Used primarily for
   * making the data files for training priors.
   */
  if (cfile != NULL) save_countvectors(cfile, hmm); 

  /* Record the null model in the HMM;
   * add prior contributions in pseudocounts and renormalize.
   */
  Plan7SetNullModel(hmm, randomseq, p1);
  if (starfile != NULL) make_star_model(hmm, starfile, pri);
  else                  P7PriorifyHMM(hmm, pri);

  /* Give the model a name; by default, the name of the alignment file
   * without any filename extension (i.e. "globins.slx" becomes "globins"
   */
  if (name == NULL) name = FileTail(seqfile, TRUE);
  Plan7SetName(hmm, name);
  free(name);
  Plan7SetComline(hmm, argc, argv);
  Plan7SetCtime(hmm);
  hmm->nseq = ainfo.nseq;
   
  /* Save new HMM to disk: open a file for appending or writing.
   */
  Plan7LSConfig(hmm);
  save_model(hmm, hmmfile, do_append, do_binary);

				/* the annotated alignment may be resaved */
  if (align_ofile != NULL) {
    char  **new_aseq;
    AINFO   new_ainfo;

    P7Traces2Alignment(dsq, ainfo.sqinfo, ainfo.wgt, ainfo.nseq, hmm->M, tr, FALSE,
		       &new_aseq, &new_ainfo);
    if ((fp = fopen(align_ofile, "w")) == NULL) {
      Warn("Failed to open alignment resave file %s; using stdout instead");
      fp = stdout;
    }
    WriteSELEX(fp, new_aseq, &new_ainfo, 50);
    if (fp != stdout) fclose(fp);
    FreeAlignment(new_aseq, &new_ainfo);
  }
  
  /* Print information for the user
   */
  printf("Constructed a hidden Markov model (length %d)\n", hmm->M);
  PrintPlan7Stats(stdout, hmm, dsq, ainfo.nseq, tr); 

  printf("\nHMM written to %s\n", hmmfile);
  if (align_ofile != NULL)
    printf("Annotated alignment written to %s\n", align_ofile); 
  if (cfile != NULL)
    printf("Counts file written to %s\n", cfile);

  /* Verbose output; show scores for each sequence
   */
  if (verbose)
    print_all_scores(stdout, hmm, &ainfo, dsq, ainfo.nseq, tr);

  /* Clean up and exit
   */
#ifdef MEMDEBUG
  malloc_chain_check(1);
#endif

  for (idx = 0; idx < ainfo.nseq; idx++) P7FreeTrace(tr[idx]);
  free(tr);
  FreePlan7(hmm);
  P7FreePrior(pri);
  FreeAlignment(aseq, &ainfo);
  Free2DArray(dsq, ainfo.nseq); 
  SqdClean();

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) {
    printf("** I have %d bytes still allocated, pal\n",
	   current_size- orig_size);
    malloc_list(2, histid1, histid2);
  }
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  return 0;
}

/* Function: save_model()
 * 
 * Purpose:  Save the new model to a file.
 * 
 * Args:     hmm       - model to save
 *           hmmfile   - file to save to (if NULL, use stdout)      
 *           do_append - TRUE to append to file
 *           do_binary - TRUE to write a binary file
 *           
 * Return:   (void)
 */          
static void
save_model(struct plan7_s *hmm, char *hmmfile, int do_append, int do_binary)
{
  FILE    *fp;
  HMMFILE *hmmfp;

  if (hmmfile == NULL)
    fp = stdout;
  else if (do_append)
    {
				/* check that it looks like an HMM file */
#ifdef REMOVED			/* This code induces an unresolved Linux/SGI NFS bug! */
      if (FileExists(hmmfile)) 
	{ 
	  hmmfp = HMMFileOpen(hmmfile, NULL);
	  if (hmmfp == NULL) {
	    Warn("%s not an HMM file; I refuse to append to it; using stdout instead\n",
		 hmmfile);
	    fp = stdout;
	  } else {
	    HMMFileClose(hmmfp);
	  }
	}
#endif

      if ((fp = fopen(hmmfile, "a")) == NULL) {
	Warn("hey, where'd your HMM file go?");
	fp = stdout;
      } 
    } 
  else 
    {
      if ((fp = fopen(hmmfile, "w")) == NULL) {
	Warn("Failed to open HMM save file %s; using stdout instead");
	fp = stdout;
      }
    }

  if (do_binary) WriteBinHMM(fp, hmm);
  else           WriteAscHMM(fp, hmm);

  if (fp != stdout) fclose(fp);
  return;
}





/* Function: print_all_scores()
 * 
 * Purpose:  For each training sequence, print its score under
 *           the final model.
 *           
 * Args:     fp   - where to print the output (usu. stdout)
 *           hmm  - newly constructed HMM, with prob's.
 *           ainfo- info with aseq
 *           dsq  - digitized unaligned training sequences.
 *           nseq - number of training sequences
 *           tr   - array of tracebacks
 *           
 * Return:   (void)                         
 */
static void
print_all_scores(FILE *fp, struct plan7_s *hmm,
		 AINFO *ainfo, char **dsq, int nseq, struct p7trace_s **tr)
{
  int idx;			/* counter for sequences */

				/* make sure model scores are ready */
  Plan7Logoddsify(hmm);
				/* header */
  fputs("**\n", fp);
  fputs("Individual training sequence scores:\n", fp);
				/* score for each sequence */
  for (idx = 0; idx < nseq; idx++) 
    {
      fprintf(fp, "%7.2f %-12s %s\n", 
	      P7TraceScore(hmm, dsq[idx], tr[idx]),
	      ainfo->sqinfo[idx].name,
	      (ainfo->sqinfo[idx].flags & SQINFO_DESC) ? 
	      ainfo->sqinfo[idx].desc : "");
      P7PrintTrace(fp, tr[idx], hmm, dsq[idx]);
    }
  fputs("\n", fp);
}



/* Function: make_star_model()
 * 
 * Purpose:  Take an HMM in counts form, and the name of one
 *           of GJM's star matrix files, and produce a star model
 *           in probability form.
 */
static void
make_star_model(struct plan7_s *hmm, char *starfile, struct p7prior_s *pri) 
{
  FILE   *fp;
  float  **mx;
  float   *pq;
  int      nq;

  if ((fp = fopen(starfile, "r")) == NULL)
    Die("Failed to open GJM's star matrix file %s", starfile);
  ReadGJMMatrices(fp, &mx, &pq, &nq); 
  fclose(fp);

  MakeStarHMM(hmm, mx, pq, nq, pri);

  FMX2Free(mx);
  free(pq);
  return;
}
  

/* Function: save_countvectors()
 * 
 * Purpose:  Save emission/transition count vectors to a file.
 *           Used for gathering the data on which to train a
 *           prior (e.g. mixture Dirichlet, star, etc.)
 *           
 *           The format of the file is one vector per line:
 *           M <f> <f> ...: 20 match emission counts in order AC..WY.
 *           I <f> <f> ...: 20 insert emission counts in order AC..WY.
 *           T <f> <f> ...: 7 transition counts in order TMM, TMI, TMD, 
 *                            TIM, TII, TDM, TDD. (see structs.h)
 *           
 * Args:     cfile  - counts file to make
 *           hmm    - counts-based HMM
 */
static void
save_countvectors(char *cfile, struct plan7_s *hmm)
{
  FILE *fp;
  int k, x;

  if ((fp = fopen(cfile, "w")) == NULL)
    Die("failed to open count vector file %s for writing", cfile);

				/* match emission vectors */
  for (k = 1; k <= hmm->M; k++)
    {
      fputs("M ", fp);
      for (x = 0; x < Alphabet_size; x++)
	fprintf(fp, "%.2f ", hmm->mat[k][x]);
      fputs("\n", fp);
    }
				/* insert emission vectors */
  for (k = 1; k < hmm->M; k++)
    {
      fputs("I ", fp);
      for (x = 0; x < Alphabet_size; x++)
	fprintf(fp, "%.2f ", hmm->ins[k][x]);
      fputs("\n", fp);
    }
				/* transition vectors */
    for (k = 1; k < hmm->M; k++)
    {
      fputs("T ", fp);
      for (x = 0; x < 7; x++)
	fprintf(fp, "%.2f ", hmm->t[k][x]);
      fputs("\n", fp);
    }

  fclose(fp);
}
