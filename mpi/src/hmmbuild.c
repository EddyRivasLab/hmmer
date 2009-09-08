/* main() for profile HMM construction from a multiple sequence alignment
 * 
 * SRE, Wed Jan  3 11:03:47 2007 [Janelia] [The Chemical Brothers]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_msaweight.h"
#include "esl_msacluster.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#define ALPHOPTS "--amino,--dna,--rna"                         /* Exclusive options for alphabet choice */
#define CONOPTS "--fast,--hand"                                /* Exclusive options for model construction                    */
#define EFFOPTS "--eent,--eclust,--eset,--enone"               /* Exclusive options for effective sequence number calculation */
#define WGTOPTS "--wgsc,--wblosum,--wpb,--wnone,--wgiven"      /* Exclusive options for relative weighting                    */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",                  1 },
  { "-n",        eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,    NULL, "name the HMM <s>",                                      1 },
  { "-o",        eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,    NULL, "direct summary output to file <f>, not stdout",         1 },
  { "-O",        eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,    NULL, "resave annotated, possibly modified MSA to file <f>",   1 },
/* Selecting the alphabet rather than autoguessing it */
  { "--amino",   eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is protein sequence data",              2 },
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is DNA sequence data",                  2 },
  { "--rna",     eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is RNA sequence data",                  2 },
/* Alternate model construction strategies */
  { "--fast",    eslARG_NONE,"default",NULL, NULL,    CONOPTS,    NULL,     NULL, "assign cols w/ >= symfrac residues as consensus",       3 },
  { "--hand",    eslARG_NONE,   FALSE, NULL, NULL,    CONOPTS,    NULL,     NULL, "manual construction (requires reference annotation)",   3 },
  { "--symfrac", eslARG_REAL,   "0.5", NULL, "0<=x<=1", NULL,   "--fast",   NULL, "sets sym fraction controlling --fast construction",     3 },
  { "--fragthresh",eslARG_REAL, "0.5", NULL, "0<=x<=1", NULL,     NULL,     NULL, "if L < x<L>, tag sequence as a fragment",               3 },
/* Alternate relative sequence weighting strategies */
  /* --wme not implemented in HMMER3 yet */
  { "--wpb",     eslARG_NONE,"default",NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff position-based weights",                      4 },
  { "--wgsc",    eslARG_NONE,   NULL,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Gerstein/Sonnhammer/Chothia tree weights",             4 },
  { "--wblosum", eslARG_NONE,   NULL,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff simple filter weights",                       4 },
  { "--wnone",   eslARG_NONE,   NULL,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "don't do any relative weighting; set all to 1",        4 },
  { "--wgiven",  eslARG_NONE,   NULL,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "use weights as given in MSA file",                     4 },
  { "--wid",     eslARG_REAL, "0.62",  NULL,"0<=x<=1",   NULL,"--wblosum",   NULL, "for --wblosum: set identity cutoff",                   4 },
/* Alternate effective sequence weighting strategies */
  { "--eent",    eslARG_NONE,"default",NULL, NULL,    EFFOPTS,    NULL,      NULL, "adjust eff seq # to achieve relative entropy target",  5 },
  { "--eclust",  eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "eff seq # is # of single linkage clusters",            5 },
  { "--enone",   eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "no effective seq # weighting: just use nseq",          5 },
  { "--eset",    eslARG_REAL,   NULL,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "set eff seq # for all models to <x>",                  5 },
  { "--ere",     eslARG_REAL,   NULL,  NULL,"x>0",       NULL, "--eent",     NULL, "for --eent: set minimum rel entropy/position to <x>",  5 },
  { "--esigma",  eslARG_REAL, "45.0",  NULL,"x>0",       NULL, "--eent",     NULL, "for --eent: set sigma param to <x>",                   5 },
  { "--eid",     eslARG_REAL, "0.62",  NULL,"0<=x<=1",   NULL,"--eclust",    NULL, "for --eclust: set fractional identity cutoff to <x>",  5 },
/* Control of E-value calibration */
  { "--EmL",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "length of sequences for MSV Gumbel mu fit",            6 },   
  { "--EmN",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "number of sequences for MSV Gumbel mu fit",            6 },   
  { "--EvL",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "length of sequences for Viterbi Gumbel mu fit",        6 },   
  { "--EvN",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "number of sequences for Viterbi Gumbel mu fit",        6 },   
  { "--EfL",     eslARG_INT,    "100", NULL,"n>0",       NULL,    NULL,      NULL, "length of sequences for Forward exp tail tau fit",     6 },   
  { "--EfN",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "number of sequences for Forward exp tail tau fit",     6 },   
  { "--Eft",     eslARG_REAL,  "0.04", NULL,"0<x<1",     NULL,    NULL,      NULL, "tail mass for Forward exponential tail tau fit",       6 },   
/* Other options */
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    "-n",  "run as an MPI parallel program",                       8 },  
#endif
  { "--informat", eslARG_STRING, NULL, NULL, NULL,      NULL,      NULL,    NULL, "assert input alifile is in format <s> (no autodetect)", 8 },
  { "--seed",     eslARG_INT,   "42", NULL, "n>=0",     NULL,      NULL,    NULL, "set RNG seed to <n> (if 0: one-time arbitrary seed)",   8 },
  { "--laplace", eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,    NULL, "use a Laplace +1 prior",                                8 },
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,    NULL, "arrest after start: for debugging MPI under gdb",       8 },  
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  FILE         *ofp;		/* output file (default is stdout) */

  char         *alifile;	/* name of the alignment file we're building HMMs from  */
  int           fmt;		/* format code for alifile */
  ESL_MSAFILE  *afp;            /* open alifile  */
  ESL_ALPHABET *abc;		/* digital alphabet */

  char         *hmmfile;        /* file to write HMM to                    */
  FILE         *hmmfp;          /* HMM output file handle                  */

  char         *postmsafile;	/* optional file to resave annotated, modified MSAs to  */
  FILE         *postmsafp;	/* open <postmsafile>, or NULL */

  P7_BG	       *bg;		/* null model                              */
  P7_PRIOR     *pri;		/* mixture Dirichlet prior for the HMM     */

  int           nali;		/* which # alignment this is in file (only valid in serial mode)   */
  int           nnamed;		/* number of alignments that had their own names */

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */
};


static char usage[]  = "[-options] <hmmfile output> <alignment file input>";
static char banner[] = "profile HMM construction from multiple sequence alignments";

static int  init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errmsg);
static int  init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errmsg);

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

static int output_header(const ESL_GETOPTS *go, const struct cfg_s *cfg);
static int output_result(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int msaidx, ESL_MSA *msa, P7_HMM *hmm, ESL_MSA *postmsa);
static int set_msa_name (const ESL_GETOPTS *go,       struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;	/* command line processing                 */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  struct cfg_s     cfg;

  /* Parse the command line
   */
  if ((go = esl_getopts_Create(options)) == NULL) esl_fatal("problem with options structure");
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
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\nOptions for selecting alphabet rather than guessing it:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      puts("\nAlternative model construction strategies:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\nAlternative relative sequence weighting strategies:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\nAlternate effective sequence weighting strategies:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\nControl of E-value calibration:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      puts("\nOther options:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 2) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      printf("\nTo see more help on other available options, do %s -h\n\n", argv[0]);
      exit(1);
    }


  /* Initialize what we can in the config structure (without knowing the alphabet yet) 
   */
  cfg.ofp         = NULL;	           /* opened in init_master_cfg() */
  cfg.alifile     = esl_opt_GetArg(go, 2);
  cfg.fmt         = eslMSAFILE_UNKNOWN;     /* autodetect alignment format by default. */ 
  cfg.afp         = NULL;	           /* created in init_master_cfg() */
  cfg.abc         = NULL;	           /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  cfg.hmmfile     = esl_opt_GetArg(go, 1); 
  cfg.hmmfp       = NULL;	           /* opened in init_master_cfg() */
  cfg.postmsafile = esl_opt_GetString(go, "-O"); /* NULL by default */
  cfg.postmsafp   = NULL;                  /* opened in init_master_cfg() */
  cfg.bg          = NULL;	           /* created in init_shared_cfg() */
  cfg.pri         = NULL;                  /* created in init_shared_cfg() */

  cfg.nali       = 0;		           /* this counter is incremented in masters */
  cfg.nnamed     = 0;		           /* 0 or 1 if a single MSA; == nali if multiple MSAs */
  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */
  cfg.do_stall   = esl_opt_GetBoolean(go, "--stall");

  if (esl_opt_IsOn(go, "--informat")) {
    cfg.fmt = esl_msa_EncodeFormat(esl_opt_GetString(go, "--informat"));
    if (cfg.fmt == eslMSAFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--informat"));
  }

  /* This is our stall point, if we need to wait until we get a
   * debugger attached to this process for debugging (especially
   * useful for MPI):
   */
  while (cfg.do_stall); 

  /* Start timing. */
  esl_stopwatch_Start(w);

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI
  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if (cfg.my_rank > 0)  mpi_worker(go, &cfg);
      else 		    mpi_master(go, &cfg);

      esl_stopwatch_Stop(w);
      esl_stopwatch_MPIReduce(w, 0, MPI_COMM_WORLD);
      MPI_Finalize();
    }
  else
#endif /*HAVE_MPI*/
    {
      serial_master(go, &cfg);
      esl_stopwatch_Stop(w);
    }

  if (cfg.my_rank == 0) {
    fputc('\n', cfg.ofp);
    esl_stopwatch_Display(cfg.ofp, w, "# CPU time: ");
  }


  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if (esl_opt_IsOn(go, "-o")) { fclose(cfg.ofp); }
    if (cfg.afp   != NULL) esl_msafile_Close(cfg.afp);
    if (cfg.abc   != NULL) esl_alphabet_Destroy(cfg.abc);
    if (cfg.hmmfp != NULL) fclose(cfg.hmmfp);
  }
  p7_bg_Destroy(cfg.bg);
  p7_prior_Destroy(cfg.pri);
  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  return 0;
}


/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Already set:
 *    cfg->hmmfile     - command line arg 1
 *    cfg->alifile     - command line arg 2
 *    cfg->postmsafile - option -O (default NULL)
 *    cfg->fmt         - format of alignment file
 * Sets: 
 *    cfg->afp       - open alignment file                
 *    cfg->abc       - digital alphabet
 *    cfg->hmmfp     - open HMM file
 *    cfg->postmsafp - open MSA resave file, or NULL
 *                   
 * Errors in the MPI master here are considered to be "recoverable",
 * in the sense that we'll try to delay output of the error message
 * until we've cleanly shut down the worker processes. Therefore
 * errors return (code, errmsg) by the ESL_FAIL mech.
 */
static int
init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errmsg)
{
  int status;

  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errmsg, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
  } else cfg->ofp = stdout;

  status = esl_msafile_Open(cfg->alifile, cfg->fmt, NULL, &(cfg->afp));
  if (status == eslENOTFOUND)    ESL_FAIL(status, errmsg, "Alignment file %s doesn't exist or is not readable\n", cfg->alifile);
  else if (status == eslEFORMAT) ESL_FAIL(status, errmsg, "Couldn't determine format of alignment %s\n", cfg->alifile);
  else if (status != eslOK)      ESL_FAIL(status, errmsg, "Alignment file open failed with error %d\n", status);

  if      (esl_opt_GetBoolean(go, "--amino"))   cfg->abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))     cfg->abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))     cfg->abc = esl_alphabet_Create(eslRNA);
  else {
    int type;
    status = esl_msafile_GuessAlphabet(cfg->afp, &type);
    if (status == eslEAMBIGUOUS)    ESL_FAIL(status, errmsg, "Failed to guess the bio alphabet used in %s.\nUse --dna, --rna, or --amino option to specify it.", cfg->alifile);
    else if (status == eslEFORMAT)  ESL_FAIL(status, errmsg, "Alignment file parse failed: %s\n", cfg->afp->errbuf);
    else if (status == eslENODATA)  ESL_FAIL(status, errmsg, "Alignment file %s is empty\n", cfg->alifile);
    else if (status != eslOK)       ESL_FAIL(status, errmsg, "Failed to read alignment file %s\n", cfg->alifile);
    cfg->abc = esl_alphabet_Create(type);
  }
  esl_msafile_SetDigital(cfg->afp, cfg->abc);

  if ((cfg->hmmfp = fopen(cfg->hmmfile, "w")) == NULL) ESL_FAIL(status, errmsg, "Failed to open HMM file %s for writing", cfg->hmmfile);

  if (cfg->postmsafile != NULL) {
    if ((cfg->postmsafp = fopen(cfg->postmsafile, "w")) == NULL) ESL_FAIL(status, errmsg, "Failed to MSA resave file %s for writing", cfg->postmsafile);
  } else cfg->postmsafp = NULL;

  output_header(go, cfg);

  /* with msa == NULL, output_result() prints the tabular results header, if needed */
  output_result(go, cfg, errmsg, 0, NULL, NULL, NULL);
  return eslOK;
}

/* init_shared_cfg() 
 * Shared initialization of cfg, after alphabet is known
 * Already set:
 *    cfg->abc
 * Sets:
 *    cfg->bg
 *    cfg->pri
 *    
 * Because this is called from an MPI worker, it cannot print; 
 * it must return error messages, not print them.
 */
static int
init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errmsg)
{
  if (esl_opt_GetBoolean(go, "--laplace"))
    {
      cfg->pri = p7_prior_CreateLaplace(cfg->abc);
    }
  else 
    {
      if      (cfg->abc->type == eslAMINO) cfg->pri = p7_prior_CreateAmino();
      else if (cfg->abc->type == eslDNA)   cfg->pri = p7_prior_CreateNucleic();  
      else if (cfg->abc->type == eslRNA)   cfg->pri = p7_prior_CreateNucleic();  
      else    ESL_FAIL(eslEINVAL, errmsg, "invalid alphabet type");
    }

  cfg->bg = p7_bg_Create(cfg->abc);
  
  if (cfg->pri == NULL) ESL_FAIL(eslEINVAL, errmsg, "alphabet initialization failed");
  if (cfg->bg  == NULL) ESL_FAIL(eslEINVAL, errmsg, "null model initialization failed");
  return eslOK;
}



/* serial_master()
 * The serial version of hmmbuild.
 * For each MSA, build an HMM and save it.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with p7_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  P7_BUILDER *bld         = NULL;
  ESL_MSA    *msa         = NULL;
  ESL_MSA    *postmsa     = NULL;
  ESL_MSA   **postmsa_ptr = (cfg->postmsafile != NULL) ? &postmsa : NULL;
  P7_HMM     *hmm         = NULL;
  char        errmsg[eslERRBUFSIZE];
  int         status;

  if ((status = init_master_cfg(go, cfg, errmsg)) != eslOK) p7_Fail(errmsg);
  if ((status = init_shared_cfg(go, cfg, errmsg)) != eslOK) p7_Fail(errmsg);
  
  if ((bld    = p7_builder_Create(go, cfg->abc))  == NULL)  p7_Fail("p7_builder_Create failed");

  cfg->nali = 0;
  while ((status = esl_msa_Read(cfg->afp, &msa)) != eslEOF)
    {
      if      (status == eslEFORMAT)  p7_Fail("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", cfg->afp->linenumber, cfg->afp->fname, cfg->afp->errbuf, cfg->afp->buf);
      else if (status != eslOK)       p7_Fail("Alignment file read unexpectedly failed with code %d\n", status);
      cfg->nali++;  

      if ((status = set_msa_name(go, cfg, errmsg, msa)) != eslOK) p7_Fail("%s\n", errmsg); /* cfg->nnamed gets incremented in this call */

                /*         bg   new-HMM trarr gm   om  */
      if ((status = p7_Builder(bld, msa, cfg->bg, &hmm, NULL, NULL, NULL, postmsa_ptr)) != eslOK) p7_Fail("build failed: %s", bld->errbuf);
      if ((status = output_result(go, cfg, errmsg, cfg->nali, msa, hmm, postmsa))       != eslOK) p7_Fail(errmsg);

      p7_hmm_Destroy(hmm);
      esl_msa_Destroy(msa);
      esl_msa_Destroy(postmsa);
    }

  p7_builder_Destroy(bld);
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of hmmbuild.
 * Follows standard pattern for a master/worker load-balanced MPI program (J1/78-79).
 * 
 * A master can only return if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * 
 * Recoverable errors include all worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers.
 * 
 * Unrecoverable errors are master-side errors that may affect MPI
 * communication, meaning we cannot count on being able to reach the
 * workers and shut them down. Unrecoverable errors result in immediate
 * p7_Fail()'s, which will cause MPI to shut down the worker processes
 * uncleanly.
 */
static void
mpi_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      xstatus       = eslOK;	/* changes from OK on recoverable error */
  int      status;
  int      have_work     = TRUE;	/* TRUE while alignments remain  */
  int      nproc_working = 0;	        /* number of worker processes working, up to nproc-1 */
  int      wi;          	        /* rank of next worker to get an alignment to work on */
  char    *buf           = NULL;	/* input/output buffer, for packed MPI messages */
  int      bn            = 0;
  ESL_MSA *msa           = NULL;
  P7_HMM  *hmm           = NULL;
  ESL_MSA **msalist      = NULL;
  ESL_MSA  *postmsa      = NULL;
  int      *msaidx       = NULL;
  char     errmsg[eslERRBUFSIZE];
  MPI_Status mpistatus; 
  int      n;
  int      pos;
  
  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */
  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errmsg)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { if ((status = init_shared_cfg(go, cfg, errmsg)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { bn = 4096; if ((buf = malloc(sizeof(char) * bn)) == NULL) { sprintf(errmsg, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((msalist = malloc(sizeof(ESL_MSA *) * cfg->nproc)) == NULL) { sprintf(errmsg, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((msaidx  = malloc(sizeof(int)       * cfg->nproc)) == NULL) { sprintf(errmsg, "allocation failed"); xstatus = eslEMEM; } }
  for (wi = 0; wi < cfg->nproc; wi++) { msalist[wi] = NULL; msaidx[wi] = 0; }
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) p7_Fail(errmsg);
  ESL_DPRINTF1(("MPI master is initialized\n"));

  /* Worker initialization:
   * Because we've already successfully initialized the master before we start
   * initializing the workers, we don't expect worker initialization to fail;
   * so we just receive a quick OK/error code reply from each worker to be sure,
   * and don't worry about an informative message. 
   */
  MPI_Bcast(&(cfg->abc->type), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (status != eslOK) p7_Fail("One or more MPI worker processes failed to initialize.");
  ESL_DPRINTF1(("%d workers are initialized\n", cfg->nproc-1));


  /* Main loop: combining load workers, send/receive, clear workers loops;
   * also, catch error states and die later, after clean shutdown of workers.
   * 
   * When a recoverable error occurs, have_work = FALSE, xstatus !=
   * eslOK, and errmsg is set to an informative message. No more
   * errmsg's can be received after the first one. We wait for all the
   * workers to clear their work units, then send them shutdown signals,
   * then finally print our errmsg and exit.
   * 
   * Unrecoverable errors just crash us out with p7_Fail().
   */
  wi = 1;
  while (have_work || nproc_working)
    {
      if (have_work) 
	{
	  if ((status = esl_msa_Read(cfg->afp, &msa)) == eslOK) 
	    {
	      cfg->nali++;  
	      ESL_DPRINTF1(("MPI master read MSA %s\n", msa->name == NULL? "" : msa->name));
	    }
	  else 
	    {
	      have_work = FALSE;
	      if      (status == eslEFORMAT)  { xstatus = eslEFORMAT; sprintf(errmsg, "Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", cfg->afp->linenumber, cfg->afp->fname, cfg->afp->errbuf, cfg->afp->buf); } 
	      else if (status != eslEOF)      { xstatus = status;     sprintf(errmsg, "Alignment file read unexpectedly failed with code %d\n", status); }
	      ESL_DPRINTF1(("MPI master has run out of MSAs (having read %d)\n", cfg->nali));
	    } 
	}

      if ((have_work && nproc_working == cfg->nproc-1) || (!have_work && nproc_working > 0))
	{
	  if (MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mpistatus) != 0) p7_Fail("mpi probe failed");
	  if (MPI_Get_count(&mpistatus, MPI_PACKED, &n)                != 0) p7_Fail("mpi get count failed");
	  wi = mpistatus.MPI_SOURCE;
	  ESL_DPRINTF1(("MPI master sees a result of %d bytes from worker %d\n", n, wi));

	  if (n > bn) {
	    if ((buf = realloc(buf, sizeof(char) * n)) == NULL) p7_Fail("reallocation failed");
	    bn = n; 
	  }
	  if (MPI_Recv(buf, bn, MPI_PACKED, wi, 0, MPI_COMM_WORLD, &mpistatus) != 0) p7_Fail("mpi recv failed");
	  ESL_DPRINTF1(("MPI master has received the buffer\n"));

	  /* If we're in a recoverable error state, we're only clearing worker results;
           * just receive them, don't unpack them or print them.
           * But if our xstatus is OK, go ahead and process the result buffer.
	   */
	  if (xstatus == eslOK)	
	    {
	      pos = 0;
	      if (MPI_Unpack(buf, bn, &pos, &xstatus, 1, MPI_INT, MPI_COMM_WORLD)     != 0)     p7_Fail("mpi unpack failed");
	      if (xstatus == eslOK) /* worker reported success. Get the HMM. */
		{
		  ESL_DPRINTF1(("MPI master sees that the result buffer contains an HMM\n"));
		  if (p7_hmm_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &(cfg->abc), &hmm) != eslOK) p7_Fail("HMM unpack failed");
		  ESL_DPRINTF1(("MPI master has unpacked the HMM\n"));

		  if (cfg->postmsafile != NULL) {
		    if (esl_msa_MPIUnpack(cfg->abc, buf, bn, &pos, MPI_COMM_WORLD, &postmsa) != eslOK) p7_Fail("postmsa unpack failed");
		  } 

		  if ((status = output_result(go, cfg, errmsg, msaidx[wi], msalist[wi], hmm, postmsa)) != eslOK) xstatus = status;

		  esl_msa_Destroy(postmsa); postmsa = NULL;
		  p7_hmm_Destroy(hmm);      hmm     = NULL;
		}
	      else	/* worker reported an error. Get the errmsg. */
		{
		  if (MPI_Unpack(buf, bn, &pos, errmsg, eslERRBUFSIZE, MPI_CHAR, MPI_COMM_WORLD) != 0) p7_Fail("mpi unpack of errmsg failed");
		  ESL_DPRINTF1(("MPI master sees that the result buffer contains an error message\n"));
		}
	    }
	  esl_msa_Destroy(msalist[wi]);
	  msalist[wi] = NULL;
	  msaidx[wi]  = 0;
	  nproc_working--;
	}

      if (have_work)
	{   
	  ESL_DPRINTF1(("MPI master is sending MSA %s to worker %d\n", msa->name == NULL ? "":msa->name, wi));
	  if (esl_msa_MPISend(msa, wi, 0, MPI_COMM_WORLD, &buf, &bn) != eslOK) p7_Fail("MPI msa send failed");
	  msalist[wi] = msa;
	  msaidx[wi]  = cfg->nali; /* 1..N for N alignments in the MSA database */
	  msa = NULL;
	  wi++;
	  nproc_working++;
	}
    }
  
  /* On success or recoverable errors:
   * Shut down workers cleanly. 
   */
  ESL_DPRINTF1(("MPI master is done. Shutting down all the workers cleanly\n"));
  for (wi = 1; wi < cfg->nproc; wi++) 
    if (esl_msa_MPISend(NULL, wi, 0, MPI_COMM_WORLD, &buf, &bn) != eslOK) p7_Fail("MPI msa send failed");
  free(buf);

  if (xstatus != eslOK) p7_Fail(errmsg);
  else                  return;
}


static void
mpi_worker(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int           xstatus = eslOK;
  int           status;
  int           type;
  P7_BUILDER   *bld         = NULL;
  ESL_MSA      *msa         = NULL;
  ESL_MSA      *postmsa     = NULL;
  ESL_MSA     **postmsa_ptr = (cfg->postmsafile != NULL) ? &postmsa : NULL;
  P7_HMM       *hmm         = NULL;
  char         *wbuf        = NULL;	/* packed send/recv buffer  */
  void         *tmp;			/* for reallocation of wbuf */
  int           wn          = 0;	/* allocation size for wbuf */
  int           sz, n;		        /* size of a packed message */
  int           pos;
  char          errmsg[eslERRBUFSIZE];

  /* After master initialization: master broadcasts its status.
   */
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return; /* master saw an error code; workers do an immediate normal shutdown. */
  ESL_DPRINTF2(("worker %d: sees that master has initialized\n", cfg->my_rank));
  
  /* Master now broadcasts worker initialization information (alphabet type) 
   * Workers returns their status post-initialization.
   * Initial allocation of wbuf must be large enough to guarantee that
   * we can pack an error result into it, because after initialization,
   * errors will be returned as packed (code, errmsg) messages.
   */
  MPI_Bcast(&type, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus == eslOK) { if ((cfg->abc = esl_alphabet_Create(type))      == NULL)    xstatus = eslEMEM; }
  if (xstatus == eslOK) { if ((status = init_shared_cfg(go, cfg, errmsg)) != eslOK)   xstatus = status;  }
  if (xstatus == eslOK) { wn = 4096;  if ((wbuf = malloc(wn * sizeof(char))) == NULL) xstatus = eslEMEM; }
  if (xstatus == eslOK) { if ((bld = p7_builder_Create(go, cfg->abc))     == NULL)    xstatus = eslEMEM; }
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); /* everyone sends xstatus back to master */
  if (xstatus != eslOK) {
    if (wbuf != NULL) free(wbuf);
    if (bld  != NULL) p7_builder_Destroy(bld);
    return; /* shutdown; we passed the error back for the master to deal with. */
  }

  ESL_DPRINTF2(("worker %d: initialized\n", cfg->my_rank));

                      /* source = 0 (master); tag = 0 */
  while (esl_msa_MPIRecv(0, 0, MPI_COMM_WORLD, cfg->abc, &wbuf, &wn, &msa) == eslOK) 
    {
      /* Build the HMM */
      ESL_DPRINTF2(("worker %d: has received MSA %s (%d columns, %d seqs)\n", cfg->my_rank, msa->name, msa->alen, msa->nseq));
      if ((status = p7_Builder(bld, msa, cfg->bg, &hmm, NULL, NULL, NULL, postmsa_ptr)) != eslOK) { strcpy(errmsg, bld->errbuf); goto ERROR; }
      ESL_DPRINTF2(("worker %d: has produced an HMM %s\n", cfg->my_rank, hmm->name));

      /* Calculate upper bound on size of sending status, HMM, and optional postmsa; make sure wbuf can hold it. */
      n = 0;
      if (MPI_Pack_size(1,    MPI_INT, MPI_COMM_WORLD, &sz) != 0)     goto ERROR;   n += sz;
      if (p7_hmm_MPIPackSize( hmm,     MPI_COMM_WORLD, &sz) != eslOK) goto ERROR;   n += sz;
      if (esl_msa_MPIPackSize(postmsa, MPI_COMM_WORLD, &sz) != eslOK) goto ERROR;   n += sz;
      if (n > wn) { ESL_RALLOC(wbuf, tmp, sizeof(char) * n); wn = n; }
      ESL_DPRINTF2(("worker %d: has calculated that HMM will pack into %d bytes\n", cfg->my_rank, n));

      /* Send status, HMM, and optional postmsa back to the master */
      pos = 0;
      if (MPI_Pack       (&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0)     goto ERROR;
      if (p7_hmm_MPIPack (hmm,                 wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK) goto ERROR;
      if (esl_msa_MPIPack(postmsa,             wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK) goto ERROR;
      MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
      ESL_DPRINTF2(("worker %d: has sent HMM to master in message of %d bytes\n", cfg->my_rank, pos));

      esl_msa_Destroy(msa);     msa     = NULL;
      esl_msa_Destroy(postmsa); postmsa = NULL;
      p7_hmm_Destroy(hmm);      hmm     = NULL;
    }

  if (wbuf != NULL) free(wbuf);
  p7_builder_Destroy(bld);
  return;

 ERROR:
  ESL_DPRINTF2(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errmsg));
  pos = 0;
  MPI_Pack(&status, 1,                MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errmsg,  eslERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
  if (wbuf != NULL) free(wbuf);
  if (msa  != NULL) esl_msa_Destroy(msa);
  if (hmm  != NULL) p7_hmm_Destroy(hmm);
  if (bld  != NULL) p7_builder_Destroy(bld);
  return;
}
#endif /*HAVE_MPI*/



static int
output_header(const ESL_GETOPTS *go, const struct cfg_s *cfg)
{
  if (cfg->my_rank > 0)  return eslOK;
  /*  if (! cfg->be_verbose) return eslOK; */
  p7_banner(cfg->ofp, go->argv[0], banner);
  
  fprintf(cfg->ofp, "# input alignment file:             %s\n", cfg->alifile);
  fprintf(cfg->ofp, "# output HMM file:                  %s\n", cfg->hmmfile);

  if (esl_opt_IsUsed(go, "-n"))          fprintf(cfg->ofp, "# name (the single) HMM:            %s\n",   esl_opt_GetString(go, "-n"));
  if (esl_opt_IsUsed(go, "-o"))          fprintf(cfg->ofp, "# output directed to file:          %s\n",   esl_opt_GetString(go, "-o"));
  if (esl_opt_IsUsed(go, "--amino"))     fprintf(cfg->ofp, "# input alignment is asserted as:   protein\n");
  if (esl_opt_IsUsed(go, "--dna"))       fprintf(cfg->ofp, "# input alignment is asserted as:   DNA\n");
  if (esl_opt_IsUsed(go, "--rna"))       fprintf(cfg->ofp, "# input alignment is asserted as:   RNA\n");
  if (esl_opt_IsUsed(go, "--fast"))      fprintf(cfg->ofp, "# model architecture construction:  fast/heuristic\n");
  if (esl_opt_IsUsed(go, "--hand"))      fprintf(cfg->ofp, "# model architecture construction:  hand-specified by RF annotation\n");
  if (esl_opt_IsUsed(go, "--symfrac"))   fprintf(cfg->ofp, "# sym fraction for model structure: %.3f\n", esl_opt_GetReal(go, "--symfrac"));
  if (esl_opt_IsUsed(go, "--fragthresh"))fprintf(cfg->ofp, "# seq called fragment if < xL    :  %.3f\n", esl_opt_GetReal(go, "--fragthresh"));
  if (esl_opt_IsUsed(go, "--wpb"))       fprintf(cfg->ofp, "# relative weighting scheme:        Henikoff PB\n");
  if (esl_opt_IsUsed(go, "--wgsc"))      fprintf(cfg->ofp, "# relative weighting scheme:        G/S/C\n");
  if (esl_opt_IsUsed(go, "--wblosum"))   fprintf(cfg->ofp, "# relative weighting scheme:        BLOSUM filter\n");
  if (esl_opt_IsUsed(go, "--wnone"))     fprintf(cfg->ofp, "# relative weighting scheme:        none\n");
  if (esl_opt_IsUsed(go, "--wid"))       fprintf(cfg->ofp, "# frac id cutoff for BLOSUM wgts:   %f\n",   esl_opt_GetReal(go, "--wid"));
  if (esl_opt_IsUsed(go, "--eent"))      fprintf(cfg->ofp, "# effective seq number scheme:      entropy weighting\n");
  if (esl_opt_IsUsed(go, "--eclust"))    fprintf(cfg->ofp, "# effective seq number scheme:      single linkage clusters\n");
  if (esl_opt_IsUsed(go, "--enone"))     fprintf(cfg->ofp, "# effective seq number scheme:      none\n");
  if (esl_opt_IsUsed(go, "--eset"))      fprintf(cfg->ofp, "# effective seq number:             set to %f\n", esl_opt_GetReal(go, "--eset"));
  if (esl_opt_IsUsed(go, "--ere") )      fprintf(cfg->ofp, "# minimum rel entropy target:       %f bits\n",   esl_opt_GetReal(go, "--ere"));
  if (esl_opt_IsUsed(go, "--esigma") )   fprintf(cfg->ofp, "# entropy target sigma parameter:   %f bits\n",   esl_opt_GetReal(go, "--esigma"));
  if (esl_opt_IsUsed(go, "--eid") )      fprintf(cfg->ofp, "# frac id cutoff for --eclust:      %f\n",        esl_opt_GetReal(go, "--eid"));
  if (esl_opt_IsUsed(go, "--EmL") )      fprintf(cfg->ofp, "# seq length for MSV Gumbel mu fit: %d\n",        esl_opt_GetInteger(go, "--EmL"));
  if (esl_opt_IsUsed(go, "--EmN") )      fprintf(cfg->ofp, "# seq number for MSV Gumbel mu fit: %d\n",        esl_opt_GetInteger(go, "--EmN"));
  if (esl_opt_IsUsed(go, "--EvL") )      fprintf(cfg->ofp, "# seq length for Vit Gumbel mu fit: %d\n",        esl_opt_GetInteger(go, "--EvL"));
  if (esl_opt_IsUsed(go, "--EvN") )      fprintf(cfg->ofp, "# seq number for Vit Gumbel mu fit: %d\n",        esl_opt_GetInteger(go, "--EvN"));
  if (esl_opt_IsUsed(go, "--EfL") )      fprintf(cfg->ofp, "# seq length for Fwd exp tau fit:   %d\n",        esl_opt_GetInteger(go, "--EfL"));
  if (esl_opt_IsUsed(go, "--EfN") )      fprintf(cfg->ofp, "# seq number for Fwd exp tau fit:   %d\n",        esl_opt_GetInteger(go, "--EfN"));
  if (esl_opt_IsUsed(go, "--Eft") )      fprintf(cfg->ofp, "# tail mass for Fwd exp tau fit:    %f\n",        esl_opt_GetReal(go, "--Eft"));
#ifdef HAVE_MPI
  if (esl_opt_IsUsed(go, "--mpi") )      fprintf(cfg->ofp, "# parallelization mode:             MPI\n");
#endif
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0) fprintf(cfg->ofp,"# random number seed:               one-time arbitrary\n");
    else                                       fprintf(cfg->ofp,"# random number seed set to:        %d\n", esl_opt_GetInteger(go, "--seed"));
  }
  if (esl_opt_IsUsed(go, "--laplace") )  fprintf(cfg->ofp, "# prior:                            Laplace +1\n");
  fprintf(cfg->ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");

  return eslOK;
}


static int
output_result(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int msaidx, ESL_MSA *msa, P7_HMM *hmm, ESL_MSA *postmsa)
{
  int status;

  /* Special case: output the tabular results header. 
   * Arranged this way to keep the two fprintf()'s close together in the code,
   * so we can keep the data and labels properly sync'ed.
   */
  if (msa == NULL)
    {
      fprintf(cfg->ofp, "#%4s %-20s %5s %5s %5s %8s %6s %s\n", " idx", "name",                 "nseq",  "alen",  "mlen",  "eff_nseq",  "re/pos",  "description");
      fprintf(cfg->ofp, "#%4s %-20s %5s %5s %5s %8s %6s %s\n", "----", "--------------------", "-----", "-----", "-----", "--------",  "------",  "-----------");
      return eslOK;
    }

  if ((status = p7_hmm_Validate(hmm, errbuf, 0.0001))       != eslOK) return status;
  if ((status = p7_hmmfile_WriteASCII(cfg->hmmfp, -1, hmm)) != eslOK) ESL_FAIL(status, errbuf, "HMM save failed");
  
	             /* #   name nseq alen M eff_nseq re/pos description*/
  fprintf(cfg->ofp, "%-5d %-20s %5d %5" PRId64 " %5d %8.2f %6.3f %s\n",
	  msaidx,
	  (msa->name != NULL) ? msa->name : "",
	  msa->nseq,
	  msa->alen,
	  hmm->M,
	  hmm->eff_nseq,
	  p7_MeanMatchRelativeEntropy(hmm, cfg->bg),
	  (msa->desc != NULL) ? msa->desc : "");
  
  if (cfg->postmsafp != NULL && postmsa != NULL) {
    esl_msa_Write(cfg->postmsafp, postmsa, eslMSAFILE_STOCKHOLM);
  }

  return eslOK;
}



/* set_msa_name() 
 * Make sure the alignment has a name; this name will
 * then be transferred to the model.
 * 
 * We can only do this for a single alignment in a file. For multi-MSA
 * files, each MSA is required to have a name already.
 *
 * Priority is:
 *      1. Use -n <name> if set, overriding any name the alignment might already have. 
 *      2. Use alignment's existing name, if non-NULL.
 *      3. Make a name, from alignment file name without path and without filename extension 
 *         (e.g. "/usr/foo/globins.slx" gets named "globins")
 * If none of these succeeds, return <eslEINVAL>.
 *         
 * If a multiple MSA database (e.g. Stockholm/Pfam), and we encounter
 * an MSA that doesn't already have a name, return <eslEINVAL> if nali > 1.
 * (We don't know we're in a multiple MSA database until we're on the second
 * alignment.)
 * 
 * If we're in MPI mode, we assume we're in a multiple MSA database,
 * even on the first alignment.
 * 
 * Because we can't tell whether we've got more than one
 * alignment 'til we're on the second one, these fatal errors
 * only happen after the first HMM has already been built.
 * Oh well.
 */
static int
set_msa_name(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, ESL_MSA *msa)
{
  char *name = NULL;
  int   status;

  if (cfg->do_mpi == FALSE && cfg->nali == 1) /* first (only?) HMM in file: */
    {
      if  (esl_opt_GetString(go, "-n") != NULL)
	{
	  if ((status = esl_msa_SetName(msa, esl_opt_GetString(go, "-n"))) != eslOK) return status;
	}
      else if (msa->name != NULL) 
	{
	  cfg->nnamed++;
	}
      else if (! cfg->afp->do_stdin)
	{
	  if ((status = esl_FileTail(cfg->afp->fname, TRUE, &name)) != eslOK) return status; /* TRUE=nosuffix */	  
	  if ((status = esl_msa_SetName(msa, name))                 != eslOK) return status;
	  free(name);
	}
      else ESL_FAIL(eslEINVAL, errbuf, "Failed to set model name: msa has no name, no msa filename, and no -n");
    }
  else 
    {
      if (esl_opt_GetString(go, "-n") != NULL) ESL_FAIL(eslEINVAL, errbuf, "Oops. Wait. You can't use -n with an alignment database.");
      else if (msa->name              != NULL) cfg->nnamed++;
      else                                     ESL_FAIL(eslEINVAL, errbuf, "Oops. Wait. I need name annotation on each alignment in a multi MSA file; failed on #%d", cfg->nali+1);

      /* special kind of failure: the *first* alignment didn't have a name, and we used the filename to
       * construct one; now that we see a second alignment, we realize this was a boo-boo*/
      if (cfg->nnamed != cfg->nali)            ESL_FAIL(eslEINVAL, errbuf, "Oops. Wait. I need name annotation on each alignment in a multi MSA file; first MSA didn't have one");
    }
  return eslOK;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
