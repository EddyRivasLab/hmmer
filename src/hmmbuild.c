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
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",              1 },
  { "-n",        eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,    NULL, "name the HMM <s>",                                  1 },
  { "-o",        eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,    NULL, "direct summary output to file <f>, not stdout",     1 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "use tabular output summary format, 1 line per HMM", 1 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "run as an MPI parallel program",                    1 },  
#endif
/* Selecting the alphabet rather than autoguessing it */
  { "--amino",   eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is protein sequence data",              2},
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is DNA sequence data",                  2},
  { "--rna",     eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is RNA sequence data",                  2},
/* Alternate model construction strategies */
  { "--fast",    eslARG_NONE,"default",NULL, NULL,    CONOPTS,    NULL,     NULL, "assign cols w/ >= symfrac residues as consensus",       3 },
  { "--hand",    eslARG_NONE,   FALSE, NULL, NULL,    CONOPTS,    NULL,     NULL, "manual construction (requires reference annotation)",   3 },
  { "--symfrac", eslARG_REAL,   "0.5", NULL, "0<=x<=1", NULL,   "--fast",   NULL, "sets sym fraction controlling --fast construction",     3 },
/* Alternate relative sequence weighting strategies */
  /* --wme not implemented in HMMER3 yet */
  { "--wgsc",    eslARG_NONE,"default",NULL, NULL,    WGTOPTS,    NULL,      NULL, "Gerstein/Sonnhammer/Chothia tree weights",         4},
  { "--wblosum", eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff simple filter weights",                   4},
  { "--wpb",     eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff position-based weights",                  4},
  { "--wnone",   eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "don't do any relative weighting; set all to 1",    4},
  { "--wgiven",  eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "use weights as given in MSA file",                 4},
  { "--pbswitch",eslARG_INT,  "1000",  NULL,"n>0",       NULL,    NULL,      NULL, "set failover to efficient PB wgts at > <n> seqs",  4},
  { "--wid",     eslARG_REAL, "0.62",  NULL,"0<=x<=1",   NULL,"--wblosum",   NULL, "for --wblosum: set identity cutoff",               4},
/* Alternate effective sequence weighting strategies */
  { "--eent",    eslARG_NONE,"default",NULL, NULL,    EFFOPTS,    NULL,      NULL, "adjust eff seq # to achieve relative entropy target", 5},
  { "--eclust",  eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "eff seq # is # of single linkage clusters",           5},
  { "--enone",   eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "no effective seq # weighting: just use nseq",         5},
  { "--eset",    eslARG_REAL,   NULL,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "set eff seq # for all models to <x>",                 5},
  { "--ere",     eslARG_REAL,   NULL,  NULL,"x>0",       NULL, "--eent",     NULL, "for --eent: set target relative entropy to <x>",      5},
  { "--eX",      eslARG_REAL,  "6.0",  NULL,"x>0",       NULL, "--eent",  "--ere", "for --eent: set minimum total rel ent param to <x>",  5},
  { "--eid",     eslARG_REAL, "0.62",  NULL,"0<=x<=1",   NULL,"--eclust",    NULL, "for --eclust: set fractional identity cutoff to <x>", 5},
/* Control of E-value calibration */
  { "--Es",      eslARG_INT,    NULL,  NULL,"n>0",       NULL,    NULL,      NULL, "set random number seed to <n>",                     6},
  { "--EvL",     eslARG_INT,    "100", NULL,"n>0",       NULL,    NULL,      NULL, "length of sequences for Viterbi Gumbel mu fit",     6},   
  { "--EvN",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "number of sequences for Viterbi Gumbel mu fit",     6},   
  { "--EfL",     eslARG_INT,    "100", NULL,"n>0",       NULL,    NULL,      NULL, "length of sequences for Forward exp tail mu fit",   6},   
  { "--EfN",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "number of sequences for Forward exp tail mu fit",   6},   
  { "--Eft",     eslARG_REAL,  "0.04", NULL,"0<x<1",     NULL,    NULL,      NULL, "tail mass for Forward exponential tail mu fit",     6},   
/* Other options */
  { "--laplace", eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,    NULL, "use a Laplace +1 prior",                            7},
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,    NULL, "arrest after start: for debugging MPI under gdb",   7},  
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

  P7_BG	       *bg;		/* null model                              */
  P7_DPRIOR    *pri;		/* mixture Dirichlet prior for the HMM     */

  int           be_verbose;	/* standard verbose output, as opposed to one-line-per-HMM summary */
  int           nali;		/* which # alignment this is in file (only valid in serial mode)   */

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */
};


static char usage[]  = "[-options] <hmmfile output> <alignment file input>";
static char banner[] = "profile HMM construction from a multiple sequence alignment";

static int  init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errmsg);
static int  init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errmsg);

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

static int process_workunit       (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***opt_tr);
static int output_result          (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int msaidx, ESL_MSA *msa, P7_HMM *hmm);

static int set_relative_weights   (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);
static int build_model            (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***opt_tr);
static int set_model_name         (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, P7_HMM *hmm);
static int set_effective_seqnumber(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, const ESL_MSA *msa, P7_HMM *hmm, const P7_BG *bg, const P7_DPRIOR *prior);
static int parameterize           (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, P7_HMM *hmm, const P7_DPRIOR *prior);
static int calibrate              (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, P7_HMM *hmm);

static double default_target_relent(const ESL_ALPHABET *abc, int M, double eX);

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
      puts("\n  where basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\n  options for selecting alphabet rather than guessing it:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      puts("\n  alternative model construction strategies:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\n  alternative relative sequence weighting strategies:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\n  alternate effective sequence weighting strategies:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\n  control of E-value calibration:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      puts("\n  other (rarely used) options:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 2) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      puts("\n  where basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      printf("\nTo see more help on other available options, do %s -h\n\n", argv[0]);
      exit(1);
    }


  /* Initialize what we can in the config structure (without knowing the alphabet yet) 
   */
  cfg.ofp        = NULL;	           /* opened in init_master_cfg() */
  cfg.alifile    = esl_opt_GetArg(go, 2);
  cfg.fmt        = eslMSAFILE_UNKNOWN;     /* autodetect alignment format by default. */ 
  cfg.afp        = NULL;	           /* created in init_master_cfg() */
  cfg.abc        = NULL;	           /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  cfg.hmmfile    = esl_opt_GetArg(go, 1); 
  cfg.hmmfp      = NULL;	           /* opened in init_master_cfg() */
  cfg.bg         = NULL;	           /* created in init_shared_cfg() */
  cfg.pri        = NULL;                   /* created in init_shared_cfg() */

  if (esl_opt_GetBoolean(go, "-1")) cfg.be_verbose = FALSE;        
  else                              cfg.be_verbose = TRUE;        
  cfg.nali       = 0;		           /* this counter is incremented in masters */
  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */
  cfg.do_stall   = esl_opt_GetBoolean(go, "--stall");


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
      cfg.be_verbose = FALSE;
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

  if (cfg.my_rank == 0) esl_stopwatch_Display(cfg.ofp, w, "# CPU time: ");


  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if (! esl_opt_IsDefault(go, "-o")) { fclose(cfg.ofp); }
    if (cfg.afp   != NULL) esl_msafile_Close(cfg.afp);
    if (cfg.abc   != NULL) esl_alphabet_Destroy(cfg.abc);
    if (cfg.hmmfp != NULL) fclose(cfg.hmmfp);
  }
  p7_bg_Destroy(cfg.bg);
  p7_dprior_Destroy(cfg.pri);
  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  return 0;
}


/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Already set:
 *    cfg->hmmfile - command line arg 1
 *    cfg->alifile - command line arg 2
 *    cfg->fmt     - format of alignment file
 * Sets: 
 *    cfg->afp     - open alignment file                
 *    cfg->abc     - digital alphabet
 *    cfg->hmmfp   - open HMM file
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

  /* with msa == NULL, output_result() prints the tabular results header, if needed */
  if (! cfg->be_verbose) output_result(go, cfg, errmsg, 0, NULL, NULL);
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
      cfg->pri = p7_dprior_CreateLaplace(cfg->abc);
    }
  else 
    {
      if      (cfg->abc->type == eslAMINO) cfg->pri = p7_dprior_CreateAmino();
      else if (cfg->abc->type == eslDNA)   cfg->pri = p7_dprior_CreateNucleic();  
      else if (cfg->abc->type == eslRNA)   cfg->pri = p7_dprior_CreateNucleic();  
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
  int      status;
  char     errmsg[eslERRBUFSIZE];
  ESL_MSA *msa = NULL;
  P7_HMM  *hmm = NULL;

  if ((status = init_master_cfg(go, cfg, errmsg)) != eslOK) p7_Fail(errmsg);
  if ((status = init_shared_cfg(go, cfg, errmsg)) != eslOK) p7_Fail(errmsg);

  cfg->nali = 0;
  while ((status = esl_msa_Read(cfg->afp, &msa)) != eslEOF)
    {
      if      (status == eslEFORMAT) p7_Fail("Alignment file parse error:\n%s\n", cfg->afp->errbuf);
      else if (status == eslEINVAL)  p7_Fail("Alignment file parse error:\n%s\n", cfg->afp->errbuf);
      else if (status != eslOK)      p7_Fail("Alignment file read failed with error code %d\n", status);
      cfg->nali++;  

      if (cfg->be_verbose) {
	if (msa->name != NULL) fprintf(cfg->ofp, "Alignment:           %s\n",  msa->name);
	else                   fprintf(cfg->ofp, "Alignment:           #%d\n", cfg->nali);
	fprintf                       (cfg->ofp, "Number of sequences: %d\n",  msa->nseq);
	fprintf                       (cfg->ofp, "Number of columns:   %" PRId64 "\n",  msa->alen);
	fputs("", cfg->ofp);
	fflush(stdout);
      }
      
      if ((status = process_workunit(go, cfg, errmsg,            msa, &hmm, NULL)) != eslOK) p7_Fail(errmsg);
      if ((status = output_result(   go, cfg, errmsg, cfg->nali, msa, hmm))        != eslOK) p7_Fail(errmsg);

      if (cfg->be_verbose) {
	fprintf(cfg->ofp, "Built a model of %d nodes.\n", hmm->M);
	fprintf(cfg->ofp, "Mean match relative entropy:  %.2f bits\n", p7_MeanMatchRelativeEntropy(hmm, cfg->bg));
	fprintf(cfg->ofp, "Mean match information:       %.2f bits\n", p7_MeanMatchInfo(hmm, cfg->bg));
      }

      p7_hmm_Destroy(hmm);
      esl_msa_Destroy(msa);
    }
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
	      if      (status == eslEFORMAT)  { xstatus = eslEFORMAT; snprintf(errmsg, eslERRBUFSIZE, "Alignment file parse error:\n%s\n", cfg->afp->errbuf); } 
	      else if (status == eslEINVAL)   { xstatus = eslEFORMAT; snprintf(errmsg, eslERRBUFSIZE, "Alignment file parse error:\n%s\n", cfg->afp->errbuf); } 
	      else if (status != eslEOF)      { xstatus = status;     snprintf(errmsg, eslERRBUFSIZE, "Alignment file read unexpectedly failed with code %d\n", status); }
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

		  if ((status = output_result(go, cfg, errmsg, msaidx[wi], msalist[wi], hmm))     != eslOK) xstatus = status;

		  p7_hmm_Destroy(hmm);
		  hmm = NULL;
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
  ESL_MSA      *msa  = NULL;
  P7_HMM       *hmm  = NULL;
  char         *wbuf = NULL;	/* packed send/recv buffer  */
  int           wn   = 0;	/* allocation size for wbuf */
  int           sz, n;		/* size of a packed message */
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
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); /* everyone sends xstatus back to master */
  if (xstatus != eslOK) {
    if (wbuf != NULL) free(wbuf);
    return; /* shutdown; we passed the error back for the master to deal with. */
  }
  ESL_DPRINTF2(("worker %d: initialized\n", cfg->my_rank));

                      /* source = 0 (master); tag = 0 */
  while (esl_msa_MPIRecv(0, 0, MPI_COMM_WORLD, cfg->abc, &wbuf, &wn, &msa) == eslOK) 
    {
      ESL_DPRINTF2(("worker %d: has received MSA %s (%d columns, %d seqs)\n", cfg->my_rank, msa->name, msa->alen, msa->nseq));
      if ((status =   process_workunit(go, cfg, errmsg, msa, &hmm, NULL)) != eslOK) goto ERROR;
      ESL_DPRINTF2(("worker %d: has produced an HMM %s\n", cfg->my_rank, hmm->name));

      n = 0;
      if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0)     goto ERROR;   n += sz;
      if (p7_hmm_MPIPackSize(hmm, MPI_COMM_WORLD, &sz)   != eslOK) goto ERROR;   n += sz;
      if (n > wn) {
	void *tmp;
	ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
	wn = n;
      }
      ESL_DPRINTF2(("worker %d: has calculated that HMM will pack into %d bytes\n", cfg->my_rank, n));

      pos = 0;
      if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD)   != 0)     goto ERROR;
      if (p7_hmm_MPIPack(hmm, wbuf, wn, &pos, MPI_COMM_WORLD)             != eslOK) goto ERROR;
      MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
      ESL_DPRINTF2(("worker %d: has sent HMM to master in message of %d bytes\n", cfg->my_rank, pos));

      esl_msa_Destroy(msa); msa = NULL;
      p7_hmm_Destroy(hmm);  hmm = NULL;
    }

  if (wbuf != NULL) free(wbuf);
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
  return;
}
#endif /*HAVE_MPI*/



/* A work unit consists of one multiple alignment, <msa>.
 * The job is to turn it into a new profile HMM, returned in <*ret_hmm>.
 */
static int
process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***opt_tr)
{
  P7_HMM *hmm = NULL;
  int status;

  if ((status =  set_relative_weights   (go, cfg, errbuf, msa))                         != eslOK) goto ERROR;
  if ((status =  build_model            (go, cfg, errbuf, msa, &hmm, opt_tr))           != eslOK) goto ERROR;
  if ((status =  set_model_name         (go, cfg, errbuf, msa, hmm))                    != eslOK) goto ERROR;
  if ((status =  set_effective_seqnumber(go, cfg, errbuf, msa, hmm, cfg->bg, cfg->pri)) != eslOK) goto ERROR;
  if ((status =  parameterize           (go, cfg, errbuf, hmm, cfg->pri))               != eslOK) goto ERROR;
  if ((status =  calibrate              (go, cfg, errbuf, hmm))                         != eslOK) goto ERROR;

  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  ESL_DPRINTF2(("worker %d: has caught an error in process_workunit\n", cfg->my_rank));
  p7_hmm_Destroy(hmm);
  p7_trace_DestroyArray(*opt_tr, msa->nseq);
  *ret_hmm = NULL;
  if (opt_tr  != NULL) *opt_tr = NULL;
  return status;
}


static int
output_result(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int msaidx, ESL_MSA *msa, P7_HMM *hmm)
{
  int status;

  /* Special case: output the tabular results header. 
   * Arranged this way to keep the two fprintf()'s close together in the code,
   * so we can keep the data and labels properly sync'ed.
   */
  if (msa == NULL && ! cfg->be_verbose) 
    {
      fprintf(cfg->ofp, "# %3s %-20s %5s %5s %5s\n", "idx", "name",                 "nseq",  "alen",  "M");
      fprintf(cfg->ofp, "#%4s %-20s %5s %5s %5s\n", "----", "--------------------", "-----", "-----", "-----");
      return eslOK;
    }

  if ((status = p7_hmm_Validate(hmm, errbuf, 0.0001))  != eslOK) return status;
  if ((status = p7_hmmfile_Write(cfg->hmmfp, hmm))     != eslOK) ESL_FAIL(status, errbuf, "HMM save failed");
  
  if (! cfg->be_verbose)	/* tabular output */
    {                    /* #   name nseq alen M */
      fprintf(cfg->ofp, "%-5d %-20s %5d %5" PRId64 " %5d\n",
	      msaidx,
	      (msa->name != NULL) ? msa->name : "",
	      msa->nseq,
	      msa->alen,
	      hmm->M);
    }
  return eslOK;
}


/* set_relative_weights():
 * Set msa->wgt vector, using user's choice of relative weighting algorithm.
 */
static int
set_relative_weights(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa)
{
  if (cfg->be_verbose) {
    fprintf(cfg->ofp, "%-40s ... ", "Relative sequence weighting");  
    fflush(cfg->ofp); 
  }

  if      (esl_opt_GetBoolean(go, "--wnone"))                  esl_vec_DSet(msa->wgt, msa->nseq, 1.);
  else if (esl_opt_GetBoolean(go, "--wgiven"))                 ;
  else if (msa->nseq >= esl_opt_GetInteger(go, "--pbswitch"))  esl_msaweight_PB(msa);
  else if (esl_opt_GetBoolean(go, "--wpb"))                    esl_msaweight_PB(msa);
  else if (esl_opt_GetBoolean(go, "--wgsc"))                   esl_msaweight_GSC(msa);
  else if (esl_opt_GetBoolean(go, "--wblosum"))                esl_msaweight_BLOSUM(msa, esl_opt_GetReal(go, "--wid"));

  if (cfg->be_verbose) fprintf(cfg->ofp, "done.\n");
  return eslOK;
}

/* build_model():
 * Given <msa>, choose HMM architecture, collect counts;
 * upon return, <*ret_hmm> is newly allocated and contains
 * relative-weighted observed counts.
 * Optionally, caller can request an array of inferred traces for
 * the <msa> too.
 */
static int
build_model(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***opt_tr)
{
  int status;

  if (cfg->be_verbose) {
    fprintf(cfg->ofp, "%-40s ... ", "Constructing model architecture"); 
    fflush(cfg->ofp);
  }

  if      (esl_opt_GetBoolean(go, "--fast")) 
    {
      status = p7_Fastmodelmaker(msa, esl_opt_GetReal(go, "--symfrac"), ret_hmm, opt_tr);
      if      (status == eslENORESULT) ESL_XFAIL(status, errbuf, "Alignment %s has no consensus columns w/ > %d%% residues - can't build a model.\n", msa->name != NULL ? msa->name : "", (int) (100 * esl_opt_GetReal(go, "--symfrac")));
      else if (status == eslEMEM)      ESL_XFAIL(status, errbuf, "Memory allocation failure in model construction.\n");
      else if (status != eslOK)        ESL_XFAIL(status, errbuf, "internal error in model construction.\n");      
    }
  else if (esl_opt_GetBoolean(go, "--hand")) 
    {
      status = p7_Handmodelmaker(msa, ret_hmm, opt_tr);
      if      (status == eslENORESULT) ESL_XFAIL(status, errbuf, "Alignment %s has no annotated consensus columns - can't build a model.\n", msa->name != NULL ? msa->name : "");
      else if (status == eslEFORMAT)   ESL_XFAIL(status, errbuf, "Alignment %s has no reference annotation line\n", msa->name != NULL ? msa->name : "");            
      else if (status == eslEMEM)      ESL_XFAIL(status, errbuf, "Memory allocation failure in model construction.\n");
      else if (status != eslOK)        ESL_XFAIL(status, errbuf, "internal error in model construction.\n");
    }

  if (cfg->be_verbose) fprintf(cfg->ofp, "done.\n");
  return eslOK;

 ERROR:
  if (cfg->be_verbose) fprintf(cfg->ofp, "FAILED.\n");
  return status;

}


/* set_model_name()
 * Give the model a name.
 * We deal with this differently depending on whether
 * we're in an alignment database or a single alignment.
 * 
 * If a single alignment, priority is:
 *      1. Use -n <name> if set.
 *      2. Use msa->name (avail in Stockholm/Pfam or SELEX formats only)
 *      3. If all else fails, use alignment file name without
 *         filename extension (e.g. "globins.slx" gets named "globins"
 *         
 * If a multiple MSA database (e.g. Stockholm/Pfam), 
 * only msa->name is applied. -n is not allowed.
 * if msa->name is unavailable, or -n was used,
 * a fatal error is thrown.
 * 
 * If we're in MPI mode, we assume we're in a multiple MSA database.
 * 
 * Because we can't tell whether we've got more than one
 * alignment 'til we're on the second one, these fatal errors
 * only happen after the first HMM has already been built.
 * Oh well.
 */
static int
set_model_name(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, P7_HMM *hmm)
{
  int status;

  if (cfg->be_verbose) {
    fprintf(cfg->ofp, "%-40s ... ", "Set model name");
    fflush(cfg->ofp);
  }

  if (cfg->do_mpi == FALSE && cfg->nali == 1)	/* first (only?) HMM in file:  */
    {
      if      (esl_opt_GetString(go, "-n") != NULL) p7_hmm_SetName(hmm, esl_opt_GetString(go, "-n"));
      else if (msa->name != NULL)                   p7_hmm_SetName(hmm, msa->name);
      else  
	{
	  char *name;
	  esl_FileTail(esl_opt_GetArg(go, 2), TRUE, &name); /* TRUE=nosuffix */
	  p7_hmm_SetName(hmm, name);
	  free(name);
	}
    }
  else				/* multi */
    {
      if (esl_opt_GetString(go, "-n") != NULL) ESL_XFAIL(eslEINVAL, errbuf, "Oops. Wait. You can't use -n with an alignment database.\n");
      else if (msa->name != NULL)              p7_hmm_SetName(hmm, msa->name);
      else                                     ESL_XFAIL(eslEINVAL, errbuf, "Oops. Wait. I need name annotation on each alignment.\n");
    }

  if (cfg->be_verbose) fprintf(cfg->ofp, "done. [%s]\n", hmm->name);
  return eslOK;

 ERROR:
  if (cfg->be_verbose) fprintf(cfg->ofp, "FAILED.\n");
  return status;
}


/* set_effective_seqnumber()
 * Incept:    SRE, Fri May 11 08:14:57 2007 [Janelia]
 *
 * <hmm> comes in with weighted observed counts. It goes out with
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
set_effective_seqnumber(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, const ESL_MSA *msa, P7_HMM *hmm, const P7_BG *bg, const P7_DPRIOR *prior)
{
  int    status = eslOK;
  double neff;

  if (cfg->be_verbose){
    fprintf(cfg->ofp, "%-40s ... ", "Set effective sequence number");
    fflush(cfg->ofp);
  }

  if      (esl_opt_GetBoolean(go, "--enone") == TRUE) 
    {
      neff = msa->nseq;
      if (cfg->be_verbose) fprintf(cfg->ofp, "done. [--enone: neff=nseq=%d]\n", msa->nseq);
    }
  else if (! esl_opt_IsDefault(go, "--eset"))
    {
      neff = esl_opt_GetReal(go, "--eset");
      if (cfg->be_verbose) fprintf(cfg->ofp, "done. [--eset: set to neff = %.2f]\n", neff);
    }
  else if (esl_opt_GetBoolean(go, "--eclust") == TRUE)
    {
      int nclust;

      status = esl_msacluster_SingleLinkage(msa, esl_opt_GetReal(go, "--eid"), NULL, NULL, &nclust);
      if      (status == eslEMEM) ESL_XFAIL(status, errbuf, "memory allocation failed");
      else if (status != eslOK)   ESL_XFAIL(status, errbuf, "single linkage clustering algorithm (at %d%% id) failed", (int)(100 * esl_opt_GetReal(go, "--eid")));

      neff = (double) nclust;
      if (cfg->be_verbose) fprintf(cfg->ofp, "done. [--eclust SLC at %.1f%%; neff = %.2f clusters]\n", 100. * esl_opt_GetReal(go, "--eid"), neff);
    }
  else if (esl_opt_GetBoolean(go, "--eent") == TRUE)
    {
      double etarget; 

      if (esl_opt_IsDefault(go, "--ere")) etarget = default_target_relent(hmm->abc, hmm->M, esl_opt_GetReal(go, "--eX"));
      else                                etarget = esl_opt_GetReal(go, "--ere");

      status = p7_EntropyWeight(hmm, bg, prior, etarget, &neff);
      if      (status == eslEMEM) ESL_XFAIL(status, errbuf, "memory allocation failed");
      else if (status != eslOK)   ESL_XFAIL(status, errbuf, "internal failure in entropy weighting algorithm");
    
      if (cfg->be_verbose) fprintf(cfg->ofp, "done. [etarget %.2f bits; neff %.2f]\n", etarget, neff);
    }
    
  hmm->eff_nseq = neff;
  p7_hmm_Scale(hmm, neff / (double) hmm->nseq);
  return eslOK;

 ERROR:
  if (cfg->be_verbose) fprintf(cfg->ofp, "FAILED.\n");
  return status;
}

/* parameterize()
 * Converts counts to probability parameters.
 */
static int
parameterize(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, P7_HMM *hmm, const P7_DPRIOR *prior)
{
  int status;

  if (cfg->be_verbose){ fprintf(cfg->ofp, "%-40s ... ", "Parameterizing"); fflush(cfg->ofp); }

  if ((status = p7_ParameterEstimation(hmm, prior)) != eslOK) ESL_XFAIL(status, errbuf, "parameter estimation failed");

  if (cfg->be_verbose) fprintf(cfg->ofp, "done.\n");
  return eslOK;

 ERROR:
  if (cfg->be_verbose) fprintf(cfg->ofp, "FAILED.\n");
  return status;
}

/* calibrate()
 * 
 * Sets the E value parameters of the model with two short simulations.
 */
static int
calibrate(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, P7_HMM *hmm)
{
  ESL_RANDOMNESS *r  = NULL;
  P7_PROFILE     *gm = NULL;
  int             vL = esl_opt_GetInteger(go, "--EvL");	/* length of random seqs for Viterbi mu sim */
  int             vN = esl_opt_GetInteger(go, "--EvN");	/* number of random seqs for Viterbi mu sim */
  int             fL = esl_opt_GetInteger(go, "--EfL");	/* length of random seqs for Forward mu sim */
  int             fN = esl_opt_GetInteger(go, "--EfN");	/* number of random seqs for Forward mu sim */
  double          ft = esl_opt_GetReal   (go, "--Eft");	/* tail mass for Forward mu sim             */
  double lambda, mu, tau;
  int    status;

  if (cfg->be_verbose) { fprintf(cfg->ofp, "%-40s ... ", "Calibrating");    fflush(cfg->ofp); }

  if (esl_opt_IsDefault(go, "--Es"))  r = esl_randomness_CreateTimeseeded();
  else                                r = esl_randomness_Create(esl_opt_GetInteger(go, "--Es"));
  if (r == NULL) ESL_XFAIL(eslEMEM, errbuf, "failed to create random number generator");

  if ((gm     = p7_profile_Create(hmm->M, cfg->abc))                  == NULL) ESL_XFAIL(eslEMEM, errbuf, "failed to allocate profile");
  if ((status = p7_ProfileConfig(hmm, cfg->bg, gm, vL, p7_LOCAL))    != eslOK) ESL_XFAIL(status,  errbuf, "failed to configure profile");
  if ((status = p7_Lambda(hmm, cfg->bg, &lambda))                    != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine lambda");
  if ((status = p7_Mu    (r, gm, cfg->bg, vL, vN, lambda, &mu))      != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine mu");
  if ((status = p7_Tau   (r, gm, cfg->bg, fL, fN, lambda, ft, &tau)) != eslOK) ESL_XFAIL(status,  errbuf, "failed to determine tau");

  hmm->evparam[p7_LAMBDA] = lambda;
  hmm->evparam[p7_MU]     = mu;
  hmm->evparam[p7_TAU]    = tau;
  hmm->flags             |= p7H_STATS;

  p7_profile_Destroy(gm);
  esl_randomness_Destroy(r);
  if (cfg->be_verbose) fprintf(cfg->ofp, "done.\n");
  return eslOK;

 ERROR:
  esl_randomness_Destroy(r);
  p7_profile_Destroy(gm);
  if (cfg->be_verbose) fprintf(cfg->ofp, "FAILED.\n");
  return status;
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


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
