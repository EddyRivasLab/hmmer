/* hmmscan: search sequence(s) against a profile HMM database
 * 
 * SRE, Mon Oct 20 08:28:05 2008 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "esl_mpi.h"
#endif 

#include "hmmer.h"

#define RNGOPTS "--Rdet,--Rseed,-Rarb"                         /* Exclusive options for controlling run-to-run variation      */

static ESL_OPTIONS options[] = {
  /* name           type          default  env  range toggles  reqs   incomp                         help                                                      docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL,  NULL,                          "show brief help on version and usage",                         1 },
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,   NULL,  NULL,  NULL,                          "direct output to file <f>, not stdout",                        1 },

  { "--hmmE",       eslARG_REAL,  "10.0", NULL, "x>0",  NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "E-value cutoff for reporting significant HMM hits",            2 },
  { "--hmmT",       eslARG_REAL,   FALSE, NULL, "x>0",  NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "bit score cutoff for reporting significant HMM hits",          2 },
  { "--domE",       eslARG_REAL,"1000.0", NULL, "x>0",  NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "E-value cutoff for reporting individual domains",              2 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, "x>0",  NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "bit score cutoff for reporting individual domains",            2 },
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL,  "--hmmE,--hmmT,--domE,--domT", "use GA gathering threshold bit score cutoffs in <hmmfile>",    2 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL,  "--hmmE,--hmmT,--domE,--domT", "use NC noise threshold bit score cutoffs in <hmmfile>",        2 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL,  "--hmmE,--hmmT,--domE,--domT", "use TC trusted threshold bit score cutoffs in <hmmfile>",      2 },
  { "--hmmZ",       eslARG_REAL,   FALSE, NULL, "x>0",  NULL,  NULL,  NULL,                          "set # of comparisons done, for E-value calculation",           2 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",  NULL,  NULL,  NULL,                          "set # of significant seqs, for domain E-value calculation",    2 },

  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, "--F1,--F2,--F3",   "Turn all heuristic filters off (less speed, more power)",      3 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,   NULL,  NULL, "--max",            "MSV threshold: promote hits w/ P <= F1",                       3 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,   NULL,  NULL, "--max",            "Vit threshold: promote hits w/ P <= F2",                       3 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,   NULL,  NULL, "--max",            "Fwd threshold: promote hits w/ P <= F3",                       3 },
  { "--nobias",     eslARG_NONE,    NULL, NULL, NULL,   NULL,  NULL, "--max",            "turn off composition bias filter",                             3 },
  { "--nonull2",    eslARG_NONE,    NULL, NULL, NULL,   NULL,  NULL,    NULL,            "turn off biased composition score corrections",                3 },
/* Control of run-to-run variation in RNG */
  { "--Rdet",       eslARG_NONE,"default",NULL, NULL,   RNGOPTS,  NULL,    NULL,         "reseed RNG to minimize run-to-run stochastic variation",       4 },
  { "--Rseed",       eslARG_INT,    NULL, NULL, NULL,   RNGOPTS,  NULL,    NULL,         "reseed RNG with fixed seed",                                   4 },
  { "--Rarb",       eslARG_NONE,    NULL, NULL, NULL,   RNGOPTS,  NULL,    NULL,         "seed RNG arbitrarily; allow run-to-run stochastic variation",  4 },
/* Other options */
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",  NULL,  NULL,  "--notextw",    "set max width of ASCII text output lines",                     5 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,      NULL,  NULL,  "--textw",      "unlimit ASCII text output line width",                         5 },
#ifdef HAVE_MPI
  //  { "--mpi",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL,  NULL,                          "run as an MPI parallel program",                               4 },
  //  { "--stall",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL,  NULL,                          "arrest after start: for debugging MPI under gdb",              4 },  
#endif 
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmm database> <query seqfile>";
static char banner[] = "search sequence(s) against a profile HMM database";

/* struct cfg_s: "Global" application configuration shared by all threads/MPI processes. */
struct cfg_s {
  /* Shared configuration in masters and workers                            */
  char            *hmmfile;     /* file to read HMM(s) from                 */
  char            *seqfile;     /* file to read sequence(s) from            */
  int              format;      /* format of seqfile                        */
  ESL_ALPHABET    *abc;         /* sequence alphabet                        */
  P7_BG           *bg;          /* null model                               */
  int              mode;	/* profile mode: e.g. p7_LOCAL              */
  P7_TOPHITS      *hitlist;	/* top-scoring sequence hits                */
  int              textw;	/* width of output ASCII text               */

  int     do_hmm_by_E;		/* TRUE to cut HMM reporting off by E       */
  double  hmmE;	                /* HMM E-value threshold                    */
  double  hmmT;	                /* HMM bit score threshold                  */
  int     do_dom_by_E;		/* TRUE to cut domain reporting off by E    */
  double  domE;	                /* domain E-value threshold                 */
  double  domT;	                /* domain bit score threshold               */
  int     model_cutoff_flag;    /* FALSE, p7H_GA, p7H_TC, or p7H_NC         */

  double  hmmZ;			/* # of HMMs searched for E-value purposes  */
  double  domZ;			/* # signific seqs for domain E-values      */
  int     hmmZ_is_fixed;	/* TRUE if hmmZ was set on cmd line         */
  int     domZ_is_fixed;	/* TRUE if domZ was set on cmd line         */

  int              do_max;	/* TRUE to run in slow/max mode             */
  double           F1;		/* MSV filter threshold                     */
  double           F2;		/* Viterbi filter threshold                 */
  double           F3;		/* uncorrected Forward filter threshold     */

  uint64_t         nseq;	/* number of sequences searched             */
  uint64_t         nmodels;	/* number of models searched                */
  uint64_t         nres;	/* # of residues searched                   */
  uint64_t         nnodes;	/* # of model nodes searched                */
  uint64_t         n_past_msv;	/* # models that pass the MSVFilter()       */
  uint64_t         n_past_vit;	/* # models that pass the ViterbiFilter()   */
  uint64_t         n_past_fwd;	/* # models that pass the ForwardFilter()   */

  /* Shared configuration for MPI */
  int              do_mpi;
  int              my_rank;
  int              nproc;
  int              do_stall;

  /* Workers only */
  ESL_RANDOMNESS  *r;		/* RNG is used for traceback sampling */
  P7_DOMAINDEF    *ddef;	/* a domain definition workbook       */

  /* Master only (i/o streams)                                          */
  P7_HMMFILE      *hfp;     /* open HMM file                            */  
  ESL_SQFILE      *sqfp;    /* open seqfile                             */
  FILE            *ofp;     /* output file for results (default stdout) */
};

static void process_commandline(int argc, char **argv, struct cfg_s *cfg, ESL_GETOPTS **ret_go);

static void init_shared_cfg(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  init_master_cfg(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  init_worker_cfg(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static void serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  output_header(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  set_pfam_bitscore_cutoffs (struct cfg_s *cfg, P7_OPROFILE *om);
static int  hmm_is_reportable         (struct cfg_s *cfg, float seq_score, float Pval);
static int  domain_is_reportable      (struct cfg_s *cfg, float dom_score, float Pval);
static void apply_thresholds_to_output(struct cfg_s *cfg, P7_TOPHITS *hitlist);

static int  output_per_model_hitlist (ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist);
static int  output_per_domain_hitlist(ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist);
static int  output_search_statistics (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_STOPWATCH *w);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go	   = NULL;                    /* command line processing  */
  ESL_STOPWATCH   *w       = esl_stopwatch_Create();  /* timing                   */
  struct cfg_s     cfg;

  /* Initializations */
  process_commandline(argc, argv, &cfg, &go);    
  init_shared_cfg(go, &cfg);                        
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */

  /* Stall here, if we need to wait for a debugger to attach to us (MPI) */
  while (cfg.do_stall); 

  /* Start timing */
  esl_stopwatch_Start(w);	                          
  
  /* Main body: hand off to serial version or MPI masters/workers, as appropriate */
#ifdef HAVE_MPI
  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      /* Initialize MPI, figure out who we are, and whether we're running
       * this show (proc 0) or working in it (procs >0).
       */
      cfg.do_mpi = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));
      if (cfg.my_rank == 0 && cfg.nproc < 2) p7_Fail("Need at least 2 MPI processes to run --mpi mode.");

#if 0
      if (cfg.my_rank > 0)   mpi_worker(go, &cfg);
      else                   mpi_master(go, &cfg);
#endif

      esl_stopwatch_Stop(w);
      esl_stopwatch_MPIReduce(w, 0, MPI_COMM_WORLD);
      MPI_Finalize();		/* both workers and masters reach this line */
    }
  else
#endif /*HAVE_MPI*/
    {      /* No MPI? Then we're just the serial master. */		
      serial_master(go, &cfg);
      esl_stopwatch_Stop(w);
    }      

  output_search_statistics(go, &cfg, w);            
  fprintf(cfg.ofp, "\n"); 

  /* Clean up and exit */
  if (cfg.my_rank == 0) 
    {	/* cleanup of master-only cfg */
      if (cfg.hfp  != NULL)               p7_hmmfile_Close(cfg.hfp);
      if (cfg.sqfp != NULL)               esl_sqfile_Close(cfg.sqfp);
      if (! esl_opt_IsDefault(go, "-o"))  fclose(cfg.ofp);
    }

  if (! cfg.do_mpi || cfg.my_rank > 0)
    {			/* cleanup of worker-only cfg */
      if (cfg.r    != NULL) esl_randomness_Destroy(cfg.r);
      if (cfg.ddef != NULL) p7_domaindef_Destroy(cfg.ddef);
    }
  p7_tophits_Destroy(cfg.hitlist);
  p7_bg_Destroy(cfg.bg);
  esl_alphabet_Destroy(cfg.abc);
  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  return eslOK;
}



/* process_commandline()
 * 
 * Processes the commandline, filling in fields in <cfg> and creating and returning
 * an <ESL_GETOPTS> options structure. The help page (hmmsearch -h) is formatted
 * here.
 */
static void
process_commandline(int argc, char **argv, struct cfg_s *cfg, ESL_GETOPTS **ret_go)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     p7_Die("Internal failure creating options object");
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere most common options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

      puts("\nOptions controlling significance thresholds for reporting:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      puts("\nOptions controlling acceleration heuristics:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 

      puts("\nOptions controlling run-to-run variation due to random number generation:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 120); 

      puts("\nOther expert options:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 2)      { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((cfg->hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <seqfile> argument on command line"); goto ERROR; }
  if ((cfg->seqfile = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <hmmfile> argument on command line"); goto ERROR; }

  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere most common options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);  
}
 

/* init_shared_cfg()
 * 
 * Initialize the parts of the <cfg> structure that are shared between
 * master and workers.
 */
static void
init_shared_cfg(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  /* cfg->hmmfile, cfg->seqfile already set in process_commandline()  */
  cfg->format   = eslSQFILE_UNKNOWN;    /* eventually, allow options to set this */
  cfg->abc      = NULL;		        /* to be initialized later */
  cfg->bg       = NULL;		        /* to be initialized later */
  cfg->mode     = p7_LOCAL;
  cfg->hitlist  = p7_tophits_Create(); /* master accumulates complete list; workers have partial lists. */

  if (esl_opt_GetBoolean(go, "--notextw")) cfg->textw = 0;
  else                                     cfg->textw = esl_opt_GetInteger(go, "--textw");

  /* Reporting thresholds */
  cfg->do_hmm_by_E       = TRUE;
  cfg->hmmE              = esl_opt_GetReal(go, "--hmmE");
  cfg->hmmT              = 0.0;
  cfg->do_dom_by_E       = TRUE;
  cfg->domE              = esl_opt_GetReal(go, "--domE");
  cfg->domT              = 0.0;
  cfg->model_cutoff_flag = 0;
  
  if (! esl_opt_IsDefault(go, "--hmmT")) 
    {
      cfg->hmmT        = esl_opt_GetReal(go, "--hmmT"); 
      cfg->do_hmm_by_E = FALSE;
    } 
  if (! esl_opt_IsDefault(go, "--domT")) 
    {
      cfg->domT        = esl_opt_GetReal(go, "--domT"); 
      cfg->do_dom_by_E = FALSE;
    }
  if (esl_opt_GetBoolean(go, "--cut_ga"))
    {
      cfg->hmmT        = cfg->domT        = 0.0;
      cfg->do_hmm_by_E = cfg->do_dom_by_E = FALSE;
      cfg->model_cutoff_flag = p7H_GA;
    }
  if (esl_opt_GetBoolean(go, "--cut_nc"))
    {
      cfg->hmmT        = cfg->domT        = 0.0;
      cfg->do_hmm_by_E = cfg->do_dom_by_E = FALSE;
      cfg->model_cutoff_flag = p7H_NC;
    }
  if (esl_opt_GetBoolean(go, "--cut_tc"))
    {
      cfg->hmmT        = cfg->domT        = 0.0;
      cfg->do_hmm_by_E = cfg->do_dom_by_E = FALSE;
      cfg->model_cutoff_flag = p7H_TC;
    }

  /* Search space sizes for E-value calculations */
  cfg->hmmZ_is_fixed = FALSE;
  cfg->domZ_is_fixed = FALSE;
  if (! esl_opt_IsDefault(go, "--hmmZ")) 
    {
      cfg->hmmZ_is_fixed = TRUE;
      cfg->hmmZ          = esl_opt_GetReal(go, "--hmmZ");
    }
  if (! esl_opt_IsDefault(go, "--domZ")) 
    {
      cfg->domZ_is_fixed = TRUE;
      cfg->domZ          = esl_opt_GetReal(go, "--domZ");
    }

  /* Heuristic filter thresholds */
  cfg->do_max = FALSE;
  cfg->F1     = esl_opt_GetReal(go, "--F1");
  cfg->F2     = esl_opt_GetReal(go, "--F2");
  cfg->F3     = esl_opt_GetReal(go, "--F3");
  if (esl_opt_GetBoolean(go, "--max")) 
    {
      cfg->do_max = TRUE;
      cfg->F1 = cfg->F2 = cfg->F3 = 1.0; 
    }


  /* Accounting as we collect results */
  cfg->nseq       = 0;
  cfg->nmodels    = 0;
  cfg->nres       = 0;
  cfg->nnodes     = 0;
  cfg->n_past_msv = 0;
  cfg->n_past_vit = 0;
  cfg->n_past_fwd = 0;

  /* These will be initialized later if --mpi is set: */
  cfg->do_mpi   = FALSE;
  cfg->my_rank  = 0;
  cfg->nproc    = 0;
  cfg->do_stall = FALSE;
  // cfg->do_stall = esl_opt_GetBoolean(go, "--stall");

  /* These are initialized later in the workers only: */
  cfg->r        = NULL;
  cfg->ddef     = NULL;

  /* These are initialized later in the master only: */
  cfg->hfp      = NULL;
  cfg->sqfp     = NULL;
  cfg->ofp      = NULL;
}


/* init_master_cfg()
 *
 * Responsible for opening the i/o streams that only a master has his
 * hands on.
 *
 * Error handling relies on these pointers being initialized to NULL
 * by the caller; this was done in init_shared_cfg().
 *                   
 * Since now we may be within MPI, errors can no longer be handled by
 * dumping a message and exiting. We have worker processes to shut
 * down, and possibly other cleanup to do, so we must try to delay
 * resolution of any error messages until after we attempt to clean
 * up. Therefore errors return (code, errmsg) by the ESL_FAIL mech.
 */
static int
init_master_cfg(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  char *filename;
  int   status;

  status = p7_hmmfile_Open(cfg->hmmfile, p7_HMMDBENV, &(cfg->hfp));
  if      (status == eslENOTFOUND) ESL_FAIL(status, errbuf, "Failed to open hmm file %s for reading.\n",     cfg->hmmfile);
  else if (status == eslEFORMAT)   ESL_FAIL(status, errbuf, "Unrecognized format, trying to open hmm file %s for reading.\n", cfg->hmmfile);
  else if (status != eslOK)        ESL_FAIL(status, errbuf, "Unexpected error %d in opening hmm file %s.\n", status, cfg->hmmfile);  
  
  status = esl_sqfile_Open(cfg->seqfile, cfg->format, NULL, &(cfg->sqfp));
  if      (status == eslENOTFOUND) ESL_FAIL(status, errbuf, "Failed to open sequence file %s for reading\n",      cfg->seqfile);
  else if (status == eslEFORMAT)   ESL_FAIL(status, errbuf, "Sequence file %s is empty or misformatted\n",        cfg->seqfile);
  else if (status == eslEINVAL)    ESL_FAIL(status, errbuf, "Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        ESL_FAIL(status, errbuf, "Unexpected error %d opening sequence file %s\n", status, cfg->seqfile);

  if (! cfg->hfp->is_pressed)      ESL_FAIL(eslEINVAL, errbuf, "Failed to open binary dbs for HMM file %s: use hmmpress first\n", cfg->hmmfile);

  filename = esl_opt_GetString(go, "-o");
  if (filename != NULL) {
    if ((cfg->ofp = fopen(filename, "w")) == NULL) 
      ESL_FAIL(eslESYS, errbuf, "Failed to open -o output file %s\n", filename);
  } else cfg->ofp = stdout;

  return eslOK;
}

/* init_worker_cfg()
 * 
 * As with the init_master_cfg, we may be within MPI now, so we
 * can't handle errors by dumping a message and exiting; 
 * we have to delay shutdown until we get MPI cleaned up. 
 * So, errors return via the ESL_FAIL mechanism.
 */
static int
init_worker_cfg(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  if ((cfg->r    = esl_randomness_CreateTimeseeded()) == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create random number generator");
  if ((cfg->ddef = p7_domaindef_Create(cfg->r))       == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create domain definition workbook");
  return eslOK;
}





static void
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  P7_OPROFILE     *om      = NULL;     /* optimized profiile                      */
  ESL_SQ          *sq      = NULL;     /* target sequence                         */
  P7_OMX          *fwd     = NULL;     /* DP matrix                               */
  P7_OMX          *bck     = NULL;     /* DP matrix                               */
  P7_OMX          *oxf     = NULL;     /* optimized DP matrix for Forward         */
  P7_OMX          *oxb     = NULL;     /* optimized DP matrix for Backward        */
  P7_HIT          *hit     = NULL;     /* ptr to the current hit output data      */
  float            usc, vfsc, fwdsc;   /* filter scores                           */
  float            filtersc;	       /* HMM null filter score                   */
  float            nullsc;             /* null model score                        */
  float            seqbias;  	
  float            seq_score;	       /* final corrected per-seq bit score       */
  float            sum_score;	       /* corrected reconstruction score for seq  */
  float            pre_score, pre2_score;
  double           P;		       /* P-value of a hit */
  int              Ld;		       /* # of residues in envelopes */
  char             errbuf[eslERRBUFSIZE];
  int              status, hstatus, sstatus;
  int              d;

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) esl_fatal(errbuf);
  if ((status = init_worker_cfg(go, cfg, errbuf)) != eslOK) esl_fatal(errbuf);

  sq = esl_sq_Create();	/* initially in text mode, until first HMM is read. */
 
  output_header(go, cfg);

  if (! cfg->hmmZ_is_fixed && cfg->hfp->is_pressed) {
    cfg->hmmZ          =  (double) cfg->hfp->ssi->nprimary;
  }

  /* We don't know the alphabet yet: cfg->abc == NULL until we read the first HMM. */

  while ( (sstatus = esl_sqio_Read(cfg->sqfp, sq)) == eslOK)
    {
      fwd = p7_omx_Create(200, sq->n, sq->n); /* initial alloc is for M=200, L; will grow as needed */
      bck = p7_omx_Create(200, sq->n, sq->n);	
      oxf = p7_omx_Create(200, 0,     sq->n); /* one-row parsing matrix, O(M+L) */
      oxb = p7_omx_Create(200, 0,     sq->n);     

      fprintf(cfg->ofp, "Query:       %s  [L=%ld]\n", sq->name, (long) sq->n);
      if (sq->acc[0]  != 0) fprintf(cfg->ofp, "Accession:   %s\n", sq->acc);
      if (sq->desc[0] != 0) fprintf(cfg->ofp, "Description: %s\n", sq->desc);
      cfg->nseq++;
      cfg->nres += sq->n;

      if (cfg->abc != NULL)	/* once we're on sequence #>2, abc is known, bg exists */
	{
	  p7_bg_SetLength(cfg->bg, sq->n);
	  p7_bg_NullOne  (cfg->bg, sq->dsq, sq->n, &nullsc);
	}
	  
      while ((hstatus = p7_oprofile_ReadMSV(cfg->hfp, &(cfg->abc), &om)) == eslOK) 
	{
	  /* One time only initializations after abc becomes known: */
	  if (cfg->bg == NULL) 	/* bg == NULL serves as a flag for the first-time init  */
	    {
	      cfg->bg = p7_bg_Create(cfg->abc);
	      if (esl_sq_Digitize(cfg->abc, sq) != eslOK) p7_Die("alphabet mismatch");
	      esl_sqfile_SetDigital(cfg->sqfp, cfg->abc);
	      p7_bg_SetLength(cfg->bg, sq->n);
	      p7_bg_NullOne  (cfg->bg, sq->dsq, sq->n, &nullsc);
	    }

	  cfg->nmodels++;
	  cfg->nnodes += om->M;
	  if (! cfg->hmmZ_is_fixed && ! cfg->hfp->is_pressed) cfg->hmmZ = cfg->nmodels;

	  /* Pfam bit score thresholds may be in use */
	  if (cfg->model_cutoff_flag && set_pfam_bitscore_cutoffs(cfg, om) != eslOK) p7_Fail("requested score cutoffs not available in model %s\n", om->name);

	  p7_omx_GrowTo(oxf, om->M, 0, sq->n); 
	  p7_oprofile_ReconfigMSVLength(om, sq->n);
	
	  /* First level filter: the MSV filter, multihit with <om> */
	  p7_MSVFilter(sq->dsq, sq->n, om, oxf, &usc);
	  seq_score = (usc - nullsc) / eslCONST_LOG2;
	  P = esl_gumbel_surv(seq_score,  om->evparam[p7_MU],  om->evparam[p7_LAMBDA]);
	  if (P > cfg->F1 && ! hmm_is_reportable(cfg, seq_score, P)) goto FINISH;

	  /* If it passes the MSV filter, read the rest of the profile */
	  p7_oprofile_ReadRest(cfg->hfp, om);
	  p7_oprofile_ReconfigRestLength(om, sq->n);

	  /* HMM filtering */
	  if (! esl_opt_GetBoolean(go, "--nobias") && !  esl_opt_GetBoolean(go, "--max"))
	    {
	      p7_bg_SetFilter(cfg->bg, om->M, om->compo); /* EXPERIMENTAL */
	      p7_bg_FilterScore(cfg->bg, sq->dsq, sq->n, &filtersc);
	      seq_score = (usc - filtersc) / eslCONST_LOG2;
	      P = esl_gumbel_surv(seq_score,  om->evparam[p7_MU],  om->evparam[p7_LAMBDA]);
	      if (P > cfg->F1 && ! hmm_is_reportable(cfg, seq_score, P)) goto FINISH;
	    }
	  else filtersc = nullsc;
	  cfg->n_past_msv++;


	  /* Second level filter: ViterbiFilter(), multihit with <om> */
	  p7_ViterbiFilter(sq->dsq, sq->n, om, oxf, &vfsc);  
	  seq_score = (vfsc-filtersc) / eslCONST_LOG2;
	  P  = esl_gumbel_surv(seq_score,  om->evparam[p7_MU], om->evparam[p7_LAMBDA]);
	  if (P > cfg->F2 && ! hmm_is_reportable(cfg, seq_score, P)) goto FINISH;
	  cfg->n_past_vit++;


	  /* Parse it with Forward and obtain its real Forward score. */
	  p7_ForwardParser(sq->dsq, sq->n, om, oxf, &fwdsc);
	  seq_score = (fwdsc-filtersc) / eslCONST_LOG2;
	  P = esl_exp_surv(seq_score,  om->evparam[p7_TAU], om->evparam[p7_LAMBDA]);
	  if (P > cfg->F3 && ! hmm_is_reportable(cfg, seq_score, P)) goto FINISH;
	  cfg->n_past_fwd++;

	  /* It's for real; do the Backwards parser pass and hand it
	   * off to domain definition logic.
	   */
	  p7_omx_GrowTo(oxb, om->M, 0, sq->n);
	  p7_BackwardParser(sq->dsq, sq->n, om, oxf, oxb, NULL);
	  status = p7_domaindef_ByPosteriorHeuristics(sq, om, oxf, oxb, fwd, bck, cfg->ddef);
	  if      (status == eslERANGE) { fprintf(stderr, "WARNING: posterior decoding failed on %s; skipping it!\n", sq->name);                       goto FINISH; }
	  else if (status != eslOK)     { fprintf(stderr, "WARNING: domain definition failed due to unknown problem on %s; skipping it!\n", sq->name); goto FINISH; }

	  /* What's the per-seq score, and is it significant enough to be reported? */
	  /* Figure out the sum of null2 corrections to be added to the null score */
	  if (! esl_opt_GetBoolean(go, "--nonull2"))
	    {
	      seqbias = esl_vec_FSum(cfg->ddef->n2sc, sq->n+1);
	      seqbias = p7_FLogsum(0.0, log(cfg->bg->omega) + seqbias);
	    }
	  else seqbias = 0.0f;
	  pre_score =  (fwdsc - nullsc) / eslCONST_LOG2; 
	  seq_score =  (fwdsc - (nullsc + seqbias)) / eslCONST_LOG2;

	  /* Calculate the "reconstruction score": estimated
	   * per-sequence score as sum of individual domains,
	   * discounting domains that aren't significant after they're
	   * null-corrected.
	   */
	  sum_score = 0.0f;
	  seqbias   = 0.0f;
	  Ld        = 0;
	  for (d = 0; d < cfg->ddef->ndom; d++) 
	    {
	      if (esl_opt_GetBoolean(go, "--nonull2")) cfg->ddef->dcl[d].domcorrection = 0.0f;
	      if (cfg->ddef->dcl[d].envsc - cfg->ddef->dcl[d].domcorrection > 0.0) 
		{
		  sum_score += cfg->ddef->dcl[d].envsc;
		  Ld        += cfg->ddef->dcl[d].jenv  - cfg->ddef->dcl[d].ienv + 1;
		  seqbias   += cfg->ddef->dcl[d].domcorrection;
		}
	    }
	  sum_score += (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3)); 
	  if (! esl_opt_GetBoolean(go, "--nonull2")) 
	    seqbias = p7_FLogsum(0.0, log(cfg->bg->omega) + seqbias);
	  pre2_score = (sum_score - nullsc) / eslCONST_LOG2;
	  sum_score  = (sum_score - (nullsc + seqbias)) / eslCONST_LOG2;

	  /* A special case: let sum_score override the seq_score when it's better, and it includes at least 1 domain */
	  if (Ld > 0 && sum_score > seq_score)
	    {
	      seq_score = sum_score;
	      pre_score = pre2_score;
	    }

	  /* Apply thresholding and determine whether to put this
	   * sequence into the hit list. E-value thresholding may
	   * only be a lower bound for now, so this list may be longer
	   * than eventually reported.
	   */
	  P =  esl_exp_surv (seq_score,  om->evparam[p7_TAU], om->evparam[p7_LAMBDA]);
	  if (hmm_is_reportable(cfg, seq_score, P))
	    {
	      /* Calculate and store the per-seq hit output information */
	      p7_tophits_CreateNextHit(cfg->hitlist, &hit);
	      if ((status  = esl_strdup(om->name, -1, &(hit->name)))  != eslOK) esl_fatal("allocation failure");
	      if ((status  = esl_strdup(om->acc,  -1, &(hit->acc)))   != eslOK) esl_fatal("allocation failure");
	      if ((status  = esl_strdup(om->desc, -1, &(hit->desc)))  != eslOK) esl_fatal("allocation failure");
	      hit->ndom       = cfg->ddef->ndom;
	      hit->nexpected  = cfg->ddef->nexpected;
	      hit->nregions   = cfg->ddef->nregions;
	      hit->nclustered = cfg->ddef->nclustered;
	      hit->noverlaps  = cfg->ddef->noverlaps;
	      hit->nenvelopes = cfg->ddef->nenvelopes;

	      hit->pre_score  = pre_score;
	      hit->pre_pvalue = esl_exp_surv (hit->pre_score,  om->evparam[p7_TAU], om->evparam[p7_LAMBDA]);

	      hit->score      = seq_score;
	      hit->pvalue     = P;
	      hit->sortkey    = -log(P);

	      hit->sum_score  = sum_score;
	      hit->sum_pvalue = esl_exp_surv (hit->sum_score,  om->evparam[p7_TAU], om->evparam[p7_LAMBDA]);

	      /* Transfer all domain coordinates (unthresholded for
               * now) with their alignment displays to the hit list,
               * associated with the sequence. Domain reporting will
               * be thresholded after complete hit list is collected,
               * because we probably need to know # of significant
               * seqs found to set domZ, and thence threshold and
               * count reported domains.
	       */
	      hit->dcl       = cfg->ddef->dcl;
	      cfg->ddef->dcl = NULL;
	      hit->best_domain = 0;
	      for (d = 0; d < hit->ndom; d++)
		{
		  Ld = hit->dcl[d].jenv - hit->dcl[d].ienv + 1;
		  hit->dcl[d].bitscore = hit->dcl[d].envsc + (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3)); 
		  if (! esl_opt_GetBoolean(go, "--nonull2")) 
		    seqbias = p7_FLogsum(0.0, log(cfg->bg->omega) + hit->dcl[d].domcorrection);
		  else  seqbias = 0.0;
		  hit->dcl[d].bitscore = (hit->dcl[d].bitscore - (nullsc + seqbias)) / eslCONST_LOG2;
		  hit->dcl[d].pvalue   = esl_exp_surv (hit->dcl[d].bitscore,  om->evparam[p7_TAU], om->evparam[p7_LAMBDA]);

		  if (hit->dcl[d].bitscore > hit->dcl[hit->best_domain].bitscore) hit->best_domain = d;
		}
	    }

	FINISH:
	  p7_domaindef_Reuse(cfg->ddef);
	  p7_oprofile_Destroy(om);
	} /* end, loop over models */
      if      (hstatus == eslEFORMAT)   p7_Fail("bad file format in HMM file %s",             cfg->hmmfile);
      else if (hstatus == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets",   cfg->hmmfile);
      else if (hstatus != eslEOF)       p7_Fail("Unexpected error in reading HMMs from %s",   cfg->hmmfile);

      /* Format the output. */
      p7_tophits_Sort(cfg->hitlist);
      apply_thresholds_to_output(cfg, cfg->hitlist);
      output_per_model_hitlist (go, cfg, cfg->hitlist);  fprintf(cfg->ofp, "\n");
      output_per_domain_hitlist(go, cfg, cfg->hitlist);  fprintf(cfg->ofp, "\n");
      fprintf(cfg->ofp, "//\n");

      p7_omx_Destroy(fwd);
      p7_omx_Destroy(bck);
      p7_omx_Destroy(oxf);
      p7_omx_Destroy(oxb);
      esl_sq_Reuse(sq);
    } /* end, loop over sequences */
  if (sstatus != eslEOF) p7_Fail("Sequence file %s has a format problem: read failed at line %d:\n%s\n",
				 cfg->seqfile, cfg->sqfp->linenumber, cfg->sqfp->errbuf);     

 
  esl_sq_Destroy(sq);
  return;
}



static int
output_header(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  if (cfg->my_rank > 0) return eslOK;
  p7_banner(cfg->ofp, go->argv[0], banner);
  
  fprintf(cfg->ofp, "# query sequence file:             %s\n", cfg->seqfile);
  fprintf(cfg->ofp, "# target HMM database:             %s\n", cfg->hmmfile);
  if (! esl_opt_IsDefault(go, "-o"))          fprintf(cfg->ofp, "# output directed to file:         %s\n",   esl_opt_GetString(go, "-o"));
  if (! esl_opt_IsDefault(go, "--hmmE"))      fprintf(cfg->ofp, "# HMM E-value threshold:        <= %g\n",   esl_opt_GetReal(go, "--hmmE"));
  if (! esl_opt_IsDefault(go, "--hmmT"))      fprintf(cfg->ofp, "# HMM bit score threshold:      <= %g\n",   esl_opt_GetReal(go, "--hmmT"));
  if (! esl_opt_IsDefault(go, "--domE"))      fprintf(cfg->ofp, "# domain E-value threshold:     <= %g\n",   esl_opt_GetReal(go, "--domE"));
  if (! esl_opt_IsDefault(go, "--domT"))      fprintf(cfg->ofp, "# domain bit score threshold:   <= %g\n",   esl_opt_GetReal(go, "--domT"));
  if (! esl_opt_IsDefault(go, "--cut_ga"))    fprintf(cfg->ofp, "# using GA bit score thresholds:   yes\n"); 
  if (! esl_opt_IsDefault(go, "--cut_nc"))    fprintf(cfg->ofp, "# using NC bit score thresholds:   yes\n");
  if (! esl_opt_IsDefault(go, "--cut_tc"))    fprintf(cfg->ofp, "# using TC bit score thresholds:   yes\n");
  if (! esl_opt_IsDefault(go, "--hmmZ"))      fprintf(cfg->ofp, "# HMM search space set to:    %.0f\n",    esl_opt_GetReal(go, "--hmmZ"));
  if (! esl_opt_IsDefault(go, "--domZ"))      fprintf(cfg->ofp, "# domain search space set to:      %.0f\n",    esl_opt_GetReal(go, "--domZ"));
  if (! esl_opt_IsDefault(go, "--max"))       fprintf(cfg->ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n");
  if (! esl_opt_IsDefault(go, "--F1"))        fprintf(cfg->ofp, "# MSV filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F1"));
  if (! esl_opt_IsDefault(go, "--F2"))        fprintf(cfg->ofp, "# Vit filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F2"));
  if (! esl_opt_IsDefault(go, "--F3"))        fprintf(cfg->ofp, "# Fwd filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F3"));
  if (! esl_opt_IsDefault(go, "--nonull2"))   fprintf(cfg->ofp, "# null2 bias corrections:          off\n");
  if (! esl_opt_IsDefault(go, "--nobias"))    fprintf(cfg->ofp, "# biased composition HMM filter:   off\n");
  if (! esl_opt_IsDefault(go, "--nonull2"))   fprintf(cfg->ofp, "# null2 bias corrections:          off\n");
  if (! esl_opt_IsDefault(go, "--Rdet") )     fprintf(cfg->ofp, "# RNG seed (run-to-run variation): reseed deterministically; minimize variation\n");
  if (! esl_opt_IsDefault(go, "--Rseed") )    fprintf(cfg->ofp, "# RNG seed (run-to-run variation): reseed to %d\n", esl_opt_GetInteger(go, "--Rseed"));
  if (! esl_opt_IsDefault(go, "--Rarb") )     fprintf(cfg->ofp, "# RNG seed (run-to-run variation): one arbitrary seed; allow run-to-run variation\n");
  if (! esl_opt_IsDefault(go, "--textw"))     fprintf(cfg->ofp, "# max ASCII text line length:      %d\n",     esl_opt_GetInteger(go, "--textw"));
  if (! esl_opt_IsDefault(go, "--notextw"))   fprintf(cfg->ofp, "# max ASCII text line length:      unlimited\n");
  fprintf(cfg->ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}


static int
set_pfam_bitscore_cutoffs(struct cfg_s *cfg, P7_OPROFILE *om)
{
  //  if (! (om->flags & cfg->model_cutoff_flag)) return eslFAIL;

  if (cfg->model_cutoff_flag == p7H_GA) { cfg->hmmT = om->cutoff[p7_GA1];  cfg->domT = om->cutoff[p7_GA2]; }
  if (cfg->model_cutoff_flag == p7H_TC) { cfg->hmmT = om->cutoff[p7_TC1];  cfg->domT = om->cutoff[p7_TC2]; }
  if (cfg->model_cutoff_flag == p7H_NC) { cfg->hmmT = om->cutoff[p7_NC1];  cfg->domT = om->cutoff[p7_NC2]; }
  return eslOK;
}


static int
hmm_is_reportable(struct cfg_s *cfg, float seq_score, float Pval)
{
  if      (  cfg->do_hmm_by_E   && Pval * cfg->hmmZ <= cfg->hmmE) return TRUE;
  else if (! cfg->do_hmm_by_E   && seq_score        >= cfg->hmmT) return TRUE;
  else return FALSE;
}


static int
domain_is_reportable(struct cfg_s *cfg, float dom_score, float Pval)
{
  if      (  cfg->do_dom_by_E   && Pval * cfg->domZ <= cfg->domE) return TRUE;
  else if (! cfg->do_dom_by_E   && dom_score        >= cfg->domT) return TRUE;
  else return FALSE;
}


/* domZ is also determined here */
static void
apply_thresholds_to_output(struct cfg_s *cfg, P7_TOPHITS *hitlist)
{
  int h, d;			/* counters over sequence hits, domains in sequences */
  
  /* First pass over HMMs, flag all reportable ones */
  hitlist->nreported = 0;
  for (h = 0; h < hitlist->N; h++)
    {
      if (hmm_is_reportable(cfg, hitlist->hit[h]->score, hitlist->hit[h]->pvalue))
	{
	  hitlist->hit[h]->is_reported = TRUE;
	  hitlist->nreported++;
	}
    }
  
  /* Now we can determined domZ, the effective search space in which additional domains are found */
  if (! cfg->domZ_is_fixed) cfg->domZ = (double) hitlist->nreported;

  /* Second pass is over domains, flagging reportable ones:
   * we always report at least the best single domain for each sequence, 
   * regardless of threshold.
   */
  for (h = 0; h < hitlist->N; h++)  
    if (hitlist->hit[h]->is_reported)
      for (d = 0; d < hitlist->hit[h]->ndom; d++)
	{
	  if (hitlist->hit[h]->best_domain == d ||
	      domain_is_reportable(cfg, hitlist->hit[h]->dcl[d].bitscore, hitlist->hit[h]->dcl[d].pvalue))
	    {
	      hitlist->hit[h]->nreported++;
	      hitlist->hit[h]->dcl[d].is_reported = TRUE;
	    }
	}
}





static int
output_per_model_hitlist(ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist)
{
  int h;
  int d;
  int namew       = ESL_MAX(8,  p7_tophits_GetMaxNameLength(hitlist));
  int descw;
  
  if (cfg->textw >  0) descw = ESL_MAX(32, cfg->textw - namew - 59);  /* 59 chars excluding desc is from the format: 22+2 +22+2 +8+2 +<name>+1 */
  else                 descw = INT_MAX;  

  fprintf(cfg->ofp, "Scores for complete sequence (score includes all domains):\n");
  fprintf(cfg->ofp, "%22s  %22s  %8s\n",                              " --- full sequence ---",   " --- best 1 domain ---",       "-#dom-");
  fprintf(cfg->ofp, "%9s %6s %5s  %9s %6s %5s  %5s %2s  %-*s %s\n", "E-value", " score", " bias", "E-value", " score",  "bias", "  exp",  "N", namew, "Model",    "Description");
  fprintf(cfg->ofp, "%9s %6s %5s  %9s %6s %5s  %5s %2s  %-*s %s\n", "-------", "------", "-----", "-------", "------", "-----", " ----", "--", namew, "--------", "-----------");

  for (h = 0; h < hitlist->N; h++)
    if (hitlist->hit[h]->is_reported)
      {
	d    = hitlist->hit[h]->best_domain;

	fprintf(cfg->ofp, "%9.2g %6.1f %5.1f  %9.2g %6.1f %5.1f  %5.1f %2d  %-*s %-.*s\n",
		hitlist->hit[h]->pvalue * (double) cfg->hmmZ,
		hitlist->hit[h]->score,
		hitlist->hit[h]->pre_score - hitlist->hit[h]->score, /* bias correction */
		hitlist->hit[h]->dcl[d].pvalue * (double) cfg->hmmZ,
		hitlist->hit[h]->dcl[d].bitscore,
		p7_FLogsum(0.0, log(cfg->bg->omega) + hitlist->hit[h]->dcl[d].domcorrection),
		hitlist->hit[h]->nexpected,
		hitlist->hit[h]->nreported,
		namew, hitlist->hit[h]->name,
		descw, (hitlist->hit[h]->desc == NULL ? "" : hitlist->hit[h]->desc));
      }
  return eslOK;
}

static int
output_per_domain_hitlist(ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist)
{
  int h, d;
  int nd;
  int namew, descw;

  fprintf(cfg->ofp, "Domain and alignment annotation for each model:\n");

  for (h = 0; h < hitlist->N; h++)
    if (hitlist->hit[h]->is_reported)
      {
	namew = strlen(hitlist->hit[h]->name);
	descw = (cfg->textw > 0 ?  ESL_MAX(32, cfg->textw - namew - 5) : INT_MAX);

	fprintf(cfg->ofp, ">> %s  %-.*s\n", hitlist->hit[h]->name, descw, (hitlist->hit[h]->desc == NULL ? "" : hitlist->hit[h]->desc));

	/* The domain table is 117 char wide:
           #  bit score    bias    E-value ind Evalue hmm from   hmm to    ali from   ali to    env from   env to    ali acc
         ---- --------- ------- ---------- ---------- -------- --------    -------- --------    -------- --------    -------
            1     123.4    23.1    9.7e-11     6.8e-9        3     1230 ..        1      492 []        2      490 .]    0.90
         1234 1234567.9 12345.7 1234567890 1234567890 12345678 12345678 .. 12345678 12345678 [] 12345678 12345678 .]    0.12
	*/
	fprintf(cfg->ofp, "  %4s %9s %7s %10s %10s %8s %8s %2s %8s %8s %2s %8s %8s %2s %7s\n",    "#", "bit score",    "bias",    "E-value", "ind Evalue", "hmm from",   "hmm to", "  ", "ali from", "ali to",   "  ", "env from",   "env to", "  ", "ali-acc");
	fprintf(cfg->ofp, "  %4s %9s %7s %10s %10s %8s %8s %2s %8s %8s %2s %8s %8s %2s %7s\n",  "---", "---------", "-------", "----------", "----------", "--------", "--------", "  ", "--------", "--------", "  ", "--------", "--------", "  ", "-------");

	nd = 0;
	for (d = 0; d < hitlist->hit[h]->ndom; d++)
	  if (hitlist->hit[h]->dcl[d].is_reported) 
	    {
	      nd++;
	      fprintf(cfg->ofp, "  %4d %9.1f %7.1f %10.2g %10.2g %8d %8d %c%c %8ld %8ld %c%c %8d %8d %c%c %7.2f\n",
		      nd,
		      hitlist->hit[h]->dcl[d].bitscore,
		      p7_FLogsum(0.0, log(cfg->bg->omega) + hitlist->hit[h]->dcl[d].domcorrection),
		      hitlist->hit[h]->dcl[d].pvalue * cfg->domZ,
		      hitlist->hit[h]->dcl[d].pvalue * cfg->hmmZ,
		      hitlist->hit[h]->dcl[d].ad->hmmfrom,
		      hitlist->hit[h]->dcl[d].ad->hmmto,
		      (hitlist->hit[h]->dcl[d].ad->hmmfrom == 1) ? '[' : '.',
		      (hitlist->hit[h]->dcl[d].ad->hmmto   == hitlist->hit[h]->dcl[d].ad->M) ? ']' : '.',
		      hitlist->hit[h]->dcl[d].ad->sqfrom,
		      hitlist->hit[h]->dcl[d].ad->sqto,
		      (hitlist->hit[h]->dcl[d].ad->sqfrom == 1) ? '[' : '.',
		      (hitlist->hit[h]->dcl[d].ad->sqto   == hitlist->hit[h]->dcl[d].ad->L) ? ']' : '.',
		      hitlist->hit[h]->dcl[d].ienv,
		      hitlist->hit[h]->dcl[d].jenv,
		      (hitlist->hit[h]->dcl[d].ienv == 1) ? '[' : '.',
		      (hitlist->hit[h]->dcl[d].jenv == hitlist->hit[h]->dcl[d].ad->L) ? ']' : '.',
		      (hitlist->hit[h]->dcl[d].oasc / (1.0 + fabs((float) (hitlist->hit[h]->dcl[d].jenv - hitlist->hit[h]->dcl[d].ienv)))));
	    }

	
	fprintf(cfg->ofp, "\n  Alignments for each domain:\n");
	nd = 0;
	for (d = 0; d < hitlist->hit[h]->ndom; d++)
	  if (hitlist->hit[h]->dcl[d].is_reported) 
	    {
	      nd++;
	      fprintf(cfg->ofp, "  == domain %d    score: %.1f bits;  conditional E-value: %.2g\n",
		      nd, 
		      hitlist->hit[h]->dcl[d].bitscore,
		      hitlist->hit[h]->dcl[d].pvalue * cfg->domZ);
	      p7_alidisplay_Print(cfg->ofp, hitlist->hit[h]->dcl[d].ad, 40, cfg->textw);
	      fprintf(cfg->ofp, "\n");
	    }
      }
  return eslOK;
}


static int 
output_search_statistics(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_STOPWATCH *w)
{
  fprintf(cfg->ofp, "Internal statistics summary:\n");
  fprintf(cfg->ofp, "----------------------------\n");
  fprintf(cfg->ofp, "Query sequence(s):           %15" PRId64 "  (%" PRId64 " residues)\n", cfg->nseq,    cfg->nres);
  fprintf(cfg->ofp, "Target HMM(s):               %15" PRId64 "  (%" PRId64 " nodes)\n",    cfg->nmodels, cfg->nnodes);

  fprintf(cfg->ofp, "Passed MSV filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n", 
	  cfg->n_past_msv,
	  (double) cfg->n_past_msv / (double) cfg->nmodels, 
	  cfg->F1 * (double) cfg->nmodels, 
	  cfg->F1);
  fprintf(cfg->ofp, "Passed Vit filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",   
	  cfg->n_past_vit,
	  (double) cfg->n_past_vit / (double) cfg->nmodels,
	  cfg->F2 * (double) cfg->nmodels,
	  cfg->F2);
  fprintf(cfg->ofp, "Passed Fwd filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",         
	  cfg->n_past_fwd, 
	  (double) cfg->n_past_fwd / (double) cfg->nmodels,
	  cfg->F3 * (double) cfg->nmodels,
	  cfg->F3);	  

  fprintf(cfg->ofp, "Initial search space (hmmZ): %15.0f  %s\n", cfg->hmmZ, cfg->hmmZ_is_fixed ? "[as set by --hmmZ on cmdline]" : "[actual number of target seqs]"); 
  fprintf(cfg->ofp, "Domain search space  (domZ): %15.0f  %s\n", cfg->domZ, cfg->domZ_is_fixed ? "[as set by --domZ on cmdline]" : "[number of seqs reported over threshold]"); 
  fprintf(cfg->ofp, "Mc/sec:                      %15.2f\n", 
	  (double) cfg->nres * (double) cfg->nnodes / (w->user * 1.0e6));
  esl_stopwatch_Display(cfg->ofp, w, "# CPU time: ");

  return eslOK;
}


