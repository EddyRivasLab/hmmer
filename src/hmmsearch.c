/* hmmsearch: search profile HMM(s) against a sequence database.
 * 
 * SRE, Thu Dec 20 07:07:25 2007 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif 

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_histogram.h"
#include "esl_stopwatch.h"
#include "esl_gumbel.h"
#include "esl_exponential.h"
#include "esl_vectorops.h"

#include "hmmer.h"


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles  reqs   incomp  help   docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "show brief help on version and usage",                         1 },
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL, "direct output to file <f>, not stdout",                        1 },

  { "--seqE",       eslARG_REAL,  "10.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "E-value cutoff for reporting significant sequence hits",       2 },
  { "--seqT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "bit score cutoff for reporting significant sequence hits",     2 },
  { "--domE",       eslARG_REAL,"1000.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "E-value cutoff for reporting individual domains",              2 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "bit score cutoff for reporting individual domains",            2 },
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--seqE,--seqT,--domE,--domT", "use GA gathering threshold bit score cutoffs in <hmmfile>",    2 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--seqE,--seqT,--domE,--domT", "use NC noise threshold bit score cutoffs in <hmmfile>",        2 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--seqE,--seqT,--domE,--domT", "use TC trusted threshold bit score cutoffs in <hmmfile>",      2 },
  { "--seqZ",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                          "set # of comparisons done, for E-value calculation",           2 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                          "set # of significant seqs, for domain E-value calculation",    2 },

  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",  3 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,      NULL,  NULL, "--max",          "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             3 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,      NULL,  NULL, "--max",          "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             3 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,      NULL,  NULL, "--max",          "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             3 },
  { "--biasfilter", eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, "--max",          "turn on composition bias filter (more speed, less power)", 3 },

  { "--nonull2",    eslARG_NONE,   NULL,  NULL, NULL,      NULL,  NULL,  NULL, "turn off biased composition score corrections",            4 },
  { "--seed",       eslARG_INT,    "42",  NULL, NULL,      NULL,  NULL,  NULL, "set random number generator seed",                         4 },  
  { "--timeseed",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "use arbitrary random number generator seed (by time())",   4 },  
#ifdef HAVE_MPI
  // { "--stall",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "arrest after start: for debugging MPI under gdb",          4 },  
  //  { "--mpi",       eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "run as an MPI parallel program",                           4 },
#endif 
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[options] <query hmmfile> <target seqfile>";
static char banner[] = "search profile HMM(s) against a sequence database";



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

  int     do_seq_by_E;		/* TRUE to cut sequence reporting off by E  */
  double  seqE;	                /* sequence E-value threshold               */
  double  seqT;	                /* sequence bit score threshold             */
  int     do_dom_by_E;		/* TRUE to cut sequence reporting off by E  */
  double  domE;	                /* domain E-value threshold                 */
  double  domT;	                /* domain bit score threshold               */
  int     model_cutoff_flag;    /* FALSE, p7H_GA, p7H_TC, or p7H_NC         */

  double  seqZ;			/* # of seqs searched for E-value purposes  */
  double  domZ;			/* # signific seqs for domain E-values      */
  int     seqZ_is_fixed;	/* TRUE if seqZ was set on cmd line         */
  int     domZ_is_fixed;	/* TRUE if domZ was set on cmd line         */

  int              do_max;	/* TRUE to run in slow/max mode             */
  double           F1;		/* MSV filter threshold                     */
  double           F2;		/* Viterbi filter threshold                 */
  double           F3;		/* uncorrected Forward filter threshold     */

  uint64_t         nseq;	/* true number of sequences searched        */
  uint64_t         nmodels;	/* number of models searched                */
  uint64_t         nres;	/* # of residues searched                   */
  uint64_t         nnodes;	/* # of model nodes searched                */
  uint64_t         n_past_msv;	/* # of seqs that pass the MSVFilter()      */
  uint64_t         n_past_vit;	/* # of seqs that pass the ViterbiFilter()  */
  uint64_t         n_past_fwd;	/* # of seqs that pass the ForwardFilter()  */

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
static void serial_master  (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  set_pfam_bitscore_cutoffs(struct cfg_s *cfg, P7_HMM *hmm);
static int  sequence_is_reportable   (struct cfg_s *cfg, float seq_score, float Pval);
static int  domain_is_reportable     (struct cfg_s *cfg, float dom_score, float Pval);
static void apply_thresholds_to_output(struct cfg_s *cfg, P7_TOPHITS *hitlist);
static int  output_header              (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  output_per_sequence_hitlist(ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist);
static int  output_per_domain_hitlist  (ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist);
static int  output_search_statistics   (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_STOPWATCH *w);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                 */
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
    {       /* No MPI? Then we're just the serial master. */
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
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

      puts("\nOptions controlling significance thresholds for reporting:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      puts("\nOptions controlling acceleration heuristics:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 

      puts("\nOther expert options:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 2)      { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((cfg->hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <seqfile> argument on command line"); goto ERROR; }
  if ((cfg->seqfile = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <hmmfile> argument on command line"); goto ERROR; }

  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere basic options are:");
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
  cfg->format   = eslSQFILE_UNKNOWN;    /* eventually, allow options to set this            */
  cfg->abc      = NULL;		        /* to be initialized later */
  cfg->bg       = NULL;		        /* to be initialized later */
  cfg->mode     = p7_LOCAL;
  cfg->hitlist  = p7_tophits_Create();  /* master accumulates complete list; workers have partial lists. */

  /* Reporting thresholds */
  cfg->do_seq_by_E       = TRUE;
  cfg->seqE              = esl_opt_GetReal(go, "--seqE");
  cfg->seqT              = 0.0;
  cfg->do_dom_by_E       = TRUE;
  cfg->domE              = esl_opt_GetReal(go, "--domE");
  cfg->domT              = 0.0;
  cfg->model_cutoff_flag = 0;
  
  if (! esl_opt_IsDefault(go, "--seqT")) 
    {
      cfg->seqT        = esl_opt_GetReal(go, "--seqT"); 
      cfg->do_seq_by_E = FALSE;
    } 
  if (! esl_opt_IsDefault(go, "--domT")) 
    {
      cfg->domT        = esl_opt_GetReal(go, "--domT"); 
      cfg->do_dom_by_E = FALSE;
    }
  if (esl_opt_GetBoolean(go, "--cut_ga"))
    {
      cfg->seqT        = cfg->domT        = 0.0;
      cfg->do_seq_by_E = cfg->do_dom_by_E = FALSE;
      cfg->model_cutoff_flag = p7H_GA;
    }
  if (esl_opt_GetBoolean(go, "--cut_nc"))
    {
      cfg->seqT        = cfg->domT        = 0.0;
      cfg->do_seq_by_E = cfg->do_dom_by_E = FALSE;
      cfg->model_cutoff_flag = p7H_NC;
    }
  if (esl_opt_GetBoolean(go, "--cut_tc"))
    {
      cfg->seqT        = cfg->domT        = 0.0;
      cfg->do_seq_by_E = cfg->do_dom_by_E = FALSE;
      cfg->model_cutoff_flag = p7H_TC;
    }

  /* Search space sizes for E-value calculations */
  cfg->seqZ_is_fixed = FALSE;
  cfg->domZ_is_fixed = FALSE;
  if (! esl_opt_IsDefault(go, "--seqZ")) 
    {
      cfg->seqZ_is_fixed = TRUE;
      cfg->seqZ          = esl_opt_GetReal(go, "--seqZ");
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
  //  cfg->do_stall = esl_opt_GetBoolean(go, "--stall");

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
 * dumping a message and exiting; we have worker processes to shut
 * down, and possibly other cleanup to do, so we must try to delay
 * resolution of any error messages until after we attempt to clean
 * up. Therefore errors return (code, errmsg) by the ESL_FAIL mech.
 */
static int
init_master_cfg(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  char *filename;
  int   status;

  status = p7_hmmfile_Open(cfg->hmmfile, NULL, &(cfg->hfp));
  if      (status == eslENOTFOUND) ESL_FAIL(status, errbuf, "Failed to open hmm file %s for reading.\n",     cfg->hmmfile);
  else if (status == eslEFORMAT)   ESL_FAIL(status, errbuf, "Unrecognized format, trying to open hmm file %s for reading.\n", cfg->hmmfile);
  else if (status != eslOK)        ESL_FAIL(status, errbuf, "Unexpected error %d in opening hmm file %s.\n", status, cfg->hmmfile);  
  
  status = esl_sqfile_Open(cfg->seqfile, cfg->format, p7_SEQDBENV, &(cfg->sqfp));
  if      (status == eslENOTFOUND) ESL_FAIL(status, errbuf, "Failed to open sequence file %s for reading\n",      cfg->seqfile);
  else if (status == eslEFORMAT)   ESL_FAIL(status, errbuf, "Sequence file %s is empty or misformatted\n",        cfg->seqfile);
  else if (status == eslEINVAL)    ESL_FAIL(status, errbuf, "Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        ESL_FAIL(status, errbuf, "Unexpected error %d opening sequence file %s\n", status, cfg->seqfile);

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
  if (esl_opt_GetBoolean(go, "--timeseed")) cfg->r = esl_randomness_CreateTimeseeded();
  else                                      cfg->r = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  if (cfg->r == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create random number generator");

  if ((cfg->ddef = p7_domaindef_Create(cfg->r)) == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create domain definition workbook");
  return eslOK;
}




static void
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  P7_HMM          *hmm     = NULL;     /* query HMM                               */
  P7_PROFILE      *gm      = NULL;     /* profile HMM                             */
  P7_OPROFILE     *om      = NULL;     /* optimized profiile                      */
  ESL_SQ          *sq      = NULL;     /* target sequence                         */
  P7_OMX          *fwd     = NULL;     /* DP matrix                               */
  P7_OMX          *bck     = NULL;     /* DP matrix                               */
  P7_OMX          *oxf     = NULL;     /* optimized DP matrix for Forward         */
  P7_OMX          *oxb     = NULL;     /* optimized DP matrix for Backward        */
  P7_HIT          *hit     = NULL;     /* ptr to the current hit output data      */
  float            usc, vfsc, fwdsc;   /* filter scores                           */
  float            filtersc;           /* HMM null filter score                   */
  float            nullsc;             /* null model score                        */
  float            seqbias;  
  float            seq_score;          /* the corrected per-seq bit score */
  float            sum_score;	       /* the corrected reconstruction score for the seq */
  float            pre_score, pre2_score; /* uncorrected bit scores for seq */
  double           P;		       /* P-value of a hit */
  double           E;		       /* bound on E-value of a hit  */
  int              Ld;		       /* # of residues in envelopes */
  char             errbuf[eslERRBUFSIZE];
  int              status, hstatus, sstatus;
  int              d;

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) esl_fatal(errbuf);
  if ((status = init_worker_cfg(go, cfg, errbuf)) != eslOK) esl_fatal(errbuf);
 
  output_header(go, cfg);

  while ((hstatus = p7_hmmfile_Read(cfg->hfp, &(cfg->abc), &hmm)) == eslOK) 
    {
      fwd = p7_omx_Create(hmm->M, 400, 400); /* initial alloc is for M, L=400; will grow as needed */
      bck = p7_omx_Create(hmm->M, 400, 400);	
      oxf = p7_omx_Create(hmm->M, 0,   400); /* one-row parsing matrix, O(M+L) */
      oxb = p7_omx_Create(hmm->M, 0,   400);     

      /* One time only initializations after abc becomes known: */
      if (cfg->bg == NULL) cfg->bg = p7_bg_Create(cfg->abc);
      if (sq      == NULL) sq      = esl_sq_CreateDigital(cfg->abc);

      /* Pfam bit score thresholds may be in use */
      if (cfg->model_cutoff_flag && set_pfam_bitscore_cutoffs(cfg, hmm) != eslOK) p7_Fail("requested score cutoffs not available in model %s\n", hmm->name);

      if (esl_opt_GetBoolean(go, "--biasfilter"))
	p7_bg_SetFilter(cfg->bg, hmm->M, hmm->compo); /* EXPERIMENTAL */
      
      gm = p7_profile_Create (hmm->M, cfg->abc);
      om = p7_oprofile_Create(hmm->M, cfg->abc);
      p7_ProfileConfig(hmm, cfg->bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; MSVFilter requires local mode */
      /* additionally this *must* be a multihit mode. unihit modes have rare numerical problems in posterior decoding. [J3/119-120] */

      p7_oprofile_Convert(gm, om);     /* <om> is now p7_LOCAL, multihit */

      fprintf(cfg->ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M);
      if (hmm->acc  != NULL) fprintf(cfg->ofp, "Accession:   %s\n", hmm->acc);
      if (hmm->desc != NULL) fprintf(cfg->ofp, "Description: %s\n", hmm->desc);
      cfg->nmodels++;
      cfg->nnodes += om->M;

      while ( (sstatus = esl_sqio_Read(cfg->sqfp, sq)) == eslOK)
	{
	  cfg->nseq++;
	  if (! cfg->seqZ_is_fixed) cfg->seqZ = (double) cfg->nseq;
	  cfg->nres += sq->n;
	  p7_oprofile_ReconfigLength(om, sq->n);
	  p7_omx_GrowTo(oxf, om->M, 0, sq->n);    /* expand the one-row omx if needed */

	  /* Null model score for this sequence.  */
	  p7_bg_SetLength(cfg->bg, sq->n);
	  p7_bg_NullOne  (cfg->bg, sq->dsq, sq->n, &nullsc);

	  /* First level filter: the MSV filter, multihit with <om> */
	  p7_MSVFilter(sq->dsq, sq->n, om, oxf, &usc);
	  seq_score = (usc - nullsc) / eslCONST_LOG2;
	  P = esl_gumbel_surv(seq_score,  hmm->evparam[p7_MU],  hmm->evparam[p7_LAMBDA]);
	  E = P * (double) cfg->seqZ;
	  if (P > cfg->F1 && ! sequence_is_reportable(cfg, seq_score, P)) { esl_sq_Reuse(sq); continue; } 

	  /* HMM filtering */
	  if (esl_opt_GetBoolean(go, "--biasfilter"))
	    {
	      p7_bg_FilterScore(cfg->bg, sq->dsq, sq->n, &filtersc);
	      seq_score = (usc - filtersc) / eslCONST_LOG2;
	      P = esl_gumbel_surv(seq_score,  hmm->evparam[p7_MU],  hmm->evparam[p7_LAMBDA]);
	      E = P * (double) cfg->seqZ;
	      if (P > cfg->F1 && ! sequence_is_reportable(cfg, seq_score, P)) { esl_sq_Reuse(sq); continue; } 
	    }
	  else filtersc = nullsc;
	  cfg->n_past_msv++;

	  //printf("Past MSV: %-20s %5d %8.4f %8.4f %8.4f %8.4f\n", sq->name, (int) sq->n, nullsc, filtersc, usc, P);

	  /* Second level filter: ViterbiFilter(), multihit with <om> */
	  if (P > cfg->F2) 		
	    {
	      p7_ViterbiFilter(sq->dsq, sq->n, om, oxf, &vfsc);  
	      seq_score = (vfsc-filtersc) / eslCONST_LOG2;
	      P  = esl_gumbel_surv(seq_score,  hmm->evparam[p7_MU],  hmm->evparam[p7_LAMBDA]);
	      E = P * (double) cfg->seqZ;
	      if (P > cfg->F2 && ! sequence_is_reportable(cfg, seq_score, P)) { esl_sq_Reuse(sq); continue; } 
	      //  printf("Past VIT: %-20s %5d %8.4f %8.4f %8.4f %8.4f\n", sq->name, (int) sq->n, nullsc, filtersc, vfsc, P);
	    }
	  cfg->n_past_vit++;

	  /* Parse it with Forward and obtain its real Forward score. */
	  p7_ForwardParser(sq->dsq, sq->n, om, oxf, &fwdsc);
	  seq_score = (fwdsc-filtersc) / eslCONST_LOG2;
	  P = esl_exp_surv(seq_score,  hmm->evparam[p7_TAU],  hmm->evparam[p7_LAMBDA]);
	  E = P * (double) cfg->seqZ;
	  if (P > cfg->F3 && ! sequence_is_reportable(cfg, seq_score, P)) { esl_sq_Reuse(sq); continue; }
	  cfg->n_past_fwd++;

	  //printf("forward = %g\n", fwdsc);
	  //  printf("Past FWD: %-20s %5d %8.4f %8.4f %8.4f %8.4f\n", sq->name, (int) sq->n, nullsc, filtersc, final_sc, P);

	  /* ok, it's for real; do the Backwards parser pass and hand
	   * it off to domain definition logic. 
	   */
	  p7_omx_GrowTo(oxb, hmm->M, 0, sq->n);
	  p7_BackwardParser(sq->dsq, sq->n, om, oxf, oxb, NULL);
	  status = p7_domaindef_ByPosteriorHeuristics(sq, om, oxf, oxb, fwd, bck, cfg->ddef);
	  if (status == eslERANGE) { 
	    fprintf(stderr, "WARNING: posterior decoding failed on %s; skipping it!\n", sq->name);
	    esl_sq_Reuse(sq);
	    continue;
	  } else if (status != eslOK) {
	    fprintf(stderr, "WARNING: domain definition failed due to unknown problem on %s; skipping it!\n", sq->name);
	    esl_sq_Reuse(sq);
	    continue;
	  }

	  /* What's the per-seq score, and is it significant enough to be reported? */
	  /* Figure out the sum of null2 corrections to be added to the null score */
	  if (! esl_opt_GetBoolean(go, "--nonull2"))
	    {
	      seqbias = esl_vec_FSum(cfg->ddef->n2sc, sq->n+1);
	      seqbias = p7_FLogsum(0.0, log(cfg->bg->omega) + seqbias);
	    }
	  else seqbias = 0.0;
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
	  P =  esl_exp_surv (seq_score,  hmm->evparam[p7_TAU], hmm->evparam[p7_LAMBDA]);
	  if (sequence_is_reportable(cfg, seq_score, P))
	    {
	      p7_tophits_CreateNextHit(cfg->hitlist, &hit);
	      if ((status  = esl_strdup(sq->name, -1, &(hit->name)))  != eslOK) esl_fatal("allocation failure");
	      if ((status  = esl_strdup(sq->acc,  -1, &(hit->acc)))   != eslOK) esl_fatal("allocation failure");
	      if ((status  = esl_strdup(sq->desc, -1, &(hit->desc)))  != eslOK) esl_fatal("allocation failure");
	      hit->ndom       = cfg->ddef->ndom;
	      hit->nexpected  = cfg->ddef->nexpected;
	      hit->nregions   = cfg->ddef->nregions;
	      hit->nclustered = cfg->ddef->nclustered;
	      hit->noverlaps  = cfg->ddef->noverlaps;
	      hit->nenvelopes = cfg->ddef->nenvelopes;

	      hit->pre_score  = pre_score;
	      hit->pre_pvalue = esl_exp_surv (hit->pre_score,  hmm->evparam[p7_TAU], hmm->evparam[p7_LAMBDA]);

	      hit->score      = seq_score;
	      hit->pvalue     = P;
	      hit->sortkey    = -log(P);

	      hit->sum_score  = sum_score;
	      hit->sum_pvalue = esl_exp_surv (hit->sum_score,  hmm->evparam[p7_TAU], hmm->evparam[p7_LAMBDA]);

	      /* Transfer all domain coordinates (unthresholded for
               * now) with their alignment displays to the hit list,
               * associated with the sequence. Domain reporting will
               * be thresholded after complete hit list is collected,
               * because we probably need to know # of significant
               * seqs found to set domZ, and thence threshold and
               * count reported domains.
	       */
	      hit->dcl         = cfg->ddef->dcl;
	      cfg->ddef->dcl   = NULL;
	      hit->best_domain = 0;
	      for (d = 0; d < hit->ndom; d++)
		{
		  Ld = hit->dcl[d].jenv - hit->dcl[d].ienv + 1;
		  hit->dcl[d].bitscore = hit->dcl[d].envsc + (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3)); 
		  if (! esl_opt_GetBoolean(go, "--nonull2")) 
		    seqbias = p7_FLogsum(0.0, log(cfg->bg->omega) + hit->dcl[d].domcorrection);
		  else
		    seqbias = 0.0;
		  hit->dcl[d].bitscore = (hit->dcl[d].bitscore - (nullsc + seqbias)) / eslCONST_LOG2;
		  hit->dcl[d].pvalue   = esl_exp_surv (hit->dcl[d].bitscore,  hmm->evparam[p7_TAU], hmm->evparam[p7_LAMBDA]);

		  if (hit->dcl[d].bitscore > hit->dcl[hit->best_domain].bitscore) hit->best_domain = d;
		}
	    }

	  esl_sq_Reuse(sq);
	  p7_domaindef_Reuse(cfg->ddef);
	}
      if (sstatus != eslEOF) p7_Fail("Sequence file %s has a format problem: read failed at line %d:\n%s\n",
				     cfg->seqfile, cfg->sqfp->linenumber, cfg->sqfp->errbuf);     

      /* Format the output. */
      p7_tophits_Sort(cfg->hitlist);
      apply_thresholds_to_output(cfg, cfg->hitlist);
      output_per_sequence_hitlist(go, cfg, cfg->hitlist);  fprintf(cfg->ofp, "\n");
      output_per_domain_hitlist  (go, cfg, cfg->hitlist);  fprintf(cfg->ofp, "\n");
      fprintf(cfg->ofp, "//\n");

      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
      p7_oprofile_Destroy(om);
    }
  if      (hstatus == eslEOD)       p7_Fail("read failed, HMM file %s may be truncated?", cfg->hmmfile);
  else if (hstatus == eslEFORMAT)   p7_Fail("bad file format in HMM file %s",             cfg->hmmfile);
  else if (hstatus == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets",   cfg->hmmfile);
  else if (hstatus != eslEOF)       p7_Fail("Unexpected error in reading HMMs from %s",   cfg->hmmfile);

  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  p7_omx_Destroy(oxf);
  p7_omx_Destroy(oxb);
  esl_sq_Destroy(sq);
  return;
}

static int
set_pfam_bitscore_cutoffs(struct cfg_s *cfg, P7_HMM *hmm)
{
  if (! (hmm->flags & cfg->model_cutoff_flag)) return eslFAIL;

  if (cfg->model_cutoff_flag == p7H_GA) { cfg->seqT = hmm->cutoff[p7_GA1];  cfg->domT = hmm->cutoff[p7_GA2]; }
  if (cfg->model_cutoff_flag == p7H_TC) { cfg->seqT = hmm->cutoff[p7_TC1];  cfg->domT = hmm->cutoff[p7_TC2]; }
  if (cfg->model_cutoff_flag == p7H_NC) { cfg->seqT = hmm->cutoff[p7_NC1];  cfg->domT = hmm->cutoff[p7_NC2]; }
  return eslOK;
}

static int
sequence_is_reportable(struct cfg_s *cfg, float seq_score, float Pval)
{
  if      (  cfg->do_seq_by_E   && Pval * cfg->seqZ <= cfg->seqE) return TRUE;
  else if (! cfg->do_seq_by_E   && seq_score        >= cfg->seqT) return TRUE;
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
  
  /* First pass over sequences, flag all reportable ones */
  hitlist->nreported = 0;
  for (h = 0; h < hitlist->N; h++)
    {
      if (sequence_is_reportable(cfg, hitlist->hit[h]->score, hitlist->hit[h]->pvalue))
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
output_header(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  if (cfg->my_rank > 0) return eslOK;
  p7_banner(cfg->ofp, go->argv[0], banner);
  
  fprintf(cfg->ofp, "# query HMM file:                  %s\n", cfg->hmmfile);
  fprintf(cfg->ofp, "# target sequence database:        %s\n", cfg->seqfile);
  if (! esl_opt_IsDefault(go, "-o"))          fprintf(cfg->ofp, "# output directed to file:      %s\n",      esl_opt_GetString(go, "-o"));
  if (! esl_opt_IsDefault(go, "--seqE"))      fprintf(cfg->ofp, "# sequence E-value threshold:   <= %g\n",   esl_opt_GetReal(go, "--seqE"));
  if (! esl_opt_IsDefault(go, "--seqT"))      fprintf(cfg->ofp, "# sequence bit score threshold: <= %g\n",   esl_opt_GetReal(go, "--seqT"));
  if (! esl_opt_IsDefault(go, "--domE"))      fprintf(cfg->ofp, "# domain E-value threshold:     <= %g\n",   esl_opt_GetReal(go, "--domE"));
  if (! esl_opt_IsDefault(go, "--domT"))      fprintf(cfg->ofp, "# domain bit score threshold:   <= %g\n",   esl_opt_GetReal(go, "--domT"));
  if (! esl_opt_IsDefault(go, "--cut_ga"))    fprintf(cfg->ofp, "# using GA bit score thresholds:   yes\n"); 
  if (! esl_opt_IsDefault(go, "--cut_nc"))    fprintf(cfg->ofp, "# using NC bit score thresholds:   yes\n");
  if (! esl_opt_IsDefault(go, "--cut_tc"))    fprintf(cfg->ofp, "# using TC bit score thresholds:   yes\n");
  if (! esl_opt_IsDefault(go, "--seqZ"))      fprintf(cfg->ofp, "# sequence search space set to:    %.0f\n",    esl_opt_GetReal(go, "--seqZ"));
  if (! esl_opt_IsDefault(go, "--domZ"))      fprintf(cfg->ofp, "# domain search space set to:      %.0f\n",    esl_opt_GetReal(go, "--domZ"));
  if (! esl_opt_IsDefault(go, "--max"))       fprintf(cfg->ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n");
  if (! esl_opt_IsDefault(go, "--F1"))        fprintf(cfg->ofp, "# MSV filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F1"));
  if (! esl_opt_IsDefault(go, "--F2"))        fprintf(cfg->ofp, "# Vit filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F2"));
  if (! esl_opt_IsDefault(go, "--F3"))        fprintf(cfg->ofp, "# Fwd filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F3"));
  if (! esl_opt_IsDefault(go, "--biasfilter"))fprintf(cfg->ofp, "# biased composition HMM filter:   on\n");
  if (! esl_opt_IsDefault(go, "--nonull2"))   fprintf(cfg->ofp, "# null2 bias corrections:          off\n");
  if (! esl_opt_IsDefault(go, "--seed"))      fprintf(cfg->ofp, "# random number seed set to:       %d\n",  esl_opt_GetInteger(go, "--seed"));
  if (! esl_opt_IsDefault(go, "--timeseed"))  fprintf(cfg->ofp, "# random number seed set to:       %ld\n",  cfg->r->seed);
  fprintf(cfg->ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}


static int
output_per_sequence_hitlist(ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist)
{
  int    h;
  int    d;
  int    namew       = ESL_MAX(8,  p7_tophits_GetMaxNameLength(hitlist));
  
  fprintf(cfg->ofp, "Scores for complete sequences (score includes all domains):\n");

  fprintf(cfg->ofp, "%25s  %25s  %13s\n", "   --- full sequence ----", "---- single best dom ----", "-- #doms --");
  fprintf(cfg->ofp, "%10s %7s %6s  %7s %6s %10s  %6s %5s  %-*s %s\n", "E-value", "  score", "  bias", "  score",   "bias",    "E-value", "   exp",     "N", namew, "Sequence", "Description");
  fprintf(cfg->ofp, "%10s %7s %6s  %7s %6s %10s  %6s %5s  %-*s %s\n", "-------", "-------", "------", "-------", "------", "----------", " -----", "-----", namew, "--------", "-----------");

  for (h = 0; h < hitlist->N; h++)
    if (hitlist->hit[h]->is_reported)
      {
	d    = hitlist->hit[h]->best_domain;

	fprintf(cfg->ofp, "%10.2g %7.1f %6.1f  %7.1f %6.1f %10.2g  %6.1f %5d  %-*s %s\n",
		hitlist->hit[h]->pvalue * (double) cfg->seqZ,
		hitlist->hit[h]->score,
		hitlist->hit[h]->pre_score - hitlist->hit[h]->score, /* bias correction */
		hitlist->hit[h]->dcl[d].bitscore,
		p7_FLogsum(0.0, log(cfg->bg->omega) + hitlist->hit[h]->dcl[d].domcorrection),
		hitlist->hit[h]->dcl[d].pvalue * (double) cfg->seqZ,
		hitlist->hit[h]->nexpected,
		hitlist->hit[h]->nreported,
		namew, 
		hitlist->hit[h]->name,
		hitlist->hit[h]->desc);
      }
  return eslOK;
}


static int
output_per_domain_hitlist(ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist)
{
  int h, d;
  int nd;

  fprintf(cfg->ofp, "Domain and alignment annotation for each sequence:\n");

  for (h = 0; h < hitlist->N; h++)
    if (hitlist->hit[h]->is_reported)
      {
	fprintf(cfg->ofp, ">> %s  %s\n", hitlist->hit[h]->name, hitlist->hit[h]->desc);

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
		      hitlist->hit[h]->dcl[d].pvalue * cfg->seqZ,
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
	      p7_alidisplay_Print(cfg->ofp, hitlist->hit[h]->dcl[d].ad, 40, -1);
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
  fprintf(cfg->ofp, "Query HMM(s):                %15ld  (%ld nodes)\n",    cfg->nmodels, cfg->nnodes);
  fprintf(cfg->ofp, "Target sequences:            %15ld  (%ld residues)\n", cfg->nseq,    cfg->nres);
  fprintf(cfg->ofp, "Passed MSV filter:           %15ld  (%.6g); expected %.1f (%.6g)\n", 
	  cfg->n_past_msv,
	  (double) cfg->n_past_msv / (double) cfg->nseq, 
	  cfg->F1 * (double) cfg->nseq, 
	  cfg->F1);
  fprintf(cfg->ofp, "Passed Vit filter:           %15ld  (%.6g); expected %.1f (%.6g)\n",   
	  cfg->n_past_vit,
	  (double) cfg->n_past_vit / (double) cfg->nseq,
	  cfg->F2 * (double) cfg->nseq,
	  cfg->F2);
  fprintf(cfg->ofp, "Passed Fwd filter:           %15ld  (%.6g); expected %.1f (%.6g)\n",         
	  cfg->n_past_fwd, 
	  (double) cfg->n_past_fwd / (double) cfg->nseq,
	  cfg->F3 * (double) cfg->nseq,
	  cfg->F3);	  

  fprintf(cfg->ofp, "Initial search space (seqZ): %15.0f  %s\n", cfg->seqZ, cfg->seqZ_is_fixed ? "[as set by --seqZ on cmdline]" : "[actual number of target seqs]"); 
  fprintf(cfg->ofp, "Domain search space  (domZ): %15.0f  %s\n", cfg->domZ, cfg->domZ_is_fixed ? "[as set by --domZ on cmdline]" : "[number of seqs reported over threshold]"); 
  fprintf(cfg->ofp, "Mc/sec:                      %15.2f\n", 
	  (double) cfg->nres * (double) cfg->nnodes / (w->user * 1.0e6));
  esl_stopwatch_Display(cfg->ofp, w, "# CPU time: ");

  return eslOK;
}


