/* hmmpfam: search sequence(s) against a profile HMM database
 * 
 * SRE, Mon Oct 20 08:28:05 2008 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif 

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


#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles  reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "show brief help on version and usage",               1 },
  { "-o",        eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL, "direct output to file <f>, not stdout",              1 },
  { "-E",        eslARG_REAL,  "10.0", NULL,"x>0",      NULL,  NULL,  NULL, "E-value cutoff for sequences",                       1 },
  { "-Z",        eslARG_INT,    FALSE, NULL,"n>0",      NULL,  NULL,  NULL, "set # of comparisons done, for E-value calculation", 1 },

  { "--F1",      eslARG_REAL,  "0.02", NULL, NULL,      NULL,  NULL,  NULL, "MSV threshold: promote hits w/ P <= F1",             2 },
  { "--F2",      eslARG_REAL,  "1e-3", NULL, NULL,      NULL,  NULL,  NULL, "Vit threshold: promote hits w/ P <= F2",             2 },
  { "--F3",      eslARG_REAL,  "1e-5", NULL, NULL,      NULL,  NULL,  NULL, "Fwd threshold: promote hits w/ P <= F3",             2 },

  { "--textw",   eslARG_INT,    "120", NULL, NULL,      NULL,  NULL,  NULL, "sets maximum ASCII text output line length",         4 },
  { "--seed",    eslARG_INT,    NULL,  NULL, NULL,      NULL,  NULL,  NULL, "random number generator seed",                       4 },  
  { "--stall",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "arrest after start: for debugging MPI under gdb",    4 },  
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "run as an MPI parallel program",                     4 },
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
static int  output_per_model_hitlist (ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist);
static int  output_per_domain_hitlist(ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist);
static int  output_alignments        (ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist);
static int  output_search_statistics (ESL_GETOPTS *go, struct cfg_s *cfg, ESL_STOPWATCH *w);

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
    {		
      /* No MPI? Then we're just the serial master. */
      serial_master(go, &cfg);
      esl_stopwatch_Stop(w);
    }      

  output_search_statistics   (go, &cfg, w);            
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
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s", go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s", go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
      //      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); /* 2= group; 2 = indentation; 80=textwidth*/
      //      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); /* 3= group; 2 = indentation; 80=textwidth*/
      //      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); /* 3= group; 2 = indentation; 80=textwidth*/
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
  cfg->format   = eslSQFILE_UNKNOWN;    /* eventually, allow options to set this            */

  /* in serial, single master has the master hitlists, of course; 
   * in MPI, workers assemble partial hitlists, send to master, 
   * master accumulates in its main hitlist. 
   */
  cfg->hitlist    = p7_tophits_Create();
  cfg->nseq       = 0;
  cfg->nmodels    = 0;
  cfg->nres       = 0;
  cfg->nnodes     = 0;
  cfg->n_past_msv = 0;
  cfg->n_past_vit = 0;
  cfg->n_past_fwd = 0;

  /* These will be initialized when we read first HMM and know our alphabet: */
  cfg->abc      = NULL;		        
  cfg->bg       = NULL;		        

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
  if (esl_opt_IsDefault(go, "--seed")) cfg->r = esl_randomness_CreateTimeseeded();
  else                                 cfg->r = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  if (cfg->r == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create random number generator");

  if ((cfg->ddef = p7_domaindef_Create(cfg->r)) == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create domain definition workbook");
  return eslOK;
}





static void
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  P7_OPROFILE     *om      = NULL;     /* optimized profiile                      */
  ESL_SQ          *sq      = NULL;     /* target sequence                         */
  P7_TRACE        *tr      = NULL;     /* trace of hmm aligned to sq              */
  P7_OMX          *fwd     = NULL;     /* DP matrix                               */
  P7_OMX          *bck     = NULL;     /* DP matrix                               */
  P7_OMX          *oxf     = NULL;     /* optimized DP matrix for Forward         */
  P7_OMX          *oxb     = NULL;     /* optimized DP matrix for Backward        */
  P7_HIT          *hit     = NULL;     /* ptr to the current hit output data      */
  double F1                = esl_opt_GetReal(go, "--F1"); /* MSVFilter threshold: must be < F1 to go on */
  double F2                = esl_opt_GetReal(go, "--F2"); /* ViterbiFilter threshold                    */
  double F3                = esl_opt_GetReal(go, "--F3"); /* ForwardFilter threshold                    */
  double Evalue_threshold  = esl_opt_GetReal(go, "-E");	  /* per-seq E-value threshold                  */
  float            usc, vfsc;          /* filter scores                           */
  float            fwdsc;
  float            final_sc;	       /* final bit score                         */
  float            nullsc;             /* null model score                        */
  float            omega  = 1.0f/256.0f;
  float            seqbias;  
  float            seq_score;          /* the per-seq bit score */
  double           P;		       /* P-value of a hit */
  double           Z;                  /* effective # of models, for P-values      */
  double           E;		       /* bound on E-value of a hit  */
  int              Ld;		       /* # of residues in envelopes */
  char             errbuf[eslERRBUFSIZE];
  int              status, hstatus, sstatus;
  int              d;

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) esl_fatal(errbuf);
  if ((status = init_worker_cfg(go, cfg, errbuf)) != eslOK) esl_fatal(errbuf);

  sq = esl_sq_Create();	/* initially in text mode, until first HMM is read. */
 
  output_header(go, cfg);

  if (! esl_opt_IsDefault(go, "-Z")) Z = (double) esl_opt_GetInteger(go, "-Z");
  else if (cfg->hfp->is_pressed)     Z = (double) cfg->hfp->ssi->nprimary;
  else                               Z = 0.0;

  /* We don't know the alphabet yet: cfg->abc == NULL until we read the first HMM.
   */

  while ( (sstatus = esl_sqio_Read(cfg->sqfp, sq)) == eslOK)
    {
      tr  = p7_trace_Create();
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

	  p7_omx_GrowTo(oxf, om->M, 0, sq->n); 
	  p7_oprofile_ReconfigLength(om, sq->n);
	
	  /* First level filter: the MSV filter, multihit with <om> */
	  p7_MSVFilter(sq->dsq, sq->n, om, oxf, &usc);
	  usc = (usc - nullsc) / eslCONST_LOG2;
	  P = esl_gumbel_surv(usc,  om->evparam[p7_MU],  om->evparam[p7_LAMBDA]);
	  if (P > F1) goto FINISH;
	  cfg->n_past_msv++;

	  /* If it passes the MSV filter, read the rest of the profile */
	  p7_oprofile_ReadRest(cfg->hfp, om);

	  /* Second level filter: ViterbiFilter(), multihit with <om> */
	  p7_ViterbiFilter(sq->dsq, sq->n, om, oxf, &vfsc);  
	  vfsc = (vfsc-nullsc) / eslCONST_LOG2;
	  P  = esl_gumbel_surv(vfsc,  om->evparam[p7_MU], om->evparam[p7_LAMBDA]);
	  if (P > F2) goto FINISH;
	  cfg->n_past_vit++;

	  /* Parse it with Forward and obtain its real Forward score. */
	  p7_ForwardParser(sq->dsq, sq->n, om, oxf, &fwdsc);
	  final_sc = (fwdsc-nullsc) / eslCONST_LOG2;
	  P = esl_exp_surv(final_sc,  om->evparam[p7_TAU], om->evparam[p7_LAMBDA]);
	  if (P > F3) goto FINISH;
	  cfg->n_past_fwd++;

	  /* It's for real; do the Backwards parser pass and hand it
	   * off to domain definition logic.
	   */
	  p7_omx_GrowTo(oxb, om->M, 0, sq->n);
	  p7_BackwardParser(sq->dsq, sq->n, om, oxf, oxb, NULL);
	  status = p7_domaindef_ByPosteriorHeuristics(sq, om, oxf, oxb, fwd, bck, cfg->ddef);
	  if      (status == eslERANGE) { fprintf(stderr, "WARNING: posterior decoding failed on %s; skipping it!\n", sq->name); goto FINISH; }
	  else if (status != eslOK)     { fprintf(stderr, "WARNING: domain definition failed due to unknown problem on %s; skipping it!\n", sq->name); goto FINISH; }

	  /* What's the per-seq score, and is it significant enough to be reported? */
	  /* Figure out the sum of null2 corrections to be added to the null score */
	  seqbias = esl_vec_FSum(cfg->ddef->n2sc, sq->n+1);
	  seqbias = p7_FLogsum(0.0, log(omega) + seqbias);
	  seq_score =  (fwdsc - (nullsc + seqbias)) / eslCONST_LOG2;
	  P         =  esl_exp_surv (seq_score,  om->evparam[p7_TAU], om->evparam[p7_LAMBDA]);
	  E = (Z > 0.0) ? P*Z : P * (double) cfg->nseq;
	    
	  if (E <= Evalue_threshold) /* P*current nseq is only a bound on the E-value */
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

	      hit->pre_score  = ( fwdsc-nullsc) / eslCONST_LOG2;
	      hit->pre_pvalue = esl_exp_surv (hit->pre_score,  om->evparam[p7_TAU], om->evparam[p7_LAMBDA]);

	      hit->score      = seq_score;
	      hit->pvalue     = P;
	      hit->sortkey    = -log(P);

	      /* Calculate the "reconstruction score": estimated per-sequence score as sum of the domains */
	      hit->sum_score = 0.0f;
	      Ld             = 0;
	      final_sc       = 0.0;
	      seqbias        = 0.0f;
	      
	      for (d = 0; d < cfg->ddef->ndom; d++) 
		{
		  if (cfg->ddef->dcl[d].envsc - cfg->ddef->dcl[d].domcorrection > 0.0) 
		    {
		      final_sc       += cfg->ddef->dcl[d].envsc;
		      seqbias        += cfg->ddef->dcl[d].domcorrection;
		      Ld             += cfg->ddef->dcl[d].jenv  - cfg->ddef->dcl[d].ienv + 1;
		    }
		}
	      final_sc       += (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3)); 
	      seqbias = p7_FLogsum(0.0, log(omega) + seqbias);
	      hit->sum_score  = (final_sc - (nullsc + seqbias)) / eslCONST_LOG2;
	      hit->sum_pvalue = esl_exp_surv (hit->sum_score,  om->evparam[p7_TAU], om->evparam[p7_LAMBDA]);

	      /* A special case: let sum_score override the seq_score when it's better, and it includes at least 1 domain */
	      if (Ld > 0 && hit->sum_score > seq_score)
		{
		  hit->score     = hit->sum_score;
		  hit->pvalue    = hit->sum_pvalue;
		  hit->sortkey   = -log(hit->sum_pvalue);
		  hit->pre_score = ( final_sc-nullsc) / eslCONST_LOG2;
		  hit->pre_pvalue = esl_exp_surv (hit->pre_score,  om->evparam[p7_TAU], om->evparam[p7_LAMBDA]);
		}

	      /* Transfer the domain coordinates and alignment display to the hit list */
	      hit->dcl       = cfg->ddef->dcl;
	      cfg->ddef->dcl = NULL;
	      for (d = 0; d < hit->ndom; d++)
		{
		  Ld = hit->dcl[d].jenv - hit->dcl[d].ienv + 1;
		  hit->dcl[d].bitscore = hit->dcl[d].envsc + (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3)); 
		  seqbias              = p7_FLogsum(0.0, log(omega) + hit->dcl[d].domcorrection);
		  hit->dcl[d].bitscore = (hit->dcl[d].bitscore - (nullsc + seqbias)) / eslCONST_LOG2;
		  hit->dcl[d].pvalue   = esl_exp_surv (hit->dcl[d].bitscore,  om->evparam[p7_TAU], om->evparam[p7_LAMBDA]);
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
      output_per_model_hitlist (go, cfg, cfg->hitlist);  fprintf(cfg->ofp, "\n");
      output_per_domain_hitlist(go, cfg, cfg->hitlist);  fprintf(cfg->ofp, "\n");
      output_alignments        (go, cfg, cfg->hitlist);  fprintf(cfg->ofp, "\n");
      fprintf(cfg->ofp, "//\n");

      p7_omx_Destroy(fwd);
      p7_omx_Destroy(bck);
      p7_omx_Destroy(oxf);
      p7_omx_Destroy(oxb);
      p7_trace_Destroy(tr);
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
  
  fprintf(cfg->ofp, "target HMM database:   %s\n", cfg->hmmfile);
  fprintf(cfg->ofp, "query sequence file:   %s\n", cfg->seqfile);
  if (! esl_opt_IsDefault(go, "-o"))     fprintf(cfg->ofp, "output directed to file:    %s\n",    esl_opt_GetString(go, "-o"));
  if (! esl_opt_IsDefault(go, "-E"))     fprintf(cfg->ofp, "sequence E-value threshold: <= %g\n", esl_opt_GetReal(go, "-E"));
  if (! esl_opt_IsDefault(go, "--F1"))   fprintf(cfg->ofp, "MSV filter P threshold:     <= %g\n", esl_opt_GetReal(go, "--F1"));
  if (! esl_opt_IsDefault(go, "--F2"))   fprintf(cfg->ofp, "Vit filter P threshold:     <= %g\n", esl_opt_GetReal(go, "--F2"));
  if (! esl_opt_IsDefault(go, "--F3"))   fprintf(cfg->ofp, "Fwd filter P threshold:     <= %g\n", esl_opt_GetReal(go, "--F3"));
  if (! esl_opt_IsDefault(go, "--seed")) fprintf(cfg->ofp, "Random generator seed:      %d\n",    esl_opt_GetInteger(go, "--seed"));
  fprintf(cfg->ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}

static int
output_per_model_hitlist(ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist)
{
  int    h;
  int    textw            = esl_opt_GetInteger(go, "--textw");
  int    namew            = ESL_MAX(8, p7_tophits_GetMaxNameLength(hitlist));
  int    descw            = textw - namew - 24;         /* 7.1f score, 10.2g Eval, 3d ndom, 4 spaces betw 5 fields: 24 char  */
  double Evalue_threshold = esl_opt_GetReal(go, "-E");
  double E;
  
  fprintf(cfg->ofp, "Scores for complete sequence (score includes all domains):\n");

  fprintf(cfg->ofp, "%7s %10s %7s %7s %10s %3s %5s %3s %3s %3s %3s %-*s %s\n", "Score", "E-value", "  bias ", " Sum sc", " Sum E-val", " N ", " exp ", "reg", "clu", " ov", "env", namew, "Sequence", "Description");
  fprintf(cfg->ofp, "%7s %10s %7s %7s %10s %3s %5s %3s %3s %3s %3s %-*s %s\n", "-----", "-------", "-------", "-------", "----------", "---", "-----", "---", "---", "---", "---", namew, "--------", "-----------");

  for (h = 0; h < hitlist->N; h++)
    {
      E = (double) cfg->nseq * hitlist->hit[h]->pvalue;

      if (E <= Evalue_threshold) 
	{
	  fprintf(cfg->ofp, "%7.1f %10.2g %7.1f %7.1f %10.2g %3d %5.1f %3d %3d %3d %3d %-*s %-.*s\n",
		  hitlist->hit[h]->score,
		  hitlist->hit[h]->pvalue * (double) cfg->nseq,
		  hitlist->hit[h]->pre_score - hitlist->hit[h]->score, /* bias correction */
		  hitlist->hit[h]->sum_score,
		  hitlist->hit[h]->sum_pvalue * (double) cfg->nseq,
		  hitlist->hit[h]->ndom,
		  hitlist->hit[h]->nexpected,
		  hitlist->hit[h]->nregions,
		  hitlist->hit[h]->nclustered,
		  hitlist->hit[h]->noverlaps,
		  hitlist->hit[h]->nenvelopes,
		  namew, hitlist->hit[h]->name,
		  descw, hitlist->hit[h]->desc);
	}
    }
  return eslOK;
}


static int
output_per_domain_hitlist(ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist)
{
  int h, d;
  int namew               = ESL_MAX(8, p7_tophits_GetMaxNameLength(hitlist));
  double Evalue_threshold = esl_opt_GetReal(go, "-E");
  double seqE;
  
  fprintf(cfg->ofp, "Scores for individual domains (scored as if only this domain occurred):\n");

  fprintf(cfg->ofp, "%-*s %10s %7s %6s %6s %7s %10s %10s\n", namew, "Sequence", " seq E-val", " Domain", "sqfrom", "sqto",     "score", " ind E-val", " dom E-val");
  fprintf(cfg->ofp, "%-*s %10s %7s %6s %6s %7s %10s %10s\n", namew, "--------", "----------", "-------", "------", "------", "-------", "----------", "----------");

  for (h = 0; h < hitlist->N; h++)
    {
      seqE = hitlist->hit[h]->pvalue * (double) cfg->nseq;
      if (seqE <= Evalue_threshold) 
	{
	  for (d = 0; d < hitlist->hit[h]->ndom; d++)
	    {
	      fprintf(cfg->ofp, "%-*s %10.2g %3d/%-3d %6d %6d %7.1f %10.2g %10.2g\n", 
		      namew, hitlist->hit[h]->name, 
		      seqE,
		      d+1, hitlist->hit[h]->ndom,
		      hitlist->hit[h]->dcl[d].ienv,
		      hitlist->hit[h]->dcl[d].jenv,
		      hitlist->hit[h]->dcl[d].bitscore,
		      hitlist->hit[h]->dcl[d].pvalue * (double) cfg->nseq,  /* E-value if this were the only domain in the seq          */
		      hitlist->hit[h]->dcl[d].pvalue * (double) (h+1));     /* E-value conditional on all seqs to this point being true */
	    }
	}
    }
  return eslOK;
}

static int
output_alignments(ESL_GETOPTS *go, struct cfg_s *cfg, P7_TOPHITS *hitlist)
{
  int    h,d;
  double seqE;
  int    textw            = esl_opt_GetInteger(go, "--textw");

  fprintf(cfg->ofp, "Alignments of each domain found in top-scoring sequences:\n");

  for (h = 0; h < hitlist->N; h++)
    {
      seqE = hitlist->hit[h]->pvalue * (double) cfg->nseq;
    
      for (d = 0; d < hitlist->hit[h]->ndom; d++)
	{
	  fprintf(cfg->ofp, "%s: domain %d of %d, from %d to %d: score %.1f, E = %10.1g\n",
		  hitlist->hit[h]->name,
		  d+1, hitlist->hit[h]->ndom,
		  hitlist->hit[h]->dcl[d].ienv,
		  hitlist->hit[h]->dcl[d].jenv,
		  hitlist->hit[h]->dcl[d].bitscore,
		  hitlist->hit[h]->dcl[d].pvalue * (double) (h+1));  
	  p7_alidisplay_Print(cfg->ofp, hitlist->hit[h]->dcl[d].ad, 40, textw);
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
  fprintf(cfg->ofp, "Query sequences:    %ld  (%ld residues)\n", cfg->nseq,    cfg->nres);
  fprintf(cfg->ofp, "Target HMM(s):      %ld  (%ld nodes)\n",    cfg->nmodels, cfg->nnodes);

  fprintf(cfg->ofp, "Passed MSV filter: %ld  (%.3f; expected %.3f)\n", 
	  cfg->n_past_msv, (double) cfg->n_past_msv / (double) cfg->nmodels, esl_opt_GetReal(go, "--F1"));
  fprintf(cfg->ofp, "Passed Vit filter: %ld  (%.4f; expected %.4f)\n",   
	  cfg->n_past_vit, (double) cfg->n_past_vit / (double) cfg->nmodels, esl_opt_GetReal(go, "--F2"));
  fprintf(cfg->ofp, "Passed Fwd filter: %ld  (%.2g; expected %.2g)\n",         
	  cfg->n_past_fwd, (double) cfg->n_past_fwd / (double) cfg->nmodels, esl_opt_GetReal(go, "--F3"));

  fprintf(cfg->ofp, "Mc/sec:            %.2f\n", 
	  (double) cfg->nres * (double) cfg->nnodes / (w->user * 1.0e6));
  esl_stopwatch_Display(cfg->ofp, w, "# CPU time: ");

  return eslOK;
}
