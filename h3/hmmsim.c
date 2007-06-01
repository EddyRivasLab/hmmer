/* main() for scoring profile HMMs against simulated sequences
 * 
 * Example:
 *  ./hmmbuild Pfam /misc/data0/databases/Pfam/Pfam-A.seed                
 *  qsub -N testrun -j y -R y -b y -cwd -V -pe lam-mpi-tight 32 'mpirun C ./mpi-hmmsim -N 10000 ./Pfam > foo.out'
 * 
 * SRE, Fri Apr 20 14:56:26 2007 [Janelia]
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
#include "esl_histogram.h"
#include "esl_vectorops.h"
#include "esl_dmatrix.h"
#include "esl_ratematrix.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"

#include "hmmer.h"

#define ALGORITHMS "--fwd,--viterbi,--island,--hybrid"
#define STYLES     "--fs,--sw,--ls,--s"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, NULL, "show brief help on version and usage",            1 },
  { "-L",        eslARG_INT,    "100", NULL, "n>0",     NULL,  NULL, NULL, "length of random target seqs",                    1 },
  { "-N",        eslARG_INT,   "1000", NULL, "n>0",     NULL,  NULL, NULL, "number of random target seqs",                    1 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, NULL, "run as an MPI parallel program",                  1 },
#endif
  { "-o",        eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,"--mpi", "output P(S>x) histogram to <f> in xy format",     2 },
  { "--xfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,"--mpi", "output bitscores as binary double vector to <f>", 2 },

  { "--fs",      eslARG_NONE,"default",NULL, NULL,    STYLES,  NULL, NULL, "multihit local alignment",                    3 },
  { "--sw",      eslARG_NONE,   FALSE, NULL, NULL,    STYLES,  NULL, NULL, "unihit local alignment",                      3 },
  { "--ls",      eslARG_NONE,   FALSE, NULL, NULL,    STYLES,  NULL, NULL, "multihit glocal alignment",                   3 },
  { "--s",       eslARG_NONE,   FALSE, NULL, NULL,    STYLES,  NULL, NULL, "unihit glocal alignment",                     3 },

  { "--viterbi", eslARG_NONE,"default",NULL, NULL, ALGORITHMS, NULL, NULL, "Score seqs with the Viterbi algorithm",       4 },
  { "--fwd",     eslARG_NONE,   FALSE, NULL, NULL, ALGORITHMS, NULL, NULL, "Score seqs with the Forward algorithm",       4 },
  { "--island",  eslARG_NONE,   FALSE, NULL, NULL, ALGORITHMS, NULL, NULL, "Score seqs with the Island algorithm",        4 },
  { "--hybrid",  eslARG_NONE,   FALSE, NULL, NULL, ALGORITHMS, NULL, NULL, "Score seqs with the Hybrid algorithm",        4 },

  { "--tmin",    eslARG_REAL,  "0.02", NULL, NULL,      NULL,  NULL, NULL, "Set lower bound tail mass for fwd,island",    5 },
  { "--tmax",    eslARG_REAL,  "0.02", NULL, NULL,      NULL,  NULL, NULL, "Set lower bound tail mass for fwd,island",    5 },
  { "--tstep",   eslARG_REAL,  "0.02", NULL, NULL,      NULL,  NULL, NULL, "Set additive step size for tmin...tmax range",5 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct cfg_s {
  char           *hmmfile;
  ESL_RANDOMNESS *r;
  ESL_ALPHABET   *abc;
  P7_BG          *bg;
  P7_HMMFILE     *hfp;
  int             ntailsettings;
  int             my_rank;
  int             nproc;
  int             do_mpi;
};

static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "collect profile HMM score distributions on random sequences";

static void serial_master   (ESL_GETOPTS *go, struct cfg_s *cfg, double *results);
#ifdef HAVE_MPI
static void mpi_master      (ESL_GETOPTS *go, struct cfg_s *cfg, double *results);
static void mpi_worker      (ESL_GETOPTS *go, struct cfg_s *cfg, double *results);
#endif 
static void process_workunit(ESL_GETOPTS *go, struct cfg_s *cfg, P7_PROFILE *gm, double *results);


int
main(int argc, char **argv)
{
  int status;
  ESL_GETOPTS     *go	   = NULL;      /* command line processing                   */
  struct cfg_s     cfg;
  double          *results = NULL;	/* mu, lambda of one or more fit, possibly for a range of tail probabilities */

  /* Process command line options.
   */
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
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\noutput options (only in serial mode, for single HMM input):");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nalternative alignment styles :");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\nalternative scoring algorithms :");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\ncontrolling range of fitted tail masses :");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  /* Initialize configuration shared across all kinds of masters
   * and workers in this .c file.
   */
  cfg.hmmfile  = esl_opt_GetArg(go, 1);
  if (cfg.hmmfile == NULL) esl_fatal("Failed to read <hmmfile> argument from command line.");
  cfg.r        = esl_randomness_CreateTimeseeded();
  cfg.abc      = esl_alphabet_Create(eslAMINO);
  cfg.bg       = p7_bg_Create(cfg.abc);
  cfg.hfp      = NULL;
  cfg.my_rank  = 0;
  cfg.nproc    = 0;
  cfg.do_mpi   = FALSE;

  cfg.ntailsettings = (int) (1 + (esl_opt_GetReal(go, "--tmax") - esl_opt_GetReal(go, "--tmin")) / esl_opt_GetReal(go, "--tstep"));
  if (cfg.ntailsettings < 1) cfg.ntailsettings = 1;
  if (esl_opt_GetBoolean(go, "--viterbi") || esl_opt_GetBoolean(go, "--hybrid")) cfg.ntailsettings = 1;

  ESL_ALLOC(results, sizeof(double) * 6 * cfg.ntailsettings);

  /* Initialize MPI, figure out who we are, and whether we're running
   * this show (proc 0) or working in it (procs >0)
   */
#ifdef HAVE_MPI
  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      cfg.do_mpi = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if (cfg.my_rank > 0)	/* we're one of the MPI workers */
	{  
	  mpi_worker(go, &cfg, results);
	}
      else 			/* we're the MPI master */
	{
	  status = p7_hmmfile_Open(cfg.hmmfile, NULL, &(cfg.hfp));
	  if (status == eslENOTFOUND) esl_fatal("Failed to open hmm file %s for reading.\n",  cfg.hmmfile);
	  else if (status != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", cfg.hmmfile);

	  mpi_master   (go, &cfg, results);
	}

      MPI_Finalize();		/* both workers and masters reach this line */
    }
  else
#endif /*HAVE_MPI*/
    {				/*  we're the serial master */
      status = p7_hmmfile_Open(cfg.hmmfile, NULL, &(cfg.hfp));
      if (status == eslENOTFOUND) esl_fatal("Failed to open hmm file %s for reading.\n",  cfg.hmmfile);
      else if (status != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", cfg.hmmfile);

      serial_master(go, &cfg, results);
    }      

  /* normal execution falls through this goto target:
   * both error and normal cases execute the same cleanup code.
   */
  status = eslOK;
 ERROR:
  if (results != NULL) free(results);
  if (cfg.hfp != NULL) p7_hmmfile_Close(cfg.hfp);
  if (cfg.bg  != NULL) p7_bg_Destroy(cfg.bg);
  if (cfg.abc != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.r   != NULL) esl_randomness_Destroy(cfg.r);
  if (go      != NULL) esl_getopts_Destroy(go);
  return status;
}

static void
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg, double *results)
{
  P7_HMM     *hmm = NULL;      /* query HMM                                 */
  P7_PROFILE *gm  = NULL;
  int     status;
  int     i;

  while ((status = p7_hmmfile_Read(cfg->hfp, &(cfg->abc), &hmm)) != eslEOF) 
    {
      if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", cfg->hmmfile);
      else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s",             cfg->hmmfile);
      else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets",   cfg->hmmfile);
      else if (status != eslOK)        esl_fatal("Unexpected error in reading HMMs from %s",   cfg->hmmfile);

      gm = p7_profile_Create(hmm->M, cfg->abc);
      hmm->gm = (P7_PROFILE *) gm;
      if      (esl_opt_GetBoolean(go, "--fs"))  p7_ProfileConfig(hmm, cfg->bg, gm, p7_LOCAL);
      else if (esl_opt_GetBoolean(go, "--sw"))  p7_ProfileConfig(hmm, cfg->bg, gm, p7_UNILOCAL);
      else if (esl_opt_GetBoolean(go, "--ls"))  p7_ProfileConfig(hmm, cfg->bg, gm, p7_GLOCAL);
      else if (esl_opt_GetBoolean(go, "--s"))   p7_ProfileConfig(hmm, cfg->bg, gm, p7_UNIGLOCAL);

      p7_ReconfigLength(gm,      esl_opt_GetInteger(go, "-L"));
      p7_bg_SetLength  (cfg->bg, esl_opt_GetInteger(go, "-L"));

      process_workunit(go, cfg, gm, results);

      for (i = 0; i < cfg->ntailsettings; i++)
	{	 /* <name> <tail prob>  <mu> <lambda> <E@10>  <mu w/ fixed lambda> <E@10> */
	  printf("%-20s  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
		 hmm->name, results[6*i], results[6*i+1], results[6*i+2], results[6*i+3], results[6*i+4], results[6*i+5]);
	  fflush(stdout);
	}

      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);      
    }
}


#ifdef HAVE_MPI
static void
mpi_master(ESL_GETOPTS *go, struct cfg_s *cfg, double *results)
{
  P7_HMM          *hmm     = NULL;      /* query HMM                                 */
  P7_HMM         **hmmlist = NULL;      /* queue of HMMs being worked on, 1..nproc-1 */
  P7_PROFILE      *gm      = NULL; 
  int              have_work;
  int              nproc_working;
  int              wi;
  int              status;
  MPI_Status       mstatus;
  int              i;

  ESL_ALLOC(hmmlist, sizeof(P7_HMM *) * cfg->nproc);
  for (wi = 0; wi < cfg->nproc; wi++) hmmlist[wi] = NULL;

  /* My design pattern for data parallelization in a master/worker model:
   * three phases: 
   *  1. load workers;
   *  2. recv result/send work loop;
   *  3. collect remaining results
   * but implemented in a single while loop to avoid redundancy.
   */
  have_work     = TRUE;
  nproc_working = 0;
  wi            = 1;
  while (have_work || nproc_working)
    {
      /* Get next work unit. */
      if ((status = p7_hmmfile_Read(cfg->hfp, &(cfg->abc), &hmm)) == eslOK) 
	{
	  gm = p7_profile_Create(hmm->M, cfg->abc);
	  hmm->gm = (P7_PROFILE *) gm;
	  if      (esl_opt_GetBoolean(go, "--fs"))  p7_ProfileConfig(hmm, cfg->bg, gm, p7_LOCAL);
	  else if (esl_opt_GetBoolean(go, "--sw"))  p7_ProfileConfig(hmm, cfg->bg, gm, p7_UNILOCAL);
	  else if (esl_opt_GetBoolean(go, "--ls"))  p7_ProfileConfig(hmm, cfg->bg, gm, p7_GLOCAL);
	  else if (esl_opt_GetBoolean(go, "--s"))   p7_ProfileConfig(hmm, cfg->bg, gm, p7_UNIGLOCAL);

	  p7_ReconfigLength(gm,      esl_opt_GetInteger(go, "-L"));
	  p7_bg_SetLength  (cfg->bg, esl_opt_GetInteger(go, "-L"));
	}
      else if (status == eslEOF)       have_work = FALSE;
      else if (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", cfg->hmmfile);
      else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s",             cfg->hmmfile);
      else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets",   cfg->hmmfile);
      else esl_fatal("Unexpected error in reading HMMs from %s", cfg->hmmfile);


      /* If we have work but no free workers, or we have no work but workers
       * are still working, then wait for a result to return from any worker.
       */
      if ( (have_work && nproc_working == cfg->nproc-1) || (! have_work && nproc_working > 0))
	{
	  MPI_Recv(results, 6*cfg->ntailsettings, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mstatus);
	  wi = mstatus.MPI_SOURCE;

	  for (i = 0; i < cfg->ntailsettings; i++)
	    {	 /* <name> <tail prob>  <mu> <lambda> <E@10>  <mu w/ fixed lambda> <E@10> */
	      printf("%-20s  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
		     hmmlist[wi]->name, results[6*i], results[6*i+1], results[6*i+2], results[6*i+3], results[6*i+4], results[6*i+5]);
	      fflush(stdout);
	    }

	  nproc_working--;
	  p7_hmm_Destroy(hmmlist[wi]);
	  hmmlist[wi] = NULL;
	}
	
      /* If we have work, assign it to a free worker;
       * else, terminate the free worker.
       */
      if (have_work) 
	{
	  p7_profile_MPISend(gm, wi);
	  hmmlist[wi] = hmm;
	  wi++;
	  nproc_working++;

	  p7_profile_Destroy(gm);
	}
    }

  /* Tell all the workers (1..nproc-1) to shut down by sending them a NULL workunit. */
  for (wi = 1; wi < cfg->nproc; wi++) p7_profile_MPISend(NULL, wi);	

  /* normal execution falls through here; both normal and error cases execute same cleanup code */
 ERROR:
  if (hmmlist != NULL) free(hmmlist);
  return;
}


/* mpi_worker()
 * The main control for an MPI worker process.
 */
static void
mpi_worker(ESL_GETOPTS *go, struct cfg_s *cfg, double *results)
{
  P7_PROFILE     *gm      = NULL;

  while (p7_profile_MPIRecv(cfg->abc, cfg->bg, &gm) == eslOK) 
    {
      process_workunit(go, cfg, gm, results);
      MPI_Send(results, 6*cfg->ntailsettings, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

      p7_profile_Destroy(gm);
    }
  return;
}
#endif /*HAVE_MPI*/


/* A work unit consists of one profile, <gm>, with <gm->abc>, <gm->bg> both set.
 * The result is the <results> array, which contains offset/mu/lambda/mu2 fits for one or
 * more tail masses.
 */
static void
process_workunit(ESL_GETOPTS *go, struct cfg_s *cfg, P7_PROFILE *gm, double *results)
{
  int             L   = esl_opt_GetInteger(go, "-L");
  int             N   = esl_opt_GetInteger(go, "-N");
  P7_GMX         *gx  = p7_gmx_Create(gm->M, L);
  ESL_HISTOGRAM  *h   = NULL;
  ESL_DSQ        *dsq = NULL;
  FILE           *histfp = NULL; 
  int             i;
  int             sc;
  int             nullsc;
  float           bitscore;
  double         *xv;
  int             n;			/* number of data points */
  double          mu, lambda;
  double          tmin  = esl_opt_GetReal(go, "--tmin");
  double          tmax  = esl_opt_GetReal(go, "--tmax");
  double          tstep = esl_opt_GetReal(go, "--tstep");
  double          tailp;
  int             status;
  double          param[2];
  double          x10;

  if (esl_opt_GetBoolean(go, "--island")) h = esl_histogram_Create(-50., 50., 0.5); /* bin the island scores; there's too many of them.  */
  else                                    h = esl_histogram_CreateFull(-50.5, 50.5, 0.5);  
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  

  /* Output options. */
  if (esl_opt_GetString(go, "-o") != NULL) 
    {
      char  *histfile = esl_opt_GetString(go, "-o");
      
      histfp   = fopen(histfile, "w");
      if (histfp == NULL) p7_Die("Failed to open output file %s\n", histfile);
    }

  /* Collect a score histogram from N random sequences of length L  */
  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(cfg->r, gm->bg->f, gm->abc->K, L, dsq);
	  
      if (esl_opt_GetBoolean(go, "--viterbi"))
	{
	  p7_GViterbi(dsq, L, gm, gx, &sc);
	  p7_bg_NullOne(gm->bg, dsq, L, &nullsc);
	  bitscore = p7_SILO2Bitscore(sc - nullsc);
	  esl_histogram_Add(h, bitscore);
	}
      else if (esl_opt_GetBoolean(go, "--fwd"))
	{
	  p7_GForward(dsq, L, gm, gx, &sc);
	  p7_bg_NullOne(gm->bg, dsq, L, &nullsc);
	  bitscore = p7_SILO2Bitscore(sc - nullsc);
	  esl_histogram_Add(h, bitscore);
	}
      else if (esl_opt_GetBoolean(go, "--island"))
	{ /* island needs mx for at least 4 rows */
	  p7_island_Viterbi(dsq, L, gm, gx, h);
	}
      else if (esl_opt_GetBoolean(go, "--hybrid"))
	{ 
	  p7_GHybrid(dsq, L, gm, gx, NULL, &sc);
	  p7_bg_NullOne(gm->bg, dsq, L, &nullsc);
	  bitscore = p7_SILO2Bitscore(sc - nullsc);
	  esl_histogram_Add(h, bitscore);
	}
    }


  /* Fit the histogram  */
  if (esl_opt_GetBoolean(go, "--viterbi") || esl_opt_GetBoolean(go, "--hybrid"))
    {
      esl_histogram_GetData(h, &xv, &n);
      esl_histogram_GetRank(h, 10, &x10);
      esl_gumbel_FitComplete(xv, n, &mu, &lambda);
      
      results[0]    = 1.0;
      results[1]    = mu;
      results[2]    = lambda;
      results[3]    = n * esl_gumbel_surv(x10, mu, lambda); /* E-value of the 10th ranked score should be ~10. */

      if (histfp != NULL) {
	param[0]      = mu;
	param[1]      = lambda;
	esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, &param);
	esl_histogram_PlotSurvival(histfp, h);
      }

      lambda = 0.693;
      esl_gumbel_FitCompleteLoc(xv, n, 0.693, &mu);
      results[4]    = mu;
      results[5]    = n * esl_gumbel_surv(x10, mu, lambda);

      if (histfp != NULL) {
	param[0]      = mu;
	param[1]      = lambda;
	esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, &param);
	esl_histogram_PlotSurvival(histfp, h);
      }

    }

  else if (esl_opt_GetBoolean(go, "--fwd"))
    {
      for (i = 0, tailp = tmin; tailp <= tmax+1e-7; tailp += tstep)
	{
	  esl_histogram_GetTailByMass(h, tailp, &xv, &n, NULL);
	  esl_histogram_GetRank(h, 10, &x10);
	  esl_exp_FitComplete(xv, n, &mu, &lambda);

	  results[i*6]   = tailp;
	  results[i*6+1] = mu;
	  results[i*6+2] = lambda;
	  results[i*6+3] = esl_opt_GetInteger(go, "-N") * tailp * esl_exp_surv(x10, mu, lambda);
	  results[i*6+4] = mu;
	  results[i*6+5] = esl_opt_GetInteger(go, "-N") * tailp * esl_exp_surv(x10, mu, 0.693);
	  i++;

	  if (histfp != NULL) {
	    param[0] = mu;
	    param[1] = lambda;
	    esl_histogram_SetExpectedTail(h, mu, tailp, &esl_exp_generic_cdf, &param);
	    esl_histogram_PlotSurvival(histfp, h);

	    param[0] = mu;
	    param[1] = 0.693; 
	    esl_histogram_SetExpectedTail(h, mu, tailp, &esl_exp_generic_cdf, &param);
	    esl_histogram_PlotSurvival(histfp, h);
	  }
	}
    }
  else if (esl_opt_GetBoolean(go, "--island"))
    {
      /* For island, we do a binned fit, because of the data volume. */
      double actual_mass;

      for (i = 0, tailp = tmin; tailp <= tmax+1e-7; tailp += tstep)
	{
	  esl_histogram_SetTailByMass(h, tailp, &actual_mass);
	  esl_exp_FitCompleteBinned(h, &mu, &lambda);
	  results[i*6]   = tailp;
	  results[i*6+1] = mu;
	  results[i*6+2] = lambda;
	  results[i*6+3] = 0.;
	  results[i*6+4] = mu;
	  results[i*6+5] = 0.;
	  i++;

	  if (histfp != NULL) {
	    param[0] = mu;
	    param[1] = lambda;
	    esl_histogram_SetExpectedTail(h, mu, actual_mass, &esl_exp_generic_cdf, &param);
	    esl_histogram_PlotSurvival(histfp, h);

	    param[0] = mu;
	    param[1] = 0.693; 
	    esl_histogram_SetExpectedTail(h, mu, actual_mass, &esl_exp_generic_cdf, &param);
	    esl_histogram_PlotSurvival(histfp, h);
	  }
	}
    }


  /* fallthrough: both normal, error cases execute same cleanup code */
 ERROR:
  if (histfp != NULL) fclose(histfp);
  if (h      != NULL) esl_histogram_Destroy(h);
  if (gx     != NULL) p7_gmx_Destroy(gx);
  if (dsq    != NULL) free(dsq);
  return;
}
