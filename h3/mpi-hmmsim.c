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
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",   1 },
  { "-L",        eslARG_INT,    "100", NULL, "n>0",     NULL,      NULL,    NULL, "length of random target seqs",           1 },
  { "-N",        eslARG_INT,   "1000", NULL, "n>0",     NULL,      NULL,    NULL, "number of random target seqs",           1 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "run as an MPI parallel program",         1 },
#endif
  { "--fs",      eslARG_NONE,"default",NULL, NULL,    STYLES,      NULL,    NULL, "multihit local alignment",               2 },
  { "--sw",      eslARG_NONE,   FALSE, NULL, NULL,    STYLES,      NULL,    NULL, "unihit local alignment",                 2 },
  { "--ls",      eslARG_NONE,   FALSE, NULL, NULL,    STYLES,      NULL,    NULL, "multihit glocal alignment",              2 },
  { "--s",       eslARG_NONE,   FALSE, NULL, NULL,    STYLES,      NULL,    NULL, "unihit glocal alignment",                2 },
  { "--viterbi", eslARG_NONE,"default",NULL, NULL, ALGORITHMS,     NULL,    NULL, "Score seqs with the Viterbi algorithm",  3 },
  { "--fwd",     eslARG_NONE,   FALSE, NULL, NULL, ALGORITHMS,     NULL,    NULL, "Score seqs with the Forward algorithm",  3 },
  { "--island",  eslARG_NONE,   FALSE, NULL, NULL, ALGORITHMS,     NULL,    NULL, "Score seqs with the Island algorithm",   3 },
  { "--hybrid",  eslARG_NONE,   FALSE, NULL, NULL, ALGORITHMS,     NULL,    NULL, "Score seqs with the Hybrid algorithm",   3 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct cfg_s {
  char           *hmmfile;
  ESL_RANDOMNESS *r;
  ESL_ALPHABET   *abc;
  P7_BG          *bg;
  P7_HMMFILE     *hfp;
  double          mintail;
  double          maxtail;
  double          tailstep;
  int             ntailsettings;
  int             my_rank;
  int             nproc;
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
  double           tailp;

  /* Process command line options.
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || 
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      esl_fatal("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 2 = indentation; 80=textwidth*/
      puts("\nalternative alignment styles :");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); /* 2 = indentation; 80=textwidth*/
      puts("\nalternative scoring algorithms :");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); /* 2 = indentation; 80=textwidth*/
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
  cfg.mintail  = 0.001;
  cfg.maxtail  = 1.0;
  cfg.tailstep = exp(0.25 * log(2.)); /* that's the fourth root of 2: so, 4 steps per doubling */
  cfg.my_rank  = 0;
  cfg.nproc    = 0;

  for (cfg.ntailsettings = 0, tailp = cfg.mintail; tailp <= cfg.maxtail; tailp *= cfg.tailstep) cfg.ntailsettings++;
  if (cfg.ntailsettings < 1) cfg.ntailsettings = 1;
  if (esl_opt_GetBoolean(go, "--viterbi") || esl_opt_GetBoolean(go, "--hybrid")) cfg.ntailsettings = 1;

  ESL_ALLOC(results, sizeof(double) * 3 * cfg.ntailsettings);

  /* Initialize MPI, figure out who we are, and whether we're running
   * this show (proc 0) or working in it (procs >0)
   */
#ifdef HAVE_MPI
  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
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
      MPI_Finalize();
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
  P7_HMM *hmm = NULL;      /* query HMM                                 */
  int     status;
  int     i;

  while ((status = p7_hmmfile_Read(cfg->hfp, &(cfg->abc), &hmm)) != eslEOF) 
    {
      if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", cfg->hmmfile);
      else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s",             cfg->hmmfile);
      else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets",   cfg->hmmfile);
      else if (status != eslOK)        esl_fatal("Unexpected error in reading HMMs from %s",   cfg->hmmfile);

      hmm->gm = p7_profile_Create(hmm->M, cfg->abc);
      if      (esl_opt_GetBoolean(go, "--fs"))  p7_ProfileConfig(hmm, cfg->bg, hmm->gm, p7_LOCAL);
      else if (esl_opt_GetBoolean(go, "--sw"))  p7_ProfileConfig(hmm, cfg->bg, hmm->gm, p7_UNILOCAL);
      else if (esl_opt_GetBoolean(go, "--ls"))  p7_ProfileConfig(hmm, cfg->bg, hmm->gm, p7_GLOCAL);
      else if (esl_opt_GetBoolean(go, "--s"))   p7_ProfileConfig(hmm, cfg->bg, hmm->gm, p7_UNIGLOCAL);

      process_workunit(go, cfg, hmm->gm, results);

      for (i = 0; i < cfg->ntailsettings; i++)
	{			                /* <name>   <tail prob>    <mu>             <lambda>  */
	  printf("%-20s  %8.4f %8.4f  %8.4f\n", hmm->name, results[3*i], results[3*i+1], results[3*i+2]);
	  fflush(stdout);
	}

      p7_profile_Destroy(hmm->gm);
      p7_hmm_Destroy(hmm);      
    }
}


#ifdef HAVE_MPI
static void
mpi_master(ESL_GETOPTS *go, struct cfg_s *cfg, double *results)
{
  P7_HMM          *hmm     = NULL;      /* query HMM                                 */
  P7_HMM         **hmmlist = NULL;      /* queue of HMMs being worked on, 1..nproc-1 */
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
	  hmm->gm = p7_profile_Create(hmm->M, cfg->abc);
	  if      (esl_opt_GetBoolean(go, "--fs"))  p7_ProfileConfig(hmm, cfg->bg, hmm->gm, p7_LOCAL);
	  else if (esl_opt_GetBoolean(go, "--sw"))  p7_ProfileConfig(hmm, cfg->bg, hmm->gm, p7_UNILOCAL);
	  else if (esl_opt_GetBoolean(go, "--ls"))  p7_ProfileConfig(hmm, cfg->bg, hmm->gm, p7_GLOCAL);
	  else if (esl_opt_GetBoolean(go, "--s"))   p7_ProfileConfig(hmm, cfg->bg, hmm->gm, p7_UNIGLOCAL);
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
	  MPI_Recv(results, 3*cfg->ntailsettings, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mstatus);
	  wi = mstatus.MPI_SOURCE;

	  for (i = 0; i < cfg->ntailsettings; i++)
	    {			                     /* <name>          <tail prob>    <mu>             <lambda>  */
	      printf("%-20s  %8.4f %8.4f  %8.4f\n", hmmlist[wi]->name, results[3*i], results[3*i+1], results[3*i+2]);
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
	  p7_profile_MPISend(hmm->gm, wi);
	  hmmlist[wi] = hmm;
	  wi++;
	  nproc_working++;

	  p7_profile_Destroy(hmm->gm);
	}
      else p7_profile_MPISend(NULL, wi);	
    }

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
      MPI_Send(results, 3*cfg->ntailsettings, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

      p7_profile_Destroy(gm);
    }
  return;
}
#endif /*HAVE_MPI*/


/* A work unit consists of one profile, <gm>, with <gm->abc>, <gm->bg> both set.
 * The result is the <results> array, which contains offset/mu/lambda fits for one or
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
  int             i;
  int             sc;
  int             nullsc;
  float           bitscore;
  double         *xv;
  int             n;			/* number of data points */
  double          mu, lambda;
  double          tailp;
  int             status;

  if (esl_opt_GetBoolean(go, "--island")) h = esl_histogram_Create(-50., 50., 0.5); /* bin the island scores; there's too many of them.  */
  else                                    h = esl_histogram_CreateFull(-50.5, 50.5, 0.5);  
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));

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
  if (esl_opt_GetBoolean(go, "--viterbi"))
    {
      esl_histogram_GetData(h, &xv, &n);
      esl_gumbel_FitComplete(xv, n, &mu, &lambda);
      
      results[0]    = 1.0;
      results[1]    = mu;
      results[2]    = lambda;
    }
  else if (esl_opt_GetBoolean(go, "--fwd"))
    {
      for (i = 0, tailp = cfg->mintail; tailp <= cfg->maxtail; tailp *= cfg->tailstep)
	{
	  esl_histogram_GetTailByMass(h, tailp, &xv, &n, NULL);
	  esl_exp_FitComplete(xv, n, &mu, &lambda);
	  results[i*3]   = tailp;
	  results[i*3+1] = mu;
	  results[i*3+2] = lambda;
	  i++;
	}
    }
  else if (esl_opt_GetBoolean(go, "--island"))
    {
      /* For island, we do a binned fit, because of the data volume. */
      double actual_mass;

      for (i = 0, tailp = cfg->mintail; tailp <= cfg->maxtail; tailp *= cfg->tailstep)
	{
	  esl_histogram_SetTailByMass(h, tailp, &actual_mass);
	  esl_exp_FitCompleteBinned(h, &mu, &lambda);
	  results[i*3]   = tailp;
	  results[i*3+1] = mu;
	  results[i*3+2] = lambda;
	  i++;
	}
    }
  else if (esl_opt_GetBoolean(go, "--hybrid"))
    {
      esl_histogram_GetData(h, &xv, &n);
      esl_gumbel_FitComplete(xv, n, &mu, &lambda);
      
      results[0]    = 1.0;
      results[1]    = mu;
      results[2]    = lambda;
    }

  /* fallthrough: both normal, error cases execute same cleanup code */
 ERROR:
  if (h   != NULL) esl_histogram_Destroy(h);
  if (gx  != NULL) p7_gmx_Destroy(gx);
  if (dsq != NULL) free(dsq);
  return;
}
