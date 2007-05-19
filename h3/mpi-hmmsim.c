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

#include "mpi.h"

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

#define ALGORITHMS "--fwd,--viterbi,--island"
enum algchoice_e { DO_VITERBI, DO_FORWARD, DO_ISLAND };

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",   1 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0",     NULL,      NULL,    NULL, "length of random target seqs",           1 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0",     NULL,      NULL,    NULL, "number of random target seqs",           1 },
  { "--viterbi", eslARG_NONE,  "TRUE", NULL, NULL, ALGORITHMS,     NULL,    NULL, "Score seqs with the Viterbi algorithm",  2 },
  { "--fwd",     eslARG_NONE,   FALSE, NULL, NULL, ALGORITHMS,     NULL,    NULL, "Score seqs with the Forward algorithm",  2 },
  { "--island",  eslARG_NONE,   FALSE, NULL, NULL, ALGORITHMS,     NULL,    NULL, "Score seqs with the Island algorithm",   2 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "mpi-hmmsim [-options] <hmmfile>";

static void mpi_worker(int N, int L, enum algchoice_e algorithm_choice, int my_rank, double mintail, double maxtail, double tailstep);

int
main(int argc, char **argv)
{
  int status;
  ESL_GETOPTS     *go	   = NULL;      /* command line processing                   */
  char            *hmmfile = NULL;	/* HMM file to read                          */
  P7_HMMFILE      *hfp     = NULL;      /* open HMM file                             */
  ESL_ALPHABET    *abc     = NULL;      /* alphabet to use                           */
  P7_HMM          *hmm     = NULL;      /* query HMM                                 */
  P7_BG           *bg      = NULL;      /* null1 model                               */
  P7_HMM         **hmmlist = NULL;      /* queue of HMMs being worked on, 1..nproc-1 */
  int              L;			/* length of random sequences                */
  int              N;			/* number of random sequences                */
  enum algchoice_e algorithm_choice;
  int              optset;
  int              my_rank;
  int              nproc;
  int              have_work;
  int              nproc_working;
  int              wi;
  double          *results;	/* mu, lambda of one or more fit, possibly for a range of tail probabilities */
  MPI_Status       mstatus;
  double           mintail, maxtail, tailstep, tailp;
  int              ntailsettings;
  int              i;

  /* Process command line options.
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    puts(usage);
    puts("\ngeneral options are:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 2 = indentation; 80=textwidth*/
    puts("\nalternative scoring algorithms :");
    esl_opt_DisplayHelp(stdout, go, 2, 2, 80); /* 2 = indentation; 80=textwidth*/
    exit(0);
  }
  L = esl_opt_GetInteger(go, "-L");
  N = esl_opt_GetInteger(go, "-N");
  if      (esl_opt_GetBoolean(go, "--viterbi")) algorithm_choice = DO_VITERBI;
  else if (esl_opt_GetBoolean(go, "--fwd"))     algorithm_choice = DO_FORWARD;
  else if (esl_opt_GetBoolean(go, "--island"))  algorithm_choice = DO_ISLAND;
  if (esl_opt_ArgNumber(go) != 1) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    exit(1);
  }
  hmmfile = esl_opt_GetArg(go, eslARG_INFILE, NULL);
  if (hmmfile == NULL) esl_fatal("Failed to read <hmmfile> argument from command line.");
  esl_getopts_Destroy(go);

  mintail = 0.001;
  maxtail = 1.0;
  tailstep = exp(0.25 * log(2.)); /* that's the fourth root of 2: so, 4 steps per doubling */

  /* Initialize MPI, figure out who we are, and whether we're running
   * this show (proc 0) or working in it (procs >0)
   */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  /* If we're a worker, go to that section of the code
   */
  if (my_rank > 0)
    {  
      mpi_worker(N, L, algorithm_choice, my_rank, mintail, maxtail, tailstep);
      MPI_Finalize();
      return 0;
    }
  /* Else, we're the master process 0; continue... */


  /* Initializations, including opening the HMM file
   */
  abc = esl_alphabet_Create(eslAMINO);

  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if (status == eslENOTFOUND) esl_fatal("Failed to open hmm file %s for reading.\n", hmmfile);
  else if (status != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", hmmfile);

  bg = p7_bg_Create(abc);
  ESL_ALLOC(hmmlist, sizeof(P7_HMM *) * nproc);
  for (wi = 0; wi < nproc; wi++) hmmlist[wi] = NULL;

  /* How many different values of the tail probability will we use, to test different fits?
   */
  for (ntailsettings = 0, tailp = mintail; tailp <= maxtail; tailp *= tailstep) ntailsettings++;
  if (ntailsettings < 1) ntailsettings = 1;
  if (algorithm_choice == DO_VITERBI) ntailsettings = 1;
  ESL_ALLOC(results, sizeof(double) * 3 * ntailsettings);

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
      if ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) == eslOK) 
	{
	  hmm->gm = p7_profile_Create(hmm->M, abc);
	  hmm->bg = bg;
	  p7_ProfileConfig(hmm, hmm->gm, p7_UNILOCAL);
	}
      else if (status == eslEOF)       have_work = FALSE;
      else if (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
      else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s", hmmfile);
      else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets", hmmfile);
      else esl_fatal("Unexpected error in reading HMMs from %s", hmmfile);


      /* If we have work but no free workers, or we have no work but workers
       * are still working, then wait for a result to return from any worker.
       */
      if ( (have_work && nproc_working == nproc-1) || (! have_work && nproc_working > 0))
	{
	  MPI_Recv(results, 3*ntailsettings, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mstatus);
	  wi = mstatus.MPI_SOURCE;

	  for (i = 0; i < ntailsettings; i++)
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


  free(hmmlist);
  p7_bg_Destroy(bg);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  MPI_Finalize();
  return 0;

 ERROR:
  if (hmmlist != NULL) free(hmmlist);
  if (bg      != NULL) p7_bg_Destroy(bg);
  if (hfp     != NULL) p7_hmmfile_Close(hfp);
  if (abc     != NULL) esl_alphabet_Destroy(abc);
  return status;
}


/* mpi_worker()
 * The main control for an MPI worker process.
 */
static void
mpi_worker(int N, int L, enum algchoice_e algorithm_choice, int my_rank, double mintail, double maxtail, double tailstep)
{
  int status;
  ESL_RANDOMNESS *r   = esl_randomness_CreateTimeseeded();
  ESL_ALPHABET   *abc = esl_alphabet_Create(eslAMINO);    
  ESL_HISTOGRAM  *h   = NULL;
  P7_BG          *bg  = p7_bg_Create(abc);
  P7_GMX         *gx  = p7_gmx_Create(200, L);
  P7_PROFILE     *gm  = NULL;
  ESL_DSQ        *dsq = NULL;
  int             i;  
  int             sc;
  int             nullsc;
  float           bitscore;
  double *xv;
  int     n;			/* number of data points */
  double  mu, lambda;
  double *results = NULL;	/* mu, lambda fits, possibly for a range of tail probabilities: up to 3*ntailsettings values */
  int     ntailsettings;
  double  tailp;

  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));

  /* How many different values of the tail probability will we use, to test different fits?
   */
  for (ntailsettings = 0, tailp = mintail; tailp <= maxtail; tailp *= tailstep)  ntailsettings++;
  if (ntailsettings < 1) ntailsettings = 1;
  ESL_ALLOC(results, sizeof(double) * 3 * ntailsettings);

  /* Main loop */
  while (p7_profile_MPIRecv(abc, bg, &gm) == eslOK) 
    {
      if (algorithm_choice == DO_ISLAND) h = esl_histogram_Create(-50., 50., 0.5); /* bin the island scores; there's too many of them.  */
      else                               h = esl_histogram_CreateFull(-50.5, 50.5, 0.5);  
      p7_gmx_GrowTo(gx, gm->M, L);

      /* Collect a score histogram from N random sequences of length L  */
      for (i = 0; i < N; i++)
	{
	  esl_rnd_xfIID(r, gm->bg->f, gm->abc->K, L, dsq);
	  
	  if (algorithm_choice == DO_VITERBI) 
	    {
	      p7_Viterbi(dsq, L, gm, gx, NULL, &sc);
	      p7_bg_NullOne(gm->bg, dsq, L, &nullsc);
	      bitscore = p7_SILO2Bitscore(sc - nullsc);
	      esl_histogram_Add(h, bitscore);
	    }
	  else if (algorithm_choice == DO_FORWARD) 
	    {
	      p7_Forward(dsq, L, gm, gx, &sc);
	      p7_bg_NullOne(gm->bg, dsq, L, &nullsc);
	      bitscore = p7_SILO2Bitscore(sc - nullsc);
	      esl_histogram_Add(h, bitscore);
	    }
	  else if (algorithm_choice == DO_ISLAND) 
	    { /* island needs mx for at least 4 rows */
	      p7_island_Viterbi(dsq, L, gm, gx, h);
	    }
	}

      /* Fit the histogram  */
      if (algorithm_choice == DO_VITERBI)
	{
	  esl_histogram_GetData(h, &xv, &n);
	  esl_gumbel_FitComplete(xv, n, &mu, &lambda);

	  results[0]    = 1.0;
	  results[1]    = mu;
	  results[2]    = lambda;
	  ntailsettings = 1;
	}
      else if (algorithm_choice == DO_FORWARD)
	{
	  for (i = 0, tailp = mintail; tailp <= maxtail; tailp *= tailstep)
	    {
	      esl_histogram_GetTailByMass(h, tailp, &xv, &n, NULL);
	      esl_exp_FitComplete(xv, n, &mu, &lambda);
	      results[i*3]   = tailp;
	      results[i*3+1] = mu;
	      results[i*3+2] = lambda;
	      i++;
	    }
	}
      else if (algorithm_choice == DO_ISLAND)
	{
	  /* For island, we do a binned fit, because of the data volume. */
	  double actual_mass;

	  for (i = 0, tailp = mintail; tailp <= maxtail; tailp *= tailstep)
	    {
	      esl_histogram_SetTailByMass(h, tailp, &actual_mass);
	      esl_exp_FitCompleteBinned(h, &mu, &lambda);
	      results[i*3]   = tailp;
	      results[i*3+1] = mu;
	      results[i*3+2] = lambda;
	      i++;
	    }
	}
      
      /* Send the result back to the master   
       * (note, in a single send, so results don't interleave)
       */
      MPI_Send(results, 3*ntailsettings, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

      esl_histogram_Destroy(h);
      p7_profile_Destroy(gm);
    }
  
  free(results);
  free(dsq);
  p7_gmx_Destroy(gx);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return;

 ERROR:
  if (results != NULL) free(results);
  if (dsq != NULL) free(dsq);
  p7_gmx_Destroy(gx);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return;
}


