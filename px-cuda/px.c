
#include "easel.h"
#include "esl_dsqdata.h"

#include "hmmer.h"

#include <pthread.h>

/* CREW   (struct crew_s)
 * Shared data amongst the threads.
 */
typedef struct crew_s {
  int               nworkers;
  struct worker_s **uw;
} CREW;


/* WORKER   (struct worker_s)
 * Data private to each worker.
 */
typedef struct worker_s {
  int          idx;    // Which worker I am; 0..nworkers-1
  pthread_t    thrid;  // My thread id
  CREW        *crew;   

  ESL_DSQDATA *dd;   // reference to reader. (threadsafe by design)
  
  P7_BG       *bg;   // clone
  P7_PROFILE  *gm;   // clone. TODO: shadow these clones instead; cloning is wasteful in memory
  P7_OPROFILE *om;   // clone

  P7_ENGINE   *eng;  // independent

  char errbuf[eslERRBUFSIZE];
  int status;
} WORKER;

static CREW *crew_Create (ESL_DSQDATA *dd, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, int n);
static int   crew_Start  (CREW *crew);
static int   crew_Finish (CREW *crew);
static void  crew_Destroy(CREW *crew);

static void *search_thread(void *p);

static CREW *
crew_Create(ESL_DSQDATA *dd, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, int n)
{
  CREW    *crew = NULL;
  int      u;
  int      status;

  ESL_ALLOC(crew, sizeof(CREW));
  crew->nworkers  = n;
  crew->uw        = NULL;

  ESL_ALLOC(crew->uw, sizeof(WORKER *) * n);
  for (u = 0; u < n; u++) crew->uw[u] = NULL;

  for (u = 0; u < n; u++)
    {
      ESL_ALLOC(crew->uw[u], sizeof(WORKER));
      crew->uw[u]->dd  = dd;   // reference
      crew->uw[u]->bg  = NULL;
      crew->uw[u]->gm  = NULL;
      crew->uw[u]->om  = NULL;
      crew->uw[u]->eng = NULL;
      
      if (u == 0) {
	crew->uw[u]->bg = bg;
	crew->uw[u]->gm = gm;
	crew->uw[u]->om = om;
      } else {
	crew->uw[u]->bg = p7_bg_Clone(bg);
	crew->uw[u]->gm = p7_profile_Clone(gm);
	crew->uw[u]->om = p7_oprofile_Create(gm->M, gm->abc);
	p7_oprofile_Convert (gm, crew->uw[u]->om);
      }
      crew->uw[u]->eng = p7_engine_Create(gm->abc, NULL, NULL, 200, 400);
    }
  return crew;

 ERROR:
  crew_Destroy(crew);
  return NULL;
}

static int
crew_Start(CREW *crew)
{
  int u;

  for (u = 0; u < crew->nworkers; u++)
    {
      crew->uw[u]->idx  = u;
      crew->uw[u]->crew = crew;
      pthread_create(&(crew->uw[u]->thrid), NULL, search_thread, crew->uw[u]);
    }
  return eslOK;
}


static int
crew_Finish(CREW *crew)
{
  int u;

  for (u = 0; u < crew->nworkers; u++)
    pthread_join( crew->uw[u]->thrid, NULL);
  return eslOK;
}

static void
crew_Destroy(CREW *crew)
{
  int u;
  if (!crew) return;
  
  if (crew->uw) {
    for (u = 0; u < crew->nworkers; u++)
      {
	if (crew->uw[u]) {
	  if (u>0) p7_bg_Destroy(crew->uw[u]->bg);
	  if (u>0) p7_profile_Destroy(crew->uw[u]->gm);
	  if (u>0) p7_oprofile_Destroy(crew->uw[u]->om);
	  p7_engine_Destroy(crew->uw[u]->eng);
	  free(crew->uw[u]);
	}
      }
    free(crew->uw);
  }
  free(crew);
}


static void *
search_thread(void *p)
{
  ESL_STOPWATCH     *w   = esl_stopwatch_Create();
  WORKER            *uw  = (WORKER *) p;  
  ESL_DSQDATA       *dd  = uw->dd;
  P7_PROFILE        *gm  = uw->gm;
  P7_OPROFILE       *om  = uw->om;
  P7_BG             *bg  = uw->bg;
  P7_ENGINE         *eng = uw->eng;
  ESL_DSQDATA_CHUNK *chu = NULL;
  int      i;
  int      status;

  //  esl_stopwatch_Start(w);

  while (( status = esl_dsqdata_Read(dd, &chu)) == eslOK)  
    {
      //     esl_stopwatch_Stop(w);
      //      printf("thread %d:  read time: %.5f sec\n", uw->idx, esl_stopwatch_GetElapsed(w));
      //esl_stopwatch_Start(w);

      for (i = 0; i < chu->N; i++)
	{
	  p7_bg_SetLength(bg, (int) chu->L[i]);            // TODO: remove need for cast
	  p7_oprofile_ReconfigLength(om, (int) chu->L[i]); //         (ditto)
	  
	  status = p7_engine_Overthruster(eng, chu->dsq[i], (int) chu->L[i], om, bg);  
	  if (status == eslFAIL) { 
	    p7_engine_Reuse(eng);
	    continue;
	  }

	  p7_profile_SetLength(gm, (int) chu->L[i]);
	  status = p7_engine_Main(eng, chu->dsq[i], (int) chu->L[i], gm); 

	  p7_engine_Reuse(eng);
	}
      //      esl_stopwatch_Stop(w);
      //printf("thread %d: chunk time: %.5f sec\n", uw->idx, esl_stopwatch_GetElapsed(w));
      //esl_stopwatch_Start(w);

      esl_dsqdata_Recycle(dd, chu);
    }
  
  esl_stopwatch_Stop(w);
  //printf("thread %d:   eof time: %.5f sec\n", uw->idx, esl_stopwatch_GetElapsed(w));
  //esl_stopwatch_Destroy(w);

  uw->errbuf[0] = '\0';
  uw->status    = eslOK;
  return (void *) uw;
}


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { "-h",        eslARG_NONE,  FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  { "-n",        eslARG_INT,     "1",  NULL, NULL,   NULL,  NULL, NULL, "set number of threads to <n>",          0 },
  { "-s",        eslARG_INT,     "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",         0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "px, the first parallel tests of H4";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_BG          *bg      = NULL;
  P7_HMM         *hmm     = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  ESL_DSQDATA    *dd      = NULL;
  CREW           *crew    = NULL;
  int             ncore   = esl_opt_GetInteger(go, "-n");
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create (hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config   (gm, hmm, bg);
  p7_oprofile_Convert (gm, om);

  p7_bg_SetFilter(bg, om->M, om->compo);

  /* Open sequence database */
  status = esl_dsqdata_Open(&abc, seqfile, ncore, &dd);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open dsqdata files:\n  %s",    dd->errbuf);
  else if (status == eslEFORMAT)   p7_Fail("Format problem in dsqdata files:\n  %s", dd->errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error in opening dsqdata (code %d)", status);

  /* Create the work crew */
  crew = crew_Create(dd, gm, om, bg, ncore);
  
  crew_Start(crew);
  crew_Finish(crew);

  crew_Destroy(crew);
  esl_dsqdata_Close(dd);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  exit(0);
}





