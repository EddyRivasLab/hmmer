/* hmmevolve: table of time dependent parameters from a eHMMER profile HMM.
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "e2_config.h"
#include "e1_rate.h"
#include "evohmmer.h"
#include "ratematrix.h"


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",          eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL,  "show brief help on version and usage",                       1 },
  { "-o",          eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,    NULL,  "send  output table to file <f>, not stdout",                 1 },
/* evomodel */
  { "--evomodel",  eslARG_STRING,"AGAX", NULL, NULL,      NULL,      NULL,    NULL,  "evolutionary model used",                                    2 },
  { "--betainf",   eslARG_REAL,  "0.69", NULL, "x>=0",    NULL,      NULL,    NULL,  "betainf = ldI/muI = beta at time infinity (if ldI<muI)",     2 },
  { "--fixtime",   eslARG_REAL,    NULL, NULL, "x>=0",    NULL,      NULL,    NULL,  "TRUE: use a fix time for the evolutionary models of a pair", 2 },
  { "--noevo",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--mx,--mxfile", "do not evolve",                                       2 },

/* options controlling wich state of model to evolve */
  { "-M",          eslARG_INT,    FALSE, NULL, "n>=0",    NULL,      NULL,    NULL,  "State to evolve",                                            2 },
/* options controlling evolutionary time range */
  { "--mintime",   eslARG_REAL,   "0.0",  NULL, "x>=0.0", NULL,      NULL,    NULL, "minimum time",                                                3 },
  { "--maxtime",   eslARG_REAL,  "20.0",  NULL, "x>=0.0", NULL,      NULL,    NULL, "maximum time",                                                3 },
/* other options */
  { "--seed",      eslARG_INT,      "0", NULL, "n>=0",    NULL,      NULL,    NULL, "set RNG seed to <n>",                                         4 },
  { "--tol",       eslARG_REAL,  "1e-2", NULL, NULL,      NULL,      NULL,    NULL, "tolerance",                                                   4 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile (single)>";
static char banner[] = "evolve transition probabilities of a eHMMER profile HMM";

static void cmdline_failure(char *argv0, char *format, ...);
static void cmdline_help(char *argv0, ESL_GETOPTS *go);

static void evolve_hmm(ESL_GETOPTS *go, FILE *ofp, HMMRATE *hmmrate, P7_HMM *hmm, ESL_ALPHABET *abc, P7_BG *bg, double tmin, double tmax, int state, char *errbuf);


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go         = NULL;             /* command line processing                 */
  ESL_ALPHABET    *abc        = NULL;             /* sequence alphabet                       */
  ESL_RANDOMNESS  *r          = NULL;             /* source of randomness                    */
  char            *hmmfile    = NULL;             /* file to read HMM(s) from                */
  P7_HMMFILE      *hfp        = NULL;             /* open hmmfile                            */
  P7_HMM          *hmm        = NULL;             /* HMM to emit from                        */
  FILE            *ofp        = NULL;	          /* output stream                           */
  HMMRATE         *hmmrate    = NULL;
  P7_BG           *bg         = NULL;
  double           tmin;
  double           tmax;
  int              nhmms      = 0;
  int              state      = -1;               /* state to plot, all by default */
  int              status;	      
  char             errbuf[eslERRBUFSIZE];

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in configuration: %s\n",       go->errbuf);
  if (esl_opt_GetBoolean(go, "-h"))                    cmdline_help   (argv[0], go);      
  if (esl_opt_ArgNumber(go) != 1)                      cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        

  if ((hmmfile = esl_opt_GetArg(go, 1)) == NULL)       cmdline_failure(argv[0], "Failed to get <hmmfile> on cmdline: %s\n", go->errbuf);

  if ( esl_opt_IsOn(go, "-o") ) {
    if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "-o"));
  } else ofp = stdout;

  /* the evolutionary model */
  ESL_ALLOC(hmmrate, sizeof(HMMRATE));
  hmmrate->emR = NULL;
  hmmrate->S   = NULL;
  hmmrate->evomodel = e1_rate_Evomodel(esl_opt_GetString(go, "--evomodel"));
  hmmrate->betainf  = esl_opt_GetReal(go, "--betainf");
  hmmrate->fixtime  = esl_opt_IsOn(go, "--fixtime")? esl_opt_GetReal(go, "--fixtime") : -1.0;
  hmmrate->tol      = esl_opt_GetReal(go, "--tol");
  hmmrate->statfp   = stdout;
  
  // time range
  tmin  = esl_opt_GetReal(go, "--mintime");
  tmax  = esl_opt_GetReal(go, "--maxtime");
  state = esl_opt_IsOn(go, "-M")? esl_opt_GetInteger(go, "-M") : -1;

  r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));

  status = p7_hmmfile_Open(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",                       status, hmmfile, errbuf);  

  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslEOF)
    {
      if      (status == eslEFORMAT)    esl_fatal("Bad file format in HMM file %s:\n%s\n",          hfp->fname, hfp->errbuf);
      else if (status == eslEINCOMPAT)  esl_fatal("HMM in %s is not in the expected %s alphabet\n", hfp->fname, esl_abc_DecodeType(abc->type));
      else if (status != eslOK)         esl_fatal("Unexpected error in reading HMMs from %s\n",     hfp->fname);
      nhmms++;
      if (nhmms > 1) break; // only one HMM

      // null model
      bg = p7_bg_Create(abc);

      /* read the substitution matrix */
      if (!esl_opt_IsOn(go, "--noevo")) {
	double *f = NULL;
	int     k;
	
	ESL_ALLOC(f, sizeof(double) * abc->K);
	for (k = 0; k < abc->K; k ++) {
	  f[k] = (double)bg->f[k];
	}
	
	hmmrate->emR = ratematrix_emrate_Create(abc, 1);
	ratematrix_emrate_Set("BLOSUM62", NULL, f, &hmmrate->emR[0], TRUE, hmmrate->tol, errbuf, FALSE);
	
	if (f) free(f);
      }

      evolve_hmm(go, ofp, hmmrate, hmm, abc, bg, tmin, tmax, state, errbuf);

      esl_alphabet_Destroy(abc); abc = NULL;
      p7_bg_Destroy(bg);         bg  = NULL;
      p7_hmm_Destroy(hmm);       hmm = NULL;
    }
  if (nhmms == 0) esl_fatal("Empty HMM file %s? No HMM data found.\n"); 

  if (hmmrate) {
    if (hmmrate->emR) ratematrix_emrate_Destroy(hmmrate->emR, 1);
    free(hmmrate);
  }
  if (esl_opt_IsOn(go, "-o")) { fclose(ofp); }
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  p7_hmmfile_Close(hfp);
  return eslOK;

 ERROR:
    if (esl_opt_IsOn(go, "-o")) { fclose(ofp); }
    if (r) esl_randomness_Destroy(r);
    if (abc) esl_alphabet_Destroy(abc);
    if (bg) p7_bg_Destroy(bg);
    if (hmm) p7_hmm_Destroy(hmm);
    if (go) esl_getopts_Destroy(go);
    if (hfp) p7_hmmfile_Close(hfp);
    return eslFAIL;
}


static void
cmdline_failure(char *argv0, char *format, ...) 
{
  va_list argp;
  printf("\nERROR: ");
  va_start(argp, format);
  vfprintf(stdout, format, argp);
  va_end(argp);
  esl_usage(stdout, argv0, usage);
  printf("\nTo see more help on available options, do %s -h\n\n", argv0);
  exit(1);
}

static void
cmdline_help(char *argv0, ESL_GETOPTS *go) 
{
  p7_banner (stdout, argv0, banner);
  esl_usage (stdout, argv0, usage);
  puts("\nCommon options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  puts("\nOptions controlling which state to evolve:");
  esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
  puts("\nOptions controlling evolutionary time range:");
  esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
  puts("\nOther options::");
  esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
  exit(0);
}

static void
evolve_hmm(ESL_GETOPTS *go, FILE *ofp, HMMRATE *hmmrate, P7_HMM *hmm, ESL_ALPHABET *abc, P7_BG *bg, double tmin, double tmax, int state, char *errbuf)
{
  P7_RATE *R  = NULL;
  double   time;
  double   tdelta = 1e-3;
  int      m;

  if (state > hmm->M) { printf("hmm state is %d but the hmm has %d states\n", state, hmm->M); return; }
      
  // make the hmm Rate
  if (p7_RateConstruct(hmm, bg, hmmrate, &R, errbuf, FALSE) != eslOK) esl_fatal("HMM RateConstruct() failed");

  fprintf(ofp, "# state  time  tMM     tMI     tMO     tII      tIM      tDD     tDM\n", state, time);
  while (time <= tmax) {
    // evolve the hmm
    p7_EvolveFromRate(hmmrate->statfp, hmm, R, bg, time, errbuf, FALSE);

    for (m = 0; m < hmm->M; m ++) {
      
      if (m == state || state < 0) {
	
	// print the transition probabilities
	fprintf(ofp, "%d %f %f %f %f %f %f %f %f\n",
		m, time,
		hmm->t[m][p7H_MM], hmm->t[m][p7H_MI], hmm->t[m][p7H_MD], 
		hmm->t[m][p7H_II], hmm->t[m][p7H_IM],
		hmm->t[m][p7H_DD], hmm->t[m][p7H_DM]);
      }
    }
    
    time += tdelta;
  }
  
  
  p7_RateDestroy(R);
  
  return;
}

