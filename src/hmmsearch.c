/* hmmsearch: search profile HMM(s) against a sequence database.
 * 
 * SRE, Thu Dec 20 07:07:25 2007 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "esl_mpi.h"
#endif 

#include "hmmer.h"

#define RNGOPTS "--Rdet,--Rseed,-Rarb"                         /* Exclusive options for controlling run-to-run variation      */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles  reqs   incomp  help   docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "show brief help on version and usage",                         1 },
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL, "direct output to file <f>, not stdout",                        1 },
  { "-A",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL, "save multiple alignment of all hits to file <s>",              1 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",        "report sequences <= this E-value threshold in output",         2 },
  { "-T",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",        "report sequences >= this score threshold in output",           2 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                                "set # of comparisons done, for E-value calculation",           2 },
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",        "report domains <= this E-value threshold in output",           2 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",        "report domains >= this score cutoff in output",                2 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                                "set # of significant seqs, for domain E-value calculation",    2 },
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "-E,-T,--domE,--domT",               "use profile's GA gathering cutoffs to set -T, --domT",         2 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "-E,-T,--domE,--domT",               "use profile's NC noise cutoffs to set -T, --domT",             2 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "-E,-T,--domE,--domT",               "use profile's TC trusted cutoffs to set -T, --domT",           2 },
  /* Control of inclusion thresholds */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include sequences <= this E-value threshold in output ali",    3 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include sequences >= this score threshold in output ali",      3 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include domains <= this E-value threshold in output ali",      3 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include domains >= this score threshold in output ali",        3 },
  { "--inc_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  "-A",  "--incE,--incT,--incdomE,--incdomT", "use profile's GA gathering cutoffs to set --incT, --incdomT",  3 },
  { "--inc_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  "-A",  "--incE,--incT,--incdomE,--incdomT", "use profile's NC noise cutoffs to set --incT, --incdomT",      3 },
  { "--inc_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  "-A",  "--incE,--incT,--incdomE,--incdomT", "use profile's TC trusted cutoffs to set --incT, --incdomT",    3 },
  /* Control of acceleration pipeline */
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",      4 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,      NULL,  NULL, "--max",          "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             4 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,      NULL,  NULL, "--max",          "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             4 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,      NULL,  NULL, "--max",          "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             4 },
  { "--nobias",     eslARG_NONE,   NULL,  NULL, NULL,      NULL,  NULL, "--max",          "turn off composition bias filter",                             4 },
  { "--nonull2",    eslARG_NONE,   NULL,  NULL, NULL,      NULL,  NULL,    NULL,          "turn off biased composition score corrections",                4 },
/* Control of run-to-run variation in RNG */
  { "--Rdet",       eslARG_NONE,"default",NULL, NULL,   RNGOPTS,  NULL,    NULL,          "reseed RNG to minimize run-to-run stochastic variation",       5 },
  { "--Rseed",       eslARG_INT,    NULL, NULL, NULL,   RNGOPTS,  NULL,    NULL,          "reseed RNG with fixed seed",                                   5 },
  { "--Rarb",       eslARG_NONE,    NULL, NULL, NULL,   RNGOPTS,  NULL,    NULL,          "seed RNG arbitrarily; allow run-to-run stochastic variation",  5 },
/* Other options */
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",  NULL,  NULL,  "--notextw",     "set max width of ASCII text output lines",                     6 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,      NULL,  NULL,  "--textw",       "unlimit ASCII text output line width",                         6 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[options] <query hmmfile> <target seqfile>";
static char banner[] = "search profile HMM(s) against a sequence database";



static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_hmmfile, char **ret_seqfile)
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

      puts("\nOptions controlling reporting thresholds:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      puts("\nOptions controlling significance (inclusion) thresholds:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 

      puts("\nOptions controlling acceleration heuristics:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      puts("\nOptions controlling run-to-run variation due to random number generation:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      puts("\nOther expert options:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                  != 2)     { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((*ret_hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <hmmfile> argument on command line"); goto ERROR; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <seqfile> argument on command line"); goto ERROR; }
  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere most common options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);  
}

static int
output_header(FILE *ofp, ESL_GETOPTS *go, char *hmmfile, char *seqfile)
{
  p7_banner(ofp, go->argv[0], banner);
  
  fprintf(ofp, "# query HMM file:                  %s\n", hmmfile);
  fprintf(ofp, "# target sequence database:        %s\n", seqfile);
  if (! esl_opt_IsDefault(go, "-o"))          fprintf(ofp, "# output directed to file:         %s\n",      esl_opt_GetString(go, "-o"));
  if (! esl_opt_IsDefault(go, "-A"))          fprintf(ofp, "# MSA of all hits saved to file:   %s\n",      esl_opt_GetString(go, "-A"));
  if (! esl_opt_IsDefault(go, "-E"))          fprintf(ofp, "# sequence E-value threshold:   <= %g\n",      esl_opt_GetReal(go, "-E"));
  if (! esl_opt_IsDefault(go, "-T"))          fprintf(ofp, "# sequence bit score threshold: >= %g\n",      esl_opt_GetReal(go, "-T"));
  if (! esl_opt_IsDefault(go, "-Z"))          fprintf(ofp, "# sequence search space set to:    %.0f\n",    esl_opt_GetReal(go, "-Z"));
  if (! esl_opt_IsDefault(go, "--domE"))      fprintf(ofp, "# domain E-value threshold:     <= %g\n",      esl_opt_GetReal(go, "--domE"));
  if (! esl_opt_IsDefault(go, "--domT"))      fprintf(ofp, "# domain bit score threshold:   >= %g\n",      esl_opt_GetReal(go, "--domT"));
  if (! esl_opt_IsDefault(go, "--domZ"))      fprintf(ofp, "# domain search space set to:      %.0f\n",    esl_opt_GetReal(go, "--domZ"));
  if (! esl_opt_IsDefault(go, "--cut_ga"))    fprintf(ofp, "# set reporting thresholds to:     GA cutoffs\n"); 
  if (! esl_opt_IsDefault(go, "--cut_nc"))    fprintf(ofp, "# set reporting thresholds to:     NC cutoffs\n"); 
  if (! esl_opt_IsDefault(go, "--cut_tc"))    fprintf(ofp, "# set reporting thresholds to:     TC cutoffs\n"); 
  if (! esl_opt_IsDefault(go, "--incE"))      fprintf(ofp, "# seq inclusion E-val thresh:   <= %g\n",      esl_opt_GetReal(go, "--incE"));
  if (! esl_opt_IsDefault(go, "--incT"))      fprintf(ofp, "# seq inclusion score thresh:   >= %g\n",      esl_opt_GetReal(go, "--incT"));
  if (! esl_opt_IsDefault(go, "--incdomE"))   fprintf(ofp, "# dom inclusion E-val thresh:   <= %g\n",      esl_opt_GetReal(go, "--incdomE"));
  if (! esl_opt_IsDefault(go, "--incdomT"))   fprintf(ofp, "# dom inclusion score thresh:   >= %g\n",      esl_opt_GetReal(go, "--incdomT"));
  if (! esl_opt_IsDefault(go, "--inc_ga"))    fprintf(ofp, "# set inclusion thresholds to:     GA cutoffs\n"); 
  if (! esl_opt_IsDefault(go, "--inc_nc"))    fprintf(ofp, "# set inclusion thresholds to:     NC cutoffs\n"); 
  if (! esl_opt_IsDefault(go, "--inc_tc"))    fprintf(ofp, "# set inclusion thresholds to:     TC cutoffs\n"); 
  if (! esl_opt_IsDefault(go, "--max"))       fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n");
  if (! esl_opt_IsDefault(go, "--F1"))        fprintf(ofp, "# MSV filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F1"));
  if (! esl_opt_IsDefault(go, "--F2"))        fprintf(ofp, "# Vit filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F2"));
  if (! esl_opt_IsDefault(go, "--F3"))        fprintf(ofp, "# Fwd filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F3"));
  if (! esl_opt_IsDefault(go, "--nobias"))    fprintf(ofp, "# biased composition HMM filter:   off\n");
  if (! esl_opt_IsDefault(go, "--nonull2"))   fprintf(ofp, "# null2 bias corrections:          off\n");
  if (! esl_opt_IsDefault(go, "--Rdet") )     fprintf(ofp, "# RNG seed (run-to-run variation): reseed deterministically; minimize variation\n");
  if (! esl_opt_IsDefault(go, "--Rseed") )    fprintf(ofp, "# RNG seed (run-to-run variation): reseed to %d\n", esl_opt_GetInteger(go, "--Rseed"));
  if (! esl_opt_IsDefault(go, "--Rarb") )     fprintf(ofp, "# RNG seed (run-to-run variation): one arbitrary seed; allow run-to-run variation\n");
  if (! esl_opt_IsDefault(go, "--textw"))     fprintf(ofp, "# max ASCII text line length:      %d\n",     esl_opt_GetInteger(go, "--textw"));
  if (! esl_opt_IsDefault(go, "--notextw"))   fprintf(ofp, "# max ASCII text line length:      unlimited\n");
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go	   = NULL;               /* command line processing                 */
  FILE            *ofp     = NULL;
  char            *hmmfile = NULL;
  P7_HMMFILE      *hfp     = NULL;
  P7_HMM          *hmm     = NULL;
  char            *dbfile  = NULL;
  int              dbfmt   = eslSQFILE_UNKNOWN;
  ESL_SQFILE      *dbfp    = NULL;
  ESL_SQ          *dbsq    = NULL;
  ESL_ALPHABET    *abc     = NULL;
  ESL_STOPWATCH   *w       = NULL;
  P7_BG           *bg      = NULL;
  int              textw;
  int              status;
  int              hstatus, sstatus;

  /* Initializations */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */
  process_commandline(argc, argv, &go, &hmmfile, &dbfile);    
  w = esl_stopwatch_Create();
  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* Open the target sequence database */
  status = esl_sqfile_Open(dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          dbfile);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            dbfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, dbfile);  

  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open hmm file %s for reading.\n",                      hmmfile);
  else if (status == eslEFORMAT)   p7_Fail("Unrecognized format, trying to open hmm file %s for reading.\n", hmmfile);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening hmm file %s.\n", status,          hmmfile);  

  /* Open the results output file */
  if (esl_opt_IsDefault(go, "-o")) ofp = stdout;
  else {
    ofp = fopen(esl_opt_GetString(go, "-o"), "w");
    if (ofp == NULL) p7_Fail("Failed to open output file %s for writing\n", esl_opt_GetString(go, "-o"));
  }

  /* Outer loop: over each query HMM in <hmmfile>. <abc> is not known 'til first HMM is read. */
  while ((hstatus = p7_hmmfile_Read(hfp, &abc, &hmm)) == eslOK) 
    {
      P7_PIPELINE     *pli     = NULL;
      P7_TOPHITS      *th      = NULL;
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;

      esl_stopwatch_Start(w);

      /* One-time initializations after alphabet <abc> becomes known */
      if (bg == NULL) 
	{
	  output_header(ofp, go, hmmfile, dbfile);
	  bg   = p7_bg_Create(abc);
	  dbsq = esl_sq_CreateDigital(abc);
	}

      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
      p7_oprofile_Convert(gm, om);                  /* <om> is now p7_LOCAL, multihit */

      /* Create processing pipeline and hit list */
      pli = p7_pipeline_Create(go, om->M, 100, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
      th  = p7_tophits_Create(); 

      fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M);
      if (hmm->acc  != NULL) fprintf(ofp, "Accession:   %s\n", hmm->acc);
      if (hmm->desc != NULL) fprintf(ofp, "Description: %s\n", hmm->desc);
      p7_pli_NewModel(pli, om, bg);

      /* Main loop: */
      while ( (sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
	{
	  p7_pli_NewSeq(pli, dbsq);
	  p7_bg_SetLength(bg, dbsq->n);
	  p7_oprofile_ReconfigLength(om, dbsq->n);
	  
	  p7_Pipeline(pli, om, bg, dbsq, th);
	  
	  esl_sq_Reuse(dbsq);
	  p7_pipeline_Reuse(pli);
	}
      if      (sstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s line %" PRId64 "):\n%s\n",
						dbfp->filename, dbfp->linenumber, dbfp->errbuf);     
      else if (sstatus != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
						sstatus, dbfp->filename);
  
      /* Print the results.  */
      p7_tophits_Sort(th);
      p7_tophits_Threshold(th, pli);
      p7_tophits_Targets(ofp, th, pli, bg, textw); fprintf(ofp, "\n\n");
      p7_tophits_Domains(ofp, th, pli, bg, textw); fprintf(ofp, "\n\n");
  
      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, pli, w);
      fprintf(ofp, "//\n");

      /* Output the results in an MSA */
      if (! esl_opt_IsDefault(go, "-A")) {
	FILE    *afp = NULL;
	ESL_MSA *msa = NULL;

	if ((afp = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL)
	  fprintf(ofp, "WARNING: failed to open alignment file %s; skipping the alignment output\n", esl_opt_GetString(go, "-A"));

	if (afp != NULL) {
	  p7_tophits_Alignment(th, abc, NULL, NULL, 0, p7_DEFAULT, &msa);
	  if (textw > 0) esl_msa_Write(afp, msa, eslMSAFILE_STOCKHOLM);
	  else           esl_msa_Write(afp, msa, eslMSAFILE_PFAM);
	  
	  fprintf(ofp, "# Alignment of all hits above threshold saved to: %s\n", esl_opt_GetString(go, "-A"));
	  
	  esl_msa_Destroy(msa);
	  fclose(afp);
	}
      }

      esl_sqfile_Position(dbfp, 0); /* rewind */
      p7_pipeline_Destroy(pli);
      p7_tophits_Destroy(th);
      p7_profile_Destroy(gm);
      p7_oprofile_Destroy(om);
      p7_hmm_Destroy(hmm);
    } /* end outer loop over query HMMs */
  if      (hstatus == eslEOD)       p7_Fail("read failed, HMM file %s may be truncated?", hmmfile);
  else if (hstatus == eslEFORMAT)   p7_Fail("bad file format in HMM file %s",             hmmfile);
  else if (hstatus == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets",   hmmfile);
  else if (hstatus != eslEOF)       p7_Fail("Unexpected error in reading HMMs from %s",   hmmfile);


  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);
  p7_bg_Destroy(bg);
  esl_sq_Destroy(dbsq);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  if (! esl_opt_IsDefault(go, "-o")) fclose(ofp);
  esl_getopts_Destroy(go);
  return eslOK;
}


 


/*****************************************************************
 * @LICENSE@
 *****************************************************************/

