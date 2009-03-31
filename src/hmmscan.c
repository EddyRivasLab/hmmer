/* hmmscan: search sequence(s) against a profile HMM database
 * 
 * SRE, Mon Oct 20 08:28:05 2008 [Janelia]
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

#include "hmmer.h"


static ESL_OPTIONS options[] = {
  /* name           type          default  env  range toggles  reqs   incomp                         help                                                      docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL,  NULL,                                   "show brief help on version and usage",                         1 },
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,   NULL,  NULL,  NULL,                                   "direct output to file <f>, not stdout",                        1 },
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
  /* Control of filter pipeline */
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, "--F1,--F2,--F3",   "Turn all heuristic filters off (less speed, more power)",      4 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,   NULL,  NULL, "--max",            "MSV threshold: promote hits w/ P <= F1",                       4 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,   NULL,  NULL, "--max",            "Vit threshold: promote hits w/ P <= F2",                       4 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,   NULL,  NULL, "--max",            "Fwd threshold: promote hits w/ P <= F3",                       4 },
  { "--nobias",     eslARG_NONE,    NULL, NULL, NULL,   NULL,  NULL, "--max",            "turn off composition bias filter",                             4 },
  { "--nonull2",    eslARG_NONE,    NULL, NULL, NULL,   NULL,  NULL,    NULL,            "turn off biased composition score corrections",                4 },
  /* Other options */
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",    NULL,  NULL,    NULL,         "set RNG seed to <n> (if 0: one-time arbitrary seed)",          6 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",  NULL,  NULL,  "--notextw",    "set max width of ASCII text output lines",                     6 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,      NULL,  NULL,  "--textw",      "unlimit ASCII text output line width",                         6 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmm database> <query seqfile>";
static char banner[] = "search sequence(s) against a profile HMM database";



/* process_commandline()
 * 
 * Processes the commandline, filling in fields in <cfg> and creating and returning
 * an <ESL_GETOPTS> options structure. The help page (hmmsearch -h) is formatted
 * here.
 */
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

      puts("\nOther expert options:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 2)      { puts("Incorrect number of command line arguments.");      goto ERROR; }
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
  
  fprintf(ofp, "# query sequence file:             %s\n", seqfile);
  fprintf(ofp, "# target HMM database:             %s\n", hmmfile);
  if (esl_opt_IsUsed(go, "-o"))          fprintf(ofp, "# output directed to file:         %s\n",   esl_opt_GetString(go, "-o"));
  if (esl_opt_IsUsed(go, "-E"))          fprintf(ofp, "# sequence E-value threshold:   <= %g\n",      esl_opt_GetReal(go, "-E"));
  if (esl_opt_IsUsed(go, "-T"))          fprintf(ofp, "# sequence bit score threshold: >= %g\n",      esl_opt_GetReal(go, "-T"));
  if (esl_opt_IsUsed(go, "-Z"))          fprintf(ofp, "# sequence search space set to:    %.0f\n",    esl_opt_GetReal(go, "-Z"));
  if (esl_opt_IsUsed(go, "--domE"))      fprintf(ofp, "# domain E-value threshold:     <= %g\n",      esl_opt_GetReal(go, "--domE"));
  if (esl_opt_IsUsed(go, "--domT"))      fprintf(ofp, "# domain bit score threshold:   >= %g\n",      esl_opt_GetReal(go, "--domT"));
  if (esl_opt_IsUsed(go, "--domZ"))      fprintf(ofp, "# domain search space set to:      %.0f\n",    esl_opt_GetReal(go, "--domZ"));
  if (esl_opt_IsUsed(go, "--cut_ga"))    fprintf(ofp, "# set reporting thresholds to:     GA cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--cut_nc"))    fprintf(ofp, "# set reporting thresholds to:     NC cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--cut_tc"))    fprintf(ofp, "# set reporting thresholds to:     TC cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--incE"))      fprintf(ofp, "# seq inclusion E-val thresh:   <= %g\n",      esl_opt_GetReal(go, "--incE"));
  if (esl_opt_IsUsed(go, "--incT"))      fprintf(ofp, "# seq inclusion score thresh:   >= %g\n",      esl_opt_GetReal(go, "--incT"));
  if (esl_opt_IsUsed(go, "--incdomE"))   fprintf(ofp, "# dom inclusion E-val thresh:   <= %g\n",      esl_opt_GetReal(go, "--incdomE"));
  if (esl_opt_IsUsed(go, "--incdomT"))   fprintf(ofp, "# dom inclusion score thresh:   >= %g\n",      esl_opt_GetReal(go, "--incdomT"));
  if (esl_opt_IsUsed(go, "--inc_ga"))    fprintf(ofp, "# set inclusion thresholds to:     GA cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--inc_nc"))    fprintf(ofp, "# set inclusion thresholds to:     NC cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--inc_tc"))    fprintf(ofp, "# set inclusion thresholds to:     TC cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--max"))       fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n");
  if (esl_opt_IsUsed(go, "--F1"))        fprintf(ofp, "# MSV filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F1"));
  if (esl_opt_IsUsed(go, "--F2"))        fprintf(ofp, "# Vit filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F2"));
  if (esl_opt_IsUsed(go, "--F3"))        fprintf(ofp, "# Fwd filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F3"));
  if (esl_opt_IsUsed(go, "--nonull2"))   fprintf(ofp, "# null2 bias corrections:          off\n");
  if (esl_opt_IsUsed(go, "--nobias"))    fprintf(ofp, "# biased composition HMM filter:   off\n");
  if (esl_opt_IsUsed(go, "--nonull2"))   fprintf(ofp, "# null2 bias corrections:          off\n");
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed")==0)fprintf(ofp, "# random number seed:              one-time arbitrary\n");
    else                                    fprintf(ofp, "# random number seed set to:       %d\n", esl_opt_GetInteger(go, "--seed"));
  }
  if (esl_opt_IsUsed(go, "--textw"))     fprintf(ofp, "# max ASCII text line length:      %d\n",     esl_opt_GetInteger(go, "--textw"));
  if (esl_opt_IsUsed(go, "--notextw"))   fprintf(ofp, "# max ASCII text line length:      unlimited\n");
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}



int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go	   = NULL;              /* command line processing                  */
  FILE            *ofp     = NULL;	        /* output file for results (default stdout) */
  char            *seqfile = NULL;              /* file to read query sequence(s) from      */
  int              seqfmt  = eslSQFILE_UNKNOWN; /* format of seqfile                        */
  ESL_SQFILE      *sqfp    = NULL;              /* open seqfile                             */
  char            *hmmfile = NULL;              /* file to read HMM(s) from                 */
  ESL_ALPHABET    *abc     = NULL;              /* sequence alphabet                        */
  P7_BG           *bg      = NULL;              /* null model                               */
  ESL_STOPWATCH   *w       = NULL;              /* timing                                   */
  ESL_SQ          *qsq     = NULL;		/* query sequence                           */
  int              textw;
  int              status;
  int              sstatus, hstatus;

  /* Initializations */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */
  process_commandline(argc, argv, &go, &hmmfile, &seqfile);    
  w = esl_stopwatch_Create();
  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* Open the query sequence database */
  status = esl_sqfile_Open(seqfile, seqfmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",      seqfile);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",        seqfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, seqfile);
  qsq = esl_sq_Create();	/* initially in text mode, until first HMM is read. */

  /* Open the results output file */
  if ( esl_opt_IsOn(go, "-o")) {
    ofp = fopen(esl_opt_GetString(go, "-o"), "w");
    if (ofp == NULL) p7_Fail("Failed to open output file %s for writing\n", esl_opt_GetString(go, "-o"));
  } else 
    ofp = stdout;
 
  /* Note: header output is deferred until we know that HMM file has been opened
   * and its alphabet matches first query sequence. That happens in a deferred 
   * initialization block, because we don't know <abc> until first HMM is seen.
   */

  /* Outside loop: over each query sequence in <seqfile>. <abc> not known when we start! */
  while ( (sstatus = esl_sqio_Read(sqfp, qsq)) == eslOK)
    {
      P7_HMMFILE      *hfp     = NULL;		/* open HMM database file                   */
      P7_PIPELINE     *pli     = NULL;		/* processing pipeline                      */
      P7_TOPHITS      *th      = NULL;        	/* top-scoring sequence hits                */
      P7_OPROFILE     *om      = NULL;		/* target profile                           */
      int              qshown  = FALSE;		/* flag for deferred query name/acc/desc    */

      esl_stopwatch_Start(w);	                          

      /* Open the target profile database */
      status = p7_hmmfile_Open(hmmfile, p7_HMMDBENV, &hfp);
      if      (status == eslENOTFOUND) p7_Fail("Failed to open hmm file %s for reading.\n",                       hmmfile);
      else if (status == eslEFORMAT)   p7_Fail("Unrecognized format, trying to open hmm file %s for reading.\n",  hmmfile);
      else if (status != eslOK)        p7_Fail("Unexpected error %d in opening hmm file %s.\n",           status, hmmfile);  
      if (! hfp->is_pressed)           p7_Fail("Failed to open binary dbs for HMM file %s: use hmmpress first\n", hmmfile);
  
      /* Create processing pipeline and hit list */
      pli = p7_pipeline_Create(go, 100, 100, p7_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
      th  = p7_tophits_Create(); 
      pli->hfp = hfp;		/* for two-stage input, pipeline needs to know about <hfp> */

      /* unlike hmmsearch (for instance) we have to defer printing the query name,acc,desc until after output_header()
       * is called, so it can't go here (where it would seem to make sense); it's below instead (using <qshown>)
       */
      p7_pli_NewSeq(pli, qsq);

      if (abc != NULL)	/* once we're on sequence #>2, abc is known, bg exists */
	p7_bg_SetLength(bg, qsq->n);
	  
      /* Main loop: over each model in the target database; collect hits in <th> list */
      while ((hstatus = p7_oprofile_ReadMSV(hfp, &abc, &om)) == eslOK) /* <abc> gets set by first HMM in file  */
	{
	  /* One time only initializations after abc becomes known: */
	  if (bg == NULL) 	/* bg == NULL serves as a flag for the first-time init  */
	    {
	      output_header(ofp, go, hmmfile, seqfile);

	      bg = p7_bg_Create(abc);
	      if (esl_sq_Digitize(abc, qsq) != eslOK) p7_Die("alphabet mismatch");
	      esl_sqfile_SetDigital(sqfp, abc);
	      p7_bg_SetLength(bg, qsq->n);
	    }

	  /* Because we deferred output_header(), we had to defer query name/acc/desc output too,
	   * else they come out in the wrong order
	   */
	  if (! qshown) 
	    {			
	      fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n);
	      if (qsq->acc[0]  != 0) fprintf(ofp, "Accession:   %s\n", qsq->acc);
	      if (qsq->desc[0] != 0) fprintf(ofp, "Description: %s\n", qsq->desc);
	      qshown = TRUE;
	    }

	  p7_pli_NewModel(pli, om, bg);
	  p7_oprofile_ReconfigLength(om, qsq->n);

	  p7_Pipeline(pli, om, bg, qsq, th);

	  p7_oprofile_Destroy(om);
	  p7_pipeline_Reuse(pli);
	}
      if      (hstatus == eslEFORMAT)   p7_Fail("bad file format in HMM file %s",             hmmfile);
      else if (hstatus == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets",   hmmfile);
      else if (hstatus != eslEOF)       p7_Fail("Unexpected error in reading HMMs from %s",   hmmfile);      

      /* Print results */
      p7_tophits_Sort(th);
      p7_tophits_Threshold(th, pli);
      p7_tophits_Targets(ofp, th, pli, textw); fprintf(ofp, "\n\n");
      p7_tophits_Domains(ofp, th, pli, textw); fprintf(ofp, "\n\n");

      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, pli, w);
      fprintf(ofp, "//\n");

      p7_hmmfile_Close(hfp);
      p7_pipeline_Destroy(pli);
      p7_tophits_Destroy(th);
      esl_sq_Reuse(qsq);
    }
  if      (sstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s line %" PRId64 "):\n%s\n",
					    sqfp->filename, sqfp->linenumber, sqfp->errbuf);     
  else if (sstatus != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					    sstatus, sqfp->filename);

  esl_sq_Destroy(qsq);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);
  p7_bg_Destroy(bg);
  if ( esl_opt_IsOn(go, "-o")) fclose(ofp);
  esl_getopts_Destroy(go);
  return eslOK;
}



/*****************************************************************
 * @LICENSE@
 *****************************************************************/

