/* shmmer: search a protein sequence against a target database with alternative score systems
 * 
 * SC, Thu Mar  4 13:36:34 EST 2010
 * Serial version computing Forward, Viterbi, S/W and Miyazawa scores
 */

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sqio.h"
#include "esl_sq.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_scorematrix.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

/* Worker information to
 * compute Forward and
 * Viterbi scores
 */
typedef struct {
  enum            p7_palg_e alg_mode;    /* FWD or VIT */
  P7_BG          *bg;
  P7_PROFILE     *gm;
  P7_OPROFILE    *om;
  P7_TOPHITS     *th;
  int             ntargets;             /* number of target sequences    */
  int             by_E;		              /* use E-value as threshold      */
  double          E;                    /* per-target E-value threshold  */
  enum            p7_evsetby_e E_setby; /* --Erank or --Edist            */
  double          T;                    /* per-target score threshold    */
  FILE           *rsfp;                 /* scores from random sequences  */
} WORKER_PINFO;

/* Worker information
 * to compute Smith-Waterman
 * and Miyazawa scores
 */
typedef struct {
	enum            p7_salg_e alg_mode;   /* SW or MIY */
  P7_SCORESYS    *sm;
  P7_TOPHITS     *th;
  int             ntargets;             /* number of target sequences   */
  int             by_E;		              /* use E-value as threshold     */
  double          E;                    /* per-target E-value threshold */
  enum            p7_evsetby_e E_setby; /* --Erank only                 */
  double          T;                    /* per-target score threshold   */
  FILE           *rsfp;                 /* scores from random sequences */
} WORKER_SINFO;

#define ALGORITHMS  "--fwd,--vit,--sw,--miy"
#define PALGORITHMS "--fwd,--vit"
#define SALGORITHMS "--sw,--miy"
#define PGAPS       "--popen,--pextend"
#define SGAPS       "--sopen,--sextend"
#define REPOPTS     "-E,-T"

static ESL_OPTIONS options[] = {
  /* name           type         default   env  range   toggles   reqs   incomp                             help                                       docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL,              "show brief help on version and usage",                         1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL,              "direct output to file <f>, not stdout",                        2 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,      NULL,  NULL,  "--textw",         "unlimit ASCII text output line width",                         2 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",  NULL,  NULL,  "--notextw",       "set max width of ASCII text output lines",                     2 },
  /* Control of scoring system */
  { "--fwd",        eslARG_NONE,"default",NULL, NULL,ALGORITHMS,  NULL,  SGAPS,             "score seqs with the Forward algorithm",                        3 },
  { "--vit",        eslARG_NONE,   FALSE, NULL, NULL,ALGORITHMS,  NULL,  SGAPS,             "score seqs with the Viterbi algorithm",                        3 },
  { "--sw",         eslARG_NONE,   FALSE, NULL, NULL,ALGORITHMS,  NULL,  PGAPS,             "score seqs with the Smith-Waterman algorithm",                 3 },
  { "--miy",        eslARG_NONE,   FALSE, NULL, NULL,ALGORITHMS,  NULL,  PGAPS,             "score seqs with the Miyazawa algorithm",                       3 },
  { "--popen",      eslARG_REAL,  "0.02", NULL, "0<=x<0.5",NULL,  NULL,  SGAPS,             "gap open probability",                                         3 }, /* pHMMER default    */
  { "--pextend",    eslARG_REAL,   "0.4", NULL, "0<=x<1",  NULL,  NULL,  SGAPS,             "gap extend probability",                                       3 }, /* pHMMER default    */
  { "--sopen",      eslARG_REAL,  "11.0", NULL, "0<=x<100",NULL,  NULL,  PGAPS,             "gap open score",                                               3 }, /* NCBI-BLAST default for BLOSUM62 */
  { "--sextend",    eslARG_REAL,   "1.0", NULL, "0<=x<100",NULL,  NULL,  PGAPS,             "gap extend score",                                             3 }, /* NCBI-BLAST default for BLOSUM62 */
  { "--mxfile",     eslARG_INFILE,  NULL, NULL, NULL,      NULL,  NULL,  NULL,              "substitution score matrix [default: BLOSUM62]",                3 },

  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",     NULL,  NULL,  NULL,              "report sequences <= this E-value threshold in output",         4 }, /* and they are sorted from low to high E-values */
  { "-T",           eslARG_REAL,   FALSE, NULL,  NULL,     NULL,  NULL,  REPOPTS,           "report sequences >= this score threshold in output",           4 }, /* and they are sorted from high to low scores   */
  /* Control of E-value calculation */
  { "--Erank",      eslARG_NONE,"default",NULL, NULL, "--Edist","--Efile",NULL,             "use rank order statistics to compute E-values",               11 },
  { "--Efile",      eslARG_INFILE,  NULL, NULL, NULL,      NULL,  NULL,  "--Edist",         "Scores from random sequences for rank order statistics",      11 }, /* CHECK reqs and incomp for this option */
  { "--Edist",      eslARG_NONE,   FALSE, NULL, NULL, "--Erank",  NULL,  SALGORITHMS,       "use fitted distribution to compute E-values",                 11 },
  { "--EvL",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,              "length of sequences for Viterbi Gumbel mu fit",               11 }, /* For Viterbi and S/W scores */
  { "--EvN",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,              "number of sequences for Viterbi Gumbel mu fit",               11 }, /* For Viterbi and S/W scores */
  { "--EfL",        eslARG_INT,    "100", NULL,"n>0",      NULL,  NULL,  NULL,              "length of sequences for Forward exp tail tau fit",            11 }, /* For Forward and (hopefully) Miy scores */
  { "--EfN",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,              "number of sequences for Forward exp tail tau fit",            11 }, /* For Forward and (hopefully) Miy scores */
  { "--Eft",        eslARG_REAL,  "0.04", NULL,"0<x<1",    NULL,  NULL,  NULL,              "tail mass for Forward exponential tail tau fit",              11 }, /* For Forward and (hopefully) Miy scores */
  /* other options */
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",    NULL,  NULL,  NULL,              "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--qformat",    eslARG_STRING,  NULL, NULL, NULL,      NULL,  NULL,  NULL,              "assert query <seqfile> is in format <s>: no autodetection",   12 },
  { "--tformat",    eslARG_STRING,  NULL, NULL, NULL,      NULL,  NULL,  NULL,              "assert target <seqdb> is in format <s>>: no autodetection",   12 },
  /* Control of reporting thresholds */
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <query seqfile> <target seqdb>";
static char banner[] = "search a protein sequence against a protein database with alternative score systems";

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * We use it to provide configuration options to serial_master. serial_ploop
 * and serial_slooop get their own WORKER_INFO objects.
 *
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 *
 * THIS IS A SERIAL IMPLEMENTATION, SO DO NOT WORRY ABOUT ABOVE TOO MUCH FOR NOW.
 */
struct cfg_s {
  char            *qfile;             /* query sequence file                             */
  char            *dbfile;            /* database file                                   */
};

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_ploop (WORKER_PINFO *info, ESL_SQFILE *dbfp);     /* computes F or V scores          */
static int  serial_sloop (WORKER_SINFO *info, ESL_SQFILE *dbfp);     /* computes S/W or Miyazawa scores */

/* process_commandline()
 * Take argc, argv, and options; parse the command line;
 * display help/usage info.
 */
static void 
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_qfile, char **ret_dbfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     esl_fatal("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n",  go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n",  go->errbuf); goto ERROR; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 120=textwidth */

      puts("\noptions directing output:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      puts("\noptions controlling scoring system:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);

      puts("\noptions controlling reporting thresholds:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);

      puts("\noptions controlling E-value calculation:");
       esl_opt_DisplayHelp(stdout, go, 11, 2, 80);

      puts("\nother expert options:");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 2)    { puts("Incorrect number of command line arguments.");    goto ERROR; }
  if ((*ret_qfile  = esl_opt_GetArg(go, 1)) == NULL) { puts("Failed to get <qfile> argument on command line"); goto ERROR; }
  if ((*ret_dbfile = esl_opt_GetArg(go, 2)) == NULL) { puts("Failed to get <seqdb> argument on command line"); goto ERROR; }

  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere basic options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);  
}

/* HEADER  (prints options set to a non-default value) */
static int
output_header(FILE *ofp, ESL_GETOPTS *go, char *qfile, char *dbfile)
{
  p7_banner(ofp, go->argv[0], banner);  /* It includes some standard hmmer message along the banner defined here */
  
  fprintf(ofp, "# query sequence file:             %s\n", qfile);
  fprintf(ofp, "# target sequence database:        %s\n", dbfile);
  if (esl_opt_IsUsed(go, "-o"))          fprintf(ofp, "# output directed to file:         %s\n",      esl_opt_GetString(go, "-o"));
  if (esl_opt_IsUsed(go, "--notextw"))   fprintf(ofp, "# max ASCII text line length:      unlimited\n");
  if (esl_opt_IsUsed(go, "--textw"))     fprintf(ofp, "# max ASCII text line length:      %d\n",             esl_opt_GetInteger(go, "--textw"));  
  if (esl_opt_IsUsed(go, "--popen"))     fprintf(ofp, "# gap open probability:            %f\n",             esl_opt_GetReal  (go, "--popen"));
  if (esl_opt_IsUsed(go, "--pextend"))   fprintf(ofp, "# gap extend probability:          %f\n",             esl_opt_GetReal  (go, "--pextend"));
  if (esl_opt_IsUsed(go, "--sopen"))     fprintf(ofp, "# gap open score:                  %f\n",             esl_opt_GetReal  (go, "--sopen"));
  if (esl_opt_IsUsed(go, "--sextend"))   fprintf(ofp, "# gap extend score:                %f\n",             esl_opt_GetReal  (go, "--sextend"));
  if (esl_opt_IsUsed(go, "--mxfile"))    fprintf(ofp, "# subst score matrix:              %s\n",             esl_opt_GetString(go, "--mxfile"));
  if (esl_opt_IsUsed(go, "--fwd"))       fprintf(ofp, "# Computing Forward scores           \n");
  if (esl_opt_IsUsed(go, "--vit"))       fprintf(ofp, "# Computing Viterbi scores           \n");
  if (esl_opt_IsUsed(go, "--sw"))        fprintf(ofp, "# Computing Smith-Waterman scores    \n");
  if (esl_opt_IsUsed(go, "--miy"))       fprintf(ofp, "# Computing Miyazawa scores          \n");
  if (esl_opt_IsUsed(go, "-E"))          fprintf(ofp, "# sequence reporting threshold:    E-value <= %g\n",  esl_opt_GetReal(go, "-E"));
  if (esl_opt_IsUsed(go, "-T"))          fprintf(ofp, "# sequence reporting threshold:    score >= %g\n",    esl_opt_GetReal(go, "-T"));
  if (esl_opt_IsUsed(go, "--Erank"))     fprintf(ofp, "# using rank order statistics to compute E-values    \n");
  if (esl_opt_IsUsed(go, "--Edist"))     fprintf(ofp, "# using fitted distribution to compute E-values      \n");
  if (esl_opt_IsUsed(go, "--EvL") )      fprintf(ofp, "# seq length, Vit Gumbel mu fit:   %d\n",     esl_opt_GetInteger(go, "--EvL"));
  if (esl_opt_IsUsed(go, "--EvN") )      fprintf(ofp, "# seq number, Vit Gumbel mu fit:   %d\n",     esl_opt_GetInteger(go, "--EvN"));
  if (esl_opt_IsUsed(go, "--EfL") )      fprintf(ofp, "# seq length, Fwd exp tau fit:     %d\n",     esl_opt_GetInteger(go, "--EfL"));
  if (esl_opt_IsUsed(go, "--EfN") )      fprintf(ofp, "# seq number, Fwd exp tau fit:     %d\n",     esl_opt_GetInteger(go, "--EfN"));
  if (esl_opt_IsUsed(go, "--Eft") )      fprintf(ofp, "# tail mass for Fwd exp tau fit:   %f\n",     esl_opt_GetReal   (go, "--Eft"));
  if (esl_opt_IsUsed(go, "--seed"))  {
      if (esl_opt_GetInteger(go, "--seed") == 0) fprintf(ofp, "# random number seed:              one-time arbitrary\n");
      else                                       fprintf(ofp, "# random number seed set to:       %d\n", esl_opt_GetInteger(go, "--seed"));
    }
  if (esl_opt_IsUsed(go, "--qformat"))   fprintf(ofp, "# query <seqfile> format asserted: %s\n",     esl_opt_GetString(go, "--qformat"));
  if (esl_opt_IsUsed(go, "--tformat"))   fprintf(ofp, "# target <seqdb> format asserted:  %s\n",     esl_opt_GetString(go, "--tformat"));

  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}

/* MAIN
 * Serial only implementation
 */
int
main(int argc, char **argv)
{
  int             status   = eslOK;

  ESL_GETOPTS     *go  = NULL;	/* command line processing */
  struct cfg_s     cfg;         /* configuration data      */

  /* Set processor specific flags */
  impl_Init();

  /* Initialize what we can in the config structure (without knowing the alphabet yet) */
  cfg.qfile      = NULL;
  cfg.dbfile     = NULL;

  /* Process command-line */
  process_commandline(argc, argv, &go, &cfg.qfile, &cfg.dbfile);    

  /* Serial */
  status = serial_master(go, &cfg);

  /* Cleanup */
  esl_getopts_Destroy(go);

  return status;
}

/* serial_master()
 * A master can only return if it's successful. All errors are handled immediately and fatally with p7_Fail().
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;             /* output file for results (default stdout)  */
	FILE            *rsfp     = NULL;               /* input file of scores from random seqs.    */
  int              qformat  = eslSQFILE_UNKNOWN;  /* format of qfile                           */
  ESL_SQFILE      *qfp      = NULL;	              /* open qfile                                */
  ESL_SQ          *qsq      = NULL;               /* query sequence                            */
  int              dbformat = eslSQFILE_UNKNOWN;  /* format of dbfile                          */
  ESL_SQFILE      *dbfp     = NULL;               /* open dbfile                               */
  ESL_ALPHABET    *abc      = NULL;               /* sequence alphabet                         */
  P7_BUILDER      *bld      = NULL;               /* HMM construction configuration            */
  ESL_STOPWATCH   *w        = NULL;               /* for timing                                */
  int              seed;
  int              nquery   = 0;                  /* number of queries processex                */
  int              textw;                         /* set max width of ASCII text output lines   */
  int              status   = eslOK;              /* general status of different function calls */
  int              rstatus  = eslOK;              /* read status                                */
  int              sstatus  = eslOK;              /* search status                              */
  int              h;
  double           evalue;

  WORKER_PINFO     *pinfo   = NULL;
  WORKER_SINFO     *sinfo   = NULL;

  /* Initialization */
  abc     = esl_alphabet_Create(eslAMINO);      /* The resulting ESL_ALPHABET object includes input map for digitalization */
  w       = esl_stopwatch_Create();
  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* Start watch */
  esl_stopwatch_Start(w);

  /* Query format */
  if (esl_opt_IsOn(go, "--qformat")) {
    qformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat")); /* Here we autodetect the format if no --qformat option is given */
    if (qformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }

  /* Target format */
  if (esl_opt_IsOn(go, "--tformat")) {
    dbformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Probabilistic score system
   */
  if (esl_opt_GetBoolean(go, "--fwd") || esl_opt_GetBoolean(go, "--vit"))
  {
  	/* Allocate builder */
  	bld = p7_builder_Create(NULL, abc);                 /* P7_BUILDER is not initialized because go = NULL */
  	if ((seed = esl_opt_GetInteger(go, "--seed")) > 0)
  	{                           /* a little wasteful - we're blowing a couple of usec by reinitializing */
  		esl_randomness_Init(bld->r, seed);
  		bld->do_reseeding = TRUE;
  	}
  	/* Initialize builder for single sequence search */
  	bld->EvL = esl_opt_GetInteger(go, "--EvL");
  	bld->EvN = esl_opt_GetInteger(go, "--EvN");
  	bld->EfL = esl_opt_GetInteger(go, "--EfL");
  	bld->EfN = esl_opt_GetInteger(go, "--EfN");
  	bld->Eft = esl_opt_GetReal   (go, "--Eft");

  	/* Set score system */
  	status = p7_builder_SetScoreSystem(bld, esl_opt_GetString(go, "--mxfile"), NULL, esl_opt_GetReal(go, "--popen"), esl_opt_GetReal(go, "--pextend"));
  	if (status != eslOK) esl_fatal("Failed to set single query seq score system:\n%s\n", bld->errbuf);

  	/* Allocate pinfo */
  	ESL_ALLOC(pinfo, sizeof(*pinfo));

  	/* Initialize pinfo */
  	if (esl_opt_GetBoolean(go, "--fwd")) pinfo->alg_mode = FWD;
  	else                                 pinfo->alg_mode = VIT;

  	pinfo->bg         = p7_bg_Create(abc);
  	pinfo->ntargets   = 0;
  	if ( esl_opt_IsOn(go,"-T") )                           /* set threshold to scores   */
  	{
  		pinfo->T = esl_opt_GetReal(go, "-T");
  		pinfo->by_E = FALSE;
  	}
  	else                                                   /* set threshold to E-values */
  		{
  		pinfo->by_E = TRUE;
  		pinfo->E    = esl_opt_GetReal(go, "-E");
  		}

  	if (esl_opt_GetBoolean(go, "--Erank")) pinfo->E_setby =  RANK_ORDER; /* note, we still compute E-values when using a score threshold */
  	else                                   pinfo->E_setby =  FIT_DIST;
  }

  /* Non-probabilistic score system
   */
  else /* --sw or --miy */
  {
  	/* Allocate sinfo */
  	ESL_ALLOC(sinfo, sizeof(*sinfo));

  	/* Set score system */
  	ESL_ALLOC(sinfo->sm, sizeof(*sinfo->sm));
   	sinfo->sm->abc = abc;                     /* set the alphabet before setting the score matrix below */

  	if (esl_opt_GetBoolean(go, "--sw"))
  	{
  		sinfo->alg_mode = SW;
    	status = p7_builder_SetSWScoreSystem(sinfo->sm, esl_opt_GetString(go, "--mxfile"), NULL, esl_opt_GetReal(go, "--sopen"), esl_opt_GetReal(go, "--sextend"));
    	if (status != eslOK) esl_fatal("Failed to set single seq score system:\n%s\n", sinfo->sm->errbuf);
  	}

  	else
  	{
  		sinfo->alg_mode = MIY;
    	status = p7_builder_SetMiyScoreSystem(sinfo->sm, esl_opt_GetString(go, "--mxfile"), NULL, esl_opt_GetReal(go, "--sopen"), esl_opt_GetReal(go, "--sextend"));
    	if (status != eslOK) esl_fatal("Failed to set single seq score system:\n%s\n", sinfo->sm->errbuf);
  	}

  	sinfo->ntargets   = 0;
   	if ( esl_opt_IsOn(go,"-T") )                    /* set threshold to scores   */
  	{
   		sinfo->T = esl_opt_GetReal(go, "-T");
   		sinfo->by_E = FALSE;
  	}
    	else                                          /* set threshold to E-values */
    		{
    		sinfo->by_E = TRUE;
    		sinfo->E    = esl_opt_GetReal(go, "-E");
    		}

    if (esl_opt_GetBoolean(go, "--Erank")) sinfo->E_setby =  RANK_ORDER;
    else                                   sinfo->E_setby =  FIT_DIST;
  }

  /* Open output files */
  if (esl_opt_IsOn(go, "-o")) { if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) esl_fatal("Failed to open output file %s for writing\n", esl_opt_GetString(go, "-o")); }

  /* Open random scores file */
  if (esl_opt_IsOn(go, "--Efile")) {if ((rsfp = fopen(esl_opt_GetString(go, "--Efile"), "r")) == NULL) esl_fatal("Failed to open random scores file %s for reading\n", esl_opt_GetString(go, "--Efile")); }

  /* Open target file (autodetect format unless given) */
  status =  esl_sqfile_OpenDigital(abc, cfg->dbfile, dbformat, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open target sequence database %s for reading\n", cfg->dbfile);
  else if (status == eslEFORMAT)   esl_fatal("Target sequence database file %s is empty or misformatted\n", cfg->dbfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        esl_fatal("Unexpected error %d opening target sequence database file %s\n", status, cfg->dbfile);

  /* Open query file (autodetect format unless given) */
  status = esl_sqfile_OpenDigital(abc, cfg->qfile, qformat, NULL, &qfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open sequence file %s for reading\n", cfg->qfile);
  else if (status == eslEFORMAT)   esl_fatal("Sequence file %s is empty or misformatted\n", cfg->qfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        esl_fatal ("Unexpected error %d opening sequence file %s\n", status, cfg->qfile);

  /* Create digital query sequence */
  qsq  = esl_sq_CreateDigital(abc);

  /* Output header (non-default options) */
  output_header(ofp, go, cfg->qfile, cfg->dbfile);

  /* Outer loop over sequence queries
   */
  while ((rstatus = esl_sqio_Read(qfp, qsq)) == eslOK)
  {
      nquery++;
      if (qsq->n == 0) continue; /* skip zero length seqs as if they aren't even present */

      /* Start watch */
      esl_stopwatch_Start(w);

      /* Seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1)
      {
      	if (! esl_sqfile_IsRewindable(dbfp)) esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile); /* NOT SURE WHAT IS GOING ON HERE!!! */
      	esl_sqfile_Position(dbfp, 0);
      }

      /* Report query */
      fprintf(ofp, "Query:  %s  [L=%ld]\n", qsq->name, (long) qsq->n);
      if (qsq->acc[0]  != '\0') fprintf(ofp, "Accession:   %s\n", qsq->acc);
      if (qsq->desc[0] != '\0') fprintf(ofp, "Description: %s\n", qsq->desc);

      /* Inner loop over sequence targets
       */
      if (esl_opt_GetBoolean(go, "--fwd") || esl_opt_GetBoolean(go, "--vit"))
      {
      	pinfo->th = p7_tophits_Create();
      	p7_SingleBuilder(bld, qsq, pinfo->bg, NULL, NULL, &pinfo->gm, &pinfo->om); /* profile is P7_LOCAL by default, This is changed later to P7_UNILOCAL */
    	  pinfo->rsfp     = rsfp;
      	sstatus = serial_ploop(pinfo, dbfp);
      }
      else /* --sw or --miy */
      {
        sinfo->th       = p7_tophits_Create();
    	  sinfo->sm->dsq  = qsq->dsq;
    	  sinfo->sm->n    = qsq->n;
    	  sinfo->rsfp     = rsfp;
    	  sstatus         = serial_sloop(sinfo, dbfp);
      }

      /* Search status */
      switch(sstatus)
      {
      case eslEFORMAT:
    	  esl_fatal("Parse failed (sequence file %s):\n%s\n", dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
    	  break;
      case eslEOF:
    	  /* do nothing */
    	  break;
      default:
    	  esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
      }

      /* Sort and print hits */
      if (esl_opt_GetBoolean(go, "--fwd") || (esl_opt_GetBoolean(go, "--vit")))
      {
      	p7_tophits_Sort(pinfo->th); /* sorting by pvalue or score */

      	for (h = 0; h < pinfo->th->N; h++)
      	{
      		evalue = pinfo->th->hit[h]->pvalue * pinfo->ntargets;

      		if ( (pinfo->by_E && evalue <= pinfo->E) || !pinfo->by_E ) /* if by score, already filtered */
      		fprintf(ofp, "Score: %9.3f   P-value: %9.3g   E-value: %9.3g   %s\n", pinfo->th->hit[h]->score, pinfo->th->hit[h]->pvalue, evalue, pinfo->th->hit[h]->name);
      	}
      }

      else /* --sw or --miy */
      {
    	  p7_tophits_Sort(sinfo->th); /* sorting by pvalue or score */

        for (h = 0; h < sinfo->th->N; h++)
        {
        	evalue = sinfo->th->hit[h]->pvalue * sinfo->ntargets;

       		if ( (pinfo->by_E && evalue <= pinfo->E) || !pinfo->by_E ) /* if by score, already filtered */
        	fprintf(ofp, "Score: %9.3f   P-value: %9.3g   E-value: %9.3g   %s\n", sinfo->th->hit[h]->score, sinfo->th->hit[h]->pvalue, evalue, sinfo->th->hit[h]->name);
        }
      }

      /* Stop watch */
      esl_stopwatch_Stop(w);

      /* End query report */
      fprintf(ofp, "//\n");

      /* Reuse*/
      esl_sq_Reuse(qsq);

      /* Cleanup */
      if (esl_opt_GetBoolean(go, "--fwd") || esl_opt_GetBoolean(go, "--vit"))
        {
        p7_profile_Destroy(pinfo->gm);
        p7_oprofile_Destroy(pinfo->om);
        p7_tophits_Destroy(pinfo->th);
        }
      else
      {
       p7_tophits_Destroy(sinfo->th);
      }

    } /* end outer loop over query sequences */

  /* Query status */
  if      (rstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",qfp->filename, esl_sqfile_GetErrorBuf(qfp));
  else if (rstatus != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",rstatus, qfp->filename);

  /* Cleanup */
  if (esl_opt_GetBoolean(go, "--fwd") || esl_opt_GetBoolean(go, "--vit"))
    {
    p7_builder_Destroy(bld);
    p7_bg_Destroy(pinfo->bg);
    free(pinfo);
    }
  else /* --sw or --miy */
    {
  	/* careful here: sinfo->sm->sq points to qsq destroyed below */
    esl_scorematrix_Destroy(sinfo->sm->S);
    if (sinfo->sm->Q != NULL) free(sinfo->sm->Q);
    free(sinfo->sm);
    free(sinfo);
    }

  esl_sqfile_Close(dbfp);
  esl_sqfile_Close(qfp);
  esl_stopwatch_Destroy(w);
  esl_sq_Destroy(qsq);
  esl_alphabet_Destroy(abc);

  if (rsfp != NULL)   fclose(rsfp);
  if (ofp  != stdout) fclose(ofp);
  return eslOK;

 ERROR:
  return eslFAIL;
}

/************************************/
/* INNER LOOP OVER SEQUENCE TARGETS */
/************************************/

static int
serial_ploop(WORKER_PINFO *info, ESL_SQFILE *dbfp)
{
	ESL_SQ         *dbsq     = NULL;
	P7_GMX         *gx       = NULL;
	P7_OMX         *ox       = NULL;
	P7_HIT         *hit      = NULL;     /* ptr to the current hit output data */
	int             status;              /* generic (reusable) function status */
	int             rstatus;             /* read status                        */
	int             cstatus;             /* computation status                 */
	float           nullsc;              /* null model score in nats (only transitions in the null model)     */
	float           rsc;                 /* raw lod sequence score in nats (only emissions in the null model) */
	float           sc;                  /* final lod sequence score in bits                                  */
	double           P;                  /* P-value of a score (in bits)         */

	dbsq = esl_sq_CreateDigital(info->gm->abc);

	/* Create dynamic programming matrices */
	gx = p7_gmx_Create(info->gm->M, 100);    /* target length of 100 is a dummy for now */
	ox = p7_omx_Create(info->om->M, 0, 100); /* target length of 100 is a dummy for now */  /* allocate memory efficient linear arrays */

	/* INNER LOOP OVER SEQUENCE TARGETS */
	while ((rstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
	{
		/* Count target sequences */
		info->ntargets++;                     /* phmmer counts all target sequences including empty ones, so I do it too */

		if (dbsq->n == 0) return eslOK;       /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */

		/* Reconfig all models using target length  */
		p7_bg_SetLength(info->bg, dbsq->n);
		p7_ReconfigUnihit(info->gm, dbsq->n);            /* Unihit model (P7_UNILOCAL) */ /* MAKE SURE THIS IS CORRECT */
		p7_oprofile_ReconfigUnihit(info->om, dbsq->n);   /* Unihit model (P7_UNILOCAL) */

		/* Reconfig dp matrices using target length */
		p7_gmx_GrowTo(gx, info->gm->M, dbsq->n);
		p7_omx_GrowTo(ox, info->om->M, 0, dbsq->n);         /* NOT SURE I NEED TO DO THIS FOR ox */

		if (info->alg_mode == FWD)  /* I may want to move the if outside the target sequence loop!!! */
		{
			cstatus = p7_ForwardParser(dbsq->dsq, dbsq->n, info->om, ox, &rsc); /* rsc is in nats */
			if (cstatus == eslERANGE) cstatus = p7_GForward(dbsq->dsq, dbsq->n, info->gm, gx, &rsc); /* if a filter overflows, failover to slow versions */
			else if (cstatus != eslOK) esl_fatal ("Failed to compute Forward scores!\n");

			/* Base null model score */
			p7_bg_NullOne(info->bg, dbsq->dsq, dbsq->n, &nullsc); /* nullsc is in nats? only scores transitions!!! */

			/* Biased composition HMM filtering
			 * DISABLED
			 */

			/* Calculate the null2-corrected per-seq score
			 * DISABLED
			 */

			/* Score (in bits) */
			sc = (rsc - nullsc) / eslCONST_LOG2;

			/* Calculate P-value */
			if (info->E_setby == RANK_ORDER) P = rank_order(info->rsfp, sc);
			else P = esl_exp_surv(sc,info->gm->evparam[p7_FTAU], info->gm->evparam[p7_FLAMBDA]);
		} /* end forward computation */

		else /* VIT */
		{
			cstatus = p7_ViterbiFilter(dbsq->dsq, dbsq->n, info->om, ox, &rsc); /* rsc is in nats */
			if (cstatus == eslERANGE) cstatus = p7_GViterbi(dbsq->dsq, dbsq->n, info->gm, gx, &rsc); /* if a filter overflows, failover to slow versions */
			else if (cstatus != eslOK) esl_fatal ("Failed to compute Viterbi scores!\n");

			/* Base null model score */
			p7_bg_NullOne(info->bg, dbsq->dsq, dbsq->n, &nullsc);

			/* Biased composition HMM filtering
			 * DISABLED
			 */

			/* Calculate the null2-corrected per-seq score
			 * DISABLED
			 */

			/* Score (in bits) */
			sc = (rsc - nullsc) / eslCONST_LOG2;

			/* Calculate P-value */
			if (info->E_setby == RANK_ORDER) P = rank_order(info->rsfp, sc);
			else P = esl_gumbel_surv(sc, info->gm->evparam[p7_VMU], info->gm->evparam[p7_VLAMBDA]);  /* ARE MU AND LAMBDA THE SAME FOR gm AND om??? they should be! */ /* GOT THE CORRECTED LAMBDA??? */
		} /* end viterbi computation */

		/* Hitlist */
		if (info->by_E == TRUE) /* threshold and sort by E-values */
		{
			if (P * info->ntargets <= info->E) /* Calculated E-value is a lower bound because we haven't yet read the full target database */
			{
				p7_tophits_CreateNextHit(info->th, &hit); /* allocates new hit structure */
				if ((status  = esl_strdup(dbsq->name, -1, &(hit->name)))  != eslOK) esl_fatal("allocation failure"); /* CHECK I AM USING THIS FUNCTION PROPERLY!!! */
				hit->score   = sc;
				hit->pvalue  = P;
				hit->sortkey = -log(P); /* per-seq output sorts on bit score if inclusion is by score  */
			}
		}

		else /* threshold and sort by scores */
		{
			if (sc >= info->T)
			{
				p7_tophits_CreateNextHit(info->th, &hit); /* allocates new hit structure */
				if ((status  = esl_strdup(dbsq->name, -1, &(hit->name)))  != eslOK) esl_fatal("allocation failure"); /* CHECK I AM USING THIS FUNCTION PROPERLY!!! */
				hit->score   = sc;
				hit->pvalue  = P;
				hit->sortkey = sc; /* per-seq output sorts on bit score if inclusion is by score  */
			}
		}

		/* Reuse */
		esl_sq_Reuse(dbsq);
		p7_gmx_Reuse(gx);
		p7_omx_Reuse(ox);

	} /* end loop over seq. targets */

	/* Cleanup */
	esl_sq_Destroy(dbsq);
	p7_gmx_Destroy(gx);
	p7_omx_Destroy(ox);

	return rstatus; /* read status */
}

static int
serial_sloop(WORKER_SINFO *info, ESL_SQFILE *dbfp)
{
	ESL_SQ   *dbsq     = NULL;
	P7_GMX   *gx       = NULL;
	P7_HIT   *hit      = NULL;     /* ptr to the current hit output data   */
	int      status;               /* generic (reusable) function status   */
	int      rstatus;              /* read status                          */
	int      cstatus;              /* computation status                   */

	float    rsc;                  /* raw sequence score                   */
	float    sc;                   /* final sequence score in bits         */
	double    P;                   /* P-value of a score (in bits)         */

	dbsq = esl_sq_CreateDigital(info->sm->abc);

	/* Create dynamic programming matrices */
	if ((gx = p7_gmx_Create(info->sm->n, 100)) == NULL) esl_fatal("Dynamic programming matrix allocation failure");; /* Length target (dummy 100 here) is the number of rows in the matrix */

	/* INNER LOOP OVER SEQUENCE TARGETS */
	while ((rstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
	{
		/* Count target sequences */
		info->ntargets++;                     /* phmmer counts all target sequences including empty ones, so I do it too */

		if (dbsq->n == 0) return eslOK;       /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */

		/* Reconfig dp matrices using target length */
		if ((status = p7_gmx_GrowTo(gx, info->sm->n, dbsq->n)) != eslOK) esl_fatal("Dynamic programming matrix reallocation failure");

		if (info->alg_mode == SW)
		{

			/* Fast sse S/W implementation here (I could use Farrar's implementation)
			 */

			/* Slow S/W implementation (it will do for now!) */
			cstatus = p7_GSmithWaterman(dbsq->dsq, dbsq->n, info->sm, gx, &rsc);
			if (cstatus != eslOK)  esl_fatal ("Failed to compute slow Smith-Waterman scores!\n");

			/* Score (units unchanged) */
			sc = rsc;

			/* Calculate P-value */
			if (info->E_setby == RANK_ORDER) P = rank_order(info->rsfp, sc);
			else esl_fatal ("P-values for Smith-Waterman scores can only be computed with rank order statistics for now!\n");
		}

		else /* MIY */
		{
			/* Fast Miyazawa implementation (not implemented)
			 */

			/* Slow Miyazawa implementation (it will do for now!) */
			cstatus = p7_GMiyazawa(dbsq->dsq, dbsq->n, info->sm, gx, &rsc);
			if (cstatus != eslOK)  esl_fatal ("Failed to compute slow Miyazawa scores!\n");

			/* Score (units unchanged) */
			sc = rsc;

			/* Calculate P-value */
			if (info->E_setby == RANK_ORDER) P = rank_order(info->rsfp, sc);
			else esl_fatal ("P-values for Miyazawa scores can only be computed with rank order statistics for now!\n");
		}

		/* Hitlist */
			if (info->by_E == TRUE) /* threshold and sort by E-values */
			{
				if (P * info->ntargets <= info->E) /* Calculated E-value is a lower bound because we haven't yet read the full target database */
				{
					p7_tophits_CreateNextHit(info->th, &hit); /* allocates new hit structure */
					if ((status  = esl_strdup(dbsq->name, -1, &(hit->name)))  != eslOK) esl_fatal("allocation failure"); /* CHECK I AM USING THIS FUNCTION PROPERLY!!! */
					hit->score   = sc;
					hit->pvalue  = P;
					hit->sortkey = -log(P); /* per-seq output sorts on bit score if inclusion is by score  */
				}
			}

			else /* threshold and sort by scores */
			{
				if (sc >= info->T)
				{
					p7_tophits_CreateNextHit(info->th, &hit); /* allocates new hit structure */
					if ((status  = esl_strdup(dbsq->name, -1, &(hit->name)))  != eslOK) esl_fatal("allocation failure"); /* CHECK I AM USING THIS FUNCTION PROPERLY!!! */
					hit->score   = sc;
					hit->pvalue  = P;
					hit->sortkey = sc; /* per-seq output sorts on bit score if inclusion is by score  */
				}
			}

		/* Reuse */
		esl_sq_Reuse(dbsq);
		p7_gmx_Reuse(gx);

	} /* end loop over seq. targets */

	/* Cleanup */
	esl_sq_Destroy(dbsq);
	p7_gmx_Destroy(gx);

	return rstatus; /* read status */
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
