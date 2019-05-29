/* Add mask line to a multiple sequence alignment, based on
 * either the position range in the alignment or the position
 * range in the model that would be produced using the alignment
 * and the same set of flags in hmmbuild
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#undef HAVE_PTHREAD

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msaweight.h"
#include "esl_msacluster.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"
#include "esl_regexp.h"

#include "hmmer.h"

typedef struct {
  P7_BG	           *bg;
  P7_BUILDER       *bld;
} WORKER_INFO;


#define ALPHOPTS "--amino,--dna,--rna"                         /* Exclusive options for alphabet choice */
#define CONOPTS "--fast,--hand"                                /* Exclusive options for model construction                    */
#define WGTOPTS "--wgsc,--wblosum,--wpb,--wnone,--wgiven"      /* Exclusive options for relative weighting                    */
#define RANGEOPTS "--modelrange,--alirange,--ali2model,--model2ali"      /* Exclusive options for relative weighting                    */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { (char *) "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, (char *) "show brief help on version and usage",                  1 },
  { (char *) "-o",        eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,    NULL, (char *) "direct summary output to file <f>, not stdout",         1 },
/* Selecting the alphabet rather than autoguessing it */
  { (char *) "--amino",   eslARG_NONE,   FALSE, NULL, NULL,  (char *)  ALPHOPTS,    NULL,     NULL, (char *) "input alignment is protein sequence data",              2 },
  { (char *) "--dna",     eslARG_NONE,   FALSE, NULL, NULL,  (char *)  ALPHOPTS,    NULL,     NULL, (char *) "input alignment is DNA sequence data",                  2 },
  {(char *)  "--rna",     eslARG_NONE,   FALSE, NULL, NULL,  (char *)  ALPHOPTS,    NULL,     NULL, (char *) "input alignment is RNA sequence data",                  2 },
/* Alternate model construction strategies */
  { (char *) "--fast",    eslARG_NONE,(char *) "default",NULL, NULL,   (char *)  CONOPTS,    NULL,     NULL, (char *) "assign cols w/ >= symfrac residues as consensus",       3 },
  { (char *) "--hand",    eslARG_NONE,   FALSE, NULL, NULL,    (char *) CONOPTS,    NULL,     NULL, (char *) "manual construction (requires reference annotation)",   3 },
  { (char *) "--symfrac", eslARG_REAL,  (char *)  "0.5", NULL,(char *)  "0<=x<=1", NULL,  (char *)  "--fast",   NULL, (char *) "sets sym fraction controlling --fast construction",     3 },
  {(char *)  "--fragthresh",eslARG_REAL,(char *) "0.5", NULL, (char *) "0<=x<=1", NULL,     NULL,     NULL, (char *) "if L <= x*alen, tag sequence as a fragment",            3 },
/* Alternate relative sequence weighting strategies */
  /* --wme not implemented in HMMER3 yet */
  {(char *)  "--wpb",     eslARG_NONE,(char *) "default",NULL, NULL,  (char *)   WGTOPTS,    NULL,      NULL, (char *) "Henikoff position-based weights",                      4 },
  { (char *) "--wgsc",    eslARG_NONE,   NULL,  NULL, NULL,  (char *)   WGTOPTS,    NULL,      NULL, (char *) "Gerstein/Sonnhammer/Chothia tree weights",             4 },
  { (char *) "--wblosum", eslARG_NONE,   NULL,  NULL, NULL,  (char *)   WGTOPTS,    NULL,      NULL, (char *) "Henikoff simple filter weights",                       4 },
  { (char *) "--wnone",   eslARG_NONE,   NULL,  NULL, NULL,   (char *)  WGTOPTS,    NULL,      NULL, (char *) "don't do any relative weighting; set all to 1",        4 },
  { (char *) "--wgiven",  eslARG_NONE,   NULL,  NULL, NULL,   (char *)  WGTOPTS,    NULL,      NULL, (char *) "use weights as given in MSA file",                     4 },
  { (char *) "--wid",     eslARG_REAL, (char *) "0.62",  NULL,(char *) "0<=x<=1",   NULL,(char *) "--wblosum",   NULL, (char *) "for --wblosum: set identity cutoff",                   4 },
/* mask ranges */
  { (char *) "--modelrange", eslARG_STRING, NULL, NULL, NULL,  NULL,    NULL,    (char *)  RANGEOPTS,  (char *) "range(s) for mask(s) in model coordinates", 5 },
  {(char *)  "--alirange",   eslARG_STRING, NULL, NULL, NULL,  NULL,    NULL,    (char *)  RANGEOPTS,  (char *) "range(s) for mask(s) in alignment coordinates", 5 },
  {(char *)  "--apendmask",    eslARG_NONE, NULL,  NULL, NULL, NULL, (char *)  WGTOPTS,         NULL,  (char *) "add to existing mask (default ignores to existing mask)",    5 },
  {(char *) "--model2ali",  eslARG_STRING, NULL, NULL, NULL,  NULL,    NULL,    (char *)  RANGEOPTS,  (char *) "print model ranges corresponding to input alignment ranges", 5 },
  { (char *) "--ali2model",  eslARG_STRING, NULL, NULL, NULL,  NULL,    NULL,    (char *)  RANGEOPTS,  (char *) "print alignment ranges corresponding to input model ranges", 5 },

/* Other options */
  { (char *) "--informat", eslARG_STRING, NULL, NULL, NULL,      NULL,      NULL,    NULL, (char *) "assert input alifile is in format <s> (no autodetect)", 8 },
  { (char *) "--seed",     eslARG_INT,   (char *) "42", NULL, (char *) "n>=0",     NULL,      NULL,    NULL, (char *) "set RNG seed to <n> (if 0: one-time arbitrary seed)",   8 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


static char usage[]  = "[-options] <msafile> <postmsafile>";
static char banner[] = "appending modelmask line to a multiple sequence alignments";

static int output_header(const ESL_GETOPTS *go, FILE *ofp, char *alifile, char *postmsafile);


static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_alifile, char **ret_postalifile)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK) { if (printf("Failed to process environment:\n%s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) { if (printf("Failed to parse command line:\n%s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK) { if (printf("Failed to parse command line:\n%s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* help format: */
  if (esl_opt_GetBoolean(go, (char *) "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);

      if (puts("\nBasic options:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);

      if (puts("\nMask range options (format:  --xxx 10-20,30-40 ) :") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);

      if (puts("\nOptions for selecting alphabet rather than guessing it:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);

      if (puts("\nAlternative model construction strategies:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);

      if (puts("\nAlternative relative sequence weighting strategies:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);

      if (puts("\nOther options:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                  > 2)    { if (puts("Incorrect number of command line arguments.")          < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_alifile     = esl_opt_GetArg(go, 1)) == NULL) { if (puts("Failed to get <msafile> argument on command line")     < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_IsUsed(go, (char *) "--alirange") || esl_opt_IsUsed(go,(char *)  "--modelrange") ) {
    if ((*ret_postalifile = esl_opt_GetArg(go, 2)) == NULL) { if (puts("Failed to get <postmsafile> argument on command line")     < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  }

  if (strcmp(*ret_alifile, "-") == 0 && ! esl_opt_IsOn(go,(char *)  "--informat"))
    { if (puts("Must specify --informat to read <alifile> from stdin ('-')") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }


  *ret_go = go;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere basic options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  printf("\nTo see more help on other available options, do:\n  %s -h\n\n", argv[0]);
  esl_getopts_Destroy(go);
  exit(1);  

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

static int
output_header(const ESL_GETOPTS *go, FILE *ofp, char *alifile, char *postmsafile)
{
  p7_banner(ofp, go->argv[0], banner);

  if (fprintf(ofp, "# input alignment file:             %s\n", alifile) < 0)     ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go,(char *)  "--alirange") || esl_opt_IsUsed(go, (char *) "--modelrange") ) {
    if (fprintf(ofp, "# output alignment file:            %s\n", postmsafile) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }

  if (esl_opt_IsUsed(go, (char *) "--alirange")   && fprintf(ofp, "# alignment range:                  %s\n",        esl_opt_GetString(go,(char *)  "--alirange"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--modelrange") && fprintf(ofp, "# model range:                      %s\n",        esl_opt_GetString(go,(char *)  "--modelrange"))< 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--apendmask")  && fprintf(ofp, "# add to existing mask:             [on]\n"                                            )< 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, (char *) "--model2ali")   && fprintf(ofp, "# ali ranges for model range:      %s\n",        esl_opt_GetString(go, (char *) "--model2ali"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--ali2model")   && fprintf(ofp, "# model ranges for ali range:      %s\n",        esl_opt_GetString(go, (char *) "--ali2model"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");


  if (esl_opt_IsUsed(go, (char *) "-o")           && fprintf(ofp, "# output directed to file:          %s\n",        esl_opt_GetString(go, (char *) "-o"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--amino")      && fprintf(ofp, "# input alignment is asserted as:   protein\n")                                        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--dna")        && fprintf(ofp, "# input alignment is asserted as:   DNA\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--rna")        && fprintf(ofp, "# input alignment is asserted as:   RNA\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--fast")       && fprintf(ofp, "# model architecture construction:  fast/heuristic\n")                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--hand")       && fprintf(ofp, "# model architecture construction:  hand-specified by RF annotation\n")                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--symfrac")    && fprintf(ofp, "# sym fraction for model structure: %.3f\n",      esl_opt_GetReal(go, (char *) "--symfrac"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--fragthresh") && fprintf(ofp, "# seq called frag if L <= x*alen:   %.3f\n",      esl_opt_GetReal(go, (char *) "--fragthresh")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--wpb")        && fprintf(ofp, "# relative weighting scheme:        Henikoff PB\n")                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--wgsc")       && fprintf(ofp, "# relative weighting scheme:        G/S/C\n")                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--wblosum")    && fprintf(ofp, "# relative weighting scheme:        BLOSUM filter\n")                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--wnone")      && fprintf(ofp, "# relative weighting scheme:        none\n")                                           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--wid")        && fprintf(ofp, "# frac id cutoff for BLOSUM wgts:   %f\n",        esl_opt_GetReal(go, (char *) "--wid"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, (char *) "--seed"))  {
    if (esl_opt_GetInteger(go, (char *) "--seed") == 0  && fprintf(ofp,"# random number seed:               one-time arbitrary\n")                        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if                              (  fprintf(ofp,"# random number seed set to:        %d\n",         esl_opt_GetInteger(go, (char *) "--seed"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }

  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}

/* lifted from esl-sfetch */
static int
parse_coord_string(const char *cstring, uint32_t *ret_start, uint32_t *ret_end)
{
  ESL_REGEXP *re = esl_regexp_Create();
  char        tok1[32];
  char        tok2[32];

  if (esl_regexp_Match(re, "^(\\d+)\\D+(\\d*)$", cstring) != eslOK) esl_fatal("-c takes arg of subseq coords <from>..<to>; %s not recognized", cstring);
  if (esl_regexp_SubmatchCopy(re, 1, tok1, 32)            != eslOK) esl_fatal("Failed to find <from> coord in %s", cstring);
  if (esl_regexp_SubmatchCopy(re, 2, tok2, 32)            != eslOK) esl_fatal("Failed to find <to> coord in %s",   cstring);

  *ret_start = atol(tok1);
  *ret_end   = (tok2[0] == '\0') ? 0 : atol(tok2);

  esl_regexp_Destroy(re);
  return eslOK;
}


int
main(int argc, char **argv)
{

  int i,j;

  ESL_GETOPTS     *go = NULL;	/* command line processing                 */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();

  int              status;
  ESL_MSA      *msa         = NULL;
  FILE         *ofp         = NULL;    /* output file (default is stdout) */
  ESL_ALPHABET *abc         = NULL;    /* digital alphabet */

  char         *alifile;               /* name of the alignment file we're building HMMs from  */
  ESL_MSAFILE  *afp         = NULL;    /* open alifile  */
  int           fmt;                   /* format code for alifile */

  char         *postmsafile;           /* optional file to resave annotated, modified MSAs to  */
  FILE         *postmsafp = NULL;      /* open <postmsafile>, or NULL */

  int           mask_range_cnt = 0;
  uint32_t      mask_starts[100];      // over-the-top allocation.
  uint32_t      mask_ends[100];

  char         *rangestr;
  char         *range;


  int     *map = NULL; /* map[i]=j,  means model position i comes from column j of the alignment; 1..alen */

  int    keep_mm;

  p7_Init();

  alifile     = NULL;
  postmsafile = NULL;

  /* Parse the command line
   */
  process_commandline(argc, argv, &go, &alifile, &postmsafile);
  keep_mm = esl_opt_IsUsed(go, (char *)  "--apendmask");

  /* Initialize what we can in the config structure (without knowing the alphabet yet).
   * Fields controlled by masters are set up in usual_master() or mpi_master()
   * Fields used by workers are set up in mpi_worker()
   */
  ofp         = NULL;
  fmt         = eslMSAFILE_UNKNOWN;    /* autodetect alignment format by default. */
  afp         = NULL;
  abc         = NULL;

  if (esl_opt_IsOn(go, (char *) "--informat")) {
    fmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, (char *) "--informat"));
    if (fmt == eslMSAFILE_UNKNOWN) p7_Fail((char *)  "%s is not a recognized input sequence file format\n", esl_opt_GetString(go, (char *) "--informat"));
  }


  /* Parse the ranges */

  if (esl_opt_IsUsed(go,(char *)  "--alirange")) {
    esl_strdup(esl_opt_GetString(go,(char *)  "--alirange"), -1, &rangestr) ;
  } else if (esl_opt_IsUsed(go, (char *) "--modelrange")) {
    esl_strdup(esl_opt_GetString(go,(char *)  "--modelrange"), -1, &rangestr) ;
  } else if (esl_opt_IsUsed(go, (char *) "--model2ali")) {
    esl_strdup(esl_opt_GetString(go,(char *)  "--model2ali"), -1, &rangestr) ;
  } else if (esl_opt_IsUsed(go, (char *) "--ali2model")) {
    esl_strdup(esl_opt_GetString(go, (char *) "--ali2model"), -1, &rangestr) ;
  } else {
    if (puts("Must specify mask range with --modelrange, --alirange, --model2ali, or --ali2model\n") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto ERROR;
  }

  while ( (status = esl_strtok(&rangestr, (char *)  ",",  &range) ) == eslOK) {
    parse_coord_string(range, mask_starts + mask_range_cnt, mask_ends + mask_range_cnt );
    mask_range_cnt++;
  }



  /* Start timing. */
  esl_stopwatch_Start(w);



  /* Open files, set alphabet.
   *   afp       - open alignment file for input
   *   abc       - alphabet expected or guessed in ali file
   *   postmsafp - open MSA output file
   *   ofp       - optional open output file, or stdout
   */
  if      (esl_opt_GetBoolean(go, (char *) "--amino"))   abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, (char *) "--dna"))     abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, (char *) "--rna"))     abc = esl_alphabet_Create(eslRNA);
  else                                          abc = NULL;
  
  status = esl_msafile_Open(&abc, alifile, NULL, fmt, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);

  if (esl_opt_IsUsed(go,(char *)  "--alirange") || esl_opt_IsUsed(go, (char *) "--modelrange") ) {
    postmsafp = fopen(postmsafile, "w");
    if (postmsafp == NULL) p7_Fail((char *)  "Failed to MSA output file %s for writing", postmsafile);
  }

  if (esl_opt_IsUsed(go, (char *) "-o")) 
    {
      ofp = fopen(esl_opt_GetString(go,(char *)  "-o"), "w");
      if (ofp == NULL) p7_Fail((char *)  "Failed to open -o output file %s\n", esl_opt_GetString(go, (char *) "-o"));
    } 
  else ofp = stdout;


  /* Looks like the i/o is set up successfully...
   * Initial output to the user
   */
  output_header(go, ofp, alifile, postmsafile);                                  /* cheery output header                                */

  /* read the alignment */
  if ((status = esl_msafile_Read(afp, &msa)) != eslOK)  esl_msafile_ReadFailure(afp, status);


  if (esl_opt_IsUsed(go,(char *)  "--alirange") || esl_opt_IsUsed(go,(char *)  "--modelrange") ) {
    /* add/modify mmline for the mask */
    if (msa->mm == NULL) {
      ESL_ALLOC(msa->mm, msa->alen);
      keep_mm = FALSE;
    }

    if (!keep_mm)
      for (i=0; i<msa->alen; i++) msa->mm[i] = '.';

  }

  // convert model coordinates to alignment coordinates, if necessary
  if (esl_opt_IsUsed(go, (char *) "--modelrange") || esl_opt_IsUsed(go, (char *) "--model2ali") || esl_opt_IsUsed(go, (char *) "--ali2model") ) {
    int      apos, idx;
    float    r;            /* weighted residue count              */
    float    totwgt;       /* weighted residue+gap count          */
    float    symfrac;

    //same as p7_builder relative_weights
    if      (esl_opt_IsOn(go,(char *)  "--wnone")  )                  { esl_vec_DSet(msa->wgt, msa->nseq, 1.); }
    else if (esl_opt_IsOn(go, (char *) "--wgiven") )                  ;
    else if (esl_opt_IsOn(go,(char *)  "--wpb")    )                  esl_msaweight_PB(msa);
    else if (esl_opt_IsOn(go, (char *) "--wgsc")   )                  esl_msaweight_GSC(msa);
    else if (esl_opt_IsOn(go,(char *)  "--wblosum"))                  esl_msaweight_BLOSUM(msa, esl_opt_GetReal(go, (char *) "--wid"));

    symfrac = esl_opt_GetReal(go, (char *) "--symfrac");

    if ((status =  esl_msa_MarkFragments_old(msa, esl_opt_GetReal(go, "--fragthresh")))           != eslOK) goto ERROR;

    // Determine weighted sym freq in each column, build a map of model mask coordinates to alignment coords
    ESL_ALLOC(map, sizeof(int)     * (msa->alen+1));
    i = 0;

    if ( esl_opt_IsOn(go, (char *) "--hand")) {
      if (msa->rf == NULL)      p7_Fail((char *) "Model file doe not contain an RF line, required for --hand.\n");
        /* Watch for off-by-one. rf is [0..alen-1]*/
       for (apos = 1; apos <= msa->alen; apos++) {
         if (!esl_abc_CIsGap(msa->abc, msa->rf[apos-1]) ) {
           map[i] = apos;
           i++;
         }
       }

    } else {

      for (apos = 1; apos <= msa->alen; apos++)
      {

          r = totwgt = 0.;
          for (idx = 0; idx < msa->nseq; idx++)
          {
            if       (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) { r += msa->wgt[idx]; totwgt += msa->wgt[idx]; }
            else if  (esl_abc_XIsGap(msa->abc,     msa->ax[idx][apos])) {                     totwgt += msa->wgt[idx]; }
            else if  (esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos])) continue;
          }

          if (r > 0. && r / totwgt >= symfrac) {
            map[i] = apos;
            i++;
          }
      }
    }


    if ( esl_opt_IsUsed(go, (char *) "--model2ali") ) {
      //print mapping
      printf ("model coordinates     alignment coordinates\n");
      for (i=0; i<mask_range_cnt; i++)
        printf ("%8d..%-8d -> %8d..%-8d\n", mask_starts[i], mask_ends[i], map[mask_starts[i]-1], map[mask_ends[i]-1]);
    } else if ( esl_opt_IsUsed(go,(char *)  "--ali2model") ) {
      //print mapping  (requires scanning the inverted map
      int alistart = 0;
      int aliend = 0;
      printf ("alignment coordinates     model coordinates\n");
      for (i=0; i<mask_range_cnt; i++) {
        //find j for ali positions
        while (map[alistart] < mask_starts[i] )
          alistart++;
        aliend = alistart;
        while (map[aliend] < mask_ends[i] )
          aliend++;

        printf ("   %8d..%-8d -> %8d..%-8d\n", map[alistart], map[aliend], alistart+1, aliend+1);
      }
    } else {
      //convert the mask coords based on map
      for (i=0; i<mask_range_cnt; i++) {
          mask_starts[i] = map[mask_starts[i]-1]; //-1 because mmline is offset by one relative to the 1-base alignment
          mask_ends[i]   = map[mask_ends[i]-1];
      }
    }
  }

  if (esl_opt_IsUsed(go,(char *)   "--alirange") || esl_opt_IsUsed(go, (char *)  "--modelrange") ) {
    //overwrite '.' with 'm' everywhere the range says to do it
    for (i=0; i<mask_range_cnt; i++)
      for (j=mask_starts[i]; j<=mask_ends[i]; j++)
        msa->mm[j-1] = 'm';

    if ((status = esl_msafile_Write(postmsafp, msa, eslMSAFILE_STOCKHOLM))  != eslOK) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  }

  esl_stopwatch_Stop(w);
  free(map);
  if (esl_opt_IsOn(go, (char *)  "-o"))  fclose(ofp);
  if (postmsafp) fclose(postmsafp);
  if (afp)       esl_msafile_Close(afp);
  if (abc)       esl_alphabet_Destroy(abc);

  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  return 0;


 ERROR:
  if (map) free(map);
  if (w)   esl_stopwatch_Destroy(w);
  if (go)  esl_getopts_Destroy(go);
   return eslFAIL;
}

