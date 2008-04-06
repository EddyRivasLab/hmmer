/* 
 * 
 * SRE, Mon Mar 31 08:59:07 2008 [Janelia]
 * SVN $Id$
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sqio.h"

static char banner[] = "concatenate sequence data into one binary file";
static char usage1[] = "[options] <seqfile> <concatdb>";

#define ALPH_OPTS "--rna,--dna,--amino" /* toggle group, alphabet type options          */

static ESL_OPTIONS options[] = {
  /* name         type           default   env range togs  reqs  incomp      help                                      docgroup */
  { "-h",         eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL,      NULL, "help; show brief info on version and usage",          1 },
  { "--informat", eslARG_STRING,  FALSE, NULL, NULL, NULL, NULL,      NULL, "specify that input file is in format <s>",            1 },
  { "--rna",      eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, ALPH_OPTS, "specify that <seqfile> contains RNA sequence",        1 },
  { "--dna",      eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, ALPH_OPTS, "specify that <seqfile> contains DNA sequence",        1 },
  { "--amino",    eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, ALPH_OPTS, "specify that <seqfile> contains protein sequence",    1 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

static void cmdline_failure(char *argv0, char *format, ...);
static void cmdline_help   (char *argv0, ESL_GETOPTS *go);

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go        = NULL;
  char           *seqfile   = NULL;
  char           *dbfile    = NULL;
  ESL_SQFILE     *sqfp      = NULL;
  FILE           *dbfp      = NULL;
  int             infmt     = eslSQFILE_UNKNOWN;
  int             alphatype = eslUNKNOWN;
  ESL_ALPHABET   *abc       = NULL;
  ESL_SQ         *sq        = NULL;
  int             status    = eslOK;
  long long       n         = 0;
  int             nseq      = 0;

  /* Parse command line */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)    cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK)    cmdline_failure(argv[0], "Error in app configuration: %s\n",   go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") )                      cmdline_help(argv[0], go);

  if (esl_opt_ArgNumber(go) != 2)                         cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");
  seqfile = esl_opt_GetArg(go, 1);
  dbfile  = esl_opt_GetArg(go, 2);

  if (esl_opt_GetString(go, "--informat") != NULL) {
    infmt = esl_sqio_FormatCode(esl_opt_GetString(go, "--informat"));
    if (infmt == eslSQFILE_UNKNOWN) esl_fatal("%s is not a valid input sequence file format for --informat"); 
  }

  /* Open the sequence file in digital mode */
  status = esl_sqfile_Open(seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file %s", seqfile);
  else if (status == eslEFORMAT)   esl_fatal("Format of seqfile %s unrecognized.", seqfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  if      (esl_opt_GetBoolean(go, "--rna"))   alphatype = eslRNA;
  else if (esl_opt_GetBoolean(go, "--dna"))   alphatype = eslDNA;
  else if (esl_opt_GetBoolean(go, "--amino")) alphatype = eslAMINO;
  else {
    status = esl_sqfile_GuessAlphabet(sqfp, &alphatype);
    if      (status == eslEAMBIGUOUS) esl_fatal("Couldn't guess alphabet from first sequence in %s", seqfile);
    else if (status == eslEFORMAT)    esl_fatal("Sequence file parse error, line %d of file %s:\n%s\n",
					       sqfp->linenumber, seqfile, sqfp->errbuf);
    else if (status == eslENODATA)    esl_fatal("Sequence file %s contains no data?", seqfile);
    else if (status != eslOK)         esl_fatal("Failed to guess alphabet (error code %d)\n", status);
  }
  abc = esl_alphabet_Create(alphatype);
  sq  = esl_sq_CreateDigital(abc);

  /* Open the output database file in binary mode */
  if ((dbfp = fopen(dbfile, "wb")) == NULL) esl_fatal("Failed to open dbfile %s for writing\n", dbfile);

  /* Read each sequence; output digital to dbfp. */
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      fwrite(sq->dsq+1, sizeof(char), sq->n, dbfp); /* dsq+1 to skip SENTINELs */
      n += sq->n;
      nseq++;
      esl_sq_Reuse(sq);
    }
  
  printf("%d residues (from %d seqs) concatenated and written (digital mode) to %s\n", 
	 n, nseq, dbfile);

  fclose(dbfp);
  esl_sqfile_Close(sqfp);
  esl_alphabet_Destroy(abc);
  esl_sq_Destroy(sq);
  esl_getopts_Destroy(go);
  return 0;
}

static void
cmdline_failure(char *argv0, char *format, ...)
{
  va_list argp;

  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  esl_usage(stdout, argv0, usage1);
  printf("\nTo see more help on available options, do %s -h\n\n", argv0);
  exit(1);
}

static void
cmdline_help(char *argv0, ESL_GETOPTS *go) 
{
  esl_banner(stdout, argv0, banner);
  esl_usage (stdout, argv0, usage1);
  puts("\n where general options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  exit(0);
}

