/* hmmemit: sample sequence(s) from a profile HMM.
 * 
 * SRE, Tue Jan  9 13:22:53 2007 [Janelia] [Verdi, Requiem]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#define MODEOPTS "--local,--unilocal,--glocal,--uniglocal"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",          eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",                1 },
  { "-c",          eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    "-p", "emit simple consensus sequence",                      1 },
  { "-o",          eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,    NULL, "send sequence output to file <f>, not stdout",        1 },
  { "-p",          eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    "-c", "sample from profile, not core model",                 1 },
  { "-N",          eslARG_INT,      "1", NULL, "n>0",     NULL,      NULL,    "-c", "number of seqs to sample",                            1 },
/* options for profile emission, with -p  */
  { "-L",          eslARG_INT,    "400", NULL, NULL,      NULL,      "-p",    "-c", "set expected length from profile to <n>",             2 },
  { "--local",     eslARG_NONE,"default",NULL, NULL,    MODEOPTS,    "-p",    "-c", "configure profile in local mode",                     2 }, 
  { "--unilocal",  eslARG_NONE,  FALSE,  NULL, NULL,    MODEOPTS,    "-p",    "-c", "configure profile in unilocal mode",                  2 }, 
  { "--glocal",    eslARG_NONE,  FALSE,  NULL, NULL,    MODEOPTS,    "-p",    "-c", "configure profile in glocal mode",                    2 }, 
  { "--uniglocal", eslARG_NONE,  FALSE,  NULL, NULL,    MODEOPTS,    "-p",    "-c", "configure profile in glocal mode",                    2 }, 
/* other options */
  { "--seed",      eslARG_INT,      "0", NULL, "n>=0",    NULL,      NULL,    NULL, "set RNG seed to <n>",                                 3 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile (single)>";
static char banner[] = "sample sequence(s) from a profile HMM";

static void cmdline_failure(char *argv0, char *format, ...);
static void cmdline_help(char *argv0, ESL_GETOPTS *go);


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go         = NULL;             /* command line processing                 */
  ESL_ALPHABET    *abc        = NULL;             /* sequence alphabet                       */
  ESL_RANDOMNESS  *r          = NULL;             /* source of randomness                    */
  char            *hmmfile    = NULL;             /* file to read HMM(s) from                */
  P7_HMMFILE      *hfp        = NULL;             /* open hmmfile                            */
  P7_HMM          *hmm        = NULL;             /* HMM to emit from                        */
  P7_PROFILE      *gm         = NULL;             /* profile HMM (scores)                    */
  P7_BG           *bg         = NULL;             /* null model                              */
  P7_TRACE        *tr         = NULL;             /* sampled trace                           */
  ESL_SQ          *sq         = NULL;             /* sampled digital sequence                */
  FILE            *ofp        = NULL;	          /* output stream                           */
  int              N          = 0;  	          /* how many sequences to emit              */
  int              do_profile = FALSE; 	          /* TRUE to emit from profile, not core     */
  int              L          = 0;	          /* expected length from a profile          */
  int              outfmt     = eslSQFILE_FASTA;
  int              mode       = p7_NO_MODE;
  int              nseq;
  int              status;	      

  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in configuration: %s\n",       go->errbuf);
  if (esl_opt_GetBoolean(go, "-h"))                    cmdline_help   (argv[0], go);      
  if (esl_opt_ArgNumber(go) != 1)                      cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        

  if ((hmmfile = esl_opt_GetArg(go, 1)) == NULL)       cmdline_failure(argv[0], "Failed to get <hmmfile> on cmdline: %s\n", go->errbuf);

  do_profile   = esl_opt_GetBoolean(go, "-p");
  N            = esl_opt_GetInteger(go, "-N");
  L            = esl_opt_GetInteger(go, "-L");

  if      (esl_opt_GetBoolean(go, "--local"))     mode = p7_LOCAL;
  else if (esl_opt_GetBoolean(go, "--unilocal"))  mode = p7_UNILOCAL;
  else if (esl_opt_GetBoolean(go, "--glocal"))    mode = p7_GLOCAL;
  else if (esl_opt_GetBoolean(go, "--uniglocal")) mode = p7_UNIGLOCAL;

  if ( esl_opt_IsOn(go, "-o"))
    {
      ofp = fopen(esl_opt_GetString(go, "-o"), "w");
      if (ofp == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "-o"));
    }
  else ofp = stdout;


  /*****************************************************************
   * Initializations, including opening the HMM file
   *****************************************************************/

  r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));

  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open HMM file %s for reading.\n",                   hmmfile);
  else if (status == eslEFORMAT)   p7_Fail("File %s does not appear to be in a recognized HMM format.\n", hmmfile);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n",       status, hmmfile);  

  else if (status != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", hmmfile);
    
  if ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslOK) {
    if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
    else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s", hmmfile);
    else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets", hmmfile);
    else                             esl_fatal("Unexpected error in reading HMMs");
  }
   
  /* init sq (need to know abc to do this */
  if (sq == NULL) sq = esl_sq_CreateDigital(abc);
  if (sq == NULL) esl_fatal("Failed to allocated sequence");

  if (esl_opt_GetBoolean(go, "-c")) 
    {
      if (p7_emit_SimpleConsensus(hmm, sq)                 != eslOK) esl_fatal("failed to create simple consensus seq");
      if (esl_sq_FormatName(sq, "%s-consensus", hmm->name) != eslOK) esl_fatal("Failed to set sequence name");
      if (esl_sqio_Write(ofp, sq, eslSQFILE_FASTA, FALSE)  != eslOK) esl_fatal("Failed to write sequence");
    }
  else
    {
      if ((tr = p7_trace_Create())              == NULL)  esl_fatal("Failed to allocate trace");
      if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");
      if ((gm = p7_profile_Create(hmm->M, abc)) == NULL)  esl_fatal("failed to create profile");

      if (p7_ProfileConfig(hmm, bg, gm, L, mode)         != eslOK) esl_fatal("failed to configure profile");
      if (p7_bg_SetLength(bg, L)                         != eslOK) esl_fatal("failed to reconfig null model length");
      if (p7_hmm_Validate    (hmm, NULL, 0.0001)         != eslOK) esl_fatal("whoops, HMM is bad!");
      if (p7_profile_Validate(gm,  NULL, 0.0001)         != eslOK) esl_fatal("whoops, profile is bad!");

      for (nseq = 1; nseq <= N; nseq++) 
	{
	  if (do_profile) {
	    status = p7_ProfileEmit(r, hmm, gm, bg, sq, tr);
	    if (status != eslOK) esl_fatal("Failed to emit sequence from hmm\n");
	  } else {
	    status = p7_CoreEmit(r, hmm, sq, tr);
	    if (status != eslOK) esl_fatal("Failed to emit sequence from hmm\n");
	  }
      
	  status = esl_sq_FormatName(sq, "%s-sample%d", hmm->name, nseq);
	  if (status != eslOK) esl_fatal("Failed to set sequence name\n");

	  status = esl_sqio_Write(ofp, sq, outfmt, FALSE);
	  if (status != eslOK) esl_fatal("Failed to write sequence\n");
	}

      p7_profile_Destroy(gm);
      p7_bg_Destroy(bg);
      p7_trace_Destroy(tr);
    }

  if (esl_opt_IsOn(go, "-o")) { fclose(ofp); }
  esl_sq_Destroy(sq);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  p7_hmmfile_Close(hfp);
  p7_hmm_Destroy(hmm);
  return eslOK;
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
  puts("\n Common options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  puts("\n Options controlling emission from profiles (with -p):");
  esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
  puts("\n Other options::");
  esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
  exit(0);
}
