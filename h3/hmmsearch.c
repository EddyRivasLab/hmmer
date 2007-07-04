/* main() for scoring a profile HMM against a sequence database.
 */

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sqio.h"
#include "esl_histogram.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",   1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "hmmsearch <query hmmfile> <target seqfile> [-options]";


int
main(int argc, char **argv)
{
  ESL_ALPHABET    *abc     = NULL;     /* sequence alphabet                       */
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                 */
  P7_HMM          *hmm     = NULL;     /* query HMM                               */
  P7_PROFILE      *gm      = NULL;     /* profile HMM                             */
  P7_BG           *bg      = NULL;     /* null model                              */
  ESL_SQ          *sq      = NULL;     /* target sequence                         */
  ESL_DSQ         *dsq     = NULL;     /* digitized target sequence               */
  P7_TRACE        *tr      = NULL;     /* trace of hmm aligned to sq              */
  P7_GMX          *mx      = NULL;     /* DP matrix                               */
  int              vsc, fsc;	       /* Viterbi, Forward scores (SILO)          */
  ESL_HISTOGRAM   *h       = esl_histogram_CreateFull(-50., 50., 0.5);

  char            *hmmfile;            /* file to read HMM(s) from                */
  char            *seqfile;	       /* file to read sequence(s) from           */
  int              format  = eslSQFILE_FASTA;
  P7_HMMFILE      *hfp     = NULL;    
  ESL_SQFILE      *sqfp    = NULL;
  int              hstatus, sstatus;



  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("Failed to parse command line: %s", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("Failed to parse command line: %s", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    puts(usage);
    puts("\noptions are:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
    return eslOK;
  }
  if (esl_opt_ArgNumber(go) != 2) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return eslFAIL;
  }
  hmmfile = esl_opt_GetArg(go, 1);
  if (hmmfile == NULL) esl_fatal("Failed to get <hmmfile> argument on command line");
  seqfile = esl_opt_GetArg(go, 2);
  if (seqfile == NULL) esl_fatal("Failed to get <seqfile> argument on command line");
  
  /*****************************************************************
   * Initializations, including opening the input files
   *****************************************************************/

  abc = esl_alphabet_Create(eslAMINO);
  tr =  p7_trace_Create(256);
  mx  = p7_gmx_Create(200, 400); /* initial alloc is for a DP problem M=200,L=400 */

  hstatus = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if (hstatus == eslENOTFOUND) esl_fatal("Failed to open hmm file %s for reading.\n", hmmfile);
  else if (hstatus != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", hmmfile);

  sstatus = esl_sqfile_Open(seqfile, format, NULL, &sqfp); /* NULL = env variable name for dir the seqfile might be in */
  if      (sstatus == eslENOTFOUND) esl_fatal("Failed to open sequence file %s for reading\n", seqfile);
  else if (sstatus == eslEFORMAT)   esl_fatal("Sequence file %s is empty or misformatted\n",   seqfile);
  else if (sstatus == eslEINVAL)    esl_fatal("Can't autodetect format of a stdin or .gz seqfile");
  else if (sstatus != eslOK)        esl_fatal("Unexpected error %d opening sequence file %s\n", sstatus, seqfile);

  sq = esl_sq_Create();

  while ( (hstatus = p7_hmmfile_Read(hfp, &abc, &hmm)) == eslOK) 
    {
      /* Configure the profile (wait 'til we see sequences to config length) */
      bg = p7_bg_Create(abc);
      gm = p7_profile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, p7_UNILOCAL);

      while ( (sstatus = esl_sqio_Read(sqfp, sq)) == eslOK)
	{
	  esl_abc_CreateDsq(abc, sq->seq, &dsq); /* replace w/ direct digital read when it's available */
	  p7_gmx_GrowTo(mx, hmm->M, sq->n); /* realloc DP matrix as needed */
	  p7_ReconfigLength(gm, sq->n);
	  p7_bg_SetLength(bg, sq->n);

	  p7_GViterbi(dsq, sq->n, gm, mx, &vsc);
	  p7_GTrace  (dsq, sq->n, gm, mx, tr);
	  p7_GForward(dsq, sq->n, gm, mx, &fsc);

	  printf("%15s  %6.2f   %6.2f\n", sq->name,
		 (((double) vsc / p7_INTSCALE) - sq->n * log(bg->p1) - log(1.-bg->p1)) / eslCONST_LOG2,
		 (((double) fsc / p7_INTSCALE) - sq->n * log(bg->p1) - log(1.-bg->p1)) / eslCONST_LOG2);

	  free(dsq);
	  esl_sq_Reuse(sq);
	  p7_trace_Reuse(tr);
	}
      if (sstatus != eslEOF) 
	esl_fatal("Sequence file %s has a format problem: read failed at line %d:\n%s\n",
		  seqfile, sqfp->linenumber, sqfp->errbuf);

      p7_profile_Destroy(gm);  
      p7_bg_Destroy(bg);
      p7_hmm_Destroy(hmm);
    }
  if      (hstatus == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
  else if (hstatus == eslEFORMAT)   esl_fatal("bad file format in HMM file %s", hmmfile);
  else if (hstatus == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets", hmmfile);
  else if (hstatus != eslEOF)       esl_fatal("Unexpected error in reading HMMs");

  p7_gmx_Destroy(mx);
  p7_trace_Destroy(tr);
  p7_hmmfile_Close(hfp);
  esl_histogram_Destroy(h);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  esl_getopts_Destroy(go);
  esl_alphabet_Destroy(abc);
  return 0;
}
