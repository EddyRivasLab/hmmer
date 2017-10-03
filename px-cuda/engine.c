#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dsqdata.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,      "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example driver for the engine";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  ESL_DSQDATA    *dd      = NULL;
  ESL_DSQDATA_CHUNK *chu  = NULL;
  P7_ENGINE      *eng     = NULL;
  int             i;
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

  /* Open sequence file */
  if (( status = esl_dsqdata_Open(&abc, seqfile, 1, &dd) ) != eslOK)
    esl_fatal("open failed");
 
  /* Create the comparison engine */
  eng   = p7_engine_Create(abc, /*params=*/NULL, /*stats=*/NULL, gm->M, /*L_hint=*/400);

  /* For each sequence in <sqfile>: */
  while ( esl_dsqdata_Read(dd, &chu) == eslOK)
    { 
      //printf("working on %s (len=%d)... ", sq->name, (int) sq->n);

      for (i = 0; i < chu->N; i++)
	{
	  p7_bg_SetLength     (bg, (int) chu->L[i]);
	  p7_oprofile_ReconfigLength(om, (int) chu->L[i]);

	  printf("seq %d %s\n", chu->i0+i, chu->name[i]);

	  status = p7_engine_Overthruster(eng, chu->dsq[i], chu->L[i], om, bg);
	  if      (status == eslFAIL) { 
	    p7_engine_Reuse(eng);
	    continue; 
	  }
	  else if (status != eslOK)   p7_Fail("overthruster failed with code %d\n", status);

	  p7_profile_SetLength(gm, (int) chu->L[i]);
	  status = p7_engine_Main(eng, chu->dsq[i], (int) chu->L[i], gm);


	  //p7_trace_DumpAnnotated(stdout, eng->tr, gm, sq->dsq);
	  p7_engine_Reuse(eng);
	}
      esl_dsqdata_Recycle(dd, chu);
    }
  
  esl_dsqdata_Close(dd);
  p7_engine_Destroy(eng);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
