#include <string.h>
#include "easel.h"
#include "esl_dsqdata.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { "-h",        eslARG_NONE,  FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",  0 },
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
  P7_ENGINE      *eng     = NULL;
  ESL_DSQDATA_CHUNK *chu = NULL;
  int             ncore   = 1;
  int  i;
  int             status;
  char outfile_name[255];
  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  strcpy(outfile_name, "uniprot_trembl_samples/");
  strcat(outfile_name, hmm->name);
  strcat(outfile_name, ".hits");
  FILE *outfile = fopen(outfile_name, "w");

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create (hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config   (gm, hmm, bg);
  p7_oprofile_Convert (gm, om);

  p7_bg_SetFilter(bg, om->M, om->compo);

  uint64_t sequence_id = 0;
  uint64_t num_hits = 0;
  /* Open sequence database */
  status = esl_dsqdata_Open(&abc, seqfile, ncore, &dd);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open dsqdata files:\n  %s",    dd->errbuf);
  else if (status == eslEFORMAT)   p7_Fail("Format problem in dsqdata files:\n  %s", dd->errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error in opening dsqdata (code %d)", status);

  eng = p7_engine_Create(abc, NULL, NULL, gm->M, 400);

  while (( status = esl_dsqdata_Read(dd, &chu)) == eslOK)  
    {
      for (i = 0; i < chu->N; i++)
	{
	  p7_bg_SetLength(bg, (int) chu->L[i]);            // TODO: remove need for cast
	  p7_oprofile_ReconfigLength(om, (int) chu->L[i]); //         (ditto)
	  
	  //	  printf("seq %d %s\n", chu->i0+i, chu->name[i]);

	  status = p7_engine_Overthruster(eng, chu->dsq[i], (int) chu->L[i], om, bg);  
	  if (status == eslFAIL) { 
	    p7_engine_Reuse(eng);
	    sequence_id++;
	    continue;
	  }

	  p7_profile_SetLength(gm, (int) chu->L[i]);
	  status = p7_engine_Main(eng, chu->dsq[i], (int) chu->L[i], gm);
	  
	  // anything that goes to main is a hit for now
	  fprintf(outfile, "%lu %s %s %s\n", sequence_id, chu->name[i], chu->acc[i], chu->desc[i]);
	  
	  p7_engine_Reuse(eng);
	  sequence_id++;
	  num_hits++;
	}
      esl_dsqdata_Recycle(dd, chu);
    }

  fprintf(outfile, "%lu hits found\n", num_hits);
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





