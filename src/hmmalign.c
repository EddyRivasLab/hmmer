/* hmmalign:  align sequences to a profile HMM
 * 
 * SRE, Wed Oct 22 10:33:15 2008 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

#define ALPHOPTS "--amino,--dna,--rna"                         /* Exclusive options for alphabet choice */

static ESL_OPTIONS options[] = {
  /* name             type        default      env  range   toggles   reqs  incomp               help                                          docgroup*/
  { "-h",          eslARG_NONE,     FALSE,     NULL, NULL,   NULL,    NULL,  NULL, "show brief help on version and usage",                        1 },
  { "-o",          eslARG_OUTFILE,   NULL,     NULL, NULL,   NULL,    NULL,  NULL, "output alignment to file <f>, not stdout",                    1 },
  { "-q",          eslARG_NONE,     FALSE,     NULL, NULL,   NULL,    NULL,  NULL, "quiet: suppress banner and informational output",             1 },

  //  { "--mapali",    eslARG_INFILE,    NULL,     NULL, NULL,   NULL,    NULL,  NULL, "include alignment in file <f> using map in HMM",              2 },
  { "--amino",     eslARG_NONE,     FALSE,     NULL, NULL, ALPHOPTS,  NULL,  NULL, "assert <seqfile>, <hmmfile> both protein: no autodetection",  2 },
  { "--dna",       eslARG_NONE,     FALSE,     NULL, NULL, ALPHOPTS,  NULL,  NULL, "assert <seqfile>, <hmmfile> both DNA: no autodetection",      2 },
  { "--rna",       eslARG_NONE,     FALSE,     NULL, NULL, ALPHOPTS,  NULL,  NULL, "assert <seqfile>, <hmmfile> both RNA: no autodetection",      2 },
  { "--informat",  eslARG_STRING,    NULL,     NULL, NULL,   NULL,    NULL,  NULL, "assert input <seqfile> is in format <s>: no autodetection",   2 },
  { "--outformat", eslARG_STRING, "Stockholm", NULL, NULL,   NULL,    NULL,  NULL, "output alignment in format <s>",                              2 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "align sequences to a profile HMM";

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
  puts("\n Less common options are:");
  esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
  puts("\nSequence input formats include: FASTA, EMBL, Genbank, Uniprot");
  puts("Alignment output formats include: Stockholm, Pfam, A2M, PSIBLAST\n");
  exit(0);
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = NULL;	/* application configuration       */
  char         *hmmfile = NULL;	/* HMM file name                   */
  char         *seqfile = NULL; /* sequence file name              */
  int           infmt   = eslSQFILE_UNKNOWN;
  int           outfmt  = eslMSAFILE_STOCKHOLM;
  P7_HMMFILE   *hfp     = NULL;	/* open HMM file                   */
  ESL_SQFILE   *sqfp    = NULL;	/* open sequence file              */
  ESL_SQ      **sq      = NULL;	/* array of sequences              */
  void         *p       = NULL;	/* tmp ptr for reallocation        */
  int           nseq    = 0;
  ESL_ALPHABET *abc     = NULL;	/* alphabet (set from the HMM file)*/
  P7_BG        *bg      = NULL;
  P7_HMM       *hmm     = NULL;
  P7_PROFILE   *gm      = NULL;
  P7_OPROFILE  *om      = NULL;
  P7_OMX       *oxf     = NULL;	/* optimized Forward matrix        */
  P7_OMX       *oxb     = NULL;	/* optimized Backward matrix       */
  P7_TRACE    **tr      = NULL;	/* array of tracebacks             */
  ESL_MSA      *msa     = NULL;	/* resulting multiple alignment    */
  int           idx;		/* counter over seqs, traces       */
  float         fwdsc;		/* Forward score                   */
  float         oasc;		/* optimal accuracy score          */
  int           status;		/* easel/hmmer return code         */

  /* Parse the command line
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in configuration: %s\n",       go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") )                   cmdline_help   (argv[0], go);
  if (esl_opt_ArgNumber(go) != 2)                      cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        

  hmmfile = esl_opt_GetArg(go, 1);
  seqfile = esl_opt_GetArg(go, 2);

  /* If caller declared an input format, decode it 
   */
  if (! esl_opt_IsDefault(go, "--informat")) {
    infmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat"));
    if (infmt == eslSQFILE_UNKNOWN) cmdline_failure(argv[0], "%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--informat"));
  }

  /* Determine output alignment file format */
  outfmt = esl_msa_EncodeFormat(esl_opt_GetString(go, "--outformat"));
  if (outfmt == eslMSAFILE_UNKNOWN)    cmdline_failure(argv[0], "%s is not a recognized output MSA file format\n", esl_opt_GetString(go, "--outformat"));


  /* If caller forced an alphabet on us, create the one the caller wants  
   */
  if      (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);

  /* Read one HMM, and make sure there's only one.
   */
  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if      (status == eslENOTFOUND) cmdline_failure(argv[0], "Failed to open HMM file %s for reading.\n",                   hmmfile);
  else if (status == eslEFORMAT)   cmdline_failure(argv[0], "File %s does not appear to be in a recognized HMM format.\n", hmmfile);
  else if (status != eslOK)        cmdline_failure(argv[0], "Unexpected error %d in opening HMM file %s.\n", status,       hmmfile);  

  status = p7_hmmfile_Read(hfp, &abc, &hmm);
  if      (status == eslEFORMAT)   cmdline_failure(argv[0], "bad file format in HMM file %s:\n       %s\n",     hfp->fname, hfp->errbuf);
  else if (status == eslEINCOMPAT) cmdline_failure(argv[0], "HMM in %s is not in the expected %s alphabet\n",   hfp->fname, esl_abc_DecodeType(abc->type));
  else if (status == eslEOF)       cmdline_failure(argv[0], "Empty HMM file %s? No HMM data found.\n",          hfp->fname);
  else if (status != eslOK)        cmdline_failure(argv[0], "Unexpected error in reading HMMs from %s\n",       hfp->fname);

  status = p7_hmmfile_Read(hfp, &abc, NULL);
  if      (status != eslEOF)       cmdline_failure(argv[0], "HMM file %s does not contain just one HMM\n", hfp->fname);
  p7_hmmfile_Close(hfp);


  /* Read digital sequences into an array.
   */
  status = esl_sqfile_OpenDigital(abc, seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          seqfile);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            seqfile);
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, seqfile);

  ESL_ALLOC(sq, sizeof(ESL_SQ *));
  sq[0] = esl_sq_CreateDigital(abc);
  nseq = 0;
  while ((status = esl_sqio_Read(sqfp, sq[nseq])) == eslOK)
    {
      nseq++;
      ESL_RALLOC(sq, p, sizeof(ESL_SQ *) * (nseq+1));
      sq[nseq] = esl_sq_CreateDigital(abc);
    }
  if      (status == eslEFORMAT) p7_Fail("Parse failed, line %ld, file %s\n%s", (long) sqfp->linenumber, sqfp->filename, sqfp->errbuf);
  else if (status != eslEOF)     p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  esl_sqfile_Close(sqfp);


  /* Remaining initializations, including trace array allocation
   */
  ESL_ALLOC(tr, sizeof(P7_TRACE *) * nseq);
  for (idx = 0; idx < nseq; idx++)
    tr[idx] = p7_trace_CreateWithPP();

  bg = p7_bg_Create(abc);
  gm = p7_profile_Create (hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);

  p7_ProfileConfig(hmm, bg, gm, sq[0]->n, p7_UNILOCAL);
  p7_oprofile_Convert(gm, om); 

  oxf = p7_omx_Create(hmm->M, sq[0]->n, sq[0]->n);
  oxb = p7_omx_Create(hmm->M, sq[0]->n, sq[0]->n);	
  

  /* Collect an OA trace for each sequence 
   */
  for (idx = 0; idx < nseq; idx++)
    {
      p7_omx_GrowTo(oxf, hmm->M, sq[idx]->n, sq[idx]->n);
      p7_omx_GrowTo(oxb, hmm->M, sq[idx]->n, sq[idx]->n);
      
      p7_oprofile_ReconfigLength(om, sq[idx]->n);

      p7_Forward (sq[idx]->dsq, sq[idx]->n, om,      oxf, &fwdsc);
      p7_Backward(sq[idx]->dsq, sq[idx]->n, om, oxf, oxb, NULL);

      status = p7_Decoding(om, oxf, oxb, oxb);      /* <oxb> is now overwritten with post probabilities     */
      if (status == eslERANGE) return eslFAIL;      /* rare: numeric overflow, generally on repetitive garbage [J3/119-212] */

      p7_OptimalAccuracy(om, oxb, oxf, &oasc);      /* <oxf> is now overwritten with OA scores              */
      p7_OATrace        (om, oxb, oxf, tr[idx]);    /* tr[idx] is now an OA traceback for seq #idx          */
 
      p7_omx_Reuse(oxf);
      p7_omx_Reuse(oxb);
    }

#if 0
  for (idx = 0; idx < nseq; idx++)
    p7_trace_Dump(stdout, tr[idx], gm, sq[idx]->dsq);
#endif

  p7_Traces2Alignment(sq, tr, nseq, om->M, p7_DEFAULT, &msa);

  esl_msa_Write(stdout, msa, outfmt);

  for (idx = 0; idx <= nseq; idx++) esl_sq_Destroy(sq[idx]);    /* including sq[nseq] because we overallocated */
  for (idx = 0; idx <  nseq; idx++) p7_trace_Destroy(tr[idx]); 
  free(sq);
  free(tr);
  esl_msa_Destroy(msa);
  p7_bg_Destroy(bg);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
  p7_omx_Destroy(oxf);
  p7_omx_Destroy(oxb);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;

 ERROR:
  return status;
}


/*****************************************************************
 * @LICENSE@
 *****************************************************************/