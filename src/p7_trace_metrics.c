#include "p7_config.h"
#include "hmmer.h"

/* 
 * 0. Complete state path, st/k/i all have to match.
 * 1. Aligned residues: if i is emitted by M/I, st/k have to match.
 * 2. Region: if i is emitted by M/I, i must be emitted by M/I in test.
 * 3. Boundary: if B/E; residue position must match.
 *
 * Collect tp,fp,fn. tp = reference, test agree
 *                   fp = test has something reference doesn't
 *                   fn = reference has something test doesn't
 */
int
p7_trace_metrics(P7_TRACE *reftr, P7_TRACE *testtr, P7_PROFILE *gm, ESL_DSQ *dsq)
{
  int zr, zt;			/* position in reference, test trace */
  int state_tp,  state_fp,  state_fn;
  int align_tp,  align_fp,  align_fn;
  int region_tp, region_fp, region_fn;
  int edge_tp,   edge_fp,   edge_fn;

  state_tp  = state_fp  = state_fn  = 0;
  align_tp  = align_fp  = align_fn  = 0;
  region_tp = region_fp = region_fn = 0;
  edge_tp   = edge_fp   = edge_fn   = 0;

  for (zr = zt = 0; zr < reftr->N; zr++)
    {
      /* Increment zt: align test to the current position of the reference  */
      if (reftr->i[zr])	 /* if reference emits i: move test to where it emitted i too */
	{
	  while (testtr->i[zt] < reftr->i[zr])
	    { /* everything we must skip over in the test is a) nonemitting and b) a false positive */

	      printf("%5s %2s %5s %5s  %c     %5d %2s %5d %5d  %c\n",
		     "", "", "", "", ' ',     
		     zt, p7_trace_DecodeStatetype(testtr->st[zt]), testtr->k[zt], testtr->i[zt], (testtr->i[zt] ? gm->abc->sym[dsq[testtr->i[zt]]] : '-'));

	      if (testtr->st[zt] == p7T_B || testtr->st[zt] == p7T_E) edge_fp++;
	      state_fp++;
	      zt++;
	    }
	}
      else 
	{ /* if reference is mute: either find a corresponding mute state in test, or halt and wait on T or next emitter */		
	  while (testtr->i[zt] == 0 && testtr->st[zt] != reftr->st[zr] && zt < testtr->N-1)
	    {  /* again, anything extra we see in the test trace is a false positive */
	      printf("%5s %2s %5s %5s  %c     %5d %2s %5d %5d  %c\n",
		     "", "", "", "", ' ',     
		     zt, p7_trace_DecodeStatetype(testtr->st[zt]), testtr->k[zt], testtr->i[zt], (testtr->i[zt] ? gm->abc->sym[dsq[testtr->i[zt]]] : '-'));

	      if (testtr->st[zt] == p7T_B || testtr->st[zt] == p7T_E) edge_fp++;
	      state_fp++;
	      zt++;
	    }
	}

      /* Now we can be in one of the following states:
       *   if reference trace emits: then test trace also emits, the
       *     same position i. (though not necessarily w/ same state, or k)
       *   
       *   if reference trace is on a mute state: test trace may emit or not.
       *      if test trace emits: then its position i is the next one that
       *          the reference trace will reach next, and it
       *          is paused here, waiting for reference to get thru its mute
       *          states and catch up.
       *      if not: then either the reference and the test are on the same
       *          mute state;
       *              or the test trace is waiting at its terminal T state.
       *
       *   corollary: if test trace is on a mute state, it must either match 
       *      the reference mute state, or it must be done (on T).
       *
       * so zr is always going to advance; but under certain conditions, zt
       * is waiting: reference trace does not emit, and test either emits or 
       * is on T.
       */
      if (reftr->i[zr] == 0 && reftr->st[zr] != testtr->st[zt])
	{
	  printf("%5d %2s %5d %5d  %c     %5s %2s %5s %5s  %c\n",
		 zr, p7_trace_DecodeStatetype(reftr->st[zr]), reftr->k[zr], reftr->i[zr], (reftr->i[zr] ? gm->abc->sym[dsq[reftr->i[zr]]] : '-'),
		 "", "", "", "", ' ');
	}
      else
	{
	  printf("%5d %2s %5d %5d  %c     %5d %2s %5d %5d  %c\n",
		 zr, p7_trace_DecodeStatetype( reftr->st[zr]),  reftr->k[zr],  reftr->i[zr], ( reftr->i[zr] ? gm->abc->sym[dsq[ reftr->i[zr]]] : '-'),
		 zt, p7_trace_DecodeStatetype(testtr->st[zt]), testtr->k[zt], testtr->i[zt], (testtr->i[zt] ? gm->abc->sym[dsq[testtr->i[zt]]] : '-'));
	}
      /* All-states comparison */
      if (reftr->st[zr] == testtr->st[zt] && reftr->i[zr]  == testtr->i[zt]  && reftr->k[zr]  == testtr->k[zt]) state_tp++;
      else {
	state_fn++;
	if (reftr->i[zr]) state_fp++; /* checking i[z] makes sure we only count each i once */
      }

      /* Alignment comparison */
      if (reftr->i[zr])		/* this makes sure we only evaluate alignment of each i once. we can go through zr loop more than once w/ testtr->i[zt] on the same i */
	{
	  if ( p7_trace_IsMain(reftr->st[zr]))
	    {
	      if ( reftr->k[zr] == testtr->k[zt] && 
		   ( (p7_trace_IsM(reftr->st[zr]) && p7_trace_IsM(testtr->st[zt])) ||
		     (p7_trace_IsI(reftr->st[zr]) && p7_trace_IsI(testtr->st[zt]))) )
		align_tp++;
	      else {
		align_fn++;	
		if (p7_trace_IsMain(testtr->st[zt])) align_fp++; /* not only did testtr fail to correctly align i, it aligned it somewhere else that's wrong (a false pos, in addition to being a false neg) */
	      }
	    }
	  else { 			/* reference emits it w/ NJC; if test puts it in model, it's a fp */
	    if (p7_trace_IsMain(testtr->st[zt])) align_fp++;
	  }
	}
      
      /* Region comparison */
      if (reftr->i[zr]) {
	if      (p7_trace_IsMain(reftr->st[zr]) && p7_trace_IsMain(testtr->st[zt])) region_tp++;
	else if (p7_trace_IsMain(reftr->st[zr]))                    	            region_fn++;
	else if (p7_trace_IsMain(testtr->st[zt]))	                            region_fp++;
      }

      /* Edge comparison */
      if (reftr->st[zr] == p7T_B) {
	if (testtr->st[zt] == p7T_B) edge_tp++;
	else                         edge_fn++;
      }
      if (reftr->st[zr] == p7T_E) {
	if (testtr->st[zt] == p7T_E) edge_tp++;
	else                         edge_fn++;
      }
      
      /* Advance zt unless it's supposed to wait */
      if (reftr->i[zr] || ( zt < testtr->N-1 && testtr->i[zt] == 0)) zt++;
    }

  printf("  which    TP   FN   FP\n");
  printf("-------- ---- ---- ----\n");
  printf("state    %4d %4d %4d\n", state_tp,  state_fn,  state_fp);
  printf("align    %4d %4d %4d\n", align_tp,  align_fn,  align_fp);
  printf("region   %4d %4d %4d\n", region_tp, region_fn, region_fp);
  printf("edge     %4d %4d %4d\n", edge_tp,   edge_fn,   edge_fp);
  return eslOK;
}

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "configured mean seq length for profile",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "example of calculating accuracy metrics for comparison of two traces";


/* Test: emit two traces from the same profile. 
 *       compare them.
 * 
 * cc -I . -I ../easel -L. -L../easel p7_trace_benchmark.c -lhmmer -leasel -lm
 */
int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  int             L       = esl_opt_GetInteger(go, "-L");
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_TRACE       *reftr   = p7_trace_Create();
  P7_REFMX       *vmx     = NULL;
  P7_TRACE       *testtr  = p7_trace_Create();
  ESL_SQ         *sq      = NULL;
  char            errbuf[eslERRBUFSIZE];
  int             i;
  int             status;


  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",                       status, hmmfile, errbuf);  

  status = p7_hmmfile_Read(hfp, &abc, &hmm);
  if      (status == eslEFORMAT)   p7_Fail("Bad file format in HMM file %s:\n%s\n",          hfp->fname, hfp->errbuf);
  else if (status == eslEINCOMPAT) p7_Fail("HMM in %s is not in the expected %s alphabet\n", hfp->fname, esl_abc_DecodeType(abc->type));
  else if (status == eslEOF)       p7_Fail("Empty HMM file %s? No HMM data found.\n",        hfp->fname);
  else if (status != eslOK)        p7_Fail("Unexpected error in reading HMMs from %s\n",     hfp->fname);
  p7_hmmfile_Close(hfp);

  bg = p7_bg_Create(abc);                
  gm = p7_profile_Create(hmm->M, abc);   

  p7_profile_ConfigUniglocal(gm, hmm, bg, L); 
  sq = esl_sq_CreateDigital(abc);

  /* Emit a sequence, and keep its trace as the "reference" */
  p7_ProfileEmit(rng, hmm, gm, bg, sq, reftr);

  /* Viterbi alignment of that sequence back to the model: Viterbi trace is the "test" */
  vmx = p7_refmx_Create(gm->M, sq->n);

  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, sq->n);
  p7_bg_SetLength(bg, sq->n);
  p7_ReferenceViterbi(sq->dsq, sq->n, gm, vmx, testtr, NULL);

  p7_trace_DumpAnnotated(stdout, reftr,  gm, sq->dsq);
  p7_trace_DumpAnnotated(stdout, testtr, gm, sq->dsq);

  /* Calculate trace accuracy metrics. */
  p7_trace_metrics(reftr, testtr, gm, sq->dsq);

  p7_refmx_Destroy(vmx);
  p7_trace_Destroy(reftr);
  p7_trace_Destroy(testtr);
  esl_sq_Destroy(sq);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
