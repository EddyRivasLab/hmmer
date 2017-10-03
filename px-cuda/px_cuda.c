#include <string.h>
#include "easel.h"
#include "esl_dsqdata.h"
#include <x86intrin.h>
#include <math.h>
#include "esl_sse.h"
#include "hmmer.h"
#include "px_cuda.h"
char * restripe_char(char *source, int source_chars_per_vector, int dest_chars_per_vector, int source_length, int *dest_length){
  char *dest;
  int dest_num_vectors;
  int source_num_vectors;
  source_num_vectors = source_length/source_chars_per_vector;
  if(source_num_vectors * source_chars_per_vector != source_length){
    source_num_vectors++;
  }
  dest_num_vectors = source_length/dest_chars_per_vector;
  if(dest_num_vectors * dest_chars_per_vector != source_length){
    dest_num_vectors++;  //round up if source length isn't a multiple of the dest vector length
  }

  dest = (char *) malloc(dest_num_vectors * dest_chars_per_vector);
  *dest_length = dest_num_vectors * dest_chars_per_vector;

  int source_pos, dest_pos;
  int i;

  for(i = 0; i < source_length; i++){
    source_pos = ((i % source_num_vectors) * source_chars_per_vector) + (i / source_num_vectors);
    dest_pos = ((i % dest_num_vectors) * dest_chars_per_vector) + (i / dest_num_vectors);
    dest[dest_pos] = source[source_pos];
  }

  // pad out the dest vector with zeroes if necessary
  for(i = source_length; i < *dest_length; i++){
      dest_pos = ((i % dest_num_vectors) * dest_chars_per_vector) + (i / dest_num_vectors);
    dest[dest_pos] = 0;
  }

  return dest;

}

int
p7_SSVFilter_shell_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *fx, float *ret_sc)
{
  int      Q          = P7_Q(om->M, p7_VWIDTH_SSE);
  __m128i  hv         = _mm_set1_epi8(-128);
  __m128i  neginfmask = _mm_insert_epi8( _mm_setzero_si128(), -128, 0);
  __m128i *dp;
  __m128i *rbv;
  __m128i  mpv;
  __m128i  sv;
  int8_t   h;
  int      i,q;
  int      status;

  char *card_rbv= NULL;
  char *card_temp_buffer = NULL;
  char *card_hv, *card_mpv;
  char *restriped_rbv;
  int card_rbv_length;
  int card_rbv_allocated = 0;
  int Q_card = P7_Q(om->M, 128);

  if (( status = p7_filtermx_Reinit(fx, om->M) ) != eslOK) goto FAILURE;
  fx->M    = om->M;
  fx->Vw   = p7_VWIDTH_SSE / sizeof(int16_t); // A hack. FILTERMX wants Vw in units of int16_t. 
  fx->type = p7F_SSVFILTER;
  dp       = (__m128i *) fx->dp;

  card_rbv_allocated = 128 * P7_Q(om->M, 128);  // compute size of buffers required for a model of width M

  cudaMalloc((void **) &card_temp_buffer, card_rbv_length);// Don't do this in production code.
  cudaMalloc((void **) &card_rbv, card_rbv_length);// Don't do this in production code.
  cudaMalloc((void **) &card_hv, 128);
  cudaMalloc((void **) &card_mpv, 128);
  <<<1, 32>>>SSV_start_cuda((unsigned int *) card_temp_buffer, (unsigned int *) card_hv, (unsigned int *) card_mpv, Q_card);
  mpv = hv;
  for (q = 0; q < Q; q++)
    dp[q] = hv;

  for (i = 1; i <= L; i++)
    {
      rbv = (__m128i *) om->rbv[dsq[i]];
      restriped_rbv = restripe_char((char *)rbv, om->V, 128, Q * p7_VWIDTH_SSE, &card_rbv_length);
      if (card_rbv_allocated < card_rbv_length){ // we need to allocate more space for the rbv vector on the card
        printf("Didn't allocate enough buffer space -- something went wrong\n");
      }

      cudaMemcpy(card_rbv, restriped_rbv, card_rbv_length, cudaMemcpyHostToDevice);
      SSV_row_cuda((unsigned int *) card_rbv, (unsigned int *) card_temp_buffer, P7_Q(om->M, 128));
      for (q = 0; q < Q; q++)
        {
          sv    = _mm_adds_epi8(mpv, rbv[q]);
          hv    = _mm_max_epi8(hv, sv);
          mpv   = dp[q];
          dp[q] = sv;
        }
      mpv = esl_sse_rightshift_int8(sv, neginfmask);
    }
  h = esl_sse_hmax_epi8(hv);

  if (h == 127)  
    { *ret_sc = eslINFINITY; return eslERANGE; }
  else if (h > -128)
    { 
      *ret_sc = ((float) h + 128.) / om->scale_b + om->tauBM - 2.0;   // 2.0 is the tauNN/tauCC "2 nat approximation"
      *ret_sc += 2.0 * logf(2.0 / (float) (L + 2));                   
      return eslOK;
    }
  else 
    {
      *ret_sc = -eslINFINITY;
      return eslOK;
    }
    
 FAILURE:
  *ret_sc = -eslINFINITY;
  return status;
}


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
    float ssv_score;

    p7_SSVFilter_shell_sse(chu->dsq[i], chu->L[i], om, eng->fx ,&ssv_score);
  	
	  
	  
	  p7_engine_Reuse(eng);
    num_hits++;
	}
      esl_dsqdata_Recycle(dd, chu);
    }
    printf("Saw %ld sequences\n", num_hits);
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





