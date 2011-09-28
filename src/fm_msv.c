/* MSV algorithm; generic (non-SIMD) version.
 * 
 * Contents:
 *   1. MSV implementation.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 * 
 * SRE, Fri Aug 15 10:38:21 2008 [Janelia]
 * SVN $Id: generic_msv.c 3582 2011-06-26 20:09:21Z wheelert $
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"


#include "hmmer.h"

/*****************************************************************
 * 1. MSV implementation.
 *****************************************************************/


int
FM_Recurse( int depth, int M,
            const FM_DATA *fmf, const FM_DATA *fmb,
            const FM_CFG *fm_cfg, const FM_HMMDATA *fm_hmmdata,
            int first, int last, float sc_threshFM,
            FM_DP_PAIR *dp_pairs, FM_INTERVAL *interval,
            int **starts, int **ends, int *hit_cnt,
            long *node_cnts, long *diag_cnts,
            char *seq
          )
{


  //bounding cutoffs
  int max_depth = 16;
  int ssv_req = 16;
  int neg_len_limit = 4;
  int consec_pos_req = 5; //6
  float score_ratio_req = 0.45; //.49

  float sc, next_score;

  int max_k;
  int c, i, k;

  FM_INTERVAL interval_new;
  interval_new.lower = -1;

  for (c=0; c< fm_cfg->meta->alph_size; c++) {//acgt
    int dppos = last;

    seq[depth-1] = fm_cfg->meta->alph[c];
    seq[depth] = '\0';

    float max_sc = 0.0;

    for (i=first; i<=last; i++) { // for each surviving diagonal from the previous round

        k = dp_pairs[i].pos + 1;  // forward
        next_score = fm_hmmdata->scores[k][c];
        sc = dp_pairs[i].score + next_score;
        if (sc > max_sc) {
          max_sc = sc;
          max_k = k;
        }

        /*
        if (sc >= sc_threshFM) { // hit threshold, so mark it as a seed and stop extending
          diags[diag_end].k = dp_pairs[i].pos - depth + 1;
          diags[diag_end].score = dp_pairs[i].score;
          diags[diag_end].interval = interval;
          diag_end++;

          (*hit_cnt)++;
          continue;
            //printf("hit %4d, depth %2d\n", *hit_cnt, depth);
        }
        */

//        if ( ) // forward condition
//          continue; // not possible to hit threshold

        if (  (sc < sc_threshFM && ( depth == ssv_req || k == M ))                                        //didn't hit threshold
           || (depth > ssv_req - 10 &&  sc + fm_hmmdata->opt_ext_fwd[k][ssv_req-depth-1] < sc_threshFM )  //can't hit threshold, even with best possible extension up to length ssv_req
           || (depth == dp_pairs[i].max_score_len + neg_len_limit)                                        //too many consecutive positions with a negative total score contribution (sort of like Xdrop)
           || (sc/depth < score_ratio_req)                                                                //score density is too low
           || (dp_pairs[i].max_consec_pos < consec_pos_req  &&                                            //a seed is expected to have at least one consec_pos_req-long run of positive-scoring matches;  if it hasn't,  (see Tue Nov 23 09:39:54 EST 2010)
                  ( (depth >= max_depth/2 &&  sc/depth < sc_threshFM/max_depth)                              // if we're close to the end of the sequence, abort -- if that end does have sufficiently long all-positive run, I'll find it on the reverse sweep
                  || depth == max_depth-consec_pos_req+1 )                                                   // if we're at least half way across the sequence, and score density is too low, abort -- if the density on the other side is high enough, I'll find it on the reverse sweep
              )
         )
        {
        //do nothing - it's been pruned

        } else if ( sc>0 && ( k < M || sc >= sc_threshFM) ) { // if either (a) score is above threshold, or it's non-negative and
                               //there are more positions in the model after this one, add it to the list of extendable diagonals
            dppos++;
            dp_pairs[dppos].pos = k;
            dp_pairs[dppos].score = sc;

            if (sc > dp_pairs[i].max_score) {
              dp_pairs[dppos].max_score = sc;
              dp_pairs[dppos].max_score_len = depth;
            } else {
              dp_pairs[dppos].max_score = dp_pairs[i].max_score;
              dp_pairs[dppos].max_score_len = dp_pairs[i].max_score_len;
            }

            dp_pairs[dppos].consec_pos =  (next_score > 0 ? dp_pairs[i].consec_pos + 1 : 0);
            dp_pairs[dppos].max_consec_pos = ESL_MAX( dp_pairs[dppos].consec_pos, dp_pairs[i].max_consec_pos);


        }
    }

/*
    if (dppos > last ){  // at least one useful extension


      node_cnts[depth-1]++;
      diag_cnts[depth-1]+=dppos-last;

      if ( depth == max_depth ) {
//      if ( max_sc  >= sc_thresh50) {
        // this is a bit aggressive - if there are multiple existing diagonals, I should
        // probably only approve the ones above thresh, and keep extending others.
        (*hit_cnt)++;
      } else {


        //node_cnts[18]++; // counting the number of interval computations

        int count;
        if (interval->lower > 0 && interval->lower <= interval->upper) {
//          count = bwt_getOccCount (meta, fmindex->occCnts_sb, fmindex->occCnts_b, fmindex->BWT, interval->lower-1, c);
          count = bwt_getOccCount_SSE (meta, fmindex->occCnts_sb, fmindex->occCnts_b, fmindex->BWT, interval->lower-1, c);
          interval_new.lower = fmindex->C[c] + count;


//          count = bwt_getOccCount (meta, fmindex->occCnts_sb, fmindex->occCnts_b, fmindex->BWT, interval->upper, c);
          count = bwt_getOccCount_SSE (meta, fmindex->occCnts_sb, fmindex->occCnts_b, fmindex->BWT, interval->upper, c);
          interval_new.upper = fmindex->C[c] + count - 1;
        }

        //printf ("%-18s : %d, %d, (%.2f, %d)\n", seq, interval_new.lower, interval_new.upper, max_sc, dppos-last );

        if ( interval_new.lower < 0 || interval_new.lower > interval_new.upper  )  //that suffix doesn't exist
          continue;


        p7_BWT_Recurse ( depth+1, M, fmf, fmb, fm_cfg, fm_hmmdata,
                    last+1, dppos, sc_threshFM, dp_pairs, &interval,
                   starts, ends, hit_cnt, node_cnts, diag_cnts
                   , seq
              );

      }

    } else {
      //printf ("internal stop: %d\n", depth);
    }
*/
  }

  return eslOK;
}


int FM_getSeeds (const P7_OPROFILE *gm, P7_GMX *gx, float sc_threshFM,
                int **starts, int** ends, int *hit_cnt,
                const FM_DATA *fmf, const FM_DATA *fmb,
                const FM_CFG *fm_cfg, const FM_HMMDATA *fm_hmmdata )
{
  FM_INTERVAL interval;
  //ESL_DSQ c;
  int i, k;
  int status;
  float sc;
  char         *seq;

  FM_DP_PAIR dp_pairs[10000]; // should always be more than enough

  long node_cnts[20];
  long diag_cnts[20];
  for (i=0; i<20; i++) {
    node_cnts[i] = diag_cnts[i] = 0;
  }



  ESL_ALLOC(seq, 50*sizeof(char));


  for (i=0; i<fm_cfg->meta->alph_size; i++) {//skip '$'
    int cnt=0;

    interval.lower = fmf->C[i];
    interval.upper  = abs(fmf->C[i+1])-1;

    if (interval.lower<0 ) //none of that character found
      continue;


    //c = i - (i <= gm->abc->K ? 1 : 0);  // shift to the easel alphabet
    seq[0] = fm_cfg->meta->alph[i];//gm->abc->sym[i];
    seq[1] = '\0';

    // fill in a DP column for the character c, (compressed so that only positive-scoring entries are kept
    for (k = 1; k < gm->M; k++) // there's no need to bother keeping an entry starting at the last position (gm->M)
    {
      sc = fm_hmmdata->scores[k][i];
      if (sc>0) { // we'll extend any positive-scoring diagonal
        dp_pairs[cnt].pos = k;
        dp_pairs[cnt].score = sc;
        dp_pairs[cnt].max_score = sc;
        dp_pairs[cnt].max_score_len = 1;
        dp_pairs[cnt].consec_pos = 1;
        dp_pairs[cnt].max_consec_pos = 1;
        cnt++;
      }

//      printf ("%2d : %2d : %.2f (scores[%d][%d])\n", 1, k, sc, k, c );

    }


    printf ("char: %d: cnt:%d\n", i, cnt);
    node_cnts[0]++;
    diag_cnts[0]+=cnt;


    FM_Recurse ( 2, gm->M, fmf, fmb, fm_cfg, fm_hmmdata,
                 0, cnt-1, sc_threshFM, dp_pairs, &interval,
                 starts, ends, hit_cnt, node_cnts, diag_cnts
                 ,seq
            );
  }

ERROR:
  ESL_EXCEPTION(eslEMEM, "Error allocating memory for hit list\n");

}

/* Function:  p7_FM_MSV()
 * Synopsis:  Finds windows with MSV scores above given threshold
 *
 * Details:   Uses FM-index to find high-scoring diagonals (seeds), then extends those
 *            seeds to maximal scoring diagonals (no gaps). Diagonals in close
 *            proximity (within a small window) are stitched together, and windows
 *            meeting the MSV scoring threshold (usually score s.t. p=0.02) are
 *            captured, and passed on to the Viterbi and Forward stages of the
 *            pipeline.
 *
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues
 *            gm      - profile (can be in any mode)
 *            gx      - DP matrix
 *            nu      - configuration: expected number of hits (use 2.0 as a default)
 *            bg      - the background model, required for translating a P-value threshold into a score threshold
 *            P       - p-value below which a region is captured as being above threshold
 *            starts  - RETURN: array of start positions for windows surrounding above-threshold areas
 *            ends    - RETURN: array of end positions for windows surrounding above-threshold areas
 *            hit_cnt - RETURN: count of entries in the above two arrays
 *
 *
 * Note:      Not worried about speed here. Based on p7_GMSV
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_FM_MSV( P7_OPROFILE *om, P7_GMX *gx, float nu, P7_BG *bg, double P, int **starts, int** ends, int *hit_cnt,
         const FM_DATA *fmf, const FM_DATA *fmb, const FM_CFG *fm_cfg, const FM_HMMDATA *fm_hmmdata)
{

  float P_fm = 0.5;
  float nullsc, sc_thresh, sc_threshFM;
  float invP, invP_FM;

  int i;
//  float      **dp    = gx->dp;
//  float       *xmx   = gx->xmx;
  float      tloop = logf((float) om->max_length / (float) (om->max_length+3));
  float      tloop_total = tloop * om->max_length;
  float      tmove = logf(     3.0f / (float) (om->max_length+3));
  float      tbmk  = logf(     2.0f / ((float) om->M * (float) (om->M+1)));
//  float    tej   = logf((nu - 1.0f) / nu);
  float      tec   = logf(1.0f / nu);
//  int          i,k;
  int          status;
  int hit_arr_size = 50; // arbitrary size - it, and the hit lists it relates to, will be increased as necessary


  /*
   * Computing the score required to let P meet the F1 prob threshold
   * In original code, converting from an MSV score S (the score getting
   * to state C) to a probability goes like this:
   *  S = XMX(L,p7G_C)
   *  usc = S + tmove + tloop_total
   *  P = f ( (usc - nullsc) / eslCONST_LOG2 , mu, lambda)
   * and we're computing the threshold score S, so reverse it:
   *  (usc - nullsc) /  eslCONST_LOG2 = inv_f( P, mu, lambda)
   *  usc = nullsc + eslCONST_LOG2 * inv_f( P, mu, lambda)
   *  S = usc - tmove - tloop_total
   *
   *  Here, I compute threshold with length model based on max_length.  Usually, the
   *  length of a window returned by this scan will be 2*max_length-1 or longer.  Doesn't
   *  really matter - in any case, both the bg and om models will change with roughly
   *  1 bit for each doubling of the length model, so they offset.
   */
  p7_bg_SetLength(bg, om->max_length);
  p7_oprofile_ReconfigMSVLength(om, om->max_length);
  p7_bg_NullOne  (bg, NULL, om->max_length, &nullsc);

  invP = esl_gumbel_invsurv(P, om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  sc_thresh =   nullsc  + (invP * eslCONST_LOG2) - tmove - tloop_total;

  invP_FM = esl_gumbel_invsurv(P_fm, om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  sc_threshFM =   nullsc  + (invP_FM * eslCONST_LOG2) - tmove - tloop_total - tmove - tbmk - tec;


  /* stuff used to keep track of hits (regions passing the p-threshold)*/
  ESL_ALLOC(*starts, hit_arr_size * sizeof(int));
  ESL_ALLOC(*ends, hit_arr_size * sizeof(int));
  (*starts)[0] = (*ends)[0] = -1;
  *hit_cnt = 0;


  FM_getSeeds(om, gx, sc_threshFM, starts, ends, hit_cnt, fmf, fmb, fm_cfg, fm_hmmdata);



  int L;  //shouldn't be here
  /*
   * merge overlapping windows, compressing list in place.
   */
  if ( *hit_cnt > 0 ) {
    int merged = 0;
    do {
      merged = 0;

      int new_hit_cnt = 0;
      for (i=1; i<*hit_cnt; i++) {
        if ((*starts)[i] <= (*ends)[new_hit_cnt]) {
          //merge windows
          merged = 1;
          if ( (*ends)[i] > (*ends)[new_hit_cnt] )
            (*ends)[new_hit_cnt] = (*ends)[i];

        } else {
          //new window
          new_hit_cnt++;
          (*starts)[new_hit_cnt] = (*starts)[i];
          (*ends)[new_hit_cnt] = (*ends)[i];
        }
      }
      *hit_cnt = new_hit_cnt + 1;
    } while (merged > 0);

    if ((*starts)[0] <  1)
      (*starts)[0] =  1;

    if ((*ends)[*hit_cnt - 1] >  L)
      (*ends)[*hit_cnt - 1] =  L;

  }
  return eslOK;
ERROR:
  ESL_EXCEPTION(eslEMEM, "Error allocating memory for hit list\n");

}
/*------------------ end, p7_FM_MSV() ------------------------*/





/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7GENERIC_MSV_BENCHMARK
/*
   gcc -g -O2      -o generic_msv_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_MSV_BENCHMARK generic_msv.c -lhmmer -leasel -lm
   icc -O3 -static -o generic_msv_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_MSV_BENCHMARK generic_msv.c -lhmmer -leasel -lm
   ./benchmark-generic-msv <hmmfile>
 */
/* As of Fri Dec 28 14:48:39 2007
 *    Viterbi  = 61.8 Mc/s
 *    Forward  =  8.6 Mc/s
 *   Backward  =  7.1 Mc/s
 *       GMSV  = 55.9 Mc/s
 * (gcc -g -O2, 3.2GHz Xeon, N=50K, L=400, M=72 RRM_1 model)
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,  "20000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic MSV";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *gx      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL);
  gx = p7_gmx_Create(gm->M, L);

  /* Baseline time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_GMSV      (dsq, L, gm, gx, 2.0, &sc);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_MSV_BENCHMARK*/
/*----------------- end, benchmark ------------------------------*/

/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7GENERIC_MSV_TESTDRIVE
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"
/* The MSV score can be validated against Viterbi (provided we trust
 * Viterbi), by creating a multihit local profile in which:
 *   1. All t_MM scores = 0
 *   2. All other core transitions = -inf
 *   3. All t_BMk entries uniformly log 2/(M(M+1))
 */
static void
utest_msv(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, P7_PROFILE *gm, int nseq, int L)
{
  P7_PROFILE *g2 = NULL;
  ESL_DSQ   *dsq = NULL;
  P7_GMX    *gx  = NULL;
  float     sc1, sc2;
  int       k, idx;

  if ((dsq    = malloc(sizeof(ESL_DSQ) *(L+2))) == NULL)  esl_fatal("malloc failed");
  if ((gx     = p7_gmx_Create(gm->M, L))        == NULL)  esl_fatal("matrix creation failed");
  if ((g2     = p7_profile_Clone(gm))           == NULL)  esl_fatal("profile clone failed");

  /* Make g2's scores appropriate for simulating the MSV algorithm in Viterbi */
  esl_vec_FSet(g2->tsc, p7P_NTRANS * g2->M, -eslINFINITY);
  for (k = 1; k <  g2->M; k++) p7P_TSC(g2, k, p7P_MM) = 0.0f;
  for (k = 0; k <  g2->M; k++) p7P_TSC(g2, k, p7P_BM) = log(2.0f / ((float) g2->M * (float) (g2->M+1)));

  for (idx = 0; idx < nseq; idx++)
    {
      if (esl_rsq_xfIID(r, bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");

      if (p7_GMSV    (dsq, L, gm, gx, 2.0, &sc1)       != eslOK) esl_fatal("MSV failed");
      if (p7_GViterbi(dsq, L, g2, gx,      &sc2)       != eslOK) esl_fatal("viterbi failed");
      if (fabs(sc1-sc2) > 0.0001) esl_fatal("MSV score not equal to Viterbi score");
    }

  p7_gmx_Destroy(gx);
  p7_profile_Destroy(g2);
  free(dsq);
  return;
}
#endif /*p7GENERIC_MSV_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/


/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
/* gcc -g -Wall -Dp7GENERIC_MSV_TESTDRIVE -I. -I../easel -L. -L../easel -o generic_msv_utest generic_msv.c -lhmmer -leasel -lm
 */
#ifdef p7GENERIC_MSV_TESTDRIVE
#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"

#include "p7_config.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "--vv",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be very verbose",                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for the generic Msv implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_BG          *bg   = NULL;
  int             M    = 100;
  int             L    = 200;
  int             nseq = 20;
  char            errbuf[eslERRBUFSIZE];

  if ((abc = esl_alphabet_Create(eslAMINO))         == NULL)  esl_fatal("failed to create alphabet");
  if (p7_hmm_Sample(r, M, abc, &hmm)                != eslOK) esl_fatal("failed to sample an HMM");
  if ((bg = p7_bg_Create(abc))                      == NULL)  esl_fatal("failed to create null model");
  if ((gm = p7_profile_Create(hmm->M, abc))         == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL)    != eslOK) esl_fatal("failed to config profile");
  if (p7_hmm_Validate    (hmm, errbuf, 0.0001)      != eslOK) esl_fatal("whoops, HMM is bad!: %s", errbuf);
  if (p7_profile_Validate(gm,  errbuf, 0.0001)      != eslOK) esl_fatal("whoops, profile is bad!: %s", errbuf);

  utest_msv(go, r, abc, bg, gm, nseq, L);

  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_MSV_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/


/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7GENERIC_MSV_EXAMPLE
/* 
   gcc -g -O2 -Dp7GENERIC_MSV_EXAMPLE -I. -I../easel -L. -L../easel -o generic_msv_example generic_msv.c -lhmmer -leasel -lm
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "--nu",      eslARG_REAL,   "2.0", NULL, NULL,  NULL,  NULL, NULL, "set nu param to <x>: expected # MSV diagonals",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of generic MSV algorithm";


int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  float           nu      = esl_opt_GetReal(go, "--nu");
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *fwd     = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  P7_TRACE       *tr      = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           sc, nullsc, seqscore, lnP;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);

  /* Allocate matrix */
  fwd = p7_gmx_Create(gm->M, sq->n);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_ReconfigLength(gm,  sq->n);
      p7_bg_SetLength(bg,    sq->n);
      p7_gmx_GrowTo(fwd, gm->M, sq->n); 

      /* Run MSV */
      p7_GMSV(sq->dsq, sq->n, gm, fwd, nu, &sc);

      /* Calculate bit score and P-value using standard null1 model*/
      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);
      seqscore = (sc - nullsc) / eslCONST_LOG2;
      lnP      =  esl_gumbel_logsurv(seqscore,  gm->evparam[p7_MMU],  gm->evparam[p7_MLAMBDA]);

      /* output suitable for direct use in profmark benchmark postprocessors:
       * <Pvalue> <bitscore> <target name> <query name>
       */
      printf("%g\t%.2f\t%s\t%s\n", exp(lnP), seqscore, sq->name, hmm->name);

      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n", sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  /* Cleanup */
  esl_sqfile_Close(sqfp); 
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_gmx_Destroy(fwd);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_MSV_EXAMPLE*/
/*-------------------- end, example -----------------------------*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
