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
#include <string.h>

/*****************************************************************
 * 1. MSV implementation.
 *****************************************************************/




/* hit_sorter(): qsort's pawn, below */
static int
hit_sorter(const void *a, const void *b)
{
    FM_DIAG *d1 = (FM_DIAG*)a;
    FM_DIAG *d2 = (FM_DIAG*)b;

    return 2 * (d1->sortkey > d2->sortkey) - 1;  // same as the test below
//    if      ( d1->sortkey > d2->sortkey) return 1;
//    else                                 return -1;

}

int
mergeSeeds(FM_DIAGLIST *seeds, int N, int msv_length) {
  int i;
  int j = 0;

  if (seeds->count == 0)
    return eslOK;

  FM_DIAG *diags = seeds->diags;
  FM_DIAG next = diags[0];

  int tmp;
  int next_is_complement;

  int curr_is_complement = next.sortkey > N;
  int compl_mult         = curr_is_complement ? 1 : -1;
  int curr_n             = next.n;
  int curr_k             = next.k;
  int curr_len           = next.length;
  int curr_end           = curr_n + curr_len - 1;
  int curr_diagval       = next.n + next.k * compl_mult;


  for( i=1; i<seeds->count; i++) {

    next = diags[i];
    next_is_complement = next.sortkey > N;
    compl_mult         = next_is_complement ? 1 : -1;

    //overlapping diagonals will share the same value of (n - k) for non-complement hits, and (n+k) for complement hits
    if (  next_is_complement == curr_is_complement              //same direction
          && ( next.n + next.k * compl_mult) == curr_diagval    //same diagonal
          && next.n + next.length < curr_n + curr_len + msv_length       //overlapping, or close to it
       ) {


      //overlapping diags; extend, if appropriate
      if (next_is_complement) {
        curr_len += curr_n - next.n;
      } else {
        tmp = next.n + next.length - 1;
        if (tmp > curr_end) {
          curr_end = tmp;
          curr_len = curr_end - curr_n + 1;
        }
      }
    } else {
      //doesn't overlap current diagonal, so store current...
      diags[j].n      = curr_n;
      diags[j].k      = curr_k;
      diags[j].length = curr_end - curr_n + 1;
      diags[j].complementarity = curr_is_complement;
      diags[j].sortkey = 0.0;

      // ... then start up a new one
      curr_n   = next.n;
      curr_k   = next.k;
      curr_len = next.length;
      curr_end = curr_n + curr_len - 1;
      curr_diagval = next.n + next.k * compl_mult;
      curr_is_complement = next_is_complement;

      j++;
    }
  }
  // store final entry
  diags[j].n      = curr_n;
  diags[j].k      = curr_k;
  diags[j].length = curr_end - curr_n + 1;
  diags[j].sortkey = 0.0;

  seeds->count = j+1;
  return eslOK;

}

uint32_t
FM_backtrackSeed(const FM_DATA *fmf, FM_CFG *fm_cfg, int i, FM_DIAG *seed) {
  int j = i;
  int len = 0;
  int c;

  while ( j != fmf->term_loc && (j & fm_cfg->maskSA)) { //go until we hit a position in the full SA that was sampled during FM index construction
    c = fm_getChar( fm_cfg->meta->alph_type, j, fmf->BWT);
    j = fm_getOccCount (fmf, fm_cfg, j-1, c);
    j += abs(fmf->C[c]);
    len++;
  }

  return len + (j==fmf->term_loc ? 0 : fmf->SA[ j >> fm_cfg->shiftSA ]) ; // len is how many backward steps we had to take to find a sampled SA position

}

int
FM_getPassingDiags(const FM_DATA *fmf, FM_CFG *fm_cfg,
            int k, int M, float sc, int depth, int fm_direction,
            int model_direction, int complementarity,
            FM_INTERVAL *interval,
            FM_DIAGLIST *seeds
            )
{

  int i;
  FM_DIAG *seed;

  //iterate over the forward interval, for each entry backtrack until hitting a sampled suffix array entry
  for (i = interval->lower;  i<= interval->upper; i++) {

    seed = fm_newSeed(seeds);
    seed->k      = k;
    seed->length = depth;
    seed->n      = FM_backtrackSeed(fmf, fm_cfg, i, seed);
    seed->complementarity = complementarity;

    /*  In both fm_forward and fm_backward cases, seed->n corresponds to the start
     * of the seed in terms of the target sequence, in a forward direction. But the meaning of
     * seed->k differs: in fm_backward, k represents the model position aligned to n, while
     * in fm_forward, k holds the model position aligned to the last character of the seed.
     * So ... shift k over in the forward case.  Doing this will mean
     */

    if(fm_direction == fm_forward) //
      seed->k += (depth - 1) * (model_direction == fm_forward ? -1 : 1) ;

    seed->sortkey =   ( complementarity == fm_complement ? fmf->N + 1 : 0)   // makes complement seeds cover a different score range than non-complements
                    + (seed->n + ( seed->k * (complementarity == fm_complement ? 1 : -1)))                                          // unique diagonal within the complement/non-complement score range
                    + ((double)seed->k/(double)(M+1))                                     // fractional part, used to sort seeds sharing a diagonal
                    ;
  }

  return eslOK;
}


int
FM_Recurse( int depth, int M, int fm_direction,
            const FM_DATA *fmf, const FM_DATA *fmb,
            FM_CFG *fm_cfg, const FM_HMMDATA *fm_hmmdata,
            int first, int last, float sc_threshFM,
            FM_DP_PAIR *dp_pairs,
            FM_INTERVAL *interval_1, FM_INTERVAL *interval_2,
            FM_DIAGLIST *seeds,
            char *seq
          )
{


  float sc, next_score;

  //int max_k;
  int c, i, k;

  FM_INTERVAL interval_1_new, interval_2_new;

  for (c=0; c< fm_cfg->meta->alph_size; c++) {//acgt
    int dppos = last;

    seq[depth-1] = fm_cfg->meta->alph[c];
    seq[depth] = '\0';

    for (i=first; i<=last; i++) { // for each surviving diagonal from the previous round

        if (dp_pairs[i].model_direction == fm_forward)
          k = dp_pairs[i].pos + 1;
        else  //fm_backward
          k = dp_pairs[i].pos - 1;

        if (dp_pairs[i].complementarity == fm_complement) {
          next_score = fm_hmmdata->s.scores_f[k][ fm_getComplement(c,fm_cfg->meta->alph_type) ];
        } else
          next_score = fm_hmmdata->s.scores_f[k][c];

        sc = dp_pairs[i].score + next_score;


        if ( sc >= sc_threshFM ) { // this is a seed I want to extend

          interval_1_new.lower = interval_1->lower;
          interval_1_new.upper = interval_1->upper;

          //fprintf(stderr, "%-18s : (%.2f)\n", seq, sc);

          if (fm_direction == fm_forward) {
            //searching for forward matches on the FM-index
            interval_2_new.lower = interval_2->lower;
            interval_2_new.upper = interval_2->upper;
            if ( interval_1_new.lower >= 0 && interval_1_new.lower <= interval_1_new.upper  )  //no use extending a non-existent string
              fm_updateIntervalForward( fmb, fm_cfg, c, &interval_1_new, &interval_2_new);

            if ( interval_1_new.lower >= 0 && interval_1_new.lower <= interval_1_new.upper  )  //no use passing a non-existent string
              FM_getPassingDiags(fmf, fm_cfg, k, M, sc, depth, fm_forward,
                                 dp_pairs[i].model_direction, dp_pairs[i].complementarity,
                                 &interval_2_new, seeds);

          } else {
            //searching for reverse matches on the FM-index
            if ( interval_1_new.lower >= 0 && interval_1_new.lower <= interval_1_new.upper  )  //no use extending a non-existent string
              fm_updateIntervalReverse( fmf, fm_cfg, c, &interval_1_new);

            if ( interval_1_new.lower >= 0 && interval_1_new.lower <= interval_1_new.upper  )  //no use passing a non-existent string
              FM_getPassingDiags(fmf, fm_cfg, k, M, sc, depth, fm_backward,
                                 dp_pairs[i].model_direction, dp_pairs[i].complementarity,
                                 &interval_1_new, seeds);

          }

        } else if (  sc <= 0                                                                                        //some other path in the string enumeration tree will do the job
            || depth == fm_cfg->max_depth                                                                            //can't extend anymore, 'cause we've reached the pruning length
            || (depth == dp_pairs[i].max_score_len + fm_cfg->neg_len_limit)                                        //too many consecutive positions with a negative total score contribution (sort of like Xdrop)
            || (sc/depth < fm_cfg->score_ratio_req)                                                                //score density is too low
            || (dp_pairs[i].max_consec_pos < fm_cfg->consec_pos_req  &&                                            //a seed is expected to have at least one run of positive-scoring matches at least length consec_pos_req;  if it hasn't,  (see Tue Nov 23 09:39:54 EST 2010)
                   ( (depth >= fm_cfg->max_depth/2 &&  sc/depth < sc_threshFM/fm_cfg->max_depth)                              // if we're close to the end of the sequence, abort -- if that end does have sufficiently long all-positive run, I'll find it on the reverse sweep
                   || depth == fm_cfg->max_depth-fm_cfg->consec_pos_req+1 )                                                   // if we're at least half way across the sequence, and score density is too low, abort -- if the density on the other side is high enough, I'll find it on the reverse sweep
               )
            || (dp_pairs[i].model_direction == fm_forward  &&
                   ( (k == M)                                                                                                          //can't extend anymore, 'cause we're at the end of the model, going forward
                  || (depth > (fm_cfg->max_depth - 10) &&  sc + fm_hmmdata->opt_ext_fwd[k][fm_cfg->max_depth-depth-1] < sc_threshFM)   //can't hit threshold, even with best possible forward extension up to length ssv_req
                  ))
            || (dp_pairs[i].model_direction == fm_backward &&
                   ( (k == 1)                                                                                                          //can't extend anymore, 'cause we're at the beginning of the model, going backwards
                  || (depth > (fm_cfg->max_depth - 10) &&  sc + fm_hmmdata->opt_ext_rev[k][fm_cfg->max_depth-depth-1] < sc_threshFM )  //can't hit threshold, even with best possible extension up to length ssv_req
                  ))
         )
        {

          //do nothing - it's been pruned

        } else { // it's possible to extend this diagonal and reach the threshold score

            dppos++;
            dp_pairs[dppos].pos = k;
            dp_pairs[dppos].score = sc;
            dp_pairs[dppos].model_direction   = dp_pairs[i].model_direction;
            dp_pairs[dppos].complementarity   = dp_pairs[i].complementarity;

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



    if ( dppos > last ){  // at least one diagonal that might reach threshold score, but hasn't yet, so extend

      interval_1_new.lower = interval_1->lower;
      interval_1_new.upper = interval_1->upper;

      if (fm_direction == fm_forward) {
        interval_2_new.lower = interval_2->lower;
        interval_2_new.upper = interval_2->upper;
        if ( interval_1_new.lower >= 0 && interval_1_new.lower <= interval_1_new.upper  )  //no use extending a non-existent string
          fm_updateIntervalForward( fmb, fm_cfg, c, &interval_1_new, &interval_2_new);

        if (  interval_1_new.lower < 0 || interval_1_new.lower > interval_1_new.upper )  //that string doesn't exist in fwd index
            continue;

        FM_Recurse(depth+1, M, fm_direction,
                  fmf, fmb, fm_cfg, fm_hmmdata,
                  last+1, dppos, sc_threshFM, dp_pairs,
                  &interval_1_new, &interval_2_new,
                  seeds,
                  seq
                  );


      } else {

        if ( interval_1_new.lower >= 0 && interval_1_new.lower <= interval_1_new.upper  )  //no use extending a non-existent string
          fm_updateIntervalReverse( fmf, fm_cfg, c, &interval_1_new);

        if (  interval_1_new.lower < 0 || interval_1_new.lower > interval_1_new.upper )  //that string doesn't exist in reverse index
            continue;

        FM_Recurse(depth+1, M, fm_direction,
                  fmf, fmb, fm_cfg, fm_hmmdata,
                  last+1, dppos, sc_threshFM, dp_pairs,
                  &interval_1_new, NULL,
                  seeds,
                  seq
                  );

      }
    } else {
      //printf ("internal stop: %d\n", depth);
    }

  }

  return eslOK;
}


int FM_getSeeds (const P7_OPROFILE *gm, P7_GMX *gx, float sc_threshFM,
                FM_DIAGLIST *seeds,
                const FM_DATA *fmf, const FM_DATA *fmb,
                FM_CFG *fm_cfg, const FM_HMMDATA *fm_hmmdata )
{
  FM_INTERVAL interval_f1, interval_f2, interval_bk;
  //ESL_DSQ c;
  int i, k;
  int status;
  float sc;
  char         *seq;

  FM_DP_PAIR *dp_pairs_fwd;
  FM_DP_PAIR *dp_pairs_rev;

  ESL_ALLOC(dp_pairs_fwd, gm->M * fm_cfg->max_depth * sizeof(FM_DP_PAIR)); // guaranteed to be enough to hold all diagonals
  ESL_ALLOC(dp_pairs_rev, gm->M * fm_cfg->max_depth * sizeof(FM_DP_PAIR));

  ESL_ALLOC(seq, 50*sizeof(char));


  for (i=0; i<fm_cfg->meta->alph_size; i++) {//skip '$'
    int fwd_cnt=0;
    int rev_cnt=0;

    interval_f1.lower = interval_f2.lower = interval_bk.lower = fmf->C[i];
    interval_f1.upper = interval_f2.upper = interval_bk.upper = abs(fmf->C[i+1])-1;

    if (interval_f1.lower<0 ) //none of that character found
      continue;


    //c = i - (i <= gm->abc->K ? 1 : 0);  // shift to the easel alphabet
    seq[0] = fm_cfg->meta->alph[i];//gm->abc->sym[i];
    seq[1] = '\0';

    // Fill in a DP column for the character c, (compressed so that only positive-scoring entries are kept)
    // There will be 4 DP columns for each character, (1) fwd-std, (2) fwd-complement, (3) rev-std, (4) rev-complement
    for (k = 1; k <= gm->M; k++) // there's no need to bother keeping an entry starting at the last position (gm->M)
    {


      sc = fm_hmmdata->s.scores_f[k][i];
      if (sc>0) { // we'll extend any positive-scoring diagonal
        if (k < gm->M-2) { // don't bother starting a forward diagonal so close to the end of the model
          //Forward pass on the FM-index
          dp_pairs_fwd[fwd_cnt].pos =             k;
          dp_pairs_fwd[fwd_cnt].score =           sc;
          dp_pairs_fwd[fwd_cnt].max_score =       sc;
          dp_pairs_fwd[fwd_cnt].max_score_len =   1;
          dp_pairs_fwd[fwd_cnt].consec_pos =      1;
          dp_pairs_fwd[fwd_cnt].max_consec_pos =  1;
          dp_pairs_fwd[fwd_cnt].complementarity = fm_nocomplement;
          dp_pairs_fwd[fwd_cnt].model_direction = fm_forward;
          fwd_cnt++;
        }
        if (k > 3) { // don't bother starting a reverse diagonal so close to the start of the model
          dp_pairs_rev[rev_cnt].pos =             k;
          dp_pairs_rev[rev_cnt].score =           sc;
          dp_pairs_rev[rev_cnt].max_score =       sc;
          dp_pairs_rev[rev_cnt].max_score_len =   1;
          dp_pairs_rev[rev_cnt].consec_pos =      1;
          dp_pairs_rev[rev_cnt].max_consec_pos =  1;
          dp_pairs_rev[rev_cnt].complementarity = fm_nocomplement;
          dp_pairs_rev[rev_cnt].model_direction = fm_backward;
          rev_cnt++;
        }
      }


      sc = fm_hmmdata->s.scores_f[k][fm_getComplement(i, fm_cfg->meta->alph_type)];
      if (sc>0) { // we'll extend any positive-scoring diagonal

        //forward on the FM, reverse on the model
        if (k > 3) { // don't bother starting a reverse diagonal so close to the start of the model
          dp_pairs_fwd[fwd_cnt].pos =             k;
          dp_pairs_fwd[fwd_cnt].score =           sc;
          dp_pairs_fwd[fwd_cnt].max_score =       sc;
          dp_pairs_fwd[fwd_cnt].max_score_len =   1;
          dp_pairs_fwd[fwd_cnt].consec_pos =      1;
          dp_pairs_fwd[fwd_cnt].max_consec_pos =  1;
          dp_pairs_fwd[fwd_cnt].complementarity = fm_complement;
          dp_pairs_fwd[fwd_cnt].model_direction = fm_backward;
          fwd_cnt++;
        }


        //reverse on the FM, forward on the model
        if (k < gm->M-2) { // don't bother starting a forward diagonal so close to the end of the model
          dp_pairs_rev[rev_cnt].pos =             k;
          dp_pairs_rev[rev_cnt].score =           sc;
          dp_pairs_rev[rev_cnt].max_score =       sc;
          dp_pairs_rev[rev_cnt].max_score_len =   1;
          dp_pairs_rev[rev_cnt].consec_pos =      1;
          dp_pairs_rev[rev_cnt].max_consec_pos =  1;
          dp_pairs_rev[rev_cnt].complementarity = fm_complement;
          dp_pairs_rev[rev_cnt].model_direction = fm_forward;
          rev_cnt++;
        }
      }


    }

    FM_Recurse ( 2, gm->M, fm_forward,
                fmf, fmb, fm_cfg, fm_hmmdata,
                 0, fwd_cnt-1, sc_threshFM, dp_pairs_fwd,
                 &interval_f1, &interval_f2,
                 seeds
                 , seq
            );


    FM_Recurse ( 2, gm->M, fm_backward,
                fmf, fmb, fm_cfg, fm_hmmdata,
                 0, rev_cnt-1, sc_threshFM, dp_pairs_rev,
                 &interval_bk, NULL,
                 seeds
                 , seq
            );

  }


  //now sort, first by direction, then N (position on database sequence), then K (model position)
  qsort(seeds->diags, seeds->count, sizeof(FM_DIAG), hit_sorter);

  //merge duplicates
  mergeSeeds(seeds, fmf->N, fm_cfg->msv_length);

  return eslOK;

ERROR:
  ESL_EXCEPTION(eslEMEM, "Error allocating memory for hit list\n");

}



/* Function:  FM_extendSeed()
 * Synopsis:  Finds windows with MSV scores above given threshold
 *
 * Details:  Starting with a short FM-based seed, tally up scores along
 *           a stretch of the corresponding diagonal close to the seed.
 *           As scores are tallied, create a window each time the score
 *           threshold is met, just like the full-DP-MSV approach does.
 *           Note that this actually wastes an opportunity to use the
 *           full-length diagonal seed to produce more constrained
 *           windows;  this is done for consistency between the two
 *           methods - a change in windowing strategy leads to differences
 *           in window size, thus subtle differences in domain envelopes,
 *           meaning different scores/e-values
 *           See ~/wheelert/notebook/2011/1101_fmindex_benchmarking/00NOTES, (Nov 7 and 8)
 *
 * Args:
 *
 * Note:
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:
 */
int
FM_extendSeed(FM_DIAG *diag, const FM_DATA *fm, const FM_HMMDATA *hmmdata, FM_CFG *cfg,
              float sc_thresh, float j_thresh, ESL_SQ  *tmp_sq, FM_WINDOWLIST *windowlist)
{

  int k,n;
  int model_start, model_end, target_start, target_end;
  int hit_start, hit_end;
  float sc;
  int c;
  uint32_t tmp_id;
  uint32_t tmp_n;
  uint32_t fm_pos;

  if (diag->complementarity == fm_complement) {
    // complementary diagonal is reverse, so K is the last model position of the diag, not the first
    model_start      = ESL_MAX(1, diag->k - diag->length - cfg->msv_length + 1 + 1 );
    model_end        = ESL_MIN(hmmdata->M, diag->k  - diag->length + cfg->msv_length ); //+1 and -1 cancel
    target_start     = diag->n - (model_end - diag->k);
    target_end       = diag->n + (diag->k - model_start);
  } else {
    model_start      = ESL_MAX(1, diag->k + diag->length - cfg->msv_length  ); //-1 and +1 cancel
    model_end        = ESL_MIN(hmmdata->M, diag->k + cfg->msv_length - 1 );
    target_start     = diag->n - (diag->k - model_start);
    target_end       = diag->n + (model_end - diag->k);
  }


  k = model_start;
  n = 1;
  sc = 0;

  tmp_id = fm_computeSequenceOffset( fm, cfg->meta, 0, target_start);
  fm_convertRange2DSQ(cfg->meta, tmp_id, target_start, target_end-target_start+1, fm->T, tmp_sq );

  if (diag->complementarity == fm_complement)
    esl_sq_ReverseComplement(tmp_sq);

  hit_start = n;

  if (sc_thresh <= j_thresh) {
    //SSV.  Just capture diags passing the sc_thresh, as is done in SIMD-MSV

    for (  ; k <= model_end; k++, n++) {

        c = tmp_sq->dsq[n];

        sc  += hmmdata->s.scores_f[k][c];

        if (sc < 0) {
          sc = 0;
          hit_start = n+1;
        }
        hit_end = n;


        if (sc > sc_thresh) {
          //grab the target position at which the threshold was reached
          if (diag->complementarity == fm_complement)
            fm_pos = target_end - hit_end + 1;
          else
            fm_pos = target_start +  hit_end - 1;

          //then convert to positions in the original sequence used to produce the db
          fm_getOriginalPosition (fm, cfg->meta, 0, hit_end-hit_start+1, fm_forward, fm_pos, &(tmp_id), &(tmp_n) );
          if (tmp_id == -1) continue; // crosses over a barrier between sequences in the digital data structure

          if (windowlist == NULL) //this is used by the MSV merging step, to figure out when the necessary score was reached.
            return tmp_n;
          else
            fm_newWindow(windowlist, tmp_id, tmp_n, fm_pos, k, hit_end-hit_start+1,  sc, diag->complementarity );
          sc = 0;
          hit_start = n+1;
        }
    }

  } else { //sc_thresh > j_thresh ... which calls for MSV
    //  Capture a diagonal if score >= j_thresh, but keep appending to that diagonal
    //  if score is increased (until sc_thresh is reached ... then start over)
    // The task of merging sub-sc_thresh diagonals is left to calling function
    int active_diag = FALSE;
    int compl_dir =  (diag->complementarity == fm_complement ? -1 : 1);

    for (  ; k <= model_end; k++, n++) {

         c = tmp_sq->dsq[n];

         sc  += hmmdata->s.scores_f[k][c];

         hit_end = n;
         if (sc < 0) {
           sc = 0;
           active_diag = FALSE;
           hit_start = n+1;
         } else if (sc > j_thresh) {
           if (active_diag) { // just tack on to the active diagonal (if score is improved)

             if (sc > windowlist->windows[windowlist->count-1].score) {
               windowlist->windows[windowlist->count-1].fm_n  += compl_dir;
               windowlist->windows[windowlist->count-1].n += compl_dir;
               windowlist->windows[windowlist->count-1].k ++;
               windowlist->windows[windowlist->count-1].score = sc;
             }

           } else { // start a new diagonal
             //grab the target position at which the threshold was reached
             if (diag->complementarity == fm_complement)
               fm_pos = target_end - hit_end + 1;
             else
               fm_pos = target_start +  hit_end - 1;

             //then convert to positions in the original sequence used to produce the db
             fm_getOriginalPosition (fm, cfg->meta, 0, hit_end-hit_start+1, fm_forward, fm_pos, &(tmp_id), &(tmp_n));
             if (tmp_id == -1) continue; // crosses over a barrier between sequences in the digital data structure

             fm_newWindow(windowlist, tmp_id, tmp_n, fm_pos, k,  hit_end-hit_start+1, sc, diag->complementarity );
             active_diag = TRUE;
           }

           if (sc > sc_thresh) { // done with that diagonal, start a new one
             sc = 0;
             active_diag = FALSE;
             hit_start = n+1;
           }

         }
     }


  }

  return eslOK;
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
p7_FM_MSV( P7_OPROFILE *om, P7_GMX *gx, float nu, P7_BG *bg, double F1,
         const FM_DATA *fmf, const FM_DATA *fmb, FM_CFG *fm_cfg, const FM_HMMDATA *fm_hmmdata,
         FM_WINDOWLIST *windowlist)
{
  float P_fm = 0.5;
  float sc_thresh, sc_threshFM, j_thresh;
  float invP, invP_FM;
  float nullsc;

  int i, j, k;
  int window_start;
  int window_end;
  int window_start_fm;
  int window_end_fm;

  ESL_SQ   *tmp_sq;
  FM_WINDOW *window;
  FM_WINDOW *prev_window;
  FM_WINDOW *prev;
  FM_WINDOW *curr;
  int active_diagonal;
  int target_len;

  FM_DIAGLIST seeds;

  float      tloop = logf((float) om->max_length / (float) (om->max_length+3));
  float      tloop_total = tloop * om->max_length;
  float      tmove = logf(     3.0f / (float) (om->max_length+3));
  float      tbmk  = logf(     2.0f / ((float) om->M * (float) (om->M+1)));
  float      tec   = logf(1.0f / nu);

  float score_modifier =  tmove + tbmk + tec + tmove + tloop_total;

  FM_METADATA *meta = fm_cfg->meta;


  fm_initSeeds(&seeds);

  /* Set false target length. This is a conservative estimate of the length of window that'll
   * soon be passed on to later phases of the pipeline;  used to recover some bits of the score
   * that we would miss if we left length parameters set to the full target length */
  p7_oprofile_ReconfigMSVLength(om, om->max_length);
  p7_bg_SetLength(bg, om->max_length);
  p7_bg_NullOne  (bg, NULL, om->max_length, &nullsc);

  tmp_sq   =  esl_sq_CreateDigital(om->abc);




  /*
   * Computing the score required to let P meet the F1 prob threshold
   * In original code, converting from an MSV score S (the score getting
   * to state C) to a probability goes like this:
   *  S = XMX(L,p7G_C)
   *  usc = S + tmove + tloop_total
   *  P = f ( (usc - nullsc) / eslCONST_LOG2 , mu, lambda)
   *  and XMX(C) was the diagonal score + tmove + tbmk + tec
   * and we're computing the threshold score S, so reverse it:
   *  (usc - nullsc) /  eslCONST_LOG2 = inv_f( P, mu, lambda)
   *  usc = nullsc + eslCONST_LOG2 * inv_f( P, mu, lambda)
   *  S = usc - tmove - tloop_total - tmove - tbmk - tec
   *
   *
   *  when I'm done gathering a diagonal, the score S will be back-converted to usc:
   *    usc = S + tmove + tbmk + tec + tmove + tloop_total
   *
   *
   *  Here, I compute threshold with length model based on max_length.  Usually, the
   *  length of a window returned by this scan will be 2*max_length-1 or longer.  Doesn't
   *  really matter - in any case, both the bg and om models will change with roughly
   *  1 bit for each doubling of the length model, so they offset.
   */
  invP = esl_gumbel_invsurv(F1, om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  sc_thresh =   (invP * eslCONST_LOG2) + nullsc - (tmove + tloop_total + tmove + tbmk + tec);
  j_thresh  =   0 - tmove - tbmk - tec; //the score required to make it worth appending a diagonal to another diagonal

  invP_FM = esl_gumbel_invsurv(P_fm, om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  sc_threshFM = ESL_MIN(fm_cfg->max_scthreshFM,  (invP_FM * eslCONST_LOG2) + nullsc - (tmove + tloop_total + tmove + tbmk + tec) ) ;


  //get diagonals that score above sc_threshFM
  FM_getSeeds(om, gx, sc_threshFM, &seeds, fmf, fmb, fm_cfg, fm_hmmdata);

  //now extend those diagonals to find ones scoring above sc_thresh (or if MSV is required, capture diagonals scoring at least j_thresh)
  for(i=0; i<seeds.count; i++) {
    FM_extendSeed( seeds.diags+i, fmf, fm_hmmdata, fm_cfg, sc_thresh, j_thresh , tmp_sq, windowlist);
  }


  if (sc_thresh >= j_thresh) {
    // MSV, so join up neighboring diags that reach sc_thresh only in combination
    // this will require repeating a pass over each diagonal that finally causes score to meet threshold
    active_diagonal = FALSE;

    for(i=0; i<windowlist->count; i++) {
      window  =  windowlist->windows + i;
      if (active_diagonal &&
          meta->seq_data[prev_window->id].id == meta->seq_data[window->id].id &&
          prev_window->complementarity == window->complementarity &&
          prev_window->n + prev_window->length + om->max_length <= window->id  // windows are close enough to bother merging
       ) {
        if ( prev_window->score + window->score < sc_thresh ) {
          // failed to hit threshold.  update cumulative score and position, and move along
          prev_window->score += window->score - j_thresh;
          prev_window->length = window->n - prev_window->n + 1; //use the length field to keep track of the range between the 1st peak and the final peak
        } else {
          // hit threshold. need to backtrack to figure out exactly when sc_thresh was actually hit, to mirror SIMD-MSV
          float req_score = sc_thresh - prev_window->score - j_thresh; // score required to push window.score + prev_window.score - j_thresh above sc_thresh
          int pos = FM_extendSeed( seeds.diags+i, fmf, fm_hmmdata, fm_cfg, req_score, 10000 , tmp_sq, NULL); // 1000 ensures that sc_thresh < j_thresh, so SSV-like extension is done

          prev_window->length =  pos - prev_window->n + 1; //use the length field to keep track of the range between the 1st peak and the final peak
          windowlist->windows[i] = *prev_window;
          active_diagonal = FALSE;
        }
      } else {
        //starting a new diagonal
        if (window->score > sc_thresh) {
          // single diagonal hit threshold
          windowlist->windows[i] = *window;
          active_diagonal = FALSE;
        } else {
          //didn't hit threshold, wait for the next one
          prev_window = window;
          active_diagonal = TRUE;
        }
      }

    }

  }


  for(i=0; i<windowlist->count; i++)
    windowlist->windows[i].score += score_modifier;


  //filter for biased composition here?

  //update window sizes
  j=0;
  for(i=0; i<windowlist->count; i++) {
    int skipped_cnt_rev = 0;
    int skipped_cnt_fwd = 0;
    window  =  windowlist->windows + i;

    //coordinates in the target sequence db
    k = window->id;
    while ( meta->seq_data[k].id == meta->seq_data[window->id].id)
      k++;

    target_len = meta->seq_data[k-1].start + meta->seq_data[k-1].length - 1;

    window_start = ESL_MAX( 1,            window->n + window->length - om->max_length) ;
    window_end   = ESL_MIN( target_len ,  window->n + window->length + om->max_length - 2);


    //Then scan backward to figure out how many characters in the fm index correspond to the window
    // (if a window spans more than one seq_data block, there will need to be a bunch of inserted Ns in the final sequence)
    k = window->id;
    while (window_start < meta->seq_data[window->id].start) {
      skipped_cnt_rev += meta->seq_data[k].start - (meta->seq_data[k-1].start + meta->seq_data[k-1].length - 1);
      k--;
    }

    //Then scan forward to figure out how many characters in the fm index correspond to the window
    k = window->id;
    while (meta->seq_data[k].start + meta->seq_data[k].length < window_end) {
      skipped_cnt_fwd += meta->seq_data[k+1].start - (meta->seq_data[k].start + meta->seq_data[k].length - 1);
      k++;
    }

    //Step back along the sequence segments, until finding the first one that contains the character that starts the sequence
    while (window_start < meta->seq_data[window->id].start)
      window->id--;


    //switch coordinates to fm_index
    window_start_fm =   window->fm_n - (window->n - window_start) + skipped_cnt_rev ;
    window_end_fm   =   window->fm_n + (window_end - window->n) - skipped_cnt_fwd;

    window->n      = window_start;
    window->length = window_end - window_start + 1;
    window->fm_n   = window_start_fm;


  }

  //filter for biased composition here?

  //Merge overlapping windows
  if (windowlist->count > 0) {
    j=1;
    for(i=1; i<windowlist->count; i++) {
      prev  =  windowlist->windows + i-1;
      curr  =  windowlist->windows + i;

      if ( prev->n + prev->length - 1 >= curr->n  &&
           prev->complementarity == curr->complementarity &&
           meta->seq_data[prev->id].id == meta->seq_data[curr->id].id

          ) { //overlap , so extend the previous window, and update its score

          prev->length = curr->n - prev->n + 1 + curr->length;

      } else {
        //no overlap, so shift the current window over to the next active slot
        windowlist->windows[j++] = windowlist->windows[i];

      }

    }
    windowlist->count = j;
  }

  return eslEOF;

//ERROR:
//  ESL_EXCEPTION(eslEMEM, "Error allocating memory for hit list\n");

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
  p7_profile_ConfigUnilocal(gm, hmm, bg, L); /* unihit local */

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
  for (k = 1; k <  g2->M; k++) P7P_TSC(g2, k, p7P_MM)  = 0.0f;
  for (k = 0; k <  g2->M; k++) P7P_TSC(g2, k, p7P_BLM) = log(2.0f / ((float) g2->M * (float) (g2->M+1)));

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

  if ((abc = esl_alphabet_Create(eslAMINO))    == NULL)  esl_fatal("failed to create alphabet");
  if (p7_hmm_Sample(r, M, abc, &hmm)           != eslOK) esl_fatal("failed to sample an HMM");
  if ((bg = p7_bg_Create(abc))                 == NULL)  esl_fatal("failed to create null model");
  if ((gm = p7_profile_Create(hmm->M, abc))    == NULL)  esl_fatal("failed to create profile");
  if (p7_profile_Config   (gm, hmm, bg)        != eslOK) esl_fatal("failed to config profile");
  if (p7_profile_SetLength(gm, L)              != eslOK) esl_fatal("failed to config profile length model");
  if (p7_hmm_Validate    (hmm, errbuf, 0.0001) != eslOK) esl_fatal("whoops, HMM is bad!: %s", errbuf);
  if (p7_profile_Validate(gm,  errbuf, 0.0001) != eslOK) esl_fatal("whoops, profile is bad!: %s", errbuf);

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
  p7_profile_ConfigLocal(gm, hmm, bg, 400);

  /* Allocate matrix */
  fwd = p7_gmx_Create(gm->M, sq->n);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_profile_SetLength(gm,  sq->n);
      p7_bg_SetLength     (bg,  sq->n);
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
