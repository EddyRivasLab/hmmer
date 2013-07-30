#include "p7_config.h"

#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_sq.h"

#include "hmmer.h"


/* hit_sorter(): qsort's pawn, below */
static int
fm_hit_sorter(const void *a, const void *b)
{
    FM_DIAG *d1 = (FM_DIAG*)a;
    FM_DIAG *d2 = (FM_DIAG*)b;

    return 2 * (d1->sortkey > d2->sortkey) - 1;  // same as the test below
//    if      ( d1->sortkey > d2->sortkey) return 1;
//    else                                 return -1;

}

static int
fm_mergeSeeds(FM_DIAGLIST *seeds, int N, int ssv_length) {
  int i;
  int j = 0;

  FM_DIAG next;
  int tmp;
  int next_is_complement;
  int curr_is_complement;
  int compl_mult;
  int curr_n;
  int curr_k;
  int curr_len;
  int curr_end;
  int curr_diagval;

  FM_DIAG *diags = seeds->diags;


  if (seeds->count == 0)
    return eslOK;

  //sort, first by direction, then N (position on database sequence), then K (model position)
  qsort(diags, seeds->count, sizeof(FM_DIAG), fm_hit_sorter);

  for (i=0; i<seeds->count; i++)
    printf("n: %d, k: %d, s: %3f\n", diags[i].n, diags[i].k, diags[i].sortkey);


  next = diags[0];

  curr_is_complement = (next.complementarity == fm_complement);
  compl_mult         = curr_is_complement ? 1 : -1;
  compl_mult         = -1;
  curr_n             = next.n;
  curr_k             = next.k;
  curr_len           = next.length;
  curr_end           = curr_n + curr_len - 1;
  curr_diagval       = next.n + next.k * compl_mult;


  for( i=1; i<seeds->count; i++) {

    next = diags[i];
    next_is_complement = (next.complementarity == fm_complement);
    //compl_mult         = next_is_complement ? 1 : -1;
    compl_mult = -1;

    //overlapping diagonals will share the same value of (n - k) for non-complement hits, and (n+k) for complement hits
    if (  next_is_complement == curr_is_complement              //same direction
          && ( next.n + next.k * compl_mult) == curr_diagval    //same diagonal
          && next.n + next.length < curr_n + curr_len + ssv_length       //overlapping, or close to it
    ) {


      //overlapping diags; extend, if appropriate
        tmp = next.n + next.length - 1;
        if (tmp > curr_end) {
          curr_end = tmp;
          curr_len = curr_end - curr_n + 1;
        }
    } else {
      //doesn't overlap current diagonal, so store current...
      diags[j].n      = curr_n;
      diags[j].k      = curr_k;
      diags[j].length = curr_end - curr_n + 1;
      diags[j].complementarity = curr_is_complement;
      diags[j].score = 0.0;

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
  diags[j].score  = 0.0;

  seeds->count = j+1;
  return eslOK;

}


static uint32_t
FM_backtrackSeed(const FM_DATA *fmf, FM_CFG *fm_cfg, int i, FM_DIAG *seed) {
  int j = i;
  int len = 0;
  int c;

  while ( j != fmf->term_loc && (j % fm_cfg->meta->freq_SA)) { //go until we hit a position in the full SA that was sampled during FM index construction
    c = fm_getChar( fm_cfg->meta->alph_type, j, fmf->BWT);
    //printf("c: %d  (j=%d)\n", c,j );
    j = fm_getOccCount (fmf, fm_cfg, j-1, c);
    //printf ("  cnt1: %d\n", j);
    j += abs(fmf->C[c]);
    //printf ("  cnt2: %d\n", j);
    len++;
  }

//  printf("stop at j: %d (%d)\n", j, fmf->term_loc);

  return len + (j==fmf->term_loc ? 0 : fmf->SA[ j / fm_cfg->meta->freq_SA ]) ; // len is how many backward steps we had to take to find a sampled SA position

}

static int
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

    if (complementarity == fm_nocomplement )
      seed->n    =  fmf->N - FM_backtrackSeed(fmf, fm_cfg, i, seed) - depth;
    else
      seed->n    =  FM_backtrackSeed(fmf, fm_cfg, i, seed) + 1;

    seed->complementarity = complementarity;

    /* seed->n corresponds to the start of the seed in terms of the
     * target sequence, in a forward direction.  seed->k holds the
     * model position at the beginning of that seed.  If model_direction
     * is fm_reverse, the ->n value is from the beginning of the revcomp,
     */
    if (model_direction == fm_forward)
      seed->k -= (depth - 1) ;


//    printf ("n: %d;  k: %d;  len: %d\n", seed->n, seed->k, seed->length);

    seed->sortkey =  (int)( complementarity == fm_complement ? fmf->N + 1 : 0)   // makes complement seeds cover a different score range than non-complements
                    +  ((int)(seed->n) - (int)(seed->k) )                        // unique diagonal within the complement/non-complement score range
                    + ((double)(seed->k)/(double)(M+1))  ;                       // fractional part, used to sort seeds sharing a diagonal
  }

  return eslOK;
}


static int
FM_Recurse( int depth, int M, int Kp, int fm_direction,
            const FM_DATA *fmf, const FM_DATA *fmb,
            FM_CFG *fm_cfg, const P7_SCOREDATA *ssvdata,
            int first, int last, float sc_threshFM,
            FM_DP_PAIR *dp_pairs,
            FM_INTERVAL *interval_1, FM_INTERVAL *interval_2,
            FM_DIAGLIST *seeds,
            char *seq
          )
{


  float sc, next_score;

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
          next_score = ssvdata->ssv_scores_f[k*Kp + fm_getComplement(c,fm_cfg->meta->alph_type) ];
        } else
          next_score = ssvdata->ssv_scores_f[k*Kp + c];

        sc = dp_pairs[i].score + next_score;

  //      printf( "%-18s : (%.2f >=? %.2f)\n", seq, sc, sc_threshFM);

        if ( sc >= sc_threshFM ) { // this is a seed I want to extend

          interval_1_new.lower = interval_1->lower;
          interval_1_new.upper = interval_1->upper;

          if (fm_direction == fm_forward) {
            if ( interval_1_new.lower >= 0 && interval_1_new.lower <= interval_1_new.upper  )  //no use extending a non-existent string
            fm_updateIntervalReverse( fmf, fm_cfg, c, &interval_1_new);

//            printf( "%-18s : (%d - %d) (%.2f)   (pass F)\n", seq, interval_1_new.lower, interval_1_new.upper, sc);

            if ( interval_1_new.lower >= 0 && interval_1_new.lower <= interval_1_new.upper  )  //no use passing a non-existent string
              FM_getPassingDiags(fmf, fm_cfg, k, M, sc, depth, fm_forward,
                                 dp_pairs[i].model_direction, dp_pairs[i].complementarity,
                                 &interval_1_new, seeds);

          } else {
            //searching for forward matches on the FM-index
            interval_2_new.lower = interval_2->lower;
            interval_2_new.upper = interval_2->upper;

            //searching for reverse matches on the FM-index
            if ( interval_1_new.lower >= 0 && interval_1_new.lower <= interval_1_new.upper  )  //no use extending a non-existent string
              fm_updateIntervalForward( fmb, fm_cfg, c, &interval_1_new, &interval_2_new);

            //printf( "%-18s : (%d - %d) (%.2f)   (pass R) (%d, %d)\n", seq, interval_1_new.lower, interval_1_new.upper, sc, interval_2_new.lower, interval_2_new.upper);


            if ( interval_2_new.lower >= 0 && interval_2_new.lower <= interval_2_new.upper  )  //no use passing a non-existent string
              FM_getPassingDiags(fmf, fm_cfg, k, M, sc, depth, fm_backward,
                                 dp_pairs[i].model_direction, dp_pairs[i].complementarity,
                                 &interval_2_new, seeds);

          }

        } else if (  sc <= 0                                                                                        //some other path in the string enumeration tree will do the job
            || depth == fm_cfg->max_depth                                                                            //can't extend anymore, 'cause we've reached the pruning length
            || (depth == dp_pairs[i].max_score_len + fm_cfg->neg_len_limit)                                        //too many consecutive positions with a negative total score contribution (sort of like Xdrop)
            || ((float)sc/(float)depth < fm_cfg->score_ratio_req)                                                                //score density is too low
            || (dp_pairs[i].max_consec_pos < fm_cfg->consec_pos_req  &&                                            //a seed is expected to have at least one run of positive-scoring matches at least length consec_pos_req;  if it hasn't,  (see Tue Nov 23 09:39:54 EST 2010)
                   ( (depth >= fm_cfg->max_depth/2 &&  sc/depth < sc_threshFM/fm_cfg->max_depth)                              // if we're close to the end of the sequence, abort -- if that end does have sufficiently long all-positive run, I'll find it on the reverse sweep
                   || depth == fm_cfg->max_depth-fm_cfg->consec_pos_req+1 )                                                   // if we're at least half way across the sequence, and score density is too low, abort -- if the density on the other side is high enough, I'll find it on the reverse sweep
               )
            || (dp_pairs[i].model_direction == fm_forward  &&
                   ( (k == M)                                                                                                          //can't extend anymore, 'cause we're at the end of the model, going forward
                  || (depth > (fm_cfg->max_depth - 10) &&  sc + ssvdata->opt_ext_fwd[k][fm_cfg->max_depth-depth-1] < sc_threshFM)   //can't hit threshold, even with best possible forward extension up to length ssv_req
                  ))
            || (dp_pairs[i].model_direction == fm_backward &&
                   ( (k == 1)                                                                                                          //can't extend anymore, 'cause we're at the beginning of the model, going backwards
                  || (depth > (fm_cfg->max_depth - 10) &&  sc + ssvdata->opt_ext_rev[k][fm_cfg->max_depth-depth-1] < sc_threshFM )  //can't hit threshold, even with best possible extension up to length ssv_req
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

        if ( interval_1_new.lower >= 0 && interval_1_new.lower <= interval_1_new.upper  )  //no use extending a non-existent string
          fm_updateIntervalReverse( fmf, fm_cfg, c, &interval_1_new);

        if (  interval_1_new.lower < 0 || interval_1_new.lower > interval_1_new.upper ) { //that string doesn't exist in fwd index
          continue;
        }
        FM_Recurse(depth+1, M, Kp, fm_direction,
                  fmf, fmb, fm_cfg, ssvdata,
                  last+1, dppos, sc_threshFM, dp_pairs,
                  &interval_1_new, NULL,
                  seeds,
                  seq
                  );


      } else {

        interval_2_new.lower = interval_2->lower;
        interval_2_new.upper = interval_2->upper;

        if ( interval_1_new.lower >= 0 && interval_1_new.lower <= interval_1_new.upper  )  //no use extending a non-existent string
          fm_updateIntervalForward( fmb, fm_cfg, c, &interval_1_new, &interval_2_new);

        if (  interval_1_new.lower < 0 || interval_1_new.lower > interval_1_new.upper ) { //that string doesn't exist in reverse index
          continue;
        }
        FM_Recurse(depth+1, M, Kp, fm_direction,
                  fmf, fmb, fm_cfg, ssvdata,
                  last+1, dppos, sc_threshFM, dp_pairs,
                  &interval_1_new, &interval_2_new,
                  seeds,
                  seq
                  );

      }
    } else {
//      printf ("internal stop: %d\n", depth);
    }

  }

  return eslOK;
}


static int FM_getSeeds (const P7_OPROFILE *gm, float sc_threshFM,
                FM_DIAGLIST *seeds,
                const FM_DATA *fmf, const FM_DATA *fmb,
                FM_CFG *fm_cfg, const P7_SCOREDATA *ssvdata )
{
  FM_INTERVAL interval_f1, interval_f2, interval_bk;
  int i, k;
  int status;
  float sc;
  char         *seq;

  FM_DP_PAIR *dp_pairs_fwd;
  FM_DP_PAIR *dp_pairs_rev;

  ESL_ALLOC(dp_pairs_fwd, gm->M * fm_cfg->max_depth * sizeof(FM_DP_PAIR)); // guaranteed to be enough to hold all diagonals
  ESL_ALLOC(dp_pairs_rev, gm->M * fm_cfg->max_depth * sizeof(FM_DP_PAIR));

  ESL_ALLOC(seq, 50*sizeof(char));


  for (i=0; i<fm_cfg->meta->alph_size; i++) {
    int fwd_cnt=0;
    int rev_cnt=0;

    interval_f1.lower = interval_f2.lower = interval_bk.lower = fmf->C[i];
    interval_f1.upper = interval_f2.upper = interval_bk.upper = abs(fmf->C[i+1])-1;

    if (interval_f1.lower<0 ) //none of that character found
      continue;


    //c = i + (i >= gm->abc->K ? 1 : 0);  // shift to the easel alphabet, where the Kth position is '-'
    seq[0] = fm_cfg->meta->alph[i];
    seq[1] = '\0';


//    printf("%-18s : %d , %d\n", seq, interval_f1.lower, interval_f1.upper);

    // Fill in a DP column for the character c, (compressed so that only positive-scoring entries are kept)
    // There will be 4 DP columns for each character, (1) fwd-std, (2) fwd-complement, (3) rev-std, (4) rev-complement
    for (k = 1; k <= gm->M; k++) // there's no need to bother keeping an entry starting at the last position (gm->M)
    {


      sc = ssvdata->ssv_scores_f[k*gm->abc->Kp + i];
      if (sc>0) { // we'll extend any positive-scoring diagonal
        /* fwd on model, fwd on FM (really, reverse on FM, but the FM is on a reversed string, so its fwd*/
/*
        if (k < gm->M-3) { // don't bother starting a forward diagonal so close to the end of the model
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
*/

        /* rev on model, rev on FM (the FM is on the unreversed string)*/
/*
        if (k > 4) { // don't bother starting a reverse diagonal so close to the start of the model
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
*/
      }


      // Now do the reverse complement
      sc = ssvdata->ssv_scores_f[k*gm->abc->Kp + fm_getComplement(i, fm_cfg->meta->alph_type)];
      if (sc>0) { // we'll extend any positive-scoring diagonal
        /* rev on model, fwd on FM (really, reverse on FM, but the FM is on a reversed string, so its fwd*/
/*
        if (k > 4) { // don't bother starting a reverse diagonal so close to the start of the model
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
*/

        /* fwd on model, rev on FM (the FM is on the unreversed string - complemented)*/

        if (k < gm->M-3) { // don't bother starting a forward diagonal so close to the end of the model
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

    FM_Recurse ( 2, gm->M, gm->abc->Kp, fm_forward,
                 fmf, fmb, fm_cfg, ssvdata,
                 0, fwd_cnt-1, sc_threshFM, dp_pairs_fwd,
                 &interval_f1, NULL,
                 seeds
                 , seq
            );

    FM_Recurse ( 2, gm->M, gm->abc->Kp, fm_backward,
                 fmf, fmb, fm_cfg, ssvdata,
                 0, rev_cnt-1, sc_threshFM, dp_pairs_rev,
                 &interval_bk, &interval_f2,
                 seeds
                 , seq
            );

  }


  //merge duplicates
  fm_mergeSeeds(seeds, fmf->N, fm_cfg->ssv_length);


  free (dp_pairs_fwd);
  free (dp_pairs_rev);
  if (seq) free(seq);
  return eslOK;

ERROR:
  ESL_EXCEPTION(eslEMEM, "Error allocating memory for hit list\n");

}

static int
FM_window_from_diag (FM_DIAG *diag, const FM_DATA *fm, const FM_METADATA *meta, P7_HMM_WINDOWLIST *windowlist) {

  int status;
  int again = TRUE;
  /*  Some of the seeds we just created may span more than one sequence in the target.
   *  Check for this, and resolve it, by trimming such an offending seed, and tacking
   *  the remainder on as a new seed
   */
  while (again) {
    uint32_t seg_id, seg_pos;
    again = FALSE;

    if (diag->complementarity == fm_nocomplement) {
      status = fm_getOriginalPosition (fm, meta, 0, diag->length, fm_forward, diag->n, &seg_id, &seg_pos);
      if (status == eslERANGE) {
        int target_end = meta->seq_data[seg_id].offset + meta->seq_data[seg_id].length - 1;
        int diag_end   = seg_pos + diag->length;
        int overext    = diag_end - target_end;

        p7_hmmwindow_new(windowlist, seg_id, seg_pos, diag->n, diag->k+diag->length-overext-1, diag->length-overext, diag->score, diag->complementarity);

        diag->n      +=  diag->length - overext;
        diag->k      +=  diag->length - overext;
        diag->length = overext;
        again        = TRUE;
      }
    } else {
      status = fm_getOriginalPosition (fm, meta, 0, diag->length, fm_backward, fm->N - diag->n + 1, &seg_id, &seg_pos);
      if (status == eslERANGE) {
        int target_end = meta->seq_data[seg_id].offset - meta->seq_data[seg_id].length + 1;
        int diag_end   = seg_pos - diag->length;
        int overext    = target_end - diag_end;

        p7_hmmwindow_new(windowlist, seg_id, seg_pos, diag->n, diag->k+diag->length-overext-1, diag->length-overext, diag->score, diag->complementarity);

        diag->n      +=  diag->length - overext;
        diag->k      +=  diag->length - overext;
        diag->length = overext;
        again        = TRUE;
      }
    }

    if (!again)
      p7_hmmwindow_new(windowlist, seg_id, seg_pos, diag->n, diag->k+diag->length-1, diag->length, diag->score, diag->complementarity);

  }

  return eslOK;

}

static int
FM_extendSeed(FM_DIAG *diag, const FM_DATA *fm, const P7_SCOREDATA *ssvdata, FM_CFG *cfg, float sc_thresh, ESL_SQ   *tmp_sq)
{
  //extend seed in both diagonal directions,
  //capturing the score

  int k,n;
  int model_start, model_end, target_start, target_end;
  int hit_start, hit_end, max_hit_start, max_hit_end;
  float sc;
  float max_sc = 0;
  int c;

  //this determines the start and end of the window that we think it's possible we'll extend to the window to (which determines the sequence we'll extract)
  model_start      = ESL_MIN( diag->k, ESL_MAX(1, diag->k + diag->length - cfg->ssv_length  )) ; //-1 and +1 cancel
  model_end        = ESL_MIN( ssvdata->M, diag->k + cfg->ssv_length - 1 );
  target_start     = diag->n - (diag->k - model_start);
  target_end       = diag->n + (model_end - diag->k);
  if (target_start < 1) {
    model_start += 1-target_start;
    target_start = 1;
  }
  if (target_end > fm->N-1) {
    model_end -= target_end - (fm->N-1);
    target_end = fm->N-1;
  }


  k = model_start;
  n = 1;
  sc = 0;

  fm_convertRange2DSQ(fm, cfg->meta, target_start, target_end-target_start+1, diag->complementarity, tmp_sq );

  //This finds the highest-scoring sub-diag in the just-determined potential diagonal range.
  hit_start = n;
  for (  ; k <= model_end; k++, n++) {
      c = tmp_sq->dsq[n];

      sc  += ssvdata->ssv_scores_f[k*tmp_sq->abc->Kp + c];

      if (sc < 0) {
        sc = 0;
        hit_start = n+1;
      }
      hit_end = n;

      if (sc > max_sc) {
          max_sc = sc;
          max_hit_start = hit_start;
          max_hit_end   = hit_end;
      }
  }

  if (diag->complementarity == fm_complement) {
    diag->n       = target_end  - max_hit_end + 1;
    diag->k       = max_hit_end;
  } else {
    diag->n       = target_start + max_hit_start - 1;
    diag->k       = model_start  + max_hit_start - 1;
  }
  diag->length  = max_hit_end - max_hit_start + 1;
  diag->score = max_sc;


  return eslOK;
}


/* Function:  p7_SSVFM_longlarget()
 * Synopsis:  Finds windows with SSV scores above given threshold
 *
 * Details:   Uses FM-index to find high-scoring diagonals (seeds), then extends those
 *            seeds to maximal scoring diagonals (no gaps). Windows meeting the SSV
 *            scoring threshold (usually score s.t. p=0.02) are captured, and passed
 *            on to the Viterbi and Forward stages of the pipeline.
 *
 * Args:      om      - optimized profile
 *            nu      - configuration: expected number of hits (use 2.0 as a default)
 *            bg      - the background model, required for translating a P-value threshold into a score threshold
 *            F1      - p-value below which a window is captured as being above threshold
 *            fmf     - data for forward traversal of the FM-index
 *            fmb     - data for backward traversal of the FM-index
 *            fm_cfg  - FM-index meta data
 *            ssvdata - compact data required for computing MSV (SSV) scores
 *            windowlist - RETURN: collection of SSV-passing windows, with meta data required for downstream stages.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_SSVFM_longlarget( P7_OPROFILE *om, float nu, P7_BG *bg, double F1,
         const FM_DATA *fmf, const FM_DATA *fmb, FM_CFG *fm_cfg, const P7_SCOREDATA *ssvdata,
         P7_HMM_WINDOWLIST *windowlist)
{
  float P_fm = 0.5;
  float sc_thresh, sc_threshFM;
  float invP, invP_FM;
  float nullsc;

  int i;

  float      tloop = logf((float) om->max_length / (float) (om->max_length+3));
  float      tloop_total = tloop * om->max_length;
  float      tmove = logf(     3.0f / (float) (om->max_length+3));
  float      tbmk  = logf(     2.0f / ((float) om->M * (float) (om->M+1)));
  float      tec   = logf(1.0f / nu);
  FM_DIAG   *diag;

  ESL_SQ   *tmp_sq;

  FM_DIAGLIST seeds;
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
   * In original code, converting from an SSV score S (the score getting
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
   *  Here, I compute threshold with length model based on max_length.  Usually, the
   *  length of a window returned by this scan will be 2*max_length-1 or longer.  Doesn't
   *  really matter - in any case, both the bg and om models will change with roughly
   *  1 bit for each doubling of the length model, so they offset.
   */
  invP = esl_gumbel_invsurv(F1, om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  sc_thresh =   (invP * eslCONST_LOG2) + nullsc - (tmove + tloop_total + tmove + tbmk + tec);

  invP_FM = esl_gumbel_invsurv(P_fm, om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  sc_threshFM = ESL_MIN(fm_cfg->max_scthreshFM,  (invP_FM * eslCONST_LOG2) + nullsc - (tmove + tloop_total + tmove + tbmk + tec) ) ;

  //get diagonals that score above sc_threshFM
  FM_getSeeds(om, sc_threshFM, &seeds, fmf, fmb, fm_cfg, ssvdata);


  //now extend those diagonals to find ones scoring above sc_thresh
  for(i=0; i<seeds.count; i++) {
    FM_extendSeed( seeds.diags+i, fmf, ssvdata, fm_cfg, sc_thresh, tmp_sq);
  }


  for(i=0; i<seeds.count; i++) {
    diag = seeds.diags+i;
    if (diag->score >= sc_thresh)
      FM_window_from_diag(diag, fmf, fm_cfg->meta, windowlist );
  }



  esl_sq_Destroy(tmp_sq);

  free(seeds.diags);
  return eslEOF;

//ERROR:
//  ESL_EXCEPTION(eslEMEM, "Error allocating memory for hit list\n");

}
/*------------------ end, FM_MSV() ------------------------*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
