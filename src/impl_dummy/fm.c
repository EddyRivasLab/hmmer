#include "hmmer.h"
#include "impl_dummy.h"


/* Function:  fm_initConfig()
 * do nothing.  This is a placeholder in the dummy implementation for
 * a function that does something in _sse and _vmx code
 */
int
fm_initConfig( FM_CFG *cfg ) {
  cfg->maskSA       =  cfg->meta->freq_SA - 1;
  cfg->shiftSA      =  cfg->meta->SA_shift;

  //bounding cutoffs
  cfg->max_depth = 16;
  cfg->neg_len_limit = 4;
  cfg->consec_pos_req = 5; //6
  cfg->score_ratio_req = 0.45; //.49

  return eslOK;
}

/* Function:  fm_destroyConfig()
 * do nothing.  This is a placeholder in the dummy implementation for
 * a function that does something in _sse and _vmx code
 */
int
fm_destroyConfig(FM_CFG *cfg ) {
  return eslOK;
}


/* Function:  fm_getOccCount()
 * Synopsis:  Compute number of occurrences of c in BWT[1..pos]
 *
 * Purpose:   Scan through the BWT to compute number of occurrence of c in BWT[0..pos].
 *
 *            First, use checkpointed occurrence counts in the arrays occCnts_sb and occCnts_b.
 *            The checkpoint used is the one closest to pos, possibly requiring that counts be
 *            subtracted from the checkpointed total
 *
 *            Counting is done by simply scanning through the BWT serially. This can be made
 *            faster with slick bit-fiddling, but speed isn't the issue here, so it's
 *            just a slow (and correct) implementation
 */
int
fm_getOccCount (const FM_DATA *fm, FM_CFG *cfg, int pos, uint8_t c) {

  int i;
  FM_METADATA *meta = cfg->meta;

  int cnt;
  const int b_pos          = (pos+1) >> meta->cnt_shift_b; //floor(pos/b_size)   : the b count element preceding pos
  const uint16_t * occCnts_b  = fm->occCnts_b;
  const uint32_t * occCnts_sb = fm->occCnts_sb;
  const int sb_pos         = (pos+1) >> meta->cnt_shift_sb; //floor(pos/sb_size) : the sb count element preceding pos


  const int cnt_mod_mask_b = meta->freq_cnt_b - 1; //used to compute the mod function
  const int b_rel_pos      = (pos+1) & cnt_mod_mask_b; // pos % b_size      : how close is pos to the boundary corresponding to b_pos
  const int up_b           = b_rel_pos>>(meta->cnt_shift_b - 1); //1 if pos is expected to be closer to the boundary of b_pos+1, 0 otherwise
  const int landmark       = ((b_pos+up_b)<<(meta->cnt_shift_b)) - 1 ;
  // get the cnt stored at the nearest checkpoint
  cnt =  FM_OCC_CNT(sb, sb_pos, c );

  if (up_b)
    cnt += FM_OCC_CNT(b, b_pos + 1, c ) ;
  else if ( b_pos !=  sb_pos * (1<<(meta->cnt_shift_sb - meta->cnt_shift_b)) )
    cnt += FM_OCC_CNT(b, b_pos, c )  ;// b_pos has cumulative counts since the prior sb_pos - if sb_pos references the same count as b_pos, it'll doublecount


  if ( landmark < fm->N || landmark == -1) {

    const uint8_t * BWT = fm->BWT;


    if (meta->alph_type == fm_DNA) {

      if (!up_b) { // count forward, adding
        for (i=1+floor(landmark/4.0); i<(pos+1)/4 ;  i++) {// floor to allow for the case of landmark = -1
          if ((BWT[i] & 0xc0)>>6 == c)  cnt++; //11000000
          if ((BWT[i] & 0x30)>>4 == c)  cnt++; //00110000
          if ((BWT[i] & 0x0c)>>2 == c)  cnt++; //00001100
          if ((BWT[i] & 0x03)    == c)  cnt++; //00000011
        }
        int remaining_cnt = pos + 1 -  i*4 ;
        if (remaining_cnt >= 1)
          if ((BWT[i] & 0xc0)>>6 == c)  cnt++; //11000000
        if (remaining_cnt >= 2)
          if ((BWT[i] & 0x30)>>4 == c)  cnt++; //00110000
        if (remaining_cnt >= 3)
          if ((BWT[i] & 0x0c)>>2 == c)  cnt++; //00001100

      } else { // count backwards, subtracting
        for (i=landmark/4; i>pos/4 ; i--) {
          if ((BWT[i] & 0xc0)>>6 == c)  cnt--; //11000000
          if ((BWT[i] & 0x30)>>4 == c)  cnt--; //00110000
          if ((BWT[i] & 0x0c)>>2 == c)  cnt--; //00001100
          if ((BWT[i] & 0x03)    == c)  cnt--; //00000011
        }
        int remaining_cnt = 3 + 4*i - pos;
        if (remaining_cnt >= 1)
          if ((BWT[i] & 0x03)    == c)  cnt--; //00000011
        if (remaining_cnt >= 2)
          if ((BWT[i] & 0x0c)>>2 == c)  cnt--; //00001100
        if (remaining_cnt >= 3)
          if ((BWT[i] & 0x30)>>4 == c)  cnt--; //00110000

      }
    } else if ( meta->alph_type == fm_DNA_full) {
      if (!up_b) { // count forward, adding
        for (i=1+floor(landmark/2.0); i<(pos+1)/2 ;  i++) {// floor to allow for the case of landmark = -1
          if ((BWT[i] & 0xf0)>>4 == c)  cnt++;
          if ((BWT[i] & 0x0f)    == c)  cnt++;
        }

        if ( !(pos & 0x1) ) {// pos is even, so there's a final singleton
          if ((BWT[i] & 0xf0)>>4 == c)  cnt++;
        }
      } else { // count backwards, subtracting
        for (i=landmark/2; i>pos/2 ; i--) {
          if ((BWT[i] & 0xf0)>>4 == c)  cnt--;  // BWT[i] contains two chars, compressed into one bit
          if ((BWT[i] & 0x0f)    == c)  cnt--;
        }
        if (!(pos & 0x1)) { // pos is even, so there's a final singleton
          if ((BWT[i] & 0x0f) == c)  cnt--;
        }
      }
    } else {
      esl_fatal("Invalid alphabet type\n");
    }
  }

  if (c==0 && pos >= fm->term_loc) { // I overcounted 'A' by one, because '$' was replaced with an 'A'
    cnt--;
  }

  return cnt;

}


/* Function:  fm_getOccCountLT()
 * Synopsis:  Compute number of occurrences of characters with value <c in BWT[1..pos]
 *
 * Purpose:   Scan through the BWT to compute number of occurrence of characters with value <c
 *            in BWT[0..pos].
 *
 *            First, use checkpointed occurrence counts in the arrays occCnts_sb and occCnts_b.
 *            The checkpoint used is the one closest to pos, possibly requiring that counts be
 *            subtracted from the checkpointed total
 *
 *            Counting is done by simply scanning through the BWT serially. This can be made
 *            faster with slick bit-fiddling, but speed isn't the issue here, so it's
 *            just a slow (and correct) implementation
 */
int
fm_getOccCountLT (const FM_DATA *fm, FM_CFG *cfg, int pos, uint8_t c, uint32_t *cnteq, uint32_t *cntlt) {


  if (c == 0 && pos >= fm->term_loc)// < 'A'?  cntlt depends on relationship of pos and the position where the '$' was replaced by 'A'
    *cntlt = 1;
  else
    *cntlt = 0;

  int i;
  FM_METADATA *meta = cfg->meta;

  //int cnt;
  const int b_pos          = (pos+1) >> meta->cnt_shift_b; //floor(pos/b_size)   : the b count element preceding pos
  const uint16_t * occCnts_b  = fm->occCnts_b;
  const uint32_t * occCnts_sb = fm->occCnts_sb;
  const int sb_pos         = (pos+1) >> meta->cnt_shift_sb; //floor(pos/sb_size) : the sb count element preceding pos


  const int cnt_mod_mask_b = meta->freq_cnt_b - 1; //used to compute the mod function
  const int b_rel_pos      = (pos+1) & cnt_mod_mask_b; // pos % b_size      : how close is pos to the boundary corresponding to b_pos
  const int up_b           = b_rel_pos>>(meta->cnt_shift_b - 1); //1 if pos is expected to be closer to the boundary of b_pos+1, 0 otherwise
  const int landmark       = ((b_pos+up_b)<<(meta->cnt_shift_b)) - 1 ;

  // get the cnt stored at the nearest checkpoint
  *cnteq = FM_OCC_CNT(sb, sb_pos, c );
  for (i=0; i<c; i++)
    *cntlt += FM_OCC_CNT(sb, sb_pos, i );

  if (up_b) {
    *cnteq += FM_OCC_CNT(b, b_pos + 1, c ) ;
    for (i=0; i<c; i++)
      *cntlt += FM_OCC_CNT(b, b_pos + 1, i ) ;
  } else if ( b_pos !=  sb_pos * (1<<(meta->cnt_shift_sb - meta->cnt_shift_b)) ) {
    *cnteq += FM_OCC_CNT(b, b_pos, c )  ;// b_pos has cumulative counts since the prior sb_pos - if sb_pos references the same count as b_pos, it'll doublecount
    for (i=0; i<c; i++)
      *cntlt += FM_OCC_CNT(b, b_pos, i ) ;
  }




  if ( landmark < fm->N || landmark == -1) {

    const uint8_t * BWT = fm->BWT;

    if (meta->alph_type == fm_DNA) {

      if (!up_b) { // count forward, adding
        for (i=1+floor(landmark/4.0); i<(pos+1)/4 ;  i++) {// floor to allow for the case of landmark = -1
          if ((BWT[i] & 0xc0)>>6 <  c)  (*cntlt)++; //11000000
          if ((BWT[i] & 0xc0)>>6 == c)  (*cnteq)++; //11000000
          if ((BWT[i] & 0x30)>>4 <  c)  (*cntlt)++; //00110000
          if ((BWT[i] & 0x30)>>4 == c)  (*cnteq)++; //00110000
          if ((BWT[i] & 0x0c)>>2 <  c)  (*cntlt)++; //00001100
          if ((BWT[i] & 0x0c)>>2 == c)  (*cnteq)++; //00001100
          if ((BWT[i] & 0x03)    <  c)  (*cntlt)++; //00000011
          if ((BWT[i] & 0x03)    == c)  (*cnteq)++; //00000011
        }
        int remaining_cnt = pos + 1 -  i*4 ;
        if (remaining_cnt >= 1) {
          if ((BWT[i] & 0xc0)>>6 <  c)  (*cntlt)++; //11000000
          if ((BWT[i] & 0xc0)>>6 == c)  (*cnteq)++; //11000000
        }
        if (remaining_cnt >= 2) {
          if ((BWT[i] & 0x30)>>4 <  c)  (*cntlt)++; //00110000
          if ((BWT[i] & 0x30)>>4 == c)  (*cnteq)++; //00110000
        }
        if (remaining_cnt >= 3) {
          if ((BWT[i] & 0x0c)>>2 <  c)  (*cntlt)++; //00001100
          if ((BWT[i] & 0x0c)>>2 == c)  (*cnteq)++; //00001100
        }

      } else { // count backwards, subtracting
        for (i=landmark/4; i>pos/4 ; i--) {
          if ((BWT[i] & 0xc0)>>6 <  c)  (*cntlt)--; //11000000
          if ((BWT[i] & 0xc0)>>6 == c)  (*cnteq)--; //11000000
          if ((BWT[i] & 0x30)>>4 <  c)  (*cntlt)--; //00110000
          if ((BWT[i] & 0x30)>>4 == c)  (*cnteq)--; //00110000
          if ((BWT[i] & 0x0c)>>2 <  c)  (*cntlt)--; //00001100
          if ((BWT[i] & 0x0c)>>2 == c)  (*cnteq)--; //00001100
          if ((BWT[i] & 0x03)    <  c)  (*cntlt)--; //00000011
          if ((BWT[i] & 0x03)    == c)  (*cnteq)--; //00000011
        }
        int remaining_cnt = 3 + 4*i - pos;
        if (remaining_cnt >= 1) {
          if ((BWT[i] & 0x03)    <  c)  (*cntlt)--; //00000011
          if ((BWT[i] & 0x03)    == c)  (*cnteq)--; //00000011
        }
        if (remaining_cnt >= 2) {
          if ((BWT[i] & 0x0c)>>2 <  c)  (*cntlt)--; //00001100
          if ((BWT[i] & 0x0c)>>2 == c)  (*cnteq)--; //00001100
        }
        if (remaining_cnt >= 3) {
          if ((BWT[i] & 0x30)>>4 <  c)  (*cntlt)--; //00110000
          if ((BWT[i] & 0x30)>>4 == c)  (*cnteq)--; //00110000
        }
      }
    } else if ( meta->alph_type == fm_DNA_full) {
      if (!up_b) { // count forward, adding
        for (i=1+floor(landmark/2.0); i<(pos+1)/2 ;  i++) {// floor to allow for the case of landmark = -1
          if ((BWT[i] & 0xf0)>>4 <  c)  (*cntlt)++;
          if ((BWT[i] & 0xf0)>>4 == c)  (*cnteq)++;
          if ((BWT[i] & 0x0f)    <  c)  (*cntlt)++;
          if ((BWT[i] & 0x0f)    == c)  (*cnteq)++;
        }

        if ( !(pos & 0x1) ) {// pos is even, so there's a final singleton
          if ((BWT[i] & 0xf0)>>4 <  c)  (*cntlt)++;
          if ((BWT[i] & 0xf0)>>4 == c)  (*cnteq)++;
        }
      } else { // count backwards, subtracting
        for (i=landmark/2; i>pos/2 ; i--) {
          if ((BWT[i] & 0xf0)>>4 <  c)  (*cntlt)--;  // BWT[i] contains two chars, compressed into one bit
          if ((BWT[i] & 0xf0)>>4 == c)  (*cnteq)--;  // BWT[i] contains two chars, compressed into one bit
          if ((BWT[i] & 0x0f)    <  c)  (*cntlt)--;
          if ((BWT[i] & 0x0f)    == c)  (*cnteq)--;
        }
        if (!(pos & 0x1)) { // pos is even, so there's a final singleton
          if ((BWT[i] & 0x0f) <  c)  (*cntlt)--;
          if ((BWT[i] & 0x0f) == c)  (*cnteq)--;
        }
      }
    } else {
      esl_fatal("Invalid alphabet type\n");
    }
  }

  if (c==0 && pos >= fm->term_loc) { // I overcounted 'A' by one, because '$' was replaced with an 'A'
    (*cnteq)--;
  }

  return eslOK;

}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
