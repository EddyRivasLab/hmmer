#include "hmmer.h"
#include "esl_vmx.h"
#include "impl_vmx.h"


int
fm_getbits_vec (vector unsigned char in, char *buf, int reverse) {
  byte_vec new;
  new.u8 = in;
  int i,j;

  for (i=0; i<16;i++) {

    for (j=0; j<8; j++) {
      if (reverse)
        buf[9*i+j] = (new.bytes[i]>>j)&0x1 ? '1' : '0';
      else
        buf[9*i+(7-j)] = (new.bytes[i]>>j)&0x1 ? '1' : '0';
    }
    buf[9*i + 8] = ' ';
  }
  buf[143] = '\0';

  return eslOK;
}

int
fm_print_vec (vector unsigned char in) {
  char str[144];
  fm_getbits_vec(in, str, 0);
  fprintf (stderr, "%s\n", str);
  return eslOK;
}


int
fm_print_vec_rev (vector unsigned char in) {
  char str[144];
  fm_getbits_vec(in, str, 1);
  fprintf (stderr, "%s\n", str);
  return eslOK;
}


/* Function:  fm_initMiscVars()
 * Purpose:   Initialize vector masks used in VMX FMindex implementation
 */
int
fm_initMiscVars( FM_MISC_VARS *misc ) {
  int status;
  int i,j;
  int trim_chunk_count;

  misc->fm_allones_v = esl_vmx_set_u8((unsigned char) 0xff);
  misc->fm_neg128_v  = esl_vmx_set_u8((int8_t) -128);
  misc->fm_zeros_v   = esl_vmx_set_u8((int8_t) 0x00);      //all zeros
  misc->fm_m0f       = esl_vmx_set_u8((int8_t) 0x0f);      //00 00 11 11
  misc->fm_one       = esl_vmx_set_u8((int8_t) 1);         //
  misc->fm_two       = esl_vmx_set_u8((int8_t) 2);         //
  misc->fm_four      = esl_vmx_set_u8((int8_t) 4);         //

  if (meta->alph_type == fm_DNA) {
    misc->fm_m01 = esl_vmx_set_u8((int8_t) 0x55);   //01 01 01 01
    misc->fm_m11 = esl_vmx_set_u8((int8_t) 0x03);  //00 00 00 11
  }
    //set up an array of vectors, one for each character in the alphabet
  misc->fm_chars_v         = NULL;
  ESL_ALLOC (misc->fm_chars_mem, misc->meta->alph_size * sizeof(vector unsigned char)  + 15 ); // +15 for manual 16-byte alignment, which matters for SIMD stuff
     misc->fm_chars_v =   (vector unsigned char *) (((unsigned long int)misc->fm_chars_mem + 15) & (~0xf));   /* align vector memory on 16-byte boundaries */

  for (i=0; i<misc->meta->alph_size; i++) {
    int8_t c = i;
    if (misc->meta->alph_type == fm_DNA) {
      // need 4 copies per byte
      c |= i<<2;
      c |= i<<4;
      c |= i<<6;
    } //else, just leave it on the right-most bits

    misc->fm_chars_v[i] = esl_vmx_set_u8(c);
  }

  /* this is a collection of masks used to clear off the left- or right- part
   *  of a register when we shouldn't be counting the whole thing
   * Incrementally chew off the 1s in chunks of 2 (for DNA) or 4 (for DNA_full)
   * from the right side, and stick each result into an element of a vector
   */
  if (misc->meta->alph_type == fm_DNA)
    trim_chunk_count = 64; //2-bit steps
  else //(meta->alph_type == fm_DNA_full)
    trim_chunk_count = 16; //8-bit steps

  //chars_per_vector = 128/meta->charBits;
  misc->fm_masks_v         = NULL;
  misc->fm_reverse_masks_v = NULL;
  ESL_ALLOC (misc->fm_masks_mem, (1+trim_chunk_count) *sizeof(vector unsigned char)  + 15 ); // +15 for manual 16-byte alignment, which matters for SIMD stuff
      misc->fm_masks_v =   (vector unsigned char *) (((unsigned long int)misc->fm_masks_mem + 15) & (~0xf));   /* align vector memory on 16-byte boundaries */

  ESL_ALLOC (misc->fm_reverse_masks_mem, (1+trim_chunk_count) *sizeof(vector unsigned char)  + 15 ); // +15 for manual 16-byte alignment, which matters for SIMD stuff
      misc->fm_reverse_masks_v =   (vector unsigned char *) (((unsigned long int)misc->fm_reverse_masks_mem + 15) & (~0xf));   /* align vector memory on 16-byte boundaries */


  {
    byte_vec arr;
    arr.u8 = fm_allones_v;

    for (i=trim_chunk_count-1; i>0; i--) {
      int byte_mask=0xff; //11 11 11 11
      int byte_i = (i-1)/(trim_chunk_count/16);
      if (meta->alph_type == fm_DNA) {
        switch (i&0x03) {
          case 1:
            byte_mask = 0xc0; //11 00 00 00
            break;
          case 2:
            byte_mask = 0xf0; //11 11 00 00
            break;
          case 3:
            byte_mask = 0xfc; //11 11 11 00
            break;
          default:
            break;
        }
      }

      arr.bytes[byte_i] = byte_mask; //chew off the appropriate number of bits
      for (j=byte_i+1; j<16; j++) {
        arr.bytes[j] = 0x0;
      }
      misc->fm_masks_v[i]                           = *(vector unsigned char*)(&(arr.u8));
      misc->fm_reverse_masks_v[trim_chunk_count-i]  = vec_andc( misc->fm_allones_v , misc->fm_masks_v[i]);   //_mm_andnot_si128(fm_masks_v[i], fm_allones_v );

    }
  }

  if (meta->alph_type == fm_DNA_full) {
    misc->fm_masks_v[16]          = misc->fm_allones_v;
    misc->fm_reverse_masks_v[16]  = misc->fm_allones_v;

  }
/*
  fprintf( stderr, "mask\n========\n");
  for (i=0; i<trim_chunk_count; i++) {
      fprintf( stderr, "%2d:", i);
            fm_print_vec(fm_masks_v[i]);
        }

  fprintf( stderr, "\n\n reverse mask\n========\n");
  for (i=0; i<trim_chunk_count; i++) {
      fprintf( stderr, "%2d:", i);
            fm_print_vec(fm_reverse_masks_v[i]);
        }
exit(1);
*/
  return eslOK;

ERROR:
  if (fm_chars_mem)         free(fm_chars_mem);
  if (fm_masks_mem)         free(fm_masks_mem);
  if (fm_reverse_masks_mem) free(fm_reverse_masks_mem);

  esl_fatal("Error allocating memory in initGlobals\n");
  return eslFAIL;
}


/* Function:  fm_destroyMiscVars()
 * Purpose:   Destroy vector masks used in SSE FMindex implementation
 */
int
fm_destroyMiscVars(FM_MISC_VARS *misc ) {
  if (misc) {
    if (misc->fm_chars_mem)         free(misc->fm_chars_mem);
    if (misc->fm_masks_mem)         free(misc->fm_masks_mem);
    if (misc->fm_reverse_masks_mem) free(misc->fm_reverse_masks_mem);
  }
  return eslOK;
}

/* Function:  fm_getOccCount()
 * Synopsis:  Compute number of occurrences of c in BWT[1..pos]
 *
 * Purpose:   Scan through the BWT to compute number of occurrence of c in BWT[0..pos],
 *            using altivec instructions to scan 16 bytes-at-a-time.
 *
 *            First, use checkpointed occurrence counts in the arrays occCnts_sb and occCnts_b.
 *            The checkpoint used is the one closest to pos, possibly requiring that counts be
 *            subtracted from the checkpointed total
 *
 *            The counting method is SIMD, loading 16 bytes (32 or 64 chars, depending on
 *            alphabet) at a time into the vector co-processors, then counting occurrences. One
 *            constraint of this approach is that occCnts_b checkpoints must be spaced at least
 *            every 32 or 64 chars (16 bytes, in pressed format), and in multiples of 64/32, so
 *            that _mm_load_si128 calls appropriately meet 16-byte-alignment requirements. That's
 *            a reasonable expectation, as spacings of 256 or more seem to give the best speed,
 *            and certainly better space-utilization.
 */
int
fm_getOccCount (FM_DATA *fm, FM_MISC_VARS *misc, int pos, uint8_t c) {

  int i;
  FM_METADATA *meta = misc->meta;

  int cnt, cnt2;
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


  if ( landmark < fm->N  || landmark == -1) {

    const uint8_t * BWT = fm->BWT;


    vector unsigned char c_v = *(misc->fm_chars_v + c);
    vector unsigned char BWT_v;
    vector unsigned char tmp_v;
    vector unsigned char tmp2_v;
    byte_vec counts;
    counts.u8 = misc->fm_neg128_v; // set to -128, offset to allow each 8-bit int to hold up to 255.
                                        // so effectively, can guarantee holding 128*16 = 2048.
                                        // Since I count from left or right, whichever is closer, this means
                                        // we can support an occ_b interval of up to 4096 with guarantee of
                                        // correctness.
    if (meta->alph_type == fm_DNA) {

      if (!up_b) { // count forward, adding
        for (i=1+floor(landmark/4.0) ; i+15<( (pos+1)/4);  i+=16) { // keep running until i begins a run that shouldn't all be counted
             BWT_v    = *(vector unsigned char*)(BWT+i);
             FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
             FM_COUNT_2BIT(tmp_v, tmp2_v, counts.i8);
        }

        int remaining_cnt = pos + 1 -  i*4 ;
        if (remaining_cnt > 0) {
             BWT_v    = *(vector unsigned char*)(BWT+i);
             FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
             tmp_v    = vec_and(tmp_v, *(misc->fm_masks_v + remaining_cnt)); // leaves only the remaining_cnt chars in the array
             FM_COUNT_2BIT(tmp_v, tmp2_v, counts.i8);
        }

      } else { // count backwards, subtracting
        for (i=(landmark/4)-15 ; i>(pos/4);  i-=16) {
             BWT_v = *(vector unsigned char*)(BWT+i);
             FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
             FM_COUNT_2BIT(tmp_v, tmp2_v, counts.i8);
        }

        int remaining_cnt = 64 - (pos + 1 - i*4);
        if (remaining_cnt > 0) {
             BWT_v = *(vector unsigned char*)(BWT+i);
             FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
             tmp_v    = vec_and(tmp_v, *(misc->fm_reverse_masks_v + remaining_cnt)); // leaves only the remaining_cnt chars in the array
             FM_COUNT_2BIT(tmp_v, tmp2_v, counts.i8);
        }
      }

    } else if ( meta->alph_type == fm_DNA_full) {

      if (!up_b) { // count forward, adding

        for (i=1+floor(landmark/2.0) ; i+15<( (pos+1)/2);  i+=16) { // keep running until i begins a run that shouldn't all be counted
             BWT_v    = *(vector unsigned char*)(BWT+i);
             FM_MATCH_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
             FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts.i8);
        }
        int remaining_cnt = pos + 1 -  i*2 ;
        if (remaining_cnt > 0) {
             BWT_v    = *(vector unsigned char*)(BWT+i);
             FM_MATCH_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
             tmp_v     = vec_and(tmp_v, *(misc->fm_masks_v + (remaining_cnt+1)/2)); // mask characters we don't want to count
             tmp2_v    = vec_and(tmp2_v, *(misc->fm_masks_v + remaining_cnt/2));
             FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts.i8);
        }

      } else { // count backwards, subtracting
        for (i=(landmark/2)-15 ; i>(pos/2);  i-=16) {
             BWT_v = *(vector unsigned char*)(BWT+i);
             FM_MATCH_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
             FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts.i8);
        }

        int remaining_cnt = 32 - (pos + 1 - i*2);
        if (remaining_cnt > 0) {
             BWT_v = *(vector unsigned char*)(BWT+i);
             FM_MATCH_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
             tmp_v     = vec_and(tmp_v, *(misc->fm_reverse_masks_v + remaining_cnt/2)); // mask characters we don't want to count
             tmp2_v    = vec_and(tmp2_v, *(misc->fm_reverse_masks_v + (remaining_cnt+1)/2));
             FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts.i8);
        }
      }
    } else {
      esl_fatal("Invalid alphabet type\n");
    }

    FM_GATHER_8BIT_COUNTS(counts.i8,counts.i32);
    cnt  +=   ( up_b == 1 ?  -1 : 1) * ( counts.ints[3] + 128*16 );

  }

  if (c==0 && pos >= fm->term_loc) { // I overcounted 'A' by one, because '$' was replaced with an 'A'
    cnt--;
  }



  return cnt ;

}


/* Function:  fm_getOccCountLT()
 * Synopsis:  Compute number of occurrences of characters with value <c in BWT[1..pos]
 *
 * Purpose:   Scan through the BWT to compute number of occurrences of characters with
 *            value <c in BWT[0..pos], using altivec instructions to scan 16 bytes-at-a-time.
 *
 *            First, use checkpointed occurrence counts in the arrays occCnts_sb and occCnts_b.
 *            The checkpoint used is the one closest to pos, possibly requiring that counts be
 *            subtracted from the checkpointed total
 *
 *            The counting method is SIMD, loading 16 bytes (32 or 64 chars, depending on
 *            alphabet) at a time into the vector co-processors, then counting occurrences. One
 *            constraint of this approach is that occCnts_b checkpoints must be spaced at least
 *            every 32 or 64 chars (16 bytes, in pressed format), and in multiples of 64/32, so
 *            that _mm_load_si128 calls appropriately meet 16-byte-alignment requirements. That's
 *            a reasonable expectation, as spacings of 256 or more seem to give the best speed,
 *            and certainly better space-utilization.
 */

int
fm_getOccCountLT (FM_DATA *fm, FM_MISC_VARS *misc, int pos, uint8_t c, uint32_t *cnteq, uint32_t *cntlt) {

  if (c == 0 && pos >= fm->term_loc)// < 'A'?  cntlt depends on relationship of pos and the position where the '$' was replaced by 'A'
    *cntlt = 1;
  else
    *cntlt = 0;

  FM_METADATA *meta = misc->meta;

  int i,j;

  int cnt, cnt2;
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



  if ( landmark < fm->N || landmark == -1 ) {

    const uint8_t * BWT = fm->BWT;


    vector unsigned char c_v; // = *(fm_chars_v + c);
    vector unsigned char BWT_v;
    vector unsigned char tmp_v;
    vector unsigned char tmp2_v;
    byte_vec counts_lt;
    counts_lt.u8 = misc->fm_neg128_v; // set to -128, offset to allow each 8-bit int to hold up to 255.
    byte_vec counts_eq;
    counts_eq.u8 = misc->fm_neg128_v; // set to -128, offset to allow each 8-bit int to hold up to 255.
                                        // so effectively, can guarantee holding 128*16 = 2048.
                                        // Since I count from left or right, whichever is closer, this means
                                        // we can support an occ_b interval of up to 4096 with guarantee of
                                        // correctness.
    if (meta->alph_type == fm_DNA) {

      if (!up_b) { // count forward, adding
        for (i=1+floor(landmark/4.0) ; i+15<( (pos+1)/4);  i+=16) { // keep running until i begins a run that shouldn't all be counted
           BWT_v    = *(vector unsigned char*)(BWT+i);
           for (j=0; j<c; j++) {
             c_v = *(misc->fm_chars_v + j);
             FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
             FM_COUNT_2BIT(tmp_v, tmp2_v, counts_lt.i8);
           }
           c_v = *(misc->fm_chars_v + c);
           FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
           FM_COUNT_2BIT(tmp_v, tmp2_v, counts_eq.i8);
        }

        int remaining_cnt = pos + 1 -  i*4 ;
        if (remaining_cnt > 0) {
          BWT_v    = *(vector unsigned char*)(BWT+i);
          for (j=0; j<c; j++) {
             c_v = *(misc->fm_chars_v + j);
             FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
             tmp_v    = vec_and(tmp_v, *(misc->fm_masks_v + remaining_cnt)); // leaves only the remaining_cnt chars in the array
             FM_COUNT_2BIT(tmp_v, tmp2_v, counts_lt.i8);
          }
          c_v = *(misc->fm_chars_v + c);
          FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
          tmp_v    = vec_and(tmp_v, *(misc->fm_masks_v + remaining_cnt)); // leaves only the remaining_cnt chars in the array
          FM_COUNT_2BIT(tmp_v, tmp2_v, counts_eq.i8);
        }

      } else { // count backwards, subtracting
        for (i=(landmark/4)-15 ; i>(pos/4);  i-=16) {
          BWT_v = *(vector unsigned char*)(BWT+i);
          for (j=0; j<c; j++) {
            c_v = *(misc->fm_chars_v + j);
            FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
            FM_COUNT_2BIT(tmp_v, tmp2_v, counts_lt.i8);
          }
          c_v = *(misc->fm_chars_v + c);
          FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
          FM_COUNT_2BIT(tmp_v, tmp2_v, counts_eq.i8);
        }

        int remaining_cnt = 64 - (pos + 1 - i*4);
        if (remaining_cnt > 0) {
          BWT_v = *(vector unsigned char*)(BWT+i);
          for (j=0; j<c; j++) {
            c_v = *(misc->fm_chars_v + j);
            FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
            tmp_v    = vec_and(tmp_v, *(misc->fm_reverse_masks_v + remaining_cnt)); // leaves only the remaining_cnt chars in the array
            FM_COUNT_2BIT(tmp_v, tmp2_v, counts_lt.i8);
          }
          c_v = *(fm_chars_v + c);
          FM_MATCH_2BIT(BWT_v, c_v, tmp_v, tmp2_v, tmp_v);
          tmp_v    = vec_and(tmp_v, *(misc->fm_reverse_masks_v + remaining_cnt)); // leaves only the remaining_cnt chars in the array
          FM_COUNT_2BIT(tmp_v, tmp2_v, counts_eq.i8);
        }
      }

    } else if ( meta->alph_type == fm_DNA_full) {

      if (!up_b) { // count forward, adding

        for (i=1+floor(landmark/2.0) ; i+15<( (pos+1)/2);  i+=16) { // keep running until i begins a run that shouldn't all be counted
          BWT_v    = *(vector unsigned char*)(BWT+i);
          FM_LT_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
          FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts_lt.i8);
          FM_MATCH_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
          FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts_eq.i8);
        }
        int remaining_cnt = pos + 1 -  i*2 ;
        if (remaining_cnt > 0) {
          BWT_v    = *(vector unsigned char*)(BWT+i);
          FM_LT_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
          tmp_v     = vec_and(tmp_v, *(misc->fm_masks_v + (remaining_cnt+1)/2)); // mask characters we don't want to count
          tmp2_v    = vec_and(tmp2_v, *(misc->fm_masks_v + remaining_cnt/2));
          FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts_lt.i8);

          FM_MATCH_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
          tmp_v     = vec_and(tmp_v, *(misc->fm_masks_v + (remaining_cnt+1)/2)); // mask characters we don't want to count
          tmp2_v    = vec_and(tmp2_v, *(misc->fm_masks_v + remaining_cnt/2));
          FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts_eq.i8);
        }

      } else { // count backwards, subtracting
        for (i=(landmark/2)-15 ; i>(pos/2);  i-=16) {
          BWT_v = *(vector unsigned char*)(BWT+i);
          FM_LT_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
          FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts_lt.i8);
          FM_MATCH_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
          FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts_eq.i8);
        }

        int remaining_cnt = 32 - (pos + 1 - i*2);
        if (remaining_cnt > 0) {
          BWT_v = *(vector unsigned char*)(BWT+i);
          FM_LT_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
          tmp_v     = vec_and(tmp_v, *(misc->fm_reverse_masks_v + remaining_cnt/2)); // mask characters we don't want to count
          tmp2_v    = vec_and(tmp2_v, *(fmisc->m_reverse_masks_v + (remaining_cnt+1)/2));
          FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts_lt.i8);

          FM_MATCH_4BIT(BWT_v, c_v, tmp_v, tmp2_v);
          tmp_v     = vec_and(tmp_v, *(misc->fm_reverse_masks_v + remaining_cnt/2)); // mask characters we don't want to count
          tmp2_v    = vec_and(tmp2_v, *(misc->fm_reverse_masks_v + (remaining_cnt+1)/2));
          FM_COUNT_4BIT( (vector signed char)tmp_v,  (vector signed char)tmp2_v, counts_eq.i8);
        }
      }
    } else {
      esl_fatal("Invalid alphabet type\n");
    }

    FM_GATHER_8BIT_COUNTS(counts_lt.i8,counts_lt.i32);
        (*cntlt)  +=   ( up_b == 1 ?  -1 : 1) * ( counts_lt.ints[3] + 128*16 );

    FM_GATHER_8BIT_COUNTS(counts_eq.i8,counts_eq.i32);
        (*cnteq)  +=   ( up_b == 1 ?  -1 : 1) * ( counts_eq.ints[3] + 128*16 );

  }

  if (c==0 && pos >= fm->term_loc) { // I overcounted 'A' by one, because '$' was replaced with an 'A'
    (*cnteq)--;
  }



  return cnt ;

}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
