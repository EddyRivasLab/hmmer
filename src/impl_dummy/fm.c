#include "hmmer.h"
#include "impl_dummy.h"



/* Function:  fm_initGlobals()
 * do nothing.  This is a placeholder in the dummy implementation for
 * a function that does something in _sse and _vmx code
 */
int
fm_initGlobals( FM_METADATA *meta  ) {
	return eslOK;
}

/* Function:  fm_destroyGlobals()
 * do nothing.  This is a placeholder in the dummy implementation for
 * a function that does something in _sse and _vmx code
 */
int
fm_destroyGlobals( ) {
    return eslOK;
}


/* Function:  bwt_getOccCount()
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
/* Function:  bwt_getOccCount_packed_dummy()
 * Synopsis:  Scan through the BWT to fully compute counts from landmark to pos
 * Incept:    TJW, Sat Dec 25 20:35:33 MST 2010 [Tucson]
 * Purpose:   Scan through the BWT to fully compute counts from landmark to pos. In this case,
 *            the BWT is assumed to be packed two chars per byte.
 */
int
fm_getOccCount (FM_METADATA *meta, FM_DATA *fm, int pos, uint8_t c) {

	int i;

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


	if ( landmark < fm->N ) {

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



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
