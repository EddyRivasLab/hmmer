#include "fm.h"
#include <sys/times.h>

static const char *optString = "o:ch?";
static const struct option longOpts[] = {
  { "out",        required_argument, NULL, 'o' },
  { "count_only", no_argument,       NULL, 'c' },

  { "help",       no_argument,       NULL, 'h' },
  { NULL,         no_argument,       NULL,   0 }
};

int occCallCnt;
//int interesting = 14853708;
int interesting = -1;



/* Global variables - initialized once in initGlobals(), to avoid recomputing them
 *
 * Purists, avert your eyes.
 */
__m128i *masks_v;
__m128i *reverse_masks_v;
__m128i *chars_v;
__m128i zeros_v;
__m128i m0f;  //00 00 11 11
__m128i m01;  //01 01 01 01
__m128i m11;  //00 00 00 11
__m128i allones_v;

int num_freq_cnts_b ;
int num_freq_cnts_sb;

int *Cvals;
int maskSA;
int shiftSA;


/* Function:  initGlobals()
 * Synopsis:  Run set of queries against an FM
 * Purpose:   Initialize a bunch of global variables, mostly vector masks, that I
 * don't want to recompute each time bwt_getOccCount() is called.
 */
int
initGlobals( FM_METADATA *meta, FM_DATA *fm  ) {
	int status;
	int i;
	int chars_per_vector;

	allones_v = _mm_set1_epi8(0xff);

	occCallCnt = 0;
	//capture abs() for C array (so I don't have to do it for every query sequence)
	Cvals = NULL;
	ESL_ALLOC(Cvals, meta->alph_size * sizeof(int) );
	for(i=0; i<meta->alph_size; i++)
		Cvals[i] = abs(fm->C[i]);

	zeros_v = _mm_set1_epi8((int8_t) 0x00);      //00 00 00 00
	m0f     = _mm_set1_epi8((int8_t) 0x0f);      //00 00 11 11

	if (meta->alph_type = fm_DNA) {
		m01 = _mm_set1_epi8((int8_t) 0x55);	 //01 01 01 01
		m11 = _mm_set1_epi8((int8_t) 0x03);  //00 00 00 11
	}

    //set up an array of vectors, one for each character in the alphabet
	chars_v         = NULL;
	ESL_ALLOC(chars_v, meta->alph_size * sizeof(__m128));
	for (i=0; i<meta->alph_size; i++) {
		int8_t c = i;
		if (meta->alph_type == fm_DNA) {
			// need 4 copies per byte
			c &= i<<2;
			c &= i<<4;
			c &= i<<6;
		} //else, just leave it on the right-most bits

		chars_v[i] = _mm_set1_epi8(c);
	}

	printf ("D\n");

	/* this is a collection of masks used to clear off the left- or right- part
	 *  of a register when we shouldn't be counting the whole thing
	 * Incrementally chew off the 1s in chunks of 2 (for DNA) or 4 (for DNA_full)
	 * from the right side, and stick each result into an element of a __m128 array
	 */
	chars_per_vector = 128/meta->charBits;
	masks_v         = NULL;
	reverse_masks_v = NULL;
	ESL_ALLOC(masks_v, chars_per_vector*sizeof(__m128));
	ESL_ALLOC(reverse_masks_v, chars_per_vector*sizeof(__m128));
	{
		byte_m128 arr;
		arr.m128 = allones_v;
		for (i=chars_per_vector-1; i>0; i--) { // don't need the 0th entry
			int byte_mask=0xff; //11 11 11 11
			if (meta->alph_type == fm_DNA) {
				switch (i&0x11) {
				  case 1:
					byte_mask ^= 0x3f; //11 00 00 00
					break;
				  case 2:
					byte_mask ^= 0x0f; //11 11 00 00
					break;
				  case 3:
					byte_mask ^= 0x03; //11 11 11 00
					break;
				  default:
			        break;
				}
			} else if (meta->alph_type == fm_DNA_full) {
				if (i&0x1)
					byte_mask ^= 0x0f; //11 11 00 00
			}
			  arr.bytes[i/meta->charBits] = byte_mask; //if it's odd, chew off the last 4 bits, else clear the whole byte
			  masks_v[i]             = *(__m128i*)(&(arr.m128));
			  reverse_masks_v[i]  = _mm_andnot_si128(masks_v[i], allones_v );
		}
	}
	printf ("E\n");
	maskSA       =  meta->freq_SA - 1;
	shiftSA      =  meta->SA_shift;

	printf ("F\n");

	return eslOK;

ERROR:
	if (masks_v)         free (masks_v);
	if (reverse_masks_v) free (reverse_masks_v);

	esl_fatal("Error allocating memory in initGlobals\n");
	//return fmFAIL;
}



/* Function:  bwt_getOccCount()
 * Synopsis:  Compute number of occurrences of c in BWT[1..pos]
 *
 * Purpose:   Scan through the BWT to compute number of occurrence of c in BWT[1..pos],
 *            using SSE to scan 16 bytes-at-a-time.
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
#ifndef FMDEBUG
inline
#endif
int
bwt_getOccCount (FM_METADATA *meta, FM_DATA *fm, int pos, uint8_t c) {

	int i;


	int cnt;
	const int b_pos          = (pos+1) >> meta->cnt_shift_b; //floor(pos/b_size)   : the b count element preceding pos
	const uint16_t * restrict occCnts_b  = fm->occCnts_b;
	const uint32_t * restrict occCnts_sb = fm->occCnts_sb;
	const int sb_pos         = (pos+1) >> meta->cnt_shift_sb; //floor(pos/sb_size) : the sb count element preceding pos


	const int cnt_mod_mask_b = meta->freq_cnt_b - 1; //used to compute the mod function
	const int b_rel_pos      = (pos+1) & cnt_mod_mask_b; // pos % b_size      : how close is pos to the boundary corresponding to b_pos
	const int up_b           = b_rel_pos>>(meta->cnt_shift_b - 1); //1 if pos is expected to be closer to the boundary of b_pos+1, 0 otherwise
	const int landmark       = ((b_pos+up_b)<<(meta->cnt_shift_b)) - 1 ;
	// get the cnt stored at the nearest checkpoint
	if ( up_b  &&  (b_pos+1) ==  (sb_pos+1) * (1<<(meta->cnt_shift_sb - meta->cnt_shift_b)) )
		cnt = FM_OCC_CNT(sb, sb_pos+1, c )  ;// b_pos+1 is referencing the same position as sb+1. I don't store those b-vectors, so just go straight to the sb vector
	else {
		cnt =  FM_OCC_CNT(sb, sb_pos, c );
		cnt += FM_OCC_CNT(b, b_pos + up_b, c ) ;
	}

	const uint8_t * BWT = fm->BWT;


	register __m128i c_v = *(chars_v + c);
	register __m128i BWT_v;
	register __m128i tmp_v;
	register __m128i tmp2_v;
	register __m128i counts_v = allones_v; // set to -128, offset to allow each 8-bit int to hold up to 255.
	                                       // so effectively, can guarantee holding 128*16 = 2048.
	                                       // Since I count from left or right, whichever is closer, this means
	                                       // we can support an occ_b interval of up to 4096 with guarantee of
	                                       // correctness.

	if (meta->alph_type == fm_DNA) {


		if (!up_b) { // count forward, adding
			for (i=1+floor(landmark/4.0) ; i+15<( (pos+1)/4);  i+=16) { // keep running until i begins a run that shouldn't all be counted
				BWT_v    = *(__m128i*)(BWT+i);
				FM_COUNTS_SSE_4PACKED(BWT_v, c_v, tmp_v, tmp2_v, counts_v);
			}
			int remaining_cnt = pos + 1 -  i*4 ;
			if (remaining_cnt > 0) {
				BWT_v    = *(__m128i*)(BWT+i);
				BWT_v    = _mm_and_si128(BWT_v, *(masks_v + remaining_cnt)); // leaves only the remaining_cnt chars in the array
				FM_COUNTS_SSE_4PACKED(BWT_v, c_v, tmp_v, tmp2_v, counts_v);
			}

		} else { // count backwards, subtracting
			for (i=(landmark/4)-15 ; i>(pos/4);  i-=16) {
			 	BWT_v = *(__m128i*)(BWT+i);
			 	FM_COUNTS_SSE_4PACKED(BWT_v, c_v, tmp_v, tmp2_v, counts_v);
			}

			int remaining_cnt = 64 - (pos + 1 - i*4);
			if (remaining_cnt > 0) {
				BWT_v = *(__m128i*)(BWT+i);
				BWT_v    = _mm_and_si128(BWT_v, *(reverse_masks_v + remaining_cnt)); // leaves only the remaining_cnt chars in the array
				FM_COUNTS_SSE_4PACKED(BWT_v, c_v, tmp_v, tmp2_v, counts_v);
			}
		}

	} else if ( meta->alph_type == fm_DNA_full) {

		if (!up_b) { // count forward, adding
			for (i=1+floor(landmark/2.0) ; i+15<( (pos+1)/2);  i+=16) { // keep running until i begins a run that shouldn't all be counted
				BWT_v    = *(__m128i*)(BWT+i);
				FM_COUNTS_SSE_2PACKED(BWT_v, c_v, tmp_v, tmp2_v, counts_v);
			}
			int remaining_cnt = pos + 1 -  i*2 ;
			if (remaining_cnt > 0) {
				BWT_v    = *(__m128i*)(BWT+i);
				BWT_v    = _mm_and_si128(BWT_v, *(masks_v + remaining_cnt)); // leaves only the remaining_cnt chars in the array
				FM_COUNTS_SSE_2PACKED(BWT_v, c_v, tmp_v, tmp2_v, counts_v);
			}

		} else { // count backwards, subtracting
			for (i=(landmark/2)-15 ; i>(pos/2);  i-=16) {
			 	BWT_v = *(__m128i*)(BWT+i);
				FM_COUNTS_SSE_2PACKED(BWT_v, c_v, tmp_v, tmp2_v, counts_v);
			}

			int remaining_cnt = 32 - (pos + 1 - i*2);
			if (remaining_cnt > 0) {
				BWT_v = *(__m128i*)(BWT+i);
				BWT_v    = _mm_and_si128(BWT_v, *(reverse_masks_v + remaining_cnt)); // leaves only the remaining_cnt chars in the array
				FM_COUNTS_SSE_2PACKED(BWT_v, c_v, tmp_v, tmp2_v, counts_v);
			}
		}
	} else {
		esl_fatal("Invalid alphabet type\n");
	}

	counts_v = _mm_xor_si128(counts_v, allones_v); //counts are stored in signed bytes, base -128. Move them to unsigned bytes
    FM_GATHER_8BIT_COUNTS(counts_v,counts_v,counts_v);

	cnt  +=   ( up_b == 1 ?  -1 : 1) * ( _mm_extract_epi16(counts_v, 0) );
	occCallCnt++;
	return cnt ;

}




/* Function:  getFMRange()
 * Synopsis:  For a given query sequence, find its interval in the FM-index
 * Purpose:   Implement Algorithm 3.6 (p17) of Firth paper (A Comparison of BWT Approaches
 *            to String Pattern Matching). All the meat is in the method of counting
 *            characters - bwt_getOccCount, which depends on compilation choices.
 */
int
getFMRange( FM_METADATA *meta, FM_DATA *fm, char * restrict query, char * restrict inv_alph,
		FM_INTERVAL *interval) {

	int count1, count2;
	int C_base;
	int i=0;
	int32_t * restrict C = fm->C;

	uint8_t c = inv_alph[query[0]];
	interval->lower = abs(C[c]);
	interval->upper  = abs(C[c+1])-1;

	while (interval->lower>0 && interval->lower <= interval->upper) {
		c = query[++i];
		if (c == '\0')  // end of query - the current range defines the hits
			break;

		c = inv_alph[c];
		C_base = Cvals[c];

		//TODO: these will often overlap - might get acceleration by merging to a single redundancy-avoiding call
		count1 = bwt_getOccCount (meta, fm, interval->lower-1, c);
		count2 = bwt_getOccCount (meta, fm, interval->upper, c);

		interval->lower = C_base + count1;
		interval->upper = C_base + count2 - 1;
	}

	return eslOK;
}



/* Function:  getChar()
 * Synopsis:  Find the character c residing at a given position in the BWT.
 * Purpose:   The returned char is used by getFMHits(), to seed a call to
 *            bwt_getOccCount().
 */
#ifndef FMDEBUG
inline
#endif
uint8_t
getChar( FM_METADATA *meta, int j, const uint8_t * restrict B ) {
	uint8_t c = -1;

	if (meta->alph_type == fm_DNA) {
        /*
         *  B[j>>2] is the byte of B in which j is found (j/4)
         *
         *  without branching, do:
         *  if (j%0 == 0) use the rightmost two bits
         *  else shift  4-(j%1) bits to the right, and use the remaining rightmost 2 bits
         *
         *  the bizarre  ((0x4 - (j&0x3))&0x3)<<1 gets the amount of shift described above
         */
		c = (B[j>>2] >> ( ((0x4 - (j&0x3))&0x3)<<1 ) ) & 0x3;
	} else if (meta->alph_type == fm_DNA_full) {
		c = (B[j>>1] >> (((j&0x1)^0x1)<<2) ) & 0xf;  //unpack the char: shift 4 bits right if it's odd, then mask off left bits in any case
	} else {
		esl_fatal("Invalid alphabet type\n");
	}

	return c;
}

/* Function:  getFMHits()
 * Synopsis:  For a given interval, identify the position in original text for each element
 *            of interval
 * Incept:    TJW, Mon Jan  3 11:35:24 EST 2011 [Janelia]
 * Purpose:   Implement Algorithm 3.7 (p17) of Firth paper (A Comparison of BWT Approaches
 *            to String Pattern Matching). Most of the meat is in the method of counting
 *            characters - bwt_getOccCount, which depends on compilation choices.
 */
#ifndef FMDEBUG
inline
#endif
int
getFMHits( FM_METADATA *meta, FM_DATA *fm, FM_INTERVAL *interval, FM_HIT *hits_ptr) {

	const uint8_t * restrict B = fm->BWT;
	const uint32_t * restrict occCnts_sb = fm->occCnts_sb;
	const uint16_t * restrict occCnts_b  = fm->occCnts_b;
	int i, j, len = 0;

	for (i = interval->lower;  i<= interval->upper; i++) {
		j = i;
		len = 0;

		while ( j & maskSA ) { //go until we hit a position in the full SA that was sampled during FM index construction
			uint8_t c = getChar( meta, j, B);
			j = bwt_getOccCount (meta, fm, j-1, c);
			j += Cvals[c] ;
			len++;
		}
		const int tmp = j >> shiftSA;
		hits_ptr[i - interval->lower].start = fm->SA[ tmp ] + len;
	}

	return eslOK;

}


/* Function:  readFM()
 * Synopsis:  Read the FM index off disk
 * Incept:    TJW, Wed Dec 22 19:05:04 MST 2010 [Tucson]
 * Purpose:   Read the FM-index as written by fmbuild.
 *            First read the metadata header, then allocate space for the full index,
 *            then read it in.
 */
int
readFM( const char *fname, FM_METADATA *meta,  FM_DATA *fm )
{
	//shortcut variables
	int *C               = NULL;


	uint16_t *occCnts_b  = NULL;
	uint32_t *occCnts_sb = NULL;
	FILE *fp             = NULL;
	int status;
	int i;
	int num_SA_samples;
	int prevC;
	int cnt;
	int bwtSize;

	if((fp = fopen(fname, "rb")) == NULL)
		esl_fatal("Cannot open file `%s': ", fname);

	//get the FM meta data
	if(fread(meta, sizeof(FM_METADATA), 1, fp) != 1)
		esl_fatal("Error reading BWT size.%s\n", " ");


	num_freq_cnts_sb = 1+ceil((float)meta->N/meta->freq_cnt_sb);
	num_freq_cnts_b  = 1+ceil((float)meta->N/meta->freq_cnt_b);
	num_SA_samples   = 1+floor((float)meta->N/meta->freq_SA);


	// allocate and read the data
	bwtSize = meta->L  * sizeof(uint8_t);

	ESL_ALLOC (fm->T, bwtSize );
	ESL_ALLOC (fm->BWT, (1 + bwtSize/16) * 16  ); // 16 stuff for ensuring 16-byte aligned end of the array, which matters for SIMD stuff
	ESL_ALLOC (fm->SA, num_SA_samples * sizeof(uint32_t));
	ESL_ALLOC (fm->C, 1+meta->alph_size * sizeof(uint32_t));
	ESL_ALLOC (fm->occCnts_b,  num_freq_cnts_b *  (meta->alph_size ) * sizeof(uint16_t)); // every freq_cnt positions, store an array of ints
	ESL_ALLOC (fm->occCnts_sb,  num_freq_cnts_sb *  (meta->alph_size ) * sizeof(uint32_t)); // every freq_cnt positions, store an array of ints

	//shortcut variables
	C          = fm->C;
	occCnts_b  = fm->occCnts_b;
	occCnts_sb = fm->occCnts_sb;

	// read FM index structures
	if( fread(fm->T, 1, (size_t)bwtSize, fp) != (size_t)bwtSize)
		esl_fatal("%s: Error reading FM T.\n", "fmsearch");
	if( fread(fm->BWT, 1, (size_t)bwtSize, fp) != (size_t)bwtSize)
		esl_fatal("%s: Error reading FM BWT.\n", "fmsearch");
	if(fread(fm->SA, sizeof(int), (size_t)num_SA_samples, fp) != (size_t)num_SA_samples)
		esl_fatal("%s: Error reading FM SA.\n", "fmsearch");
	if(fread(occCnts_b, (meta->alph_size) * sizeof(uint16_t), (size_t)num_freq_cnts_b, fp) != (size_t)num_freq_cnts_b)
		esl_fatal("%s: Error reading FM Occ_b.\n", "fmsearch");
	if(fread(occCnts_sb, (meta->alph_size) * sizeof(uint32_t), (size_t)num_freq_cnts_sb, fp) != (size_t)num_freq_cnts_sb)
		esl_fatal("%s: Error reading FM Occ_sb.\n", "fmsearch");

	/*compute the first position of each letter in the alphabet in a sorted list
	* (with an extra value to simplify lookup of the last position for the last letter).
	* Negative values indicate that there are zero of that character in T, can be
	* used to establish the end of the prior range*/
	C[0] = 0;
	for (i=0; i<meta->alph_size; i++) {
	  prevC = abs(C[i]);

	  cnt = FM_OCC_CNT( sb, num_freq_cnts_sb-1, i);

	  if (cnt==0) {// none of this character
		  C[i+1] = prevC;
		  C[i] *= -1; // use negative to indicate that there's no character of this type, the number gives the end point of the previous
	  } else {
		  C[i+1] = prevC + cnt;
	  }
	}
	C[meta->alph_size] *= -1;



/*
	for (i=0; i<num_freq_cnts_b; i++) {
		  printf("b: %3d : ", i);
		  int c;
		  for (int c=0; c<meta->alph_size; c++) {
			  //FM_OCC_CNT_INCHAR( i, c)
			  printf(" %4d ",   FM_OCC_CNT_INCHAR(i,c) );
		  }
		  printf("\n");

	}

	for (int c=0; c<=meta->alph_size; c++) {
		printf ("C[%d] = %d\n" , c, C[c]);
	}
	exit(1);
*/


	return eslOK;

ERROR:
    if (fm->T)           free (fm->T);
    if (fm->BWT)         free (fm->BWT);
	if (fm->SA)          free (fm->SA);
	if (fm->C)           free (fm->C);
	if (fm->occCnts_b)   free (fm->occCnts_b);
	if (fm->occCnts_sb)  free (fm->occCnts_sb);
	esl_fatal("Error allocating memory in %s\n", "readFM");

}


/* Function:  display_usage()
 * Synopsis:  Print usage instructions
 * Incept:    TJW, Mon Jan 31 13:54:01 MST 2010
 */
void
display_usage () {

	fprintf (stderr,
			"Usage: fmsearch [optional args] file_fmindex file_queryseqs\n\n"
			"--out filename  (-o) \n"
			"   File to which fmsearch writes misses and the positions of\n"
			"   hits. Default is to not write them out (useful for timing\n"
			"   experiments).  If filename is 'stdout', then output\n"
			"   will go to stdout, instead of a file.\n"
			"\n"
			"--count_only  (-c)\n"
			"   Only compute the count of hits for a query, not the location\n"
			"   of each hit."
			"\n"
			);
	exit(0);
}


/* Function:  main()
 * Synopsis:  Run set of queries against an FM
 * Incept:    TJW, Fri Dec 24 21:30:51 MST 2010 [Tucson]
 * Purpose:   Read in a FM and a file of query sequences.
 *            For each query, find matching FM interval, then collect positions in
 *            the original text T for the corresponding occurences. These positions
 *            are 0-based (so first character is position 0).
 */
int
main(int argc,  char *argv[]) {

	void* tmp; // used for RALLOC calls
	clock_t t1, t2;
	struct tms ts1, ts2;
	//start timer
	t1 = times(&ts1);

	const char *fname_fm      = NULL;
	const char *fname_queries = NULL;
	char *inv_alph            = NULL;
	char *alph                = NULL;
	FM_HIT *hits              = NULL;
	char *line                = NULL;
	int status        = eslOK;
	int hit_cnt       = 0;
	int hit_indiv_cnt = 0;
	int miss_cnt      = 0;
	int hit_num       = 0;
	int hits_size     = 0;

	int count_only    = 0;

	FM_METADATA meta;
	FM_DATA fm;
	FM_INTERVAL interval;
	FILE* fp = NULL;
	FILE* out = NULL;

	int opt = 0;
	int longIndex = 0;

	opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
	while( opt != -1 ) {
	  switch( opt ) {
		  case 'o':
			  if (strcmp (optarg, "stdout") == 0)
				  out = stdout;
			  else
				  out = fopen(optarg,"w");
			  break;
		  case 'c':
			  count_only = 1;
			  break;
		  case 'h':   /* fall-through is intentional */
		  case '?':
			  display_usage();
			  break;
		  default:
			  fprintf(stderr, "Unknown flag: '%s'\n", longOpts[longIndex].name );
			  exit(1);
			  break;
	  }
	  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
	}


	if (argc - optind != 2) {
	  fprintf (stderr, "Must specify file_fmindex and file_queryseqs; "
			  "for more details, run with -h flag\n");
	  exit(1);
	}

	fname_fm = argv[optind];
	fname_queries = argv[optind+1];

	printf("A\n");


	readFM ( fname_fm, &meta, &fm );

	printf("B\n");
	initGlobals(&meta, &fm);

	printf("C\n");
		exit(1);

	createAlphabet(meta.alph_type, &alph, &inv_alph, &(meta.alph_size), NULL); // don't override charBits




	fp = fopen(fname_queries,"r");
	if (fp == NULL)
		esl_fatal("Unable to open file %s\n", fname_queries);

	ESL_ALLOC(line, FM_MAX_LINE * sizeof(char));

	hits_size = 200;
	ESL_RALLOC(hits, tmp, hits_size * sizeof(FM_HIT));

	while(fgets(line, FM_MAX_LINE, fp) ) {
		int qlen=0;
		while (line[qlen] != '\0' && line[qlen] != '\n')  qlen++;
		if (line[qlen] == '\n')  line[qlen] = '\0';

		getFMRange(&meta, &fm, line, inv_alph, &interval);

		if (interval.lower>0 && interval.lower <= interval.upper) {

			hit_num =  interval.upper - interval.lower + 1;

			if (hit_num > hits_size) {
				ESL_RALLOC(hits, tmp, hit_num * sizeof(FM_HIT));
				hits_size = hit_num;
			}

			if (!count_only)
				getFMHits(&meta, &fm, &interval, hits);

			hit_cnt++;
			hit_indiv_cnt += hit_num;


			//PRINTOUT ( "HIT:  %s  (%d, %d)\n", line, interval.lower, interval.upper);
			int i;
			for (i = 0; i< hit_num; i++) {
				int pos_fwd = meta.N - hits[i].start - qlen - 1; // flip the sequence reversal
				//PRINTOUT ( "\t%d\n", pos_fwd);
			}
		} else {
			//PRINTOUT ( "MISS: %s\n", line);
			miss_cnt++;
		}


	}

    free (hits);
    fclose(fp);

    // compute and print the elapsed time in millisec
    t2 = times(&ts2);
    double clk_ticks = sysconf(_SC_CLK_TCK);
    double elapsedTime = (t2-t1)/clk_ticks;
    //elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to s
    double throughput = occCallCnt/elapsedTime;

    fprintf (stderr, "hit: %-10d  (%d)\n", hit_cnt, hit_indiv_cnt);
    fprintf (stderr, "miss:%-10d\n", miss_cnt);
    fprintf (stderr, "run time:  %.2f seconds\n", elapsedTime);
    fprintf (stderr, "occ calls: %12s\n", commaprint(occCallCnt));
    fprintf (stderr, "occ/sec:   %12s\n", commaprint(throughput));


	exit(eslOK);


ERROR:
	printf ("failure allocating memory for hits\n");
	exit(status);


}





/*****************************************************************
 * @LICENSE@
 *****************************************************************/
