#include "divsufsort.h"
#include "fm.h"


//static const char *optString = "a:b:f:deupsrnl:h?";
static const char *optString = "a:b:s:del:h?";
static const struct option longOpts[] = {
  { "alph",       required_argument, NULL, 'a' },
  { "bin_length", required_argument, NULL, 'b' },
  { "sa_freq",    required_argument, NULL, 's' },
  { "help",       no_argument,       NULL, 'h' },
  { NULL,         no_argument,       NULL,   0 }
};

void
display_usage () {

	fprintf (stderr,
			"Usage: fmbuild [optional args] infile outfile \n\n"
			"--alph (dna|dna_full|amino|alpha)\n"
			"   Select alphabet for FMindex:\n"
			"   * dna:      ACGT (default)\n"
			"   * dna_full: ACGTRYMKSWHBVDN\n"
			"   * amino:    ACDEFGHIKLMNPQRSTVWYBJZOUX\n"
			"\n"
			"--bin_length len   (-b) \n"
			"   Length of bins for occurence counting; shorter bins\n"
			"   typically lead to shorter run times. Must be power\n"
			"   of 2, between 128 and 4096 (default 512)\n"
			"\n"
			"--sa_freq len  (-s)\n"
			"   Number of characters between sampled suffix array values\n"
			"   Must be power of 2, typically 8 to 32 (default 32)\n"
			"\n"
			"--help\n"
			"   Produce this message\n"
			"\n");
	exit(0);
}

/* Function:  main()
 * Synopsis:  break input sequence set into chunks, for each one building the
 *            Burrows-Wheeler transform and corresponding FM-index. Maintain requisite
 *            meta data.
 * Notes:     Currently depends on the divsufsort-lite code of Yuta Mori, though this
 *            could easily be replaced.
 */
int
main(int argc, char *argv[]) {
	FILE *fp             = NULL;
	uint8_t *T           = NULL;
	uint8_t *BWT         = NULL;
	int *SA              = NULL; //what I write will be 32-bit ints, but I need to keep this as int so it'll work with libdivsufsort
	uint32_t *occCnts_sb = NULL; // same indexing as above
	uint32_t *cnts_sb    = NULL;
	uint16_t *occCnts_b  = NULL; // this is logically a 2D array, but will be indexed as occ_cnts[alph_size*index + char]  (instead of occ_cnts[index][char])
	uint16_t *cnts_b     = NULL;
	FM_METADATA *meta    = NULL;
	char *inv_alph       = NULL;
	char *alph           = NULL;



	long i,j, c, k;
	int status;
	int result;
	int opt = 0;
	int longIndex = 0;

	int num_freq_cnts_sb ;
	int num_freq_cnts_b ;
	int num_SA_samples ;

	FM_ALLOC (meta, sizeof(FM_METADATA));
	meta->alph_type   = fm_DNA;
	meta->freq_SA     = 8;
	meta->freq_cnt_b  = 512;
	meta->freq_cnt_sb = pow(2,16); //65536 - that's the # values in a short


	opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
	while( opt != -1 ) {
	  switch( opt ) {
		  case 'a':
			  if (strcmp(optarg, "dna")==0)
				  meta->alph_type = fm_DNA;
			  else if (strcmp(optarg, "dna_full")==0)
				  meta->alph_type = fm_DNA_full;
			  else if (strcmp(optarg, "amino")==0)
				  meta->alph_type = fm_AMINO;
			  else
				  FM_FAIL("Unknown alphabet type. Try 'dna', 'dna_full', or 'amino'\n%s", "");
			  break;
		  case 'b':
			  meta->freq_cnt_b = atoi(optarg);
			  if ( meta->freq_cnt_b < 32 || meta->freq_cnt_b >4096 ||  (meta->freq_cnt_b & (meta->freq_cnt_b - 1)) /* test power of 2*/ ) {
				  fprintf (stderr, "bin_length must be a power of 2, at least 128, and at most 4096\n");
				  exit(1);
			  }
			  break;
		  case 's':
			  meta->freq_SA = atoi(optarg);
			  if ( (meta->freq_SA & (meta->freq_SA - 1)) /* test power of 2*/ ) {
				  fprintf (stderr, "SA_freq must be a power of 2\n");
				  exit(1);
			  }
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
	  fprintf (stderr, "must specify infile and outfile; for more details, run with -h flag\n");
	  exit(1);
	}


	const char *fname_in = argv[optind];
	const char *fname_out = argv[optind+1];


	/* Open a file for reading. */
	if((fp = fopen(fname_in, "rb")) == NULL) {
	FM_FAIL("%s: Cannot open file `%s': ", argv[0], fname_in);
	}

	/* Get the file size. */
	//TODO: this is where the loop should start, just use a sequence block
	if(fseek(fp, 0, SEEK_END) == 0) {
	  meta->N = ftell(fp);
	rewind(fp);
	if(meta->N < 0)
	  FM_FAIL("%s: Cannot ftell `%s': ", argv[0], fname_in);
	} else {
	FM_FAIL( "%s: Cannot fseek `%s': ", argv[0], fname_in);
	}

	meta->SA_shift = log2(meta->freq_SA);
	meta->cnt_shift_b = log2(meta->freq_cnt_b);

	num_freq_cnts_b  = 1+ceil((float)meta->N/meta->freq_cnt_b);
	num_SA_samples   = 1+floor((float)meta->N/meta->freq_SA);

	meta->cnt_shift_sb = log2(meta->freq_cnt_sb);
	num_freq_cnts_sb = 1+ceil((float)meta->N/meta->freq_cnt_sb);


	//getInverseAlphabet
	createAlphabet(meta->alph_type, &alph, &inv_alph, &(meta->alph_size), &(meta->charBits));



	/* Allocate n bytes for BWT and 4n bytes for SA*/
	FM_ALLOC (T, meta->N * sizeof(uint8_t));
	FM_ALLOC (BWT, meta->N * sizeof(uint8_t));
	FM_ALLOC (SA, meta->N * sizeof(int));
	FM_ALLOC (occCnts_sb, num_freq_cnts_sb *  meta->alph_size * sizeof(uint32_t)); // every freq_cnt_sb positions, store an array of ints
	FM_ALLOC (cnts_sb,    meta->alph_size * sizeof(uint32_t));
	FM_ALLOC (occCnts_b,  num_freq_cnts_b *  meta->alph_size * sizeof(uint16_t)); // every freq_cnt_b positions, store an array of 8-byte ints
	FM_ALLOC (cnts_b,     meta->alph_size * sizeof(uint16_t));

	if((T == NULL)  || (SA==NULL) || (BWT==NULL) || (cnts_b==NULL) || (occCnts_b==NULL) || (cnts_sb==NULL) || (occCnts_sb==NULL) ) {
	FM_FAIL( "%s: Cannot allocate memory.\n", argv[0]);
	}


	/* Read n bytes of data. */
	if(fread(T, sizeof(unsigned char), (size_t)meta->N, fp) != (size_t)meta->N) {
	FM_FAIL("%s: %s '%s'", argv[0],
			((ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in"),
			argv[1]);
	}
	fclose(fp);




	/* Reverse the sequence in place, for suitable BWT traversal.
	*  Convert input sequence to alphabet of 1..k for alphabet of size k.
	*  (a) collapsing upper/lower case for appropriate sorting.
	*  (b) reserving 0 for '$', which must be lexicographically smallest
	*/
	T[meta->N-1] = 0; // by default, this is LF;  convert it to 0, essentially '$'

	for(i=0; i < (meta->N-1)/2; ++i) {
	  const uint8_t tmpc1 = inv_alph[T[meta->N-2-i]];
	  const uint8_t tmpc2 = inv_alph[T[i]];

	  if (tmpc1 == 255 || tmpc2 == 255)
		  fprintf (stderr, "requested alphabet doesn't match input test: T[%ld] = (%d), T[%ld] = (%d)\n",
				  meta->N-2-i, T[meta->N-2-i], i, T[i]);

	  T[i] = tmpc1;
	  T[meta->N-2-i] = tmpc2;

	}
	if ( ! (meta->N & 0x1) )  //odd # chars, need to fix the middle one
	  T[meta->N/2 - 1 ] = inv_alph[T[meta->N/2 - 1 ]];

	/* Construct the Suffix Array */
	result = divsufsort(T, SA, meta->N);
	if ( result < 0 )
	FM_FAIL("%s: Error building BWT.\n", argv[0]);


	/* Construct the BWT, SA landmarks, and FM-index */
	for (c=0; c<meta->alph_size; c++) {
	  cnts_sb[c] = 0;
	  cnts_b[c] = 0;
	  FM_OCC_CNT(sb, 0, c ) = 0;
	  FM_OCC_CNT(b, 0, c ) = 0;
	}

	BWT[0] =  SA[0]==0 ? 0 /* '$' */ : T[ SA[0]-1] ;
	cnts_sb[BWT[0]]++;
	cnts_b[BWT[0]]++;


	//Scan through BWT to build the FM index structures
	for(j=1; j < meta->N; ++j) {
	  BWT[j] =  SA[j]==0 ? 0 /* '$' */ : T[ SA[j]-1] ;

	  if ( !(j % meta->freq_SA) )
		  SA[ j>>meta->SA_shift ] = ( SA[j] == meta->N - 1 ? -1 : SA[j] ) ; // handle the wrap-around '$'

	  cnts_sb[BWT[j]]++;
	  cnts_b[BWT[j]]++;

	  const long joffset = j+1;
	  if ( !(  joffset % meta->freq_cnt_b) ) {  // (j+1)%freq_cnt_b==0  , i.e. every freq_cnt_bth position, noting that it's a zero-based count

		  for (c=0; c<meta->alph_size; c++)
			  FM_OCC_CNT(b, (joffset>>meta->cnt_shift_b), c ) = cnts_b[c];

		  if ( !(joffset % meta->freq_cnt_sb) ) {  // j%freq_cnt_sb==0
			  for (c=0; c<meta->alph_size; c++) {
				  FM_OCC_CNT(sb, (joffset>>meta->cnt_shift_sb), c ) = cnts_sb[c];
				  cnts_b[c] = 0;
			  }
		  }

	  }
	}

	//wrap up the counting;
	for (c=0; c<meta->alph_size; c++) {
	  FM_OCC_CNT(b, num_freq_cnts_b-1, c ) = cnts_b[c];
	  FM_OCC_CNT(sb, num_freq_cnts_sb-1, c ) = cnts_sb[c];
	}

	  /* Convert BWT to packed versions if appropriate. */
	if (meta->alph_type == fm_DNA) {
		 //4 chars per byte.  Counting will be done based on quadruples 0..3; 4..7; 8..11; etc.
		  for(i=0; i < meta->N-3; i+=4)
			  BWT[i>>2] = BWT[i]<<6 | BWT[i+1]<<4 | BWT[i+2]<<2 | BWT[i+3];
		  if (i>=meta->N-3)  BWT[i>>2] =  BWT[i]<<6;
		  if (i>=meta->N-2)  BWT[i>>2] =  BWT[i+1]<<4;
		  if (i==meta->N-1)  BWT[i>>2] =  BWT[i+2]<<2;
	} else if (meta->alph_type == fm_DNA_full) {
		//2 chars per byte.  Counting will be done based on quadruples 0..3; 4..7; 8..11; etc.
		  for(i=0; i < meta->N-1; i+=2)
			  BWT[i>>1] = BWT[i]<<4 | BWT[i+1];
		  if (i==meta->N-1)
			  BWT[i>>1] =  BWT[i]<<4 ;
	}


	/* Open a file for writing. */
	if((fp = fopen(fname_out, "wb")) == NULL)
	FM_FAIL( "%s: Cannot open file `%s': ", argv[0], fname_out);
	/* Write the FM-index meta data */
	if(fwrite(meta, sizeof(FM_METADATA), 1, fp) != 1)
	  FM_FAIL( "%s: Error writing BWT.\n", argv[0]);
	if(fwrite(BWT, sizeof(uint8_t), (size_t)((1+meta->N)/(8/meta->charBits)), fp) != (size_t)((1+meta->N)/(8/meta->charBits)))
	  FM_FAIL( "%s: Error writing FM index.\n", argv[0]);
	if(fwrite(SA, sizeof(int), (size_t)num_SA_samples, fp) != (size_t)num_SA_samples)
	  FM_FAIL( "%s: Error writing FM index.\n", argv[0]);
	if(fwrite(occCnts_b, sizeof(uint16_t)*(meta->alph_size), (size_t)num_freq_cnts_b, fp) != (size_t)num_freq_cnts_b)
	  FM_FAIL( "%s: Error writing FM index.\n", argv[0]);
	if(fwrite(occCnts_sb, sizeof(uint32_t)*(meta->alph_size), (size_t)num_freq_cnts_sb, fp) != (size_t)num_freq_cnts_sb)
	  FM_FAIL( "%s: Error writing FM index.\n", argv[0]);

	fclose(fp);
	free(T);
	free(BWT);
	free(SA);
	free(occCnts_b);
	free(cnts_b);
	free(occCnts_sb);
	free(cnts_sb);
	free(meta);
	free(inv_alph);
	free(alph);



	return (fmOK);


ERROR:
	/* Deallocate memory. */
	if (fp)         fclose(fp);
	if (T)          free(T);
	if (BWT)        free(BWT);
	if (SA)         free(SA);
	if (occCnts_b)  free(occCnts_b);
	if (cnts_b)     free(cnts_b);
	if (occCnts_sb) free(occCnts_sb);
	if (cnts_sb)    free(cnts_sb);
	if (meta)       free(meta);
	if (inv_alph)   free(inv_alph);
	if (alph)       free(alph);


	fprintf (stderr, "failure during memory allocation\n");

  	exit(status);

}
