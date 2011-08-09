#include "divsufsort.h"
#include "hmmer.h"
#include "easel.h"

//#define PRINTBWT 1
//#define PRINTOCC 1

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                      1 },

  { "--alph",       eslARG_STRING,     "dna", NULL, NULL,    NULL,  NULL,  NULL,             "alphabet [dna,dna_full,amino]",                             2 },
  { "--bin_length", eslARG_INT,        "256", NULL, NULL,    NULL,  NULL,  NULL,             "bin length (power of 2;  32<=b<=4096)",                     2 },
  { "--sa_freq",    eslARG_INT,        "8",   NULL, NULL,    NULL,  NULL,  NULL,             "suffix array sample rate (power of 2)",                     2 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[options] <seqfile> <fmfile>";
static char banner[] = "build a HMMER binary-formatted database from an input sequence file";


static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_seqfile, char **ret_fmfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     p7_Die("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nBasic options:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/

      puts("\nSpecial options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); /* 2= group; 2 = indentation; 120=textwidth*/

      exit(0);
  }

  if (esl_opt_ArgNumber(go)                  != 2)     { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <seqfile> argument on command line"); goto ERROR; }
  if ((*ret_fmfile  = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <fmfile> argument on command line");   goto ERROR; }

  /* Validate any attempted use of stdin streams */
  if (esl_strcmp(*ret_seqfile, "-") == 0 && esl_strcmp(*ret_fmfile, "-") == 0) {
    puts("Either <seqfile> or <fmfile> may be '-' (to read from stdin), but not both.");
    goto ERROR;
  }

  *ret_go = go;
  return;

 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere basic options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);
}

static int
output_header(FILE *ofp, const ESL_GETOPTS *go, char *seqfile, char *fmfile)
{
  p7_banner(ofp, go->argv[0], banner);

  fprintf(ofp, "# input sequence file:                     %s\n", seqfile);
  fprintf(ofp, "# output binary-formatted HMMER database:  %s\n", fmfile);
  fprintf(ofp, "# alphabet     :                           %s\n", esl_opt_GetString(go, "--alph"));
  fprintf(ofp, "# bin_length   :                           %d\n", esl_opt_GetInteger(go, "--bin_length"));
  fprintf(ofp, "# suffix array sample rate:                %d\n", esl_opt_GetInteger(go, "--sa_freq"));
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
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
	uint32_t *SAsamp     = NULL;
	uint32_t *occCnts_sb = NULL; // same indexing as above
	uint32_t *cnts_sb    = NULL;
	uint16_t *occCnts_b  = NULL; // this is logically a 2D array, but will be indexed as occ_cnts[alph_size*index + char]  (instead of occ_cnts[index][char])
	uint16_t *cnts_b     = NULL;
	FM_METADATA *meta    = NULL;
	char *inv_alph       = NULL;
	char *alph           = NULL;

        clock_t t1, t2;
        struct tms ts1, ts2;

	long i,j, c;
	int status = eslOK;
	int result;
	int opt = 0;
	int longIndex = 0;

	int num_freq_cnts_sb ;
	int num_freq_cnts_b ;
	int num_SA_samples ;


	int chars_per_byte;
	char *fname_in = NULL;
	char *fname_out= NULL;

	ESL_GETOPTS     *go  = NULL;    /* command line processing                 */

	ESL_ALLOC (meta, sizeof(FM_METADATA));
	meta->alph_type   = fm_DNA;
	meta->freq_SA     = 8;
	meta->freq_cnt_b  = 256;
	meta->freq_cnt_sb = pow(2,16); //65536 - that's the # values in a short

	process_commandline(argc, argv, &go, &fname_in, &fname_out);

	if (esl_opt_IsOn(go, "--alph")) { alph    = esl_opt_GetString(go, "--alph") ; }
	if ( esl_strcmp(alph, "dna")==0)
	  meta->alph_type = fm_DNA;
	else if (esl_strcmp(alph, "dna_full")==0)
	  meta->alph_type = fm_DNA_full;
	else if (esl_strcmp(alph, "amino")==0)
	  meta->alph_type = fm_AMINO;
	else
	  esl_fatal("Unknown alphabet type. Try 'dna', 'dna_full', or 'amino'\n%s", "");

    alph = NULL;

	if (esl_opt_IsOn(go, "--bin_length")) meta->freq_cnt_b = esl_opt_GetInteger(go, "--bin_length");
	if ( meta->freq_cnt_b < 32 || meta->freq_cnt_b >4096 ||  (meta->freq_cnt_b & (meta->freq_cnt_b - 1))  ) // test power of 2
	  esl_fatal("bin_length must be a power of 2, at least 128, and at most 4096\n");

	if (esl_opt_IsOn(go, "--sa_freq")) meta->freq_SA = esl_opt_GetInteger(go, "--sa_freq");
	if ( (meta->freq_SA & (meta->freq_SA - 1))  )  // test power of 2
		esl_fatal ("SA_freq must be a power of 2\n");


        //start timer
        t1 = times(&ts1);

	/* Open a file for reading. */
	if((fp = fopen(fname_in, "rb")) == NULL) {
	esl_fatal("%s: Cannot open file `%s': ", argv[0], fname_in);
	}

	/* Get the file size. */
	//TODO: this is where the loop should start, just use a sequence block
	if(fseek(fp, 0, SEEK_END) == 0) {
	  meta->N = ftell(fp);
	  rewind(fp);
	  if(meta->N < 0)
	    esl_fatal("%s: Cannot ftell `%s': ", argv[0], fname_in);
	} else {
	  esl_fatal( "%s: Cannot fseek `%s': ", argv[0], fname_in);
	}


	output_header(stdout, go, fname_in, fname_out);

	meta->SA_shift = log2(meta->freq_SA);
	meta->cnt_shift_b = log2(meta->freq_cnt_b);

	num_freq_cnts_b  = 1+ceil((float)meta->N/meta->freq_cnt_b);
	num_SA_samples   = 1+floor((float)meta->N/meta->freq_SA);

	meta->cnt_shift_sb = log2(meta->freq_cnt_sb);
	num_freq_cnts_sb = 1+ceil((float)meta->N/meta->freq_cnt_sb);


	//getInverseAlphabet
	fm_createAlphabet(meta->alph_type, &alph, &inv_alph, &(meta->alph_size), &(meta->charBits));



    //shift inv_alph up one, to make space for '$' at 0
	for (i=0; i<256; i++)
		if ( inv_alph[i] >= 0)
			inv_alph[i]++;


	/* Allocate n bytes for BWT and 4n bytes for SA*/
	ESL_ALLOC (T, meta->N * sizeof(uint8_t));
	ESL_ALLOC (BWT, meta->N * sizeof(uint8_t));
	ESL_ALLOC (SA, meta->N * sizeof(int));
	ESL_ALLOC (SAsamp, num_SA_samples * sizeof(uint32_t));
	ESL_ALLOC (occCnts_sb, num_freq_cnts_sb *  meta->alph_size * sizeof(uint32_t)); // every freq_cnt_sb positions, store an array of ints
	ESL_ALLOC (cnts_sb,    meta->alph_size * sizeof(uint32_t));
	ESL_ALLOC (occCnts_b,  num_freq_cnts_b *  meta->alph_size * sizeof(uint16_t)); // every freq_cnt_b positions, store an array of 8-byte ints
	ESL_ALLOC (cnts_b,     meta->alph_size * sizeof(uint16_t));

	if((T == NULL)  || (SA==NULL) || (BWT==NULL) || (cnts_b==NULL) || (occCnts_b==NULL) || (cnts_sb==NULL) || (occCnts_sb==NULL) ) {
	esl_fatal( "%s: Cannot allocate memory.\n", argv[0]);
	}


	/* Read n bytes of data. */
	if(fread(T, sizeof(unsigned char), (size_t)meta->N, fp) != (size_t)meta->N) {
	esl_fatal("%s: %s '%s'", argv[0],
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
		  esl_fatal("requested alphabet doesn't match input test: T[%ld] = (%d), T[%ld] = (%d)\n",
				  meta->N-2-i, T[meta->N-2-i], i, T[i]);

	  T[i] = tmpc1;
	  T[meta->N-2-i] = tmpc2;

	}
	if ( ! (meta->N & 0x1) )  //odd # chars, need to fix the middle one
	  T[meta->N/2 - 1 ] = inv_alph[T[meta->N/2 - 1 ]];

	/* Construct the Suffix Array */
	result = divsufsort(T, SA, meta->N);
	if ( result < 0 )
	  esl_fatal("%s: Error building BWT.\n", argv[0]);


	/* Construct the BWT, SA landmarks, and FM-index */
	for (c=0; c<meta->alph_size; c++) {
	  cnts_sb[c] = 0;
	  cnts_b[c] = 0;
	  FM_OCC_CNT(sb, 0, c ) = 0;
	  FM_OCC_CNT(b, 0, c ) = 0;
	}

	BWT[0] =  SA[0]==0 ? 0 /* '$' */ : T[ SA[0]-1]-1 ; //move values down so 'a'=0...'t'=3; store 'a' in place of '$'
	cnts_sb[BWT[0]]++;
	cnts_b[BWT[0]]++;

	//Scan through SA to build the BWT and FM index structures
	for(j=1; j < meta->N; ++j) {
	  if (SA[j]==0) { //'$'
		  meta->term_loc = j;
		  BWT[j] =  0; //store 'a' in place of '$'
	  } else {
	      BWT[j] =  T[ SA[j]-1]-1 ; //move values down so 'a'=0...'t'=3;
	  }

	  //sample the SA
	  if ( !(j % meta->freq_SA) )
		  SAsamp[ j>>meta->SA_shift ] = ( SA[j] == meta->N - 1 ? -1 : SA[j] ) ; // handle the wrap-around '$'

	  cnts_sb[BWT[j]]++;
	  cnts_b[BWT[j]]++;

	  const long joffset = j+1;
	  if ( !(  joffset % meta->freq_cnt_b) ) {  // (j+1)%freq_cnt_b==0  , i.e. every freq_cnt_bth position, noting that it's a zero-based count

		  for (c=0; c<meta->alph_size; c++) {
			  FM_OCC_CNT(b, (joffset>>meta->cnt_shift_b), c ) = cnts_b[c];
		  }

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


#ifdef PRINTOCC
	printf ("SB\n---\n");
	for(j=0; j < num_freq_cnts_sb; ++j) {
		int sum = 0;
		printf("%ld  ", j);
		for (c=0; c<meta->alph_size; c++) {
			printf("%8d  ", FM_OCC_CNT(sb, j, c ));
			sum+=FM_OCC_CNT(sb, j, c );
		}
		printf("%8d\n", sum);
	}

	printf ("\nB\n---\n");
	for(j=0; j < num_freq_cnts_b; ++j) {
		int sum = 0;
		printf("%ld  ", j);
		for (c=0; c<meta->alph_size; c++) {
			printf("%8d  ", FM_OCC_CNT(b, j, c ));
			sum+=FM_OCC_CNT(b, j, c );
		}
		printf("%8d\n", sum);
	}


#endif

#ifdef PRINTBWT
	for(j=0; j < meta->N; ++j)
	  printf ("%d", BWT[j]);
#endif

	/* Convert BWT and T to packed versions if appropriate. */
	if (meta->alph_type == fm_DNA) {
		 //4 chars per byte.  Counting will be done based on quadruples 0..3; 4..7; 8..11; etc.
		  for(i=0; i < meta->N-3; i+=4) {
			  BWT[i>>2] = BWT[i]<<6 | BWT[i+1]<<4 | BWT[i+2]<<2 | BWT[i+3];
			    T[i>>2] =   T[i]<<6 |   T[i+1]<<4 |   T[i+2]<<2 | T[i+3];
		  }
		  if (i>=meta->N-3) {
			  BWT[i>>2] =  BWT[i]<<6;
			    T[i>>2] =    T[i]<<6;
		  }
		  if (i>=meta->N-2) {
			  BWT[i>>2] =  BWT[i+1]<<4;
			    T[i>>2] =    T[i+1]<<4;
		  }
		  if (i==meta->N-1)  {
			  BWT[i>>2] =  BWT[i+2]<<2;
			    T[i>>2] =    T[i+2]<<2;
		  }
	} else if (meta->alph_type == fm_DNA_full) {
		//2 chars per byte.  Counting will be done based on quadruples 0..3; 4..7; 8..11; etc.
		  for(i=0; i < meta->N-1; i+=2) {
			  BWT[i>>1] = BWT[i]<<4 | BWT[i+1];
			    T[i>>1] =   T[i]<<4 |   T[i+1];
		  }
		  if (i==meta->N-1) {
			  BWT[i>>1] =  BWT[i]<<4 ;
			    T[i>>1] =    T[i]<<4 ;
		  }
	}

	chars_per_byte = 8/meta->charBits;
	meta->L = (size_t)((chars_per_byte-1+meta->N)/chars_per_byte);


	/* Open a file for writing. */
	if((fp = fopen(fname_out, "wb")) == NULL)
	  esl_fatal( "%s: Cannot open file `%s': ", argv[0], fname_out);
	/* Write the FM-index meta data */
	if(fwrite(meta, sizeof(FM_METADATA), 1, fp) != 1)
	  esl_fatal( "%s: Error writing meta data for FM index.\n", argv[0]);

	if(fwrite(T, sizeof(uint8_t), meta->L, fp) != meta->L)
	  esl_fatal( "%s: Error writing T in FM index.\n", argv[0]);
	if(fwrite(BWT, sizeof(uint8_t), meta->L, fp) != meta->L)
	  esl_fatal( "%s: Error writing BWT in FM index.\n", argv[0]);
	if(fwrite(SAsamp, sizeof(uint32_t), (size_t)num_SA_samples, fp) != (size_t)num_SA_samples)
	  esl_fatal( "%s: Error writing SA in FM index.\n", argv[0]);
	if(fwrite(occCnts_b, sizeof(uint16_t)*(meta->alph_size), (size_t)num_freq_cnts_b, fp) != (size_t)num_freq_cnts_b)
	  esl_fatal( "%s: Error writing occCnts_b in FM index.\n", argv[0]);
	if(fwrite(occCnts_sb, sizeof(uint32_t)*(meta->alph_size), (size_t)num_freq_cnts_sb, fp) != (size_t)num_freq_cnts_sb)
	  esl_fatal( "%s: Error writing occCnts_sb in FM index.\n", argv[0]);

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

    esl_getopts_Destroy(go);

    // compute and print the elapsed time in millisec
    t2 = times(&ts2);
    {
      double clk_ticks = sysconf(_SC_CLK_TCK);
      double elapsedTime = (t2-t1)/clk_ticks;

      fprintf (stderr, "run time:  %.2f seconds\n", elapsedTime);
    }


	return (eslOK);


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
