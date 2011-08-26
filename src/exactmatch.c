#include "hmmer.h"
#include <sys/times.h>

#include "easel.h"
#include <string.h>

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,      FALSE, NULL, NULL,    NULL,  NULL,  NULL,          "show brief help on version and usage",                      1 },

  { "--out",       eslARG_STRING,     "none", NULL, NULL,    NULL,  NULL,  NULL,          "save list of hits to file <s>  ('-' writes to stdout)",     2 },
  { "--count_only", eslARG_NONE,      FALSE, NULL, NULL,    NULL,  NULL,  NULL,          "compute just counts, not locations",                        2 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


static char usage[]  = "[options] <fmfile> <qfile>";
static char banner[] = "Find all instances of each <qfile> sequence in the database represented by <fmfile>";

int occCallCnt;

int num_freq_cnts_b ;
int num_freq_cnts_sb;

int32_t *Cvals;
int maskSA;
int shiftSA;


static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_fmfile, char **ret_qfile)
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
  if ((*ret_fmfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <fmfile> argument on command line"); goto ERROR; }
  if ((*ret_qfile  = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <qfile> argument on command line");   goto ERROR; }

  /* Validate any attempted use of stdin streams */
  if (esl_strcmp(*ret_fmfile, "-") == 0 && esl_strcmp(*ret_qfile, "-") == 0) {
    puts("Either <fmfile> or <qfile> may be '-' (to read from stdin), but not both.");
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
output_header(FILE *ofp, const ESL_GETOPTS *go, char *fmfile, char *qfile, FM_METADATA *meta)
{
  p7_banner(ofp, go->argv[0], banner);

  fprintf(ofp, "# input binary-formatted HMMER database:   %s\n", fmfile);
  fprintf(ofp, "# input file of query sequences:           %s\n", qfile);

  if (esl_opt_IsUsed(go, "--out")) {
  fprintf(ofp, "# output file containing list of hits:     ");
	  char *outfile = esl_opt_GetString(go, "--out");
	  if (esl_strcmp(outfile, "-"))
		  fprintf(ofp, "stdout\n");
	  else
		  fprintf(ofp, "%s\n", outfile);
  }

  if (esl_opt_IsUsed(go, "--count_only"))
	  fprintf(ofp, "# output only counts, not hit locations\n");

  char *alph;
  if (meta->alph_type == fm_DNA)
	  alph = "dna";
  else if (meta->alph_type == fm_DNA_full)
	  alph = "dna_full";
  else if (meta->alph_type == fm_AMINO)
  	  alph = "amino";
  fprintf(ofp, "# alphabet     :                           %s\n", alph);

  fprintf(ofp, "# bin_length   :                           %d\n", meta->freq_cnt_b);
  fprintf(ofp, "# suffix array sample rate:                %d\n", meta->freq_SA);
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}



//see: http://c-faq.com/stdio/commaprint.html
char *
commaprint(unsigned long n) {
        static int comma = ',';
        static char retbuf[30];
        char *p = &retbuf[sizeof(retbuf)-1];
        int i = 0;
        *p = '\0';

        do {
                if(i%3 == 0 && i != 0)
                        *--p = comma;
                *--p = '0' + n % 10;
                n /= 10;
                i++;
        } while(n != 0);

        return p;
}

/* Function:  getFMRange()
 * Synopsis:  For a given query sequence, find its interval in the FM-index
 * Purpose:   Implement Algorithm 3.6 (p17) of Firth paper (A Comparison of BWT Approaches
 *            to String Pattern Matching). All the meat is in the method of counting
 *            characters - bwt_getOccCount, which depends on compilation choices.
 */
int
getFMRange( FM_METADATA *meta, FM_DATA *fm, char *query, char *inv_alph, FM_INTERVAL *interval) {

	int count1, count2;
	int i=0;

	uint8_t c = inv_alph[(int)query[0]];
	interval->lower = Cvals[c] + (c==0?1:0);
	interval->upper  = Cvals[c+1]-1;

	while (interval->lower>0 && interval->lower <= interval->upper) {
		c = query[++i];
		if (c == '\0')  // end of query - the current range defines the hits
			break;

		c = inv_alph[c];

		//TODO: these will often overlap - might get acceleration by merging to a single redundancy-avoiding call
		count1 = fm_getOccCount (meta, fm, interval->lower-1, c);
		count2 = fm_getOccCount (meta, fm, interval->upper, c);

// 		   printf("%d %d %d -> %d,%d\n", c, interval->lower-1, interval->upper, count1, count2);

		interval->lower = Cvals[c] + count1;
		interval->upper = Cvals[c] + count2 - 1;

		occCallCnt+=2;
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
getChar( FM_METADATA *meta, int j, const uint8_t *B ) {
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
 * Purpose:   Implement Algorithm 3.7 (p17) of Firth paper (A Comparison of BWT Approaches
 *            to String Pattern Matching). Most of the meat is in the method of counting
 *            characters - bwt_getOccCount, which depends on compilation choices.
 */
#ifndef FMDEBUG
inline
#endif
int
getFMHits( FM_METADATA *meta, FM_DATA *fm, FM_INTERVAL *interval, FM_HIT *hits_ptr) {

	const uint8_t *B = fm->BWT;
	int i, j, len = 0;

	for (i = interval->lower;  i<= interval->upper; i++) {
		j = i;
		len = 0;

		while ( j & maskSA ) { //go until we hit a position in the full SA that was sampled during FM index construction
			uint8_t c = getChar( meta, j, B);
			j = fm_getOccCount (meta, fm, j-1, c);
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
	int compressed_bytes;
	int chars_per_byte;

	if((fp = fopen(fname, "rb")) == NULL)
		esl_fatal("Cannot open file `%s': ", fname);

	//get the FM meta data
	if(fread(meta, sizeof(FM_METADATA), 1, fp) != 1)
		esl_fatal("Error reading BWT size.%s\n", " ");


	num_freq_cnts_sb = 1+ceil((float)fm->N/meta->freq_cnt_sb);
	num_freq_cnts_b  = 1+ceil((float)fm->N/meta->freq_cnt_b);
	num_SA_samples   = 1+floor((float)fm->N/meta->freq_SA);


	chars_per_byte = 8/meta->charBits;
	compressed_bytes = 	((chars_per_byte-1+fm->N)/chars_per_byte);

	// allocate and read the data
	bwtSize = compressed_bytes  * sizeof(uint8_t);
        ESL_ALLOC (fm->T, bwtSize );
	ESL_ALLOC (fm->BWT_mem, bwtSize + 31 ); // +31 for manual 16-byte alignment  ( typically only need +15, but this allows offset in memory, plus offset in case of <16 bytes of characters at the end)
        fm->BWT = 	(uint8_t *) (((unsigned long int)fm->BWT_mem + 15) & (~0xf));   // align vector memory on 16-byte boundaries
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
	C[0] = 1;


	Cvals = NULL; 	//capture abs() for C array (so I don't have to do it for every query sequence)
	ESL_ALLOC(Cvals, (meta->alph_size+1) * sizeof(int32_t) );
	for(i=0; i<=meta->alph_size; i++)
		Cvals[i] = abs(fm->C[i]);



	return eslOK;

ERROR:
    if (fm->T)           free (fm->T);
    if (fm->BWT_mem)     free (fm->BWT_mem);
	if (fm->SA)          free (fm->SA);
	if (fm->C)           free (fm->C);
	if (fm->occCnts_b)   free (fm->occCnts_b);
	if (fm->occCnts_sb)  free (fm->occCnts_sb);
	esl_fatal("Error allocating memory in %s\n", "readFM");
    return eslFAIL;
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

	char *fname_fm      = NULL;
	char *fname_queries = NULL;
	char *inv_alph      = NULL;
	char *alph          = NULL;
	FM_HIT *hits        = NULL;
	char *line          = NULL;
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
    char *outname = NULL;

	ESL_GETOPTS     *go  = NULL;    /* command line processing                 */
	process_commandline(argc, argv, &go, &fname_fm, &fname_queries);



	if (esl_opt_IsOn(go, "--out")) {
		outname = esl_opt_GetString(go, "--out");
		if ( esl_strcmp ("-", outname) == 0 ) {
			out = stdout;
			outname = "stdout";
		} else {
			out = fopen(optarg,"w");
		}
	}

	if (esl_opt_IsOn(go, "--count_only"))
		count_only = 1;

	occCallCnt = 0;


	readFM ( fname_fm, &meta, &fm );

	output_header(stdout, go, fname_fm, fname_queries, &meta);


	/* initialize a few global variables, then call initGlobals
	 * to do architecture-specific initialization
	 */
	maskSA       =  meta.freq_SA - 1;
	shiftSA      =  meta.SA_shift;
	fm_initGlobals(&meta, &fm);


	fm_createAlphabet(meta.alph_type, &alph, &inv_alph, &(meta.alph_size), NULL); // don't override charBits


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

			//if (!count_only)
			//	getFMHits(&meta, &fm, &interval, hits);

			hit_cnt++;
			hit_indiv_cnt += hit_num;

//			printf ("%s : %d\n", line, interval.upper-interval.lower+1);


			//PRINTOUT ( "HIT:  %s  (%d, %d)\n", line, interval.lower, interval.upper);
			//int i;
			//for (i = 0; i< hit_num; i++) {
			//	int pos_fwd = meta.N - hits[i].start - qlen - 1; // flip the sequence reversal
				//PRINTOUT ( "\t%d\n", pos_fwd);
			//}
		} else {
//			printf ("no hit\n\n");
			//PRINTOUT ( "MISS: %s\n", line);
			miss_cnt++;
		}


	}

    free (hits);
    fclose(fp);
    fm_destroyGlobals();

    // compute and print the elapsed time in millisec
    t2 = times(&ts2);
    {
      double clk_ticks = sysconf(_SC_CLK_TCK);
      double elapsedTime = (t2-t1)/clk_ticks;
      double throughput = occCallCnt/elapsedTime;

      fprintf (stderr, "hit: %-10d  (%d)\n", hit_cnt, hit_indiv_cnt);
      fprintf (stderr, "miss:%-10d\n", miss_cnt);
      fprintf (stderr, "run time:  %.2f seconds\n", elapsedTime);
      fprintf (stderr, "occ calls: %12s\n", commaprint(occCallCnt));
      fprintf (stderr, "occ/sec:   %12s\n", commaprint(throughput));
    }

	exit(eslOK);


ERROR:
	printf ("failure allocating memory for hits\n");
	exit(status);


}





/*****************************************************************
 * @LICENSE@
 *****************************************************************/
