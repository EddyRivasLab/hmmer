#include "divsufsort.h"
#include "hmmer.h"
#include "easel.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_mem.h"
#include <string.h>

//#define PRINTBWT 1
//#define PRINTOCC 1

#define FM_BLOCK_COUNT 100000 //max number of SQ objects in a block
#define FM_BLOCK_OVERLAP 100000 //100 Kbases of overlap, at most, between adjascent FM-index blocks

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  NULL,       "show brief help on version and usage",                      1 },

  { "--informat",   eslARG_STRING,     FALSE, NULL, NULL, NULL, NULL,      NULL,        "specify that input file is in format <s>",                  2 },
  { "--alph",       eslARG_STRING,     "dna", NULL, NULL,    NULL,  NULL,  NULL,        "alphabet [dna,dna_full,amino]",                             2 },
  { "--bin_length", eslARG_INT,        "256", NULL, NULL,    NULL,  NULL,  NULL,        "bin length (power of 2;  32<=b<=4096)",                     2 },
  { "--sa_freq",    eslARG_INT,        "8",   NULL, NULL,    NULL,  NULL,  NULL,        "suffix array sample rate (power of 2)",                     2 },
  { "--block_size", eslARG_INT,        "15",  NULL, NULL,    NULL,  NULL,  NULL,        "input sequence broken into chunks this size (Mbases)",      2 },


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



int
ensure_meta_alloc (FM_METADATA *meta, int numseqs, int *allocedseqs) {
    int status = eslOK;
	int to_alloc = *allocedseqs;
	if (numseqs == *allocedseqs) {
		to_alloc *= 2;
	}

	if (numseqs == 0 || numseqs == *allocedseqs) { // either first allocation, or increase in size
		ESL_REALLOC (meta->name_lengths, to_alloc * sizeof(uint16_t*));
		ESL_REALLOC (meta->names,        to_alloc * sizeof(char*));
		ESL_REALLOC (meta->starts,       to_alloc * sizeof(uint32_t*));
		ESL_REALLOC (meta->lengths,      to_alloc * sizeof(uint32_t*));
		ESL_REALLOC (meta->ids,          to_alloc * sizeof(uint32_t*));
		if (meta->names == NULL || meta->starts == NULL || meta->lengths == NULL || meta->ids == NULL)
			esl_fatal("unable to allocate memory to store FM meta data\n");
	}

	*allocedseqs = to_alloc;

	return eslOK;

ERROR:
	return status;
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

	int chars_per_byte;
	int num_freq_cnts_sb ;
	int num_freq_cnts_b ;
	int num_SA_samples ;

	int             infmt     = eslSQFILE_UNKNOWN;
	int             alphatype = eslUNKNOWN;
	ESL_ALPHABET   *abc       = NULL;
	ESL_SQ         *sq        = NULL;
	ESL_SQFILE     *sqfp      = NULL;

    ESL_SQ       *tmpsq;
    ESL_SQ_BLOCK *block;

	char *fname_in = NULL;
	char *fname_out= NULL;
	int block_size = 15000000;
	int sq_cnt = 0;
	int use_tmpsq = 0;
	int reported_N = 0;
	int block_length;
	int max_block_size;

	int numblocks = 0;
	int allocedblocks = 10;
	int numseqs;
	int allocedseqs = 1000;

	int compressed_bytes;
    int term_loc;
    uint64_t offset = 0;

	ESL_GETOPTS     *go  = NULL;    /* command line processing                 */

	ESL_ALLOC (meta, sizeof(FM_METADATA));
	if (meta == NULL)
		esl_fatal("unable to allocate memory to store FM meta data\n");

	meta->alph_type   = fm_DNA;
	meta->freq_SA     = 8;
	meta->freq_cnt_b  = 256;
	meta->freq_cnt_sb = pow(2,16); //65536 - that's the # values in a short

	ESL_ALLOC (meta->block_offsets,  allocedblocks * sizeof(uint16_t*));

	status = ensure_meta_alloc(meta, 0, &allocedseqs);
	if (status != eslOK) {
        esl_exception(eslEMEM, __FILE__, __LINE__, "failure allocting memory");
        goto ERROR;
	}

	process_commandline(argc, argv, &go, &fname_in, &fname_out);

	if (esl_opt_IsOn(go, "--alph")) { alph    = esl_opt_GetString(go, "--alph") ; }
	if ( esl_strcmp(alph, "dna")==0) {
	  meta->alph_type = fm_DNA;
	  alphatype = eslDNA;
	} else if (esl_strcmp(alph, "dna_full")==0) {
	  meta->alph_type = fm_DNA_full;
	  alphatype = eslDNA;
	} else if (esl_strcmp(alph, "amino")==0) {
	  meta->alph_type = fm_AMINO;
	  alphatype = eslAMINO;
	} else {
	  esl_fatal("Unknown alphabet type. Try 'dna', 'dna_full', or 'amino'\n%s", "");
	}
    alph = NULL;

	if (esl_opt_IsOn(go, "--bin_length")) meta->freq_cnt_b = esl_opt_GetInteger(go, "--bin_length");
	if ( meta->freq_cnt_b < 32 || meta->freq_cnt_b >4096 ||  (meta->freq_cnt_b & (meta->freq_cnt_b - 1))  ) // test power of 2
	  esl_fatal("bin_length must be a power of 2, at least 128, and at most 4096\n");

	if (esl_opt_IsOn(go, "--sa_freq")) meta->freq_SA = esl_opt_GetInteger(go, "--sa_freq");
	if ( (meta->freq_SA & (meta->freq_SA - 1))  )  // test power of 2
		esl_fatal ("SA_freq must be a power of 2\n");


	if (esl_opt_IsOn(go, "--block_size")) block_size = 1000000 * esl_opt_GetInteger(go, "--block_size");
	if ( block_size <=0  )
		esl_fatal ("block_size must be a positive number\n");

    //start timer
    t1 = times(&ts1);



	output_header(stdout, go, fname_in, fname_out);

	meta->SA_shift = log2(meta->freq_SA);
	meta->cnt_shift_b = log2(meta->freq_cnt_b);
	meta->cnt_shift_sb = log2(meta->freq_cnt_sb);


	//getInverseAlphabet
	fm_createAlphabet(meta->alph_type, &alph, &inv_alph, &(meta->alph_size), &(meta->charBits));

    //shift inv_alph up one, to make space for '$' at 0
	for (i=0; i<256; i++)
		if ( inv_alph[i] >= 0)
			inv_alph[i]++;



	output_header(stdout, go, fname_in, fname_out);

	meta->SA_shift = log2(meta->freq_SA);
	meta->cnt_shift_b = log2(meta->freq_cnt_b);
	meta->cnt_shift_sb = log2(meta->freq_cnt_sb);


	//getInverseAlphabet
	fm_createAlphabet(meta->alph_type, &alph, &inv_alph, &(meta->alph_size), &(meta->charBits));

    //shift inv_alph up one, to make space for '$' at 0
	for (i=0; i<256; i++)
		if ( inv_alph[i] >= 0)
			inv_alph[i]++;

    if (esl_opt_GetString(go, "--informat") != NULL) {
      infmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat"));
      if (infmt == eslSQFILE_UNKNOWN) esl_fatal("%s is not a valid input sequence file format for --informat");
    }

    status = esl_sqfile_Open(fname_in, infmt, NULL, &sqfp);
    if      (status == eslENOTFOUND) esl_fatal("No such file %s", fname_in);
    else if (status == eslEFORMAT)   esl_fatal("Format of seqfile %s unrecognized.", fname_in);
    else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

    abc     = esl_alphabet_Create(alphatype);
    sq      = esl_sq_CreateDigital(abc);
    tmpsq   =  esl_sq_CreateDigital(abc);

    esl_sqfile_SetDigital(sqfp, abc);
    block = esl_sq_CreateDigitalBlock(FM_BLOCK_COUNT, abc);

    max_block_size = FM_BLOCK_OVERLAP+block_size;

	/* Allocate BWT, Text, SA, and FM-index data structures, allowing storage of maximally large sequence*/
	ESL_ALLOC (T, max_block_size * sizeof(uint8_t));
	ESL_ALLOC (BWT, max_block_size * sizeof(uint8_t));
	ESL_ALLOC (SA, max_block_size * sizeof(int));
	ESL_ALLOC (SAsamp,     (1+max_block_size/meta->freq_SA) * sizeof(uint32_t));
	ESL_ALLOC (occCnts_sb, (1+ceil((float)max_block_size/meta->freq_cnt_sb)) *  meta->alph_size * sizeof(uint32_t)); // every freq_cnt_sb positions, store an array of ints
	ESL_ALLOC (cnts_sb,    meta->alph_size * sizeof(uint32_t));
	ESL_ALLOC (occCnts_b,  ( 1+ceil((float)max_block_size/meta->freq_cnt_b)) *  meta->alph_size * sizeof(uint16_t)); // every freq_cnt_b positions, store an array of 8-byte ints
	ESL_ALLOC (cnts_b,     meta->alph_size * sizeof(uint16_t));

	if((T == NULL)  || (BWT == NULL)  || (SA==NULL) || (SAsamp==NULL) || (BWT==NULL) || (cnts_b==NULL) || (occCnts_b==NULL) || (cnts_sb==NULL) || (occCnts_sb==NULL) ) {
	esl_fatal( "%s: Cannot allocate memory.\n", argv[0]);
	}


	// Open a file for writing.
	if((fp = fopen(fname_out, "wb")) == NULL)
	  esl_fatal( "%s: Cannot open file `%s': ", argv[0], fname_out);

//	skip 64 bits, which will be used to store the offset of the metadata,
//	once I know where that will be
    if(fseek(fp, sizeof(uint64_t), 0) != 0)
    	esl_fatal("error writing FM-index meta data");

    numblocks = 0;
	numseqs = 0;

    /* Main loop: */
    while (status == eslOK ) {

        //reset block as an empty vessel
        for (i=0; i<block->count; i++)
            esl_sq_Reuse(block->list + i);

        if (use_tmpsq) {
            esl_sq_Copy(tmpsq , block->list);
            block->complete = FALSE;  //this lets ReadBlock know that it needs to append to a small bit of previously-read seqeunce
            block->list->C = FM_BLOCK_OVERLAP; // overload the ->C value, which ReadBlock uses to determine how much
                                                   // overlap should be retained in the ReadWindow step
        }


        status = esl_sqio_ReadBlock(sqfp, block, block_size, TRUE);
        if (status == eslEOF)
        	continue;
        if (status != eslOK) {
            esl_exception(eslEMEM, __FILE__, __LINE__, "failure reading sequence block");
            goto ERROR;
        }


        if (block->complete || block->count == 0) {
            use_tmpsq = FALSE;
        } else {
            /* The final sequence on the block was a probably-incomplete window of the active sequence.
             * We won't use it on this pass through block-construction, but will save it for the
             * next pass.
             */
            esl_sq_Copy(block->list + (block->count - 1) , tmpsq);
            use_tmpsq = TRUE;
        }

        block->first_seqidx = sq_cnt;
        sq_cnt += block->count - (use_tmpsq ? 1 : 0);// if there's an incomplete sequence read into the block wait to count it until it's complete.



    	/* Read dseqs from block into text element T.
    	*  Convert the dsq from esl-alphabet to fm-alphabet (1..k for alphabet of size k).
    	*  (a) collapsing upper/lower case for appropriate sorting.
    	*  (b) reserving 0 for '$', which must be lexicographically smallest
    	*      (these will later be shifted to 0-based alphabet, once SA has been built)
    	*
    	*/
        block_length = 0;
        for (i=0; i<block->count; i++) {

        	ensure_meta_alloc(meta, numseqs, &allocedseqs);
        	if (status != eslOK) {
                esl_exception(eslEMEM, __FILE__, __LINE__, "failure allocting memory");
                goto ERROR;
        	}
			meta->name_lengths[numseqs] = strlen(block->list[i].name) ;
			ESL_REALLOC (meta->names[numseqs], meta->name_lengths[numseqs] * sizeof(char));
        	if (meta->names[numseqs] == NULL)
        		esl_fatal("unable to allocate memory to store FM meta data(3)\n");

        	//meta data
        	esl_memstrcpy(block->list[i].name, meta->name_lengths[numseqs], meta->names[numseqs]);
        	meta->ids[numseqs]   = block->first_seqidx + i;
        	meta->starts[numseqs] = block->list[i].start;


        	for (j=1; j<=block->list[i].n; j++) {
        		c = abc->sym[block->list[i].dsq[j]];
        		if ( meta->alph_type == fm_DNA) {
        			if (inv_alph[c] == -1) {
        				if (!reported_N) {
							printf ("You have selected alph_type 'dna', but your sequence database contains ambiguity code 'N'\n");
							printf ("the database will be built, but hits including these positions will not be found.\n");
							printf ("Use alph_type 'dna_full' if this bothers you.\n");
							reported_N = 1;
        				}

        				numseqs++;

        				//start a new block
        				ensure_meta_alloc(meta, numseqs, &allocedseqs);
        				if (status != eslOK) {
        			        esl_exception(eslEMEM, __FILE__, __LINE__, "failure allocting memory");
        			        goto ERROR;
        				}
        				meta->name_lengths[numseqs] = strlen(block->list[i].name) ;
        				ESL_REALLOC (meta->names[numseqs], meta->name_lengths[numseqs] * sizeof(char));
        	        	if ( meta->names[numseqs] == NULL)
        	        		esl_fatal("unable to allocate memory to store FM meta data(3)\n");

        	        	//meta data
        	        	esl_memstrcpy(block->list[i].name, meta->name_lengths[numseqs], meta->names[numseqs]);
        	        	meta->ids[numseqs]   = block->first_seqidx + i ;
        	        	meta->starts[numseqs] = block->list[i].start + j;
        			}
        		} else if (inv_alph[c] == -1) {
    	    		esl_fatal("requested alphabet doesn't match input text\n");
        		}

    	    	T[block_length] = inv_alph[c];

    	    	block_length++;
    	    	meta->lengths[numseqs]++;

        	}
        	numseqs++;
        }
    	T[block_length] = 0; // last character 0 is effectively '$' for suffix array


    	num_freq_cnts_b  = 1+ceil((float)block_length/meta->freq_cnt_b);
    	num_freq_cnts_sb = 1+ceil((float)block_length/meta->freq_cnt_sb);
    	num_SA_samples   = 1+floor((float)block_length/meta->freq_SA);

    	meta->seq_count = numseqs;



        //###########################################
    	//     build Suffix Array, then FM-index
    	//###########################################

    	chars_per_byte = 8/meta->charBits;
    	compressed_bytes = 	((chars_per_byte-1+block_length)/chars_per_byte);

    	// Construct the Suffix Array
    	status = divsufsort(T, SA, block_length);
    	if ( status < 0 )
    	  esl_fatal("%s: Error building BWT.\n", argv[0]);


    	// Construct the BWT, SA landmarks, and FM-index
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
    	for(j=1; j < block_length; ++j) {
    	  if (SA[j]==0) { //'$'
    		  term_loc = j;
    		  BWT[j] =  0; //store 'a' in place of '$'
    	  } else {
    	      BWT[j] =  T[ SA[j]-1]-1 ; //move values down so 'a'=0...'t'=3;
    	  }

    	  //sample the SA
    	  if ( !(j % meta->freq_SA) )
    		  SAsamp[ j>>meta->SA_shift ] = ( SA[j] == block_length - 1 ? -1 : SA[j] ) ; // handle the wrap-around '$'

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



    	// Convert BWT and T to packed versions if appropriate.
    	if (meta->alph_type == fm_DNA) {
    		 //4 chars per byte.  Counting will be done based on quadruples 0..3; 4..7; 8..11; etc.
    		  for(i=0; i < block_length-3; i+=4) {
    			  BWT[i>>2] = BWT[i]<<6 | BWT[i+1]<<4 | BWT[i+2]<<2 | BWT[i+3];
    			    T[i>>2] =   T[i]<<6 |   T[i+1]<<4 |   T[i+2]<<2 | T[i+3];
    		  }
    		  if (i>=block_length-3) {
    			  BWT[i>>2] =  BWT[i]<<6;
    			    T[i>>2] =    T[i]<<6;
    		  }
    		  if (i>=block_length-2) {
    			  BWT[i>>2] =  BWT[i+1]<<4;
    			    T[i>>2] =    T[i+1]<<4;
    		  }
    		  if (i==block_length-1)  {
    			  BWT[i>>2] =  BWT[i+2]<<2;
    			    T[i>>2] =    T[i+2]<<2;
    		  }
    	} else if (meta->alph_type == fm_DNA_full) {
    		//2 chars per byte.  Counting will be done based on quadruples 0..3; 4..7; 8..11; etc.
    		  for(i=0; i < block_length-1; i+=2) {
    			  BWT[i>>1] = BWT[i]<<4 | BWT[i+1];
    			    T[i>>1] =   T[i]<<4 |   T[i+1];
    		  }
    		  if (i==block_length-1) {
    			  BWT[i>>1] =  BWT[i]<<4 ;
    			    T[i>>1] =    T[i]<<4 ;
    		  }
    	}


    	// Write the FM-index meta data
    	if(fwrite(&block_length, sizeof(uint32_t), 1, fp) !=  1)
    		esl_fatal( "%s: Error writing block_length in FM index.\n", argv[0]);
    	if(fwrite(&term_loc, sizeof(uint32_t), 1, fp) !=  1)
    		esl_fatal( "%s: Error writing terminal location in FM index.\n", argv[0]);

    	if(fwrite(T, sizeof(uint8_t), compressed_bytes, fp) != compressed_bytes)
    	  esl_fatal( "%s: Error writing T in FM index.\n", argv[0]);
    	if(fwrite(BWT, sizeof(uint8_t), compressed_bytes, fp) != compressed_bytes)
    	  esl_fatal( "%s: Error writing BWT in FM index.\n", argv[0]);
    	if(fwrite(SAsamp, sizeof(uint32_t), (size_t)num_SA_samples, fp) != (size_t)num_SA_samples)
    	  esl_fatal( "%s: Error writing SA in FM index.\n", argv[0]);
    	if(fwrite(occCnts_b, sizeof(uint16_t)*(meta->alph_size), (size_t)num_freq_cnts_b, fp) != (size_t)num_freq_cnts_b)
    	  esl_fatal( "%s: Error writing occCnts_b in FM index.\n", argv[0]);
    	if(fwrite(occCnts_sb, sizeof(uint32_t)*(meta->alph_size), (size_t)num_freq_cnts_sb, fp) != (size_t)num_freq_cnts_sb)
    	  esl_fatal( "%s: Error writing occCnts_sb in FM index.\n", argv[0]);


    	if (numblocks == allocedblocks) {
    		allocedblocks *= 2;
    		ESL_REALLOC(meta->block_offsets, allocedblocks);
    	}
    	meta->block_offsets[numblocks] = offset;

    	offset +=  2*sizeof(uint32_t)
    	         + 2*sizeof(uint8_t)*compressed_bytes
    	         + sizeof(uint32_t)*num_SA_samples
    	         + sizeof(uint16_t)*(meta->alph_size)*num_freq_cnts_b
    	         + sizeof(uint32_t)*(meta->alph_size)*num_freq_cnts_sb ;

    	numblocks++;

    }


	if(fwrite(meta, sizeof(FM_METADATA), 1, fp) != 1)
	  esl_fatal( "%s: Error writing meta data for FM index.\n", argv[0]);

	//go back to beginning of file, and store the location at which the meta data was written
	fseek(fp, 0, 0);
	if(fwrite(&offset, sizeof(uint64_t), 1, fp) != 1)
	  esl_fatal( "%s: Error writing meta data for FM index.\n", argv[0]);






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
