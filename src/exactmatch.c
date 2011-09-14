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
output_header(FM_METADATA *meta, FILE *ofp, const ESL_GETOPTS *go, char *fmfile, char *qfile)
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

/* Function:  getSARangeForward()
 * Synopsis:  For a given query sequence, find its interval in the FM-index, using forward search
 * Purpose:   Implement Algorithm 4 (i371) of Simpson (Bioinformatics 2010). I's the forward
 *            search on a bi-directional BWT, as described by Lam 2009.
 *            All the meat is in the method of counting characters - bwt_getOccCount, which
 *            depends on compilation choices.
 *
 *            Note: it seems odd, but is correct, that the fm-index passed in to this function
 *            is the backward index corresponding to the forward index into which I want to
 *            do a forward search
 */
int
getSARangeForward(FM_DATA *fm, FM_MISC_VARS *misc, char *query, char *inv_alph, FM_INTERVAL *interval) {

  uint32_t occLT_l, occLT_u, occ_l, occ_u;
  int i=0;
  int lower_b, upper_b;

  uint8_t c = inv_alph[(int)query[0]];
  interval->lower  = lower_b = abs(fm->C[c]);
  interval->upper  = upper_b = abs(fm->C[c+1])-1;

  while (lower_b>0 && lower_b <= upper_b) {
    c = query[++i];
    if (c == '\0')  // end of query - the current range defines the hits
      break;

    c = inv_alph[c];

    fm_getOccCountLT (fm, misc, lower_b-1, c, &occ_l, &occLT_l);
    fm_getOccCountLT (fm, misc, upper_b,   c, &occ_u, &occLT_u);

    interval->lower += (occLT_u - occLT_l);
    interval->upper = interval->lower + (occ_u - occ_l) - 1;

    lower_b = abs(fm->C[c]) + occ_l;
    upper_b = abs(fm->C[c]) + occ_u - 1;

    misc->occCallCnt+=2;
  }

  return eslOK;
}


/* Function:  getSARangeReverse()
 * Synopsis:  For a given query sequence, find its interval in the FM-index, using backward search
 * Purpose:   Implement Algorithm 3.6 (p17) of Firth paper (A Comparison of BWT Approaches
 *            to String Pattern Matching). This is what Simpson and Lam call "Reverse Search".
 *            All the meat is in the method of counting characters - bwt_getOccCount, which
 *            depends on compilation choices.
 */
int
getSARangeReverse( FM_DATA *fm, FM_MISC_VARS *misc, char *query, char *inv_alph, FM_INTERVAL *interval) {

  int count1, count2;
  int i=0;

  uint8_t c = inv_alph[(int)query[0]];
  interval->lower  = abs(fm->C[c]);
  interval->upper  = abs(fm->C[c+1])-1;

  while (interval->lower>0 && interval->lower <= interval->upper) {
    c = query[++i];
    if (c == '\0')  // end of query - the current range defines the hits
      break;

    c = inv_alph[c];

    //TODO: counting in these calls will often overlap
      // - might get acceleration by merging to a single redundancy-avoiding call
    count1 = fm_getOccCount (fm, misc, interval->lower-1, c);
    count2 = fm_getOccCount (fm, misc, interval->upper, c);

    interval->lower = abs(fm->C[c]) + count1;
    interval->upper = abs(fm->C[c]) + count2 - 1;

    misc->occCallCnt+=2;
  }

  return eslOK;
}



/* Function:  getChar()
 * Synopsis:  Find the character c residing at a given position in the BWT.
 * Purpose:   The returned char is used by getFMHits(), to seed a call to
 *            bwt_getOccCount().
 */
//#ifndef FMDEBUG
//inline
//#endif
uint8_t
getChar(uint8_t alph_type, int j, const uint8_t *B ) {
  uint8_t c = -1;

  if (alph_type == fm_DNA) {
    /*
     *  B[j>>2] is the byte of B in which j is found (j/4)
     *
     *  Let j' be the final two bits of j (j&0x2)
     *  The char bits are the two starting at position 2*j'.
     *  Without branching, grab them by shifting B[j>>2] right 6-2*j' bits,
     *  then masking to keep the final two bits
     */
    c = (B[j>>2] >> ( 0x6 - ((j&0x3)<<1) ) & 0x3);
  } else if (alph_type == fm_DNA_full) {
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
//#ifndef FMDEBUG
//inline
//#endif
int
getFMHits( FM_DATA *fm, FM_MISC_VARS *misc, FM_INTERVAL *interval, int block_id, int hit_offset, FM_HIT *hits_ptr, int fm_direction) {

  int i, j, len = 0;

  for (i = interval->lower;  i<= interval->upper; i++) {
    j = i;
    len = 0;

    while ( j & misc->maskSA ) { //go until we hit a position in the full SA that was sampled during FM index construction
      uint8_t c = getChar( misc->meta->alph_type, j, fm->BWT);
      j = fm_getOccCount (fm, misc, j-1, c);
      j += abs(fm->C[c]);
      len++;
    }
    const int tmp = j >> misc->shiftSA;
    hits_ptr[hit_offset + i - interval->lower].block = block_id;
    hits_ptr[hit_offset + i - interval->lower].start = fm->SA[ tmp ] + len;
    hits_ptr[hit_offset + i - interval->lower].direction = fm_direction;
  }

  return eslOK;

}

/* Function:  freeFM()
 * Synopsis:  release the memory required to store an individual FM-index
 */
void
freeFM ( FM_DATA *fm, int freeSA)
{
  if (fm->T)            free (fm->T);
  if (fm->BWT_mem)      free (fm->BWT_mem);
  if (fm->C)            free (fm->C);
  if (fm->occCnts_b)    free (fm->occCnts_b);
  if (fm->occCnts_sb)   free (fm->occCnts_sb);

  if (freeSA && fm->SA) free (fm->SA);
}

/* Function:  readFM()
 * Synopsis:  Read the FM index off disk
 * Purpose:   Read the FM-index as written by fmbuild.
 *            First read the metadata header, then allocate space for the full index,
 *            then read it in.
 */
int
readFM( FILE *fp, FM_DATA *fm, FM_METADATA *meta, int getAll )
{
  //shortcut variables
  int *C               = NULL;

  int status;
  int i;

  uint16_t *occCnts_b  = NULL;  //convenience variables, used to simplify macro calls
  uint32_t *occCnts_sb = NULL;

  int compressed_bytes;
  int num_freq_cnts_b;
  int num_freq_cnts_sb;
  int num_SA_samples;
  int prevC;
  int cnt;
  int chars_per_byte = 8/meta->charBits;


  if(fread(&(fm->N), sizeof(uint32_t), 1, fp) !=  1)
    esl_fatal( "%s: Error reading block_length in FM index.\n", __FILE__);
  if(fread(&(fm->term_loc), sizeof(uint32_t), 1, fp) !=  1)
    esl_fatal( "%s: Error reading terminal location in FM index.\n", __FILE__);
  if(fread(&(fm->seq_offset), sizeof(uint16_t), 1, fp) !=  1)
    esl_fatal( "%s: Error reading seq_offset in FM index.\n", __FILE__);


  compressed_bytes =   ((chars_per_byte-1+fm->N)/chars_per_byte);
  num_freq_cnts_b  = 1+ceil((float)fm->N/meta->freq_cnt_b);
  num_freq_cnts_sb = 1+ceil((float)fm->N/meta->freq_cnt_sb);
  num_SA_samples   = 1+floor((float)fm->N/meta->freq_SA);

  // allocate space, then read the data
  if (getAll) ESL_ALLOC (fm->T, compressed_bytes );
  ESL_ALLOC (fm->BWT_mem, compressed_bytes + 31 ); // +31 for manual 16-byte alignment  ( typically only need +15, but this allows offset in memory, plus offset in case of <16 bytes of characters at the end)
     fm->BWT =   (uint8_t *) (((unsigned long int)fm->BWT_mem + 15) & (~0xf));   // align vector memory on 16-byte boundaries
  if (getAll) ESL_ALLOC (fm->SA, num_SA_samples * sizeof(uint32_t));
  ESL_ALLOC (fm->C, 1+meta->alph_size * sizeof(uint32_t));
  ESL_ALLOC (fm->occCnts_b,  num_freq_cnts_b *  (meta->alph_size ) * sizeof(uint16_t)); // every freq_cnt positions, store an array of ints
  ESL_ALLOC (fm->occCnts_sb,  num_freq_cnts_sb *  (meta->alph_size ) * sizeof(uint32_t)); // every freq_cnt positions, store an array of ints


  if(getAll && fread(fm->T, sizeof(uint8_t), compressed_bytes, fp) != compressed_bytes)
    esl_fatal( "%s: Error reading T in FM index.\n", __FILE__);
  if(fread(fm->BWT, sizeof(uint8_t), compressed_bytes, fp) != compressed_bytes)
    esl_fatal( "%s: Error reading BWT in FM index.\n", __FILE__);
  if(getAll && fread(fm->SA, sizeof(uint32_t), (size_t)num_SA_samples, fp) != (size_t)num_SA_samples)
    esl_fatal( "%s: Error reading SA in FM index.\n", __FILE__);

  if(fread(fm->occCnts_b, sizeof(uint16_t)*(meta->alph_size), (size_t)num_freq_cnts_b, fp) != (size_t)num_freq_cnts_b)
    esl_fatal( "%s: Error reading occCnts_b in FM index.\n", __FILE__);
  if(fread(fm->occCnts_sb, sizeof(uint32_t)*(meta->alph_size), (size_t)num_freq_cnts_sb, fp) != (size_t)num_freq_cnts_sb)
    esl_fatal( "%s: Error reading occCnts_sb in FM index.\n", __FILE__);


  //shortcut variables
  C          = fm->C;
  occCnts_b  = fm->occCnts_b;
  occCnts_sb = fm->occCnts_sb;


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

  return eslOK;

ERROR:
  freeFM(fm, getAll);
  esl_fatal("Error allocating memory in %s\n", "readFM");
  return eslFAIL;
}


/* Function:  readFMmeta()
 * Synopsis:  Read metadata from disk for the set of FM-indexes stored in a HMMER binary file
 *
 * Input: file pointer to binary file
 * Output: return filled meta struct
 */
int
readFMmeta( FILE *fp, FM_METADATA *meta)
{
  int status;
  int i;

  if( fread(&(meta->fwd_only),     sizeof(uint8_t),  1, fp) != 1 ||
      fread(&(meta->alph_type),    sizeof(uint8_t),  1, fp) != 1 ||
      fread(&(meta->alph_size),    sizeof(uint8_t),  1, fp) != 1 ||
      fread(&(meta->charBits),     sizeof(uint8_t),  1, fp) != 1 ||
      fread(&(meta->freq_SA),      sizeof(uint32_t), 1, fp) != 1 ||
      fread(&(meta->freq_cnt_sb),  sizeof(uint32_t), 1, fp) != 1 ||
      fread(&(meta->freq_cnt_b),   sizeof(uint32_t), 1, fp) != 1 ||
      fread(&(meta->SA_shift),     sizeof(uint8_t),  1, fp) != 1 ||
      fread(&(meta->cnt_shift_sb), sizeof(uint8_t),  1, fp) != 1 ||
      fread(&(meta->cnt_shift_b),  sizeof(uint8_t),  1, fp) != 1 ||
      fread(&(meta->block_count),  sizeof(uint16_t), 1, fp) != 1 ||
      fread(&(meta->seq_count),    sizeof(uint32_t), 1, fp) != 1
  )
    esl_fatal( "%s: Error reading meta data for FM index.\n", __FILE__);



  ESL_ALLOC (meta->seq_data,  meta->seq_count   * sizeof(FM_SEQDATA));
  if (meta->seq_data == NULL  )
    esl_fatal("unable to allocate memory to store FM meta data\n");

//  if(    fread(misc->meta->block_sizes,   sizeof(uint64_t), misc->meta->block_count, fp) != misc->meta->block_count )
//    esl_fatal( "%s: Error reading meta data for FM index.\n", __FILE__);


  for (i=0; i<meta->seq_count; i++) {
    if( fread(&(meta->seq_data[i].id),          sizeof(uint32_t),  1, fp) != 1 ||
        fread(&(meta->seq_data[i].start),       sizeof(uint32_t),  1, fp) != 1 ||
        fread(&(meta->seq_data[i].length),      sizeof(uint32_t),  1, fp) != 1 ||
        fread(&(meta->seq_data[i].offset),      sizeof(uint32_t),  1, fp) != 1 ||
        fread(&(meta->seq_data[i].name_length), sizeof(uint16_t),  1, fp) != 1
        )
      esl_fatal( "%s: Error reading meta data for FM index.\n", __FILE__);

    ESL_ALLOC (meta->seq_data[i].name, (1+meta->seq_data[i].name_length) * sizeof(char));

    if( fread(meta->seq_data[i].name,  sizeof(char), meta->seq_data[i].name_length+1  , fp) !=  meta->seq_data[i].name_length+1 )
      esl_fatal( "%s: Error reading meta data for FM index.\n", __FILE__);

  }

  return eslOK;

ERROR:

  if (meta->seq_data) {
    for (i=0; i<meta->seq_count; i++)
      free(meta->seq_data[i].name);
    free(meta->seq_data);
  }
  free(meta);

   esl_fatal("Error allocating memory in %s\n", "readFM");
   return eslFAIL;
}


/* Function:  computeSequenceOffset()
 * Synopsis:  search in the meta->seq_data array for the sequence id corresponding to the requested position
 *
 * Input: file pointer to binary file
 * Output: return filled meta struct
 */
//inline
uint32_t
computeSequenceOffset (FM_DATA *fms, FM_METADATA *meta, int block, int pos) {

  uint32_t lo = fms[block].seq_offset;
  uint32_t hi  = (block == meta->block_count-1 ? meta->seq_count : fms[block+1].seq_offset) - 1;
  uint32_t mid;

  /*  //linear scan
  for (mid=lo+1; i<=hi; i++) {
    if (meta->seq_data[i].offset > pos) // the position of interest belongs to the previous sequence
      break;
  }
  return i-1;
    */

  //binary search, first handling edge cases
  if (lo==hi)                           return lo;
  if (meta->seq_data[hi].offset <= pos) return hi;

  while (1) {
    mid = (lo + hi + 1) / 2;  /* round up */
    if      (meta->seq_data[mid].offset < pos) lo = mid; /* too far left  */
    else if (meta->seq_data[mid-1].offset > pos) hi = mid; /* too far right */
    else break;                 /* found it */
  }
  return mid-1;

}


/* hit_sorter(): qsort's pawn, below */
static int
hit_sorter(const void *a, const void *b)
{
  FM_HIT *h1 = (FM_HIT*)a;
  FM_HIT *h2 = (FM_HIT*)b;

  if      (h1->sortkey > h2->sortkey) return  1;
  else if (h1->sortkey < h2->sortkey) return -1;
  else {
    if  (h1->start > h2->start) return  1;
    else                        return -1;
  }
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
    int i;
  int count_only    = 0;

  FM_DATA *fmsf;
  FM_DATA *fmsb;
  FM_INTERVAL interval;
  FILE* fp_fm = NULL;
  FILE* fp = NULL;
  FILE* out = NULL;
  char *outname = NULL;

  ESL_GETOPTS     *go  = NULL;    /* command line processing                 */
  FM_MISC_VARS *misc;
  FM_METADATA *meta;

  //start timer
  t1 = times(&ts1);

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




  if((fp_fm = fopen(fname_fm, "rb")) == NULL)
    esl_fatal("Cannot open file `%s': ", fname_fm);


  ESL_ALLOC(misc, sizeof(FM_MISC_VARS));
  ESL_ALLOC(misc->meta, sizeof(FM_METADATA));
  meta = misc->meta;
  misc->occCallCnt = 0;

  readFMmeta( fp_fm, misc->meta);

  //read in FM-index bl.ocks
  ESL_ALLOC(fmsf, misc->meta->block_count * sizeof(FM_DATA) );
  for (i=0; i<meta->block_count; i++)
    readFM( fp_fm, fmsf+i, misc->meta, 1 );

  if (!meta->fwd_only) {
    ESL_ALLOC(fmsb, meta->block_count * sizeof(FM_DATA) );
    for (i=0; i<meta->block_count; i++)
      readFM( fp_fm, fmsb+i, misc->meta, 0 );
  }
  fclose(fp_fm);

  output_header(misc->meta, stdout, go, fname_fm, fname_queries);


  /* initialize a few global variables, then call initGlobals
   * to do architecture-specific initialization
   */
  misc->maskSA       =  meta->freq_SA - 1;
  misc->shiftSA      =  meta->SA_shift;
  fm_initMiscVars(misc);


  fm_createAlphabet(meta->alph_type, &alph, &inv_alph, &(meta->alph_size), NULL); // don't override charBits


  fp = fopen(fname_queries,"r");
  if (fp == NULL)
    esl_fatal("Unable to open file %s\n", fname_queries);

  ESL_ALLOC(line, FM_MAX_LINE * sizeof(char));

  hits_size = 200;
  ESL_ALLOC(hits, hits_size * sizeof(FM_HIT));


  while(fgets(line, FM_MAX_LINE, fp) ) {
    int qlen=0;
    while (line[qlen] != '\0' && line[qlen] != '\n')  qlen++;
    if (line[qlen] == '\n')  line[qlen] = '\0';

    hit_num = 0;

    for (i=0; i<meta->block_count; i++) {

      getSARangeForward(fmsb+i, misc, line, inv_alph, &interval);// yes, use the backward fm to produce a forward search on the forward fm
      if (interval.lower>0 && interval.lower <= interval.upper) {
        int new_hit_num =  interval.upper - interval.lower + 1;
        hit_num += new_hit_num;
        if (!count_only) {
          if (hit_num > hits_size) {
            hits_size = 2*hit_num;
            ESL_RALLOC(hits, tmp, hits_size * sizeof(FM_HIT));
          }
          getFMHits(fmsb+i, misc, &interval, i, hit_num-new_hit_num, hits, fm_forward);
        }
      }


      /* find reverse hits, using backward search on the forward FM*/
      if (!meta->fwd_only) {
        getSARangeReverse(fmsf+i, misc, line, inv_alph, &interval);
        if (interval.lower>0 && interval.lower <= interval.upper) {
          int new_hit_num =  interval.upper - interval.lower + 1;
          hit_num += new_hit_num;
          if (!count_only) {
            if (hit_num > hits_size) {
              hits_size = 2*hit_num;
              ESL_RALLOC(hits, tmp, hits_size * sizeof(FM_HIT));
            }
            getFMHits(fmsf+i, misc, &interval, i, hit_num-new_hit_num, hits, fm_backward);
          }
        }
      }
    }


    if (hit_num > 0) {
      hit_cnt++;

      if (count_only) {
        hit_indiv_cnt += hit_num;
      } else {
        //printf ("HIT:  %s (%d)\n", line, hit_num);
        //printf ("HIT:  %s  (%d, %d)\n", line, interval.lower, interval.upper );
        //for each hit, identify the sequence id and position within that sequence
        for (i = 0; i< hit_num; i++) {
          int block = hits[i].block;
          int seq_offset = computeSequenceOffset( fmsf, meta, block, hits[i].start);
          int pos =  ( hits[i].start - meta->seq_data[ seq_offset ].offset) + meta->seq_data[ seq_offset ].start - 1;
          //reuse hit variables.  Now "block" has the index into the matching sequence (in meta), and "start" has the pos within that sequence
          hits[i].block   = seq_offset;
          hits[i].start   = pos;
          hits[i].sortkey = meta->seq_data[ seq_offset ].id;
        }

        //now sort according the the sequence_id corresponding to that seq_offset
        qsort(hits, hit_num, sizeof(FM_HIT), hit_sorter);

        hit_indiv_cnt++;
        for (i = 1; i< hit_num; i++) {
          if (meta->seq_data[ hits[i].block ].id != meta->seq_data[ hits[i-1].block ].id ||
              hits[i].start != hits[i-1].start )
          {
            //printf ( "\t%10s (%d)\n",meta->seq_data[ hits[i].block ].name, hits[i].start);
            hit_indiv_cnt++;
          }
        }
        //printf ("\n");
      }
    } else {
      miss_cnt++;
    }


  }

  for (i=0; i<meta->block_count; i++) {
    freeFM( fmsb+i, 1 );
    if (!meta->fwd_only)
      freeFM( fmsf+i, 0 );
  }

  for (i=0; i<meta->seq_count; i++)
    free (meta->seq_data[i].name);

  free (meta->seq_data);
  free (meta);
  free (hits);
  free (line);

  fclose(fp);
  fm_destroyMiscVars(misc);
  free(misc);

  // compute and print the elapsed time in millisec
  t2 = times(&ts2);
  {
    double clk_ticks = sysconf(_SC_CLK_TCK);
    double elapsedTime = (t2-t1)/clk_ticks;
    double throughput = misc->occCallCnt/elapsedTime;

    fprintf (stderr, "hit: %-10d  (%d)\n", hit_cnt, hit_indiv_cnt);
    fprintf (stderr, "miss:%-10d\n", miss_cnt);
    fprintf (stderr, "run time:  %.2f seconds\n", elapsedTime);
    fprintf (stderr, "occ calls: %12s\n", commaprint(misc->occCallCnt));
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
