#include "hmmer.h"
//#include "easel.h"
//#include <string.h>



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
getSARangeForward(FM_DATA *fm, FM_CFG *cfg, char *query, char *inv_alph, FM_INTERVAL *interval)
{

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

    fm_getOccCountLT (fm, cfg, lower_b-1, c, &occ_l, &occLT_l);
    fm_getOccCountLT (fm, cfg, upper_b,   c, &occ_u, &occLT_u);

    interval->lower += (occLT_u - occLT_l);
    interval->upper = interval->lower + (occ_u - occ_l) - 1;

    lower_b = abs(fm->C[c]) + occ_l;
    upper_b = abs(fm->C[c]) + occ_u - 1;

    cfg->occCallCnt+=2;
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
getSARangeReverse( FM_DATA *fm, FM_CFG *cfg, char *query, char *inv_alph, FM_INTERVAL *interval)
{

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
    count1 = fm_getOccCount (fm, cfg, interval->lower-1, c);
    count2 = fm_getOccCount (fm, cfg, interval->upper, c);

    interval->lower = abs(fm->C[c]) + count1;
    interval->upper = abs(fm->C[c]) + count2 - 1;

    cfg->occCallCnt+=2;
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
getChar(uint8_t alph_type, int j, const uint8_t *B )
{
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


/* Function:  freeFM()
 * Synopsis:  release the memory required to store an individual FM-index
 */
void
freeFM ( FM_DATA *fm, int isMainFM)
{

  if (fm->BWT_mem)      free (fm->BWT_mem);
  if (fm->C)            free (fm->C);
  if (fm->occCnts_b)    free (fm->occCnts_b);
  if (fm->occCnts_sb)   free (fm->occCnts_sb);

  if (isMainFM && fm->T)  free (fm->T);
  if (isMainFM && fm->SA) free (fm->SA);
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

  if( fread(&(meta->fwd_only),     sizeof(meta->fwd_only),     1, fp) != 1 ||
      fread(&(meta->alph_type),    sizeof(meta->alph_type),    1, fp) != 1 ||
      fread(&(meta->alph_size),    sizeof(meta->alph_size),    1, fp) != 1 ||
      fread(&(meta->charBits),     sizeof(meta->charBits),     1, fp) != 1 ||
      fread(&(meta->freq_SA),      sizeof(meta->freq_SA),      1, fp) != 1 ||
      fread(&(meta->freq_cnt_sb),  sizeof(meta->freq_cnt_sb),  1, fp) != 1 ||
      fread(&(meta->freq_cnt_b),   sizeof(meta->freq_cnt_b),   1, fp) != 1 ||
      fread(&(meta->SA_shift),     sizeof(meta->SA_shift),     1, fp) != 1 ||
      fread(&(meta->cnt_shift_sb), sizeof(meta->cnt_shift_sb), 1, fp) != 1 ||
      fread(&(meta->cnt_shift_b),  sizeof(meta->cnt_shift_b),  1, fp) != 1 ||
      fread(&(meta->block_count),  sizeof(meta->block_count),  1, fp) != 1 ||
      fread(&(meta->seq_count),    sizeof(meta->seq_count),    1, fp) != 1
  )
    esl_fatal( "%s: Error reading meta data for FM index.\n", __FILE__);


  ESL_ALLOC (meta->seq_data,  meta->seq_count   * sizeof(FM_SEQDATA));
  if (meta->seq_data == NULL  )
    esl_fatal("unable to allocate memory to store FM meta data\n");


  for (i=0; i<meta->seq_count; i++) {
    if( fread(&(meta->seq_data[i].id),          sizeof(meta->seq_data[i].id),          1, fp) != 1 ||
        fread(&(meta->seq_data[i].start),       sizeof(meta->seq_data[i].start),       1, fp) != 1 ||
        fread(&(meta->seq_data[i].length),      sizeof(meta->seq_data[i].length),      1, fp) != 1 ||
        fread(&(meta->seq_data[i].offset),      sizeof(meta->seq_data[i].offset),      1, fp) != 1 ||
        fread(&(meta->seq_data[i].name_length), sizeof(meta->seq_data[i].name_length), 1, fp) != 1
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
 * Synopsis:  Search in the meta->seq_data array for the sequence id corresponding to the
 *            requested position. The matching entry is the one with the largest index i
 *            such that seq_data[i].offset < pos
 *
 *
 * Input: file pointer to binary file
 * Output: return filled meta struct
 */
//inline
uint32_t
computeSequenceOffset (FM_DATA *fms, FM_METADATA *meta, int block, int pos)
{

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
    if      (meta->seq_data[mid].offset <= pos) lo = mid; /* too far left  */
    else if (meta->seq_data[mid-1].offset > pos) hi = mid; /* too far right */
    else return mid-1;                 /* found it */
  }


}


