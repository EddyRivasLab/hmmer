/* Architecture-independent functions used in FM-index range compuations
 * and list management.
 *
 * Contents:
 *   1. List management
 *   2. Interval / range computation
 *   3. Functions related to the original sequence
 *   4. FM data initialization, configuration, and reading from file
 *
 *  # 5. Unit tests.
 *  # 6. Test driver.
 *  # 7. Copyright and license.
 *
 */
#include "hmmer.h"




/* Function:  fm_initSeeds()
 *
 * Synopsis:  initialize the object used to store a list of seed diagonals
 *
 * Returns:   eslEMEM in event of allocation failure, otherwise eslOK
 */
int
fm_initSeeds (FM_DIAGLIST *list) {
  int status;
  list->size = 1000;
  ESL_ALLOC(list->diags, list->size * sizeof(FM_DIAG));
  list->count = 0;

  return eslOK;

ERROR:
  return eslEMEM;

}

/* Function:  fm_newSeed()
 *
 * Synopsis:  return a pointer to the next seed element on the list,
 *            increasing the size of the list, if necessary.
 *
 * Returns:   NULL in event of allocation failure, otherwise pointer to
 *            the next seed diagonal
 */
FM_DIAG *
fm_newSeed (FM_DIAGLIST *list) {
  int status;

  if (list->count == list->size) {
    list->size *= 4;
    ESL_REALLOC(list->diags, list->size * sizeof(FM_DIAG));
  }
  list->count++;
  return list->diags + (list->count - 1);

ERROR:
  return NULL;
}

/* Function:  fm_initWindows()
 *
 * Synopsis:  initialize the object used to store a list of windows around seeds
 *
 * Returns:   eslEMEM in event of allocation failure, otherwise eslOK
 */
int
fm_initWindows (P7_MSV_WINDOWLIST *list) {
  int status;
  list->size = 1000;
  ESL_ALLOC(list->windows, list->size * sizeof(P7_MSV_WINDOW));
  list->count = 0;

  return eslOK;

ERROR:
  return eslEMEM;

}

/* Function:  fm_newWindow()
 *
 * Synopsis:  return a pointer to the next window element on the list,
 *            increasing the size of the list, if necessary.
 *
 * Purpose:   accepts <id>, <pos>, <fm_pos>, <k>, <length>, <score>,
 *            and <complementarity>, assigns those to the next window
 *            element, then returns it.
 *
 * Returns:   NULL in event of allocation failure, otherwise pointer to
 *            the next seed diagonal
 */
P7_MSV_WINDOW *
fm_newWindow (P7_MSV_WINDOWLIST *list, uint32_t id, uint32_t pos, uint32_t fm_pos, uint16_t k, uint32_t length, float score, uint8_t complementarity) {
  int status;
  P7_MSV_WINDOW *window;

  if (list->count == list->size) {
    list->size *= 4;
    ESL_REALLOC(list->windows, list->size * sizeof(P7_MSV_WINDOW));
  }
  window = list->windows + list->count;

  window->id      = id;
  window->n       = pos;
  window->fm_n    = fm_pos;
  window->k       = k;
  window->length  = length;
  window->score   = score;
  window->complementarity  = complementarity;

  list->count++;

  return window;

ERROR:
  return NULL;
}



/*********************************************************************
 *# 2. Interval / range computation
 *********************************************************************/

/* Function:  fm_updateIntervalForward()
 *
 * Synopsis:  return a pointer to the next window element on the list,
 *            increasing the size of the list, if necessary.
 *
 * Purpose:   Implement Algorithm 4 (i371) of Simpson (Bioinformatics 2010)
 *
 * Returns:   eslOK
 */
int
fm_updateIntervalForward( const FM_DATA *fm, FM_CFG *cfg, char c, FM_INTERVAL *interval_bk, FM_INTERVAL *interval_f) {
  uint32_t occLT_l, occLT_u, occ_l, occ_u;

  fm_getOccCountLT (fm, cfg, interval_bk->lower - 1, c, &occ_l, &occLT_l);
  fm_getOccCountLT (fm, cfg, interval_bk->upper,     c, &occ_u, &occLT_u);

  interval_f->lower += (occLT_u - occLT_l);
  interval_f->upper = interval_f->lower + (occ_u - occ_l) - 1;

  interval_bk->lower = abs(fm->C[(int)c]) + occ_l;
  interval_bk->upper = abs(fm->C[(int)c]) + occ_u - 1;

  return eslOK;
}

/* Function:  getSARangeForward()
 * Synopsis:  For a given query sequence, find its interval in the FM-index, using forward search
 * Purpose:   Implement Algorithm 4 (i371) of Simpson (Bioinformatics 2010). It's the forward
 *            search on a bi-directional BWT, as described by Lam 2009.
 *            All the meat is in the method of counting characters - bwt_getOccCount, which
 *            depends on compilation choices.
 *
 *            Note: it seems odd, but is correct, that the fm-index passed in to this function
 *            is the backward index corresponding to the forward index into which I want to
 *            do a forward search
 */
int
fm_getSARangeForward( const FM_DATA *fm, FM_CFG *cfg, char *query, char *inv_alph, FM_INTERVAL *interval)
{

  int i=0;
  FM_INTERVAL interval_bk;

  uint8_t c = inv_alph[(int)query[0]];
  interval->lower  = interval_bk.lower = abs(fm->C[c]);
  interval->upper  = interval_bk.upper = abs(fm->C[c+1])-1;

  while (interval_bk.lower>0 && interval_bk.lower <= interval_bk.upper) {
    c = query[++i];
    if (c == '\0')  // end of query - the current range defines the hits
      break;

    c = inv_alph[c];
    fm_updateIntervalForward( fm, cfg, c, &interval_bk, interval);
    cfg->occCallCnt+=2;
  }

  return eslOK;
}


int
fm_updateIntervalReverse( const FM_DATA *fm, FM_CFG *cfg, char c, FM_INTERVAL *interval) {
  int count1, count2;
  //TODO: counting in these calls will often overlap
    // - might get acceleration by merging to a single redundancy-avoiding call
  count1 = fm_getOccCount (fm, cfg, interval->lower-1, c);
  count2 = fm_getOccCount (fm, cfg, interval->upper, c);

  interval->lower = abs(fm->C[(int)c]) + count1;
  interval->upper = abs(fm->C[(int)c]) + count2 - 1;

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
fm_getSARangeReverse( const FM_DATA *fm, FM_CFG *cfg, char *query, char *inv_alph, FM_INTERVAL *interval)
{


  int i=0;

  char c = inv_alph[(int)query[0]];
  interval->lower  = abs(fm->C[(int)c]);
  interval->upper  = abs(fm->C[(int)c+1])-1;

//  fprintf (stderr, "%d : %12d..%12d\n", c, interval->lower, interval->upper );

  while (interval->lower>0 && interval->lower <= interval->upper) {
    c = query[++i];
    if (c == '\0')  // end of query - the current range defines the hits
      break;

    c = inv_alph[(int)c];

    fm_updateIntervalReverse(fm, cfg, c, interval);

//    fprintf (stderr, "%d : %12d..%12d\n", c, interval->lower, interval->upper );

    cfg->occCallCnt+=2;
  }

  return eslOK;
}



/*********************************************************************
 *# 3. Functions related to the original sequence
 *********************************************************************/


/* Function:  getChar()
 * Synopsis:  Find the character c residing at a given position in the BWT.
 * Purpose:   The returned char is used by getFMHits(), to seed a call to
 *            bwt_getOccCount().
 */
//#ifndef FMDEBUG
//inline
//#endif
uint8_t
fm_getChar(uint8_t alph_type, int j, const uint8_t *B )
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


/* Function:  getChar()
 * Synopsis:  Find the character c residing at a given position in the BWT.
 * Purpose:   The returned char is used by getFMHits(), to seed a call to
 *            bwt_getOccCount().
 */
//#ifndef FMDEBUG
//inline
//#endif
int
fm_convertRange2DSQ(FM_METADATA *meta, int id, int first, int length, const uint8_t *B, ESL_SQ *sq )
{
  int i;
  uint8_t c;

  esl_sq_GrowTo(sq, length);
  sq->n = length;

  if (meta->alph_type == fm_DNA || meta->alph_type == fm_RNA) {
    /*
     *  B[j>>2] is the byte of B in which j is found (j/4)
     *
     *  Let j' be the final two bits of j (j&0x2)
     *  The char bits are the two starting at position 2*j'.
     *  Without branching, grab them by shifting B[j>>2] right 6-2*j' bits,
     *  then masking to keep the final two bits
     */
    for (i = first; i<= first+length-1; i++)
      sq->dsq[i-first+1] = (B[i>>2] >> ( 0x6 - ((i&0x3)<<1) ) & 0x3);

    sq->dsq[length+1] = eslDSQ_SENTINEL;

  } else if (meta->alph_type == fm_DNA_full || meta->alph_type == fm_RNA_full ) {
    for (i = first; i<= first+length-1; i++) {
      c = (B[i>>1] >> (((i&0x1)^0x1)<<2) ) & 0xf;  //unpack the char: shift 4 bits right if it's odd, then mask off left bits in any case
      sq->dsq[i-first+1] = c + (c < 4 ? 0 : 1);
    }
    sq->dsq[length+1] = eslDSQ_SENTINEL;
  } else {
    esl_fatal("Invalid alphabet type\n");
  }

  return eslOK;
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
fm_computeSequenceOffset (const FM_DATA *fms, FM_METADATA *meta, int block, int pos)
{

  uint32_t lo = fms[block].seq_offset;
  uint32_t hi  = lo + fms[block].seq_cnt - 1;
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


/* Function:  FM_getOriginalPosition()
 * Synopsis:  find
 * Purpose:   Given:
 *            fms       - an array of FM-indexes
 *            meta      - the fm metadata
 *            fm_id     - index of the fm-index in which a hit is sought
 *            length    - length of the hit in question
 *            direction - direction of the hit in question
 *            fm_pos    - position in the fm-index
 *
 *            Returns
 *            *segment_id - index of the sequence segment captured in the FM-index
 *            *seg_pos    - position in the original sequence, as compressed in the FM binary data structure
 */
int
fm_getOriginalPosition (const FM_DATA *fms, FM_METADATA *meta, int fm_id, int length, int direction, uint32_t fm_pos,
                  uint32_t *segment_id, uint32_t *seg_pos
) {
  *segment_id = fm_computeSequenceOffset( fms, meta, fm_id, fm_pos);
  *seg_pos    =  ( fm_pos - meta->seq_data[ *segment_id ].offset) + meta->seq_data[ *segment_id ].start - 1;

  //verify that the hit doesn't extend beyond the bounds of the target sequence
  if (direction == fm_forward) {
    if (*seg_pos + length > meta->seq_data[ *segment_id ].start + meta->seq_data[ *segment_id ].length ) {
      *segment_id  = *seg_pos  = -1;  // goes into the next sequence, so it should be ignored
      return eslOK;
    }
  } else { //backward
    if ((int)*seg_pos - length + 1 < 0 ) {
      *segment_id  = *seg_pos  = -1; // goes into the previous sequence, so it should be ignored
      return eslOK;
    }
  }


  return eslOK;
}

/*********************************************************************
 *# 4. FM data initialization, configuration, and reading from file
 *********************************************************************/

int
fm_initConfigGeneric( FM_CFG *cfg, ESL_GETOPTS *go ) {

  cfg->maskSA       =  cfg->meta->freq_SA - 1;
  cfg->shiftSA      =  cfg->meta->SA_shift;

  cfg->msv_length      = (go ? esl_opt_GetInteger(go, "--fm_msv_length") : -1);
  cfg->max_depth       = (go ? esl_opt_GetInteger(go, "--fm_max_depth") :  -1);
  cfg->neg_len_limit   = (go ? esl_opt_GetInteger(go, "--fm_max_neg_len") : -1);
  cfg->consec_pos_req  = (go ? esl_opt_GetInteger(go, "--fm_req_pos") : -1);
  cfg->score_ratio_req = (go ? esl_opt_GetReal(go, "--fm_sc_ratio") : -1.0);
  cfg->max_scthreshFM  = (go ? esl_opt_GetReal(go, "--fm_max_scthresh") : -1.0);


/*

  //bounding cutoffs
  cfg->max_scthreshFM  = 10.5;

  cfg->max_depth       = 18;
  cfg->neg_len_limit   = 4;
  cfg->consec_pos_req  = 4;
  cfg->score_ratio_req = 0.40;
  cfg->msv_length      = 50;

*/
  return eslOK;
}

/* Function:  freeFM()
 * Synopsis:  release the memory required to store an individual FM-index
 */
void
fm_freeFM ( FM_DATA *fm, int isMainFM)
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
fm_readFM( FM_DATA *fm, FM_METADATA *meta, int getAll )
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


  if(fread(&(fm->N), sizeof(fm->N), 1, meta->fp) !=  1)
    esl_fatal( "%s: Error reading block_length in FM index.\n", __FILE__);
  if(fread(&(fm->term_loc), sizeof(fm->term_loc), 1, meta->fp) !=  1)
    esl_fatal( "%s: Error reading terminal location in FM index.\n", __FILE__);
  if(fread(&(fm->seq_offset), sizeof(fm->seq_offset), 1, meta->fp) !=  1)
    esl_fatal( "%s: Error reading seq_offset in FM index.\n", __FILE__);
  if(fread(&(fm->overlap), sizeof(fm->overlap), 1, meta->fp) !=  1)
    esl_fatal( "%s: Error reading overlap in FM index.\n", __FILE__);
  if(fread(&(fm->seq_cnt), sizeof(fm->seq_cnt), 1, meta->fp) !=  1)
    esl_fatal( "%s: Error reading seq_cnt in FM index.\n", __FILE__);


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


  if(getAll && fread(fm->T, sizeof(uint8_t), compressed_bytes, meta->fp) != compressed_bytes)
    esl_fatal( "%s: Error reading T in FM index.\n", __FILE__);
  if(fread(fm->BWT, sizeof(uint8_t), compressed_bytes, meta->fp) != compressed_bytes)
    esl_fatal( "%s: Error reading BWT in FM index.\n", __FILE__);
  if(getAll && fread(fm->SA, sizeof(uint32_t), (size_t)num_SA_samples, meta->fp) != (size_t)num_SA_samples)
    esl_fatal( "%s: Error reading SA in FM index.\n", __FILE__);

  if(fread(fm->occCnts_b, sizeof(uint16_t)*(meta->alph_size), (size_t)num_freq_cnts_b, meta->fp) != (size_t)num_freq_cnts_b)
    esl_fatal( "%s: Error reading occCnts_b in FM index.\n", __FILE__);
  if(fread(fm->occCnts_sb, sizeof(uint32_t)*(meta->alph_size), (size_t)num_freq_cnts_sb, meta->fp) != (size_t)num_freq_cnts_sb)
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
  fm_freeFM(fm, getAll);
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
fm_readFMmeta( FM_METADATA *meta)
{
  int status;
  int i;

  if( fread(&(meta->fwd_only),     sizeof(meta->fwd_only),     1, meta->fp) != 1 ||
      fread(&(meta->alph_type),    sizeof(meta->alph_type),    1, meta->fp) != 1 ||
      fread(&(meta->alph_size),    sizeof(meta->alph_size),    1, meta->fp) != 1 ||
      fread(&(meta->charBits),     sizeof(meta->charBits),     1, meta->fp) != 1 ||
      fread(&(meta->freq_SA),      sizeof(meta->freq_SA),      1, meta->fp) != 1 ||
      fread(&(meta->freq_cnt_sb),  sizeof(meta->freq_cnt_sb),  1, meta->fp) != 1 ||
      fread(&(meta->freq_cnt_b),   sizeof(meta->freq_cnt_b),   1, meta->fp) != 1 ||
      fread(&(meta->SA_shift),     sizeof(meta->SA_shift),     1, meta->fp) != 1 ||
      fread(&(meta->cnt_shift_sb), sizeof(meta->cnt_shift_sb), 1, meta->fp) != 1 ||
      fread(&(meta->cnt_shift_b),  sizeof(meta->cnt_shift_b),  1, meta->fp) != 1 ||
      fread(&(meta->block_count),  sizeof(meta->block_count),  1, meta->fp) != 1 ||
      fread(&(meta->seq_count),    sizeof(meta->seq_count),    1, meta->fp) != 1 ||
      fread(&(meta->char_count),   sizeof(meta->char_count),   1, meta->fp) != 1
  )
    esl_fatal( "%s: Error reading meta data for FM index.\n", __FILE__);


  ESL_ALLOC (meta->seq_data,  meta->seq_count   * sizeof(FM_SEQDATA));
  if (meta->seq_data == NULL  )
    esl_fatal("unable to allocate memory to store FM meta data\n");


  for (i=0; i<meta->seq_count; i++) {
    if( fread(&(meta->seq_data[i].id),           sizeof(meta->seq_data[i].id),           1, meta->fp) != 1 ||
        fread(&(meta->seq_data[i].start),        sizeof(meta->seq_data[i].start),        1, meta->fp) != 1 ||
        fread(&(meta->seq_data[i].length),       sizeof(meta->seq_data[i].length),       1, meta->fp) != 1 ||
        fread(&(meta->seq_data[i].offset),       sizeof(meta->seq_data[i].offset),       1, meta->fp) != 1 ||
        fread(&(meta->seq_data[i].name_length),  sizeof(meta->seq_data[i].name_length),  1, meta->fp) != 1 ||
        fread(&(meta->seq_data[i].acc_length),   sizeof(meta->seq_data[i].acc_length),   1, meta->fp) != 1 ||
        fread(&(meta->seq_data[i].source_length),sizeof(meta->seq_data[i].source_length),1, meta->fp) != 1 ||
        fread(&(meta->seq_data[i].desc_length),  sizeof(meta->seq_data[i].desc_length),  1, meta->fp) != 1
        )
      esl_fatal( "%s: Error reading meta data for FM index.\n", __FILE__);

    ESL_ALLOC (meta->seq_data[i].name,  (1+meta->seq_data[i].name_length)   * sizeof(char));
    ESL_ALLOC (meta->seq_data[i].acc,   (1+meta->seq_data[i].acc_length)    * sizeof(char));
    ESL_ALLOC (meta->seq_data[i].source,(1+meta->seq_data[i].source_length) * sizeof(char));
    ESL_ALLOC (meta->seq_data[i].desc,  (1+meta->seq_data[i].desc_length)   * sizeof(char));

    if(
        fread(meta->seq_data[i].name,   sizeof(char), meta->seq_data[i].name_length+1   , meta->fp) !=  meta->seq_data[i].name_length+1  ||
        fread(meta->seq_data[i].acc,    sizeof(char), meta->seq_data[i].acc_length+1    , meta->fp) !=  meta->seq_data[i].acc_length+1  ||
        fread(meta->seq_data[i].source, sizeof(char), meta->seq_data[i].source_length+1 , meta->fp) !=  meta->seq_data[i].source_length+1  ||
        fread(meta->seq_data[i].desc,   sizeof(char), meta->seq_data[i].desc_length+1   , meta->fp) !=  meta->seq_data[i].desc_length+1 )
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





/* Function:  fm_configAlloc()
 * Synopsis:  Allocate a 16-byte-aligned <FM_CFG> model object, and it's FM_METADATA
 *
 */
int
fm_configAlloc(void **mem, FM_CFG **cfg)
{
  int status;

  if (mem == NULL || cfg == NULL)
    esl_fatal("null pointer when allocating FM configuration\n");

  *mem = NULL;
  *cfg = NULL;

  ESL_ALLOC(*mem, sizeof(FM_CFG)+ 15 );
    *cfg =   (FM_CFG *) (((unsigned long int)(*mem) + 15) & (~0xf));   /* align vector memory on 16-byte boundaries */

  ESL_ALLOC((*cfg)->meta, sizeof(FM_METADATA));

  return eslOK;

ERROR:
  if (*cfg != NULL)
   if ((*cfg)->meta != NULL) free ((*cfg)->meta);

  if (*mem != NULL)
    free (*mem);

  return eslEMEM;
}


/************************************************************
 * @LICENSE@
 *
 * SVN $Id: fm_general.c 3784 2011-12-07 21:51:25Z wheelert $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/src/hmmer/trunk/src/fm_general.c $
 ************************************************************/

