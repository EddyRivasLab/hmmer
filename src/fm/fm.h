#ifndef P7_FM_INCLUDED
#define P7_FM_INCLUDED

#include "p7_config.h"
#include <stdio.h"

#define FM_MAX_LINE 256

/* Structure the 2D occ array into a single array.  "type" is either b or sb.
 * Note that one extra count value is required by RLE, one 4-byte int for
 * each superblock count vector, and one 2-byte short for each block count
 * vector. This is small overhead, even for a small alphabet like dna.
 */
#define FM_OCC_CNT( type, i, c)  ( occCnts_##type[(meta->alph_size)*(i) + (c)])

enum fm_alphabettypes_e {
  fm_DNA        = 0,  //acgt,  2 bit
  fm_DNA_full   = 1,  //includes ambiguity codes, 4 bit
  fm_RNA        = 2,  //acgu,  2 bit
  fm_RNA_full   = 3,  //includes ambiguity codes, 4 bit
  fm_AMINO      = 4,  // 5 bit
};

enum fm_direction_e {
  fm_forward    = 0,
  fm_backward   = 1,
};

enum fm_complementarity_e {
  fm_nocomplement    = 0,
  fm_complement   = 1,
};

typedef struct fm_interval_s {
  int   lower;
  int   upper;
} FM_INTERVAL;

typedef struct fm_hit_s {
  uint32_t  start;
  uint32_t  block;
  int       direction;
  int       length;
  int       sortkey;
} FM_HIT;


typedef struct fm_seqdata_s {
  uint32_t id;
  uint32_t start;
  uint32_t length;
  uint32_t offset;
  uint16_t name_length;
  uint16_t source_length;
  uint16_t acc_length;
  uint16_t desc_length;
  char     *name;
  char     *source;
  char     *acc;
  char     *desc;
} FM_SEQDATA;


typedef struct fm_metadata_s {
  uint8_t  fwd_only;
  uint8_t  alph_type;
  uint8_t  alph_size;
  uint8_t  charBits;
  uint32_t freq_SA; //frequency with which SA is sampled
  uint32_t freq_cnt_sb; //frequency with which full cumulative counts are captured
  uint32_t freq_cnt_b; //frequency with which intermittent counts are captured
  uint8_t  SA_shift;
  uint8_t  cnt_shift_sb;
  uint8_t  cnt_shift_b;
  uint16_t block_count;
  uint32_t seq_count;
  uint64_t char_count; //total count of characters including those in and out of the alphabet
  char     *alph;
  char     *inv_alph;
  FILE       *fp;
  FM_SEQDATA *seq_data;
} FM_METADATA;



typedef struct fm_data_s {
  uint32_t N; //length of text
  uint32_t term_loc; // location in the BWT at which the '$' char is found (replaced in the sequence with 'a')
  uint32_t seq_offset;
  uint32_t overlap; // number of bases at the beginning that overlap the FM-index for the preceding block
  uint16_t seq_cnt;
  uint8_t  *T;  //text corresponding to the BWT
  uint8_t  *BWT_mem;
  uint8_t  *BWT;
  uint32_t *SA; // sampled suffix array
  int32_t  *C; //the first position of each letter of the alphabet if all of T is sorted.  (signed, as I use that to keep tract of presence/absence)
  uint32_t *occCnts_sb;
  uint16_t *occCnts_b;
} FM_DATA;

typedef struct fm_dp_pair_s {
  uint16_t    pos;  // position of the diagonal in the model.
  float       score;
  float       max_score;
  uint8_t     max_score_len; // how long was the diagonal when the maximum observed score was seen?
  uint8_t     consec_pos;
  uint8_t     max_consec_pos;
  uint8_t     model_direction;
  uint8_t     complementarity;
} FM_DP_PAIR;


typedef struct fm_diag_s {
  uint32_t    n;  //position of the database sequence at which the diagonal starts
  double       sortkey;
  uint16_t    k;  //position of the model at which the diagonal starts
  uint16_t    length;
  uint8_t     complementarity;
} FM_DIAG;

typedef struct fm_diaglist_s {
  FM_DIAG   *diags;
  int       count;
  int       size;
} FM_DIAGLIST;

#endif /*P7_FM_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

