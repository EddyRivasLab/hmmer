/* H4_HIT: info about a profile/seq comparison saved to top hits list
 * 
 * 
 */


/* 
 * 
 * MEMORY-CRITICAL STRUCTURE. Do not add fields without review/consultation.
 */
typedef struct {
  int32_t    profile_id;
  int64_t    seq_id;
  float      seq_score;
  float      logP;

  H4_DOMAIN *dom;
  int        D;
  float      D_exp;

  uint8_t    flags;
} H4_HIT;

#define h4H_IS_INCLUDED  (1<<0)
#define h4H_IS_REPORTED	 (1<<1)		
// NWO: will probably need to bring IS_NEW, IS_DROPPED, IS_DUPLICATE flags over??
