/* H4_DOMAIN: information about one domain in a stored hit
 *
 * Under construction. 
 * NWO: combining base/p7_domain, base/p7_envelope, base/p7_alidisplay
 * and using new magic capsules.
 */


/* 
 * 
 * MEMORY-CRITICAL STRUCTURE. Do not add fields without review/consultation.
 */
typedef struct {         // currently 37B + len(aliz) + len(ppz)
  int32_t   i0, k0;      // anchor for this domain         (1..L, 1..M) (-1 if unused)
  int32_t   ia, ib;      // envelope/ali start/end on sequence   (1..L)
  int32_t   ka, kb;      // envelope/ali start/end on profile    (1..M) (1,M if glocal)
  int32_t   oa, ob;      // outer envelope start/end on sequence (1..L)
  float     envscore;    // envelope score in bits
  char     *aliz;        // alignment magic capsule, \0-term
  char     *ppz;         // ali accuracy magic capsule, \0-term
  uint8_t   flags;       // various flags; see definitions below.
} H4_DOMAIN;


#define h4D_ENVSC_APPROX (1<<0)
#define h4D_IS_GLOCAL    (1<<1)
#define h4D_IS_INCLUDED  (1<<2)
#define h4D_IS_REPORTED  (1<<3)

