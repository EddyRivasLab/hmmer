/* H4_MODE: "algorithm-dependent" parameters of a HMMER4 profile
 * 
 * See also:
 *   h4_profile : the profile itself; "model-dependent" parameters and metadata.
 */
#ifndef h4MODE_INCLUDED
#define h4MODE_INCLUDED

#include <h4_config.h>

/* Constants defining fixed sizes.
 * (These are to avoid cryptic bare numbers like '5' in the code.) 
 */
#define h4_NX      5     // number of 'special' states w/ params, ENJCB
#define h4_NXT     2     // each special state has a LOOP and MOVE transition.


/* H4_MODE
 * 
 * We need to reconfigure seqlength-dependent parameters of the
 * profile for each target sequence. We pull these changeable
 * parameters out into a separate structure, so we can treat
 * a H4_PROFILE as const in big (parallelized) searches. 
 */
typedef struct {
  float   xf[h4_NX][h4_NXT];    // probabilities [ENJCB][LOOP,MOVE] (i.e. in FB filter)
  float   xsc[h4_NX][h4_NXT];   // log2 bit scores (i.e. in reference DP)
  int16_t xw[h4_NX][h4_NXT];    // 16-bit scaled/offset log2 bit scores (i.e. in VF)

  float nullsc;                 // null1 score correction for length L
  int   L;                      // current configured target seq length                      unset=-1
  float nj;                     // expected # of J's: 0.0 = unihit; 1.0 = standard multihit. unset=-1.
  float pglocal;                // B->G probability
} H4_MODE;


/* Indices for five (h4_NX) special state types, each with two
 * (h4_NXT) transition parameters.
 */
#define h4_LOOP  0         //       xsc[][LOOP] is:     xsc[][MOVE] is:
#define h4_MOVE  1         //    -------------------    ---------------
#define h4_E     0         //         E -> J                E -> C
#define h4_N     1         //         N -> N                N -> B
#define h4_J     2         //         J -> J                J -> B
#define h4_C     3         //         C -> C                C -> T
#define h4_B     4         //         B -> L                B -> G


extern H4_MODE *h4_mode_Create      (void);
extern H4_MODE *h4_mode_Clone(const H4_MODE *mo);
extern int      h4_mode_Copy (const H4_MODE *mo, H4_MODE *m2);
extern int      h4_mode_SetCustom   (H4_MODE *mo, int L, float nj, float pglocal);
extern int      h4_mode_SetDefault  (H4_MODE *mo);
extern int      h4_mode_SetLocal    (H4_MODE *mo);
extern int      h4_mode_SetGlocal   (H4_MODE *mo);
extern int      h4_mode_SetUnihit   (H4_MODE *mo);
extern int      h4_mode_SetUnilocal (H4_MODE *mo);
extern int      h4_mode_SetUniglocal(H4_MODE *mo);
extern int      h4_mode_SetGlobal   (H4_MODE *mo);        
extern int      h4_mode_SetLength   (H4_MODE *mo, int L);
extern void     h4_mode_Destroy     (H4_MODE *mo);

extern int      h4_mode_Dump(FILE *fp, const H4_MODE *mo);
extern int      h4_mode_SameAsSSV(const H4_MODE *mo, H4_MODE **ret_xmo);
extern int      h4_mode_SameAsVF (const H4_MODE *mo, H4_MODE **ret_xmo);

#endif //h4MODE_INCLUDED
