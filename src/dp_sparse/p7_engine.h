/* The comparison engine.
 * Given a profile and a sequence; calculates scores, anchors, envelopes, and alignments.
 */
#ifndef p7ENGINE_INCLUDED
#define p7ENGINE_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"

#include "base/p7_anchors.h"
#include "base/p7_anchorhash.h"
#include "base/p7_bg.h"
#include "base/p7_envelopes.h"
#include "base/p7_trace.h"

#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/p7_filtermx.h"
#include "dp_vector/p7_oprofile.h"

#include "dp_sparse/p7_sparsemx.h"

#include "search/p7_mpas.h"


/* P7_ENGINE_PARAMS 
 * Configuration/control settings for the Engine.
 */
typedef struct p7_engine_params_s {
  uint32_t rng_seed;           // random number generator seed. >0 for specific, reproducible seed; 0=random seed. Default = 42.
  int      rng_reproducible;   // TRUE to reseed RNG at every comparison, enabling reproducible results.
  int      sparsify_ramlimit;  // Memory redline for checkpointed decoding, in MB. Default = p7_SPARSIFY_RAMLIMIT [p7_config.h]
  float    sparsify_thresh;    // (i,k) supercell included in sparsemask if pp>this probability, 0<=x<1. Default = p7_SPARSIFY_THRESH [p7_config.h]
  int      do_biasfilter;      // TRUE to use ad hoc "bias filter" after MSV/SSV step

  P7_MPAS_PARAMS *mpas_params;  // optional config/control parameters for MPAS algorithm; or NULL for defaults

} P7_ENGINE_PARAMS;


/* P7_ENGINE_STATS
 * Statistics collection for the Engine.
 */
typedef struct p7_engine_stats_s {
  int n_past_msv;
  int n_past_bias;
  int n_ran_vit;
  int n_past_vit;
  int n_past_fwd;

} P7_ENGINE_STATS;

/* P7_ENGINE
 * The Engine.
 */
typedef struct p7_engine_s {
  ESL_RANDOMNESS *rng;    // Random number generator; used for MPAS sampling

  P7_FILTERMX    *fx;     // one-row vectorized DP for MSV, Vit filters. O(M) mem
  P7_CHECKPTMX   *cx;     // Checkpointed vector local F/B/D matrix.     O(M \sqrt L) mem.  (p7_SPARSIFY_RAMLIMIT = 128M)
  P7_SPARSEMASK  *sm;     // Sparse mask.                                O(L) mem. 

  int       used_main;    // if TRUE, main engine was used and these structures need to be Reuse()'d:
  P7_SPARSEMX    *sxf;    // Sparse Forward matrix                       O(L)
  P7_SPARSEMX    *sxd;    // Sparse Decoding (also briefly Backward)     O(L)
  P7_SPARSEMX    *asf;    // ASC Sparse Forward mx
  P7_SPARSEMX    *asb;    // ASC Sparse Backward mx  (could be optimized away, but current impl can't decode while overwriting <asb>)
  P7_SPARSEMX    *asd;    // ASC Sparse Decoding mx
  P7_ANCHORS     *vanch;  // Initial anchor set implied by the Viterbi parse
  P7_ANCHORS     *anch;   // Anchor set optimized by MPAS
  P7_ANCHORHASH  *ahash;  // used in MPAS algorithm, for checking for anchorsets already seen
  P7_ENVELOPES   *env;    // inferred envelope
  P7_TRACE       *tr;     // inferred alignment
  
  float          *wrkM;   // O(M) temporary workspace. (stochastic trace; sparse null2)
  float          *wrkKp;  // O(Kp) temporary workspace (sparse null2)

  float           nullsc; // null raw score
  float           biassc; // ad hoc "bias filter" score, acts as modified null
  float           mfsc;   // MSV raw score
  float           vfsc;   // ViterbiFilter raw score, complete sequence
  float           ffsc;   // ForwardFilter raw score, complete sequence
  float           vsc;    // sparse Viterbi score
  float           fsc;    // sparse Forward score
  float           asc_f;  // ASC Forward score, s^A_f

  float           F1;
  float           F2;
  float           F3;

  P7_ENGINE_PARAMS *params; // config/control parameters for the Engine
  P7_ENGINE_STATS  *stats;  // optional stats collection for the Engine, or NULL

} P7_ENGINE;

extern P7_ENGINE_PARAMS *p7_engine_params_Create (P7_MPAS_PARAMS *mpas_params);
extern void              p7_engine_params_Destroy(P7_ENGINE_PARAMS *prm);

extern P7_ENGINE_STATS  *p7_engine_stats_Create(void);
extern void              p7_engine_stats_Destroy(P7_ENGINE_STATS *prm);

extern P7_ENGINE *p7_engine_Create (const ESL_ALPHABET *abc, P7_ENGINE_PARAMS *prm, P7_ENGINE_STATS *stats, int M_hint, int L_hint);
extern int        p7_engine_Reuse  (P7_ENGINE *eng);
extern void       p7_engine_Destroy(P7_ENGINE *eng);

extern int p7_engine_Overthruster(P7_ENGINE *eng, ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_BG *bg);
extern int p7_engine_Main        (P7_ENGINE *eng, ESL_DSQ *dsq, int L, P7_PROFILE  *gm);

#endif /*p7ENGINE_INCLUDED*/

