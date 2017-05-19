/* P7_OPROFILE: a search profile in vectorized form.
 * 
 * Independent of vector ISA. (Do not add any ISA-specific code.)
 * See notes in p7_oprofile.md.
 */
#ifndef p7OPROFILE_INCLUDED
#define p7OPROFILE_INCLUDED
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"

#include "base/general.h"
#include "base/p7_bg.h"
#include "base/p7_hmm.h"
#include "base/p7_profile.h"

#include "dp_vector/simdvec.h"

#define p7O_NXSTATES  4    // special states stored: ENJC
#define p7O_NXTRANS   2    // special states all have 2 transitions: move, loop
#define p7O_NTRANS    8    // 7 core transitions + BMk entry

enum p7o_xstates_e      { p7O_E    = 0, p7O_N    = 1,  p7O_J  = 2,  p7O_C  = 3 };
enum p7o_xtransitions_e { p7O_MOVE = 0, p7O_LOOP = 1 };
enum p7o_tsc_e          { p7O_BM   = 0, p7O_MM   = 1,  p7O_IM = 2,  p7O_DM = 3, p7O_MD   = 4, p7O_MI   = 5,  p7O_II = 6,  p7O_DD = 7 };
// Don't change the above order. Some routines assume p7O_DD is last; e.g. p7_oprofile_tqz_from_y(). Some routines do loops from p7O_BM..p7O_II.


typedef struct p7_oprofile_s {
  /* SSVFilter uses scaled int8_t scores: e.g. 16x per 128b vector                   */
  int8_t  **rbv;                /* match scores [x=0..Kp-1][q=0..Qb-1]               */
  int8_t   *rbv_mem;            /* one aligned allocation that <rbv> ptrs point into */
  float     tauBM;              /* constant B->Mk score:    log 2/M(M+1)             */  
  float     scale_b;            /* typically 3 / log2: scores scale to 1/3 bits      */

  /* ViterbiFilter uses scaled int16_t scores: e.g. 8x per 128b vector               */
  int16_t **rwv;                /* match scores [x=0..Kp-1][q=0..Qw-1]               */
  int16_t  *twv;                /* transition score blocks [8*Qw]                    */
  int16_t  *rwv_mem;            /* aligned allocation for <rwv>                      */
  int16_t  *twv_mem;            /* aligned allocation for <twv>                      */

  int16_t   xw[p7O_NXSTATES][p7O_NXTRANS]; /* ENJC state transition costs            */
  float     scale_w;            /* score units: typically 500 / log(2), 1/500 bits   */
  int16_t   base_w;             /* offset of int16 scores: typically +12000          */
  int16_t   ddbound_w;          /* threshold precalculated for lazy DD evaluation    */

  /* Forward, Backward use single-precision floats: e.g. 4x per 128b vector          */
  float  **rfv;                 /* match scores [x=0..Kp-1][q=0..Qf-1]               */
  float  *tfv;                  /* transition probability blocks [8*Qf]              */
  float  *rfv_mem;              /* aligned allocation for <rfv>                      */
  float  *tfv_mem;              /* aligned allocation for <tfv>                      */
  
  float    xf[p7O_NXSTATES][p7O_NXTRANS]; /* ENJC transition costs                   */

  /* Disk offset information for hmmpfam's fast model retrieval                      */
  off_t  offs[p7_NOFFSETS];     /* p7_{MFP}OFFSET, or -1                             */

  /* Disk offset bookkeeping for h3f:                                                */
  off_t  roff;                  /* record offset (start of record); -1 if none       */
  off_t  eoff;                  /* offset to last byte of record; -1 if unknown      */

  /* Information, annotation copied from parent profile:                             */
  char  *name;                  /* unique name of model                              */
  char  *acc;                   /* unique accession of model, or NULL                */
  char  *desc;                  /* brief (1-line) description of model, or NULL      */
  char  *rf;                    /* reference line           1..M; *rf='\0' = unused  */
  char  *mm;                    /* modelmask line           1..M; *mm='\0' = unused  */
  char  *cs;                    /* consensus structure line 1..M, *cs='\0' = unused  */
  char  *consensus;             /* consensus residues for ali display, 1..M          */
  float  evparam[p7_NEVPARAM];  /* parameters for determining E-values, or UNSET     */
  float  cutoff[p7_NCUTOFFS];   /* per-seq/per-dom bit cutoffs, or UNSET             */
  float  compo[p7_MAXABET];     /* per-model HMM filter composition, or UNSET        */
  const ESL_ALPHABET *abc;      /* copy of ptr to alphabet information               */

  /* Information about current configuration, size, allocation                       */
  int    L;                     /* current configured target seq length              */
  int    M;                     /* model length                                      */
  int    V;                     /* vectors contains V bytes (V/2 int16; V/4 float)   */
  int    max_length;            /* upper bound on emitted sequence length            */
  int    allocM;                /* maximum model length currently allocated for      */
  int    allocQb;               /* P7_Q(allocM,p7_VMAX_SSE): alloc size for rb       */
  int    allocQw;               /* P7_Q(allocM,p7_VMAX_VF):  alloc size for tw, rw   */
  int    allocQf;               /* P7_Q(allocM,p7_VMAX_FB):  alloc size for tf, rf   */

  int    mode;                  /* currently must be p7_LOCAL | p7_UNILOCAL          */
  float  nj;                    /* expected # of J's: 0 or 1, uni vs. multihit       */
  int    is_shadow;             /* TRUE if this <om> shadows another: ptrs are refs  */
} P7_OPROFILE;

typedef struct {
  int            count;       /* number of <P7_OPROFILE> objects in the block */
  int            listSize;    /* maximum number elements in the list          */
  P7_OPROFILE  **list;        /* array of <P7_OPROFILE> objects               */
} P7_OM_BLOCK;


// The definition of p7_oprofile_FGetEmission has moved into the p7_alidisplay_<simd> files to make it easier to
// provide SIMD ISA-specific versions of this inlined function.  The original code is retained here for archival purposes
/* retrieve match odds ratio [k][x]
 * this gets used in p7_alidisplay.c, when we're deciding if a residue is conserved or not */
/*static inline float 
p7_oprofile_FGetEmission(const P7_OPROFILE *om, int k, int x)
{
  union { __m128 v; float p[4]; } u;
  int   Q = P7_NVF(om->M);
  int   q = ((k-1) % Q);
  int   r = (k-1)/Q;
  u.v = om->rfv[x][q];
  return u.p[r];
} */

extern P7_OPROFILE *p7_oprofile_Create(int allocM, const ESL_ALPHABET *abc);
extern int          p7_oprofile_IsLocal(const P7_OPROFILE *om);
extern size_t       p7_oprofile_Sizeof (const P7_OPROFILE *om);
extern P7_OPROFILE *p7_oprofile_Shadow (const P7_OPROFILE *om);
extern void         p7_oprofile_Destroy(P7_OPROFILE *om);

extern int          p7_oprofile_Convert(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_ReconfigLength  (P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMultihit(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigUnihit  (P7_OPROFILE *om, int L);

extern int p7_oprofile_y_from_tk (int t, int k,        int Q, int V);
extern int p7_oprofile_y_from_tqz(int t, int q, int z, int Q, int V);
extern int p7_oprofile_k_from_tqz(int t, int q, int z, int Q, int V);
extern int p7_oprofile_tqz_from_y(int y, int Q, int V, int *ret_t, int *ret_q, int *ret_z);

extern int          p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om);
extern int          p7_oprofile_Sample(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, const P7_BG *bg, int M, int L,
				       P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **ret_om);
extern int          p7_oprofile_Compare(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);

#endif /*p7OPROFILE_INCLUDED*/
