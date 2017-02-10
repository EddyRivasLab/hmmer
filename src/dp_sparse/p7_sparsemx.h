/* Data structures used by sparse dynamic programming.
 * 
 * P7_SPARSEMASK defined the sparse cell mask: which i,k cells are to
 * be included in a DP matrix.
 * 
 * P7_SPARSEMX is a sparse DP matrix, as constrained by a P7_SPARSEMASK.
 * 
 * Contents:
 *    1. P7_SPARSEMASK structure declaration
 *    2. P7_SPARSEMX structure declaration, indexing constants
 *    3. Function declarations from p7_sparsemx.c
 */
#ifndef p7SPARSEMX_INCLUDED
#define p7SPARSEMX_INCLUDED
#include "p7_config.h"

#include "esl_random.h"

#include "base/p7_trace.h"
#include "dp_reference/p7_refmx.h"

/*****************************************************************
 * 1. The P7_SPARSEMASK structure
 *****************************************************************/

/* Maximum width of SIMD vectors in any compiled implementation.
 * (Current widest = AVX-512.)  
 * Perhaps this should be in dp_vector/simdvec.h, or even p7_config.h
 * We could check for what implementations we support.
 */
#define p7_MAXVB  64       // in bytes
#define p7_MAXVW  32       // in words (int16)
#define p7_MAXVF  16       // in floats (or int32)

#define p7_VDEFAULT 4      // default/placeholder V to use in sparsemask_Create(), _Reinit() calls outside vector code

struct p7_sparsemask_seg_s {
    int ia;                // <seg[s]> = ia,ib pairs for each segment <s>, s=1..nseg. 
    int ib;                //    also, seg[0] and seg[nseg+1] contain sentinel values. 
};

typedef struct {
  int      L;		  // sequence has 1..L residues 
  int      M;		  // profile has 1..M positions 

  /* Dimensions that depend on SIMD vector size, in the interface to vectorized fb filter */
  int      Q;		  // number of striped vectors in fwdfilter; width of each striped segment; size of "slots"
  int      V;             // width of a SIMD vector in floats (4 for SSE, for example); # of "slots"

  struct p7_sparsemask_seg_s *seg;         
  int    **k;		  // k[0,1..L] = ptrs into kmem, rows of sparse k indices; k[0]=NULL; k[i]=NULL if n[i]=0 
  int     *n;		  // number of cells included on each row; n[0]=0; n[i] <= M 
  int     *kmem;	  // memory that k[] are pointing into, storing k indices of included cells, kmem[0..ncells-1]
  int      S;             // number of sparsified segments 
  int      nrow;          // number of included rows; \sum_{i=1}^{L} \delta(n[i]) 
  int64_t  ncells;        // number of included supercells; \sum_{i=1}^{L} n[i]        
  int      ralloc;        // k[] is allocated for ralloc rows; L+1 <= ralloc 
  int64_t  kalloc;        // kmem[] is allocated for kalloc cells; ncells <= kalloc 
  int      salloc;        // seg[] is allocated for salloc ia,ib pairs; nseg+2 <= salloc, +2 because of sentinels at 0,nseg+1

  /* "Slots" are used to convert striped vectors in f/b filter into correct M..1 cell index order in <kmem>; see note [3] */
  /* These array sizes (and last_k, below) must be able to hold up to the maximum vector width in floats: AVX-512, 16 floats. */
  int  *s[p7_MAXVF];      // slot pointers s[0..V-1] into <kmem>, for temporary storage of a striped vector row's sparse cells 
  int   sn[p7_MAXVF];     // number of sparse cells stored so far in each slot; sn[0..V-1]

  /* These are used for argument validation; the construction API is a little counterintuitive because i,q run in reverse */
  int   last_i;           // i of last StartRow(i) call; rows must be added in L..1 order; initialized to L+1 
  int   last_k[p7_MAXVF]; // k of last cell added in slot r; cells must be added in M..1 order; initialized to M+1 or -1 

  /* memory allocation profiling, statistics */
  int   n_krealloc;	  // number of times we reallocated <kmem> in this instance of the structure 
  int   n_rrealloc;	  // ditto for <k>, <n> 
  int   n_srealloc;	  // ditto for <seg>      
} P7_SPARSEMASK;


/*****************************************************************
 * 2. The P7_SPARSEMX structure
 *****************************************************************/

#define p7S_NSCELLS 6
#define p7S_ML 0
#define p7S_MG 1
#define p7S_IL 2
#define p7S_IG 3
#define p7S_DL 4
#define p7S_DG 5

#define p7S_NXCELLS 9
#define p7S_E  0
#define p7S_N  1
#define p7S_J  2
#define p7S_B  3
#define p7S_L  4
#define p7S_G  5
#define p7S_C  6
#define p7S_JJ 7	/* in decoding (only) we separately decode J occupancy vs JJ emission */
#define p7S_CC 8	/* ditto for C */


/* the same data structure gets used in several DP contexts, which
 * have different boundary conditions/sentinel values; for example, in
 * a Decoding matrix, impossible cells are 0.0, whereas in a Viterbi,
 * Forward, or Backward matrix, impossible cells are -eslINFINITY. The
 * <type> field gets set by each algorithm implementation, and this
 * gets used for example by a _Validate() call.
 * 
 * Some of these codes must be sync'ed with p7_refmx.h. Some unit tests
 * compare reference, sparse matrices, including their type codes.
 */
#define p7S_UNSET         0  // = p7R_UNSET
#define p7S_FORWARD       1  // = p7R_FORWARD
#define p7S_BACKWARD      2  // = p7R_BACKWARD
#define p7S_DECODING      3  // = p7R_DECODING
#define p7S_VITERBI       4  // = p7R_VITERBI
#define p7S_AEC_ALIGN     5  // = p7R_AEC_ALIGN
#define p7S_ASC_FWD       6
#define p7S_ASC_BCK       7
#define p7S_ASC_DECODE    8
#define p7S_ENVSCORE      9   // deprecated. Will remove, when sparse_envscore.c goes away.
#define p7S_MASSTRACE     10  //  ... ditto, for sparse_masstrace.c

typedef struct {
  float  *dp;		// main DP supercells. sm->ncells <= dalloc. each supercell contains p7S_NSCELLS values. 
  float  *xmx;		// special DP supercells. there are <sm->nrow>+<sm->S> of these, each w/ p7S_NXCELLS values. 

  int64_t dalloc;	// current <dp> allocation, denominated in supercells, (each p7S_NSCELLS wide) 
  int     xalloc;	// current <xmx> allocation, denominated in total rows (each p7S_NXCELLS wide; xalloc >= nrow+nseg) 

  const P7_SPARSEMASK *sm;

  int     type;		// p7S_UNSET | p7R_VITERBI | p7R_FORWARD... etc 
} P7_SPARSEMX;

/*****************************************************************
 * 3. Function declarations from p7_sparsemx.c, sparse_fwdback.c
 *****************************************************************/

/* P7_SPARSEMASK object management */
extern P7_SPARSEMASK *p7_sparsemask_Create   (int M, int L, int V);
extern int            p7_sparsemask_Reinit   (P7_SPARSEMASK *sm, int M, int L, int V);
extern size_t         p7_sparsemask_Sizeof   (const P7_SPARSEMASK *sm);
extern size_t         p7_sparsemask_MinSizeof(const P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Reuse    (P7_SPARSEMASK *sm);
extern void           p7_sparsemask_Destroy  (P7_SPARSEMASK *sm);

/* P7_SPARSEMASK construction API */
extern int            p7_sparsemask_Add      (P7_SPARSEMASK *sm, int q, int r);
extern int            p7_sparsemask_StartRow (P7_SPARSEMASK *sm, int i);
extern int            p7_sparsemask_FinishRow(P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Finish   (P7_SPARSEMASK *sm);
extern int            p7_sparsemask_AddAll   (P7_SPARSEMASK *sm);

/* P7_SPARSEMASK debugging tools */
extern int            p7_sparsemask_Dump        (FILE *ofp, P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Compare     (const P7_SPARSEMASK *sm1, const P7_SPARSEMASK *sm2);
extern int            p7_sparsemask_Validate    (const P7_SPARSEMASK *sm, char *errbuf);
extern int            p7_sparsemask_SetFromTrace(P7_SPARSEMASK *sm, ESL_RANDOMNESS *rng, const P7_TRACE *tr);

/* P7_SPARSEMX object management */
extern P7_SPARSEMX   *p7_sparsemx_Create   (P7_SPARSEMASK *sm);
extern int            p7_sparsemx_Reinit   (P7_SPARSEMX *sx, const P7_SPARSEMASK *sm);
extern int            p7_sparsemx_Zero     (P7_SPARSEMX *sx);
extern size_t         p7_sparsemx_Sizeof   (const P7_SPARSEMX *sx);
extern size_t         p7_sparsemx_MinSizeof(const P7_SPARSEMASK *sm); // not a typo: yes, it takes a SPARSEMASK, not a SPARSEMX
extern float          p7_sparsemx_GetSpecial(const P7_SPARSEMX *sx, int i, int s);
extern int            p7_sparsemx_Reuse    (P7_SPARSEMX *sx);
extern void           p7_sparsemx_Destroy  (P7_SPARSEMX *sx);

/* Extracting information from a sparse DP matrix */
extern int   p7_sparsemx_TracePostprobs(const P7_SPARSEMX *sxd, P7_TRACE *tr);
extern int   p7_sparsemx_CountTrace(const P7_TRACE *tr, P7_SPARSEMX *sxd);
extern int   p7_sparsemx_ExpectedDomains(const P7_SPARSEMX *sxd, int iae, int ibe, float *ret_ndom_expected);

/* P7_SPARSEMX debugging tools */
extern char *p7_sparsemx_DecodeState(int type);
extern char *p7_sparsemx_DecodeSpecial(int type);
extern int   p7_sparsemx_Dump(FILE *ofp, P7_SPARSEMX *sx);
extern int   p7_sparsemx_DumpWindow(FILE *ofp, const P7_SPARSEMX *sx, int i1, int i2, int ka, int kb);
extern int   p7_sparsemx_Copy2Reference(const P7_SPARSEMX *sx, P7_REFMX *rx);
extern int   p7_sparsemx_Compare(const P7_SPARSEMX *sx1, const P7_SPARSEMX *sx2, float tol);
extern int   p7_sparsemx_CompareReference(const P7_SPARSEMX *sx, const P7_REFMX *rx, float tol);
extern int   p7_sparsemx_CompareReferenceAsBound(const P7_SPARSEMX *sx, const P7_REFMX *rx, float tol);
extern int   p7_sparsemx_CompareDecoding(const P7_SPARSEMX *sxe, const P7_SPARSEMX *sxa, float tol);
extern int   p7_sparsemx_Validate(const P7_SPARSEMX *sx, char *errbuf);



extern int   p7_sparsemx_PlotDomainInference(FILE *ofp, const P7_SPARSEMX *sxd, int ia, int ib, const P7_TRACE *tr);

#endif /*p7SPARSEMX_INCLUDED*/ 
