/* Data structures used by sparse dynamic programming.
 * 
 * P7_SPARSEMASK defined the sparse cell mask: which i,k cells are to
 * be included in a DP matrix.
 * 
 * P7_SPARSEMX is a sparse DP matrix, as constrained by a P7_SPARSEMASK.
 * 
 * Contents:
 *    1. P7_SPARSEMASK structure declaration, indexing constants
 *    2. P7_SPARSEMX structure declaration, indexing constants
 *    3. Function declarations from p7_sparsemx.c
 *    4. Notes:
 *        [1] On the layout of P7_SPARSEMASK: why kmem[] is in reverse order during construction
 *        [2] On phases of construction of P7_SPARSEMASK: why k[], i[] aren't set until the end
 *        [3] On sorting striped indices; why four "slots" are used, then contiguated.
 *    5. Copyright and license information
 */
#ifndef p7SPARSEMX_INCLUDED
#define p7SPARSEMX_INCLUDED

#include "p7_config.h"

#include "esl_random.h"

#include "base/p7_trace.h"
#include "dp_reference/p7_refmx.h"
#include "dp_vector/simdvec.h"
 #include "hardware/hardware.h"
/*****************************************************************
 * 1. The P7_SPARSEMASK structure
 *****************************************************************/
typedef struct p7_sparsemask_seg_s {
    int ia;               // <seg[s]> = ia,ib pairs for each segment <s>, s=1..nseg. 
    int ib;               //    also, seg[0] and seg[nseg+1] contain sentinel values. 
  }       p7_sparsemask_seg_s;

typedef struct {
  int      L;		// sequence has 1..L residues 
  int      M;		// profile has 1..M positions 
  
  int      Q;		// number of striped vectors in fwdfilter; width of each striped segment; width of "slots"
  SIMD_TYPE simd;
 #ifdef HAVE_SSE2 
  p7_sparsemask_seg_s *seg;         
  int    **k;		// k[0,1..L] = ptrs into kmem, rows of sparse k indices; k[0]=NULL; k[i]=NULL if n[i]=0 
  int     *n;		// number of cells included on each row; n[0]=0; n[i] <= M 
  int     *kmem;	// memory that k[] are pointing into, storing k indices of included cells, kmem[0..ncells-1]
  int      S;     // number of sparsified segments 
  int      nrow;        // number of included rows; \sum_{i=1}^{L} \delta(n[i]) 
  int64_t  ncells;  // number of included supercells; \sum_{i=1}^{L} n[i]        
  int      ralloc;  // k[] is allocated for ralloc rows; L+1 <= ralloc 
  int64_t  kalloc;  // kmem[] is allocated for kalloc cells; ncells <= kalloc 
  int      salloc;  // seg[] is allocated for salloc ia,ib pairs; nseg+2 <= salloc, +2 because of sentinels at 0,nseg+1

  /* "Slots" are used to convert striped vectors in f/b filter into correct M..1 cell index order in <kmem>; see note [3] */
  int  *s[p7_VNF];  // slot pointers s[0..3] into <kmem>, for temporary storage of a striped vector row's sparse cells 
  int   sn[p7_VNF]; // number of sparse cells stored so far in each slot; sn[0..3]
#endif

#ifdef HAVE_AVX2  
  int      Q_AVX;   // number of striped vectors in fwdfilter; width of each striped segment; width of "slots"

  p7_sparsemask_seg_s *seg_AVX;         
  int    **k_AVX;   // k[0,1..L] = ptrs into kmem, rows of sparse k indices; k[0]=NULL; k[i]=NULL if n[i]=0 
  int     *n_AVX;   // number of cells included on each row; n[0]=0; n[i] <= M 
  int     *kmem_AVX;  // memory that k[] are pointing into, storing k indices of included cells, kmem[0..ncells-1]
  int      S_AVX;     // number of sparsified segments 
  int      nrow_AVX;        // number of included rows; \sum_{i=1}^{L} \delta(n[i]) 
  int64_t  ncells_AVX;  // number of included supercells; \sum_{i=1}^{L} n[i]        

  int      ralloc_AVX;  // k[] is allocated for ralloc rows; L+1 <= ralloc 
  int64_t  kalloc_AVX;  // kmem[] is allocated for kalloc cells; ncells <= kalloc 
  int      salloc_AVX;  // seg[] is allocated for salloc ia,ib pairs; nseg+2 <= salloc, +2 because of sentinels at 0,nseg+1

  /* "Slots" are used to convert striped vectors in f/b filter into correct M..1 cell index order in <kmem>; see note [3] */
  int  *s_AVX[p7_VNF_AVX];  // slot pointers s[0..7] into <kmem>, for temporary storage of a striped vector row's sparse cells 
  int   sn_AVX[p7_VNF_AVX]; // number of sparse cells stored so far in each slot; sn[0..7]
#endif

#ifdef HAVE_AVX512  
  int      Q_AVX_512;   // number of striped vectors in fwdfilter; width of each striped segment; width of "slots"

  struct p7_sparsemask_seg_s *seg_AVX_512;         
  int    **k_AVX_512;   // k[0,1..L] = ptrs into kmem, rows of sparse k indices; k[0]=NULL; k[i]=NULL if n[i]=0 
  int     *n_AVX_512;   // number of cells included on each row; n[0]=0; n[i] <= M 
  int     *kmem_AVX_512;  // memory that k[] are pointing into, storing k indices of included cells, kmem[0..ncells-1]
  int      S_AVX_512;     // number of sparsified segments 
  int      nrow_AVX_512;        // number of included rows; \sum_{i=1}^{L} \delta(n[i]) 
  int64_t  ncells_AVX_512;  // number of included supercells; \sum_{i=1}^{L} n[i]        

  int      ralloc_AVX_512;  // k[] is allocated for ralloc rows; L+1 <= ralloc 
  int64_t  kalloc_AVX_512;  // kmem[] is allocated for kalloc cells; ncells <= kalloc 
  int      salloc_AVX_512;  // seg[] is allocated for salloc ia,ib pairs; nseg+2 <= salloc, +2 because of sentinels at 0,nseg+1

  /* "Slots" are used to convert striped vectors in f/b filter into correct M..1 cell index order in <kmem>; see note [3] */
  int  *s_AVX_512[p7_VNF_AVX_512];  // slot pointers s[0..3] into <kmem>, for temporary storage of a striped vector row's sparse cells 
  int   sn_AVX_512[p7_VNF_AVX_512]; // number of sparse cells stored so far in each slot; sn[0..15]
#endif
  
  

  /* These are used for argument validation; the construction API is a little counterintuitive because i,q run in reverse */
 #ifdef HAVE_SSE2
   int   last_i;    // i of last StartRow(i) call; rows must be added in L..1 order; initialized to L+1 
  int   last_k[p7_VNF]; // k of last cell added in slot r; cells must be added in M..1 order; initialized to M+1 or -1 
#endif

#ifdef HAVE_AVX2
   int   last_i_AVX;    // i of last StartRow(i) call; rows must be added in L..1 order; initialized to L+1 
  int   last_k_AVX[p7_VNF_AVX]; // k of last cell added in slot r; cells must be added in M..1 order; initialized to M+1 or -1 
#endif

#ifdef HAVE_AVX512
   int   last_i_AVX_512;    // i of last StartRow(i) call; rows must be added in L..1 order; initialized to L+1 
  int   last_k_AVX_512[p7_VNF_AVX_512]; // k of last cell added in slot r; cells must be added in M..1 order; initialized to M+1 or -1 
#endif
  /* memory allocation profiling, statistics */
  int   n_krealloc;	// number of times we reallocated <kmem> in this instance of the structure 
  int   n_rrealloc;	// ditto for <k>, <n> 
  int   n_srealloc;	// ditto for <seg>      
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
extern P7_SPARSEMASK *p7_sparsemask_Create   (int M, int L, SIMD_TYPE simd);
extern P7_SPARSEMASK *p7_sparsemask_Create_sse   (int M, int L);
extern P7_SPARSEMASK *p7_sparsemask_Create_avx   (int M, int L);
extern P7_SPARSEMASK *p7_sparsemask_Create_avx512   (int M, int L);

extern int            p7_sparsemask_Reinit   (P7_SPARSEMASK *sm, int M, int L);
extern int            p7_sparsemask_Reinit_sse   (P7_SPARSEMASK *sm, int M, int L);
extern int            p7_sparsemask_Reinit_avx  (P7_SPARSEMASK *sm, int M, int L);
extern int            p7_sparsemask_Reinit_avx512   (P7_SPARSEMASK *sm, int M, int L);

extern size_t         p7_sparsemask_Sizeof   (const P7_SPARSEMASK *sm);
extern size_t         p7_sparsemask_Sizeof_sse   (const P7_SPARSEMASK *sm);
extern size_t         p7_sparsemask_Sizeof_avx   (const P7_SPARSEMASK *sm);
extern size_t         p7_sparsemask_Sizeof_avx512  (const P7_SPARSEMASK *sm);

extern size_t         p7_sparsemask_MinSizeof(const P7_SPARSEMASK *sm);
extern size_t         p7_sparsemask_MinSizeof_sse(const P7_SPARSEMASK *sm);
extern size_t         p7_sparsemask_MinSizeof_avx(const P7_SPARSEMASK *sm);
extern size_t         p7_sparsemask_MinSizeof_avx512(const P7_SPARSEMASK *sm);

extern int            p7_sparsemask_Reuse    (P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Reuse_sse    (P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Reuse_avx    (P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Reuse_avx512    (P7_SPARSEMASK *sm);

extern void           p7_sparsemask_Destroy  (P7_SPARSEMASK *sm);
extern void           p7_sparsemask_Destroy_sse  (P7_SPARSEMASK *sm);
extern void           p7_sparsemask_Destroy_avx  (P7_SPARSEMASK *sm);
extern void           p7_sparsemask_Destroy_avx512  (P7_SPARSEMASK *sm);

/* P7_SPARSEMASK construction API */
extern int            p7_sparsemask_AddAll   (P7_SPARSEMASK *sm);
extern int            p7_sparsemask_StartRow (P7_SPARSEMASK *sm, int i);
extern int            p7_sparsemask_Add      (P7_SPARSEMASK *sm, int q, int r);
extern int            p7_sparsemask_StartRow_sse (P7_SPARSEMASK *sm, int i);
extern int            p7_sparsemask_Add_sse      (P7_SPARSEMASK *sm, int q, int r);

extern int            p7_sparsemask_FinishRow(P7_SPARSEMASK *sm);
extern int            p7_sparsemask_FinishRow_sse(P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Finish_sse  (P7_SPARSEMASK *sm);
extern int            p7_sparsemask_StartRow_avx (P7_SPARSEMASK *sm, int i);
extern int            p7_sparsemask_Add_avx      (P7_SPARSEMASK *sm, int q, int r);
extern int            p7_sparsemask_FinishRow_avx(P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Finish_avx   (P7_SPARSEMASK *sm);

extern int            p7_sparsemask_StartRow_avx512 (P7_SPARSEMASK *sm, int i);
extern int            p7_sparsemask_Add_avx512      (P7_SPARSEMASK *sm, int q, int r);
extern int            p7_sparsemask_FinishRow_avx512(P7_SPARSEMASK *sm);
extern int            p7_sparsemask_FinishRow_avx512   (P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Finish_avx512   (P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Finish   (P7_SPARSEMASK *sm);

/* P7_SPARSEMASK debugging tools */
extern int            p7_sparsemask_Dump(FILE *ofp, P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Compare(const P7_SPARSEMASK *sm1, const P7_SPARSEMASK *sm2);
extern int            p7_sparsemask_Validate(const P7_SPARSEMASK *sm, char *errbuf);
extern int            p7_sparsemask_SetFromTrace(P7_SPARSEMASK *sm, ESL_RANDOMNESS *rng, const P7_TRACE *tr);

extern int            p7_sparsemask_Dump_sse(FILE *ofp, P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Compare_sse(const P7_SPARSEMASK *sm1, const P7_SPARSEMASK *sm2);
extern int            p7_sparsemask_Validate_sse(const P7_SPARSEMASK *sm, char *errbuf);
extern int            p7_sparsemask_SetFromTrace_sse(P7_SPARSEMASK *sm, ESL_RANDOMNESS *rng, const P7_TRACE *tr);

extern int            p7_sparsemask_Dump_avx(FILE *ofp, P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Compare_avx(const P7_SPARSEMASK *sm1, const P7_SPARSEMASK *sm2);
extern int            p7_sparsemask_Validate_avx(const P7_SPARSEMASK *sm, char *errbuf);
extern int            p7_sparsemask_SetFromTrace_avx(P7_SPARSEMASK *sm, ESL_RANDOMNESS *rng, const P7_TRACE *tr);

extern int            p7_sparsemask_Dump_avx512(FILE *ofp, P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Compare_avx512(const P7_SPARSEMASK *sm1, const P7_SPARSEMASK *sm2);
extern int            p7_sparsemask_Validate_avx512(const P7_SPARSEMASK *sm, char *errbuf);
extern int            p7_sparsemask_SetFromTrace_avx512(P7_SPARSEMASK *sm, ESL_RANDOMNESS *rng, const P7_TRACE *tr);

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
extern char *p7_sparsemx_DecodeState  (int type);
extern char *p7_sparsemx_DecodeSpecial(int type);
extern int   p7_sparsemx_Dump(FILE *ofp, P7_SPARSEMX *sx);
extern int   p7_sparsemx_DumpWindow(FILE *ofp, const P7_SPARSEMX *sx, int i1, int i2, int ka, int kb);
extern int   p7_sparsemx_Copy2Reference  (const P7_SPARSEMX *sx, P7_REFMX *rx);
extern int   p7_sparsemx_Compare(const P7_SPARSEMX *sx1, const P7_SPARSEMX *sx2, float tol);
extern int   p7_sparsemx_CompareReference       (const P7_SPARSEMX *sx, const P7_REFMX *rx, float tol);
extern int   p7_sparsemx_CompareReferenceAsBound(const P7_SPARSEMX *sx, const P7_REFMX *rx, float tol);
extern int   p7_sparsemx_CompareDecoding        (const P7_SPARSEMX *sxe, const P7_SPARSEMX *sxa, float tol);
extern int   p7_sparsemx_Validate(const P7_SPARSEMX *sx, char *errbuf);
extern int   p7_sparsemx_Validate_sse(const P7_SPARSEMX *sx, char *errbuf);
extern int   p7_sparsemx_Validate_avx(const P7_SPARSEMX *sx, char *errbuf);
extern int   p7_sparsemx_Validate_avx512(const P7_SPARSEMX *sx, char *errbuf);
extern int   p7_sparsemx_PlotDomainInference(FILE *ofp, const P7_SPARSEMX *sxd, int ia, int ib, const P7_TRACE *tr);


/*****************************************************************
 * 4. Notes
 *****************************************************************/

/* [1] Notes on the layout of the P7_SPARSEMASK structure.
 * 
 *     Structure is designed for how the vectorized forward/backward
 *     filter works. The filter does posterior decoding on a backwards
 *     pass from L..1, in linear memory, discarding decoded rows
 *     immediately. Therefore sparse cells are found in reverse order,
 *     at least w.r.t. rows i. However, we don't know how many sparse
 *     cells there are, and may need to reallocate k index memory
 *     (<kmem>) during collection. This makes it problematic to store
 *     sparse cells in 1..L order. But we want them so, and we want
 *     them packed and contiguous. 
 *    
 *     So: during construction, indices are stored in reverse order;
 *     when construction is done, we reverse the entire kmem[] array.
 *     
 *     Benchmarks disagree on whether this has an impact. In
 *     fwdfilter_benchmark, running Caudal_act or Patched against
 *     random sequences, wall time shows a ~5% hit; but gprof
 *     profiling claims a negligible hit [SRE:J10/39]. Even if
 *     it's as much as 5%, it's worth it in code clarity.

 *     So, during construction, indices are
 *     stored in reverse order in <kmem>. For example, for a sparse
 *     matrix
 *          k= 1 2 3 4
 *         i=1 o o o .
 *           2 . . o .
 *           3 . . o o
 *           4 . . . .
 *     kmem gets built up in reverse, but n[i] counters are in the correct direction:
 *        kmem[]:   4 3 3 3 2 1
 *        n[]:      0 3 1 2 0
 *        ncells:   6     
 *     Then in _Finish() we reverse kmem, and set k[] ptrs:
 *        kmem[]:   1 2 3 3 3 4
 *        n[]:      0 3 1 2 0
 *        k[]:      NULL kmem kmem+3 kmem+4 NULL
 *     i.e.
 *        kmem:     [ 1 2 3 ] [ 3 ] [ 3 4 ]
 *               x    ^        ^     ^       x
 *               k[0] k[1]    k[2]   k[3]    k[4]   
 *            n[0]=0  n[1]=3  n[2]=1 n[3]=2  n[4]=0

 *     To traverse sparse cells in forward order:
 *        for (i = 1; i <= L; i++)
 *          for (z = 0; z < n[i]; z++)
 *            k = k[i][z];
 *     To traverse them in backwards order:
 *        for (i = L; i >= 1; i--)
 *           for (z = n[i]-1; z >= 0, z--)
 *              k = k[i][z];
 *
 * [2] Notes on how the P7_SPARSEMASK structure is allocated and built.
 * 
 *     Constraints in play:
 *       - minimize reallocation and maximize data contiguity
 *       - fb filter collects rows in reverse order i=L..1,k=M..1, then reverses the whole kmem[]
 *       - fb filter is vectorized and its k indices are striped, out of order

 *     Some fields are set at creation (or reinit) time; some during
 *     collection of sparse cells; and some when collection is
 *     finished; call these "create", "collect", and "finish" time.
 *     
 *     L,M,ralloc are set at creation/reinit time, as the size of the
 *     DP problem. We need to know L at create time because k[] and
 *     n[] need to be allocated 0,1..L; on reinit, if L+1 > ralloc, we
 *     reallocate them. 
 *     
 *     During collect time, we set <n[]> as we append to <kmem> array in
 *     reverse order, counting <ncells>. If <kmem> must be
 *     reallocated, <kalloc> will change.
 *     
 *     At finish time, we reverse kmem[], then set everything
 *     else. Row pointers <k[]> are set, using <n[]>. Segments are
 *     determined, reallocating <seg[]> and resetting <salloc> if
 *     needed, setting <seg[]>, <nseg>, and <nrow>.
 *     
 *     
 *     
 * [3] Sorting striped indices 
 *     Because the rows in the f/b filter are striped, sparse cells
 *     are found out of order, and need to be sorted M..1 for storage
 *     in <kmem>. We can avoid an nlog(n) sort by temporarily storing
 *     the Q independent striped segments in separate "slots", and 
 *     contiguating the slots afterwards. The slots are set up at the
 *     end of <kmem>. This dictates the api for collecting cells on 
 *     each row:
 *
 *         p7_sparsemask_StartRow(sm, i)  
 *           sets up slot ptrs s[] in kmem; reallocates kmem if needed;
 *           zeros slot counts sn[].
 *           
 *         for each sparse cell k, accessed in reverse order M..1 in each slot, slots in any order:
 *           p7_sparsemask_Add(sm, i, k, slot)
 *             appends k index to slot vector s[slot], increments slot count sn[slot]
 *             
 *         p7_sparsemask_FinishRow(sm, i)
 *           Concats slots onto kmem, increments ncells. 
 *           Slots are now invalid until the next StartRow().
 *           
 *      These only need to be called on a row with one or more sparse
 *      cells.  If no sparse cells are added, FinishRow() need not be
 *      called, but it's harmless to call it anyway.
 *      
 *      What this looks like in kmem, for V=4 slots and Q=3 cells max per slot,
 *      after we set it up with StartRow():
 *        kmem = [ 0 .. ncells-1] [ . . . ] [ . . . ] [ . . . ] [ . . . ]
 *                                  ^         ^         ^         ^
 *                                  s[3]      s[2]      s[1]      s[0]
 *                                 sn[3]=0   sn[2]=0   sn[1]=0   sn[0]=0
 *                                 
 *      As we're filling the slots with Add() it might look something like:                           
 *        kmem = [ 0 .. ncells-1] [ 11 10 . ] [ . . . ] [ 5 4 . ] [ 1 . . ]
 *                                  ^           ^         ^         ^
 *                                  s[3]        s[2]      s[1]      s[0]
 *                                 sn[3]=2     sn[2]=0   sn[1]=2   sn[0]=0
 *        
 *      and when we collapse it with FinishRow():               
 *        kmem = [ 0 ..  11 10 5 4 1 ]
 *      with ncells is incremented by 5. Remember, kmem[] is collected in reverse
 *      order during collection.
 */


#endif /*p7SPARSEMX_INCLUDED*/ 
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
