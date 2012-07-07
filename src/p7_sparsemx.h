
#ifndef P7_SPARSEMX_INCLUDED
#define P7_SPARSEMX_INCLUDED

#include "p7_refmx.h"


typedef struct {
  int      L;		/* sequence has 1..L residues */
  int      M;		/* profile has 1..M positions */
  int      Q;		/* number of striped vectors in fwdfilter; width of striped segments */
  int      W;		/* maximum number of sparse cells per row; <=M; 1/p for posterior probability p per cell */
  int      Ws;		/* maximum # of elements needed in each "slot" for temporary storage from striped vectors */

  int     *i;    	/* ia,ib pairs for each segment. i[0..2*nseg-1] */
  int    **k;		/* k[0,1..L] = ptrs into kmem, rows of sparse k indices; k[0]=NULL; k[i]=NULL if n[i]=0 */
  int     *n;		/* number of cells included on each row; n[0]=0; n[i] <= M */
  int     *kmem;	/* memory that k[] are pointing into, storing k indices of included cells */

  int      nseg;	/* number of sparsified segments */
  int      nrow;        /* number of included rows; \sum_{i=1}^{L} \delta(n[i]) */
  int64_t  ncells;	/* number of included cells; \sum_{i=1}^{L} n[i]        */

  int      ralloc;	/* k[] is allocated for ralloc rows; L+1 <= ralloc */
  int64_t  kalloc;	/* kmem[] is allocated for kalloc cells; ncells <= kalloc */
  int      ialloc;	/* i[] is allocated for up to 2*ialloc coords ia,ib; nseg <= ialloc */

  /* "Slots" are used to convert striped vectors in f/b filter into correct M..1 cell index order in <kmem>; see note [3] */
  int  *s[p7_VNF];	/* slot pointers s[0..Q-1] into <kmem>, for temporary storage of a striped vector row's sparse cells */
  int   sn[p7_VNF];	/* number of sparse cells stored so far in each slot; sn[0..p7_VNF-1] */

  /* These are used for argument validation; the API is a little counterintuitive */
  int   last_i;  	/* i of last StartRow(i) call; rows must be added in L..1 order; initialized to L+1 */
  int   last_k[p7_VNF]; /* k of last cell added in slot r; cells must be added in M..1 order; initialized to M+1 or -1 */

  /* memory allocation profiling, statistics */
  int   n_krealloc;	/* number of times we reallocated <kmem> in this instance of the structure */
  int   n_irealloc;	/* ditto for <i>      */
  int   n_rrealloc;	/* ditto for <k>, <n> */
} P7_SPARSEMASK;


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
 */
#define p7S_UNSET     0
#define p7S_FORWARD   1
#define p7S_BACKWARD  2
#define p7S_DECODING  3
#define p7S_ALIGNMENT 4
#define p7S_VITERBI   5


typedef struct {
  float  *dp;		/* main DP supercells. sm->ncells <= dalloc. each supercell contains p7S_NSCELLS values. */
  float  *xmx;		/* special DP supercells. there are <sm->nrow>+<sm->nseg> of these, each w/ p7S_NXCELLS values. */

  int64_t dalloc;	/* current <dp> allocation, denominated in supercells, (each p7S_NSCELLS wide) */
  int     xalloc;	/* current <xmx> allocation, denominated in total rows (each p7S_NXCELLS wide; xalloc >= nrow+nseg) */

  P7_SPARSEMASK *sm;

  int     type;		/* p7S_UNSET | p7R_VITERBI | p7R_FORWARD... etc */
} P7_SPARSEMX;

extern P7_SPARSEMASK *p7_sparsemask_Create(int M, int L, float pthresh);
extern int            p7_sparsemask_Reinit (P7_SPARSEMASK *sm, int M, int L, float pthresh);
extern int            p7_sparsemask_Reuse  (P7_SPARSEMASK *sm);
extern void           p7_sparsemask_Destroy(P7_SPARSEMASK *sm);

extern int            p7_sparsemask_AddAll   (P7_SPARSEMASK *sm);
extern int            p7_sparsemask_StartRow (P7_SPARSEMASK *sm, int i);
extern int            p7_sparsemask_Add      (P7_SPARSEMASK *sm, int q, int r);
extern int            p7_sparsemask_FinishRow(P7_SPARSEMASK *sm);
extern int            p7_sparsemask_Finish   (P7_SPARSEMASK *sm);

extern int            p7_sparsemask_Dump(FILE *ofp, P7_SPARSEMASK *sm);

extern P7_SPARSEMX   *p7_sparsemx_Create(P7_SPARSEMASK *sm);
extern int            p7_sparsemx_Reinit (P7_SPARSEMX *sx, P7_SPARSEMASK *sm);
extern int            p7_sparsemx_Reuse  (P7_SPARSEMX *sx);
extern void           p7_sparsemx_Destroy(P7_SPARSEMX *sx);

extern char *p7_sparsemx_DecodeSpecial(int type);
extern int   p7_sparsemx_Dump(FILE *ofp, P7_SPARSEMX *sx);
extern int   p7_sparsemx_DumpWindow(FILE *ofp, P7_SPARSEMX *sx, int i1, int i2, int ka, int kb);
extern int   p7_sparsemx_Copy2Reference(P7_SPARSEMX *sx, P7_REFMX *rx);

#endif /*P7_SPARSEMX_INCLUDED*/


/* [1] Notes on the layout of the P7_SPARSEMASK structure.
 * 
 *     Structure is designed for how the vectorized forward/backward
 *     filter works. The filter does posterior decoding on a backwards
 *     pass from L..1, in linear memory, discarding decoded rows
 *     immediately. Therefore sparse cells are found in reverse order,
 *     at least w.r.t. rows i. However, we don't know how many sparse
 *     cells there are, and may need to reallocate k index memory
 *     (<kmem>) during collection. This makes it problematic to store
 *     sparse cells in 1..L order.
 *     
 *     In hope of optimizing contiguity, sparse cell indices are
 *     stored in reverse order in <kmem>. For example, for a sparse
 *     matrix
 *          k= 1 2 3 4
 *         i=1 o o o .
 *           2 . . o .
 *           3 . . o o
 *           4 . . . .
 *     the fields associated with sparse cells are:
 *        kmem[]:   4 3 3 3 2 1
 *        n[]:      0 3 1 2 0
 *        k[]:      NULL kmem+3 kmem+2 kmem NULL
 *        ncells:   6
 *     i.e.
 *        kmem:     [ 4  3 ] [ 3 ] [ 3 2 1 ]
 *               x    ^        ^     ^       x
 *               k[4] k[3]     k[2]  k[1]    k[0]   
 *            n[4]=0  n[3]=2  n[2]=1 n[1]=3  n[0]=0
 *            
 *     To traverse sparse cells in forward order:
 *        for (i = 1; i <= L; i++)
 *          for (z = n[i]-1; z; z--)
 *            k = k[i][z];
 *     To traverse them in backwards order:
 *        for (i = L; i >= 1; i--)
 *           for (z = 0; z < n[i], z++)
 *              k = k[i][z];
 *
 *     We're sort of hoping that the z--/z++ order of evaluation is
 *     enough of a system hint to prefetch memory in the correct
 *     traversal direction along kmem.
 *     
 * [2] Notes on how the P7_SPARSEMASK structure is allocated and built.
 * 
 *     Constraints in play:
 *       - minimize reallocation and maximize data contiguity
 *       - fb filter collects rows in reverse order i=L..1,k=M..1
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
 *     At collect time, we set <n[]> as we append to <kmem> array in
 *     reverse order, counting <ncells>. If <kmem> must be
 *     reallocated, <kalloc> will change.
 *     
 *     At finish time, we set everything else. Row pointers <k[]> are
 *     set, using <n[]>. Segments are determined, reallocating <i[]>
 *     and resetting <ialloc> if needed, setting <i[]>, <nseg>, and
 *     <nrow>.
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
 *         p7_sparsemask_FihishRow(sm, i)
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
 *      with ncells is incremented by 5.
 */


 
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
