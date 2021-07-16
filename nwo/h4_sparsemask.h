/* H4_SPARSEMASK: which cells are to be included in sparse DP
 */
#ifndef h4SPARSEMASK_INCLUDED
#define h4SPARSEMASK_INCLUDED
#include "h4_config.h"

#include "esl_random.h"

#include "h4_path.h"

#include "simdvec.h"

struct h4_sparsemask_seg_s {
    int ia;                // <seg[s]> = ia,ib pairs for each segment <s>, s=1..nseg. 
    int ib;                //    also, seg[0] and seg[S+1] contain sentinel values. 
};

typedef struct {
  int      L;		  // sequence has 1..L residues 
  int      M;		  // profile has 1..M positions 

  /* Dimensions that depend on SIMD vector size, in the interface to vectorized fb filter */
  int      V;             // width of a SIMD vector in floats (4 for SSE, for example); # of "slots"
  int      Q;		  // number of striped vectors in fwdfilter; width of each striped segment; size of "slots"

  struct h4_sparsemask_seg_s *seg;  // seg[1..S] segment locations
  int    **k;		  // k[0,1..L] = ptrs into kmem, rows of sparse k indices; k[0]=NULL; k[i]=NULL if n[i]=0 
  int     *n;		  // number of cells included on each row; n[0]=0; n[i] <= M 
  int     *kmem;	  // memory that k[] are pointing into, storing k indices of included cells, kmem[0..ncells-1]
  int      S;             // number of sparsified segments 
  int      nrow;          // number of included rows; \sum_{i=1}^{L} \delta(n[i]) 
  int64_t  ncells;        // number of included supercells; \sum_{i=1}^{L} n[i]        
  int      ralloc;        // k[] is allocated for ralloc rows; L+1 <= ralloc 
  int64_t  kalloc;        // kmem[] is allocated for kalloc cells; ncells <= kalloc 
  int      salloc;        // seg[] is allocated for salloc ia,ib pairs; S+2 <= salloc, +2 because of sentinels at 0,S+1

  /* "Slots" are used to convert striped vectors in f/b filter into correct M..1 cell index order in <kmem>; see note [3] */
  /* These array sizes (and last_k, below) must be able to hold up to the maximum vector width in floats: AVX-512, 16 floats. */
  int  *s[h4_VMAX_FB];   // slot pointers s[0..V-1] into <kmem>, for temporary storage of a striped vector row's sparse cells 
  int   sn[h4_VMAX_FB];  // number of sparse cells stored so far in each slot; sn[0..V-1]

  /* These are used for argument validation; the construction API is a little counterintuitive because i,q run in reverse */
  int   last_i;             // i of last StartRow(i) call; rows must be added in L..1 order; initialized to L+1 
  int   last_k[h4_VMAX_FB]; // k of last cell added in slot r; cells must be added in M..1 order; initialized to M+1 or -1 

  /* memory allocation profiling, statistics */
  int   n_krealloc;	  // number of times we reallocated <kmem> in this instance of the structure 
  int   n_rrealloc;	  // ditto for <k>, <n> 
  int   n_srealloc;	  // ditto for <seg>      
} H4_SPARSEMASK;

extern H4_SPARSEMASK *h4_sparsemask_Create(int M, int L);
extern int            h4_sparsemask_Reinit(H4_SPARSEMASK *sm, int M, int L);
extern size_t         h4_sparsemask_Sizeof   (const H4_SPARSEMASK *sm);
extern size_t         h4_sparsemask_MinSizeof(const H4_SPARSEMASK *sm);
extern int            h4_sparsemask_Reuse  (H4_SPARSEMASK *sm);
extern void           h4_sparsemask_Destroy(H4_SPARSEMASK *sm);

extern int            h4_sparsemask_StartRow (H4_SPARSEMASK *sm, int i);
extern int            h4_sparsemask_Add      (H4_SPARSEMASK *sm, int q, int r);
extern int            h4_sparsemask_FinishRow(H4_SPARSEMASK *sm);
extern int            h4_sparsemask_Finish   (H4_SPARSEMASK *sm);
extern int            h4_sparsemask_AddAll   (H4_SPARSEMASK *sm);

extern int            h4_sparsemask_Dump(FILE *ofp, H4_SPARSEMASK *sm);
extern int            h4_sparsemask_Compare(const H4_SPARSEMASK *sm1, const H4_SPARSEMASK *sm2);
extern int            h4_sparsemask_Validate(const H4_SPARSEMASK *sm, char *errbuf);
extern int            h4_sparsemask_SetFromTrace(H4_SPARSEMASK *sm, ESL_RANDOMNESS *rng, const H4_PATH *pi);


#endif // h4SPARSEMASK_INCLUDED
