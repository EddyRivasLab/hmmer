/* h4_sparsemx: DP matrix for production implementation
 */
#ifndef h4SPARSEMX_INCLUDED
#define h4SPARSEMX_INCLUDED
#include "h4_config.h"

#include "h4_refmx.h"
#include "h4_sparsemask.h"


typedef struct {
  float  *dp;		// main DP supercells. sm->ncells <= dalloc. each supercell contains h4S_NSCELLS values. 
  float  *xmx;		// special DP supercells. there are <sm->nrow>+<sm->S> of these, each w/ h4S_NXCELLS values. 

  int64_t dalloc;	// current <dp> allocation, denominated in supercells, (each h4S_NSCELLS wide) 
  int     xalloc;	// current <xmx> allocation, denominated in total rows (each h4S_NXCELLS wide; xalloc >= nrow+nseg) 

  const H4_SPARSEMASK *sm;

  int     type;		// h4S_UNSET | h4S_VITERBI | h4S_FORWARD... etc 
} H4_SPARSEMX;



/* One "supercell" contains 6 contiguous values for profile states.
 * Same indexing as in H4_REFMX.
 */
#define h4S_NSCELLS 6
#define h4S_ML 0
#define h4S_MG 1
#define h4S_IL 2
#define h4S_IG 3
#define h4S_DL 4
#define h4S_DG 5

/* "Special" states for each DP row, for states outside the core. 
 * Same indexing as in H4_REFMX.
 */
#define h4S_NXCELLS 9
#define h4S_E  0
#define h4S_N  1
#define h4S_J  2
#define h4S_B  3
#define h4S_L  4
#define h4S_G  5
#define h4S_C  6
#define h4S_JJ 7	/* in decoding (only) we separately decode J occupancy vs JJ emission */
#define h4S_CC 8	/* ditto for C */

/* The same data structure gets used in several DP contexts, which
 * have different boundary conditions/sentinel values; for example, in
 * a Decoding matrix, impossible cells are 0.0, whereas in a Viterbi,
 * Forward, or Backward matrix, impossible cells are -eslINFINITY. The
 * <type> field gets set by each algorithm implementation, and this
 * gets used for example by a _Validate() call.
 * 
 * Some of these codes must be sync'ed with p7_refmx.h. Some unit tests
 * compare reference, sparse matrices, including their type codes.
 */
#define h4S_UNSET         0  // = h4R_UNSET
#define h4S_FORWARD       1  // = h4R_FORWARD
#define h4S_BACKWARD      2  // = h4R_BACKWARD
#define h4S_DECODING      3  // = h4R_DECODING
#define h4S_VITERBI       4  // = h4R_VITERBI


/* H4_SPARSEMX */
extern H4_SPARSEMX *h4_sparsemx_Create(const H4_SPARSEMASK *sm);
extern int          h4_sparsemx_Reinit(H4_SPARSEMX *sx, const H4_SPARSEMASK *sm);
extern size_t       h4_sparsemx_Sizeof(const H4_SPARSEMX *sx);
extern size_t       h4_sparsemx_MinSizeof(const H4_SPARSEMASK *sm);
extern void         h4_sparsemx_Destroy(H4_SPARSEMX *sx);

/* debugging tools */
extern char *h4_sparsemx_DecodeState(int type);
extern char *h4_sparsemx_DecodeSpecial(int type);
extern int   h4_sparsemx_Dump(FILE *ofp, H4_SPARSEMX *sx);
extern int   h4_sparsemx_DumpWindow(FILE *ofp, const H4_SPARSEMX *sx, int ia, int ib, int ka, int kb);
extern int   h4_sparsemx_CompareReference(const H4_SPARSEMX *sx, const H4_REFMX *rx, float tol);
extern int   h4_sparsemx_CompareReferenceAsBound(const H4_SPARSEMX *sx, const H4_REFMX *rx, float tol);

extern int   h4_sparsemx_Validate(const H4_SPARSEMX *sx, char *errbuf);

#endif //h4SPARSEMX_INCLUDED
