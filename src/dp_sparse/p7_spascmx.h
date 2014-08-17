/* There is no P7_SPASCMX structure. Rather, the p7_spascmx code
 * provides for using P7_SPARSEMX and P7_SPARSEMASK in anchor set
 * constrained DP calculations.
 */
#ifndef P7_SPASCMX_INCLUDED
#define P7_SPASCMX_INCLUDED

#include "p7_config.h"

#include "base/p7_anchors.h"
#include "base/p7_profile.h"
#include "dp_sparse/p7_sparsemx.h"


/* To create a sparse ASC matrix: see p7_sparsemx_Create()
 * To reuse it:                   see p7_sparsemx_Reuse()
 * To destroy it:                 see p7_sparsemx_Destroy()
 *
 * Only _Reinit() needs to be specialized, because that's where we
 * need to know the exact # of cells in the sparse ASC matrix.
 */

extern int    p7_spascmx_Resize(P7_SPARSEMX *asx, const P7_SPARSEMASK *sm, const P7_ANCHOR *anch, int D);
extern size_t p7_spascmx_MinSizeof(const P7_SPARSEMASK *sm, const P7_ANCHOR *anch, int D, int64_t *opt_dalloc, int *opt_xalloc);
extern int    p7_spascmx_Dump(FILE *fp, const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D);
extern int    p7_spascmx_CompareReference(const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D, const P7_REFMX *rxu, const P7_REFMX *rxd, float tol);
extern int    p7_spascmx_Validate(const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D, char *errbuf);




/*****************************************************************
 * 2. Footnotes
 ***************************************************************** 
 *
 * [1] ON SPECIAL CASE VALUES IN SPARSE ASC MATRIX
 * 
 * The _Validate() routine checks for the following patterns of
 * unreachable values in the Fwd, Bck, and Decoding matrices:
 *
 *  FORWARD:
 *          [M I D]..[M I D]..[M I D]  [M I D]  [M I D]..[M I D]..[M I D]
 *  k:         1        k      k0-1      k0      k0+1       k        M        
 *    ------------------------------------------------------------------
 *    u1-1 | * * *    * * *    * * * |  <= UP sector
 *    u1   | . * *    . * .    . * . |
 *    ..   | . . *   {. . .}   . . . |
 *    u2   | . . *    . . .    . . . |
 *    ------------------------------------------------------------------
 *    i0   |                         |  . * *    * * .    * * .    * * .
 *    i0+1 |         DOWN sector =>  |  * . *    . * *    . * .    . * .  
 *    ..   |                         |  * . *    . . *   {. . .}   . * .
 *    d2   |                         |  * . *    . . *    . . .    . * . 
 *
 *  
 *  BACKWARD:
 *          [M I D]..[M I D]..[M I D]  [M I D]  [M I D]..[M I D]..[M I D]
 *  k:         1        k      k0-1      k0      k0+1       k        M
 *    ------------------------------------------------------------------
 *    u1-1 | * * *    * * *    * * * |  <= UP sector
 *    u1   | . . l    . . .    . . * |
 *    ..   | . . l   {. . .}   . . * |
 *    u2   | . * l    . * .    . . . |
 *    ------------------------------------------------------------------
 *    i0   |                         |  . . .    . . .    . . .    . * .
 *    ..   |         DOWN sector =>  |  . . .    . . .   {. . .}   . * .  
 *    d2   |                         |  . * .    . * .    . * .    . * . 
 *
 * 
 *  DECODING:
 *          [M I D]..[M I D]..[M I D]  [M I D]  [M I D]..[M I D]..[M I D]
 *  k:         1        k      k0-1      k0      k0+1       k        M      
 *    ----------------------------------------------------------------------
 *    u1-1 | 0 0 g    0 0 g    0 0 ? |  <= UP sector                        
 *    u1   | . 0 g    . 0 .    . 0 ? |                                      
 *    ..   | . . g   {. . .}   . . 0 |                                      
 *    u2   | . 0 g    . 0 .    . . . |                                      
 *    ------------------------------------------------------------------
 *    i0   |                         |  . 0 0    0 0 .    0 0 .    0 0 .
 *    i0+1 |         DOWN sector =>  |  0 . 0    . 0 0    . 0 .    . 0 .
 *    ..   |                         |  0 . 0    . . 0   {. . .}   . 0 .
 *    d2   |                         |  0 0 0    . 0 0    . 0 .    . 0 .
 * 
 * Rationale for UP sector unreachable special cases:
 *  u1-1 (prev to 1st row) not in UP sector, but has special needs; see [1].
 *  M states are always reachable.
 *  In top row u1, IL/IG unreached in Fwd, because no M/I on row u1-1.
 *  In bot row u2, for all k<k0-1, IL/IG unreached in Bck, because no row u2+1;
 *     IL/IG(k0-1) do have path to anchor cell. 
 *  In 1st col k=1, DL1 unreachable in Fwd, because no M/D at k=0; 
 *     DG1 unreachable in F/B because of wing retraction;
 *     DG1 valid in decoding because of G->DG1..DGk-1->MGk wing unfold
 *  In last col k0-1, for all i<u2, DL/DG unreachable in Bck, because no k+1;
 *     DL/DG at u2,k0-1 do have path to anchor cell.
 *  If 0 rows in UP: u1-1 special case only: and k0-1 case is 0 0 g! (hence ? mark)
 *  If 1 row in UP:  u1-1, u1 : and k0-1 case is 0 0 g!
 *  If 2 rows in UP: u1-1 u1, u2
 *
 * Rationale for DOWN sector unreachable special cases:
 *  In anchor cell i0,k0 (top left corner of DOWN sector) only MG/ML are used.
 *  For rest of top row i0, ML/MG, IL/IG unreachable in Fwd, because no i-1 row
 *  For rest of left col k0, ML/MG, DL/DG unreachable in Fwd, because no k-1 col (and no LG->Mk entry in DOWN)
 *  In bot row d2, IL/IG unreachable in Bck, because no i+1 row
 *  In 2nd col k0+1, DL/DG unreachable in Fwd, because M/D unreachable in k-1 col 
 *  In 2nd row i0+1, k>k0, IL/IG unreachable in Fwd, because M/I unreachable on i0 row except anchor i0,k0
 *  In last col k=M, IL/IGm state never exists.
 *  If 1 row in DOWN:  i0 
 *  If 2 rows in DOWN: i0, d2
 *  If 3 rows in DOWN: i0, i0+1, d2
 *  
 *****************************************************************
 *
 * [2] ON STORAGE OF UP SECTOR BOUNDARY ROWS i1-1
 *    
 * The "i1-1" row is the row previous to the first row of the UP
 * sector. This may either be a row ia(g)-1 preceding a sparse segment
 * (where no sparse supercell storage exists; including row 0 case),
 * or an anchor row i0 (where sparse cells may exist because of an
 * arguable design decision). (The first case is the first UP sector
 * in a segment, where we're U only; the second case is all other UP
 * sectors in the segment, where we're D-U.)  No emitting state (M/I)
 * can be reached on this row, and neither can the DL states, but the
 * DG states need consideration, as follows.
 *    
 * The {LG} states that can transition to ML/MG(i1) are on i1-1
 * row. DG states on row i1-1 can be reached by G->DG1..DGk-1->MGk
 * wing unfolding. We only unfold in decoding, so these DG states are
 * unreachable in F/B. We don't store any main supercells for row 0
 * nor for rows ia(g)-1, so we do no unfolding of DG probability
 * there. BUT: because we decided to keep sparse UP supercells on
 * anchor rows i0(d-1), that means we may have supercells on these
 * i1-1 prev rows. Thus for completeness, we decode MGs by wing
 * unfolding in any such cells.
 *     
 * This design decision could be made a different way. Since we never
 * care about DG decoding on rows 0 or ia(g)-1, why care about DG
 * decoding when the UP boundary is i0(d-1)? We could simply not allow
 * any sparse supercells other than the i0,k0 DOWN anchor cell on
 * anchor rows i0. We'd save some supercell storage that way (though
 * probably not much, because anchor cell i0,k0 is usually highly
 * probable by construction, so other supercells on same row are
 * unlikely to appear in sparse mask). No good reason for not doing it
 * this way... I didn't realize it was an issue until late, and I
 * already have it coded with these supercells included, and I'm
 * loathe to code it differently (and possibly discover unexpected
 * kinks).
 */


#endif /*P7_SPASCMX_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
