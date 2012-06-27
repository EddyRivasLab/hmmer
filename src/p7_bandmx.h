/* P7_BANDMX is the banded DP matrix used for banded Forward/Backward,
 * posterior decoding, and optimal accuracy alignment, using a dual-mode
 * local/glocal model.
 * 
 * Contents:
 *    1. The P7_BANDMX structure
 *    2. Function declarations
 *    x. Copyright and license information
 */
#ifndef P7_BANDMX_INCLUDED
#define P7_BANDMX_INCLUDED

#include "p7_gbands.h"

#define p7B_NSCELLS 6
#define p7B_ML 0
#define p7B_MG 1
#define p7B_IL 2
#define p7B_IG 3
#define p7B_DL 4
#define p7B_DG 5

#define p7B_NXCELLS 9
#define p7B_E  0
#define p7B_N  1
#define p7B_J  2
#define p7B_B  3
#define p7B_L  4
#define p7B_G  5
#define p7B_C  6
#define p7B_JJ 7	/* in decoding (only) we separately decode J occupancy vs JJ emission */
#define p7B_CC 8	/* ditto for C */

/* The P7_BANDMX object.
 */
typedef struct {
  float     *dp;      /* main DP supercells; there are <bnd->ncell> of these, each containing p7B_NSCELLS values     */
  float     *xmx;     /* special DP supercells; there are <bnd->nrow+bnd->nseg> of these, each w/ p7B_NXCELLS values */
  
  int64_t    dalloc;  /* current <dp> allocation, denominated in "supercells" (each p7B_NSCELLS wide) */
  int        xalloc;  /* current <xmx> allocation, denominated in banded rows (each p7B_NXCELLS wide) */

  P7_GBANDS *bnd;     /* reference copy; caller remains responsible for free'ing its bands */
} P7_BANDMX;


extern P7_BANDMX *p7_bandmx_Create (P7_GBANDS *bnd);
extern size_t     p7_bandmx_Sizeof (P7_BANDMX *bmx);
extern int        p7_bandmx_Reinit (P7_BANDMX *bmx, P7_GBANDS *bnd);
extern int        p7_bandmx_Reuse  (P7_BANDMX *bmx);
extern void       p7_bandmx_Destroy(P7_BANDMX *bmx);

extern char *p7_bandmx_DecodeSpecial(int type);
extern int   p7_bandmx_Dump      (FILE *ofp, P7_BANDMX *bmx);
extern int   p7_bandmx_DumpWindow(FILE *ofp, P7_BANDMX *bmx, int istart, int iend, int kstart, int kend);

/* from banded_fwdback.c */
extern int p7_BandedForward (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_BANDMX *bmf, float *opt_sc);
extern int p7_BandedBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_BANDMX *bmb, float *opt_sc);
extern int p7_BandedDecoding(const P7_PROFILE *gm, const P7_BANDMX *bmf, P7_BANDMX *bmb, P7_BANDMX *bmd);
extern int p7_BandedAlign   (const P7_PROFILE *gm, float gamma, const P7_BANDMX *bmd, P7_BANDMX *bma, P7_TRACE *tr, float *opt_gain);
#endif /*P7_BANDMX_INCLUDED*/

/*****************************************************************
 * 3. Notes
 *****************************************************************/

/* See notes in p7_gbands.h for layout of bands, and idioms for
 * traversing bands (same idioms apply for traversing a banded
 * matrix).
 * 
 * For each banded segment i..j, we store an extra xmx row i-1.  This
 * allows us to do posterior decoding on {NJBLGC,JJ,CC}(i-1), values
 * that we may want. For example, the posterior probability that a
 * domain starts at residue i is in B(i-1); specifically that a local
 * vs glocal alignment starts at i is in L(i-1), G(i-1) respectively;
 * and the single per-domain score for one domain within i..j can be
 * estimated as a function of differencing N/C/J(i-1) against N/C/J(j).
 * xref J9/128-130.
 */

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
