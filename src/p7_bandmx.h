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

#define p7B_NXCELLS 7
#define p7B_E  0
#define p7B_N  1
#define p7B_J  2
#define p7B_B  3
#define p7B_L  4
#define p7B_G  5
#define p7B_C  6

#define p7B_FORWARD  1		
#define p7B_BACKWARD 2		
#define p7B_DECODING 3
#define p7B_ALIGN    4

/* The P7_BANDMX object.
 *                          dp1       dp2
 * p7_BandedForward()       Fwd        -
 * p7_BandedBackward()      Fwd       Bck
 * p7_BandedDecoding()      Fwd       PP
 * p7_BandedAlignment()     OA        PP
 */
typedef struct {
  float     *dp1,  *dp2;	/* main DP cells, for two banded matrices <dp1>, <dp2> */
  float     *xmx1, *xmx2;	/* special DP cells, ditto                             */
  
  int64_t    dalloc;		/* current <dp1,dp2> allocation, denominated in "supercells" (each p7B_NSCELLS wide) */
  int        xalloc;		/* current <xmx1,xmx2> allocation, denominated in banded rows (each p7B_NXCELLS wide) */

  P7_GBANDS *bnd;	        /* reference copy; caller remains responsible for free'ing its bands */
} P7_BANDMX;


extern P7_BANDMX *p7_bandmx_Create (P7_GBANDS *bnd);
extern size_t     p7_bandmx_Sizeof (P7_BANDMX *bmx);
extern int        p7_bandmx_GrowTo (P7_BANDMX *bmx, P7_GBANDS *bnd);
extern int        p7_bandmx_Reuse  (P7_BANDMX *bmx);
extern void       p7_bandmx_Destroy(P7_BANDMX *bmx);


extern char *p7_bandmx_DecodeSpecial(int type);
extern int   p7_bandmx_Dump(FILE *ofp, P7_BANDMX *bmx, int which);

#endif /*P7_BANDMX_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
