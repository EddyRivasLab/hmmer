#ifndef DEFAULTSTRUCTSH_INCLUDED
#define DEFAULTSTRUCTSH_INCLUDED

struct logodds_s {
  /* Note: We have to provide a definition for this structure, but
   *       the default implementation doesn't use it, so we will
   *	   just leave it blank. - CRS 10 June 2005
   */
};

/* Declaration of Plan7 dynamic programming matrix structure.
 *
 * Note:  This was originally defined in structs.h, but this structure is
 *        really implementation dependent, so I moved it here.  Since the
 *        default implementation doesn't use the _mem pointers, I 
 *        commented those parts out.
 *          - CRS 21 June 2005
 */
struct dpmatrix_s {
  int **xmx;			/* special scores [0.1..N][BECJN]     */
  int **mmx;			/* match scores [0.1..N][0.1..M]      */
  int **imx;			/* insert scores [0.1..N][0.1..M-1.M] */
  int **dmx;			/* delete scores [0.1..N][0.1..M-1.M] */

  /* Hidden ptrs where the real memory is kept; this trick was
   * introduced by Erik Lindahl with the Altivec port; it's used to
   * align xmx, etc. on 16-byte boundaries for cache optimization.
   */
  void *xmx_mem, *mmx_mem, *imx_mem, *dmx_mem;

  int *  workspace;      /* Workspace for altivec (aligned ptr)    */
  int *  workspace_mem;  /* Actual allocated pointer for workspace */
  
  /* The other trick brought in w/ the Lindahl Altivec port; dp matrix
   * is retained and grown, rather than reallocated for every HMM or sequence.
   * Keep track of current allocated-for size in rows (sequence length N)
   * and columns (HMM length M). Also keep track of pad sizes: how much
   * we should overallocate rows or columns when we reallocate. If pad = 0,
   * then we're not growable in this dimension.
   */
  int maxN;			/* alloc'ed for seq of length N; N+1 rows */
  int maxM;			/* alloc'ed for HMM of length M; M+1 cols */

  int padN;			/* extra pad in sequence length/rows */
  int padM;			/* extra pad in HMM length/columns   */
};

#endif
