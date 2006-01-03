#ifndef ALTIVECSTRUCTSH_INCLUDED
#define ALTIVECSTRUCTSH_INCLUDED


/*
 * Structure: logodds_s
 * Date:      CRS, 10 June 2005
 *
 * Purpose:   Used to customize the layout of the data used by the dynamic
 *            programming algorithms.
 */
struct logodds_s {
  /* Note: The Altivec implementation needs the logodds scores aligned on
   *       16-byte boundaries.  In the original code, this was accomplished
   *       by adding _mem pointers, which point to  where the real memory 
   *       is alloc'ed and free'd, as opposed to where it is accessed.  
   *       We move those pointers here, to accomplish this under the new
   *       architecture.
   */
  int  **tsc;                   /* transition scores     [0.6][1.M-1]       -*/
  int  **msc;                   /* match emission scores [0.MAXCODE-1][1.M] -*/
  int  **isc;                   /* ins emission scores [0.MAXCODE-1][1.M-1] -*/
  int    xsc[4][2];             /* N,E,C,J transitions                      -*/
  int   *bsc;                   /* begin transitions     [1.M]              -*/
  int   *esc;			/* end transitions       [1.M]              -*/
  int  *tsc_mem, *msc_mem, *isc_mem, *bsc_mem, *esc_mem; 
};

/*
 * Structure: cust_dpmatrix_s
 * Date:      CRS, 13 July 2005 
 *
 * Purpose:   Stores the calculations of the dynamic programming algorithms.
 *
 */
typedef struct {
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
} cust_dpmatrix_s;


#endif /*ALTIVECSTRUCTSH_INCLUDED*/
