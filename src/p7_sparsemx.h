
#ifndef P7_SPARSEMX_INCLUDED
#define P7_SPARSEMX_INCLUDED

typedef struct {
  int      L;		/* sequence has 1..L residues */
  int      M;		/* profile has 1..M positions */
  
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

typedef struct {
  float  *dp;		/* main DP supercells. sm->ncells <= dalloc. each supercell contains p7S_NSCELLS values. */
  float  *xmx;		/* special DP supercells. there are <sm->nrow>+<sm->nseg> of these, each w/ p7S_NXCELLS values. */

  int64_t dalloc;	/* current <dp> allocation, denominated in supercells, (each p7S_NSCELLS wide) */
  int     xalloc;

  P7_SPARSEMASK *sm;
} P7_SPARSEMX;

extern P7_SPARSEMASK *p7_sparsemask_CreateFull(int M, int L);
extern void           p7_sparsemask_Destroy(P7_SPARSEMASK *sm);

extern P7_SPARSEMX   *p7_sparsemx_Create(P7_SPARSEMASK *sm);
extern void           p7_sparsemx_Destroy(P7_SPARSEMX *sx);

extern char *p7_sparsemx_DecodeSpecial(int type);
extern int   p7_sparsemx_Dump(FILE *ofp, P7_SPARSEMX *sx);
extern int   p7_sparsemx_DumpWindow(FILE *ofp, P7_SPARSEMX *sx, int i1, int i2, int ka, int kb);

#endif /*P7_SPARSEMX_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
