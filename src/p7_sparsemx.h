typedef struct {
  int      L;		/* sequence has 1..L residues */
  int      M;		/* profile has 1..M positions */
  
  int    **kin;		/* kin[0,1..L] = ptrs into kmem, rows of sparse k indices; kin[0]=NULL; kin[i]=NULL if nin[i]=0 */
  int     *nin;		/* number of cells included on each row; nin[0]=0; nin[i] <= M */
  int     *kmem;	/* memory that kin[] are pointing into, storing k indices of included cells */

  int      ralloc;	/* kin[] is allocated for ralloc rows; L+1 <= ralloc */
  int64_t  kalloc;	/* kmem[] is allocated for kalloc cells; ncells <= kalloc */

  int64_t  ncells;	/* number of included cells; \sum_{i=0}^{L} nin[i] */
} P7_SPARSEMASK;

