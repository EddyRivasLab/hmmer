#ifndef P7_MASSTRACE_INCLUDED
#define P7_MASSTRACE_INCLUDED


typedef struct {
  float *kmass;	   /* Cumulative mass for ka,kb endpoints on either side of k0. k<=k0: P(ka <= k). k>=k0: P(kb >= k). [k0] = 1. [0]=[M+1]=0. */
  float *imass;	   /* Optional; can be <NULL>. If present, analogous to kmass, for ia,ib sequence endpoints, 0.1..L.L+1 */

  /* The i0,k0,L,M fields are only set when kmass/imass are valid. Otherwise (after Create(), or Reuse()) they're 0. */
  int    i0;	   /* position of anchor on the sequence, 1..L; or 0 */
  int    k0;	   /* position of anchor on the model, 1..M; or 0 */
  int    st0;	   /* state type of anchor */
  int    M;	   /* length of profile; or 0 */
  int    L;	   /* length of seq; or 0 */

  int    kalloc;   /* current allocated size of kmass;  >= M+2 */
  int    ialloc;   /* current allocated size of imass;  >= L+2, if imass is used (non-NULL) */
} P7_MASSTRACE;


/* Recursion in mass trace algorithm works by propagating posterior
 * alignment mass, apportioned by probability of a possible traceback
 * path.  When no paths are possible, the F/B score on a DP cell is
 * -inf, and all paths into that cell were -inf; the mass trace
 * recursion then has rho = 0, pathsc = -inf, and totsc = -inf, which
 * gives us a term of 0 * exp(-inf-inf), which would give us NaN. We
 * have to guard against this; that's what the macro does.
 */
#define P7_MASSTRACE_INCREMENT(rho, pathsc, totsc) ( ((rho) == 0. || (totsc) == -eslINFINITY) ? 0. : (rho) * expf( (pathsc) - (totsc) ))

extern P7_MASSTRACE *p7_masstrace_Create    (int M_hint, int L_hint);
extern P7_MASSTRACE *p7_masstrace_CreateSlim(int M_hint, int L_hint);
extern int           p7_masstrace_GrowTo (P7_MASSTRACE *mt, int M, int L);
extern int           p7_masstrace_Zero   (P7_MASSTRACE *mt, int M, int L);
extern int           p7_masstrace_Reuse  (P7_MASSTRACE *mt);
extern void          p7_masstrace_Destroy(P7_MASSTRACE *mt);

extern int p7_masstrace_CountTrace(const P7_TRACE *tr, int i0, int k0, int st0, P7_MASSTRACE *mt, int *ntr);
extern int p7_masstrace_FinishCount(P7_MASSTRACE *mt, int ntr);

extern int   p7_masstrace_Dump     (FILE *ofp, P7_MASSTRACE *mt);
extern int   p7_masstrace_PlotImass(FILE *ofp, P7_MASSTRACE *mt);
extern int   p7_masstrace_PlotKmass(FILE *ofp, P7_MASSTRACE *mt);
extern int   p7_masstrace_Compare(const P7_MASSTRACE *mte, const P7_MASSTRACE *mta, float tol);
extern float p7_masstrace_GetMaxAbsDiff(const P7_MASSTRACE *mte, const P7_MASSTRACE *mta);
extern int   p7_masstrace_Validate(const P7_MASSTRACE *mt, char *errbuf);


#endif /*P7_MASSTRACE_INCLUDED*/

/* 
 * A "domain" is one pass through the model, a subpath that starts
 * with B and ends with E.
 * Let:
 *   ka = first {MD}k in a domain; e.g. B->L->MLka or B->G->{MD}G1.  (ka=1 for glocal path)
 *   kb = last {MDk} in a domain; e.g. {MD}Lkb->E or {DM}Gm->E.      (kb=M for glocal) 
 *   ia = first emitted residue x_i in domain
 *   ib = last emitted residue x_i in domain
 *   
 *   i0,k0 = anchor of the mass trace; we assign a probability mass of
 *           1.0 that domain includes this point. The "anchor" is generally
 *           the point of highest posterior probability in some domain in the
 *           optimal overall parse of a target seq.
 *           
 * The mass trace algorithm calculates cumulative posterior
 * probability masses for the possible locations of ka,kb,ia,ib:
 * 
 * P(ka <= k) : probability that domain alignment extends at least as far left as start position k on the model;
 * P(kb >= k) : analogous, for the end position
 * P(ia <= i) : analogous, for start position on sequence
 * P(ib >= i) : analogous, for end position on sequence.
 * 
 * P(ka <= k) and P(kb >= k) are nonoverlapping (P(ka <= k) = 1.0 for
 * k>=k0; P(kb >= k) = 1.0 for k <=k0) so we store them both in kmass:
 *   for k<=k0: kmass[k] = P(ka<=k)
 *   for k>=k0: kmass[k] = P(kb>=k)
 * Edge cases: kmass[k0]=1. kmass[0]=kmass[M+1]=0.
 * Analogous for imass[i], cumulative distributions on each side of i0.
 * 
 * Probability that ka start =k: kmass[k]-kmass[k-1], defined for k<=k0  (this is why we store a kmass[0]=0 sentinel)
 * Probability that kb end   =k: kmass[k]-kmass[k+1], defined for k>=k0  ( and the kmass[M+1]=0 sentinel )
 * and analogously for ia,ib starts.
 * 
 * The "envelope" is defined by four start/end coords on the seq and model:
 *   iae = \min_{i \leq i0} P(ia \leq i) \geq massthresh
 *   ibe = \max_{i \geq i0} P(ib \geq i) \geq massthresh
 *   
 *   kae = \min_{k \leq k0} P(ka \leq k) \leq massthresh
 *   kbe = \max_{k \geq k0} P(kb \geq k) \leq massthresh
 * i.e. the widest start/end points that contain at least a
 * posterior probability of <massthresh>. (The probability that
 * the domain starts at some i < iae is less than <massthresh>.)

 * Normally (in production code), to save time, we do not calculate or
 * store the complete cumulative distributions. The Up and Down
 * recursions may terminate at some row i>1 (Up) or i<L (Down).  In
 * this production mode, <imass> storage is not needed at all; the Up
 * and Down algorithms only need to store one current mass for current
 * row i, recursing upwards and downward in rows i until this current
 * rowmass ceases to satisfy the within-envelope threshold
 * conditions. <kmass> is needed as tmp workspace to accumulate path
 * probabilities on each row i.  Because not all i may be calculated,
 * now kmass[k] is only a bound: for k < k0: kmass[k] <= P(ka<=k) for
 * k > k0: kmass[k] <= P(kb>=k). 
 * 
 * (Note that it's not quite as simple as calculating rows i until iae
 * or ibe has been identified; more rows i<iae and/or i>ibe may also need
 * to be calculated until kae,kbe have also been identified.)
 *   
 * In some cases (development code, unit tests, statistics gathering)
 * we may want completely calculated cumulative distributions, perhaps
 * because we're going to plot them, or because we're going to compare
 * them to the same cumulative distributions generated by a second
 * independent means (as in unit testing), as opposed to merely
 * testing the envelope mass threshold condition. Having a non-NULL
 * <imass> signals to the Up() and Down() recursions that they need to
 * calculate all rows i, so imass[] and kmass[] are complete.
 * 
 * The <P7_MASSTRACE> object stores the <imass> and <kmass>
 * distributions for use by the Up() and Down() recursions. 
 * 
 * If <imass> is <NULL> the recursions may truncate their calculations
 * after identifying the envelope coords; <kmass> is best thought of
 * as a temporary workspace, and the object's function is primarily to
 * manage the allocation of <kmass>, minimizing reallocation by
 * reusing it where possible.
 * 
 * If <imass> is <non-NULL>, the recursions will complete their
 * calculations to the full set of rows 1..L, in order to collect the
 * complete <imass> and <kmass> cumulative distributions.  The
 * object's function is not only to manage the allocations of <imass>
 * and <kmass>, but also to store these distributions for subsequent
 * analysis. In unit tests, alternative means of collecting these
 * distributions (for example by a stochastic trace ensemble) are
 * used, and P7_MASSTRACE objects are compared.
 */
