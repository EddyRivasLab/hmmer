/* P7_PROFILE: a scoring profile, and its implicit model.
 */
#ifndef p7PROFILE_INCLUDED
#define p7PROFILE_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include "esl_alphabet.h"

#include "base/general.h"

/* Indices for six special state types x that have transition parameters
 * in gm->xsc[x][y], where all of them have two choices y.
 */
#define p7P_NXSTATES 6
#define p7P_NXTRANS  2
#define p7P_LOOP 0              /* gm->xsc[x][p7P_LOOP]    gm->xsc[x][p7P_MOVE] */        
#define p7P_MOVE 1              /* -------------------     -------------------- */
#define p7P_E  0            	/*     E->J                       E->C          */
#define p7P_N  1		/*     N->N                       N->B          */
#define p7P_J  2                /*     J->J                       J->B          */
#define p7P_C  3		/*     C->C                       C->T          */
#define p7P_B  4		/*     B->L                       B->G          */
#define p7P_G  5		/*     G->M1                      G->D1         */

/* Indices x for main model transition scores gm->tsc[k][x] 
 * Order is optimized for dynamic programming algorithms.
 */
#define p7P_NTRANS 10
#define p7P_MM   0
#define p7P_IM   1
#define p7P_DM   2
#define p7P_LM   3	/* local submodel entry L->Mk; stored off-by-one, tsc[k-1][LM] = L->Mk */
#define p7P_GM   4	/* wing-retracted glocal submodel entry G->Mk;    tsc[k-1][GM] = G->Mk */
#define p7P_MD   5  
#define p7P_DD   6
#define p7P_MI   7
#define p7P_II   8
#define p7P_DGE  9	/* wing-retracted glocal exit, DD component, tDGEk = Dk+1..Dm->E, 0.0 for k=M-1,M */

/* Indices for residue emission score vectors */
#define p7P_NR   2
#define p7P_M    0
#define p7P_I    1

/* A brief cheat sheet on the boundary conditions / edge cases of the
 * profile's transition and emission score arrays:
 * 
 *               k = 0                             k = 1                                   k = M-1                             k = M
 * tsc: [ MM IM DM LM GM MD DD MI II DGE ] [ MM IM DM LM GM MD DD MI II DGE ] ... [ MM IM DM LM GM MD DD MI II DGE ] [ MM IM DM LM GM MD DD MI II DGE ] 
 *         *  *  *  .  .  *  *  *  *   *      .  .  .  .  .  .  .  .  .   .          .  .  .  .  .  .  .  .  .   0      *  *  *  *  *  0  0  *  *   0
 *                  ^  ^ {LG}->Mk are stored off-by-one. These are {LG}->M1, for example
 *                  
 *        k=0 -inf: p7_profile_Create(). D_0 does not exist. I_0 exists in core HMM, but is removed from profile.
 *        k=M -inf: p7_profile_ConfigCustom().
 *        DGE_M 0:  modelconfig.c::set_glocal_exit().
 *        DGE_M-1 0: modelconfig.c::set_glocal_exit().
 *        DGE_0:  That would be the D1->Dm->E path, which is valid, but mute; we elide it by leaving DGE0 to -inf
 *                  
 *              k = 0     k = 1      k = M-1    k = M
 * rsc[x]:     [ M  I ] [ M  I ] ... [ M  I ] [ M  I ] 
 *  [gap]:       *  *     *  *         *  *     *  *
 *  [missing]:   *  *     *  *         *  *     *  * 
 *  [residues]:  *  *     .  0         .  0     .  *
 * 
 *        k=0 -inf: p7_profile_Create(). M_0, I_0 are nonexistent.
 *        all k -inf for gap, missing residue: p7_profile_Create()
 *        insert emissions (1..M-1) hardwired to 0: modelconfig.c::p7_profile_ConfigCustom()
 *        k=M insert emission -inf: modelconfig.c::p7_profile_ConfigCustom()
 */
typedef struct p7_profile_s {
  /* Model parameters:                                                               */
  int     M;		/* number of nodes in the model                              */
  float  *tsc;          /* transitions  [0.1..M][0..p7P_NTRANS-1], hand-indexed      */
  float **rsc;          /* emissions [0..Kp-1][0.1..M][p7P_NR], hand-indexed         */
  float   xsc[p7P_NXSTATES][p7P_NXTRANS]; /* special transitions [ENJCBG][LOOP,MOVE] */

  /* Memory allocation:                                                              */
  int     allocM;	/* max # of nodes allocated in this structure                */

  /* Configuration: length model, multi vs unihit, local vs glocal:                  */
  int     L;		/* current configured target seq length           (unset:-1) */
  float   nj;           /* exp # of J's; 0.0=unihit 1.0=standard multihit (unset:-1) */
  float   pglocal;	/* base B->G; 0.0=local; 0.5=dual; 1.0=glocal     (unset:-1) */

  /* Annotation copied from parent HMM:                                                   */
  char  *name;			/* unique name of model                                   */
  char  *acc;			/* unique accession of model, or NULL                     */
  char  *desc;                  /* brief (1-line) description of model, or NULL           */
  char  *rf;                    /* reference line from alignment 1..M; *rf=0 means unused */
  char  *mm;                    /* modelmask line           1..M; *ref=0: unused          */
  char  *cs;                    /* consensus structure line 1..M, *cs=0 means unused      */
  char  *consensus;		/* consensus residues to display in alignments, 1..M      */
  float  evparam[p7_NEVPARAM]; 	/* parameters for determining E-values, or UNSET          */
  float  cutoff[p7_NCUTOFFS]; 	/* per-seq/per-domain bit score cutoffs, or UNSET         */
  float  compo[p7_MAXABET];	/* per-model HMM filter composition, or UNSET             */

  /* Disk offset information supporting fast model retrieval:                             */
  off_t  offs[p7_NOFFSETS];     /* p7_{MFP}OFFSET, or -1                                  */
  off_t  roff;                  /* record offset (start of record); -1 if none            */
  off_t  eoff;                  /* offset to last byte of record; -1 if unknown           */

  /* Derived information:                                                                 */
  int     max_length;	/* calc'ed upper bound on emitted seq length (nhmmer) (unset:-1)  */

  /* Associated objects:                                                                  */
  const ESL_ALPHABET *abc;	/* copy of pointer to appropriate alphabet                */
} P7_PROFILE;

/* Convenience macros for accessing transition, emission scores */
/* _LM,GM are specially stored off-by-one: [k-1][p7P_{LG}M] is score for *entering* at Mk */
#define P7P_TSC(gm, k, s) ((gm)->tsc[(k) * p7P_NTRANS + (s)])
#define P7P_MSC(gm, k, x) ((gm)->rsc[x][(k) * p7P_NR + p7P_M])
#define P7P_ISC(gm, k, x) ((gm)->rsc[x][(k) * p7P_NR + p7P_I])



extern P7_PROFILE *p7_profile_Create(int M, const ESL_ALPHABET *abc);
extern P7_PROFILE *p7_profile_Clone(const P7_PROFILE *gm);
extern int         p7_profile_Copy(const P7_PROFILE *src, P7_PROFILE *dst);
extern int         p7_profile_SetNullEmissions(P7_PROFILE *gm);
extern int         p7_profile_Reuse(P7_PROFILE *gm);
extern size_t      p7_profile_Sizeof(P7_PROFILE *gm);
extern void        p7_profile_Destroy(P7_PROFILE *gm);
extern int         p7_profile_IsLocal   (const P7_PROFILE *gm);
extern int         p7_profile_IsGlocal  (const P7_PROFILE *gm);
extern int         p7_profile_IsMultihit(const P7_PROFILE *gm);
extern float       p7_profile_GetT(const P7_PROFILE *gm, char st1, int k1, char st2, int k2);
extern int         p7_profile_Dump(FILE *fp, P7_PROFILE *gm);
extern int         p7_profile_Validate(const P7_PROFILE *gm, char *errbuf, float tol);
extern char       *p7_profile_DecodeT(int tidx);
extern int         p7_profile_GetMutePathLogProb(const P7_PROFILE *gm, double *ret_mute_lnp);
extern int         p7_profile_Compare(P7_PROFILE *gm1, P7_PROFILE *gm2, float tol);


/* Whereas the core HMM is a model of a single global alignment to a
 * homologous domain, a profile is a model of local and glocal homologous
 * alignments embedded in a longer sequence. The profile is constructed
 * from the core HMM by adding nine additional states: 
 *
 *   S->N->B-->L ... -->E->C->T
 *          \_>G ... _/ |
 *         ^            |
 *         |____ J _____|
 *
 * These nine states control the "algorithm-dependent" configuration
 * of the profile and of our alignment/scoring algorithms (Viterbi or 
 * Forward). There are three components to this configuration:
 *     1. The length model:
 *          N,C,J transitions control the expected length of
 *          nonhomologous sequence segments. The length model is
 *          typically reset for each new target sequence length L.
 *     2. Multihit vs. unihit
 *          The E transitions control whether the model is allowed
 *          to loop around for multiple homology segments per
 *          target sequence. Usually we're multihit.
 *     3. Search mode: local vs. glocal paths
 *          The B->L, B->G probabilities control whether the model
 *          only allows local alignment, only allows glocal alignment,
 *          or a mixture of the two. Usually we're "dual-mode"
 *          with 50/50 probability of local vs. glocal.
 *         
 * For efficiency, two of the nine states (S,T) only appear in our
 * formal state diagrams. They are not explicitly represented in the
 * P7_PROFILE parameters that DP calculations use:
 *     1. The S->N probability is always 1.0, so our DP codes 
 *        initialize on the N state and don't store anything for S.
 *     2. The C->T log probability is added on to any final 
 *        lod score. We know the T state only has nonzero probability
 *        at position i=L in the target, so DP algorithms don't
 *        bother to store it in their matrix.
 *        
 * The profile is missing some states of the core HMM:
 *     1. The D1 state is removed by "wing retraction".
 *     2. The I0 and Im states aren't present.
 * Transition costs from the I0/Im states are initialized to -inf.
 * Transition costs from the D1 state are present in the P7_PROFILE
 * data structure, but unused by DP algorithms; see the "wing 
 * retraction" section below for an explanation.
 *     
 * Some alignment/scoring DP implementations may implicitly override the
 * length model, multihit mode, or the search mode (so long as they
 * document it). Examples:
 *     1. For the length model, we sometimes make the approximation
 *          that there are ~L NN/CC/JJ transitions (assuming the
 *          target sequence length is large enough that the sum 
 *          of homologous segment lengths is negligible). This
 *          allows what the code calls the "3 nat" or "2 nat"
 *          approximation: L \log \frac{L}{L+3} ~ 3.0 (multihit mode),
 *          L \log \frac{L}{L+2} ~ 2.0 (unihit mode).
 *     2. To override multihit/unihit mode: Multihit vs. unihit mode
 *          also affects the length model parameters, so the
 *          algorithm must not only replace tEC/tEJ parameters, but
 *          must also use appropriate N,C,J parameters. The MSV
 *          filter (impl_xxx/msvfilter.c) is an example where the
 *          algorithm is always multihit, and where NN,CC,JJ 
 *          parameters are handled implicitly w/ the 3-nat approximation.
 *     3. To override search mode and use only one path (local or
 *          glocal): The algorithm ignores the tBL/tBG log probability
 *          and evaluates the path it wants as if it had B->L or B->G
 *          set to 1.0.
 *          
 * The full (dual-mode, local/glocal) model is only implemented in
 * P7_PROFILE, not in a vectorized P7_OPROFILE. Vector implementations
 * are always and only in local alignment mode, because of numeric
 * range limitations. Vector implementations are only used for fast
 * filtering before a final banded calculation with the full model.
 *          
 * Wing retraction:        
 *   For the local path, all algorithms must use the L->Mk "uniform"
 *   entry distribution, and assume that all Mk->E and Dk->E are
 *   1.0. The L->Mk entry distribution supersedes all paths involving
 *   terminal deletion. This is the "implicit probabilistic model"
 *   of local alignment [Eddy08].  The D_1 and I_M states don't exist
 *   on the local path, and algorithms must treat them as zero probability.
 * 
 *   For the glocal path, all algorithms must use a G->Mk "wing
 *   retracted" entry distribution <p7P_GM>, and some algorithms may
 *   optionally use a DGk+1->E wing retracted exit distribution
 *   <p7P_DGE>. Wing retracted entry G->Mk is required, because this is
 *   the way we remove a mute cycle (B->G->D1..Dm->E->J->B) from the
 *   profile: by removing the D1 state and neglecting the probability
 *   of the mute cycle in all DP calculations.
 *      
 *   Wing retracted entry, see modelconfig.c::set_glocal_entry()
 *      tGM1 = log t(G->M1) 
 *      tGMk = log t(G->D1) + \sum_j={1..k-2} log t(Dj->Dj+1) + log t(Dk-1->Mk)
 *      stored off-by-one: tGMk is stored at TSC(p7P_GM, k-1) in profile structure.
 *      
 *   Wing retracted exit, see modelconfig.c::set_glocal_exit()
 *      tDGkE = log t(Dk+1->...Dm->E)
 *            = \sum_j={k+1..m-1} log t(Dj->Dj+1)    (recall that Dm->E = 1.0)
 *      valid for k=0..M: 
 *      boundary conditions:
 *      TSC(DGE,M) = TSC(DGE,M-1) = 0    
 *   
 *   Wing retracted exits are used in sparse DP calculations: see
 *   sparse_fwdback.c. A DP calculation may also use Mm->E and Dm->E
 *   (both 1.0) explicitly, rather than using the DGE wing retracted
 *   exit.  The reference implementation does this, for example (see
 *   reference_fwdback.c or reference_viterbi.c).
 *
 *   Wing retraction is only used internally in DP implementations.
 *   Any resulting tracebacks (alignments) still represent a
 *   G->D1...Dk-1->Mk entry explicitly; you'll see this in traceback
 *   code, where it adds k-1 delete states to a trace whenever it sees
 *   a G->Mk entry.  Thus we may still need the D1->{DM}
 *   parameterization, to call p7_trace_Score() for example; and so,
 *   even though the D1 state is implicitly removed by all DP
 *   implementations, on both the local and glocal paths, its
 *   parameterization is still present in the P7_PROFILE.
 */

#endif /*p7PROFILE_INCLUDED*/

