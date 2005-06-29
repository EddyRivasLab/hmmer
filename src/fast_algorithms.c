/* fast_algorithms.c
 * Optimized routines to replace slower implementations in core_algorithms.c.
 * 
 * The routines in core_algorithms.c are designed for clarity
 * and maintainability, not for speed. Implementations here
 * are designed for speed, not clarity. If you're trying to 
 * understand the code, or optimize for a specific platform,
 * you are probably better off looking at core_algorithms.c.
 * 
 * P7Viterbi() is the key function to target optimization to.
 * The implementation in core_algorithms.c is currently ifdef'ed 
 * out of the code. The implementation that is used by default 
 * is here, in fast_algorithms.c. A third implementation, from
 * Erik Lindahl at Stanford, is Mac/Altivec specific.
 * 
 * Which implementation is used is controlled by ifdef's. The
 * default implementation uses a fast implementation of 
 * P7Viterbi() from here. Other options (mutually exclusive):
 * 
 * -DSLOW
 *   enable original core_algorithms.c code: slower than default,
 *   but might be easier to follow, for someone trying
 *   to understand the DP code.
 * -DALTIVEC
 *   enable Erik Lindahl's Altivec code for Macintosh OSX
 *
 * SRE, Sun Nov 10 08:54:48 2002 [AA 3080, Denver to StL]
 * SVN $Id$
 */

#include "config.h"
#include "squidconf.h"

#include "squid.h"

#include "plan7.h"
#include "structs.h"
#include "funcs.h"

#if (defined __GNUC__) && (defined __APPLE__)
#include <ppc_intrinsics.h>
#endif

/*
 * Note: This function is no longer needed, since the architecture for
 *       specifing customized implementations has changed.  Thus, I 
 *	 am commenting this code out.  - CRS 21 June 2005
 *
#ifdef ALTIVEC

/\*################################################################
 * The Altivec port, for Macintosh PowerPC.
 * Erik Lindahl, Stanford University, 2002.
 * 
 * Replaces the following functions:
 *   AllocPlan7Body()      plan7.c               (data alignment on 16-byte boundaries)
 *   CreatePlan7Matrix()   core_algorithms.c     (data alignment on 16-byte boundaries)
 *   ResizePlan7Matrix()   core_algorithms.c     (data alignment on 16-byte boundaries)
 *   P7Viterbi()           core_algorithms.c     (total recode, w/ Altivec instructions)
 ################################################################*\/    

void
AllocPlan7Body(struct plan7_s *hmm, int M) 
{
  int k, x;

  hmm->M = M;

  hmm->rf     = MallocOrDie ((M+2) * sizeof(char));
  hmm->cs     = MallocOrDie ((M+2) * sizeof(char));
  hmm->ca     = MallocOrDie ((M+2) * sizeof(char));
  hmm->map    = MallocOrDie ((M+1) * sizeof(int));

  hmm->t      = MallocOrDie (M     *           sizeof(float *));
  hmm->tsc    = MallocOrDie (7     *           sizeof(int *));
  hmm->mat    = MallocOrDie ((M+1) *           sizeof(float *));
  hmm->ins    = MallocOrDie (M     *           sizeof(float *));
  hmm->msc    = MallocOrDie (MAXCODE   *       sizeof(int *));
  hmm->isc    = MallocOrDie (MAXCODE   *       sizeof(int *)); 
  
  hmm->t[0]   = MallocOrDie ((7*M)     *       sizeof(float));
  /\* Allocate extra memory so tsc[TMM,TIM,TDM,TMD,TDD] start on the 
   * 16-byte cache boundary, and tsc[TMI,TII] start
   * 12 bytes offset from the boundary. 
   *\/
  hmm->tsc_mem = MallocOrDie (((7*(M+16)))  *   sizeof(int));
  hmm->mat[0] = MallocOrDie ((MAXABET*(M+1)) * sizeof(float));
  hmm->ins[0] = MallocOrDie ((MAXABET*M) *     sizeof(float));
  /\* Allocate extra mem. to make sure all members of msc,isc start
   * on 12-byte offsets from cache boundary.
   *\/
  hmm->msc_mem = MallocOrDie ((MAXCODE*(M+1+16)) * sizeof(int));
  hmm->isc_mem = MallocOrDie ((MAXCODE*(M+16)) *   sizeof(int));

  /\* note allocation strategy for important 2D arrays -- trying
   * to keep locality as much as possible, cache efficiency etc.
   *\/
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * MAXABET;
    if (k < M) {
      hmm->ins[k] = hmm->ins[0] + k * MAXABET;
      hmm->t[k]   = hmm->t[0]   + k * 7;
    }
  }
  
  /\* align tsc pointers *\/
  hmm->tsc[TMM] = (int *) (((((size_t) hmm->tsc_mem) + 15) & (~0xf)));
  hmm->tsc[TMI] = (int *) (((((size_t) hmm->tsc_mem) + (M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  hmm->tsc[TMD] = (int *) (((((size_t) hmm->tsc_mem) + 2*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TIM] = (int *) (((((size_t) hmm->tsc_mem) + 3*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TII] = (int *) (((((size_t) hmm->tsc_mem) + 4*(M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  hmm->tsc[TDM] = (int *) (((((size_t) hmm->tsc_mem) + 5*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TDD] = (int *) (((((size_t) hmm->tsc_mem) + 6*(M+12)*sizeof(int) + 15) & (~0xf)));

  for (x = 0; x < MAXCODE; x++) {
    hmm->msc[x] = (int *) (((((size_t)hmm->msc_mem) + x*(M+1+12)*sizeof(int) + 15) & (~0xf)) + 12);
    hmm->isc[x] = (int *) (((((size_t)hmm->isc_mem) + x*(M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  }
  /\* tsc[0] is used as a boundary condition sometimes [Viterbi()],
   * so set to -inf always.
   *\/
  for (x = 0; x < 7; x++)
    hmm->tsc[x][0] = -INFTY;

  hmm->begin  = MallocOrDie  ((M+1) * sizeof(float));
  hmm->bsc_mem= MallocOrDie  ((M+1+12) * sizeof(int));
  hmm->end    = MallocOrDie  ((M+1) * sizeof(float));
  hmm->esc_mem= MallocOrDie  ((M+1+12) * sizeof(int));

  hmm->bsc = (int *) (((((size_t) hmm->bsc_mem) + 15) & (~0xf)) + 12);
  hmm->esc = (int *) (((((size_t) hmm->esc_mem) + 15) & (~0xf)) + 12);
  
  return;
}  

#endif /\*the ALTIVEC port*\/
*
*
*/


/************************************************************
 * @LICENSE@
 ************************************************************/

