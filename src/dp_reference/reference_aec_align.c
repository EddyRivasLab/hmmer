/* Reference implementation of anchor/envelope constrained (AEC) alignment;
 * by maximum expected gain (MEG).
 * 
 * Reference implementation is for development and testing. Not used
 * in HMMER's main executables. Production code uses sparse dynamic
 * programming; see src/dp_sparse.
 * 
 * Contents:
 *    1. AEC MEG alignment
 *    2. Choice functions for traceback
 *    3. Traceback routine
 *    4. Footnotes
 *    5. Unit tests.
 *    6. Test driver
 *    7. Example
 */
#include <p7_config.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "base/p7_envelopes.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_aec_align.h"

static int reference_aec_trace(const P7_PROFILE *gm, P7_ENVELOPES *env, const P7_REFMX *mx, P7_TRACE *tr);


/*****************************************************************
 * 1. AEC MEG alignment, fill.
 *****************************************************************/

/* Function:  p7_reference_AEC_Align()
 * Synopsis:  Maximum expected gain, anchor/envelope-constrained alignment.
 *
 * Purpose:   Calculates an anchor/envelope-constrained (AEC) maximum expected
 *            gain (MEG) alignment, given profile <gm>, envelope data <env>,
 *            and ASC posterior decoding matrix UP/DOWN pair <apu>/<apd>.
 *
 *            Caller provides <mx> as space for the DP
 *            calculation. This will be reallocated if needed, so it
 *            can be of any allocated size.
 *
 * Args:      gm   : profile
 *            env  : envelopes 
 *            apu  : ASC Decoding UP matrix
 *            apd  : ASC Decoding DOWN matrix
 *            mx   : MEG alignment matrix space
 *            tr   : RETURN: MEG path
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_reference_AEC_Align(const P7_PROFILE *gm, P7_ENVELOPES *env, const P7_REFMX *apu, const P7_REFMX *apd, P7_REFMX *mx, P7_TRACE *tr)
{
  float *tsc = gm->tsc;      // TSC() access macro now works
  int    M   = env->M;
  int    L   = env->L;
  int    d,i,k;
  int    ia,i0,k0,ib;        // tmp vars for simplifying notation.
  int    is_glocal;
  float  xX;                 // AEC only needs one special, and doesn't need to store it.
  int    status;

  /* Contract checks / argument validation */
  ESL_DASSERT1(( gm->M == M ));
  ESL_DASSERT1(( apu->M == M && apu->L == L && apu->type == p7R_ASC_DECODE_UP));
  ESL_DASSERT1(( apd->M == M && apd->L == L && apd->type == p7R_ASC_DECODE_DOWN));
  
  /* Reallocation if needed */
  if ( (status = p7_refmx_GrowTo   (mx, M, L))                != eslOK) return status;
  if ( (status = p7_refmx_SetType  (mx, M, L, p7R_AEC_ALIGN)) != eslOK) return status;
  if ( (status = p7_refmx_SetValues(mx, -eslINFINITY))        != eslOK) return status; 

  xX = 0.;
  for (d = 1; d <= env->D; d++)
    {
      /* Simplify notation with some tmp vars for the coords of this domain <d>. */
      ia        = env->arr[d].ia; 
      i0        = env->arr[d].i0;
      k0        = env->arr[d].k0;
      ib        = env->arr[d].ib;
      is_glocal = (env->arr[d].flags & p7E_IS_GLOCAL);
      
      /* Initialization of row ia(d), if UP sector exists (usually does) */
      if (ia < i0)
	{
	  for (k = 1; k < k0; k++)
	    {
	      P7R_MX(mx,ia,k,p7R_ML) = P7R_MX(apu,ia,k,p7R_ML) + P7R_MX(apu,ia,k,p7R_MG) + 
	 	                       (is_glocal ? P7_DELTAT(xX, TSC(p7P_GM, k-1)) :
					            P7_DELTAT(xX, TSC(p7P_LM, k-1)));
	      P7R_MX(mx,ia,k,p7R_IL) = -eslINFINITY;              
	      P7R_MX(mx,ia,k,p7R_DL) = ESL_MAX( P7_DELTAT(P7R_MX(mx,ia,k-1,p7R_ML), TSC(p7P_MD,k-1)),
						P7_DELTAT(P7R_MX(mx,ia,k-1,p7R_DL), TSC(p7P_DD,k-1)));
	    }
	}

      /* Remaining rows of UP sector (if any) */
      for (i = ia+1; i < i0; i++)
	for (k = 1; k < k0; k++)
	  {
	    P7R_MX(mx,i,k,p7R_ML) = P7R_MX(apu,i,k,p7R_ML) + P7R_MX(apu,i,k,p7R_MG) + 
	                            ESL_MAX( ESL_MAX( P7_DELTAT( P7R_MX(mx,i-1,k-1,p7R_ML), TSC(p7P_MM,k-1)),
				                      P7_DELTAT( P7R_MX(mx,i-1,k-1,p7R_IL), TSC(p7P_IM,k-1))),
		                                      P7_DELTAT( P7R_MX(mx,i-1,k-1,p7R_DL), TSC(p7P_DM,k-1)));
	    P7R_MX(mx,i,k,p7R_IL) = P7R_MX(apu,i,k,p7R_IL) + P7R_MX(apu,i,k,p7R_IG) + 
	                            ESL_MAX( P7_DELTAT( P7R_MX(mx,i-1,k,p7R_ML), TSC(p7P_MI,k)),
				             P7_DELTAT( P7R_MX(mx,i-1,k,p7R_IL), TSC(p7P_II,k)));
	    P7R_MX(mx,i,k,p7R_DL) = ESL_MAX( P7_DELTAT( P7R_MX(mx,i,k-1,p7R_ML), TSC(p7P_MD,k-1)),
					     P7_DELTAT( P7R_MX(mx,i,k-1,p7R_DL), TSC(p7P_DD,k-1)));
	  }
	      
      /* Anchor cell i0(d),k0(d); usually initialized from last
       * supercell of UP sector, but rarely UP sector may not
       * exist. 
       */
      if (ia < i0) P7R_MX(mx,i0,k0,p7R_ML) = P7R_MX(apd,i0,k0,p7R_ML) + P7R_MX(apd,i0,k0,p7R_MG) + 
		                             ESL_MAX( ESL_MAX( P7_DELTAT( P7R_MX(mx,i0-1,k0-1,p7R_ML), TSC(p7P_MM,k0-1)),
				                               P7_DELTAT( P7R_MX(mx,i0-1,k0-1,p7R_IL), TSC(p7P_IM,k0-1))),
		                                               P7_DELTAT( P7R_MX(mx,i0-1,k0-1,p7R_DL), TSC(p7P_DM,k0-1)));
      else         P7R_MX(mx,i0,k0,p7R_ML) = P7R_MX(apd,i0,k0,p7R_ML) + P7R_MX(apd,i0,k0,p7R_MG) +
	 	                             (is_glocal ? P7_DELTAT(xX, TSC(p7P_GM, k0-1)) :
					                  P7_DELTAT(xX, TSC(p7P_LM, k0-1)));
      P7R_MX(mx,i0,k0,p7R_IL) = -eslINFINITY;
      P7R_MX(mx,i0,k0,p7R_DL) = -eslINFINITY;

      /* Remainder of the top DOWN row, i0 for k>k0 */
      for (k = k0+1; k <= M; k++)
	{
	  P7R_MX(mx,i0,k,p7R_ML) = -eslINFINITY;
	  P7R_MX(mx,i0,k,p7R_IL) = -eslINFINITY;
	  P7R_MX(mx,i0,k,p7R_DL) = ESL_MAX( P7_DELTAT( P7R_MX(mx,i0,k-1,p7R_ML), TSC(p7P_MD,k-1)),
				            P7_DELTAT( P7R_MX(mx,i0,k-1,p7R_DL), TSC(p7P_DD,k-1)));
	}

      /* DOWN sector recursion */
      for (i = i0+1; i <= ib; i++)
	for (k = k0; k <= M; k++)
	  {
	    P7R_MX(mx,i,k,p7R_ML) = P7R_MX(apd,i,k,p7R_ML) + P7R_MX(apd,i,k,p7R_MG) + 
	                            ESL_MAX( ESL_MAX( P7_DELTAT( P7R_MX(mx,i-1,k-1,p7R_ML), TSC(p7P_MM,k-1)),
				                      P7_DELTAT( P7R_MX(mx,i-1,k-1,p7R_IL), TSC(p7P_IM,k-1))),
		                                      P7_DELTAT( P7R_MX(mx,i-1,k-1,p7R_DL), TSC(p7P_DM,k-1)));
	    P7R_MX(mx,i,k,p7R_IL) = P7R_MX(apd,i,k,p7R_IL) + P7R_MX(apd,i,k,p7R_IG) + 
	                            ESL_MAX( P7_DELTAT( P7R_MX(mx,i-1,k,p7R_ML), TSC(p7P_MI,k)),
				             P7_DELTAT( P7R_MX(mx,i-1,k,p7R_IL), TSC(p7P_II,k)));
	    P7R_MX(mx,i,k,p7R_DL) = ESL_MAX( P7_DELTAT( P7R_MX(mx,i,k-1,p7R_ML), TSC(p7P_MD,k-1)),
					     P7_DELTAT( P7R_MX(mx,i,k-1,p7R_DL), TSC(p7P_DD,k-1)));
	  }

      /* Termination: exits from row ib. */
      if (is_glocal)
	xX = ESL_MAX(P7R_MX(mx,ib,M,p7R_ML), 
		     P7R_MX(mx,ib,M,p7R_DL));
      else           
	for (k = k0; k <= M; k++) {
	  xX = ESL_MAX(xX, P7R_MX(mx,ib,k,p7R_ML));
	  xX = ESL_MAX(xX, P7R_MX(mx,ib,k,p7R_DL));
	}
    } /* end loop over domains d */

  return reference_aec_trace(gm, env, mx, tr);
}
/*---------------- end, AEC/MEG alignment fill ------------------*/




/*****************************************************************
 * 2. Choice functions for AEC/MEG traceback.
 *****************************************************************/

/* Style here is taken from reference_trace.c. The choice functions
 * are abstracted, so the same traceback engine works on any
 * optimization criterion for the AEC path. (Even though the only
 * currently implemented optimization is MEG.)
 */

static inline int
reference_aec_select_m(const P7_PROFILE *gm, const P7_REFMX *mx, int is_glocal, int i, int k)
{
  int          gstate[3] = { p7T_MG, p7T_IG, p7T_DG };  
  int          lstate[3] = { p7T_ML, p7T_IL, p7T_DL };
  float        path[3];

  path[0] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_ML), P7P_TSC(gm, k-1, p7P_MM));
  path[1] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_IL), P7P_TSC(gm, k-1, p7P_IM));
  path[2] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_DL), P7P_TSC(gm, k-1, p7P_DM));
  
  return (is_glocal ? gstate[esl_vec_FArgMax(path, 3)] : lstate[esl_vec_FArgMax(path, 3)]);
}

static inline int
reference_aec_select_i(const P7_PROFILE *gm, const P7_REFMX *mx, int is_glocal, int i, int k)
{
  float        path[2];
  path[0] = P7_DELTAT( P7R_MX(mx, i-1, k, p7R_ML), P7P_TSC(gm, k, p7P_MI));
  path[1] = P7_DELTAT( P7R_MX(mx, i-1, k, p7R_IL), P7P_TSC(gm, k, p7P_II));
  
  if (is_glocal)  return ( path[0] >= path[1] ? p7T_MG : p7T_IG);
  else            return ( path[0] >= path[1] ? p7T_ML : p7T_IL);
}

static inline int
reference_aec_select_d(const P7_PROFILE *gm, const P7_REFMX *mx, int is_glocal, int i, int k)
{
  float        path[2];
  path[0] = P7_DELTAT( P7R_MX(mx, i, k-1, p7R_ML), P7P_TSC(gm, k-1, p7P_MD));
  path[1] = P7_DELTAT( P7R_MX(mx, i, k-1, p7R_DL), P7P_TSC(gm, k-1, p7P_DD));
  
  if (is_glocal)  return ( path[0] >= path[1] ? p7T_MG : p7T_DG);
  else            return ( path[0] >= path[1] ? p7T_ML : p7T_DL);
}

static inline int
reference_aec_select_e(const P7_PROFILE *gm, const P7_REFMX *mx, int is_glocal, int ib, int k0, int *ret_k)
{
  float max  = -eslINFINITY;
  int   kmax = -1;
  int   smax = -1;  
  int   k;

  if (is_glocal) 
    {
      kmax = gm->M;
      smax = ( P7R_MX(mx,ib,gm->M,p7R_ML) >= P7R_MX(mx,ib,gm->M,p7R_DL) ? p7T_MG : p7T_DG);
    }
  else 
    {
      for (k = k0; k <= gm->M; k++)
	{
	  if (P7R_MX(mx,ib,k,p7R_ML) > max) { max = P7R_MX(mx,ib,k,p7R_ML); smax = p7T_ML; kmax = k; }
	  if (P7R_MX(mx,ib,k,p7R_DL) > max) { max = P7R_MX(mx,ib,k,p7R_DL); smax = p7T_DL; kmax = k; }
	}
    }
  *ret_k = kmax;
  return   smax;
}


/*****************************************************************
 * 3. Traceback routine.
 *****************************************************************/

static int
reference_aec_trace(const P7_PROFILE *gm, P7_ENVELOPES *env, const P7_REFMX *mx, P7_TRACE *tr)
{
  int i,k,d;
  int sprv, scur;
  int status;

  /* Contract checks, argument validation */
  ESL_DASSERT1(( mx->M == gm->M && gm->M == env->M ));
  ESL_DASSERT1(( mx->L == env->L ));
  ESL_DASSERT1(( mx->type == p7R_AEC_ALIGN ));

  if   ((status = p7_trace_Append(tr, p7T_T, 0, 0)) != eslOK) return status;
  if   ((status = p7_trace_Append(tr, p7T_C, 0, 0)) != eslOK) return status; 

  for (d = env->D; d >= 1; d--)
    {
      /* Residues ib(d)+1 .. ia(d+1)-1 are assigned to C|J.
       * Sentinel at ia(D+1) = L+1, so no special case needed for ia(d+1)-1 at d=D.
       */
      for (i = env->arr[d+1].ia-1; i > env->arr[d].ib; i--)
	if ((status = p7_trace_Append(tr, (d == env->D ? p7T_C : p7T_J), 0, i)) != eslOK) return status;
      if ((status = p7_trace_Append(tr, p7T_E, 0, i)) != eslOK) return status;
      scur = p7T_E;
      k    = gm->M+1;
      // i now ib(d); scur now E; k is M+1 awaiting traceback from E

      /* Weird but true, we don't need to do DOWN and UP sectors separately here.
       * We're guaranteed that the path will pass thru the anchor i0,k0,M.
       * From there, we know we'll connect to i-1,k-1 supercell.
       * (Sparse, though, can't get away with this, and must do DOWN, UP separately
       * because of its more complicated traversal of matrix supercells.)
       */
      while (i > env->arr[d].ia || (scur != p7T_ML && scur != p7T_MG)) // that is: until we reach M on first row ia(d)...
	{
	  switch (scur) {
	  case p7T_ML: case p7T_MG: sprv = reference_aec_select_m(gm, mx, (env->arr[d].flags & p7E_IS_GLOCAL), i, k);        i--; k--; break;
	  case p7T_IL: case p7T_IG: sprv = reference_aec_select_i(gm, mx, (env->arr[d].flags & p7E_IS_GLOCAL), i, k);        i--;      break;
	  case p7T_DL: case p7T_DG: sprv = reference_aec_select_d(gm, mx, (env->arr[d].flags & p7E_IS_GLOCAL), i, k);             k--; break;
	  case p7T_E:               sprv = reference_aec_select_e(gm, mx, (env->arr[d].flags & p7E_IS_GLOCAL), i, env->arr[d].k0, &k); break;
	  default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in reference AEC traceback");
	  }

	  if (scur == p7T_E) env->arr[d].kb = k;

	  if ( (status = p7_trace_Append(tr, sprv, k, i)) != eslOK) return status;
	  scur = sprv;
	}

      /* Now we're on ia(d) and scur = ML|MG, k=ka(d).
       * Glocal alignments must do a left wing unfolding.
       */
      if (env->arr[d].flags & p7E_IS_GLOCAL)
	{
	  for (; k > 1; k--) 
	    if ( (status = p7_trace_Append(tr, p7T_DG, k-1, i))  != eslOK) return status;
	  sprv = p7T_G;
	}
      else sprv = p7T_L;
      env->arr[d].ka = k;

      if ( (status = p7_trace_Append(tr, sprv,                     0, i)) != eslOK) return status;
      if ( (status = p7_trace_Append(tr, p7T_B,                    0, i)) != eslOK) return status;
      if ( (status = p7_trace_Append(tr, (d == 1 ? p7T_N : p7T_J), 0, i)) != eslOK) return status;
      // i is now ia(d), and we're on B. 
    }
      
  for (i = env->arr[1].ia-1; i >= 1; i--)
    if ((status = p7_trace_Append(tr, p7T_N, 0, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_S, 0, 0))   != eslOK) return status;
  
  tr->M = env->M;
  tr->L = env->L;
  return p7_trace_Reverse(tr);
}



/*****************************************************************
 * 4. Footnotes
 *****************************************************************/

/* 
 * Because we fit the calculation into one matrix, we can't 
 * initialize an ia-1 top row for the up matrices, because 
 * we may be adjacent to an ib row of an immediately preceding
 * domain.
 */


/*****************************************************************
 * 5. Unit tests
 *****************************************************************/
#ifdef p7REFERENCE_AEC_ALIGN_TESTDRIVE
#include "hmmer.h"

/* "crashtestdummy" test
 * Compare randomly selected profile to sequences sampled
 * from that profile. 
 * 
 * The test is not very stringent, because we don't know the
 * true alignment. This is just a test that nothing obviously
 * horrible happens, like a crash, or obviously invalid data
 * structures.
 * 
 * Extended from reference_envelopes.c unit tests.
 *
 * We test:
 *    1. Coordinates of each envelope/alignment are coherent:
 *       1 <= oa <= ia <= i0 <= ib <= ob <= L
 *       1 <= ka <= k0 <= kb <= M
 *       
 *    2. Envelopes do not overlap (assuming default threshold of
 *       0.5 when defining them):
 *         ia(d) > ib(d-1)  for d = 1..D-1
 *       (Outer envelopes, in contrast, can overlap.)
 *       
 *    3. envsc(d) <= asc_sc <= fwdsc.
 *    
 *    4. If D=1 (single domain) in both the generated trace
 *       and the inferred envelopes, and the domain coords in 
 *       the trace are encompassed by the outer envelope,
 *       then envsc(d) >= generated trace score.
 *
 *    5. MEG trace passes validation.
 *
 *    6. Indexed trace domains agree with envelope data.
 */
static void
utest_crashtestdummy(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int N)
{
  char             msg[] = "reference_aec_align:: crash test dummy failed";
  ESL_SQ          *sq    = esl_sq_CreateDigital(abc);
  P7_BG           *bg    = p7_bg_Create(abc);
  P7_HMM          *hmm   = NULL;
  P7_PROFILE      *gm    = p7_profile_Create(M, abc);
  P7_TRACE        *gtr   = p7_trace_Create();            // generated trace
  P7_TRACE        *vtr   = p7_trace_Create();            // Viterbi trace
  P7_TRACE        *tr    = p7_trace_Create();            // MEG trace
  P7_REFMX        *rxf   = p7_refmx_Create(M, 20);       // Fwd, Vit ~~> ASC Decode UP
  P7_REFMX        *rxd   = p7_refmx_Create(M, 20);       // Bck, Decode ~~> ASC Decode DOWN
  P7_REFMX        *afu   = p7_refmx_Create(M, 20);       // ASC Fwd UP
  P7_REFMX        *afd   = p7_refmx_Create(M, 20);       // ASC Fwd DOWN
  P7_REFMX        *apu   = rxf;                          // for 'clarity' we use two names for this mx
  P7_REFMX        *apd   = rxd;                          //   ... and this one too.
  float           *wrk   = NULL;
  P7_ANCHORS      *anch  = p7_anchors_Create();
  P7_ANCHORHASH   *ah    = p7_anchorhash_Create();
  P7_ENVELOPES    *env   = p7_envelopes_Create();
  float            tol   = 0.001;
  char             errbuf[eslERRBUFSIZE];
  float  gsc, fsc, asc;
  int    idx;
  int    d;
  
  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);

  for (idx = 0; idx < N; idx++)
    {
      /* Emit sequence from model, using an arbitrary length model of <M>;
       * restrict the emitted sequence length to 6M, arbitrarily, to 
       * keep it down to something reasonable.
       */
      if ( p7_profile_SetLength(gm, M) != eslOK) esl_fatal(msg);
      do {
	esl_sq_Reuse(sq);
	if (p7_ProfileEmit(rng, hmm, gm, bg, sq, gtr) != eslOK) esl_fatal(msg);
      } while (sq->n > M * 6); 
      if (p7_trace_Index   (gtr)                      != eslOK) esl_fatal(msg);
      if (p7_trace_Score   (gtr, sq->dsq, gm, &gsc)   != eslOK) esl_fatal(msg);

      /* Reset the length model to the actual length sq->n, then
       * put it through the domain postprocessing analysis pipeline
       */
      if ( p7_profile_SetLength(gm, sq->n)                          != eslOK) esl_fatal(msg);
     
      /* First pass analysis */
      if ( p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf, vtr, NULL) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd,      NULL) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd)  != eslOK) esl_fatal(msg);

      /* Anchor determination (MPAS algorithm) */
      if ( p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, vtr, &wrk, ah,
				afu, afd, anch, &asc, NULL, NULL)  != eslOK) esl_fatal(msg);

      /* Reuse rxf,rxd as apu, apd; finish ASC analysis with Backward, Decoding */
      p7_refmx_Reuse(apu);  p7_refmx_Reuse(apd);
      if ( p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch->a, anch->D, apu, apd, NULL)               != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch->a, anch->D, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(msg);

      /* Envelope calculation */
      if ( p7_reference_Envelopes (sq->dsq, sq->n, gm, anch->a, anch->D, apu, apd, afu, afd, env) != eslOK) esl_fatal(msg);

      //p7_envelopes_Dump(stdout, env);

      p7_refmx_Reuse(afu);
      if ( p7_reference_AEC_Align(gm, env, apu, apd, afu, tr) != eslOK) esl_fatal(msg);

      //p7_refmx_Dump(stdout, afu);
      //p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);

      /* Test 1. Coords of each domain are coherent */
      if (anch->D != env->D) esl_fatal(msg);
      for (d = 1; d <= anch->D; d++)
	{
	  if (! (1 <= env->arr[d].oa &&
		 env->arr[d].oa  <= env->arr[d].ia   &&
		 env->arr[d].ia  <= env->arr[d].i0   &&
		 env->arr[d].i0  <= env->arr[d].ib   &&
		 env->arr[d].ib  <= env->arr[d].ob   &&
		 env->arr[d].ob  <= sq->n)) esl_fatal(msg);
	  if (! (1 <= env->arr[d].ka &&
		 env->arr[d].ka <= env->arr[d].k0 &&
		 env->arr[d].k0 <= env->arr[d].kb &&
		 env->arr[d].kb <= gm->M)) esl_fatal(msg);
	}

      /* Test 2. Envelopes do not overlap. */
      for (d = 1; d <= anch->D; d++)
	if (! (env->arr[d].ia > env->arr[d-1].ib)) esl_fatal(msg);  // at d=1, env->arr[d-1=0].ib is a sentinel, 0

      /* Test 3. envsc(d) <= asc_sc <= fwdsc */
      for (d = 1; d <= anch->D; d++)
	if (! (env->arr[d].env_sc <= asc+tol && asc <= fsc+tol)) esl_fatal(msg);

      /* Test 4, only on D=1 case with generated trace's domain 
       * encompassed by the outer envelope 
       */
      if (gtr->ndom == 1 &&  anch->D   == 1 && 
	  gtr->sqfrom[0] >= env->arr[1].oa &&
	  gtr->sqto[0]   <= env->arr[1].ob)
	if (! ( env->arr[1].env_sc >= gsc)) esl_fatal(msg);

      /* Test 5. MEG trace passes validation */
      if (p7_trace_Validate(tr,  abc, sq->dsq, errbuf)  != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
      if (p7_trace_Validate(gtr, abc, sq->dsq, errbuf)  != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
      if (p7_trace_Validate(vtr, abc, sq->dsq, errbuf)  != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);

      /* Test 6. Indexed domains agree with envelope data */
      p7_trace_Index(tr);
      if (tr->ndom != env->D) esl_fatal(msg);
      for (d = 1; d <= env->D; d++)                          // domain numbering in traces is 0..ndom-1, off by one from 1..D in anchors, envelopes
	if (! ( tr->sqfrom[d-1]  == env->arr[d].ia  &&
		tr->sqto[d-1]    == env->arr[d].ib  &&
		tr->hmmfrom[d-1] == env->arr[d].ka  &&
		tr->hmmto[d-1]   == env->arr[d].kb)) esl_fatal(msg);


      p7_envelopes_Reuse(env);
      p7_anchors_Reuse(anch);
      p7_anchorhash_Reuse(ah);
      p7_refmx_Reuse(rxf); p7_refmx_Reuse(rxd);
      p7_refmx_Reuse(afu); p7_refmx_Reuse(afd);
      p7_trace_Reuse(gtr); p7_trace_Reuse(vtr); p7_trace_Reuse(tr);
      esl_sq_Reuse(sq);
    }
      
  if (wrk) free(wrk);
  p7_envelopes_Destroy(env);
  p7_anchors_Destroy(anch);
  p7_anchorhash_Destroy(ah);
  p7_refmx_Destroy(afu); p7_refmx_Destroy(afd);
  p7_refmx_Destroy(rxf); p7_refmx_Destroy(rxd);
  p7_trace_Destroy(vtr); p7_trace_Destroy(gtr); p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_sq_Destroy(sq);
}

/* "singlemulti" test.
 * 
 * Use p7_modelsample_SinglePathedASC() to create a
 * profile/sequence/anchorset comparison that has only a single
 * possible path when anchor set constrained. Now the expected true
 * path and envelope(s) are known, from that single path, so we can compare.
 * 
 * In order to guarantee only one possible path, while allowing
 * multiple domains, the profile is limited to glocal-only.
 *
 * This is an extension of a very similar unit test from
 * reference_envelopes.c.
 * 
 * We test:
 *     1. The MEG trace is identical to the generated path.
 *     2. Trace and envelopes agree on number of domains.
 *     3. For each domain, oa==ia, ob==ib, and these coords
 *        agree with the trace.
 *     4. In the case of a single domain (D=1), the envelope
 *        score == the trace score.
 */
static void
utest_singlemulti(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int N)
{
  char          msg[] = "reference_aec_align:: singlemulti unit test failed";
  P7_BG        *bg    = p7_bg_Create(abc);
  P7_HMM       *hmm   = NULL;
  P7_PROFILE   *gm    = NULL;
  ESL_DSQ      *dsq   = NULL;
  int           L;
  P7_TRACE     *gtr   = NULL;
  P7_TRACE     *tr    = p7_trace_Create();
  P7_ANCHOR    *anch  = NULL;
  int           D;
  P7_REFMX     *afu   = p7_refmx_Create(M, 20);
  P7_REFMX     *afd   = p7_refmx_Create(M, 20);
  P7_REFMX     *apu   = p7_refmx_Create(M, 20);
  P7_REFMX     *apd   = p7_refmx_Create(M, 20);
  P7_ENVELOPES *env   = p7_envelopes_Create();
  float         gsc;
  float         tol   = 0.001;
  int           idx;
  int           d;
  
  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_SinglePathedASC(rng, M, bg, &hmm, &gm, &dsq, &L, &gtr, &anch, &D, &gsc) != eslOK) esl_fatal(msg);

      if ( p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd, NULL)               != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCBackward(dsq, L, gm, anch, D, apu, apd, NULL)               != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCDecoding(dsq, L, gm, anch, D, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(msg);

      if ( p7_reference_Envelopes (dsq, L, gm, anch, D, apu, apd, afu, afd, env)      != eslOK) esl_fatal(msg);

      p7_refmx_Reuse(afu);
      if ( p7_reference_AEC_Align( gm, env, apu, apd, afu, tr) != eslOK) esl_fatal(msg);
      if ( p7_trace_Index(tr)                                  != eslOK) esl_fatal(msg);

      //p7_refmx_Dump(stdout, afu);
      //p7_trace_DumpAnnotated(stdout, gtr, gm, dsq);
      //p7_trace_DumpAnnotated(stdout,  tr, gm, dsq);

      /* Test 1. The traces are identical. */
      if ( p7_trace_Compare(gtr, tr, 0.0) != eslOK) esl_fatal(msg);

      /* Test 2. Domain #'s agree */
      if (! (gtr->ndom == D && env->D == D)) esl_fatal(msg);
      if (! (tr->ndom  == D))                esl_fatal(msg);

      /* Test 3. Envelope coords (and outer env coords) match trace */
      for (d = 1; d <= D; d++)
	{
	  if (! (env->arr[d].ia == gtr->sqfrom[d-1] &&   // beware: domains in trace are 0..ndom-1, off by one from env's 1..D
		 env->arr[d].ia ==  tr->sqfrom[d-1] &&
		 env->arr[d].ia == env->arr[d].oa)) esl_fatal(msg);

	  if (! (env->arr[d].ib == gtr->sqto[d-1]   &&
		 env->arr[d].ib ==  tr->sqto[d-1]   &&
		 env->arr[d].ib == env->arr[d].ob)) esl_fatal(msg);

	  if (! (env->arr[d].ka == gtr->hmmfrom[d-1] &&
		 env->arr[d].ka ==  tr->hmmfrom[d-1])) esl_fatal(msg);

	  if (! (env->arr[d].kb == gtr->hmmto[d-1] &&
		 env->arr[d].kb ==  tr->hmmto[d-1])) esl_fatal(msg);
	}

      /* Test 4. If D == 1, envelope score == trace score. */
      if (D == 1 &&  esl_FCompare_old(env->arr[1].env_sc, gsc, tol) != eslOK) esl_fatal(msg);

      p7_envelopes_Reuse(env);
      p7_refmx_Reuse(afu); p7_refmx_Reuse(afd);
      p7_refmx_Reuse(apu); p7_refmx_Reuse(apd);
      p7_trace_Reuse(tr);

      free(dsq);
      free(anch);
      p7_trace_Destroy(gtr);
      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
    }
     
  p7_trace_Destroy(tr);
  p7_envelopes_Destroy(env);
  p7_refmx_Destroy(afu); p7_refmx_Destroy(afd);
  p7_refmx_Destroy(apu); p7_refmx_Destroy(apd);
  p7_bg_Destroy(bg);
}


#endif /*p7REFERENCE_AEC_ALIGN_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/


/*****************************************************************
 * 6. Test driver
 *****************************************************************/
#ifdef p7REFERENCE_AEC_ALIGN_TESTDRIVE

#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for AEC/MEG alignment inference";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  int             M    = 20;

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_crashtestdummy(rng, M, abc, 10);
  utest_singlemulti   (rng, M, abc, 10);

  fprintf(stderr, "#  status = ok\n");

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}

#endif /*p7REFERENCE_AEC_ALIGN_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/



/*****************************************************************
 * 7. Example
 *****************************************************************/
#ifdef p7REFERENCE_AEC_ALIGN_EXAMPLE
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of ASC MEG alignment step";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_ANCHORS     *anch    = p7_anchors_Create();
  P7_ANCHORHASH  *ah      = p7_anchorhash_Create();
  P7_REFMX       *rxf     = NULL;
  P7_REFMX       *rxd     = NULL;
  P7_REFMX       *afu     = NULL;
  P7_REFMX       *afd     = NULL;
  P7_REFMX       *apu     = NULL;
  P7_REFMX       *apd     = NULL;
  P7_TRACE       *tr      = NULL;
  float          *wrk     = NULL;
  P7_ENVELOPES   *env     = p7_envelopes_Create();
  float           fsc, vsc, asc, asc_b;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Read one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);
  esl_sqfile_Close(sqfp);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);

  /* Set the profile and null model's target length models */
  p7_bg_SetLength     (bg, sq->n);
  p7_profile_SetLength(gm, sq->n);

  /* Allocate DP matrices and tracebacks */
  rxf = p7_refmx_Create(gm->M, sq->n);
  rxd = p7_refmx_Create(gm->M, sq->n);
  tr  = p7_trace_Create();
  afu = p7_refmx_Create(gm->M, sq->n);
  afd = p7_refmx_Create(gm->M, sq->n);

  /* First pass analysis */
  p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf,  tr, &vsc);
  p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc);
  p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd, NULL);   
  p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd);   

  /* MPAS algorithm gets us an anchor set */
  p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, tr, &wrk, ah,
		       afu, afd, anch, &asc, NULL, NULL);

  
  /* We no longer need rxf and rxd. 
   * Use their space for apu/apd pair, which will briefly
   * hold ASC Backward matrices, then get used for ASC Decoding.
   */
  apu = rxf; p7_refmx_Reuse(apu);
  apd = rxd; p7_refmx_Reuse(apd);

  p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch->a, anch->D, apu, apd, &asc_b);
  
  /* ASC Decoding takes afu/afd and abu/abd as input;
   * overwrites abu/abd with decoding matrices
   */
  p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch->a, anch->D, afu, afd, apu, apd, apu, apd);

  /* Envelope calculation needs to get four matrices:
   * ASC Decoding pair, apu/apd, and it will leave these constant;
   * ASC Forward pair,  afu/afd, and it will overwrite these.
   */
  p7_reference_Envelopes(sq->dsq, sq->n, gm, anch->a, anch->D, apu, apd, afu, afd, env);

  /* MEG alignment step uses afu as its matrix; apu/apd decoding matrices as its input */
  p7_refmx_Reuse(afu);
  p7_trace_Reuse(tr);
  p7_reference_AEC_Align(gm, env, apu, apd, afu, tr);

  //p7_refmx_Dump(stdout, afu);
  p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);
  p7_envelopes_Dump(stdout, env);

  p7_envelopes_Destroy(env);
  p7_anchors_Destroy(anch);
  p7_anchorhash_Destroy(ah);
  if (wrk) free(wrk);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(rxd);
  p7_refmx_Destroy(rxf);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_AEC_ALIGN_EXAMPLE*/

