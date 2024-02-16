/* Reference implementation of anchorset/envelope-constrained (AEC) alignment.
 *
 * This is the reference implementation, used for testing. It is not used in HMMER's
 * main programs. The production code uses sparse AEC DP.
 *
 * Contents:
 *   1. AEC alignment, DP fill
 *   2. Choice functions for AEC traceback.
 *   3. AEC traceback
 *   4. Experiment: coverage stats, ensemble vs Viterbi inference [was decoding-hunt]
 *   5. Experiment 2: details on one profile/sequence comparison  [was decoding-study]
 *   6. Unit tests
 *   7. Test driver
 *   8. Example
 */
#include <h4_config.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "h4_envset.h"
#include "h4_path.h"
#include "h4_profile.h"
#include "h4_refmx.h"

static int reference_aec_trace(const H4_PROFILE *hmm, const H4_REFMX *mx, H4_ENVSET *env, H4_PATH *pi);

/*****************************************************************
 * 1. AEC alignment, fill
 *****************************************************************/

/* Function:  h4_reference_aec_Align()
 * Synopsis:  Anchor/envelope-constrained (AEC) alignment
 * Incept:    SRE, Fri 09 Jul 2021
 *
 * Purpose:   Given profile <hmm>, envelopes <env>, and ASC posterior
 *            decoding matrices <apu>/<apd>, calculate an AEC path
 *            for the sequence (including all domains) in a DP matrix <mx>
 *            and return the inferred path in <pi>, and set the <ka>/<kb>
 *            profile bounds in <env>.
 *
 *            Caller provides an allocated <mx> for a DP matrix space
 *            and an allocated <pi> for the path; they are reused and
 *            reallocated as needed.
 *
 *            Sequence doesn't need to be provided, because everything
 *            we need to know about it is in <env> and <apu>/<apd>.
 *            
 * Args:      input:
 *            hmm - profile
 *            apu - ASC posterior decoding UP matrix
 *            apd -  ... and DOWN matrix
 *
 *            input/output (object data modified):
 *            env - envelope data; ka/kb endpoints assigned here
 *            
 *            output (allocated object provided):
 *            mx  - DP matrix
 *            pi  - optimal AEC path for entire sequence
 *            
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on reallocation failure
 */
int
h4_reference_aec_Align(const H4_PROFILE *hmm, const H4_REFMX *apu, const H4_REFMX *apd, H4_ENVSET *env, H4_REFMX *mx, H4_PATH *pi)
{
  int   M  = env->M;
  int   L  = env->L;
  float xX = 0.;
  int   ia, i0, k0, ib, is_glocal;
  int   d, i, k;
  int   status;

  /* input validation */
  ESL_DASSERT1(( hmm->M == M && apu->M == M && apd->M == M));
  ESL_DASSERT1(( apu->L == L && apd->L == L ));
  ESL_DASSERT1(( apu->type == h4R_ASC_DECODE_UP   ));
  ESL_DASSERT1(( apd->type == h4R_ASC_DECODE_DOWN ));

  /* reallocation/reuse, as needed */
  if (( status = h4_refmx_GrowTo(mx, apu->M, apu->L))                 != eslOK) return status;
  if (( status = h4_refmx_SetType(mx, apu->M, apu->L, h4R_AEC_ALIGN)) != eslOK) return status;
  if (( status = h4_path_Reuse(pi))                                   != eslOK) return status;

  for (d = 1; d <= env->D; d++)
    {
      /* simplify notation below with some tmp vars for this envelope [d] */
      ia = env->e[d].ia;
      i0 = env->e[d].i0;
      k0 = env->e[d].k0;
      ib = env->e[d].ib;
      is_glocal = (env->e[d].flags & h4E_IS_GLOCAL);

      /* Initialization of row ia(d), if UP sector exists */
      if (ia < i0)
	{
          H4R_MX(mx,ia,0,h4R_ML) = H4R_MX(mx,ia,0,h4R_IL) = H4R_MX(mx,ia,0,h4R_DL) = -eslINFINITY;
	  for (k = 1; k < k0; k++)
	    {
	      H4R_MX(mx,ia,k,h4R_ML) = H4R_MX(apu,ia,k,h4R_ML) + H4R_MX(apu,ia,k,h4R_MG) + 
                                       (is_glocal ? H4_DELTAT(xX, hmm->tsc[k-1][h4_GM]) :
                                                    H4_DELTAT(xX, hmm->tsc[k-1][h4_LM]));
	      H4R_MX(mx,ia,k,h4R_IL) = (is_glocal ? H4R_MX(apu,ia,k,h4R_IL) + H4R_MX(apu,ia,k,h4R_IG) + H4_DELTAT(xX, hmm->tsc[k-1][h4_GI]):
                                                    -eslINFINITY);
	      H4R_MX(mx,ia,k,h4R_DL) = ESL_MAX( ESL_MAX(H4_DELTAT(H4R_MX(mx,ia,k-1,h4R_ML), hmm->tsc[k-1][h4_MD]),
						        H4_DELTAT(H4R_MX(mx,ia,k-1,h4R_IL), hmm->tsc[k-1][h4_ID])),
                                                        H4_DELTAT(H4R_MX(mx,ia,k-1,h4R_DL), hmm->tsc[k-1][h4_DD]));
            }
        }

      /* Remaining rows of UP sector (if any) */
      for (i = ia+1; i < i0; i++)
        {
          H4R_MX(mx,i,0,h4R_ML) = H4R_MX(mx,i,0,h4R_IL) = H4R_MX(mx,i,0,h4R_DL) = -eslINFINITY;
          for (k = 1; k < k0; k++)
            {
              H4R_MX(mx,i,k,h4R_ML) = H4R_MX(apu,i,k,h4R_ML) + H4R_MX(apu,i,k,h4R_MG) + 
                                      ESL_MAX( ESL_MAX( H4_DELTAT( H4R_MX(mx,i-1,k-1,h4R_ML), hmm->tsc[k-1][h4_MM]),
                                                        H4_DELTAT( H4R_MX(mx,i-1,k-1,h4R_IL), hmm->tsc[k-1][h4_IM])),
                                                        H4_DELTAT( H4R_MX(mx,i-1,k-1,h4R_DL), hmm->tsc[k-1][h4_DM]));
              H4R_MX(mx,i,k,h4R_IL) = H4R_MX(apu,i,k,h4R_IL) + H4R_MX(apu,i,k,h4R_IG) + 
                                      ESL_MAX( ESL_MAX( H4_DELTAT( H4R_MX(mx,i-1,k,h4R_ML), hmm->tsc[k][h4_MI]),
                                                        H4_DELTAT( H4R_MX(mx,i-1,k,h4R_IL), hmm->tsc[k][h4_II])),
                                                        H4_DELTAT( H4R_MX(mx,i-1,k,h4R_DL), hmm->tsc[k][h4_DI]));
              H4R_MX(mx,i,k,h4R_DL) = ESL_MAX( ESL_MAX( H4_DELTAT( H4R_MX(mx,i,k-1,h4R_ML), hmm->tsc[k-1][h4_MD]),
                                                        H4_DELTAT( H4R_MX(mx,i,k-1,h4R_IL), hmm->tsc[k-1][h4_ID])),
                                                        H4_DELTAT( H4R_MX(mx,i,k-1,h4R_DL), hmm->tsc[k-1][h4_DD]));
            }
        }

      /* Anchor cell i0(d),k0(d); usually initialized from last
       * supercell of UP sector, but rarely UP sector may not
       * exist. 
       */
      if (ia < i0) H4R_MX(mx,i0,k0,h4R_ML) = H4R_MX(apd,i0,k0,h4R_ML) + H4R_MX(apd,i0,k0,h4R_MG) + 
                                             ESL_MAX( ESL_MAX( H4_DELTAT( H4R_MX(mx,i0-1,k0-1,h4R_ML), hmm->tsc[k0-1][h4_MM]),
				                               H4_DELTAT( H4R_MX(mx,i0-1,k0-1,h4R_IL), hmm->tsc[k0-1][h4_IM])),
                                                               H4_DELTAT( H4R_MX(mx,i0-1,k0-1,h4R_DL), hmm->tsc[k0-1][h4_DM]));
      else         H4R_MX(mx,i0,k0,h4R_ML) = H4R_MX(apd,i0,k0,h4R_ML) + H4R_MX(apd,i0,k0,h4R_MG) +
                                             (is_glocal ? H4_DELTAT(xX, hmm->tsc[k0-1][h4_GM]) : 
                                                          H4_DELTAT(xX, hmm->tsc[k0-1][h4_LM]));
      H4R_MX(mx,i0,k0,h4R_IL) = -eslINFINITY;
      H4R_MX(mx,i0,k0,h4R_DL) = -eslINFINITY;

      /* Remainder of the top DOWN row, i0 for k>k0 */
      for (k = k0+1; k <= M; k++)
	{
	  H4R_MX(mx,i0,k,h4R_ML) = -eslINFINITY;
	  H4R_MX(mx,i0,k,h4R_IL) = -eslINFINITY;
	  H4R_MX(mx,i0,k,h4R_DL) = ESL_MAX( ESL_MAX( H4_DELTAT( H4R_MX(mx,i0,k-1,h4R_ML), hmm->tsc[k-1][h4_MD]),
				                     H4_DELTAT( H4R_MX(mx,i0,k-1,h4R_IL), hmm->tsc[k-1][h4_ID])),
                                                     H4_DELTAT( H4R_MX(mx,i0,k-1,h4R_DL), hmm->tsc[k-1][h4_DD]));
        }

      /* DOWN sector recursion */
      for (i = i0+1; i <= ib; i++)
        {
          H4R_MX(mx,i,k0-1,h4R_ML) = H4R_MX(mx,i,k0-1,h4R_IL) = H4R_MX(mx,i,k0-1,h4R_DL) = -eslINFINITY;  
          for (k = k0; k <= M; k++)
            {
              H4R_MX(mx,i,k,h4R_ML) = H4R_MX(apd,i,k,h4R_ML) + H4R_MX(apd,i,k,h4R_MG) + 
                                      ESL_MAX( ESL_MAX( H4_DELTAT( H4R_MX(mx,i-1,k-1,h4R_ML), hmm->tsc[k-1][h4_MM]),
                                                        H4_DELTAT( H4R_MX(mx,i-1,k-1,h4R_IL), hmm->tsc[k-1][h4_IM])),
                                                        H4_DELTAT( H4R_MX(mx,i-1,k-1,h4R_DL), hmm->tsc[k-1][h4_DM]));
              H4R_MX(mx,i,k,h4R_IL) = H4R_MX(apd,i,k,h4R_IL) + H4R_MX(apd,i,k,h4R_IG) + 
                                      ESL_MAX( ESL_MAX( H4_DELTAT( H4R_MX(mx,i-1,k,h4R_ML), hmm->tsc[k][h4_MI]),
                                                        H4_DELTAT( H4R_MX(mx,i-1,k,h4R_IL), hmm->tsc[k][h4_II])),
                                                        H4_DELTAT( H4R_MX(mx,i-1,k,h4R_DL), hmm->tsc[k][h4_DI]));
              H4R_MX(mx,i,k,h4R_DL) = ESL_MAX( ESL_MAX( H4_DELTAT( H4R_MX(mx,i,k-1,h4R_ML), hmm->tsc[k-1][h4_MD]),
                                                        H4_DELTAT( H4R_MX(mx,i,k-1,h4R_IL), hmm->tsc[k-1][h4_ID])),
                                                        H4_DELTAT( H4R_MX(mx,i,k-1,h4R_DL), hmm->tsc[k-1][h4_DD]));
            }
        }

      /* Termination: exits from row ib. */
      if (is_glocal)
	xX = ESL_MAX(H4R_MX(mx,ib,M,h4R_ML), 
		     H4R_MX(mx,ib,M,h4R_DL));
      else           
	for (k = k0; k <= M; k++) {
	  xX = ESL_MAX(xX, H4R_MX(mx,ib,k,h4R_ML));
	  xX = ESL_MAX(xX, H4R_MX(mx,ib,k,h4R_DL));
	}
    }

  return reference_aec_trace(hmm, mx, env, pi);  // traceback also records ka/kb endpoints
}


/*****************************************************************
 * 2. Choice functions for AEC traceback.
 *****************************************************************/

/* Style here is taken from reference_dp.c. The choice functions
 * are abstracted, so the same traceback engine works on any
 * optimization criterion for the AEC path. (Even though the only
 * currently implemented optimization is MEG.)
 */
static inline int
reference_aec_select_m(const H4_PROFILE *hmm, const H4_REFMX *mx, int is_glocal, int i, int k)
{
  float        path[3];

  path[0] = H4_DELTAT( H4R_MX(mx, i-1, k-1, h4R_ML), hmm->tsc[k-1][h4_MM]);
  path[1] = H4_DELTAT( H4R_MX(mx, i-1, k-1, h4R_IL), hmm->tsc[k-1][h4_IM]);
  path[2] = H4_DELTAT( H4R_MX(mx, i-1, k-1, h4R_DL), hmm->tsc[k-1][h4_DM]);
  
  return ((is_glocal ? h4P_MG : h4P_ML) + esl_vec_FArgMax(path, 3)); // chummy assumption of M-I-D order in h4P
}

static inline int
reference_aec_select_i(const H4_PROFILE *hmm, const H4_REFMX *mx, int is_glocal, int i, int k)
{
  float        path[3];

  path[0] = H4_DELTAT( H4R_MX(mx, i-1, k, h4R_ML), hmm->tsc[k][h4_MI]);
  path[1] = H4_DELTAT( H4R_MX(mx, i-1, k, h4R_IL), hmm->tsc[k][h4_II]);
  path[2] = H4_DELTAT( H4R_MX(mx, i-1, k, h4R_DL), hmm->tsc[k][h4_DI]);
  
  return ((is_glocal ? h4P_MG : h4P_ML) + esl_vec_FArgMax(path, 3)); 
}

static inline int
reference_aec_select_d(const H4_PROFILE *hmm, const H4_REFMX *mx, int is_glocal, int i, int k)
{
  float        path[3];

  path[0] = H4_DELTAT( H4R_MX(mx, i, k-1, h4R_ML), hmm->tsc[k-1][h4_MD]);
  path[1] = H4_DELTAT( H4R_MX(mx, i, k-1, h4R_IL), hmm->tsc[k-1][h4_ID]);
  path[2] = H4_DELTAT( H4R_MX(mx, i, k-1, h4R_DL), hmm->tsc[k-1][h4_DD]);

  return ((is_glocal ? h4P_MG : h4P_ML) + esl_vec_FArgMax(path, 3)); 
}

static inline int
reference_aec_select_e(const H4_PROFILE *hmm, const H4_REFMX *mx, int is_glocal, int ib, int k0, int *ret_k)
{
  float max  = -eslINFINITY;
  int   kmax = -1;
  int   smax = -1;  
  int   k;

  if (is_glocal) 
    {
      kmax = hmm->M;
      smax = ( H4R_MX(mx,ib,hmm->M,h4R_ML) >= H4R_MX(mx,ib,hmm->M,h4R_DL) ? h4P_MG : h4P_DG);
    }
  else 
    {
      for (k = k0; k <= hmm->M; k++)
	{
	  if (H4R_MX(mx,ib,k,h4R_ML) > max) { max = H4R_MX(mx,ib,k,h4R_ML); smax = h4P_ML; kmax = k; }
	  if (H4R_MX(mx,ib,k,h4R_DL) > max) { max = H4R_MX(mx,ib,k,h4R_DL); smax = h4P_DL; kmax = k; }
	}
    }
  *ret_k = kmax;
  return   smax;
}

/*****************************************************************
 * 3. AEC traceback 
 *****************************************************************/

static int
reference_aec_trace(const H4_PROFILE *hmm, const H4_REFMX *mx, H4_ENVSET *env, H4_PATH *pi)
{
  int i,k,d,s;
  int is_glocal;
  int status;

  /* argument validation */
  ESL_DASSERT1(( mx->M == hmm->M && hmm->M == env->M ));
  ESL_DASSERT1(( mx->L == env->L ));
  ESL_DASSERT1(( mx->type == h4R_AEC_ALIGN ));

  for (d = env->D; d >= 1; d--)
    {
      /* Residues ib(d)+1 .. ia(d+1)-1 are assigned to C|J.
       * That's ia(d+1)-ib(d)-1 residues; +1 state for emit-on-transition.
       * Sentinel at ia(D+1) = L+1, so no special case needed for ia(d+1)-1 at d=D.
       */
      if ((status = h4_path_AppendElement(pi, (d == env->D ? h4P_C : h4P_J), env->e[d+1].ia - env->e[d].ib)) != eslOK) return status;
      s = h4P_E;
      k = hmm->M+1;
      i = env->e[d].ib;
      is_glocal = env->e[d].flags & h4E_IS_GLOCAL;

      /* Weird but true, we don't need to do DOWN and UP sectors separately here.
       * We're guaranteed that the path will pass thru the anchor i0,k0,M.
       * From the anchor supercell, we know we'll connect to i-1,k-1 supercell.
       * From anywhere else on anchor row i0, can only be in D, looking back at k-1.
       * Same on row ia; can only be in D, looking at k-1, until we reach the 
       * alignment start, which can be MG|ML|IG.
       */
      do
	{
	  switch (s) {
	  case h4P_ML: case h4P_MG: s = reference_aec_select_m(hmm, mx, is_glocal, i, k);       i--; k--; break;
	  case h4P_IL: case h4P_IG: s = reference_aec_select_i(hmm, mx, is_glocal, i, k);       i--;      break;
	  case h4P_DL: case h4P_DG: s = reference_aec_select_d(hmm, mx, is_glocal, i, k);            k--; break;
	  case h4P_E:               s = reference_aec_select_e(hmm, mx, is_glocal, i, env->e[d].k0, &k);  env->e[d].kb = k; break; // record kb endpoint
	  default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in reference AEC traceback");
	  }
	  if ( (status = h4_path_Append(pi, s)) != eslOK) return status;
	} while (i > env->e[d].ia || h4_path_IsD(s)); // i.e. stop on MG|ML|IG on row ia

      /* Glocal alignments must do left wing unfolding of k or k-1 delete states, depending whether we're in IG or MG */
      if (is_glocal && k > 1)
        {
          if ((status = h4_path_AppendElement(pi, h4P_DG, (s == h4P_MG ? k-1 : k))) != eslOK) return status;
          k = 1;
        }

      h4_path_AppendElement(pi, (is_glocal? h4P_G : h4P_L), k);
      env->e[d].ka = k; // record ka endpoint
      // i is now ia(d)
    }
      
  /* Remaining i-1 residues are in N */
  if ((status = h4_path_AppendElement(pi, h4P_N, i)) != eslOK) return status;
  return h4_path_Reverse(pi);
}


/***************************************************************** 
 * 4. Experiment: coverage stats, ensemble vs Viterbi inference
 *****************************************************************/
#ifdef h4REFERENCE_AEC_EXPERIMENT

#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "general.h"
#include "h4_anchorhash.h"
#include "h4_anchorset.h"
#include "h4_envset.h"
#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"
#include "h4_refmx.h"

#include "reference_dp.h"
#include "reference_asc.h"
#include "reference_mpas.h"
#include "reference_envelopes.h"
#include "reference_aec.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                   docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show brief help summary",                   0 },
  { "--local",    eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "use local mode, not dual glocal/local",     0 },
  { "--seed",     eslARG_INT,     "0", NULL, NULL,  NULL,   NULL,  NULL, "set random number generator seed",          0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show HMMER version number",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 2, argc, argv, "coverage statistics, probabilistic inference vs. Viterbi", "[-options] <hmmfile> <seqfile>");
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  int             infmt   = eslSQFILE_UNKNOWN;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  ESL_ALPHABET   *abc     = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = NULL;
  H4_PATH        *pi      = h4_path_Create();
  H4_REFMX       *rxf     = NULL;
  H4_REFMX       *rxd     = NULL;
  H4_REFMX       *afu     = NULL;
  H4_REFMX       *afd     = NULL;
  H4_REFMX       *apu     = NULL;
  H4_REFMX       *apd     = NULL;
  H4_REFMX       *aec     = NULL;
  float          *wrk     = NULL;
  H4_ANCHORSET   *anch    = h4_anchorset_Create(0,0,0);
  H4_ANCHORHASH  *ah      = h4_anchorhash_Create();
  H4_ENVSET      *env     = h4_envset_Create(0,0,0);
  float          vsc, fsc, bsc, asc;
  int            vdom, vfull, vcov;
  int            pdom, pfull, pcov;
  int            dcov50, dcov90;
  int            status;

  if (esl_opt_GetBoolean(go, "--local")) h4_mode_SetLocal(mo);

  status = h4_hmmfile_Open(hmmfile, NULL, &hfp);
  if (status != eslOK) esl_fatal("Error: failed to open %s for reading profile HMM(s)\n%s\n", strcmp(hmmfile, "-") == 0? "<stdin>" : hmmfile, hfp->errmsg);

  status = h4_hmmfile_Read(hfp, &abc, &hmm);
  if      (status == eslEFORMAT)   esl_fatal("Parse failed, bad profile HMM file format in %s:\n   %s",  strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status == eslEOF)       esl_fatal("Empty input? No profile HMM found in %s\n",                strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status != eslOK)        esl_fatal("Unexpected error reading profile HMM from %s (code %d)\n", strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, status);

  status = esl_sqfile_OpenDigital(abc, seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open %s for reading sequences\n",             strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile);
  else if (status == eslEFORMAT)   esl_fatal("Format of sequence file %s unrecognized\n",             strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile);
  else if (status != eslOK)        esl_fatal("Unexpected error opening sequence file %s (code %d)\n", strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile, status);

  sq = esl_sq_CreateDigital(abc);

  rxf = h4_refmx_Create(hmm->M, 400);
  rxd = h4_refmx_Create(hmm->M, 400);
  afu = h4_refmx_Create(hmm->M, 400);
  afd = h4_refmx_Create(hmm->M, 400);
  apu = rxf;  // reusing space to conserve memory
  apd = rxd;
  aec = afu;

  while ((status = esl_sqio_Read(sqfp, sq)) != eslEOF)
    {
      if      (status == eslEFORMAT) esl_fatal("Failed to parse sequence format of %s\n%s\n",           strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile, esl_sqfile_GetErrorBuf(sqfp));
      else if (status != eslOK)      esl_fatal("Unexpected error reading sequence from %s (code %d)\n", strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile, status);

      h4_mode_SetLength(mo, sq->n);
      h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxf, pi, &vsc);
      h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,     &fsc);
      h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxd,     &bsc);   
      h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxd, rxd);   

      vdom   = h4_path_GetDomainCount(pi);
      vfull  = h4_path_GetFullDomainCount(pi, hmm->M);
      vcov   = h4_path_GetHomologyCoverage(pi);
      dcov50 = h4_refmx_GetDecodedDomainCoverage(rxd, 0.50);
      dcov90 = h4_refmx_GetDecodedDomainCoverage(rxd, 0.90);

      h4_reference_MPAS(rng, sq->dsq, sq->n, hmm, mo, rxf, rxd, pi, &wrk, ah,
                        afu, afd, anch, &asc, NULL, NULL);

      h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch,           apu, apd, NULL);
      h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, apu, apd, apu, apd);
      h4_reference_Envelopes   (sq->dsq, sq->n, hmm, mo, anch, apu, apd, afu, afd, env);
      h4_reference_aec_Align(hmm, apu, apd, env, aec, pi);

      pdom  = h4_path_GetDomainCount(pi);
      pfull = h4_path_GetFullDomainCount(pi, hmm->M);
      pcov  = h4_path_GetHomologyCoverage(pi);

      printf("%-30s %9.2f %9.2f %9.2f %4d %4d %4d %4d %5d %5d %5d %5d\n",
             sq->name,
             vsc, fsc, asc,
             vdom,  pdom,
             vfull, pfull,
             vcov,  pcov, dcov50, dcov90);
      
      esl_sq_Reuse(sq);
    }

  esl_sqfile_Close(sqfp);
  h4_hmmfile_Close(hfp);
  
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  esl_sq_Destroy(sq);
  h4_path_Destroy(pi);
  h4_refmx_Destroy(rxf);   h4_refmx_Destroy(rxd);
  h4_refmx_Destroy(afu);   h4_refmx_Destroy(afd);
  h4_anchorset_Destroy(anch);
  h4_anchorhash_Destroy(ah);
  h4_envset_Destroy(env);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
}

#endif // h4REFERENCE_AEC_EXPERIMENT

/*****************************************************************
 * 5. Experiment 2: details on one profile/seq comparison 
 *****************************************************************/
#ifdef h4REFERENCE_AEC_EXPERIMENT2

#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "general.h"
#include "h4_anchorhash.h"
#include "h4_anchorset.h"
#include "h4_envset.h"
#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"
#include "h4_refmx.h"

#include "reference_dp.h"
#include "reference_asc.h"
#include "reference_mpas.h"
#include "reference_envelopes.h"
#include "reference_aec.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                   docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show brief help summary",                            1 },
  { "--local",    eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "use local mode, not dual glocal/local",              1 },
  { "--epath",    eslARG_OUTFILE,NULL, NULL, NULL,  NULL,   NULL,  NULL, "save ensemble verbose path dump to file <f>",        1 },
  { "--vpath",    eslARG_OUTFILE,NULL, NULL, NULL,  NULL,   NULL,  NULL, "save Viterbi verbose path dump to file <f>",         1 },

  { "--ia",       eslARG_INT,    NULL, NULL, NULL,  NULL,   NULL,  NULL, "when plotting, start sequence at i=ia not 1",        2 },
  { "--ib",       eslARG_INT,    NULL, NULL, NULL,  NULL,   NULL,  NULL, "when plotting, end sequence at at i=ib not L",       2 },
  { "--ka",       eslARG_INT,    NULL, NULL, NULL,  NULL,   NULL,  NULL, "when plotting, start profile at k=ka not 1",         2 },
  { "--kb",       eslARG_INT,    NULL, NULL, NULL,  NULL,   NULL,  NULL, "when plotting, end profile at at k=kb not M",        2 },
  { "--diplot",   eslARG_OUTFILE,NULL, NULL, NULL,  NULL,   NULL,  NULL, "save domain inference xy plots to file <f>",         2 },
  { "--heatmap",  eslARG_OUTFILE,NULL, NULL, NULL,  NULL,   NULL,  NULL, "save PostScript heatmap of decoding matrix to <f>",  2 },
  { "--ppmx",     eslARG_OUTFILE,NULL, NULL, NULL,  NULL,   NULL,  NULL, "save posterior decoding matrix to <f>",              2 },

  { "--seed",     eslARG_INT,     "0", NULL, NULL,  NULL,   NULL,  NULL, "set random number generator seed",                   3 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show HMMER version number",                          3 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static ESL_GETOPTS *
process_cmdline(const ESL_OPTIONS *options, int argc, char **argv)
{
  ESL_GETOPTS *go            = esl_getopts_Create(options);
  char        *lastslash     = strrchr(argv[0], '/');
  char        *cmdname       = (lastslash? lastslash+1 : argv[0]);
  char         usage[]       = "[-options] <hmmfile> <seqfile>";
  char         description[] = "coverage statistics, ensemble inference vs. Viterbi";
  
  h4_Init();

  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK) 
    {
      esl_printf("Failed to parse command line: %s\n", go->errbuf);
      esl_printf("Usage:\n  %s %s\n", cmdname, usage);
      esl_printf("\nTo see more help on available options, do `%s -h`\n\n", cmdname);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      esl_printf("%s : %s\n", cmdname, description);
      esl_printf("\nUsage:\n  %s %s\n", cmdname, usage);
      
      esl_printf("\nOptions:\n");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      esl_printf("\noptions for dumping plots and data; possibly windowed ia..ib/ka..kb\n");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      esl_printf("\nother options:\n");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 2)
    {
      esl_printf("Incorrect number of command line arguments.\n");
      esl_printf("Usage:\n  %s %s\n", cmdname, usage);
      esl_printf("\nTo see more help on available options, do `%s -h`\n\n", cmdname);
      exit(1);
    }
  return go;
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = process_cmdline(options, argc, argv);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  int             infmt   = eslSQFILE_UNKNOWN;
  char           *outfile = NULL;
  FILE           *ofp     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  ESL_ALPHABET   *abc     = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = NULL;
  H4_PATH        *pi      = h4_path_Create();
  H4_REFMX       *rxf     = NULL;
  H4_REFMX       *rxd     = NULL;
  H4_REFMX       *afu     = NULL;
  H4_REFMX       *afd     = NULL;
  H4_REFMX       *apu     = NULL;
  H4_REFMX       *apd     = NULL;
  H4_REFMX       *aec     = NULL;
  float          *wrk     = NULL;
  H4_ANCHORSET   *anch    = h4_anchorset_Create(0,0,0);
  H4_ANCHORHASH  *ah      = h4_anchorhash_Create();
  H4_ENVSET      *env     = h4_envset_Create(0,0,0);
  H4_MPAS_PARAMS *prm     = h4_mpas_params_Create();
  H4_MPAS_STATS  *stats   = h4_mpas_stats_Create();
  int             ia      = (esl_opt_IsOn(go, "--ia") ? esl_opt_GetInteger(go, "--ia") : 0);
  int             ib      = (esl_opt_IsOn(go, "--ib") ? esl_opt_GetInteger(go, "--ib") : 0);
  int             ka      = (esl_opt_IsOn(go, "--ka") ? esl_opt_GetInteger(go, "--ka") : 0);
  int             kb      = (esl_opt_IsOn(go, "--kb") ? esl_opt_GetInteger(go, "--kb") : 0);
  float           vsc, fsc, bsc, asc;
  int             d,D;
  int             domia, domib, domka, domkb;
  int             status;
  
  if (esl_opt_GetBoolean(go, "--local")) h4_mode_SetLocal(mo);
  prm->be_verbose     = TRUE;

  status = h4_hmmfile_Open(hmmfile, NULL, &hfp);
  if (status != eslOK) esl_fatal("Error: failed to open %s for reading profile HMM(s)\n%s\n", strcmp(hmmfile, "-") == 0? "<stdin>" : hmmfile, hfp->errmsg);

  status = h4_hmmfile_Read(hfp, &abc, &hmm);
  if      (status == eslEFORMAT)   esl_fatal("Parse failed, bad profile HMM file format in %s:\n   %s",  strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status == eslEOF)       esl_fatal("Empty input? No profile HMM found in %s\n",                strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status != eslOK)        esl_fatal("Unexpected error reading profile HMM from %s (code %d)\n", strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, status);

  status = esl_sqfile_OpenDigital(abc, seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open %s for reading sequences\n",             strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile);
  else if (status == eslEFORMAT)   esl_fatal("Format of sequence file %s unrecognized\n",             strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile);
  else if (status != eslOK)        esl_fatal("Unexpected error opening sequence file %s (code %d)\n", strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile, status);

  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) esl_fatal("Failed to parse sequence format of %s\n%s\n",           strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile, esl_sqfile_GetErrorBuf(sqfp));
  else if (status == eslEOF)     esl_fatal("Sequence input from %s was empty\n",                    strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile);
  else if (status != eslOK)      esl_fatal("Unexpected error reading sequence from %s (code %d)\n", strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile, status);

  h4_mode_SetLength(mo, sq->n);

  rxf = h4_refmx_Create(hmm->M, sq->n);
  rxd = h4_refmx_Create(hmm->M, sq->n);
  afu = h4_refmx_Create(hmm->M, sq->n);
  afd = h4_refmx_Create(hmm->M, sq->n);
  apu = rxf;  // reusing space to conserve memory
  apd = rxd;
  aec = afu;

  h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxf, pi, &vsc);
  h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,     &fsc);
  h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxd,     &bsc);   
  h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxd, rxd);   

  printf("# Domain structure according to Viterbi parse:\n");
  D = h4_path_GetDomainCount(pi);
  for (d = 1; d <= D; d++)
    {
      h4_path_FetchDomainBounds(pi, d, &domia, &domib, &domka, &domkb);
      printf("(     %c%4d..%4d%c     ) ", domka==1 ? '[' : '.', domia, domib, domkb==hmm->M ? ']' : '.');
    }
  printf("\n");
  for (d = 1; d <= D; d++)
    {
      h4_path_FetchDomainBounds(pi, d, &domia, &domib, &domka, &domkb);
      printf("      %c%4d..%4d%c       ", domka==1 ? '[' : '.', domka, domkb, domkb==hmm->M ? ']' : '.');
    }
  printf("\n\n");

  if ((outfile = esl_opt_GetString(go, "--ppmx")) != NULL)
    {
      if ((ofp = fopen(outfile, "w")) == NULL)
        esl_fatal("Failed to open decoding matrix dump output file %s\n", outfile);

      h4_refmx_DumpWindow(ofp, rxd,  (ia ? ia : 1), (ib ? ib : sq->n), (ka ? ka : 1), (kb ? kb : hmm->M));
      fclose(ofp);
    }

  if (( outfile = esl_opt_GetString(go, "--vpath")) != NULL)
    {
      if ((ofp = fopen(outfile, "w")) == NULL)
        esl_fatal("Failed to open viterbi path output file %s\n", outfile);
      h4_path_DumpVerbose(ofp, pi, hmm, mo, sq->dsq);
      fclose(ofp);
    }

  if (( outfile = esl_opt_GetString(go, "--diplot")) != NULL)
    {
      if ((ofp = fopen(outfile, "w")) == NULL)
        esl_fatal("Failed to open posterior decoding xy plots output file %s\n", outfile);
      h4_refmx_PlotDomainInference(ofp, rxd, pi, (ia ? ia : 1), (ib ? ib : sq->n));
      fclose(ofp);
    }

  if (( outfile = esl_opt_GetString(go, "--heatmap")) != NULL)
    {
      if ((ofp = fopen(outfile, "w")) == NULL)
        esl_fatal("Failed to open heatmap Postscript output file %s\n", outfile);
      h4_refmx_PlotHeatMap(ofp, rxd, (ia ? ia : 1), (ib ? ib : sq->n), (ka ? ka : 1), (kb ? kb : hmm->M) );
      fclose(ofp);
    }

  h4_reference_MPAS(rng, sq->dsq, sq->n, hmm, mo, rxf, rxd, pi, &wrk, ah,
                    afu, afd, anch, &asc, prm, stats);

  h4_mpas_stats_Dump(stdout, stats);

  h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch,           apu, apd, NULL);
  h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, apu, apd, apu, apd);
  h4_reference_Envelopes   (sq->dsq, sq->n, hmm, mo, anch, apu, apd, afu, afd, env);
  h4_reference_aec_Align(hmm, apu, apd, env, aec, pi);

  printf("# Domain structure according to ensemble inference:\n");
  for (d = 1; d <= env->D; d++)
    printf("(%4d %c%4d..%4d%c %4d) ", env->e[d].oa, env->e[d].ka==1 ? '[' : '.', env->e[d].ia, env->e[d].ib, env->e[d].kb==hmm->M ? ']' : '.', env->e[d].ob);
  printf("\n");
  for (d = 1; d <= env->D; d++)
    printf("      %c%4d..%4d%c       ", env->e[d].ka==1 ? '[' : '.', env->e[d].ka, env->e[d].kb, env->e[d].kb==hmm->M ? ']' : '.');
  printf("\n\n");


  if (( outfile = esl_opt_GetString(go, "--epath")) != NULL)
    {
      if ((ofp = fopen(outfile, "w")) == NULL)
        esl_fatal("Failed to open ensemble path output file %s\n", outfile);
      h4_path_DumpVerbose(ofp, pi, hmm, mo, sq->dsq);
      fclose(ofp);
    }

  esl_sqfile_Close(sqfp);
  h4_hmmfile_Close(hfp);
  
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  esl_sq_Destroy(sq);
  h4_path_Destroy(pi);
  h4_refmx_Destroy(rxf);   h4_refmx_Destroy(rxd);
  h4_refmx_Destroy(afu);   h4_refmx_Destroy(afd);
  h4_anchorset_Destroy(anch);
  h4_anchorhash_Destroy(ah);
  h4_envset_Destroy(env);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif // h4REFERENCE_AEC_EXPERIMENT2


/*****************************************************************
 * 6. Unit tests
 *****************************************************************/
#ifdef h4REFERENCE_AEC_TESTDRIVE

#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "h4_mpas.h"
#include "h4_path.h"
#include "h4_anchorset.h"
#include "h4_anchorhash.h"

#include "emit.h"
#include "modelsample.h"
#include "reference_aec.h"
#include "reference_asc.h"
#include "reference_dp.h"
#include "reference_envelopes.h"
#include "reference_mpas.h"

static void
utest_generation(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int N)
{
  char             msg[] = "reference_aec:: generation unit test failed";
  H4_PROFILE      *hmm   = NULL;
  H4_MODE         *mo    = h4_mode_Create();
  ESL_SQ          *sq    = esl_sq_CreateDigital(abc);
  H4_PATH         *gpi   = h4_path_Create();              // generated path
  H4_PATH         *pi    = h4_path_Create();              
  H4_PATH         *aecpi = h4_path_Create();              // AEC inferred path
  H4_REFMX        *rxf   = h4_refmx_Create(M, 20);        // Vit, Fwd ... then ASC Decode UP
  H4_REFMX        *rxd   = h4_refmx_Create(M, 20);        // Bck, Decode ... then ASC Decode DOWN
  H4_REFMX        *afu   = h4_refmx_Create(M, 20);        // ASC Forward UP
  H4_REFMX        *afd   = h4_refmx_Create(M, 20);        // ASC Forward DOWN
  H4_REFMX        *apu   = rxf;                           // for "clarity", we use two names for one shared matrix
  H4_REFMX        *apd   = rxd;                           //  ... this one too
  H4_REFMX        *aec   = afu;
  H4_ANCHORSET    *anch  = h4_anchorset_Create(0,0,0);
  H4_ANCHORHASH   *ah    = h4_anchorhash_Create();
  H4_ENVSET       *env   = h4_envset_Create(0,0,0); 
  float           *wrk   = NULL;                          // workspace needed (and managed) by MPAS
  float            tol   = 0.001;
  int              idx;
  float            gsc, vsc, fsc, asc, pisc;
  int              ia,ib;
  int              d;
  char             errbuf[eslERRBUFSIZE];

  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);

  for (idx = 1; idx <= N; idx++)
    {
      /* Emit sequence from model, using an arbitrary length model of <M>,
       * and restrict the emitted sequence length to 6M, arbitrarily, to 
       * keep it down to something reasonable. Then reset length model
       * to the actual length.
       */
      if ( h4_mode_SetLength(mo, M)                   != eslOK) esl_fatal(msg);
      do {
        if ( h4_emit(rng, hmm, mo, sq, gpi)           != eslOK) esl_fatal(msg);
      } while (sq->n > M * 6); 
      if ( h4_mode_SetLength(mo, sq->n)               != eslOK) esl_fatal(msg);
      if ( h4_path_Score(gpi, sq->dsq, hmm, mo, &gsc) != eslOK) esl_fatal(msg);

      /* First pass analysis */
      if ( h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxf, pi,  &vsc) != eslOK) esl_fatal(msg);   // use <rxf> matrix temporarily
      if ( h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc) != eslOK) esl_fatal(msg);   // overwrites Viterbi matrix with Forward in <rxf>
      if ( h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxd,      NULL) != eslOK) esl_fatal(msg);   // use <rxd> matrix temporarily
      if ( h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxd, rxd)  != eslOK) esl_fatal(msg);   // overwrites Backward with Decoding in <rxd>

      /* Determine most probable anchor set, followed by ASC decoding */
      if ( h4_reference_MPAS(rng, sq->dsq, sq->n, hmm, mo, rxf, rxd, pi, &wrk, ah, afu, afd, anch, &asc, NULL, NULL) != eslOK) esl_fatal(msg);
      if ( h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch, apu, apd, NULL)               != eslOK) esl_fatal(msg);
      if ( h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(msg);

      /* Envelope inference */
      if ( h4_reference_Envelopes(sq->dsq, sq->n, hmm, mo, anch, apu, apd, afu, afd, env) != eslOK) esl_fatal(msg);

      /* AEC alignment and its score */
      if ( h4_reference_aec_Align(hmm, apu, apd, env, aec, aecpi) != eslOK) esl_fatal(msg);
      if ( h4_path_Score(aecpi, sq->dsq, hmm, mo, &pisc)          != eslOK) esl_fatal(msg);

      /* Test 1. AEC path is valid */
      if ( h4_path_Validate(aecpi, hmm->M, sq->n, errbuf) != eslOK)
        esl_fatal("%s:\n  %s\n", msg, errbuf);

      /* Test 2. 
       *   AEC path score <= Viterbi; 
       *   AEC path score <= ASC <= Forward;
       *   Viterbi <= Forward;
       *   Generated path score <= Viterbi;
       *   envsc[d] <= ASC score for all d
       */
      if (pisc > vsc+tol) esl_fatal(msg);
      if (pisc > asc+tol) esl_fatal(msg);
      if (pisc > fsc+tol) esl_fatal(msg);
      if (asc  > fsc+tol) esl_fatal(msg);
      if (vsc  > fsc+tol) esl_fatal(msg);
      if (gsc  > vsc+tol) esl_fatal(msg);

      for (d = 1; d <= env->D; d++)
        if (env->e[d].env_sc > asc+tol) esl_fatal(msg);

      /* Test 3. Coords for each domain are coherent. */
      if (anch->D != env->D) esl_fatal(msg);
      for (d = 1; d <= anch->D; d++)
        {
          if (! (            1 <= env->e[d].oa &&
                  env->e[d].oa <= env->e[d].ia &&
                  env->e[d].ia <= env->e[d].i0 &&
                  env->e[d].i0 <= env->e[d].ib &&
                  env->e[d].ib <= env->e[d].ob &&
                  env->e[d].ob <= sq->n))
            esl_fatal(msg);
        }

      /* Test 4. Envelopes do not overlap (in ia/ib coords. Outer oa/ob coords can.) */
      for (d = 1; d <= env->D; d++)
        if (! (env->e[d].ia > env->e[d-1].ib)) esl_fatal(msg);  // boundary condition e[0].ib = 0

      /* Test 5. For D=1 case (only), if outer envelope encompasses
       *         generated trace's domain, env_sc > trace score.
       */
      if (h4_path_GetDomainCount(gpi) == 1 && anch->D == 1)
        {
          h4_path_FetchDomainBounds(gpi, 1, &ia, &ib, /*ka=*/NULL, /*kb=*/NULL);
          if (env->e[1].oa <= ia && env->e[1].ob >= ib)
            {
              if (! (env->e[1].env_sc >= gsc)) esl_fatal(msg);
            }
        }
    } 

  free(wrk);
  h4_envset_Destroy(env);
  h4_anchorhash_Destroy(ah);
  h4_anchorset_Destroy(anch);
  h4_refmx_Destroy(afu);  h4_refmx_Destroy(afd);    
  h4_refmx_Destroy(rxf);  h4_refmx_Destroy(rxd);    
  h4_path_Destroy(pi);    h4_path_Destroy(gpi);   h4_path_Destroy(aecpi);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
}

static void
utest_singlemulti(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int N)
{
  char msg[] = "reference_aec:: singlemulti unit test failed";
  ESL_SQ       *sq   = NULL;
  H4_PROFILE   *hmm  = NULL;
  H4_MODE      *mo   = NULL;
  H4_ANCHORSET *anch = NULL;
  H4_PATH      *gpi  = NULL;
  H4_REFMX     *afu  = h4_refmx_Create(M, 20);
  H4_REFMX     *afd  = h4_refmx_Create(M, 20);
  H4_REFMX     *apu  = h4_refmx_Create(M, 20);
  H4_REFMX     *apd  = h4_refmx_Create(M, 20);
  H4_ENVSET    *env  = h4_envset_Create(0,0,0); 
  H4_REFMX     *aec  = afu;
  H4_PATH      *api  = h4_path_Create();
  float         gsc, asc, psc;
  int           ia,ib,ka,kb;
  int           idx;
  int           d;

  for (idx = 0; idx < N; idx++)
    {
      if ( h4_modelsample_SinglePathASC(rng, abc, M, &hmm, &sq, &anch, &mo, &gpi, &gsc) != eslOK) esl_fatal(msg);

      if ( h4_reference_asc_Forward (sq->dsq, sq->n, hmm, mo, anch, afu, afd, &asc)               != eslOK) esl_fatal(msg);
      if ( h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch, apu, apd, NULL)               != eslOK) esl_fatal(msg);
      if ( h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(msg);

      if ( h4_reference_Envelopes(sq->dsq, sq->n, hmm, mo, anch, apu, apd, afu, afd, env)         != eslOK) esl_fatal(msg);

      if ( h4_reference_aec_Align(hmm, apu, apd, env, aec, api) != eslOK) esl_fatal(msg);
      if ( h4_path_Score(api, sq->dsq, hmm, mo, &psc)           != eslOK) esl_fatal(msg);

      /* Test 1. Domain #'s agree */
      if ( h4_path_GetDomainCount(gpi) != anch->D) esl_fatal(msg);
      if ( env->D                      != anch->D) esl_fatal(msg);
      if ( h4_path_GetDomainCount(api) != anch->D) esl_fatal(msg);

      /* Test 2. Inferred path is identical to the unique generated path. */
      if ( h4_path_Compare(gpi, api)   != eslOK)   esl_fatal(msg);

      /* Test 3. Envelope coords (even outer ones) match generated path */
      for (d = 1; d <= anch->D; d++)
        {
          if ( h4_path_FetchDomainBounds(gpi, d, &ia, &ib, &ka, &kb) != eslOK) esl_fatal(msg);

          if (env->e[d].oa != ia) esl_fatal(msg);
          if (env->e[d].ia != ia) esl_fatal(msg);
          if (env->e[d].ib != ib) esl_fatal(msg);
          if (env->e[d].ob != ib) esl_fatal(msg);
          if (env->e[d].ka != ka) esl_fatal(msg);
          if (env->e[d].kb != kb) esl_fatal(msg);
        }

      /* Test 4. If D==1, env score == path score */
      if (anch->D == 1 && esl_FCompare(gsc, env->e[1].env_sc, /*r_tol*/ 0.0, /*a_tol*/ 0.001) != eslOK) esl_fatal(msg);

      h4_path_Destroy(gpi);       gpi  = NULL;
      h4_anchorset_Destroy(anch); anch = NULL;
      h4_mode_Destroy(mo);        mo   = NULL;
      h4_profile_Destroy(hmm);    hmm  = NULL;
      esl_sq_Destroy(sq);         sq   = NULL;
    }

  h4_refmx_Destroy(afu); h4_refmx_Destroy(afd);
  h4_refmx_Destroy(apu); h4_refmx_Destroy(apd);
  h4_envset_Destroy(env);
  h4_path_Destroy(api);
}
  


#endif //h4REFERENCE_AEC_TESTDRIVE


/*****************************************************************
 * 7. Test driver
 *****************************************************************/
#ifdef h4REFERENCE_AEC_TESTDRIVE

#include "easel.h"
#include "esl_getopts.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                   docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show brief help summary",                   0 },
  { "-s",         eslARG_INT,     "0", NULL, NULL,  NULL,   NULL,  NULL, "set random number generator seed",          0 },
  { "-M",         eslARG_INT,    "10", NULL, NULL,  NULL,   NULL,  NULL, "set test profile length",                   0 },
  { "-N",         eslARG_INT,   "100", NULL, NULL,  NULL,   NULL,  NULL, "number of profile/seq comparisons to test", 0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show HMMER version number",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = h4_CreateDefaultApp(options, 0, argc, argv, "test driver for reference_aec", "[-options]");
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc = esl_alphabet_Create(eslAMINO);
  int             M   = esl_opt_GetInteger(go, "-M");
  int             N   = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng)); 

  utest_generation (rng, M, abc, N);
  utest_singlemulti(rng, M, abc, N);

  fprintf(stderr, "#  status   = ok\n");

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif //h4REFERENCE_AEC_TESTDRIVE


/*****************************************************************
 * 8. Example
 *****************************************************************/
#ifdef h4REFERENCE_AEC_EXAMPLE
#include <h4_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_anchorhash.h"
#include "h4_anchorset.h"
#include "h4_envset.h"
#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"
#include "h4_refmx.h"

#include "general.h"
#include "reference_aec.h"
#include "reference_asc.h"
#include "reference_dp.h"
#include "reference_envelopes.h"
#include "reference_mpas.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show HMMER version info",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of reference AEC alignment";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 2, argc, argv, banner, usage);  
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  H4_REFMX       *rxf     = NULL;
  H4_REFMX       *rxd     = NULL;
  H4_REFMX       *afu     = NULL;
  H4_REFMX       *afd     = NULL;
  H4_REFMX       *apu     = NULL;
  H4_REFMX       *apd     = NULL;
  H4_REFMX       *aec     = NULL;
  H4_PATH        *pi      = h4_path_Create();
  H4_ANCHORSET   *anch    = h4_anchorset_Create(0,0,0);
  H4_ANCHORHASH  *ah      = h4_anchorhash_Create();
  H4_ENVSET      *env     = h4_envset_Create(0,0,0);
  float          *wrk     = NULL;
  float           fsc, vsc, asc, pisc;
  char            errbuf[eslERRBUFSIZE];
  int             status;

  /* Read one profile */
  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);

  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);
 
  /* Read one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      esl_fatal("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  /* Set target length model */
  h4_mode_SetLength(mo, sq->n);

  /* Allocate four DP matrices */
  rxf = h4_refmx_Create(hmm->M, sq->n);
  rxd = h4_refmx_Create(hmm->M, sq->n);
  afu = h4_refmx_Create(hmm->M, sq->n);
  afd = h4_refmx_Create(hmm->M, sq->n);
  apu = rxf; // reusing space to conserve memory
  apd = rxd;
  aec = afu; 

  /* First pass analysis */
  h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxf, pi,  &vsc);
  h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc);
  h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxd,      NULL);   
  h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxd, rxd);   
  // now rxf = Forward matrix; rxd = Decoding

  /* Determine most probable anchor set */
  h4_reference_MPAS(rng, sq->dsq, sq->n, hmm, mo, rxf, rxd, pi, &wrk, ah,
                    afu, afd, anch, &asc, NULL, NULL);
  /* now we don't need rxf/rxd anymore. we can reuse them as apu/apd for ASC Backwards and Decoding. */

  /* ASC decoding */
  h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch, apu, apd, NULL);
  h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, apu, apd, apu, apd);

  /* Envelope determination */
  h4_reference_Envelopes(sq->dsq, sq->n, hmm, mo, anch, apu, apd, afu, afd, env);
  /* now we don't need afu/afd, and can reuse them for AEC */

  /* AEC alignment */
  h4_reference_aec_Align(hmm, apu, apd, env, aec, pi);

  if ((status = h4_path_Validate(pi, hmm->M, sq->n, errbuf)) != eslOK)
    esl_fatal("path validation failed:\n  %s\n", errbuf);
  
  h4_path_Score(pi, sq->dsq, hmm, mo, &pisc);

  printf("Viterbi sc = %.2f\n", vsc);
  printf("Forward sc = %.2f\n", fsc);
  printf("ASC sc     = %.2f\n", asc);
  printf("AEC ali sc = %.2f\n", pisc);

  h4_anchorset_Dump(stdout, anch);
  h4_envset_Dump(stdout, env);
  h4_path_Dump(stdout, pi);

  free(wrk);
  h4_refmx_Destroy(rxf);   h4_refmx_Destroy(rxd);
  h4_refmx_Destroy(afu);   h4_refmx_Destroy(afd);
  h4_anchorset_Destroy(anch);
  h4_anchorhash_Destroy(ah);
  h4_envset_Destroy(env);
  esl_sqfile_Close(sqfp);
  h4_hmmfile_Close(hfp);
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  h4_path_Destroy(pi);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif //REFERENCE_AEC_EXAMPLE
