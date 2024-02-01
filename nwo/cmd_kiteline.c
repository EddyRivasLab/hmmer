#include <h4_config.h>

#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_cpu.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_subcmd.h"

#include "h4_anchorhash.h"
#include "h4_anchorset.h"
#include "h4_checkptmx.h"
#include "h4_envset.h"
#include "h4_filtermx.h"
#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_mpas.h"
#include "h4_profile.h"
#include "h4_refmx.h"
#include "h4_sparsemask.h"

#include "general.h"
#include "calibrate.h"
#include "fbfilter.h"
#include "ssvfilter.h"
#include "vitfilter.h"
#include "sparse_dp.h"
#include "zigar.h"

#include "reference_dp.h"
#include "reference_asc.h"
#include "reference_mpas.h"
#include "reference_envelopes.h"
#include "reference_aec.h"

static ESL_OPTIONS kiteline_options[] = {
   /* name             type      default                         env  range      toggles reqs incomp  help                                                docgroup*/
  { "-h",           eslARG_NONE,  FALSE,                        NULL, NULL,      NULL, NULL, NULL,     "show brief help on version and usage",                     1 },
  { "--seed",       eslARG_INT,     "0",                        NULL, "n>=0",    NULL, NULL, NULL,     "set random number seed to <n>",                            1 },
};


int
h4_cmd_kiteline(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_subcmd_CreateDefaultApp(topcmd, sub, kiteline_options, argc, argv);
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
  H4_FILTERMX    *fx      = NULL;
  H4_CHECKPTMX   *cpx     = NULL;
  // H4_SPARSEMASK  *sm      = NULL;
  // H4_SPARSEMX    *sxf     = h4_sparsemx_Create(NULL);
  // H4_SPARSEMX    *sxd     = h4_sparsemx_Create(NULL);
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
  double          lambda, sfmu, vfmu, fftau;
  float           sfsc, vfsc, ffsc, vsc, fsc, bsc, asc;
  double          sfP,  vfP,  ffP,  fwdP;
  int             d, z;
  char           *zali    = NULL;
  int             status;

  h4_Init();

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

  /* Calibrate the model.
   * Eventually, this will already be done in hmmbuild and stored in the model.
   */
  h4_lambda (     hmm, &lambda);
  h4_ssv_mu (rng, hmm, /*L=*/200, /*N=*/200, lambda,       &sfmu);
  h4_vit_mu (rng, hmm, /*L=*/200, /*N=*/200, lambda,       &vfmu);
  h4_fwd_tau(rng, hmm, /*L=*/100, /*N=*/200, lambda, 0.04, &fftau);

  fx  = h4_filtermx_Create(hmm->M);
  cpx = h4_checkptmx_Create(hmm->M, 400, ESL_MBYTES(32));
  //  sm  = h4_sparsemask_Create(hmm->M, 400);
  rxf = h4_refmx_Create(hmm->M, 400);
  rxd = h4_refmx_Create(hmm->M, 400);
  afu = h4_refmx_Create(hmm->M, 400);
  afd = h4_refmx_Create(hmm->M, 400);
  apu = rxf;  // reusing space to conserve memory
  apd = rxd;
  aec = afu;

  while ((status = esl_sqio_Read(sqfp, sq)) != eslEOF)
    {
      if      (status == eslEFORMAT) esl_fatal("Error: failed to parse sequence format of %s\n%s\n",    strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile, esl_sqfile_GetErrorBuf(sqfp));
      else if (status != eslOK)      esl_fatal("Unexpected error reading sequence from %s (code %d)\n", strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile, status);

      h4_mode_SetLength(mo, sq->n);

      h4_ssvfilter(sq->dsq, sq->n, hmm, mo,      &sfsc);
      sfP  = esl_gumbel_surv(sfsc, sfmu,  lambda);
      if (sfP > 0.02) goto NO_HIT;

      h4_vitfilter(sq->dsq, sq->n, hmm, mo, fx,  &vfsc);
      vfP  = esl_gumbel_surv(vfsc, vfmu,  lambda);
      if (vfP > 0.001) goto NO_HIT;

      h4_fwdfilter(sq->dsq, sq->n, hmm, mo, cpx, &ffsc);
      ffP  = esl_exp_surv   (ffsc, fftau, lambda);
      if (ffP > 1e-5) goto NO_HIT;

      // h4_bckfilter(sq->dsq, sq->n, hmm, mo, cpx, sm, h4SPARSIFY_THRESH);

      // h4_sparse_Viterbi (sq->dsq, sq->n, hmm, mo, sm, sxf, pi, &vsc);
      // h4_sparse_Forward (sq->dsq, sq->n, hmm, mo, sm, sxf,     &fsc);
      // h4_sparse_Backward(sq->dsq, sq->n, hmm, mo, sm, sxd,     &bsc);
      // h4_sparse_Decoding(sq->dsq, sq->n, hmm, mo, sxf, sxd, sxd);

      /* First pass analysis */
      h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxf, pi, &vsc);
      h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,     &fsc);
      h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxd,     &bsc);   
      h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxd, rxd);   

      fwdP  = esl_exp_surv(fsc, fftau-1.0, lambda);  // the -1.0 is for pG=0.5

      /* MPAS to identify domains and their anchors */
      h4_reference_MPAS(rng, sq->dsq, sq->n, hmm, mo, rxf, rxd, pi, &wrk, ah,
                        afu, afd, anch, &asc, NULL, NULL);

      /* ASC decoding */
      h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch, apu, apd, NULL);
      h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, apu, apd, apu, apd);

      /* Envelope determination */
      h4_reference_Envelopes(sq->dsq, sq->n, hmm, mo, anch, apu, apd, afu, afd, env);

      /* AEC alignment */
      h4_reference_aec_Align(hmm, apu, apd, env, aec, pi);

      //printf("%-30s %10.2f %10.2g %10.2f %10.2g %10.2f %10.2g %10.2f %10.2f %10.2f %10.2g %4d\n", sq->name, sfsc, sfP, vfsc, vfP, ffsc, ffP, vsc, fsc, bsc, fwdP, anch->D);

      z = 0;
      for (d = 1; d <= env->D; d++)
        {
          while (pi->st[z] != h4P_L && pi->st[z] != h4P_G) z++;
          h4_zigar_Encode(pi, z, &zali);
          z++;
          
          printf("%-30s %10.2f %10.2g %4d %4d %6d %6d %6d %6d %6d %6d %6d %6d %10.2f %s\n",
                 sq->name,
                 fsc, fwdP,
                 d, env->D,
                 env->e[d].ia, env->e[d].ib, env->e[d].oa, env->e[d].ob, (int) sq->n,
                 env->e[d].ka, env->e[d].kb, hmm->M,
                 env->e[d].env_sc,
                 zali);

          free(zali);
        }

    NO_HIT:
      //h4_sparsemask_Reuse(sm);
      esl_sq_Reuse(sq);
    }
         
  esl_sqfile_Close(sqfp);
  h4_hmmfile_Close(hfp);

  h4_anchorhash_Destroy(ah);
  h4_anchorset_Destroy(anch);
  h4_refmx_Destroy(rxf);
  h4_refmx_Destroy(rxd);
  h4_refmx_Destroy(afu);
  h4_refmx_Destroy(afd);  // ... apu, apd, aec point to rxf, rxd, afu and don't need to be free'd
  h4_path_Destroy(pi);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  h4_filtermx_Destroy(fx);
  h4_checkptmx_Destroy(cpx);
  //h4_sparsemask_Destroy(sm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}

