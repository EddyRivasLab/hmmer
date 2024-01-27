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

#include "h4_checkptmx.h"
#include "h4_filtermx.h"
#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_sparsemask.h"

#include "general.h"
#include "calibrate.h"
#include "fbfilter.h"
#include "ssvfilter.h"
#include "vitfilter.h"
#include "sparse_dp.h"


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
  H4_SPARSEMASK  *sm      = NULL;
  H4_SPARSEMX    *sxf     = h4_sparsemx_Create(NULL);
  H4_SPARSEMX    *sxd     = h4_sparsemx_Create(NULL);
  H4_PATH        *vpi     = h4_path_Create();
  double          lambda, sfmu, vfmu, fftau;
  float           sfsc, vfsc, ffsc, vsc, fsc, bsc;
  double          sfP,  vfP,  ffP;
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
  sm  = h4_sparsemask_Create(hmm->M, 400);

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

      h4_bckfilter(sq->dsq, sq->n, hmm, mo, cpx, sm, h4SPARSIFY_THRESH);

      h4_sparse_Viterbi (sq->dsq, sq->n, hmm, mo, sm, sxf, vpi, &vsc);
      h4_sparse_Forward (sq->dsq, sq->n, hmm, mo, sm, sxf,      &fsc);
      h4_sparse_Backward(sq->dsq, sq->n, hmm, mo, sm, sxd,      &bsc);
      h4_sparse_Decoding(sq->dsq, sq->n, hmm, mo, sxf, sxd, sxd);

      printf("%-30s %10.2f %10.2g %10.2f %10.2g %10.2f %10.2g %10.2f %10.2f %10.2f\n", sq->name, sfsc, sfP, vfsc, vfP, ffsc, ffP, vsc, fsc, bsc);

    NO_HIT:
      h4_sparsemask_Reuse(sm);
      esl_sq_Reuse(sq);
    }
         

  esl_randomness_Destroy(rng);
  esl_sqfile_Close(sqfp);
  h4_hmmfile_Close(hfp);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  h4_filtermx_Destroy(fx);
  h4_checkptmx_Destroy(cpx);
  h4_sparsemask_Destroy(sm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}

