
extern int p7_Calibrate(P7_HMM *hmm, P7_BUILDER *cfg_b, ESL_RANDOMNESS **byp_rng, P7_BG **byp_bg, P7_PROFILE **byp_gm, P7_OPROFILE **byp_om);
extern int p7_Lambda(P7_HMM *hmm, P7_BG *bg, double *ret_lambda);
extern int p7_MSVMu     (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda,               double *ret_mmu);
extern int p7_ViterbiMu (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda,               double *ret_vmu);
extern int p7_Tau       (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda, double tailp, double *ret_tau);
