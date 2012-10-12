extern int p7_profile_Config         (P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg);
extern int p7_profile_ConfigLocal    (P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L);
extern int p7_profile_ConfigUnilocal (P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L);
extern int p7_profile_ConfigGlocal   (P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L);
extern int p7_profile_ConfigUniglocal(P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L);
extern int p7_profile_ConfigCustom   (P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L, float nj, float pglocal);
extern int p7_profile_SetLength      (P7_PROFILE *gm, int L);

