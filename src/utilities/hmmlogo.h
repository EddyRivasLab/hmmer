extern float hmmlogo_maxHeight (P7_BG *bg);
extern int hmmlogo_emissionHeightsDivRelent (P7_HMM *hmm, P7_BG *bg, float *rel_ents, float **heights );
extern int hmmlogo_posScoreHeightsDivRelent (P7_HMM *hmm, P7_BG *bg, float *rel_ents, float **heights );
extern int hmmlogo_ScoreHeights (P7_HMM *hmm, P7_BG *bg, float **heights );
extern int hmmlogo_IndelValues (P7_HMM *hmm, float *insert_P, float *insert_expL, float *delete_P );
