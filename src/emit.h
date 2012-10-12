extern int p7_CoreEmit   (ESL_RANDOMNESS *r, const P7_HMM *hmm,                                        ESL_SQ *sq, P7_TRACE *tr);
extern int p7_ProfileEmit(ESL_RANDOMNESS *r, const P7_HMM *hmm, const P7_PROFILE *gm, const P7_BG *bg, ESL_SQ *sq, P7_TRACE *tr);
extern int p7_emit_SimpleConsensus(const P7_HMM *hmm, ESL_SQ *sq);
extern int p7_emit_FancyConsensus (const P7_HMM *hmm, float min_lower, float min_upper, ESL_SQ *sq);
